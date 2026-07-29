#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: combine_final_tables.R <input_dir> <output_dir> [max_gap_bp] [bed_panel]")
}

input_dir  <- args[1]
output_dir <- args[2]
max_gap_bp <- if (length(args) >= 3) as.integer(args[3]) else 500000L
bed_panel  <- if (length(args) >= 4 && nchar(args[4]) > 0) args[4] else NULL
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# --- Annotated files (*_final_CNVs_annotated.tsv) ---
annotated_files <- list.files(input_dir, pattern = "_final_CNVs_annotated\\.tsv$",
                              full.names = TRUE, recursive = FALSE)
if (length(annotated_files) == 0) {
  cat("WARNING: no *_final_CNVs_annotated.tsv files found in", input_dir, "\n")
  annotated_all <- data.table()
} else {
  cat(sprintf("Annotated files found: %d\n", length(annotated_files)))
  annotated_list <- lapply(annotated_files, function(f) {
    dt <- tryCatch(fread(f, sep = "\t", header = TRUE),
                   error = function(e) {
                     cat("  SKIP:", basename(f), "-", conditionMessage(e), "\n")
                     NULL
                   })
    if (!is.null(dt) && nrow(dt) > 0) dt[, source_sample := basename(f)]
    dt
  })
  annotated_all <- rbindlist(Filter(Negate(is.null), annotated_list), fill = TRUE)
  cat(sprintf("  -> %d rows after binding\n", nrow(annotated_all)))
}

# --- Merged/smoothed files (*_smoothed_variants.tsv) ---
merged_files <- list.files(input_dir, pattern = "_smoothed_variants\\.tsv$",
                           full.names = TRUE, recursive = FALSE)
if (length(merged_files) == 0) {
  cat("WARNING: no *_smoothed_variants.tsv files found in", input_dir, "\n")
  merged_all <- data.table()
} else {
  cat(sprintf("Merged files found: %d\n", length(merged_files)))
  merged_list <- lapply(merged_files, function(f) {
    tryCatch(fread(f, sep = "\t", header = TRUE),
             error = function(e) {
               cat("  SKIP:", basename(f), "-", conditionMessage(e), "\n")
               NULL
             })
  })
  merged_all <- rbindlist(Filter(Negate(is.null), merged_list), fill = TRUE)
  cat(sprintf("  -> %d rows after binding\n", nrow(merged_all)))
}

# --- Normalise column names before join ---
if (nrow(annotated_all) > 0) {
  if ("CHR"  %in% names(annotated_all)) setnames(annotated_all, "CHR",  "CHROM")
  if ("STOP" %in% names(annotated_all)) setnames(annotated_all, "STOP", "END")
  if ("source_sample" %in% names(annotated_all)) {
    annotated_all[, sample := sub("_final_CNVs_annotated\\.tsv$", "", source_sample)]
    annotated_all[, source_sample := NULL]
  }
}

# --- Write per-type bound tables ---
fwrite(merged_all,    file.path(output_dir, "all_samples_merged.tsv"),    sep = "\t")
fwrite(annotated_all, file.path(output_dir, "all_samples_annotated.tsv"), sep = "\t")
cat("Wrote all_samples_merged.tsv and all_samples_annotated.tsv\n")

# --- Merge by sample + CHROM, START, END ---
if (nrow(merged_all) == 0 || nrow(annotated_all) == 0) {
  cat("WARNING: one or both tables empty — combined output will equal the non-empty table\n")
  combined <- if (nrow(merged_all) > 0) merged_all else annotated_all
} else {
  combined <- merge(merged_all, annotated_all,
                    by = c("sample", "CHROM", "START", "END"),
                    all.x = TRUE)
  cat(sprintf("Combined: %d rows x %d cols\n", nrow(combined), ncol(combined)))
}

fwrite(combined, file.path(output_dir, "all_samples_combined.tsv"), sep = "\t")
cat("Wrote all_samples_combined.tsv\n")

# ---------------------------------------------------------------------------
# Per-exon breakdown for consensus-path regions
# Consensus candidates use UNION coordinates (large when gene-level callers are
# involved).  Re-intersecting with the capture BED gives per-exon entries that
# match the granularity produced by merging_and_smoothing.R for the voting path,
# so smooth_combined() behaves identically for both paths.
# Only applied when 'consensus_score' column is present (consensus path only).
# ---------------------------------------------------------------------------
if (!is.null(bed_panel) && file.exists(bed_panel) &&
    nrow(combined) > 0 &&
    "consensus_score" %in% names(combined)) {

  is_cons <- !is.na(combined[["consensus_score"]])

  if (any(is_cons)) {
    cat(sprintf("Per-exon breakdown: %d consensus regions using BED panel %s\n",
                sum(is_cons), basename(bed_panel)))

    bed <- tryCatch(
      fread(bed_panel, header = FALSE, sep = "\t"),
      error = function(e) { cat("  WARNING: could not read BED panel:", conditionMessage(e), "\n"); NULL }
    )

    if (!is.null(bed) && ncol(bed) >= 3) {
      bed_col_names <- c("b_CHROM", "b_START", "b_END",
                         if (ncol(bed) >= 4) "b_NAME" else NULL)
      setnames(bed, seq_along(bed_col_names), bed_col_names)
      bed[, b_CHROM := gsub("^chr", "", as.character(b_CHROM))]
      bed[, b_START := as.integer(b_START) + 1L]   # BED 0-based → 1-based
      bed[, b_END   := as.integer(b_END)]
      setkey(bed, b_CHROM, b_START, b_END)

      cons_dt  <- combined[ is_cons]
      other_dt <- combined[!is_cons]

      cons_dt[, CHROM := as.character(CHROM)]
      setkey(cons_dt, CHROM, START, END)

      ov <- tryCatch(
        foverlaps(cons_dt, bed,
                  by.x = c("CHROM", "START", "END"),
                  by.y = c("b_CHROM", "b_START", "b_END"),
                  type = "any", nomatch = NA),
        error = function(e) { cat("  WARNING: foverlaps failed:", conditionMessage(e), "\n"); NULL }
      )

      if (!is.null(ov)) {
        matched   <- ov[!is.na(b_START)]
        unmatched <- ov[ is.na(b_START)]

        if (nrow(matched) > 0) {
          matched[, START := b_START]
          matched[, END   := b_END]
          cat(sprintf("  -> %d per-exon entries (from %d consensus regions)\n",
                      nrow(matched), sum(is_cons)))
        }

        drop_cols <- intersect(c("b_CHROM","b_START","b_END","b_NAME"), names(ov))
        for (col in drop_cols) {
          if (col %in% names(matched))   matched[,   (col) := NULL]
          if (col %in% names(unmatched)) unmatched[, (col) := NULL]
        }

        combined <- rbindlist(list(other_dt, matched, unmatched), fill = TRUE)
        cat(sprintf("  -> %d total entries after breakdown\n", nrow(combined)))
      }
    }
  }
}

# ---------------------------------------------------------------------------
# Smoothing: merge adjacent segments per sample, same type, gap <= max_gap_bp
# ---------------------------------------------------------------------------
smooth_combined <- function(dt, max_gap_bp = 500000L) {
  if (nrow(dt) == 0) return(data.table())

  na_str <- function(x) {
    if (is.null(x) || length(x) == 0 || (length(x) == 1 && is.na(x))) "" else as.character(x)
  }

  split_field <- function(x, sep) {
    s <- na_str(x)
    if (s == "") character(0) else trimws(strsplit(s, sep, fixed = TRUE)[[1]])
  }

  dt <- dt[order(sample, CHROM, START)]

  segments <- list()
  cur <- NULL

  for (i in seq_len(nrow(dt))) {
    r <- dt[i]

    same_seg <- !is.null(cur) &&
      identical(cur$sample, r[["sample"]]) &&
      identical(cur$CHROM,  r[["CHROM"]])  &&
      identical(cur$type,   na_str(r[["consensus_type"]])) &&
      !is.na(r[["START"]]) &&
      (r[["START"]] - cur$END) <= max_gap_bp

    if (same_seg) {
      cur$END       <- max(cur$END, r[["END"]])
      cur$n_targets <- cur$n_targets + 1L
      cur$n_callers <- max(cur$n_callers, r[["n_callers"]], na.rm = TRUE)

      cur$callers <- paste(sort(unique(c(
        split_field(cur$callers, ";"),
        split_field(r[["callers"]], ";")
      ))), collapse = ";")

      cur$genes <- paste(unique(c(
        split_field(cur$genes, ", "),
        split_field(r[["GFT_genes"]], ", ")
      )), collapse = ", ")

      lbl <- na_str(r[["cn_label"]])
      if (lbl != "")
        cur$cn_labels <- if (cur$cn_labels == "") lbl
                         else paste(cur$cn_labels, lbl, sep = "|")

      cur$classifications <- paste(sort(unique(c(
        split_field(cur$classifications, ";"),
        split_field(r[["Classification"]], ";")
      ))), collapse = ";")

      sc <- suppressWarnings(as.numeric(na_str(r[["Total score"]])))
      if (length(sc) == 1L && !is.na(sc))
        cur$max_total_score <- max(cur$max_total_score, sc, na.rm = TRUE)

      cur$target_names <- paste(unique(c(
        split_field(cur$target_names, ";"),
        split_field(r[["target_names"]], ";")
      )), collapse = ";")

      cur$BED_gene_name <- paste(unique(c(
        split_field(cur$BED_gene_name, ";"),
        split_field(r[["BED_gene_name"]], ";")
      )), collapse = ";")

      cur$gene_biotype <- paste(unique(c(
        split_field(cur$gene_biotype, ";"),
        split_field(r[["gene_biotype"]], ";")
      )), collapse = ";")

      cur$dosage_sensitive_genes <- paste(unique(c(
        split_field(cur$dosage_sensitive_genes, ", "),
        split_field(r[["Known or predicted dosage-sensitive genes"]], ", ")
      )), collapse = ", ")

    } else {
      if (!is.null(cur)) segments <- c(segments, list(cur))
      cur <- list(
        sample                  = r[["sample"]],
        CHROM                   = r[["CHROM"]],
        START                   = r[["START"]],
        END                     = r[["END"]],
        type                    = na_str(r[["consensus_type"]]),
        n_targets               = 1L,
        n_callers               = r[["n_callers"]],
        callers                 = na_str(r[["callers"]]),
        genes                   = na_str(r[["GFT_genes"]]),
        cn_labels               = na_str(r[["cn_label"]]),
        classifications         = na_str(r[["Classification"]]),
        max_total_score         = suppressWarnings(as.numeric(na_str(r[["Total score"]]))),
        target_names            = na_str(r[["target_names"]]),
        BED_gene_name           = na_str(r[["BED_gene_name"]]),
        gene_biotype            = na_str(r[["gene_biotype"]]),
        dosage_sensitive_genes  = na_str(r[["Known or predicted dosage-sensitive genes"]])
      )
    }
  }
  if (!is.null(cur)) segments <- c(segments, list(cur))

  if (length(segments) == 0) return(data.table())
  rbindlist(lapply(segments, as.data.table), fill = TRUE)
}

if (nrow(combined) > 0) {
  cat(sprintf("Smoothing combined data (max_gap_bp = %d bp)...\n", max_gap_bp))
  smoothed <- smooth_combined(combined, max_gap_bp = max_gap_bp)
  cat(sprintf("  -> %d smoothed segments\n", nrow(smoothed)))
  fwrite(smoothed, file.path(output_dir, "all_samples_smoothed.tsv"), sep = "\t")
  cat("Wrote all_samples_smoothed.tsv\n")
} else {
  smoothed <- data.table()
  fwrite(smoothed, file.path(output_dir, "all_samples_smoothed.tsv"), sep = "\t")
}

# ---------------------------------------------------------------------------
# HTML final report
# ---------------------------------------------------------------------------

# Render a base-R plot expression to a PNG and return an <img> tag with the
# image embedded as a base64 data URI. PNG bypasses all SVG font-metric issues.
plot_to_img <- function(expr, width = 8, height = 5, res = 96) {
  tmp <- tempfile(fileext = ".png")
  on.exit(unlink(tmp), add = TRUE)
  grDevices::png(tmp, width = round(width * res), height = round(height * res),
                 res = res, bg = "white")
  tryCatch(
    force(expr),
    error = function(e) {
      graphics::plot.new()
      graphics::text(0.5, 0.5, paste("Plot error:", conditionMessage(e)),
                     cex = 0.85, col = "#dc2626")
    }
  )
  grDevices::dev.off()
  b64 <- base64enc::base64encode(tmp)
  paste0('<img src="data:image/png;base64,', b64,
         '" style="width:100%;height:auto;display:block">')
}

generate_html_report <- function(dt, output_path) {
  esc <- function(x) {
    x <- ifelse(is.na(x), "", as.character(x))
    x <- gsub("&", "&amp;", x, fixed = TRUE)
    x <- gsub("<", "&lt;",  x, fixed = TRUE)
    x <- gsub(">", "&gt;",  x, fixed = TRUE)
    x <- gsub('"', "&quot;", x, fixed = TRUE)
    x
  }
  fmt_bp <- function(x) {
    x <- suppressWarnings(as.numeric(x))
    ifelse(is.na(x), "",
      ifelse(x >= 1e6, paste0(round(x / 1e6, 2), " Mb"),
      ifelse(x >= 1e3, paste0(round(x / 1e3,  1), " kb"),
             paste0(x, " bp"))))
  }

  n_samples  <- if (nrow(dt) > 0 && "sample" %in% names(dt)) length(unique(dt$sample)) else 0L
  n_variants <- nrow(dt)

  # --- type breakdown cards ---
  type_html <- ""
  if (n_variants > 0 && "type" %in% names(dt)) {
    tc <- sort(table(dt[["type"]]), decreasing = TRUE)
    for (nm in names(tc)) {
      cls <- if (grepl("(?i)del", nm, perl = TRUE)) "bdel"
             else if (grepl("(?i)dup", nm, perl = TRUE)) "bdup"
             else "both"
      type_html <- paste0(type_html,
        '<div class="si"><span class="badge ', cls, '">', esc(nm), '</span>',
        '<span class="sn">', tc[[nm]], '</span></div>\n')
    }
  }

  # --- top 10 genes ---
  top_genes_html <- "<em>No gene data</em>"
  if (n_variants > 0 && "genes" %in% names(dt)) {
    gv <- unlist(strsplit(paste(na.omit(dt$genes), collapse = ","), "[,;]+"))
    gv <- trimws(gv[nchar(trimws(gv)) > 0 & trimws(gv) != "NA"])
    if (length(gv) > 0) {
      gt <- sort(table(gv), decreasing = TRUE)
      top10 <- head(gt, 10)
      rows <- paste(sapply(seq_along(top10), function(i)
        sprintf("<tr><td>%s</td><td>%d</td></tr>", esc(names(top10)[i]), top10[[i]])),
        collapse = "\n")
      top_genes_html <- paste0(
        '<table class="mt"><thead><tr><th>Gene</th><th>Count</th></tr></thead><tbody>',
        rows, '</tbody></table>')
    }
  }

  # --- top 10 dosage-sensitive genes ---
  dosage_html <- "<em>No data</em>"
  if (n_variants > 0 && "dosage_sensitive_genes" %in% names(dt)) {
    dv <- unlist(strsplit(paste(na.omit(dt$dosage_sensitive_genes), collapse = ", "), ",\\s*"))
    dv <- trimws(dv[nchar(trimws(dv)) > 0 & trimws(dv) != "NA"])
    if (length(dv) > 0) {
      dt2 <- sort(table(dv), decreasing = TRUE)
      top10 <- head(dt2, 10)
      rows <- paste(sapply(seq_along(top10), function(i)
        sprintf("<tr><td>%s</td><td>%d</td></tr>", esc(names(top10)[i]), top10[[i]])),
        collapse = "\n")
      dosage_html <- paste0(
        '<table class="mt"><thead><tr><th>Gene</th><th>Count</th></tr></thead><tbody>',
        rows, '</tbody></table>')
    }
  }

  # ---------------------------------------------------------------------------
  # PNG plots (embedded inline as base64 data URIs)
  # ---------------------------------------------------------------------------
  plots_html <- ""
  if (n_variants > 0) {

    CX_MAIN <- 1.3
    CX_LAB  <- 1.15
    CX_AXIS <- 1.05

    STD_W   <- 8
    STD_H   <- 6
    STD_RES <- 130

    # Plot 1 – full-width vertical bar: samples on x-axis, counts on y-axis
    p1_h   <- min(10, max(5, n_samples * 0.18 + 2))
    p1_cx  <- max(0.45, min(0.85, 14 / n_samples * 1.8))
    p1_bot <- max(6, max(nchar(sort(unique(dt$sample)))) * 0.52)
    svg1 <- plot_to_img({
      types_u   <- sort(unique(dt$type))
      samples_u <- sort(unique(dt$sample))
      mat <- matrix(0L, nrow = length(types_u), ncol = length(samples_u),
                    dimnames = list(types_u, samples_u))
      by_st <- dt[, .N, by = .(sample, type)]
      for (k in seq_len(nrow(by_st)))
        mat[by_st$type[k], by_st$sample[k]] <- by_st$N[k]
      bar_cols <- sapply(rownames(mat), function(t)
        if (grepl("(?i)del", t, perl = TRUE)) "#dc2626"
        else if (grepl("(?i)dup", t, perl = TRUE)) "#2563eb"
        else "#9ca3af")
      par(mar = c(p1_bot, 5.5, 3.5, 1))
      barplot(mat, horiz = FALSE, las = 2, col = bar_cols, border = NA,
              ylab = "Number of variants", main = "Variants per sample",
              cex.names = p1_cx, cex.axis = CX_AXIS, cex.lab = CX_LAB,
              cex.main = CX_MAIN,
              legend.text = rownames(mat),
              args.legend = list(x = "topright", bty = "n", cex = 1.0))
    }, width = 14, height = p1_h, res = 150)

    # Plot 2 – chromosome distribution: vertical bars, chromosomes on x-axis
    svg2 <- plot_to_img({
      chr_order <- c(paste0("chr", c(as.character(1:22), "X", "Y", "M")),
                     c(as.character(1:22), "X", "Y", "MT"))
      ch_counts <- table(dt$CHROM)
      ord  <- chr_order[chr_order %in% names(ch_counts)]
      rest <- sort(setdiff(names(ch_counts), ord))
      ch_counts <- ch_counts[c(ord, rest)]
      alt_col <- rep(c("#2563eb", "#3b82f6"), length.out = length(ch_counts))
      bot <- max(5, max(nchar(names(ch_counts))) * 0.52)
      par(mar = c(bot, 5.5, 3.5, 1))
      barplot(ch_counts, las = 2, col = alt_col, border = NA,
              ylab = "Number of variants",
              main = "Variants per chromosome",
              cex.names = 1.0, cex.axis = CX_AXIS,
              cex.lab = CX_LAB, cex.main = CX_MAIN)
    }, width = STD_W, height = STD_H, res = STD_RES)

    # Plot 3 – classification consequence pie chart
    svg3 <- plot_to_img({
      if ("classifications" %in% names(dt)) {
        cv <- unlist(strsplit(paste(na.omit(dt$classifications), collapse = ";"), ";"))
        cv <- trimws(cv[nchar(trimws(cv)) > 0 & trimws(cv) != "NA"])
      } else cv <- character(0)
      if (length(cv) > 0) {
        ct  <- sort(table(cv), decreasing = TRUE)
        pal <- c("#dc2626", "#f97316", "#eab308", "#22c55e", "#2563eb",
                 "#7c3aed", "#0891b2", "#db2777", "#65a30d", "#6366f1")
        cols <- pal[(seq_along(ct) - 1L) %% length(pal) + 1L]
        pct  <- round(100 * ct / sum(ct), 1)
        lbls <- paste0(names(ct), "\n", pct, "%")
        par(mar = c(2, 2, 3.5, 2))
        pie(ct, labels = lbls, col = cols, border = "white",
            main = "Consequence distribution",
            cex = 0.95, cex.main = CX_MAIN)
      } else {
        plot.new()
        text(0.5, 0.5, "No classification data", cex = 1.3, col = "grey60")
      }
    }, width = STD_W, height = STD_H, res = STD_RES)

    # Plot 4 – top 15 affected genes (horizontal bar)
    svg4 <- plot_to_img({
      gv <- unlist(strsplit(paste(na.omit(dt$genes), collapse = ","), "[,;]+"))
      gv <- trimws(gv[nchar(trimws(gv)) > 0 & trimws(gv) != "NA"])
      if (length(gv) > 0) {
        top15    <- tail(sort(table(gv)), 15)
        left_mar <- max(8, max(nchar(names(top15))) * 0.72)
        par(mar = c(5, left_mar, 3.5, 2))
        barplot(top15, horiz = TRUE, las = 1, col = "#7c3aed", border = NA,
                xlab = "Number of variants",
                main = "Top 15 affected genes",
                cex.names = 0.95, cex.axis = CX_AXIS,
                cex.lab = CX_LAB, cex.main = CX_MAIN)
      } else {
        plot.new()
        text(0.5, 0.5, "No gene data", cex = 1.3, col = "grey60")
      }
    }, width = STD_W, height = STD_H, res = STD_RES)

    # Plot 5 – caller contributions (horizontal bar)
    svg5 <- plot_to_img({
      cv <- unlist(strsplit(paste(na.omit(dt$callers), collapse = ";"), ";"))
      cv <- trimws(cv[nchar(trimws(cv)) > 0 & trimws(cv) != "NA"])
      if (length(cv) > 0) {
        ct  <- sort(table(cv), decreasing = TRUE)
        pal <- c("#2563eb","#dc2626","#059669","#d97706","#7c3aed",
                 "#0891b2","#db2777","#65a30d","#ea580c","#4f46e5")
        cols     <- pal[(seq_along(ct) - 1L) %% length(pal) + 1L]
        left_mar <- max(8, max(nchar(names(ct))) * 0.72)
        par(mar = c(5, left_mar, 3.5, 2))
        barplot(ct, horiz = TRUE, las = 1, col = cols, border = NA,
                xlab = "Number of variants",
                main = "Callers contributing to variants",
                cex.names = 0.95, cex.axis = CX_AXIS,
                cex.lab = CX_LAB, cex.main = CX_MAIN)
      } else {
        plot.new()
        text(0.5, 0.5, "No caller data", cex = 1.3, col = "grey60")
      }
    }, width = STD_W, height = STD_H, res = STD_RES)

    plots_html <- paste0(
      '<div class="plots-grid">\n',
      '<div class="plot-card full"><h3>Variants per sample</h3>',       svg1, '</div>\n',
      '<div class="plot-card"><h3>Chromosome distribution</h3>',        svg2, '</div>\n',
      '<div class="plot-card"><h3>Consequence distribution</h3>',       svg3, '</div>\n',
      '<div class="plot-card"><h3>Top 15 affected genes</h3>',          svg4, '</div>\n',
      '<div class="plot-card"><h3>Caller contributions</h3>',           svg5, '</div>\n',
      '</div>\n'
    )
  }

  # --- main table ---
  view_cols    <- c("sample", "CHROM", "START", "END", "type",
                    "genes", "cn_labels", "classifications",
                    "target_names", "dosage_sensitive_genes")
  view_headers <- c("Sample", "Chr", "Start", "End", "Type",
                    "Genes", "CN Labels", "Classification",
                    "Target names", "Dosage-sensitive genes")
  header_html <- paste(sapply(seq_along(view_headers), function(i)
    paste0('<th onclick="srt(this)">', view_headers[i],
           ' <span class="si2">&#x21C5;</span></th>')),
    collapse = "")

  data_rows_html <- ""
  if (n_variants > 0) {
    data_rows_html <- paste(sapply(seq_len(nrow(dt)), function(i) {
      r  <- dt[i]
      tp <- tolower(trimws(as.character(r[["type"]])))
      rc <- if (grepl("del", tp)) ' class="rdel"' else if (grepl("dup", tp)) ' class="rdup"' else ""
      sz <- ""
      if (!is.na(r[["START"]]) && !is.na(r[["END"]]))
        sz <- fmt_bp(as.numeric(r[["END"]]) - as.numeric(r[["START"]]))
      cells <- paste(sapply(view_cols, function(col) {
        val <- if (col %in% names(r)) esc(r[[col]]) else ""
        if (col == "END" && sz != "")
          val <- paste0(val, ' <small>(', sz, ')</small>')
        paste0("<td>", val, "</td>")
      }), collapse = "")
      paste0("<tr", rc, ">", cells, "</tr>")
    }), collapse = "\n")
  }

  rep_date <- format(Sys.time(), "%Y-%m-%d %H:%M")

  css <- '
<style>
*,*::before,*::after{box-sizing:border-box;margin:0;padding:0}
:root{--bg:#f5f7fa;--card:#fff;--bd:#dde1e7;--tx:#1a1d21;--mu:#6b7280;
  --del:#dc2626;--dup:#2563eb;--dbg:#fff0f0;--ubg:#eff4ff;--acc:#059669}
@media(prefers-color-scheme:dark){:root{--bg:#111317;--card:#1e2126;--bd:#2d3139;
  --tx:#d1d5db;--mu:#9ca3af;--dbg:#2a1515;--ubg:#141a2e}}
body{background:var(--bg);color:var(--tx);font-family:system-ui,sans-serif;font-size:13px}
a{color:var(--acc)}
header{background:var(--card);border-bottom:2px solid var(--bd);padding:18px 28px;
  display:flex;align-items:center;justify-content:space-between}
header h1{font-size:1.25rem;font-weight:700}
header .sub{color:var(--mu);font-size:.82rem;margin-top:2px}
.wrap{padding:22px 28px;max-width:1700px;margin:0 auto}
.cards{display:grid;grid-template-columns:repeat(auto-fit,minmax(170px,1fr));gap:14px;margin-bottom:20px}
.card{background:var(--card);border:1px solid var(--bd);border-radius:8px;padding:14px 18px}
.card h3{font-size:.72rem;text-transform:uppercase;letter-spacing:.05em;color:var(--mu);margin-bottom:6px}
.card .big{font-size:2rem;font-weight:700;line-height:1}
.si{display:flex;align-items:center;justify-content:space-between;margin:3px 0}
.sn{font-weight:700}
.badge{display:inline-block;padding:1px 7px;border-radius:10px;font-size:.78rem;font-weight:600}
.bdel{background:#fee2e2;color:var(--del)}.bdup{background:#dbeafe;color:var(--dup)}
.both{background:#e5e7eb;color:var(--mu)}
.panels{display:grid;grid-template-columns:1fr 1fr;gap:14px;margin-bottom:20px}
.panel{background:var(--card);border:1px solid var(--bd);border-radius:8px;padding:14px}
.panel h2{font-size:.88rem;font-weight:600;margin-bottom:10px;color:var(--tx)}
.mt{width:100%;border-collapse:collapse;font-size:.82rem}
.mt th,.mt td{padding:5px 8px;border-bottom:1px solid var(--bd);text-align:left}
.mt th{color:var(--mu);font-weight:600;font-size:.72rem;text-transform:uppercase}
.plots-grid{display:grid;grid-template-columns:repeat(2,1fr);
  gap:16px;margin-bottom:20px}
.plot-card{background:var(--card);border:1px solid var(--bd);border-radius:8px;
  padding:14px 16px;overflow:hidden}
.plot-card.full{grid-column:1/-1}
.plot-card h3{font-size:.72rem;font-weight:600;color:var(--mu);text-transform:uppercase;
  letter-spacing:.05em;margin-bottom:8px}
.plot-card img{width:100%;height:auto;display:block;border-radius:4px}
.section-label{font-size:.9rem;font-weight:700;color:var(--tx);margin:4px 0 12px}
.section{background:var(--card);border:1px solid var(--bd);border-radius:8px;overflow:hidden}
.sh{display:flex;align-items:center;justify-content:space-between;padding:12px 16px;border-bottom:1px solid var(--bd)}
.sh h2{font-size:.88rem;font-weight:600}
.srch{padding:5px 10px;border:1px solid var(--bd);border-radius:6px;background:var(--bg);
  color:var(--tx);font-size:.82rem;width:260px;outline:none}
.srch:focus{border-color:var(--acc)}
.tw{overflow-x:auto;max-height:72vh;overflow-y:auto}
table.main{width:100%;border-collapse:collapse;font-size:.80rem}
table.main thead tr{position:sticky;top:0;z-index:1}
table.main th{background:var(--bg);color:var(--mu);font-size:.70rem;text-transform:uppercase;
  letter-spacing:.04em;padding:7px 9px;text-align:left;cursor:pointer;white-space:nowrap;
  border-bottom:2px solid var(--bd);user-select:none}
table.main th:hover{color:var(--tx)}
table.main td{padding:6px 9px;border-bottom:1px solid var(--bd);vertical-align:top;
  max-width:260px;word-break:break-word}
.rdel{background:var(--dbg)}.rdup{background:var(--ubg)}
.hidden{display:none}
.si2{font-size:.65rem;opacity:.45}
th.asc .si2::after{content:" ▲";opacity:.8}th.desc .si2::after{content:" ▼";opacity:.8}
small{color:var(--mu);font-size:.78em}
.cnt{color:var(--mu);font-size:.78rem;padding:6px 16px}
@media(max-width:900px){.panels{grid-template-columns:1fr}
  header{flex-direction:column;align-items:flex-start;gap:6px}}
</style>'

  js <- '
<script>
function filterTable(q){
  q=q.toLowerCase();
  var rows=document.querySelectorAll("#tb tr");
  var n=0;
  rows.forEach(function(r){
    var hide=q!==""&&!r.textContent.toLowerCase().includes(q);
    r.classList.toggle("hidden",hide);
    if(!hide)n++;
  });
  document.getElementById("rowcnt").textContent=n+" variants shown";
}
var _dir={};
function srt(th){
  var col=th.cellIndex;
  var tbody=document.getElementById("tb");
  var rows=Array.from(tbody.querySelectorAll("tr"));
  var asc=_dir[col]!==true;_dir[col]=asc;
  document.querySelectorAll("table.main th").forEach(function(h){h.classList.remove("asc","desc")});
  th.classList.add(asc?"asc":"desc");
  rows.sort(function(a,b){
    var av=a.cells[col]?a.cells[col].textContent.trim():"";
    var bv=b.cells[col]?b.cells[col].textContent.trim():"";
    var an=parseFloat(av),bn=parseFloat(bv);
    if(!isNaN(an)&&!isNaN(bn))return asc?an-bn:bn-an;
    return asc?av.localeCompare(bv):bv.localeCompare(av);
  });
  rows.forEach(function(r){tbody.appendChild(r)});
}
</script>'

  body <- paste0(
    '<!DOCTYPE html>\n<html lang="en">\n<head>\n<meta charset="UTF-8">\n',
    '<meta name="viewport" content="width=device-width,initial-scale=1">\n',
    '<title>CNV Final Report</title>\n', css, '\n</head>\n<body>\n',
    '<header><div><h1>CNV Final Report</h1>',
    '<div class="sub">Generated: ', rep_date, '</div></div>',
    '<div class="sub">Bronco CNV Pipeline</div></header>\n',
    '<div class="wrap">\n',
    '<div class="cards">',
    '<div class="card"><h3>Samples in cohort</h3><div class="big">', n_samples, '</div></div>',
    '<div class="card"><h3>Total variants in cohort</h3><div class="big">', n_variants, '</div></div>',
    '<div class="card"><h3>By type</h3>', type_html, '</div>',
    '</div>\n',
    ## variants table — shown first
    '<div class="section">',
    '<div class="sh"><h2>Final variants table</h2>',
    '<input class="srch" type="search" placeholder="Search&hellip;" oninput="filterTable(this.value)"></div>',
    '<div class="cnt" id="rowcnt">', n_variants, ' variants shown</div>',
    '<div class="tw"><table class="main">',
    '<thead><tr>', header_html, '</tr></thead>',
    '<tbody id="tb">', data_rows_html, '</tbody>',
    '</table></div></div>\n',
    ## gene summary panels
    '<div class="panels" style="margin-top:20px">',
    '<div class="panel"><h2>Top genes</h2>', top_genes_html, '</div>',
    '<div class="panel"><h2>Top dosage-sensitive genes (with reported consequence)</h2>', dosage_html, '</div>',
    '</div>\n',
    ## plots
    if (plots_html != "") paste0('<p class="section-label">Plots</p>\n', plots_html) else "",
    '</div>\n', js, '\n</body>\n</html>\n'
  )

  writeLines(body, output_path)
  invisible(output_path)
}

generate_html_report(smoothed, file.path(output_dir, "final_report.html"))
cat("Wrote final_report.html\n")
