library(data.table)
library(vcfR)

extract_info_value <- function(info_vec, key) {
  pattern <- paste0("(^|;)", key, "=([^;]+)")
  out <- sub(pattern, "\\2", info_vec)
  out[!grepl(pattern, info_vec)] <- NA_character_
  out
}

normalize_svtype <- function(x) {
  x <- toupper(as.character(x))
  x[x %in% c("DELETION", "DEL", "<DEL>")] <- "DEL"
  x[x %in% c("DUPLICATION", "DUP", "<DUP>")] <- "DUP"
  x
}

read_all_tsv <- function(tsv_files) {

  file_info <- data.table(filepath = tsv_files)
  file_info[, filename := basename(filepath)]
  file_info <- file_info[grepl("_normalized\\.tsv$", filename)]
  file_info[, base := sub("_normalized\\.tsv$", "", filename)]

  valid_callers <- c("cnvkit", "cnMOPS", "ExomeDepth", "panelcnMOPS", "gatk", "jabcontool", "xhmm", "freec", "conifer")

  file_info[, varcaller := NA_character_]
  for (vc in valid_callers) {
    file_info[grepl(paste0("_", vc, "$"), base), varcaller := vc]
  }

  pattern_callers <- paste(valid_callers, collapse = "|")
  file_info[!is.na(varcaller),
            sample := sub(paste0("_(", pattern_callers, ")$"), "", base)]

  file_info <- file_info[!is.na(varcaller)]

  if (nrow(file_info) == 0) {
    return(data.table())
  }

  all_calls <- rbindlist(
    lapply(seq_len(nrow(file_info)), function(i) {
      dt <- fread(file_info$filepath[i], sep = "\t", header = TRUE)
      if (!nrow(dt)) return(NULL)
      dt[, sample := file_info$sample[i]]
      dt[, varcaller := file_info$varcaller[i]]
      dt
    }),
    fill = TRUE
  )

  all_calls
}

read_all_vcf <- function(vcf_files) {

  if (length(vcf_files) == 0) {
    return(data.table())
  }

  file_info <- data.table(filepath = vcf_files)
  file_info[, filename := basename(filepath)]
  file_info <- file_info[grepl("\\.vcf(\\.gz)?$", filename)]
  file_info[, base := sub("\\.vcf(\\.gz)?$", "", filename)]

  valid_callers <- c("gatk")

  file_info[, varcaller := NA_character_]
  for (vc in valid_callers) {
    file_info[grepl(paste0("_", vc, "$"), base), varcaller := vc]
  }

  pattern_callers <- paste(valid_callers, collapse = "|")
  file_info[!is.na(varcaller),
            sample := sub(paste0("_(", pattern_callers, ")$"), "", base)]

  file_info <- file_info[!is.na(varcaller)]

  if (nrow(file_info) == 0) {
    return(data.table())
  }

  vcf_calls <- rbindlist(
    lapply(seq_len(nrow(file_info)), function(i) {
      v <- read.vcfR(file_info$filepath[i], verbose = FALSE)
      fix <- as.data.table(getFIX(v))

      cn_mat <- extract.gt(v, element = "CN", as.numeric = TRUE)

      cn_vec <- as.numeric(cn_mat[, 1])

      dt <- data.table(
        CHROM = as.character(fix$CHROM),
        START = suppressWarnings(as.numeric(fix$POS)),
        END = fix[, as.numeric(sub(".*_", "", ID))],
        CN_ESTIMATE = cn_vec,
        sample = file_info$sample[i],
        varcaller = file_info$varcaller[i]
      )

      dt[is.na(END), END := START]
      dt <- dt[CN_ESTIMATE != 2,]
      dt
    }),
    fill = TRUE
  )
  vcf_calls
}

listofCHR <- c(
  "chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
  "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18",
  "chr19","chr20","chr21","chr22","chrX","chrY"
)

listofCHR2 <- c(
  "1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16",
  "17","18","19","20","21","22","X","Y"
)


run_all <- function(args){

    bed_file <- args[1]
    sample_name <- args[2]
    min_number_callers <- as.numeric(args[3])
    threshold_size <- as.numeric(args[4])

    tsv_files <- list.files(path = "./", pattern = "_normalized\\.tsv$", full.names = TRUE)
    tsv_files <- tsv_files[!grepl("final_CNV_probs_jabcontool_normalized\\.tsv$", tsv_files)]
    tsv_calls <- read_all_tsv(tsv_files)

    vcf_files <- list.files(path = "./", pattern = "\\.vcf(\\.gz)?$", full.names = TRUE)
    vcf_calls <- read_all_vcf(vcf_files)
    if (nrow(vcf_calls) > 0) {
      vcf_calls <- vcf_calls[CN_ESTIMATE != 2, TYPE := ifelse(CN_ESTIMATE < 2, "DEL", "DUP")]
    }

  if (file.exists("final_CNV_probs_jabcontool_normalized.tsv")) {
    jabcontool_calls <- fread("final_CNV_probs_jabcontool_normalized.tsv")
    jabcontool_calls[, varcaller := "jabcontool"]

  if ("SAMPLE" %in% colnames(jabcontool_calls)) {
    jabcontool_calls[, SAMPLE := sub("\\.region_coverage\\.tsv$", "", SAMPLE)]
    jabcontool_calls[, SAMPLE := sub("\\.bam$", "", SAMPLE)]  # por si acaso
    setnames(jabcontool_calls, "SAMPLE", "sample")
  }
}
    all_calls <- rbindlist(
      list(tsv_calls, vcf_calls, jabcontool_calls),
      fill = TRUE,
      use.names = TRUE
    )

    targets_raw <- fread(file = bed_file, sep = "\t", header = FALSE)
    if (ncol(targets_raw) < 3) {
    stop("BED file MUST contain at least three columns including: CHROM, START, END")
     }

     targets <- copy(targets_raw[, 1:3])
     setnames(targets, c("CHROM", "START", "END"))
     targets[, CHROM := as.character(CHROM)]
     targets[, START := as.numeric(START)]
     targets[, END := as.numeric(END)]

     if (grepl("^chr", targets$CHROM[1])) {
      genome_style <- "UCSC"
      targets <- targets[CHROM %in% listofCHR]
    } else {
      genome_style <- "ENSEMBL"
      targets <- targets[CHROM %in% listofCHR2]
    }
    targets[, target_id := paste0(CHROM, ":", START, "-", END)]

    if (ncol(targets_raw) >= 4) {
      targets[, target_name := as.character(targets_raw[[4]])]
    } else {
      targets[, target_name := target_id]
    }

    setDT(targets)
    setkey(targets, CHROM, START, END)

    required_cols <- c("CHROM", "START", "END", "CN_ESTIMATE", "TYPE", "sample", "varcaller")
    missing_cols <- setdiff(required_cols, colnames(all_calls))

    all_calls[, CHROM := as.character(CHROM)]
    all_calls[, START := as.numeric(START)]
    all_calls[, END := as.numeric(END)]
    all_calls[, CN_ESTIMATE := suppressWarnings(as.numeric(CN_ESTIMATE))]
    if ("ExomeDepth" %in% all_calls$varcaller) {
      all_calls[varcaller == "ExomeDepth" & !is.na(CN_ESTIMATE),
                CN_ESTIMATE := CN_ESTIMATE * 2]
    }
    all_calls[, TYPE := as.character(TYPE)]
    all_calls[toupper(TYPE) %in% c("DELETION", "DEL"), TYPE := "DEL"]
    all_calls[toupper(TYPE) %in% c("DUPLICATION", "DUP"), TYPE := "DUP"]

    if (genome_style == "UCSC") {
      all_calls[, CHROM := ifelse(grepl("^chr", CHROM), CHROM, paste0("chr", CHROM))]
      all_calls <- all_calls[CHROM %in% listofCHR]
    } else {
      all_calls[, CHROM := gsub("^chr", "", CHROM)]
      all_calls <- all_calls[CHROM %in% listofCHR2]
    }

    all_calls <- all_calls[(END - START) <= threshold_size]

    all_calls[, CALL_CHROM := CHROM]
    all_calls[, CALL_START := START]
    all_calls[, CALL_END := END]

    all_calls <- all_calls[!is.na(CHROM) & !is.na(START) & !is.na(END) & START <= END]

    setDT(all_calls)
    setkey(all_calls, CHROM, START, END)

    targets <- targets[!is.na(CHROM) & !is.na(START) & !is.na(END) & START <= END]
    setkey(targets, CHROM, START, END)

    ov <- foverlaps(
      x = targets,
      y = all_calls,
      by.x = c("CHROM", "START", "END"),
      type = "any",
      mult = "all",
      nomatch = NULL
    )

    if (!nrow(ov)) {
      stop("There were no overlaps between the target regions on the BED file and the CNV calls")
    }

    ov_summary <- ov[, .(
      CN_ESTIMATE = median(CN_ESTIMATE, na.rm = TRUE),
      TYPE = {
        t <- unique(na.omit(TYPE))
        if (length(t) == 1) {
          t
        } else if (length(t) == 0) {
          NA_character_
        } else {
          "CONFLICT"
        }
      },
      caller_coords = paste(
        unique(paste0(CALL_CHROM, ":", CALL_START, "-", CALL_END)),
        collapse = ";"
      )
    ), by = .(sample, varcaller, target_id)]

    target_meta <- unique(targets[, .(target_id, CHROM, START, END, target_name)])

    ov_summary <- merge(
      ov_summary,
      target_meta,
      by = "target_id",
      all.x = TRUE
    )

    wide_types <- dcast(
      ov_summary,
      sample + target_id + CHROM + START + END + target_name ~ varcaller,
      value.var = "TYPE"
    )

    wide_coords <- dcast(
      ov_summary,
      sample + target_id ~ varcaller,
      value.var = "caller_coords"
    )

    coord_cols <- setdiff(colnames(wide_coords), c("sample", "target_id"))
    if (length(coord_cols) > 0) {
      setnames(wide_coords, coord_cols, paste0(coord_cols, "_coords"))
    }

    wide_cn <- dcast(
      ov_summary,
      sample + target_id ~ varcaller,
      value.var = "CN_ESTIMATE"
    )

    cn_cols <- setdiff(colnames(wide_cn), c("sample", "target_id"))
    if (length(cn_cols) > 0) {
      setnames(wide_cn, cn_cols, paste0(cn_cols, "_CN"))
    }

    wide_calls <- merge(
      wide_types,
      wide_coords,
      by = c("sample", "target_id"),
      all.x = TRUE
    )

    wide_calls <- merge(
      wide_calls,
      wide_cn,
      by = c("sample", "target_id"),
      all.x = TRUE
    )

    type_cols <- intersect(
      unique(all_calls$varcaller),
      colnames(wide_calls)
    )

    min_number_callers <- min(min_number_callers, length(type_cols))

    for (col in type_cols) {
      wide_calls[is.na(get(col)), (col) := "0"]
    }

    wide_calls[, Count_Detected := rowSums(as.matrix(.SD) != "0"), .SDcols = type_cols]

    dir.create("merged_variants", showWarnings = FALSE)
    write.table(
      wide_calls,
      file  = paste0("merged_variants/", sample_name, "_raw_callers_matrix.tsv"),
      sep   = "\t",
      quote = FALSE,
      row.names = FALSE,
      col.names = TRUE
    )

    wide_calls <- wide_calls[Count_Detected >= min_number_callers]
    wide_calls[, n_DEL := rowSums(.SD == "DEL", na.rm = TRUE), .SDcols = type_cols]
    wide_calls[, n_DUP := rowSums(.SD == "DUP", na.rm = TRUE), .SDcols = type_cols]
    wide_calls[, consensus_type := NA_character_]
    wide_calls[n_DEL >= 2 & n_DUP < 2, consensus_type := "DEL"]
    wide_calls[n_DUP >= 2 & n_DEL < 2, consensus_type := "DUP"]
    wide_calls[n_DEL >= 2 & n_DUP >= 2, consensus_type := "CONFLICT"]
    coord_cols_present <- intersect(paste0(type_cols, "_coords"), colnames(wide_calls))
    cn_cols_present <- intersect(paste0(type_cols, "_CN"), colnames(wide_calls))


    for (col in coord_cols_present) {
      wide_calls[is.na(get(col)), (col) := "0"]
    }

    for (col in cn_cols_present) {
      wide_calls[is.na(get(col)), (col) := 0]
    }

    final_calls <- wide_calls[consensus_type %in% c("DEL", "DUP")]

    final_bed <- final_calls[, .(
      CHROM,
      START,
      END,
      consensus_type
    )]

    final_tsv <- final_calls
    setorder(final_tsv, CHROM, START)

    write.table(
      final_bed,
      file = paste0("merged_variants/", sample_name, "_merged_target_consensus.bed"),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE,
      col.names = FALSE
    )

    write.table(
      final_tsv,
      file = paste0("merged_variants/", sample_name, "_merged_target_consensus.tsv"),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE,
      col.names = TRUE
    )

    # === Smoothed independent variant calls ===
    # Merges overlapping target hits into single variant calls using the actual
    # caller coordinates: START = earliest start across all callers, END = farthest
    # end across all callers, so the smoothed call spans the full reported extent.
    if (nrow(final_calls) > 0 && length(coord_cols_present) > 0) {

        parse_bounds <- function(row_vec, type, max_span = 500000) {
            valid <- row_vec[!is.na(row_vec) & row_vec != "0" & nchar(as.character(row_vec)) > 0]
            if (!length(valid)) return(NA_real_)
            txt <- paste(valid, collapse = ";")
            # coords format: "CHROM:START-END" — capture the ":digits-digits" part
            m <- regmatches(txt, gregexpr(":([0-9]+)-([0-9]+)", txt))[[1]]
            if (!length(m)) return(NA_real_)
            parts  <- strsplit(sub("^:", "", m), "-")
            starts <- as.numeric(sapply(parts, `[`, 1))
            ends   <- as.numeric(sapply(parts, `[`, 2))
            # drop any individual coord that spans more than max_span (likely an artifact)
            keep   <- (ends - starts) <= max_span
            starts <- starts[keep]
            ends   <- ends[keep]
            if (!length(starts)) return(NA_real_)
            if (type == "start") min(starts, na.rm = TRUE) else max(ends, na.rm = TRUE)
        }

        caller_type_cols <- intersect(type_cols, colnames(final_calls))

        fc_coord_data <- final_calls[, coord_cols_present, with = FALSE]
        call_starts   <- apply(fc_coord_data, 1, parse_bounds, type = "start")
        call_ends     <- apply(fc_coord_data, 1, parse_bounds, type = "end")

        smoothed <- cbind(
            final_calls[, .(sample, CHROM, target_name, consensus_type,
                             target_START = START, target_END = END)],
            final_calls[, c(caller_type_cols, coord_cols_present), with = FALSE]
        )
        smoothed[, call_start := fifelse(is.na(call_starts), target_START, as.numeric(call_starts))]
        smoothed[, call_end   := fifelse(is.na(call_ends),   target_END,   as.numeric(call_ends))]

        setorder(smoothed, sample, CHROM, consensus_type, call_start)

        # Build a key of non-trivial caller coord strings per row for exact-match merging.
        # Two rows belong to the same variant only if at least one caller reported the
        # identical coordinate string for both (same call mapped to multiple targets).
        # This prevents nearby but distinct variants from being merged just because their
        # caller coordinates happen to overlap.
        smoothed[, coord_key := {
            apply(.SD, 1, function(row) {
                vals <- row[!is.na(row) & row != "0" & nchar(as.character(row)) > 0]
                if (!length(vals)) "" else paste(sort(vals), collapse = "|")
            })
        }, .SDcols = coord_cols_present]

        smoothed[, var_group := {
            n <- .N
            if (n == 0L) return(integer(0))
            g       <- integer(n)
            cur_grp <- 1L
            g[1L]   <- 1L
            # group_keys: union of all coord strings seen so far per group
            group_keys <- list()
            group_keys[[1L]] <- strsplit(coord_key[1L], "\\|")[[1]]
            if (n > 1L) {
                for (i in 2L:n) {
                    row_keys <- strsplit(coord_key[i], "\\|")[[1]]
                    # find the most recent group whose coord set overlaps this row
                    matched_grp <- NA_integer_
                    for (grp in cur_grp:1L) {
                        if (length(intersect(row_keys, group_keys[[grp]])) > 0L) {
                            matched_grp <- grp
                            break
                        }
                    }
                    if (!is.na(matched_grp)) {
                        g[i] <- matched_grp
                        group_keys[[matched_grp]] <- union(group_keys[[matched_grp]], row_keys)
                    } else {
                        cur_grp <- cur_grp + 1L
                        g[i]    <- cur_grp
                        group_keys[[cur_grp]] <- row_keys
                    }
                }
            }
            g
        }, by = .(sample, CHROM, consensus_type)]

        smoothed[, coord_key := NULL]

        # Aggregate each merged group into one variant row
        smooth_output <- smoothed[, {
            active <- caller_type_cols[
                vapply(caller_type_cols,
                       function(col) any(get(col) != "0", na.rm = TRUE),
                       logical(1))
            ]
            .(
                START        = min(call_start),
                END          = max(call_end),
                n_targets    = .N,
                n_callers    = length(active),
                callers      = paste(active, collapse = ";"),
                target_names = paste(sort(unique(target_name)), collapse = ";")
            )
        }, by = .(sample, CHROM, consensus_type, var_group)]

        smooth_output[, var_group := NULL]

        # Deduplicate: rows with identical coords after merging are the same variant
        smooth_output <- unique(smooth_output, by = c("sample", "CHROM", "START", "END", "consensus_type"))

        setorder(smooth_output, sample, CHROM, START)

        write.table(
            smooth_output,
            file      = paste0("merged_variants/", sample_name, "_smoothed_variants.tsv"),
            sep       = "\t",
            quote     = FALSE,
            row.names = FALSE,
            col.names = TRUE
        )
    }

}

    args <- commandArgs(trailingOnly = T)
    run_all(args)