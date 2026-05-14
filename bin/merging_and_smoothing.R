suppressPackageStartupMessages({
  library(GenomicRanges)
  library(rtracklayer)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(VariantAnnotation)
  library(optparse)
})

bf_to_prob <- function(bf) {
  bf <- as.numeric(bf)
  bf / (bf + 1)
}

cn_from_log2 <- function(log2_ratio) {
  round(2 * (2 ^ as.numeric(log2_ratio)))
}

cn_from_ratio <- function(ratio) {
  round(2 * as.numeric(ratio))
}

cn_string_to_int <- function(cn_str) {
  as.integer(str_replace(cn_str, "CN", ""))
}

direction_from_cn <- function(cn) {
  ifelse(cn > 2, 1L, ifelse(cn < 2, -1L, 0L))
}

filter_by_size <- function(gr, max_size, tool_name) {
  if (length(gr) == 0 || is.null(max_size)) return(gr)
  sizes  <- width(gr)
  keep   <- sizes <= max_size
  n_drop <- sum(!keep)
  if (n_drop > 0) {
    cat(sprintf("    [size filter] %s: %d call(s) removed(s) (> %s nt)\n",
                tool_name, n_drop,
                formatC(max_size, format = "d", big.mark = ",")))
  }
  gr[keep]
}

parse_cnvkit <- function(file, max_size = NULL) {
  dt <- fread(file)
  cn <- cn_from_log2(dt$FOLD_CHANGE_LOG)
  filter_by_size(
    GRanges(
      seqnames   = as.character(dt$CHROM),
      ranges     = IRanges(start = dt$POS, end = dt$END),
      tool       = "cnvkit",
      direction  = ifelse(dt$SVTYPE == "DUP", 1L, -1L),
      magnitude  = abs(as.numeric(dt$FOLD_CHANGE_LOG)),  
      confidence = pmin(as.numeric(dt$PROBES) / 50, 1),  
      cn_estimate = cn,
      qual_flag  = 0L,
      call_id    = paste0("cnvkit_", seq_len(nrow(dt)))
    ),
    max_size, "cnvkit"
  )
}

parse_exomedepth <- function(file, max_size = NULL) {
  dt <- fread(file)
  filter_by_size(
    GRanges(
      seqnames   = as.character(dt$chromosome),
      ranges     = IRanges(start = dt$start, end = dt$end),
      tool       = "exomedepth",
      direction  = ifelse(dt$type == "duplication", 1L, -1L),
      magnitude  = abs(log2(pmax(as.numeric(dt$reads.ratio), 0.01))),
      confidence = bf_to_prob(dt$BF),
      cn_estimate = cn_from_ratio(dt$reads.ratio),
      qual_flag  = 0L,
      call_id    = paste0("exomedepth_", seq_len(nrow(dt)))
    ),
    max_size, "exomedepth"
  )
}

parse_cnmops <- function(file, max_size = NULL) {
  dt <- fread(file)
  cn <- cn_string_to_int(dt$CN)
  filter_by_size(
    GRanges(
      seqnames   = as.character(dt$chromosome),
      ranges     = IRanges(start = dt$start, end = dt$end),
      tool       = "cnmops",
      direction  = direction_from_cn(cn),
      magnitude  = abs(cn - 2) / 2,
      confidence = ifelse(abs(as.numeric(dt$mean)) > 1, 0.9, 0.5),
      cn_estimate = cn,
      qual_flag  = 0L,
      call_id    = paste0("cnmops_", seq_len(nrow(dt)))
    ),
    max_size, "cnmops"
  )
}

parse_panelcnmops <- function(file, max_size = NULL) {
  dt <- fread(file)
  cn <- cn_string_to_int(dt$CN)
  qual_flag <- ifelse(str_detect(as.character(dt$lowQual), "lowQual"), 1L, 0L)
  rc_norm_ratio <- as.numeric(dt$RC.norm) / pmax(as.numeric(dt$medRC.norm), 1)
  filter_by_size(
    GRanges(
      seqnames   = as.character(dt$chromosome),
      ranges     = IRanges(start = dt$start, end = dt$end),
      tool       = "panelcnmops",
      direction  = direction_from_cn(cn),
      magnitude  = abs(log2(pmax(rc_norm_ratio, 0.01))),
      confidence = ifelse(qual_flag == 1, 0.3, 0.9),
      cn_estimate = cn,
      qual_flag  = qual_flag,
      call_id    = paste0("panelcnmops_", seq_len(nrow(dt)))
    ),
    max_size, "panelcnmops"
  )
}

parse_freec <- function(file, max_size = NULL) {
  dt <- fread(file)
  cn <- as.integer(dt$CN)
  filter_by_size(
    GRanges(
      seqnames   = as.character(dt$chromosome),
      ranges     = IRanges(start = dt$start, end = dt$end),
      tool       = "freec",
      direction  = direction_from_cn(cn),
      magnitude  = abs(cn - 2) / 2,
      confidence = 0.7,
      cn_estimate = cn,
      qual_flag  = 0L,
      call_id    = paste0("freec_", seq_len(nrow(dt)))
    ),
    max_size, "freec"
  )
}

parse_xhmm <- function(file, max_size = NULL) {
  dt <- fread(file)
  cn <- as.integer(dt$CN)
  filter_by_size(
    GRanges(
      seqnames   = as.character(dt$chromosome),
      ranges     = IRanges(start = dt$start, end = dt$end),
      tool       = "xhmm",
      direction  = direction_from_cn(cn),
      magnitude  = abs(cn - 2) / 2,
      confidence = 0.7,
      cn_estimate = cn,
      qual_flag  = 0L,
      call_id    = paste0("xhmm_", seq_len(nrow(dt)))
    ),
    max_size, "xhmm"
  )
}

parse_conifer <- function(file, max_size = NULL) {
  dt <- fread(file)
  if (nrow(dt) == 0) {
    return(GRanges(
      tool = character(0), direction = integer(0),
      magnitude = numeric(0), confidence = numeric(0),
      cn_estimate = integer(0), qual_flag = integer(0),
      call_id = character(0)
    ))
  }
  cn <- as.integer(dt$CN)
  filter_by_size(
    GRanges(
      seqnames    = as.character(dt$chromosome),
      ranges      = IRanges(start = dt$start, end = dt$end),
      tool        = "conifer",
      direction   = direction_from_cn(cn),
      magnitude   = abs(cn - 2) / 2,
      confidence  = 0.7,
      cn_estimate = cn,
      qual_flag   = 0L,
      call_id     = paste0("conifer_", seq_len(nrow(dt)))
    ),
    max_size, "conifer"
  )
}

parse_jabcontool <- function(file, max_size = NULL, sample_id = NULL) {
  dt <- fread(file)
  if (!is.null(sample_id) && "sample" %in% names(dt)) {
    expected <- paste0(sample_id, ".region_coverage.tsv")
    dt <- dt[sample == expected]
    if (nrow(dt) == 0) {
      cat(sprintf("    JabConTool: sample '%s' not found en %s\n",
                  sample_id, basename(file)))
      return(GRanges(
        tool = character(0), direction = integer(0),
        magnitude = numeric(0), confidence = numeric(0),
        cn_estimate = integer(0), qual_flag = integer(0),
        call_id = character(0)
      ))
    }
    cat(sprintf("    JabConTool: %d rows for sample '%s'\n", nrow(dt), sample_id))
  }
  
  region_rows <- !is.na(dt$cov) & dt$cov != "" &
    (is.na(dt$pos) | dt$pos == "")
  dt <- dt[region_rows]
  if (nrow(dt) == 0) {
    return(GRanges(
      tool = character(0), direction = integer(0),
      magnitude = numeric(0), confidence = numeric(0),
      cn_estimate = integer(0), qual_flag = integer(0),
      call_id = character(0)
    ))
  }
  
  cn_pred <- as.integer(dt$cn_pred)
  is_cnv  <- cn_pred != 2
  dt      <- dt[is_cnv]
  cn_pred <- cn_pred[is_cnv]
  if (nrow(dt) == 0) {
    return(GRanges(
      tool = character(0), direction = integer(0),
      magnitude = numeric(0), confidence = numeric(0),
      cn_estimate = integer(0), qual_flag = integer(0),
      call_id = character(0)
    ))
  }
  
  ll_cols <- c("0", "1", "2", "3", "4", "5")
  ll_mat  <- as.matrix(dt[, ..ll_cols])
  storage.mode(ll_mat) <- "numeric"
  
  ll_pred   <- ll_mat[cbind(seq_len(nrow(dt)), cn_pred + 1L)] 
  ll_normal <- ll_mat[, 3]                                      
  margin    <- ll_pred - ll_normal                              
  
  confidence <- 1 / (1 + exp(-margin))
  
  filter_by_size(
    GRanges(
      seqnames    = as.character(dt$chr),
      ranges      = IRanges(start = dt$start, end = dt$end),
      tool        = "jabcontool",
      direction   = direction_from_cn(cn_pred),
      magnitude   = abs(cn_pred - 2) / 2,
      confidence  = confidence,
      cn_estimate = cn_pred,
      qual_flag   = 0L,
      call_id     = paste0("jabcontool_", seq_len(nrow(dt)))
    ),
    max_size, "jabcontool"
  )
}

parse_gatk_gcnv <- function(file, max_size = NULL) {
  if (!file.exists(file)) {
    return(GRanges(
      tool = character(0), direction = integer(0),
      magnitude = numeric(0), confidence = numeric(0),
      cn_estimate = integer(0), qual_flag = integer(0),
      call_id = character(0)
    ))
  }
  
  vcf <- readVcf(file)
  if (length(vcf) == 0) {
    return(GRanges(
      tool = character(0), direction = integer(0),
      magnitude = numeric(0), confidence = numeric(0),
      cn_estimate = integer(0), qual_flag = integer(0),
      call_id = character(0)
    ))
  }
  
  rr <- rowRanges(vcf)
  end_vec <- info(vcf)$END
  if (length(end_vec) == 0 || all(is.na(end_vec))) {
    end_vec <- end(rr)
  }
  
  geno_cn <- geno(vcf)$CN
  geno_qa <- if ("QA" %in% names(geno(vcf))) geno(vcf)$QA else NULL
  geno_qs <- if ("QS" %in% names(geno(vcf))) geno(vcf)$QS else NULL
  
  cn <- as.integer(geno_cn[, 1])
  qa <- if (!is.null(geno_qa)) as.numeric(geno_qa[, 1]) else rep(NA_real_, length(cn))
  qs <- if (!is.null(geno_qs)) as.numeric(geno_qs[, 1]) else rep(NA_real_, length(cn))
  
  is_cnv <- !is.na(cn) & cn != 2
  if (!any(is_cnv)) {
    return(GRanges(
      tool = character(0), direction = integer(0),
      magnitude = numeric(0), confidence = numeric(0),
      cn_estimate = integer(0), qual_flag = integer(0),
      call_id = character(0)
    ))
  }
  
  rr_f      <- rr[is_cnv]
  cn_f      <- cn[is_cnv]
  qa_f      <- qa[is_cnv]
  qs_f      <- qs[is_cnv]
  end_vec_f <- end_vec[is_cnv]
  
  q_score   <- ifelse(!is.na(qs_f), qs_f, qa_f)
  confidence <- ifelse(is.na(q_score),
                       0.7,
                       1 - 10^(-q_score / 10))   
  confidence <- pmin(pmax(confidence, 0), 1)    
  
  chrs <- as.character(seqnames(rr_f))
  chrs <- sub("^chr", "", chrs)
  
  filter_by_size(
    GRanges(
      seqnames    = chrs,
      ranges      = IRanges(start = start(rr_f), end = as.integer(end_vec_f)),
      tool        = "gatk_gcnv",
      direction   = direction_from_cn(cn_f),
      magnitude   = abs(cn_f - 2) / 2,
      confidence  = confidence,
      cn_estimate = cn_f,
      qual_flag   = 0L,
      call_id     = paste0("gatk_gcnv_", seq_along(cn_f))
    ),
    max_size, "gatk_gcnv"
  )
}


load_capture_exons <- function(bed_file) {
  bed <- fread(bed_file, header = FALSE)
  if (ncol(bed) >= 4) {
    setnames(bed, 1:4, c("chr", "start", "end", "exon_id"))
  } else {
    setnames(bed, 1:3, c("chr", "start", "end"))
    bed[, exon_id := paste0("target_", seq_len(.N))]
  }
  GRanges(
    seqnames = as.character(bed$chr),
    ranges   = IRanges(start = bed$start + 1, end = bed$end),
    exon_id  = bed$exon_id
  )
}

load_cds_from_gtf <- function(gtf_file, capture_gr = NULL) {
  
  cat("    Reading GTF (all sources, only CDS feature)...\n")
  
  awk_filter <- paste0("awk ", shQuote('$3 == "CDS"'))
  cmd <- paste0("grep -v ", shQuote("^#"), " ", shQuote(gtf_file),
                " | ", awk_filter)
  cds_raw <- fread(cmd = cmd, header = FALSE, sep = "\t",
                   col.names = c("chr","source","feature","start","end",
                                 "score","strand","frame","attributes"))
  cat(sprintf("    %d CDS features loaded \n", nrow(cds_raw)))
  
  cds_raw[, gene_name := {
    attrs    <- attributes
    has_name <- grepl('gene_name "', attrs, fixed = TRUE)
    result   <- character(length(attrs))
    if (any(has_name))
      result[has_name] <- gsub('.*gene_name "([^"]+)".*', "\\1",
                               attrs[has_name], perl = TRUE)
    if (any(!has_name))
      result[!has_name] <- gsub('.*gene_id "([^"]+)".*', "\\1",
                                attrs[!has_name], perl = TRUE)
    result
  }]
  
  cds_raw[, chr_clean := as.character(chr)]
  
  cds_raw[, exon_num := gsub('.*exon_number "([^"]+)".*', "\\1",
                             attributes, perl = TRUE)]
  
  cds_raw[, exon_id := paste0(gene_name, ":", chr_clean, ":",
                              start, "-", end, ":exon", exon_num)]
  
  cds_unique <- unique(cds_raw, by = c("chr_clean", "start", "end"))
  cat(sprintf("    %d unique CDS after coordinate deduplication \n", nrow(cds_unique)))
  
  gr <- GRanges(
    seqnames  = cds_unique$chr_clean,
    ranges    = IRanges(start = cds_unique$start, end = cds_unique$end),
    strand    = cds_unique$strand,
    exon_id   = cds_unique$exon_id,
    gene_name = cds_unique$gene_name
  )
  
  if (!is.null(capture_gr)) {
    hits <- findOverlaps(gr, capture_gr, ignore.strand = TRUE)
    keep <- unique(queryHits(hits))
    gr   <- gr[keep]
    cat(sprintf("    %d CDS overlap with the regions of interest\n", length(gr)))
  }
  
  gr
}

project_tool_to_exons <- function(tool_calls_gr, exons_gr, tool_name) {
  if (length(tool_calls_gr) == 0) {
    empty <- data.table(
      exon_id              = character(0),
      n_calls              = integer(0),
      direction            = integer(0),
      direction_consistent = integer(0),
      magnitude_max        = numeric(0),
      magnitude_mean       = numeric(0),
      confidence_max       = numeric(0),
      confidence_mean      = numeric(0),
      cn_estimate          = integer(0),
      qual_flag            = integer(0),
      call_size_log10      = numeric(0),
      exon_coverage        = numeric(0)
    )
    feat_cols <- setdiff(names(empty), "exon_id")
    setnames(empty, feat_cols, paste0(tool_name, "_", feat_cols))
    empty[, paste0(tool_name, "_called") := integer(0)]
    return(empty)
  }
  
  hits <- findOverlaps(exons_gr, tool_calls_gr)
  
  if (length(hits) == 0) {
    empty <- data.table(
      exon_id              = character(0),
      n_calls              = integer(0),
      direction            = integer(0),
      direction_consistent = integer(0),
      magnitude_max        = numeric(0),
      magnitude_mean       = numeric(0),
      confidence_max       = numeric(0),
      confidence_mean      = numeric(0),
      cn_estimate          = integer(0),
      qual_flag            = integer(0),
      call_size_log10      = numeric(0),
      exon_coverage        = numeric(0)
    )
    feat_cols <- setdiff(names(empty), "exon_id")
    setnames(empty, feat_cols, paste0(tool_name, "_", feat_cols))
    empty[, paste0(tool_name, "_called") := integer(0)]
    return(empty)
  }
  
  exon_idx <- queryHits(hits)
  call_idx <- subjectHits(hits)
  
  exon_starts    <- start(exons_gr)[exon_idx]
  exon_ends      <- end(exons_gr)[exon_idx]
  call_starts    <- start(tool_calls_gr)[call_idx]
  call_ends      <- end(tool_calls_gr)[call_idx]
  exon_widths    <- exon_ends - exon_starts + 1
  call_widths    <- call_ends - call_starts + 1
  overlap_widths <- pmax(0, pmin(exon_ends, call_ends) - pmax(exon_starts, call_starts) + 1)
  
  hits_dt <- data.table(
    exon_id       = exons_gr$exon_id[exon_idx],
    direction     = tool_calls_gr$direction[call_idx],
    magnitude     = tool_calls_gr$magnitude[call_idx],
    confidence    = tool_calls_gr$confidence[call_idx],
    cn_estimate   = tool_calls_gr$cn_estimate[call_idx],
    qual_flag     = tool_calls_gr$qual_flag[call_idx],
    call_size     = call_widths,
    exon_coverage = overlap_widths / pmax(exon_widths, 1),
    call_id       = tool_calls_gr$call_id[call_idx],
    call_chr      = as.character(seqnames(tool_calls_gr)[call_idx]),
    call_start    = call_starts,
    call_end      = call_ends
  )
  
  agg <- hits_dt[, {
    best <- which.max(confidence)
    .(
      n_calls              = .N,
      direction            = direction[best],
      direction_consistent = as.integer(length(unique(direction)) == 1),
      magnitude_max        = max(magnitude, na.rm = TRUE),
      magnitude_mean       = mean(magnitude, na.rm = TRUE),
      confidence_max       = max(confidence, na.rm = TRUE),
      confidence_mean      = mean(confidence, na.rm = TRUE),
      cn_estimate          = as.integer(median(cn_estimate, na.rm = TRUE)),
      qual_flag            = max(qual_flag, na.rm = TRUE),
      call_size_log10      = log10(max(call_size, na.rm = TRUE)),
      exon_coverage        = max(exon_coverage, na.rm = TRUE),
      call_chr             = call_chr[best],
      call_start           = call_start[best],
      call_end             = call_end[best]
    )
  }, by = exon_id]
  
  feat_cols <- setdiff(names(agg), "exon_id")
  setnames(agg, feat_cols, paste0(tool_name, "_", feat_cols))
  
  agg[, paste0(tool_name, "_called") := 1L]
  
  agg
}

build_per_exon_matrix <- function(tools_gr_list, exons_gr) {
    base <- data.table(
    exon_id   = exons_gr$exon_id,
    chr       = as.character(seqnames(exons_gr)),
    start     = start(exons_gr),
    end       = end(exons_gr),
    exon_size = width(exons_gr)
  )
  
  for (tool_name in names(tools_gr_list)) {
    cat(sprintf("  Projecting %s...\n", tool_name))
    proj <- project_tool_to_exons(tools_gr_list[[tool_name]],
                                  exons_gr, tool_name)
    base <- merge(base, proj, by = "exon_id", all.x = TRUE)
  }
  
  tool_names <- names(tools_gr_list)
  for (tool_name in tool_names) {
    cols_called    <- paste0(tool_name, "_called")
    cols_dir       <- paste0(tool_name, "_direction")
    cols_consist   <- paste0(tool_name, "_direction_consistent")
    cols_numeric   <- paste0(tool_name, c("_n_calls", "_magnitude_max",
                                          "_magnitude_mean", "_confidence_max",
                                          "_confidence_mean", "_qual_flag",
                                          "_call_size_log10", "_exon_coverage"))
    cols_cn        <- paste0(tool_name, "_cn_estimate")
    
    if (cols_called  %in% names(base)) base[is.na(get(cols_called)), (cols_called) := 0L]
    if (cols_dir     %in% names(base)) base[is.na(get(cols_dir)),    (cols_dir)    := 0L]
    if (cols_consist %in% names(base)) base[is.na(get(cols_consist)),(cols_consist):= 1L]
    if (cols_cn      %in% names(base)) base[is.na(get(cols_cn)),     (cols_cn)     := 2L]
    for (cn in cols_numeric) {
      if (cn %in% names(base)) base[is.na(get(cn)), (cn) := 0]
    }
  }
  
  called_cols <- grep("_called$", names(base), value = TRUE)
  dir_cols    <- grep("_direction$", names(base), value = TRUE)
  conf_cols   <- grep("_confidence_max$", names(base), value = TRUE)
  
  base[, n_tools_called := rowSums(.SD), .SDcols = called_cols]
  base[, n_tools_DEL := rowSums(.SD == -1), .SDcols = dir_cols]
  base[, n_tools_DUP := rowSums(.SD ==  1), .SDcols = dir_cols]
  base[, direction_net := n_tools_DUP - n_tools_DEL]
  base[, direction_consensus := ifelse(n_tools_DEL > n_tools_DUP, -1L,
                                       ifelse(n_tools_DUP > n_tools_DEL,  1L, 0L))]
  base[, has_conflict := as.integer(n_tools_DEL > 0 & n_tools_DUP > 0)]
  base[, sum_confidence := rowSums(.SD), .SDcols = conf_cols]
  
  base
}

assign_labels_from_gt <- function(per_exon_dt, gt_bed_file, exons_gr) {
  gt_bed <- fread(gt_bed_file, header = FALSE)
  setnames(gt_bed, 1:4, c("chr", "start", "end", "type"))
  
  gt_gr <- GRanges(
    seqnames = as.character(gt_bed$chr),
    ranges   = IRanges(start = gt_bed$start + 1, end = gt_bed$end),
    type     = gt_bed$type
  )
  
  hits <- findOverlaps(exons_gr, gt_gr)
  exon_label_dt <- data.table(
    exon_id   = exons_gr$exon_id[queryHits(hits)],
    gt_type   = gt_gr$type[subjectHits(hits)]
  )
  
  agg_labels <- exon_label_dt[, .(
    label = 1L,
    gt_direction = ifelse("DUP" %in% toupper(gt_type), 1L,
                          ifelse("DEL" %in% toupper(gt_type) | 
                                   "deletion" %in% gt_type, -1L, 0L))
  ), by = exon_id]
  
  per_exon_dt <- merge(per_exon_dt, agg_labels, by = "exon_id", all.x = TRUE)
  per_exon_dt[is.na(label), label := 0L]
  per_exon_dt[is.na(gt_direction), gt_direction := 0L]
  per_exon_dt
}

run_all <- function(input_dir,
                         capture_bed,
                         output_dir,
                         sample_id       = "sample1",
                         gt_bed          = NULL,
                         jabcontool_file = NULL,
                         gatk_vcf_file   = NULL,
                         file_overrides  = list(),
                         min_callers     = 2L,
                         smooth_gap_bp   = 500000L,
                         max_call_size   = 500000L,
                         gtf_file        = NULL) {


  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  cat("==> 1. Loading reference...\n")
  if (!is.null(gtf_file) && file.exists(gtf_file)) {
    cat("    GTF detected: using exon CDS as analysis unit\n")
    capture_gr <- load_capture_exons(capture_bed)   
    exons_gr   <- load_cds_from_gtf(gtf_file, capture_gr = capture_gr)
    ref_label  <- "CDS exons"
  } else {
    cat("    BED-only mode: using kit target regions as analysis unit\n")
    exons_gr  <- load_capture_exons(capture_bed)
    ref_label <- "target regions"
  }
  cat(sprintf("    %d %s loaded as reference unit\n",
              length(exons_gr), ref_label))
  
  cat("==> 2. Loading calls per tool...\n")
  
  files <- list(
    cnvkit       = file.path(input_dir, paste0(sample_id, "_cnvkit.tsv")),
    exomedepth   = file.path(input_dir, paste0(sample_id, "_ExomeDepth.tsv")),
    cnmops       = file.path(input_dir, paste0(sample_id, "_cnMOPS.tsv")),
    panelcnmops  = file.path(input_dir, paste0(sample_id, "_panelcnMOPS.tsv")),
    freec        = file.path(input_dir, paste0(sample_id, "_freec.tsv")),
    xhmm         = file.path(input_dir, paste0(sample_id, "_xhmm.tsv")),
    conifer      = file.path(input_dir, paste0(sample_id, "_conifer.tsv")),
    jabcontool   = if (!is.null(jabcontool_file)) jabcontool_file
    else file.path(input_dir, "jabcontool_test.tsv"),
    gatk_gcnv    = if (!is.null(gatk_vcf_file)) gatk_vcf_file
    else file.path(input_dir, paste0(sample_id, "_gatk.vcf.gz"))
  )
  
  for (k in names(file_overrides)) {
    files[[k]] <- file_overrides[[k]]
  }
  
  parsers <- list(
    cnvkit       = parse_cnvkit,
    exomedepth   = parse_exomedepth,
    cnmops       = parse_cnmops,
    panelcnmops  = parse_panelcnmops,
    freec        = parse_freec,
    xhmm         = parse_xhmm,
    conifer      = parse_conifer,
    jabcontool   = parse_jabcontool,
    gatk_gcnv    = parse_gatk_gcnv
  )
  
  tools_gr_list <- list()
  for (tool in names(files)) {
    if (file.exists(files[[tool]])) {
      cat(sprintf("    %s... ", tool))
      if (tool == "jabcontool") {
        tools_gr_list[[tool]] <- parsers[[tool]](files[[tool]],
                                                 max_size  = max_call_size,
                                                 sample_id = sample_id)
      } else {
        tools_gr_list[[tool]] <- parsers[[tool]](files[[tool]], max_size = max_call_size)
      }
      cat(sprintf("%d calls\n", length(tools_gr_list[[tool]])))
    } else {
      cat(sprintf("    %s SKIPPED (does not exist %s)\n", tool, files[[tool]]))
    }
  }
  
  cat("==> 3. Building per-exon matrix...\n")
  per_exon_dt <- build_per_exon_matrix(tools_gr_list, exons_gr)
  
  if (!is.null(gt_bed) && file.exists(gt_bed)) {
    cat("==> 4. Assigning ground truth labels...\n")
    per_exon_dt <- assign_labels_from_gt(per_exon_dt, gt_bed, exons_gr)
    cat(sprintf("    %d positive exons / %d total (%.2f%%)\n",
                sum(per_exon_dt$label),
                nrow(per_exon_dt),
                100 * mean(per_exon_dt$label)))
  } else {
    cat("==> 4. No ground truth found: only features will be created (prediction mode)\n")
  }
  
  cat("==> 5. Writing outputs...\n")
  
  out_matrix <- file.path(output_dir, paste0(sample_id, "_per_exon_matrix.tsv"))
  fwrite(per_exon_dt, out_matrix, sep = "\t")
  cat(sprintf("    per_exon_matrix : %s (%d x %d)\n",
              out_matrix, nrow(per_exon_dt), ncol(per_exon_dt)))
  
  out_consensus     <- file.path(output_dir, paste0(sample_id, "_merged_target_consensus.tsv"))
  out_consensus_bed <- file.path(output_dir, paste0(sample_id, "_merged_target_consensus.bed"))
  make_merged_target_consensus(per_exon_dt, tool_names = names(tools_gr_list),
                               sample_id = sample_id, out_file = out_consensus,
                               bed_file = out_consensus_bed, min_callers = min_callers)
  cat(sprintf("    merged_target   : %s\n", out_consensus))
  cat(sprintf("    merged_target   : %s\n", out_consensus_bed))

  out_smooth     <- file.path(output_dir, paste0(sample_id, "_smoothed_variants.tsv"))
  out_smooth_bed <- file.path(output_dir, paste0(sample_id, "_smoothed_variants.bed"))
  make_smoothed_variants(per_exon_dt, sample_id = sample_id,
                         out_file = out_smooth, bed_file = out_smooth_bed,
                         min_callers = min_callers, max_gap_bp = smooth_gap_bp)
  cat(sprintf("    smoothed_vars   : %s\n", out_smooth))
  cat(sprintf("    smoothed_vars   : %s\n", out_smooth_bed))

  invisible(per_exon_dt)
}

make_merged_target_consensus <- function(per_exon_dt, tool_names, sample_id,
                                         out_file, bed_file = NULL, min_callers = 2) {
  called_cols  <- intersect(paste0(tool_names, "_called"),  names(per_exon_dt))
  active_tools <- sub("_called$", "", called_cols)

  dt <- per_exon_dt[n_tools_called >= min_callers & end > start]
  if (nrow(dt) == 0) {
    fwrite(data.table(), out_file, sep = "\t")
    if (!is.null(bed_file))
      fwrite(data.table(CHR = character(), START = integer(), END = integer(), TYPE = character()),
             bed_file, sep = "\t", col.names = FALSE)
    cat(sprintf("       0 exons with >= %d callers\n", min_callers))
    return(invisible(NULL))
  }

  for (tool in active_tools) {
    col_c  <- paste0(tool, "_called")
    col_d  <- paste0(tool, "_direction")
    col_e  <- paste0(tool, "_cn_estimate")
    col_cs <- paste0(tool, "_call_start")   
    col_ce <- paste0(tool, "_call_end")
    col_cc <- paste0(tool, "_call_chr")
    
    dt[, (tool) := ifelse(get(col_c)==1 & get(col_d)==-1, "DEL",
                          ifelse(get(col_c)==1 & get(col_d)== 1, "DUP", 0))]
    
    has_native <- all(c(col_cc, col_cs, col_ce) %in% names(dt))
    dt[, paste0(tool,"_coords") := ifelse(
      get(col_c) == 1,
      if (has_native)
        paste0(get(col_cc), ":", get(col_cs), "-", get(col_ce))
      else
        paste0(chr, ":", start, "-", end),
      0
    )]
    
    dt[, paste0(tool,"_CN") :=
         ifelse(get(col_c)==1 & col_e %in% names(dt), get(col_e), 0)]
  }
  
  dt[, consensus_type := ifelse(n_tools_DEL > n_tools_DUP, "DEL",
                                ifelse(n_tools_DUP > n_tools_DEL, "DUP", "CONFLICT"))]
  
  out <- data.table(
    sample      = sample_id,
    target_id   = paste0(dt$chr,":",dt$start,"-",dt$end),
    CHROM       = dt$chr,
    START       = dt$start,
    END         = dt$end,
    target_name = dt$exon_id
  )
  for (tool in active_tools)                   out[, (tool)                    := dt[[tool]]]
  for (tool in active_tools)                   out[, paste0(tool,"_coords")    := dt[[paste0(tool,"_coords")]]]
  for (tool in active_tools)                   out[, paste0(tool,"_CN")        := dt[[paste0(tool,"_CN")]]]
  out[, Count_Detected := dt$n_tools_called]
  out[, n_DEL          := dt$n_tools_DEL]
  out[, n_DUP          := dt$n_tools_DUP]
  out[, consensus_type := dt$consensus_type]
  
  fwrite(out, out_file, sep = "\t")
  if (!is.null(bed_file))
    fwrite(out[, .(CHR = CHROM, START, END, TYPE = consensus_type)], bed_file, sep = "\t", col.names = FALSE)
  cat(sprintf("       %d exons written in merged_target_consensus\n", nrow(out)))
}

make_smoothed_variants <- function(per_exon_dt, sample_id, out_file,
                                   bed_file = NULL, min_callers = 2, max_gap_bp = 500000) {
  dt <- per_exon_dt[n_tools_called >= min_callers &
                      direction_consensus != 0 &
                      has_conflict == 0 &
                      end > start]
  if (nrow(dt) == 0) {
    fwrite(data.table(), out_file, sep = "\t")
    if (!is.null(bed_file))
      fwrite(data.table(CHR = character(), START = integer(), END = integer(), TYPE = character()),
             bed_file, sep = "\t", col.names = FALSE)
    cat(sprintf("       0 CNVs in smoothed_variants with current filters\n"))
    return(invisible(NULL))
  }
  dt <- dt[order(chr, start)]
  called_cols  <- grep("_called$", names(dt), value = TRUE)
  active_tools <- sub("_called$", "", called_cols)
  
  cnvs    <- list()
  current <- NULL
  for (i in seq_len(nrow(dt))) {
    row   <- dt[i]
    dir   <- row$direction_consensus
    nc    <- as.integer(row$n_tools_called)
    callers_i <- paste(active_tools[sapply(active_tools, function(t) {
      col <- paste0(t, "_called")
      col %in% names(row) && row[[col]] == 1
    })], collapse = ";")
    
    if (!is.null(current) &&
        current$chr == row$chr &&
        current$dir == dir &&
        (row$start - current$end) <= max_gap_bp) {
      current$end    <- as.integer(row$end)
      current$n_targets <- current$n_targets + 1L
      current$n_callers <- max(current$n_callers, nc)
      merged_set <- union(strsplit(current$callers, ";")[[1]],
                          strsplit(callers_i, ";")[[1]])
      current$callers <- paste(sort(merged_set), collapse = ";")
      current$target_names <- paste(current$target_names, row$exon_id, sep = ";")
    } else {
      if (!is.null(current)) cnvs <- c(cnvs, list(current))
      current <- list(chr = as.character(row$chr), start = as.integer(row$start),
                      end = as.integer(row$end), dir = dir,
                      type = ifelse(dir==1,"DUP","DEL"),
                      n_targets = 1L, n_callers = nc,
                      callers = callers_i,
                      target_names = as.character(row$exon_id))
    }
  }
  if (!is.null(current)) cnvs <- c(cnvs, list(current))
  
  out <- rbindlist(lapply(cnvs, function(x) data.table(
    sample = sample_id, CHROM = x$chr, consensus_type = x$type,
    START = x$start, END = x$end, n_targets = x$n_targets,
    n_callers = x$n_callers, callers = x$callers, target_names = x$target_names
  )))
  fwrite(out, out_file, sep = "\t")
  if (!is.null(bed_file))
    fwrite(out[, .(CHR = CHROM, START, END, TYPE = consensus_type)], bed_file, sep = "\t", col.names = FALSE)
  cat(sprintf("       %d CNVs in smoothed_variants\n", nrow(out)))
}

option_list <- list(
  make_option("--input_dir",       type = "character", help = "Directory with per-tool TSV/VCF files"),
  make_option("--capture_bed",     type = "character", help = "Capture regions BED file"),
  make_option("--output_dir",      type = "character", help = "Output directory"),
  make_option("--sample_id",       type = "character", default = "sample1",  help = "Sample identifier [default: %default]"),
  make_option("--gt_bed",          type = "character", default = NULL,        help = "Ground-truth BED (optional)"),
  make_option("--jabcontool_file", type = "character", default = NULL,        help = "JabConTool TSV file (optional)"),
  make_option("--gatk_vcf_file",   type = "character", default = NULL,        help = "GATK gCNV VCF/VCF.GZ file (optional)"),
  make_option("--min_callers",     type = "integer",   default = 2L,          help = "Minimum callers required [default: %default]"),
  make_option("--smooth_gap_bp",   type = "integer",   default = 500000L,     help = "Maximum gap (bp) for smoothing [default: %default]"),
  make_option("--max_call_size",   type = "integer",   default = 500000L,     help = "Maximum call size (bp) [default: %default]"),
  make_option("--gtf_file",        type = "character", default = NULL,        help = "GTF annotation file (optional, activates CDS mode)")
)

opt <- parse_args(OptionParser(option_list = option_list))

for (required in c("input_dir", "capture_bed", "output_dir")) {
  if (is.null(opt[[required]])) {
    stop(sprintf("Required flag --%s is missing.", required))
  }
}

run_all(
  input_dir       = opt$input_dir,
  capture_bed     = opt$capture_bed,
  output_dir      = opt$output_dir,
  sample_id       = opt$sample_id,
  gt_bed          = opt$gt_bed,
  jabcontool_file = opt$jabcontool_file,
  gatk_vcf_file   = opt$gatk_vcf_file,
  min_callers     = opt$min_callers,
  smooth_gap_bp   = opt$smooth_gap_bp,
  max_call_size   = opt$max_call_size,
  gtf_file        = opt$gtf_file
)
