#!/usr/bin/env Rscript
# apply_qc.R
#
# Per-sample QC metric calculation and threshold application.
# Reads the decontX H5, GTF, and pre-computed scDblFinder results,
# then calculates per-cell metrics and applies hard + soft thresholds.
#
# Inputs:  H5 + GTF + scdblfinder_<sampleId>.csv.gz + threshold params
# Outputs: qc_metrics_<sampleId>.csv.gz + qc_summary_<sampleId>.csv

suppressPackageStartupMessages({
  library(magrittr)
  library(rhdf5)
  library(Matrix)
  library(data.table)
  library(strex)
})

# ---------------------------------------------------------------------------
# Read H5 (stacked S+U+A)
# ---------------------------------------------------------------------------
.get_h5_mx <- function(af_mat_f) {
  h5   <- H5Fopen(af_mat_f, flags = "H5F_ACC_RDONLY")
  mat  <- sparseMatrix(
    i    = as.vector(h5$matrix$indices + 1L),
    p    = as.vector(h5$matrix$indptr),
    x    = as.vector(h5$matrix$data),
    repr = "C",
    dims = h5$matrix$shape
  )
  colnames(mat) <- as.character(h5$matrix$barcodes)
  rownames(mat) <- as.character(h5$matrix$features$name)
  H5Fclose(h5)
  mat
}

# ---------------------------------------------------------------------------
# Parse GTF — returns ensembl_id / symbol / gene_type table
# ---------------------------------------------------------------------------
.parse_gtf_annotations <- function(gtf_f) {
  decomp_cmd <- if (grepl("\\.gz$", gtf_f)) {
    paste0("zcat ", shQuote(gtf_f))
  } else {
    paste0("cat ", shQuote(gtf_f))
  }
  gtf_lines <- fread(
    cmd       = paste0(decomp_cmd, " | grep -v '^#'"),
    sep       = "\t",
    header    = FALSE,
    select    = c(3, 9),
    col.names = c("feature", "attributes")
  )
  gtf_lines <- gtf_lines[feature == "gene"]

  gtf_lines[, ensembl_id := str_match(attributes, 'gene_id "([^"]+)"')[, 2]]
  gtf_lines[, symbol     := str_match(attributes, 'gene_name "([^"]+)"')[, 2]]
  gtf_lines[, gene_type  := str_match(attributes, 'gene_type "([^"]+)"')[, 2]]

  na_type <- is.na(gtf_lines$gene_type)
  if (any(na_type)) {
    gtf_lines[na_type, gene_type := str_match(attributes, 'gene_biotype "([^"]+)"')[, 2]]
  }

  gtf_lines[, attributes := NULL]
  gtf_lines[, feature    := NULL]
  unique(gtf_lines, by = "ensembl_id")
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
run_apply_qc <- function() {

  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 16) {
    stop("Usage: apply_qc.R <sampleId> <h5_file> <genome_gtf> <dbl_csv> ",
         "<metadata_csv> ",
         "<hard_min_counts> <hard_min_feats> <hard_max_mito> ",
         "<min_counts> <min_feats> <max_mito> <min_mito> ",
         "<max_splice> <min_splice> ",
         "<mad_filter> <mad_nmads>")
  }

  current_sample_id <- args[1]
  h5_file         <- args[2]
  gtf_file        <- args[3]
  dbl_csv         <- args[4]
  metadata_csv    <- args[5]
  hard_min_counts <- as.numeric(args[6])
  hard_min_feats  <- as.numeric(args[7])
  hard_max_mito   <- as.numeric(args[8])
  min_counts      <- as.numeric(args[9])
  min_feats       <- as.numeric(args[10])
  max_mito        <- as.numeric(args[11])
  min_mito        <- as.numeric(args[12])
  max_splice      <- as.numeric(args[13])
  min_splice      <- as.numeric(args[14])
  mad_filter      <- tolower(args[15]) == "true"
  mad_nmads       <- as.numeric(args[16])

  message("=== APPLY_QC: ", current_sample_id, " ===")
  message("Hard thresholds: min_counts=", hard_min_counts,
          " min_feats=", hard_min_feats, " max_mito=", hard_max_mito)
  message("Soft thresholds: counts=[", min_counts, ",Inf]",
          " feats=[", min_feats, ",Inf]",
          " mito=[", min_mito, ",", max_mito, "]",
          " splice=[", min_splice, ",", max_splice, "]")
  message("MAD filter: ", mad_filter, if (mad_filter) paste0(" (nmads=", mad_nmads, ")") else "")

  # -------------------------------------------------------------------------
  # 1. Load filtered H5
  # -------------------------------------------------------------------------
  mat <- .get_h5_mx(h5_file)
  message("Loaded H5: ", nrow(mat), " features x ", ncol(mat), " barcodes")

  # -------------------------------------------------------------------------
  # 2. Split S / U / A and compute totals
  # -------------------------------------------------------------------------
  splice_ns       <- c("S", "U", "A")
  usa_ls          <- lapply(splice_ns, function(s) {
    idx <- grep(paste0("_", s, "$"), rownames(mat))
    mat[idx, , drop = FALSE]
  }) %>% setNames(splice_ns)

  proper_gs       <- sub("_[SUA]$", "", rownames(usa_ls[["S"]]))
  total_spliced   <- Matrix::colSums(usa_ls[["S"]])
  total_unspliced <- Matrix::colSums(usa_ls[["U"]])

  counts_mat      <- Reduce("+", usa_ls)
  rownames(counts_mat) <- proper_gs

  # -------------------------------------------------------------------------
  # 3. Gene annotations — exclude rRNA/Mt_rRNA, identify mitochondrial genes
  # -------------------------------------------------------------------------
  gene_annots  <- .parse_gtf_annotations(gtf_file)
  rrna_ids     <- gene_annots[gene_type == "rRNA",    ensembl_id]
  mt_rrna_ids  <- gene_annots[gene_type == "Mt_rRNA", ensembl_id]
  exclude_gs   <- rownames(counts_mat) %in% c(rrna_ids, mt_rrna_ids)
  message("Excluding ", sum(exclude_gs), " rRNA / Mt_rRNA genes")
  counts_mat   <- counts_mat[!exclude_gs, ]

  mito_symbols <- gene_annots[grepl("^MT-", symbol), ensembl_id]
  mito_gs      <- rownames(counts_mat) %in% mito_symbols
  mito_sum     <- Matrix::colSums(counts_mat[mito_gs, , drop = FALSE])
  message("Mitochondrial genes found: ", sum(mito_gs))

  # -------------------------------------------------------------------------
  # 4. Per-cell QC metrics
  # -------------------------------------------------------------------------
  cell_counts  <- Matrix::colSums(counts_mat)
  cell_feats   <- Matrix::colSums(counts_mat > 0)
  mito_pct     <- mito_sum / cell_counts
  splice_frac  <- (total_spliced + 1) / (total_spliced + total_unspliced + 2)

  qc_dt <- data.table(
    cell_id         = colnames(counts_mat),
    sample_id       = current_sample_id,
    sum             = as.numeric(cell_counts),
    detected        = as.integer(cell_feats),
    mito_sum        = as.numeric(mito_sum),
    mito_pct        = mito_pct,
    total_spliced   = as.numeric(total_spliced),
    total_unspliced = as.numeric(total_unspliced),
    log_counts      = log10(cell_counts + 1),
    log_feats       = log10(cell_feats + 1),
    logit_mito      = qlogis((mito_sum + 1) / (cell_counts + 2)),
    logit_spliced   = qlogis(splice_frac)
  )

  qc_dt  <- qc_dt[sum > 0]
  n_total <- nrow(qc_dt)
  message("Total cells after zero removal: ", n_total)

  # -------------------------------------------------------------------------
  # 5. Read doublet results and merge
  # -------------------------------------------------------------------------
  dbl_dt    <- fread(dbl_csv)
  dbl_merge <- dbl_dt[, .(cell_id, scdbl_class)]
  qc_dt     <- merge(qc_dt, dbl_merge, by = "cell_id", all.x = TRUE)
  qc_dt[is.na(scdbl_class), scdbl_class := "singlet"]
  qc_dt[, is_singlet := scdbl_class == "singlet"]

  # -------------------------------------------------------------------------
  # 6. Attach sample-level metadata if available
  # -------------------------------------------------------------------------
  if (basename(metadata_csv) != "NO_FILE") {
    meta_dt <- fread(metadata_csv)
    if (!"sample_id" %in% names(meta_dt)) {
      stop("Normalized metadata file is missing a sample_id column")
    }
    sample_meta <- meta_dt[meta_dt[["sample_id"]] == current_sample_id]
    if (nrow(sample_meta) != 1) {
      available_ids <- sort(unique(meta_dt[["sample_id"]]))
      stop(
        "Expected exactly one metadata row for sample ", current_sample_id,
        "; found ", nrow(sample_meta),
        ". Available sample_id values include: ",
        paste(utils::head(available_ids, 10), collapse = ", "),
        if (length(available_ids) > 10) " ..." else ""
      )
    }
    qc_dt <- merge(qc_dt, sample_meta, by = "sample_id", all.x = TRUE, sort = FALSE)
  }

  # -------------------------------------------------------------------------
  # 7. Apply hard and soft QC thresholds
  # -------------------------------------------------------------------------
  qc_dt[, keep_hard :=
    is_singlet &
    log_counts >= log10(hard_min_counts + 1) &
    log_feats  >= log10(hard_min_feats + 1) &
    logit_mito <  qlogis(hard_max_mito)
  ]

  qc_dt[, keep_flat :=
    keep_hard &
    log_counts    >= log10(min_counts + 1) &
    log_feats     >= log10(min_feats + 1) &
    logit_mito    >  qlogis(min_mito) &
    logit_mito    <  qlogis(max_mito) &
    logit_spliced >  qlogis(min_splice) &
    logit_spliced <  qlogis(max_splice)
  ]

  n_flat <- sum(qc_dt$keep_flat)
  message("Cells passing flat thresholds: ", n_flat, " / ", n_total,
          " (", round(100 * n_flat / n_total, 1), "%)")

  # -------------------------------------------------------------------------
  # 7b. MAD adaptive filtering — applied to flat-passing cells only
  # -------------------------------------------------------------------------
  .mad_bounds <- function(x, nmads) {
    m <- median(x, na.rm = TRUE)
    d <- mad(x, na.rm = TRUE)
    c(lo = m - nmads * d, hi = m + nmads * d)
  }

  if (mad_filter) {
    flat_pass <- qc_dt[keep_flat == TRUE]

    b_counts <- .mad_bounds(flat_pass$log_counts,    mad_nmads)
    b_feats  <- .mad_bounds(flat_pass$log_feats,     mad_nmads)
    b_mito   <- .mad_bounds(flat_pass$logit_mito,    mad_nmads)
    b_splice <- .mad_bounds(flat_pass$logit_spliced, mad_nmads)

    qc_dt[, mad_lo_log_counts    := b_counts["lo"]]
    qc_dt[, mad_hi_log_counts    := b_counts["hi"]]
    qc_dt[, mad_lo_log_feats     := b_feats["lo"]]
    qc_dt[, mad_hi_log_feats     := b_feats["hi"]]
    qc_dt[, mad_lo_logit_mito    := b_mito["lo"]]
    qc_dt[, mad_hi_logit_mito    := b_mito["hi"]]
    qc_dt[, mad_lo_logit_spliced := b_splice["lo"]]
    qc_dt[, mad_hi_logit_spliced := b_splice["hi"]]

    qc_dt[, keep :=
      keep_flat &
      log_counts    >= b_counts["lo"] & log_counts    <= b_counts["hi"] &
      log_feats     >= b_feats["lo"]  & log_feats     <= b_feats["hi"] &
      logit_mito    >= b_mito["lo"]   & logit_mito    <= b_mito["hi"] &
      logit_spliced >= b_splice["lo"] & logit_spliced <= b_splice["hi"]
    ]

    n_mad_rm_counts <- sum(flat_pass$log_counts    < b_counts["lo"] | flat_pass$log_counts    > b_counts["hi"])
    n_mad_rm_feats  <- sum(flat_pass$log_feats     < b_feats["lo"]  | flat_pass$log_feats     > b_feats["hi"])
    n_mad_rm_mito   <- sum(flat_pass$logit_mito    < b_mito["lo"]   | flat_pass$logit_mito    > b_mito["hi"])
    n_mad_rm_splice <- sum(flat_pass$logit_spliced < b_splice["lo"] | flat_pass$logit_spliced > b_splice["hi"])

    n_kept <- sum(qc_dt$keep)
    message("Cells passing MAD filter: ", n_kept, " / ", n_flat,
            " (removed ", n_flat - n_kept, ")")
    message("  MAD removals: counts=", n_mad_rm_counts, " feats=", n_mad_rm_feats,
            " mito=", n_mad_rm_mito, " splice=", n_mad_rm_splice)
  } else {
    qc_dt[, keep := keep_flat]
    qc_dt[, c("mad_lo_log_counts", "mad_hi_log_counts",
               "mad_lo_log_feats",  "mad_hi_log_feats",
               "mad_lo_logit_mito", "mad_hi_logit_mito",
               "mad_lo_logit_spliced", "mad_hi_logit_spliced") := NA_real_]
    n_kept          <- n_flat
    n_mad_rm_counts <- NA_integer_
    n_mad_rm_feats  <- NA_integer_
    n_mad_rm_mito   <- NA_integer_
    n_mad_rm_splice <- NA_integer_
  }

  message("Cells passing final QC: ", n_kept, " / ", n_total,
          " (", round(100 * n_kept / n_total, 1), "%)")

  # -------------------------------------------------------------------------
  # 8. Write outputs
  # -------------------------------------------------------------------------
  qc_out <- sprintf("qc_metrics_%s.csv.gz", current_sample_id)
  fwrite(qc_dt, qc_out)
  message("Written QC metrics: ", qc_out)

  summary_dt <- data.table(
    sample_id       = current_sample_id,
    n_total         = n_total,
    n_singlets      = sum(qc_dt$is_singlet),
    n_doublets      = sum(!qc_dt$is_singlet),
    n_pass_hard     = sum(qc_dt$keep_hard),
    n_pass_flat     = n_flat,
    n_pass_qc       = n_kept,
    n_excluded      = n_total - n_kept,
    pct_excluded    = round(100 * (1 - n_kept / n_total), 1),
    mad_filter      = mad_filter,
    mad_nmads       = if (mad_filter) mad_nmads else NA_real_,
    n_mad_removed   = if (mad_filter) n_flat - n_kept else NA_integer_,
    n_mad_rm_counts = n_mad_rm_counts,
    n_mad_rm_feats  = n_mad_rm_feats,
    n_mad_rm_mito   = n_mad_rm_mito,
    n_mad_rm_splice = n_mad_rm_splice,
    median_counts   = round(median(qc_dt[keep == TRUE, sum]), 0),
    median_feats    = round(median(qc_dt[keep == TRUE, detected]), 0),
    median_mito_pct = round(median(qc_dt[keep == TRUE, mito_pct]) * 100, 2)
  )
  fwrite(summary_dt, sprintf("qc_summary_%s.csv", current_sample_id))
  message("Written summary: ", sprintf("qc_summary_%s.csv", current_sample_id))

  message("=== APPLY_QC done ===")
}

run_apply_qc()
