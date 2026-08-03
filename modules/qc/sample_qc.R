#!/usr/bin/env Rscript
# sample_qc.R
#
# Per-sample cell-level QC mirroring scprocess's SampleQC.R.
#
# Takes the decontX-filtered H5 (stacked S+U+A format) and genome GTF,
# calculates per-cell QC metrics, runs scDblFinder, applies hard + soft
# thresholds, and saves:
#   - qc_metrics_<sampleId>.csv.gz     (per-cell QC with keep/exclude flags)
#   - qc_summary_<sampleId>.csv        (summary stats for the report)
#   - scdblfinder_<sampleId>.csv.gz    (doublet detection results)

suppressPackageStartupMessages({
  library(magrittr)
  library(rhdf5)
  library(Matrix)
  library(data.table)
  library(SingleCellExperiment)
  library(scDblFinder)
  library(strex)
})

# ---------------------------------------------------------------------------
# Read H5 (stacked S+U+A) — mirrors scprocess .get_h5_mx
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
# Parse GTF to get gene annotations (ensembl_id, symbol, gene_type)
# ---------------------------------------------------------------------------
.parse_gtf_annotations <- function(gtf_f) {
  # read only gene-level entries (skip comment lines that start with #)
  # handle both plain and gzipped GTFs
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

  # extract fields from attributes column
  gtf_lines[, ensembl_id := str_match(attributes, 'gene_id "([^"]+)"')[, 2]]
  gtf_lines[, symbol     := str_match(attributes, 'gene_name "([^"]+)"')[, 2]]
  gtf_lines[, gene_type  := str_match(attributes, 'gene_type "([^"]+)"')[, 2]]

  # fall back: if gene_type missing, try gene_biotype
  na_type <- is.na(gtf_lines$gene_type)
  if (any(na_type)) {
    gtf_lines[na_type, gene_type := str_match(attributes, 'gene_biotype "([^"]+)"')[, 2]]
  }

  gtf_lines[, attributes := NULL]
  gtf_lines[, feature    := NULL]

  # de-duplicate (some GTFs have overlapping entries)
  gtf_lines <- unique(gtf_lines, by = "ensembl_id")

  return(gtf_lines)
}

# ---------------------------------------------------------------------------
# Main QC function
# ---------------------------------------------------------------------------
run_sample_qc <- function() {

  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 12) {
    stop("Usage: sample_qc.R <sampleId> <h5_file> <genome_gtf> ",
         "<hard_min_counts> <hard_min_feats> <hard_max_mito> ",
         "<min_counts> <min_feats> <max_mito> <min_mito> ",
         "<max_splice> <min_splice>")
  }

  sample_id        <- args[1]
  h5_file          <- args[2]
  gtf_file         <- args[3]
  hard_min_counts  <- as.numeric(args[4])
  hard_min_feats   <- as.numeric(args[5])
  hard_max_mito    <- as.numeric(args[6])
  min_counts       <- as.numeric(args[7])
  min_feats        <- as.numeric(args[8])
  max_mito         <- as.numeric(args[9])
  min_mito         <- as.numeric(args[10])
  max_splice       <- as.numeric(args[11])
  min_splice       <- as.numeric(args[12])

  message("=== SAMPLE_QC: ", sample_id, " ===")
  message("H5:              ", h5_file)
  message("GTF:             ", gtf_file)
  message("Hard thresholds: min_counts=", hard_min_counts,
          " min_feats=", hard_min_feats, " max_mito=", hard_max_mito)
  message("Soft thresholds: counts=[", min_counts, ",Inf]",
          " feats=[", min_feats, ",Inf]",
          " mito=[", min_mito, ",", max_mito, "]",
          " splice=[", min_splice, ",", max_splice, "]")

  # -------------------------------------------------------------------------
  # 1. Load filtered H5 (stacked S+U+A from decontX)
  # -------------------------------------------------------------------------
  mat <- .get_h5_mx(h5_file)
  message("Loaded H5: ", nrow(mat), " features x ", ncol(mat), " barcodes")

  # -------------------------------------------------------------------------
  # 2. Split S / U / A and compute totals
  # -------------------------------------------------------------------------
  splice_ns  <- c("S", "U", "A")
  usa_ls     <- lapply(splice_ns, function(s) {
    idx <- grep(paste0("_", s, "$"), rownames(mat))
    mat[idx, , drop = FALSE]
  }) %>% setNames(splice_ns)

  # proper gene names (strip _S/_U/_A suffix)
  proper_gs  <- sub("_[SUA]$", "", rownames(usa_ls[["S"]]))

  # per-cell spliced / unspliced totals
  total_spliced   <- Matrix::colSums(usa_ls[["S"]])
  total_unspliced <- Matrix::colSums(usa_ls[["U"]])

  # combined counts matrix (S + U + A)
  counts_mat <- Reduce("+", usa_ls)
  rownames(counts_mat) <- proper_gs

  # -------------------------------------------------------------------------
  # 3. Gene annotations — identify rRNA and mitochondrial genes
  # -------------------------------------------------------------------------
  gene_annots <- .parse_gtf_annotations(gtf_file)

  # rRNA and Mt_rRNA genes to exclude
  rrna_ids    <- gene_annots[gene_type == "rRNA", ensembl_id]
  mt_rrna_ids <- gene_annots[gene_type == "Mt_rRNA", ensembl_id]
  exclude_gs  <- rownames(counts_mat) %in% c(rrna_ids, mt_rrna_ids)
  message("Excluding ", sum(exclude_gs), " rRNA / Mt_rRNA genes")
  counts_mat  <- counts_mat[!exclude_gs, ]

  # mitochondrial genes (by symbol prefix MT-)
  mito_symbols <- gene_annots[grepl("^MT-", symbol), ensembl_id]
  mito_gs      <- rownames(counts_mat) %in% mito_symbols
  mito_sum     <- Matrix::colSums(counts_mat[mito_gs, , drop = FALSE])
  message("Mitochondrial genes found: ", sum(mito_gs))

  # -------------------------------------------------------------------------
  # 4. Per-cell QC metrics (mirrors scprocess make_qc_dt)
  # -------------------------------------------------------------------------
  cell_counts  <- Matrix::colSums(counts_mat)
  cell_feats   <- Matrix::colSums(counts_mat > 0)
  mito_pct     <- mito_sum / cell_counts
  splice_frac  <- (total_spliced + 1) / (total_spliced + total_unspliced + 2)

  qc_dt <- data.table(
    cell_id        = colnames(counts_mat),
    sample_id      = sample_id,
    sum            = as.numeric(cell_counts),
    detected       = as.integer(cell_feats),
    mito_sum       = as.numeric(mito_sum),
    mito_pct       = mito_pct,
    total_spliced  = as.numeric(total_spliced),
    total_unspliced = as.numeric(total_unspliced),
    log_counts     = log10(cell_counts + 1),
    log_feats      = log10(cell_feats + 1),
  logit_mito     = qlogis((mito_sum + 1) / (cell_counts + 2)),
    logit_spliced  = qlogis(splice_frac)
  )

  # exclude any zero-count cells
  qc_dt <- qc_dt[sum > 0]
  n_total <- nrow(qc_dt)
  message("Total cells after zero removal: ", n_total)

  # -------------------------------------------------------------------------
  # 5. Doublet detection with scDblFinder
  #    Use combined S+U+A counts (mirrors scprocess SampleQC.R which passes
  #    the full counts_mat to scDblFinder)
  # -------------------------------------------------------------------------
  message("Running scDblFinder on combined S+U+A counts (",
          nrow(counts_mat), " genes) ...")
  keep_for_dbl <- qc_dt$detected >= 100
  n_dbl_input  <- sum(keep_for_dbl)

  if (n_dbl_input >= 100) {
    dbl_cells   <- qc_dt[keep_for_dbl, cell_id]
    dbl_mat     <- counts_mat[, dbl_cells, drop = FALSE]
    dbl_sce     <- SingleCellExperiment(assays = list(counts = dbl_mat))
    dbl_res     <- scDblFinder(dbl_sce, returnType = "table",
                               multiSampleMode = "singleModel", verbose = FALSE)
    dbl_dt      <- as.data.table(dbl_res, keep.rownames = "cell_id")
    dbl_dt      <- dbl_dt[type == "real"]
    dbl_dt[, sample_id := sample_id]
    setnames(dbl_dt, "class", "scdbl_class")
    message("scDblFinder: ", sum(dbl_dt$scdbl_class == "doublet"), " doublets / ",
            nrow(dbl_dt), " cells")
  } else {
    message("Too few cells (", n_dbl_input, ") for scDblFinder — skipping doublet detection")
    dbl_dt <- data.table(
      cell_id     = qc_dt$cell_id,
      sample_id   = sample_id,
      scdbl_class = "singlet",
      scdbl_score = NA_real_
    )
  }

  # save doublet results
  dbl_out <- sprintf("scdblfinder_%s.csv.gz", sample_id)
  fwrite(dbl_dt, dbl_out)
  message("Written doublet results: ", dbl_out)

  # -------------------------------------------------------------------------
  # 6. Merge doublet info and apply QC thresholds
  # -------------------------------------------------------------------------
  # join doublet class to qc_dt
  dbl_merge <- dbl_dt[, .(cell_id, scdbl_class)]
  qc_dt     <- merge(qc_dt, dbl_merge, by = "cell_id", all.x = TRUE)
  qc_dt[is.na(scdbl_class), scdbl_class := "singlet"]

  # restrict to singlets for QC
  qc_dt[, is_singlet := scdbl_class == "singlet"]

  # hard baseline thresholds (very permissive — removes obvious junk)
  qc_dt[, keep_hard :=
    is_singlet &
    log_counts >= log10(hard_min_counts + 1) &
    log_feats  >= log10(hard_min_feats + 1) &
    logit_mito <  qlogis(hard_max_mito)
  ]

  # soft per-sample thresholds (more stringent)
  qc_dt[, keep :=
    keep_hard &
    log_counts    >= log10(min_counts + 1) &
    log_feats     >= log10(min_feats + 1) &
    logit_mito    >  qlogis(min_mito) &
    logit_mito    <  qlogis(max_mito) &
    logit_spliced >  qlogis(min_splice) &
    logit_spliced <  qlogis(max_splice)
  ]

  n_kept <- sum(qc_dt$keep)
  message("Cells passing QC: ", n_kept, " / ", n_total,
          " (", round(100 * n_kept / n_total, 1), "%)")

  # save qc metrics
  qc_out <- sprintf("qc_metrics_%s.csv.gz", sample_id)
  fwrite(qc_dt, qc_out)
  message("Written QC metrics: ", qc_out)

  # -------------------------------------------------------------------------
  # 7. Summary CSV for the report
  # -------------------------------------------------------------------------
  summary_dt <- data.table(
    sample_id        = sample_id,
    n_total          = n_total,
    n_singlets       = sum(qc_dt$is_singlet),
    n_doublets       = sum(!qc_dt$is_singlet),
    n_pass_hard      = sum(qc_dt$keep_hard),
    n_pass_qc        = n_kept,
    n_excluded       = n_total - n_kept,
    pct_excluded     = round(100 * (1 - n_kept / n_total), 1),
    median_counts    = round(median(qc_dt[keep == TRUE, sum]), 0),
    median_feats     = round(median(qc_dt[keep == TRUE, detected]), 0),
    median_mito_pct  = round(median(qc_dt[keep == TRUE, mito_pct]) * 100, 2)
  )
  fwrite(summary_dt, sprintf("qc_summary_%s.csv", sample_id))
  message("Written summary: ", sprintf("qc_summary_%s.csv", sample_id))

  message("=== SAMPLE_QC done ===")
}

run_sample_qc()
