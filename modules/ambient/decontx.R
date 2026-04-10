#!/usr/bin/env Rscript
# decontx.R
#
# Loads the stacked S+U+A H5 from BARCODE_ESTIMATION and the knee CSV,
# calls cells / empties, runs decontX, and saves:
#   - filt_counts_<sampleId>.h5              (decontX-corrected cell matrix, 10x v3)
#   - cell_barcodes_<sampleId>.csv           (called cell barcodes, no header)
#   - dcx_params_<sampleId>.csv             (contamination + cluster per barcode)
#   - barcodes_qc_metrics_<sampleId>.csv.gz (pre_S/U/A and post_S/U/A per barcode)
#   - dcx_summary_<sampleId>.csv            (per-sample summary for the report)
#   - ambient_params_<sampleId>.env         (KEY=VALUE for Nextflow env() outputs)

suppressPackageStartupMessages({
  library(magrittr)
  library(rhdf5)
  library(DropletUtils)
  library(data.table)
  library(Matrix)
  library(decontX)
  library(strex)
})

# ---------------------------------------------------------------------------
# Read H5 directly with rhdf5 — mirrors scprocess .get_h5_mx(sel_s = '')
# Avoids DropletUtils read10xCounts appending "-1" to barcodes and
# DelayedArray colname-assignment failures on large matrices.
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

run_decontx <- function() {

  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 3) {
    stop("Usage: decontx.R <sampleId> <h5_file> <knee_csv>")
  }

  sample_id <- args[1]
  h5_file   <- args[2]
  knee_csv  <- args[3]

  message("Sample:   ", sample_id)
  message("H5:       ", h5_file)
  message("Knee CSV: ", knee_csv)

  # -------------------------------------------------------------------------
  # Load H5 (stacked S+U+A matrix, 10x v3, written by barcode_estimation.R)
  # -------------------------------------------------------------------------
  af_mat <- .get_h5_mx(h5_file)
  message("Loaded H5: ", nrow(af_mat), " features x ", ncol(af_mat), " barcodes")

  # -------------------------------------------------------------------------
  # Load knee data
  # -------------------------------------------------------------------------
  knee_dt <- fread(knee_csv)

  expected_cells  <- unique(knee_dt$expected_cells)
  total_included  <- unique(knee_dt$total_droplets_included)

  if (length(expected_cells) != 1 || length(total_included) != 1) {
    stop("expected_cells / total_droplets_included must be unique scalars in knee CSV")
  }
  message("expected_cells:          ", expected_cells)
  message("total_droplets_included: ", total_included)

  # -------------------------------------------------------------------------
  # Call cells and empties
  # -------------------------------------------------------------------------
  # Cells: top expected_cells barcodes by rank (rank 1 = highest UMIs)
  cell_bcs  <- knee_dt[order(rank), barcode][seq_len(expected_cells)]
  cell_bcs  <- cell_bcs[cell_bcs %in% colnames(af_mat)]

  # Empties: barcodes flagged in the empty plateau by barcode_estimation.R
  empty_bcs <- knee_dt[in_empty_plateau == TRUE, barcode]
  empty_bcs <- empty_bcs[empty_bcs %in% colnames(af_mat)]

  message("Cells:   ", length(cell_bcs))
  message("Empties: ", length(empty_bcs))

  if (length(cell_bcs) == 0)  stop("No cell barcodes found in H5 after filtering")
  if (length(empty_bcs) == 0) stop("No empty barcodes found in H5 after filtering")

  cell_mat  <- af_mat[, cell_bcs,  drop = FALSE]
  empty_mat <- af_mat[, empty_bcs, drop = FALSE]

  # -------------------------------------------------------------------------
  # Run decontX
  # -------------------------------------------------------------------------
  message("Running decontX ...")
  dcx_res   <- decontX(x = cell_mat, background = empty_mat)

  clean_mat <- round(dcx_res$decontXcounts)
  keep_idx  <- which(Matrix::colSums(clean_mat) > 0)
  message("Barcodes with non-zero post-decontX counts: ", length(keep_idx))
  clean_mat <- clean_mat[, keep_idx, drop = FALSE]

  # -------------------------------------------------------------------------
  # Helper: sum S / U / A rows by suffix
  # -------------------------------------------------------------------------
  .sum_usa <- function(mat, prefix) {
    suffixes <- c("S", "U", "A")
    sums_ls  <- lapply(suffixes, function(s) {
      idx <- str_detect(rownames(mat), paste0("_", s, "$"))
      as.numeric(Matrix::colSums(mat[idx, , drop = FALSE]))
    })
    dt <- as.data.table(setNames(sums_ls, paste0(prefix, "_", suffixes)))
    dt[, barcode := colnames(mat)]
    setcolorder(dt, c("barcode", paste0(prefix, "_", suffixes)))
    dt
  }

  # -------------------------------------------------------------------------
  # Save filtered H5
  # -------------------------------------------------------------------------
  h5_out <- sprintf("filt_counts_%s.h5", sample_id)
  write10xCounts(h5_out, clean_mat, version = "3", overwrite = TRUE)
  message("Written filtered H5: ", h5_out)

  # -------------------------------------------------------------------------
  # Save cell barcodes CSV (no header — mirrors scprocess bcs_f format)
  # -------------------------------------------------------------------------
  bcs_out <- sprintf("cell_barcodes_%s.csv", sample_id)
  fwrite(data.table(barcode = colnames(clean_mat)), bcs_out,
         col.names = FALSE, quote = FALSE)
  message("Written cell barcodes: ", bcs_out)

  # -------------------------------------------------------------------------
  # Save dcx_params CSV (contamination + cluster per barcode)
  # -------------------------------------------------------------------------
  dcx_out <- sprintf("dcx_params_%s.csv", sample_id)
  dcx_dt  <- data.table(
    barcode           = colnames(clean_mat),
    pct_contamination = dcx_res$contamination[keep_idx],
    dcx_cluster       = sprintf("cl%02d", dcx_res$z[keep_idx])
  )
  fwrite(dcx_dt, dcx_out)
  message("Written dcx params: ", dcx_out)

  # -------------------------------------------------------------------------
  # Save barcodes_qc_metrics CSV.gz
  # pre_S/U/A = counts from original af_mat (all barcodes)
  # post_S/U/A = counts from decontX output (cells only; NA for discarded bcs)
  # -------------------------------------------------------------------------
  pre_dt  <- .sum_usa(af_mat,    "pre")
  post_dt <- .sum_usa(clean_mat, "post")
  qc_dt   <- merge(pre_dt, post_dt, by = "barcode", all.x = TRUE)

  # Keep only barcodes with >= 10 total pre-count (mirrors get_usa_dt)
  qc_dt[, pre_total := pre_S + pre_U + pre_A]
  qc_dt   <- qc_dt[pre_total >= 10][, pre_total := NULL]

  qc_out  <- sprintf("barcodes_qc_metrics_%s.csv.gz", sample_id)
  fwrite(qc_dt, qc_out)
  message("Written QC metrics: ", qc_out)

  # -------------------------------------------------------------------------
  # Save per-sample summary CSV for the report
  # -------------------------------------------------------------------------
  mean_contam <- round(mean(dcx_res$contamination[keep_idx]) * 100, 2)
  summary_dt  <- data.table(
    sample_id               = sample_id,
    n_cells                 = length(keep_idx),
    n_empties               = length(empty_bcs),
    expected_cells          = as.integer(expected_cells),
    total_droplets_included = as.integer(total_included),
    mean_pct_contamination  = mean_contam
  )
  fwrite(summary_dt, sprintf("dcx_summary_%s.csv", sample_id))
  message("Written summary: ", sprintf("dcx_summary_%s.csv", sample_id))

  # -------------------------------------------------------------------------
  # Write ambient_params.env for Nextflow env() outputs
  # -------------------------------------------------------------------------
  env_out <- sprintf("ambient_params_%s.env", sample_id)
  writeLines(c(
    sprintf("DCX_EXPECTED_CELLS=%s",    as.integer(expected_cells)),
    sprintf("DCX_TOTAL_INCLUDED=%s",    as.integer(total_included)),
    sprintf("DCX_MEAN_CONTAMINATION=%s", mean_contam)
  ), env_out)
  message("Written ambient params: ", env_out)
}

run_decontx()
