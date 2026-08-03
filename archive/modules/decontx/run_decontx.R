#!/usr/bin/env Rscript
# run_decontx.R
#
# Ambient RNA removal using decontX, mirroring scprocess/scripts/ambient.R.
#
# Cell calling follows scprocess's barcodeRanks method:
#   - Cells  : top expected_cells barcodes by rank (from knee_data)
#   - Empties: barcodes flagged in_empty_plateau == TRUE (from knee_data)
#   No emptyDrops is run — it is too slow on unfiltered matrices.
#
# Inputs:
#   sample_id     — sample identifier
#   h5_file       — unfiltered total (S+U+A) count matrix from BARCODE_ESTIMATION
#   knee_data_csv — per-barcode CSV from BARCODE_ESTIMATION:
#                   columns: barcode, rank, total, expected_cells, in_empty_plateau, ...
#   ncores        — parallelism (currently unused; kept for interface compatibility)
#   h5_out        — output path for ambient-corrected cell matrix
#   barcodes_out  — output path for retained barcode list
#   contamination_out — output path for per-barcode contamination metrics
#
# Usage:
#   Rscript run_decontx.R <sample_id> <h5_file> <knee_data_csv> <ncores>
#                         <h5_out> <barcodes_out> <contamination_out>

suppressPackageStartupMessages({
    library(DropletUtils)
    library(rhdf5)
    library(decontX)
    library(Matrix)
    library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 8) {
    stop(paste(
        "Usage: run_decontx.R <sample_id> <h5_file> <knee_data_csv> <ncores>",
        "<h5_out> <barcodes_out> <contamination_out> <usa_metrics_out>"
    ))
}

sample_id         <- args[1]
h5_file           <- args[2]
knee_data_csv     <- args[3]
ncores            <- as.integer(args[4])
h5_out            <- args[5]
barcodes_out      <- args[6]
contamination_out <- args[7]
usa_metrics_out   <- args[8]

message("=== DECONTX ===")
message("Sample:        ", sample_id)
message("H5 input:      ", h5_file)
message("Knee data:     ", knee_data_csv)
message("Cores:         ", ncores)

# ---------------------------------------------------------------------------
# 1. Load knee data (produced by BARCODE_ESTIMATION, sorted by rank)
# ---------------------------------------------------------------------------
message("Loading knee data ...")
knee_dt <- fread(knee_data_csv)
setorder(knee_dt, rank)

expected_cells <- knee_dt$expected_cells[1]
message("Expected cells (from knee):  ", expected_cells)
message("In-empty-plateau barcodes:   ", sum(knee_dt$in_empty_plateau))

# ---------------------------------------------------------------------------
# 2. Cell calling — barcodeRanks method (mirrors scprocess)
#    Cells  = top expected_cells barcodes (already sorted by rank)
#    Empties = barcodes in the empty plateau (in_empty_plateau == TRUE)
# ---------------------------------------------------------------------------
cell_bcs  <- knee_dt$barcode[seq_len(expected_cells)]
empty_bcs <- knee_dt[in_empty_plateau == TRUE, barcode]

message("Cells called:  ", length(cell_bcs))
message("Empty droplets:", length(empty_bcs))

if (length(cell_bcs) == 0) {
    stop("No cells called — check knee_data_csv and expected_cells.")
}
if (length(empty_bcs) == 0) {
    stop("No empty-plateau barcodes found — check in_empty_plateau column in knee_data_csv.")
}

# ---------------------------------------------------------------------------
# 3. Load unfiltered count matrix and extract cell / empty submatrices
# ---------------------------------------------------------------------------
message("Loading count matrix ...")
h5_obj <- H5Fopen(h5_file, flags = "H5F_ACC_RDONLY")
mat <- sparseMatrix(
    i    = as.vector(h5_obj$matrix$indices + 1),
    p    = as.vector(h5_obj$matrix$indptr),
    x    = as.vector(h5_obj$matrix$data),
    repr = "C",
    dims = h5_obj$matrix$shape
)
colnames(mat) <- h5_obj$matrix$barcodes
rownames(mat) <- as.character(h5_obj$matrix$features$name)
H5Fclose(h5_obj)
message("Barcodes in H5: ", ncol(mat))

# Intersect with barcodes actually present in the H5
cell_bcs_present  <- intersect(cell_bcs,  colnames(mat))
empty_bcs_present <- intersect(empty_bcs, colnames(mat))

if (length(cell_bcs_present) == 0) {
    stop("None of the cell barcodes from knee_data are present in the H5 matrix.")
}

message("Cell barcodes matched in H5:  ", length(cell_bcs_present))
message("Empty barcodes matched in H5: ", length(empty_bcs_present))

cell_mat  <- mat[, cell_bcs_present,  drop = FALSE]
empty_mat <- mat[, empty_bcs_present, drop = FALSE]

# ---------------------------------------------------------------------------
# 4. decontX — ambient correction using the empty-plateau as background
#    varGenes = 2000 matches scprocess default
# ---------------------------------------------------------------------------
message("Running decontX ...")
dcx_res   <- decontX(x = cell_mat, background = empty_mat, varGenes = 2000)

clean_mat <- round(dcx_res$decontXcounts)
keep_idx  <- which(Matrix::colSums(clean_mat) != 0)
message("Keeping ", length(keep_idx), " barcodes with non-zero counts in decontaminated matrix")
clean_mat <- clean_mat[, keep_idx]

# ---------------------------------------------------------------------------
# 5. Write outputs
# ---------------------------------------------------------------------------
message("Writing filtered H5: ", h5_out)
write10xCounts(h5_out, clean_mat, version = "3", overwrite = TRUE)

message("Writing barcode list: ", barcodes_out)
cell_bcs_dt <- data.table(barcode = colnames(clean_mat))
fwrite(cell_bcs_dt, barcodes_out, col.names = FALSE, quote = FALSE)

message("Writing contamination metrics: ", contamination_out)
dcx_clust  <- factor(sprintf("cl%02d", dcx_res$z))
dcx_contam <- dcx_res$contamination

contam_dt <- data.table(
    barcode           = colnames(clean_mat),
    pct.contamination = dcx_contam[keep_idx],
    dcx_cluster       = dcx_clust[keep_idx]
)
fwrite(contam_dt, contamination_out)

message("Median contamination: ", round(median(contam_dt$pct.contamination) * 100, 2), "%")

# ---------------------------------------------------------------------------
# 6. Write per-barcode total-count metrics (pre + post decontX)
#    Covers barcodes up to total_droplets_included — the same set the ambient
#    report plots (cells + empties used as background + excluded barcodes).
#    post_total is NA for barcodes not retained as cells.
#
#    NOTE: the H5 from BARCODE_ESTIMATION is already a collapsed S+U+A total
#    (gene names have no -S/-U/-A suffix).  We therefore store pre_total and
#    post_total (plain colSums) rather than separate S/U/A columns.
# ---------------------------------------------------------------------------
message("Writing USA metrics: ", usa_metrics_out)

total_droplets_included <- knee_dt$total_droplets_included[1]

# Pre-decontX: top total_droplets_included barcodes by rank
top_bcs         <- knee_dt$barcode[seq_len(min(total_droplets_included, nrow(knee_dt)))]
top_bcs_present <- intersect(top_bcs, colnames(mat))
pre_dt <- data.table(
    barcode   = top_bcs_present,
    pre_total = as.numeric(Matrix::colSums(mat[, top_bcs_present, drop = FALSE]))
)

# Post-decontX: only retained cell barcodes (NA elsewhere)
post_dt <- data.table(
    barcode    = colnames(clean_mat),
    post_total = as.numeric(Matrix::colSums(clean_mat))
)

# Merge: all top barcodes, post_total NA for non-cells
usa_dt <- merge(pre_dt, post_dt, by = "barcode", all.x = TRUE)
fwrite(usa_dt, usa_metrics_out)

message("USA metrics rows: ", nrow(usa_dt))
message("=== DONE ===")
