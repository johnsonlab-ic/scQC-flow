#!/usr/bin/env Rscript
# knee_calling.R
#
# "Classic" knee/shin cell calling — replicates the original (pre-CellSweep)
# decontX-era run exactly, minus the decontX() correction call itself. Mirrors
# scprocess's own selection exactly (dev/_archive/scprocess/scripts/ambient.R
# `call_cells_and_empties()`, and pseudobulk_and_empties.R `.get_one_empty_pb()`
# — both use `in_empty_plateau` as the actual ambient/background barcode set;
# `total_droplets_included` is NEVER used for selection anywhere in scprocess,
# only as a cosmetic bound in one report's dot-colouring):
#   cell     = top `expected_cells` barcodes by rank   (from knee_data)
#   ambient  = in_empty_plateau == TRUE                (from knee_data;
#              this is the CellSweep background)
#   excluded = everything else (never touched by anything downstream)
#
# Emits the same downstream contract as the GMM / EmptyDrops callers so
# QC / HVG / integration / annotation / AMBIENT_DE / CellSweep run unchanged.
#
# Usage: knee_calling.R <sampleId> <h5_file> <knee_data_csv>

suppressPackageStartupMessages({
  library(rhdf5); library(DropletUtils); library(data.table); library(Matrix)
})

a <- commandArgs(trailingOnly = TRUE)        # <sampleId> <h5_file> <knee_csv>
sid <- a[1]; h5_file <- a[2]; knee_csv <- a[3]
message("Knee/shin calling: ", sid)

.get_h5_mx <- function(f) {
  h <- H5Fopen(f, flags = "H5F_ACC_RDONLY")
  m <- sparseMatrix(i = as.vector(h$matrix$indices + 1L), p = as.vector(h$matrix$indptr),
                    x = as.vector(h$matrix$data), repr = "C", dims = h$matrix$shape)
  colnames(m) <- as.character(h$matrix$barcodes)
  rownames(m) <- as.character(h$matrix$features$name)
  H5Fclose(h); m
}

knee <- fread(knee_csv)
setorder(knee, rank)
expected_cells <- as.integer(knee$expected_cells[1])
total_included <- as.integer(knee$total_droplets_included[1])
message("  expected_cells (knee): ", expected_cells,
        " | in_empty_plateau: ", sum(knee$in_empty_plateau))

# Selected by rank ORDER, not a `rank <= N` threshold: barcodeRanks() uses
# averaged/fractional ranks for tied UMI totals, so a threshold comparison
# can select a different count than `order(rank)[seq_len(N)]` (this is what
# both decontx.R and scprocess's ambient.R used).
cell_bcs_top <- knee$barcode[seq_len(expected_cells)]
knee[, is_cell  := barcode %in% cell_bcs_top]
knee[, is_empty := in_empty_plateau == TRUE]
knee[, population := fifelse(is_cell, "cell", fifelse(is_empty, "ambient", "excluded"))]
knee[, splice_pct := as.numeric(spliced) / pmax(as.numeric(total), 1)]

cell_bcs  <- knee[is_cell  == TRUE, barcode]
empty_bcs <- knee[is_empty == TRUE, barcode]
if (length(cell_bcs)  == 0) stop("No cells called (check expected_cells)")
if (length(empty_bcs) == 0) stop("No empty-plateau barcodes (check in_empty_plateau)")
message(sprintf("  cells=%d  ambient=%d  excluded=%d",
                length(cell_bcs), length(empty_bcs),
                sum(knee$population == "excluded")))

af  <- .get_h5_mx(h5_file)                   # stacked S/U/A features x barcodes
cell_bcs <- cell_bcs[cell_bcs %in% colnames(af)]
write10xCounts(sprintf("filt_counts_%s.h5", sid), af[, cell_bcs, drop = FALSE],
               version = "3", overwrite = TRUE)
fwrite(data.table(barcode = cell_bcs),  sprintf("cell_barcodes_%s.csv", sid),
       col.names = FALSE, quote = FALSE)
fwrite(data.table(barcode = empty_bcs), sprintf("empty_barcodes_%s.csv", sid),
       col.names = FALSE, quote = FALSE)
fwrite(knee[, .(barcode, rank, total, spliced, splice_pct, population, is_cell, is_empty)],
       sprintf("knee_labels_%s.csv.gz", sid))
fwrite(data.table(sample_id = sid, n_total = nrow(knee),
                  n_cell = length(cell_bcs), n_empty = length(empty_bcs),
                  n_excluded = sum(knee$population == "excluded"),
                  expected_cells = expected_cells,
                  total_droplets_included = total_included,
                  knee1 = knee$knee1[1], shin1 = knee$shin1[1],
                  knee2 = knee$knee2[1], shin2 = knee$shin2[1]),
       sprintf("knee_summary_%s.csv", sid))
writeLines(c(sprintf("KNEE_N_CELL=%d", length(cell_bcs)),
             sprintf("KNEE_N_EMPTY=%d", length(empty_bcs))),
           sprintf("cell_calling_%s.env", sid))
message("Done: ", sid)
