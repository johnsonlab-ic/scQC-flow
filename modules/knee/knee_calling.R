#!/usr/bin/env Rscript
# knee_calling.R
#
# "Classic" knee/shin cell calling — the barcode-rank selection that originally
# fed decontX (decontX itself is dropped on this branch). Mirrors the archived
# decontX cell selection exactly:
#   is_cell  = top `expected_cells` barcodes by rank   (from knee_data)
#   is_empty = in_empty_plateau == TRUE                (from knee_data)
#   (remaining barcodes => 'ambiguous', neither cell nor empty)
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
message("  expected_cells (knee): ", expected_cells,
        " | in_empty_plateau: ", sum(knee$in_empty_plateau))

knee[, is_cell  := rank <= expected_cells]
knee[, is_empty := in_empty_plateau == TRUE]
knee[, population := fifelse(is_cell, "cell", fifelse(is_empty, "empty", "ambiguous"))]
knee[, splice_pct := as.numeric(spliced) / pmax(as.numeric(total), 1)]

cell_bcs  <- knee[is_cell  == TRUE, barcode]
empty_bcs <- knee[is_empty == TRUE, barcode]
if (length(cell_bcs)  == 0) stop("No cells called (check expected_cells)")
if (length(empty_bcs) == 0) stop("No empty-plateau barcodes (check in_empty_plateau)")
message(sprintf("  cells=%d  empty=%d  ambiguous=%d",
                length(cell_bcs), length(empty_bcs),
                sum(knee$population == "ambiguous")))

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
                  n_ambiguous = sum(knee$population == "ambiguous"),
                  expected_cells = expected_cells,
                  knee1 = knee$knee1[1], shin1 = knee$shin1[1],
                  knee2 = knee$knee2[1], shin2 = knee$shin2[1]),
       sprintf("knee_summary_%s.csv", sid))
writeLines(c(sprintf("KNEE_N_CELL=%d", length(cell_bcs)),
             sprintf("KNEE_N_EMPTY=%d", length(empty_bcs))),
           sprintf("cell_calling_%s.env", sid))
message("Done: ", sid)
