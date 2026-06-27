#!/usr/bin/env Rscript
# emptydrops_calling.R
#
# CellSweep-recommended cell calling: identify cells with DropletUtils::emptyDrops
# (no GMM, no splice logic). Emits the same downstream contract as CELL_CALLING so
# QC / HVG / integration / AMBIENT_DE / CellSweep run unchanged.
#
#   is_cell  = emptyDrops FDR <= FDR_CUT  (real cells, profile deviates from ambient)
#   is_empty = total UMI < LOWER          (the ambient pool emptyDrops tests against)
#   (lower <= UMI, not significant => 'ambiguous', neither cell nor empty)
#
# Outputs (Nextflow cwd):
#   filt_counts_<id>.h5             stacked S+U+A matrix subset to is_cell
#   cell_barcodes_<id>.csv          is_cell barcodes (no header)
#   empty_barcodes_<id>.csv         is_empty barcodes (ambient profile input)
#   emptydrops_labels_<id>.csv.gz   per-candidate FDR + population (UMI >= LOWER)
#   emptydrops_summary_<id>.csv     per-sample counts + params
#   cell_calling_<id>.env           KEY=VALUE for Nextflow env() outputs

suppressPackageStartupMessages({
  library(rhdf5); library(DropletUtils); library(data.table); library(Matrix)
})
LOWER   <- 100      # emptyDrops ambient threshold (also the is_empty cutoff)
FDR_CUT <- 0.01

a <- commandArgs(trailingOnly = TRUE)        # <sampleId> <h5_file>
sid <- a[1]; h5_file <- a[2]
message("EmptyDrops calling: ", sid)

.get_h5_mx <- function(f) {
  h <- H5Fopen(f, flags = "H5F_ACC_RDONLY")
  m <- sparseMatrix(i = as.vector(h$matrix$indices + 1L), p = as.vector(h$matrix$indptr),
                    x = as.vector(h$matrix$data), repr = "C", dims = h$matrix$shape)
  colnames(m) <- as.character(h$matrix$barcodes)
  rownames(m) <- as.character(h$matrix$features$name)
  H5Fclose(h); m
}

af   <- .get_h5_mx(h5_file)                  # features (stacked S/U/A) x barcodes
feat <- rownames(af); base <- sub("_(S|U|A)$", "", feat); genes <- unique(base)
G    <- sparseMatrix(i = match(base, genes), j = seq_along(base), x = 1,
                     dims = c(length(genes), length(base)))
summed  <- G %*% af; rownames(summed) <- genes   # gene-level counts for emptyDrops
total   <- Matrix::colSums(summed)
spliced <- Matrix::colSums(af[grepl("_S$", feat), , drop = FALSE])

set.seed(100)
ed <- emptyDrops(summed, lower = LOWER)
is_cell <- !is.na(ed$FDR) & ed$FDR <= FDR_CUT
bc <- colnames(af)
dt <- data.table(barcode = bc, total = as.numeric(total), spliced = as.numeric(spliced),
                 splice_pct = as.numeric(spliced / pmax(total, 1)),
                 FDR = ed$FDR, is_cell = is_cell, is_empty = (total < LOWER))
dt[, population := fifelse(is_cell, "cell", fifelse(total < LOWER, "empty", "ambiguous"))]
message(sprintf("  cells=%d  empty=%d  ambiguous=%d",
                sum(dt$is_cell), sum(dt$is_empty), sum(dt$population == "ambiguous")))

cell_bcs <- dt[is_cell == TRUE, barcode]
if (length(cell_bcs) == 0) stop("emptyDrops called no cells")
write10xCounts(sprintf("filt_counts_%s.h5", sid), af[, cell_bcs, drop = FALSE],
               version = "3", overwrite = TRUE)
fwrite(data.table(barcode = cell_bcs), sprintf("cell_barcodes_%s.csv", sid),
       col.names = FALSE, quote = FALSE)
fwrite(data.table(barcode = dt[is_empty == TRUE, barcode]),
       sprintf("empty_barcodes_%s.csv", sid), col.names = FALSE, quote = FALSE)
fwrite(dt[total >= LOWER, .(barcode, total, spliced, splice_pct, FDR, population, is_cell, is_empty)],
       sprintf("emptydrops_labels_%s.csv.gz", sid))

fwrite(data.table(sample_id = sid, n_total = nrow(dt), n_cell = sum(dt$is_cell),
                  n_empty = sum(dt$is_empty), n_ambiguous = sum(dt$population == "ambiguous"),
                  n_is_cell = sum(dt$is_cell), n_is_empty = sum(dt$is_empty),
                  lower = LOWER, fdr_cut = FDR_CUT),
       sprintf("emptydrops_summary_%s.csv", sid))
writeLines(c(sprintf("ED_N_CELL=%d", sum(dt$is_cell)),
             sprintf("ED_N_EMPTY=%d", sum(dt$is_empty))),
           sprintf("cell_calling_%s.env", sid))
message("Done: ", sid)
