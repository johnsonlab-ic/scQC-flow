#!/usr/bin/env Rscript
# alevinfry_to_h5.R
# Convert simpleaf/alevin-fry quantification output to a 10x-format h5 file.
#
# Usage: Rscript alevinfry_to_h5.R <sample_id> <quant_dir> <h5_out>
#
#   quant_dir  — the {sampleId}_af_quant directory produced by SIMPLEAF_QUANT.
#                Contains af_quant/ (what loadFry reads) and af_map/.
#   h5_out     — output .h5 path (10x v3 format, unfiltered droplets).
#
# The output h5 contains total counts (S + U + A), all non-zero barcodes.
# This is compatible with CellBender, Seurat, and other downstream tools.

suppressPackageStartupMessages({
    library(fishpond)
    library(SingleCellExperiment)
    library(DropletUtils)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
    stop("Usage: alevinfry_to_h5.R <sample_id> <quant_dir> <h5_out>")
}

sample_id <- args[1]
quant_dir <- args[2]
h5_out    <- args[3]

fry_dir <- file.path(quant_dir, "af_quant")
if (!dir.exists(fry_dir)) {
    stop(paste("af_quant directory not found inside quant_dir:", fry_dir))
}

message("Sample:  ", sample_id)
message("fry_dir: ", fry_dir)

# outputFormat = "scRNA" loads the spliced + unspliced + ambiguous (S+U+A)
# counts summed into a single 'counts' assay — the total count matrix.
sce <- loadFry(fry_dir, outputFormat = "scRNA")
mat <- counts(sce)

message("Barcodes before zero-count filter: ", ncol(mat))
mat <- mat[, colSums(mat) > 0]
message("Barcodes after zero-count filter:  ", ncol(mat))
message("Genes: ", nrow(mat))

message("Writing 10x-format h5: ", h5_out)
write10xCounts(h5_out, mat, version = "3", overwrite = TRUE)
message("Done.")
