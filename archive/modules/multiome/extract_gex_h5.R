#!/usr/bin/env Rscript
# extract_gex_h5.R
# Extract Gene Expression modality from multiome H5 and write as single-modality H5
# This prepares the raw counts for CellBender which expects single-modality data

suppressPackageStartupMessages({
    library(argparse)
    library(Seurat)
    library(DropletUtils)
})

# Parse command line arguments
parser <- ArgumentParser(description = "Extract Gene Expression from multiome H5")
parser$add_argument("--input_h5", type = "character", required = TRUE,
                    help = "Input multiome raw H5 file")
parser$add_argument("--output_h5", type = "character", required = TRUE,
                    help = "Output Gene Expression only H5 file")
parser$add_argument("--sample_name", type = "character", required = TRUE,
                    help = "Sample name for logging")

args <- parser$parse_args()

cat("=== Extracting Gene Expression from Multiome H5 ===\n")
cat("Sample:", args$sample_name, "\n")
cat("Input H5:", args$input_h5, "\n")
cat("Output H5:", args$output_h5, "\n\n")

# Read the multiome H5 - returns a list with Gene Expression and Peaks
cat("Reading multiome H5 file...\n")
counts <- Read10X_h5(args$input_h5)

# Check if it's a list (multiome) or single matrix
if (is.list(counts)) {
    cat("Detected multiome data with modalities:", paste(names(counts), collapse = ", "), "\n")
    
    # Extract Gene Expression modality
    if ("Gene Expression" %in% names(counts)) {
        gex_counts <- counts$`Gene Expression`
        cat("Extracted Gene Expression modality\n")
    } else {
        stop("No 'Gene Expression' modality found in H5 file")
    }
} else {
    cat("Single modality data detected, using as-is\n")
    gex_counts <- counts
}

cat("Gene Expression matrix dimensions:", nrow(gex_counts), "genes x", ncol(gex_counts), "barcodes\n")

# Write out as 10X-compatible H5 using DropletUtils
cat("Writing Gene Expression H5 file...\n")
write10xCounts(
    path = args$output_h5,
    x = gex_counts,
    type = "HDF5",
    genome = "GRCh38",
    version = "3",
    overwrite = TRUE
)

cat("Successfully wrote Gene Expression H5 to:", args$output_h5, "\n")
cat("=== Extraction complete ===\n")
