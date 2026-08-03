#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(argparse)
  library(data.table)
})

parser <- ArgumentParser(description = "Format annotation metadata for export")
parser$add_argument("--input_csv", type = "character", required = TRUE,
  help = "Input annotation cell-label CSV.GZ")
parser$add_argument("--prefix", type = "character", required = TRUE,
  help = "Export column prefix")
parser$add_argument("--output_csv", type = "character", required = TRUE,
  help = "Output export metadata CSV.GZ")

opts <- parser$parse_args()

prefix <- gsub("[^A-Za-z0-9_]", "_", opts$prefix)
dt <- fread(opts$input_csv)
out_dt <- data.table(cell_id = dt$cell_id)
out_dt[[paste0(prefix, "_annotation")]] <- dt$label
out_dt[[paste0(prefix, "_annotation_score")]] <- dt$label_score
if ("n_markers" %in% names(dt)) {
  out_dt[[paste0(prefix, "_annotation_n_markers")]] <- dt$n_markers
}

fwrite(out_dt, opts$output_csv)