#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

option_list <- list(
  make_option("--input_csv", type = "character", help = "Input annotation cell-label CSV.GZ"),
  make_option("--prefix", type = "character", help = "Export column prefix"),
  make_option("--output_csv", type = "character", help = "Output export metadata CSV.GZ")
)

parser <- OptionParser(option_list = option_list)
opts <- parse_args(parser)

required <- c("input_csv", "prefix", "output_csv")
missing <- required[vapply(required, function(key) is.null(opts[[key]]) || !nzchar(opts[[key]]), logical(1))]
if (length(missing) > 0) {
  print_help(parser)
  stop(sprintf("Missing required arguments: %s", paste(missing, collapse = ", ")), call. = FALSE)
}

prefix <- gsub("[^A-Za-z0-9_]", "_", opts$prefix)
dt <- fread(opts$input_csv)
out_dt <- data.table(cell_id = dt$cell_id)
out_dt[[paste0(prefix, "_annotation")]] <- dt$label
out_dt[[paste0(prefix, "_annotation_score")]] <- dt$label_score
if ("n_markers" %in% names(dt)) {
  out_dt[[paste0(prefix, "_annotation_n_markers")]] <- dt$n_markers
}

fwrite(out_dt, opts$output_csv)