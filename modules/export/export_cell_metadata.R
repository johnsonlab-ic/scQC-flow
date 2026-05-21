#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

parse_args <- function(args) {
  if (length(args) %% 2 != 0) {
    stop("Arguments must be provided as --key value pairs")
  }
  out <- list()
  idx <- seq(1, length(args), by = 2)
  for (i in idx) {
    key <- sub("^--", "", args[[i]])
    out[[key]] <- args[[i + 1]]
  }
  out
}

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  required <- c("qc_pattern", "integration_csv", "annotation_pattern", "utils_r", "output_all_csv", "output_clean_csv")
  missing <- setdiff(required, names(args))
  if (length(missing) > 0) {
    stop(sprintf("Missing required arguments: %s", paste(missing, collapse = ", ")))
  }

  source(args$utils_r)

  all_dir <- dirname(args$output_all_csv)
  clean_dir <- dirname(args$output_clean_csv)
  dir.create(all_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(clean_dir, recursive = TRUE, showWarnings = FALSE)

  metadata_dt <- load_export_metadata(
    args$integration_csv,
    args$qc_pattern,
    args$annotation_pattern
  )

  fwrite(metadata_dt, args$output_all_csv)
  message(sprintf("Wrote unified cell metadata: %s (%d rows)", args$output_all_csv, nrow(metadata_dt)))

  has_doublet_cols <- all(c("is_dbl", "in_dbl_cl") %in% names(metadata_dt))
  if (!has_doublet_cols) {
    message("Skipping clean metadata export: doublet columns not found in merged metadata")
    return(invisible(NULL))
  }

  clean_dt <- metadata_dt[is_dbl == FALSE & in_dbl_cl == FALSE]
  fwrite(clean_dt, args$output_clean_csv)
  message(sprintf("Wrote clean unified cell metadata: %s (%d rows)", args$output_clean_csv, nrow(clean_dt)))
}

main()