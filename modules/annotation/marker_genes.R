#!/usr/bin/env Rscript

source("annotation_utils.R")

run_annotation_markers <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 10) {
    stop(
      paste(
        "Usage: marker_genes.R <integration_csv> <genome_gtf> <marker_csv>",
        "<sel_res> <min_cl_size> <min_cells>",
        "<out_markers> <out_logcpms> <out_panel> <h5_pattern>"
      )
    )
  }

  integration_f <- args[1]
  genome_gtf <- args[2]
  marker_csv <- args[3]
  sel_res <- args[4]
  min_cl_size <- as.integer(args[5])
  min_cells <- as.integer(args[6])
  out_markers <- args[7]
  out_logcpms <- args[8]
  out_panel <- args[9]
  h5_pattern <- args[10]

  h5_files <- Sys.glob(h5_pattern)
  assert_that(length(h5_files) > 0, msg = sprintf("No filtered H5 files matched '%s'", h5_pattern))

  message("=== ANNOTATION ===")
  message("Selected resolution: ", sel_res)
  message("Filtered H5 files: ", length(h5_files))

  biotypes_dt <- parse_gtf_annotations(genome_gtf)
  ann_dt <- load_annotation_cells(integration_f, sel_res, min_cl_size)
  pb_obj <- build_pseudobulk_from_h5s(h5_files, ann_dt, biotypes_dt)
  logcpms_dt <- make_logcpms_all(pb_obj, min_cells = min_cells)
  marker_dt <- calc_find_markers_pseudobulk(logcpms_dt, pb_obj$row_dt)
  panel_dt <- normalize_marker_panel(marker_csv, biotypes_dt, marker_dt)

  fwrite(marker_dt, out_markers)
  fwrite(logcpms_dt, out_logcpms)
  fwrite(panel_dt, out_panel)

  message("Wrote marker statistics: ", out_markers)
  message("Wrote logCPM summaries: ", out_logcpms)
  message("Wrote processed marker panel: ", out_panel)
  message("=== ANNOTATION done ===")
}

run_annotation_markers()