#!/usr/bin/env Rscript

source("annotation_utils.R")

run_zoom_markers <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 8) {
    stop(
      paste(
        "Usage: zoom_markers.R <integration_csv> <genome_gtf> <sel_res>",
        "<min_cl_size> <min_cells> <out_markers> <out_logcpms> <h5_pattern>"
      )
    )
  }

  integration_f <- args[1]
  genome_gtf <- args[2]
  sel_res <- args[3]
  min_cl_size <- as.integer(args[4])
  min_cells <- as.integer(args[5])
  out_markers <- args[6]
  out_logcpms <- args[7]
  h5_pattern <- args[8]

  h5_files <- Sys.glob(h5_pattern)
  assert_that(length(h5_files) > 0, msg = sprintf("No filtered H5 files matched '%s'", h5_pattern))

  message("=== ZOOM MARKERS ===")
  message("Selected resolution: ", sel_res)
  message("Filtered H5 files: ", length(h5_files))

  biotypes_dt <- parse_gtf_annotations(genome_gtf)
  ann_dt <- load_annotation_cells(integration_f, sel_res, min_cl_size)
  pb_obj <- build_pseudobulk_from_h5s(h5_files, ann_dt, biotypes_dt)
  prep_obj <- prepare_cluster_matrices(pb_obj, min_cells = min_cells)
  logcpms_dt <- make_logcpms_for_genes(prep_obj, pb_obj$row_dt, gene_ids = pb_obj$row_dt$gene_id)
  marker_dt <- calc_find_markers_pseudobulk(prep_obj, pb_obj$row_dt)

  fwrite(marker_dt, out_markers)
  fwrite(logcpms_dt, out_logcpms)

  message("Wrote marker statistics: ", out_markers)
  message("Wrote logCPM summaries: ", out_logcpms)
  message("=== ZOOM MARKERS done ===")
}

run_zoom_markers()