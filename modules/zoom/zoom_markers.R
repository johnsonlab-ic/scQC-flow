#!/usr/bin/env Rscript

source("annotation_utils.R")

run_zoom_markers <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 8) {
    stop(
      paste(
        "Usage: zoom_markers.R <integration_csv> <genome_gtf> <sel_res>",
        "<min_cl_size> <min_cells> <out_markers> <out_logcpms> <h5_pattern>",
        "[n_cores]"
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
  n_cores <- if (length(args) >= 9L) as.integer(args[9]) else 1L
  if (is.na(n_cores) || n_cores < 1L) n_cores <- 1L

  h5_files <- Sys.glob(h5_pattern)
  assert_that(length(h5_files) > 0, msg = sprintf("No filtered H5 files matched '%s'", h5_pattern))

  message("=== ZOOM MARKERS ===")
  message("Selected resolution: ", sel_res)
  message("Filtered H5 files: ", length(h5_files))
  message("Workers: ", n_cores)

  biotypes_dt <- parse_gtf_annotations(genome_gtf)
  ann_dt <- load_annotation_cells(integration_f, sel_res, min_cl_size)
  pb_obj <- build_pseudobulk_from_h5s(h5_files, ann_dt, biotypes_dt, n_cores = n_cores)
  prep_obj <- prepare_cluster_matrices(pb_obj, min_cells = min_cells)
  marker_dt <- calc_find_markers_pseudobulk(prep_obj, pb_obj$row_dt)

  # Compute logCPMs only for the top markers that will be plotted in the report.
  # Use a generous threshold (FDR < 0.1, top 50 per cluster) to cover any
  # threshold tuning at report render time, avoiding full-gene computation.
  top_for_cpms <- get_top_markers(marker_dt, fdr_cut = 0.1, n_top = 50, max_zero_p = 1.0)
  logcpm_gene_ids <- unique(top_for_cpms$gene_id)
  message("Computing logCPMs for ", length(logcpm_gene_ids), " marker genes (not all ",
          nrow(pb_obj$row_dt), " genes)")
  logcpms_dt <- make_logcpms_for_genes(prep_obj, pb_obj$row_dt, gene_ids = logcpm_gene_ids)

  fwrite(marker_dt, out_markers)
  fwrite(logcpms_dt, out_logcpms)

  message("Wrote marker statistics: ", out_markers)
  message("Wrote logCPM summaries: ", out_logcpms)
  message("=== ZOOM MARKERS done ===")
}

run_zoom_markers()