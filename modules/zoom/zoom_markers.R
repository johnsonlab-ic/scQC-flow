#!/usr/bin/env Rscript

source("annotation_utils.R")

run_zoom_markers <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 9) {
    stop(
      paste(
        "Usage: zoom_markers.R <integration_csv> <genome_gtf> <sel_res>",
        "<min_cl_size> <min_cells> <out_markers> <out_logcpms> <out_marker_expr> <h5_pattern>",
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
  out_marker_expr <- args[8]
  h5_pattern <- args[9]
  n_cores <- if (length(args) >= 10L) as.integer(args[10]) else 1L
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

  cl_col <- paste0("RNA_snn_res.", sel_res)
  int_dt <- fread(integration_f)
  assert_that(cl_col %in% names(int_dt), msg = sprintf("Integration output is missing '%s'", cl_col))
  umap_dt <- merge(
    ann_dt[, .(sample_id, cell_id, cluster)],
    int_dt[, .(cell_id, UMAP1, UMAP2)],
    by = "cell_id",
    all.x = FALSE,
    sort = FALSE
  )

  # Compute logCPMs only for the top markers that will be plotted in the report.
  # Use a generous threshold (FDR < 0.1, top 50 per cluster) to cover any
  # threshold tuning at report render time, avoiding full-gene computation.
  top_for_cpms <- get_top_markers(marker_dt, fdr_cut = 0.1, n_top = 50, max_zero_p = 1.0)
  logcpm_gene_ids <- unique(top_for_cpms$gene_id)
  message("Computing logCPMs for ", length(logcpm_gene_ids), " marker genes (not all ",
          nrow(pb_obj$row_dt), " genes)")
  logcpms_dt <- make_logcpms_for_genes(prep_obj, pb_obj$row_dt, gene_ids = logcpm_gene_ids)

  top_sel_dt <- unique(top_for_cpms[, .(label = as.character(cluster), symbol, gene_id)])
  top_marker_expr_dt <- if (nrow(top_sel_dt) > 0) {
    load_h5_marker_expression_multi(
      h5_files,
      sel_dt_list = list(top = top_sel_dt),
      int_dt = umap_dt
    )$top
  } else {
    data.table()
  }

  fwrite(marker_dt, out_markers)
  fwrite(logcpms_dt, out_logcpms)
  saveRDS(top_marker_expr_dt, out_marker_expr)

  message("Wrote marker statistics: ", out_markers)
  message("Wrote logCPM summaries: ", out_logcpms)
  message("Wrote top-marker expression cache: ", out_marker_expr)
  message("=== ZOOM MARKERS done ===")
}

run_zoom_markers()