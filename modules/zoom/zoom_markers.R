#!/usr/bin/env Rscript

source("annotation_utils.R")

run_zoom_markers <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 9) {
    stop(
      paste(
        "Usage: zoom_markers.R <integration_csv> <genome_gtf> <sel_res>",
        "<min_cl_size> <min_cells> <out_markers> <out_logcpms> <out_marker_expr> <h5_pattern>",
        "[n_cores] [marker_top_n] [marker_min_cpm] [marker_fdr_cut] [marker_max_zero_p] [marker_not_ok_re]"
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
  marker_top_n <- if (length(args) >= 11L) as.integer(args[11]) else 10L
  if (is.na(marker_top_n) || marker_top_n < 1L) marker_top_n <- 10L
  marker_min_cpm <- if (length(args) >= 12L) as.numeric(args[12]) else 50
  if (is.na(marker_min_cpm) || marker_min_cpm < 0) marker_min_cpm <- 50
  marker_fdr_cut <- if (length(args) >= 13L) as.numeric(args[13]) else 0.01
  if (is.na(marker_fdr_cut) || marker_fdr_cut < 0) marker_fdr_cut <- 0.01
  marker_max_zero_p <- if (length(args) >= 14L) as.numeric(args[14]) else 0.5
  if (is.na(marker_max_zero_p) || marker_max_zero_p < 0 || marker_max_zero_p > 1) marker_max_zero_p <- 0.5
  marker_not_ok_re <- if (length(args) >= 15L) args[15] else "(lincRNA|lncRNA|pseudogene|antisense)"
  marker_top_n_plot <- min(marker_top_n, 5L)

  h5_files <- Sys.glob(h5_pattern)
  assert_that(length(h5_files) > 0, msg = sprintf("No filtered H5 files matched '%s'", h5_pattern))

  message("=== ZOOM MARKERS ===")
  message("Selected resolution: ", sel_res)
  message("Filtered H5 files: ", length(h5_files))
  message("Workers: ", n_cores)
  message("Marker cache top-n: ", marker_top_n_plot)
  message("Marker cache min CPM: ", marker_min_cpm)
  message("Marker cache FDR cut: ", marker_fdr_cut)
  message("Marker cache max zero proportion: ", marker_max_zero_p)

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

  # Compute logCPMs and expression only for the exact marker subset the report
  # will facet over, so report-side plotting cannot request genes absent from
  # the cached zoom marker tables.
  marker_cache_dt <- marker_dt
  if (nzchar(marker_not_ok_re)) {
    marker_cache_dt <- marker_cache_dt[grepl(marker_not_ok_re, gene_type, ignore.case = TRUE, perl = TRUE, invert = TRUE)]
  }
  marker_cache_dt <- marker_cache_dt[logcpm.sel >= log(marker_min_cpm + 1)]
  top_for_cpms <- get_top_markers(
    marker_cache_dt,
    fdr_cut = marker_fdr_cut,
    n_top = marker_top_n_plot,
    max_zero_p = marker_max_zero_p
  )
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