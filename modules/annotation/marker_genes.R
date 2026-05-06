#!/usr/bin/env Rscript

source("annotation_utils.R")

run_annotation_markers <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 14) {
    stop(
      paste(
        "Usage: marker_genes.R <integration_csv> <genome_gtf> <marker_csv>",
        "<sel_res> <min_cl_size> <min_cells>",
        "<out_markers> <out_logcpms> <out_panel> <out_marker_expr> <out_cell_labels>",
        "<out_pseudobulk> <n_cores> <h5_pattern>"
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
  out_marker_expr <- args[10]
  out_cell_labels <- args[11]
  out_pseudobulk <- args[12]
  n_cores <- as.integer(args[13])
  h5_pattern <- args[14]

  if (is.na(n_cores) || n_cores < 1L) {
    n_cores <- 1L
  }

  h5_files <- Sys.glob(h5_pattern)
  assert_that(length(h5_files) > 0, msg = sprintf("No filtered H5 files matched '%s'", h5_pattern))

  message("=== ANNOTATION ===")
  message("Selected resolution: ", sel_res)
  message("Filtered H5 files: ", length(h5_files))
  message("Workers: ", n_cores)

  biotypes_dt <- parse_gtf_annotations(genome_gtf)
  ann_dt <- load_annotation_cells(integration_f, sel_res, min_cl_size)

  if (file.exists(out_pseudobulk)) {
    message("Loading pseudobulk checkpoint: ", out_pseudobulk)
    pb_obj <- readRDS(out_pseudobulk)
  } else {
    pb_obj <- build_pseudobulk_from_h5s(h5_files, ann_dt, biotypes_dt, n_cores = n_cores)
    saveRDS(pb_obj, out_pseudobulk, compress = FALSE)
    message("Wrote pseudobulk checkpoint: ", out_pseudobulk)
  }

  prep_obj <- prepare_cluster_matrices(pb_obj, min_cells = min_cells)
  marker_dt <- calc_find_markers_pseudobulk(prep_obj, pb_obj$row_dt)
  panel_dt <- normalize_marker_panel(marker_csv, biotypes_dt, marker_dt)

  cluster_label_dt <- assign_cluster_labels_from_panel(marker_dt, panel_dt, min_cpm = 50)

  cl_col <- paste0("RNA_snn_res.", sel_res)
  int_dt <- fread(integration_f)
  assert_that(cl_col %in% names(int_dt), msg = sprintf("Integration output is missing '%s'", cl_col))
  umap_dt <- int_dt[!is.na(get(cl_col)) & sample_id != "", .(
    sample_id,
    cell_id,
    cluster = get(cl_col),
    UMAP1,
    UMAP2
  )]
  cl_ord <- order_cluster_labels(unique(umap_dt$cluster))
  umap_dt[, cluster := factor(as.character(cluster), levels = cl_ord)]
  cluster_label_dt <- merge(data.table(cluster = cl_ord), cluster_label_dt, by = "cluster", all.x = TRUE, sort = FALSE)
  cluster_label_dt[is.na(label), `:=`(label = "Unassigned", label_score = NA_real_, n_markers = 0L)]
  cell_labels_dt <- merge(umap_dt[, .(sample_id, cell_id, cluster)], cluster_label_dt[, .(cluster, label, label_score, n_markers)], by = "cluster", all.x = TRUE, sort = FALSE)

  top_mkrs_dt <- marker_dt[
    logFC > 0 &
    logcpm.sel >= log(50 + 1) &
    !grepl("(lincRNA|lncRNA|pseudogene|antisense)", gene_type, perl = TRUE)
  ] |> get_top_markers(fdr_cut = 0.01, n_top = 10, max_zero_p = 0.5)

  cpms_gene_ids <- unique(c(panel_dt$gene_id, top_mkrs_dt$gene_id))
  logcpms_dt <- make_logcpms_for_genes(prep_obj, pb_obj$row_dt, gene_ids = cpms_gene_ids)

  expr_ls <- load_h5_marker_expression_multi(
    h5_files,
    sel_dt_list = list(panel = panel_dt),
    int_dt = umap_dt
  )
  marker_expr_dt <- expr_ls$panel

  fwrite(marker_dt, out_markers)
  fwrite(logcpms_dt, out_logcpms)
  fwrite(panel_dt, out_panel)
  saveRDS(marker_expr_dt, out_marker_expr)
  fwrite(cell_labels_dt, out_cell_labels)

  message("Wrote marker statistics: ", out_markers)
  message("Wrote logCPM summaries: ", out_logcpms)
  message("Wrote processed marker panel: ", out_panel)
  message("Wrote marker-expression cache: ", out_marker_expr)
  message("Wrote per-cell annotation labels: ", out_cell_labels)
  message("=== ANNOTATION done ===")
}

run_annotation_markers()