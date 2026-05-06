#!/usr/bin/env Rscript

source("annotation_utils.R")

run_prep_report_inputs <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 4) {
    stop(
      paste(
        "Usage: prep_report_inputs.R <integration_csv> <marker_stats_csv>",
        "<out_top_marker_expr> <h5_pattern>"
      )
    )
  }

  integration_f <- args[1]
  marker_stats_f <- args[2]
  out_top_marker_expr <- args[3]
  h5_pattern <- args[4]

  sel_res <- Sys.getenv("ANNOTATION_SEL_RES", "0.2")
  min_cpm_mkr <- as.numeric(Sys.getenv("ANNOTATION_MIN_CPM_MKR", "50"))
  not_ok_re <- Sys.getenv("ANNOTATION_NOT_OK_RE", "(lincRNA|lncRNA|pseudogene|antisense)")
  top_n <- as.integer(Sys.getenv("ANNOTATION_TOP_N", "10"))
  fdr_cut <- as.numeric(Sys.getenv("ANNOTATION_FDR_CUT", "0.01"))
  max_zero_p <- as.numeric(Sys.getenv("ANNOTATION_MAX_ZERO_P", "0.5"))

  h5_files <- Sys.glob(h5_pattern)
  if (length(h5_files) == 0) {
    saveRDS(data.table(), out_top_marker_expr)
    message("No filtered H5 inputs matched pattern; wrote empty top-marker cache.")
    return(invisible(NULL))
  }

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
  if (nrow(umap_dt) == 0) {
    saveRDS(data.table(), out_top_marker_expr)
    message("No cells available at selected annotation resolution; wrote empty top-marker cache.")
    return(invisible(NULL))
  }

  mkrs_dt <- fread(marker_stats_f)
  top_min_dt <- mkrs_dt[str_detect(gene_type, not_ok_re, negate = TRUE)][
    logcpm.sel >= log(min_cpm_mkr + 1)
  ] |> get_top_markers(fdr_cut = fdr_cut, n_top = top_n, max_zero_p = max_zero_p)

  top_sel_dt <- unique(top_min_dt[, .(label = as.character(cluster), symbol, gene_id)])
  if (nrow(top_sel_dt) == 0) {
    saveRDS(data.table(), out_top_marker_expr)
    message("No top markers passed filters; wrote empty top-marker cache.")
    return(invisible(NULL))
  }

  top_marker_expr_dt <- load_h5_marker_expression_multi(
    h5_files,
    sel_dt_list = list(top = top_sel_dt),
    int_dt = umap_dt
  )$top

  saveRDS(top_marker_expr_dt, out_top_marker_expr)
  message("Wrote top-marker expression cache: ", out_top_marker_expr)
}

run_prep_report_inputs()
