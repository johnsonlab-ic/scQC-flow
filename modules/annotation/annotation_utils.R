suppressPackageStartupMessages({
  library(assertthat)
  library(circlize)
  library(ComplexHeatmap)
  library(data.table)
  library(edgeR)
  library(ggplot2)
  library(ggrastr)
  library(grid)
  library(Matrix)
  library(patchwork)
  library(rhdf5)
  library(scales)
  library(scater)
  library(stringr)
})

nice_cols <- c(
  '#DC050C', '#FB8072', '#1965B0', '#7BAFDE', '#882E72', '#B17BA6',
  '#FF7F00', '#FDB462', '#E7298A', '#E78AC3', '#33A02C', '#B2DF8A'
)

parse_gtf_annotations <- function(gtf_f) {
  decomp_cmd <- if (grepl("\\.gz$", gtf_f)) {
    paste0("zcat ", shQuote(gtf_f))
  } else {
    paste0("cat ", shQuote(gtf_f))
  }

  gtf_lines <- fread(
    cmd = paste0(decomp_cmd, " | grep -v '^#'"),
    sep = "\t",
    header = FALSE,
    select = c(3, 9),
    col.names = c("feature", "attributes")
  )
  gtf_lines <- gtf_lines[feature == "gene"]

  gtf_lines[, gene_id := str_match(attributes, 'gene_id "([^"]+)"')[, 2]]
  gtf_lines[, symbol := str_match(attributes, 'gene_name "([^"]+)"')[, 2]]
  gtf_lines[, gene_type := str_match(attributes, 'gene_type "([^"]+)"')[, 2]]
  na_type <- is.na(gtf_lines$gene_type)
  if (any(na_type)) {
    gtf_lines[na_type, gene_type := str_match(attributes, 'gene_biotype "([^"]+)"')[, 2]]
  }

  out_dt <- unique(gtf_lines[, .(gene_id, symbol, gene_type)], by = "gene_id")
  out_dt[is.na(symbol) | symbol == "", symbol := gene_id]
  out_dt[is.na(gene_type) | gene_type == "", gene_type := "unknown"]
  out_dt
}

read_sparse_h5_matrix <- function(h5_f) {
  h5 <- H5Fopen(h5_f, flags = "H5F_ACC_RDONLY")
  on.exit(H5Fclose(h5))

  mat <- sparseMatrix(
    i = as.vector(h5$matrix$indices + 1L),
    p = as.vector(h5$matrix$indptr),
    x = as.vector(h5$matrix$data),
    repr = "C",
    dims = h5$matrix$shape
  )
  colnames(mat) <- as.character(h5$matrix$barcodes)
  rownames(mat) <- as.character(h5$matrix$features$name)
  mat
}

sum_sua_counts <- function(stacked_mat) {
  suffixes <- c("S", "U", "A")
  mats <- lapply(suffixes, function(suffix) {
    idx <- grep(paste0("_", suffix, "$"), rownames(stacked_mat))
    stacked_mat[idx, , drop = FALSE]
  })
  names(mats) <- suffixes

  if (any(vapply(mats, nrow, integer(1)) == 0L)) {
    stop("Expected stacked S/U/A rows in filtered H5 matrix")
  }

  counts_mat <- Reduce(`+`, mats)
  rownames(counts_mat) <- sub("_[SUA]$", "", rownames(mats[[1]]))
  counts_mat
}

infer_sample_id_from_h5 <- function(h5_f) {
  sample_id <- basename(h5_f)
  sample_id <- sub("^filt_counts_", "", sample_id)
  sub("\\.h5$", "", sample_id)
}

order_cluster_labels <- function(cluster_ids) {
  cluster_ids <- unique(as.character(cluster_ids))
  if (length(cluster_ids) == 0) {
    return(character())
  }
  if (all(grepl("^cl[0-9]+$", cluster_ids))) {
    return(cluster_ids[order(as.integer(sub("^cl", "", cluster_ids)))])
  }
  sort(cluster_ids)
}

load_annotation_cells <- function(integration_f, sel_res, min_cl_size) {
  cl_var <- paste0("RNA_snn_res.", sel_res)
  int_dt <- fread(integration_f)
  assert_that(cl_var %in% names(int_dt), msg = sprintf("Integration output is missing '%s'", cl_var))

  keep_idx <- !is.na(int_dt[[cl_var]]) & !is.na(int_dt$sample_id) & int_dt$sample_id != ""
  if ("is_dbl" %in% names(int_dt)) {
    keep_idx <- keep_idx & !as.logical(int_dt$is_dbl)
  }
  if ("in_dbl_cl" %in% names(int_dt)) {
    keep_idx <- keep_idx & !as.logical(int_dt$in_dbl_cl)
  }

  ann_dt <- int_dt[keep_idx, .(sample_id, cell_id, cluster = get(cl_var))]
  cl_counts <- ann_dt[, .N, by = cluster][N >= min_cl_size]
  keep_clusters <- order_cluster_labels(cl_counts$cluster)
  ann_dt <- ann_dt[cluster %in% keep_clusters]
  ann_dt[, cluster := factor(as.character(cluster), levels = keep_clusters)]
  ann_dt
}

build_pseudobulk_from_h5s <- function(h5_files, ann_dt, biotypes_dt, n_cores = 1L) {
  sample_ids <- unique(ann_dt$sample_id)
  cluster_ids <- levels(ann_dt$cluster)

  h5_sample_ids <- vapply(h5_files, infer_sample_id_from_h5, character(1))
  h5_map <- setNames(h5_files, h5_sample_ids)
  assert_that(all(sample_ids %in% names(h5_map)), msg = "Missing filtered H5 files for one or more integrated samples")

  tasks <- lapply(sample_ids, function(sample_id_val) {
    list(
      sample_id = sample_id_val,
      h5_f = h5_map[[sample_id_val]],
      sample_cells = ann_dt[sample_id == sample_id_val, .(cell_id, cluster)]
    )
  })
  result_dir <- tempfile("annotation_pb_")
  dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(result_dir, recursive = TRUE, force = TRUE), add = TRUE)

  workers <- max(1L, min(as.integer(n_cores), length(sample_ids)))
  build_one_sample <- function(task) {
    counts_mat <- read_sparse_h5_matrix(task$h5_f) |> sum_sua_counts()

    sample_id_val <- task$sample_id
    sample_cells <- task$sample_cells
    global_ids <- paste(sample_id_val, colnames(counts_mat), sep = "_")
    matched <- match(global_ids, sample_cells$cell_id)
    keep_idx <- which(!is.na(matched))
    assert_that(length(keep_idx) > 0, msg = sprintf("No integrated cells matched filtered H5 for sample %s", sample_id_val))

    counts_sel <- counts_mat[, keep_idx, drop = FALSE]
    sample_clusters <- as.character(sample_cells$cluster[matched[keep_idx]])

    cluster_fac <- factor(sample_clusters, levels = cluster_ids)
    cluster_model <- sparseMatrix(
      i = seq_along(cluster_fac),
      j = as.integer(cluster_fac),
      x = 1,
      dims = c(length(cluster_fac), length(cluster_ids)),
      dimnames = list(NULL, cluster_ids)
    )

    cluster_sums <- counts_sel %*% cluster_model
    colnames(cluster_sums) <- cluster_ids
    n_cells <- Matrix::colSums(cluster_model)

    sample_result <- list(
      sample_id = sample_id_val,
      gene_ids = rownames(counts_sel),
      cluster_sums = cluster_sums,
      n_cells = as.numeric(n_cells)
    )

    out_f <- file.path(result_dir, sprintf("%s.rds", sample_id_val))
    saveRDS(sample_result, out_f, compress = FALSE)
    list(sample_id = sample_id_val, result_f = out_f)
  }

  if (workers == 1L) {
    sample_refs <- lapply(tasks, build_one_sample)
  } else {
    sample_refs <- parallel::mclapply(tasks, build_one_sample, mc.cores = workers)
  }

  sample_results <- lapply(sample_refs, function(ref) {
    assert_that(file.exists(ref$result_f), msg = sprintf("Missing pseudobulk temp result for sample %s", ref$sample_id))
    readRDS(ref$result_f)
  })

  gene_ids <- sample_results[[1]]$gene_ids
  for (res in sample_results) {
    assert_that(identical(gene_ids, res$gene_ids), msg = "Gene IDs differ across filtered H5 files")
  }
  sample_lu <- setNames(sample_results, vapply(sample_results, function(x) x$sample_id, character(1)))

  n_cells_dt <- rbindlist(lapply(sample_results, function(res) {
    data.table(
      sample_id = res$sample_id,
      cluster = cluster_ids,
      n_cells = as.integer(res$n_cells)
    )
  }))

  cluster_mats <- lapply(cluster_ids, function(cluster_id) {
    mats <- lapply(sample_ids, function(sample_id_val) {
      res <- sample_lu[[sample_id_val]]
      x <- res$cluster_sums[, cluster_id, drop = FALSE]
      colnames(x) <- sample_id_val
      x
    })
    do.call(cbind, mats)
  })
  names(cluster_mats) <- cluster_ids

  row_dt <- merge(
    data.table(gene_id = gene_ids),
    biotypes_dt,
    by = "gene_id",
    all.x = TRUE,
    sort = FALSE
  )
  row_dt[is.na(symbol) | symbol == "", symbol := gene_id]
  row_dt[is.na(gene_type) | gene_type == "", gene_type := "unknown"]

  list(cluster_mats = cluster_mats, n_cells_dt = n_cells_dt, row_dt = row_dt)
}

prepare_cluster_matrices <- function(pb_obj, min_cells = 10) {
  out <- list()

  for (cluster_id_val in names(pb_obj$cluster_mats)) {
    x <- pb_obj$cluster_mats[[cluster_id_val]]
    sample_ids <- colnames(x)
    cell_dt <- pb_obj$n_cells_dt[cluster == cluster_id_val]
    n_cells_vec <- cell_dt[match(sample_ids, sample_id), n_cells]

    use_idx <- n_cells_vec >= min_cells
    if (!any(use_idx)) {
      next
    }

    x <- x[, use_idx, drop = FALSE]
    sample_ids <- sample_ids[use_idx]
    n_cells_vec <- n_cells_vec[use_idx]

    if (ncol(x) > 1) {
      out_idx <- scater::isOutlier(Matrix::colSums(x), log = TRUE, type = "lower", nmads = 3)
      x <- x[, !out_idx, drop = FALSE]
      sample_ids <- sample_ids[!out_idx]
      n_cells_vec <- n_cells_vec[!out_idx]
    }

    if (ncol(x) == 0) {
      next
    }

    dge_obj <- edgeR::DGEList(x, remove.zeros = TRUE)
    dge_obj <- edgeR::calcNormFactors(dge_obj, method = "TMMwsp")
    lib_sizes <- dge_obj$samples$lib.size * dge_obj$samples$norm.factors

    out[[cluster_id_val]] <- list(
      counts = x,
      sample_ids = sample_ids,
      n_cells = n_cells_vec,
      lib_sizes = lib_sizes
    )
  }

  assert_that(length(out) > 0, msg = "No pseudobulk samples passed annotation minimum-cell filtering")
  out
}

make_logcpms_for_genes <- function(prep_obj, rows_dt, gene_ids, pseudo_count = 10) {
  gene_ids <- unique(gene_ids)
  if (length(gene_ids) == 0) {
    return(data.table())
  }

  out_ls <- lapply(names(prep_obj), function(cluster_id_val) {
    x <- prep_obj[[cluster_id_val]]$counts
    keep_idx <- which(rownames(x) %in% gene_ids)
    if (length(keep_idx) == 0) {
      return(NULL)
    }

    x <- x[keep_idx, , drop = FALSE]
    sm <- summary(x)
    dt <- data.table(
      gene_id = rownames(x)[sm$i],
      sample_id = colnames(x)[sm$j],
      count = as.numeric(sm$x)
    )

    # Keep explicit zeros for selected genes/samples so plotting behaves consistently.
    full_dt <- CJ(gene_id = rownames(x), sample_id = prep_obj[[cluster_id_val]]$sample_ids)
    dt <- merge(full_dt, dt, by = c("gene_id", "sample_id"), all.x = TRUE)
    dt[is.na(count), count := 0]
    dt <- merge(
      dt,
      data.table(
        sample_id = prep_obj[[cluster_id_val]]$sample_ids,
        lib_size = prep_obj[[cluster_id_val]]$lib_sizes,
        n_cells = prep_obj[[cluster_id_val]]$n_cells
      ),
      by = "sample_id",
      all.x = TRUE
    )
    dt[, `:=`(
      cluster = cluster_id_val,
      logcpm = log(count / lib_size * 1e6 + pseudo_count)
    )]
    dt
  })

  logcpms_dt <- rbindlist(Filter(Negate(is.null), out_ls), use.names = TRUE)
  merge(logcpms_dt, rows_dt, by = "gene_id", all.x = TRUE)
}

calc_find_markers_pseudobulk <- function(prep_obj, rows_dt) {
  cl_ls <- names(prep_obj)

  count_ls <- lapply(cl_ls, function(cluster_id) {
    x <- prep_obj[[cluster_id]]$counts
    colnames(x) <- sprintf("%s-%s", cluster_id, colnames(x))
    x
  })
  x_full <- do.call(cbind, count_ls)

  des_all <- data.table(
    col_lab = colnames(x_full),
    cluster = str_extract(colnames(x_full), "^[^-]+"),
    sample_id = str_extract(colnames(x_full), "(?<=-).+")
  )

  keep_gs <- edgeR::filterByExpr(edgeR::DGEList(x_full), group = des_all$cluster, min.count = 1)
  x_full <- x_full[keep_gs, , drop = FALSE]

  dge <- edgeR::DGEList(x_full, remove.zeros = TRUE)
  dge <- edgeR::calcNormFactors(dge, method = "TMMwsp")
  dge <- edgeR::estimateDisp(dge)
  logcpm_mat <- edgeR::cpm(dge, log = TRUE, prior.count = 10)

  mkrs_ls <- vector("list", length(cl_ls))
  zeros_ls <- vector("list", length(cl_ls))
  relcpms_ls <- vector("list", length(cl_ls))

  for (ii in seq_along(cl_ls)) {
    cluster_id <- cl_ls[[ii]]
    is_cluster <- des_all$cluster == cluster_id

    this_design <- model.matrix(~ is_cluster)
    fit <- edgeR::glmQLFit(dge, design = this_design)
    fit <- edgeR::glmTreat(fit, coef = "is_clusterTRUE")
    top_dt <- topTags(fit, n = Inf, sort.by = "none") |>
      as.data.frame() |>
      as.data.table(keep.rownames = "gene_id")
    top_dt[, cluster := cluster_id]
    mkrs_ls[[ii]] <- top_dt

    this_x <- x_full[, is_cluster, drop = FALSE]
    zeros_ls[[ii]] <- data.table(
      cluster = cluster_id,
      gene_id = rownames(this_x),
      n_zero = Matrix::rowSums(this_x == 0),
      n_cl = ncol(this_x)
    )

    relcpms_ls[[ii]] <- data.table(
      cluster = cluster_id,
      gene_id = rownames(logcpm_mat),
      logcpm.sel = rowMeans(logcpm_mat[, is_cluster, drop = FALSE]),
      logcpm.other = rowMeans(logcpm_mat[, !is_cluster, drop = FALSE])
    )
  }

  mkrs_dt <- rbindlist(mkrs_ls, use.names = TRUE)
  zeros_dt <- rbindlist(zeros_ls, use.names = TRUE)
  rel_cpms <- rbindlist(relcpms_ls, use.names = TRUE)

  mkrs_dt <- merge(mkrs_dt, zeros_dt, by = c("cluster", "gene_id"))
  mkrs_dt <- merge(mkrs_dt, rel_cpms, by = c("cluster", "gene_id"))
  mkrs_dt <- merge(mkrs_dt, rows_dt, by = "gene_id")
  mkrs_dt[order(cluster, PValue)]
}

normalize_marker_panel <- function(marker_f, biotypes_dt, marker_calcs_dt) {
  panel_dt <- fread(marker_f)

  if (!"label" %in% names(panel_dt) && all(c("cell_type", "marker_gene") %in% names(panel_dt))) {
    setnames(panel_dt, c("cell_type", "marker_gene"), c("label", "symbol"))
  }
  if ("ensembl_id" %in% names(panel_dt) && !"gene_id" %in% names(panel_dt)) {
    setnames(panel_dt, "ensembl_id", "gene_id")
  }

  assert_that("label" %in% names(panel_dt), msg = "Marker CSV must contain 'label' or 'cell_type'")
  assert_that(any(c("symbol", "gene_id") %in% names(panel_dt)), msg = "Marker CSV must contain 'symbol', 'gene_id', or 'marker_gene'")

  panel_dt[, input_order := seq_len(.N)]
  label_levels <- unique(panel_dt$label)

  if (!"symbol" %in% names(panel_dt)) {
    panel_dt <- merge(panel_dt, biotypes_dt[, .(gene_id, symbol)], by = "gene_id", all.x = TRUE)
  }
  if (!"gene_id" %in% names(panel_dt)) {
    gene_lu <- marker_calcs_dt[symbol %in% panel_dt$symbol,
      .SD[which.max(logcpm.sel)], by = symbol][, .(symbol, gene_id)]
    panel_dt <- merge(panel_dt, gene_lu, by = "symbol", all.x = TRUE)
  }

  dup_symbols <- panel_dt[, .N, by = symbol][N > 1 & !is.na(symbol)]
  assert_that(nrow(dup_symbols) == 0, msg = "Marker genes must be unique across labels within one CSV file")

  panel_dt <- panel_dt[!is.na(symbol) & !is.na(gene_id)]
  panel_dt <- panel_dt[gene_id %in% marker_calcs_dt$gene_id]
  assert_that(uniqueN(panel_dt$gene_id) >= 2, msg = "Fewer than two provided marker genes were found in the marker calculations")

  panel_dt[, label := factor(label, levels = label_levels)]
  setorder(panel_dt, input_order)
  panel_dt[, input_order := NULL]
  panel_dt[, .(label, symbol, gene_id)]
}

get_top_markers <- function(input_mkrs, fdr_cut = 0.01, n_top = 10, max_zero_p = 0.5) {
  top_tmp <- input_mkrs[order(cluster, -abs(logFC))]
  rbind(
    top_tmp[(FDR < fdr_cut) & (logFC > 0) & (n_zero / n_cl <= max_zero_p)],
    top_tmp[(FDR < fdr_cut) & (logFC < 0)],
    top_tmp[FDR >= fdr_cut & ((n_zero / n_cl <= max_zero_p) | (logFC < 0))]
  )[, .SD[seq_len(min(.N, n_top))], by = cluster][, prop_zero := n_zero / n_cl]
}

plot_top_marker_genes <- function(sel_cl, top_mkrs_dt, logcpms_all, cl_order = NULL, pseudo_count = 10) {
  sel_mkrs <- top_mkrs_dt[cluster == sel_cl]
  if (nrow(sel_mkrs) == 0) {
    return(ggplot() + theme_void() + labs(title = sprintf("No top markers for %s", sel_cl)))
  }

  fdr_show <- 0.01
  sel_mkrs[, symbol_lab := ifelse(
    FDR < fdr_show,
    sprintf("%s (FDR = %.1e)", symbol, FDR),
    sprintf("%s (FDR = %.3f)", symbol, FDR)
  )]

  plot_dt <- merge(logcpms_all, sel_mkrs[, .(gene_id, symbol_lab)], by = "gene_id")
  if (nrow(plot_dt) == 0) {
    return(ggplot() + theme_void() + labs(title = sprintf("No plottable top markers for %s", sel_cl)))
  }
  plot_dt[, is_sel := ifelse(cluster == sel_cl, "test cluster", "other")]
  if (is.null(cl_order)) {
    all_levels <- c(sel_cl, setdiff(order_cluster_labels(unique(plot_dt$cluster)), sel_cl))
    plot_dt[, cluster := factor(cluster, levels = all_levels)]
  } else {
    plot_dt[, cluster := factor(cluster, levels = cl_order)]
  }
  plot_dt[, symbol_lab := factor(symbol_lab, levels = unique(sel_mkrs$symbol_lab))]

  ggplot(plot_dt) +
    aes(x = cluster, y = logcpm, fill = is_sel, size = n_cells) +
    geom_point(position = position_jitter(width = 0.18, height = 0), shape = 21, colour = "black", alpha = 0.8) +
    scale_fill_manual(values = c(`test cluster` = "#1b1b1b", other = "#c7c7c7")) +
    scale_size(range = c(1.5, 6)) +
    facet_wrap(~ symbol_lab, scales = "free_y", ncol = 2) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    labs(
      x = "cluster",
      y = "log CPM",
      fill = "cluster status",
      size = "# cells",
      title = sprintf("Top %d marker genes for %s", nrow(sel_mkrs), sel_cl)
    )
}

plot_selected_marker_panel <- function(sel_panel_dt, logcpms_all, cl_order = NULL,
                                       highlight_clusters = character(), title = NULL) {
  sel_panel_dt <- unique(sel_panel_dt[, .(label, symbol, gene_id)])
  if (nrow(sel_panel_dt) == 0) {
    plot_title <- if (is.null(title)) "No marker genes available" else title
    return(ggplot() + theme_void() + labs(title = plot_title))
  }

  plot_dt <- merge(
    logcpms_all,
    sel_panel_dt[, .(gene_id, panel_symbol = symbol)],
    by = "gene_id"
  )
  if (nrow(plot_dt) == 0) {
    plot_title <- if (is.null(title)) "No marker genes matched the pseudobulk table" else title
    return(ggplot() + theme_void() + labs(title = plot_title))
  }
  plot_dt[, symbol := panel_symbol]
  plot_dt[, panel_symbol := NULL]

  plot_dt[, is_sel := if (length(highlight_clusters) > 0) {
    ifelse(cluster %in% highlight_clusters, "matched cluster", "other")
  } else {
    "all clusters"
  }]

  if (is.null(cl_order)) {
    plot_dt[, cluster := factor(cluster, levels = order_cluster_labels(unique(plot_dt$cluster)))]
  } else {
    plot_dt[, cluster := factor(cluster, levels = cl_order)]
  }
  plot_dt[, symbol := factor(symbol, levels = unique(sel_panel_dt$symbol))]

  fill_vals <- c(
    "matched cluster" = "#1b1b1b",
    other = "#c7c7c7",
    "all clusters" = "#1b1b1b"
  )

  ggplot(plot_dt) +
    aes(x = cluster, y = logcpm, fill = is_sel, size = n_cells) +
    geom_point(position = position_jitter(width = 0.18, height = 0), shape = 21, colour = "black", alpha = 0.8) +
    scale_fill_manual(values = fill_vals[unique(plot_dt$is_sel)], drop = FALSE) +
    scale_size(range = c(1.5, 6)) +
    facet_wrap(~ symbol, scales = "free_y", ncol = 2) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    labs(
      x = "cluster",
      y = "log CPM",
      fill = "cluster status",
      size = "# cells",
      title = if (is.null(title)) as.character(unique(sel_panel_dt$label)[1]) else title
    )
}

make_top_marker_heatmap_panel <- function(top_mkrs_dt, cl_order = NULL) {
  if (nrow(top_mkrs_dt) == 0) {
    return(data.table(label = factor(), symbol = character(), gene_id = character()))
  }

  panel_dt <- copy(top_mkrs_dt)
  if (!is.null(cl_order)) {
    panel_dt[, cluster := factor(as.character(cluster), levels = cl_order)]
  }
  panel_dt[, abs_logFC := abs(logFC)]
  setorder(panel_dt, cluster, -abs_logFC, FDR)
  panel_dt <- panel_dt[, .SD[1], by = gene_id]
  panel_dt[, abs_logFC := NULL]
  label_levels <- if (is.null(cl_order)) unique(as.character(panel_dt$cluster)) else cl_order
  panel_dt[, label := factor(as.character(cluster), levels = label_levels)]
  panel_dt[, .(label, symbol, gene_id)]
}

assign_cluster_labels_from_panel <- function(mkrs_dt, panel_dt, min_cpm = 10) {
  score_dt <- merge(
    mkrs_dt[gene_id %in% panel_dt$gene_id],
    unique(panel_dt[, .(gene_id, label)]),
    by = "gene_id"
  )

  if (nrow(score_dt) == 0) {
    return(data.table(cluster = unique(mkrs_dt$cluster), label = "Unassigned", label_score = NA_real_, n_markers = 0L))
  }

  score_dt <- score_dt[, .(
    label_score = mean(ifelse(logcpm.sel >= log(min_cpm + 1), pmax(logFC, 0), 0), na.rm = TRUE),
    n_markers = uniqueN(gene_id),
    n_positive = sum((logcpm.sel >= log(min_cpm + 1)) & (logFC > 0), na.rm = TRUE)
  ), by = .(cluster, label)]

  setorder(score_dt, cluster, -label_score, -n_positive, label)
  out_dt <- score_dt[, .SD[1], by = cluster]
  out_dt[is.na(label_score) | label_score <= 0, label := "Unassigned"]
  out_dt
}

plot_heatmap_of_selected_genes <- function(mkrs_dt, panel_dt, max_fc = 3, min_cpm = 10,
                                           pseudocount = 10) {
  mkrs_sel <- mkrs_dt[gene_id %in% panel_dt$gene_id]
  if (nrow(mkrs_sel) == 0) {
    stop("No selected marker genes were found in the marker statistics table")
  }

  max_cpms <- mkrs_sel[, .(max_logcpm = max(logcpm.sel, na.rm = TRUE)), by = symbol]
  keep_symbols <- max_cpms[max_logcpm >= log(min_cpm + pseudocount), symbol]
  if (length(keep_symbols) == 0) {
    keep_symbols <- unique(mkrs_sel$symbol)
  }
  mkrs_sel <- mkrs_sel[symbol %in% keep_symbols]

  symbol_order <- panel_dt[symbol %in% keep_symbols, unique(symbol)]
  log2fc_mat <- dcast(
    mkrs_sel[, .(symbol, cluster, log2fc = logFC)],
    symbol ~ cluster,
    value.var = "log2fc",
    fill = 0
  ) |> as.matrix(rownames = "symbol")
  log2fc_mat <- log2fc_mat[symbol_order, , drop = FALSE]

  signif_mat <- dcast(
    mkrs_sel[, .(
      symbol,
      cluster,
      signif = ifelse(
        logcpm.sel < log(min_cpm + 1),
        "",
        ifelse((FDR > 0.05) | (logFC < 0), "", ifelse(FDR > 0.01, "*", ifelse(FDR > 0.001, "**", "***")))
      )
    )],
    symbol ~ cluster,
    value.var = "signif",
    fill = ""
  ) |> as.matrix(rownames = "symbol")
  signif_mat <- signif_mat[symbol_order, , drop = FALSE]

  col_groups <- panel_dt[match(rownames(log2fc_mat), symbol), as.character(label)]
  col_groups <- factor(col_groups, levels = unique(as.character(panel_dt$label)))
  label_levels <- levels(col_groups)
  label_cols <- setNames(grDevices::hcl.colors(max(length(label_levels), 1), palette = "Set 2")[seq_along(label_levels)], label_levels)
  col_fun <- circlize::colorRamp2(c(-max_fc, 0, max_fc), c("#2166ac", "#f7f7f7", "#b2182b"))
  top_annot <- HeatmapAnnotation(
    label = col_groups,
    col = list(label = label_cols),
    show_annotation_name = FALSE,
    simple_anno_size = unit(3, "mm")
  )

  heatmap_mat <- t(log2fc_mat)
  signif_mat <- t(signif_mat)

  Heatmap(
    heatmap_mat,
    name = "log2FC",
    col = col_fun,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    column_split = col_groups,
    cluster_column_slices = FALSE,
    column_gap = unit(2, "mm"),
    top_annotation = top_annot,
    cell_fun = function(j, i, x, y, width, height, fill) {
      lab <- signif_mat[i, j]
      if (nzchar(lab)) {
        grid.text(lab, x, y, gp = gpar(fontsize = 8))
      }
    },
    column_names_gp = gpar(fontsize = 8),
    row_names_gp = gpar(fontsize = 8),
    column_names_side = "top",
    row_names_side = "left",
    row_title = NULL,
    column_title = NULL,
    heatmap_legend_param = list(at = c(-max_fc, 0, max_fc), title = "log2FC")
  )
}

plot_umap_cluster <- function(umap_dt, clust_dt, name) {
  stopifnot(
    all(c("UMAP1", "UMAP2") %in% names(umap_dt)),
    "cell_id" %in% names(umap_dt),
    "cell_id" %in% names(clust_dt),
    "cluster" %in% names(clust_dt)
  )

  plot_dt <- merge(umap_dt, clust_dt, by = "cell_id", all.x = TRUE)[, .(
    UMAP1 = rescale(UMAP1, to = c(0.05, 0.95)),
    UMAP2 = rescale(UMAP2, to = c(0.05, 0.95)),
    cluster
  )]
  if (is.factor(plot_dt$cluster)) {
    plot_dt[, cluster := factor(cluster, levels = levels(cluster))]
  } else {
    plot_dt[, cluster := factor(cluster, levels = unique(as.character(cluster)))]
  }
  plot_dt <- plot_dt[sample(.N, .N)]

  cl_labels <- plot_dt[, .(N = .N), by = cluster]
  cl_labels <- cl_labels[match(levels(plot_dt$cluster), cluster)]
  cl_labels[, cl_label := sprintf("%s (%d)", cluster, signif(N, 2))]
  label_lu <- setNames(cl_labels$cl_label, cl_labels$cluster)

  cl_cols <- rep(nice_cols, times = 10)[seq_along(label_lu)]
  names(cl_cols) <- label_lu
  n_col <- 4
  n_rows_lgd <- ceiling(length(label_lu) / n_col)

  ggplot(plot_dt) +
    aes(x = UMAP1, y = UMAP2, colour = label_lu[as.character(cluster)]) +
    geom_point_rast(size = 0.1, raster.dpi = 150) +
    scale_colour_manual(values = cl_cols,
      guide = guide_legend(override.aes = list(size = 3), nrow = n_rows_lgd)) +
    scale_x_continuous(breaks = pretty_breaks(), limits = c(0, 1)) +
    scale_y_continuous(breaks = pretty_breaks(), limits = c(0, 1)) +
    theme_bw() +
    theme(panel.grid = element_blank(), aspect.ratio = 1,
      legend.title.position = "left", legend.position = "bottom",
      axis.ticks = element_blank(), axis.text = element_blank()) +
    labs(colour = name)
}

load_h5_marker_expression_multi <- function(h5_files, sel_dt_list, int_dt, pseudo_count = 10) {
  if (length(sel_dt_list) == 0) {
    return(list())
  }

  clean_sel <- lapply(sel_dt_list, function(sel_dt) {
    if (is.null(sel_dt) || nrow(sel_dt) == 0) {
      return(data.table(label = character(), symbol = character(), gene_id = character()))
    }
    unique(sel_dt[, .(label, symbol, gene_id)])
  })

  marker_ids <- unique(unlist(lapply(clean_sel, function(x) x$gene_id), use.names = FALSE))
  out_names <- names(clean_sel)
  if (is.null(out_names) || any(out_names == "")) {
    out_names <- paste0("set", seq_along(clean_sel))
    names(clean_sel) <- out_names
  }

  if (length(marker_ids) == 0 || length(h5_files) == 0) {
    empty <- data.table()
    return(setNames(replicate(length(clean_sel), empty, simplify = FALSE), names(clean_sel)))
  }

  cell_meta <- unique(int_dt[, .(sample_id, cell_id, UMAP1, UMAP2)])
  out_ls <- setNames(lapply(seq_along(clean_sel), function(...) vector("list", length(h5_files))), names(clean_sel))

  for (ii in seq_along(h5_files)) {
    h5_f <- h5_files[[ii]]
    sample_id_val <- infer_sample_id_from_h5(h5_f)
    sample_meta <- cell_meta[sample_id == sample_id_val]
    if (nrow(sample_meta) == 0) {
      next
    }

    counts_mat <- sum_sua_counts(read_sparse_h5_matrix(h5_f))
    row_idx <- which(rownames(counts_mat) %in% marker_ids)
    if (length(row_idx) == 0) {
      next
    }

    global_ids <- paste(sample_id_val, colnames(counts_mat), sep = "_")
    keep_idx <- which(global_ids %in% sample_meta$cell_id)
    if (length(keep_idx) == 0) {
      next
    }

    lib_sizes <- Matrix::colSums(counts_mat[, keep_idx, drop = FALSE])
    expr_mat <- as.matrix(counts_mat[row_idx, keep_idx, drop = FALSE])
    expr_mat <- sweep(expr_mat, 2, pmax(lib_sizes, 1), "/") * 1e6
    expr_mat <- log(expr_mat + pseudo_count)
    colnames(expr_mat) <- global_ids[keep_idx]

    expr_base_dt <- as.data.table(expr_mat, keep.rownames = "gene_id")
    expr_base_dt <- melt(expr_base_dt, id.vars = "gene_id", variable.name = "cell_id", value.name = "logcount")
    expr_base_dt <- merge(expr_base_dt, sample_meta[, .(cell_id, UMAP1, UMAP2)], by = "cell_id", all.x = TRUE)

    for (nm in names(clean_sel)) {
      sel_dt <- clean_sel[[nm]]
      if (nrow(sel_dt) == 0) {
        next
      }
      tmp_dt <- merge(expr_base_dt, sel_dt, by = "gene_id", all = FALSE, allow.cartesian = TRUE)
      out_ls[[nm]][[ii]] <- tmp_dt
    }
  }

  out_dt <- lapply(names(clean_sel), function(nm) {
    dt <- rbindlist(Filter(Negate(is.null), out_ls[[nm]]), use.names = TRUE)
    if (nrow(dt) == 0) {
      return(dt)
    }
    dt[, expr := logcount]
    dt[, symbol := factor(symbol, levels = unique(clean_sel[[nm]]$symbol))]
    dt
  })
  names(out_dt) <- names(clean_sel)
  out_dt
}

load_h5_marker_expression <- function(h5_files, sel_dt, int_dt, pseudo_count = 10) {
  load_h5_marker_expression_multi(
    h5_files,
    sel_dt_list = list(default = sel_dt),
    int_dt = int_dt,
    pseudo_count = pseudo_count
  )$default
}

plot_expression_umap_pair <- function(expr_dt, meta_dt,
                                      cluster_col = "cluster",
                                      gene_name = NULL,
                                      cluster_name = NULL) {
  stopifnot(
    all(c("cell_id", "expr") %in% names(expr_dt)),
    all(c("cell_id", "UMAP1", "UMAP2", cluster_col) %in% names(meta_dt))
  )

  expr_dt <- copy(expr_dt)[, .(expr = max(expr, na.rm = TRUE)), by = cell_id]
  plot_dt <- merge(meta_dt, expr_dt, by = "cell_id", all.x = TRUE)
  plot_dt[is.na(expr), expr := 0]

  cluster_vals <- meta_dt[[cluster_col]]
  cluster_levels <- if (is.factor(cluster_vals)) {
    levels(cluster_vals)
  } else {
    order_cluster_labels(unique(cluster_vals))
  }
  plot_dt <- plot_dt[!is.na(get(cluster_col))]
  plot_dt[, cluster_plot := factor(as.character(get(cluster_col)), levels = cluster_levels)]

  label_dt <- plot_dt[, .(
    UMAP1 = median(UMAP1),
    UMAP2 = median(UMAP2)
  ), by = cluster_plot][!is.na(cluster_plot)]

  cl_cols <- rep(nice_cols, times = ceiling(max(1, length(cluster_levels)) / length(nice_cols)) + 1L)[seq_along(cluster_levels)]
  names(cl_cols) <- cluster_levels
  cluster_title <- if (is.null(cluster_name)) cluster_col else cluster_name
  gene_title <- if (is.null(gene_name) || gene_name == "") "marker" else gene_name

  base_theme <- theme_bw() +
    theme(
      panel.grid = element_blank(),
      aspect.ratio = 1,
      legend.title.position = "left",
      legend.position = "bottom",
      axis.ticks = element_blank(),
      axis.text = element_blank()
    )

  g_cluster <- ggplot(plot_dt[sample(.N, .N)]) +
    aes(x = UMAP1, y = UMAP2, colour = cluster_plot) +
    geom_point_rast(size = 0.1, raster.dpi = 150) +
    geom_text(
      data = label_dt,
      aes(x = UMAP1, y = UMAP2, label = cluster_plot),
      inherit.aes = FALSE,
      colour = "black",
      size = 3
    ) +
    scale_colour_manual(values = cl_cols, guide = "none") +
    base_theme +
    labs(title = cluster_title)

  g_expr <- ggplot(plot_dt[order(expr)]) +
    aes(x = UMAP1, y = UMAP2, colour = expr) +
    geom_point_rast(size = 0.1, raster.dpi = 150) +
    scale_colour_viridis_c(
      option = "magma",
      breaks = pretty_breaks(),
      guide = guide_colourbar(title.position = "top"),
      na.value = "grey90"
    ) +
    base_theme +
    labs(title = sprintf("%s expression", gene_title), colour = "log CPM")

  g_cluster + g_expr + plot_layout(widths = c(1, 1))
}

downsample_feature_umap_cells <- function(plot_dt, expr_cols, max_cells = 30000L, keep_top_frac = 0.2) {
  stopifnot(all(expr_cols %in% names(plot_dt)))

  if (nrow(plot_dt) <= max_cells) {
    return(plot_dt)
  }

  keep_top_n <- min(as.integer(ceiling(max_cells * keep_top_frac)), nrow(plot_dt), max_cells)
  expr_mat <- as.matrix(plot_dt[, ..expr_cols])
  max_expr <- do.call(pmax, c(split(expr_mat, col(expr_mat)), list(na.rm = TRUE)))
  max_expr[!is.finite(max_expr)] <- 0

  top_idx <- order(max_expr, decreasing = TRUE)[seq_len(keep_top_n)]
  if (keep_top_n >= max_cells) {
    keep_idx <- top_idx
  } else {
    rest_idx <- setdiff(seq_len(nrow(plot_dt)), top_idx)
    rest_n <- min(length(rest_idx), max_cells - keep_top_n)
    keep_idx <- c(top_idx, if (rest_n > 0) sample(rest_idx, rest_n) else integer())
  }

  plot_dt[sort(keep_idx)]
}

plot_top_marker_umap_facet <- function(sel_cl, top_mkrs_dt, top_expr_dt, int_umap_dt, ncol = 4,
                                       max_cells = 30000L, raster_dpi = 150) {
  sel_mkrs <- top_mkrs_dt[cluster == sel_cl, .(gene_id, symbol)]
  if (nrow(sel_mkrs) == 0 || nrow(top_expr_dt) == 0) return(NULL)

  gene_levels <- unique(sel_mkrs$symbol)
  base_dt <- int_umap_dt[, .(cell_id, UMAP1, UMAP2)]

  sel_expr <- top_expr_dt[
    as.character(label) == as.character(sel_cl) & symbol %in% gene_levels,
    .(cell_id, symbol, expr)
  ]
  if (nrow(sel_expr) == 0) return(NULL)

  expr_wide <- dcast(sel_expr, cell_id ~ symbol, value.var = "expr", fill = 0)
  plot_dt <- merge(base_dt, expr_wide, by = "cell_id", all.x = TRUE)
  present_genes <- gene_levels[gene_levels %in% names(plot_dt)]
  if (length(present_genes) == 0) return(NULL)
  for (g in present_genes) set(plot_dt, which(is.na(plot_dt[[g]])), g, 0)
  plot_dt <- downsample_feature_umap_cells(plot_dt, present_genes, max_cells = max_cells)

  plot_dt_long <- melt(
    plot_dt,
    id.vars    = c("cell_id", "UMAP1", "UMAP2"),
    measure.vars = present_genes,
    variable.name = "symbol",
    value.name    = "expr"
  )
  plot_dt_long[, symbol := factor(as.character(symbol), levels = present_genes)]

  ggplot(plot_dt_long[order(expr)], aes(UMAP1, UMAP2, colour = expr)) +
    geom_point_rast(size = 0.15, stroke = 0, raster.dpi = raster_dpi) +
    scale_colour_viridis_c(option = "magma", guide = "none") +
    facet_wrap(~ symbol, ncol = ncol) +
    theme_void() +
    theme(
      strip.text   = element_text(size = 7, margin = margin(b = 2)),
      panel.spacing = unit(2, "pt")
    )
}