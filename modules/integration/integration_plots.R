# integration_plots.R — plotting helper functions for integration_report.qmd
# Sourced into the Quarto report; not run directly.
# Replicates scprocess/code/integration.R plotting functions exactly.

suppressPackageStartupMessages({
  library(data.table)
  library(forcats)
  library(ggh4x)
  library(ggplot2)
  library(ggplot.multistats)
  library(patchwork)
  library(scales)
  library(ggrepel)
})

# ---------------------------------------------------------------------------
# Colour palette (matches scprocess nice_cols)
# ---------------------------------------------------------------------------
nice_cols <- c(
  '#DC050C', '#FB8072', '#1965B0', '#7BAFDE', '#882E72', '#B17BA6',
  '#FF7F00', '#FDB462', '#E7298A', '#E78AC3', '#33A02C', '#B2DF8A'
)

order_cluster_labels <- function(cluster_ids) {
  cluster_ids <- unique(as.character(cluster_ids))
  if (length(cluster_ids) == 0) {
    return(character())
  }
  if (all(grepl('^cl[0-9]+$', cluster_ids))) {
    return(cluster_ids[order(as.integer(sub('^cl', '', cluster_ids)))])
  }
  sort(cluster_ids)
}

make_plot_palette <- function(levels_vec) {
  levels_vec <- unique(as.character(levels_vec))
  if (length(levels_vec) == 0) {
    return(setNames(character(), character()))
  }
  palette_vals <- rep(nice_cols, length.out = length(levels_vec))
  names(palette_vals) <- levels_vec
  palette_vals
}

cluster_resolution_label <- function(cl_col) {
  res_val <- sub('^RNA_snn_res\\.', '', cl_col)
  sub('^leiden_', '', res_val)
}

# ---------------------------------------------------------------------------
# UMAP density (binned 2D histogram, log10 colour scale)
# Matches scprocess plot_umap_density exactly.
# ---------------------------------------------------------------------------
plot_umap_density <- function(input_dt) {
  umap_dt <- copy(input_dt)[, .(
    UMAP1 = rescale(UMAP1, to = c(0.05, 0.95)),
    UMAP2 = rescale(UMAP2, to = c(0.05, 0.95))
  )]

  ggplot(umap_dt) + aes(x = UMAP1, y = UMAP2) +
    geom_bin2d(bins = 50) +
    scale_fill_distiller(palette = "RdBu", trans = "log10") +
    scale_x_continuous(breaks = pretty_breaks(), limits = c(0, 1)) +
    scale_y_continuous(breaks = pretty_breaks(), limits = c(0, 1)) +
    theme_bw() +
    theme(panel.grid = element_blank(),
      legend.title.position = "left", legend.position = "bottom",
      axis.ticks = element_blank(), axis.text = element_blank(),
      aspect.ratio = 1)
}

# ---------------------------------------------------------------------------
# Doublet UMAP: binned hex showing mean doublet proportion per bin
# Matches scprocess plot_umap_doublets exactly.
# ---------------------------------------------------------------------------
plot_umap_doublets <- function(int_dt) {

  # Two-pass output: use dbl_UMAP1/2 from pass 1 and is_dbl column
  if (all(c("dbl_UMAP1", "dbl_UMAP2", "is_dbl") %in% names(int_dt))) {
    dbl_dt <- copy(int_dt)[!is.na(dbl_UMAP1), .(
      UMAP1  = rescale(dbl_UMAP1, to = c(0.05, 0.95)),
      UMAP2  = rescale(dbl_UMAP2, to = c(0.05, 0.95)),
      is_dbl = as.numeric(is_dbl)
    )]
  } else {
    # Fallback for legacy format
    stopifnot(all(c("UMAP1", "UMAP2", "is_doublet") %in% names(int_dt)))
    dbl_dt <- copy(int_dt)[, .(
      UMAP1  = rescale(UMAP1, to = c(0.05, 0.95)),
      UMAP2  = rescale(UMAP2, to = c(0.05, 0.95)),
      is_dbl = as.numeric(is_doublet)
    )]
  }

  ggplot(dbl_dt) +
    aes(x = UMAP1, y = UMAP2, z = is_dbl) +
    stat_summary_hex(fun = 'mean', bins = 30) +
    scale_fill_distiller(palette = 'PiYG', limits = c(0, 1),
      breaks = pretty_breaks()) +
    labs(fill = 'mean doublets\nby bin') +
    scale_x_continuous(breaks = pretty_breaks()) +
    scale_y_continuous(breaks = pretty_breaks()) +
    theme_bw() +
    theme(panel.grid = element_blank(), axis.ticks = element_blank(),
      axis.text = element_blank(), aspect.ratio = 1)
}

# ---------------------------------------------------------------------------
# Doublet cluster scatter: log10(cluster size) vs doublet pct
# Matches scprocess plot_doublet_clusters exactly.
# ---------------------------------------------------------------------------
plot_doublet_clusters <- function(int_dt, dbl_cl_prop, cl_col = NULL) {
  # Two-pass format: dbl_cluster and is_dbl are precomputed per cell in int_dt
  if (all(c("dbl_cluster", "is_dbl") %in% names(int_dt))) {
    dbl_clusts <- int_dt[!is.na(dbl_cluster), .(
      dbl_cl_size = .N,
      dbl_cl_prop = sum(as.logical(is_dbl), na.rm = TRUE) / .N * 100
    ), by = .(dbl_cluster)]
  } else {
    # Legacy fallback: cluster column from pass-2 + is_doublet flag
    stopifnot(!is.null(cl_col))
    cl_col <- resolve_cluster_column(int_dt, cl_col)
    stopifnot(all(c(cl_col, "is_doublet") %in% names(int_dt)))
    dbl_clusts <- int_dt[, .(
      dbl_cl_size = .N,
      dbl_cl_prop = sum(is_doublet) / .N * 100
    ), by = cl_col]
    setnames(dbl_clusts, cl_col, "dbl_cluster")
  }

  log_brks <- log10(c(1e1, 2e1, 5e1, 1e2, 2e2, 5e2, 1e3, 2e3, 5e3, 1e4, 2e4, 5e4))
  log_labs <- c("10", "20", "50", "100", "200", "500",
    "1k", "2k", "5k", "10k", "20k", "50k")

  ggplot(dbl_clusts) +
    aes(x = log10(dbl_cl_size), y = dbl_cl_prop) +
    geom_hline(yintercept = dbl_cl_prop * 100, linetype = "dashed", colour = 'black') +
    geom_point() +
    scale_x_continuous(breaks = log_brks, labels = log_labs) +
    scale_y_continuous(breaks = pretty_breaks()) +
    expand_limits(y = c(0, 100)) +
    theme_classic() +
    labs(x = "# cells in cluster", y = "doublet pct. of cluster")
}

# ---------------------------------------------------------------------------
# UMAP coloured by Leiden cluster assignment
# Matches scprocess plot_umap_cluster exactly.
# ---------------------------------------------------------------------------
plot_umap_cluster <- function(umap_dt, clust_dt, name) {
  stopifnot(
    all(c('UMAP1', 'UMAP2') %in% names(umap_dt)),
    'cell_id' %in% names(umap_dt),
    'cell_id' %in% names(clust_dt),
    'cluster' %in% names(clust_dt)
  )

  plot_dt <- merge(umap_dt, clust_dt, by = 'cell_id', all.x = TRUE)[, .(
    UMAP1   = rescale(UMAP1, to = c(0.05, 0.95)),
    UMAP2   = rescale(UMAP2, to = c(0.05, 0.95)),
    cluster,
    centroid_group = if ('centroid_group' %in% names(clust_dt)) centroid_group else cluster,
    centroid_label = if ('centroid_label' %in% names(clust_dt)) centroid_label else cluster
  )]

  if (is.factor(clust_dt$cluster)) {
    cluster_levels <- levels(clust_dt$cluster)
  } else {
    cluster_levels <- order_cluster_labels(plot_dt$cluster)
  }
  plot_dt[, cluster := factor(as.character(cluster), levels = cluster_levels)]
  plot_dt[, centroid_group := as.character(centroid_group)]
  plot_dt[, centroid_label := as.character(centroid_label)]
  plot_dt <- plot_dt[sample(.N, .N)]

  label_dt <- plot_dt[
    !is.na(centroid_group) & nzchar(centroid_group),
    .(
      UMAP1 = median(UMAP1),
      UMAP2 = median(UMAP2),
      centroid_label = centroid_label[1]
    ),
    by = centroid_group
  ]

  # define colours
  cl_cols <- make_plot_palette(levels(plot_dt$cluster))
  n_col      <- 4
  n_rows_lgd <- ceiling(length(cl_cols) / n_col)

  ggplot(plot_dt) +
    aes(x = UMAP1, y = UMAP2, colour = cluster) +
    geom_point(size = 0.1) +
    geom_text_repel(
      data = label_dt,
      aes(x = UMAP1, y = UMAP2, label = centroid_label),
      inherit.aes = FALSE,
      size = 3,
      seed = 1,
      min.segment.length = 0,
      max.overlaps = Inf,
      box.padding = 0.3,
      point.padding = 0.15,
      segment.color = 'grey55'
    ) +
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

plot_cluster_composition_bars <- function(int_dt, cl_cols, split_var) {
  stopifnot(length(cl_cols) > 0, split_var %in% names(int_dt))

  plot_ls <- lapply(cl_cols, function(cl_col) {
    stopifnot(cl_col %in% names(int_dt))
    this_dt <- copy(int_dt)[
      !is.na(get(cl_col)) & !is.na(get(split_var)) & as.character(get(split_var)) != '',
      .(cluster = as.character(get(cl_col)), split = as.character(get(split_var)))
    ]
    if (nrow(this_dt) == 0) {
      return(NULL)
    }
    this_dt[, cluster := factor(cluster, levels = order_cluster_labels(cluster))]
    this_dt[, resolution := factor(cluster_resolution_label(cl_col), levels = vapply(cl_cols, cluster_resolution_label, character(1)))]
    this_dt[, .N, by = .(resolution, cluster, split)][, prop := N / sum(N), by = .(resolution, cluster)]
  })

  plot_dt <- rbindlist(plot_ls, use.names = TRUE, fill = TRUE)
  if (nrow(plot_dt) == 0) {
    return(
      ggplot() +
        annotate('text', x = 0, y = 0, label = sprintf('No cells had usable values for %s', split_var)) +
        theme_void()
    )
  }

  split_levels <- plot_dt[, .(N = sum(N)), by = split][order(-N), split]
  plot_dt[, split := factor(split, levels = split_levels)]
  split_cols <- setNames(
    rep(grDevices::hcl.colors(max(length(split_levels), 3), palette = 'Set 3'), length.out = length(split_levels)),
    split_levels
  )

  ggplot(plot_dt, aes(x = cluster, y = prop, fill = split)) +
    geom_col(width = 0.82) +
    facet_wrap(~ resolution, ncol = 1, scales = 'free_x') +
    scale_fill_manual(values = split_cols, drop = FALSE) +
    scale_y_continuous(labels = label_percent(accuracy = 1), expand = expansion(mult = c(0, 0.02))) +
    theme_bw() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      strip.background = element_rect(fill = 'white'),
      legend.position = 'bottom',
      legend.title.position = 'left'
    ) +
    guides(fill = guide_legend(ncol = max(1, ceiling(length(split_levels) / 12)))) +
    labs(x = 'cluster', y = 'proportion', fill = split_var)
}

# ---------------------------------------------------------------------------
# UMAP coloured by sample
# ---------------------------------------------------------------------------
plot_umap_sample <- function(int_dt) {
  stopifnot(all(c("UMAP1", "UMAP2", "sample_id") %in% names(int_dt)))

  plot_dt <- copy(int_dt)[, .(
    UMAP1     = rescale(UMAP1, to = c(0.05, 0.95)),
    UMAP2     = rescale(UMAP2, to = c(0.05, 0.95)),
    sample_id
  )][sample(.N, .N)]

  n_samples <- uniqueN(plot_dt$sample_id)

  ggplot(plot_dt, aes(x = UMAP1, y = UMAP2, colour = sample_id)) +
    geom_point(size = 0.1) +
    scale_x_continuous(breaks = pretty_breaks(), limits = c(0, 1)) +
    scale_y_continuous(breaks = pretty_breaks(), limits = c(0, 1)) +
    guides(colour = guide_legend(
      override.aes = list(size = 3),
      ncol = max(1, ceiling(n_samples / 20))
    )) +
    theme_bw() +
    theme(panel.grid = element_blank(), aspect.ratio = 1,
      legend.title.position = "left", legend.position = "bottom",
      axis.ticks = element_blank(), axis.text = element_blank()) +
    labs(colour = "sample")
}

# ---------------------------------------------------------------------------
# UMAP coloured by a single metadata / categorical variable
# ---------------------------------------------------------------------------
plot_umap_by_var <- function(int_dt, var_name) {
  stopifnot(all(c("UMAP1", "UMAP2", var_name) %in% names(int_dt)))

  plot_dt <- copy(int_dt)[!is.na(get(var_name)), c("UMAP1", "UMAP2", var_name), with = FALSE]
  plot_dt[, UMAP1 := rescale(UMAP1, to = c(0.05, 0.95))]
  plot_dt[, UMAP2 := rescale(UMAP2, to = c(0.05, 0.95))]
  plot_dt <- plot_dt[sample(.N, .N)]

  if (is.numeric(plot_dt[[var_name]])) {
    ggplot(plot_dt, aes(x = UMAP1, y = UMAP2, colour = get(var_name))) +
      geom_point(size = 0.1) +
      scale_x_continuous(breaks = pretty_breaks(), limits = c(0, 1)) +
      scale_y_continuous(breaks = pretty_breaks(), limits = c(0, 1)) +
      scale_colour_distiller(palette = "YlGnBu", direction = 1) +
      theme_bw() +
      theme(panel.grid = element_blank(), aspect.ratio = 1,
        legend.title.position = "left", legend.position = "bottom",
        axis.ticks = element_blank(), axis.text = element_blank()) +
      labs(colour = var_name)
  } else {
    ggplot(plot_dt, aes(x = UMAP1, y = UMAP2, colour = as.factor(get(var_name)))) +
      geom_point(size = 0.1) +
      scale_x_continuous(breaks = pretty_breaks(), limits = c(0, 1)) +
      scale_y_continuous(breaks = pretty_breaks(), limits = c(0, 1)) +
      guides(colour = guide_legend(override.aes = list(size = 3))) +
      theme_bw() +
      theme(panel.grid = element_blank(), aspect.ratio = 1,
        legend.title.position = "left", legend.position = "bottom",
        axis.ticks = element_blank(), axis.text = element_blank()) +
      labs(colour = var_name)
  }
}

# ---------------------------------------------------------------------------
# Resolve cluster column name (leiden_X vs RNA_snn_res.X)
# ---------------------------------------------------------------------------
resolve_cluster_column <- function(int_dt, cl_col) {
  if (cl_col %in% names(int_dt)) return(cl_col)
  if (grepl('^leiden_', cl_col)) {
    alt_col <- sub('^leiden_', 'RNA_snn_res.', cl_col)
    if (alt_col %in% names(int_dt)) return(alt_col)
  }
  if (grepl('^RNA_snn_res\\.', cl_col)) {
    alt_col <- sub('^RNA_snn_res\\.', 'leiden_', cl_col)
    if (alt_col %in% names(int_dt)) return(alt_col)
  }
  stop(sprintf('Cluster column not found: %s', cl_col))
}

# ---------------------------------------------------------------------------
# Cluster mixing diagnostics across samples
# Matches scprocess plot_cluster_entropies exactly.
# ---------------------------------------------------------------------------
plot_cluster_entropies <- function(input_dt, batch_var = 'sample_id', what = c('norm', 'raw')) {
  what <- match.arg(what)
  stopifnot(all(c('cluster', batch_var) %in% names(input_dt)))

  ns_dt <- copy(input_dt)[, .N, by = c(batch_var, 'cluster')][
    , p_sample := N / sum(N), by = batch_var
  ][
    , p_cluster := N / sum(N), by = cluster
  ][
    , p_cl_norm := p_sample / sum(p_sample), by = cluster
  ]

  entropy_dt <- ns_dt[, .(
    h_cl_raw     = -sum(p_cluster * log2(p_cluster), na.rm = TRUE),
    h_cl_norm    = -sum(p_cl_norm * log2(p_cl_norm), na.rm = TRUE),
    max_pct_raw  = 100 * max(p_cluster),
    max_pct_norm = 100 * max(p_cl_norm),
    N            = sum(N)
  ), by = cluster]
  labels_dt <- entropy_dt[order(cluster)]

  cl_ls   <- sort(unique(entropy_dt$cluster))
  cl_cols <- nice_cols[seq_along(cl_ls)]
  names(cl_cols) <- cl_ls

  ggplot(entropy_dt) +
    aes_string(x = paste0('h_cl_', what), y = paste0('max_pct_', what)) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "grey") +
    geom_text_repel(data = labels_dt, aes(label = cluster), size = 3,
      min.segment.length = 0, max.overlaps = Inf, box.padding = 0.5) +
    geom_point(shape = 21, aes(size = sqrt(N), fill = cluster)) +
    scale_x_continuous(breaks = pretty_breaks(n = 3)) +
    scale_y_continuous(breaks = pretty_breaks(n = 3)) +
    scale_fill_manual(values = cl_cols, guide = "none") +
    expand_limits(x = 0, y = 0) +
    scale_size(
      range  = c(1, 8),
      breaks = sqrt(c(2e2, 5e2, 1e3, 2e3, 5e3, 1e4, 2e4, 5e4)),
      labels = c('200', '500', '1k', '2k', '5k', '10k', '20k', '50k')
    ) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    labs(
      x    = 'entropy',
      y    = sprintf('max. pct. of one %s', batch_var),
      size = 'total # cells'
    )
}

# ---------------------------------------------------------------------------
# QC metric distributions per cluster
# Matches scprocess plot_cluster_qc_distns exactly.
# ---------------------------------------------------------------------------
plot_cluster_qc_distns <- function(qc_melt, clust_dt, name, min_cl_size = 1e2) {
  stopifnot(all(c('cell_id', 'cluster') %in% names(clust_dt)))

  if (is.numeric(name)) {
    # exclude tiny clusters
    cl_n_dt  <- clust_dt[, .N, by = cluster]
    cls_tiny <- cl_n_dt$N < min_cl_size
    if (any(cls_tiny)) {
      clust_dt <- clust_dt[cluster %in% cl_n_dt[!cls_tiny]$cluster]
    }
    clust_dt[, cluster := fct_infreq(as.factor(cluster))]
    name <- sprintf('res = %s', name)
  } else {
    clust_dt[, cluster := fct_infreq(as.factor(cluster))]
  }

  if (nrow(clust_dt) == 0) {
    return(
      ggplot() +
        annotate("text", x = 0, y = 0, label = sprintf("%s: no clusters passed QC plot filters", name)) +
        theme_void()
    )
  }

  plot_dt <- merge(qc_melt, clust_dt, by = 'cell_id')

  if (nrow(plot_dt) == 0 || !any(!is.na(plot_dt$qc_full))) {
    return(
      ggplot() +
        annotate("text", x = 0, y = 0, label = sprintf("%s: no QC rows matched integrated cells", name)) +
        theme_void()
    )
  }

  cl_cols <- rep(nice_cols, times = 10)[seq_along(levels(plot_dt$cluster))]
  names(cl_cols) <- levels(plot_dt$cluster)

  log_brks   <- log10(c(1e1, 2e1, 5e1, 1e2, 2e2, 5e2, 1e3, 2e3, 5e3, 1e4, 2e4, 5e4, 1e5, 2e5, 5e5))
  log_labs   <- c("10", "20", "50", "100", "200", "500",
    "1k", "2k", "5k", "10k", "20k", "50k", "100k", "200k", "500k")
  logit_brks <- qlogis(c(1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 0.10, 0.30, 0.50, 0.70, 0.90, 0.97, 0.99))
  logit_labs <- c("0.01%", "0.03%", "0.1%", "0.3%", "1%", "3%", "10%", "30%",
    "50%", "70%", "90%", "97%", "99%")

  ggplot(plot_dt) + aes(x = cluster, y = qc_val, fill = cluster) +
    geom_violin(kernel = 'rectangular', adjust = 0.5, scale = 'width', width = 0.8, colour = NA) +
    scale_fill_manual(values = cl_cols, guide = "none") +
    facet_grid(qc_full ~ ., scales = 'free_y') +
    facetted_pos_scales(
      y = list(
        qc_full == "no. of UMIs"  ~ scale_y_continuous(breaks = log_brks, labels = log_labs),
        qc_full == "no. of genes" ~ scale_y_continuous(breaks = log_brks, labels = log_labs),
        qc_full == "mito pct."    ~ scale_y_continuous(breaks = logit_brks, labels = logit_labs),
        qc_full == "spliced pct." ~ scale_y_continuous(breaks = logit_brks, labels = logit_labs)
      )
    ) +
    theme_bw() +
    theme(
      axis.text.x      = element_text(angle = 90, hjust = 1, vjust = 0.5),
      panel.grid       = element_blank(),
      strip.background = element_rect(fill = 'white')
    ) +
    labs(x = name, y = "QC metric")
}

# ---------------------------------------------------------------------------
# Cell counts + percentage per cluster (returns a data.table for kable)
# ---------------------------------------------------------------------------
cluster_counts_tbl <- function(int_dt, cl_col) {
  cl_col <- resolve_cluster_column(int_dt, cl_col)
  stopifnot(cl_col %in% names(int_dt))

  tbl <- int_dt[, .N, by = cl_col][order(-N)]
  setnames(tbl, c("cluster", "n_cells"))
  tbl[, pct := round(n_cells / sum(n_cells) * 100, 1)]
  tbl
}
