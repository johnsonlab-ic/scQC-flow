# integration_plots.R — plotting helper functions for integration_report.qmd
# Sourced into the Quarto report; not run directly.
# Adapted from scprocess/scripts/integration.R.

suppressPackageStartupMessages({
  library(data.table)
  library(forcats)
  library(ggh4x)
  library(ggplot2)
  library(scales)
  library(ggrepel)
})

# ---------------------------------------------------------------------------
# Shared UMAP theme
# ---------------------------------------------------------------------------
umap_theme <- theme_bw(base_size = 11) +
  theme(
    panel.grid   = element_blank(),
    axis.ticks   = element_blank(),
    axis.text    = element_blank(),
    aspect.ratio = 1,
    legend.position = "right"
  )

# ---------------------------------------------------------------------------
# UMAP coloured by sample
# ---------------------------------------------------------------------------
plot_umap_sample <- function(int_dt) {
  stopifnot(all(c("UMAP1", "UMAP2", "sample_id") %in% names(int_dt)))

  n_samples <- int_dt[, uniqueN(sample_id)]

  ggplot(int_dt[sample(nrow(int_dt))],
         aes(x = UMAP1, y = UMAP2, colour = sample_id)) +
    geom_point(size = 0.3, alpha = 0.5) +
    guides(colour = guide_legend(
      override.aes = list(size = 3, alpha = 1),
      ncol = max(1, ceiling(n_samples / 20))
    )) +
    labs(colour = "sample", title = "UMAP coloured by sample") +
    umap_theme
}

# ---------------------------------------------------------------------------
# UMAP coloured by a single metadata / categorical variable
# ---------------------------------------------------------------------------
plot_umap_by_var <- function(int_dt, var_name) {
  stopifnot(all(c("UMAP1", "UMAP2", var_name) %in% names(int_dt)))

  ggplot(int_dt[sample(nrow(int_dt))],
         aes(x = UMAP1, y = UMAP2,
             colour = as.factor(get(var_name)))) +
    geom_point(size = 0.3, alpha = 0.5) +
    guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    labs(colour = var_name,
         title  = paste("UMAP coloured by", var_name)) +
    umap_theme
}

resolve_cluster_column <- function(int_dt, cl_col) {
  if (cl_col %in% names(int_dt)) {
    return(cl_col)
  }

  if (grepl('^leiden_', cl_col)) {
    alt_col <- sub('^leiden_', 'RNA_snn_res.', cl_col)
    if (alt_col %in% names(int_dt)) {
      return(alt_col)
    }
  }

  if (grepl('^RNA_snn_res\\.', cl_col)) {
    alt_col <- sub('^RNA_snn_res\\.', 'leiden_', cl_col)
    if (alt_col %in% names(int_dt)) {
      return(alt_col)
    }
  }

  stop(sprintf('Cluster column not found: %s', cl_col))
}

# ---------------------------------------------------------------------------
# UMAP coloured by Leiden cluster assignment
# cl_col: column name, e.g. "leiden_0.5"
# ---------------------------------------------------------------------------
plot_umap_cluster <- function(int_dt, cl_col) {
  cl_col <- resolve_cluster_column(int_dt, cl_col)
  stopifnot(all(c("UMAP1", "UMAP2", cl_col) %in% names(int_dt)))

  res_val <- sub('^RNA_snn_res\\.', '', cl_col)
  res_val <- sub('^leiden_', '', res_val)
  n_cl    <- int_dt[, uniqueN(get(cl_col))]

  ggplot(int_dt[sample(nrow(int_dt))],
         aes(x = UMAP1, y = UMAP2,
             colour = as.factor(get(cl_col)))) +
    geom_point(size = 0.3, alpha = 0.5) +
    guides(colour = guide_legend(
      override.aes = list(size = 3, alpha = 1),
      ncol = max(1, ceiling(n_cl / 20))
    )) +
    labs(colour  = "cluster",
         title   = sprintf("UMAP — Leiden res %.1f (%d clusters)",
                           as.numeric(res_val), n_cl)) +
    umap_theme
}

# ---------------------------------------------------------------------------
# UMAP cell density (2-D histogram, log10 colour scale)
# Adapted from scprocess plot_umap_density().
# ---------------------------------------------------------------------------
plot_umap_density <- function(int_dt) {
  stopifnot(all(c("UMAP1", "UMAP2") %in% names(int_dt)))

  ggplot(int_dt, aes(x = UMAP1, y = UMAP2)) +
    geom_bin2d(bins = 60) +
    scale_fill_distiller(palette = "RdBu", trans = "log10") +
    labs(fill = "cells\n(log10)", title = "UMAP density") +
    umap_theme
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

# ---------------------------------------------------------------------------
# Cluster mixing diagnostics across samples
# ---------------------------------------------------------------------------
plot_cluster_entropies <- function(input_dt, batch_var = 'sample_id', what = c('norm', 'raw')) {
  what <- match.arg(what)
  stopifnot(all(c('cluster', batch_var) %in% names(input_dt)))

  ns_dt <- input_dt[, .N, by = c(batch_var, 'cluster')][
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

  cl_ls   <- sort(unique(entropy_dt$cluster))
  cl_cols <- rep(c(
    '#DC050C', '#FB8072', '#1965B0', '#7BAFDE', '#882E72', '#B17BA6',
    '#FF7F00', '#FDB462', '#E7298A', '#E78AC3', '#33A02C', '#B2DF8A'
  ), length.out = length(cl_ls))
  names(cl_cols) <- cl_ls

  ggplot(entropy_dt) +
    aes_string(x = paste0('h_cl_', what), y = paste0('max_pct_', what)) +
    geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, colour = 'grey50') +
    geom_point(aes(size = sqrt(N), fill = cluster), shape = 21) +
    geom_text_repel(aes(label = cluster), size = 3, max.overlaps = Inf) +
    scale_fill_manual(values = cl_cols, guide = 'none') +
    scale_size(
      range  = c(1, 8),
      breaks = sqrt(c(2e2, 5e2, 1e3, 2e3, 5e3, 1e4, 2e4, 5e4)),
      labels = c('200', '500', '1k', '2k', '5k', '10k', '20k', '50k')
    ) +
    expand_limits(x = 0, y = 0) +
    theme_bw(base_size = 11) +
    theme(panel.grid = element_blank()) +
    labs(
      x    = 'entropy',
      y    = sprintf('max. pct. of one %s', batch_var),
      size = 'total # cells'
    )
}

# ---------------------------------------------------------------------------
# QC metric distributions per cluster
# ---------------------------------------------------------------------------
plot_cluster_qc_distns <- function(qc_melt, clust_dt, name) {
  stopifnot(all(c('cell_id', 'cluster') %in% names(clust_dt)))
  plot_dt <- merge(qc_melt, clust_dt, by = 'cell_id')
  plot_dt[, cluster := fct_infreq(as.factor(cluster))]

  cl_levels <- levels(plot_dt$cluster)
  cl_cols <- rep(c(
    '#DC050C', '#FB8072', '#1965B0', '#7BAFDE', '#882E72', '#B17BA6',
    '#FF7F00', '#FDB462', '#E7298A', '#E78AC3', '#33A02C', '#B2DF8A'
  ), length.out = length(cl_levels))
  names(cl_cols) <- cl_levels

  log_brks   <- log10(c(1e1, 2e1, 5e1, 1e2, 2e2, 5e2, 1e3, 2e3, 5e3, 1e4, 2e4, 5e4, 1e5, 2e5, 5e5))
  log_labs   <- c('10', '20', '50', '100', '200', '500', '1k', '2k', '5k', '10k', '20k', '50k', '100k', '200k', '500k')
  logit_brks <- qlogis(c(1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 0.10, 0.30, 0.50, 0.70, 0.90, 0.97, 0.99))
  logit_labs <- c('0.01%', '0.03%', '0.1%', '0.3%', '1%', '3%', '10%', '30%', '50%', '70%', '90%', '97%', '99%')

  ggplot(plot_dt) +
    aes(x = cluster, y = qc_val, fill = cluster) +
    geom_violin(colour = NA, kernel = 'rectangular', adjust = 0.5, scale = 'width', width = 0.8) +
    scale_fill_manual(values = cl_cols, guide = 'none') +
    facet_grid(qc_full ~ ., scales = 'free_y') +
    facetted_pos_scales(
      y = list(
        qc_full == 'no. of UMIs'  ~ scale_y_continuous(breaks = log_brks, labels = log_labs),
        qc_full == 'no. of genes' ~ scale_y_continuous(breaks = log_brks, labels = log_labs),
        qc_full == 'mito pct.'    ~ scale_y_continuous(breaks = logit_brks, labels = logit_labs),
        qc_full == 'spliced pct.' ~ scale_y_continuous(breaks = logit_brks, labels = logit_labs)
      )
    ) +
    theme_bw(base_size = 11) +
    theme(
      axis.text.x      = element_text(angle = 90, hjust = 1, vjust = 0.5),
      panel.grid       = element_blank(),
      strip.background = element_rect(fill = 'white')
    ) +
    labs(x = name, y = 'QC metric')
}
