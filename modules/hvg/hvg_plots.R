# hvg_plots.R — plotting helper functions for hvg_report.qmd
# Sourced into the Quarto report; not run directly.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

# ---------------------------------------------------------------------------
# Mean–variance scatter
# Coloured by HVG status (HVG = blue, non-HVG = grey), matching mlm_hvgs
# framing: normalized variance on y-axis, log10(mean) on x-axis.
# ---------------------------------------------------------------------------
plot_mean_variance <- function(hvg_dt) {
  stopifnot(all(c("means", "variances_norm", "highly_variable") %in% names(hvg_dt)))

  n_selected <- sum(hvg_dt$highly_variable, na.rm = TRUE)
  n_total    <- nrow(hvg_dt)

  col_vals    <- c(hvg = "#1965B0", other = "grey70")
  status_labs <- c(hvg = "highly variable gene", other = "other")

  plot_dt <- copy(hvg_dt)
  plot_dt[, status := ifelse(highly_variable, "hvg", "other")]
  plot_dt[, status := factor(status, levels = names(status_labs))]

  ggplot(plot_dt[order(-status)],
         aes(x = log10(means + 1e-4), y = variances_norm, colour = status, alpha = status)) +
    geom_hline(yintercept = 0, linewidth = 0.1, colour = "grey20", alpha = 0.5) +
    geom_point(size = 0.3) +
    scale_colour_manual(values = col_vals, breaks = names(status_labs),
                        labels = status_labs, drop = FALSE,
                        guide  = guide_legend(override.aes = list(alpha = 1, size = 2))) +
    scale_alpha_manual(values = c(hvg = 0.8, other = 0.3), guide = "none") +
    scale_x_continuous(breaks = pretty_breaks()) +
    scale_y_continuous(breaks = pretty_breaks()) +
    theme_classic(base_size = 14) +
    labs(
      x      = "log10(mean expression + 0.0001)",
      y      = "mean standardized var. across samples",
      colour = "HVG classification\nof gene",
      title  = sprintf("HVG selection: %d / %d genes selected", n_selected, n_total)
    )
}

# ---------------------------------------------------------------------------
# HVG rank histogram
# Distribution of scanpy HVG ranks for the selected genes.
# ---------------------------------------------------------------------------
plot_hvg_ranks <- function(hvg_dt) {
  stopifnot("hvg_rank" %in% names(hvg_dt))

  hvg_only <- hvg_dt[highly_variable == TRUE]

  ggplot(hvg_only, aes(x = hvg_rank)) +
    geom_histogram(bins = 50, fill = "#1965B0", colour = "white") +
    scale_x_continuous(breaks = pretty_breaks()) +
    scale_y_continuous(breaks = pretty_breaks()) +
    theme_classic(base_size = 14) +
    labs(
      x     = "HVG rank",
      y     = "count",
      title = sprintf("Distribution of HVG ranks (n = %d)", nrow(hvg_only))
    )
}

# ---------------------------------------------------------------------------
# Batch prevalence bar chart
# For each value of highly_variable_nbatches, the number of selected HVGs.
# Analogous to mlm_hvgs section 2: shows how "sticky" each HVG is across
# batches — genes called in many batches are robust HVGs.
# ---------------------------------------------------------------------------
plot_hvg_nbatches <- function(hvg_dt) {
  stopifnot("highly_variable_nbatches" %in% names(hvg_dt))

  nb_dt <- hvg_dt[highly_variable == TRUE & !is.na(highly_variable_nbatches),
                  .N, by = highly_variable_nbatches][order(highly_variable_nbatches)]

  ggplot(nb_dt, aes(x = highly_variable_nbatches, y = N)) +
    geom_col(fill = "#1965B0") +
    scale_x_continuous(breaks = pretty(nb_dt$highly_variable_nbatches)) +
    scale_y_continuous(breaks = pretty_breaks()) +
    theme_classic(base_size = 14) +
    labs(
      x     = "number of batches where selected as HVG",
      y     = "count of genes",
      title = "Batch prevalence of selected HVGs"
    )
}

# ---------------------------------------------------------------------------
# Top-genes horizontal barplot
# Shows the top n_top genes ranked by `rank_by` column as a horizontal bar
# chart with gene labels — adapted analogue of the ambient-profile heatmaps
# in mlm_hvgs (var / mean orderings).
# ---------------------------------------------------------------------------
plot_top_genes_bar <- function(hvg_dt, rank_by = c("variances_norm", "means"),
                               n_top = 40) {
  rank_by <- match.arg(rank_by)
  stopifnot(rank_by %in% names(hvg_dt), "highly_variable" %in% names(hvg_dt))

  plot_dt <- copy(hvg_dt)[order(-get(rank_by))][seq_len(min(n_top, .N))]
  plot_dt[, gene_id := factor(gene_id, levels = rev(gene_id))]
  plot_dt[, hvg_lab := ifelse(highly_variable, "HVG", "non-HVG")]

  x_lab <- switch(rank_by,
    variances_norm = "normalised variance",
    means          = "mean expression (log-normalised scale)"
  )
  title_str <- switch(rank_by,
    variances_norm = sprintf("Top %d genes by normalised variance", min(n_top, nrow(plot_dt))),
    means          = sprintf("Top %d genes by mean expression",     min(n_top, nrow(plot_dt)))
  )

  ggplot(plot_dt, aes(x = get(rank_by), y = gene_id, fill = hvg_lab)) +
    geom_col() +
    scale_fill_manual(values = c("HVG" = "#1965B0", "non-HVG" = "grey70"),
                      guide  = guide_legend(reverse = TRUE)) +
    scale_x_continuous(breaks = pretty_breaks(), expand = expansion(mult = c(0, 0.05))) +
    theme_classic(base_size = 11) +
    theme(axis.text.y = element_text(size = 7)) +
    labs(x = x_lab, y = NULL, fill = NULL, title = title_str)
}
