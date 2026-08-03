# hvg_plots.R — plotting helpers for hvg_report.qmd
# Adapted from scprocess/code/hvgs.R to exactly replicate mlm_hvgs.Rmd output.
# Sourced into the Quarto report; not run directly.

suppressPackageStartupMessages({
  library(magrittr)
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(scales)
  library(SummarizedExperiment)
  library(DESeq2)
  library(ComplexHeatmap)
  library(circlize)
  library(viridis)
  library(strex)
})

# ---------------------------------------------------------------------------
# Compute DESeq2 VST on pseudobulk empty counts, and annotate rowData
# with EdgeR logFC/pval/padj. Mirrors scprocess calc_vst_obj().
# ---------------------------------------------------------------------------
calc_vst_obj <- function(pb, edger_dt) {
  empty_mat <- assay(pb)
  suppressMessages({
    dds <- DESeqDataSetFromMatrix(
      countData = empty_mat,
      colData   = data.frame(dummy = rep(1, ncol(empty_mat))),
      design    = ~ 1
    )
  })

  vst_obj <- tryCatch({
    tryCatch(
      vst(dds, blind = TRUE),
      error = function(e) varianceStabilizingTransformation(dds, blind = TRUE)
    )
  }, error = function(e) {
    suppressMessages({
      dds2 <- DESeqDataSetFromMatrix(
        countData = empty_mat + 1L,
        colData   = data.frame(dummy = rep(1, ncol(empty_mat))),
        design    = ~ 1
      )
    })
    vst(dds2, blind = TRUE)
  })

  # Add EdgeR logFC / pval to rowData
  edger_tmp <- copy(edger_dt)[, .(
    gene_id, log2fc.empty = logFC, pval.empty = PValue, padj.empty = FDR
  )]
  row_dt  <- rowData(vst_obj) %>%
    as.data.frame %>%
    as.data.table(keep.rownames = "gene_id")
  row_dt  <- merge(row_dt, edger_tmp, by = "gene_id", all.x = TRUE) %>%
    setkey("gene_id") %>%
    .[ rownames(vst_obj) ]
  rowData(vst_obj) <- as(row_dt, "DataFrame")

  vst_obj
}

# ---------------------------------------------------------------------------
# Scatter: normalised HVG variance vs ambient log2FC.
# Mirrors scprocess plot_hvg_stats_vs_empty_log2fc() exactly.
# ---------------------------------------------------------------------------
plot_hvg_stats_vs_empty_log2fc <- function(hvgs_dt, edger_dt, n_top = 20) {
  col_vals    <- c(hvg = "#1965B0", dirty = "#DC050C", boring = "grey")
  status_labs <- c(
    hvg   = "highly variable gene",
    dirty = "\"ambient\" gene",
    boring = "other"
  )

  plot_dt <- merge(hvgs_dt, edger_dt, by = "gene_id", all.x = TRUE) %>%
    .[is.na(FDR), is_ambient := FALSE] %>%
    .[is.na(FDR), logFC      := 0   ] %>%
    .[, .(
      gene_id, log2fc = logFC, padj = FDR,
      hv_n = highly_variable_nbatches, mean_var = variances_norm,
      is_hvg = highly_variable, is_ambient
    )] %>%
    .[, status := ifelse(is_hvg, "hvg", ifelse(is_ambient, "dirty", "boring")) %>%
        factor(levels = names(status_labs))]

  labels_dt <- plot_dt[status != "boring"] %>%
    .[order(status, -mean_var)] %>%
    .[, .SD[seq_len(min(.N, n_top))], by = status] %>%
    .[, symbol := gene_id %>% str_extract(".+(?=_ENS)")]

  ggplot(plot_dt[order(-status)]) +
    aes(x = log2fc, y = mean_var, colour = status) +
    geom_hline(yintercept = 0, linewidth = 0.1, colour = "grey20", alpha = 0.5) +
    geom_vline(xintercept = 0, linewidth = 0.1, colour = "grey20", alpha = 0.5) +
    geom_point(size = 0.2, alpha = 0.5, show.legend = TRUE) +
    geom_label_repel(
      data          = labels_dt,
      aes(label     = symbol),
      size          = 2,
      max.overlaps  = Inf,
      show.legend   = FALSE,
      label.padding = 0.1
    ) +
    scale_x_continuous(breaks = pretty_breaks()) +
    scale_y_continuous(breaks = pretty_breaks()) +
    scale_colour_manual(
      values = col_vals,
      breaks = names(status_labs),
      labels = status_labs,
      drop   = FALSE,
      guide  = guide_legend(override.aes = list(alpha = 1, size = 2))
    ) +
    theme_classic(base_size = 14) +
    theme(plot.caption = element_text(vjust = -0.5)) +
    labs(
      x      = "log2fc of \"empty\" drops vs all cells",
      y      = "mean standardized var. across samples",
      colour = "HVG classification\nof gene"
    )
}

# ---------------------------------------------------------------------------
# Volcano: ambient gene identification at varying CPM thresholds.
# Mirrors scprocess plot_ambient_gene_calculations() exactly.
# ---------------------------------------------------------------------------
plot_ambient_gene_calculations <- function(edger_dt, max_padj = 0.01, n_top = 10,
                                           min_cpm_empty = 0) {
  pval_1   <- edger_dt[FDR <= max_padj]$PValue %>% max(na.rm = TRUE)
  pval_2   <- edger_dt[FDR  > max_padj]$PValue %>% min(na.rm = TRUE)
  max_pval <- mean(c(pval_1, pval_2))

  tmp_dt <- if (min_cpm_empty > 0) {
    edger_dt[mean_logcpm.empties > log(min_cpm_empty + 1)]
  } else {
    copy(edger_dt)
  }

  labels_dt <- rbind(
    tmp_dt[logFC > 0][order(-logFC )][seq_len(min(.N, n_top))],
    tmp_dt[logFC > 0][order( PValue)][seq_len(min(.N, n_top))]
  ) %>% unique %>%
    .[, symbol := gene_id %>% str_extract(".+(?=_ENS)")]

  ggplot(tmp_dt) +
    aes(x = logFC, y = -log10(FDR)) +
    geom_hline(yintercept = -log10(max_pval), colour = "grey", linetype = "dashed") +
    geom_point(size = 0.1) +
    geom_label_repel(data = labels_dt, aes(label = symbol), size = 3, max.overlaps = Inf) +
    scale_x_continuous(breaks = pretty_breaks()) +
    scale_y_continuous(breaks = pretty_breaks()) +
    theme_classic(base_size = 14) +
    theme(plot.caption = element_text(vjust = -0.5)) +
    labs(
      x = "log2fc of \"empty\" drops vs all cells",
      y = "-log10( BH-adjusted p-value of \"empty\" drops vs all cells )"
    )
}

# ---------------------------------------------------------------------------
# Heatmap: top ambient genes across pseudobulk empty profiles.
# Mirrors scprocess plot_heatmap_of_ambient_profiles().
# Bug fix vs original: uses the vst_obj argument (not a global).
# ---------------------------------------------------------------------------
plot_heatmap_of_ambient_profiles <- function(vst_obj,
  top_var = c("var", "mean", "log2fc.empty", "pval.empty"), n_top = 40) {

  top_var <- match.arg(top_var)
  vst_mat <- assay(vst_obj)

  if (top_var == "var") {
    top_vals <- apply(vst_mat, 1, var) %>% setNames(rownames(vst_mat))
  } else if (top_var == "mean") {
    top_vals <- rowMeans(vst_mat) %>% setNames(rownames(vst_mat))
  } else if (top_var == "log2fc.empty") {
    top_vals <- rowData(vst_obj)$log2fc.empty %>% setNames(rownames(vst_obj))
  } else {  # pval.empty
    rows_dt  <- rowData(vst_obj) %>%
      as.data.frame %>%
      as.data.table(keep.rownames = "gene_id")
    rows_dt  <- rows_dt[!is.na(log2fc.empty) & log2fc.empty > 0]
    top_vals <- (rows_dt$pval.empty * -1) %>% setNames(rows_dt$gene_id)
  }

  top_gs   <- sort(top_vals, na.last = NA) %>% tail(n_top) %>% names
  plot_mat <- vst_mat[top_gs, , drop = FALSE]
  row_labs <- top_gs %>% str_extract(".+(?=_ENS)")

  mat_min  <- max(5, floor(min(plot_mat, na.rm = TRUE)))
  mat_max  <- ceiling(max(plot_mat, na.rm = TRUE))
  col_fun  <- colorRamp2(seq(mat_min, mat_max, length.out = 100), viridis(100))

  Heatmap(
    matrix               = plot_mat,
    col                  = col_fun,
    row_labels           = row_labs,
    row_names_gp         = gpar(fontsize = 8),
    column_names_gp      = gpar(fontsize = 8),
    heatmap_legend_param = list(title = "log2cpm-like\nvalues from\nDESeq2::vst"),
    row_names_side       = "left",
    column_names_side    = "top",
    na_col               = "grey"
  )
}
