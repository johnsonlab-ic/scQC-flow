# qc_plots.R — plotting and summary helper functions for qc_report.qmd
# Ported from scprocess code/qc.R

suppressPackageStartupMessages({
  library(magrittr)
  library(data.table)
  library(forcats)
  library(assertthat)
  library(stringr)
  library(ggplot2)
  library(scales)
  library(patchwork)
  library(ggh4x)
  library(UpSetR)
})

# ---------------------------------------------------------------------------
# Axis breaks and labels for transformed metrics
# ---------------------------------------------------------------------------
n_brks      <- c(1e1, 2e1, 5e1, 1e2, 2e2, 5e2,
                 1e3, 2e3, 5e3, 1e4, 2e4, 5e4, 1e5) %>% log10
n_labs      <- c("10", "20", "50", "100", "200", "500",
                 "1k", "2k", "5k", "10k", "20k", "50k", "100k")
log_brks    <- c(1e1, 3e1, 1e2, 3e2, 1e3, 3e3, 1e4, 3e4, 1e5, 3e5) %>% log10
log_labs    <- c("10", "30", "100", "300", "1k", "3k", "10k", "30k", "100k", "300k")
mito_brks   <- c(1e-5, 1e-4, 1e-3, 1e-2, 0.10, 0.50, 0.90, 0.99, 0.999) %>% qlogis
mito_labs   <- c("0.001%", "0.01%", "0.1%", "1%", "10%", "50%", "90%", "99%", "99.9%")
splice_brks <- c(1e-4, 1e-3, 1e-2, 0.10, 0.50, 0.90, 0.99, 0.999) %>% qlogis
splice_labs <- c("0.01%", "0.1%", "1%", "10%", "50%", "90%", "99%", "99.9%")

# ---------------------------------------------------------------------------
# Build cuts data.table from threshold parameters
# ---------------------------------------------------------------------------
make_cuts_dt <- function(b_lvls, min_counts, min_feats,
                         min_mito, max_mito, min_splice, max_splice) {
  rows <- lapply(b_lvls, function(b) {
    data.table(
      batch_var         = b,
      log_counts_min    = log10(min_counts + 1),
      log_feats_min     = log10(min_feats + 1),
      logit_mito_min    = qlogis(min_mito),
      logit_mito_max    = qlogis(max_mito),
      logit_spliced_min = qlogis(min_splice),
      logit_spliced_max = qlogis(max_splice),
      n_cells_min       = -Inf
    )
  }) %>% rbindlist

  cuts_dt <- rows %>%
    melt(id = "batch_var", variable.name = "cut_var", value.name = "cut_point") %>%
    .[, qc_var := str_extract(cut_var, "^(.+)(?=_(min|max))")] %>%
    .[, minmax := str_extract(cut_var, "(min|max)")] %>%
    dcast(batch_var + qc_var ~ minmax, value.var = "cut_point") %>%
    .[is.na(min), min := -Inf] %>%
    .[is.na(max), max := Inf]

  return(cuts_dt)
}

# ---------------------------------------------------------------------------
# Build hard-threshold cuts data.table (defines pass zone for hard thresholds)
# ---------------------------------------------------------------------------
make_hard_cuts_dt <- function(b_lvls, hard_min_counts, hard_min_feats, hard_max_mito) {
  rows <- lapply(b_lvls, function(b) {
    data.table(
      batch_var      = b,
      log_counts_min = log10(hard_min_counts + 1),
      log_feats_min  = log10(hard_min_feats + 1),
      logit_mito_max = qlogis(hard_max_mito)
    )
  }) %>% rbindlist

  hard_cuts_dt <- rows %>%
    melt(id = "batch_var", variable.name = "cut_var", value.name = "cut_point") %>%
    .[, qc_var := str_extract(cut_var, "^(.+)(?=_(min|max))")] %>%
    .[, minmax := str_extract(cut_var, "(min|max)")] %>%
    dcast(batch_var + qc_var ~ minmax, value.var = "cut_point") %>%
    .[is.na(min), min := -Inf] %>%
    .[is.na(max), max := Inf]

  return(hard_cuts_dt)
}

# ---------------------------------------------------------------------------
# Violin + marginal cell-count plot (mirrors scprocess plot_qc_ranges_marginals)
# ---------------------------------------------------------------------------
plot_qc_ranges_marginals <- function(qc_input, b_lvls, qc_names, qc_lu, cuts_dt,
                                     hard_cuts_dt = NULL,
                                     batch_label = "sample_id") {
  tmp_names <- intersect(qc_names, colnames(qc_input))
  qc_melt   <- copy(qc_input) %>%
    melt(measure = tmp_names, val = "qc_val", var = "qc_var") %>%
    .[, qc_full := qc_lu[as.character(qc_var)]] %>%
    .[, qc_var  := factor(qc_var, levels = tmp_names)] %>%
    .[, qc_full := fct_reorder(qc_full, as.integer(qc_var))]

  hlines_dt <- cuts_dt %>% copy %>%
    .[, qc_full := qc_lu[qc_var]] %>%
    .[, qc_var  := factor(qc_var, levels = qc_names)] %>%
    .[, qc_full := fct_reorder(qc_full, as.integer(qc_var))]

  # medians per sample
  qc_meds <- qc_melt %>%
    .[, .(log10_N = log10(.N),
          q50     = median(qc_val, na.rm = TRUE)),
      by = c("batch_var", "qc_var", "qc_full")]

  # n_cells bar data
  n_dt <- qc_meds[, .(batch_var, `no. of cells` = log10_N)] %>% unique %>%
    melt.data.table(id = "batch_var", var = "var", val = "value") %>%
    .[, batch_var := factor(batch_var, levels = rev(b_lvls))]

  n_lims   <- c(min(n_dt$value) + log10(0.5), max(n_dt$value) + log10(2))
  qc_melt  <- qc_melt[, batch_var := factor(batch_var, levels = rev(b_lvls))]
  qc_meds  <- qc_meds[, batch_var := factor(batch_var, levels = rev(b_lvls))]
  hlines_dt <- hlines_dt[, batch_var := factor(batch_var, levels = rev(b_lvls))]

  # Prepare hard exclusion zone data (lower and upper strips per variable)
  hard_lower_dt <- NULL
  hard_upper_dt <- NULL
  if (!is.null(hard_cuts_dt)) {
    hard_hlines_dt <- hard_cuts_dt %>% copy %>%
      .[, qc_full := qc_lu[qc_var]] %>%
      .[, qc_var  := factor(qc_var, levels = qc_names)] %>%
      .[, qc_full := fct_reorder(qc_full, as.integer(qc_var))] %>%
      .[, batch_var := factor(batch_var, levels = rev(b_lvls))]
    hard_lower_dt <- hard_hlines_dt[is.finite(min)]
    hard_upper_dt <- hard_hlines_dt[is.finite(max)]
  }

  # cell-count panel
  g_n <- ggplot() +
    geom_rect(data = hlines_dt[qc_var == "n_cells"],
      aes(xmin = as.integer(batch_var) - 0.4, xmax = as.integer(batch_var) + 0.4,
          ymin = min, ymax = max),
      fill = "grey80", colour = NA, alpha = 0.2) +
    geom_point(data = n_dt, aes(y = value, x = as.integer(batch_var)),
      size = 4, shape = 21, fill = "grey40") +
    scale_x_continuous(breaks = seq_len(length(b_lvls)), labels = rev(b_lvls)) +
    facet_grid(. ~ var, scales = "free", space = "free_y") +
    scale_y_continuous(breaks = n_brks, labels = n_labs) +
    expand_limits(y = n_lims) +
    coord_flip(xlim = c(0.5, length(b_lvls) + 0.5), expand = FALSE) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          strip.text.y = element_blank()) +
    labs(y = NULL, x = batch_label)

  # violin panels — hard exclusion zones (darker) drawn first, soft pass zone on top
  g_violin <- ggplot() +
    (if (!is.null(hard_lower_dt) && nrow(hard_lower_dt) > 0)
      geom_rect(data = hard_lower_dt,
        aes(xmin = as.integer(batch_var) - 0.4, xmax = as.integer(batch_var) + 0.4,
            ymin = -Inf, ymax = min),
        fill = "grey55", colour = NA, alpha = 0.3)
    else NULL) +
    (if (!is.null(hard_upper_dt) && nrow(hard_upper_dt) > 0)
      geom_rect(data = hard_upper_dt,
        aes(xmin = as.integer(batch_var) - 0.4, xmax = as.integer(batch_var) + 0.4,
            ymin = max, ymax = Inf),
        fill = "grey55", colour = NA, alpha = 0.3)
    else NULL) +
    geom_rect(data = hlines_dt[qc_var != "n_cells"],
      aes(xmin = as.integer(batch_var) - 0.4, xmax = as.integer(batch_var) + 0.4,
          ymin = min, ymax = max),
      fill = "grey80", colour = NA, alpha = 0.2) +
    geom_violin(data = qc_melt[!is.na(qc_val)],
      aes(x = batch_var, y = qc_val), colour = NA, fill = "grey40",
      kernel = "rectangular", adjust = 0.1, scale = "width", width = 0.8) +
    facet_grid(. ~ qc_full, scales = "free", space = "free_y") +
    scale_x_discrete(breaks = levels(qc_meds$batch_var), drop = FALSE) +
    facetted_pos_scales(y = list(
      qc_full == "no. of UMIs"  ~ scale_y_continuous(breaks = log_brks, labels = log_labs),
      qc_full == "no. of genes" ~ scale_y_continuous(breaks = log_brks, labels = log_labs),
      qc_full == "mito. pct."   ~ scale_y_continuous(breaks = mito_brks, labels = mito_labs),
      qc_full == "spliced pct." ~ scale_y_continuous(breaks = splice_brks, labels = splice_labs)
    )) +
    coord_flip() +
    theme_classic() +
    theme(axis.ticks.y = element_blank(),
          axis.text.y  = element_blank(),
          axis.text.x  = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    labs(x = NULL, y = NULL)

  g <- g_n + g_violin + plot_layout(widths = c(1, 5))
  return(g)
}

# ---------------------------------------------------------------------------
# Cell-level pairwise scatter plots (mirrors scprocess plot_qc_metric_scatter)
# ---------------------------------------------------------------------------
plot_qc_metric_scatter <- function(dt, qc_names, qc_lu, cuts_one, name) {
  qc_names <- intersect(qc_names, colnames(dt))
  melt_dt  <- dt %>%
    melt(measure = qc_names, value.name = "qc_val", variable.name = "qc_var") %>%
    .[, qc_full := qc_lu[as.character(qc_var)]] %>%
    .[, qc_var  := factor(qc_var, levels = qc_names)] %>%
    .[, qc_full := fct_reorder(qc_full, as.integer(qc_var))]

  plot_dt <- merge(
    melt_dt[, .(cell_id, qc_x = qc_full, val_x = qc_val)],
    melt_dt[, .(cell_id, qc_y = qc_full, val_y = qc_val)],
    by = "cell_id", allow.cartesian = TRUE
  ) %>% .[as.integer(qc_x) > as.integer(qc_y)]

  rects_dt <- merge(
    cuts_one[qc_var %in% qc_names] %>%
      .[, .(dummy = "d", qc_var.y = qc_var, min.y = min, max.y = max)],
    cuts_one[qc_var %in% qc_names] %>%
      .[, .(dummy = "d", qc_var.x = qc_var, min.x = min, max.x = max)],
    by = "dummy", allow.cartesian = TRUE
  ) %>%
    .[, qc_x     := qc_lu[qc_var.x]] %>%
    .[, qc_var.x := factor(qc_var.x, levels = qc_names)] %>%
    .[, qc_x     := fct_reorder(qc_x, as.integer(qc_var.x))] %>%
    .[, qc_y     := qc_lu[qc_var.y]] %>%
    .[, qc_var.y := factor(qc_var.y, levels = qc_names)] %>%
    .[, qc_y     := fct_reorder(qc_y, as.integer(qc_var.y))] %>%
    .[as.integer(qc_x) > as.integer(qc_y)]

  g <- ggplot() +
    geom_rect(data = rects_dt,
      aes(xmin = min.x, xmax = max.x, ymin = min.y, ymax = max.y),
      fill = "grey80", colour = NA, alpha = 0.4) +
    geom_bin2d(data = plot_dt, aes(x = val_x, y = val_y), bins = 40) +
    scale_fill_distiller(palette = "RdBu", trans = "log10") +
    facet_grid(qc_y ~ qc_x, scales = "free") +
    theme_classic() +
    theme(panel.grid = element_blank(),
          strip.background = element_rect(fill = "white")) +
    facetted_pos_scales(
      x = list(
        qc_x == "no. of UMIs"  ~ scale_x_continuous(breaks = log_brks, labels = log_labs),
        qc_x == "no. of genes" ~ scale_x_continuous(breaks = log_brks, labels = log_labs),
        qc_x == "mito. pct."   ~ scale_x_continuous(breaks = mito_brks, labels = mito_labs),
        qc_x == "spliced pct." ~ scale_x_continuous(breaks = splice_brks, labels = splice_labs)
      ),
      y = list(
        qc_y == "no. of UMIs"  ~ scale_y_continuous(breaks = log_brks, labels = log_labs),
        qc_y == "no. of genes" ~ scale_y_continuous(breaks = log_brks, labels = log_labs),
        qc_y == "mito. pct."   ~ scale_y_continuous(breaks = mito_brks, labels = mito_labs),
        qc_y == "spliced pct." ~ scale_y_continuous(breaks = splice_brks, labels = splice_labs)
      )
    ) +
    labs(x = "QC metric 1", y = "QC metric 2", fill = "no. of cells", title = name)

  return(g)
}

# ---------------------------------------------------------------------------
# Per-sample MAD threshold violin plot
# ---------------------------------------------------------------------------
plot_mad_thresholds <- function(qc_input, b_lvls, qc_lu, nmads) {
  qc_names_mad <- c("log_counts", "log_feats", "logit_mito", "logit_spliced")

  flat_dt <- qc_input[keep_flat == TRUE]
  flat_dt[, batch_var := factor(batch_var, levels = rev(b_lvls))]

  bound_cols <- c("batch_var",
    paste0("mad_lo_", qc_names_mad),
    paste0("mad_hi_", qc_names_mad))
  bounds_wide <- unique(flat_dt[, bound_cols, with = FALSE])
  bounds_wide[, batch_pos := as.integer(factor(batch_var, levels = rev(b_lvls)))]

  bounds_lo <- melt(bounds_wide,
    id.vars = c("batch_var", "batch_pos"),
    measure.vars = paste0("mad_lo_", qc_names_mad),
    variable.name = "metric", value.name = "bound_val")
  bounds_lo[, qc_var     := sub("^mad_lo_", "", metric)]
  bounds_lo[, bound_type := "lower"]

  bounds_hi <- melt(bounds_wide,
    id.vars = c("batch_var", "batch_pos"),
    measure.vars = paste0("mad_hi_", qc_names_mad),
    variable.name = "metric", value.name = "bound_val")
  bounds_hi[, qc_var     := sub("^mad_hi_", "", metric)]
  bounds_hi[, bound_type := "upper"]

  bounds_long <- rbind(
    bounds_lo[, .(batch_var, batch_pos, qc_var, bound_val, bound_type)],
    bounds_hi[, .(batch_var, batch_pos, qc_var, bound_val, bound_type)]
  )
  bounds_long[, qc_full := qc_lu[qc_var]]
  bounds_long[, qc_var  := factor(qc_var, levels = qc_names_mad)]
  bounds_long[, qc_full := fct_reorder(qc_full, as.integer(qc_var))]

  flat_melt <- melt(flat_dt,
    measure.vars = qc_names_mad, variable.name = "qc_var", value.name = "qc_val")
  flat_melt[, qc_full   := qc_lu[as.character(qc_var)]]
  flat_melt[, qc_var    := factor(qc_var, levels = qc_names_mad)]
  flat_melt[, qc_full   := fct_reorder(qc_full, as.integer(qc_var))]
  flat_melt[, batch_pos := as.integer(batch_var)]

  ggplot() +
    geom_violin(data = flat_melt[!is.na(qc_val)],
      aes(x = batch_pos, y = qc_val, group = batch_pos),
      colour = NA, fill = "grey40",
      kernel = "rectangular", adjust = 0.1, scale = "width", width = 0.8) +
    geom_segment(data = bounds_long,
      aes(x = batch_pos - 0.35, xend = batch_pos + 0.35,
          y = bound_val, yend = bound_val,
          colour = bound_type),
      linewidth = 1.0) +
    scale_colour_manual(
      values = c(lower = "#D73027", upper = "#FC8D59"),
      labels = c(lower = sprintf("lower (median − %g×MAD)", nmads),
                 upper = sprintf("upper (median + %g×MAD)", nmads))) +
    scale_x_continuous(
      breaks = seq_len(length(b_lvls)),
      labels = rev(b_lvls)) +
    facet_grid(. ~ qc_full, scales = "free") +
    facetted_pos_scales(y = list(
      qc_full == "no. of UMIs"  ~ scale_y_continuous(breaks = log_brks, labels = log_labs),
      qc_full == "no. of genes" ~ scale_y_continuous(breaks = log_brks, labels = log_labs),
      qc_full == "mito. pct."   ~ scale_y_continuous(breaks = mito_brks, labels = mito_labs),
      qc_full == "spliced pct." ~ scale_y_continuous(breaks = splice_brks, labels = splice_labs)
    )) +
    coord_flip() +
    theme_classic() +
    theme(
      axis.text.y     = element_text(size = 7),
      axis.text.x     = element_text(angle = 90, hjust = 1, vjust = 0.5),
      legend.position = "bottom"
    ) +
    labs(x = NULL, y = NULL, colour = NULL) +
    guides(colour = guide_legend(override.aes = list(linewidth = 2)))
}

# ---------------------------------------------------------------------------
# Per-metric MAD removal breakdown (for summary table)
# ---------------------------------------------------------------------------
calc_mad_exclusions <- function(qc_dt) {
  flat_dt <- qc_dt[keep_flat == TRUE]

  flat_dt[, `:=`(
    mad_fail_counts = log_counts    < mad_lo_log_counts    | log_counts    > mad_hi_log_counts,
    mad_fail_feats  = log_feats     < mad_lo_log_feats     | log_feats     > mad_hi_log_feats,
    mad_fail_mito   = logit_mito    < mad_lo_logit_mito    | logit_mito    > mad_hi_logit_mito,
    mad_fail_splice = logit_spliced < mad_lo_logit_spliced | logit_spliced > mad_hi_logit_spliced
  )]

  out <- flat_dt[, .(
    `N after flat`    = .N,
    `N removed by MAD`= sum(keep_flat & !keep),
    `% MAD (UMIs)`    = round(100 * mean(mad_fail_counts,  na.rm = TRUE), 1),
    `% MAD (genes)`   = round(100 * mean(mad_fail_feats,   na.rm = TRUE), 1),
    `% MAD (mito)`    = round(100 * mean(mad_fail_mito,    na.rm = TRUE), 1),
    `% MAD (splice)`  = round(100 * mean(mad_fail_splice,  na.rm = TRUE), 1)
  ), by = batch_var]

  setnames(out, "batch_var", "sample_id")
  return(out)
}

# ---------------------------------------------------------------------------
# QC summary table (mirrors scprocess calc_qc_summary)
# ---------------------------------------------------------------------------
calc_qc_summary <- function(qc_dt, kept_dt, cuts_dt, qc_lu, mad_filter = FALSE) {
  tbl_tmp <- merge(
    qc_dt[, .(n_pre_QC = .N), by = .(batch_var)],
    kept_dt[, .(n_post_QC = .N), by = .(batch_var)],
    by = "batch_var", all = TRUE
  ) %>%
    .[is.na(n_post_QC), n_post_QC := 0] %>%
    .[, n_excluded   := n_pre_QC - n_post_QC] %>%
    .[, pct_excluded := round(100 * (1 - n_post_QC / n_pre_QC), 1)] %>%
    .[, batch_excluded := n_post_QC == 0]

  # per-threshold exclusion breakdown
  exclude_dt <- .calc_exclusions_by_qc_var(qc_dt, cuts_dt, qc_lu)

  tbl_pretty <- tbl_tmp %>%
    .[, .(batch_var, excluded = batch_excluded,
          `N pre-QC` = n_pre_QC, `N post-QC` = n_post_QC,
          `N excluded` = n_excluded, `pct. excluded` = pct_excluded)] %>%
    merge(exclude_dt, by = "batch_var") %>%
    setnames("batch_var", "sample_id") %>%
    .[order(-`pct. excluded`, -`N post-QC`)]

  if (mad_filter) {
    mad_cols   <- calc_mad_exclusions(qc_dt)
    tbl_pretty <- merge(tbl_pretty, mad_cols, by = "sample_id", all.x = TRUE)
  }

  return(tbl_pretty)
}

.calc_exclusions_by_qc_var <- function(qc_dt, cuts_dt, qc_lu) {
  cut_vars    <- intersect(cuts_dt$qc_var, colnames(qc_dt))
  qc_tmp      <- qc_dt[, c("batch_var", "cell_id", cut_vars), with = FALSE]
  all_batches <- unique(qc_tmp$batch_var)

  exclude_tmp <- lapply(cut_vars, function(vv) {
    lapply(all_batches, function(bb) {
      cuts_tmp <- cuts_dt[(batch_var == bb) & (qc_var == vv)]
      spec     <- c(cuts_tmp$min, cuts_tmp$max)
      qc_tmp[batch_var == bb] %>%
        .[, .(.N, N_keep = sum(get(vv) %between% spec)), by = .(batch_var)] %>%
        .[, qc_var    := vv] %>%
        .[, col_title := sprintf("pct. excluded by %s", qc_lu[[vv]])]
    }) %>% rbindlist
  }) %>% rbindlist %>%
    .[, pct_exc   := round((N - N_keep) / N * 100, 1)] %>%
    .[, qc_var    := factor(qc_var, levels = names(qc_lu))] %>%
    .[, col_title := fct_reorder(col_title, as.integer(qc_var))]

  exclude_dt <- exclude_tmp %>%
    dcast(batch_var ~ col_title, value.var = "pct_exc")

  return(exclude_dt)
}

# ---------------------------------------------------------------------------
# UpSet plot of exclusion reasons (mirrors scprocess plot_upset_of_exclusions)
# ---------------------------------------------------------------------------
plot_upset_of_exclusions <- function(qc_tmp, qc_names, qc_lu, cuts_dt) {
  var_lu  <- c(
    log_counts    = "umis",
    log_feats     = "genes",
    logit_mito    = "mito",
    logit_spliced = "spliced"
  )
  eps <- 1e-10

  # Doublets and hard-threshold failures — captured before restricting to hard-pass
  extra_ls <- list()
  if ("scdbl_class" %in% colnames(qc_tmp)) {
    dbl_cells <- qc_tmp[scdbl_class == "doublet"]$cell_id
    if (length(dbl_cells) > 0) extra_ls[["doublet"]] <- dbl_cells
    hard_fail <- qc_tmp[keep_hard == FALSE & scdbl_class != "doublet"]$cell_id
  } else {
    hard_fail <- qc_tmp[keep_hard == FALSE]$cell_id
  }
  if (length(hard_fail) > 0) extra_ls[["hard_fail"]] <- hard_fail

  # Soft threshold exclusions among hard-passing singlets only
  hard_pass <- qc_tmp[keep_hard == TRUE]

  cuts_tmp <- cuts_dt[batch_var == unique(qc_tmp$batch_var)] %>%
    .[qc_var != "n_cells"]
  cut_vars <- unique(cuts_tmp$qc_var)

  tmp_ls <- lapply(cut_vars, function(cc) {
    cut_spec  <- c(cuts_tmp[qc_var == cc]$min, cuts_tmp[qc_var == cc]$max)
    exc_cells <- list(
      low  = hard_pass[get(cc) + eps < cut_spec[1]]$cell_id,
      high = hard_pass[get(cc) - eps > cut_spec[2]]$cell_id
    )
    names(exc_cells) <- paste(names(exc_cells), var_lu[[cc]], sep = "_")
    exc_cells[sapply(exc_cells, length) > 0]
  }) %>% do.call(c, .)

  ok_cells <- qc_tmp[keep == TRUE]$cell_id
  upset_ls <- c(extra_ls, tmp_ls, list(passed_qc = ok_cells))

  upset_dt <- names(upset_ls) %>%
    lapply(function(nn) data.table(set = nn, cell_id = upset_ls[[nn]])) %>%
    rbindlist %>%
    .[, dummy := 1] %>%
    dcast(cell_id ~ set, value.var = "dummy", fill = 0)

  upset(upset_dt, sets = colnames(upset_dt)[-1], order.by = "freq",
        mb.ratio = c(0.7, 0.3), sets.bar.color = "#FB8072")
}
