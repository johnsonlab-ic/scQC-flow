# barcode_v2_plots.R

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

call_reason_levels <- c(
  "retained_upper_knee",
  "called_emptydrops",
  "excluded_emptydrops",
  "excluded_low_count"
)

call_reason_labels <- c(
  retained_upper_knee = "obvious good",
  called_emptydrops = "called by EmptyDrops",
  excluded_emptydrops = "rejected by EmptyDrops",
  excluded_low_count = "obvious bad"
)

call_reason_cols <- c(
  retained_upper_knee = "#2166ac",
  called_emptydrops = "#1b9e77",
  excluded_emptydrops = "#d95f02",
  excluded_low_count = "#bdbdbd"
)

splice_breaks <- qlogis(c(0.01, 0.03, 0.10, 0.30, 0.50, 0.70, 0.90, 0.97, 0.99))
splice_labels <- c("1%", "3%", "10%", "30%", "50%", "70%", "90%", "97%", "99%")


plot_barcode_rank_regions_v2 <- function(audit_dt) {
  plot_dt <- copy(audit_dt)[total > 0]
  plot_dt[, call_reason := factor(call_reason, levels = call_reason_levels)]

  ggplot(plot_dt, aes(x = rank, y = total, colour = call_reason, shape = final_call)) +
    geom_point(size = 0.5, alpha = 0.7, stroke = 0) +
    geom_hline(yintercept = unique(plot_dt$knee1), linetype = "dashed", colour = "#2166ac") +
    geom_hline(yintercept = unique(plot_dt$low_count_threshold), linetype = "dashed", colour = "#d95f02") +
    geom_vline(xintercept = unique(plot_dt$total_droplets_included), linetype = "dotted", colour = "#636363") +
    scale_x_log10() +
    scale_y_log10() +
    scale_colour_manual(values = call_reason_cols, labels = call_reason_labels, drop = FALSE) +
    scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), labels = c(`TRUE` = "called cell", `FALSE` = "not called")) +
    theme_classic(base_size = 11) +
    labs(x = "barcode rank", y = "total UMIs", colour = NULL, shape = NULL)
}


plot_emptydrops_ambiguous_v2 <- function(audit_dt, ed_fdr) {
  plot_dt <- copy(audit_dt)[tested_by_emptydrops == TRUE & !is.na(ED_FDR)]
  if (nrow(plot_dt) == 0) {
    return(
      ggplot() +
        annotate("text", x = 0, y = 0, label = "No ambiguous droplets were tested by EmptyDrops") +
        theme_void()
    )
  }

  plot_dt[, decision := ifelse(final_call, "called", "rejected")]

  ggplot(plot_dt, aes(x = total, y = pmax(-log10(ED_FDR), 0), colour = decision)) +
    geom_point(size = 0.7, alpha = 0.75) +
    geom_hline(yintercept = -log10(ed_fdr), linetype = "dashed", colour = "#525252") +
    scale_x_log10() +
    scale_colour_manual(values = c(called = "#1b9e77", rejected = "#d95f02")) +
    theme_classic(base_size = 11) +
    labs(x = "total UMIs", y = "-log10(EmptyDrops FDR)", colour = NULL)
}


plot_splice_diagnostics_v2 <- function(audit_dt) {
  plot_dt <- copy(audit_dt)[call_reason %in% call_reason_levels & is.finite(logit_spliced_pct)]
  plot_dt[, call_reason := factor(call_reason, levels = call_reason_levels)]

  ggplot(plot_dt, aes(x = call_reason, y = logit_spliced_pct, fill = call_reason)) +
    geom_violin(scale = "width", colour = NA, width = 0.9, alpha = 0.9) +
    scale_fill_manual(values = call_reason_cols, labels = call_reason_labels, drop = FALSE) +
    scale_y_continuous(breaks = splice_breaks, labels = splice_labels) +
    theme_classic(base_size = 11) +
    theme(axis.text.x = element_text(angle = 20, hjust = 1), legend.position = "none") +
    labs(x = NULL, y = "spliced percentage")
}