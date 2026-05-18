#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(argparse)
  library(data.table)
})

parser <- ArgumentParser(description = "Combine per-sample annotation outputs for one method")
parser$add_argument("--method_id", type = "character", required = TRUE,
  help = "Method identifier")
parser$add_argument("--cells_pattern", type = "character", required = TRUE,
  help = "Glob for per-sample annotation cell CSVs")
parser$add_argument("--export_pattern", type = "character", required = TRUE,
  help = "Glob for per-sample export metadata CSVs")
parser$add_argument("--output_cells_csv", type = "character", required = TRUE,
  help = "Combined annotation cells CSV.GZ")
parser$add_argument("--output_cluster_csv", type = "character", required = TRUE,
  help = "Combined cluster summary CSV.GZ")
parser$add_argument("--output_export_csv", type = "character", required = TRUE,
  help = "Combined export metadata CSV.GZ")

args <- parser$parse_args()

make_cluster_prediction_summary <- function(cells_dt) {
  cluster_cols <- grep("^(RNA_snn_res\\.|leiden_)", names(cells_dt), value = TRUE)
  if (length(cluster_cols) == 0) {
    return(data.table(
      cluster_col = character(),
      cluster_resolution = character(),
      cluster = character(),
      n_cells = integer(),
      top_label = character(),
      cs = numeric(),
      mean_score = numeric()
    ))
  }

  rbindlist(lapply(cluster_cols, function(cluster_col) {
    dt <- cells_dt[!is.na(get(cluster_col)) & nzchar(as.character(get(cluster_col))), .(
      cluster = as.character(get(cluster_col)),
      label = as.character(label),
      score = as.numeric(score)
    )]
    if (nrow(dt) == 0) {
      return(NULL)
    }

    counts_dt <- dt[, .N, by = .(cluster, label)]
    totals_dt <- dt[, .(
      n_cells = .N,
      mean_score = mean(score, na.rm = TRUE)
    ), by = cluster]

    summary_dt <- merge(counts_dt, totals_dt, by = "cluster", all.x = TRUE, sort = FALSE)
    summary_dt[, cs := N / n_cells]
    setorder(summary_dt, cluster, -cs, -N, label)
    summary_dt <- summary_dt[, .SD[1], by = cluster]
    summary_dt[, `:=`(
      cluster_col = cluster_col,
      cluster_resolution = sub("^(RNA_snn_res\\.|leiden_)", "", cluster_col),
      top_label = label
    )]
    summary_dt[, .(cluster_col, cluster_resolution, cluster, n_cells, top_label, cs, mean_score)]
  }), use.names = TRUE, fill = TRUE)
}

cells_paths <- Sys.glob(args$cells_pattern)
export_paths <- Sys.glob(args$export_pattern)

if (length(cells_paths) == 0) {
  stop(sprintf("No annotation cell CSVs matched pattern '%s'", args$cells_pattern), call. = FALSE)
}
if (length(export_paths) == 0) {
  stop(sprintf("No annotation export CSVs matched pattern '%s'", args$export_pattern), call. = FALSE)
}

cells_dt <- rbindlist(lapply(cells_paths, fread), use.names = TRUE, fill = TRUE)
cells_dt[, cell_id := as.character(cell_id)]
cells_dt <- unique(cells_dt, by = "cell_id")

export_dt <- rbindlist(lapply(export_paths, fread), use.names = TRUE, fill = TRUE)
export_dt[, cell_id := as.character(cell_id)]
export_dt <- unique(export_dt, by = "cell_id")

cluster_dt <- make_cluster_prediction_summary(cells_dt)

fwrite(cells_dt, args$output_cells_csv)
fwrite(cluster_dt, args$output_cluster_csv)
fwrite(export_dt, args$output_export_csv)

message(sprintf(
  "Combined %d per-sample annotation outputs for %s",
  length(cells_paths),
  args$method_id
))