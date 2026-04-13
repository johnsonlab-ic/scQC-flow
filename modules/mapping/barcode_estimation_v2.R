#!/usr/bin/env Rscript
# barcode_estimation_v2.R
#
# Standalone barcode caller for alevin-fry outputs.
#
# Outputs:
#   - af_counts_mat_v2.h5               raw stacked S+U+A matrix (10x v3 format)
#   - barcode_audit_v2_<id>.csv.gz      per-barcode audit table
#   - barcode_summary_v2_<id>.csv       per-sample decision summary
#   - cell_barcodes_v2_<id>.csv         final called barcodes
#   - ambient_params_v2_<id>.env        future-compatible thresholds / counts

suppressPackageStartupMessages({
  library(argparse)
  library(magrittr)
  library(fishpond)
  library(SingleCellExperiment)
  library(DropletUtils)
  library(BiocParallel)
  library(data.table)
  library(parallel)
  library(strex)
})


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

parser <- ArgumentParser(description = "Standalone barcode caller for alevin-fry outputs")
parser$add_argument("--sample_id", type = "character", required = TRUE)
parser$add_argument("--quant_dir", type = "character", required = TRUE)
parser$add_argument("--h5_out", type = "character", default = "af_counts_mat_v2.h5")
parser$add_argument("--audit_out", type = "character", default = "barcode_audit_v2.csv.gz")
parser$add_argument("--summary_out", type = "character", default = "barcode_summary_v2.csv")
parser$add_argument("--barcodes_out", type = "character", default = "cell_barcodes_v2.csv")
parser$add_argument("--ambient_env_out", type = "character", default = "ambient_params_v2.env")
parser$add_argument("--min_umis_empty", type = "integer", default = 5)
parser$add_argument("--niters", type = "integer", default = 10000)
parser$add_argument("--ed_fdr", type = "double", default = 0.001)
parser$add_argument("--splice_context", type = "character", default = "snrna")
parser$add_argument("--low_count_strategy", type = "character", default = "shin2")

args <- parser$parse_args()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

run_barcode_estimation_v2 <- function(args) {

  sample_id          <- args$sample_id
  quant_dir          <- args$quant_dir
  h5_out             <- args$h5_out
  audit_out          <- args$audit_out
  summary_out        <- args$summary_out
  barcodes_out       <- args$barcodes_out
  ambient_env_out    <- args$ambient_env_out
  min_umis_empty     <- args$min_umis_empty
  niters             <- args$niters
  ed_fdr             <- args$ed_fdr
  splice_context     <- tolower(args$splice_context)
  low_count_strategy <- args$low_count_strategy

  if (!(splice_context %in% c("snrna", "scrna", "auto"))) {
    stop("splice_context must be one of: snrna, scrna, auto")
  }
  if (!(low_count_strategy %in% c("knee2", "shin2"))) {
    stop("low_count_strategy must be one of: knee2, shin2")
  }

  message("=== BARCODE_ESTIMATION_V2 ===")
  message("Sample:          ", sample_id)
  message("Quant dir:       ", quant_dir)
  message("min_umis_empty:  ", min_umis_empty)
  message("EmptyDrops FDR:  ", ed_fdr)
  message("EmptyDrops iter: ", niters)
  message("Splice context:  ", splice_context)

  sce <- loadFry(
    quant_dir,
    outputFormat = list(S = c("S"), U = c("U"), A = c("A"))
  )

  total_counts <- assay(sce, "S") + assay(sce, "U") + assay(sce, "A")
  n_barcodes_input <- ncol(total_counts)

  mat <- assayNames(sce) %>% lapply(function(mode_name) {
    this_mat <- assay(sce, mode_name)
    rownames(this_mat) <- paste0(rownames(this_mat), "_", mode_name)
    this_mat
  }) %>% do.call(rbind, .)

  keep_bcs <- Matrix::colSums(mat) > 0
  mat <- mat[, keep_bcs, drop = FALSE]
  total_counts <- total_counts[, keep_bcs, drop = FALSE]
  sce <- sce[, keep_bcs]
  n_zero_removed <- n_barcodes_input - ncol(total_counts)

  message("Input barcodes:      ", n_barcodes_input)
  message("Zero-count removed:  ", n_zero_removed)
  message("Non-zero barcodes:   ", ncol(total_counts))

  write10xCounts(h5_out, mat, version = "3", overwrite = TRUE)
  message("Written raw H5: ", h5_out)

  barcode_extra_dt <- data.table(
    barcode          = colnames(total_counts),
    detected_genes   = Matrix::colSums(total_counts > 0),
    spliced          = Matrix::colSums(assay(sce, "S")),
    unspliced        = Matrix::colSums(assay(sce, "U")),
    ambiguous_counts = Matrix::colSums(assay(sce, "A"))
  )
  barcode_extra_dt[, spliced_pct := (spliced + 1) / (spliced + unspliced + 2)]
  barcode_extra_dt[, logit_spliced_pct := qlogis(pmin(pmax(spliced_pct, 1e-6), 1 - 1e-6))]

  thresh_ls <- calc_ambient_params_v2(
    split_mat = mat,
    run = sample_id,
    min_umis_empty = min_umis_empty,
    low_count_threshold = low_count_strategy,
    run_var = "sample_id"
  )

  audit_dt <- merge(thresh_ls$ranks_dt, barcode_extra_dt, by = "barcode", all.x = TRUE)
  setorder(audit_dt, rank)

  audit_dt[, decision_zone := fifelse(
    total >= knee1,
    "obvious_good",
    fifelse(total < low_count_threshold, "obvious_bad", "ambiguous")
  )]
  audit_dt[, tested_by_emptydrops := decision_zone == "ambiguous"]

  knee1_val <- as.integer(unique(audit_dt$knee1))
  empty_bcs <- audit_dt[in_empty_plateau == TRUE, barcode]
  empty_idx <- which(colnames(total_counts) %in% empty_bcs)
  if (length(empty_idx) == 0) {
    empty_idx <- which(audit_dt$decision_zone == "obvious_bad")
  }
  if (length(empty_idx) == 0) {
    fallback_n <- min(100L, nrow(audit_dt))
    empty_idx <- order(Matrix::colSums(total_counts))[seq_len(fallback_n)]
  }

  n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1"))
  bpparam <- MulticoreParam(workers = max(1L, n_cores), progressbar = FALSE)

  message("Running EmptyDrops with ", length(empty_idx), " background barcodes and retain=", knee1_val)
  ed_res <- emptyDrops(
    m = total_counts,
    niters = niters,
    BPPARAM = bpparam,
    known.empty = empty_idx,
    retain = knee1_val
  )

  ed_dt <- as.data.table(as.data.frame(ed_res), keep.rownames = TRUE)
  setnames(ed_dt, "rn", "barcode")
  setnames(
    ed_dt,
    c("Total", "LogProb", "PValue", "Limited", "FDR"),
    c("ED_Total", "ED_LogProb", "ED_PValue", "ED_Limited", "ED_FDR")
  )

  audit_dt <- merge(
    audit_dt,
    ed_dt[, .(barcode, ED_Total, ED_LogProb, ED_PValue, ED_Limited, ED_FDR)],
    by = "barcode",
    all.x = TRUE
  )
  setorder(audit_dt, rank)

  audit_dt[, final_call := decision_zone == "obvious_good" |
    (tested_by_emptydrops & !is.na(ED_FDR) & ED_FDR <= ed_fdr)]

  audit_dt[, call_reason := fcase(
    decision_zone == "obvious_good", "retained_upper_knee",
    decision_zone == "obvious_bad", "excluded_low_count",
    tested_by_emptydrops & final_call, "called_emptydrops",
    tested_by_emptydrops & !final_call, "excluded_emptydrops",
    default = "unclassified"
  )]

  audit_dt <- annotate_splice_flags(audit_dt, splice_context)

  summary_dt <- make_summary_dt(
    audit_dt = audit_dt,
    sample_id = sample_id,
    total_barcodes_input = n_barcodes_input,
    zero_total_removed = n_zero_removed,
    ed_fdr = ed_fdr,
    niters = niters,
    splice_context = splice_context,
    low_count_strategy = low_count_strategy,
    min_umis_empty = min_umis_empty
  )

  fwrite(audit_dt, file = audit_out, compress = "gzip")
  fwrite(summary_dt, file = summary_out)
  fwrite(audit_dt[final_call == TRUE, .(barcode)], file = barcodes_out)

  writeLines(c(
    sprintf("RUN=%s", sample_id),
    sprintf("CB_EXPECTED_CELLS_V2=%s", summary_dt$final_cells),
    sprintf("CB_LOW_COUNT_THRESHOLD_V2=%s", summary_dt$low_count_threshold),
    sprintf("KNEE1_V2=%s", summary_dt$knee1),
    sprintf("SHIN1_V2=%s", summary_dt$shin1),
    sprintf("KNEE2_V2=%s", summary_dt$knee2),
    sprintf("SHIN2_V2=%s", summary_dt$shin2),
    sprintf("ED_FDR_V2=%s", summary_dt$ed_fdr_threshold)
  ), ambient_env_out)

  message("Written audit table: ", audit_out)
  message("Written summary:     ", summary_out)
  message("Written barcodes:    ", barcodes_out)
  message("Written env:         ", ambient_env_out)
  message("=== BARCODE_ESTIMATION_V2 done ===")
}


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

annotate_splice_flags <- function(audit_dt, splice_context) {
  audit_dt[, splice_flag := fifelse(final_call, "expected_for_context", "not_called")]

  if (splice_context == "auto") {
    audit_dt[final_call == TRUE, splice_flag := "not_evaluated"]
    return(audit_dt)
  }

  called_pct <- audit_dt[final_call == TRUE & is.finite(spliced_pct), spliced_pct]
  if (length(called_pct) < 5) {
    audit_dt[final_call == TRUE, splice_flag := "not_evaluated"]
    return(audit_dt)
  }

  q1 <- unname(quantile(called_pct, 0.25, na.rm = TRUE))
  q3 <- unname(quantile(called_pct, 0.75, na.rm = TRUE))
  iqr_val <- q3 - q1
  if (!is.finite(iqr_val) || iqr_val == 0) {
    audit_dt[final_call == TRUE, splice_flag := "not_evaluated"]
    return(audit_dt)
  }

  if (splice_context == "snrna") {
    upper_cut <- q3 + 1.5 * iqr_val
    audit_dt[final_call == TRUE, splice_flag := fifelse(
      spliced_pct > upper_cut,
      "context_outlier",
      "expected_for_context"
    )]
  } else {
    lower_cut <- q1 - 1.5 * iqr_val
    audit_dt[final_call == TRUE, splice_flag := fifelse(
      spliced_pct < lower_cut,
      "context_outlier",
      "expected_for_context"
    )]
  }

  audit_dt
}


make_summary_dt <- function(audit_dt, sample_id, total_barcodes_input, zero_total_removed,
                            ed_fdr, niters, splice_context, low_count_strategy,
                            min_umis_empty) {
  data.table(
    sample_id = sample_id,
    total_barcodes_input = total_barcodes_input,
    zero_total_removed = zero_total_removed,
    total_droplets = nrow(audit_dt),
    obvious_good = sum(audit_dt$decision_zone == "obvious_good"),
    ambiguous_tested = sum(audit_dt$tested_by_emptydrops),
    obvious_bad = sum(audit_dt$decision_zone == "obvious_bad"),
    empty_plateau_barcodes = sum(audit_dt$in_empty_plateau, na.rm = TRUE),
    retained_upper_knee = sum(audit_dt$call_reason == "retained_upper_knee"),
    called_by_emptydrops = sum(audit_dt$call_reason == "called_emptydrops"),
    excluded_low_count = sum(audit_dt$call_reason == "excluded_low_count"),
    excluded_emptydrops = sum(audit_dt$call_reason == "excluded_emptydrops"),
    final_cells = sum(audit_dt$final_call),
    splice_context_outliers = sum(audit_dt$splice_flag == "context_outlier", na.rm = TRUE),
    knee1 = as.integer(unique(audit_dt$knee1)[1]),
    shin1 = as.integer(unique(audit_dt$shin1)[1]),
    knee2 = as.integer(unique(audit_dt$knee2)[1]),
    shin2 = as.integer(unique(audit_dt$shin2)[1]),
    total_droplets_included = as.integer(unique(audit_dt$total_droplets_included)[1]),
    low_count_threshold = as.integer(unique(audit_dt$low_count_threshold)[1]),
    ed_fdr_threshold = ed_fdr,
    min_umis_empty = min_umis_empty,
    ed_niters = niters,
    low_count_strategy = low_count_strategy,
    splice_context = splice_context,
    median_spliced_pct_called = median(audit_dt[final_call == TRUE, spliced_pct], na.rm = TRUE),
    median_spliced_pct_ambiguous = median(audit_dt[tested_by_emptydrops == TRUE, spliced_pct], na.rm = TRUE),
    median_spliced_pct_obvious_bad = median(audit_dt[decision_zone == "obvious_bad", spliced_pct], na.rm = TRUE)
  )
}


calc_ambient_params_v2 <- function(split_mat, run,
                                   min_umis_empty = 5,
                                   min_umis_cells = NULL,
                                   rank_empty_plateau = NULL,
                                   low_count_threshold = "shin2",
                                   expected_cells = NA,
                                   total_included = NA,
                                   run_var = "sample_id",
                                   knee1 = NA, shin1 = NA,
                                   knee2 = NA, shin2 = NA) {
  if (is.character(low_count_threshold)) {
    if (!(low_count_threshold %in% c("knee2", "shin2"))) {
      stop('low_count_threshold must be "knee2", "shin2", or an integer')
    }
  }

  knee1_ls <- .get_knee_and_shin_1_v2(split_mat, min_umis_cells, knee1, shin1, knee2)
  knee2_ls <- .get_knee_and_shin_2_v2(
    split_mat,
    knee1_ls$ranks_dt,
    rank_empty_plateau,
    min_umis_empty,
    knee1_ls$shin1_x,
    knee2,
    shin2
  )
  params_ls <- .get_params_ls_v2(knee1_ls, knee2_ls, low_count_threshold,
    expected_cells, total_included)

  ranks_dt <- knee1_ls$ranks_dt %>%
    .[, (run_var) := run] %>%
    .[, `:=`(
      knee1 = knee1_ls$sel_knee["knee"],
      shin1 = knee1_ls$sel_knee["shin"],
      knee2 = knee2_ls$sel_knee["knee"],
      shin2 = knee2_ls$sel_knee["shin"],
      total_droplets_included = params_ls$total_included,
      low_count_threshold = params_ls$lc,
      expected_cells = params_ls$expected_cells
    )]

  ranks_dt <- .get_empty_plateau_v2(
    knee_df = ranks_dt,
    shin1 = knee1_ls$sel_knee["shin"],
    total_included = params_ls$total_included,
    knee2 = knee2_ls$sel_knee["knee"]
  )

  list(
    ranks_dt = ranks_dt,
    params_ls = params_ls,
    knee1_ls = knee1_ls,
    knee2_ls = knee2_ls
  )
}


.get_empty_plateau_v2 <- function(knee_df, shin1, total_included, knee2) {
  shin1_idx <- which.min(abs(knee_df$total - shin1))[1]
  shin1_x <- knee_df[shin1_idx, rank]

  empty_start <- copy(knee_df)[, n := .I] %>%
    .[rank %between% c(shin1_x, total_included), n] %>%
    log10() %>% mean() %>% `^`(10, .)

  empty_end <- copy(knee_df)[total == knee2, unique(rank)]

  knee_df[, in_empty_plateau := fifelse(rank %between% c(empty_start, empty_end), TRUE, FALSE)]
  knee_df
}


.get_knee_and_shin_1_v2 <- function(split_mat, min_umis_cells,
                                    knee1 = NA, shin1 = NA, knee2 = NA) {
  if (all(sapply(c(knee1, shin1, knee2), function(param) !is.na(param)))) {
    low <- median(c(knee2, shin1))
    ranks_obj <- barcodeRanks(split_mat, lower = low)
    sel_knee <- c(shin = shin1, knee = knee1)

  } else if (!is.null(min_umis_cells)) {
    ranks_obj <- barcodeRanks(split_mat, lower = min_umis_cells)
    sel_knee <- c(
      shin = as.integer(round(as.numeric(as.character(metadata(ranks_obj)$inflection)))),
      knee = as.integer(round(as.numeric(as.character(metadata(ranks_obj)$knee))))
    )

  } else {
    ranks_ls <- lapply(seq(1000, 100, by = -100),
      function(lower_cut) barcodeRanks(split_mat, lower = lower_cut))
    ks_strs <- lapply(ranks_ls, function(x) {
      paste0(
        as.integer(round(as.numeric(as.character(metadata(x)$inflection)))),
        "_",
        as.integer(round(as.numeric(as.character(metadata(x)$knee))))
      )
    }) %>% unlist()

    votes_tbl <- table(ks_strs)
    sel_cut <- names(votes_tbl)[which.max(votes_tbl)]
    sel_knee <- strsplit(sel_cut, "_")[[1]] %>% as.integer() %>% setNames(c("shin", "knee"))
    ranks_obj <- ranks_ls[[which(ks_strs == sel_cut)[1]]]
  }

  ranks_dt <- ranks_obj %>%
    as.data.frame() %>%
    as.data.table(keep.rownames = TRUE) %>%
    setnames("rn", "barcode") %>%
    .[order(rank)]

  shin1_idx <- which.min(abs(ranks_dt$total - sel_knee[1]))[1]
  shin1_x <- ranks_dt[shin1_idx, rank]

  list(ranks_dt = ranks_dt, sel_knee = sel_knee, shin1_x = shin1_x)
}


.get_knee_and_shin_2_v2 <- function(split_mat, ranks_dt, rank_empty_plateau,
                                    min_umis_empty, shin1_x, knee2, shin2) {
  if (all(sapply(c(knee2, shin2), function(param) !is.na(param)))) {
    shin2_corr <- ranks_dt[which.min(abs(ranks_dt$total - shin2))[1], total]
    knee2_corr <- ranks_dt[which.min(abs(ranks_dt$total - knee2))[1], total]
    sel_knee <- c(shin = shin2_corr, knee = knee2_corr)

  } else if (!is.null(rank_empty_plateau)) {
    ranks_smol <- ranks_dt[rank > rank_empty_plateau, barcode]
    ranks_obj <- barcodeRanks(split_mat[, ranks_smol], lower = min_umis_empty)
    sel_knee <- c(
      shin = as.integer(round(as.numeric(as.character(metadata(ranks_obj)$inflection)))),
      knee = as.integer(round(as.numeric(as.character(metadata(ranks_obj)$knee))))
    )

  } else {
    cuts <- .calc_small_knee_cuts_ls_v2(ranks_dt, min_umis_empty, shin1_x)
    ranks_ls <- lapply(cuts, function(this_cut) {
      smol <- ranks_dt[rank > this_cut, barcode]
      obj <- barcodeRanks(split_mat[, smol], lower = min_umis_empty)
      sh2 <- as.integer(round(as.numeric(as.character(metadata(obj)$inflection))))
      kn2 <- as.integer(round(as.numeric(as.character(metadata(obj)$knee))))
      paste0(sh2, "_", kn2)
    }) %>% unlist()

    shin_tbl <- str_before_first(ranks_ls, "_") %>% table()
    sel_i <- names(shin_tbl)[which.max(shin_tbl)]

    match_idx <- grepl(paste0("^", sel_i, "_"), ranks_ls)
    match_ks <- str_after_last(ranks_ls[match_idx], "_") %>% as.numeric()
    med_val <- median(match_ks)
    sel_k <- match_ks[which.min(abs(match_ks - med_val))[1]]

    sel_knee <- c(shin = as.integer(sel_i), knee = as.integer(sel_k))
  }

  knee2_x <- ranks_dt[total == sel_knee["knee"], rank] %>% unique()
  list(sel_knee = sel_knee, knee2_x = knee2_x)
}


.calc_small_knee_cuts_ls_v2 <- function(ranks_dt, min_umis_empty, shin_x) {
  last <- tail(unique(ranks_dt[total > min_umis_empty, total]), n = 3)[1]
  last_x <- ranks_dt[total == last, rank][1]

  middle <- copy(ranks_dt) %>%
    .[, n := .I] %>%
    .[rank %between% c(shin_x, last_x), n] %>%
    log10() %>% mean() %>% `^`(10, .)

  10^seq(log10(shin_x), log10(middle), length.out = 10)
}


.get_params_ls_v2 <- function(knee1_ls, knee2_ls, low_count_threshold,
                              expected_cells = NA, total_included = NA) {
  ranks_dt <- knee1_ls$ranks_dt
  shin1_x <- knee1_ls$shin1_x
  knee2_x <- knee2_ls$knee2_x

  if (is.na(expected_cells)) {
    expected_cells <- shin1_x
  }

  if (is.na(total_included)) {
    total_included <- copy(ranks_dt) %>%
      .[, n := .I] %>%
      .[rank %between% c(shin1_x, knee2_x), n] %>%
      log10() %>% mean() %>% `^`(10, .) %>% round()
  }

  if (is.character(low_count_threshold)) {
    lc <- if (low_count_threshold == "knee2") knee2_ls$sel_knee["knee"]
    else knee2_ls$sel_knee["shin"]
  } else {
    expected_x <- ranks_dt[which.min(abs(rank - expected_cells)), total]
    total_x <- ranks_dt[which.min(abs(rank - total_included)), total]
    if (!((low_count_threshold < expected_x) & (low_count_threshold < total_x))) {
      stop("low_count_threshold must be less than expected_cells and total_droplets_included")
    }
    lc <- low_count_threshold
  }

  list(expected_cells = expected_cells, total_included = total_included, lc = lc)
}


run_barcode_estimation_v2(args)