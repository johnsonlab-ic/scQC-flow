#!/usr/bin/env Rscript
# barcode_estimation_v2.R
#
# Standalone barcode caller for alevin-fry outputs using EmptyDrops.
#
# Knee/ambient parameters (knee1, shin1, knee2, shin2) are estimated with the
# same voting approach used in scprocess (barcode_estimation.R).  EmptyDrops
# is then called with:
#   lower  = shin2  → barcodes at or below shin2 define the ambient profile
#                     and are NOT tested (avoiding FDR inflation across ~1.9M tests)
#   retain = knee1  → barcodes above knee1 are auto-retained without testing
#
# Only the barcodes strictly between shin2 and knee1 (~100–200 k) are actually
# tested and FDR-corrected.
#
# Outputs
#   af_counts_mat_v2.h5              raw stacked S+U+A matrix (10x v3 format)
#   barcode_audit_v2_<id>.csv.gz     per-barcode audit table
#   barcode_summary_v2_<id>.csv      per-sample decision summary
#   cell_barcodes_v2_<id>.csv        final called barcodes
#   ambient_params_v2_<id>.env       thresholds for downstream (decontX, QC)

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
parser$add_argument("--sample_id",          type = "character", required = TRUE)
parser$add_argument("--quant_dir",          type = "character", required = TRUE)
parser$add_argument("--h5_out",             type = "character", default = "af_counts_mat_v2.h5")
parser$add_argument("--audit_out",          type = "character", default = "barcode_audit_v2.csv.gz")
parser$add_argument("--summary_out",        type = "character", default = "barcode_summary_v2.csv")
parser$add_argument("--barcodes_out",       type = "character", default = "cell_barcodes_v2.csv")
parser$add_argument("--ambient_env_out",    type = "character", default = "ambient_params_v2.env")
parser$add_argument("--min_umis_empty",     type = "integer",   default = 5)
parser$add_argument("--niters",             type = "integer",   default = 10000)
parser$add_argument("--ed_fdr",             type = "double",    default = 0.001)
parser$add_argument("--splice_context",     type = "character", default = "snrna")
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

  # -------------------------------------------------------------------------
  # Load alevin-fry output
  # -------------------------------------------------------------------------

  sce <- loadFry(
    quant_dir,
    outputFormat = list(S = c("S"), U = c("U"), A = c("A"))
  )

  total_counts     <- assay(sce, "S") + assay(sce, "U") + assay(sce, "A")
  n_barcodes_input <- ncol(total_counts)

  # Stacked S+U+A matrix (for H5 output used by downstream QC)
  mat <- assayNames(sce) %>% lapply(function(mode_name) {
    m             <- assay(sce, mode_name)
    rownames(m)   <- paste0(rownames(m), "_", mode_name)
    m
  }) %>% do.call(rbind, .)

  keep_bcs       <- Matrix::colSums(mat) > 0
  mat            <- mat[, keep_bcs, drop = FALSE]
  total_counts   <- total_counts[, keep_bcs, drop = FALSE]
  sce            <- sce[, keep_bcs]
  n_zero_removed <- n_barcodes_input - ncol(total_counts)

  message("Input barcodes:      ", n_barcodes_input)
  message("Zero-count removed:  ", n_zero_removed)
  message("Non-zero barcodes:   ", ncol(total_counts))

  write10xCounts(h5_out, mat, version = "3", overwrite = TRUE)
  message("Written raw H5: ", h5_out)

  # -------------------------------------------------------------------------
  # Per-barcode splice stats
  # -------------------------------------------------------------------------

  barcode_extra_dt <- data.table(
    barcode          = colnames(total_counts),
    detected_genes   = Matrix::colSums(total_counts > 0),
    spliced          = Matrix::colSums(assay(sce, "S")),
    unspliced        = Matrix::colSums(assay(sce, "U")),
    ambiguous_counts = Matrix::colSums(assay(sce, "A"))
  )
  barcode_extra_dt[, spliced_pct := (spliced + 1) / (spliced + unspliced + 2)]
  barcode_extra_dt[, logit_spliced_pct := qlogis(spliced_pct)]

  # -------------------------------------------------------------------------
  # Knee / ambient parameter estimation (mirrors scprocess calc_ambient_params)
  # -------------------------------------------------------------------------

  thresh_dt <- calc_ambient_params(
    split_mat           = mat,
    run                 = sample_id,
    min_umis_empty      = min_umis_empty,
    low_count_threshold = low_count_strategy,
    run_var             = "sample_id"
  )

  knee1_val  <- as.integer(unique(thresh_dt$knee1)[1])
  shin1_val  <- as.integer(unique(thresh_dt$shin1)[1])
  knee2_val  <- as.integer(unique(thresh_dt$knee2)[1])
  shin2_val  <- as.integer(unique(thresh_dt$shin2)[1])
  lc_val     <- as.integer(unique(thresh_dt$low_count_threshold)[1])
  retain_val <- knee1_val

  message("knee1: ", knee1_val, "  shin1: ", shin1_val)
  message("knee2: ", knee2_val, "  shin2: ", shin2_val)
  message("EmptyDrops lower (ambient floor): ", lc_val)
  message("EmptyDrops retain (auto-call floor): ", retain_val)

  # Safety check: lower must be strictly below retain
  if (lc_val >= retain_val) {
    warning(sprintf(
      "low_count_threshold (%d) >= knee1 (%d) — forcing retain to knee1 + 1.",
      lc_val, retain_val
    ))
    retain_val <- lc_val + 1L
  }

  # -------------------------------------------------------------------------
  # Decision zones (for audit / reporting)
  # -------------------------------------------------------------------------

  audit_dt <- merge(thresh_dt, barcode_extra_dt, by = "barcode", all.x = TRUE)
  setorder(audit_dt, rank)

  audit_dt[, decision_zone := fifelse(
    total > retain_val,  "obvious_good",
    fifelse(total <= lc_val, "obvious_bad", "ambiguous")
  )]
  audit_dt[, tested_by_emptydrops := decision_zone == "ambiguous"]

  message("Decision zones:")
  message("  obvious_good (total > ", retain_val, "): ", sum(audit_dt$decision_zone == "obvious_good"))
  message("  ambiguous    (", lc_val, " < total <= ", retain_val, "): ", sum(audit_dt$decision_zone == "ambiguous"))
  message("  obvious_bad  (total <= ", lc_val, "): ", sum(audit_dt$decision_zone == "obvious_bad"))

  # -------------------------------------------------------------------------
  # EmptyDrops
  #
  # Key design: use lower = lc_val (NOT known.empty).
  # When known.empty is provided, DropletUtils ignores lower for ambient
  # estimation and tests ALL barcodes — FDR over 1.9M tests kills every call.
  # With lower = lc_val, the ~1.85M barcodes at or below lc_val define the
  # ambient profile and are excluded from testing entirely.  Only the
  # ~100-200k barcodes in the ambiguous zone are tested and FDR-corrected.
  # -------------------------------------------------------------------------

  n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1"))
  bpparam <- MulticoreParam(workers = max(1L, n_cores), progressbar = FALSE)

  message("Running EmptyDrops: lower=", lc_val, ", retain=", retain_val, ", niters=", niters)

  ed_res <- emptyDrops(
    m       = total_counts,
    lower   = lc_val,
    retain  = retain_val,
    niters  = niters,
    BPPARAM = bpparam
  )

  ed_dt <- as.data.table(as.data.frame(ed_res), keep.rownames = TRUE)
  setnames(ed_dt, "rn", "barcode")
  setnames(
    ed_dt,
    c("Total", "LogProb", "PValue", "Limited", "FDR"),
    c("ED_Total", "ED_LogProb", "ED_PValue", "ED_Limited", "ED_FDR")
  )

  n_tested  <- sum(!is.na(ed_dt$ED_PValue))
  n_limited <- sum(ed_dt$ED_Limited, na.rm = TRUE)
  message("EmptyDrops: ", n_tested, " barcodes tested, ", n_limited,
          " Limited (", round(100 * n_limited / max(n_tested, 1L), 1), "%)",
          if (n_limited / max(n_tested, 1L) > 0.5) " — consider increasing --niters" else "")

  # -------------------------------------------------------------------------
  # Final calls
  # -------------------------------------------------------------------------

  audit_dt <- merge(
    audit_dt,
    ed_dt[, .(barcode, ED_Total, ED_LogProb, ED_PValue, ED_Limited, ED_FDR)],
    by     = "barcode",
    all.x  = TRUE
  )
  setorder(audit_dt, rank)

  # A barcode is a cell if it is above retain (auto-retained, FDR=NA from ED)
  # OR it is in the ambiguous zone and passes the FDR threshold.
  audit_dt[, final_call := decision_zone == "obvious_good" |
    (tested_by_emptydrops & !is.na(ED_FDR) & ED_FDR <= ed_fdr)]

  audit_dt[, call_reason := fcase(
    decision_zone == "obvious_good",           "retained_upper_knee",
    decision_zone == "obvious_bad",            "excluded_low_count",
    tested_by_emptydrops & final_call,         "called_emptydrops",
    tested_by_emptydrops & !final_call,        "excluded_emptydrops",
    default = "unclassified"
  )]

  audit_dt <- annotate_splice_flags(audit_dt, splice_context)

  message("Final calls:")
  message("  retained_upper_knee: ", sum(audit_dt$call_reason == "retained_upper_knee"))
  message("  called_emptydrops:   ", sum(audit_dt$call_reason == "called_emptydrops"))
  message("  total cells:         ", sum(audit_dt$final_call))

  # -------------------------------------------------------------------------
  # Summary table
  # -------------------------------------------------------------------------

  summary_dt <- data.table(
    sample_id               = sample_id,
    total_barcodes_input    = n_barcodes_input,
    zero_total_removed      = n_zero_removed,
    total_droplets          = nrow(audit_dt),
    obvious_good            = sum(audit_dt$decision_zone == "obvious_good"),
    ambiguous_tested        = sum(audit_dt$tested_by_emptydrops),
    obvious_bad             = sum(audit_dt$decision_zone == "obvious_bad"),
    empty_plateau_barcodes  = sum(audit_dt$in_empty_plateau, na.rm = TRUE),
    retained_upper_knee     = sum(audit_dt$call_reason == "retained_upper_knee"),
    called_by_emptydrops    = sum(audit_dt$call_reason == "called_emptydrops"),
    excluded_low_count      = sum(audit_dt$call_reason == "excluded_low_count"),
    excluded_emptydrops     = sum(audit_dt$call_reason == "excluded_emptydrops"),
    final_cells             = sum(audit_dt$final_call),
    splice_context_outliers = sum(audit_dt$splice_flag == "context_outlier", na.rm = TRUE),
    knee1                   = knee1_val,
    shin1                   = shin1_val,
    knee2                   = knee2_val,
    shin2                   = shin2_val,
    total_droplets_included = as.integer(unique(thresh_dt$total_droplets_included)[1]),
    low_count_threshold     = lc_val,
    ed_retain_threshold     = retain_val,
    ed_fdr_threshold        = ed_fdr,
    min_umis_empty          = min_umis_empty,
    ed_niters               = niters,
    low_count_strategy      = low_count_strategy,
    splice_context          = splice_context,
    median_spliced_pct_called   = median(audit_dt[final_call == TRUE,         spliced_pct], na.rm = TRUE),
    median_spliced_pct_ambient  = median(audit_dt[decision_zone == "obvious_bad", spliced_pct], na.rm = TRUE)
  )

  # -------------------------------------------------------------------------
  # Write outputs
  # -------------------------------------------------------------------------

  fwrite(audit_dt, file = audit_out, compress = "gzip")
  fwrite(summary_dt, file = summary_out)
  fwrite(audit_dt[final_call == TRUE, .(barcode)], file = barcodes_out)

  writeLines(c(
    sprintf("RUN=%s",                       sample_id),
    sprintf("CB_EXPECTED_CELLS_V2=%s",      summary_dt$final_cells),
    sprintf("CB_LOW_COUNT_THRESHOLD_V2=%s", lc_val),
    sprintf("KNEE1_V2=%s",                  knee1_val),
    sprintf("SHIN1_V2=%s",                  shin1_val),
    sprintf("KNEE2_V2=%s",                  knee2_val),
    sprintf("SHIN2_V2=%s",                  shin2_val),
    sprintf("ED_FDR_V2=%s",                 ed_fdr)
  ), ambient_env_out)

  message("Written audit table: ", audit_out)
  message("Written summary:     ", summary_out)
  message("Written barcodes:    ", barcodes_out)
  message("Written env:         ", ambient_env_out)
  message("=== BARCODE_ESTIMATION_V2 done ===")
}


# ---------------------------------------------------------------------------
# Helper: splice context flags
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

  q1      <- unname(quantile(called_pct, 0.25, na.rm = TRUE))
  q3      <- unname(quantile(called_pct, 0.75, na.rm = TRUE))
  iqr_val <- q3 - q1

  if (!is.finite(iqr_val) || iqr_val == 0) {
    audit_dt[final_call == TRUE, splice_flag := "not_evaluated"]
    return(audit_dt)
  }

  if (splice_context == "snrna") {
    upper_cut <- q3 + 1.5 * iqr_val
    audit_dt[final_call == TRUE, splice_flag := fifelse(
      spliced_pct > upper_cut, "context_outlier", "expected_for_context"
    )]
  } else {
    lower_cut <- q1 - 1.5 * iqr_val
    audit_dt[final_call == TRUE, splice_flag := fifelse(
      spliced_pct < lower_cut, "context_outlier", "expected_for_context"
    )]
  }

  audit_dt
}


# ---------------------------------------------------------------------------
# Knee / ambient parameter functions (verbatim from scprocess barcode_estimation.R)
# ---------------------------------------------------------------------------

calc_ambient_params <- function(split_mat, run,
                                min_umis_empty      = 5,
                                min_umis_cells      = NULL,
                                rank_empty_plateau  = NULL,
                                low_count_threshold = "shin2",
                                expected_cells      = NA,
                                total_included      = NA,
                                run_var             = "sample_id",
                                knee1 = NA, shin1 = NA,
                                knee2 = NA, shin2 = NA) {
  if (is.character(low_count_threshold)) {
    if (!(low_count_threshold %in% c("knee2", "shin2"))) {
      stop('low_count_threshold must be "knee2", "shin2", or an integer')
    }
  }

  knee1_ls  <- .get_knee_and_shin_1(split_mat, min_umis_cells, knee1, shin1, knee2)
  knee2_ls  <- .get_knee_and_shin_2(split_mat, knee1_ls$ranks_dt,
                                    rank_empty_plateau, min_umis_empty,
                                    knee1_ls$shin1_x, knee2, shin2)
  params_ls <- .get_params_ls(knee1_ls, knee2_ls, low_count_threshold,
                               expected_cells, total_included)

  bender_ps <- knee1_ls$ranks_dt %>%
    .[, (run_var) := run] %>%
    .[, `:=`(
      knee1                   = knee1_ls$sel_knee["knee"],
      shin1                   = knee1_ls$sel_knee["shin"],
      knee2                   = knee2_ls$sel_knee["knee"],
      shin2                   = knee2_ls$sel_knee["shin"],
      total_droplets_included = params_ls$total_included,
      low_count_threshold     = params_ls$lc,
      expected_cells          = params_ls$expected_cells
    )]

  bender_ps <- .get_empty_plateau(
    knee_df        = bender_ps,
    shin1          = knee1_ls$sel_knee["shin"],
    total_included = params_ls$total_included,
    knee2          = knee2_ls$sel_knee["knee"]
  )

  return(bender_ps)
}


.get_empty_plateau <- function(knee_df, shin1, total_included, knee2) {
  shin1_idx   <- which.min(abs(knee_df$total - shin1))[1]
  shin1_x     <- knee_df[shin1_idx, rank]

  empty_start <- copy(knee_df)[, n := .I] %>%
    .[rank %between% c(shin1_x, total_included), n] %>%
    log10() %>% mean() %>% (function(x) 10^x)()

  empty_end <- copy(knee_df)[total == knee2, unique(rank)]

  knee_df[, in_empty_plateau :=
    fifelse(rank %between% c(empty_start, empty_end), TRUE, FALSE)]
  return(knee_df)
}


.get_knee_and_shin_1 <- function(split_mat, min_umis_cells,
                                  knee1 = NA, shin1 = NA, knee2 = NA) {
  if (all(sapply(c(knee1, shin1, knee2), function(p) !is.na(p)))) {
    low       <- median(c(knee2, shin1))
    ranks_obj <- barcodeRanks(split_mat, lower = low)
    sel_knee  <- c(shin = shin1, knee = knee1)

  } else if (!is.null(min_umis_cells)) {
    ranks_obj <- barcodeRanks(split_mat, lower = min_umis_cells)
    sel_knee  <- c(
      shin = as.integer(round(as.numeric(as.character(metadata(ranks_obj)$inflection)))),
      knee = as.integer(round(as.numeric(as.character(metadata(ranks_obj)$knee))))
    )

  } else {
    ranks_ls <- lapply(seq(1000, 100, by = -100),
                       function(x) barcodeRanks(split_mat, lower = x))
    ks_strs  <- lapply(ranks_ls, function(x)
      paste0(
        as.integer(round(as.numeric(as.character(metadata(x)$inflection)))),
        "_",
        as.integer(round(as.numeric(as.character(metadata(x)$knee))))
      )) %>% unlist()

    votes_tbl <- table(ks_strs)
    sel_cut   <- names(votes_tbl)[which.max(votes_tbl)]
    sel_knee  <- strsplit(sel_cut, "_")[[1]] %>% as.integer() %>%
      setNames(c("shin", "knee"))
    ranks_obj <- ranks_ls[[which(ks_strs == sel_cut)[1]]]
  }

  ranks_dt  <- ranks_obj %>% as.data.frame() %>%
    as.data.table(keep.rownames = TRUE) %>%
    setnames("rn", "barcode") %>%
    .[order(rank)]

  shin1_idx <- which.min(abs(ranks_dt$total - sel_knee[1]))[1]
  shin1_x   <- ranks_dt[shin1_idx, rank]

  list(ranks_dt = ranks_dt, sel_knee = sel_knee, shin1_x = shin1_x)
}


.get_knee_and_shin_2 <- function(split_mat, ranks_dt, rank_empty_plateau,
                                  min_umis_empty, shin1_x, knee2, shin2) {
  if (all(sapply(c(knee2, shin2), function(p) !is.na(p)))) {
    shin2_corr <- ranks_dt[which.min(abs(ranks_dt$total - shin2))[1], total]
    knee2_corr <- ranks_dt[which.min(abs(ranks_dt$total - knee2))[1], total]
    sel_knee   <- c(shin = shin2_corr, knee = knee2_corr)

  } else if (!is.null(rank_empty_plateau)) {
    ranks_smol <- ranks_dt[rank > rank_empty_plateau, barcode]
    ranks_obj  <- barcodeRanks(split_mat[, ranks_smol], lower = min_umis_empty)
    sel_knee   <- c(
      shin = as.integer(round(as.numeric(as.character(metadata(ranks_obj)$inflection)))),
      knee = as.integer(round(as.numeric(as.character(metadata(ranks_obj)$knee))))
    )

  } else {
    cuts     <- .calc_small_knee_cuts_ls(ranks_dt, min_umis_empty, shin1_x)
    ranks_ls <- lapply(cuts, function(this_cut) {
      smol <- ranks_dt[rank > this_cut, barcode]
      obj  <- barcodeRanks(split_mat[, smol], lower = min_umis_empty)
      sh2  <- as.integer(round(as.numeric(as.character(metadata(obj)$inflection))))
      kn2  <- as.integer(round(as.numeric(as.character(metadata(obj)$knee))))
      paste0(sh2, "_", kn2)
    }) %>% unlist()

    shin_tbl  <- str_before_first(ranks_ls, "_") %>% table()
    sel_i     <- names(shin_tbl)[which.max(shin_tbl)]

    match_idx <- grepl(paste0("^", sel_i, "_"), ranks_ls)
    match_ks  <- str_after_last(ranks_ls[match_idx], "_") %>% as.numeric()
    med_val   <- median(match_ks)
    sel_k     <- match_ks[which.min(abs(match_ks - med_val))[1]]

    sel_knee  <- c(shin = as.integer(sel_i), knee = as.integer(sel_k))
  }

  knee2_x <- ranks_dt[total == sel_knee["knee"], rank] %>% unique()
  list(sel_knee = sel_knee, knee2_x = knee2_x)
}


.calc_small_knee_cuts_ls <- function(ranks_dt, min_umis_empty, shin_x) {
  last   <- tail(unique(ranks_dt[total > min_umis_empty, total]), n = 3)[1]
  last_x <- ranks_dt[total == last, rank][1]

  middle <- copy(ranks_dt) %>%
    .[, n := .I] %>%
    .[rank %between% c(shin_x, last_x), n] %>%
    log10() %>% mean() %>% `^`(10, .)

  10^seq(log10(shin_x), log10(middle), length.out = 10)
}


.get_params_ls <- function(knee1_ls, knee2_ls, low_count_threshold,
                            expected_cells = NA, total_included = NA) {
  ranks_dt <- knee1_ls$ranks_dt
  shin1_x  <- knee1_ls$shin1_x
  knee2_x  <- knee2_ls$knee2_x

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
          else                                 knee2_ls$sel_knee["shin"]
  } else {
    expected_x <- ranks_dt[which.min(abs(rank - expected_cells)), total]
    total_x    <- ranks_dt[which.min(abs(rank - total_included)), total]
    if (!((low_count_threshold < expected_x) & (low_count_threshold < total_x))) {
      stop("low_count_threshold must be less than expected_cells and total_droplets_included")
    }
    lc <- low_count_threshold
  }

  list(expected_cells = expected_cells, total_included = total_included, lc = lc)
}


run_barcode_estimation_v2(args)
