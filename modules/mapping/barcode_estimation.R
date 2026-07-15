#!/usr/bin/env Rscript
# barcode_estimation.R
#
# Loads an alevin-fry quant directory and produces:
#   - af_counts_mat.h5          (S+U+A stacked matrix, 10x v3 format)
#   - knee_plot_data_<id>.csv   (per-barcode ranks, knee params, splice stats)
#   - ambient_params_<id>.env   (CellBender/decontX priors as KEY=VALUE lines)
#
# Logic mirrors scprocess scripts/mapping.R exactly.

suppressPackageStartupMessages({
  library(magrittr)
  library(fishpond)
  library(SingleCellExperiment)
  library(DropletUtils)
  library(data.table)
  library(parallel)
  library(strex)
  library(jsonlite)
})

# ---------------------------------------------------------------------------
# Args
# ---------------------------------------------------------------------------

run_barcode_estimation <- function() {

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6) {
  stop("Usage: barcode_estimation.R <sampleId> <quantDir> <h5Out> <kneeOut> <ambientEnvOut> <statsOut>")
}

sample_id  <- args[1]
quant_dir  <- args[2]   # top-level simpleaf output dir: contains af_map/ and af_quant/
h5_out     <- args[3]
knee_out   <- args[4]
ambient_env_out <- args[5]
stats_out  <- args[6]

message("Sample:   ", sample_id)
message("Quant:    ", quant_dir)
message("H5 out:   ", h5_out)
message("Knee out: ", knee_out)
message("Params:   ", ambient_env_out)
message("Stats out:", stats_out)

# ---------------------------------------------------------------------------
# Load alevin-fry output (S, U, A modes)
# ---------------------------------------------------------------------------

sce <- loadFry(
  file.path(quant_dir, "af_quant"),
  outputFormat = list(S = c("S"), U = c("U"), A = c("A"))
)

# Stack S, U, A into a single matrix with prefixed row names
mat <- assayNames(sce) %>% lapply(function(n) {
  m            <- assay(sce, n)
  rownames(m)  <- paste0(rownames(m), "_", n)
  m
}) %>% do.call(rbind, .)

# Drop all-zero barcodes
mat <- mat[, colSums(mat) > 0]
message("Barcodes retained (total > 0): ", ncol(mat))

# Save H5
write10xCounts(h5_out, mat, version = "3", overwrite = TRUE)
message("Written H5: ", h5_out)

# Per-barcode spliced / unspliced counts (needed for knee data)
splice_dt <- data.table(
  barcode   = colnames(sce),
  spliced   = colSums(assay(sce, "S")),
  unspliced = colSums(assay(sce, "U"))
)

# ---------------------------------------------------------------------------
# Knee / ambient parameter estimation  (mirrors calc_ambient_params)
# ---------------------------------------------------------------------------

bender_ps <- calc_ambient_params(
  split_mat = mat,
  run       = sample_id,
  run_var   = "sample_id"
)

# Merge in splice stats
bender_ps <- merge(bender_ps, splice_dt, by = "barcode") %>%
  .[order(rank)]

# Save knee data CSV (uncompressed — report process can gzip if desired)
fwrite(bender_ps, file = knee_out)
message("Written knee data: ", knee_out)

# Save ambient params as shell-style KEY=VALUE lines
amb <- list(
  run                      = sample_id,
  cb_total_droplets_included = as.integer(unique(bender_ps$total_droplets_included)),
  cb_expected_cells        = as.integer(unique(bender_ps$expected_cells)),
  cb_low_count_threshold   = as.integer(unique(bender_ps$low_count_threshold)),
  knee1                    = as.integer(unique(bender_ps$knee1)),
  shin1                    = as.integer(unique(bender_ps$shin1)),
  knee2                    = as.integer(unique(bender_ps$knee2)),
  shin2                    = as.integer(unique(bender_ps$shin2))
)
writeLines(c(
  sprintf("RUN=%s", amb$run),
  sprintf("CB_TOTAL_DROPLETS_INCLUDED=%s", amb$cb_total_droplets_included),
  sprintf("CB_EXPECTED_CELLS=%s", amb$cb_expected_cells),
  sprintf("CB_LOW_COUNT_THRESHOLD=%s", amb$cb_low_count_threshold),
  sprintf("KNEE1=%s", amb$knee1),
  sprintf("SHIN1=%s", amb$shin1),
  sprintf("KNEE2=%s", amb$knee2),
  sprintf("SHIN2=%s", amb$shin2)
), ambient_env_out)
message("Written ambient params: ", ambient_env_out)

# ---------------------------------------------------------------------------
# Read-level + per-barcode stats from piscem/alevin-fry
# (mirrors dev/sn_report_review_2026-07-03/scripts/extract_alevinfry_stats.R)
# ---------------------------------------------------------------------------

mi <- fromJSON(file.path(quant_dir, "af_map", "map_info.json"))
fd <- fread(file.path(quant_dir, "af_quant", "featureDump.txt"))

total_mapped <- sum(fd$MappedReads)
total_dedup  <- sum(fd$DeduplicatedReads)

# UMI floor (200) matches the hard "definite empty" floor used elsewhere in
# the pipeline -- keeps the above-floor summary from being swamped by the
# millions of background barcodes with 1-2 genes.
fd_above <- fd[DeduplicatedReads >= 200]

stats_dt <- data.table(
  sample_id              = sample_id,
  num_reads              = as.integer(mi$num_reads),
  num_mapped             = as.integer(mi$num_mapped),
  percent_mapped         = as.numeric(mi$percent_mapped),
  runtime_seconds        = as.numeric(mi$runtime_seconds),
  n_barcodes_detected    = nrow(fd),
  n_barcodes_above_floor = nrow(fd_above),
  saturation             = 1 - total_dedup / total_mapped,
  median_genes_per_bc    = median(fd_above$NumGenesExpressed, na.rm = TRUE),
  median_mean_by_max     = median(fd_above$MeanByMax, na.rm = TRUE)
)
fwrite(stats_dt, file = stats_out)
message("Written alevin-fry stats: ", stats_out)
}


# ===========================================================================
# Functions (mirrors scprocess mapping.R)
# ===========================================================================

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

  knee1_ls <- .get_knee_and_shin_1(split_mat, min_umis_cells, knee1, shin1, knee2)
  knee2_ls <- .get_knee_and_shin_2(split_mat, knee1_ls$ranks_dt,
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
  shin1_idx  <- which.min(abs(knee_df$total - shin1))[1]
  shin1_x    <- knee_df[shin1_idx, rank]

  # Plateau starts immediately where the cell cutoff ends (shin1_x), not at
  # the geometric mean of [shin1_x, total_included]. The geometric mean left
  # an unclaimed rank band (neither cell nor ambient) between the cell
  # cutoff and the ambient plateau on every sample -- see
  # dev/sn_report_review_2026-07-03/ for the investigation.
  empty_start <- shin1_x

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
      smol      <- ranks_dt[rank > this_cut, barcode]
      obj       <- barcodeRanks(split_mat[, smol], lower = min_umis_empty)
      sh2       <- as.integer(round(as.numeric(as.character(metadata(obj)$inflection))))
      kn2       <- as.integer(round(as.numeric(as.character(metadata(obj)$knee))))
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


run_barcode_estimation()
