#!/usr/bin/env Rscript
# cell_calling.R
#
# Splice-aware per-sample cell selection. Replaces decontX's cell-calling role.
# Reads the stacked S+U+A H5 + knee CSV from BARCODE_ESTIMATION and classifies
# barcodes with an anchored per-sample Gaussian mixture in
#   x = (log10(UMI + ps), logit(splice%))
# then enforces two HARD biological guardrails as overrides on the GMM call.
#
# Classifier: anchored 3-component GMM (ambient / damaged / nuclei), seeded from
# fixed threshold quadrants for component IDENTITY only (not the final call),
# fit (mclust::me, model VVV = free covariance) on a stratified subsample of the
# cell-candidate space (total >= EMPTY_FLOOR), scored analytically for ALL
# candidates. See gmm_population_classify.R for the diagnostic version.
#
# Hard guardrails (override the GMM):
#   1. total UMI < EMPTY_FLOOR (200)        -> empty   (never a nucleus)
#   2. splice% > damaged_cut (adaptive)     -> damaged (high splice at ANY UMI)
#      GMM-damaged with splice <= cut       -> reassigned to argmax{ambient,nuclei}
#
# Populations: empty / ambient / damaged / cell(=nuclei) / ambiguous(low posterior)
#   is_empty = empty U ambient (U damaged if --include_damaged)
#   is_cell  = cell U ambiguous            (FRINGE -> is_cell; see TODO in MEMORY.md:
#              compare is_cell vs is_empty vs separate-ambiguous routing later)
#
# Outputs (Nextflow cwd):
#   filt_counts_<id>.h5             stacked S+U+A matrix subset to is_cell (UNCORRECTED)
#   cell_barcodes_<id>.csv          is_cell barcodes (no header) — decontX contract
#   empty_barcodes_<id>.csv         is_empty barcodes (no header) — ambient profile input
#   cell_calling_labels_<id>.csv.gz per-barcode posteriors + population + flags
#   cell_calling_summary_<id>.csv   per-sample counts + fitted cuts
#   cell_calling_gmm_<id>.rds       fitted GMM params + cuts (for the report)
#   cell_calling_<id>.env           KEY=VALUE for Nextflow env() outputs

suppressPackageStartupMessages({
  library(rhdf5)
  library(DropletUtils)
  library(data.table)
  library(Matrix)
  library(mclust)
  library(mvtnorm)
})

# Hard guardrails / defaults ---------------------------------------------------
EMPTY_FLOOR   <- 200      # total UMI below this is ALWAYS empty (guardrail 1)
PSCOUNT       <- 10
DAMAGED_NMADS <- 3.0      # damaged splice cut = median + NMADS*MAD of nuclei logit-splice
DAMAGED_FLOOR <- 0.35     # damaged splice% guardrail floor (never below)
DAMAGED_CAP   <- 0.70
TAU           <- 0.80     # max posterior below this => ambiguous (fringe)
SUB_CAP       <- 60000L   # barcodes used to FIT the GMM (scoring uses all candidates)
MIN_SEED      <- 50L      # drop a GMM component if its seed quadrant is smaller
GMM_PRIOR_DOF <- 2000     # strength of the covariance-regularising prior (see fit below)

# Seed quadrants — EM INITIALISATION ONLY (anchors component identity) ----------
# TWO-component scheme: the GMM separates only ambient vs nuclei (2 ellipses).
# 'damaged' and 'empty' are handled by hard rules, not GMM components.
SEED_AMBIENT_HI <- 1000
SEED_INTACT_UMI <- 2000
SEED_INTACT_SP  <- 0.30
POPS <- c("ambient", "nuclei")

.get_h5_mx <- function(af_mat_f) {
  h5  <- H5Fopen(af_mat_f, flags = "H5F_ACC_RDONLY")
  mat <- sparseMatrix(
    i = as.vector(h5$matrix$indices + 1L), p = as.vector(h5$matrix$indptr),
    x = as.vector(h5$matrix$data), repr = "C", dims = h5$matrix$shape)
  colnames(mat) <- as.character(h5$matrix$barcodes)
  rownames(mat) <- as.character(h5$matrix$features$name)
  H5Fclose(h5)
  mat
}

run_cell_calling <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 3) {
    stop("Usage: cell_calling.R <sampleId> <h5_file> <knee_csv> [include_damaged=FALSE]")
  }
  sample_id       <- args[1]
  h5_file         <- args[2]
  knee_csv        <- args[3]
  include_damaged <- length(args) >= 4 && toupper(args[4]) %in% c("TRUE", "T", "1", "YES")
  alpha_csv       <- if (length(args) >= 5) args[5] else ""   # optional prior CellSweep alpha_hat
  alpha_max       <- if (length(args) >= 6) as.numeric(args[6]) else 0.5
  message("Sample: ", sample_id, " | include_damaged in is_empty: ", include_damaged)

  dt <- fread(knee_csv)
  dt <- dt[total > 0]
  dt[, splice_pct := spliced / total]
  dt[, x := log10(total + PSCOUNT)]
  dt[, y := qlogis((spliced + 1) / (total + 2))]

  # Guardrail 1: hard empty floor — below this never sees the GMM ---------------
  dt[, population := NA_character_]
  dt[total < EMPTY_FLOOR, population := "empty"]
  cand <- dt[total >= EMPTY_FLOOR]
  if (nrow(cand) < 100) stop("Too few candidate barcodes above EMPTY_FLOOR")

  # Seed quadrants (init only) -------------------------------------------------
  cand[, seed := fcase(
    total > SEED_INTACT_UMI & splice_pct < SEED_INTACT_SP,         "nuclei",
    total <= SEED_AMBIENT_HI,                                      "ambient",
    default =                                                      "other"
  )]
  seed_means <- cand[seed %in% POPS, .(mx = mean(x), my = mean(y), n = .N), by = seed]
  use_pops   <- POPS[POPS %in% seed_means[n >= MIN_SEED, seed]]
  G          <- length(use_pops)
  if (G < 2) stop("Fewer than 2 usable GMM components")
  sm <- as.matrix(seed_means[match(use_pops, seed), .(mx, my)])

  # Fit anchored GMM on a stratified subsample ---------------------------------
  set.seed(42)
  fit_idx <- cand[, .(i = .I[sample(.N, min(.N, SUB_CAP))]), by = seed]$i
  if (length(fit_idx) > SUB_CAP) fit_idx <- sample(fit_idx, SUB_CAP)
  Xfit   <- as.matrix(cand[fit_idx, .(x, y)])
  d2seed <- sapply(seq_len(G), function(j) (Xfit[, 1] - sm[j, 1])^2 + (Xfit[, 2] - sm[j, 2])^2)
  z_init <- matrix(0, nrow(Xfit), G)
  z_init[cbind(seq_len(nrow(Xfit)), max.col(-d2seed))] <- 1
  # Conjugate (inverse-Wishart) prior on the component covariances. An unconstrained
  # VVV fit lets the nuclei ellipse balloon on quality-heterogeneous samples (splice-
  # axis sd ~0.7 vs ~0.3 on clean ones), which (a) spuriously collapses the ambient/
  # nuclei separability metric even when the component means are well separated and
  # (b) over-calls nuclei. The prior shrinks covariances toward
  # diag(0.3^2 log10-UMI, 0.4^2 logit-splice); dof sets its strength (weak vs the ~60k
  # fit points, but enough to rein in the broad nuclei cloud). Clean samples are
  # essentially unchanged; degraded samples get tighter, more conservative nuclei calls.
  gmm_prior <- priorControl(scale = diag(c(0.09, 0.16)) * GMM_PRIOR_DOF,
                            dof = GMM_PRIOR_DOF, shrinkage = 0)
  fit <- me(modelName = "VVV", data = Xfit, z = z_init, prior = gmm_prior)
  if (is.null(fit$parameters$mean)) stop("GMM EM failed")
  par <- fit$parameters
  for (j in seq_len(G))
    message(sprintf("  %-8s x=%.2f y=%.2f pi=%.3f",
                    use_pops[j], par$mean[1, j], par$mean[2, j], par$pro[j]))

  # --- Separability of ambient vs nuclei (sample-quality signal) --------------
  # Bhattacharyya distance between the two fitted 2-D Gaussians; higher = better
  # separated. Collapses on degraded/soup-dominated samples (see dev/cellcalling).
  .bhattacharyya <- function(mu1, S1, mu2, S2) {
    dmu <- mu1 - mu2
    S   <- (S1 + S2) / 2
    0.125 * as.numeric(t(dmu) %*% solve(S) %*% dmu) +
      0.5 * log(det(S) / sqrt(det(S1) * det(S2)))
  }
  sep_bhatt <- NA_real_; sep_splice_gap <- NA_real_; sep_umi_gap <- NA_real_
  ai <- match("ambient", use_pops); ni <- match("nuclei", use_pops)
  if (!is.na(ai) && !is.na(ni)) {
    sep_bhatt      <- .bhattacharyya(par$mean[, ai], par$variance$sigma[, , ai],
                                     par$mean[, ni], par$variance$sigma[, , ni])
    sep_splice_gap <- abs(par$mean[2, ai] - par$mean[2, ni])   # logit-splice axis
    sep_umi_gap    <- abs(par$mean[1, ai] - par$mean[1, ni])   # log10-UMI axis
    message(sprintf("  separability: bhattacharyya=%.3f splice_gap=%.3f umi_gap=%.3f",
                    sep_bhatt, sep_splice_gap, sep_umi_gap))
  }

  # Score ALL candidates analytically -----------------------------------------
  Xall <- as.matrix(cand[, .(x, y)])
  dens <- sapply(seq_len(G), function(j)
    par$pro[j] * dmvnorm(Xall, mean = par$mean[, j], sigma = par$variance$sigma[, , j]))
  if (is.null(dim(dens))) dens <- matrix(dens, ncol = G)
  rs   <- rowSums(dens); bad <- rs <= 0 | !is.finite(rs)
  post <- dens / rs; post[bad, ] <- 1 / G
  cls  <- max.col(post)
  conf <- post[cbind(seq_len(nrow(post)), cls)]; conf[bad] <- 0
  ent  <- -rowSums(ifelse(post > 0, post * log(post), 0))
  gmm_class <- use_pops[cls]
  for (j in seq_len(G)) cand[, paste0("p_", use_pops[j]) := post[, j]]
  cand[, c("gmm_class", "confidence", "entropy") := .(gmm_class, conf, ent)]

  # --- damaged splice cut: median + NMADS*MAD of nuclei-component logit-splice ----
  nuc_pool <- cand[gmm_class == "nuclei" & splice_pct < SEED_INTACT_SP]
  if (nrow(nuc_pool) < 50) nuc_pool <- cand[gmm_class == "nuclei"]
  damaged_cut <- if (nrow(nuc_pool) >= 50)
    plogis(median(nuc_pool$y) + DAMAGED_NMADS * mad(nuc_pool$y)) else 0.50
  damaged_cut <- min(max(damaged_cut, DAMAGED_FLOOR), DAMAGED_CAP)
  # cell floor = right edge of the ambient ellipse (98th pctile of GMM-ambient UMI).
  # The ONE hard UMI rule: anything to the LEFT of ambient (lower UMI) is ambient/empty.
  amb_umi    <- cand[gmm_class == "ambient", total]
  cell_floor <- if (length(amb_umi) >= 50) as.numeric(quantile(amb_umi, 0.98)) else EMPTY_FLOOR
  cell_floor <- max(cell_floor, EMPTY_FLOOR)
  message(sprintf("  cell floor = %.0f UMI | damaged splice cut = %.3f | GMM: %s",
                  cell_floor, damaged_cut, paste(use_pops, collapse="/")))

  # --- Final classification ----------------------------------------------------
  #   empty   : UMI < 200            (set above, never sees the GMM)
  #   ambient : UMI < cell_floor     (left of ambient's right edge -> ambient/empty)
  #   damaged : splice% > damaged_cut
  #   else    : GMM 2-component argmax -> ambient | nuclei   (nuclei = the cells)
  cand[, population := fcase(
    total < cell_floor,        "ambient",
    splice_pct > damaged_cut,  "damaged",
    default =                  gmm_class
  )]

  dt <- rbind(dt[total < EMPTY_FLOOR], cand, fill = TRUE)
  dt[, is_empty := population %in% c("empty", "ambient") |
                   (include_damaged & population == "damaged")]
  dt[, is_cell  := population == "nuclei"]

  # Optional post-hoc CellSweep alpha filter (TEST): drop is_cell with alpha_hat >= alpha_max.
  if (nzchar(alpha_csv) && file.exists(alpha_csv)) {
    af <- fread(alpha_csv)
    hi <- af[alpha_hat >= alpha_max, barcode]
    n0 <- sum(dt$is_cell)
    dt[barcode %in% hi, is_cell := FALSE]
    message(sprintf("  alpha filter (keep <%.2f): is_cell %d -> %d (dropped %d high-alpha)",
                    alpha_max, n0, sum(dt$is_cell), n0 - sum(dt$is_cell)))
  }

  tab <- dt[, .(n = .N, pct = round(100 * .N / nrow(dt), 1)), by = population][
    order(match(population, c("empty", "ambient", "damaged", "nuclei")))]
  print(tab)
  message("is_cell: ", sum(dt$is_cell), " | is_empty: ", sum(dt$is_empty))

  # Subset stacked matrix to is_cell (UNCORRECTED) -----------------------------
  af_mat   <- .get_h5_mx(h5_file)
  cell_bcs <- dt[is_cell == TRUE, barcode]
  cell_bcs <- cell_bcs[cell_bcs %in% colnames(af_mat)]
  if (length(cell_bcs) == 0) stop("No is_cell barcodes present in H5")
  write10xCounts(sprintf("filt_counts_%s.h5", sample_id),
                 af_mat[, cell_bcs, drop = FALSE], version = "3", overwrite = TRUE)
  fwrite(data.table(barcode = cell_bcs),
         sprintf("cell_barcodes_%s.csv", sample_id), col.names = FALSE, quote = FALSE)
  fwrite(data.table(barcode = dt[is_empty == TRUE, barcode]),
         sprintf("empty_barcodes_%s.csv", sample_id), col.names = FALSE, quote = FALSE)

  # per-barcode labels: GMM posteriors + argmax kept for the report's GMM-prediction plot
  pcols <- intersect(paste0("p_", use_pops), names(dt))
  fwrite(dt[, c("barcode", "rank", "total", "spliced", "splice_pct", pcols,
                "gmm_class", "population", "is_cell", "is_empty"), with = FALSE],
         sprintf("cell_calling_labels_%s.csv.gz", sample_id))

  .n <- function(p) { v <- tab[population == p, n]; if (length(v)) v else 0L }
  summ <- data.table(
    sample_id = sample_id, n_total = nrow(dt),
    n_empty = .n("empty"), n_ambient = .n("ambient"), n_damaged = .n("damaged"),
    n_nuclei = .n("nuclei"),
    n_is_empty = sum(dt$is_empty), n_is_cell = sum(dt$is_cell),
    cell_floor_umi = as.integer(round(cell_floor)),
    damaged_splice_cut = round(damaged_cut, 4),
    bhattacharyya = round(sep_bhatt, 3),
    splice_gap_logit = round(sep_splice_gap, 3),
    umi_gap_log10 = round(sep_umi_gap, 3),
    include_damaged_in_empty = include_damaged, gmm_G = G,
    gmm_components = paste(use_pops, collapse = "|"))
  fwrite(summ, sprintf("cell_calling_summary_%s.csv", sample_id))

  saveRDS(list(sample_id = sample_id, use_pops = use_pops,
               mu = par$mean, sigma = par$variance$sigma, pro = par$pro,
               cell_floor = cell_floor, damaged_cut = damaged_cut,
               empty_floor = EMPTY_FLOOR),
          sprintf("cell_calling_gmm_%s.rds", sample_id))

  writeLines(c(
    sprintf("CC_N_CELL=%d", sum(dt$is_cell)),
    sprintf("CC_N_EMPTY=%d", sum(dt$is_empty)),
    sprintf("CC_DAMAGED_SPLICE_CUT=%.4f", damaged_cut),
    sprintf("CC_BHATTACHARYYA=%.4f", sep_bhatt),
    sprintf("CC_SPLICE_GAP=%.4f", sep_splice_gap)),
    sprintf("cell_calling_%s.env", sample_id))
  message("Done: ", sample_id)
}

run_cell_calling()
