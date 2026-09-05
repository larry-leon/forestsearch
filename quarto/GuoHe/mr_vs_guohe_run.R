# mr_vs_guohe_run.R
#
# Stage 2 driver for the MR-vs-Guo&He comparison
# (dev/tasks/TASK_mr_vs_guohe_2026-09-04.md + addendum A; compute GO given on
# the Stage 1 projection: 16 cells x 2000 replicates, ~20 core-h).
#
# WHAT RUNS PER REPLICATE (paired by the replication's own seed):
#   1. Regenerate the replicate's data bit-for-bit from the stored seed base.
#   2. Recompute the naive column with the replication's own functions and
#      assert identical() against the stored bundle row (Gate 2 requires this
#      to hold on EVERY replicate of EVERY cell).
#   3. Read the Guo & He per-replicate estimate/bound off the stored columns
#      (debiased = r*_bias + theta, lower = theta - r*_dist).  NO Guo & He
#      re-run (Stage 0/D5).  Sec 5.1 note: the bundles do not store the
#      engine's selected index; theta uses the recomputed argmax, which equals
#      the engine's selection structurally -- the engine ranks the identical
#      full-sample coefficient vector (replication V3a/V3b: engine naive ==
#      independent coxph max to 1e-10, selection == independent argmax) and
#      its min_events = 5 guard cannot trigger at ~100 events per subgroup;
#      in Sec 5.2, where the agreement IS stored, sel_agree = 2000/2000 in
#      all six cells.
#   4. Run MR (fs_mr_inference: maxeff, unrestricted admission, draws = 5000,
#      centred Poisson, IJ interval; seed = base + m + 700000) on the
#      replication's selected candidate, plus the addendum-A record
#      (Sigma-hat, A6 mass, M_eff at 2e5 MC draws, implied tie residual).
#
# PRECONDITION (stop-on-failure, runs before any cell, always): 5 paired
# replicates each of t35_beta2_03 and t6_k12 -- regenerated naive column
# identical() to the stored bundle values and argmax matching the stored
# selection.  Any mismatch aborts the run.
#
# Output: one mr_vs_guohe_<id>.rds per cell (skip-if-exists; --force), with
# per-replicate results, gate tallies, seeds, elapsed, sessionInfo.
#
# Usage:
#   Rscript quarto/GuoHe/mr_vs_guohe_run.R --cores=120
#   Rscript quarto/GuoHe/mr_vs_guohe_run.R --precheck-only
#   Rscript quarto/GuoHe/mr_vs_guohe_run.R --cells=t7_beta2_00 --force

suppressMessages(library(survival))
suppressMessages(library(forestsearch))

.mv2_dir <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(normalizePath(f[1])) else normalizePath(getwd())
})
source(file.path(.mv2_dir, "mr_vs_guohe_sim.R"))  # adapter (sources both replication harnesses)

RNGkind("Mersenne-Twister", "Inversion", "Rejection")

MV2_M0_DRAWS <- 2e5        # addendum A1.4 Monte Carlo size (GO as projected)
MV2_M0_SEED  <- 20260904L  # fixed; M_eff is then deterministic given R-hat

# ---- arguments -------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
flag <- function(nm) any(args == paste0("--", nm))
opt <- function(nm, default) {
  hit <- grep(paste0("^--", nm, "="), args, value = TRUE)
  if (!length(hit)) default else sub(paste0("^--", nm, "="), "", hit[1])
}
n_cores <- as.integer(opt("cores", "120"))
force <- flag("force")
precheck_only <- flag("precheck-only")
cells_opt <- opt("cells", "")

MV2_CELLS <- c(sprintf("t35_beta2_%02d", 0:5),
               sprintf("t6_k%02d", c(2L, 6L, 10L, 12L)),
               sprintf("t7_beta2_%02d", 0:5))
if (nzchar(cells_opt)) MV2_CELLS <- strsplit(cells_opt, ",")[[1]]

# ---- per-replicate workers -------------------------------------------------

.mv2_gh_cols <- function(row, theta) {
  out <- list()
  for (i in 1:4) {
    out[[sprintf("gh_r%d_deb", i)]] <- row[[sprintf("r%d_bias", i)]] + theta
    out[[sprintf("gh_r%d_low", i)]] <- theta - row[[sprintf("r%d_dist", i)]]
    out[[sprintf("gh_r%d_cov", i)]] <- row[[sprintf("r%d_cover", i)]]
  }
  out
}

.mv2_mr_cols <- function(mr, arec, theta, t_mr, t_arec) {
  est <- log(mr$debiased$est)
  lo1 <- log(mr$debiased$lower_1s)
  c(list(
    mr_est = est, mr_se_ij = mr$debiased$se_ij, mr_se_wald = mr$debiased$se_wald,
    mr_ij_source = mr$debiased$ij_source, mr_ij_draws = mr$debiased$ij_draws,
    mr_lower_1s = lo1,
    mr_lower_2s = log(mr$debiased$lower), mr_upper_2s = log(mr$debiased$upper),
    mr_cover = as.integer(lo1 <= theta),
    mr_bias_sel = mr$selection_bias, mr_bias_fix = mr$fixed_bias,
    mr_selection_rate = mr$selection_rate, mr_mean_r = mr$mean_r,
    mr_na = as.integer(!(is.finite(est) && is.finite(mr$debiased$se_ij) && is.finite(lo1)))),
    list(
      p_hat_H = arec$p_hat_H,
      p_top1 = arec$p_hat_top3[1], p_top2 = arec$p_hat_top3[2], p_top3 = arec$p_hat_top3[3],
      p_lab1 = arec$p_hat_top3_labels[1], p_lab2 = arec$p_hat_top3_labels[2],
      p_lab3 = arec$p_hat_top3_labels[3],
      Sigma_HH = arec$Sigma_HH, A6_mass = arec$A6_mass, A6_mass_std = arec$A6_mass_std,
      m0_hat = arec$m0_hat, m0_mc_se = arec$m0_mc_se, M_eff = arec$M_eff,
      tie_resid_implied = arec$tie_resid_implied,
      t_mr_s = t_mr, t_arec_s = t_arec))
}

mv2_rep_51 <- function(id, m, row, beta, n, base) {
  df <- mv_gh51_regen(id, m)
  k <- length(beta)
  fits <- gh_subgroup_fits(df, beta)
  nv <- gh_naive(fits, beta)
  naive_ok <- identical(nv$sel, row$naive_sel) &&
    identical(nv$beta_s, row$naive_beta_s) &&
    identical(nv$cover, row$naive_cover) &&
    identical(nv$dist, row$naive_dist) &&
    identical(nv$bias, row$naive_bias)
  theta <- nv$beta_s
  ci <- mv_cand_idx_51(df, k)
  t0 <- proc.time()[["elapsed"]]
  asm <- mv_assemble(df, ci, mv_spec51)
  mr <- mv_mr(df, ci, names(ci)[nv$sel], mv_spec51, seed = base + m + MV_SEED_MR)
  t1 <- proc.time()[["elapsed"]]
  arec <- mv_a_record(asm, mr, ndraw_m0 = MV2_M0_DRAWS, seed_m0 = MV2_M0_SEED)
  t2 <- proc.time()[["elapsed"]]
  as.data.frame(c(list(
    id = id, m = m, seed_data = base + m, seed_mr = base + m + MV_SEED_MR,
    k_family = k, theta = theta, sel = nv$sel,
    sel_agree_mr = as.integer(which.max(asm$beta_hat) == nv$sel),
    naive_ok = as.integer(naive_ok),
    naive_est = nv$point, naive_se = fits$se[nv$sel],
    naive_lower = nv$lower, naive_cover = nv$cover),
    .mv2_gh_cols(row, theta),
    .mv2_mr_cols(mr, arec, theta, t1 - t0, t2 - t1)),
    stringsAsFactors = FALSE)
}

mv2_rep_52 <- function(id, m, row, base, tru) {
  cand <- mv_gh52_regen(id, m)
  fits <- gh52_subgroup_fits(cand$df, cand)
  nv <- gh52_naive(fits, cand, tru)
  naive_ok <- identical(nv$point, row$naive_point) &&
    identical(nv$lower, row$naive_lower) &&
    identical(nv$c_hat, row$c_hat_naive) &&
    identical(nv$gamma_s, row$gamma_s_naive) &&
    identical(nv$cover, row$naive_cover)
  sel <- which.max(replace(fits$est, !is.finite(fits$est), -Inf))
  sel_ok <- identical(cand$cuts[sel], row$c_hat_gh) &&
    identical(gh52_truth_at(tru, cand$cuts[sel]), row$gamma_s)
  theta <- row$gamma_s
  ci <- mv_cand_idx_52(cand)
  t0 <- proc.time()[["elapsed"]]
  asm <- mv_assemble(cand$df, ci, mv_spec52)
  mr <- mv_mr(cand$df, ci, cand$names[sel], mv_spec52, seed = base + m + MV_SEED_MR)
  t1 <- proc.time()[["elapsed"]]
  arec <- mv_a_record(asm, mr, ndraw_m0 = MV2_M0_DRAWS, seed_m0 = MV2_M0_SEED)
  t2 <- proc.time()[["elapsed"]]
  as.data.frame(c(list(
    id = id, m = m, seed_data = base + m, seed_mr = base + m + MV_SEED_MR,
    k_family = length(cand$cuts), theta = theta, sel = sel,
    c_hat = cand$cuts[sel], n_sel = row$n_sel,
    sel_agree_mr = as.integer(which.max(asm$beta_hat) == sel),
    naive_ok = as.integer(naive_ok && sel_ok),
    naive_est = nv$point, naive_se = (nv$point - nv$lower) / stats::qnorm(0.95),
    naive_lower = nv$lower, naive_cover = nv$cover,
    gamma_s_naive = nv$gamma_s),
    .mv2_gh_cols(row, theta),
    .mv2_mr_cols(mr, arec, theta, t1 - t0, t2 - t1)),
    stringsAsFactors = FALSE)
}

# ---- precondition: paired naive identity on 5 + 5 replicates ---------------
# Stop-on-failure BEFORE any full cell, per the compute-GO message.

mv2_precheck <- function() {
  cat("=== PRECONDITION: naive identity, t35_beta2_03 and t6_k12, m = 1..5 ===\n")
  all_ok <- TRUE
  for (id in c("t35_beta2_03", "t6_k12")) {
    bun <- readRDS(file.path(.mv2_dir, paste0("guohe_repro_", id, ".rds")))
    sc <- mv_gh51_scenario(id)
    base <- mv_gh51_base(id)
    stopifnot(bun$seed_base == base)
    for (m in 1:5) {
      df <- mv_gh51_regen(id, m)
      fits <- gh_subgroup_fits(df, sc$beta)
      nv <- gh_naive(fits, sc$beta)
      row <- bun$results[m, ]
      ok_n <- identical(nv$beta_s, row$naive_beta_s) &&
        identical(nv$cover, row$naive_cover) &&
        identical(nv$dist, row$naive_dist) &&
        identical(nv$bias, row$naive_bias)
      ok_s <- identical(nv$sel, row$naive_sel)
      cat(sprintf("  %s m=%d  naive %s  argmax %s (sel %d, stored %d)\n",
                  id, m, if (ok_n) "identical" else "MISMATCH",
                  if (ok_s) "match" else "MISMATCH", nv$sel, row$naive_sel))
      all_ok <- all_ok && ok_n && ok_s
    }
  }
  if (!all_ok) stop("PRECONDITION FAILED -- stopping before any Stage 2 compute.",
                    call. = FALSE)
  cat("PRECONDITION PASS\n\n")
}

# ---- run -------------------------------------------------------------------

mv2_run_cell <- function(id) {
  f_out <- file.path(.mv2_dir, paste0("mr_vs_guohe_", id, ".rds"))
  if (file.exists(f_out) && !force) {
    cat(sprintf("[skip] %s (exists)\n", id))
    return(invisible(NULL))
  }
  sec52 <- grepl("^t7_", id)
  bun <- readRDS(file.path(.mv2_dir, paste0("guohe_repro_", id, ".rds")))
  base <- if (sec52) mv_gh52_base(id) else mv_gh51_base(id)
  stopifnot(bun$seed_base == base, bun$n_rep_used == nrow(bun$results))
  n_rep <- nrow(bun$results)
  tru <- NULL
  sc <- NULL
  if (sec52) {
    b2 <- as.integer(sub("^t7_beta2_", "", id))
    tru <- readRDS(file.path(.mv2_dir, sprintf("guohe_sec52_truth_beta2_%02d.rds", b2)))
  } else sc <- mv_gh51_scenario(id)

  t0 <- proc.time()[["elapsed"]]
  cat(sprintf("[run ] %s  %d reps ...\n", id, n_rep)); utils::flush.console()
  rows <- parallel::mclapply(
    seq_len(n_rep),
    function(m) {
      try(if (sec52) mv2_rep_52(id, m, bun$results[m, ], base, tru)
          else mv2_rep_51(id, m, bun$results[m, ], sc$beta, sc$n, base),
          silent = TRUE)
    },
    mc.cores = n_cores, mc.preschedule = FALSE)
  bad <- vapply(rows, function(z) inherits(z, "try-error"), logical(1))
  if (any(bad)) {
    warning(sum(bad), " replicate(s) errored in ", id, ". First: ",
            as.character(rows[[which(bad)[1]]]))
  }
  res <- do.call(rbind, rows[!bad])
  el <- proc.time()[["elapsed"]] - t0
  gate2 <- list(n_rep_expected = n_rep, n_rep_done = nrow(res),
                n_errored = sum(bad),
                naive_mismatch = sum(res$naive_ok == 0L),
                mr_na = sum(res$mr_na),
                sel_agree_mr = sum(res$sel_agree_mr))
  saveRDS(list(
    id = id, section = if (sec52) "5.2" else "5.1",
    source_bundle = paste0("guohe_repro_", id, ".rds"),
    draws = MV_DRAWS, multiplier = MV_MULTIPLIER,
    m0_draws = MV2_M0_DRAWS, m0_seed = MV2_M0_SEED,
    seed_base = base, mr_seed_offset = MV_SEED_MR,
    gate2 = gate2, elapsed_sec = el,
    sessionInfo = utils::capture.output(utils::sessionInfo()),
    results = res), f_out)
  cat(sprintf("[done] %s  %d/%d reps in %.1f min  naive_mismatch %d  mr_na %d  -> %s\n",
              id, gate2$n_rep_done, n_rep, el / 60, gate2$naive_mismatch,
              gate2$mr_na, basename(f_out)))
  utils::flush.console()
  invisible(NULL)
}

mv2_precheck()
if (precheck_only) {
  cat("(--precheck-only: stopping here.)\n")
} else {
  cat(sprintf("Stage 2 run: %d cells, cores = %d, draws = %d (%s), m0 draws = %g\n\n",
              length(MV2_CELLS), n_cores, MV_DRAWS, MV_MULTIPLIER, MV2_M0_DRAWS))
  for (id in MV2_CELLS) mv2_run_cell(id)
  cat("\n=== GATE 2 TALLY ===\n")
  for (id in MV2_CELLS) {
    f <- file.path(.mv2_dir, paste0("mr_vs_guohe_", id, ".rds"))
    if (!file.exists(f)) { cat(sprintf("%-16s MISSING\n", id)); next }
    g <- readRDS(f)$gate2
    cat(sprintf("%-16s reps %d/%d  errored %d  naive_mismatch %d  mr_na %d  mr_argmax_agree %d\n",
                id, g$n_rep_done, g$n_rep_expected, g$n_errored,
                g$naive_mismatch, g$mr_na, g$sel_agree_mr))
  }
}
