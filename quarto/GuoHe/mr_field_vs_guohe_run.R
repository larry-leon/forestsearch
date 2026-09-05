# mr_field_vs_guohe_run.R
#
# Stage 2 driver for the MR (field) three-way comparison
# (dev/tasks/TASK_mr_field_vs_guohe_2026-09-05.md; E1-E5 at defaults).
#
# Per replicate, paired on the SAME seeds as the 2026-09-04 run:
#   1. Regenerate the data bit-for-bit; recompute the naive column and assert
#      identical() against the replication bundle (as on 2026-09-04).
#   2. ONE gate call with ci_method = "field" (MR seed base + m + 700000L,
#      unchanged; field draws under the derived seed + 900000L inside the
#      gate).  The debiased element is IJ-computed exactly as under "ij", so
#      the same call yields the current MR row.
#   3. Pairing proof for the new run: MR-current est / se_ij / bias_sel /
#      bias_fix and p_hat_H must be identical() to the 2026-09-04
#      mr_vs_guohe_<id>.rds row (same seed -> same Xi -> same values).
#   4. Record the field fields and the two MR (field) coverage indicators
#      (one-sided at q95, primary per E4; two-sided quantile supplementary).
#   5. NO M_eff pass (E5): the addendum-A columns (p_hat_*, Sigma_HH, A6_*,
#      m0_*, M_eff, tie_resid_implied) are joined from the 2026-09-04 bundle
#      by (id, m) at cell level, with seed_data equality asserted row-wise.
#
# Output: mr_field_vs_guohe_<id>.rds per cell (skip-if-exists; --force), with
# per-replicate results, gate tallies, seeds/settings, elapsed, sessionInfo.
# The 2026-09-04 bundles are read-only inputs and are not modified.
#
# Usage:
#   Rscript quarto/GuoHe/mr_field_vs_guohe_run.R --cores=120
#   Rscript quarto/GuoHe/mr_field_vs_guohe_run.R --cells=t7_beta2_00 --force

suppressMessages(library(survival))
suppressMessages(library(forestsearch))

.mf2_dir <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(normalizePath(f[1])) else normalizePath(getwd())
})
source(file.path(.mf2_dir, "mr_vs_guohe_sim.R"))

RNGkind("Mersenne-Twister", "Inversion", "Rejection")

MF_R_OUT <- 1000L   # E2 defaults
MF_R_IN  <- 500L

args <- commandArgs(trailingOnly = TRUE)
flag <- function(nm) any(args == paste0("--", nm))
opt <- function(nm, default) {
  hit <- grep(paste0("^--", nm, "="), args, value = TRUE)
  if (!length(hit)) default else sub(paste0("^--", nm, "="), "", hit[1])
}
n_cores <- as.integer(opt("cores", "120"))
force <- flag("force")
cells_opt <- opt("cells", "")

MF_CELLS <- c(sprintf("t35_beta2_%02d", 0:5),
              sprintf("t6_k%02d", c(2L, 6L, 10L, 12L)),
              sprintf("t7_beta2_%02d", 0:5))
if (nzchar(cells_opt)) MF_CELLS <- strsplit(cells_opt, ",")[[1]]

# Addendum-A columns joined from the 2026-09-04 bundles (E5)
MF_JOIN_COLS <- c("sel_agree_mr", "p_hat_H", "p_top1", "p_top2", "p_top3",
                  "p_lab1", "p_lab2", "p_lab3", "Sigma_HH", "A6_mass",
                  "A6_mass_std", "m0_hat", "m0_mc_se", "M_eff",
                  "tie_resid_implied")

.mf_gh_cols <- function(row, theta) {
  out <- list()
  for (i in 1:4) {
    out[[sprintf("gh_r%d_deb", i)]] <- row[[sprintf("r%d_bias", i)]] + theta
    out[[sprintf("gh_r%d_low", i)]] <- theta - row[[sprintf("r%d_dist", i)]]
    out[[sprintf("gh_r%d_cov", i)]] <- row[[sprintf("r%d_cover", i)]]
  }
  out
}

.mf_mr_cols <- function(mr, row_c, theta, t_mr) {
  est <- log(mr$debiased$est)
  lo1 <- log(mr$debiased$lower_1s)
  # pairing proof vs the 2026-09-04 row (same seed -> identical values)
  cur_ok <- identical(est, row_c$mr_est) &&
    identical(mr$debiased$se_ij, row_c$mr_se_ij) &&
    identical(mr$selection_bias, row_c$mr_bias_sel) &&
    identical(mr$fixed_bias, row_c$mr_bias_fix) &&
    identical(unname(mr$reselection$p_hat[mr$selected_index]), row_c$p_hat_H)
  f <- mr$field
  e2 <- log(f$est2); flo1 <- log(f$lower_1s)
  flo2 <- log(f$lower_2s); fup2 <- log(f$upper_2s)
  list(
    cur_ok = as.integer(cur_ok),
    mr_est = est, mr_se_ij = mr$debiased$se_ij, mr_se_wald = mr$debiased$se_wald,
    mr_ij_source = mr$debiased$ij_source, mr_ij_draws = mr$debiased$ij_draws,
    mr_lower_1s = lo1,
    mr_lower_2s = log(mr$debiased$lower), mr_upper_2s = log(mr$debiased$upper),
    mr_cover = as.integer(lo1 <= theta),
    mr_bias_sel = mr$selection_bias, mr_bias_fix = mr$fixed_bias,
    mr_selection_rate = mr$selection_rate, mr_mean_r = mr$mean_r,
    mr_na = as.integer(!(is.finite(est) && is.finite(mr$debiased$se_ij) && is.finite(lo1))),
    fld_lambda_mean = f$lambda_mean, fld_lambda_sd = f$lambda_sd,
    fld_q05 = f$q05, fld_q25 = f$q25, fld_q50 = f$q50, fld_q75 = f$q75,
    fld_q95 = f$q95, fld_q025 = f$q025, fld_q975 = f$q975,
    fld_n_out_used = f$n_out_used, fld_n_in_used_mean = f$n_in_used_mean,
    fld_est2 = e2, fld_lower_1s = flo1,
    fld_lower_2s = flo2, fld_upper_2s = fup2,
    fld_lower_se = log(f$lower_se), fld_upper_se = log(f$upper_se),
    fld_cover_1s = as.integer(flo1 <= theta),
    fld_cover_2s = as.integer(flo2 <= theta && theta <= fup2),
    fld_na = as.integer(!(is.finite(e2) && is.finite(flo1) && is.finite(f$lambda_sd))),
    t_field_s = f$timing_seconds, t_mr_s = t_mr)
}

mf_rep_51 <- function(id, m, row_r, row_c, beta, n, base) {
  df <- mv_gh51_regen(id, m)
  k <- length(beta)
  fits <- gh_subgroup_fits(df, beta)
  nv <- gh_naive(fits, beta)
  naive_ok <- identical(nv$sel, row_r$naive_sel) &&
    identical(nv$beta_s, row_r$naive_beta_s) &&
    identical(nv$cover, row_r$naive_cover) &&
    identical(nv$dist, row_r$naive_dist) &&
    identical(nv$bias, row_r$naive_bias)
  theta <- nv$beta_s
  t0 <- proc.time()[["elapsed"]]
  mr <- mv_mr(df, mv_cand_idx_51(df, k), paste0("S", nv$sel), mv_spec51,
              seed = base + m + MV_SEED_MR, ci_method = "field",
              field_R_out = MF_R_OUT, field_R_in = MF_R_IN)
  t1 <- proc.time()[["elapsed"]]
  as.data.frame(c(list(
    id = id, m = m, seed_data = base + m, seed_mr = base + m + MV_SEED_MR,
    k_family = k, theta = theta, sel = nv$sel,
    naive_ok = as.integer(naive_ok),
    naive_est = nv$point, naive_se = fits$se[nv$sel],
    naive_lower = nv$lower, naive_cover = nv$cover),
    .mf_gh_cols(row_r, theta),
    .mf_mr_cols(mr, row_c, theta, t1 - t0)),
    stringsAsFactors = FALSE)
}

mf_rep_52 <- function(id, m, row_r, row_c, base, tru) {
  cand <- mv_gh52_regen(id, m)
  fits <- gh52_subgroup_fits(cand$df, cand)
  nv <- gh52_naive(fits, cand, tru)
  naive_ok <- identical(nv$point, row_r$naive_point) &&
    identical(nv$lower, row_r$naive_lower) &&
    identical(nv$c_hat, row_r$c_hat_naive) &&
    identical(nv$gamma_s, row_r$gamma_s_naive) &&
    identical(nv$cover, row_r$naive_cover)
  sel <- which.max(replace(fits$est, !is.finite(fits$est), -Inf))
  sel_ok <- identical(cand$cuts[sel], row_r$c_hat_gh) &&
    identical(gh52_truth_at(tru, cand$cuts[sel]), row_r$gamma_s)
  theta <- row_r$gamma_s
  t0 <- proc.time()[["elapsed"]]
  mr <- mv_mr(cand$df, mv_cand_idx_52(cand), cand$names[sel], mv_spec52,
              seed = base + m + MV_SEED_MR, ci_method = "field",
              field_R_out = MF_R_OUT, field_R_in = MF_R_IN)
  t1 <- proc.time()[["elapsed"]]
  as.data.frame(c(list(
    id = id, m = m, seed_data = base + m, seed_mr = base + m + MV_SEED_MR,
    k_family = length(cand$cuts), theta = theta, sel = sel,
    c_hat = cand$cuts[sel], n_sel = row_r$n_sel,
    naive_ok = as.integer(naive_ok && sel_ok),
    naive_est = nv$point, naive_se = (nv$point - nv$lower) / stats::qnorm(0.95),
    naive_lower = nv$lower, naive_cover = nv$cover,
    gamma_s_naive = nv$gamma_s),
    .mf_gh_cols(row_r, theta),
    .mf_mr_cols(mr, row_c, theta, t1 - t0)),
    stringsAsFactors = FALSE)
}

mf_run_cell <- function(id) {
  f_out <- file.path(.mf2_dir, paste0("mr_field_vs_guohe_", id, ".rds"))
  if (file.exists(f_out) && !force) {
    cat(sprintf("[skip] %s (exists)\n", id))
    return(invisible(NULL))
  }
  sec52 <- grepl("^t7_", id)
  rep_bun <- readRDS(file.path(.mf2_dir, paste0("guohe_repro_", id, ".rds")))
  cmp_bun <- readRDS(file.path(.mf2_dir, paste0("mr_vs_guohe_", id, ".rds")))
  base <- if (sec52) mv_gh52_base(id) else mv_gh51_base(id)
  stopifnot(rep_bun$seed_base == base, cmp_bun$seed_base == base,
            nrow(cmp_bun$results) == nrow(rep_bun$results),
            identical(cmp_bun$results$m, seq_len(nrow(cmp_bun$results))))
  n_rep <- nrow(rep_bun$results)
  tru <- NULL; sc <- NULL
  if (sec52) {
    b2 <- as.integer(sub("^t7_beta2_", "", id))
    tru <- readRDS(file.path(.mf2_dir, sprintf("guohe_sec52_truth_beta2_%02d.rds", b2)))
  } else sc <- mv_gh51_scenario(id)

  t0 <- proc.time()[["elapsed"]]
  cat(sprintf("[run ] %s  %d reps ...\n", id, n_rep)); utils::flush.console()
  rows <- parallel::mclapply(
    seq_len(n_rep),
    function(m) {
      try(if (sec52) mf_rep_52(id, m, rep_bun$results[m, ], cmp_bun$results[m, ], base, tru)
          else mf_rep_51(id, m, rep_bun$results[m, ], cmp_bun$results[m, ],
                         sc$beta, sc$n, base),
          silent = TRUE)
    },
    mc.cores = n_cores, mc.preschedule = FALSE)
  bad <- vapply(rows, function(z) inherits(z, "try-error"), logical(1))
  if (any(bad)) {
    warning(sum(bad), " replicate(s) errored in ", id, ". First: ",
            as.character(rows[[which(bad)[1]]]))
  }
  res <- do.call(rbind, rows[!bad])
  # E5 join: addendum-A columns from the 2026-09-04 bundle by (id, m);
  # row order is m for both (asserted above), seed equality asserted here.
  old <- cmp_bun$results[res$m, ]
  stopifnot(identical(old$seed_data, res$seed_data))
  res <- cbind(res, old[, MF_JOIN_COLS])
  el <- proc.time()[["elapsed"]] - t0
  gate2 <- list(n_rep_expected = n_rep, n_rep_done = nrow(res),
                n_errored = sum(bad),
                naive_mismatch = sum(res$naive_ok == 0L),
                cur_mismatch = sum(res$cur_ok == 0L),
                mr_na = sum(res$mr_na), fld_na = sum(res$fld_na))
  saveRDS(list(
    id = id, section = if (sec52) "5.2" else "5.1",
    source_bundle = paste0("guohe_repro_", id, ".rds"),
    join_bundle = paste0("mr_vs_guohe_", id, ".rds"),
    join_cols = MF_JOIN_COLS,
    draws = MV_DRAWS, multiplier = MV_MULTIPLIER,
    field_R_out = MF_R_OUT, field_R_in = MF_R_IN, field_seed_offset = 900000L,
    seed_base = base, mr_seed_offset = MV_SEED_MR,
    gate2 = gate2, elapsed_sec = el,
    sessionInfo = utils::capture.output(utils::sessionInfo()),
    results = res), f_out)
  cat(sprintf("[done] %s  %d/%d reps in %.1f min  naive_mm %d  cur_mm %d  mr_na %d  fld_na %d  -> %s\n",
              id, gate2$n_rep_done, n_rep, el / 60, gate2$naive_mismatch,
              gate2$cur_mismatch, gate2$mr_na, gate2$fld_na, basename(f_out)))
  utils::flush.console()
  invisible(NULL)
}

cat(sprintf("MR (field) Stage 2 run: %d cells, cores = %d, MR draws = %d (%s), field R_out/R_in = %d/%d\n\n",
            length(MF_CELLS), n_cores, MV_DRAWS, MV_MULTIPLIER, MF_R_OUT, MF_R_IN))
for (id in MF_CELLS) mf_run_cell(id)

cat("\n=== GATE 2 TALLY ===\n")
for (id in MF_CELLS) {
  f <- file.path(.mf2_dir, paste0("mr_field_vs_guohe_", id, ".rds"))
  if (!file.exists(f)) { cat(sprintf("%-16s MISSING\n", id)); next }
  g <- readRDS(f)$gate2
  cat(sprintf("%-16s reps %d/%d  errored %d  naive_mm %d  cur_mm %d  mr_na %d  fld_na %d\n",
              id, g$n_rep_done, g$n_rep_expected, g$n_errored,
              g$naive_mismatch, g$cur_mismatch, g$mr_na, g$fld_na))
}
