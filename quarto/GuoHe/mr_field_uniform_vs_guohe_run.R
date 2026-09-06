# mr_field_uniform_vs_guohe_run.R
#
# Stage 2b driver for the uniform (kappa) calibration check set
# (dev/tasks/TASK_mr_field_uniform_2026-09-05.md, H5 as adjudicated at
# Gate 1: t35 beta2 = 0.3/0.5, t6 k = 10, t7 beta2 = 0.0; campaign tag
# "uniform").  Derived from mr_field_vs_guohe_run.R with three deltas:
#
#   1. The gate call adds field_uniform = TRUE (kappa sweep under the gate's
#      derived seed + 910000L; the field's own draws are untouched, so every
#      pre-existing column must reproduce the 2026-09-05 field bundle).
#   2. Gate 2b pairing proof: MR-current est/se_ij/bias terms AND the field
#      block (est2, lambda mean/sd, all seven quantiles, both intervals,
#      draw-usage counts) are asserted identical() against the committed
#      mr_field_vs_guohe_<id>.rds row -- "all non-uniform columns identical
#      to the 2026-09-05 field bundles".
#   3. Output stem mr_field_uniform_vs_guohe_<id>.rds (skip-if-exists;
#      --force).  The committed bundles are read-only inputs.
#
# Usage:
#   Rscript quarto/GuoHe/mr_field_uniform_vs_guohe_run.R --cores=100
#   Rscript quarto/GuoHe/mr_field_uniform_vs_guohe_run.R --cells=t7_beta2_00 --force

suppressMessages(library(survival))
suppressMessages(library(forestsearch))

.mu_dir <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(normalizePath(f[1])) else normalizePath(getwd())
})
source(file.path(.mu_dir, "mr_vs_guohe_sim.R"))

RNGkind("Mersenne-Twister", "Inversion", "Rejection")

MU_R_OUT <- 1000L   # E2 defaults, unchanged from the field arc
MU_R_IN  <- 500L

args <- commandArgs(trailingOnly = TRUE)
flag <- function(nm) any(args == paste0("--", nm))
opt <- function(nm, default) {
  hit <- grep(paste0("^--", nm, "="), args, value = TRUE)
  if (!length(hit)) default else sub(paste0("^--", nm, "="), "", hit[1])
}
n_cores <- as.integer(opt("cores", "100"))
force <- flag("force")
cells_opt <- opt("cells", "")

# H5 (Gate 1 adjudication): the four check cells.
MU_CELLS <- c("t35_beta2_03", "t35_beta2_05", "t6_k10", "t7_beta2_00")
if (nzchar(cells_opt)) MU_CELLS <- strsplit(cells_opt, ",")[[1]]

# Field-bundle columns that must reproduce exactly (same seeds; the kappa
# sweep runs after the field stream).
MU_FIELD_ID_COLS <- c("mr_est", "mr_se_ij", "mr_bias_sel", "mr_bias_fix",
                      "fld_est2", "fld_lambda_mean", "fld_lambda_sd",
                      "fld_q05", "fld_q25", "fld_q50", "fld_q75", "fld_q95",
                      "fld_q025", "fld_q975", "fld_lower_1s", "fld_lower_2s",
                      "fld_upper_2s", "fld_n_out_used")

.mu_uni_cols <- function(mr, theta) {
  u <- mr$field$uniform
  if (is.null(u) || !is.null(u$note))
    return(list(uni_na = 1L, uni_kappa = NA_real_, uni_kappa_mcse = NA_real_,
                uni_M = NA_integer_, uni_mass = NA_real_, uni_minC1 = NA_real_,
                uni_lower_2u = NA_real_, uni_upper_2u = NA_real_,
                uni_cover_2u = NA_integer_, uni_secs = NA_real_))
  lo <- log(u$lower_2u); hi <- log(u$upper_2u)
  list(uni_na = 0L, uni_kappa = u$kappa, uni_kappa_mcse = u$kappa_mcse,
       uni_M = u$M, uni_mass = u$mass_covered, uni_minC1 = u$minC1,
       uni_lower_2u = lo, uni_upper_2u = hi,
       uni_cover_2u = as.integer(lo <= theta && theta <= hi),
       uni_secs = u$timing_seconds)
}

.mu_field_row <- function(mr, theta) {
  f <- mr$field
  list(mr_est = log(mr$debiased$est), mr_se_ij = mr$debiased$se_ij,
       mr_bias_sel = mr$selection_bias, mr_bias_fix = mr$fixed_bias,
       fld_est2 = log(f$est2), fld_lambda_mean = f$lambda_mean,
       fld_lambda_sd = f$lambda_sd,
       fld_q05 = f$q05, fld_q25 = f$q25, fld_q50 = f$q50, fld_q75 = f$q75,
       fld_q95 = f$q95, fld_q025 = f$q025, fld_q975 = f$q975,
       fld_lower_1s = log(f$lower_1s), fld_lower_2s = log(f$lower_2s),
       fld_upper_2s = log(f$upper_2s), fld_n_out_used = f$n_out_used,
       fld_cover_1s = as.integer(log(f$lower_1s) <= theta),
       fld_cover_2s = as.integer(log(f$lower_2s) <= theta &&
                                   theta <= log(f$upper_2s)))
}

mu_rep <- function(id, m, row_f, base, sec52, sc, tru) {
  if (sec52) {
    cand <- mv_gh52_regen(id, m)
    fits <- gh52_subgroup_fits(cand$df, cand)
    sel <- which.max(replace(fits$est, !is.finite(fits$est), -Inf))
    theta <- row_f$theta
    t0 <- proc.time()[["elapsed"]]
    mr <- mv_mr(cand$df, mv_cand_idx_52(cand), cand$names[sel], mv_spec52,
                seed = base + m + MV_SEED_MR, ci_method = "field",
                field_R_out = MU_R_OUT, field_R_in = MU_R_IN,
                field_uniform = TRUE)
    t1 <- proc.time()[["elapsed"]]
    kfam <- length(cand$cuts)
  } else {
    df <- mv_gh51_regen(id, m)
    k <- length(sc$beta)
    fits <- gh_subgroup_fits(df, sc$beta)
    nv <- gh_naive(fits, sc$beta)
    sel <- nv$sel
    theta <- row_f$theta
    t0 <- proc.time()[["elapsed"]]
    mr <- mv_mr(df, mv_cand_idx_51(df, k), paste0("S", sel), mv_spec51,
                seed = base + m + MV_SEED_MR, ci_method = "field",
                field_R_out = MU_R_OUT, field_R_in = MU_R_IN,
                field_uniform = TRUE)
    t1 <- proc.time()[["elapsed"]]
    kfam <- k
  }
  fr <- .mu_field_row(mr, theta)
  # Gate 2b pairing proof: every pre-existing column identical() to the
  # committed 2026-09-05 field-bundle row.
  pair_ok <- all(vapply(MU_FIELD_ID_COLS,
                        function(cc) identical(fr[[cc]], row_f[[cc]]),
                        logical(1)))
  as.data.frame(c(list(id = id, m = m, seed_data = base + m,
                       seed_mr = base + m + MV_SEED_MR,
                       k_family = kfam, theta = theta, sel = sel,
                       pair_ok = as.integer(pair_ok)),
                  fr, .mu_uni_cols(mr, theta),
                  list(t_gate_s = t1 - t0)),
                stringsAsFactors = FALSE)
}

mu_run_cell <- function(id) {
  f_out <- file.path(.mu_dir, paste0("mr_field_uniform_vs_guohe_", id, ".rds"))
  if (file.exists(f_out) && !force) {
    cat(sprintf("[skip] %s (exists)\n", id)); return(invisible(NULL))
  }
  sec52 <- grepl("^t7_", id)
  fld_bun <- readRDS(file.path(.mu_dir, paste0("mr_field_vs_guohe_", id, ".rds")))
  base <- if (sec52) mv_gh52_base(id) else mv_gh51_base(id)
  stopifnot(fld_bun$seed_base == base,
            identical(fld_bun$results$m, seq_len(nrow(fld_bun$results))))
  n_rep <- nrow(fld_bun$results)
  tru <- NULL; sc <- NULL
  if (sec52) {
    b2 <- as.integer(sub("^t7_beta2_", "", id))
    tru <- readRDS(file.path(.mu_dir, sprintf("guohe_sec52_truth_beta2_%02d.rds", b2)))
  } else sc <- mv_gh51_scenario(id)

  t0 <- proc.time()[["elapsed"]]
  cat(sprintf("[run ] %s  %d reps ...\n", id, n_rep)); utils::flush.console()
  rows <- parallel::mclapply(
    seq_len(n_rep),
    function(m) try(mu_rep(id, m, fld_bun$results[m, ], base, sec52, sc, tru),
                    silent = TRUE),
    mc.cores = n_cores, mc.preschedule = FALSE)
  bad <- vapply(rows, function(z) inherits(z, "try-error"), logical(1))
  if (any(bad))
    warning(sum(bad), " replicate(s) errored in ", id, ". First: ",
            as.character(rows[[which(bad)[1]]]))
  res <- do.call(rbind, rows[!bad])
  el <- proc.time()[["elapsed"]] - t0
  gate2 <- list(n_rep_expected = n_rep, n_rep_done = nrow(res),
                n_errored = sum(bad), pair_mismatch = sum(res$pair_ok == 0L),
                uni_na = sum(res$uni_na))
  saveRDS(list(
    id = id, section = if (sec52) "5.2" else "5.1", campaign = "uniform",
    source_bundle = paste0("mr_field_vs_guohe_", id, ".rds"),
    pair_cols = MU_FIELD_ID_COLS,
    draws = MV_DRAWS, multiplier = MV_MULTIPLIER,
    field_R_out = MU_R_OUT, field_R_in = MU_R_IN,
    uniform_seed_offset = 910000L,
    uniform_settings = list(delta_grid = seq(0, 4, by = 0.5), mass = 0.99,
                            M_cap = 40L, R_rep = 300L, R_out = 300L,
                            R_in = 150L, kappa_grid = "seq(1, 3, by = 0.01)"),
    seed_base = base, mr_seed_offset = MV_SEED_MR,
    gate2 = gate2, elapsed_sec = el,
    sessionInfo = utils::capture.output(utils::sessionInfo()),
    results = res), f_out)
  cat(sprintf("[done] %s  %d/%d reps in %.1f min  pair_mm %d  uni_na %d  -> %s\n",
              id, gate2$n_rep_done, n_rep, el / 60, gate2$pair_mismatch,
              gate2$uni_na, basename(f_out)))
  utils::flush.console()
  invisible(NULL)
}

cat(sprintf("Uniform Stage 2b run: %d cells, cores = %d, field R_out/R_in = %d/%d, kappa sweep at adjudicated H2/H3 (M_cap 40, grid [1,3])\n\n",
            length(MU_CELLS), n_cores, MU_R_OUT, MU_R_IN))
for (id in MU_CELLS) mu_run_cell(id)

cat("\n=== GATE 2b TALLY ===\n")
for (id in MU_CELLS) {
  f <- file.path(.mu_dir, paste0("mr_field_uniform_vs_guohe_", id, ".rds"))
  if (!file.exists(f)) { cat(sprintf("%-16s MISSING\n", id)); next }
  g <- readRDS(f)$gate2
  cat(sprintf("%-16s reps %d/%d  errored %d  pair_mm %d  uni_na %d\n",
              id, g$n_rep_done, g$n_rep_expected, g$n_errored,
              g$pair_mismatch, g$uni_na))
}
