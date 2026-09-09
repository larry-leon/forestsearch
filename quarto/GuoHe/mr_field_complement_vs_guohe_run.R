# mr_field_complement_vs_guohe_run.R
#
# T1 driver for Tier 1 of the fs-post-selection supplement: the complement and
# joint constructions on t7 (dev/tasks/claude_cc_task_guohe_supplement_2026-09-09.md
# section 4, as amended by ..._v2.md A3/A4).
#
# TRANSPLANT of mr_field_vs_guohe_run.R (50d6b641). Named-line changes only:
#   1. Default cells: t7 only (--cells/--force/skip-if-exists scaffolding kept).
#   2. The gate call is a DIRECT fs_mr_inference() call whose argument-assembly
#      block is copied verbatim from mv_mr() (mr_vs_guohe_sim.R:132-151), plus
#      exactly four additions: field_complement = TRUE, include_complement = TRUE,
#      field_scale_complement = "selected", return_reselection = TRUE.
#      mr_vs_guohe_sim.R itself is NOT modified (A3) -- it has no
#      field_scale_complement formal and is a committed input to two campaigns.
#   3. Emitted columns extended: complement upper bounds at 0.95 and 0.975,
#      harm-side lower bounds at 0.95 and 0.975, the joint both-correct
#      indicator, complement diagnostics, and p-hat(H-hat).
#   4. --pilot scaffolding transplanted from guohe_sec52_run.R:70, 81-84, 144,
#      208-225 (A4).
#   5. Output mr_field_complement_vs_guohe_<id>.rds.
#
# The 0.975 pair and the both-correct indicator are the engine's own, per A2:
# .fs_mr_field_joint (R/fs_mr_inference.R:1195-1220) returns bonf_lower_H,
# bonf_upper_Hc and bonf_joint_prob at gamma = alpha/2 = 0.025, wired for both
# scalings at :1141-1144. Nothing under R/ is changed or needed.
#
# Truth for the complement is EXACTLY 0: guohe_sec52_truth.R:75 sets
# b <- ifelse(w <= GH52_C_LO, beta2, 0) with GH52_C_LO = 30, and the cutpoint
# grid starts at 30, so every H-hat^c = {W > c-hat} lies inside {W > 30} where
# the true effect is identically zero. Complement coverage is therefore
# upper >= 0 with no dilution term.
#
# Seeds are unchanged from the stored bundles (seed_data = base + m;
# seed_mr = base + m + MV_SEED_MR; field draws at the derived seed + 900000L
# inside the gate), so every stored-comparable column must reproduce exactly.
#
# Usage:
#   Rscript quarto/GuoHe/mr_field_complement_vs_guohe_run.R --probe --cells=t7_beta2_00,t7_beta2_05 --reps=3
#   Rscript quarto/GuoHe/mr_field_complement_vs_guohe_run.R --pilot --cores=10
#   Rscript quarto/GuoHe/mr_field_complement_vs_guohe_run.R --cores=10
#   Rscript quarto/GuoHe/mr_field_complement_vs_guohe_run.R --cells=t7_beta2_00 --force

suppressMessages(library(survival))
suppressMessages(library(forestsearch))

.mfc_dir <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(normalizePath(f[1])) else normalizePath(getwd())
})
source(file.path(.mfc_dir, "mr_vs_guohe_sim.R"))

RNGkind("Mersenne-Twister", "Inversion", "Rejection")

MF_R_OUT <- 1000L   # E2 defaults
MF_R_IN  <- 500L

args <- commandArgs(trailingOnly = TRUE)
flag <- function(nm) any(args == paste0("--", nm))
opt <- function(nm, default) {
  hit <- grep(paste0("^--", nm, "="), args, value = TRUE)
  if (!length(hit)) default else sub(paste0("^--", nm, "="), "", hit[1])
}
n_cores <- as.integer(opt("cores", "10"))
force <- flag("force")
cells_opt <- opt("cells", "")
pilot <- flag("pilot")       # A4 transplant: guohe_sec52_run.R:70
probe <- flag("probe")       # Stage 1 identity probe (A3), prints only

# Named-line change 1: t7 only (v1 section 4).
MFC_CELLS <- sprintf("t7_beta2_%02d", 0:5)
if (nzchar(cells_opt)) MFC_CELLS <- strsplit(cells_opt, ",")[[1]]

# A4 transplant: guohe_sec52_run.R:81-84 -- pilot is the maximal-bias cell at
# reps = 20, output suffixed _pilot, projection printed, production not run.
n_rep_opt <- NA_integer_
if (pilot) {
  MFC_CELLS <- "t7_beta2_00"
  n_rep_opt <- as.integer(opt("reps", "20"))
} else if (nzchar(opt("reps", ""))) {
  n_rep_opt <- as.integer(opt("reps", ""))
}

# Addendum-A columns joined from the 2026-09-04 bundles (E5), unchanged
MF_JOIN_COLS <- c("sel_agree_mr", "p_hat_H", "p_top1", "p_top2", "p_top3",
                  "p_lab1", "p_lab2", "p_lab3", "Sigma_HH", "A6_mass",
                  "A6_mass_std", "m0_hat", "m0_mc_se", "M_eff",
                  "tie_resid_implied")

.mfc_gh_cols <- function(row, theta) {
  out <- list()
  for (i in 1:4) {
    out[[sprintf("gh_r%d_deb", i)]] <- row[[sprintf("r%d_bias", i)]] + theta
    out[[sprintf("gh_r%d_low", i)]] <- theta - row[[sprintf("r%d_dist", i)]]
    out[[sprintf("gh_r%d_cov", i)]] <- row[[sprintf("r%d_cover", i)]]
  }
  out
}

# ---- the gate call (A3) ----------------------------------------------------
# The argument list below is COPIED VERBATIM from mv_mr()'s executing block,
# quarto/GuoHe/mr_vs_guohe_sim.R:132-151, with exactly four arguments added
# (marked A3). mv_mr() cannot be used directly: it has no
# field_scale_complement formal, and modifying it is forbidden by A3.
mfc_mr <- function(df, cands, sel_label, spec, draws = MV_DRAWS,
                   multiplier = MV_MULTIPLIER, seed = NULL,
                   ci_method = "field", field_R_out = MF_R_OUT,
                   field_R_in = MF_R_IN,
                   field_uniform = FALSE, ij_residual = "two_term") {
  args <- list(
    df = df, candidates = cands, spec = spec,
    selected_members = cands[[sel_label]],
    admission = list(effect_floor = NULL, consistency = NULL),
    reselection = "maxeff",
    draws = draws, multiplier = multiplier,
    ci_method = ci_method, seed = seed, return_reselection = TRUE,
    include_complement = TRUE,                    # A3 addition (was FALSE)
    # winner-only IJ variants (TASK_complement_refinements_2026-09-06);
    # "two_term" keeps the reported interval byte-identical.
    ij_residual = ij_residual)
  if (identical(ci_method, "field"))
    args <- c(args, list(field_R_out = field_R_out, field_R_in = field_R_in,
                         # kappa(Sigma-hat) sweep (TASK_mr_field_uniform_2026-09-05);
                         # FALSE keeps the 2026-09-05 field output byte-identical.
                         field_uniform = field_uniform,
                         # complement field (TASK_mr_field_complement_2026-09-06);
                         field_complement = TRUE,              # A3 addition
                         field_scale_complement = "selected")) # A3 addition
  do.call(forestsearch:::fs_mr_inference, args)
}

.nz <- function(x) if (is.null(x) || !length(x)) NA_real_ else as.numeric(x)[1]
.lg <- function(x) { v <- .nz(x); if (is.na(v)) NA_real_ else log(v) }

.mfc_mr_cols <- function(mr, row_c, theta, t_mr) {
  est <- log(mr$debiased$est)
  lo1 <- log(mr$debiased$lower_1s)
  # pairing proof vs the 2026-09-04 row (same seed -> identical values)
  cur_ok <- identical(est, row_c$mr_est) &&
    identical(mr$debiased$se_ij, row_c$mr_se_ij) &&
    identical(mr$selection_bias, row_c$mr_bias_sel) &&
    identical(mr$fixed_bias, row_c$mr_bias_fix) &&
    identical(unname(mr$reselection$p_hat[mr$selected_index]), row_c$p_hat_H)
  f  <- mr$field
  fc <- f$complement
  jt <- f$joint
  js <- f$joint_s
  gc <- mr$complement
  e2 <- log(f$est2); flo1 <- log(f$lower_1s)
  flo2 <- log(f$lower_2s); fup2 <- log(f$upper_2s)

  # --- the two 0.975 one-sided bounds and the both-correct indicator (A2) ---
  # Engine's own Bonferroni pair at gamma = alpha/2 = 0.025. The harm side is
  # identical in `joint` and `joint_s` (only the complement draws differ), so
  # bonf_lower_H is quoted once.
  h_lo_975  <- .lg(jt$bonf_lower_H)
  c_up_975  <- .lg(jt$bonf_upper_Hc)
  c_up_975s <- .lg(js$bonf_upper_Hc)
  h_lo_975s <- .lg(js$bonf_lower_H)
  # 0.95 one-sided: harm lower (field) and complement upper (field / field-s).
  h_lo_95   <- flo1
  c_up_95   <- .lg(fc$upper_1s)
  c_up_95s  <- .lg(fc$upper_1s_s)

  p_hat_H_now <- unname(mr$reselection$p_hat[mr$selected_index])

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

    # ---- NEW: harm side at both levels -------------------------------------
    h_lower_95 = h_lo_95,
    h_lower_975 = h_lo_975,
    h_cover_95 = as.integer(is.finite(h_lo_95) && h_lo_95 <= theta),
    h_cover_975 = as.integer(is.finite(h_lo_975) && h_lo_975 <= theta),

    # ---- NEW: complement upper at both levels, unscaled and field-s --------
    # Truth for the complement is exactly 0 (see header), so coverage is
    # upper >= 0.
    c_upper_95 = c_up_95, c_upper_975 = c_up_975,
    c_upper_95_s = c_up_95s, c_upper_975_s = c_up_975s,
    c_cover_95 = as.integer(is.finite(c_up_95) && c_up_95 >= 0),
    c_cover_975 = as.integer(is.finite(c_up_975) && c_up_975 >= 0),
    c_cover_95_s = as.integer(is.finite(c_up_95s) && c_up_95s >= 0),
    c_cover_975_s = as.integer(is.finite(c_up_975s) && c_up_975s >= 0),

    # ---- NEW: joint both-correct indicators (Bonferroni pair) --------------
    joint_cover = as.integer(is.finite(h_lo_975) && is.finite(c_up_975) &&
                               h_lo_975 <= theta && c_up_975 >= 0),
    joint_cover_s = as.integer(is.finite(h_lo_975s) && is.finite(c_up_975s) &&
                                 h_lo_975s <= theta && c_up_975s >= 0),
    joint_gamma = .nz(jt$gamma), joint_gamma_s = .nz(js$gamma),
    joint_prob = .nz(jt$joint_prob), joint_prob_s = .nz(js$joint_prob),
    joint_bonf_prob = .nz(jt$bonf_joint_prob),
    joint_bonf_prob_s = .nz(js$bonf_joint_prob),
    joint_corr = .nz(jt$corr), joint_corr_s = .nz(js$corr),
    joint_n_draws = .nz(jt$n_joint_draws),

    # ---- NEW: complement diagnostics ---------------------------------------
    cfld_lambda_mean = .nz(fc$lambda_mean), cfld_lambda_sd = .nz(fc$lambda_sd),
    cfld_se_field = .nz(fc$se_field), cfld_se_field_s = .nz(fc$se_field_s),
    cfld_lambda_mean_s = .nz(fc$lambda_mean_s),
    cfld_est2 = .lg(fc$est2), cfld_est2_s = .lg(fc$est2_s),
    cfld_lower_1s = .lg(fc$lower_1s), cfld_lower_1s_s = .lg(fc$lower_1s_s),
    cfld_lower_2s = .lg(fc$lower_2s), cfld_upper_2s = .lg(fc$upper_2s),
    cfld_lower_2s_s = .lg(fc$lower_2s_s), cfld_upper_2s_s = .lg(fc$upper_2s_s),
    cfld_n_in_used_mean = .nz(fc$n_in_used_mean),
    cgate_naive_est = .lg(gc$naive$est),
    cgate_est = .lg(gc$debiased$est),
    cgate_se_ij = .nz(gc$debiased$se_ij),
    cgate_se_wald = .nz(gc$debiased$se_wald),
    cgate_bias_sel = .nz(gc$selection_bias),
    cgate_bias_fix = .nz(gc$fixed_bias),
    cgate_n = .nz(gc$n),
    c_na = as.integer(!(is.finite(c_up_95) && is.finite(c_up_975) &&
                          is.finite(c_up_95s) && is.finite(c_up_975s))),

    # ---- NEW: p-hat(H-hat) recomputed in THIS run --------------------------
    p_hat_H_now = p_hat_H_now,

    t_field_s = f$timing_seconds, t_mr_s = t_mr)
}

mfc_rep_51 <- function(id, m, row_r, row_c, beta, n, base) {
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
  mr <- mfc_mr(df, mv_cand_idx_51(df, k), paste0("S", nv$sel), mv_spec51,
               seed = base + m + MV_SEED_MR)
  t1 <- proc.time()[["elapsed"]]
  as.data.frame(c(list(
    id = id, m = m, seed_data = base + m, seed_mr = base + m + MV_SEED_MR,
    k_family = k, theta = theta, sel = nv$sel,
    naive_ok = as.integer(naive_ok),
    naive_est = nv$point, naive_se = fits$se[nv$sel],
    naive_lower = nv$lower, naive_cover = nv$cover),
    .mfc_gh_cols(row_r, theta),
    .mfc_mr_cols(mr, row_c, theta, t1 - t0)),
    stringsAsFactors = FALSE)
}

mfc_rep_52 <- function(id, m, row_r, row_c, base, tru) {
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
  # N3 (v4): the top-1 minus top-2 oriented-score gap at selection. orient = +1
  # on t7, so the oriented score is fits$est and the selection is which.max().
  # Recorded per replicate so any selection mismatch can be reported with the
  # margin it crossed. Measured floor over 106 selections is 1.438e-04 against
  # a worst float deviation of 4.44e-16 -- a ratio of 3.2e11.
  .e <- fits$est[is.finite(fits$est)]
  sel_gap <- if (length(.e) >= 2L) { .s <- sort(.e, decreasing = TRUE); .s[1] - .s[2] } else NA_real_
  theta <- row_r$gamma_s
  t0 <- proc.time()[["elapsed"]]
  mr <- mfc_mr(cand$df, mv_cand_idx_52(cand), cand$names[sel], mv_spec52,
               seed = base + m + MV_SEED_MR)
  t1 <- proc.time()[["elapsed"]]
  as.data.frame(c(list(
    id = id, m = m, seed_data = base + m, seed_mr = base + m + MV_SEED_MR,
    k_family = length(cand$cuts), theta = theta, sel = sel,
    c_hat = cand$cuts[sel], n_sel = row_r$n_sel, sel_gap = sel_gap,
    naive_ok = as.integer(naive_ok && sel_ok),
    naive_est = nv$point, naive_se = (nv$point - nv$lower) / stats::qnorm(0.95),
    naive_lower = nv$lower, naive_cover = nv$cover,
    gamma_s_naive = nv$gamma_s),
    .mfc_gh_cols(row_r, theta),
    .mfc_mr_cols(mr, row_c, theta, t1 - t0)),
    stringsAsFactors = FALSE)
}

# ---- Stage 1 identity probe (A3) -------------------------------------------
# The hand-assembled fs_mr_inference() call must reproduce the stored path
# exactly. Columns compared against mr_field_vs_guohe_<id>.rds by identical().
MFC_PROBE_NAIVE <- c("naive_est", "naive_se", "naive_lower", "naive_cover")
MFC_PROBE_IJ    <- c("mr_est", "mr_se_ij", "mr_se_wald", "mr_lower_1s",
                     "mr_lower_2s", "mr_upper_2s", "mr_cover",
                     "mr_bias_sel", "mr_bias_fix", "mr_selection_rate",
                     "mr_mean_r")
MFC_PROBE_FIELD <- c("fld_lambda_mean", "fld_lambda_sd", "fld_q05", "fld_q25",
                     "fld_q50", "fld_q75", "fld_q95", "fld_q025", "fld_q975",
                     "fld_n_out_used", "fld_n_in_used_mean",
                     "fld_est2", "fld_lower_1s", "fld_lower_2s", "fld_upper_2s",
                     "fld_lower_se", "fld_upper_se", "fld_cover_1s",
                     "fld_cover_2s")

# ---- N1 pairing standard (v4; replaces A3's plain identical()) -------------
# Integer / selection / flag / seed columns: identical() REQUIRED, every
# replicate; any failure is a STOP. Floating-point columns: all.equal() at
# 1e-8, with the worst absolute and relative deviation and the bit-identical
# fraction recorded per cell. The cross-machine comparison is a provenance
# measurement, not a gate: stored bundles are x86_64 / R 4.6.1 / reference
# BLAS, this Mac is arm64 / R 4.5.2 / Accelerate.
MFC_DISCRETE <- c("sel", "c_hat", "n_sel", "seed_data", "seed_mr", "m",
                  "k_family", "naive_cover", "mr_cover", "fld_cover_1s",
                  "fld_cover_2s", "fld_n_out_used", "mr_ij_source",
                  "mr_ij_draws",
                  paste0("gh_r", 1:4, "_cov"))
# RECLASSIFIED (Larry, 2026-09-09): `theta` and `gamma_s_naive` moved from the
# discrete class to the float class. N1's discrete class exists to prove the
# SELECTION is unchanged, not to constrain computed continuous values. The
# selection keys -- c_hat, c_hat_naive, sel, n_sel and the seeds -- are verified
# identical at 24,000/24,000. The truth lookups are NOT labels: gh52_truth_at()
# is stats::approx(truth$c_grid, y, xout = c_hat, rule = 2)$y
# (guohe_sec52_truth.R:314), i.e. linear interpolation, so its return is
# computed arithmetic and belongs to the float class. Measured residual across
# all six cells: 24/24,000 values differ (0.10%), worst |diff| 5.55e-17,
# all.equal at 1e-8 TRUE throughout. At beta2 = 0 the truth curve is identically
# zero, so interpolation returns exact 0 and that cell shows 0/2000 -- which is
# why the six-replicate probe never surfaced this.
MFC_FLOAT <- c(setdiff(MFC_PROBE_NAIVE, MFC_DISCRETE),
               setdiff(MFC_PROBE_IJ, MFC_DISCRETE),
               setdiff(MFC_PROBE_FIELD, MFC_DISCRETE),
               "theta", "gamma_s_naive")
MFC_TOL <- 1e-8

# Compare one recomputed cell against the stored bundle under N1.
# Returns the discrete failures (STOP-worthy) and the float deviation summary.
.mfc_pair_n1 <- function(res, st) {
  disc_bad <- character(0)
  for (nm in intersect(MFC_DISCRETE, intersect(names(res), names(st)))) {
    a <- res[[nm]]; b <- st[[nm]]
    if (!identical(a, b)) {
      first <- which(!mapply(identical, as.list(a), as.list(b)))[1]
      disc_bad <- c(disc_bad, sprintf("%s (first differing replicate m = %s)",
                                      nm, res$m[first]))
    }
  }
  worst_abs <- 0; worst_abs_at <- ""; worst_rel <- 0; worst_rel_at <- ""
  n_val <- 0L; n_id <- 0L; ae_bad <- character(0)
  for (nm in intersect(MFC_FLOAT, intersect(names(res), names(st)))) {
    a <- as.numeric(res[[nm]]); b <- as.numeric(st[[nm]])
    ok <- is.finite(a) & is.finite(b)
    if (!any(ok)) next
    n_val <- n_val + sum(ok); n_id <- n_id + sum(a[ok] == b[ok] &
                                                  mapply(identical, as.list(a[ok]), as.list(b[ok])))
    ad <- abs(a[ok] - b[ok])
    # near-zero values scored on absolute difference only (mr_mean_r is ~0 by
    # construction, so a relative figure there is meaningless)
    rl <- ifelse(abs(b[ok]) > 1e-10, ad / abs(b[ok]), ad)
    if (max(ad) > worst_abs) {
      worst_abs <- max(ad)
      worst_abs_at <- sprintf("%s / m=%s", nm, res$m[ok][which.max(ad)])
    }
    if (max(rl) > worst_rel) {
      worst_rel <- max(rl)
      worst_rel_at <- sprintf("%s / m=%s", nm, res$m[ok][which.max(rl)])
    }
    if (!isTRUE(all.equal(a[ok], b[ok], tolerance = MFC_TOL)))
      ae_bad <- c(ae_bad, nm)
  }
  list(disc_bad = disc_bad, ae_bad = ae_bad,
       worst_abs = worst_abs, worst_abs_at = worst_abs_at,
       worst_rel = worst_rel, worst_rel_at = worst_rel_at,
       n_val = n_val, n_identical = n_id,
       frac_identical = if (n_val) n_id / n_val else NA_real_)
}

mfc_probe_cell <- function(id, n_probe) {
  f_old <- file.path(.mfc_dir, paste0("mr_field_vs_guohe_", id, ".rds"))
  if (!file.exists(f_old)) stop("stored bundle absent: ", f_old)
  old <- readRDS(f_old)$results
  rep_bun <- readRDS(file.path(.mfc_dir, paste0("guohe_repro_", id, ".rds")))
  cmp_bun <- readRDS(file.path(.mfc_dir, paste0("mr_vs_guohe_", id, ".rds")))
  base <- mv_gh52_base(id)
  b2 <- as.integer(sub("^t7_beta2_", "", id))
  tru <- readRDS(file.path(.mfc_dir, sprintf("guohe_sec52_truth_beta2_%02d.rds", b2)))
  cat(sprintf("\n[probe] %s  %d replicates\n", id, n_probe))
  rows <- lapply(seq_len(n_probe), function(m)
    mfc_rep_52(id, m, rep_bun$results[m, ], cmp_bun$results[m, ], base, tru))
  res <- do.call(rbind, rows)
  stopifnot(identical(old$seed_data[seq_len(n_probe)], res$seed_data))
  groups <- list(naive = MFC_PROBE_NAIVE, IJ = MFC_PROBE_IJ, field = MFC_PROBE_FIELD)
  bad <- character(0)
  for (g in names(groups)) {
    for (nm in groups[[g]]) {
      a <- res[[nm]]; b <- old[[nm]][seq_len(n_probe)]
      if (!identical(a, b)) {
        first <- which(!mapply(identical, as.list(a), as.list(b)))[1]
        bad <- c(bad, sprintf("%s/%s (first differing replicate m = %s; new %.17g vs stored %.17g)",
                              g, nm, first, as.numeric(a[first]), as.numeric(b[first])))
      }
    }
  }
  cat(sprintf("  naive columns : %s\n", if (!length(grep("^naive/", bad))) "IDENTICAL" else "MISMATCH"))
  cat(sprintf("  IJ columns    : %s\n", if (!length(grep("^IJ/", bad))) "IDENTICAL" else "MISMATCH"))
  cat(sprintf("  field columns : %s\n", if (!length(grep("^field/", bad))) "IDENTICAL" else "MISMATCH"))
  cat(sprintf("  naive_ok      : %d/%d   cur_ok: %d/%d\n",
              sum(res$naive_ok), nrow(res), sum(res$cur_ok), nrow(res)))
  cat(sprintf("  new columns   : c_upper_95 %s | c_upper_975 %s | joint_cover %s\n",
              paste(sprintf("%.4f", res$c_upper_95), collapse = ", "),
              paste(sprintf("%.4f", res$c_upper_975), collapse = ", "),
              paste(res$joint_cover, collapse = ", ")))
  if (length(bad)) {
    cat("  MISMATCHES:\n"); for (b in bad) cat("    ", b, "\n")
  }
  invisible(list(id = id, n = n_probe, bad = bad, res = res))
}

mfc_run_cell <- function(id) {
  f_out <- file.path(.mfc_dir, paste0("mr_field_complement_vs_guohe_", id,
                                      if (pilot) "_pilot" else "", ".rds"))
  if (file.exists(f_out) && !force) {
    cat(sprintf("[skip] %s (exists)\n", id))
    return(invisible(NULL))
  }
  sec52 <- grepl("^t7_", id)
  rep_bun <- readRDS(file.path(.mfc_dir, paste0("guohe_repro_", id, ".rds")))
  cmp_bun <- readRDS(file.path(.mfc_dir, paste0("mr_vs_guohe_", id, ".rds")))
  base <- if (sec52) mv_gh52_base(id) else mv_gh51_base(id)
  stopifnot(rep_bun$seed_base == base, cmp_bun$seed_base == base,
            nrow(cmp_bun$results) == nrow(rep_bun$results),
            identical(cmp_bun$results$m, seq_len(nrow(cmp_bun$results))))
  n_rep <- nrow(rep_bun$results)
  if (!is.na(n_rep_opt)) n_rep <- min(n_rep, n_rep_opt)
  tru <- NULL; sc <- NULL
  if (sec52) {
    b2 <- as.integer(sub("^t7_beta2_", "", id))
    tru <- readRDS(file.path(.mfc_dir, sprintf("guohe_sec52_truth_beta2_%02d.rds", b2)))
  } else sc <- mv_gh51_scenario(id)

  t0 <- proc.time()[["elapsed"]]
  cat(sprintf("[run ] %s  %d reps ...\n", id, n_rep)); utils::flush.console()
  rows <- parallel::mclapply(
    seq_len(n_rep),
    function(m) {
      try(if (sec52) mfc_rep_52(id, m, rep_bun$results[m, ], cmp_bun$results[m, ], base, tru)
          else mfc_rep_51(id, m, rep_bun$results[m, ], cmp_bun$results[m, ],
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
  # Production pairing check against the stored field bundle under the N1
  # standard (v4): discrete columns identical() -- any failure is a STOP;
  # floats at all.equal 1e-8, with the deviation summary recorded per cell.
  f_old <- file.path(.mfc_dir, paste0("mr_field_vs_guohe_", id, ".rds"))
  pair <- NULL; sel_mismatch <- NA_integer_; sel_mismatch_detail <- NULL
  if (file.exists(f_old)) {
    st <- readRDS(f_old)$results[res$m, ]
    pair <- .mfc_pair_n1(res, st)
    # N3: per-replicate selected-cutpoint tally against the stored row.
    smm <- which(res$c_hat != st$c_hat)
    sel_mismatch <- length(smm)
    if (sel_mismatch)
      sel_mismatch_detail <- data.frame(
        id = id, m = res$m[smm], c_hat_new = res$c_hat[smm],
        c_hat_stored = st$c_hat[smm], sel_gap = res$sel_gap[smm])
  }
  el <- proc.time()[["elapsed"]] - t0
  gate2 <- list(n_rep_expected = n_rep, n_rep_done = nrow(res),
                n_errored = sum(bad),
                naive_mismatch = sum(res$naive_ok == 0L),
                cur_mismatch = sum(res$cur_ok == 0L),
                # N1
                discrete_mismatch = if (is.null(pair)) NA_integer_ else length(pair$disc_bad),
                discrete_detail = if (is.null(pair)) NULL else pair$disc_bad,
                float_allequal_fail = if (is.null(pair)) NA_integer_ else length(pair$ae_bad),
                worst_abs = if (is.null(pair)) NA_real_ else pair$worst_abs,
                worst_abs_at = if (is.null(pair)) NA_character_ else pair$worst_abs_at,
                worst_rel = if (is.null(pair)) NA_real_ else pair$worst_rel,
                worst_rel_at = if (is.null(pair)) NA_character_ else pair$worst_rel_at,
                frac_identical = if (is.null(pair)) NA_real_ else pair$frac_identical,
                n_float_values = if (is.null(pair)) NA_integer_ else pair$n_val,
                # N3
                sel_mismatch = sel_mismatch,
                sel_mismatch_detail = sel_mismatch_detail,
                sel_gap_min = min(res$sel_gap, na.rm = TRUE),
                sel_gap_p01 = unname(quantile(res$sel_gap, 0.01, na.rm = TRUE, type = 7)),
                sel_gap_median = median(res$sel_gap, na.rm = TRUE),
                mr_na = sum(res$mr_na), fld_na = sum(res$fld_na),
                c_na = sum(res$c_na))
  # N1 / N3 STOP conditions, enforced at the cell boundary
  if (isTRUE(gate2$discrete_mismatch > 0L))
    stop("N1 STOP in ", id, ": discrete/selection/seed column(s) not identical(): ",
         paste(gate2$discrete_detail, collapse = "; "), call. = FALSE)
  if (isTRUE(gate2$float_allequal_fail > 0L))
    stop("N1 STOP in ", id, ": float column(s) outside all.equal tolerance ",
         MFC_TOL, ": ", paste(pair$ae_bad, collapse = ", "), call. = FALSE)
  if (isTRUE(gate2$sel_mismatch > 0L))
    stop("N3 STOP in ", id, ": ", gate2$sel_mismatch,
         " selected-cutpoint mismatch(es) vs the stored bundle. Detail: ",
         paste(sprintf("m=%s new=%.6f stored=%.6f gap=%.6g",
                       sel_mismatch_detail$m, sel_mismatch_detail$c_hat_new,
                       sel_mismatch_detail$c_hat_stored,
                       sel_mismatch_detail$sel_gap), collapse = "; "),
         call. = FALSE)
  saveRDS(list(
    id = id, section = if (sec52) "5.2" else "5.1",
    source_bundle = paste0("guohe_repro_", id, ".rds"),
    join_bundle = paste0("mr_vs_guohe_", id, ".rds"),
    pair_bundle = paste0("mr_field_vs_guohe_", id, ".rds"),
    join_cols = MF_JOIN_COLS,
    draws = MV_DRAWS, multiplier = MV_MULTIPLIER,
    field_R_out = MF_R_OUT, field_R_in = MF_R_IN, field_seed_offset = 900000L,
    field_complement = TRUE, include_complement = TRUE,
    field_scale_complement = "selected", return_reselection = TRUE,
    complement_truth = 0,
    pair_standard = list(rule = "v4 N1", discrete = "identical()",
                         float = "all.equal", tolerance = MFC_TOL,
                         discrete_cols = MFC_DISCRETE, float_cols = MFC_FLOAT),
    seed_base = base, mr_seed_offset = MV_SEED_MR, pilot = pilot,
    gate2 = gate2, elapsed_sec = el,
    sessionInfo = utils::capture.output(utils::sessionInfo()),
    results = res), f_out)
  cat(sprintf("[done] %s  %d/%d reps in %.1f min  naive_mm %d  cur_mm %d  disc_mm %s  sel_mm %s  worst_abs %.3g  frac_id %.3f  mr_na %d  fld_na %d  c_na %d  -> %s\n",
              id, gate2$n_rep_done, n_rep, el / 60, gate2$naive_mismatch,
              gate2$cur_mismatch, gate2$discrete_mismatch, gate2$sel_mismatch,
              gate2$worst_abs, gate2$frac_identical, gate2$mr_na,
              gate2$fld_na, gate2$c_na, basename(f_out)))
  utils::flush.console()
  invisible(res)
}

cat(sprintf("MR (field + complement) T1 run: %d cells, cores = %d, MR draws = %d (%s), field R_out/R_in = %d/%d%s\n\n",
            length(MFC_CELLS), n_cores, MV_DRAWS, MV_MULTIPLIER, MF_R_OUT, MF_R_IN,
            if (pilot) "  [PILOT]" else if (probe) "  [PROBE]" else ""))

if (probe) {
  n_probe <- if (!is.na(n_rep_opt)) n_rep_opt else 3L
  out <- lapply(MFC_CELLS, function(id) mfc_probe_cell(id, n_probe))
  nbad <- sum(vapply(out, function(z) length(z$bad), integer(1)))
  cat(sprintf("\n=== STAGE 1 PROBE: %s ===\n",
              if (nbad == 0L) "PASS (all compared columns identical)" else
                sprintf("FAIL (%d column(s) differ)", nbad)))
  quit(save = "no", status = if (nbad == 0L) 0L else 1L)
}

all_res <- list()
for (id in MFC_CELLS) all_res[[id]] <- mfc_run_cell(id)

# A4 transplant: guohe_sec52_run.R:214-225 -- projection, then STOP for the owner.
if (pilot) {
  r <- all_res[[MFC_CELLS[1]]]
  if (!is.null(r)) {
    per_rep <- mean(r$t_mr_s + r$t_field_s, na.rm = TRUE)
    total_core_h <- per_rep * 6 * 2000 / 3600
    cat("\n==== PILOT PROJECTION ====\n")
    cat(sprintf("  mean per-replicate (serial, gate + field): %8.2f s\n", per_rep))
    cat(sprintf("  full study 6 x 2000 reps, single core    : %8.1f core-h\n", total_core_h))
    cat(sprintf("  projected wall-clock at %3d cores        : %8.1f min\n",
                n_cores, total_core_h / n_cores * 60))
    cat(sprintf("  Gate 1a envelope: <= 40 core-h AND <= 90 min Mac wall -> %s\n",
                if (total_core_h <= 40 && total_core_h / n_cores * 60 <= 90)
                  "WITHIN" else "EXCEEDED"))
    cat("\nREPORT THIS PROJECTION TO THE OWNER BEFORE LAUNCHING PRODUCTION.\n")
  }
} else {
  cat("\n=== GATE 2 TALLY ===\n")
  for (id in MFC_CELLS) {
    f <- file.path(.mfc_dir, paste0("mr_field_complement_vs_guohe_", id, ".rds"))
    if (!file.exists(f)) { cat(sprintf("%-16s MISSING\n", id)); next }
    g <- readRDS(f)$gate2
    cat(sprintf("%-16s reps %d/%d  err %d  naive_mm %d  cur_mm %d  disc_mm %s  sel_mm %s  ae_fail %s  worst_abs %.3g  worst_rel %.3g  frac_id %.3f  mr_na %d  fld_na %d  c_na %d\n",
                id, g$n_rep_done, g$n_rep_expected, g$n_errored,
                g$naive_mismatch, g$cur_mismatch, g$discrete_mismatch,
                g$sel_mismatch, g$float_allequal_fail, g$worst_abs,
                g$worst_rel, g$frac_identical, g$mr_na, g$fld_na, g$c_na))
  }
}
