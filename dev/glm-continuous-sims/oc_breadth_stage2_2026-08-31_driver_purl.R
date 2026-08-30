# knitr::purl() of quarto/simulations/actg175/continuous/sim_fs_maxeffCons_mr_md120_knoise0_n500_batch_1_1000.qmd,
# truncated before the summary-prep chunk (setup-knobs, build-dgm, machinery, run-batch only): the
# driver executed without a render.  Run from the driver directory:
#   FS_STAGE2_RUN=comparability Rscript <this file>   (run 1, effect.threshold = 30)
#   FS_STAGE2_RUN=direct        Rscript <this file>   (run 2, effect.threshold = c1*)
## ----setup-knobs--------------------------------------------------------------
#| message: false
.t0_doc <- Sys.time()   # start of document wall-clock (total reported at the end)
# forestsearch must be INSTALLED (devtools::install()), not load_all(): the
# doFuture workers spawn separate R processes that only see the installed
# package.
library(forestsearch)
suppressPackageStartupMessages({
  library(doFuture); library(foreach); library(data.table)
})

# ── Parallelization topology ────────────────────────────────────────────────
#   "sims"  -> replicates fan out across workers (search + bootstrap sequential
#              within each). Best when n_sims >= n_workers.
#   "boots" -> replicates sequential; the search AND the nb_boots bootstrap fan
#              out across workers WITHIN each replicate. Best for FB runs.
# Bootstrap resample indices are seeded on the main process, so the two modes
# differ in wall-clock only, not in detection/estimates.
# ALIGNED TO THE SURVIVAL TEMPLATE (top-priority fix).  Two changes, both
# measured against the rendered batch_1_50 report:
#
#   1. parallel_mode "sims" -> "boots".  Under "sims" the inner plan is
#      SEQUENTIAL, so `forestsearch_bootstrap_dofuture()` ran the entire
#      nb_boots loop single-threaded inside each replicate -- the document's own
#      topology note above says "boots" is "Best for FB runs", but the knob said
#      otherwise.  Measured: FB 390 s at B=20 => 19.5 s PER BOOTSTRAP, serial.
#      The survival template runs "boots" and spreads its 300 bootstraps across
#      the pool (132 s per replicate, i.e. 0.44 s/bootstrap wall).
#
#   2. n_workers: detectCores() - 8 gave 120 on this 128-core host, but only
#      n_sims replicates ever receive work under "sims" -- so 100 workers sat
#      idle while multisession still spawned 120 R processes, each loading
#      forestsearch.  That is the k = 1.62 overhead the report printed (651 s
#      wall against ~400 s of modelled work).  The survival formula sizes the
#      pool off PHYSICAL cores at 90%, which is also what the package's own
#      .default_parallel_workers() does (80% of physical).
#
# Projected effect at production (B = 300, 50 replicates): "sims" leaves the
# bootstrap serial at ~1.6 h per replicate with most cores idle; "boots" spreads
# each replicate's 300 bootstraps across the pool, ~1 min per replicate.
parallel_mode <- "sims"
n_workers     <- ceiling(0.90 * max(1L, parallel::detectCores(logical = FALSE) - 1L))

# ── Smoke-test counts (the only knobs to change for production) ─────────────
n_sims       <- 1000L     # QUICK RUN for review        (production: 50L per batch)
diag_mode    <- FALSE   # TRUE caps to 5 reps and relays worker errors verbatim
nb_boots     <- NULL     # PRESENTATION REVIEW ONLY -- below any FB validity
                        # floor (manuscript uses 300); do not quote FB numbers
                        # from this render.            (production: 300L)
mr_draws     <- 5000L   # MR multiplier draws
k_random_noise <- 0     # add k standard-normal noise confounders (0 = none)

# ── Batching across renders (seed-disjoint pieces; pool via combine) ────────
# Every per-replicate seed is seed_base + sim_id (DGM, search, bootstrap, and
# noise), so disjoint sim_id ranges across renders pool to exactly one big run.
run_mode     <- "batch"   # "batch" (run + save one piece) | "combine" (merge)
sim_id_start <- 1L
seed_base    <- 8316951L

# ── The identifier: engine + focus (these ARE the estimator; they feed the stem)
subgroup_method <- "consistency"
sg_focus        <- "maxeffCons"  # highest effect meeting the consistency rate

# ── Scenario knob and sample size ───────────────────────────────────────────
target_md_harm <- -120    # STAGE 2 (cc_task_oc_breadth_stage2_2026-08-31): the
                          # forecast rung q = 120; the DGM is built DIRECTLY at
                          # the Stage 1 k_inter (below), not calibrated.
                          # (raw scale; oriented +40). NULL-ANALOGUE CELL (D3):
                          # set this to the calibrator's no-interaction effect.
n_sample       <- 500L    #                             (production: 500L)

# NULL / ALTERNATIVE CELL SWITCH (alignment item D3).  Defined HERE, not with
# the DGM, because the output stem below encodes the cell -- otherwise a null
# run would overwrite the alternative bundle of the same name.  See @sec-dgm for
# what it does.  Pair with target_hr_harm = 1.0 on the survival side.
null_cell      <- FALSE

# STAGE 2 -- the locked forecast (commit c9cb0ca2) supplies k_inter and c1*.
# Run 1 (comparability): effect.threshold = 30, the driver's own.
# Run 2 (direct):        FS_STAGE2_RUN=direct -> effect.threshold = c1* at full
#                        precision and the output stem gains "_c1star".
# Everything else is identical between the two runs, including the seed table.
stage2_forecast <- readRDS(file.path("..", "..", "..", "..", "dev", "glm-continuous-sims",
                                     "oc_breadth_ladder_2026-08-30_forecast120.rds"))
stopifnot(stage2_forecast$q == 120, abs(stage2_forecast$k_inter + 93.7447641240) < 1e-9)
stage2_direct <- identical(Sys.getenv("FS_STAGE2_RUN"), "direct")

method_tag <- if (identical(subgroup_method, "consistency")) "fs" else subgroup_method
# One source of truth: fs_focus_tag() holds the (subgroup_method, sg_focus)
# -> stem tag map for every driver.  It was pasted into each one and drifted;
# the method-blind copies tagged DINA/GRF runs with the consistency name for a
# rule those engines never ran.
focus_tag <- forestsearch::fs_focus_tag(subgroup_method, sg_focus)
if (!identical(focus_tag, sg_focus))
  cat(sprintf("NOTE: sg_focus '%s' is an alias on the %s path; output stem uses focus tag '%s'.
",
              sg_focus, subgroup_method, focus_tag))
rds_stem <- sprintf("%s_%s_mr_md%02d_knoise%d_n%d",
                    method_tag, focus_tag, abs(target_md_harm),
                    k_random_noise, n_sample)
# Quick-run renders are stem-tagged so their bundle can NEVER pool with real
# batches: combine globs match the untagged stem, and the tag also renames the
# results directory.  Set FALSE for production batches.
quickrun <- FALSE
if (isTRUE(null_cell)) rds_stem <- sub("_md[0-9]+_", "_mdnull_", rds_stem)
if (isTRUE(quickrun)) rds_stem <- paste0(rds_stem, "_quickrun")
if (isTRUE(stage2_direct)) rds_stem <- paste0(rds_stem, "_c1star")
results_dir <- file.path("mr_md_harm", sprintf("%s_s%d_d%d",
                         rds_stem, n_sims, mr_draws))
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# ── Search configuration (committed MD values; resample identification) ─────
consistency_method <- "resample"
# MR's de-biasing rule is derived by forestsearch() from sg_focus — not set
# here. Override only via mr_inference_args(reselection=) (advanced use).
selection_rule      <- "neighborhood"  # band rule for effMaxSG/effMinSG arms;
effect_neighborhood <- 0.10            # inert under maxcons-family foci; pinned
stop_threshold      <- NULL            # pinned so the resolved setting is visible
md_threshold   <- if (isTRUE(stage2_direct)) stage2_forecast$c1_star else 30
                        # effect.threshold, oriented scale: 30 (run 1, the
                        # driver's own) or c1* = 135.7411624608 read from the
                        # locked forecast (run 2, the direct run).
md_consistency <- 10    # consistency.threshold. OPEN QUESTION D2: sits below
                        # both CATE values, hence non-discriminating.
pconsistency   <- 0.90
fs_splits <- 400L; maxk <- 2L; n_min <- 60L; d0_min <- 12L; d1_min <- 12L

# ── Candidate cuts on the true-region covariates (J-quantile grids) ──────────
# THIS WAS MISSING, and the omission mattered more than it looks.  With no cut
# arguments supplied at all, forestsearch() used its DEFAULT cut set -- q25,
# median, q75 per continuous covariate -- and on this DGM those defaults contain
# BOTH true boundaries EXACTLY:
#
#   age     default cuts 29.0 / 34.0 / 40.0    <- the median IS the true cut 34
#   preanti default cuts  0.0 / 136.0 / 744.5  <- q75 IS the true cut 744.5
#
# That is not coincidence: actg_age_cut and actg_preanti_cut were themselves set
# to those quantiles.  The search was therefore being handed the exact answer,
# and could reach Jaccard 1.000 against the true region.
#
# The survival template deliberately avoids this.  Its true z1 cut is er's q25
# (= 8), which is also a default cut, and `conf.cont_jcuts = list(er = 10)`
# REPLACES er's defaults with a J-quantile grid (0, 3, 9, 17, 30, ...) that does
# NOT contain 8.  The search must approximate the boundary rather than be given
# it.  This restores that design symmetry on the continuous side.
#
# Cuts are placed at the k/(J+1) empirical quantiles, k = 1..J
# (get_FSdata_helpers.R cut_var_jq), and jcuts REPLACE the default set for the
# named variables (get_fsdata.R: names are stripped from conf.cont_medians).
# Measured consequence on this DGM:
#
#   default cuts : best achievable Jaccard 1.000 (age > 34.0 & preanti <= 744.5)
#   J = 10       : best achievable Jaccard 0.882 (age > 33.0 & preanti <= 677.9)
#
# CAVEAT on preanti.  Over 40% of patients are antiretroviral-naive
# (preanti = 0), so the lower half of its J-grid collapses:
#   0, 0, 0, 0, 55.8, 214.1, 432, 677.9, 872.8, 1070.5
# Four of the ten cuts are the same "preanti <= 0" condition, so the effective
# resolution is ~7 distinct cuts, concentrated above the median.  Raising J for
# preanti buys resolution only in the upper tail.
#
# Set to NULL to revert to default cuts (and to the exactly-representable truth).
fs_conf.cont_jcuts <- list(age = 10, preanti = 10)
use_lasso <- FALSE; use_dina <- FALSE; use_grf <- FALSE; use_twostage <- TRUE
is_rct <- TRUE; vi_grf_min <- -0.2
adverse_outcome <- FALSE   # higher CD4 change = better; harm = negative MD

# ── DGM knobs (NOW transplanted verbatim from the committed harness) ────────
# mr_coverage_sweep_md_harm.qmd:102-114.  Do NOT recalibrate: cal_target_md =
# -40 converges at machine precision with exactly these values.
#
# The as-received document had the SURVIVAL template's DGM parameterisation
# here, which cannot calibrate on this DGM -- see the commit message.  The
# harm region is pre-binarised into z1/z2 and the cuts are plain factor
# levels; `type = "less_eq"` is not a cutpoint type the package implements
# (generate_aft_dgm_helpers.R:750-785 accepts "greater", "multiple", "custom"
# only), and the k_inter search interval is [0, 120], not [-3, 0].
actg_arms        <- c(1L, 3L)
actg_treat_arm   <- 1L
actg_age_cut     <- 34
actg_preanti_cut <- 744.5
dgm_n_super <- 5000L
cal_k_grid_range <- c(0, 120); cal_grid_step <- 2
dgm_factor_vars   <- paste0("z", 1:12)
dgm_subgroup_vars <- c("z1", "z2")
dgm_subgroup_cuts <- list(z1 = 1L, z2 = 1L)
# The ANALYSIS covariates stay the raw ACTG175 variables -- they are NOT the
# DGM's z-factors.  age and preanti enter continuous, so the search has to
# recover the z1/z2 region through its own quartile cuts; that is the design,
# and tying analysis_binary_vars to dgm_factor_vars (as received) would have
# handed the search the answer pre-binarised.
analysis_continuous_vars <- c("age","preanti","wtkg","karnof","cd40","cd80")
analysis_binary_vars     <- c("hemo","homo","drugs","race","gender","symptom")

# ── str2: the literature's proxy for the true partner (OFF by default) ───────
# `preanti` (days of prior ART) and `str2` (naive vs experienced) are the SAME
# construct at two resolutions, and the ACTG175 HTE literature turns on exactly
# that distinction: the standard 12-covariate Lu et al. (2013) candidate set
# carries the BINARY str2 and not the continuous preanti, which is why prior ITR
# methods could not recover cuts at intermediate prior-therapy exposure -- the
# gap this design exploits.
#
# Measured on the analysis data (arms 1/3, n = 1083):
#   cor(str2, preanti)                  =  0.678
#   cor(str2, 1{preanti <= 744.5})      = -0.488
#   every str2 == 0 patient (451/451) satisfies preanti <= 744.5
# so {str2 == 0} is a strict SUBSET of the true partner condition -- a genuine
# coarsened substitute, and the direct analogue of the age-for-meno proxy in the
# GBSG cell.  Nothing else in the pool comes close (next best |r| = 0.15).
#
# CHANGES THE SEARCH.  Adding str2 enlarges the candidate family, so detection
# rates and realized rules WILL move relative to the previously committed runs;
# this is not a presentation-only knob.  Turned ON deliberately: with it FALSE
# str2 cannot be selected and the identification proxy row reads 0% by
# construction rather than as a finding.  Set FALSE to reproduce the older
# configuration exactly.
include_str2 <- TRUE
if (isTRUE(include_str2))
  analysis_binary_vars <- c(analysis_binary_vars, "str2")
outcome_name <- "y_sim"; treat_name <- "treat_sim"
id_name <- "id"; harm_col <- "flag_harm"
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || all(is.na(a))) b else a
run_fb <- !is.null(nb_boots) && nb_boots >= 1L
cat(sprintf("Batch %d..%d | n=%d | FB=%s(nb=%s) | MR draws=%d | mode=%s x %d workers\n",
            sim_id_start, sim_id_start + n_sims - 1L, n_sample,
            run_fb, nb_boots %||% 0L, mr_draws, parallel_mode, n_workers))


## ----build-dgm----------------------------------------------------------------
actg_df <- subset(speff2trial::ACTG175, arms %in% actg_arms)
actg_df$id <- seq_len(nrow(actg_df))
# TREATMENT IS SWITCHED (ddI = 1), verbatim from the committed harness
# (mr_coverage_sweep_md_harm.qmd:157-158).  The -40 calibration is defined
# against THIS coding; the as-received `treat <- (arms == 1)` is its negation
# and would have calibrated a different cell under the same name.
actg_df$treat_orig <- ifelse(actg_df$arms == actg_treat_arm, 1L, 0L)
actg_df$treat      <- 1L - actg_df$treat_orig
actg_df$cd4_change <- actg_df$cd420 - actg_df$cd40
actg_df <- actg_df[!is.na(actg_df$cd420), ]
actg_df$z1  <- as.factor(ifelse(actg_df$age > actg_age_cut, 1L, 0L))
actg_df$z2  <- as.factor(ifelse(actg_df$preanti <= actg_preanti_cut, 1L, 0L))
actg_df$z3  <- as.factor(ifelse(actg_df$wtkg <= 75, 1L, 0L))
actg_df$z4  <- as.factor(ifelse(actg_df$karnof <= median(actg_df$karnof), 1L, 0L))
actg_df$z5  <- as.factor(ifelse(actg_df$cd40 <= median(actg_df$cd40), 1L, 0L))
actg_df$z6  <- as.factor(ifelse(actg_df$cd80 <= median(actg_df$cd80), 1L, 0L))
actg_df$z7  <- as.factor(actg_df$hemo);  actg_df$z8  <- as.factor(actg_df$homo)
actg_df$z9  <- as.factor(actg_df$drugs); actg_df$z10 <- as.factor(actg_df$race)
actg_df$z11 <- as.factor(actg_df$gender);actg_df$z12 <- as.factor(actg_df$symptom)
for (v in analysis_binary_vars) actg_df[[v]] <- as.factor(actg_df[[v]])
confounders_analysis <- c(analysis_continuous_vars, analysis_binary_vars)

# NULL / ALTERNATIVE CELL PAIRING (alignment item D3).
# The survival template's cell is selected by `target_hr_harm` -- 1.0 is its
# NULL, 1.5/2.0 its alternatives.  This document had no equivalent switch: it
# could only build the alternative, so a survival null was being compared
# against a continuous alternative.  That mismatch is not cosmetic; detection
# was 100% here against 60% there, and because the FB bootstrap skips the whole
# re-selection block on a resample that finds nothing
# (bootstrap_analysis_dofuture.R:626), it also drove much of the per-bootstrap
# cost gap.
#
# `model = "null"` forces beta_inter <- 0 inside generate_glm_dgm()
# (generate_glm_dgm.R:347), so there is no treatment-by-subgroup interaction and
# effect_Q == effect_Qc == effect_ITT.  No calibration is needed or possible --
# there is no target to hit -- so the calibrator is bypassed rather than asked
# for a target of zero.
# (`null_cell` is set in the setup chunk, since the output stem depends on it.)
dgm <- if (isTRUE(null_cell)) {
  generate_glm_dgm(
    data = actg_df, factor_vars = dgm_factor_vars,
    outcome_var = "cd4_change", treatment_var = "treat",
    outcome_type = "continuous", effect_measure = "MD",
    subgroup_vars = dgm_subgroup_vars, subgroup_cuts = dgm_subgroup_cuts,
    model = "null", n_super = dgm_n_super, seed = seed_base, verbose = FALSE)
} else {
  # STAGE 2: the exact linear route Stage 1 used for the rung (build_direct in
  # oc_breadth_ladder_2026-08-30.R) -- every argument of the MD40 build's direct
  # form, with k_inter read from the locked forecast.
  generate_glm_dgm(
    data = actg_df, factor_vars = dgm_factor_vars,
    outcome_var = "cd4_change", treatment_var = "treat",
    outcome_type = "continuous", effect_measure = "MD",
    subgroup_vars = dgm_subgroup_vars, subgroup_cuts = dgm_subgroup_cuts,
    model = "alt", k_treat = 1, k_inter = stage2_forecast$k_inter,
    adverse_outcome = FALSE, n_super = dgm_n_super, seed = seed_base,
    verbose = FALSE)
}

# Exact scoring frame: for continuous, df_super passes through unchanged —
# beta(Hhat) is an exact finite-population mean, zero Monte Carlo error.
eval_df <- fs_build_eval_frame(dgm, outcome_type = "continuous")

# Computed sampling scale, replacing the measured anchor.  fs_dgm_scale()
# enumerates sigma^2 + V_g[mu_0] + C_g[mu_0, tau] + V_g[tau]/2 on df_super, so
# the scale is a constant of the fixture rather than a number read off one run.
# The anchor it replaces (13.6786 at n = 1000) implied a bracket of 16,119
# against a residual variance of 16,256 -- below the noise floor, which on the
# true region is impossible since the bracket there is sigma^2 + V_Q[mu_0].
dgm_scale <- if (isTRUE(null_cell)) NULL else fs_dgm_scale(dgm)

# truth is ASSEMBLED from the calibrated object, as in the committed harness
# (mr_coverage_sweep_md_harm.qmd:187-200).  There is no `dgm$truth` element, so
# the as-received `dgm$truth %||% dgm` fell through to dgm itself, whose
# $effect_Q is also absent -- which made the stopifnot below read
# target_md_harm on both sides and assert -40 == -40.  A vacuous guard in the
# exact position where the calibration is supposed to be verified.
inQ        <- eval_df[[harm_col]] == 1L
truth <- list(effect_Q     = dgm$hazard_ratios$harm_subgroup,
              effect_Qc    = dgm$hazard_ratios$no_harm_subgroup,
              beta_inter   = dgm$model_params$beta_inter,
              prevalence_Q = mean(inQ))
# Calibration check applies to the ALTERNATIVE cell only; under model = "null"
# there is no target, and the meaningful assertion is that the null DGM is
# what it claims: Q is empty, effect_Q is undefined on an empty region, and
# the planted interaction is absent.
if (isTRUE(null_cell)) {
  stopifnot(sum(inQ) == 0L,
            is.na(truth$effect_Q),
            isTRUE(all.equal(truth$beta_inter, 0)))
  cat("NULL CELL: Q empty; effect_Q undefined (NA); beta_inter == 0 (no treatment-by-subgroup interaction).\n")
} else {
  stopifnot(abs(truth$effect_Q - target_md_harm) < 1e-8)
}

# ── FULL REFERENCE-TARGET SET, and why the obvious construction is vacuous ───
#
# The survival template carries three population references: the overall causal
# effect, the CDE theta-ddagger, and the marginal theta-dagger.  Those are
# genuinely DIFFERENT numbers there because the hazard ratio is NON-collapsible:
# cde_H = mean(exp(theta_1))/mean(exp(theta_0)) (a potential-outcome average)
# and marg_H = dgm$hr_H_true (a Cox FIT) do not coincide.
#
# TRAP, and the reason this block is written by hand.  The apparent continuous
# analogue -- calling fs_betaHhat_theta_dagger_check(outcome_type = "continuous")
# -- would produce a VACUOUS check.  On the continuous path that function
# delegates to .fs_region_effect(), which for continuous/count delegates to
# compute_aor() = mean(mu1[S]) - mean(mu0[S]) (betaHhat_truth.R, oc_analyses.R).
# The calibrator built effect_Q with the IDENTICAL arithmetic
# (generate_glm_dgm.R:578-583).  So theta-dagger and effect_Q would be the same
# expression evaluated twice, and asserting their equality proves nothing --
# the same failure mode as the `-40 == -40` guard noted above.
#
# A real check needs the two quantities computed by DIFFERENT ROUTES:
#
#   effect_Q / effect_Qc  STRUCTURAL, exact.  Potential-outcome average over the
#                         super-population; zero Monte Carlo error.
#   marg_H   / marg_Hc    FITTED.  lm(y_sim ~ treat_sim) on ONE realized draw of
#                         the whole super-population -- realized outcomes and
#                         realized randomization, so it carries sampling noise.
#
# Under an identity link the MD is collapsible, so these agree IN EXPECTATION.
# Reporting the gap against its own Monte Carlo scale therefore VERIFIES
# collapsibility on this DGM instead of assuming it.  For a non-collapsible
# measure the same two numbers would separate -- which is exactly what the
# survival cell shows, and why both are worth carrying.
truth$effect_ITT <- dgm$hazard_ratios$overall   # ITT analogue of hr_causal

marg_seed <- 20260628L    # fixed; this frame is for the CHECK only, never for
                          # beta(Hhat) scoring, which stays exact on df_super.
marg_frame <- simulate_from_glm_dgm(dgm, n = nrow(dgm$df_super),
                                    replace = FALSE, seed = marg_seed)
.fit_md <- function(d) {
  if (!nrow(d) || length(unique(d[[treat_name]])) < 2L)
    return(c(est = NA_real_, se = NA_real_))
  fit <- stats::lm(stats::reformulate(treat_name, response = outcome_name),
                   data = d)
  co <- summary(fit)$coefficients[treat_name, ]
  c(est = unname(co[["Estimate"]]), se = unname(co[["Std. Error"]]))
}
.mH  <- .fit_md(marg_frame[marg_frame[[harm_col]] == 1L, , drop = FALSE])
.mHc <- .fit_md(marg_frame[marg_frame[[harm_col]] == 0L, , drop = FALSE])
truth$marg_H     <- unname(.mH[["est"]]);  truth$marg_H_se  <- unname(.mH[["se"]])
truth$marg_Hc    <- unname(.mHc[["est"]]); truth$marg_Hc_se <- unname(.mHc[["se"]])

# Collapsibility read-out.  RAW scale throughout (negative = harm), matching
# effect_Q and betaHhat_H.  z is the gap in units of the fitted SE: |z| well
# under ~3 is agreement to Monte Carlo error.  Reported, NOT asserted -- a hard
# stopifnot here would fire on ordinary sampling noise.
.z_H  <- (truth$marg_H  - truth$effect_Q)  / truth$marg_H_se
.z_Hc <- (truth$marg_Hc - truth$effect_Qc) / truth$marg_Hc_se
cat(sprintf("DGM: effect_Q = %.10f | effect_Qc = %.10f | effect_ITT = %.10f | prevalence(Q) = %.4f | N_super = %d\n",
            truth$effect_Q, truth$effect_Qc, truth$effect_ITT,
            truth$prevalence_Q, nrow(eval_df)))
cat(sprintf("COLLAPSIBILITY CHECK (structural vs fitted, raw scale):\n"))
cat(sprintf("  Q  : structural %+.6f | fitted %+.6f (SE %.6f) | gap %+.6f | z = %+.2f\n",
            truth$effect_Q, truth$marg_H, truth$marg_H_se,
            truth$marg_H - truth$effect_Q, .z_H))
cat(sprintf("  Qc : structural %+.6f | fitted %+.6f (SE %.6f) | gap %+.6f | z = %+.2f\n",
            truth$effect_Qc, truth$marg_Hc, truth$marg_Hc_se,
            truth$marg_Hc - truth$effect_Qc, .z_Hc))
if (any(abs(c(.z_H, .z_Hc)) > 4, na.rm = TRUE))
  warning("collapsibility check: structural and fitted targets differ by more ",
          "than 4 SE; investigate the simulator or the orientation before ",
          "reading any downstream table.", call. = FALSE)


## ----machinery----------------------------------------------------------------
.classify <- function(sg_hat, tru) {
  out <- c(sens = NA_real_, spec = NA_real_, ppv = NA_real_, npv = NA_real_)
  if (is.null(sg_hat) || length(sg_hat) != length(tru)) return(out)
  sg_hat <- as.integer(sg_hat); tru <- as.integer(tru)
  tp <- sum(sg_hat==1 & tru==1); fp <- sum(sg_hat==1 & tru==0)
  fn <- sum(sg_hat==0 & tru==1); tn <- sum(sg_hat==0 & tru==0)
  out["sens"] <- if (tp+fn>0) tp/(tp+fn) else NA_real_
  out["spec"] <- if (tn+fp>0) tn/(tn+fp) else NA_real_
  out["ppv"]  <- if (tp+fp>0) tp/(tp+fp) else NA_real_
  out["npv"]  <- if (tn+fn>0) tn/(tn+fn) else NA_real_
  out
}

# One record; status distinguishes CONFIG-ERROR / NO-DETECTION / DETECTED.
.na_record <- function(sim_id) data.frame(
  sim_id = sim_id, status = NA_character_, detected = NA_integer_,
  mr_ok = NA_integer_, err_msg = NA_character_, mr_msg = NA_character_,
  n_sel = NA_integer_, n_harm = NA_integer_, n_true = NA_integer_,
  sg_def = NA_character_, covs = NA_character_,
  betaHhat_H = NA_real_, betaHhat_Hc = NA_real_,
  fb_secs = NA_real_, fit_mr_secs = NA_real_, fb_err = NA_character_,
  fb_src1 = NA_real_, fb_src2 = NA_real_, fb_nres = NA_integer_,
  # SCHEMA CHANGE: the complement (Hc) arm is now carried for EVERY estimator,
  # matching the survival template.  Previously only MR had one, so three of the
  # four estimators silently dropped Hhat^c -- and in the FB case the package had
  # already COMPUTED it (boot$Hc_estimates) and the recorder discarded it.
  # NOT schema-compatible with bundles written before this change: those lack the
  # or_Hc_/nv_Hc_/fb_Hc_ columns, so rbind() in combine mode will reject them.
  # Re-run any batch you intend to pool with new ones.
  or_H_est = NA_real_, or_H_lo = NA_real_, or_H_hi = NA_real_, or_H_se = NA_real_,
  or_Hc_est = NA_real_, or_Hc_lo = NA_real_, or_Hc_hi = NA_real_, or_Hc_se = NA_real_,
  nv_H_est = NA_real_, nv_H_lo = NA_real_, nv_H_hi = NA_real_, nv_H_se = NA_real_,
  nv_Hc_est = NA_real_, nv_Hc_lo = NA_real_, nv_Hc_hi = NA_real_, nv_Hc_se = NA_real_,
  fb_H_est = NA_real_, fb_H_lo = NA_real_, fb_H_hi = NA_real_, fb_H_se = NA_real_,
  fb_Hc_est = NA_real_, fb_Hc_lo = NA_real_, fb_Hc_hi = NA_real_, fb_Hc_se = NA_real_,
  mr_H_est = NA_real_, mr_H_lo = NA_real_, mr_H_hi = NA_real_, mr_H_se_ij = NA_real_,
  mr_Hc_est = NA_real_, mr_Hc_lo = NA_real_, mr_Hc_hi = NA_real_, mr_Hc_se_ij = NA_real_,
  ij_source = NA_character_,
  sens = NA_real_, spec = NA_real_, ppv = NA_real_, npv = NA_real_,
  stringsAsFactors = FALSE)

# Oracle: difference in means refit on a TRUE region (identity link).
# Now takes the membership vector so the SAME code serves Q and its complement
# (the survival template fits .cox_hr_ci() on both arms); previously this was
# hard-wired to flag_harm == 1 and the complement oracle did not exist.
.oracle_md_on <- function(df, keep) {
  if (!any(keep)) return(c(est=NA_real_, lo=NA_real_, hi=NA_real_, se=NA_real_))
  d <- df[keep, , drop = FALSE]
  if (length(unique(d[[treat_name]])) < 2L)
    return(c(est=NA_real_, lo=NA_real_, hi=NA_real_, se=NA_real_))
  fit <- stats::lm(stats::reformulate(treat_name, response = outcome_name),
                   data = d)
  co <- summary(fit)$coefficients[treat_name, ]
  est <- unname(co["Estimate"]); se <- unname(co["Std. Error"])
  s <- if (isTRUE(adverse_outcome)) 1 else -1        # oriented: + = harm
  # ORDER MATTERS.  This previously returned the vector reordered as
  # [c(1,3,2,4)] -- i.e. names in the order est, hi, lo, se -- and the recorder
  # assigned it POSITIONALLY, so or_*_lo received the upper bound and or_*_hi
  # the lower one.  Every oracle interval came out inverted: mean CI length
  # -3.92 x SE (reported as -78.45 with SE 20.01) and coverage identically 0.
  # The reorder served no purpose; it is removed, the redundant sort with it,
  # and the recorder now assigns BY NAME so a future reordering cannot repeat
  # this.  With s = -1 the bounds are already ordered: lo = -est - 1.96 se
  # < hi = -est + 1.96 se, because se > 0.
  c(est = s * est,
    lo  = s * est - 1.96 * se,
    hi  = s * est + 1.96 * se,
    se  = se)
}
.oracle_md    <- function(df) .oracle_md_on(df, df[[harm_col]] == 1L)
.oracle_md_c  <- function(df) .oracle_md_on(df, df[[harm_col]] == 0L)

# Inner (within-replicate) plan.  MUST be sequential under parallel_mode =
# "sims": the FB call at the bottom of record_replicate() passes inner_parallel
# UNCONDITIONALLY, so with replicates already fanned out across n_workers an
# unguarded multisession inner plan would nest n_workers x n_workers processes.
# Under "boots" the replicate loop is sequential, so the inner plan is where all
# the parallelism lives and multisession is correct. Same pattern as the
# committed sweep (mr_coverage_sweep_md_harm.qmd:141-143).
inner_parallel <- if (identical(parallel_mode, "sims")) {
  list(plan = "sequential")
} else {
  list(plan = "multisession", workers = n_workers)
}

record_replicate <- function(sim_id) {
  rec  <- .na_record(sim_id)
  sd_i <- seed_base + sim_id
  # Pin the RNG kind: simulate_from_glm_dgm() calls set.seed(seed) with
  # kind = NULL, inheriting the ambient generator (addendum D.2).
  RNGkind("L'Ecuyer-CMRG"); set.seed(sd_i)

  df <- tryCatch(simulate_from_glm_dgm(dgm, n = n_sample, seed = sd_i),
                 error = function(e) { rec$err_msg <<- conditionMessage(e); NULL })
  if (is.null(df)) { rec$status <- "CONFIG-ERROR"; return(rec) }
  df[[id_name]] <- df[[id_name]] %||% seq_len(nrow(df))
  if (k_random_noise > 0) {
    set.seed(sd_i + 10^7)
    for (j in seq_len(k_random_noise))
      df[[paste0("noise", j)]] <- stats::rnorm(nrow(df))
  }
  confs <- c(confounders_analysis,
             if (k_random_noise > 0) paste0("noise", seq_len(k_random_noise)))
  rec$n_true <- sum(df[[harm_col]] == 1L)
  .oH  <- .oracle_md(df);   .oHc <- .oracle_md_c(df)
  rec$or_H_est  <- .oH[["est"]];  rec$or_H_lo  <- .oH[["lo"]]
  rec$or_H_hi   <- .oH[["hi"]];   rec$or_H_se  <- .oH[["se"]]
  rec$or_Hc_est <- .oHc[["est"]]; rec$or_Hc_lo <- .oHc[["lo"]]
  rec$or_Hc_hi  <- .oHc[["hi"]];  rec$or_Hc_se <- .oHc[["se"]]

  t0 <- proc.time()[3]; fs.est <- NULL
  msgs <- tryCatch(utils::capture.output(
    fs.est <- suppressWarnings(forestsearch(
      df.analysis = df, confounders.name = confs,
      outcome.name = outcome_name, treat.name = treat_name, id.name = id_name,
      outcome_type = "continuous", effect_measure = "MD",
      effect.threshold = md_threshold, consistency.threshold = md_consistency,
      pconsistency.threshold = pconsistency, fs.splits = fs_splits,
      n.min = n_min, d0.min = d0_min, d1.min = d1_min, maxk = maxk,
      vi.grf.min = vi_grf_min, sg_focus = sg_focus,
      selection_rule = selection_rule,
      effect_neighborhood = effect_neighborhood,
      stop_threshold = stop_threshold,
      consistency_method = consistency_method,
      conf.cont_jcuts = fs_conf.cont_jcuts,
      use_lasso = use_lasso, use_dina = use_dina, use_grf = use_grf,
      use_twostage = use_twostage, is.RCT = is_rct,
      adverse_outcome = adverse_outcome,
      details = FALSE, quiet = FALSE, seedit = sd_i,
      # ALIGNED: pass inner_parallel UNCONDITIONALLY, as the survival template
      # does.  The conditional handed forestsearch() `parallel_args = NULL`
      # under "sims", and NULL is NOT the formal default -- an explicitly
      # supplied NULL has length 0, so forestsearch_main.R:1450
      # (`if (length(parallel_args) > 0)`) SKIPS plan setup entirely and the
      # search inherits whatever plan is ambient in the worker.  That happened
      # to be safe (future nests to sequential), but it left the search's
      # topology implicit and unlike the survival path.  inner_parallel already
      # resolves to list(plan = "sequential") under "sims", so this states the
      # same intent explicitly and makes "boots" hand the search the pool.
      parallel_args = inner_parallel,
      mr_inference = TRUE,
      mr_inference_args = list(ci_method = "ij", draws = mr_draws,
                               include_complement = TRUE))),
    type = "message"),
    error = function(e) {
      if (isTRUE(diag_mode)) warning(sprintf("FS-ERR sim %d: %s",
                                             sim_id, conditionMessage(e)))
      rec$err_msg <<- conditionMessage(e); character(0) })
  rec$fit_mr_secs <- proc.time()[3] - t0
  mr_lines <- grep("multiplier resampling|mr_inference", msgs, value = TRUE)
  if (length(mr_lines)) rec$mr_msg <- paste(mr_lines, collapse = " | ")
  if (!is.na(rec$err_msg)) { rec$status <- "CONFIG-ERROR"; return(rec) }
  if (is.null(fs.est) || is.null(fs.est$sg.harm)) {
    rec$status <- "NO-DETECTION"; rec$detected <- 0L; return(rec)
  }
  rec$status <- "DETECTED"; rec$detected <- 1L

  # The realized rule, as a single string.  fs.est$sg.harm is a CHARACTER
  # VECTOR of rule terms, not a list with a $definition element, so the
  # as-received `att$sg_def %||% fs.est$sg.harm$definition` gave NA_character_.
  # sg_def is the key fs_attach_betaHhat() joins on, so an NA here silently
  # routes every detected replicate into the undetected branch and hands it the
  # ITT complement as its target.  Collapse verbatim as both references do.
  rec$sg_def <- paste(fs.est$sg.harm, collapse = " & ")
  # `covs` (variable names only, for the identification figures) stays NA as
  # received: the survival template fills it with rule_covs(), which is not a
  # package export -- it is sourced from a helper that must sit in the render
  # working directory, and there is no such file here.  Putting the full rule
  # string in a column that downstream code reads as a variable list would be
  # worse than leaving it empty, so it is left empty and noted instead.

  # Classification of Hhat against the true harm flag, on the ANALYSIS sample.
  # fs.est$sg.harm$flag and fs.est$df_flag do not exist either; the per-subject
  # 0/1 membership vector is grp.consistency$sg.harm.id (survival
  # template:529-532), aligned to the rows of df.  Left as-received this block
  # would have written NA to n_harm/sens/spec/ppv/npv on every replicate --
  # readout 4 of the batch, gone, and gone quietly.
  sgv <- fs.est$grp.consistency$sg.harm.id
  if (!is.null(sgv) && length(sgv) == nrow(df)) {
    rec$n_harm <- sum(sgv == 1L, na.rm = TRUE)
    rec[c("sens","spec","ppv","npv")] <- as.list(.classify(sgv, df[[harm_col]]))
  }
  # Naive and MR columns live under fs.est$mr_inference, in the naive/debiased/
  # complement sub-lists -- NOT at the top level of fs.est and not under an
  # `fs.est$mr` element (there is none).  The as-received `%||%` chain over
  # fs.est[[nm]] therefore resolved to the NA seed on every replicate: a silent
  # all-NA MR arm that mr_ok, keyed on is.na(mr_H_est), would have reported as
  # a genuine MR failure.  Field paths taken verbatim from
  # mr_coverage_sweep_md_harm.qmd:365-382 and the survival template:506-526.
  #
  # mr_ok is keyed on the PRESENCE of the gate, not on a non-NA estimate, so an
  # MR that ran and legitimately returned NA is not confused with one that
  # never ran.  Detection is already recorded above and is independent of it.
  g <- fs.est$mr_inference
  rec$mr_ok <- as.integer(!is.null(g))
  if (!is.null(g)) {
    rec$n_sel     <- g$n_selected %||% NA_integer_
    rec$ij_source <- g$debiased$ij_source %||% NA_character_
    rec$nv_H_est <- g$naive$est    %||% NA_real_
    rec$nv_H_lo  <- g$naive$lower  %||% NA_real_
    rec$nv_H_hi  <- g$naive$upper  %||% NA_real_
    rec$nv_H_se  <- g$debiased$se_wald %||% NA_real_
    rec$mr_H_est <- g$debiased$est   %||% NA_real_
    rec$mr_H_lo  <- g$debiased$lower %||% NA_real_
    rec$mr_H_hi  <- g$debiased$upper %||% NA_real_
    rec$mr_H_se_ij <- g$debiased$se_ij %||% NA_real_
    if (is.list(g$complement) && !is.null(g$complement$debiased)) {
      # NAIVE complement, added with the Hc arm: the gate already returns it
      # (g$complement$naive) and the survival template records it as nv_Hc_*.
      rec$nv_Hc_est   <- g$complement$naive$est    %||% NA_real_
      rec$nv_Hc_lo    <- g$complement$naive$lower  %||% NA_real_
      rec$nv_Hc_hi    <- g$complement$naive$upper  %||% NA_real_
      rec$nv_Hc_se    <- g$complement$debiased$se_wald %||% NA_real_
      rec$mr_Hc_est   <- g$complement$debiased$est   %||% NA_real_
      rec$mr_Hc_lo    <- g$complement$debiased$lower %||% NA_real_
      rec$mr_Hc_hi    <- g$complement$debiased$upper %||% NA_real_
      rec$mr_Hc_se_ij <- g$complement$debiased$se_ij %||% NA_real_
    }
  }

  # FB: full bootstrap of the identical fitted procedure (Eq. 7, H2).
  if (run_fb) {
    tb <- proc.time()[3]
    # SILENCED (summary-upgrade spec A2): the function's own verbosity knob is
    # `details`, already FALSE by default and pinned here, but its progress
    # lines are gated on the ORIGINAL call's quiet flag
    # (bootstrap_dofuture_main.R:455-467), which this document sets FALSE --
    # so stdout is additionally captured.  The fb_err tryCatch is preserved
    # exactly as committed; on error fb stays NULL.
    fb <- NULL
    tryCatch(
      invisible(utils::capture.output(
        fb <- forestsearch_bootstrap_dofuture(
          fs.est, nb_boots = nb_boots, seed = sd_i,
          details = FALSE,
          parallel_args = inner_parallel),
        type = "output")),
      error = function(e) { rec$fb_err <<- conditionMessage(e); NULL })
    rec$fb_secs <- proc.time()[3] - tb
    # forestsearch_bootstrap_dofuture() returns the bias-corrected estimates in
    # `H_estimates`, a one-row data.table with columns H0/H1/H2 and, per
    # correction order, sdH*/H*_lower/H*_upper (bootstrap_calculations_helpers.R
    # get_dfRes()).  H2 is Leon Eq.7 -- both bias sources -- which is the FB
    # comparator this document is built around.  There is no top-level fb_H_est
    # and no fb$H2 list, so the as-received chain left all four columns NA.
    #
    # SCALE.  get_dfRes() is invoked with est.loghr = FALSE for identity-scale
    # measures (bootstrap_dofuture_main.R:557-563 lists MD among them), so H2
    # and sdH2 are already on the MD scale.  The survival template's
    # .se_beta_hr(sdH2, H2) = sdH2/H2 is the delta-method step from an HR SD to
    # a log-HR SE and is deliberately NOT carried over: on an identity scale
    # sdH2 IS the standard error, and dividing by H2 would be meaningless.
    #
    # The per-column %||% NA_real_ is retained from the template for the reason
    # given there: get_dfRes() omits H2 when the second-order adjustment is
    # unavailable, and a bare NULL assignment would DELETE the column from rec.
    # rbind() then errors on the short record and one malformed replicate takes
    # down the whole batch -- neither run path survives it.
    #
    # ORIENTATION.  FB returns the RAW outcome scale (negative = harm); the MR
    # gate returns the ORIENTED scale (positive = harm).  Verified directly on
    # replicate 1: the FB uncorrected estimate and the gate's naive estimate
    # are EXACT negatives,
    #
    #   H0            = -103.1366726933
    #   g$naive$est   =  103.1366726933
    #   |H0 + naive|  =  0.000e+00
    #
    # so this is a sign flip, not a re-estimation -- the same relationship, and
    # the same hazard, that commit 1830ca92 documents between the estimator
    # columns and betaHhat/md_trial_H.  fb_* is recorded ORIENTED here so that
    # it matches its siblings or_H_*, nv_H_* and mr_H_* in this schema; left
    # raw it would read as a sign error against them in every downstream table.
    # Interval bounds are negated AND swapped (a lower bound on -X is minus the
    # upper bound on X); sdH2 is a standard deviation and is unaffected.
    if (!is.null(fb) && !is.null(fb$H_estimates)) {
      he <- fb$H_estimates
      .orient <- if (isTRUE(adverse_outcome)) 1 else -1
      rec$fb_H_est <- .orient * (he$H2 %||% NA_real_)
      fb_ci <- sort(.orient * c(he$H2_lower %||% NA_real_,
                                he$H2_upper %||% NA_real_), na.last = TRUE)
      rec$fb_H_lo  <- fb_ci[1]
      rec$fb_H_hi  <- fb_ci[2]
      rec$fb_H_se  <- he$sdH2 %||% NA_real_   # scale-free, NOT negated

      # ---- Two-source decomposition of the Eq.7 correction ------------------
      # RAW scale (not oriented), so these compose additively with H0/H2 as
      # written in the source rather than needing a sign map.
      #
      #   fb_src1 = H_star     - H_obs      (resampling noise on the ORIGINAL
      #                                      subgroup: same rule, bootstrap data)
      #   fb_src2 = Hstar_star - Hstar_obs  (RE-SELECTION: the resample's own
      #                                      subgroup, scored on its own data
      #                                      vs on the original data)
      #
      # forestsearch_bootstrap_dofuture() does NOT return the four constituents
      # -- H_obs, H_star, Hstar_star, Hstar_obs are locals inside the resample
      # loop (bootstrap_analysis_dofuture.R:488, :697, :712) and only the two
      # combinations survive into boot$results. They are nonetheless EXACTLY
      # recoverable, because the two combinations are linearly independent in
      # the two unknowns:
      #
      #   H_biasadj_1 = H_obs - (Hstar_star - Hstar_obs)
      #   H_biasadj_2 = 2*H_obs - (H_star + Hstar_star - Hstar_obs)
      #
      #   => fb_src2 = Hstar_star - Hstar_obs = H_obs - H_biasadj_1
      #   => fb_src1 = H_star     - H_obs     = H_biasadj_1 - H_biasadj_2
      #
      # H_obs is the original-sample fit on the original subgroup, which is
      # exactly H0 on the identity scale (est.loghr = FALSE for MD, so ci_est()
      # is a pass-through; verified in V_fb_adjudication.md section V1, where
      # H0 = -103.1366726933 is the exact negative of the gate's naive
      # estimate). Averaged over resamples, na.rm because a failed draw sets
      # both adjustments NA (bootstrap_analysis_dofuture.R:729-731).
      br <- fb$results
      if (!is.null(br) && all(c("H_biasadj_1","H_biasadj_2") %in% names(br))) {
        b1 <- br$H_biasadj_1; b2 <- br$H_biasadj_2
        ok <- is.finite(b1) & is.finite(b2)
        rec$fb_nres <- sum(ok)
        if (any(ok)) {
          rec$fb_src1 <- mean(b1[ok] - b2[ok])
          rec$fb_src2 <- mean((he$H0 %||% NA_real_) - b1[ok])
        }
      }
    }

    # FB COMPLEMENT.  forestsearch_bootstrap_dofuture() returns Hc_estimates on
    # the SAME call, with the same H0/H1/H2 + sdH*/H*_lower/H*_upper layout
    # (get_dfRes() is invoked once per arm, bootstrap_dofuture_main.R:573-596).
    # It was already being COMPUTED here and then discarded; the survival
    # template records it as fb_Hc_*.  Oriented and bound-swapped exactly as the
    # H arm above, for the same reason.  No src decomposition on this arm: the
    # two bias sources are recorded for the harm block only, as before.
    if (!is.null(fb) && !is.null(fb$Hc_estimates)) {
      hce <- fb$Hc_estimates
      .orient_c <- if (isTRUE(adverse_outcome)) 1 else -1
      rec$fb_Hc_est <- .orient_c * (hce$H2 %||% NA_real_)
      fbc_ci <- sort(.orient_c * c(hce$H2_lower %||% NA_real_,
                                   hce$H2_upper %||% NA_real_), na.last = TRUE)
      rec$fb_Hc_lo  <- fbc_ci[1]
      rec$fb_Hc_hi  <- fbc_ci[2]
      rec$fb_Hc_se  <- hce$sdH2 %||% NA_real_
    }
  }
  rec
}


## ----run-batch----------------------------------------------------------------
#| output: false
#| message: false
# Nothing from the replicate loop prints to the document (summary-upgrade spec
# A1/A3): the batch-means kable now lives in the summary layer below, and this
# chunk's only side effect is the saved bundle.
#
# run_mode = "batch" runs and saves ONE seed-disjoint piece; run_mode =
# "combine" runs nothing and merges the pieces already on disk.  The combine
# branch is ported from the survival template, which had it and this document
# did not -- so the md40 batches written so far (res_1_50, res_51_100) could not
# be pooled by the document that produced them.
combine_glob  <- file.path(results_dir, sprintf("%s_res_*.rds", rds_stem))
combine_files <- NULL    # optional explicit character vector; NULL -> the glob
save_combined <- TRUE
combined_path <- NULL    # NULL -> "<results_dir>/<stem>_combined_<min>_<max>.rds",
                         # deliberately OUTSIDE combine_glob so a re-combine
                         # never re-ingests its own pooled output.

if (identical(run_mode, "batch")) {

sim_ids <- seq(sim_id_start, length.out = if (isTRUE(diag_mode)) 5L else n_sims)
t_all <- proc.time()[3]
if (identical(parallel_mode, "sims")) {
  registerDoFuture(); future::plan(future::multisession, workers = n_workers)
  results <- foreach(i = sim_ids, .combine = rbind,
                     .options.future = list(seed = TRUE)) %dofuture%
    record_replicate(i)
  future::plan(future::sequential)
} else {
  results <- do.call(rbind, lapply(sim_ids, record_replicate))
}
elapsed <- proc.time()[3] - t_all   # run-loop wall-clock, for the timing section

# Exact beta(Hhat) at the realized rules, on the fixed evaluation frame.
# fs_attach_betaHhat() is a POST-LOOP, whole-frame call: its signature is
# (results, frame, focus, ...) -- a results data.frame carrying an sg_def
# column, keyed to the rules it finds there -- not a per-replicate call on an
# fs.est object.  The as-received document called it inside record_replicate()
# as fs_attach_betaHhat(fs.est, eval_df = eval_df, outcome_type = ...), which
# fails three ways at once: `eval_df` is not a formal (the function takes no
# ..., so R raises "unused argument" and the tryCatch-free call aborts the
# replicate), `focus` has no default and was not supplied, and `fs.est` is not
# a results frame.  Call signature and argument values follow the committed
# mr_coverage_sweep_md_harm.qmd:411-413 exactly -- same DGM, same outcome type,
# same focus.
#
# Undetected and CONFIG-ERROR rows carry sg_def = NA and take the no-subgroup
# record (nH_eval = 0, betaHhat_Hc = the ITT effect), so they remain part of
# the Hhat^c denominator instead of dropping out.
results <- fs_attach_betaHhat(
  results, eval_df, focus = "harm", outcome_type = "continuous",
  effect_measure = "MD")

out <- file.path(results_dir, sprintf("%s_res_%d_%d.rds", rds_stem,
                 min(sim_ids), max(sim_ids)))
.payload <- list(results = results, truth = truth,
             meta = list(n_sample = n_sample, n_sims = length(sim_ids),
                         nb_boots = nb_boots, mr_draws = mr_draws,
                         sg_focus = sg_focus, selection_rule = selection_rule,
                         consistency_method = consistency_method,
                         subgroup_method = subgroup_method,
                         target_md_harm = target_md_harm,
                         effect_threshold = md_threshold,
                         consistency_threshold = md_consistency,
                         adverse_outcome = adverse_outcome,
                         seed_base = seed_base, sim_id_start = sim_id_start,
                         parallel_mode = parallel_mode, n_workers = n_workers,
                         pkg_version = as.character(utils::packageVersion("forestsearch")),
                         built_at = Sys.time()))

# Scale and OC summary travel WITH the payload.  Computed here, where the DGM
# and the results are both in hand, so a document's tables and its payload
# cannot disagree and the analytic documents read numbers instead of
# transcribing literals.
if (!is.null(dgm_scale)) .payload$scale <- dgm_scale
.payload$oc    <- fs_mr_oc_summary(.payload)
saveRDS(.payload[intersect(c("results", "truth", "scale", "oc", "meta"),
                           names(.payload))], out)

} else if (identical(run_mode, "combine")) {
  # ---- Merge seed-disjoint batches ----------------------------------------
  # Ported from the survival template's combine branch, with the same three
  # refusals: batches that disagree on a count-invariant knob are not poolable;
  # a repeated sim_id means the seed ranges overlapped, so the pool would
  # double-count; and the pooled bundle is written outside combine_glob.
  files <- if (!is.null(combine_files)) combine_files else Sys.glob(combine_glob)
  files <- files[file.exists(files)]
  if (!length(files))
    stop("combine mode: no batch files found (combine_glob = '",
         combine_glob, "').", call. = FALSE)
  bundles <- lapply(files, readRDS)

  .meta_get <- function(b, k) b$meta[[k]] %||% NA
  for (k in c("n_sample", "nb_boots", "mr_draws", "subgroup_method",
              "sg_focus", "target_md_harm", "seed_base")) {
    vals <- unique(lapply(bundles, .meta_get, k = k))
    if (length(vals) > 1L)
      stop(sprintf("combine mode: batches disagree on meta$%s (%s).",
                   k, paste(unlist(vals), collapse = ", ")), call. = FALSE)
  }

  truths <- lapply(bundles, `[[`, "truth")
  if (length(truths) > 1L &&
      !all(vapply(truths[-1], identical, logical(1), truths[[1]])))
    warning("combine mode: truth targets differ across batches; using the first.",
            call. = FALSE)
  truth <- truths[[1]]

  # Schema guard: the Hc arm was added to .na_record after the first md40
  # batches were written, so an older bundle has FEWER columns and rbind()
  # would abort the whole render with a column-count error.  Name the offender.
  nms <- lapply(bundles, function(b) names(b$results))
  if (length(unique(nms)) > 1L) {
    ref <- nms[[1]]
    bad <- which(!vapply(nms, identical, logical(1), ref))
    stop(sprintf(paste("combine mode: batch files have different column sets;",
                       "%s differ(s) from %s.\n  Most likely they predate the",
                       "Hc-arm schema change -- re-run them before pooling.",
                       "\n  Missing here: %s"),
                 paste(basename(files)[bad], collapse = ", "),
                 basename(files)[1],
                 paste(setdiff(ref, nms[[bad[1]]]), collapse = ", ")),
         call. = FALSE)
  }

  results <- do.call(rbind, lapply(bundles, `[[`, "results"))
  dup <- unique(results$sim_id[duplicated(results$sim_id)])
  if (length(dup))
    stop(sprintf(paste("combine mode: %d duplicated sim_id across batches",
                       "(e.g. %s); the seed ranges overlap, so the pool is NOT",
                       "independent."), length(dup),
                 paste(utils::head(sort(dup), 5L), collapse = ", ")),
         call. = FALSE)

  # Make the report self-describing from the FILES, not the setup knobs.
  n_sample <- bundles[[1]]$meta$n_sample %||% n_sample
  nb_boots <- bundles[[1]]$meta$nb_boots %||% nb_boots
  mr_draws <- bundles[[1]]$meta$mr_draws %||% mr_draws
  run_fb   <- !is.null(nb_boots) && nb_boots >= 1L
  sim_ids  <- sort(results$sim_id)
  n_sims   <- length(sim_ids)
  elapsed  <- NA_real_    # no single wall-clock when pooling; timing prints "-"

  rng <- range(results$sim_id)
  cat(sprintf("Combined %d batch file(s) -> %d rows, sim_id %d-%d.\n",
              length(files), nrow(results), rng[1], rng[2]))
  cat("  files: ", paste(basename(files), collapse = ", "), "\n", sep = "")

  if (isTRUE(save_combined)) {
    cpath <- combined_path %||%
      file.path(results_dir, sprintf("%s_combined_%d_%d.rds",
                                     rds_stem, rng[1], rng[2]))
    .cpayload <- list(results = results, truth = truth,
                 meta = list(n_sample = n_sample, n_sims = n_sims,
                             nb_boots = nb_boots, mr_draws = mr_draws,
                             sg_focus = bundles[[1]]$meta$sg_focus %||% NA_character_,
                             subgroup_method = bundles[[1]]$meta$subgroup_method %||% NA_character_,
                             target_md_harm = bundles[[1]]$meta$target_md_harm %||% NA_real_,
                             adverse_outcome = bundles[[1]]$meta$adverse_outcome %||% NA,
                             seed_base = bundles[[1]]$meta$seed_base %||% NA_integer_,
                             sim_id_start = rng[1], sim_id_end = rng[2],
                             n_batches = length(files),
                             source_files = basename(files),
                             pkg_version = bundles[[1]]$meta$pkg_version %||% NA_character_,
                             built_at = Sys.time()))
    # Pooled bundles carry scale and oc too; omitting them here would make a
    # combined payload silently poorer than the batches it merges.
    if (!is.null(dgm_scale)) .cpayload$scale <- dgm_scale
    .cpayload$oc    <- fs_mr_oc_summary(.cpayload)
    saveRDS(.cpayload[intersect(c("results", "truth", "scale", "oc", "meta"),
                              names(.cpayload))], cpath)
    cat(sprintf("Saved pooled bundle (%d rows, %d batches) to %s\n",
                nrow(results), length(files), cpath))
  }

} else {
  stop("run_mode must be \"batch\" or \"combine\"; got \"", run_mode, "\".",
       call. = FALSE)
}


