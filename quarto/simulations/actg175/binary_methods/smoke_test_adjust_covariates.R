# =============================================================================
# Smoke tests: adjust_covariates (Cox covariate adjustment)
# -----------------------------------------------------------------------------
# Exercises the new `adjust_covariates` argument end-to-end across:
#   Layer 1 (candidate search), Layer 2 (consistency), Layer 3 (bootstrap).
#
# Run AFTER:  devtools::document(); devtools::install()
# (install, not load_all, so multisession/callr workers see the changes)
#
# Dataset: survival::gbsg (always available). hormon = treatment (0/1),
# rfstime/status = outcome, with several prognostic covariates.
#
# Design choice: confounders.name (candidate-defining pool) is kept disjoint
# from the adjustment covariates (meno, age) so the demonstration is clean and
# free of defining-variable overlap. strata()/linear terms both shown.
# =============================================================================

library(forestsearch)
library(survival)

# ---- Data prep --------------------------------------------------------------
gbsg <- survival::gbsg
df <- within(gbsg, {
  tte   <- rfstime / 30.4375      # days -> months
  event <- status
  treat <- hormon
  meno  <- factor(meno)           # stratification factor (NOT a confounder here)
})

confs <- c("size", "grade", "nodes", "pgr", "er")   # candidate-defining pool

# Common, deliberately small/fast settings for a smoke test
fs_common <- list(
  df.analysis      = df,
  outcome.name     = "tte",
  event.name       = "event",
  treat.name       = "treat",
  id.name          = "pid",    # gbsg's subject id column (NOT the default "id")
  confounders.name = confs,
  hr.threshold     = 1.10,
  hr.consistency   = 1.0,
  pconsistency.threshold = 0.80,
  fs.splits        = 200,      # small for speed; bump up for real runs
  maxk             = 2,
  use_lasso        = FALSE,    # disable pruning so candidate pool is stable
  use_grf          = FALSE,    # skip GRF screening to keep it fast/deterministic
  parallel_args    = list(plan = "multisession", workers = 2,
                          show_message = FALSE),
  seedit           = 8316951,
  details          = FALSE
)

run_fs <- function(label, ...) {
  cat("\n----", label, "----\n")
  t0 <- proc.time()[3]
  fs <- tryCatch(
    do.call(forestsearch, modifyList(fs_common, list(...))),
    error = function(e) { cat("  ERROR:", conditionMessage(e), "\n"); NULL }
  )
  dt <- round(proc.time()[3] - t0, 1)
  if (!is.null(fs)) {
    sg <- fs$sg.harm
    cat("  ran OK in", dt, "s | harm subgroup:",
        if (is.null(sg)) "none found" else paste(sg, collapse = " & "), "\n")
  }
  invisible(fs)
}

# =============================================================================
# 0. INSTANT primitive check (no full pipeline) -- run this first
#    Confirms the adjusted Cox formula + scorer agree with manual coxph.
# =============================================================================
cat("=== 0. Primitive sanity (instant) ===\n")
sub <- transform(df, Y = tte, Event = event, Treat = treat)
hr_unadj  <- get_split_hr_fast(sub)
hr_strata <- get_split_hr_fast(sub, adjust_covariates = "strata(meno)")
hr_manual <- exp(coef(coxph(Surv(Y, Event) ~ Treat + strata(meno),
                            data = sub))[["Treat"]])
cat(sprintf("  get_split_hr_fast strata=%.4f vs manual=%.4f  -> %s\n",
            hr_strata, hr_manual,
            if (abs(hr_strata - hr_manual) < 1e-6) "MATCH" else "MISMATCH"))
cat(sprintf("  unadjusted HR=%.4f, strata(meno) HR=%.4f\n", hr_unadj, hr_strata))
cat("  build_cox_formula(strata): ",
    deparse(build_cox_formula("tte", "event", "treat",
                              adjust_covariates = c("strata(meno)", "age"))), "\n")

# =============================================================================
# 1. BASELINE: unadjusted (adjust_covariates = NULL == previous behavior)
# =============================================================================
fs0 <- run_fs("1. Unadjusted (baseline)")                       # default NULL

# =============================================================================
# 2. STRATIFIED adjustment: stratify baseline hazard by meno
# =============================================================================
fs1 <- run_fs("2. adjust_covariates = strata(meno)",
              adjust_covariates = "strata(meno)")

# =============================================================================
# 3. LINEAR adjustment: include age as a linear covariate
# =============================================================================
fs2 <- run_fs("3. adjust_covariates = age (linear)",
              adjust_covariates = "age")

# =============================================================================
# 4. MIXED: stratify by meno AND adjust linearly for age
# =============================================================================
fs3 <- run_fs("4. adjust_covariates = c(strata(meno), age)",
              adjust_covariates = c("strata(meno)", "age"))

# =============================================================================
# 5. GUARD: adjust_covariates + propensity-score adjustment is rejected
#    (these are mutually exclusive). Expect a clear ERROR, not a crash.
# =============================================================================
fs_guard <- run_fs("5. Guard: adjust_covariates + ps_method (expect ERROR)",
                   adjust_covariates = "age",
                   is.RCT = FALSE, ps_method = "logistic")

# =============================================================================
# 6. ERROR HANDLING: unknown / reserved column names (expect ERRORs)
# =============================================================================
fs_e1 <- run_fs("6a. Unknown column (expect ERROR)",
                adjust_covariates = "not_a_column")
fs_e2 <- run_fs("6b. Reserved name (expect ERROR)",
                adjust_covariates = "strata(Treat)")

# =============================================================================
# 7. BOOTSTRAP (Layer 3): bias correction on an ADJUSTED fit
#    Slower -- only run if Layer 1/2 above produced a subgroup (fs1$sg.harm).
#    Confirms the reported HR is on the same adjusted scale as validation.
# =============================================================================
if (!is.null(fs1) && !is.null(fs1$sg.harm)) {
  cat("\n---- 7. Bootstrap on adjusted fit (Layer 3) ----\n")
  fs1_bc <- tryCatch(
    forestsearch_bootstrap_dofuture(
      fs.est        = fs1,
      nb_boots      = 50,                 # tiny for smoke; use >=1000 for real
      parallel_args = list(plan = "multisession", workers = 2)
    ),
    error = function(e) { cat("  ERROR:", conditionMessage(e), "\n"); NULL }
  )
  if (!is.null(fs1_bc)) cat("  bootstrap ran OK; class:", class(fs1_bc)[1], "\n")
} else {
  cat("\n---- 7. Bootstrap skipped (no harm subgroup from step 2) ----\n")
}

# =============================================================================
# EXPECTED OUTCOMES
#   0      : MATCH printed; primitive correct.
#   1-4    : "ran OK"; subgroup may differ (or appear/disappear) across
#            adjustment modes -- that is the feature working, not a bug.
#   5,6a,6b: "ERROR: ..." printed (guards firing as intended).
#   7      : "bootstrap ran OK" if a subgroup was found in step 2.
# =============================================================================
cat("\n=== smoke tests complete ===\n")
