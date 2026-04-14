#!/usr/bin/env Rscript
# ===================================================================
# qc_null_classification.R
#
# Quick check: verify flag_harm = 0 under null for both survival
# and GLM DGMs, and that classification metrics follow León et al.
# (2024, Section 3): sens(H) = NA, ppv(H) = 0.
#
# Run after: devtools::document() -> devtools::install() -> restart R
# Expected: 6/6 checks pass
# ===================================================================

library(forestsearch)
library(survival)
library(speff2trial)
library(data.table)

pass <- 0L; fail <- 0L
check <- function(cond, label) {
  if (cond) {
    cat(sprintf("  PASS: %s\n", label))
    pass <<- pass + 1L
  } else {
    cat(sprintf("  FAIL: %s\n", label))
    fail <<- fail + 1L
  }
}

# ── 1. Survival null DGM ─────────────────────────────────────────
cat("--- Survival Null DGM ---\n")
dgm_s <- setup_gbsg_dgm(model = "null", k_treat = 1.0, verbose = FALSE)
check(sum(dgm_s$df_super$flag_harm) == 0,
  sprintf("flag_harm = 0 (got %d)", sum(dgm_s$df_super$flag_harm)))

# ── 2. GLM null DGM ─────────────────────────────────────────────
cat("\n--- GLM Null DGM ---\n")
data(ACTG175)
actg <- ACTG175[ACTG175$arms %in% c(1, 3), ]
actg$y_binary <- as.integer(actg$cd420 > actg$cd40)
actg$treat <- ifelse(actg$arms == 1, 0L, 1L)
actg$z1 <- as.factor(as.integer(actg$age > 34))
actg$z2 <- as.factor(as.integer(actg$cd80 > 1292))
set.seed(99)
for (i in 3:6) actg[[paste0("z", i)]] <- as.factor(sample(0:1, nrow(actg), TRUE))

dgm_g <- generate_glm_dgm(
  data = actg, factor_vars = paste0("z", 1:6),
  outcome_var = "y_binary", treatment_var = "treat",
  outcome_type = "binary", effect_measure = "OR",
  subgroup_vars = c("z1", "z2"), subgroup_cuts = list(z1 = 1L, z2 = 1L),
  model = "null", n_super = 5000L, seed = 42L
)

check(sum(dgm_g$df_super$flag_harm) == 0,
  sprintf("flag_harm = 0 (got %d)", sum(dgm_g$df_super$flag_harm)))

check(identical(dgm_g$model, "null"),
  sprintf("dgm$model = 'null' (got '%s')", dgm_g$model))

# ── 3. End-to-end: run_simulation_analysis on GLM null ───────────
cat("\n--- run_simulation_analysis (GLM null) ---\n")
r <- run_simulation_analysis(
  sim_id = 1L, dgm = dgm_g, n_sample = 800L,
  confounders_base = paste0("z", 1:6),
  n_add_noise = 0L,
  run_fs = TRUE, run_fs_grf = FALSE, run_grf = FALSE,
  fs_params = list(
    outcome.name = "y_sim", event.name = "y_sim",
    treat.name = "treat_sim", id.name = "id",
    outcome_type = "binary", effect_measure = "OR",
    use_lasso = FALSE, use_grf = TRUE,
    hr.threshold = 1.667, hr.consistency = 1.25,
    pconsistency.threshold = 0.80,
    fs.splits = 100L, n.min = 40L, d0.min = 8L, d1.min = 8L,
    maxk = 2L, vi.grf.min = -0.2,
    parallel_args = list(plan = "sequential", workers = 1,
                          show_message = FALSE)
  ),
  verbose = FALSE
)

cat(sprintf("  FS: any.H=%d, sens=%s, ppv=%s\n",
    r$any.H[r$analysis == "FS"],
    ifelse(is.na(r$sens[r$analysis == "FS"]), "NA",
           sprintf("%.3f", r$sens[r$analysis == "FS"])),
    ifelse(is.na(r$ppv[r$analysis == "FS"]), "NA",
           sprintf("%.3f", r$ppv[r$analysis == "FS"]))))

if (r$any.H[r$analysis == "FS"] == 1) {
  check(is.na(r$sens[r$analysis == "FS"]),
    "sens = NA under null (Q = emptyset)")
  check(r$ppv[r$analysis == "FS"] == 0,
    "ppv = 0 under null")
} else {
  check(is.na(r$sens[r$analysis == "FS"]),
    "sens = NA (no detection)")
  check(is.na(r$ppv[r$analysis == "FS"]),
    "ppv = NA (no detection)")
}

# ── 4. Sanity: Alt DGM still works ──────────────────────────────
cat("\n--- GLM Alt DGM (sanity) ---\n")
dgm_a <- generate_glm_dgm(
  data = actg, factor_vars = paste0("z", 1:6),
  outcome_var = "y_binary", treatment_var = "treat",
  outcome_type = "binary", effect_measure = "OR",
  subgroup_vars = c("z1", "z2"), subgroup_cuts = list(z1 = 1L, z2 = 1L),
  model = "alt", k_inter = 1.15, n_super = 5000L, seed = 42L
)
prev <- round(100 * mean(dgm_a$df_super$flag_harm), 1)
check(prev > 5,
  sprintf("Alt flag_harm prevalence = %.1f%% (should be > 5%%)", prev))

# ── Summary ─────────────────────────────────────────────────────
cat(sprintf("\n=== %d passed, %d failed ===\n", pass, fail))
if (fail > 0) stop(sprintf("QC FAILED: %d checks failed", fail))
cat("Null classification fix verified.\n")
