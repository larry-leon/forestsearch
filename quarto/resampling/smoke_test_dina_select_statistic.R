# ===========================================================================
# Minimal smoke test: select_statistic = "effect" for DINA and GRF
#   + Tier-1 bootstrap reference (the Layer-2 gap)
# ---------------------------------------------------------------------------
# DINA (subgroup_method = "dina", dina_args$select_statistic):
#   "dina"   -> winner = argmax DINA subgroup-mean tau-hat (legacy)
#   "effect" -> winner = argmax inferential Cox HR over the same family
# GRF  (subgroup_method = "grf", grf_selection = "frontier", grf_select_statistic):
#   "dr"     -> winner = frontier pick on the DR score (legacy)
#   "effect" -> winner = frontier pick re-ranked on the inferential Cox HR
# In both, "effect" is scored with the SAME estimator the Tier-2 gate de-biases.
# Then per mode: Tier-2 de-biased gate vs Tier-1 bootstrap H2 (Leon 2024 Eq. 7,
# BOTH bias sources); the gap (Tier-2 debiased - bootstrap H2) is the
# family-construction (Layer-2) optimism the single-fit gate cannot see.
#
# Self-contained: simulates its own planted-harm Cox DGM, no external helpers.
# REINSTALL FIRST so the new keys are recognised:  devtools::install()
# ===========================================================================

suppressMessages({
  library(forestsearch)
  library(survival)
  library(future)
})

# --- Setup knobs -----------------------------------------------------------
NB_BOOTS      <- 300L                                   # Tier-1 boots (>=100 avoids the unreliable-CI warning; prod 200-500)
BOOT_WORKERS  <- min(16L, parallel::detectCores())      # bootstrap parallelism
RUN_BOOTSTRAP <- TRUE                                    # FALSE to skip the (slower) Tier-1 arm; GRF boots dominate runtime

# --- 1. Planted-harm survival DGM ------------------------------------------
# Harm subgroup = low ER ({er <= 6}): treatment INCREASES hazard there
# (HR ~ 2.5), roughly neutral elsewhere.  Treatment is randomized.
set.seed(20260610)
n <- 1000L
er        <- round(rexp(n, rate = 1 / 20))
age       <- round(rnorm(n, 60, 12))
biomarker <- rbinom(n, 1, 0.5)
treat     <- rbinom(n, 1, 0.5)

harm    <- as.integer(er <= 6)
b_treat <- ifelse(harm == 1L, log(2.5), log(0.95))
lp      <- 0.015 * (age - 60) + 0.20 * biomarker + treat * b_treat
ttime   <- rexp(n, rate = 0.10 * exp(lp))
ctime   <- rexp(n, rate = 0.05)
tte     <- pmin(ttime, ctime)
event   <- as.integer(ttime <= ctime)

dat  <- data.frame(id = seq_len(n), tte = tte, event = event,
                   treat = treat, er = er, age = age, biomarker = biomarker)
covs <- c("er", "age", "biomarker")

cat(sprintf("Planted harm subgroup {er <= 6}: n = %d (%.1f%%) | events = %d/%d\n",
            sum(harm), 100 * mean(harm), sum(event), n))

# --- 2. Drivers + shared summaries -----------------------------------------
run_dina <- function(stat) {
  forestsearch(
    df.analysis = dat, outcome.name = "tte", event.name = "event",
    treat.name = "treat", id.name = "id", confounders.name = covs,
    outcome_type = "survival", subgroup_method = "dina",
    hr.threshold = 1.25, n.min = 60, mr_inference = TRUE,
    dina_args = list(select_statistic = stat),
    parallel_args = list(plan = "sequential", show_message = FALSE),
    details = FALSE, quiet = TRUE)
}

run_grf <- function(stat) {
  forestsearch(
    df.analysis = dat, outcome.name = "tte", event.name = "event",
    treat.name = "treat", id.name = "id", confounders.name = covs,
    outcome_type = "survival", subgroup_method = "grf",
    grf_selection = "frontier", grf_select_statistic = stat,
    hr.threshold = 1.25, n.min = 60, mr_inference = TRUE,
    parallel_args = list(plan = "sequential", show_message = FALSE),
    details = FALSE, quiet = TRUE)
}

summarize <- function(fs, tag) {
  cat("\n==", tag, "==\n")
  if (is.null(fs$sg.harm)) { cat("  no subgroup selected\n"); return(invisible(FALSE)) }
  cat("  selected sg.harm:", paste(fs$sg.harm, collapse = " & "), "\n")
  g <- fs$mr_inference
  if (is.null(g)) { cat("  mr_inference = NULL (gate did not run)\n"); return(invisible(FALSE)) }
  cat(sprintf("  n_family = %d   n_selected = %d   measure = %s\n",
              g$n_family, g$n_selected, g$measure))
  cat(sprintf("  naive    est = %.3f  [%.3f, %.3f]\n",
              g$naive$est, g$naive$lower, g$naive$upper))
  cat(sprintf("  debiased est = %.3f  [%.3f, %.3f]  (se_ij = %.3f)\n",
              g$debiased$est, g$debiased$lower, g$debiased$upper, g$debiased$se_ij))
  cat(sprintf("  selection_bias = %.4f   harm_flag = %s\n",
              g$selection_bias, g$harm_flag))
  invisible(is.finite(g$debiased$est))
}

# Cross-mode comparison (native-statistic vs effect) for a method
cross_mode <- function(fs_a, fs_b, lab_a, lab_b, method) {
  pass <- function(fs) !is.null(fs$mr_inference) && is.finite(fs$mr_inference$debiased$est)
  la <- if (!is.null(fs_a$sg.harm)) paste(fs_a$sg.harm, collapse = " & ") else NA
  lb <- if (!is.null(fs_b$sg.harm)) paste(fs_b$sg.harm, collapse = " & ") else NA
  cat(sprintf("\n[%s] PASS (finite debiased CI): %s = %s | %s = %s\n",
              method, lab_a, isTRUE(pass(fs_a)), lab_b, isTRUE(pass(fs_b))))
  cat(sprintf("  subgroup identical across modes: %s\n", identical(la, lb)))
  if (!identical(la, lb))
    cat(sprintf("    %-6s -> %s\n    %-6s -> %s\n", lab_a, la, lab_b, lb))
  if (pass(fs_a) && pass(fs_b))
    cat(sprintf("  naive HR: %s = %.3f | %s = %.3f  (effect >= native expected when winners differ)\n",
                lab_a, fs_a$mr_inference$naive$est, lab_b, fs_b$mr_inference$naive$est))
  invisible(NULL)
}

# --- 3. DINA: dina vs effect -----------------------------------------------
cat("\n", strrep("#", 60), "\n  DINA  (subgroup_method = 'dina')\n", strrep("#", 60), "\n", sep = "")
fs_dina <- run_dina("dina");   summarize(fs_dina, "select_statistic = 'dina'")
fs_deff <- run_dina("effect"); summarize(fs_deff, "select_statistic = 'effect'")
cross_mode(fs_dina, fs_deff, "dina", "effect", "DINA")

# --- 4. GRF (frontier): dr vs effect ---------------------------------------
cat("\n", strrep("#", 60), "\n  GRF   (subgroup_method = 'grf', grf_selection = 'frontier')\n", strrep("#", 60), "\n", sep = "")
fs_gdr  <- run_grf("dr");     summarize(fs_gdr,  "grf_select_statistic = 'dr'")
fs_geff <- run_grf("effect"); summarize(fs_geff, "grf_select_statistic = 'effect'")
cross_mode(fs_gdr, fs_geff, "dr", "effect", "GRF")

# --- 5. Tier-1 bootstrap reference: the Layer-2 gap ------------------------
# forestsearch_bootstrap_dofuture() re-runs the FULL selection per resample
# (inheriting this fit's call via args_call_all, incl. the select_statistic),
# so H2 captures Layer-1 (winner selection) AND Layer-2 (family construction) --
# the Tier-1 reference the single-fit Tier-2 gate only approximates.
boot_compare <- function(fs, tag) {
  if (is.null(fs$sg.harm) || is.null(fs$mr_inference)) {
    cat(sprintf("\n== %-22s : skipped (no subgroup / no gate)\n", tag)); return(invisible(NULL))
  }
  g <- fs$mr_inference
  boot <- tryCatch(
    forestsearch_bootstrap_dofuture(
      fs.est = fs, nb_boots = NB_BOOTS, seed = 8316951L,
      show_three = FALSE, details = FALSE,
      parallel_args = list(plan = "multisession", workers = BOOT_WORKERS,
                           show_message = FALSE)),
    error = function(e) { cat("  bootstrap error:", conditionMessage(e), "\n"); NULL })
  plan("sequential")
  if (is.null(boot) || is.null(boot$H_estimates)) {
    cat(sprintf("\n== %-22s : bootstrap returned no H_estimates\n", tag)); return(invisible(NULL))
  }
  he  <- boot$H_estimates                          # H2 = Leon 2024 Eq. 7
  t1  <- as.numeric(he$H2); t1l <- as.numeric(he$H2_lower); t1u <- as.numeric(he$H2_upper)
  h0  <- as.numeric(he$H0)                         # bootstrap naive (== gate naive)
  gap <- g$debiased$est - t1

  cat(sprintf("\n== %-22s : Tier-2 gate vs Tier-1 bootstrap (NB=%d) ==\n", tag, NB_BOOTS))
  cat(sprintf("  subgroup            : %s\n", paste(fs$sg.harm, collapse = " & ")))
  cat(sprintf("  naive HR            : gate %.3f | boot H0 %.3f  (scale check: %s)\n",
              g$naive$est, h0,
              if (is.finite(h0) && abs(log(g$naive$est) - log(h0)) < 1e-2) "OK" else "DIFF"))
  cat(sprintf("  Tier-2 de-biased HR : %.3f  [%.3f, %.3f]\n",
              g$debiased$est, g$debiased$lower, g$debiased$upper))
  cat(sprintf("  Tier-1 bootstrap H2 : %.3f  [%.3f, %.3f]\n", t1, t1l, t1u))
  cat(sprintf("  Layer-2 gap (Tier-2 debiased - bootstrap H2) = %+.3f\n", gap))
  invisible(c(t2 = g$debiased$est, t1 = t1, gap = gap))
}

if (RUN_BOOTSTRAP) {
  cat("\n", strrep("=", 60), "\n", sep = "")
  cat("Tier-1 bootstrap reference -- the Layer-2 (family-construction) gap\n")
  cat("(GRF bootstraps refit a forest per resample and dominate the runtime)\n")
  boot_compare(fs_dina, "DINA / dina")
  boot_compare(fs_deff, "DINA / effect")
  boot_compare(fs_gdr,  "GRF  / dr")
  boot_compare(fs_geff, "GRF  / effect")
  cat("\nReading: the gap is the family-construction optimism the single-fit\n")
  cat("Tier-2 gate cannot see. For the 'effect' modes (faithful Layer-1) the\n")
  cat("gap is the honest residual the full bootstrap still corrects; a large\n")
  cat("gap means the bootstrap stays the reference, a small gap means Tier-2\n")
  cat("already captures most of the optimism.\n")
}
