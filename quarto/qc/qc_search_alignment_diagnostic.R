#!/usr/bin/env Rscript
# =============================================================================
# QC: Search Configuration Diagnostic -- All Scenarios
#
# Tests interpret_search_config() and verifies alignment between
# the diagnostic output and actual ForestSearch detection for every
# scenario in the Outcome Orientation Guide.
#
# Requires: forestsearch installed (devtools::load_all() or CRAN)
# =============================================================================

library(forestsearch)
library(survival)

# # Overlay fixed files (if running from the outputs directory)
# for (f in c("interpret_search_config.R", "forestsearch_main.R")) {
#   if (file.exists(f)) source(f)
# }

cat("=============================================================\n")
cat("  QC: Search Alignment Diagnostic -- All Scenarios\n")
cat("=============================================================\n\n")

pass <- 0L; fail <- 0L

check <- function(label, cond, msg = "") {
  if (isTRUE(cond)) {
    pass <<- pass + 1L; cat(sprintf("  [PASS] %s\n", label))
  } else {
    fail <<- fail + 1L; cat(sprintf("  [FAIL] %s %s\n", label, msg))
  }
}


# --- Shared synthetic data -----------------------------------------------
set.seed(42)
N <- 800
syn <- data.frame(
  id    = seq_len(N),
  treat = rep(0:1, each = N / 2),
  z1    = as.factor(sample(0:1, N, TRUE)),
  z2    = as.factor(sample(0:1, N, TRUE)),
  z3    = as.factor(sample(0:1, N, TRUE))
)
in_Q <- syn$z1 == 1 & syn$z2 == 1
cat(sprintf("Q prevalence: %.1f%%\n\n", 100 * mean(in_Q)))


# --- Helper: run FS and return found/not-found ---------------------------
run_fs <- function(df, outcome, event, treat, otype, emeasure,
                   threshold, consist, ...) {
  result <- suppressWarnings(forestsearch(
    df.analysis = df,
    confounders.name = c("z1", "z2", "z3"),
    outcome.name = outcome,
    event.name = event,
    treat.name = treat,
    id.name = "id",
    outcome_type = otype,
    effect_measure = emeasure,
    effect.threshold = threshold,
    consistency.threshold = consist,
    pconsistency.threshold = 0.80,
    fs.splits = 200L,
    n.min = 40L,
    d0.min = 8L,
    d1.min = 8L,
    maxk = 2L,
    use_lasso = FALSE,
    use_grf = TRUE,
    is.RCT = TRUE,
    seedit = 42L,
    details = TRUE,            # triggers alignment diagnostic
    quiet = FALSE,
    parallel_args = list(plan = "sequential", workers = 1,
                         show_message = FALSE),
    ...
  ))
  found <- !is.null(result) && !is.null(result$grp.consistency)
  vars <- character(0)
  if (found && !is.null(result$grp.consistency$sg.harm)) {
    # sg.harm is a character vector like c("{z2}", "{z1}")
    # Strip braces to get variable names
    vars <- gsub("[{}]", "", result$grp.consistency$sg.harm)
  }
  list(found = found, vars = vars)
}


# ===========================================================================
cat("--- 1. Diagnostic Output -- Binary Scenarios --------------\n")
# ===========================================================================

# 1a: Binary adverse + harm (direct)
diag_1a <- interpret_search_config(
  "binary", "OR", adverse_outcome = TRUE,
  effect_threshold = log(1.5), consistency_threshold = log(1.2),
  use_lasso = FALSE, use_grf = TRUE,
  outcome.name = "y_event", event.name = "y_event",
  treat.name = "treat", quiet = FALSE
)
check("Diag 1a: binary harm/adverse -- no warnings",
      diag_1a$n_warnings == 0)

# 1b: Binary positive + benefit (switched, adverse_outcome=TRUE default)
diag_1b <- interpret_search_config(
  "binary", "OR", adverse_outcome = TRUE,
  effect_threshold = log(1.667), consistency_threshold = log(1.25),
  use_lasso = FALSE, use_grf = TRUE,
  outcome.name = "y_improve", event.name = "y_improve",
  treat.name = "treat_sw", quiet = FALSE
)
check("Diag 1b: binary benefit/positive (switched) -- no warnings",
      diag_1b$n_warnings == 0)

# 1c: Binary positive + adverse_outcome=FALSE (the broken pattern)
diag_1c <- interpret_search_config(
  "binary", "OR", adverse_outcome = FALSE,
  effect_threshold = log(1.667), consistency_threshold = log(1.25),
  use_lasso = FALSE, use_grf = TRUE,
  outcome.name = "y_improve", event.name = "y_improve",
  treat.name = "treat_sw", quiet = FALSE
)
check("Diag 1c: binary benefit + adverse_outcome=FALSE -- WARNING issued",
      diag_1c$n_warnings > 0)
check("Diag 1c: warning mentions estimator flip",
      any(grepl("flip|COMPLEMENT|1-Y", diag_1c$warnings)))

cat("\n")


# ===========================================================================
cat("--- 2. Diagnostic Output -- Continuous Scenarios ----------\n")
# ===========================================================================

# 2a: Continuous harm/adverse (direct)
diag_2a <- interpret_search_config(
  "continuous", "MD", adverse_outcome = FALSE,
  effect_threshold = 30, consistency_threshold = 10,
  use_lasso = TRUE, use_grf = TRUE,
  outcome.name = "cd4_change", event.name = "cd4_change",
  treat.name = "treat", quiet = FALSE
)
check("Diag 2a: continuous default (adverse_outcome=FALSE) -- no warnings",
      diag_2a$n_warnings == 0)

# 2b: Continuous with adverse_outcome=TRUE (e.g., pain score)
diag_2b <- interpret_search_config(
  "continuous", "MD", adverse_outcome = TRUE,
  effect_threshold = 5, consistency_threshold = 2,
  use_lasso = TRUE, use_grf = TRUE,
  outcome.name = "pain_score", event.name = "pain_score",
  treat.name = "treat", quiet = FALSE
)
check("Diag 2b: continuous adverse_outcome=TRUE -- alignment note issued",
      diag_2b$n_warnings > 0)
check("Diag 2b: note mentions GRF flip",
      any(grepl("GRF|1-pain_score", diag_2b$warnings)))

cat("\n")


# ===========================================================================
cat("--- 3. Diagnostic Output -- Count Scenarios ----------------\n")
# ===========================================================================

diag_3a <- interpret_search_config(
  "count", "IRR", adverse_outcome = TRUE,
  effect_threshold = log(1.25), consistency_threshold = log(1.0),
  use_lasso = TRUE, use_grf = TRUE,
  outcome.name = "status", event.name = "status",
  treat.name = "hormon", offset.name = "time_months",
  quiet = FALSE
)
check("Diag 3a: count/IRR with offset -- no warnings",
      diag_3a$n_warnings == 0)
check("Diag 3a: mentions offset",
      grepl("time_months", diag_3a$search_finds))

cat("\n")


# ===========================================================================
cat("--- 4. Diagnostic Output -- Survival Scenarios ------------\n")
# ===========================================================================

# 4a: Survival harm/adverse (standard)
diag_4a <- interpret_search_config(
  "survival", "HR", adverse_outcome = TRUE,
  effect_threshold = log(1.25), consistency_threshold = log(1.0),
  use_lasso = TRUE, use_grf = TRUE,
  outcome.name = "rfstime", event.name = "status",
  treat.name = "hormon", quiet = FALSE
)
check("Diag 4a: survival harm/adverse -- no warnings",
      diag_4a$n_warnings == 0)
check("Diag 4a: mentions causal_survival_forest",
      grepl("causal_survival_forest", diag_4a$grf_note))

# 4b: Survival with GRF only (no LASSO) -- should warn
diag_4b <- interpret_search_config(
  "survival", "HR", adverse_outcome = TRUE,
  effect_threshold = log(1.25), consistency_threshold = log(1.0),
  use_lasso = FALSE, use_grf = TRUE,
  outcome.name = "rfstime", event.name = "status",
  treat.name = "hormon", quiet = FALSE
)
check("Diag 4b: survival GRF-only -- alignment note issued",
      diag_4b$n_warnings > 0)

cat("\n")


# ===========================================================================
cat("--- 5. End-to-End: Harm + Adverse Binary -----------------\n")
# ===========================================================================

# Strong signal: treatment increases P(event) for Q
p_event <- ifelse(syn$treat == 1 & in_Q, 0.80,
            ifelse(syn$treat == 1, 0.35, 0.35))
syn$y_event <- rbinom(N, 1, p_event)

r5 <- run_fs(syn, "y_event", "y_event", "treat", "binary", "OR", 2.0, 1.5)
check("E2E harm+adverse binary: detected", r5$found)
if (r5$found)
  check("E2E harm+adverse binary: z1 in subgroup", "z1" %in% r5$vars)

cat("\n")


# ===========================================================================
cat("--- 6. E2E: Benefit + Positive Binary (Switched) ---------\n")
# ===========================================================================

# Combo (treat_orig=1) helps Q improve more
p_improve <- ifelse(syn$treat == 0 & in_Q, 0.80,
              ifelse(syn$treat == 0, 0.45,
               ifelse(in_Q, 0.40, 0.45)))
syn$y_improve <- rbinom(N, 1, p_improve)
syn$treat_sw <- 1L - syn$treat

r6 <- run_fs(syn, "y_improve", "y_improve", "treat_sw", "binary", "OR", 1.5, 1.2)
check("E2E benefit+positive binary: detected", r6$found)
if (r6$found)
  check("E2E benefit+positive binary: z1 in subgroup", "z1" %in% r6$vars)

cat("\n")


# ===========================================================================
cat("--- 7. E2E: Harm + Positive Continuous (Negate Y) --------\n")
# ===========================================================================

# Treatment reduces beneficial outcome for Q
syn$y_good <- 50 + 5 * syn$treat + rnorm(N, sd = 10)
syn$y_good[in_Q & syn$treat == 1] <- syn$y_good[in_Q & syn$treat == 1] - 30

# Negate
syn$y_bad <- -syn$y_good

r7 <- run_fs(syn, "y_bad", "y_bad", "treat", "continuous", "MD", 15, 5)
check("E2E harm+positive cont (negated): detected", r7$found)

cat("\n")


# ===========================================================================
cat("--- 8. E2E: Benefit + Positive Continuous (Switched) -----\n")
# ===========================================================================

# Combo helps Q on continuous outcome
syn$y_cd4 <- 30 + 3 * (1 - syn$treat) + rnorm(N, sd = 12)
syn$y_cd4[in_Q & syn$treat == 0] <- syn$y_cd4[in_Q & syn$treat == 0] + 25

r8 <- run_fs(syn, "y_cd4", "y_cd4", "treat_sw", "continuous", "MD", 15, 5)
check("E2E benefit+positive cont (switched): detected", r8$found)
if (r8$found)
  check("E2E benefit+positive cont: z1 in subgroup", "z1" %in% r8$vars)

cat("\n")


# ===========================================================================
cat("--- 9. E2E: Harm + Adverse Count (Direct) ----------------\n")
# ===========================================================================

# Treatment increases event rate for Q
lambda <- exp(0.5 + 0.2 * syn$treat)
lambda[in_Q & syn$treat == 1] <- lambda[in_Q & syn$treat == 1] * 3
syn$y_ct <- rpois(N, lambda)
syn$fu <- rep(1, N)

r9 <- run_fs(syn, "y_ct", "y_ct", "treat", "count", "IRR", 2.0, 1.5,
             offset.name = "fu")
check("E2E harm+adverse count: detected", r9$found)

cat("\n")


# ===========================================================================
cat("--- 10. E2E: Harm + Adverse Survival (Direct) ------------\n")
# ===========================================================================

# Generate survival data with known subgroup harm
set.seed(42)
syn_surv <- syn[, c("id", "treat", "z1", "z2", "z3")]

# Baseline hazard ≈ exponential with rate 0.05
# Treatment increases hazard for Q (harm)
base_rate <- 0.05
rate <- base_rate * exp(0.3 * syn_surv$treat)
rate[in_Q & syn_surv$treat == 1] <- rate[in_Q & syn_surv$treat == 1] * 2.5

syn_surv$time <- rexp(N, rate)
syn_surv$event <- 1L  # no censoring for simplicity
# Truncate at 50
syn_surv$event[syn_surv$time > 50] <- 0L
syn_surv$time <- pmin(syn_surv$time, 50)

# Verify
cox_Q <- coxph(Surv(time, event) ~ treat,
               data = syn_surv[in_Q, ])
cox_Qc <- coxph(Surv(time, event) ~ treat,
                data = syn_surv[!in_Q, ])
cat(sprintf("  True HR(Q): %.2f, HR(Qc): %.2f\n",
    exp(coef(cox_Q)), exp(coef(cox_Qc))))

r10 <- suppressWarnings(forestsearch(
  df.analysis = syn_surv,
  confounders.name = c("z1", "z2", "z3"),
  outcome.name = "time",
  event.name = "event",
  treat.name = "treat",
  id.name = "id",
  hr.threshold = 1.5,
  hr.consistency = 1.2,
  pconsistency.threshold = 0.80,
  fs.splits = 200L,
  n.min = 40L,
  d0.min = 8L,
  d1.min = 8L,
  maxk = 2L,
  use_lasso = TRUE,
  use_grf = TRUE,
  is.RCT = TRUE,
  seedit = 42L,
  details = TRUE,
  quiet = FALSE,
  parallel_args = list(plan = "sequential", workers = 1,
                       show_message = FALSE)
))

r10_found <- !is.null(r10) && !is.null(r10$grp.consistency)
check("E2E survival harm+adverse: detected", r10_found)
if (r10_found) {
  sg10 <- gsub("[{}]", "", r10$grp.consistency$sg.harm)
  check("E2E survival: z1 in subgroup", "z1" %in% sg10)
}

cat("\n")


# ===========================================================================
cat("--- 11. E2E: Benefit + Adverse Survival (Switched) -------\n")
# ===========================================================================

# Same data but switch treatment: find where original treatment HELPS
syn_surv$treat_sw <- 1L - syn_surv$treat

r11 <- suppressWarnings(forestsearch(
  df.analysis = syn_surv,
  confounders.name = c("z1", "z2", "z3"),
  outcome.name = "time",
  event.name = "event",
  treat.name = "treat_sw",       # switched!
  id.name = "id",
  hr.threshold = 1.5,
  hr.consistency = 1.2,
  pconsistency.threshold = 0.80,
  fs.splits = 200L,
  n.min = 40L,
  d0.min = 8L,
  d1.min = 8L,
  maxk = 2L,
  use_lasso = TRUE,
  use_grf = TRUE,
  is.RCT = TRUE,
  seedit = 42L,
  quiet = TRUE,
  parallel_args = list(plan = "sequential", workers = 1,
                       show_message = FALSE)
))

# Under switched treatment, Q is where "control" (now treat=1) hurts
# → equivalent to original treatment helping Qc
# This should NOT find the same Q (it should find Qc or nothing)
r11_found <- !is.null(r11) && !is.null(r11$grp.consistency)
cat(sprintf("  Switched survival: found=%s\n", r11_found))
if (r11_found) {
  sg11 <- gsub("[{}]", "", r11$grp.consistency$sg.harm)
  cat(sprintf("  Variables: %s\n", paste(sg11, collapse = ", ")))
}
# Note: detection here depends on the specific data; the key test is
# that the diagnostic correctly describes the switched search.

cat("\n")


# ===========================================================================
cat("--- 12. E2E: Null -- No Subgroup (All Types) --------------\n")
# ===========================================================================

# Binary null
syn$y_null_b <- rbinom(N, 1, 0.5)
r12b <- run_fs(syn, "y_null_b", "y_null_b", "treat", "binary", "OR", 1.5, 1.2)
cat(sprintf("  Binary null: found=%s (should usually be FALSE)\n", r12b$found))

# Continuous null
syn$y_null_c <- rnorm(N, 50, 10)
r12c <- run_fs(syn, "y_null_c", "y_null_c", "treat", "continuous", "MD", 10, 5)
cat(sprintf("  Continuous null: found=%s (should usually be FALSE)\n", r12c$found))

cat("\n")


# ===========================================================================
cat("=============================================================\n")
cat(sprintf("  RESULTS: %d passed, %d failed\n", pass, fail))
cat("=============================================================\n")
if (fail > 0) {
  cat("\n*** FAILURES -- review above ***\n")
} else {
  cat("\n=== All checks passed ===\n")
}
