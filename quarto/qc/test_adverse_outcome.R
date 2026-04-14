# =============================================================================
# QC: adverse_outcome flip in make_effect_estimator
# =============================================================================
# Tests the Y -> 1-Y flip for beneficial binary outcomes.
# Run after devtools::load_all() or source the relevant files.

cat("\n")
cat("=============================================================\n")
cat("  QC: adverse_outcome in effect estimator\n")
cat("=============================================================\n\n")

library(speff2trial)
library(forestsearch)
# Source the files under test
#source("/mnt/user-data/outputs/glm_effect_estimators.R")

n_pass <- 0L
n_fail <- 0L

pass <- function(name) {
  n_pass <<- n_pass + 1L
  cat(sprintf("  PASS: %s\n", name))
}

fail <- function(name, msg = "") {
  n_fail <<- n_fail + 1L
  cat(sprintf("  FAIL: %s %s\n", name, msg))
}

check <- function(name, expr) {
  result <- tryCatch({
    force(expr)
    TRUE
  }, error = function(e) {
    fail(name, paste("--", e$message))
    FALSE
  })
  if (result) pass(name)
  invisible(result)
}


# =============================================================================
# Synthetic data
# =============================================================================

set.seed(42)
n <- 500
df <- data.frame(
  id    = 1:n,
  treat = rbinom(n, 1, 0.5)
)

# Beneficial outcome: treatment INCREASES improvement
p_improve_ctrl <- 0.40
p_improve_trt  <- 0.65
df$y_improve <- rbinom(n, 1,
  ifelse(df$treat == 1, p_improve_trt, p_improve_ctrl))

# Adverse outcome: treatment INCREASES failure (= 1 - improvement)
df$y_failure <- 1L - df$y_improve

cat("--- Synthetic data ---\n")
cat(sprintf("  N = %d, treat = %.1f%%\n", n, 100 * mean(df$treat)))
cat(sprintf("  P(improve|treat=1) = %.3f, P(improve|treat=0) = %.3f\n",
    mean(df$y_improve[df$treat == 1]),
    mean(df$y_improve[df$treat == 0])))
cat(sprintf("  P(failure|treat=1) = %.3f, P(failure|treat=0) = %.3f\n",
    mean(df$y_failure[df$treat == 1]),
    mean(df$y_failure[df$treat == 0])))


# =============================================================================
# GROUP 1: OR consistency -- OR(1-Y) should equal 1/OR(Y)
# =============================================================================

cat("\n--- Group 1: OR direction consistency ---\n")

fn_adverse <- make_effect_estimator(
  "binary", "treat", "y_improve",
  effect_measure = "OR", adverse_outcome = TRUE
)
fn_beneficial <- make_effect_estimator(
  "binary", "treat", "y_improve",
  effect_measure = "OR", adverse_outcome = FALSE
)

res_adv <- fn_adverse(df)
res_ben <- fn_beneficial(df)

or_raw   <- exp(res_adv$estimate)
or_flip  <- exp(res_ben$estimate)

cat(sprintf("  OR(Y, adverse=TRUE):  %.4f (treatment helps improvement)\n", or_raw))
cat(sprintf("  OR(Y, adverse=FALSE): %.4f (on failure scale)\n", or_flip))
cat(sprintf("  Product: %.4f (should be ~1.0)\n", or_raw * or_flip))

check("1.1 OR(adverse=TRUE) * OR(adverse=FALSE) ~ 1.0", {
  stopifnot(abs(or_raw * or_flip - 1.0) < 0.01)
})

check("1.2 OR(adverse=TRUE) > 1 when treatment helps", {
  stopifnot(or_raw > 1.0)
})

check("1.3 OR(adverse=FALSE) < 1 when treatment helps (on failure scale)", {
  stopifnot(or_flip < 1.0)
})

# Cross-check: OR on explicit y_failure should match adverse_outcome=FALSE
fn_explicit <- make_effect_estimator(
  "binary", "treat", "y_failure",
  effect_measure = "OR", adverse_outcome = TRUE
)
res_explicit <- fn_explicit(df)
or_explicit <- exp(res_explicit$estimate)

check("1.4 OR(adverse=FALSE on improve) = OR(adverse=TRUE on failure)", {
  stopifnot(abs(or_flip - or_explicit) < 0.001)
})


# =============================================================================
# GROUP 2: RD consistency -- RD(1-Y) should equal -RD(Y)
# =============================================================================

cat("\n--- Group 2: RD direction consistency ---\n")

fn_rd_adv <- make_effect_estimator(
  "binary", "treat", "y_improve",
  effect_measure = "RD", adverse_outcome = TRUE
)
fn_rd_ben <- make_effect_estimator(
  "binary", "treat", "y_improve",
  effect_measure = "RD", adverse_outcome = FALSE
)

res_rd_adv <- fn_rd_adv(df)
res_rd_ben <- fn_rd_ben(df)

cat(sprintf("  RD(adverse=TRUE):  %.4f (positive = treatment helps)\n",
    res_rd_adv$estimate))
cat(sprintf("  RD(adverse=FALSE): %.4f (positive = treatment hurts on failure)\n",
    res_rd_ben$estimate))
cat(sprintf("  Sum: %.6f (should be ~0.0)\n",
    res_rd_adv$estimate + res_rd_ben$estimate))

check("2.1 RD(adverse=TRUE) + RD(adverse=FALSE) ~ 0", {
  stopifnot(abs(res_rd_adv$estimate + res_rd_ben$estimate) < 0.001)
})

check("2.2 RD(adverse=TRUE) > 0 when treatment helps improvement", {
  stopifnot(res_rd_adv$estimate > 0)
})

check("2.3 RD(adverse=FALSE) < 0 when treatment helps (failure decreases)", {
  stopifnot(res_rd_ben$estimate < 0)
})


# =============================================================================
# GROUP 3: RR consistency
# =============================================================================

cat("\n--- Group 3: RR direction consistency ---\n")

fn_rr_adv <- make_effect_estimator(
  "binary", "treat", "y_improve",
  effect_measure = "RR", adverse_outcome = TRUE
)
fn_rr_ben <- make_effect_estimator(
  "binary", "treat", "y_improve",
  effect_measure = "RR", adverse_outcome = FALSE
)

res_rr_adv <- fn_rr_adv(df)
res_rr_ben <- fn_rr_ben(df)

rr_raw  <- exp(res_rr_adv$estimate)
rr_flip <- exp(res_rr_ben$estimate)

cat(sprintf("  RR(adverse=TRUE):  %.4f (treatment helps improvement)\n", rr_raw))
cat(sprintf("  RR(adverse=FALSE): %.4f (on failure scale)\n", rr_flip))

check("3.1 RR(adverse=TRUE) > 1 when treatment helps", {
  stopifnot(rr_raw > 1.0)
})

check("3.2 RR(adverse=FALSE) < 1 when treatment helps (failure decreases)", {
  stopifnot(rr_flip < 1.0)
})


# =============================================================================
# GROUP 4: adverse_outcome = TRUE is unchanged (backward compat)
# =============================================================================

cat("\n--- Group 4: backward compatibility ---\n")

# Default (adverse_outcome = TRUE) should produce identical results
# to an estimator without the parameter
fn_default <- make_effect_estimator(
  "binary", "treat", "y_failure",
  effect_measure = "OR"
)
fn_explicit_true <- make_effect_estimator(
  "binary", "treat", "y_failure",
  effect_measure = "OR", adverse_outcome = TRUE
)

res_def <- fn_default(df)
res_exp <- fn_explicit_true(df)

check("4.1 default adverse_outcome = TRUE produces same OR", {
  stopifnot(abs(res_def$estimate - res_exp$estimate) < 1e-10)
})

check("4.2 default adverse_outcome = TRUE produces same SE", {
  stopifnot(abs(res_def$se - res_exp$se) < 1e-10)
})

# Same for RD
fn_def_rd <- make_effect_estimator("binary", "treat", "y_failure",
  effect_measure = "RD")
fn_exp_rd <- make_effect_estimator("binary", "treat", "y_failure",
  effect_measure = "RD", adverse_outcome = TRUE)

check("4.3 default RD identical to adverse_outcome = TRUE", {
  res_a <- fn_def_rd(df)
  res_b <- fn_exp_rd(df)
  stopifnot(abs(res_a$estimate - res_b$estimate) < 1e-10)
})


# =============================================================================
# GROUP 5: Non-binary outcome types unaffected
# =============================================================================

cat("\n--- Group 5: non-binary types unaffected ---\n")

check("5.1 survival estimator ignores adverse_outcome", {
  # Should not error (adverse_outcome is not used for survival)
  fn <- make_effect_estimator(
    "survival", "treat", "time",
    event.name = "event", adverse_outcome = FALSE
  )
  stopifnot(is.function(fn))
})

check("5.2 continuous estimator ignores adverse_outcome", {
  df$y_cont <- rnorm(n)
  fn <- make_effect_estimator(
    "continuous", "treat", "y_cont",
    adverse_outcome = FALSE
  )
  res <- fn(df)
  stopifnot(!is.null(res$estimate))
})


# =============================================================================
# GROUP 6: Edge cases
# =============================================================================

cat("\n--- Group 6: Edge cases ---\n")

# 6.1 All Y = 1 (all improve, no failure)
check("6.1 all Y = 1: adverse_outcome=FALSE handles gracefully", {
  df_all1 <- data.frame(treat = c(0,0,0,1,1,1), y = c(1,1,1,1,1,1))
  fn <- make_effect_estimator("binary", "treat", "y",
    effect_measure = "RD", adverse_outcome = FALSE)
  res <- fn(df_all1)
  # After flip, all Y = 0. RD should be 0.
  stopifnot(abs(res$estimate) < 1e-10)
})

# 6.2 All Y = 0 (all fail, no improvement)
check("6.2 all Y = 0: adverse_outcome=FALSE handles gracefully", {
  df_all0 <- data.frame(treat = c(0,0,0,1,1,1), y = c(0,0,0,0,0,0))
  fn <- make_effect_estimator("binary", "treat", "y",
    effect_measure = "RD", adverse_outcome = FALSE)
  res <- fn(df_all0)
  # After flip, all Y = 1. RD should be 0.
  stopifnot(abs(res$estimate) < 1e-10)
})

# 6.3 Small sample (n = 10)
check("6.3 small sample (n=10) works with flip", {
  df_small <- data.frame(
    treat = c(0,0,0,0,0, 1,1,1,1,1),
    y     = c(0,1,0,0,1, 1,1,1,0,1)
  )
  fn <- make_effect_estimator("binary", "treat", "y",
    effect_measure = "OR", adverse_outcome = FALSE)
  res <- fn(df_small)
  stopifnot(!is.na(res$estimate), res$converged)
})

# 6.4 Perfect separation (all ctrl fail, all trt succeed)
check("6.4 perfect separation handled", {
  df_sep <- data.frame(
    treat = c(0,0,0,0,0, 1,1,1,1,1),
    y     = c(0,0,0,0,0, 1,1,1,1,1)
  )
  fn <- make_effect_estimator("binary", "treat", "y",
    effect_measure = "OR", adverse_outcome = FALSE)
  res <- fn(df_sep)
  # After flip: ctrl all 1, trt all 0 -> perfect separation the other way
  # Should still return a result (may use fallback)
  stopifnot(!is.null(res$estimate))
})

# 6.5 Repeated calls on same estimator (closure stability)
check("6.5 closure is stable across repeated calls", {
  fn <- make_effect_estimator("binary", "treat", "y_improve",
    effect_measure = "OR", adverse_outcome = FALSE)
  r1 <- fn(df)
  r2 <- fn(df)
  stopifnot(abs(r1$estimate - r2$estimate) < 1e-10)
})

# 6.6 Original data is NOT mutated by flip
check("6.6 original data not mutated", {
  df_copy <- data.frame(treat = df$treat, y_improve = df$y_improve)
  original_y <- df_copy$y_improve
  fn <- make_effect_estimator("binary", "treat", "y_improve",
    effect_measure = "OR", adverse_outcome = FALSE)
  res <- fn(df_copy)
  stopifnot(identical(df_copy$y_improve, original_y))
})


# =============================================================================
# GROUP 7: Harm vs benefit search scenarios
# =============================================================================

cat("\n--- Group 7: Four-scenario framework ---\n")

# Create data with known direction:
# treat = ddI (switched), control = combo
# In Q (z1==1): ddI HURTS improvement (combo benefits)
# In Qc: no differential effect

set.seed(123)
n2 <- 1000
df2 <- data.frame(
  id = 1:n2,
  z1 = rbinom(n2, 1, 0.35),
  treat = rbinom(n2, 1, 0.5)
)

# Improvement probabilities
p_ctrl <- 0.60
p_trt_qc <- 0.60  # no effect in Qc
p_trt_q  <- 0.35  # ddI hurts Q (less improvement)

df2$p <- ifelse(df2$treat == 1,
  ifelse(df2$z1 == 1, p_trt_q, p_trt_qc),
  p_ctrl
)
df2$y_improve <- rbinom(n2, 1, df2$p)
df2$y_failure <- 1L - df2$y_improve

# (A) adverse_outcome=TRUE on failure outcome, no switch needed for harm search
check("7.1 Scenario A: adverse=TRUE on failure, harm search", {
  fn <- make_effect_estimator("binary", "treat", "y_failure",
    effect_measure = "OR", adverse_outcome = TRUE)
  res_q <- fn(df2[df2$z1 == 1, ])
  or_q <- exp(res_q$estimate)
  # ddI increases failure for Q -> OR(failure) > 1
  stopifnot(or_q > 1.0)
})

# (B) adverse_outcome=TRUE on failure outcome, benefit search (switch treat)
check("7.2 Scenario B: adverse=TRUE on failure, benefit search (switched)", {
  df2_switch <- df2
  df2_switch$treat <- 1L - df2_switch$treat  # switch roles
  fn <- make_effect_estimator("binary", "treat", "y_failure",
    effect_measure = "OR", adverse_outcome = TRUE)
  res_q <- fn(df2_switch[df2_switch$z1 == 1, ])
  or_q <- exp(res_q$estimate)
  # After switching, "treat" is combo. combo DECREASES failure for Q
  # -> OR(failure, switched) < 1
  stopifnot(or_q < 1.0)
})

# (C) adverse_outcome=FALSE on improvement, benefit search (switch treat)
check("7.3 Scenario C: adverse=FALSE on improvement, benefit search (switched)", {
  df2_switch <- df2
  df2_switch$treat <- 1L - df2_switch$treat
  fn <- make_effect_estimator("binary", "treat", "y_improve",
    effect_measure = "OR", adverse_outcome = FALSE)
  res_q <- fn(df2_switch[df2_switch$z1 == 1, ])
  or_q <- exp(res_q$estimate)
  # adverse_outcome=FALSE flips to failure internally
  # After switching, "treat" is combo. On failure scale: combo decreases
  # failure -> OR(failure, switched) < 1
  # Same as scenario B
  stopifnot(or_q < 1.0)
})

# (D) adverse_outcome=FALSE on improvement, harm search (no switch)
check("7.4 Scenario D: adverse=FALSE on improvement, harm search", {
  fn <- make_effect_estimator("binary", "treat", "y_improve",
    effect_measure = "OR", adverse_outcome = FALSE)
  res_q <- fn(df2[df2$z1 == 1, ])
  or_q <- exp(res_q$estimate)
  # On failure scale: ddI increases failure for Q -> OR > 1
  stopifnot(or_q > 1.0)
})

# (E) Verify C and A produce same OR in Q (both on failure scale)
check("7.5 Scenario C OR = Scenario A OR (same failure scale)", {
  fn_A <- make_effect_estimator("binary", "treat", "y_failure",
    effect_measure = "OR", adverse_outcome = TRUE)
  fn_C <- make_effect_estimator("binary", "treat", "y_improve",
    effect_measure = "OR", adverse_outcome = FALSE)
  or_A <- exp(fn_A(df2[df2$z1 == 1, ])$estimate)
  or_C <- exp(fn_C(df2[df2$z1 == 1, ])$estimate)
  # Both compute OR on failure scale -> should be identical
  stopifnot(abs(or_A - or_C) < 0.001)
})


# =============================================================================
# GROUP 8: ACTG175 realistic test
# =============================================================================

cat("\n--- Group 8: ACTG175 realistic data ---\n")

actg_df <- subset(ACTG175, arms %in% c(1, 3))
actg_df$treat_orig <- ifelse(actg_df$arms == 1, 1L, 0L)
actg_df$treat <- 1L - actg_df$treat_orig  # switched
actg_df$y_improve <- as.integer(actg_df$cd420 > actg_df$cd40)
actg_df <- actg_df[!is.na(actg_df$cd420), ]
actg_df$z1 <- as.factor(ifelse(actg_df$age > 34, 1L, 0L))
actg_df$z2 <- as.factor(ifelse(actg_df$preanti <= 744.5, 1L, 0L))
actg_df$flag_Q <- as.integer(actg_df$z1 == 1 & actg_df$z2 == 1)

# With adverse_outcome=FALSE, OR should be on failure scale
fn_flip <- make_effect_estimator("binary", "treat", "y_improve",
  effect_measure = "OR", adverse_outcome = FALSE)

# Q subgroup: combo benefits -> ddI hurts -> OR(failure, Q, switched) > 1
res_Q <- fn_flip(actg_df[actg_df$flag_Q == 1, ])
or_Q <- exp(res_Q$estimate)

# Qc: no strong effect expected
res_Qc <- fn_flip(actg_df[actg_df$flag_Q == 0, ])
or_Qc <- exp(res_Qc$estimate)

cat(sprintf("  ACTG175 (switched, adverse_outcome=FALSE):\n"))
cat(sprintf("    OR(Q, failure scale):  %.3f\n", or_Q))
cat(sprintf("    OR(Qc, failure scale): %.3f\n", or_Qc))

# Cross-check with explicit failure outcome
actg_df$y_failure <- 1L - actg_df$y_improve
fn_explicit <- make_effect_estimator("binary", "treat", "y_failure",
  effect_measure = "OR", adverse_outcome = TRUE)
or_Q_explicit <- exp(fn_explicit(actg_df[actg_df$flag_Q == 1, ])$estimate)

check("8.1 ACTG175: flip OR matches explicit failure OR", {
  stopifnot(abs(or_Q - or_Q_explicit) < 0.001)
})

check("8.2 ACTG175: OR(Q) > 1 (ddI hurts Q on failure scale)", {
  # In Q under switched treatment: ddI has lower improvement = higher failure
  # OR(failure) should be > 1
  stopifnot(or_Q > 1.0)
})

check("8.3 ACTG175: OR(Q) > OR(Qc) (Q is the harm subgroup)", {
  stopifnot(or_Q > or_Qc)
})

# ITT check
res_ITT <- fn_flip(actg_df)
or_ITT <- exp(res_ITT$estimate)
cat(sprintf("    OR(ITT, failure scale): %.3f\n", or_ITT))

check("8.4 ACTG175: ITT OR reasonable", {
  stopifnot(or_ITT > 0.5, or_ITT < 3.0)
})


# =============================================================================
# SUMMARY
# =============================================================================

cat("\n=============================================================\n")
cat(sprintf("  %d PASSED, %d FAILED, %d total\n",
    n_pass, n_fail, n_pass + n_fail))
if (n_fail == 0) {
  cat("  All tests passed.\n")
} else {
  cat("  *** FAILURES DETECTED ***\n")
}
cat("=============================================================\n\n")

if (n_fail > 0) stop(sprintf("%d test(s) failed", n_fail))
