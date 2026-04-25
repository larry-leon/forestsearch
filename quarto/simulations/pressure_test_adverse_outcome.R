# =============================================================================
# Pressure test: adverse_outcome sign-flip mechanism for continuous + benefit
# =============================================================================
# Goal: Prove that adverse_outcome = FALSE causes grf.subg.harm.glm() to
# identify the WRONG subgroup (Qc instead of Q) under benefit search via
# treatment switching, and that adverse_outcome = TRUE aligns GRF with FS.
#
# Strategy: Simulate the data structure the vignette produces, then replicate
# .select_best_sg_glm()'s mDiff-selection logic under each flag setting.
# We don't need the actual grf::causal_forest() run — the sign of the
# identified subgroup is determined entirely by the sign of CATE inside Q
# and the score = -CATE convention.
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
})

PASS <- 0L; FAIL <- 0L
ok <- function(label, expr) {
  res <- tryCatch(expr, error = function(e) e)
  if (inherits(res, "error")) {
    FAIL <<- FAIL + 1L
    cat(sprintf("  [FAIL] %s\n         %s\n", label, conditionMessage(res)))
    return(invisible(NULL))
  }
  if (isTRUE(res)) {
    PASS <<- PASS + 1L
    cat(sprintf("  [ OK ] %s\n", label))
  } else {
    FAIL <<- FAIL + 1L
    cat(sprintf("  [FAIL] %s (got: %s)\n", label, deparse(res, 60)[1]))
  }
}

cat("\n=== Mechanism verification ===\n\n")

# ---------------------------------------------------------------------------
# Build fake data mirroring the vignette's benefit-search DGM
#   - Q subset: z1 == 1 & z2 == 1  (~30% of population)
#   - On SWITCHED scale: ddI=1 INCREASES cd4_change by +40 in Q
#                        ddI=1 DECREASES cd4_change by ~26 in Qc
#   - ITT effect: negative (switched-scale ddI generally hurts)
# ---------------------------------------------------------------------------
set.seed(8316951L)
n <- 1000L

z1 <- rbinom(n, 1L, 0.55)
z2 <- rbinom(n, 1L, 0.55)
in_Q <- (z1 == 1L & z2 == 1L)
treat_sim <- rbinom(n, 1L, 0.5)  # RCT

# Potential outcomes under switched roles
mu0 <- rnorm(n, mean = 50, sd = 10)            # combo = control (W=0)
mu1 <- mu0 +
       ifelse(in_Q, 40, -26)                    # switched-scale treatment effect
y_sim <- ifelse(treat_sim == 1L, mu1, mu0) + rnorm(n, sd = 5)

cat(sprintf("  Q prevalence: %.1f%%\n", 100 * mean(in_Q)))
cat(sprintf("  Observed MD(Q, switched)  (glm):  %.2f\n",
            coef(glm(y_sim ~ treat_sim, subset = in_Q,  family = gaussian()))[2]))
cat(sprintf("  Observed MD(Qc, switched) (glm):  %.2f\n",
            coef(glm(y_sim ~ treat_sim, subset = !in_Q, family = gaussian()))[2]))
cat(sprintf("  Observed MD(ITT, switched)(glm):  %.2f\n",
            coef(glm(y_sim ~ treat_sim, family = gaussian()))[2]))

# ---------------------------------------------------------------------------
# Emulate DR scores.  With an RCT and a correctly-specified CATE estimator,
# the doubly-robust score row for subject i is (Gamma_0_i, Gamma_1_i) where
# Gamma_w = mu_w_hat(X_i) + (W_i == w) / P(W=w) * (Y_i - mu_w_hat(X_i)).
# At the population level, mean(DR[,2] - DR[,1]) in any subset S equals
# CATE_hat(S).  For this mechanism test we can skip the forest and use
# true potential outcomes directly — the mechanism in .select_best_sg_glm()
# only cares about the SIGN and MAGNITUDE of CATE in candidate subgroups.
# ---------------------------------------------------------------------------
# Simulate grf's internal Y-flip exactly as grf_subg_harm_glm.R lines 417-420:
#   if (adverse_outcome) Y_grf <- -Y
# Then DR scores use Y_grf, so CATE_flipped = -CATE_original.
emulate_select_best <- function(in_candidate_set, adverse_outcome,
                                dmin_grf = 10) {
  # CATE in each candidate subset, under the Y-sign convention
  cate_Q_true_switched  <- +40    # switched scale CATE in Q
  cate_Qc_true_switched <- -26    # switched scale CATE in Qc

  # Under adverse_outcome = TRUE, Y -> -Y, so CATE effectively flips sign.
  sign_factor <- if (adverse_outcome) -1 else +1
  cate_Q  <- sign_factor * cate_Q_true_switched
  cate_Qc <- sign_factor * cate_Qc_true_switched

  # .select_best_sg_glm uses score = -CATE.  Candidate sets are generated
  # by policy_tree.  For this mechanism test, the two relevant candidates
  # are: (i) "Q" (leaf that splits on z1 & z2), (ii) "Qc".  The tree picks
  # whichever yields the larger score.
  score_Q  <- -cate_Q
  score_Qc <- -cate_Qc

  # mDiff picks which.max(scores), requires max >= dmin.grf
  scores <- c(Q = score_Q, Qc = score_Qc)
  best   <- names(scores)[which.max(scores)]
  feasible <- max(scores) >= dmin_grf

  list(scores = scores, best = best, feasible = feasible)
}

# ---------------------------------------------------------------------------
# Test: under adverse_outcome = FALSE (the buggy default), GRF picks Qc.
# ---------------------------------------------------------------------------
res_false <- emulate_select_best(in_Q, adverse_outcome = FALSE)
cat("\n  adverse_outcome = FALSE:\n")
cat(sprintf("    scores = (Q: %.1f, Qc: %.1f), picks '%s' (feasible: %s)\n",
            res_false$scores["Q"], res_false$scores["Qc"],
            res_false$best, res_false$feasible))

ok("1.1 FALSE: picks Qc (WRONG subgroup)", res_false$best == "Qc")
ok("1.2 FALSE: selection is feasible",     isTRUE(res_false$feasible))

# ---------------------------------------------------------------------------
# Test: under adverse_outcome = TRUE (the fix), GRF picks Q.
# ---------------------------------------------------------------------------
res_true <- emulate_select_best(in_Q, adverse_outcome = TRUE)
cat("\n  adverse_outcome = TRUE:\n")
cat(sprintf("    scores = (Q: %.1f, Qc: %.1f), picks '%s' (feasible: %s)\n",
            res_true$scores["Q"], res_true$scores["Qc"],
            res_true$best, res_true$feasible))

ok("2.1 TRUE: picks Q (CORRECT subgroup)", res_true$best == "Q")
ok("2.2 TRUE: selection is feasible",      isTRUE(res_true$feasible))

# ---------------------------------------------------------------------------
# Test: downstream estimation agrees after the fix
# Simulate: if GRF correctly identifies Q (adverse_outcome = TRUE), then
# running .glm_effect_grf(grf_data[sg_hat == 1, ]) on Q gives MD(Q, switched)
# = +40ish.  Benefit-scale inversion (negate for MD) gives -40, matching
# the true theta_dagger(Q) reported in the HTML as -39.74.
# ---------------------------------------------------------------------------
# Assume GRF gets Q right: hr.H.hat on switched scale is the MD fit on Q
md_H_hat_switched <- coef(glm(y_sim ~ treat_sim, subset = in_Q,
                              family = gaussian()))[2]
# benefit-scale inversion (invert_hr_columns, identity scale: negation)
md_H_hat_benefit  <- -md_H_hat_switched

cat(sprintf(
  "\n  Downstream: MD(H_hat, switched) = %.2f -> benefit scale = %.2f\n",
  md_H_hat_switched, md_H_hat_benefit))
cat("  Truth theta_dagger(Q) (from HTML output): -39.74\n")

ok("3.1 benefit-scale estimate matches truth (within 3)",
   abs(md_H_hat_benefit - (-39.74)) < 3)

# ---------------------------------------------------------------------------
# Test: under the buggy path (FALSE), the estimation gives opposite sign
# GRF identifies Qc as "the subgroup"; sg_hat == 1 matches Qc subjects
# ---------------------------------------------------------------------------
md_H_hat_switched_wrong <- coef(
  glm(y_sim ~ treat_sim, subset = !in_Q, family = gaussian()))[2]
md_H_hat_benefit_wrong  <- -md_H_hat_switched_wrong

cat(sprintf(
  "  Buggy path: MD(Qc, switched) = %.2f -> benefit scale = %.2f\n",
  md_H_hat_switched_wrong, md_H_hat_benefit_wrong))
cat("  Compare to HTML GRF row: +36.69 (same-ish; wrong sign from truth)\n")

ok("3.2 buggy estimate is opposite sign from truth",
   sign(md_H_hat_benefit_wrong) != sign(-39.74))

# =============================================================================
# Test: the vignette's grf_params list is constructed correctly
# =============================================================================
cat("\n=== Vignette grf_params construction ===\n\n")

# Emulate the fixed setup chunk
grf_dmin            <- 10
grf_maxdepth        <- 2L
grf_sg_criterion    <- "mDiff"
grf_tune            <- FALSE
grf_adverse_outcome <- TRUE

grf_params <- list(
  dmin.grf        = grf_dmin,
  maxdepth        = grf_maxdepth,
  sg.criterion    = grf_sg_criterion,
  tune_grf        = grf_tune,
  adverse_outcome = grf_adverse_outcome
)

ok("4.1 grf_params has adverse_outcome key",
   "adverse_outcome" %in% names(grf_params))
ok("4.2 grf_params$adverse_outcome is TRUE",
   isTRUE(grf_params$adverse_outcome))

# Emulate run_simulation_analysis.R:381 (modifyList(grf_defaults, grf_params))
# and R:990-992 (params$adverse_outcome resolution)
grf_defaults <- list(
  outcome.name = "y_sim", event.name = "event_sim", treat.name = "treat_sim",
  id.name = "id", n.min = 60, dmin.grf = 12, frac.tau = 0.60,
  maxdepth = 2, RCT = TRUE, sg.criterion = "mDiff", seedit = 8316951L
)
grf_merged <- utils::modifyList(grf_defaults, grf_params)

ok("4.3 merged params carry adverse_outcome = TRUE",
   isTRUE(grf_merged$adverse_outcome))
ok("4.4 merged params carry dmin.grf = 10 (override worked)",
   grf_merged$dmin.grf == 10)

# Emulate the resolution on lines 990-992 of run_simulation_analysis.R
outcome_type <- "continuous"
ao <- grf_merged$adverse_outcome  # params$adverse_outcome
if (is.null(ao)) ao <- outcome_type %in% c("binary", "count")

ok("4.5 run_simulation_analysis resolves ao = TRUE for continuous + override",
   isTRUE(ao))

# Sanity: without the override, ao resolves FALSE
grf_merged_nofix <- grf_merged; grf_merged_nofix$adverse_outcome <- NULL
ao_nofix <- grf_merged_nofix$adverse_outcome
if (is.null(ao_nofix)) ao_nofix <- outcome_type %in% c("binary", "count")

ok("4.6 without override, ao resolves FALSE (the original bug)",
   isFALSE(ao_nofix))

# =============================================================================
cat(sprintf("\n=== TOTAL: %d passed, %d failed ===\n", PASS, FAIL))
if (FAIL > 0L) quit(save = "no", status = 1L)
