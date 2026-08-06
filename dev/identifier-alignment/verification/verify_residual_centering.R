# verify_residual_centering.R
#
# V7 -- the residual-centering condition underlying the finite-B correction.
#
# @sec-wager's bridge asserts that the manuscript's rbar2 = B^-1 sum_b r_b^2
# "plays the role of v-hat", where Wager-Hastie-Efron's v-hat is the CENTERED
# quantity B^-1 sum_b (t_b - tbar)^2.  The two coincide iff mean(r_b) = 0 over
# the draws the IJ actually uses.
#
# This script checks that condition on the package's real code path.  The IJ
# computation is neither re-implemented nor copied: line 22 binds
# forestsearch:::.fs_mr_ij_var (R/fs_mr_inference.R:195-206) directly, so what
# is exercised IS the shipped function.  The script therefore REQUIRES
# forestsearch to be installed; everything else in it is base R.
#
# It previously source()d a verbatim snapshot of that function, and was base-R
# only.  The snapshot was accurate when written and would have stayed accurate
# only by luck: a second copy of a package function is free to drift from the
# original, and V7 would then have gone on passing against a computation the
# package no longer uses -- a quantity defined once and reconstructed at a
# second site, which is the failure the surrounding audit exists to catch.
# Binding the real function makes "what is exercised is the shipped
# computation" true permanently rather than only at the moment of writing.
#
# Reproduced from the source at lines 434-479.  NOTE: this is the behaviour as
# it was WHEN THIS SCRIPT WAS WRITTEN, i.e. the defect being diagnosed:
#   sel_bias[b]    <- P[s, b]  where a winner exists, else NA
#   selection_bias <- mean(sel_bias, na.rm = TRUE)   # over WINNER draws
#   fixed_bias     <- mean(P[sel, ])                 # over ALL draws  <-- the
#                                                    #   mismatched denominator
#   r_H            <- (selection_bias + fb) - sel_bias - P[sel, ]
#   ijH            <- .fs_mr_ij_var(Xi, r_H, which(is.finite(sel_bias)))
# The package has since been repaired: fixed_bias now averages over the winner
# draws too, so both terms share one denominator and mean(r) is zero there.
# Neither arm below matches the shipped code any more -- both hold fixed_bias at
# all-B -- and that is deliberate: this script is the record of the DIAGNOSIS,
# and its "winners" arm is what the code did before the repair.

.fs_mr_ij_var <- forestsearch:::.fs_mr_ij_var   # the SHIPPED function, not a copy

set.seed(20260805)

report <- function(label, target, estimate, tol) {
  ok <- abs(estimate - target) <= tol
  cat(sprintf("  %-44s target %12.8f  est %12.8f  %s\n",
              label, target, estimate, if (ok) "PASS" else "FAIL"))
  invisible(ok)
}

# --- Build a family with a genuine no-winner subset -------------------------
make_run <- function(n = 250L, draws = 4000L, n_cand = 3L,
                     bh = c(0.30, 0.26, 0.22), floor_t = 0.42, sel = 1L) {
  db <- matrix(rnorm(n * n_cand), n, n_cand) / sqrt(n)
  Xi <- matrix(rpois(n * draws, 1) - 1, n, draws)
  P <- crossprod(db, Xi)
  beta_star <- bh + P
  admitted <- beta_star >= floor_t
  sel_bias <- rep(NA_real_, draws)
  winner <- rep(NA_integer_, draws)
  for (b in seq_len(draws)) {
    cand <- which(admitted[, b])
    if (!length(cand)) next
    s <- cand[which.max(beta_star[cand, b])]
    sel_bias[b] <- P[s, b]; winner[b] <- s
  }
  list(Xi = Xi, P = P, sel_bias = sel_bias, sel = sel, n = n, draws = draws)
}

analyse <- function(run, convention = c("winners", "B")) {
  convention <- match.arg(convention)
  P <- run$P; sel_bias <- run$sel_bias; sel <- run$sel; draws <- run$draws
  ok <- which(is.finite(sel_bias))
  selection_bias <- if (convention == "winners") {
    mean(sel_bias, na.rm = TRUE)                  # the shipped line
  } else {
    sum(sel_bias, na.rm = TRUE) / draws           # Eq. 9 as written
  }
  fixed_bias <- mean(P[sel, ])                    # over ALL draws, both ways
  r_H <- (selection_bias + fixed_bias) - sel_bias - P[sel, ]
  ij <- .fs_mr_ij_var(run$Xi, r_H, ok)            # the package's function
  rb <- r_H[ok]
  # Wager's v-hat is the CENTERED second moment; the code uses the raw one.
  hat_V_wager <- ij$tilde_V - (run$n / length(ok)) * mean((rb - mean(rb))^2)
  list(selection_bias = selection_bias, rbar = mean(rb),
       tilde_V = ij$tilde_V, hat_V_code = ij$hat_V, hat_V_wager = hat_V_wager,
       p = length(ok) / draws, ok = ok, fixed_bias = fixed_bias, P = P, sel = sel)
}

cat("== V7a: selection_rate < 1 (floor binds; some draws admit nobody) ==\n")
run_a <- make_run()
w <- analyse(run_a, "winners")
b <- analyse(run_a, "B")
cat(sprintf("  selection_rate p = %.4f\n", w$p))
cat(sprintf("  mean(r_b) over used draws: winners %+.6e   by-B %+.6e\n",
            w$rbar, b$rbar))
# The exact identity for the shipped convention.
ok <- w$ok
pred_rbar <- w$fixed_bias - mean(run_a$P[run_a$sel, ok])
report("mean(r) identity: fb - mean_ok(D_H)", pred_rbar, w$rbar, 1e-12)
report("V_tilde invariant to convention", w$tilde_V, b$tilde_V, 1e-12)
cat(sprintf("  hat_V as coded : winners %.6f   by-B %.6f   (differ)\n",
            w$hat_V_code, b$hat_V_code))
report("hat_V with Wager's centered v-hat", w$hat_V_wager, b$hat_V_wager, 1e-12)
cat("\n")

cat("== V7b: selection_rate == 1 (no floor; every draw has a winner) ==\n")
run_b <- make_run(floor_t = -Inf)
w1 <- analyse(run_b, "winners")
b1 <- analyse(run_b, "B")
cat(sprintf("  selection_rate p = %.4f\n", w1$p))
report("mean(r_b) is exactly zero", 0, w1$rbar, 1e-12)
report("both conventions coincide", w1$selection_bias, b1$selection_bias, 1e-12)
report("code hat_V == Wager hat_V", w1$hat_V_wager, w1$hat_V_code, 1e-12)
cat("  -> at p = 1 the bridge in @sec-wager is exact and F13 is vacuous.\n\n")

cat("== V7c: how far the raw second moment departs, as p falls ==\n")
cat(sprintf("  %-8s %-10s %-14s %-14s %-14s\n",
            "floor", "p", "mean(r)", "hat_V code", "hat_V Wager"))
for (ft in c(-Inf, 0.30, 0.42, 0.55, 0.70)) {
  r <- make_run(floor_t = ft)
  a <- analyse(r, "winners")
  cat(sprintf("  %-8s %-10.4f %+-14.6e %-14.6f %-14.6f\n",
              format(ft), a$p, a$rbar, a$hat_V_code, a$hat_V_wager))
}
cat("\n  Note: the code/Wager gap at fixed convention is (n/B_ok)*mean(r)^2,\n")
cat("  which is second order in mean(r); the CONVENTION gap in hat_V is\n")
cat("  first order in Delta and is what V6 measured -- V6 is in the sibling\n")
cat("  script verify_closedform_fixture.R, not this one.\n")
