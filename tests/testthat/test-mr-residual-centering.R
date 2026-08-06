# =============================================================================
# Contract: the IJ residual is centered over the draws the IJ uses
# =============================================================================
# THE ASSERTION THAT WOULD HAVE CAUGHT F13.
#
# The defect was never a denominator on its own.  It was that the two bias terms
# carried DIFFERENT denominators: bias_sel over the draws that produced a winner,
# bias_fix over all B.  The residual of Eq. 12,
#
#     r_b = (bias_sel + bias_fix) - D_{H*_b}(b) - D_H(b),
#
# then mixed differently normalised quantities, so mean(r_b) was not zero -- and
# mean(r) = 0 is exactly the condition Eq. 13's uncentered rbar2 needs in order
# to equal Wager-Hastie-Efron's centered v-hat.
#
# Two repairs are coherent and BOTH center to machine precision, so centering
# does not choose between them; they differ in estimand:
#
#   (iii) Eq. 9 literal -- both terms over all B, D := 0 on a no-winner draw,
#         IJ over all B.
#   (iv)  CONDITIONAL   -- both terms over the winner draws, IJ over the same.
#
# The package implements (iv): the reported analysis exists only because a
# subgroup was identified, so the correction is estimated on the draws where a
# selection occurred.  This test asserts the property that holds under it.
#
# The two conventions COINCIDE at selection_rate == 1, so a fixture must make
# the admission floor bind on a real minority of draws or it asserts nothing.
#
# This mirrors the arithmetic of fs_mr_inference.R rather than calling it,
# because the property is about how the pieces combine and the fixture needs a
# controlled selection rate.  .fs_mr_ij_var() is the package's own.

test_that("the residual is exactly centered over the draws the IJ uses", {
  set.seed(20260805)
  n <- 200L; draws <- 3000L; n_cand <- 3L; sel <- 1L
  bh <- c(0.30, 0.26, 0.22)
  floor_t <- 0.42                      # binds on a substantial minority

  db <- matrix(rnorm(n * n_cand), n, n_cand) / sqrt(n)
  Xi <- matrix(rpois(n * draws, 1) - 1, n, draws)
  P  <- crossprod(db, Xi)
  beta_star <- bh + P

  sel_bias <- rep(NA_real_, draws)
  for (b in seq_len(draws)) {
    cand <- which(beta_star[, b] >= floor_t)
    if (!length(cand)) next
    sel_bias[b] <- P[cand[which.max(beta_star[cand, b])], b]
  }

  selection_rate <- mean(!is.na(sel_bias))
  # The fixture must exercise the case the change is about.
  expect_gt(selection_rate, 0.5)
  expect_lt(selection_rate, 0.95)

  # Convention (iv): both bias terms over the SAME draws, and the IJ over those.
  ok_H           <- which(is.finite(sel_bias))
  selection_bias <- mean(sel_bias, na.rm = TRUE)
  fixed_bias     <- mean(P[sel, ok_H])
  r_H            <- (selection_bias + fixed_bias) - sel_bias - P[sel, ]

  # THE CONTRACT.
  expect_lt(abs(mean(r_H[ok_H])), 1e-10)
  expect_false(anyNA(r_H[ok_H]))       # defined on every draw the IJ uses

  # And the consequence: with mean(r) = 0 the coded finite-B correction IS
  # Wager's centered one, so Eq. 13 needs no amendment.
  ij <- forestsearch:::.fs_mr_ij_var(Xi, r_H, ok_H)
  rb <- r_H[ok_H]
  hat_V_wager <- ij$tilde_V - (n / length(ok_H)) * mean((rb - mean(rb))^2)
  expect_equal(ij$hat_V, hat_V_wager, tolerance = 1e-12)
  expect_identical(ij$B_ok, length(ok_H))

  # The OLD convention -- bias_fix over all B while bias_sel is over winners --
  # fails this, which is what makes the test meaningful.
  fixed_old <- mean(P[sel, ])
  r_old     <- (selection_bias + fixed_old) - sel_bias - P[sel, ]
  expect_gt(abs(mean(r_old[ok_H])), 1e-6)

  # The mismatch is exactly the difference between the two denominators.
  expect_equal(mean(r_old[ok_H]), fixed_old - mean(P[sel, ok_H]),
               tolerance = 1e-12)
})


test_that("at selection_rate == 1 the conditional and all-B forms coincide", {
  set.seed(11L)
  n <- 150L; draws <- 800L; sel <- 1L
  db <- matrix(rnorm(n * 3L), n, 3L) / sqrt(n)
  Xi <- matrix(rpois(n * draws, 1) - 1, n, draws)
  P  <- crossprod(db, Xi)
  beta_star <- c(0.30, 0.26, 0.22) + P

  # No floor: every draw admits a candidate.
  sel_bias <- vapply(seq_len(draws),
                     function(b) P[which.max(beta_star[, b]), b], numeric(1))
  expect_identical(mean(!is.na(sel_bias)), 1)

  ok_H <- which(is.finite(sel_bias))
  expect_identical(ok_H, seq_len(draws))
  # With every draw a winner draw the two denominators are the same set, so the
  # conditional repair is a no-op -- which is why rate-1 analyses are unaffected.
  expect_identical(mean(P[sel, ok_H]), mean(P[sel, ]))

  r <- (mean(sel_bias) + mean(P[sel, ok_H])) - sel_bias - P[sel, ]
  expect_lt(abs(mean(r)), 1e-10)
})
