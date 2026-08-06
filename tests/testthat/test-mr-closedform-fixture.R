# Closed-form fixture: the static half.
#
# Carries ONLY the Level 1 algebraic identities and the S3 IJ ratio, per
# dev/identifier-alignment/SPEC_closedform_fixture.md. These are deterministic
# properties that should fail in CI. The Level 2 and Level 3 Monte Carlo
# comparisons live in dev/closedform-fixture/closedform_fixture.qmd, which is
# rendered and read rather than run in a suite that would go flaky.
#
# The mean difference is the measure under test because its score is linear in
# (alpha, beta): Lemma 1's curvature term vanishes identically, so every
# per-candidate object is closed form. It is also the only measure fitted with
# lm(), the path 411a448 repaired -- a sigma_D assertion here fails on contact
# if that path breaks again, rather than degrading silently to n_family = 0.

SPEC_MD <- list(outcome_type = "continuous", effect_measure = "MD",
                treat.name = "trt", outcome.name = "y", event.name = NULL,
                offset.name = NULL, adjust_covariates = NULL,
                adverse_outcome = TRUE)

# Deterministic centered scores -- constructed, not sampled, so every arm mean
# and arm variance is fixed before any fit runs.
.fx_z <- function(m) { z <- stats::qnorm((seq_len(m) - 0.5) / m); z - mean(z) }

.fx_blocks <- function(n_per_arm, delta, s) {
  out <- NULL
  for (b in seq_along(delta)) {
    zb <- .fx_z(n_per_arm)
    out <- rbind(out, data.frame(
      trt = rep(c(1L, 0L), each = n_per_arm),
      y   = c(delta[b] + s[b] * zb, s[b] * zb)))
  }
  out$id <- seq_len(nrow(out))
  k <- 2 * n_per_arm
  list(df = out,
       rows = lapply(seq_along(delta), function(b) ((b - 1) * k + 1):(b * k)))
}

# Closed forms. NOTE the (n_a - 1) denominators: .dfbeta_glm() returns the
# one-step influence (X'X)^-1 x_i e_i / (1 - h_i), and for a two-group model
# h_i = 1/n_a, so db_i = e_i / (n_a - 1). The Neyman s^2/n_a form differs by
# the leverage factor and is NOT what the package computes.
.fx_cf_beta <- function(df, idx) {
  d <- df[idx, ]; mean(d$y[d$trt == 1L]) - mean(d$y[d$trt == 0L])
}
.fx_cf_db <- function(df, idx) {
  d <- df[idx, ]
  y1 <- d$y[d$trt == 1L]; y0 <- d$y[d$trt == 0L]
  n1 <- length(y1); n0 <- length(y0)
  db <- numeric(nrow(d))
  db[d$trt == 1L] <-  (y1 - mean(y1)) / (n1 - 1)
  db[d$trt == 0L] <- -(y0 - mean(y0)) / (n0 - 1)
  db
}
.fx_cf_sigma <- function(df, idx) sqrt(sum(.fx_cf_db(df, idx)^2))

NOFLOOR <- list(effect_floor = NULL, consistency = NULL)


test_that("Level 1: beta_hat is the arm-mean difference, exactly", {
  bl <- .fx_blocks(40, delta = c(1.2, 0.4), s = c(1.0, 1.3))
  df <- bl$df
  cands <- list(g1 = c(bl$rows[[1]], bl$rows[[2]][1:20]),
                g2 = c(bl$rows[[2]], bl$rows[[1]][1:20]))
  A <- .fs_mr_assemble(df, cands, SPEC_MD)

  expect_equal(length(A$beta_hat), 2L)
  for (g in seq_along(cands))
    expect_equal(A$beta_hat[g], .fx_cf_beta(df, cands[[g]]), tolerance = 1e-10)
})


test_that("Level 1: sigma_D is the one-step influence norm, exactly", {
  bl <- .fx_blocks(40, delta = c(1.2, 0.4), s = c(1.0, 1.3))
  df <- bl$df
  cands <- list(g1 = c(bl$rows[[1]], bl$rows[[2]][1:20]),
                g2 = c(bl$rows[[2]], bl$rows[[1]][1:20]))
  A <- .fs_mr_assemble(df, cands, SPEC_MD)

  for (g in seq_along(cands))
    expect_equal(A$sigma_D[g], .fx_cf_sigma(df, cands[[g]]), tolerance = 1e-10)

  # The Neyman s^2/n form is NOT the package's quantity. Asserted so that a
  # future change TO the Neyman form is caught rather than absorbed.
  d <- df[cands[[1]], ]
  y1 <- d$y[d$trt == 1L]; y0 <- d$y[d$trt == 0L]
  neyman <- sqrt(sum((y1 - mean(y1))^2) / length(y1)^2 +
                 sum((y0 - mean(y0))^2) / length(y0)^2)
  expect_false(isTRUE(all.equal(A$sigma_D[1], neyman, tolerance = 1e-6)))
})


test_that("Level 1: B columns are e_i/(n_a - 1) and sum to zero within arm", {
  bl <- .fx_blocks(40, delta = c(1.2, 0.4), s = c(1.0, 1.3))
  df <- bl$df
  cands <- list(g1 = c(bl$rows[[1]], bl$rows[[2]][1:20]),
                g2 = c(bl$rows[[2]], bl$rows[[1]][1:20]))
  A <- .fs_mr_assemble(df, cands, SPEC_MD)

  for (g in seq_along(cands)) {
    idx <- cands[[g]]
    expect_equal(A$B[idx, g], .fx_cf_db(df, idx), tolerance = 1e-10)
    # zero outside the candidate
    expect_equal(A$B[-idx, g], rep(0, nrow(df) - length(idx)), tolerance = 1e-12)
    # within-arm sums vanish
    tr <- df$trt[idx]
    expect_equal(sum(A$B[idx, g][tr == 1L]), 0, tolerance = 1e-10)
    expect_equal(sum(A$B[idx, g][tr == 0L]), 0, tolerance = 1e-10)
  }
})


test_that("Level 1: Sigma = crossprod(B) matches the membership computation", {
  bl <- .fx_blocks(40, delta = c(1.2, 0.4), s = c(1.0, 1.3))
  df <- bl$df
  cands <- list(g1 = c(bl$rows[[1]], bl$rows[[2]][1:20]),
                g2 = c(bl$rows[[2]], bl$rows[[1]][1:20]))
  A <- .fs_mr_assemble(df, cands, SPEC_MD)
  Sig <- crossprod(A$B)

  Bcf <- matrix(0, nrow(df), length(cands))
  for (g in seq_along(cands)) Bcf[cands[[g]], g] <- .fx_cf_db(df, cands[[g]])
  expect_equal(unname(Sig), unname(crossprod(Bcf)), tolerance = 1e-10)

  # sigma_D is the square root of the diagonal -- same object, not a parallel one
  expect_equal(unname(A$sigma_D), unname(sqrt(diag(Sig))), tolerance = 1e-12)

  # Disjoint memberships give exactly zero covariance.
  bl2 <- .fx_blocks(40, delta = c(1.0, 1.0), s = c(1.0, 1.0))
  A2 <- .fs_mr_assemble(bl2$df, list(g1 = bl2$rows[[1]], g2 = bl2$rows[[2]]),
                        SPEC_MD)
  expect_equal(crossprod(A2$B)[1, 2], 0, tolerance = 1e-14)
})


test_that("Level 1 reaches the top level: se_wald equals sigma_D of the winner", {
  # Guards the "asserting at the leaf" failure mode: proves fs_mr_inference()
  # actually calls the assembler and that the value reaches the returned object.
  bl <- .fx_blocks(40, delta = c(1.2, 0.4), s = c(1.0, 1.3))
  df <- bl$df
  cands <- list(g1 = c(bl$rows[[1]], bl$rows[[2]][1:20]),
                g2 = c(bl$rows[[2]], bl$rows[[1]][1:20]))
  A <- .fs_mr_assemble(df, cands, SPEC_MD)

  res <- fs_mr_inference(df, cands, SPEC_MD, selected_members = cands[[1]],
                         admission = NOFLOOR, reselection = "maxeff",
                         draws = 2000L, multiplier = "gaussian", seed = 42L)

  expect_equal(res$debiased$se_wald, A$sigma_D[1], tolerance = 1e-10)
  expect_equal(res$debiased$se_wald, .fx_cf_sigma(df, cands[[1]]), tolerance = 1e-10)
  expect_equal(res$n_family, 2L)
  expect_equal(res$selection_rate, 1)
})


test_that("S3: var_ij / sigma2_D approaches 4 as B grows", {
  # THE assertion that certifies both residual terms reach the variance. An
  # implementation dropping fixed_bias returns 1 here and passes everything
  # else in this fixture.
  bl <- .fx_blocks(60, delta = 1.0, s = 1.0)
  df <- bl$df; idx <- bl$rows[[1]]
  A <- .fs_mr_assemble(df, list(g1 = idx), SPEC_MD)
  s2 <- A$sigma_D[1]^2

  ratios <- vapply(c(5000L, 20000L, 80000L), function(B) {
    r <- fs_mr_inference(df, list(g1 = idx), SPEC_MD, selected_members = idx,
                         admission = NOFLOOR, reselection = "maxeff",
                         draws = B, multiplier = "gaussian", seed = 7L)
    expect_equal(r$selection_rate, 1)
    r$debiased$var_ij / s2
  }, numeric(1))

  # Approaches 4 from below; the largest B must be close, and must not be 1.
  expect_gt(ratios[3], 3.8)
  expect_lt(ratios[3], 4.2)
  expect_gt(ratios[3], ratios[1])      # moving toward 4, not away
  expect_gt(min(ratios), 3.0)          # decisively not the drop-fixed_bias value of 1
})


test_that("Level 3: mean_r is zero where the admission floor binds", {
  # THE assertion that would have caught F13 in CI rather than by a bisect three
  # commits later. rbar = 0 is the property whose absence WAS the defect: the two
  # bias terms carried different denominators, so the residual mixed differently
  # normalised quantities.
  #
  # Two blocks, so the candidate is a strict subset and the complement is
  # non-empty -- otherwise mean_r_c degrades to NA and asserts nothing.
  bl <- .fx_blocks(60, delta = c(1.0, 0.2), s = c(1.0, 1.0))
  df <- bl$df; idx <- bl$rows[[1]]
  A <- .fs_mr_assemble(df, list(g1 = idx), SPEC_MD)

  res <- fs_mr_inference(
    df, list(g1 = idx), SPEC_MD, selected_members = idx,
    admission = list(effect_floor = A$beta_hat[1] + 0.4307 * A$sigma_D[1],
                     consistency = NULL),
    reselection = "maxeff", draws = 20000L, multiplier = "gaussian",
    seed = 7L, include_complement = TRUE)

  # RATE FIRST. A change that inadvertently drives the rate to 1 must fail
  # loudly here rather than pass the residual check vacuously: at rate 1 the two
  # candidate draw sets coincide and centering holds trivially.
  expect_gt(res$selection_rate, 0.3)
  expect_lt(res$selection_rate, 0.9)

  expect_lt(abs(res$mean_r), 1e-10)
  expect_false(is.na(res$mean_r_c))
  expect_lt(abs(res$mean_r_c), 1e-10)

  # The old arrangement -- fixed_bias over all B, selection_bias over the winner
  # draws -- would NOT satisfy this. Without asserting the counterfactual is
  # non-zero, a correct implementation is indistinguishable from one where the
  # quantity happens to be small.
  rbar_old <- res$fixed_bias * res$selection_rate - res$fixed_bias
  expect_gt(abs(rbar_old), 1e-3)
})


test_that("Level 3: mean_r_c is NA when no complement is fit", {
  # A single candidate spanning the whole frame has an empty complement.
  # NA is missing data, not a value; it must not be substituted with zero.
  bl <- .fx_blocks(60, delta = 1.0, s = 1.0)
  df <- bl$df; idx <- bl$rows[[1]]

  res <- fs_mr_inference(df, list(g1 = idx), SPEC_MD, selected_members = idx,
                         admission = NOFLOOR, reselection = "maxeff",
                         draws = 2000L, multiplier = "gaussian",
                         seed = 7L, include_complement = TRUE)
  expect_true(is.na(res$mean_r_c))
  expect_lt(abs(res$mean_r), 1e-10)
})


test_that("Level 3 wiring: tilde_V is invariant to a constant shift in r", {
  # .fs_mr_ij_var() row-centers the multiplier matrix, so adding a constant to
  # every r_b cannot move tilde_V. Binds the SHIPPED function, not a copy.
  set.seed(99)
  N <- 40L; Bd <- 500L
  Xi <- matrix(stats::rnorm(N * Bd), N, Bd)
  r  <- stats::rnorm(Bd)
  ok <- seq_len(Bd)

  base  <- .fs_mr_ij_var(Xi, r, ok)
  for (shift in c(-3.5, 0.75, 10)) {
    sh <- .fs_mr_ij_var(Xi, r + shift, ok)
    expect_equal(sh$tilde_V, base$tilde_V, tolerance = 1e-12)
  }

  # Restricted to a strict subset of draws, the invariance still holds --
  # this is the selection_rate < 1 case.
  ok2 <- sort(sample(Bd, 180L))
  b2  <- .fs_mr_ij_var(Xi, r, ok2)
  s2  <- .fs_mr_ij_var(Xi, r + 4.25, ok2)
  expect_equal(s2$tilde_V, b2$tilde_V, tolerance = 1e-12)
  expect_equal(b2$B_ok, 180L)
})
