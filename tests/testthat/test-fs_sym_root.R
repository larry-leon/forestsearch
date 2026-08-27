# Tests for fs_sym_root().
#
# All deterministic (fixed seeds, no Monte Carlo over the assertions), so
# nothing here is flaky on CRAN.  The continuity test is the point of the
# function's existence: the asymmetric root V D^{1/2} satisfies the same
# covariance identity yet is NOT a continuous function of its input on a
# rank-deficient matrix with degenerate eigenvalues -- the case the candidate
# covariances of the analytic documents present by construction (complement
# pairs), and the mechanism behind the 94fd4dad figure movement.  A test that
# only shows the fix works is weaker than one that shows the defect was real,
# so the asymmetric root's failure is asserted alongside the symmetric root's
# success on the same input.

# ---- fixtures ---------------------------------------------------------------

orth <- function(n, seed) {
  set.seed(seed)
  qr.Q(qr(matrix(rnorm(n * n), n)))
}

# Full-rank 16 x 16, distinct eigenvalues 1..16.
Q_full <- orth(16, 461)
S_full <- Q_full %*% diag(as.numeric(1:16)) %*% t(Q_full)
S_full <- (S_full + t(S_full)) / 2

# Rank 7 of 16 -- the shape that motivated the promotion -- with repeated
# nonzero eigenvalues, so the spectrum carries degenerate subspaces both at
# zero and away from it.
ev_def <- c(3, 3, 2, 1, 1, 0.5, 0.25, rep(0, 9))
Q_def  <- orth(16, 4642)
S_def  <- Q_def %*% diag(ev_def) %*% t(Q_def)
S_def  <- (S_def + t(S_def)) / 2

# The asymmetric root the function replaces: same covariance identity,
# basis-dependent.
asym_root <- function(S, scale = 2) {
  eS <- eigen((S + t(S)) / 2, symmetric = TRUE)
  eS$vectors %*% diag(sqrt(pmax(scale * eS$values, 0)), nrow = nrow(S))
}

# ---- covariance identity ----------------------------------------------------

test_that("R %*% t(R) reproduces scale * S on a full-rank input", {
  R <- fs_sym_root(S_full, scale = 2)
  expect_lt(max(abs(R %*% t(R) - 2 * S_full)), 1e-12)
})

test_that("R %*% t(R) reproduces scale * S on a rank-7-of-16 input", {
  expect_equal(sum(eigen(S_def, symmetric = TRUE,
                         only.values = TRUE)$values > 1e-8), 7L)
  R <- fs_sym_root(S_def, scale = 2)
  expect_lt(max(abs(R %*% t(R) - 2 * S_def)), 1e-12)
})

# ---- symmetry ---------------------------------------------------------------

test_that("the returned root is symmetric", {
  expect_lt(max(abs(fs_sym_root(S_full) - t(fs_sym_root(S_full)))), 1e-14)
  expect_lt(max(abs(fs_sym_root(S_def)  - t(fs_sym_root(S_def)))),  1e-14)
})

# ---- continuity: the property that motivated the function -------------------

test_that("a 1e-12 perturbation moves the symmetric root by less than 1e-6,
           while the asymmetric root V D^{1/2} fails this on the same input", {
  set.seed(97)
  E <- matrix(rnorm(256), 16)
  P <- (E + t(E)) / 2
  P <- P / max(abs(eigen(P, symmetric = TRUE, only.values = TRUE)$values))
  S_pert <- S_def + 1e-12 * P

  d_sym <- max(abs(fs_sym_root(S_pert) - fs_sym_root(S_def)))
  expect_lt(d_sym, 1e-6)

  # The defect was real: the basis rotates inside the degenerate subspaces,
  # so the asymmetric root jumps by order one under the same perturbation.
  d_asym <- max(abs(asym_root(S_pert) - asym_root(S_def)))
  expect_gt(d_asym, 1e-6)
})

# ---- reproducibility across arithmetically equivalent routes ----------------

test_that("arithmetically equivalent constructions give bit-identical roots", {
  set.seed(31)
  A  <- matrix(rnorm(256), 16)                 # deliberately asymmetric
  S0 <- S_def + 1e-9 * A
  # Route 1 vs its transpose: (S0 + t(S0))/2 is bitwise symmetric in the
  # argument order, so the roots must be identical() -- not merely close.
  expect_identical(fs_sym_root(S0), fs_sym_root(t(S0)))
  # Route 2: pre-symmetrised input.  ((M + t(M))/2 == M bitwise for M built
  # as (S0 + t(S0))/2, so the root must again be bit-identical.
  M <- (S0 + t(S0)) / 2
  expect_identical(fs_sym_root(S0), fs_sym_root(M))
})

# ---- scale ------------------------------------------------------------------

test_that("scale is honoured", {
  R1 <- fs_sym_root(S_full, scale = 1)
  R2 <- fs_sym_root(S_full, scale = 2)
  expect_lt(max(abs(R1 %*% t(R1) - S_full)), 1e-12)
  expect_lt(max(abs(R2 %*% t(R2) - 2 * S_full)), 1e-12)
  expect_equal(R2, sqrt(2) * R1, tolerance = 1e-12)
})

# ---- clamping ---------------------------------------------------------------

test_that("a slightly negative eigenvalue is clamped to zero without error", {
  Qn <- orth(16, 8)
  Sn <- Qn %*% diag(c(as.numeric(15:1), -1e-15)) %*% t(Qn)
  Sn <- (Sn + t(Sn)) / 2
  expect_silent(R <- fs_sym_root(Sn, scale = 2))
  expect_false(any(is.nan(R)))
  expect_false(any(is.na(R)))
  # Reconstruction matches the clamped spectrum: the negative direction
  # contributes exactly zero.
  Sc <- Qn %*% diag(c(as.numeric(15:1), 0)) %*% t(Qn)
  expect_lt(max(abs(R %*% t(R) - 2 * (Sc + t(Sc)) / 2)), 1e-12)
})
