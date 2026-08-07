# Design checks for the GLM/continuous simulation workstream.
# Three questions the handoff needs answered with numbers, not assertions.

set.seed(11)

# ---------------------------------------------------------------------------
# Q1. How sharp is the B->infinity oracle relative to the pipeline's own B,
#     as M and the correlation structure vary?  If the oracle is not decisively
#     sharper, the error decomposition is not worth building.
# ---------------------------------------------------------------------------
oracle_bias <- function(bh, Sigma, nb) {
  M <- length(bh); L <- chol(Sigma)
  Z <- matrix(rnorm(nb * M), nb, M) %*% L
  w <- max.col(sweep(Z, 2, bh, "+"))
  mean(Z[cbind(seq_len(nb), w)])
}

make_sigma <- function(M, rho) {
  S <- matrix(rho, M, M); diag(S) <- 1; S
}

cat("== Q1: oracle precision vs pipeline-scale, B = 20,000 ==\n")
cat(sprintf("  %-4s %-6s %12s %12s %10s\n", "M", "rho", "oracle SE", "pipe SD", "ratio"))
for (M in c(2L, 5L, 20L)) {
  for (rho in c(0.0, 0.7)) {
    bh <- seq(0.30, 0.20, length.out = M)
    S  <- make_sigma(M, rho)
    o  <- replicate(10, oracle_bias(bh, S, 2e6))
    p  <- replicate(10, oracle_bias(bh, S, 2e4))
    cat(sprintf("  %-4d %-6.1f %12.2e %12.2e %10.0fx\n",
                M, rho, sd(o)/sqrt(10), sd(p), sd(p)/(sd(o)/sqrt(10))))
  }
}

# ---------------------------------------------------------------------------
# Q2. Is the conditional estimand beta(H-hat) analytically available for a
#     Gaussian DGM with a rectangle subgroup?  For MD it is a difference of
#     conditional means, so YES in closed form -- check against simulation.
# ---------------------------------------------------------------------------
cat("\n== Q2: conditional estimand for a rectangle subgroup, Gaussian DGM ==\n")
# Y = a + tau * 1{X1 > c} * A + eps ; subgroup H = {X1 > c}
# beta(H) = tau exactly; beta(H') for a MIS-specified rectangle is a mixture.
tau <- 0.8; cut_true <- 0.5; n <- 200000L
X1 <- runif(n); A <- rbinom(n, 1, 0.5)
Y  <- 1 + tau * (X1 > cut_true) * A + rnorm(n)

for (cut_used in c(0.5, 0.3, 0.7)) {
  H <- X1 > cut_used
  emp <- mean(Y[H & A == 1]) - mean(Y[H & A == 0])
  # closed form: tau * P(X1 > cut_true | X1 > cut_used)
  cf <- tau * (1 - max(cut_true, cut_used)) / (1 - cut_used)
  cat(sprintf("  cut %.1f  empirical %.4f   closed form %.4f   diff %.4f\n",
              cut_used, emp, cf, abs(emp - cf)))
}

# ---------------------------------------------------------------------------
# Q3. Does the exact-Sigma property survive a REAL fitted family, or only a
#     constructed one?  Sigma = crossprod(B) is exact for MD by construction,
#     so the question is whether it stays exact once memberships overlap.
# ---------------------------------------------------------------------------
cat("\n== Q3: Sigma = crossprod(B) with overlapping memberships (MD) ==\n")
nn <- 400L
y  <- rnorm(nn); a <- rbinom(nn, 1, 0.5)
mk_db <- function(idx) {
  d <- numeric(nn)
  i1 <- idx[a[idx] == 1]; i0 <- idx[a[idx] == 0]
  d[i1] <-  (y[i1] - mean(y[i1])) / (length(i1) - 1)
  d[i0] <- -(y[i0] - mean(y[i0])) / (length(i0) - 1)
  d
}
g1 <- 1:250; g2 <- 150:400; g3 <- c(1:100, 300:400)
B  <- cbind(mk_db(g1), mk_db(g2), mk_db(g3))
S_analytic <- crossprod(B)

Xi <- matrix(rnorm(nn * 4e5), nn, 4e5)
S_draws <- tcrossprod(crossprod(B, Xi)) / 4e5
cat(sprintf("  max |Sigma_analytic - Sigma_draws| = %.2e   (max |Sigma| = %.2e)\n",
            max(abs(S_analytic - S_draws)), max(abs(S_analytic))))
cat("  overlaps: g1/g2 = 101 subjects, g1/g3 = 100, g2/g3 = 101\n")
