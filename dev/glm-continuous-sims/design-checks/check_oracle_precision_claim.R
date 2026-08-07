# Does oracle_design.R's Q1 measure anything about M or rho?
#
# Q1 reports  sd(p) / (sd(o) / sqrt(10))  where o uses nb = 2e6 and p uses
# nb = 2e4, and reads the resulting 24--32x as evidence that the oracle is
# "24--32x sharper ... across family sizes and correlation structures".
#
# Hypothesis: that ratio is an arithmetic consequence of the draw budgets,
#   sqrt(10 * 2e6 / 2e4) = sqrt(1000) = 31.62,
# identical for every (M, rho), and the 24--32 spread is noise in an sd
# estimated from 10 observations (~24% relative uncertainty).
#
# Three tests below. A, C confirm or refute the hypothesis; B reports the
# quantity Q1 should have reported.

set.seed(11)

oracle_bias <- function(bh, Sigma, nb) {
  M <- length(bh)
  L <- chol(Sigma)
  Z <- matrix(rnorm(nb * M), nb, M) %*% L
  w <- max.col(sweep(Z, 2, bh, "+"))
  mean(Z[cbind(seq_len(nb), w)])
}

make_sigma <- function(M, rho) {
  S <- matrix(rho, M, M)
  diag(S) <- 1
  S
}

grid <- expand.grid(M = c(2L, 5L, 20L), rho = c(0.0, 0.7))

# ---------------------------------------------------------------------------
# A. Matched draw budget. If the oracle is intrinsically sharper than the
#    pipeline-scale estimate, sd_o / sd_p < 1 here. If the Q1 ratio was an
#    artifact of nb, this is 1.0 for every configuration.
# ---------------------------------------------------------------------------
n_rep <- 200L
nb_matched <- 2e4

cat("== A: matched budget, nb = 2e4 for both, ", n_rep, " reps ==\n", sep = "")
cat(sprintf("  %-4s %-6s %12s %12s %10s\n", "M", "rho", "sd(oracle)", "sd(pipe)", "ratio"))
for (i in seq_len(nrow(grid))) {
  M <- grid$M[i]
  rho <- grid$rho[i]
  bh <- seq(0.30, 0.20, length.out = M)
  S <- make_sigma(M, rho)
  o <- replicate(n_rep, oracle_bias(bh, S, nb_matched))
  p <- replicate(n_rep, oracle_bias(bh, S, nb_matched))
  cat(sprintf(
    "  %-4d %-6.1f %12.2e %12.2e %10.2f\n",
    M, rho, sd(o), sd(p), sd(o) / sd(p)
  ))
}

# ---------------------------------------------------------------------------
# B. The quantity Q1 obscured: does the selection bias itself, and the
#    per-draw variability of estimating it, actually depend on M and rho?
#    This is what a precision argument needs and what the ratio table hid.
# ---------------------------------------------------------------------------
cat("\n== B: selection bias and its per-draw sd, by (M, rho) ==\n")
cat(sprintf("  %-4s %-6s %12s %12s\n", "M", "rho", "E[bias]", "sd at 2e4"))
for (i in seq_len(nrow(grid))) {
  M <- grid$M[i]
  rho <- grid$rho[i]
  bh <- seq(0.30, 0.20, length.out = M)
  S <- make_sigma(M, rho)
  big <- replicate(20L, oracle_bias(bh, S, 2e5))
  small <- replicate(200L, oracle_bias(bh, S, 2e4))
  cat(sprintf("  %-4d %-6.1f %12.4f %12.2e\n", M, rho, mean(big), sd(small)))
}

# ---------------------------------------------------------------------------
# C. Scaling. sd should fall as nb^(-1/2). If so, the Q1 ratio is fully
#    determined by the budgets and predicted by sqrt(10 * nb_o / nb_p).
# ---------------------------------------------------------------------------
cat("\n== C: sd vs nb at M = 5, rho = 0.7; predicted slope -0.5 ==\n")
bh <- seq(0.30, 0.20, length.out = 5L)
S <- make_sigma(5L, 0.7)
nbs <- c(2.5e3, 5e3, 1e4, 2e4, 4e4)
sds <- vapply(nbs, function(nb) sd(replicate(150L, oracle_bias(bh, S, nb))), numeric(1))
for (j in seq_along(nbs)) {
  cat(sprintf("  nb = %7.0f   sd = %.2e\n", nbs[j], sds[j]))
}
cat(sprintf("  fitted log-log slope = %.3f\n", coef(lm(log(sds) ~ log(nbs)))[2]))
cat(sprintf("  Q1's reported ratio, predicted as sqrt(10 * 2e6 / 2e4) = %.2f\n",
            sqrt(10 * 2e6 / 2e4)))
