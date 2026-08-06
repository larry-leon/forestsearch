# verify_closedform_fixture.R
#
# Numerical verification of the closed-form results in
# theory_derivations_closedform_section.qmd, plus two results not yet in
# that document:
#
#   (V4) the M = 1 degenerate IJ limit  V_tilde -> 4 * sigma2_D, and
#   (V5) the closed-form discriminator between the divide-by-B and
#        divide-by-winners conventions for bias_sel under an active
#        admission floor (audit finding F13).
#
# Base R only. Each block implements the *pipeline* formula directly
# (multiplier draws, influence dot products, the IJ of @eq-ij) and
# compares it to the closed form. Nothing is simulated from the target
# law where the construction can be simulated instead.

set.seed(20260805)

report <- function(label, target, estimate, tol) {
  ok <- abs(estimate - target) <= tol
  cat(sprintf(
    "  %-46s target %10.6f   est %10.6f   |diff| %8.6f   tol %8.6f   %s\n",
    label, target, estimate, abs(estimate - target), tol,
    if (ok) "PASS" else "FAIL"
  ))
  invisible(ok)
}

# ---------------------------------------------------------------------------
# V1. Exact MVN law of the perturbations (@eq-exact-mvn), and
# V2. bias_sel = tau * dnorm(d / tau) for two admitted candidates
#     (@eq-bias-two), Gaussian multipliers, five configurations.
# ---------------------------------------------------------------------------
cat("== V1 + V2: exact law and two-candidate selection bias ==\n")

run_two_candidate <- function(n, overlap, scale2, d_over_tau, b_draws = 4e5,
                              chunk = 2e3) {
  # Influence vectors with controlled overlap: shared block of size
  # round(overlap * n) contributes to both candidates.
  n_shared <- round(overlap * n)
  base_1 <- rnorm(n)
  base_2 <- rnorm(n)
  db1 <- base_1 / sqrt(n)
  db2 <- numeric(n)
  db2[seq_len(n_shared)] <- base_1[seq_len(n_shared)]
  db2[seq(n_shared + 1L, n)] <- base_2[seq(n_shared + 1L, n)]
  db2 <- scale2 * db2 / sqrt(n)

  s11 <- sum(db1^2)
  s22 <- sum(db2^2)
  s12 <- sum(db1 * db2)
  tau <- sqrt(s11 + s22 - 2 * s12)
  d <- d_over_tau * tau

  # Multiplier draws in chunks; selection on perturbed effects
  # (both candidates admitted on every draw: no floor).
  sum_winner <- 0
  m1 <- m2 <- 0
  c11 <- c22 <- c12 <- 0
  b_done <- 0
  while (b_done < b_draws) {
    b_now <- min(chunk, b_draws - b_done)
    g <- matrix(rnorm(b_now * n), b_now, n)
    d1 <- g %*% db1
    d2 <- g %*% db2
    win1 <- (d + d1) > d2
    sum_winner <- sum_winner + sum(ifelse(win1, d1, d2))
    m1 <- m1 + sum(d1); m2 <- m2 + sum(d2)
    c11 <- c11 + sum(d1^2); c22 <- c22 + sum(d2^2); c12 <- c12 + sum(d1 * d2)
    b_done <- b_done + b_now
  }
  bias_mc <- sum_winner / b_draws
  bias_cf <- tau * dnorm(d / tau)

  # MC tolerance: sd of the winner's perturbation is <= max(s11, s22)^{1/2}
  # up to the selection shift; 4 sd of the mean is a generous gate.
  tol <- 4 * sqrt(max(s11, s22) / b_draws)
  report(
    sprintf("bias_sel  (ovl %.2f, sc %.1f, d/tau %.1f)",
            overlap, scale2, d_over_tau),
    bias_cf, bias_mc, tol
  )
  # Second-moment check of @eq-exact-mvn on the same draws.
  report("   Sigma_12 from draws", s12, c12 / b_draws,
         4 * sqrt(s11 * s22 / b_draws))
  invisible(NULL)
}

run_two_candidate(n = 400, overlap = 0.00, scale2 = 1.0, d_over_tau = 0.0)
run_two_candidate(n = 400, overlap = 0.50, scale2 = 1.0, d_over_tau = 0.0)
run_two_candidate(n = 400, overlap = 0.50, scale2 = 1.0, d_over_tau = 1.0)
run_two_candidate(n = 400, overlap = 0.80, scale2 = 0.6, d_over_tau = 0.5)
run_two_candidate(n = 400, overlap = 0.90, scale2 = 1.0, d_over_tau = 2.0)

# Tie case sanity constant: orthogonal, common scale sigma -> sigma / sqrt(pi).
cat(sprintf("  tie-case constant check: 1/sqrt(pi) = %.6f\n\n", 1 / sqrt(pi)))

# ---------------------------------------------------------------------------
# V3. E[max] constants for M orthogonal standard normals (@eq-bias-orthog):
#     integral of M * z * dnorm(z) * pnorm(z)^(M-1).
# ---------------------------------------------------------------------------
cat("== V3: E[max] of M iid standard normals ==\n")
emax <- function(m) {
  integrate(function(z) m * z * dnorm(z) * pnorm(z)^(m - 1),
            -Inf, Inf, rel.tol = 1e-10)$value
}
report("E[max], M = 2  vs 1/sqrt(pi)", 1 / sqrt(pi), emax(2), 1e-8)
report("E[max], M = 3  vs 3/(2 sqrt(pi))", 1.5 / sqrt(pi), emax(3), 1e-8)
report("E[max], M = 4  vs 1.02938", 1.02938, emax(4), 5e-6)
cat("\n")

# ---------------------------------------------------------------------------
# V4. The M = 1 degenerate IJ limit.
#
# Claim (derived this session, previously an open sketch):
# with a single candidate, H*_b = H on every draw, so the residual
# @eq-residual carries the same perturbation twice,
#     r_b = 2 * (mean(D) - D_b),
# and as B -> infinity the IJ covariances tend to -2 * db_i, whence
#     V_tilde -> 4 * sigma2_D
# for any centered unit-variance multiplier law. Implemented literally
# from @eq-biases, @eq-residual, @eq-ij with Poisson and Gaussian weights.
# ---------------------------------------------------------------------------
cat("== V4: M = 1 IJ limit, V_tilde / sigma2_D -> 4 ==\n")

run_m1_ij <- function(n, b_draws, law = c("poisson", "gaussian")) {
  law <- match.arg(law)
  db <- rnorm(n) / sqrt(n)
  sigma2_d <- sum(db^2)

  k <- switch(law,
    poisson  = matrix(rpois(b_draws * n, 1), b_draws, n),
    gaussian = 1 + matrix(rnorm(b_draws * n), b_draws, n)
  )
  d_b <- as.vector((k - 1) %*% db)          # @eq-perturb at the single g
  bias_sel <- mean(d_b)                     # winner is g on every draw
  bias_fix <- mean(d_b)                     # same term at M = 1
  r_b <- (bias_sel + bias_fix) - d_b - d_b  # @eq-residual

  k_centered <- sweep(k, 2, colMeans(k))
  cov_i <- as.vector(crossprod(k_centered, r_b)) / b_draws  # @eq-ij
  v_tilde <- sum(cov_i^2)
  v_hat <- v_tilde - (n / b_draws) * mean(r_b^2)

  cat(sprintf(
    "  law %-8s n %4d  B %6d   V_tilde/sigma2_D = %.4f   V_hat/sigma2_D = %.4f\n",
    law, n, b_draws, v_tilde / sigma2_d, v_hat / sigma2_d
  ))
  invisible(c(v_tilde / sigma2_d, v_hat / sigma2_d))
}

run_m1_ij(n = 300, b_draws = 5e3,  law = "poisson")
run_m1_ij(n = 300, b_draws = 2e4,  law = "poisson")
run_m1_ij(n = 300, b_draws = 8e4,  law = "poisson")
run_m1_ij(n = 300, b_draws = 2e4,  law = "gaussian")
cat("  (ratio should approach 4 in B; V_hat removes the finite-B excess)\n\n")

# ---------------------------------------------------------------------------
# V5. F13 discriminator: M = 1 with an active admission floor.
#
# Admission per draw on the perturbed effect: D_b >= c, c = t - beta_hat.
# Under Gaussian multipliers D ~ N(0, sigma2) exactly, so
#   divide-by-B      : E[D 1{D >= c}]  = sigma * dnorm(c / sigma)
#   divide-by-winners: E[D | D >= c]   = sigma * dnorm(c/sigma) / (1 - pnorm(c/sigma))
# The two conventions have *different* closed forms; a fixture pinned to
# one of them detects the other.
# ---------------------------------------------------------------------------
cat("== V5: F13 discriminator under an active floor (M = 1) ==\n")
sigma <- 0.8
c_floor <- 0.5 * sigma          # floor binds on ~31% of draws
b_draws <- 2e6
d_draw <- rnorm(b_draws, 0, sigma)
admit <- d_draw >= c_floor

est_by_b <- mean(d_draw * admit)
est_by_w <- mean(d_draw[admit])
cf_by_b <- sigma * dnorm(c_floor / sigma)
cf_by_w <- cf_by_b / pnorm(c_floor / sigma, lower.tail = FALSE)

report("divide-by-B       vs sigma*phi(c/sigma)", cf_by_b, est_by_b,
       4 * sigma / sqrt(b_draws))
report("divide-by-winners vs phi/Phi-bar ratio", cf_by_w, est_by_w,
       4 * sigma / sqrt(sum(admit)))
cat(sprintf("  separation of the two conventions: %.4f vs %.4f\n",
            cf_by_b, cf_by_w))
cat(sprintf("  inflation factor 1/Phi-bar(c/sigma) = %.4f\n\n", cf_by_w / cf_by_b))

# ---------------------------------------------------------------------------
# V6. F13: what the convention choice moves, and what it provably cannot.
#
# In @eq-residual only the SCALAR bias_sel differs between conventions, so
# switching convention shifts every r_b by the same constant
#     Delta = bias_sel^(winners) - bias_sel^(B) = bias_sel^(winners) * (1 - p),
# with p the selection rate. Three exact consequences:
#
#   (a) V_tilde is INVARIANT. cov_i = B^-1 sum_b (K_bi - Kbar_i) r_b, and the
#       row-centering annihilates any constant since sum_b (K_bi - Kbar_i) = 0.
#   (b) V_hat is NOT invariant: the finite-B correction carries mean(r^2),
#       which shifts by 2*Delta*rbar + Delta^2.
#   (c) The point estimate shifts by exactly -Delta.
#
# Verified below on synthetic P/Xi with a genuine no-winner subset. The
# residual is assembled per @eq-residual, the IJ per @eq-ij-matrix.
# ---------------------------------------------------------------------------
cat("== V6: F13 convention -- invariance of V_tilde, sensitivity of V_hat ==\n")

n <- 250L
b_tot <- 4000L
n_cand <- 3L

db_mat <- matrix(rnorm(n * n_cand), n, n_cand) / sqrt(n)   # B_eff
xi <- matrix(rpois(n * b_tot, 1) - 1, n, b_tot)            # centered Poisson
p_mat <- crossprod(db_mat, xi)                             # P = B_eff' Xi
sel <- 1L                                                  # the observed winner

# Winner per draw among admitted; a floor makes some draws admit nobody.
bh <- c(0.30, 0.26, 0.22)
beta_star <- bh + p_mat
floor_t <- 0.42
admitted <- beta_star >= floor_t
has_win <- colSums(admitted) > 0L
winner <- rep(NA_integer_, b_tot)
for (b in which(has_win)) {
  cand <- which(admitted[, b])
  winner[b] <- cand[which.max(beta_star[cand, b])]
}
sel_bias_vec <- rep(NA_real_, b_tot)
sel_bias_vec[has_win] <- p_mat[cbind(winner[has_win], which(has_win))]

p_rate <- mean(has_win)
bias_sel_w <- mean(sel_bias_vec, na.rm = TRUE)   # divide-by-winners
bias_sel_b <- sum(sel_bias_vec, na.rm = TRUE) / b_tot  # divide-by-B
fixed_bias <- mean(p_mat[sel, ])                 # divide-by-B, both conventions
delta <- bias_sel_w - bias_sel_b

cat(sprintf("  selection_rate p = %.4f   bias_sel(w) = %.5f   bias_sel(B) = %.5f\n",
            p_rate, bias_sel_w, bias_sel_b))
cat(sprintf("  Delta = bias_sel(w) - bias_sel(B) = %.5f   (identity check %.5f)\n",
            delta, bias_sel_w * (1 - p_rate)))

ij_from <- function(bias_sel_used) {
  r_all <- (bias_sel_used + fixed_bias) - sel_bias_vec - p_mat[sel, ]
  ok <- which(has_win)
  rb <- r_all[ok]
  xk <- xi[, ok, drop = FALSE]
  xc <- xk - rowMeans(xk)
  cov_i <- as.numeric(xc %*% rb) / length(ok)
  tilde_v <- sum(cov_i^2)
  hat_v <- tilde_v - (n / length(ok)) * mean(rb^2)
  c(tilde_V = tilde_v, hat_V = hat_v, rbar = mean(rb))
}

ij_w <- ij_from(bias_sel_w)
ij_b <- ij_from(bias_sel_b)

report("V_tilde invariant across conventions", ij_w[["tilde_V"]],
       ij_b[["tilde_V"]], 1e-12)
cat(sprintf("  V_hat  by-winners %.6e   by-B %.6e   shift %.6e\n",
            ij_w[["hat_V"]], ij_b[["hat_V"]],
            ij_b[["hat_V"]] - ij_w[["hat_V"]]))
# Predicted V_hat shift: -(n/B_ok) * (2*(-Delta)*rbar_w + Delta^2)
n_ok <- sum(has_win)
pred_shift <- -(n / n_ok) * (-2 * delta * ij_w[["rbar"]] + delta^2)
report("   V_hat shift vs predicted", pred_shift,
       ij_b[["hat_V"]] - ij_w[["hat_V"]], 1e-10)
cat(sprintf("  point estimate shifts by exactly -Delta = %+.5f\n", -delta))
