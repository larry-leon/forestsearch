# Self-contained validation of the consistency-resampling linearization.
# Base-R Breslow Cox for Surv(y,d) ~ z (single binary covariate), so it runs
# without the survival package. The TRUE consistency rate refits the SAME Cox
# on each random half; the approximation uses beta_hat +/- D with
# D = sum_i G_i * dfbeta_i. This isolates the quality of the linearization.

# ---- base-R Breslow Cox, single binary covariate -------------------------
cox_bin <- function(y, d, z) {
  ord <- order(y)
  y <- y[ord]; d <- d[ord]; z <- z[ord]
  et <- sort(unique(y[d == 1]))
  if (length(et) == 0) return(NULL)
  # risk-set counts at each event time
  Y1 <- sapply(et, function(t) sum(z[y >= t] == 1))
  Y0 <- sapply(et, function(t) sum(z[y >= t] == 0))
  s1 <- sapply(et, function(t) sum(z[y == t & d == 1]))      # sum z among failures
  dk <- sapply(et, function(t) sum(y == t & d == 1))         # #failures
  Ufun <- function(b) {
    eb <- exp(b); zbar <- (Y1 * eb) / (Y1 * eb + Y0)
    sum(s1 - dk * zbar)
  }
  if (sum(dk) < 2 || all(Y1 == 0) || all(Y0 == 0)) return(NULL)
  rt <- tryCatch(uniroot(Ufun, interval = c(-15, 15), extendInt = "yes", tol = 1e-9),
                 error = function(e) NULL)
  if (is.null(rt)) return(NULL)
  b <- rt$root; eb <- exp(b)
  zbar <- (Y1 * eb) / (Y1 * eb + Y0)
  info <- sum(dk * zbar * (1 - zbar))
  if (!is.finite(b) || !is.finite(info) || info <= 0) return(NULL)
  dLam0 <- dk / (Y1 * eb + Y0)                               # Breslow baseline
  # per-subject score residuals L_i (sum to ~0 at b_hat)
  L <- numeric(length(y))
  for (i in seq_along(y)) {
    idx <- which(et <= y[i])
    term2 <- exp(b * z[i]) * sum((z[i] - zbar[idx]) * dLam0[idx])
    L[i] <- d[i] * (z[i] - zbar[which.min(abs(et - y[i]))]) * (d[i] == 1) - term2
  }
  # cleaner score-residual: d_i(z_i - zbar(y_i)) - term2, with zbar at y_i if event
  zbar_at <- approx(et, zbar, xout = y, rule = 2)$y
  L <- d * (z - zbar_at) - sapply(seq_along(y), function(i) {
    idx <- which(et <= y[i]); exp(b * z[i]) * sum((z[i] - zbar[idx]) * dLam0[idx])
  })
  list(beta = b, info = info, dfbeta = L / info, n = length(y), d = sum(d))
}

# ---- TRUE consistency rate: repeated Bernoulli split + refit --------------
true_rate <- function(y, d, z, hr.c = 1.0, R = 1000, min_rows = 5, min_ev = 2, seed = 1) {
  set.seed(seed); thr <- log(hr.c); ns <- 0L; nv <- 0L; n <- length(y)
  for (r in seq_len(R)) {
    A <- rbinom(n, 1, 0.5) == 1
    if (sum(A) < min_rows || sum(!A) < min_rows) next
    if (sum(d[A]) < min_ev || sum(d[!A]) < min_ev) next
    fa <- cox_bin(y[A], d[A], z[A]); fb <- cox_bin(y[!A], d[!A], z[!A])
    if (is.null(fa) || is.null(fb)) next
    nv <- nv + 1L
    if (fa$beta > thr && fb$beta > thr) ns <- ns + 1L
  }
  if (nv > 0) ns / nv else NA_real_
}

# ---- approximation --------------------------------------------------------
approx_rate <- function(fit, hr.c = 1.0, draws = 4000, seed = 1) {
  thr <- log(hr.c); delta <- fit$beta - thr
  closed <- max(0, 2 * pnorm(delta / (sqrt(sum(fit$dfbeta^2)))) - 1)
  set.seed(seed); n <- fit$n
  G <- matrix(2 * rbinom(draws * n, 1, 0.5) - 1, nrow = draws)   # Rademacher = split
  D <- as.numeric(G %*% fit$dfbeta)
  mc <- mean((fit$beta - abs(D)) > thr)
  c(beta = fit$beta, sigma_D = sqrt(sum(fit$dfbeta^2)),
    rate_true = NA, rate_closed = closed, rate_mc = mc)
}

# ---- simulate a subgroup with target marginal HR --------------------------
sim_sub <- function(n, hr, cens = 0.45, seed = 1) {
  set.seed(seed)
  z <- rep(0:1, length.out = n)
  lam <- 0.1 * ifelse(z == 1, hr, 1)            # exponential hazards
  T <- rexp(n, lam)
  C <- rexp(n, rate = lam * cens / (1 - cens))  # tune censoring
  y <- pmin(T, C); d <- as.integer(T <= C)
  data.frame(y = y, d = d, z = z)
}

# ---- run a grid ------------------------------------------------------------
cat(sprintf("%5s %5s %8s %8s | %9s %9s %9s | %9s %9s\n",
            "n", "d", "HR", "sigmaD", "true", "closed", "mc", "err.clsd", "err.mc"))
cat(strrep("-", 86), "\n")
grid <- expand.grid(n = c(60, 100, 200, 400), hr = c(1.0, 1.5, 2.0, 3.0))
res <- data.frame()
for (k in seq_len(nrow(grid))) {
  df <- sim_sub(grid$n[k], grid$hr[k], seed = 100 + k)
  fit <- cox_bin(df$y, df$d, df$z)
  if (is.null(fit)) next
  tr <- true_rate(df$y, df$d, df$z, hr.c = 1.0, R = 1500, seed = 7)
  ar <- approx_rate(fit, hr.c = 1.0, draws = 4000, seed = 7)
  cat(sprintf("%5d %5d %8.2f %8.3f | %9.3f %9.3f %9.3f | %+9.3f %+9.3f\n",
              fit$n, fit$d, grid$hr[k], ar["sigma_D"], tr,
              ar["rate_closed"], ar["rate_mc"],
              ar["rate_closed"] - tr, ar["rate_mc"] - tr))
  res <- rbind(res, data.frame(n = fit$n, d = fit$d, hr = grid$hr[k],
                               true = tr, closed = ar["rate_closed"], mc = ar["rate_mc"]))
}
cat(strrep("-", 86), "\n")
cat(sprintf("Mean |error| closed = %.4f , mc = %.4f\n",
            mean(abs(res$closed - res$true)), mean(abs(res$mc - res$true))))
cat(sprintf("Max  |error| closed = %.4f , mc = %.4f\n",
            max(abs(res$closed - res$true)), max(abs(res$mc - res$true))))
