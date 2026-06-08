# Extended self-contained validation: unadjusted AND covariate-adjusted Cox.
# General Breslow Cox (Newton-Raphson, reverse-cumsum risk sets) for p covariates,
# returning beta, p x p information, and the n x p dfbeta influence matrix.
# The adjusted consistency rate refits the SAME multi-covariate Cox on each
# random half (re-estimating the nuisance coefficient too); the approximation
# uses only the TREATMENT column of dfbeta:  D = sum_i G_i * dfbeta[i, "treat"].

# ---- general Breslow Cox -------------------------------------------------
coxmv <- function(y, d, X, maxit = 50, tol = 1e-9) {
  X <- as.matrix(X); n <- nrow(X); p <- ncol(X)
  ord <- order(y); y <- y[ord]; d <- d[ord]; X <- X[ord, , drop = FALSE]
  if (sum(d) < 2) return(NULL)
  rcs <- function(v) rev(cumsum(rev(v)))          # reverse cumulative sum
  first_idx <- match(y, y)                        # Breslow risk-set start (ties)
  ev <- which(d == 1)
  et_pos <- unique(first_idx[ev])                 # risk index of each event time
  dk <- as.numeric(table(first_idx[ev])[as.character(et_pos)])  # #events / time

  b <- rep(0, p)
  for (it in seq_len(maxit)) {
    w <- as.vector(exp(X %*% b))
    S0 <- rcs(w)
    S1 <- apply(X * w, 2, rcs)                     # n x p
    if (p == 1) S1 <- matrix(S1, ncol = 1)
    # reverse-cumsum of cross products for S2
    S2 <- array(0, c(n, p, p))
    for (a in 1:p) for (bb in a:p) {
      m <- rcs(w * X[, a] * X[, bb]); S2[, a, bb] <- m; S2[, bb, a] <- m
    }
    U <- rep(0, p); I <- matrix(0, p, p)
    for (jj in seq_along(et_pos)) {
      ri <- et_pos[jj]; s0 <- S0[ri]; xbar <- S1[ri, ] / s0
      # Breslow: each of the dk[jj] failures shares this risk set
      sfail <- colSums(X[first_idx == ri & d == 1, , drop = FALSE])
      U <- U + (sfail - dk[jj] * xbar)
      s2 <- S2[ri, , ] / s0
      I <- I + dk[jj] * (s2 - outer(xbar, xbar))
    }
    step <- tryCatch(solve(I, U), error = function(e) NULL)
    if (is.null(step)) return(NULL)
    b <- b + step
    if (max(abs(step)) < tol) break
  }
  if (any(!is.finite(b))) return(NULL)
  Iinv <- tryCatch(solve(I), error = function(e) NULL)
  if (is.null(Iinv)) return(NULL)

  # score residuals L_i (n x p), sum_i L_i = U(b_hat) ~ 0
  w <- as.vector(exp(X %*% b)); S0 <- rcs(w)
  S1 <- apply(X * w, 2, rcs); if (p == 1) S1 <- matrix(S1, ncol = 1)
  xbar_et <- S1[et_pos, , drop = FALSE] / S0[et_pos]      # nev x p
  a_k <- dk / S0[et_pos]                                  # Breslow dLambda0
  cumA <- cumsum(a_k)                                     # over event times (asc)
  cumAx <- apply(xbar_et * a_k, 2, cumsum); if (p == 1) cumAx <- matrix(cumAx, ncol = 1)
  et_times <- y[et_pos]
  nk <- findInterval(y, et_times)                         # #event times <= y_i
  L <- matrix(0, n, p)
  xbar_at <- matrix(0, n, p)
  pos <- match(y, et_times); haspos <- !is.na(pos)
  xbar_at[haspos, ] <- xbar_et[pos[haspos], ]
  for (i in seq_len(n)) {
    k <- nk[i]
    term2 <- if (k > 0) w[i] * (X[i, ] * cumA[k] - cumAx[k, ]) else rep(0, p)
    L[i, ] <- d[i] * (X[i, ] - xbar_at[i, ]) - term2
  }
  dfbeta <- L %*% Iinv                                    # n x p
  colnames(dfbeta) <- colnames(X)
  list(beta = b, info = I, dfbeta = dfbeta, n = n, d = sum(d))
}

# ---- true consistency rate: repeated split + refit (treat = column 1) -----
true_rate_mv <- function(y, d, X, hr.c = 1.0, R = 800, min_rows = 5, min_ev = 2, seed = 1) {
  set.seed(seed); thr <- log(hr.c); ns <- 0L; nv <- 0L; n <- length(y)
  X <- as.matrix(X)
  for (r in seq_len(R)) {
    A <- rbinom(n, 1, 0.5) == 1
    if (sum(A) < min_rows || sum(!A) < min_rows) next
    if (sum(d[A]) < min_ev || sum(d[!A]) < min_ev) next
    fa <- coxmv(y[A], d[A], X[A, , drop = FALSE])
    fb <- coxmv(y[!A], d[!A], X[!A, , drop = FALSE])
    if (is.null(fa) || is.null(fb)) next
    nv <- nv + 1L
    if (fa$beta[1] > thr && fb$beta[1] > thr) ns <- ns + 1L
  }
  if (nv > 0) ns / nv else NA_real_
}

# ---- approximation (treatment column only) --------------------------------
approx_mv <- function(fit, hr.c = 1.0, draws = 4000, seed = 1) {
  thr <- log(hr.c); dfb <- fit$dfbeta[, 1]                # treatment column
  sigmaD <- sqrt(sum(dfb^2)); delta <- fit$beta[1] - thr
  closed <- max(0, 2 * pnorm(delta / sigmaD) - 1)
  set.seed(seed); G <- matrix(2 * rbinom(draws * fit$n, 1, 0.5) - 1, nrow = draws)
  D <- as.numeric(G %*% dfb)
  mc <- mean((fit$beta[1] - abs(D)) > thr)
  c(beta = fit$beta[1], sigma_D = sigmaD, rate_closed = closed, rate_mc = mc)
}

# ---- simulate subgroup with treatment HR + a prognostic nuisance covariate-
sim_adj <- function(n, hr, beta_nuis = 0.7, cens = 0.45, seed = 1) {
  set.seed(seed)
  z  <- rep(0:1, length.out = n)               # treatment
  xn <- rnorm(n)                                # prognostic nuisance
  lam <- 0.1 * exp(log(hr) * z + beta_nuis * xn)
  T <- rexp(n, lam); C <- rexp(n, rate = lam * cens / (1 - cens))
  list(y = pmin(T, C), d = as.integer(T <= C), X = cbind(treat = z, xn = xn))
}

run_grid <- function(label, adjusted) {
  cat("\n===", label, "===\n")
  cat(sprintf("%5s %5s %6s %7s | %8s %8s %8s | %8s %8s\n",
              "n", "d", "HR", "sigmaD", "true", "closed", "mc", "e.clsd", "e.mc"))
  cat(strrep("-", 78), "\n")
  grid <- expand.grid(n = c(100, 200, 400), hr = c(1.5, 2.0))
  errs_c <- errs_m <- numeric(0)
  for (k in seq_len(nrow(grid))) {
    s <- sim_adj(grid$n[k], grid$hr[k], seed = 200 + k)
    X <- if (adjusted) s$X else s$X[, 1, drop = FALSE]
    fit <- coxmv(s$y, s$d, X); if (is.null(fit)) next
    tr <- true_rate_mv(s$y, s$d, X, hr.c = 1.0, R = 800, seed = 11)
    ar <- approx_mv(fit, hr.c = 1.0, draws = 4000, seed = 11)
    cat(sprintf("%5d %5d %6.2f %7.3f | %8.3f %8.3f %8.3f | %+8.3f %+8.3f\n",
                fit$n, fit$d, grid$hr[k], ar["sigma_D"], tr,
                ar["rate_closed"], ar["rate_mc"],
                ar["rate_closed"] - tr, ar["rate_mc"] - tr))
    errs_c <- c(errs_c, ar["rate_closed"] - tr); errs_m <- c(errs_m, ar["rate_mc"] - tr)
  }
  cat(sprintf("Mean |err| closed=%.4f  mc=%.4f ; Max |err| closed=%.4f  mc=%.4f\n",
              mean(abs(errs_c)), mean(abs(errs_m)), max(abs(errs_c)), max(abs(errs_m))))
}

# sanity: unadjusted via the general Cox should mirror the single-covariate run
run_grid("UNADJUSTED (treatment only)", adjusted = FALSE)
run_grid("ADJUSTED (treatment + prognostic nuisance, nuisance re-fit per half)",
         adjusted = TRUE)
