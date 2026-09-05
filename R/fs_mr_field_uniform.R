# fs_mr_field_uniform.R
#
# Uniform calibration of the field interval: the kappa(Sigma-hat) sweep
# (TASK_mr_field_uniform_2026-09-05; method acknowledged by Larry as a
# proposal for evaluation).  Called by fs_mr_inference() after the field
# block when `field_uniform = TRUE`; also callable directly on a synthetic
# Sigma for the limit-experiment reference checks (dev/tasks/limit_sweep.py,
# limit_sweep_K.py travel with the task document).
#
# Everything here CONSUMES the gate's objects and selection semantics; it
# alters neither.  Draws are made under the caller-supplied seed (the gate
# derives `seed + 910000L`, mirroring the field block's `seed + 900000L`),
# AFTER the field's own stream is fully consumed, so all prior output --
# "ij", "wald", and the field block itself -- is byte-identical whether or
# not the sweep runs.

#' Uniform (kappa) calibration of the field interval
#'
#' Computes, at analysis time and from the trial's own influence structure,
#' the smallest widening factor `kappa` in `[1, 2]` such that the two-sided
#' field interval `c + kappa * (q - c)` attains `1 - alpha` coverage
#' uniformly over a winner-profile protection family: hypothetical true
#' fields equal to the gate's shrunk field `w` with the winner's entry set to
#' `max_{g != winner} w_g + delta * sigma_sel`, for `delta` on `delta_grid`
#' (in SE units of the winner).  Nothing is tuned to any simulation design.
#'
#' The computation restricts to a mass-carrying candidate set (the smallest
#' `M <= M_cap` candidates by re-selection frequency whose cumulative share
#' of the re-selection mass is at least `mass`, the winner always included),
#' then per `delta` replays `R_rep` hypothetical trials with Gaussian
#' multipliers of covariance exactly `Sigma_hat` (via `crossprod(B, xi)`, or
#' a Cholesky root of a supplied `Sigma`), applies the gate's own selection
#' map (thresholded argmax under `reselection = "maxeff"`, vectorized; any
#' other rule via the per-draw `.fs_mr_select` path), drops trials with no
#' winner (conditioning on detection, as in the real analysis), and runs the
#' field procedure on each hypothetical trial to obtain the coverage
#' profiles C1(delta) (one-sided) and C2(delta; kappa) (two-sided, widened).
#'
#' The guarantee, as documented for the gate: the one-sided bound is
#' uniformly valid; the two-sided interval widened by the returned `kappa`
#' is uniformly valid over the winner-profile family; the plain quantile
#' interval is approximate.
#'
#' @param B Optional n x K influence matrix (the gate's `asm$B`); multiplier
#'   perturbations are `crossprod(B, rnorm(n))`, covariance `t(B) %*% B`.
#' @param Sigma Optional K x K covariance, used when `B` is `NULL` (the
#'   reference checks); perturbations via its Cholesky root.
#' @param w Numeric K: the gate's shrunk field (`beta_hat` with the winner's
#'   entry replaced by the de-biased estimate), working scale.
#' @param sel Integer: the winner's index in `w`.
#' @param sigma_sel The winner's SE on the working scale (`sdv[sel]`).
#' @param p_hat Numeric K: re-selection frequencies from the gate's
#'   multiplier pass (`tabulate(winner)/draws`).
#' @param t_g Numeric K admission floors, or `NULL` for no admission filter
#'   (the gate's resolved `t_g`, restricted with the candidate set).
#' @param reselection,sz,effect_neighborhood,selection_rule,log_scale The
#'   gate's selection-map configuration; `zcons_c`, when non-`NULL`, is the
#'   consistency centre `c_cons` and `sdv` the per-candidate SEs so the
#'   `maxcons` score can be reconstructed on the restricted family.
#' @param sdv Numeric K per-candidate SEs (needed only for `maxcons`).
#' @param zcons_c Consistency centre `c_cons`, or `NULL`.
#' @param delta_grid Winner-separation grid in SE units (H2 default 0-4 by 0.5).
#' @param mass,M_cap Mass-carrying reduction controls (H3 defaults 0.99, 12).
#' @param R_rep,R_out,R_in Monte Carlo sizes per delta (H4 defaults 300/300/150).
#' @param alpha Two-sided miscoverage target (0.05).
#' @param kappa_grid Candidate widening factors (default 1 to 2 by 0.01).
#' @param seed Optional integer; drawn once at entry.
#' @return List: `kappa` (kappa*), `kappa_mcse` (Monte Carlo SE of kappa*
#'   from the binding profile's binomial error and the local slope of
#'   min-delta C2 in kappa), `M`, `mass_covered`, `keep` (indices), `minC1`,
#'   `C1` (profile over `delta_grid`), `C2_k1` (profile at kappa = 1),
#'   `C2_kstar` (profile at kappa*), `kappa_grid`, `C2_min` (min-over-delta
#'   profile over `kappa_grid`), `n_kept` (trials with a winner per delta),
#'   `delta_grid`, `timing_seconds`.
#' @keywords internal
fs_mr_field_uniform <- function(B = NULL, Sigma = NULL, w, sel, sigma_sel,
                                p_hat, t_g = NULL,
                                reselection = "maxeff", sz = NULL,
                                effect_neighborhood = 0.10,
                                selection_rule = "neighborhood",
                                log_scale = TRUE,
                                sdv = NULL, zcons_c = NULL,
                                delta_grid = seq(0, 4, by = 0.5),
                                mass = 0.99, M_cap = 12L,
                                R_rep = 300L, R_out = 300L, R_in = 150L,
                                alpha = 0.05,
                                kappa_grid = seq(1, 2, by = 0.01),
                                seed = NULL) {
  t0 <- proc.time()
  if (!is.null(seed)) set.seed(as.integer(seed))
  K <- length(w)
  stopifnot(sel >= 1L, sel <= K, length(p_hat) == K)

  # -- 1. Mass-carrying set: smallest M by p_hat with cumulative share >= mass,
  #       capped at M_cap, winner always included. --------------------------
  tot <- sum(p_hat)
  ord <- order(p_hat, decreasing = TRUE)
  cum <- cumsum(p_hat[ord])
  M0  <- if (tot > 0) which(cum >= mass * tot)[1] else K
  keep <- ord[seq_len(min(max(M0, 1L), M_cap, K))]
  if (!(sel %in% keep)) keep <- c(keep[seq_len(length(keep) - 1L)], sel)
  keep <- sort(unique(keep))
  M <- length(keep)
  mass_covered <- if (tot > 0) sum(p_hat[keep]) / tot else NA_real_
  j0 <- match(sel, keep)

  w_r  <- w[keep]
  tg_r <- if (!is.null(t_g)) t_g[keep] else NULL
  sz_r <- if (!is.null(sz)) sz[keep] else rep(1L, M)
  sdv_r <- if (!is.null(sdv)) sdv[keep] else NULL

  # -- Draw machinery: covariance Sigma_hat restricted to `keep`, exactly. --
  if (!is.null(B)) {
    Br <- B[, keep, drop = FALSE]
    n  <- nrow(Br)
    draw <- function(R) crossprod(Br, matrix(stats::rnorm(n * R), n, R))  # M x R
  } else if (!is.null(Sigma)) {
    Lr <- t(chol(Sigma[keep, keep, drop = FALSE]))
    draw <- function(R) Lr %*% matrix(stats::rnorm(M * R), M, R)          # M x R
  } else stop("fs_mr_field_uniform: one of B or Sigma is required")

  # -- Selection map on the restricted family, the gate's own semantics. ----
  # Vectorized thresholded argmax for "maxeff" (identical to
  # .admit() + .fs_mr_select(rule = "maxeff"): passers = which(v >= t_g),
  # winner = passer with maximal v; which.max and max.col(ties = "first")
  # break ties identically).  Any other rule goes through .fs_mr_select per
  # draw, exactly as the field block's sel_one does.
  fast <- identical(reselection, "maxeff")
  sel_mat <- if (fast) {
    function(V) {                        # V: M x R -> integer R (NA = no winner)
      if (!is.null(tg_r)) V <- replace(V, V < tg_r, -Inf)
      g <- max.col(t(V), ties.method = "first")
      g[!is.finite(V[cbind(g, seq_len(ncol(V)))])] <- NA_integer_
      g
    }
  } else {
    zc <- if (!is.null(zcons_c) && !is.null(sdv_r))
      function(v) (v - zcons_c) / sdv_r else function(v) NULL
    one <- function(v) {
      pass <- if (!is.null(tg_r)) which(v >= tg_r) else seq_along(v)
      if (!length(pass)) return(NA_integer_)
      .fs_mr_select(v, zc(v), sz_r, pass, reselection,
                    effect_neighborhood, selection_rule, log_scale)
    }
    function(V) vapply(seq_len(ncol(V)), function(j) one(V[, j]), integer(1))
  }

  nde <- length(delta_grid)
  C1 <- rep(NA_real_, nde)
  n_kept <- integer(nde)
  cov2_hits <- matrix(0L, nde, length(kappa_grid))   # counts over kept trials
  ii_in <- seq_len(R_in)

  for (d in seq_len(nde)) {
    mu <- w_r
    if (M >= 2L) mu[j0] <- max(w_r[-j0]) + delta_grid[d] * sigma_sel

    Zw <- draw(R_rep)                     # zeta_r, M x R_rep
    W  <- mu + Zw
    G  <- sel_mat(W)
    kept <- which(!is.na(G))
    n_kept[d] <- length(kept)
    if (!length(kept)) next

    gidx <- rep(seq_len(R_out), each = R_in)

    Lam  <- rep(NA_real_, length(kept))
    q025 <- q95 <- q975 <- cbar <- rep(NA_real_, length(kept))

    for (i in seq_along(kept)) {
      r <- kept[i]
      # All perturbations FRESH PER TRIAL (the reference draws outer
      # perturbations per trial, limit_sweep_K.py; per-trial inner batches
      # additionally make the trials independent, so per-delta coverage has
      # plain binomial Monte Carlo error -- shared batches leave a common
      # non-averaging noise term that a K = 1 check exposes immediately).
      # Inner batches are still shared WITHIN the trial (across outer draws),
      # the reference economy.
      Zi1 <- draw(R_in)                 # m-hat(W_r) inner batch
      Zi2 <- draw(R_in)                 # field-pass inner batch
      Zo  <- draw(R_out)                # field outer draws
      Zi2_big <- Zi2[, rep(ii_in, times = R_out), drop = FALSE]  # M x (R_out*R_in)
      # m-hat(W_r) by the inner draws.
      Gi <- sel_mat(W[, r] + Zi1)
      oki <- which(!is.na(Gi))
      if (!length(oki)) next
      mW <- mean(Zi1[cbind(Gi[oki], oki)])
      Lam[i] <- Zw[G[r], r] - mW

      # The field procedure on this hypothetical trial.
      wr <- W[, r]; wr[G[r]] <- W[G[r], r] - mW
      V  <- wr + Zo
      Gs <- sel_mat(V)
      oks <- which(!is.na(Gs))
      if (length(oks) < 2L) next
      # Inner m-hat for every outer draw in one vectorized selection call,
      # then grouped means by outer draw via rowsum (no per-draw scans).
      Gin <- sel_mat(V[, gidx, drop = FALSE] + Zi2_big)
      zin <- Zi2_big[cbind(Gin, seq_along(Gin))]      # NA where no inner winner
      ok_in <- !is.na(zin)
      sums <- rowsum(ifelse(ok_in, zin, 0), gidx)     # R_out x 1, ordered 1:R_out
      cnts <- rowsum(as.numeric(ok_in), gidx)
      mV_all <- ifelse(cnts > 0, sums / cnts, NA_real_)
      Ls <- Zo[cbind(Gs[oks], oks)] - mV_all[oks]
      Ls <- Ls[is.finite(Ls)]
      if (length(Ls) < 2L) next
      qs <- stats::quantile(Ls, c(.025, .95, .975), names = FALSE, type = 7)
      q025[i] <- qs[1]; q95[i] <- qs[2]; q975[i] <- qs[3]; cbar[i] <- mean(Ls)
    }

    okc <- which(is.finite(Lam) & is.finite(q95))
    n_kept[d] <- length(okc)
    if (!length(okc)) next
    C1[d] <- mean(Lam[okc] <= q95[okc])
    for (k in seq_along(kappa_grid)) {
      kap <- kappa_grid[k]
      lo <- cbar[okc] + kap * (q025[okc] - cbar[okc])
      hi <- cbar[okc] + kap * (q975[okc] - cbar[okc])
      cov2_hits[d, k] <- sum(Lam[okc] >= lo & Lam[okc] <= hi)
    }
  }

  C2 <- sweep(cov2_hits, 1, pmax(n_kept, 1L), "/")
  C2[n_kept == 0L, ] <- NA_real_
  C2_min <- apply(C2, 2, min, na.rm = TRUE)
  hit <- which(C2_min >= 1 - alpha)
  kappa <- if (length(hit)) kappa_grid[hit[1]] else max(kappa_grid)

  # Monte Carlo SE of kappa*: binomial SE of the binding profile at kappa*,
  # divided by the local slope of the min-over-delta C2 in kappa.
  kidx <- match(kappa, kappa_grid)
  dbind <- which.min(C2[, kidx])
  seC <- sqrt(C2_min[kidx] * (1 - C2_min[kidx]) / max(n_kept[dbind], 1L))
  lo_k <- max(1L, kidx - 5L); hi_k <- min(length(kappa_grid), kidx + 5L)
  slope <- if (hi_k > lo_k)
    (C2_min[hi_k] - C2_min[lo_k]) / (kappa_grid[hi_k] - kappa_grid[lo_k]) else NA_real_
  kappa_mcse <- if (is.finite(slope) && slope > 0) seC / slope else NA_real_

  k1i <- match(1, kappa_grid)
  list(kappa = kappa, kappa_mcse = kappa_mcse,
       M = M, mass_covered = mass_covered, keep = keep,
       minC1 = suppressWarnings(min(C1, na.rm = TRUE)),
       C1 = C1,
       C2_k1 = if (!is.na(k1i)) C2[, k1i] else rep(NA_real_, nde),
       C2_kstar = C2[, kidx],
       kappa_grid = kappa_grid, C2_min = C2_min, C2 = C2,
       n_kept = n_kept, delta_grid = delta_grid,
       R_rep = as.integer(R_rep), R_out = as.integer(R_out),
       R_in = as.integer(R_in),
       timing_seconds = as.numeric((proc.time() - t0)["elapsed"]))
}
