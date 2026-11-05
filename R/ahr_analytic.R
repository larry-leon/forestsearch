## =============================================================================
## ahr_analytic.R
##
## The analytic layer of the MRCT consistency work: for any patient set of a
## Weibull AFT DGM built by inject_mrct_structure(), analysed at a calendar
## time with staggered (uniform) entry, the exact limit of the two-sample Cox
## estimate (the Struthers-Kalbfleisch root), its information and robust
## variance, the expected events, the correlation between nested sets'
## estimates, the AHR and CDE of the set, and the normal-approximation
## probabilities of the regional-consistency trigger at a trial size n.
##
## Moved from mrct-consistency (docs/theory/scripts/sk_limit.R,
## docs/theory/scripts/ahr_setup.R, design layer; sims/regional_entry_design.R)
## in the per-patient-entry formulation of sims/regional_entry_design.R: every
## patient carries their own entry factor Hw_i(u) and window density pW_i(u);
## with no entry delay every Hw_i is the common Hw and this is the library of
## sk_limit.R line for line. Nothing is drawn.
##
## Notation (techreport docs/manuscript/techreport/analytic-layer.qmd, Chapter
## 6; theory note docs/theory/AHR_mrct_dgm.qmd): patient i has hazard
## lambda0(s) exp(theta_a,i) under arm a, Lambda0(s) = (s e^-mu)^(1/tau),
## theta_1 = theta_0 + psi. Integrals are in operational time u = Lambda0(s);
## for a set S the entry-weighted at-risk and event functions are
##   A_H(u) = mean_S Hw_i(u) exp(-u e^theta_1i), a_H = mean_S Hw_i e^theta_1i exp(-u e^theta_1i)
## (B_H, b_H with theta_0), A_P, B_P likewise with pW_i, and the limiting score is
##   U(beta) = int { a_H - r (a_H + b_H) } du,  r = e^beta A_H / (e^beta A_H + B_H).
## Quadrature: composite Gauss-Legendre on geometric panels in u with breaks at
## each entry group's kink u(x - a - R) and end u(x - a).
## =============================================================================


## -----------------------------------------------------------------------------
## 0. Quadrature
## -----------------------------------------------------------------------------

#' Gauss-Legendre nodes and the matrix of cumulative integrals to the nodes
#' @noRd
.ahr_gl <- function(n) {
  i <- seq_len(n - 1); bb <- i / sqrt(4 * i^2 - 1)
  J <- matrix(0, n, n); J[cbind(i, i + 1)] <- bb; J[cbind(i + 1, i)] <- bb
  e <- eigen(J, symmetric = TRUE); o <- order(e$values)
  x <- e$values[o]; w <- 2 * e$vectors[1, o]^2
  ## Legendre P_0..P_n at the nodes
  P <- matrix(1, n, n + 1); P[, 2] <- x
  for (m in 1:(n - 1)) P[, m + 2] <- ((2 * m + 1) * x * P[, m + 1] - m * P[, m]) / (m + 1)
  ## C[j, k]: integral from -1 to x_j of the degree-(n-1) interpolant through f(x_k)
  IP <- matrix(0, n, n); IP[, 1] <- x + 1
  for (m in 1:(n - 1)) IP[, m + 1] <- (P[, m + 2] - P[, m]) / (2 * m + 1)
  coefm <- t(P[, 1:n] * w) * (2 * (0:(n - 1)) + 1) / 2
  list(x = x, w = w, C = IP %*% coefm)
}

#' Baseline cumulative hazard and its inverse
#' @noRd
.ahr_Lam0 <- function(s, pop) (s * exp(-pop$mu))^(1 / pop$tau)
#' @noRd
.ahr_s <- function(u, pop) exp(pop$mu) * u^pop$tau

#' The quadrature grid in u, with each entry group's entry factor
#'
#' Finite horizon: geometric panels from ufirst = 1e-3 / emax to the library's
#' kink u(horizon - R), one panel to u(horizon), plus every delay group's kink
#' and end as extra breaks (identical breaks dropped, so no delay gives the
#' library's grid exactly). Infinite horizon: geometric panels to 36 / emin, no
#' censoring (Hw = 1, pW = 0 for every group).
#' @noRd
.ahr_grid <- function(pop, horizon, nper = 16, per_decade = 4, emax, emin) {
  gl <- .ahr_gl(nper)
  ufirst <- 1e-3 / emax
  delays <- pop$delays
  if (is.finite(horizon)) {
    if (any(horizon - delays - pop$R <= 0))
      stop("horizon - entry delay - max_entry must be positive for every ",
           "entry group (horizon ", horizon, ").", call. = FALSE)
    ukink <- .ahr_Lam0(horizon - pop$R, pop); uend <- .ahr_Lam0(horizon, pop)
    k <- ceiling(per_decade * log10(ukink / ufirst))
    br_lib <- c(0, ufirst * (ukink / ufirst)^((0:k) / k), uend)
    ukink_g <- .ahr_Lam0(horizon - delays - pop$R, pop)
    uend_g <- .ahr_Lam0(horizon - delays, pop)
    br <- sort(unique(c(br_lib, ukink_g, uend_g)))
  } else {
    uend <- 36 / emin
    k <- ceiling(per_decade * log10(uend / ufirst))
    br <- c(0, ufirst * (uend / ufirst)^((0:k) / k))
    ukink_g <- rep(Inf, length(delays)); uend_g <- rep(Inf, length(delays))
  }
  np <- length(br) - 1; lo <- br[-length(br)]; hw <- diff(br) / 2
  u <- as.vector(t(outer(lo + hw, rep(1, nper)) + outer(hw, gl$x)))
  w <- as.vector(t(outer(hw, gl$w)))
  s <- .ahr_s(u, pop)
  if (is.finite(horizon)) {
    dens <- pop$tau * exp(pop$mu) * u^(pop$tau - 1) / pop$R
    Hw <- lapply(seq_along(delays), function(j)
      pmin(1, pmax(horizon - delays[j] - s, 0) / pop$R))
    pW <- lapply(seq_along(delays), function(j)
      ifelse(u > ukink_g[j] & u < uend_g[j], dens, 0))
  } else {
    Hw <- lapply(delays, function(a) rep(1, length(u)))
    pW <- lapply(delays, function(a) rep(0, length(u)))
  }
  list(u = u, w = w, br = br, np = np, nper = nper, C = gl$C, gw = gl$w, hw = hw,
       horizon = horizon, Hw = Hw, pW = pW, ukink = ukink_g, uend = uend_g)
}

#' Cumulative integral of node values f from 0 to each node
#' @noRd
.ahr_cum <- function(f, g) {
  F <- matrix(f, g$np, g$nper, byrow = TRUE)
  within <- (F %*% t(g$C)) * g$hw
  tot <- as.vector(F %*% g$gw) * g$hw
  start <- c(0, cumsum(tot))[seq_len(g$np)]
  as.vector(t(within + start))
}

#' Laplace-transform sums of a group: sum exp(-u e), sum e exp(-u e)
#' @noRd
.ahr_lt <- function(e, u) {
  S <- numeric(length(u)); D <- numeric(length(u))
  for (ch in split(seq_along(e), ceiling(seq_along(e) / 4000))) {
    M <- exp(-outer(e[ch], u)); S <- S + colSums(M); D <- D + colSums(M * e[ch])
  }
  list(S = S, D = D)
}

.ahr_fin <- function(v, to = 0) { v[!is.finite(v)] <- to; v }


## -----------------------------------------------------------------------------
## 1. The population
## -----------------------------------------------------------------------------

#' The population of the analytic layer
#'
#' Builds, from a Weibull AFT DGM made by [inject_mrct_structure()], the
#' per-patient quantities the analytic layer works on: each super-population
#' patient's log hazard under control (\eqn{\theta_0}) and log hazard ratio
#' (\eqn{\psi}, so \eqn{\theta_1 = \theta_0 + \psi}), the H* (`flag_harm`)
#' and region indicators, the Weibull `mu` and `tau`, and each patient's entry
#' window: entry at \eqn{a_i + E}, \eqn{E \sim U(0, R)}, analysis at calendar
#' time `analysis_time`.
#'
#' The DGM's `theta_0`, `theta_1` and `loghr_po` columns are checked to be
#' consistent with its linear predictors and with each other (as the
#' mrct-consistency `ahr_pop()` did); the layer reads the hazards from them and
#' uses `dgm$model_params$mu` as it stands, so apply [scale_dgm_hazard()]
#' before building the population.
#'
#' @section Limits:
#' Weibull AFT DGMs from [inject_mrct_structure()] only; entry uniform on
#' \eqn{(a_i, a_i + R)} for each entry group; administrative censoring only
#' (the random censoring model of the DGM, e.g. `cens_adjust`, is not in the
#' layer); correlations exact for nested sets. Each distinct delay is one
#' quadrature group with its own breaks, so a delay vector with many distinct
#' values is slow. Derivations: the technical report
#' (`docs/manuscript/techreport/analytic-layer.qmd`, Chapter 6, of the
#' mrct-consistency repository) and the theory notes
#' (`docs/theory/AHR_mrct_dgm.qmd`, `docs/theory/AHR_complete.qmd`).
#'
#' @param dgm An `aft_dgm_flex` object made by [inject_mrct_structure()].
#' @param analysis_time Calendar time of the analysis, in the DGM's time units.
#' @param max_entry Length R of each patient's uniform entry window.
#' @param entry_delay `NULL` (no delay), a non-negative numeric vector of
#'   length `nrow(dgm$df_super)` (each patient's delay \eqn{a_i}), or
#'   `list(var = <column>, delay = a)`: delay `a` for patients with
#'   `dgm$df_super[[var]] == 1`, none otherwise (e.g.
#'   `list(var = "z_region", delay = 6)`: the region enrols 6 months later).
#'   The analytic counterpart of `simulate_from_dgm(entry_delay = )`.
#' @return An object of class `ahr_population`: a list with `theta0`, `psi`,
#'   `band` (H*), `region`, `delay`, `mu`, `tau`, `x` (`analysis_time`), `R`
#'   (`max_entry`), `N`, the distinct `delays`, and an environment holding
#'   the quadrature grids and Laplace-sum caches (keyed by horizon).
#' @seealso [ahr_limit()], [ahr_design()], [ahr_solve()], [scale_dgm_hazard()]
#' @examples
#' \dontrun{
#' pop <- ahr_population(dgm, analysis_time = 40, max_entry = 11)
#' pop6 <- ahr_population(dgm, 40, 11, entry_delay = list(var = "z_region", delay = 6))
#' }
#' @export
ahr_population <- function(dgm, analysis_time, max_entry, entry_delay = NULL) {
  if (!inherits(dgm, c("aft_dgm_flex", "aft_dgm")) || is.null(dgm$mrct))
    stop("dgm must be an aft_dgm_flex object made by inject_mrct_structure().",
         call. = FALSE)
  ds <- dgm$df_super; mp <- dgm$model_params
  need <- c("theta_0", "theta_1", "lin_pred_0", "lin_pred_1", "loghr_po", "flag_harm")
  miss <- setdiff(need, names(ds))
  if (length(miss))
    stop("dgm$df_super lacks: ", paste(miss, collapse = ", "), ".", call. = FALSE)
  zr <- dgm$mrct$region$z_name
  if (is.null(zr) || !zr %in% names(ds))
    stop("dgm$df_super lacks the region indicator dgm$mrct$region$z_name.",
         call. = FALSE)
  if (!isTRUE(all.equal(ds$theta_0, -ds$lin_pred_0 / mp$tau)) ||
      !isTRUE(all.equal(ds$theta_1, -ds$lin_pred_1 / mp$tau)) ||
      !isTRUE(all.equal(ds$loghr_po, ds$theta_1 - ds$theta_0)))
    stop("dgm$df_super theta_0, theta_1 and loghr_po are not consistent with ",
         "the linear predictors and tau.", call. = FALSE)
  for (nm in c("analysis_time", "max_entry")) {
    v <- get(nm)
    if (!is.numeric(v) || length(v) != 1L || !is.finite(v) || v <= 0)
      stop(nm, " must be one finite positive number.", call. = FALSE)
  }
  N <- nrow(ds)
  delay <- if (is.null(entry_delay)) {
    rep(0, N)
  } else if (is.list(entry_delay)) {
    if (!all(c("var", "delay") %in% names(entry_delay)) ||
        !entry_delay$var %in% names(ds))
      stop("entry_delay must be list(var = <column of dgm$df_super>, delay = a).",
           call. = FALSE)
    a <- entry_delay$delay
    if (!is.numeric(a) || length(a) != 1L || !is.finite(a) || a < 0)
      stop("entry_delay$delay must be one non-negative number.", call. = FALSE)
    a * as.numeric(ds[[entry_delay$var]] == 1)
  } else {
    if (!is.numeric(entry_delay) || length(entry_delay) != N ||
        any(!is.finite(entry_delay)) || any(entry_delay < 0))
      stop("A numeric entry_delay must have length nrow(dgm$df_super) = ", N,
           " with finite non-negative values.", call. = FALSE)
    as.numeric(entry_delay)
  }
  delays <- sort(unique(delay))
  if (any(analysis_time - delays - max_entry <= 0))
    stop("analysis_time - entry delay - max_entry must be positive.", call. = FALSE)
  structure(list(theta0 = ds$theta_0, psi = ds$loghr_po, band = ds$flag_harm == 1,
                 region = ds[[zr]] == 1, delay = delay, delays = delays,
                 group = match(delay, delays), mu = unname(mp$mu), tau = mp$tau,
                 cell = paste0(ifelse(ds$flag_harm == 1, "band", "nonband"), "_",
                               ifelse(ds[[zr]] == 1, "region", "nonregion")),
                 x = analysis_time, R = max_entry, N = N, cache = new.env()),
            class = "ahr_population")
}

#' @rdname ahr_population
#' @param x An `ahr_population`.
#' @param ... Unused.
#' @export
print.ahr_population <- function(x, ...) {
  cat(sprintf("ahr_population: %d patients (region %.1f%%, H* %.1f%%), analysis at %g, entry U(0, %g)\n",
              x$N, 100 * mean(x$region), 100 * mean(x$band), x$x, x$R))
  if (any(x$delays > 0))
    cat(sprintf("  entry delays: %s\n", paste(format(x$delays), collapse = ", ")))
  invisible(x)
}


## -----------------------------------------------------------------------------
## 2. Sets and cells
## -----------------------------------------------------------------------------

## The standard sets as unions of the cells H* x region.
.ahr_set_cells <- list(
  all = c("band_region", "band_nonregion", "nonband_region", "nonband_nonregion"),
  nonregion = c("band_nonregion", "nonband_nonregion"),
  region = c("band_region", "nonband_region"),
  region_Gstar = "nonband_region",
  Hstar = c("band_region", "band_nonregion"),
  Hstar_region = "band_region",
  Hstar_nonregion = "band_nonregion")

.ahr_cell <- function(pop) pop$cell

.ahr_set_index <- function(pop, set) {
  if (is.character(set)) {
    if (length(set) != 1L || !set %in% names(.ahr_set_cells))
      stop("set must be one of ", paste(names(.ahr_set_cells), collapse = ", "),
           ", or a logical vector.", call. = FALSE)
    return(.ahr_cell(pop) %in% .ahr_set_cells[[set]])
  }
  if (!is.logical(set) || length(set) != pop$N || anyNA(set))
    stop("A logical set must have length pop$N = ", pop$N, " without NA.",
         call. = FALSE)
  if (!any(set)) stop("The set is empty.", call. = FALSE)
  set
}

## AHR and CDE of a set, closed form (CDE = the hazard ratio of the set at s -> 0+).
.ahr_causal <- function(pop, S, psi) {
  e0 <- exp(pop$theta0[S])
  c(ahr = exp(mean(psi[S])), cde = sum(e0 * exp(psi[S])) / sum(e0),
    band_share = mean(pop$band[S]), omega0 = sum(e0[pop$band[S]]) / sum(e0))
}

.ahr_psi <- function(pop, b, h, psi, b_R = b, h_R = h) {
  if (!is.null(psi)) {
    if (!is.numeric(psi) || length(psi) != pop$N || any(!is.finite(psi)))
      stop("psi must be a finite numeric vector of length pop$N.", call. = FALSE)
    return(psi)
  }
  if (is.null(b) || is.null(h)) return(pop$psi)
  for (v in list(b, h))
    if (!is.numeric(v) || length(v) != 1L || !is.finite(v) || v <= 0)
      stop("b and h must be positive numbers.", call. = FALSE)
  for (v in list(b_R, h_R))
    if (!is.numeric(v) || length(v) != 1L || !is.finite(v) || v <= 0)
      stop("b_R and h_R must be positive numbers.", call. = FALSE)
  psi <- log(b) + log(h / b) * pop$band
  if (b_R != b || h_R != h) {
    ## four cells H* x region: log b, log h outside the region, log b_R, log h_R in it
    psi[pop$region & !pop$band] <- log(b_R)
    psi[pop$region & pop$band] <- log(h_R)
  }
  psi
}

## The treated-arm multiplier of a standard cell at (b, h, b_R, h_R).
.ahr_cell_mult <- function(cell, b, h, b_R, h_R) {
  band <- startsWith(cell, "band_"); reg <- endsWith(cell, "_region")
  if (band) { if (reg) h_R else h } else { if (reg) b_R else b }
}


## -----------------------------------------------------------------------------
## 3. The core
## -----------------------------------------------------------------------------

## parts: one per type group k (band / nonband, or the whole set), each with the
## entry-weighted sums AH, aH, BH, bH, AP, BP over the group's patients at the
## nodes, its patient count n and its log-HR psi (NA if not constant).
.ahr_core <- function(parts, g, variance = TRUE) {
  sumf <- function(nm) Reduce(`+`, lapply(parts, `[[`, nm))
  n <- sum(vapply(parts, `[[`, 0, "n"))
  AH <- sumf("AH") / n; aH <- sumf("aH") / n; BH <- sumf("BH") / n; bH <- sumf("bH") / n
  AP <- sumf("AP") / n; BP <- sumf("BP") / n
  w <- g$w
  fin <- .ahr_fin
  rfun <- function(beta) fin(exp(beta) * AH / (exp(beta) * AH + BH), 0.5)
  score <- function(beta) sum(w * (aH - rfun(beta) * (aH + bH)))
  beta <- stats::uniroot(score, c(-10, 10), tol = 1e-13)$root
  dN <- w * (aH + bH) / 2                      # expected events per patient at the nodes (1:1)
  r <- rfun(beta)
  d1 <- sum(w * aH); d0 <- sum(w * bH)         # P(event in window | arm)
  i1 <- sum(dN * r * (1 - r))                  # information per patient, 1:1
  delta <- fin(log(AH / BH))                   # pooled log risk-set imbalance
  out <- list(beta = beta, hr = exp(beta), d1 = d1, d0 = d0, d = (d1 + d0) / 2,
              i1 = i1, n_pop = n, u = g$u, dN = dN, r = r, delta = delta,
              loghr_t = fin(log((aH / AH) / (bH / BH))),
              AH = AH, aH = aH, BH = BH, bH = bH, AP = AP, BP = BP)
  out$types <- do.call(rbind, lapply(names(parts), function(k) {
    p <- parts[[k]]; dNk <- w * (p$aH + p$bH) / (2 * n)
    dk <- sum(dNk)
    own <- if (is.na(p$psi)) rep(NA_real_, length(dNk)) else fin(log(p$aH / p$bH) - p$psi)
    data.frame(type = k, n = p$n, psi = p$psi, events = dk, share = dk / sum(dN),
               events_exp = sum(w * p$aH) / n, events_ctl = sum(w * p$bH) / n,
               m1_own = sum(dNk * own), m2_own = sum(dNk * own^2),
               m1_pool = sum(dNk * delta), m2_pool = sum(dNk * delta^2),
               stringsAsFactors = FALSE)
  }))
  rownames(out$types) <- NULL
  if (!variance) return(out)
  ## robust (Lin-Wei) variance, independent arms: per-patient influence
  ## phi = 1{T <= W} (z - r(T)) - g_z(min(T, W)), g_z(t) = e^{beta z} int_0^t (z - r) dLambda*
  q <- fin((aH + bH) / (exp(beta) * AH + BH))
  g1 <- .ahr_cum(exp(beta) * (1 - r) * q, g); g0 <- .ahr_cum(-r * q, g)
  f1 <- (1 - r) - g1; f0 <- -r - g0
  E1 <- sum(w * aH * f1^2) + sum(w * AP * g1^2)
  E0 <- sum(w * bH * f0^2) + sum(w * BP * g0^2)
  m1 <- sum(w * aH * f1) - sum(w * AP * g1)
  m0 <- sum(w * bH * f0) - sum(w * BP * g0)
  out$f1 <- f1; out$f0 <- f0; out$g1 <- g1; out$g0 <- g0
  out$sandwich <- (E1 + E0) / 2                # Bernoulli(1/2) allocation
  out$mean_phi <- c(m1, m0)
  out$var_model <- 1 / i1
  out$var_robust <- out$sandwich / i1^2
  out$var_robust_fixed <- (out$sandwich - m1^2) / i1^2   # allocation fixed at n/2
  out
}

## One part's entry-weighted arrays from per-entry-group Laplace sums.
## l1, l0: lists over entry groups (index j into g$Hw) of .ahr_lt() results.
.ahr_part <- function(l1, l0, js, g, n, psi) {
  p <- list(AH = 0, aH = 0, BH = 0, bH = 0, AP = 0, BP = 0, n = n, psi = psi)
  for (i in seq_along(js)) {
    H <- g$Hw[[js[i]]]; P <- g$pW[[js[i]]]
    p$AH <- p$AH + H * l1[[i]]$S; p$aH <- p$aH + H * l1[[i]]$D
    p$BH <- p$BH + H * l0[[i]]$S; p$bH <- p$bH + H * l0[[i]]$D
    p$AP <- p$AP + P * l1[[i]]$S; p$BP <- p$BP + P * l0[[i]]$S
  }
  p
}

## The cached grid of a horizon (emax, emin as the library's ahr_cache()).
.ahr_cache_grid <- function(pop, horizon, hmax = 30, bmin = 0.2) {
  key <- sprintf("grid|%s", format(horizon, digits = 15))
  if (is.null(pop$cache[[key]]))
    pop$cache[[key]] <- .ahr_grid(pop, horizon, emax = exp(max(pop$theta0)) * hmax,
                                  emin = exp(min(pop$theta0)) * bmin)
  pop$cache[[key]]
}

## Laplace sums of a (cell, entry group) at multiplier mult, cached per horizon.
.ahr_cache_lt <- function(pop, horizon, cell, j, mult) {
  key <- sprintf("lt|%s|%s|%d|%.12g", format(horizon, digits = 15), cell, j, mult)
  if (is.null(pop$cache[[key]])) {
    g <- .ahr_cache_grid(pop, horizon)
    ix <- which(.ahr_cell(pop) == cell & pop$group == j)
    pop$cache[[key]] <- .ahr_lt(exp(pop$theta0[ix]) * mult, g$u)
  }
  pop$cache[[key]]
}

## A standard set at (b, h) on the cached grid (the library's sk_fast()).
.ahr_fast <- function(pop, set, b, h, horizon, variance, b_R = b, h_R = h) {
  g <- .ahr_cache_grid(pop, horizon)
  cell_all <- .ahr_cell(pop)
  cells <- .ahr_set_cells[[set]]
  parts <- list()
  for (ty in c("band", "nonband")) {
    ck <- cells[startsWith(cells, paste0(ty, "_"))]
    if (!length(ck)) next
    l1 <- list(); l0 <- list(); js <- integer(0); n <- 0; mults <- numeric(0)
    for (k in ck) for (j in seq_along(pop$delays)) {
      nk <- sum(cell_all == k & pop$group == j)
      if (!nk) next
      mult <- .ahr_cell_mult(k, b, h, b_R, h_R)
      l1[[length(l1) + 1]] <- .ahr_cache_lt(pop, horizon, k, j, mult)
      l0[[length(l0) + 1]] <- .ahr_cache_lt(pop, horizon, k, j, 1)
      js <- c(js, j); n <- n + nk; mults <- c(mults, mult)
    }
    ## a type's log-HR is NA when its cells carry different multipliers
    if (n) parts[[ty]] <- .ahr_part(l1, l0, js, g, n,
                                    if (all(mults == mults[1])) log(mults[1]) else NA_real_)
  }
  out <- .ahr_core(parts, g, variance)
  out$grid <- g
  out
}

## Any set, any psi, on the set's own grid (the library's sk_limit()).
.ahr_direct <- function(pop, S, psi, two, horizon, variance) {
  rng <- range(c(pop$theta0[S], pop$theta0[S] + psi[S]))
  g <- .ahr_grid(pop, horizon, emax = exp(rng[2]), emin = exp(rng[1]))
  grp <- if (two) ifelse(pop$band, "band", "nonband") else rep("set", pop$N)
  parts <- lapply(split(which(S), grp[S]), function(ix) {
    js <- sort(unique(pop$group[ix]))
    l1 <- list(); l0 <- list()
    for (j in js) {
      ij <- ix[pop$group[ix] == j]
      e0 <- exp(pop$theta0[ij]); e1 <- e0 * exp(psi[ij])
      l1[[length(l1) + 1]] <- .ahr_lt(e1, g$u); l0[[length(l0) + 1]] <- .ahr_lt(e0, g$u)
    }
    .ahr_part(l1, l0, js, g, length(ix),
              if (two && all(psi[ix] == psi[ix[1]])) psi[ix[1]] else NA_real_)
  })
  out <- .ahr_core(parts, g, variance)
  out$grid <- g
  out
}

## Dispatch: the cached route for a standard set at (b, h) within the cache's
## range (0.2 <= b, h <= 30), else the direct route. force_fast: always the
## cached route (ahr_design(), whose correlations need one grid for all sets).
.ahr_exact <- function(pop, set, b, h, psi, horizon, variance, force_fast = FALSE,
                       b_R = b, h_R = h) {
  S <- .ahr_set_index(pop, set)
  psi_v <- .ahr_psi(pop, b, h, psi, b_R, h_R)
  two <- is.null(psi) && !is.null(b) && !is.null(h)
  fast <- two && is.character(set) &&
    (force_fast || (min(b, h, b_R, h_R) >= 0.2 && max(b, h, b_R, h_R) <= 30))
  out <- if (fast) .ahr_fast(pop, set, b, h, horizon, variance, b_R, h_R) else
    .ahr_direct(pop, S, psi_v, two, horizon, variance)
  out$omega <- if (two && "band" %in% out$types$type)
    out$types$share[out$types$type == "band"] else if (two) 0 else NA_real_
  out$S <- S
  out$psi_v <- psi_v
  out
}


## -----------------------------------------------------------------------------
## 4. Exports: ahr_limit(), ahr_design(), ahr_solve()
## -----------------------------------------------------------------------------

#' The exact limit of the two-sample Cox estimate on a patient set
#'
#' For a set of patients of an [ahr_population()], the limit \eqn{\beta^*} of
#' the log hazard ratio estimated by a Cox fit of treatment alone on that set
#' (1:1 randomisation, analysis at the population's calendar time with each
#' patient's uniform entry window, administrative censoring only): the root of
#' the limiting score (Struthers and Kalbfleisch), with the information, the
#' model-based and robust (Lin-Wei) variances per patient, the event
#' probabilities, the type table and \eqn{\omega} (the H* share of events),
#' and the set's AHR, CDE and its margHR at infinity (no censoring). Computed
#' by quadrature; nothing is drawn.
#'
#' A standard set at given (`b`, `h`) with \eqn{0.2 \le b, h \le 30} uses the
#' population's cached grid and Laplace sums (the route of [ahr_design()]);
#' any other call (a logical set, `psi`, or the DGM's own effects) builds the
#' set's own grid.
#'
#' @inheritSection ahr_population Limits
#' @param pop An [ahr_population()].
#' @param set A set name: `"all"`, `"nonregion"`, `"region"`,
#'   `"region_Gstar"` (region outside H*), `"Hstar"`, `"Hstar_region"`,
#'   `"Hstar_nonregion"`; or a logical vector of length `pop$N`.
#' @param b,h Hazard ratios outside and inside H*: the treated log hazard is
#'   \eqn{\theta_0 + \log b + \log(h/b) 1\{H^*\}}. Both `NULL` (and `psi`
#'   `NULL`): the DGM's own per-patient log hazard ratios (`loghr_po`).
#' @param psi `NULL`, or a per-patient log hazard ratio vector of length
#'   `pop$N` (overrides `b`, `h`, `b_R`, `h_R`).
#' @param horizon Calendar time of the analysis; `NULL` (default) is
#'   `pop$x`; `Inf` gives the uncensored limit.
#' @param variance Logical: compute the robust variance (default `TRUE`).
#' @param b_R,h_R Hazard ratios inside the region, outside and inside H*:
#'   with `b` and `h` they give the four-cell treated log hazard
#'   \eqn{\theta_0 + \log b} (outside H*, outside the region),
#'   \eqn{\log h} (H*, outside), \eqn{\log b_R} (outside H*, region),
#'   \eqn{\log h_R} (H*, region). `NULL` (default) means `b_R = b` and
#'   `h_R = h`, the two-cell model. Used only with `b` and `h`. When they
#'   differ from `b`, `h`, a type group whose cells carry different hazard
#'   ratios has `psi = NA` in `types`.
#' @return A list: `beta` (the limit) and `hr` (margHR at the horizon),
#'   `i1` (information per patient), `var_model`, `var_robust`,
#'   `var_robust_fixed` (allocation fixed at n/2) (variances per patient, so
#'   divide by the set's size), `sandwich`, `d1`, `d0`, `d` (event
#'   probabilities by arm and at 1:1, i.e. expected events per patient),
#'   `types` (per type group: patients, log-HR, events, event share, own and
#'   pooled imbalance moments), `omega`, `ahr`, `cde`, `band_share`, `omega0`,
#'   `marghr_inf`, `n_pop`, `share` (\eqn{n_S/N}), `horizon`, and the node
#'   arrays and `grid` used by [ahr_design()].
#' @references Struthers, C. A. and Kalbfleisch, J. D. (1986). Misspecified
#'   proportional hazards models. *Biometrika* 73, 363-369.
#' @seealso [ahr_population()], [ahr_design()]
#' @examples
#' \dontrun{
#' pop <- ahr_population(dgm, 40, 11)
#' e <- ahr_limit(pop, "region", b = 0.55, h = 2)
#' c(margHR40 = e$hr, AHR = e$ahr, CDE = e$cde, margHR_inf = e$marghr_inf)
#' }
#' @export
ahr_limit <- function(pop, set, b = NULL, h = NULL, psi = NULL, horizon = NULL,
                      variance = TRUE, b_R = NULL, h_R = NULL) {
  if (!inherits(pop, "ahr_population"))
    stop("pop must be an ahr_population().", call. = FALSE)
  if (xor(is.null(b), is.null(h)))
    stop("Give both b and h, or neither.", call. = FALSE)
  if ((!is.null(b_R) || !is.null(h_R)) && is.null(b))
    stop("b_R and h_R need b and h.", call. = FALSE)
  if (is.null(b_R)) b_R <- b
  if (is.null(h_R)) h_R <- h
  if (is.null(horizon)) horizon <- pop$x
  if (!is.numeric(horizon) || length(horizon) != 1L || horizon <= 0)
    stop("horizon must be one positive number (or Inf).", call. = FALSE)
  out <- .ahr_exact(pop, set, b, h, psi, horizon, variance, b_R = b_R, h_R = h_R)
  cz <- .ahr_causal(pop, out$S, out$psi_v)
  out$ahr <- cz[["ahr"]]; out$cde <- cz[["cde"]]
  out$band_share <- cz[["band_share"]]; out$omega0 <- cz[["omega0"]]
  out$marghr_inf <- if (is.infinite(horizon)) out$hr else
    .ahr_exact(pop, set, b, h, psi, Inf, FALSE, b_R = b_R, h_R = h_R)$hr
  out$share <- out$n_pop / pop$N
  out$horizon <- horizon
  out$S <- NULL; out$psi_v <- NULL
  out
}

## Asymptotic correlation of the estimators of two nested sets (small within big),
## same grid: Corr = C sqrt(pi) / sqrt(Sigma_big Sigma_small).
.ahr_corr <- function(big, small, w) {
  C <- (sum(w * small$aH * big$f1 * small$f1) + sum(w * small$AP * big$g1 * small$g1) +
        sum(w * small$bH * big$f0 * small$f0) + sum(w * small$BP * big$g0 * small$g0)) / 2
  pi <- small$n_pop / big$n_pop
  c(exact = C * sqrt(pi) / sqrt(big$sandwich * small$sandwich),
    se_ratio = sqrt(pi * small$i1 / big$i1), pi = pi)
}

#' Design quantities and trigger probabilities at a trial size
#'
#' At trial size `n`, for the standard sets of an [ahr_population()] with
#' hazard ratio `b` outside H* and `h` inside: each set's limit and its
#' model-based and robust SEs (using \eqn{E[1/n_S]},
#' \eqn{n_S \sim Bin(n, \pi_S)}), expected events, the exact and SE-ratio
#' correlations between nested sets' estimates, and, from the trivariate
#' normal approximation (estimates \eqn{\sim N(\beta^*, SE_{robust}^2)},
#' exact correlations, the Wald limit of the ITT with the model SE):
#' \itemize{
#'   \item `p_itt_ub_lt1`: P(ITT upper 95\% limit < 1);
#'   \item `p_region_itt_gt1`: P(region ITT estimate > 1);
#'   \item `p_region_Gstar_lt1`: P(region(G*) estimate < 1);
#'   \item `p_trigger`: P(ITT ub < 1 and region ITT > 1);
#'   \item `p_trigger_and_Gstar_lt1`: P(trigger and region(G*) < 1), the
#'     oracle joint.
#' }
#' Probabilities use `mvtnorm::pmvnorm(algorithm = Miwa(steps = 1024))`.
#'
#' @inheritSection ahr_population Limits
#' @param pop An [ahr_population()].
#' @param b,h Hazard ratios outside and inside H*. All sets use the
#'   population's cached grid (built for \eqn{0.2 \le b, h \le 30}; outside
#'   that range the grid's first panel is coarser than the library's).
#' @param n Trial size (both arms).
#' @param sets Standard set names to report (see [ahr_limit()]); the trigger
#'   sets `"all"`, `"region"` and `"region_Gstar"` are always computed.
#' @param b_R,h_R Hazard ratios inside the region, outside and inside H*
#'   (the four-cell model of [ahr_limit()]); default `b_R = b`, `h_R = h`,
#'   the two-cell model.
#' @return A list: `sets` (a data.frame per set: `set`, `share`, `beta`,
#'   `marghr` (at the analysis time), `se_model`, `se_robust`, `se_d4`
#'   (\eqn{\sqrt{4/events}}), `info_ratio`, `events` (expected, of `n`),
#'   `ahr`, `cde`), `corr` (a data.frame per nested pair: `big`, `small`,
#'   `exact`, `se_ratio`, `pi`), `prob` (the named probabilities above), and
#'   `b`, `h`, `n`, `b_R`, `h_R`.
#' @seealso [ahr_limit()], [ahr_solve()]
#' @examples
#' \dontrun{
#' pop <- ahr_population(dgm, 40, 11)
#' ahr_design(pop, b = 0.55, h = 2, n = 1000)$prob
#' }
#' @importFrom mvtnorm pmvnorm Miwa
#' @importFrom stats pnorm qnorm uniroot
#' @export
ahr_design <- function(pop, b, h, n, sets = c("all", "region", "region_Gstar"),
                       b_R = b, h_R = h) {
  if (!inherits(pop, "ahr_population"))
    stop("pop must be an ahr_population().", call. = FALSE)
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) || n < 1)
    stop("n must be one positive number.", call. = FALSE)
  if (!is.character(sets) || !all(sets %in% names(.ahr_set_cells)))
    stop("sets must be among ", paste(names(.ahr_set_cells), collapse = ", "), ".",
         call. = FALSE)
  trig <- c("all", "region", "region_Gstar")
  nm <- unique(c(sets, trig))
  e <- lapply(stats::setNames(nm, nm), function(st)
    .ahr_exact(pop, st, b, h, NULL, pop$x, TRUE, force_fast = TRUE, b_R = b_R, h_R = h_R))
  w <- e[[1]]$grid$w
  q <- lapply(e, function(x) {
    pi <- x$n_pop / pop$N
    inv_n <- (1 + (1 - pi) / (n * pi)) / (n * pi)      # E[1 / n_S], n_S ~ Bin(n, pi)
    list(beta = x$beta, pi = pi, se_model = sqrt(x$var_model * inv_n),
         se_robust = sqrt(x$var_robust * inv_n), se_d4 = sqrt(4 / x$d * inv_n),
         info_ratio = x$i1 / (x$d / 4), events = n * pi * x$d)
  })
  tab <- do.call(rbind, lapply(sets, function(st) {
    cz <- .ahr_causal(pop, e[[st]]$S, e[[st]]$psi_v)
    data.frame(set = st, share = q[[st]]$pi, beta = q[[st]]$beta, marghr = exp(q[[st]]$beta),
               se_model = q[[st]]$se_model, se_robust = q[[st]]$se_robust,
               se_d4 = q[[st]]$se_d4, info_ratio = q[[st]]$info_ratio,
               events = q[[st]]$events, ahr = cz[["ahr"]], cde = cz[["cde"]],
               stringsAsFactors = FALSE)
  }))
  ## nested pairs (small within big) among the computed sets
  prs <- list()
  for (bg in nm) for (sm in nm) if (bg != sm &&
      all(.ahr_set_cells[[sm]] %in% .ahr_set_cells[[bg]]))
    prs[[length(prs) + 1]] <- c(big = bg, small = sm)
  corr <- do.call(rbind, lapply(prs, function(p) {
    cc <- .ahr_corr(e[[p[["big"]]]], e[[p[["small"]]]], w)
    data.frame(big = p[["big"]], small = p[["small"]], exact = cc[["exact"]],
               se_ratio = cc[["se_ratio"]], pi = cc[["pi"]], stringsAsFactors = FALSE)
  }))
  rho <- function(bg, sm) corr$exact[corr$big == bg & corr$small == sm]
  ## probabilities (the library's design_prob(), exact correlations)
  z975 <- stats::qnorm(0.975)
  mu <- c(q$all$beta, q$region$beta, q$region_Gstar$beta)
  sd <- c(q$all$se_robust, q$region$se_robust, q$region_Gstar$se_robust)
  R <- diag(3)
  R[1, 2] <- R[2, 1] <- rho("all", "region")
  R[1, 3] <- R[3, 1] <- rho("all", "region_Gstar")
  R[2, 3] <- R[3, 2] <- rho("region", "region_Gstar")
  cutA <- -z975 * q$all$se_model
  zA <- (cutA - mu[1]) / sd[1]; zB <- (log(1) - mu[2]) / sd[2]; zC <- (0 - mu[3]) / sd[3]
  pm <- function(lo, up, idx) as.numeric(mvtnorm::pmvnorm(
    lower = lo, upper = up, corr = R[idx, idx], algorithm = mvtnorm::Miwa(steps = 1024)))
  prob <- c(p_itt_ub_lt1 = stats::pnorm(zA), p_region_itt_gt1 = 1 - stats::pnorm(zB),
            p_region_Gstar_lt1 = stats::pnorm(zC),
            p_trigger = pm(c(-Inf, zB), c(zA, Inf), 1:2),
            p_trigger_and_Gstar_lt1 = pm(c(-Inf, zB, -Inf), c(zA, Inf, zC), 1:3))
  list(sets = tab, corr = corr, prob = prob, b = b, h = h, n = n, b_R = b_R, h_R = h_R)
}

#' Solve P(trigger) = target in one free parameter
#'
#' Scans `f` over a grid of the free parameter, brackets every crossing of
#' `target`, and refines each with [stats::uniroot()]: the pattern of the
#' mrct-consistency boundary and regional-entry scripts. `f` maps the free
#' parameter to the probability (typically
#' `ahr_design(...)$prob[["p_trigger"]]`); when the parameter enters the build
#' (the region prognostic log-HR `r`, the entry delay), `f` builds its own
#' population. Roots are only found where the scan brackets them: two
#' crossings between adjacent scan points are missed.
#'
#' @param f Function of the free parameter (first argument) returning one
#'   number.
#' @param target The value to solve for.
#' @param interval Length-2 range of the free parameter.
#' @param ... Further arguments passed to `f`.
#' @param n_scan Number of equally spaced scan points on `interval`
#'   (default 25); ignored when `scan` is given.
#' @param scan Optional explicit scan points (sorted, within `interval`).
#' @param tol Tolerance of [stats::uniroot()] in the free parameter.
#' @return A list: `roots` (sorted), `n_roots`, `values` (`f` at the roots),
#'   `scan` (a data.frame of the scan points `x` and `value`), and
#'   `evaluations` (calls to `f`).
#' @seealso [ahr_design()]
#' @examples
#' \dontrun{
#' f <- function(a) ahr_design(ahr_population(dgm, 40, 11,
#'   entry_delay = list(var = "z_region", delay = a)), 0.55, 2, 1000)$prob[["p_trigger"]]
#' ahr_solve(f, 0.25, c(0, 24), scan = 0:24)
#' }
#' @export
ahr_solve <- function(f, target, interval, ..., n_scan = 25L, scan = NULL, tol = 1e-8) {
  if (!is.function(f)) stop("f must be a function.", call. = FALSE)
  if (!is.numeric(interval) || length(interval) != 2L || !(interval[1] < interval[2]))
    stop("interval must be c(lower, upper) with lower < upper.", call. = FALSE)
  xs <- as.numeric(if (is.null(scan)) seq(interval[1], interval[2], length.out = n_scan) else sort(scan))
  if (length(xs) < 2L) stop("The scan needs at least two points.", call. = FALSE)
  ne <- 0L
  ff <- function(x) {
    ne <<- ne + 1L
    v <- f(x, ...)
    if (!is.numeric(v) || length(v) != 1L || is.na(v))
      stop("f must return one number (got NA or other at ", x, ").", call. = FALSE)
    v - target
  }
  d <- vapply(xs, ff, 0)
  roots <- xs[d == 0]
  k <- length(d)
  j <- which(sign(d[-k]) * sign(d[-1]) < 0)
  for (jj in j) roots <- c(roots, stats::uniroot(ff, xs[jj + 0:1], f.lower = d[jj],
                                                f.upper = d[jj + 1], tol = tol)$root)
  roots <- sort(roots)
  vals <- vapply(roots, function(x) f(x, ...), 0)
  list(roots = roots, n_roots = length(roots), values = vals,
       scan = data.frame(x = xs, value = d + target), evaluations = ne + length(roots))
}


## -----------------------------------------------------------------------------
## 5. The approximation ladder (internal; from sk_limit.R's ahr_ladder())
## -----------------------------------------------------------------------------

## Type form: Slide 9 over types with the band's event share, one- and two-term
## corrections with the type's own and the pooled imbalance; time form: on the
## composition-induced profile theta_S(u); events form: the score read as
## D1 / D = E_N expit(beta + delta), expanded about the mean imbalance.
.ahr_ladder <- function(e) {
  expit <- function(z) 1 / (1 + exp(-z)); logit <- function(p) log(p / (1 - p))
  ty <- e$types; d <- sum(ty$events); p <- expit(ty$psi); v <- p * (1 - p); pk <- ty$events / d
  uroot <- function(G) tryCatch(stats::uniroot(G, c(1e-9, 1 - 1e-9), tol = 1e-13)$root,
                                error = function(err) NA_real_)
  P0 <- sum(pk * p)
  m1o <- sum(v * ty$m1_own); m1p <- sum(ty$m1_pool)
  m2o <- sum(v * (1 - 2 * p) * ty$m2_own); m2p <- sum(ty$m2_pool)
  P1 <- P0 + (m1o - P0 * (1 - P0) * m1p) / d
  G2 <- function(P) d * (P0 - P) + m1o - P * (1 - P) * m1p +
    0.5 * (m2o - P * (1 - P) * (1 - 2 * P) * m2p)
  P2 <- uroot(G2)
  dbar_own <- ty$m1_own / ty$events; dbar_pool <- ty$m1_pool / ty$events; dbar <- m1p / d
  dec <- c(mismatch = sum(pk * v * (dbar_own - dbar_pool)),
           cov = sum(pk * (v - sum(pk * v)) * (dbar_pool - dbar)),
           penalty = -sum(pk * (p - P0)^2) * dbar)
  pt <- expit(e$loghr_t); vt <- pt * (1 - pt); dN <- e$dN; dl <- e$delta
  T0 <- sum(pt * dN) / d
  T1 <- T0 + sum((vt - T0 * (1 - T0)) * dl * dN) / d
  T2 <- uroot(function(P) sum((pt - P) * dN) + sum((vt - P * (1 - P)) * dl * dN) +
                0.5 * sum((vt * (1 - 2 * pt) - P * (1 - P) * (1 - 2 * P)) * dl^2 * dN))
  Pb <- e$d1 / (e$d1 + e$d0); vard <- sum(dl^2 * dN) / d - dbar^2
  E0 <- logit(Pb); E1 <- E0 - dbar
  pc <- uroot(function(q) q + 0.5 * q * (1 - q) * (1 - 2 * q) * vard - Pb)
  E2 <- logit(pc) - dbar
  list(beta = c(exact = e$beta, slide9_types = logit(P0), type_1term = logit(P1),
                type_2term = logit(P2), slide9_time = logit(T0), time_1term = logit(T1),
                time_2term = logit(T2), events_ratio = E0, events_1term = E1,
                events_2term = E2),
       decomposition = c(dP = P1 - P0, dec), dbar = dbar, var_delta = vard,
       dbar_own = stats::setNames(dbar_own, ty$type),
       dbar_pool = stats::setNames(dbar_pool, ty$type),
       omega = if ("band" %in% ty$type) pk[ty$type == "band"] else 0, P0 = P0)
}
