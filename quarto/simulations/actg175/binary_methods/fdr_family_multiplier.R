# =============================================================================
# fdr_family_multiplier.R
#
# PROTOTYPE / SKETCH -- family-level false-discovery (type-I) approximation for
# forest-search subgroup identification, by a multiplier bootstrap over the whole
# candidate family.  This is the computational companion to Section 2.2 /
# Equation (eq-fdr) of the resampling-theory note: it estimates
#
#     Pr( UNION_g { beta_g >= c_screen  AND  consistency_g >= p_star } | H0 )
#
# the probability that the search flags *any* of the S overlapping candidates
# under a null configuration -- a union over a strongly correlated family that a
# per-candidate formula (Leon et al. 2024, eq. 3) cannot capture.
#
# The correlation is supplied by the per-subject influence (dfbeta) contributions:
# all candidates are perturbed by the SAME multiplier vector propagated through a
# shared N x S influence matrix, so overlap among candidates (shared subjects) is
# respected exactly, with one matrix multiply per batch of draws and no refitting.
#
# Base R only (stats).  Functions are written package-ready (roxygen) so they can
# graduate into forestsearch; the demo at the bottom is dev-only.
#
# Usage:
#   source("fdr_family_multiplier.R")
#   run_fdr_family_demo()              # self-contained simulated illustration
# =============================================================================


#' Draw mean-zero, unit-variance multipliers
#'
#' @param n,draws Matrix dimensions (subjects x draws).
#' @param type One of `"rademacher"` (the 50/50 split), `"gaussian"`, or
#'   `"poisson"` (centred Poisson, `Pois(1) - 1`).
#' @return An `n x draws` numeric matrix.
#' @keywords internal
.draw_multipliers <- function(n, draws, type = c("rademacher", "gaussian", "poisson")) {
  type <- match.arg(type)
  m <- n * draws
  g <- switch(
    type,
    rademacher = sample(c(-1, 1), m, replace = TRUE),
    gaussian   = stats::rnorm(m),
    poisson    = stats::rpois(m, 1) - 1   # mean 0, variance 1
  )
  matrix(g, nrow = n, ncol = draws)
}


#' Assemble the N x S influence matrix from per-candidate dfbeta pieces
#'
#' Each candidate subgroup contributes its treatment-coefficient `dfbeta` values
#' at the rows (subjects) it contains; all other rows are zero.  The resulting
#' column `g` is the influence representation of candidate `g`'s effect estimate,
#' and `crossprod()` of two columns is the cross-candidate covariance that drives
#' the family-level calculation.
#'
#' @param candidates A list with one element per candidate subgroup; element `g`
#'   is `list(idx = <integer subject indices into 1:N>, dfbeta = <numeric, same
#'   length as idx>)`.  In `forestsearch` these are exactly the pieces returned by
#'   the consistency engine (`residuals(coxfit, type = "dfbeta")[, treat]` for
#'   Cox, or `stats::dfbeta(glmfit)[, treat]` for GLM), restricted to the
#'   candidate's members.
#' @param N Total number of subjects.
#' @param sparse Logical; if `TRUE` and the `Matrix` package is available, return
#'   a sparse `dgCMatrix` (recommended for many candidates). Default `FALSE`.
#' @return An `N x S` (dense or sparse) influence matrix.
#' @export
assemble_influence <- function(candidates, N, sparse = FALSE) {
  S <- length(candidates)
  if (S == 0L) stop("`candidates` is empty.")
  if (sparse && requireNamespace("Matrix", quietly = TRUE)) {
    i_list <- integer(0); j_list <- integer(0); x_list <- numeric(0)
    for (g in seq_len(S)) {
      cg <- candidates[[g]]
      if (length(cg$idx) != length(cg$dfbeta))
        stop(sprintf("candidate %d: `idx` and `dfbeta` length mismatch.", g))
      i_list <- c(i_list, cg$idx)
      j_list <- c(j_list, rep.int(g, length(cg$idx)))
      x_list <- c(x_list, cg$dfbeta)
    }
    return(Matrix::sparseMatrix(i = i_list, j = j_list, x = x_list, dims = c(N, S)))
  }
  B <- matrix(0, nrow = N, ncol = S)
  for (g in seq_len(S)) {
    cg <- candidates[[g]]
    if (length(cg$idx) != length(cg$dfbeta))
      stop(sprintf("candidate %d: `idx` and `dfbeta` length mismatch.", g))
    B[cg$idx, g] <- cg$dfbeta
  }
  B
}


#' Family-level false-discovery approximation via multiplier bootstrap
#'
#' Approximates the search-level false-discovery (type-I) probability of the
#' forest-search screening-and-consistency rule by a multiplier bootstrap over the
#' entire candidate family, reusing the per-subject influence contributions
#' instead of refitting on simulated null data sets.
#'
#' A candidate `g` is *flagged* when its (perturbed) effect estimate clears both
#' the screening threshold and the closed-form consistency cutoff,
#' \deqn{\hat\beta_g \ge c_{\mathrm{screen}} \quad\text{and}\quad
#'       2\Phi((\hat\beta_g - c_{\mathrm{cons}})/\sigma_{D,g}) - 1 \ge p^\star,}
#' which is equivalent to \eqn{\hat\beta_g \ge t_g} with
#' \eqn{t_g = \max(c_{\mathrm{screen}},\ c_{\mathrm{cons}} + z\,\sigma_{D,g})} and
#' \eqn{z = \Phi^{-1}((1 + p^\star)/2)}.
#'
#' For each draw a single multiplier vector \eqn{\xi^{(b)}} is propagated to every
#' candidate through the shared influence matrix, giving
#' \eqn{\hat\beta_g^{*(b)} = \beta_{0,g} + \sum_i \xi_i^{(b)}\,\mathrm{db}_{g,i}}
#' (column `g` of `influence`).  The family discovery indicator is
#' \eqn{\bigl[\exists g:\ \hat\beta_g^{*(b)} \ge t_g\bigr]}, and the returned
#' `fdr` is its average over draws.
#'
#' @param influence `N x S` influence matrix (see [assemble_influence()]); column
#'   `g` holds candidate `g`'s treatment-coefficient dfbeta, zero outside `g`.
#' @param beta_hat Length-`S` observed candidate effect estimates on the
#'   natural-parameter (log) scale. Used only to report the observed flag set and,
#'   when `center = "observed"`, as the bootstrap centre.
#' @param sigma_D Length-`S` split-perturbation SDs. `NULL` (default) uses
#'   `sqrt(colSums(influence^2))`, the robust SE of each candidate coefficient.
#' @param c_screen Screening threshold on the comparison scale, e.g. `log(1.25)`.
#' @param c_consistency Consistency threshold on the comparison scale, e.g.
#'   `log(1.0) = 0`. Default `0`.
#' @param p_star Consistency-rate cutoff. Default `0.90`.
#' @param beta_null Scalar or length-`S` null configuration (log scale) at which
#'   the family is centred for the type-I/FDR calculation, e.g. the log ITT effect
#'   for a uniform-benefit null, or `0` for the boundary null. Default `0`.
#' @param center `"null"` (default) centres the bootstrap at `beta_null` (type-I /
#'   FDR); `"observed"` centres at `beta_hat` (sampling variability / stability of
#'   the discovery set around the observed configuration).
#' @param draws Number of multiplier draws. Default `2000`.
#' @param multiplier `"rademacher"` (the 50/50 split), `"gaussian"`, or
#'   `"poisson"`. Default `"rademacher"`.
#' @param seed Optional integer seed.
#'
#' @return A list with: `fdr` (family discovery probability), `e_num_flagged`
#'   (expected number of candidates flagged per draw), `num_flagged_quantiles`,
#'   `marginal` (per-candidate closed-form pass probability at the null),
#'   `union_independent` (independence union bound, for contrast), `thresholds`
#'   (`t_g`), `n_observed_flagged`, and metadata.
#' @export
fdr_family_multiplier <- function(influence,
                                  beta_hat,
                                  sigma_D = NULL,
                                  c_screen,
                                  c_consistency = 0,
                                  p_star = 0.90,
                                  beta_null = 0,
                                  center = c("null", "observed"),
                                  draws = 2000L,
                                  multiplier = c("rademacher", "gaussian", "poisson"),
                                  seed = NULL) {
  center     <- match.arg(center)
  multiplier <- match.arg(multiplier)
  if (!is.null(seed)) set.seed(seed)

  influence <- as.matrix(influence)
  N <- nrow(influence); S <- ncol(influence)
  if (length(beta_hat) != S) stop("`beta_hat` must have length ncol(influence).")
  if (is.null(sigma_D)) sigma_D <- sqrt(colSums(influence^2))
  if (length(sigma_D) != S) stop("`sigma_D` must have length ncol(influence).")
  if (any(sigma_D <= 0)) stop("`sigma_D` must be positive for every candidate.")
  if (length(beta_null) == 1L) beta_null <- rep(beta_null, S)
  if (length(beta_null) != S) stop("`beta_null` must be scalar or length ncol(influence).")

  # Per-candidate flag threshold on the (perturbed) estimate scale.
  z   <- stats::qnorm((1 + p_star) / 2)
  t_g <- pmax(c_screen, c_consistency + z * sigma_D)        # length S

  centre <- if (center == "null") beta_null else beta_hat   # length S

  # Multiplier bootstrap: one draw -> all candidates at once.
  #   P[g, b] = sum_i influence[i, g] * Xi[i, b]   (S x draws);  Var = sigma_D[g]^2.
  Xi <- .draw_multipliers(N, draws, multiplier)             # N x draws
  P  <- crossprod(influence, Xi)                            # S x draws
  beta_star <- centre + P                                   # recycles `centre` down columns

  flagged   <- beta_star >= t_g                             # recycles `t_g` down columns
  num_flag  <- colSums(flagged)
  discovery <- num_flag > 0L

  # Per-candidate marginal pass probability at the null (closed form) and the
  # independence union bound, for contrast with the correlated family answer.
  marginal          <- stats::pnorm((t_g - beta_null) / sigma_D, lower.tail = FALSE)
  union_independent <- 1 - prod(1 - marginal)

  list(
    fdr                   = mean(discovery),
    e_num_flagged         = mean(num_flag),
    num_flagged_quantiles = stats::quantile(num_flag, c(.5, .9, .95, .99, 1)),
    marginal              = marginal,
    union_independent     = union_independent,
    thresholds            = t_g,
    n_observed_flagged    = sum(beta_hat >= t_g),
    center                = center,
    multiplier            = multiplier,
    draws                 = as.integer(draws),
    n_subjects            = N,
    n_candidates          = S
  )
}


# =============================================================================
# DEV-ONLY DEMO  (not exported; illustrates the routine on simulated null data)
# =============================================================================

#' Self-contained demonstration on a simulated null (no harm subgroup)
#'
#' Generates a randomized trial with a *homogeneous* protective treatment effect
#' (so no subgroup is truly harmful), builds a family of overlapping one- and
#' two-factor candidate subgroups, fits a treatment-only logistic model in each to
#' obtain the log-OR and its dfbeta, and runs [fdr_family_multiplier()].  Prints
#' the correlated family false-discovery probability alongside the per-candidate
#' marginal and the independence union bound to show (i) the multiplicity inflation
#' relative to any single candidate and (ii) the reduction relative to a naive
#' independence assumption once candidate overlap is accounted for.
#'
#' @param N Trial size. @param draws Multiplier draws. @param seed RNG seed.
#' @keywords internal
run_fdr_family_demo <- function(N = 1000L, draws = 4000L, seed = 1L) {
  set.seed(seed)
  trt <- rbinom(N, 1, 0.5)
  z1 <- rbinom(N, 1, 0.5); z2 <- rbinom(N, 1, 0.5)
  x1 <- rnorm(N); x2 <- rnorm(N); x3 <- rnorm(N)

  # Homogeneous protective effect: treat log-OR = -0.3 for everyone (NULL: no harm).
  lp <- -0.2 - 0.3 * trt + 0.4 * z1 - 0.3 * z2 + 0.25 * x1
  y  <- rbinom(N, 1, plogis(lp))
  dat <- data.frame(y, trt, z1, z2,
                    x1lo = as.integer(x1 <= median(x1)),
                    x2lo = as.integer(x2 <= median(x2)),
                    x3lo = as.integer(x3 <= median(x3)))

  facs <- list(
    "z1=1"     = dat$z1 == 1,   "z1=0"     = dat$z1 == 0,
    "z2=1"     = dat$z2 == 1,   "z2=0"     = dat$z2 == 0,
    "x1<=med"  = dat$x1lo == 1, "x1>med"   = dat$x1lo == 0,
    "x2<=med"  = dat$x2lo == 1, "x2>med"   = dat$x2lo == 0,
    "x3<=med"  = dat$x3lo == 1, "x3>med"   = dat$x3lo == 0
  )
  fn <- names(facs)

  # Single- and two-factor candidate definitions.
  defs <- list()
  for (a in seq_along(facs)) defs[[length(defs) + 1L]] <- list(name = fn[a], mem = facs[[a]])
  for (a in seq_along(facs)) for (b in seq_along(facs)) if (b > a) {
    defs[[length(defs) + 1L]] <- list(name = paste(fn[a], "&", fn[b]),
                                      mem  = facs[[a]] & facs[[b]])
  }

  min_n <- 60L
  candidates <- list(); beta_hat <- numeric(0); keep <- character(0)
  for (cd in defs) {
    idx <- which(cd$mem)
    if (length(idx) < min_n) next
    d <- dat[idx, , drop = FALSE]
    if (length(unique(d$trt)) < 2L) next
    ev1 <- sum(d$y[d$trt == 1]); ev0 <- sum(d$y[d$trt == 0])
    if (ev1 %in% c(0L, sum(d$trt == 1)) || ev0 %in% c(0L, sum(d$trt == 0))) next  # avoid separation
    fit <- try(suppressWarnings(glm(y ~ trt, family = binomial(), data = d)), silent = TRUE)
    if (inherits(fit, "try-error")) next
    cf <- stats::coef(fit)[["trt"]]
    if (is.na(cf)) next
    db <- try(stats::dfbeta(fit)[, "trt"], silent = TRUE)
    if (inherits(db, "try-error") || anyNA(db)) next
    candidates[[length(candidates) + 1L]] <- list(idx = idx, dfbeta = as.numeric(db))
    beta_hat <- c(beta_hat, as.numeric(cf)); keep <- c(keep, cd$name)
  }

  B   <- assemble_influence(candidates, N)
  itt <- stats::coef(glm(y ~ trt, family = binomial(), data = dat))[["trt"]]

  res <- fdr_family_multiplier(
    influence     = B,
    beta_hat      = beta_hat,
    c_screen      = log(1.25),
    c_consistency = log(1.0),
    p_star        = 0.90,
    beta_null     = itt,          # uniform-benefit null centre (log ITT effect)
    center        = "null",
    draws         = draws,
    multiplier    = "rademacher",
    seed          = seed
  )

  cat("---- Family-level FDR (multiplier bootstrap) : simulated NULL ----\n")
  cat(sprintf("Subjects N                          : %d\n", N))
  cat(sprintf("Candidate subgroups evaluated S     : %d\n", length(candidates)))
  cat(sprintf("ITT log-OR (null centre)            : %+.3f  (OR %.2f)\n", itt, exp(itt)))
  cat(sprintf("Multiplier / draws                  : %s / %d\n", res$multiplier, res$draws))
  cat("-----------------------------------------------------------------\n")
  cat(sprintf("Max per-candidate marginal flag prob: %.3f\n", max(res$marginal)))
  cat(sprintf("Independence union bound            : %.3f\n", res$union_independent))
  cat(sprintf("Family FDR (correlated, this tool)  : %.3f   <-- search-level type-I\n", res$fdr))
  cat(sprintf("Expected # candidates flagged/draw  : %.2f\n", res$e_num_flagged))
  cat(sprintf("Candidates flagged in OBSERVED data : %d\n", res$n_observed_flagged))
  cat("-----------------------------------------------------------------\n")
  cat("Note: family FDR exceeds any single-candidate marginal (multiplicity)\n",
      "but is below the independence bound (overlapping candidates are\n",
      "positively correlated, so flags cluster rather than multiply).\n", sep = "")
  invisible(res)
}

# Run the demo when executed via `Rscript fdr_family_multiplier.R`.
if (identical(environment(), globalenv()) && !interactive() &&
    sys.nframe() == 0L && length(commandArgs(trailingOnly = TRUE)) == 0L) {
  run_fdr_family_demo()
}
