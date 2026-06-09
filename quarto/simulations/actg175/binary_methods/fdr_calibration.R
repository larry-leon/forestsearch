# =============================================================================
# fdr_calibration.R
#
# PROTOTYPE / SKETCH -- calibrate the forest-search screening / consistency
# thresholds to a TARGET search-level false-discovery (type-I) rate, using the
# family-level multiplier bootstrap of fdr_family_multiplier().
#
# Motivation.  With c_screen = log(1.25) and c_consistency = log(1.0) the
# *per-candidate* type-I is controlled (~5%, Leon et al. 2024 eq. 3), but the
# *search-level* rate -- the chance that ANY of the many overlapping candidates
# is flagged -- can be far higher (e.g. ~0.44 for the FS arm in the ACTG175
# binary/OR null simulation).  That union over a correlated family is exactly the
# object fdr_family_multiplier() approximates, so we can read the search-level
# FDR off ONE null candidate family and sweep thresholds to hit a target without
# re-running the 1000-simulation null loop.
#
# Pipeline:
#   1. enumerate_candidates()    build the one/two-factor candidate family from
#                                the SAME binary factor indicators the search uses
#   2. build_influence_matrix()  re-derive each candidate's treatment dfbeta via
#                                the engine's own .consistency_*_pieces() (so the
#                                influence contributions match the real run)
#   3. calibrate_fdr_grid()      sweep (c_screen, c_consistency, p_star) and read
#                                off the family FDR for each setting
#
# Dependencies: fdr_family_multiplier.R (sibling) and, for the dfbeta, either an
# installed `forestsearch` (uses forestsearch:::.consistency_*_pieces) or
# consistency_resample.R on the path.  Base R + stats otherwise.
#
# Usage:
#   source("fdr_calibration.R")
#   run_fdr_calibration_demo()
# =============================================================================

if (!exists("fdr_family_multiplier")) {
  .fm <- getOption("fdr_family_src", "fdr_family_multiplier.R")
  if (file.exists(.fm)) source(.fm) else
    stop("fdr_family_multiplier() not found; place fdr_family_multiplier.R on the path.")
}


#' Resolve the engine's consistency-pieces functions (installed pkg or dev file)
#' @keywords internal
.resolve_pieces_fns <- function() {
  if (requireNamespace("forestsearch", quietly = TRUE) &&
      exists(".consistency_glm_pieces", envir = asNamespace("forestsearch"), inherits = FALSE)) {
    ns <- asNamespace("forestsearch")
    return(list(glm = get(".consistency_glm_pieces", ns),
                cox = get(".consistency_cox_pieces", ns)))
  }
  if (!exists(".consistency_glm_pieces", inherits = TRUE)) {
    src <- getOption("fdr_consistency_src", "consistency_resample.R")
    if (file.exists(src)) source(src) else
      stop("Install forestsearch, or put consistency_resample.R on the path ",
           "(or set options(fdr_consistency_src = ...)).")
  }
  list(glm = get(".consistency_glm_pieces"), cox = get(".consistency_cox_pieces"))
}

#' Call a pieces function with only the arguments it accepts (robust to signature)
#' @keywords internal
.call_pieces <- function(fn, dfg, args) {
  fa   <- names(formals(fn))
  args <- args[names(args) %in% fa]
  do.call(fn, c(list(dfg), args))
}


#' Enumerate one- and two-factor candidate subgroups from factor indicators
#'
#' Mirrors Steps 1--2 of the forest search: single-factor subgroups and their
#' two-factor intersections, restricted to a minimum size (and optionally a
#' minimum event count per arm).  Pass the SAME binary factor-level indicators the
#' search used (one logical vector per factor level, length `nrow(data)`); the
#' result is a list of integer row-index vectors, one per candidate.
#'
#' @param data Analysis data frame.
#' @param factors Named list of logical vectors (binary factor-level indicators).
#' @param maxk Maximum factors per candidate (1 or 2). Default 2.
#' @param min_n Minimum candidate size. Default 60.
#' @param min_events_per_arm Optional minimum events in each treatment arm.
#' @param treat.name,outcome.name,event.name Column names used for the event
#'   guard (outcome.name for GLM, event.name for survival).
#' @return Named list of integer index vectors (candidate memberships).
#' @export
enumerate_candidates <- function(data, factors, maxk = 2L, min_n = 60L,
                                 min_events_per_arm = NULL, treat.name = NULL,
                                 outcome.name = NULL, event.name = NULL) {
  stopifnot(is.list(factors), length(factors) >= 1L)
  fn <- names(factors); if (is.null(fn)) fn <- paste0("f", seq_along(factors))

  ev_ok <- function(idx) {
    if (length(idx) < min_n) return(FALSE)
    if (is.null(min_events_per_arm)) return(TRUE)
    z  <- data[[treat.name]][idx]
    yv <- if (!is.null(event.name)) data[[event.name]][idx] else data[[outcome.name]][idx]
    sum(yv[z == 1], na.rm = TRUE) >= min_events_per_arm &&
      sum(yv[z == 0], na.rm = TRUE) >= min_events_per_arm
  }

  cands <- list(); nms <- character(0)
  for (a in seq_along(factors)) {
    idx <- which(factors[[a]])
    if (ev_ok(idx)) { cands[[length(cands) + 1L]] <- idx; nms <- c(nms, fn[a]) }
  }
  if (maxk >= 2L) for (a in seq_along(factors)) for (b in seq_along(factors)) if (b > a) {
    idx <- which(factors[[a]] & factors[[b]])
    if (ev_ok(idx)) { cands[[length(cands) + 1L]] <- idx; nms <- c(nms, paste(fn[a], fn[b], sep = " & ")) }
  }
  names(cands) <- nms
  cands
}


#' Build the influence matrix for a candidate family
#'
#' For each candidate subgroup, re-derives the treatment-coefficient `dfbeta`,
#' point estimate, and split-perturbation SD through the engine's own
#' `.consistency_cox_pieces()` / `.consistency_glm_pieces()`, so the influence
#' contributions match the real consistency calculation.  Candidates whose model
#' fails to converge (or that yield a non-finite SD) are dropped.
#'
#' @param data Analysis data frame.
#' @param candidates Named list of integer row-index vectors (see
#'   [enumerate_candidates()]).
#' @param spec List describing the effect model: `outcome_type`
#'   (`"survival"`/`"binary"`/`"continuous"`/`"count"`), `effect_measure`,
#'   `treat.name`, `outcome.name`, `offset.name`, `adjust_covariates`,
#'   `adverse_outcome`, and for survival `tte.name`/`event.name`.
#' @return List with `B` (N x S' influence matrix, usable candidates only),
#'   `beta_hat`, `sigma_D`, `names`, `n_evaluated`, `n_usable`.
#' @export
build_influence_matrix <- function(data, candidates, spec, N = nrow(data)) {
  fns    <- .resolve_pieces_fns()
  is_cox <- identical(spec$outcome_type, "survival")
  S <- length(candidates)
  B <- matrix(0, nrow = N, ncol = S)
  beta_hat <- rep(NA_real_, S); sigma_D <- rep(NA_real_, S); ok <- logical(S)

  glm_args <- list(outcome_type = spec$outcome_type, effect_measure = spec$effect_measure,
                   treat.name = spec$treat.name, outcome.name = spec$outcome.name,
                   offset.name = spec$offset.name, adjust_covariates = spec$adjust_covariates,
                   adverse_outcome = isTRUE(spec$adverse_outcome))
  cox_args <- list(tte.name = spec$tte.name, event.name = spec$event.name,
                   treat.name = spec$treat.name, adjust_covariates = spec$adjust_covariates)

  for (g in seq_len(S)) {
    idx <- candidates[[g]]
    dfg <- data[idx, , drop = FALSE]
    pieces <- tryCatch(
      .call_pieces(if (is_cox) fns$cox else fns$glm, dfg, if (is_cox) cox_args else glm_args),
      error = function(e) NULL)
    if (is.null(pieces) || !is.finite(pieces$sigma_D) || pieces$sigma_D <= 0) next
    if (length(pieces$dfbeta) != length(idx)) next
    B[idx, g] <- pieces$dfbeta
    beta_hat[g] <- pieces$beta_hat; sigma_D[g] <- pieces$sigma_D; ok[g] <- TRUE
  }

  list(B = B[, ok, drop = FALSE], beta_hat = beta_hat[ok], sigma_D = sigma_D[ok],
       names = if (!is.null(names(candidates))) names(candidates)[ok] else NULL,
       n_evaluated = S, n_usable = sum(ok))
}


#' Sweep screening / consistency thresholds and report the family FDR
#'
#' @param infl Output of [build_influence_matrix()].
#' @param beta_null Null configuration (scalar log effect or length-S' vector),
#'   e.g. the log ITT effect for a uniform-benefit null.
#' @param grid Data frame with columns `c_screen` and `p_star` (and optionally
#'   `c_consistency`; default `log(1.0)`).
#' @param multiplier,draws,seed Passed to [fdr_family_multiplier()].
#' @return `grid` augmented with `fdr`, `union_indep`, `e_flagged`, `obs_flagged`.
#' @export
calibrate_fdr_grid <- function(infl, beta_null, grid,
                               multiplier = "rademacher", draws = 4000L, seed = 1L) {
  if (is.null(grid$c_consistency)) grid$c_consistency <- log(1.0)
  out <- grid
  out$fdr <- NA_real_; out$union_indep <- NA_real_
  out$e_flagged <- NA_real_; out$obs_flagged <- NA_integer_
  for (i in seq_len(nrow(grid))) {
    r <- fdr_family_multiplier(
      influence = infl$B, beta_hat = infl$beta_hat, sigma_D = infl$sigma_D,
      c_screen = grid$c_screen[i], c_consistency = grid$c_consistency[i],
      p_star = grid$p_star[i], beta_null = beta_null, center = "null",
      draws = draws, multiplier = multiplier, seed = seed)
    out$fdr[i] <- r$fdr; out$union_indep[i] <- r$union_independent
    out$e_flagged[i] <- r$e_num_flagged; out$obs_flagged[i] <- r$n_observed_flagged
  }
  out
}


# =============================================================================
# DEV-ONLY DEMO -- ACTG175-style binary/OR null (homogeneous protective effect)
# =============================================================================

#' Self-contained calibration demonstration on a simulated null
#'
#' Mimics the ACTG175 binary/OR setting: a uniformly protective treatment
#' (OR ~ 0.66, no harm subgroup), six continuous covariates cut at q1/median/q3
#' plus two binary factors, giving a few hundred overlapping candidates. Builds
#' the influence matrix and sweeps the screening OR and consistency level,
#' printing the family false-discovery rate for each setting.
#'
#' @keywords internal
run_fdr_calibration_demo <- function(N = 1000L, draws = 3000L, seed = 7L) {
  set.seed(seed)
  trt <- rbinom(N, 1, 0.5)
  X <- as.data.frame(matrix(rnorm(N * 6), N, 6)); names(X) <- paste0("x", 1:6)
  b1 <- rbinom(N, 1, 0.5); b2 <- rbinom(N, 1, 0.4)
  lp <- -0.30 - 0.42 * trt + 0.30 * X$x1 - 0.20 * X$x2 + 0.25 * b1   # OR(trt) = exp(-0.42) ~ 0.66
  y  <- rbinom(N, 1, plogis(lp))
  dat <- data.frame(y, trt, X, b1, b2)

  factors <- list()
  for (v in paste0("x", 1:6)) {
    q <- stats::quantile(dat[[v]], c(.25, .5, .75))
    factors[[paste0(v, "<=q1")]]  <- dat[[v]] <= q[1]
    factors[[paste0(v, "<=med")]] <- dat[[v]] <= q[2]
    factors[[paste0(v, "<=q3")]]  <- dat[[v]] <= q[3]
  }
  factors[["b1=1"]] <- dat$b1 == 1; factors[["b1=0"]] <- dat$b1 == 0
  factors[["b2=1"]] <- dat$b2 == 1; factors[["b2=0"]] <- dat$b2 == 0

  cands <- enumerate_candidates(dat, factors, maxk = 2L, min_n = 60L,
                                min_events_per_arm = 5L, treat.name = "trt", outcome.name = "y")
  spec <- list(outcome_type = "binary", effect_measure = "OR", treat.name = "trt",
               outcome.name = "y", offset.name = NULL, adjust_covariates = NULL,
               adverse_outcome = TRUE)
  infl <- build_influence_matrix(dat, cands, spec)
  itt  <- stats::coef(glm(y ~ trt, family = binomial(), data = dat))[["trt"]]

  grid <- expand.grid(c_screen = log(c(1.25, 1.50, 1.75)),
                      p_star   = c(0.90, 0.95, 0.975))
  tab <- calibrate_fdr_grid(infl, beta_null = itt, grid = grid, draws = draws, seed = seed)
  tab <- tab[order(tab$c_screen, tab$p_star), ]
  tab$OR_screen <- round(exp(tab$c_screen), 2)

  cat("---- Threshold calibration via family multiplier bootstrap : simulated NULL ----\n")
  cat(sprintf("Subjects N                : %d\n", N))
  cat(sprintf("Candidate family          : %d enumerated, %d usable\n",
              infl$n_evaluated, infl$n_usable))
  cat(sprintf("ITT log-OR (null centre)  : %+.3f  (OR %.2f)\n", itt, exp(itt)))
  cat(sprintf("Multiplier / draws        : rademacher / %d\n", draws))
  cat("--------------------------------------------------------------------------------\n")
  print(tab[, c("OR_screen", "p_star", "fdr", "union_indep", "e_flagged")], row.names = FALSE)
  cat("--------------------------------------------------------------------------------\n")
  cat("Read down a column: tightening screening OR and/or the consistency level p_star\n")
  cat("lowers the search-level family FDR. Pick the setting meeting your target rate.\n")
  invisible(list(infl = infl, table = tab, itt = itt))
}

if (identical(environment(), globalenv()) && !interactive() &&
    sys.nframe() == 0L && length(commandArgs(trailingOnly = TRUE)) == 0L) {
  run_fdr_calibration_demo()
}
