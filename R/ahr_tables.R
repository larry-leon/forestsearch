#' Tables of the MRCT analytic layer, over cells and sets
#'
#' Tidy data frames from the existing analytic functions, without quadrature of their
#' own: `ahr_truths_table()` loops over [ahr_limit()] and `ahr_design_table()` over
#' [ahr_design()].
#'
#' `ahr_design_table()`'s `r1_` to `r7_` columns are the seven design probabilities of
#' a regional-consistency analysis. With \eqn{A} = \{ITT ub < 1.0\}, \eqn{B} = \{region
#' ITT > 1.0\} and \eqn{C^*} = \{region(\eqn{G^*}) < 1.0\}: (1) \eqn{A}, (2) \eqn{B},
#' (3) consistency \eqn{A \setminus B}, (4) non-consistency (the trigger) \eqn{A \cap
#' B}, (5) oracle consistency \eqn{A \cap B \cap C^*}, (6) oracle non-consistency, and
#' (7) consistency or oracle consistency. The function stops unless (3) + (4) = (1),
#' (5) + (6) = (4) and (3) + (5) = (7) at 1e-12.
#'
#' @param x The cell or cells: an [ahr_population()] result, a DGM from
#'   [inject_mrct_structure()], or a named list of either. The list's names label the
#'   cells.
#' @param sets The sets, as [ahr_limit()] takes them: a character vector of set names
#'   (`"all"`, `"nonregion"`, `"region"`, `"region_Gstar"`, `"Hstar"`, `"Hstar_region"`,
#'   `"Hstar_nonregion"`), or a list of those and logical vectors of length `pop$N`. A
#'   name on an element labels the set, and a logical vector needs one.
#' @param b,h Hazard ratios outside and inside \eqn{H^*}: one number, or a vector named
#'   by cell. In `ahr_truths_table()`, both `NULL` (the default) use the DGM's own
#'   per-patient log hazard ratios.
#' @param analysis_time,max_entry Passed to [ahr_population()] when `x` holds DGMs, and
#'   not used for an [ahr_population()] result, which carries its own.
#' @param n The trial size.
#' @param ... Passed to [ahr_limit()] (for instance `horizon`, `b_R`, `h_R`) or
#'   [ahr_design()] (`b_R`, `h_R`).
#'
#' @return `ahr_truths_table()`: a plain data frame with one row per cell and set:
#'   `cell`, `set`, `horizon` (the analysis time), `marghr` (the marginal hazard ratio at
#'   the horizon), `marghr_inf` (with no censoring), `ahr`, `cde`.
#'
#'   `ahr_design_table()`: a plain data frame with one row per cell: `cell`, `n`, `b`,
#'   `h`; the five probabilities of [ahr_design()] under their names there
#'   (`p_itt_ub_lt1`, `p_region_itt_gt1`, `p_region_Gstar_lt1`, `p_trigger`,
#'   `p_trigger_and_Gstar_lt1`); `p_Gstar_lt1_given_trigger`, their ratio; and the seven
#'   rows above, `r1_itt_ub_lt1`, `r2_region_itt_gt1`, `r3_consistency`,
#'   `r4_non_consistency`, `r5_oracle_consistency`, `r6_oracle_non_consistency`,
#'   `r7_consistency_or_oracle`.
#'
#' @seealso [ahr_population()], [ahr_limit()], [ahr_design()].
#'
#' @examples
#' set.seed(1)
#' n <- 400
#' d <- data.frame(age = rnorm(n, 50, 10), lab = rnorm(n, 300, 100),
#'                 f1 = rbinom(n, 1, 0.4), treat = rbinom(n, 1, 0.5))
#' d$tte <- rexp(n, 0.05 * exp(-0.5 * d$treat))
#' d$event <- as.integer(d$tte < 40)
#' d$tte <- pmin(d$tte, 40)
#' dgm <- suppressWarnings(inject_mrct_structure(
#'   d, outcome_var = "tte", event_var = "event", treatment_var = "treat",
#'   continuous_vars = c("age", "lab"), factor_vars = "f1", x_pred = "lab",
#'   spline_spec = list(knot = 300, zeta = 800, log_hrs = rep(log(0.55), 3)),
#'   region = list(prevalence = 0.25, or_pred = 1, x1_vars = "f1", or_x1 = 3,
#'                 loghr = -log(3), band = list(cut = 230, or = 1)),
#'   x3 = list(vars = "lab", cuts = list(lab = 230), loghr = log(3 / 0.55)),
#'   n_super = 2000, expand = "bootstrap", seed = 7L,
#'   continuous_vars_cens = c("age", "lab"), factor_vars_cens = "f1"))
#' pop <- ahr_population(dgm, analysis_time = 30, max_entry = 6)
#' ahr_truths_table(pop, sets = c("all", "region", "region_Gstar"), b = 0.55, h = 3)
#' ahr_design_table(pop, b = 0.55, h = 3, n = 800)
#'
#' @name ahr_tables
NULL


#' @rdname ahr_tables
#' @export
ahr_truths_table <- function(x, sets = c("all", "nonregion", "region", "region_Gstar"),
                             b = NULL, h = NULL, analysis_time = NULL, max_entry = NULL, ...) {
  pops <- .ahr_tables_pops(x, analysis_time, max_entry)
  sets_named <- .ahr_tables_sets(sets)
  result <- do.call(rbind, lapply(names(pops), function(k) {
    do.call(rbind, lapply(names(sets_named), function(s) {
      e <- ahr_limit(pops[[k]], sets_named[[s]], b = .ahr_tables_value(b, k, "b"),
                     h = .ahr_tables_value(h, k, "h"), variance = FALSE, ...)
      data.frame(cell = k, set = s, horizon = e$horizon, marghr = e$hr, marghr_inf = e$marghr_inf,
                 ahr = e$ahr, cde = e$cde, stringsAsFactors = FALSE)
    }))
  }))
  rownames(result) <- NULL
  result
}


#' @rdname ahr_tables
#' @export
ahr_design_table <- function(x, b, h, n, analysis_time = NULL, max_entry = NULL, ...) {
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) || n <= 0)
    stop("n must be one positive number, the trial size.", call. = FALSE)
  if (missing(b) || missing(h) || is.null(b) || is.null(h))
    stop("b and h are required: ahr_design() needs the hazard ratios.", call. = FALSE)
  pops <- .ahr_tables_pops(x, analysis_time, max_entry)
  near <- function(a, b) abs(a - b) <= 1e-12 * max(1, abs(b))
  result <- do.call(rbind, lapply(names(pops), function(k) {
    b_k <- .ahr_tables_value(b, k, "b"); h_k <- .ahr_tables_value(h, k, "h")
    p <- ahr_design(pops[[k]], b = b_k, h = h_k, n = n, ...)$prob
    pA <- p[["p_itt_ub_lt1"]]; pB <- p[["p_region_itt_gt1"]]
    pAB <- p[["p_trigger"]]; pABC <- p[["p_trigger_and_Gstar_lt1"]]
    rows <- c(pA, pB, pA - pAB, pAB, pABC, pAB - pABC, (pA - pAB) + pABC)
    stopifnot(near(rows[3] + rows[4], rows[1]), near(rows[5] + rows[6], rows[4]),
              near(rows[3] + rows[5], rows[7]))
    data.frame(cell = k, n = n, b = b_k, h = h_k,
               p_itt_ub_lt1 = pA, p_region_itt_gt1 = pB,
               p_region_Gstar_lt1 = p[["p_region_Gstar_lt1"]], p_trigger = pAB,
               p_trigger_and_Gstar_lt1 = pABC, p_Gstar_lt1_given_trigger = pABC / pAB,
               r1_itt_ub_lt1 = rows[1], r2_region_itt_gt1 = rows[2], r3_consistency = rows[3],
               r4_non_consistency = rows[4], r5_oracle_consistency = rows[5],
               r6_oracle_non_consistency = rows[6], r7_consistency_or_oracle = rows[7],
               stringsAsFactors = FALSE)
  }))
  rownames(result) <- NULL
  result
}


#' The cells as a named list of ahr_population objects
#' @noRd
.ahr_tables_pops <- function(x, analysis_time, max_entry) {
  is_cell <- function(v) inherits(v, c("ahr_population", "aft_dgm_flex", "aft_dgm"))
  cells <- if (is_cell(x)) list(cell = x) else x
  if (!is.list(cells) || !length(cells) || !all(vapply(cells, is_cell, NA)))
    stop("x must be an ahr_population() result, a DGM, or a list of either.", call. = FALSE)
  if (is.null(names(cells)) || any(!nzchar(names(cells))) || anyDuplicated(names(cells)))
    stop("A list x needs unique names: they label the cells.", call. = FALSE)
  lapply(cells, function(v) {
    if (inherits(v, "ahr_population")) return(v)
    if (is.null(analysis_time) || is.null(max_entry))
      stop("analysis_time and max_entry are required when x holds DGMs.", call. = FALSE)
    ahr_population(v, analysis_time = analysis_time, max_entry = max_entry)
  })
}


#' The sets as a named list, as ahr_limit() takes them
#' @noRd
.ahr_tables_sets <- function(sets) {
  sets <- as.list(sets)
  if (!length(sets)) stop("sets must name at least one set.", call. = FALSE)
  nm <- names(sets); if (is.null(nm)) nm <- rep("", length(sets))
  for (i in seq_along(sets)) {
    if (nzchar(nm[i])) next
    if (!is.character(sets[[i]]))
      stop("A logical vector in sets needs a name: it labels the set.", call. = FALSE)
    nm[i] <- sets[[i]]
  }
  if (anyDuplicated(nm)) stop("The sets must have unique labels.", call. = FALSE)
  stats::setNames(sets, nm)
}


#' A hazard ratio for one cell: one number, or a vector named by cell
#' @noRd
.ahr_tables_value <- function(v, cell, arg) {
  if (is.null(v) || length(v) == 1L) return(if (is.null(v)) NULL else unname(v))
  if (is.null(names(v)) || !cell %in% names(v))
    stop(arg, " must be one number or a vector named by cell; cell '", cell, "' is not in it.",
         call. = FALSE)
  unname(v[[cell]])
}
