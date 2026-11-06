#' Naive, bootstrap and multiplier estimates, side by side, and their cost
#'
#' Plain data frames for the comparison of a [forestsearch()] analysis's naive
#' estimates with its full-bootstrap and multiplier-resampling corrections, method by
#' method; the caller styles them.
#'
#' **What the comparison shows.** The table compares displays; it does not adjudicate
#' one quantity. The bootstrap column carries the point estimate only: at a small `B`
#' that point is an illustration, and its interval is not a result. For families
#' generated from a model (the DINA and GRF paths), the full bootstrap, which re-runs the
#' search in every replicate, and multiplier resampling, which perturbs one fixed
#' family, target different estimands.
#'
#' `mr_cost_table()`'s `mr_seconds` is the multiplier stage as the fit records it
#' (`mr_inference$timing_seconds`), and the projection is linear in `B` at the worker
#' count the bootstrap ran on.
#'
#' @param fits A named list of [forestsearch()] fits with multiplier resampling and the
#'   complement de-biased. The names label the methods.
#' @param boots A list of [forestsearch_bootstrap_dofuture()] results with the same
#'   names.
#' @param est.scale `"hr"` (the default), or `"1/hr"`, which inverts every ratio (limits
#'   swap) and is accepted for ratio measures only.
#' @param digits Decimals of the estimates. Default 2.
#' @param boot_seconds The wall-clock seconds of each bootstrap, a vector named by
#'   method. `NULL` (the default) reads what each result stores, the time of its
#'   replicate loop (`attr(boot$results, "timing")$total_minutes`), which leaves out the
#'   call's setup and summaries.
#' @param project_B The number of replicates to project the cost to. Default 1000.
#'
#' @return `mr_method_comparison()`: a plain data frame with two rows per method, `H` and
#'   `Hc`: `method`, `region`, `n` (patients, with their percent), `naive` (estimate with
#'   its robust 95% interval), `bootstrap` (the point estimate), and `multiplier` (the
#'   two-term de-biased estimate with its IJ 95% interval). The interval reads `NA` on a
#'   fit whose stored limits are not the two-term IJ interval; see [mr_values()].
#'
#'   `mr_cost_table()`: a plain data frame with one row per method: `method`,
#'   `mr_seconds`, `bootstrap_B`, `bootstrap_seconds`, `seconds_per_replicate`, and
#'   `projected_B_<project_B>`.
#'
#' @seealso [mr_values()], [mr_estimates_table()], [forestsearch_bootstrap_dofuture()].
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' n <- 500
#' d <- data.frame(id = seq_len(n), treat = rbinom(n, 1, 0.5), x1 = rnorm(n), x2 = rnorm(n))
#' d$time <- rexp(n, 0.1 * ifelse(d$treat == 1, ifelse(d$x1 < -0.6, 3, 0.45), 1))
#' d$event <- as.integer(d$time < 30)
#' d$time <- pmin(d$time, 30)
#' fit <- forestsearch(
#'   df.analysis = d, outcome.name = "time", event.name = "event", treat.name = "treat",
#'   id.name = "id", confounders.name = c("x1", "x2"), use_grf = FALSE, use_lasso = FALSE,
#'   maxk = 1, fs.splits = 200, sg_focus = "hrMaxSG", mr_inference = TRUE,
#'   mr_inference_args = list(draws = 200, field_R_out = 50, field_R_in = 25),
#'   parallel_args = list(plan = "sequential"), quiet = TRUE, details = FALSE,
#'   plot.sg = FALSE, show_candidate_summary = FALSE)
#' boot <- forestsearch_bootstrap_dofuture(fit, nb_boots = 2,
#'                                         parallel_args = list(plan = "sequential"))
#' mr_method_comparison(list(FS = fit), list(FS = boot))
#' mr_cost_table(list(FS = fit), list(FS = boot))
#' }
#'
#' @name mr_method_tables
NULL


#' @rdname mr_method_tables
#' @export
mr_method_comparison <- function(fits, boots, est.scale = "hr", digits = 2) {
  .mr_method_check(fits, boots)
  est.scale <- match.arg(est.scale, c("hr", "1/hr"))
  inv <- est.scale == "1/hr"
  tr <- function(x) if (inv) 1 / x else x
  f1 <- function(x) sprintf("%.*f", digits, tr(x))
  ci <- function(e, l, u) {
    sprintf("%s (%s, %s)", f1(e), if (inv) f1(u) else f1(l), if (inv) f1(l) else f1(u))
  }
  result <- do.call(rbind, lapply(names(fits), function(method) {
    fit <- fits[[method]]; boot <- boots[[method]]
    if (inv && !isTRUE(fit$mr_inference$log_scale))
      stop("est.scale = \"1/hr\" applies to ratio measures only; method '", method,
           "' is on an identity scale.", call. = FALSE)
    v <- mr_values(fit)
    H <- v[v$row == "H", ]; Hc <- v[v$row == "Hc", ]
    if (!nrow(Hc) || is.na(Hc$n))
      stop("Method '", method, "': the fit has no de-biased complement; fit with ",
           "mr_inference_args = list(include_complement = TRUE).", call. = FALSE)
    point <- c(.mr_method_point(boot, "H_estimates", method), .mr_method_point(boot, "Hc_estimates", method))
    # H2 is stored on the fit's own est.scale; mr_inference is always on the HR scale
    stored <- fit$args_call_all$est.scale
    if (identical(stored, "1/hr")) point <- 1 / point
    pct <- function(k) sprintf("%d (%d%%)", k, round(100 * k / (H$n + Hc$n)))
    data.frame(
      method = method, region = c("H", "Hc"),
      n = c(pct(H$n), pct(Hc$n)),
      naive = c(ci(H$naive_est, H$naive_lower, H$naive_upper),
                ci(Hc$naive_est, Hc$naive_lower, Hc$naive_upper)),
      bootstrap = f1(point),
      multiplier = c(ci(H$debiased_est, H$ij_lower, H$ij_upper),
                     ci(Hc$debiased_est, Hc$ij_lower, Hc$ij_upper)),
      stringsAsFactors = FALSE)
  }))
  rownames(result) <- NULL
  result
}


#' @rdname mr_method_tables
#' @export
mr_cost_table <- function(fits, boots, boot_seconds = NULL, project_B = 1000L) {
  .mr_method_check(fits, boots)
  if (!is.numeric(project_B) || length(project_B) != 1L || !is.finite(project_B) || project_B <= 0)
    stop("project_B must be one positive number.", call. = FALSE)
  methods <- names(fits)
  if (is.null(boot_seconds)) {
    boot_seconds <- vapply(methods, function(m) {
      tm <- attr(boots[[m]]$results, "timing")$total_minutes
      if (is.null(tm))
        stop("Method '", m, "': the bootstrap result stores no wall time; give boot_seconds.",
             call. = FALSE)
      60 * tm
    }, 0)
  } else if (is.null(names(boot_seconds)) || !all(methods %in% names(boot_seconds))) {
    stop("boot_seconds must be a vector named by method.", call. = FALSE)
  }
  cost <- data.frame(
    method = methods,
    mr_seconds = vapply(methods, function(m) as.numeric(fits[[m]]$mr_inference$timing_seconds), 0),
    bootstrap_B = vapply(methods, function(m) as.numeric(boots[[m]]$nb_boots), 0),
    bootstrap_seconds = as.numeric(boot_seconds[methods]),
    stringsAsFactors = FALSE)
  cost$seconds_per_replicate <- cost$bootstrap_seconds / cost$bootstrap_B
  cost[[paste0("projected_B_", project_B)]] <- round(cost$seconds_per_replicate * project_B)
  rownames(cost) <- NULL
  cost
}


#' The two lists: named alike, and each fit with its multiplier resampling
#' @noRd
.mr_method_check <- function(fits, boots) {
  if (!is.list(fits) || is.null(names(fits)) || !length(fits) || any(!nzchar(names(fits))) ||
      anyDuplicated(names(fits)))
    stop("fits must be a named list of forestsearch() fits; the names label the methods.",
         call. = FALSE)
  if (!is.list(boots) || !all(names(fits) %in% names(boots)))
    stop("boots must be a list with the names of fits: ", paste(names(fits), collapse = ", "),
         ".", call. = FALSE)
  for (m in names(fits))
    if (is.null(fits[[m]][["mr_inference"]]))
      stop("Method '", m, "': the fit has no mr_inference results; run ",
           "forestsearch(mr_inference = TRUE).", call. = FALSE)
  invisible(TRUE)
}


#' The bootstrap's corrected point estimate for one region
#' @noRd
.mr_method_point <- function(boot, element, method) {
  v <- boot[[element]][["H2"]]
  if (length(v) != 1L)
    stop("Method '", method, "': the bootstrap result has no ", element, "$H2.", call. = FALSE)
  as.numeric(v)
}
