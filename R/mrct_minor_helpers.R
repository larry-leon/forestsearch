#' Read a shipped record of the MRCT consistency campaign
#'
#' The records behind the MRCT consistency articles ship with the package, in
#' `extdata/mrct-consistency`, with their provenance in its `PROVENANCE.md`. From a source
#' tree the local copy is used; otherwise the installed one.
#'
#' @param file The file's name, for instance `"rep5_trial.csv"`.
#' @param ... Passed to [utils::read.csv()].
#'
#' @return A data frame, read with `stringsAsFactors = FALSE` and `check.names = FALSE`.
#'
#' @seealso [mrct_case_study()], [mrct_trial_analysis()].
#'
#' @examples
#' trial <- mrct_record("rep5_trial.csv")
#' dim(trial)
#' mrct_record("mr_field_S2_rep5.csv")$sg_found
#'
#' @export
mrct_record <- function(file, ...) {
  if (!is.character(file) || length(file) != 1L || !nzchar(file))
    stop("file must be one file name.", call. = FALSE)
  dd <- .mrct_record_dir()
  path <- file.path(dd, file)
  if (!file.exists(path))
    stop("There is no record '", file, "' in ", dd, ".", call. = FALSE)
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE, ...)
}


#' The directory of the shipped records: the source tree's, else the installed one
#' @noRd
.mrct_record_dir <- function(candidates = c(file.path("..", "..", "inst", "extdata", "mrct-consistency"),
                                            file.path("inst", "extdata", "mrct-consistency"),
                                            system.file("extdata", "mrct-consistency", package = "forestsearch"))) {
  dd <- candidates[nzchar(candidates) & dir.exists(candidates)][1]
  if (is.na(dd))
    stop("The MRCT records directory (extdata/mrct-consistency) was not found; ",
         "is forestsearch installed in full?", call. = FALSE)
  if (!file.exists(file.path(dd, "PROVENANCE.md")))
    stop("The records directory ", dd, " has no PROVENANCE.md, so its records cannot be ",
         "trusted as the shipped ones.", call. = FALSE)
  dd
}


#' Membership of a forestsearch fit's declared subgroup
#'
#' Applies the fit's declared rule to the rows of `data`, by [get_dfpred()].
#'
#' @param fit A [forestsearch()] fit.
#' @param data A data frame with the columns the rule names, such as the fit's whole
#'   trial.
#'
#' @return A logical vector, one value per row: `TRUE` in the declared subgroup
#'   \eqn{\hat H}. All `FALSE` when nothing was declared.
#'
#' @seealso [get_dfpred()], [forestsearch()].
#'
#' @examples
#' set.seed(1)
#' n <- 500
#' d <- data.frame(id = seq_len(n), treat = rbinom(n, 1, 0.5), x1 = rnorm(n), x2 = rnorm(n))
#' d$time <- rexp(n, 0.1 * ifelse(d$treat == 1, ifelse(d$x1 < -0.6, 3, 0.45), 1))
#' d$event <- as.integer(d$time < 30)
#' d$time <- pmin(d$time, 30)
#' fit <- forestsearch(
#'   df.analysis = d, outcome.name = "time", event.name = "event", treat.name = "treat",
#'   id.name = "id", confounders.name = c("x1", "x2"), use_grf = FALSE, use_lasso = FALSE,
#'   maxk = 1, fs.splits = 200, sg_focus = "hrMaxSG",
#'   parallel_args = list(plan = "sequential"), quiet = TRUE, details = FALSE,
#'   plot.sg = FALSE, show_candidate_summary = FALSE)
#' fit$sg.harm
#' table(fs_in_subgroup(fit, d))
#'
#' @export
fs_in_subgroup <- function(fit, data) {
  if (!inherits(fit, "forestsearch")) stop("fit must be a forestsearch() fit.", call. = FALSE)
  if (!is.data.frame(data)) stop("data must be a data frame.", call. = FALSE)
  if (is.null(fit[["sg.harm"]])) return(rep(FALSE, nrow(data)))
  get_dfpred(as.data.frame(data), sg.harm = fit[["sg.harm"]], version = 2)$treat.recommend == 0
}


#' Build an MRCT data-generating mechanism, keeping its warnings
#'
#' A thin wrapper: [inject_mrct_structure()], then [scale_dgm_hazard()] when `m != 1`.
#' The build's warnings (usually the censoring model's fit) are muffled and kept in the
#' `"build_warnings"` attribute.
#'
#' @param data The seed data: `seed_data` of [inject_mrct_structure()].
#' @param ... Passed to [inject_mrct_structure()].
#' @param m The hazard multiplier of [scale_dgm_hazard()]. Default 1, no scaling.
#'
#' @return The DGM, with the attribute `"build_warnings"`: a character vector,
#'   `character(0)` for a clean build.
#'
#' @seealso [inject_mrct_structure()], [scale_dgm_hazard()].
#'
#' @examples
#' set.seed(1)
#' n <- 400
#' d <- data.frame(age = rnorm(n, 50, 10), lab = rnorm(n, 300, 100),
#'                 f1 = rbinom(n, 1, 0.4), treat = rbinom(n, 1, 0.5))
#' d$tte <- rexp(n, 0.05 * exp(-0.5 * d$treat))
#' d$event <- as.integer(d$tte < 40)
#' d$tte <- pmin(d$tte, 40)
#' dgm <- mrct_build_dgm(
#'   d, outcome_var = "tte", event_var = "event", treatment_var = "treat",
#'   continuous_vars = c("age", "lab"), factor_vars = "f1", x_pred = "lab",
#'   spline_spec = list(knot = 300, zeta = 800, log_hrs = rep(log(0.55), 3)),
#'   region = list(prevalence = 0.25, or_pred = 1, x1_vars = "f1", or_x1 = 3,
#'                 loghr = -log(3), band = list(cut = 230, or = 1)),
#'   n_super = 2000, expand = "bootstrap", seed = 7L,
#'   continuous_vars_cens = c("age", "lab"), factor_vars_cens = "f1",
#'   m = 3)
#' length(attr(dgm, "build_warnings"))
#'
#' @export
mrct_build_dgm <- function(data, ..., m = 1) {
  if (!is.numeric(m) || length(m) != 1L || !is.finite(m) || m <= 0)
    stop("m must be one positive number.", call. = FALSE)
  w <- character(0)
  dgm <- withCallingHandlers(
    inject_mrct_structure(data, ...),
    warning = function(cnd) { w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleWarning") })
  if (m != 1) dgm <- scale_dgm_hazard(dgm, m = m)
  attr(dgm, "build_warnings") <- w
  dgm
}
