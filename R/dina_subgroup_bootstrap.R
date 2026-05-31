# File: R/dina_subgroup_bootstrap.R
# Part of the forestsearch package.
#
# Bootstrap inference for the subgroup discovered by dina_subgroup().
# Reports two conditional-on-signature effect CIs for the original-data
# subgroup S*:
#   (1) DINA effect_ci -- percentile CI on a*'beta_b with a* fixed from
#       S* (the BLP-analog); and
#   (2) within-subgroup standard-model CI -- a plain Cox/GLM treatment
#       contrast refit within the fixed signature on each resample
#       (see dina_subgroup_refit()).
# Neither CI adjusts for the data-driven selection of the signature;
# selection stability is summarized separately by selection_frequency.
# Full design rationale is in the exported function's roxygen below.


# ---------------------------------------------------------------------------
# Main exported function
# ---------------------------------------------------------------------------

#' Selection-adjusted bootstrap CI for dina_subgroup()
#'
#' Performs selection-adjusted inference for the harm subgroup discovered
#' by `dina_subgroup()`.  For each of `n_boot` bootstrap resamples of the
#' input data frame, the function refits DINA via `dina()` and re-runs
#' the subgroup search via `dina_subgroup()`, then reports three classes
#' of output:
#'
#'   * An **effect CI** on `a*^T beta_b`, where `a* = (1, x_bar_S)` is
#'     held FIXED from the original-data subgroup `S*` and `beta_b` is
#'     each iteration's DINA coefficient vector.  This is the
#'     bootstrap analog of the Wald CI returned by
#'     `dina_subgroup()` -- conditional on the discovered subgroup
#'     definition, with width reflecting sampling variance of
#'     `beta-hat`.
#'   * **Stability CIs** on the bootstrap-selected `n_subgroup` and the
#'     cut threshold(s), restricted to iterations that selected the same
#'     cut-set (covariate(s) and direction(s)) as the original-data
#'     subgroup.  For a depth-2 conjunction the threshold CI is reported
#'     per covariate.
#'   * A **selection-frequency table** counting how often each cut-set
#'     was chosen.  The structural-selection stability diagnostic.
#'
#' Each bootstrap iteration calls `dina()` (single-pass, sandwich
#' variance) under the hood, not `dina_bagged()`.
#'
#' If you have already run `dina()` and `dina_subgroup()` on the original
#' data, pass those objects via `fit` and/or `sg` so the original-data
#' point estimate reuses them rather than being recomputed here.  This
#' guarantees the reported point estimate matches your standalone result
#' exactly; otherwise the internal refit uses its own cross-fitting fold
#' assignment, which at a boundary-case `m_diff` can disagree with the
#' upstream fit on whether a subgroup qualifies.
#'
#' @section Why the bootstrap target is the *fixed-subgroup* effect:
#' A naive bootstrap of `dina_subgroup()`'s `mean_tau_hat` is
#' uninformative.  Under the default `sg_focus = "maxSG"` the search
#' selects the *largest* subgroup whose mean tau-hat exceeds `m_diff`,
#' so the chosen subgroup lies at the boundary by construction: adding
#' one more patient would push the mean below `m_diff`.  Each bootstrap
#' iteration's selected subgroup is similarly pinned, so the bootstrap
#' distribution of `mean_tau_hat` clusters tightly above `m_diff`
#' and the resulting CI quantifies threshold-grid granularity
#' rather than effect uncertainty.  (Under the band foci
#' `"effMaxSG"` / `"effMinSG"` the per-iteration subgroup is instead
#' anchored to the maximum effect rather than the `m_diff` boundary,
#' but the same caveat applies: `mean_tau_hat` is a selection-pinned
#' quantity, not an effect estimate.)  The fixed-subgroup linear
#' functional `a*^T beta_b` instead targets the clinically
#' meaningful estimand "treatment effect for patients whose
#' covariates lie in the discovered subgroup".  Because DINA's
#' `beta-hat` is cross-fit (and hence asymptotically unbiased even
#' for data-driven contrasts), the percentile distribution of
#' `a*^T beta_b` gives a non-parametric CI that complements the
#' sandwich Wald CI of `dina_subgroup()`.
#'
#' @section Two complementary effect estimates:
#' When `refit = TRUE` (default), the function reports two treatment
#' effects for the discovered subgroup, on the same natural-parameter
#' scale:
#' \itemize{
#'   \item the \strong{DINA effect} (`effect_ci`): the BLP-analog
#'     `a*^T beta`, a projection of the cross-fitted CATE surface
#'     averaged over the subgroup; and
#'   \item the \strong{within-subgroup standard-model effect}
#'     (`refit_effect_ci`): the treatment coefficient from a plain Cox
#'     model (or GLM) fit on the subgroup, refit within the FIXED
#'     discovered signature on each resample.  This is the contrast a
#'     clinical-trial reader usually expects.
#' }
#' The two coincide only under correct linearity; reporting both is
#' informative.  Crucially, \emph{both} CIs are conditional on the
#' discovered signature treated as pre-specified -- neither corrects for
#' the data-driven selection of the signature.  A selection-adjusted
#' (de-biased) interval requires the more involved bias-corrected
#' bootstrap of Leon et al. (2024, Section 3.2) and is not provided
#' here; selection \emph{stability} is instead summarized by
#' `selection_frequency`.
#'
#' @section Performance:
#' With `parallel = "boots"`, each bootstrap iteration runs in a worker
#' via `foreach::foreach() %dofuture%`.  Set `future::plan()` to a
#' parallel backend (e.g., `future::multisession`) before calling, and
#' reset to `"sequential"` afterwards.  Workers run from the installed
#' package; `devtools::install()` is required after code edits before
#' the parallel branch will see them.
#'
#' @param df data frame containing the outcome, treatment, and
#'   covariate columns.
#' @param outcome character(1); name of the response column in `df`.
#'   For `family = "cox"` this is the event/censoring time column.
#' @param treatment character(1); name of the binary treatment
#'   indicator column.
#' @param covariates character vector; names of the covariate columns
#'   to search over.  All must be numeric.
#' @param family one of `"gaussian"`, `"binomial"`, `"poisson"`, `"cox"`.
#' @param status character(1) or `NULL`; for `family = "cox"` only,
#'   the event indicator column name.
#' @param m_diff scalar harm threshold on the natural-parameter scale.
#' @param n_min positive integer; minimum subgroup size.  Default `60L`.
#' @param direction one of `"both"` (default), `"left"`, `"right"`.
#' @param max_depth integer, `1L` (default) or `2L`, forwarded to
#'   `dina_subgroup()` for both the point estimate and every bootstrap
#'   iteration.  `2L` allows the selected subgroup to be an AND-conjunction
#'   of two covariates; see [dina_subgroup()].  When `sg` is supplied, the
#'   per-iteration search follows the depth of `sg` (its `max_depth`).
#' @param grid_probs numeric vector of probabilities in `[0, 1]`; the
#'   per-covariate quantile grid for depth-2 pair thresholds, forwarded to
#'   `dina_subgroup()`.  Only consulted when the effective depth is `2`.
#'   Default interior deciles `seq(0.1, 0.9, 0.1)`.
#' @param sg_focus character; subgroup selection criterion, forwarded to
#'   `dina_subgroup()` for both the original-data point estimate and every
#'   bootstrap iteration.  One of `"maxSG"` (default), `"minSG"`, `"eff"`,
#'   `"effMaxSG"`, `"effMinSG"` (canonical `"hr"`, `"hrMaxSG"`,
#'   `"hrMinSG"` also accepted).  See [dina_subgroup()].
#' @param effect_neighborhood numeric in `[0, 1)`; relative tolerance for
#'   the `"effMaxSG"` / `"effMinSG"` effect band, forwarded to
#'   `dina_subgroup()`.  Default `0.10`.  Ignored for the non-band foci.
#' @param alpha confidence level for the percentile CI.  Default `0.05`
#'   (percentile endpoints at the `alpha/2` and `1 - alpha/2` quantiles).
#' @param n_boot positive integer; number of bootstrap iterations.
#'   Default `200L`.  Warns if `n_boot < 100`.
#' @param parallel one of `"none"` (default) or `"boots"`.  When
#'   `"boots"`, bootstrap iterations dispatch via
#'   `foreach::foreach() %dofuture%` against the active
#'   `future::plan()`.
#' @param seed optional integer seed for reproducibility of the
#'   bootstrap indices.  Default `NULL`.
#' @param verbose logical; if `TRUE`, print a progress message every
#'   10 iterations (serial mode) or a start summary (parallel mode).
#'   Default `FALSE`.
#' @param fit optional precomputed DINA fit (`"dina"` or `"dina_bagged"`)
#'   on the original `df`.  When supplied (and `sg` is `NULL`), the
#'   function reuses it for the original-data point estimate instead of
#'   refitting, then runs only the subgroup search.  Its `family` must
#'   match `family`.  Default `NULL`.
#' @param sg optional precomputed `"dina_subgroup"` result on the original
#'   `df`.  When supplied, it is used directly as the point estimate with
#'   no refit and no re-search; `fit` is then ignored.  Its stored
#'   `m_diff`, `n_min`, `alpha`, `family`, and row count must match the
#'   corresponding arguments / `df`, or an error is raised.  Supplying the
#'   upstream `sg` is the recommended way to guarantee the bootstrap's
#'   point estimate is identical to a standalone `dina_subgroup()` result.
#'   Default `NULL`.
#' @param refit logical; if \code{TRUE} (default), also fit the
#'   \emph{standard} within-subgroup treatment-effect model (Cox for
#'   survival, GLM otherwise) via \code{\link{dina_subgroup_refit}} on
#'   the original data, and bootstrap it within the fixed discovered
#'   signature.  Reported alongside the DINA effect.  Skipped (with a
#'   warning) if the original-data subgroup is not found or the
#'   within-subgroup fit errors.
#' @param refit_strata \code{NULL} (default) or a character vector of
#'   column names entered as Cox \code{strata()} terms in the
#'   within-subgroup model.  Cox family only.
#' @param refit_confounders adjustment set for the within-subgroup
#'   model: \code{"none"} (default, unadjusted), \code{NULL}
#'   (automatic: the DINA \code{covariates} with the subgroup-defining
#'   covariate omitted), or a character vector of column names.  See
#'   \code{\link{dina_subgroup_refit}}.
#' @param ... further arguments forwarded to `dina()` for each
#'   per-iteration DINA fit (e.g., `n_folds`, `cens_type`,
#'   `cens_params`, `propensity_method`, `baseline_method`).
#'
#' @return An object of class `"dina_subgroup_bootstrap"`, a list with
#'   components:
#'   \describe{
#'     \item{point}{the `"dina_subgroup"` object computed on the
#'       original (unbootstrapped) data.}
#'     \item{effect_ci}{named length-2 numeric vector (`lower`,
#'       `upper`) giving the bootstrap percentile CI on the
#'       fixed-subgroup linear functional `a*^T beta_b`, where
#'       `a* = (1, x_bar_S)` is the covariate-mean vector of the
#'       original-data subgroup and `beta_b` is each bootstrap
#'       iteration's DINA coefficient vector.  This is the
#'       meaningful effect-uncertainty CI; width reflects sampling
#'       variance of `beta-hat`.}
#'     \item{effect_dist}{numeric vector of length `n_boot`
#'       containing `a*^T beta_b` for each iteration (`NA` for
#'       iterations whose DINA fit failed).}
#'     \item{refit}{the \code{"dina_subgroup_refit"} object (standard
#'       within-subgroup model on the original data), or `NULL` if
#'       `refit = FALSE`, no subgroup was found, or the fit failed.}
#'     \item{refit_effect_ci}{named length-2 numeric vector (`lower`,
#'       `upper`); bootstrap percentile CI on the within-subgroup
#'       standard-model treatment effect, refit within the FIXED
#'       discovered signature on each resample.  Conditional on the
#'       signature; comparable to `effect_ci`.  `(NA, NA)` if `refit`
#'       was not performed.}
#'     \item{refit_effect_dist}{numeric vector of length `n_boot` of
#'       per-iteration within-signature standard-model treatment
#'       effects (`NA` where not computed or the fit failed).}
#'     \item{n_subgroup_ci}{named length-2 numeric vector; percentile CI
#'       for the bootstrap-selected subgroup size, restricted to
#'       iterations matching the original-data **cut-set** (same
#'       covariate(s) and direction(s), via a threshold-independent key).}
#'     \item{threshold_ci}{boundary-location stability, restricted to the
#'       same matched iterations.  For a single-covariate subgroup, a
#'       named length-2 numeric vector.  For a depth-2 conjunction, a
#'       matrix with one row per covariate (rownames = covariate) and
#'       columns `lower`/`upper`, since the two thresholds are not
#'       comparable on a common scale.}
#'     \item{n_modal_match}{integer; number of bootstrap iterations whose
#'       selection matched the original-data cut-set key.  Denominator for
#'       the stability CIs.}
#'     \item{selection_frequency}{a one-dimensional `table` of canonical
#'       cut-set keys (e.g. `"nodes >="`, or `"age >= & nodes >="` for a
#'       conjunction) across found iterations.  The structural-selection
#'       stability diagnostic, generalizing the depth-1
#'       `(covariate, direction)` table.}
#'     \item{boot_results}{data frame with one row per bootstrap
#'       iteration; columns `found`, `key` (canonical cut-set key),
#'       `depth` (1 or 2), `n_subgroup`, `mean_tau_hat`, `failed`, plus
#'       back-compatible `covariate`, `direction`, `threshold` columns
#'       (the scalar cut for single-covariate selections, `NA` for depth-2
#'       conjunctions -- whose per-iteration cuts are summarized by `key`).
#'       The `mean_tau_hat` column records the bootstrap-selected
#'       subgroup's boundary value (pinned near `m_diff` by the search
#'       rule) and is retained for diagnostic transparency rather than as
#'       a meaningful effect estimate.}
#'     \item{n_boot, n_boot_found, n_boot_failed}{convergence
#'       diagnostics.}
#'     \item{alpha, m_diff, n_min, sg_focus, effect_neighborhood,
#'       family}{echoed inputs.}
#'     \item{call}{the matched call.}
#'   }
#'
#' @seealso `dina_subgroup()` for the underlying threshold search;
#'   `dina()` for the per-iteration DINA fit;
#'   \code{\link{dina_subgroup_refit}} for the standard within-subgroup
#'   model reported when `refit = TRUE`.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' n <- 400
#' df_demo <- data.frame(
#'   w  = stats::rbinom(n, 1, 0.5),
#'   x1 = stats::runif(n, -1, 1),
#'   x2 = stats::runif(n, -1, 1)
#' )
#' tau_x      <- 0.4 + 1.2 * df_demo$x1
#' df_demo$y  <- 0.5 * df_demo$x1 + df_demo$w * tau_x + stats::rnorm(n)
#'
#' bs <- dina_subgroup_bootstrap(
#'   df         = df_demo,
#'   outcome    = "y",
#'   treatment  = "w",
#'   covariates = c("x1", "x2"),
#'   family     = "gaussian",
#'   m_diff     = 0.5,
#'   n_boot     = 50L,
#'   seed       = 1L
#' )
#' print(bs)
#' }
#'
#' @importFrom foreach foreach
#' @importFrom future nbrOfWorkers
#' @export
dina_subgroup_bootstrap <- function(df,
                                    outcome,
                                    treatment,
                                    covariates,
                                    family = c("gaussian", "binomial",
                                               "poisson", "cox"),
                                    status = NULL,
                                    m_diff,
                                    n_min = 60L,
                                    direction = c("both", "left", "right"),
                                    max_depth = 1L,
                                    grid_probs = seq(0.1, 0.9, by = 0.1),
                                    sg_focus = "maxSG",
                                    effect_neighborhood = 0.10,
                                    alpha = 0.05,
                                    n_boot = 200L,
                                    parallel = c("none", "boots"),
                                    seed = NULL,
                                    verbose = FALSE,
                                    fit = NULL,
                                    sg = NULL,
                                    refit = TRUE,
                                    refit_strata = NULL,
                                    refit_confounders = "none",
                                    ...) {

  family    <- match.arg(family)
  parallel  <- match.arg(parallel)
  direction <- match.arg(direction)
  call      <- match.call()

  # ---- Argument validation ----------------------------------------------
  if (missing(m_diff) || length(m_diff) != 1L || !is.numeric(m_diff) ||
      !is.finite(m_diff)) {
    stop("`m_diff` must be a single finite numeric value.")
  }
  n_boot <- as.integer(n_boot)
  if (length(n_boot) != 1L || is.na(n_boot) || n_boot < 2L) {
    stop("`n_boot` must be a single integer >= 2.")
  }
  if (n_boot < 100L) {
    warning("`n_boot` < 100; bootstrap percentile CI will be unstable. ",
            "Consider increasing `n_boot`.", call. = FALSE)
  }
  if (length(alpha) != 1L || !is.numeric(alpha) ||
      alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be a single numeric in (0, 1).")
  }

  # Normalize the GLM-natural sg_focus vocabulary to the canonical form
  # (matches dina_subgroup()), then whitelist.  The per-iteration search
  # and the original-data point estimate both use the canonical value.
  sg_focus <- .normalize_sg_focus(sg_focus)
  valid_sg_focus <- c("hr", "maxSG", "minSG", "hrMaxSG", "hrMinSG")
  if (!is.character(sg_focus) || length(sg_focus) != 1L ||
      !sg_focus %in% valid_sg_focus) {
    stop("`sg_focus` must be one of \"maxSG\", \"minSG\", \"eff\", ",
         "\"effMaxSG\", \"effMinSG\" (canonical forms \"hr\", ",
         "\"hrMaxSG\", \"hrMinSG\" are also accepted).")
  }
  .validate_effect_neighborhood(effect_neighborhood)

  max_depth <- as.integer(max_depth)
  if (length(max_depth) != 1L || is.na(max_depth) ||
      !max_depth %in% c(1L, 2L)) {
    stop("`max_depth` must be 1 or 2.")
  }
  if (max_depth == 2L &&
      (!is.numeric(grid_probs) || length(grid_probs) < 1L ||
       anyNA(grid_probs) || any(grid_probs < 0) || any(grid_probs > 1))) {
    stop("`grid_probs` must be a numeric vector in [0, 1] with no NAs.")
  }

  if (parallel == "boots") {
    n_workers <- future::nbrOfWorkers()
    if (n_workers <= 1L) {
      warning(
        "`parallel = \"boots\"` requested but `future::plan()` reports ",
        "only 1 worker; bootstraps will run serially. ",
        "Call e.g. `future::plan(future::multisession, workers = N)` ",
        "before this function to enable actual parallelism.",
        call. = FALSE
      )
    } else if (verbose) {
      message(sprintf(
        "dina_subgroup_bootstrap: running %d bootstraps across %d workers",
        n_boot, n_workers
      ))
    }
  }

  dina_args <- list(...)
  # Prevent the user's `seed` (if passed via ...) from being applied to
  # every per-iteration DINA fit identically.  The bootstrap-index seed
  # is governed by the explicit `seed` argument above; the per-iteration
  # fits are intentionally unseeded so each resample gets its own
  # cross-fitting fold assignment.
  dina_args[["seed"]] <- NULL

  # ---- Point estimate on the original data ------------------------------
  # Obtain the original-data point estimate by one of three routes, in
  # priority order:
  #
  #   1. `sg` supplied  -> use it directly (no refit, no re-search).
  #   2. `fit` supplied -> reuse the fit; run the (cheap) subgroup search.
  #   3. neither        -> refit DINA on the original data, forwarding the
  #                        explicit `seed` so the point fit is reproducible.
  #
  # Passing the upstream `fit` and/or `sg` makes the bootstrap's point
  # estimate IDENTICAL to the caller's standalone result by construction.
  # This is the recommended usage: it removes the cross-fitting
  # fold-assignment mismatch that otherwise arises when the point estimate
  # is refit here with a different (or absent) seed than the caller used,
  # which at a boundary-case `m_diff` can flip the original-data subgroup
  # between found and not-found.
  if (!is.null(sg)) {
    if (!inherits(sg, "dina_subgroup")) {
      stop("`sg` must be a \"dina_subgroup\" object as returned by ",
           "dina_subgroup().")
    }
    if (!is.null(sg$n_total) && !identical(as.integer(sg$n_total),
                                           as.integer(nrow(df)))) {
      stop("`sg` was computed on ", sg$n_total, " rows but `df` has ",
           nrow(df), " rows; `sg` must be computed on the same `df`.")
    }
    if (isTRUE(sg$found) && length(sg$mask) != nrow(df)) {
      stop("`sg$mask` length (", length(sg$mask), ") does not match ",
           "nrow(df) (", nrow(df), "); `sg` must be computed on the ",
           "same `df`.")
    }
    if (!isTRUE(all.equal(sg$m_diff, m_diff))) {
      stop("`sg$m_diff` (", format(sg$m_diff), ") does not match the ",
           "`m_diff` argument (", format(m_diff), ").")
    }
    if (!identical(as.integer(sg$n_min), as.integer(n_min))) {
      stop("`sg$n_min` (", sg$n_min, ") does not match the `n_min` ",
           "argument (", n_min, ").")
    }
    if (!isTRUE(all.equal(sg$alpha, alpha))) {
      stop("`sg$alpha` (", format(sg$alpha), ") does not match the ",
           "`alpha` argument (", format(alpha), ").")
    }
    if (!is.null(sg$family) && !identical(sg$family, family)) {
      stop("`sg$family` (", sg$family, ") does not match the `family` ",
           "argument (", family, ").")
    }
    # The supplied point estimate must have used the same selection rule
    # the per-iteration search will use, or the original-data and
    # bootstrap subgroups would be defined by different criteria.
    if (!is.null(sg$sg_focus) && !identical(sg$sg_focus, sg_focus)) {
      stop("`sg$sg_focus` (", sg$sg_focus, ") does not match the resolved ",
           "`sg_focus` (", sg_focus, ").")
    }
    if (sg_focus %in% c("hrMaxSG", "hrMinSG") &&
        !is.null(sg$effect_neighborhood) &&
        !isTRUE(all.equal(sg$effect_neighborhood, effect_neighborhood))) {
      stop("`sg$effect_neighborhood` (", format(sg$effect_neighborhood),
           ") does not match the `effect_neighborhood` argument (",
           format(effect_neighborhood), ").")
    }
    sg_point <- sg
  } else if (!is.null(fit)) {
    if (!inherits(fit, "dina")) {
      stop("`fit` must be a \"dina\" (or \"dina_bagged\") object as ",
           "returned by dina() / dina_bagged().")
    }
    if (!is.null(fit$family) && !identical(fit$family, family)) {
      stop("`fit$family` (", fit$family, ") does not match the `family` ",
           "argument (", family, ").")
    }
    sg_point <- dina_subgroup(
      fit = fit, df = df, covariates = covariates,
      m_diff = m_diff, n_min = n_min,
      direction = direction, max_depth = max_depth, grid_probs = grid_probs,
      sg_focus = sg_focus,
      effect_neighborhood = effect_neighborhood, alpha = alpha
    )
  } else {
    fit_point <- do.call(dina, c(
      list(df = df, outcome = outcome, treatment = treatment,
           covariates = covariates, family = family, status = status,
           seed = seed),
      dina_args
    ))
    sg_point <- dina_subgroup(
      fit = fit_point, df = df, covariates = covariates,
      m_diff = m_diff, n_min = n_min,
      direction = direction, max_depth = max_depth, grid_probs = grid_probs,
      sg_focus = sg_focus,
      effect_neighborhood = effect_neighborhood, alpha = alpha
    )
  }

  # Per-iteration selection must use the same depth as the point estimate.
  # When `sg` was supplied it may have been computed at max_depth = 2, so
  # follow its depth; otherwise the point estimate above used `max_depth`.
  # `grid_probs` is not stored on the subgroup object, so it always comes
  # from the argument (default deciles).
  sel_max_depth <- if (!is.null(sg_point$max_depth))
                     as.integer(sg_point$max_depth) else max_depth
  sel_grid_probs <- grid_probs
  sg_is_depth2   <- isTRUE(sg_point$found) && length(sg_point$covariate) > 1L

  # ---- Within-subgroup standard-model point estimate --------------------
  # The DINA effect above is the BLP-analog a*^T beta.  Clinical reporting
  # usually also wants the *standard* treatment-effect model fit within the
  # discovered subgroup (a plain Cox model for survival, a GLM otherwise).
  # Compute it on the original data via dina_subgroup_refit(); the fixed
  # signature and resolved adjustment set are reused per bootstrap
  # iteration below so the bootstrap CI is conditional on the SAME
  # signature -- coherent with, and directly comparable to, the DINA
  # effect_ci.  Neither CI is selection-adjusted.
  do_refit          <- isTRUE(refit) && isTRUE(sg_point$found)
  if (do_refit && sg_is_depth2) {
    warning("`refit = TRUE` is not yet supported for depth-2 (conjunction) ",
            "subgroups; skipping the within-subgroup standard-model refit. ",
            "The DINA effect_ci and stability diagnostics are unaffected.",
            call. = FALSE)
    do_refit <- FALSE
  }
  refit_point       <- NULL
  refit_conf_used   <- character(0)
  refit_signature   <- NULL

  if (do_refit) {
    refit_signature <- list(
      covariate = sg_point$covariate,
      direction = sg_point$direction,
      threshold = sg_point$threshold
    )
    refit_conf_used <- .dina_resolve_confounders(
      confounders  = refit_confounders,
      covariates   = covariates,
      sg_covariate = sg_point$covariate
    )
    refit_point <- tryCatch(
      dina_subgroup_refit(
        sg = sg_point, df = df, treatment = treatment, outcome = outcome,
        covariates = covariates, status = status,
        strata = refit_strata, confounders = refit_confounders,
        alpha = alpha
      ),
      error = function(e) {
        warning("Within-subgroup standard-model point estimate failed (",
                conditionMessage(e),
                "); skipping the within-subgroup refit.", call. = FALSE)
        NULL
      }
    )
    if (is.null(refit_point)) do_refit <- FALSE
  }

  # ---- Bootstrap indices ------------------------------------------------
  n <- nrow(df)
  boot <- .dina_make_bootstrap_indices(n, n_boot, seed = seed)

  # ---- Bootstrap loop ---------------------------------------------------
  if (parallel == "boots") {
    b <- NULL  # silence R CMD check "no visible binding"
    boot_list <- foreach::foreach(
      b = seq_len(n_boot),
      .options.future = list(
        packages = "forestsearch",
        seed     = TRUE
      )
    ) %dofuture% {
      .dina_safe_one_bootstrap(
        idx = boot$indices[, b],
        df = df, outcome = outcome, treatment = treatment,
        covariates = covariates, family = family, status = status,
        m_diff = m_diff, n_min = n_min, direction = direction,
        max_depth = sel_max_depth, grid_probs = sel_grid_probs,
        sg_focus = sg_focus, effect_neighborhood = effect_neighborhood,
        alpha = alpha, dina_args = dina_args,
        refit = do_refit, refit_signature = refit_signature,
        refit_confounders = refit_conf_used, refit_strata = refit_strata
      )
    }
  } else {
    boot_list <- lapply(seq_len(n_boot), function(b) {
      if (verbose && b %% 10L == 0L) {
        message(sprintf("dina_subgroup_bootstrap: iteration %d of %d",
                        b, n_boot))
      }
      .dina_safe_one_bootstrap(
        idx = boot$indices[, b],
        df = df, outcome = outcome, treatment = treatment,
        covariates = covariates, family = family, status = status,
        m_diff = m_diff, n_min = n_min, direction = direction,
        max_depth = sel_max_depth, grid_probs = sel_grid_probs,
        sg_focus = sg_focus, effect_neighborhood = effect_neighborhood,
        alpha = alpha, dina_args = dina_args,
        refit = do_refit, refit_signature = refit_signature,
        refit_confounders = refit_conf_used, refit_strata = refit_strata
      )
    })
  }

  # ---- Aggregate --------------------------------------------------------
  # Per-iteration bootstrap-selected subgroup statistics (for stability
  # diagnostics).  The effect-uncertainty target -- the fixed-subgroup
  # linear functional a*^T beta_b -- is computed separately below using
  # each iteration's full coefficient vector and the FIXED a* from the
  # original-data subgroup.
  #
  # Each iteration's selected cut(s) are carried in `boot_cuts` (a list of
  # covariate/direction/threshold vectors, length 1 for a single cut, 2 for
  # a conjunction), with a scalar canonical `key` for structural matching.
  boot_cuts <- lapply(boot_list, function(r) {
    if (is.null(r$cuts)) {
      list(covariate = character(0), direction = character(0),
           threshold = numeric(0))
    } else r$cuts
  })

  boot_df <- data.frame(
    found        = vapply(boot_list, function(r) {
                    if (is.na(r$found)) NA else as.logical(r$found)
                   }, logical(1L)),
    key          = vapply(boot_list, function(r) {
                    if (is.null(r$key)) NA_character_ else r$key
                   }, character(1L)),
    depth        = vapply(boot_list, function(r) {
                    if (is.null(r$depth)) NA_integer_ else as.integer(r$depth)
                   }, integer(1L)),
    n_subgroup   = vapply(boot_list, `[[`, integer(1L), "n_subgroup"),
    mean_tau_hat = vapply(boot_list, `[[`, numeric(1L), "mean_tau_hat"),
    failed       = vapply(boot_list, `[[`, logical(1L), "failed"),
    stringsAsFactors = FALSE
  )
  # Back-compatible scalar cut columns: populated for single-covariate
  # selections (identical to the legacy depth-1 output); NA for depth-2
  # conjunctions, whose cuts live in `key` and `boot_cuts`.
  boot_df$covariate <- vapply(boot_cuts, function(c)
    if (length(c$covariate) == 1L) c$covariate else NA_character_,
    character(1L))
  boot_df$direction <- vapply(boot_cuts, function(c)
    if (length(c$direction) == 1L) c$direction else NA_character_,
    character(1L))
  boot_df$threshold <- vapply(boot_cuts, function(c)
    if (length(c$threshold) == 1L) c$threshold else NA_real_, numeric(1L))

  n_boot_failed <- sum(boot_df$failed)
  fit_mask      <- !boot_df$failed
  found_mask    <- !is.na(boot_df$found) & boot_df$found
  n_boot_found  <- sum(found_mask)

  # ---- Fixed-subgroup effect distribution -------------------------------
  # The point-estimate subgroup S* and its covariate-mean vector a* are
  # held FIXED from the original-data fit.  For each bootstrap iteration
  # whose DINA fit succeeded, compute T_b = a*^T beta_b.  Per the design
  # rationale (see dina_inference_design_notes.md, Section 4), this is the
  # meaningful effect-uncertainty target -- width reflects sampling
  # variance of beta-hat, not boundary slack.
  effect_dist <- rep(NA_real_, n_boot)
  effect_ci   <- c(lower = NA_real_, upper = NA_real_)

  if (isTRUE(sg_point$found)) {
    X_all <- .dina_extract_X_from_df(df, covariates)
    x_bar_S <- colMeans(X_all[sg_point$mask, , drop = FALSE])
    a_star  <- c(1, x_bar_S)

    # Stack per-iteration beta vectors into a B x (d+1) matrix; rows for
    # failed iterations remain NA.  Matrix-vector product yields a_star^T
    # beta_b for every iteration in one shot.
    beta_mat <- do.call(rbind, lapply(boot_list, `[[`, "beta"))
    effect_dist <- as.numeric(beta_mat %*% a_star)

    effect_finite <- effect_dist[is.finite(effect_dist)]
    if (length(effect_finite) >= 2L) {
      qs <- stats::quantile(
        effect_finite, probs = c(alpha / 2, 1 - alpha / 2),
        names = FALSE, na.rm = TRUE
      )
      effect_ci <- c(lower = qs[1L], upper = qs[2L])
    }
  }

  # ---- Within-subgroup standard-model effect distribution ---------------
  # Each iteration refit the standard model (Cox / GLM) within the FIXED
  # original signature applied to that resample; refit_effect is the
  # treatment coefficient (NA if the within-signature fit was not done or
  # failed for that resample, e.g. an empty arm).  Percentile CI over the
  # finite values gives a within-signature bootstrap CI directly
  # comparable to the DINA effect_ci above.  Conditional on the signature;
  # not selection-adjusted.
  refit_effect_dist <- vapply(boot_list, `[[`, numeric(1L), "refit_effect")
  refit_effect_ci   <- c(lower = NA_real_, upper = NA_real_)
  if (do_refit) {
    refit_finite <- refit_effect_dist[is.finite(refit_effect_dist)]
    if (length(refit_finite) >= 2L) {
      rqs <- stats::quantile(
        refit_finite, probs = c(alpha / 2, 1 - alpha / 2),
        names = FALSE, na.rm = TRUE
      )
      refit_effect_ci <- c(lower = rqs[1L], upper = rqs[2L])
    }
  }

  # ---- Stability CIs on subgroup structure ------------------------------
  # Restricted to iterations that selected the SAME structural cut-set as
  # the original-data subgroup -- same covariate(s) and direction(s),
  # regardless of thresholds -- via the canonical key.  This generalizes
  # the depth-1 "(covariate, direction)" match to conjunctions.  Threshold
  # CIs are reported per component (one per covariate in the matched
  # cut-set), since thresholds on different covariates are not comparable;
  # for a single-covariate subgroup this is the usual length-2 CI.
  key_point <- if (isTRUE(sg_point$found))
                 .dina_cutset_key(sg_point$covariate, sg_point$direction)
               else NA_character_

  n_subgroup_ci <- c(lower = NA_real_, upper = NA_real_)
  threshold_ci  <- c(lower = NA_real_, upper = NA_real_)
  n_modal_match <- NA_integer_

  if (isTRUE(sg_point$found) && n_boot_found > 0L) {
    modal_match <- found_mask & !is.na(boot_df$key) &
      boot_df$key == key_point
    n_modal_match <- as.integer(sum(modal_match))

    if (n_modal_match >= 2L) {
      n_sg_qs <- stats::quantile(
        boot_df$n_subgroup[modal_match],
        probs = c(alpha / 2, 1 - alpha / 2),
        names = FALSE, na.rm = TRUE
      )
      n_subgroup_ci <- c(lower = n_sg_qs[1L], upper = n_sg_qs[2L])

      # Align each matched iteration's thresholds to the point subgroup's
      # covariates (matched iterations share the same covariate set, so the
      # match() never misses).
      covs_pt <- sg_point$covariate
      matched <- which(modal_match)
      thr_mat <- matrix(NA_real_, nrow = length(matched),
                        ncol = length(covs_pt))
      for (r in seq_along(matched)) {
        cc <- boot_cuts[[matched[r]]]
        thr_mat[r, ] <- cc$threshold[match(covs_pt, cc$covariate)]
      }
      per_comp <- t(vapply(seq_len(ncol(thr_mat)), function(j) {
        stats::quantile(thr_mat[, j], probs = c(alpha / 2, 1 - alpha / 2),
                        names = FALSE, na.rm = TRUE)
      }, numeric(2L)))
      if (length(covs_pt) == 1L) {
        threshold_ci <- c(lower = per_comp[1L, 1L], upper = per_comp[1L, 2L])
      } else {
        dimnames(per_comp) <- list(covs_pt, c("lower", "upper"))
        threshold_ci <- per_comp
      }
    }
  }

  # ---- Selection-frequency diagnostic -----------------------------------
  # Counts of the canonical cut-set key across found iterations -- the
  # depth-agnostic generalization of the (covariate, direction) table.
  selection_frequency <- if (n_boot_found > 0L) {
    table(selection = boot_df$key[found_mask])
  } else {
    NULL
  }

  if (n_boot_failed > 0L) {
    warning(sprintf("%d of %d bootstrap iterations errored and were excluded.",
                    n_boot_failed, n_boot), call. = FALSE)
  }
  if (sum(fit_mask) < n_boot %/% 2L) {
    warning(sprintf(
      "Only %d of %d bootstrap iterations produced a usable DINA fit; ",
      sum(fit_mask), n_boot
    ), "effect CI is based on a small sample.", call. = FALSE)
  }

  out <- list(
    point                = sg_point,
    effect_ci            = effect_ci,
    effect_dist          = effect_dist,
    refit                = refit_point,
    refit_effect_ci      = refit_effect_ci,
    refit_effect_dist    = refit_effect_dist,
    n_subgroup_ci        = n_subgroup_ci,
    threshold_ci         = threshold_ci,
    n_modal_match        = n_modal_match,
    selection_frequency  = selection_frequency,
    boot_results         = boot_df,
    n_boot               = n_boot,
    n_boot_found         = n_boot_found,
    n_boot_failed        = n_boot_failed,
    alpha                = alpha,
    m_diff               = m_diff,
    n_min                = n_min,
    sg_focus             = sg_focus,
    effect_neighborhood  = effect_neighborhood,
    family               = family,
    call                 = call
  )
  class(out) <- "dina_subgroup_bootstrap"
  out
}


# ---------------------------------------------------------------------------
# Internal helper: canonical cut-set key for selection stability
# ---------------------------------------------------------------------------

#' Canonical, threshold-independent key for a (possibly multi-covariate)
#' subgroup selection.
#'
#' Two selections share a key iff they use the same set of covariates with
#' the same per-covariate directions, regardless of thresholds (which are
#' the quantity whose stability is being measured) and regardless of cut
#' order.  Depth-1 keys look like \code{"nodes >="}; depth-2 keys like
#' \code{"age >= & nodes >="} (component cuts sorted for order-independence).
#' Returns \code{NA_character_} for an empty selection.
#'
#' @noRd
.dina_cutset_key <- function(covariate, direction) {
  if (length(covariate) == 0L) return(NA_character_)
  op <- ifelse(direction == "left", "<=", ">=")
  paste(sort(paste0(covariate, " ", op)), collapse = " & ")
}


# ---------------------------------------------------------------------------
# Internal helper: one bootstrap iteration with error handling
# ---------------------------------------------------------------------------

#' Run one bootstrap iteration with error handling, returning a uniform
#' list result.
#'
#' Two independent pieces of work per resample, each in its own tryCatch
#' so one failing does not lose the other:
#'   (1) refit DINA + run dina_subgroup() -> beta_b and the resample's
#'       selected-subgroup stats;
#'   (2) when `refit` is TRUE, apply the FIXED original signature
#'       (`refit_signature`) to the resample and fit the standard
#'       within-subgroup model via .dina_fit_subgroup_model(), returning
#'       its treatment coefficient.
#' Failures in either surface as NA-valued fields (and, for (1),
#' `failed = TRUE`) rather than propagating as conditions.  Both parallel
#' and serial code paths consume the same list contract downstream.
#'
#' @noRd
.dina_safe_one_bootstrap <- function(idx, df, outcome, treatment, covariates,
                                     family, status,
                                     m_diff, n_min, direction,
                                     max_depth = 1L,
                                     grid_probs = seq(0.1, 0.9, by = 0.1),
                                     sg_focus = "maxSG",
                                     effect_neighborhood = 0.10, alpha,
                                     dina_args,
                                     refit = FALSE, refit_signature = NULL,
                                     refit_confounders = character(0),
                                     refit_strata = NULL) {
  d <- length(covariates)
  na_beta <- rep(NA_real_, d + 1L)

  # df_b is plain row-indexing; it cannot fail.
  df_b <- df[idx, , drop = FALSE]

  # ---- (1) DINA refit + subgroup search (independent tryCatch) ----------
  dina_part <- tryCatch({
    fit_b <- do.call(dina, c(
      list(df = df_b, outcome = outcome, treatment = treatment,
           covariates = covariates, family = family, status = status),
      dina_args
    ))
    beta_b <- as.numeric(stats::coef(fit_b))
    sg_b <- dina_subgroup(
      fit = fit_b, df = df_b, covariates = covariates,
      m_diff = m_diff, n_min = n_min,
      direction = direction, max_depth = max_depth, grid_probs = grid_probs,
      sg_focus = sg_focus,
      effect_neighborhood = effect_neighborhood, alpha = alpha
    )
    list(beta = beta_b, sg_b = sg_b, failed = FALSE, message = "")
  }, error = function(e) {
    list(beta = na_beta, sg_b = NULL, failed = TRUE,
         message = conditionMessage(e))
  })

  # ---- (2) Within-signature standard model (independent tryCatch) -------
  # Apply the fixed original signature to this resample, then fit the
  # standard Cox/GLM model on that subset.  Independent of the DINA refit
  # above, so a DINA failure does not suppress this and vice versa.
  refit_effect <- NA_real_
  if (isTRUE(refit) && !is.null(refit_signature)) {
    refit_effect <- tryCatch({
      mask_b <- .dina_signature_mask(
        df_b, refit_signature$covariate,
        refit_signature$direction, refit_signature$threshold
      )
      sub_b <- df_b[mask_b, , drop = FALSE]
      fm <- .dina_fit_subgroup_model(
        df_sub = sub_b, treatment = treatment, outcome = outcome,
        status = status, family = family,
        confounders = refit_confounders, strata = refit_strata
      )
      fm$effect
    }, error = function(e) NA_real_)
  }

  # ---- Assemble the uniform per-iteration contract ----------------------
  sg_b  <- dina_part$sg_b
  found <- if (isTRUE(dina_part$failed)) NA else isTRUE(sg_b$found)

  if (isTRUE(found)) {
    list(
      beta         = dina_part$beta,
      found        = TRUE,
      cuts         = list(covariate = sg_b$covariate,
                          direction = sg_b$direction,
                          threshold = as.numeric(sg_b$threshold)),
      key          = .dina_cutset_key(sg_b$covariate, sg_b$direction),
      depth        = length(sg_b$covariate),
      n_subgroup   = as.integer(sg_b$n_subgroup),
      mean_tau_hat = as.numeric(sg_b$mean_tau_hat),
      refit_effect = refit_effect,
      failed       = FALSE,
      message      = ""
    )
  } else if (isFALSE(found)) {
    list(
      beta         = dina_part$beta,
      found        = FALSE,
      cuts         = list(covariate = character(0), direction = character(0),
                          threshold = numeric(0)),
      key          = NA_character_,
      depth        = NA_integer_,
      n_subgroup   = NA_integer_,
      mean_tau_hat = NA_real_,
      refit_effect = refit_effect,
      failed       = FALSE,
      message      = ""
    )
  } else {
    list(
      beta         = na_beta,
      found        = NA,
      cuts         = list(covariate = character(0), direction = character(0),
                          threshold = numeric(0)),
      key          = NA_character_,
      depth        = NA_integer_,
      n_subgroup   = NA_integer_,
      mean_tau_hat = NA_real_,
      refit_effect = refit_effect,
      failed       = TRUE,
      message      = dina_part$message
    )
  }
}


# ---------------------------------------------------------------------------
# S3 methods
# ---------------------------------------------------------------------------

#' Print a `dina_subgroup_bootstrap` result.
#'
#' @param x a `"dina_subgroup_bootstrap"` object.
#' @param digits number of digits for numeric summary.
#' @param max_select maximum number of rows of the selection-frequency
#'   table to display.  Default `5L`.
#' @param ... unused.
#' @return invisibly returns `x`.
#' @export
print.dina_subgroup_bootstrap <- function(x,
                                          digits = max(3L,
                                                       getOption("digits") - 3L),
                                          max_select = 5L,
                                          ...) {
  cat("Bootstrap inference for dina_subgroup()\n")
  cat("  Family:           ", x$family, "\n", sep = "")
  sgf <- if (is.null(x$sg_focus)) "maxSG" else x$sg_focus
  focus_label <- switch(
    sgf,
    hr      = "eff (most extreme effect)",
    maxSG   = "maxSG (largest qualifying subgroup)",
    minSG   = "minSG (smallest qualifying subgroup)",
    hrMaxSG = "effMaxSG (largest within effect band)",
    hrMinSG = "effMinSG (smallest within effect band)",
    sgf
  )
  cat("  Focus:            ", focus_label, "\n", sep = "")
  if (isTRUE(sgf %in% c("hrMaxSG", "hrMinSG"))) {
    en <- if (is.null(x$effect_neighborhood)) 0.10 else x$effect_neighborhood
    cat("  Effect band:      within ",
        format(100 * en, digits = digits),
        "% of the maximum qualifying effect\n", sep = "")
  }
  cat("  m_diff:           ", format(x$m_diff, digits = digits), "\n", sep = "")
  cat("  n_min:            ", x$n_min, "\n", sep = "")
  cat("  Bootstrap iters:  ", x$n_boot_found, " found / ",
      sum(!is.na(x$boot_results$found) & !x$boot_results$found),
      " not-found / ", x$n_boot_failed, " failed (of ",
      x$n_boot, " total)\n", sep = "")

  cat("\nPoint estimate (original data):\n")
  if (isTRUE(x$point$found)) {
    ops <- ifelse(x$point$direction == "left", "<=", ">=")
    sig <- paste(paste0(x$point$covariate, " ", ops, " ",
                        format(x$point$threshold, digits = digits)),
                 collapse = " & ")
    cat("  Signature:    ", sig, "\n", sep = "")
    cat("  n_subgroup:   ", x$point$n_subgroup, " / ",
        x$point$n_total, "\n", sep = "")
    cat("  mean tau-hat: ", format(x$point$mean_tau_hat, digits = digits),
        "\n", sep = "")
  } else {
    cat("  No qualifying subgroup on the original data.\n")
  }

  ci_pct <- format(100 * (1 - x$alpha), digits = 3L)

  cat("\nFixed-subgroup effect CI (", ci_pct,
      " pct percentile):  [",
      format(x$effect_ci[["lower"]], digits = digits), ", ",
      format(x$effect_ci[["upper"]], digits = digits), "]\n", sep = "")
  cat("  Target:       a*^T beta_b with a* = (1, x_bar_S) fixed from\n")
  cat("                the original-data subgroup S*.\n")
  cat("  Interpretation: bootstrap percentile CI for the effect in\n")
  cat("                the discovered subgroup; complements the Wald\n")
  cat("                CI returned by dina_subgroup().\n")

  # ---- Within-subgroup standard model -----------------------------------
  if (!is.null(x$refit)) {
    rf <- x$refit
    cat("\nWithin-subgroup standard model (", rf$effect_scale, "):\n",
        sep = "")
    adj <- if (length(rf$confounders_used) == 0L) "unadjusted"
           else paste(rf$confounders_used, collapse = ", ")
    cat("  Model:        standard ",
        if (rf$family == "cox") "Cox" else "GLM",
        " within ", rf$signature, " (", adj, ")\n", sep = "")
    cat("  Point effect: ", format(rf$effect, digits = digits),
        "  (Wald CI [", format(rf$ci[["lower"]], digits = digits), ", ",
        format(rf$ci[["upper"]], digits = digits), "])\n", sep = "")
    cat("  Bootstrap CI (", ci_pct, " pct percentile):  [",
        format(x$refit_effect_ci[["lower"]], digits = digits), ", ",
        format(x$refit_effect_ci[["upper"]], digits = digits), "]\n",
        sep = "")
    if (isTRUE(rf$ratio_scale)) {
      ratio_lab <- sub("^log-", "", rf$effect_scale)
      cat("  ", ratio_lab, " (point):  ",
          format(exp(rf$effect), digits = digits),
          "   bootstrap CI:  [",
          format(exp(x$refit_effect_ci[["lower"]]), digits = digits), ", ",
          format(exp(x$refit_effect_ci[["upper"]]), digits = digits), "]\n",
          sep = "")
    }
    cat("  Interpretation: standard within-subgroup treatment contrast,\n")
    cat("                conditional on the discovered signature.\n")
  }

  if (!is.na(x$n_modal_match) && x$n_modal_match >= 2L) {
    cat("\nStructural stability CIs (", ci_pct,
        " pct percentile, ", x$n_modal_match,
        " iterations matching the modal cut-set):\n", sep = "")
    cat("  n_subgroup:   [",
        format(x$n_subgroup_ci[["lower"]], digits = digits), ", ",
        format(x$n_subgroup_ci[["upper"]], digits = digits), "]\n", sep = "")
    tci <- x$threshold_ci
    if (is.matrix(tci)) {
      for (k in seq_len(nrow(tci))) {
        cat("  threshold[", rownames(tci)[k], "]: [",
            format(tci[k, "lower"], digits = digits), ", ",
            format(tci[k, "upper"], digits = digits), "]\n", sep = "")
      }
    } else {
      cat("  threshold:    [",
          format(tci[["lower"]], digits = digits), ", ",
          format(tci[["upper"]], digits = digits), "]\n", sep = "")
    }
  }

  if (!is.null(x$selection_frequency)) {
    sel_long <- as.data.frame(x$selection_frequency,
                              stringsAsFactors = FALSE)
    sel_long <- sel_long[sel_long$Freq > 0L, , drop = FALSE]
    sel_long <- sel_long[order(-sel_long$Freq), , drop = FALSE]
    n_show <- min(max_select, nrow(sel_long))
    cat("\nTop ", n_show,
        " subgroup selections across found iterations:\n", sep = "")
    for (i in seq_len(n_show)) {
      cat(sprintf("  %s:  %d of %d  (%.0f pct)\n",
                  sel_long$selection[i],
                  sel_long$Freq[i], x$n_boot_found,
                  100 * sel_long$Freq[i] / x$n_boot_found))
    }
    if (nrow(sel_long) > n_show) {
      cat(sprintf("  ... %d more selections.\n",
                  nrow(sel_long) - n_show))
    }
  }

  cat("\nNote: both effect CIs (DINA and within-subgroup standard model)\n")
  cat("are conditional on the discovered signature treated as\n")
  cat("pre-specified; neither adjusts for signature selection.\n")

  invisible(x)
}
