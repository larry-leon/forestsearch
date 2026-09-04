# run_subgroup_sims.R -----------------------------------------------------
# [delivery sentinel: p45r1-c394ff46]
# Package-level runner for the extreme-subgroups simulation study.
#
# The per-iteration computation is lifted verbatim from the committed
# vignettes (quarto/extreme_subgroups/fixed_random, commit 002f4f37):
# `run_one_sim()` + `cox_sR()` + the `%dofuture%` loop in chunks
# `sim-helpers` / `uniform-sims`.  Generalization happens only at three
# seams -- the simulate_from_dgm() argument set (resample vs fixed), the
# random-benchmark specification, and the per-subgroup fit function --
# so that with the default arguments the returned matrices are
# bit-identical to those vignettes (verified by dev/accept_phase1_bitident.R).
#
# The returned object is byte-compatible with the RDS "payload" written
# by the vignettes: saveRDS(run_subgroup_sims(...), path) produces a file
# the design-comparison memo loads unchanged.
#
# Phase 4.1 (GLM-continuous dispatch): run_subgroup_sims() also accepts a
# "glm_dgm" (generate_glm_dgm()); simulation routes to
# simulate_from_glm_dgm() through the sim_fun/sim_args seam below, and the
# default per-subgroup fit becomes subgroup_glm() -- a delegator over the
# package's native effect-estimator layer (make_effect_estimator()).  With
# a survival DGM every code path, argument default, and constructed object
# is unchanged (byte-identical results; dev/accept_phase41_*).

#' Random-benchmark subgroup specification
#'
#' Describes the size-matched random benchmark subgroups regenerated in
#' every simulated trial. With the default `nested = TRUE`, a single
#' index draw of `max(sizes)` patients is taken per trial and the smaller
#' benchmarks are its prefixes (`random15` \eqn{\subset} `random20`
#' \eqn{\subset} `random40` \eqn{\subset} `random60`), exactly as in the
#' extreme-subgroups vignettes. Membership columns are named
#' `paste0(prefix, sizes)` and take values 0/1.
#'
#' @param sizes Integer vector of benchmark sizes. Order determines the
#'   column-creation order only; the draw uses `max(sizes)`.
#' @param nested Logical. `TRUE` (default): one draw, prefix nesting --
#'   the vignette scheme. `FALSE`: an independent draw per size, in the
#'   order given (consumes more RNG; not the vignette scheme).
#' @param prefix Column-name prefix, default `"random"`.
#'
#' @return An object of class `"benchmark_spec"`.
#' @seealso [run_subgroup_sims()]
#' @export
#' @examples
#' benchmark_spec()
#' benchmark_spec(sizes = c(100, 50), prefix = "rnd")
benchmark_spec <- function(sizes = c(60L, 40L, 20L, 15L),
                           nested = TRUE,
                           prefix = "random") {
  sizes <- as.integer(sizes)
  if (length(sizes) < 1L || anyNA(sizes) || any(sizes < 1L)) {
    stop("`sizes` must be positive integers.")
  }
  if (anyDuplicated(sizes)) stop("`sizes` must be distinct.")
  structure(
    list(sizes = sizes, nested = isTRUE(nested), prefix = as.character(prefix)),
    class = c("benchmark_spec", "list")
  )
}

#' Cox per-subgroup fit function
#'
#' Returns a function `data -> c(HR, UB)` fitting `formula` by
#' [survival::coxph()] and extracting the first row of
#' `summary(fit)$conf.int`: the hazard-ratio point estimate
#' (`"exp(coef)"`) and the upper 95% bound (`"upper .95"`). Estimation
#' failure or an empty coefficient table returns `c(NA, NA)`. This
#' reproduces the vignettes' `cox_sR()` exactly when called with
#' `survival::Surv(y_sim, event_sim) ~ treat_sim + survival::strata(grade)`.
#'
#' Custom fitters may be passed to [run_subgroup_sims()] directly: any
#' function of a data frame returning a length-2 numeric
#' `(estimate, upper bound)` is accepted. Fitters run on parallel
#' workers, so they should reference non-base functions by
#' `pkg::fun` or rely on packages listed in `future_packages`.
#'
#' @section Formula environment (`lean`):
#' A formula captures the environment where it was created. In a script
#' or a knitr chunk that is the global environment, and `future`'s
#' globals machinery gives formulas special treatment: a fit whose
#' formula points at a populated environment imposes a multi-second
#' per-future cost on parallel runs (measured ~2-4 s x ~116 futures in
#' the Phase 2 timing diagnostic), independent of the serialized size.
#' With `lean = TRUE` (default) the formula is re-homed into a fresh,
#' empty environment parented on the survival namespace, cutting the
#' costly chain while keeping [survival::Surv()], [survival::strata()],
#' and base/stats symbols resolvable. Consequence: the formula must
#' reference only data columns and functions reachable from the
#' survival namespace (use `pkg::` prefixes for anything else). Set
#' `lean = FALSE` to retain the calling environment -- restoring
#' arbitrary-symbol resolution at that parallel cost. Independently of
#' `lean`, the returned closure is stripped of source references and
#' given a minimal environment holding only the formula (keep-source
#' installs otherwise attach ~250 KB of package-source baggage to every
#' worker shipment).
#'
#' @param formula A model formula for [survival::coxph()], evaluated
#'   against each subgroup's data.
#' @param lean Re-home the formula environment as described above
#'   (default `TRUE`).
#'
#' @return A function of one argument (`data`) returning
#'   `c(hr, ub)`; carries the formula in `attr(, "formula")`.
#' @export
#' @examples
#' fit <- subgroup_cox(survival::Surv(y_sim, event_sim) ~ treat_sim)
subgroup_cox <- function(formula, lean = TRUE) {
  stopifnot(inherits(formula, "formula"))
  if (isTRUE(lean)) {
    environment(formula) <- new.env(parent = asNamespace("survival"))
  }
  f <- function(data) {
    fit <- tryCatch(
      survival::coxph(formula, data = data),
      error = function(e) NULL
    )
    if (is.null(fit) || nrow(summary(fit)$conf.int) == 0) {
      return(c(NA_real_, NA_real_))
    }
    ci <- summary(fit)$conf.int[1, ]
    c(ci["exp(coef)"], ci["upper .95"])
  }
  attr(f, "formula") <- formula
  # Regardless of `lean`, ship a minimal artifact: strip source
  # references (keep-source installs -- devtools' default -- attach a
  # srcref whose srcfile drags the package source, ~250 KB measured, and
  # the execution frame's self-binding of `f` would re-drag it), and
  # rebuild the closure environment to hold exactly `formula`, parented
  # on the survival namespace.  Caller-symbol resolution in the formula
  # is governed by the FORMULA's environment (the `lean` switch above),
  # not the closure's, so this is behavior-neutral.
  f <- utils::removeSource(f)
  env_f <- new.env(parent = asNamespace("survival"))
  assign("formula", formula, envir = env_f)
  environment(f) <- env_f
  attr(f, "formula") <- formula
  f
}

#' GLM per-subgroup fit function
#'
#' Returns a function `data -> c(estimate, UB)` for non-survival outcomes
#' by delegating estimation to the package's native effect-estimator
#' layer: the closure from [make_effect_estimator()] (the same code the
#' [forestsearch()] search, bootstrap, and consistency loops use),
#' including the bit-identical unadjusted fast path
#' (`.make_lm_estimator_fast()`) under the same routing predicate
#' `forestsearch()` applies. The estimator's internal-scale
#' `estimate`/`se` are mapped to the wrapper contract as
#' `upper = estimate + qnorm(1 - (1 - level)/2) * se`, back-transformed
#' (`exp()`) for ratio measures -- the package's own display convention
#' -- so ratio-measure matrices land on the natural scale exactly as
#' [subgroup_cox()]'s `exp(coef)` / `upper .95` do. Estimation failure,
#' non-convergence, `NA` estimate/SE (e.g. a single-arm subgroup, where
#' the treatment coefficient is aliased), or a non-finite result (e.g. a
#' zero-residual-df fit) returns `c(NA, NA)`.
#'
#' @section Direction convention:
#' Estimates follow the estimator layer's harm-positive normalization:
#' with `adverse_outcome = FALSE` (higher outcome = better) the outcome
#' is negated internally so that `estimate > 0` consistently indicates
#' treatment harm, aligning with the survival convention where HR > 1 =
#' harm. Consequently the null value is always 0 (identity measures) or
#' 1 (ratio measures), whatever the raw outcome's direction.
#'
#' @section Supported settings:
#' Phase 4.5 validates `outcome_type = "continuous"` (MD),
#' `outcome_type = "binary"` (`"OR"` -- the binary default -- and
#' `"RD"`), and `outcome_type = "count"` (`"IRR"` -- the count default
#' -- and `"IRD"`; `offset.name` is required, see below). Binary
#' `"RR"` (estimator-supported; wrapper parity pending) still stops
#' with an informative error; the estimator layer
#' ([make_effect_estimator()]) already supports it and the gate lifts
#' with its parity phase -- no interface change.
#'
#' @section Serialization (`lean`):
#' The returned closure is stripped of source references and given a
#' minimal environment (the estimator closure and two scalars), and the
#' estimator closure itself is rebuilt onto a fresh environment holding
#' only its forced argument values (parented on the package namespace,
#' which serializes by reference) with its own source references
#' stripped. This severs shared srcfile environments (keep-source
#' installs) and the constructor's frame chain from the serialized
#' graph, so the fitter ships to parallel workers in the low
#' single-digit KB regardless of session state. The provenance formula
#' attached in `attr(, "formula")` (built from the column names, e.g.
#' `y_sim ~ treat_sim`) is re-homed into an empty stats-parented
#' environment when `lean = TRUE`, keeping the parallel-shipment guard
#' in [run_subgroup_sims()] meaningful.
#'
#' @param outcome.name,treat.name Column names of the outcome and the
#'   0/1 treatment indicator; defaults match [simulate_from_glm_dgm()]'s
#'   `y_sim` / `treat_sim`.
#' @param outcome_type `"continuous"`, `"binary"`, or `"count"` (all
#'   validated; count requires `offset.name`).
#' @param effect_measure `NULL` resolves to the wrapper default for the
#'   outcome type: `"MD"` (continuous), `"OR"` (binary -- the
#'   [generate_glm_dgm()] default and the calibration target; note that
#'   [make_effect_estimator()]'s own binary default is `"RD"`), or
#'   `"IRR"` (count -- all layers agree). Binary also accepts `"RD"`;
#'   count also accepts `"IRD"`; see the supported-settings section for
#'   the gated measures.
#' @param adverse_outcome Higher outcome = harm (`TRUE`, default). See
#'   the direction convention above.
#' @param robust_se,offset.name,adjust_covariates Passed to
#'   [make_effect_estimator()]. `robust_se` is consulted only where the
#'   estimator layer uses sandwich-style SEs (the count IRR path, the
#'   RR modified-Poisson fallback, and other rate estimators); the
#'   unadjusted OR path, the RD tiers, the unadjusted IRD delta chain,
#'   and MD use model-based/analytic SEs, so the flag is a no-op for
#'   those measures (asserted in the phase gates). `offset.name` is
#'   REQUIRED for `outcome_type = "count"`: it names a strictly
#'   positive exposure / follow-up time column in the analysis data
#'   (the estimator adds an internal `.log_offset` column, clamping
#'   non-positive values at `.Machine$double.eps`); for plain counts
#'   supply a unit column of 1s -- `log(1) = 0` reproduces the
#'   no-offset Poisson model exactly. Supplying `adjust_covariates`
#'   disables the unadjusted fast path (matching `forestsearch()`'s
#'   routing).
#' @param level Confidence level for the upper bound (default `0.95`).
#' @param lean Re-home the provenance formula's environment (default
#'   `TRUE`).
#'
#' @return A function of one argument (`data`) returning
#'   `c(estimate, ub)`; carries the provenance formula in
#'   `attr(, "formula")` and scale metadata in `attr(, "effect")` -- a
#'   list with `measure`, `outcome_type`, `log_scale`, `null_value`,
#'   `adverse_outcome`, `est_label`, `ub_label`, `est_thresholds`,
#'   `ub_thresholds`, `level` -- which [run_subgroup_sims()] stamps onto
#'   GLM results for scale-aware summaries and plots. A custom fitter
#'   can opt into the same behavior by attaching an `"effect"` attribute
#'   of this shape.
#' @seealso [subgroup_cox()], [run_subgroup_sims()],
#'   [make_effect_estimator()], [generate_glm_dgm()]
#' @export
#' @examples
#' \dontrun{
#' fit <- subgroup_glm()                      # y_sim ~ treat_sim, MD
#' fit(simulate_from_glm_dgm(dgm, n = 500, seed = 1))
#'
#' # Binary (Phase 4.4): OR by default, RD available
#' fit_or <- subgroup_glm(outcome_type = "binary")
#' fit_rd <- subgroup_glm(outcome_type = "binary", effect_measure = "RD")
#'
#' # Count (Phase 4.5): IRR by default, IRD available; offset required
#' fit_irr <- subgroup_glm(outcome_type = "count", offset.name = "t_exp")
#' fit_ird <- subgroup_glm(outcome_type = "count", offset.name = "t_exp",
#'                         effect_measure = "IRD")
#' }
subgroup_glm <- function(outcome.name      = "y_sim",
                         treat.name        = "treat_sim",
                         outcome_type      = "continuous",
                         effect_measure    = NULL,
                         adverse_outcome   = TRUE,
                         robust_se         = TRUE,
                         offset.name       = NULL,
                         adjust_covariates = NULL,
                         level             = 0.95,
                         lean              = TRUE) {
  outcome_type <- match.arg(outcome_type,
                            c("continuous", "binary", "count"))
  if (outcome_type == "count") {
    # Count (Phase 4.5). NULL resolves to "IRR" -- the estimator-layer
    # and generate_glm_dgm() defaults agree. Validated: IRR, IRD.
    # offset.name is required at construction: make_effect_estimator()
    # already errors without it for rate measures; this wrapper-level
    # stop pre-empts that terser message and names the unit-exposure
    # convention for plain counts.
    if (is.null(effect_measure)) effect_measure <- "IRR"
    if (!effect_measure %in% c("IRR", "IRD")) {
      stop("subgroup_glm(): effect_measure = '", effect_measure,
           "' is not available for outcome_type = \"count\" ",
           "(Phase 4.5 validates \"IRR\" and \"IRD\").", call. = FALSE)
    }
    if (is.null(offset.name)) {
      stop("subgroup_glm(): outcome_type = \"count\" requires ",
           "offset.name (an exposure / follow-up time column). For ",
           "plain counts without exposure, supply a unit column of 1s ",
           "-- log(1) = 0 reproduces the no-offset Poisson model ",
           "exactly.", call. = FALSE)
    }
  }
  if (outcome_type == "continuous") {
    if (is.null(effect_measure)) effect_measure <- "MD"
    if (!identical(effect_measure, "MD")) {
      stop("subgroup_glm(): effect_measure = '", effect_measure,
           "' is not available for outcome_type = \"continuous\" ",
           "(the continuous estimator is MD).", call. = FALSE)
    }
  } else if (outcome_type == "binary") {
    # Binary (Phase 4.4). NULL resolves to "OR" -- the
    # generate_glm_dgm() default and the calibration target -- rather
    # than make_effect_estimator()'s "RD"; the divergence is documented
    # in the Rd. Validated: OR, RD. Gated with distinct messages: RR
    # (estimator-supported, wrapper parity pending) and the rate
    # measures IRR/IRD (Phase 4.5 count-arc territory).
    if (is.null(effect_measure)) effect_measure <- "OR"
    if (effect_measure %in% c("IRR", "IRD")) {
      stop("subgroup_glm(): effect_measure = '", effect_measure,
           "' is a rate measure (requires offset.name) and arrives ",
           "with the Phase 4.5 count arc, not via ",
           "outcome_type = \"binary\".", call. = FALSE)
    }
    if (identical(effect_measure, "RR")) {
      stop("subgroup_glm(): effect_measure = \"RR\" is supported by ",
           "the estimator layer but not yet wrapper-validated ",
           "(log-binomial -> modified-Poisson chain); Phase 4.4 ",
           "validates \"OR\" and \"RD\".", call. = FALSE)
    }
    if (!effect_measure %in% c("OR", "RD")) {
      stop("subgroup_glm(): effect_measure = '", effect_measure,
           "' is not available for outcome_type = \"binary\" ",
           "(Phase 4.4 validates \"OR\" and \"RD\").", call. = FALSE)
    }
  }
  stopifnot(is.numeric(level), length(level) == 1L,
            level > 0, level < 1)

  est_fn <- make_effect_estimator(
    outcome_type      = outcome_type,
    treat.name        = treat.name,
    outcome.name      = outcome.name,
    offset.name       = offset.name,
    effect_measure    = effect_measure,
    adverse_outcome   = adverse_outcome,
    robust_se         = robust_se,
    adjust_covariates = adjust_covariates
  )
  # Unadjusted fast path -- the same predicate forestsearch() uses at its
  # construction site (the wrapper exposes no PS adjustment, so that
  # clause drops): swap in the bit-identical lm.fit closure.
  if (outcome_type == "continuous" &&
      length(.fs_adjust_terms(adjust_covariates)) == 0L) {
    est_fn <- .make_lm_estimator_fast(
      treat.name      = treat.name,
      outcome.name    = outcome.name,
      adverse_outcome = adverse_outcome
    )
  }

  # Phase 4.1r2 serialization hygiene: rebuild the estimator closure onto
  # a fresh minimal environment (its forced argument values only; parent
  # is the package namespace, which serializes by reference) and strip
  # source references. This severs every session-mutable carrier from the
  # shipped fitter's serialized graph -- shared srcfile environments left
  # by keep-source installs, and the constructor's execution-frame chain
  # -- fixing the G7-F pathology where a freshly built fitter serialized
  # at ~3.3 MB inside a loaded session (34 KB clean, 1.4 KB with both
  # carriers cut; see dev/diag_phase41_serialization.R). After this block
  # the fitter's ordinary-environment content is exactly: est_fn's value
  # bindings, z, log_scale, and an empty formula environment.
  est_env <- new.env(parent = asNamespace("forestsearch"))
  for (nm in setdiff(ls(environment(est_fn), all.names = TRUE), "...")) {
    assign(nm, get(nm, envir = environment(est_fn)), envir = est_env)
  }
  environment(est_fn) <- est_env
  est_fn <- utils::removeSource(est_fn)

  z         <- stats::qnorm(1 - (1 - level) / 2)
  log_scale <- is_log_scale(effect_measure_to_label(effect_measure))

  f <- function(data) {
    r <- tryCatch(est_fn(data), error = function(e) NULL)
    if (is.null(r) || !isTRUE(r$converged) ||
        is.na(r$estimate) || is.na(r$se)) {
      return(c(NA_real_, NA_real_))
    }
    est <- r$estimate
    up  <- est + z * r$se
    out <- if (log_scale) c(exp(est), exp(up)) else c(est, up)
    if (!all(is.finite(out))) return(c(NA_real_, NA_real_))
    out
  }

  fml <- .build_adjusted_formula(outcome.name, treat.name,
                                 adjust_covariates)
  if (isTRUE(lean)) {
    environment(fml) <- new.env(parent = asNamespace("stats"))
  }
  f <- utils::removeSource(f)
  env_f <- new.env(parent = asNamespace("forestsearch"))
  assign("est_fn", est_fn, envir = env_f)
  assign("z", z, envir = env_f)
  assign("log_scale", log_scale, envir = env_f)
  environment(f) <- env_f
  attr(f, "formula") <- fml
  attr(f, "effect") <- list(
    measure         = effect_measure,
    outcome_type    = outcome_type,
    log_scale       = log_scale,
    null_value      = if (log_scale) 1 else 0,
    adverse_outcome = isTRUE(adverse_outcome),
    est_label       = effect_measure,
    ub_label        = paste0("UB(", effect_measure, ")"),
    # Scale-aware stamps (Phase 4.4): ratio measures inherit the HR
    # legacy tails c(0.5, 1) / c(2, 3) -- summary() then reuses the
    # exact legacy header strings and the high-risk panel works out of
    # the box. Identity measures keep the 4.1 stamps (no universal
    # tails exist); the identity branch evaluates to the pre-4.4
    # values, so continuous output is unchanged.
    est_thresholds  = if (log_scale) c(0.5, 1) else c(NA_real_, 0),
    ub_thresholds   = if (log_scale) c(2, 3) else c(NA_real_, NA_real_),
    level           = level
  )
  f
}

#' Validate subgroup definitions against a data frame
#'
#' Evaluates every subgroup `id` expression against `data` (augmented
#' with `flag_itt = 1L`, the scalar `cutpoints` as columns, and 0-valued
#' benchmark membership columns) and fails with a named, collected error
#' if any expression errors, is non-logical, or has the wrong length.
#' Run automatically at the top of [run_subgroup_sims()] so that a typo
#' in a subgroup definition stops the study in milliseconds instead of
#' surfacing as a silent all-`NA` column after the simulation loop.
#'
#' @param subgroups List of subgroup definitions; each element a list
#'   with character scalars `id` (an R expression over columns), `name`,
#'   and `grp`.
#' @param data A data frame carrying the covariates the `id` expressions
#'   reference (e.g. `dgm$df_source` or `dgm$df_super`).
#' @param cutpoints Named list of scalars exposed as columns (e.g.
#'   `list(age_med = 53)`), matching the runner's `cutpoints` argument.
#' @param benchmarks A [benchmark_spec()] (or `NULL`) whose membership
#'   columns are stubbed as `0L` during validation.
#'
#' @return Invisibly `TRUE`; errors otherwise.
#' @export
validate_subgroups <- function(subgroups, data,
                               cutpoints = list(),
                               benchmarks = benchmark_spec()) {
  if (!is.list(subgroups) || length(subgroups) == 0L) {
    stop("`subgroups` must be a non-empty list.")
  }
  for (fld in c("id", "name", "grp")) {
    bad <- vapply(subgroups, function(s) {
      !is.list(s) || !is.character(s[[fld]]) || length(s[[fld]]) != 1L ||
        is.na(s[[fld]])
    }, logical(1L))
    if (any(bad)) {
      stop("subgroup element(s) ", paste(which(bad), collapse = ", "),
           " lack a character scalar `", fld, "`.")
    }
  }
  nms <- vapply(subgroups, `[[`, character(1L), "name")
  if (anyDuplicated(nms)) {
    stop("duplicate subgroup names: ",
         paste(unique(nms[duplicated(nms)]), collapse = ", "))
  }

  probe <- as.data.frame(data)
  probe$flag_itt <- 1L
  for (nm in names(cutpoints)) probe[[nm]] <- cutpoints[[nm]]
  if (!is.null(benchmarks)) {
    for (s in benchmarks$sizes) probe[[paste0(benchmarks$prefix, s)]] <- 0L
  }

  problems <- character(0)
  for (sg in subgroups) {
    v <- tryCatch(
      eval(parse(text = sg$id), envir = probe, enclos = parent.frame()),
      error = function(e) e
    )
    if (inherits(v, "error")) {
      problems <- c(problems, sprintf("'%s' (%s): %s",
                                      sg$name, sg$id, conditionMessage(v)))
    } else if (!is.logical(v) || length(v) != nrow(probe)) {
      problems <- c(problems, sprintf(
        "'%s' (%s): evaluates to %s of length %d, expected logical of length %d",
        sg$name, sg$id, typeof(v), length(v), nrow(probe)))
    }
  }
  if (length(problems)) {
    stop("invalid subgroup definition(s):\n  ",
         paste(problems, collapse = "\n  "))
  }
  invisible(TRUE)
}

# Internal: one simulated trial -> per-subgroup (estimate, UB, N) rows.
# Body order and semantics mirror the vignettes' run_one_sim() exactly:
# simulate -> flag_itt -> cut-point columns -> benchmark seed/draw ->
# subgroup subset()/fit loop with the < min_n skip rule.  The simulate
# seam is outcome-agnostic: the runner supplies `sim_fun`
# (simulate_from_dgm or simulate_from_glm_dgm) and a COMPLETE
# design-specific `sim_args` (survival resample: n/analysis_time/
# max_entry/cens_adjust; survival fixed: baseline/analysis_time/
# max_entry/cens_adjust; glm: n); only `seed` is appended here.  All
# arguments are named, so on the survival path the composed call is
# identical to the pre-seam form.
.run_subgroup_sims_one <- function(ss, dgm, sim_fun, sim_args,
                                   seed_base, rand_seed_offset,
                                   subgroups, cut_points, benchmarks,
                                   min_n, fit) {

  # 1. Draw one synthetic trial from the DGM.
  df_s <- do.call(sim_fun, c(
    list(dgm = dgm),
    sim_args,
    list(seed = seed_base + ss)
  ))

  # 2. Expose flag_itt and scalar cut-points as columns
  df_s$flag_itt <- 1L
  for (nm in names(cut_points)) df_s[[nm]] <- cut_points[[nm]]

  # 3. Regenerate random-benchmark subgroups for this trial
  if (!is.null(benchmarks)) {
    set.seed(seed_base + ss + rand_seed_offset)
    n_s <- nrow(df_s)
    if (benchmarks$nested) {
      r_idx <- sample.int(n_s, min(max(benchmarks$sizes), n_s),
                          replace = FALSE)
      for (s in benchmarks$sizes) {
        df_s[[paste0(benchmarks$prefix, s)]] <-
          as.integer(seq_len(n_s) %in% r_idx[seq_len(min(s, length(r_idx)))])
      }
    } else {
      for (s in benchmarks$sizes) {
        idx_s <- sample.int(n_s, min(s, n_s), replace = FALSE)
        df_s[[paste0(benchmarks$prefix, s)]] <-
          as.integer(seq_len(n_s) %in% idx_s)
      }
    }
  }

  # 4. Fit in every pre-defined subgroup
  sg_ids <- vapply(subgroups, `[[`, character(1L), "id")
  n_sg   <- length(subgroups)
  hr_row <- rep(NA_real_, n_sg)
  ub_row <- rep(NA_real_, n_sg)
  n_row  <- rep(NA_real_, n_sg)

  for (gg in seq_len(n_sg)) {
    df_sg <- tryCatch(
      subset(df_s, eval(parse(text = sg_ids[[gg]]))),
      error = function(e) NULL
    )
    if (is.null(df_sg) || nrow(df_sg) < min_n) next
    n_row[gg]  <- nrow(df_sg)
    r          <- fit(df_sg)
    hr_row[gg] <- r[1]
    ub_row[gg] <- r[2]
  }

  list(hr = hr_row, ub = ub_row, n = n_row)
}

#' Run the extreme-subgroups simulation study
#'
#' Draws `n_sims` synthetic trials from a calibrated DGM via
#' [simulate_from_dgm()], fits `fit` in every pre-defined subgroup of
#' every trial, and returns the raw result matrices plus full
#' provenance. With the default `fit`, `benchmarks`, seeds, and skip
#' rule, the matrices are bit-identical to those produced by the
#' `extreme_subgroups` vignettes for the matching `baseline`, and the
#' returned object is byte-compatible with their RDS payload:
#' `saveRDS(sims, path)` yields a file the design-comparison memo loads
#' unchanged.
#'
#' Per-iteration seeding is explicit (trial `ss` uses
#' `seed_base + ss` for the DGM draw and
#' `seed_base + ss + rand_seed_offset` for benchmark membership), so
#' results are independent of the parallel schedule: any `workers`
#' setting, including sequential, gives identical matrices.
#'
#' GLM dispatch: when `dgm` inherits `"glm_dgm"` ([generate_glm_dgm()]),
#' trials are drawn by [simulate_from_glm_dgm()] (`n` is required; the
#' survival-only `analysis_time` / `max_entry` / `cens_adjust` must not
#' be supplied, and `baseline = "fixed"` is unavailable), the default
#' `fit` becomes [subgroup_glm()] constructed from the DGM's own
#' `outcome_type` and `effect_measure` (binary DGMs also pass their
#' `adverse_outcome`; a continuous DGM resolves to the same fitter as
#' before), and the result additionally carries
#' the fitter's `effect` scale metadata with `design = "glm"`. The
#' `sim_hrs` / `sim_ubs` matrices then hold that fitter's
#' `(estimate, upper bound)` pairs -- the field names are retained for
#' structural compatibility, exactly as `hr.threshold` serves
#' generic-threshold duty for GLM outcomes in [forestsearch()].
#'
#' @param dgm A DGM from [generate_aft_dgm_flex()] (survival) or
#'   [generate_glm_dgm()] (GLM). For `baseline = "fixed"` (survival
#'   only) it must carry `df_source`.
#' @param subgroups List of subgroup definitions, each a list with
#'   character scalars `id` (an expression over trial columns, exposed
#'   cut-points, `flag_itt`, and benchmark columns), `name`, and `grp`.
#' @param n_sims Number of simulated trials.
#' @param fit Per-subgroup analysis function `data -> c(estimate, upper)`;
#'   see [subgroup_cox()]. Default is an unstratified Cox fit on the
#'   simulated columns (the vignettes pass a grade-stratified formula);
#'   for a `"glm_dgm"` the default is [subgroup_glm()] constructed from
#'   the DGM's `outcome_type` / `effect_measure` (and, for binary, its
#'   `adverse_outcome`).
#' @param baseline `"resample"` (random-X: draw `n` patients from the
#'   super-population per trial) or `"fixed"` (conditional-on-X: every
#'   trial is the `df_source` panel). Survival only; a `"glm_dgm"`
#'   always resamples.
#' @param n Patients per trial. Required for `"resample"` and for a
#'   `"glm_dgm"`; must be `NULL` for `"fixed"`, where N is
#'   `nrow(dgm$df_source)` by construction.
#' @param analysis_time,max_entry,cens_adjust Passed to
#'   [simulate_from_dgm()] (administrative censoring time, recruitment
#'   window, calibrated censoring shift). Survival-only: supplying any
#'   of them with a `"glm_dgm"` is an error.
#' @param cutpoints Named list of scalars exposed as columns of each
#'   simulated trial so `id` expressions like `"age <= age_med"`
#'   resolve.
#' @param benchmarks A [benchmark_spec()], or `NULL` to skip benchmark
#'   generation entirely (no RNG is consumed for it).
#' @param min_n Subgroups with fewer rows are skipped for that trial
#'   (row stays `NA`), default `5L` as in the vignettes.
#' @param workers Parallel workers. `NULL` (default) reuses a
#'   caller-established non-sequential [future::plan()] when one exists
#'   -- no worker-pool spawn or teardown, the mode the vignettes use
#'   after their `parallel-setup` chunk -- and otherwise creates a
#'   multisession plan on `ceiling(0.90 * availableCores())`, restored
#'   on exit. A numeric value forces a plan for this call (restored on
#'   exit): a fraction in (0, 1) of `availableCores()`, an integer
#'   count, or `1` for sequential. Results are identical under every
#'   setting (seeding is per-iteration); only wall time differs.
#'   `sim_config` records `workers`, `plan_reused`, and `t_plan_secs`.
#' @param seed_base,rand_seed_offset Seed scheme; defaults (`0L`,
#'   `1e6L`) reproduce the vignettes.
#' @param hr_true,k_treat Optional provenance values stored verbatim in
#'   the result (the vignettes store `0.70` and the calibrated
#'   `k_treat`). `hr_true` is the true effect on the fitter's estimate
#'   scale: a hazard ratio for Cox fits, a mean difference for
#'   [subgroup_glm()] fits.
#' @param future_packages Packages loaded on workers, default
#'   `c("forestsearch", "survival")` as in the vignettes.
#' @param validate Run [validate_subgroups()] first (default `TRUE`).
#' @param verbose Print progress lines.
#'
#' @return An object of class `"subgroup_sims"`: a named list with
#'   `design`, `n_sims`, matrices `sim_hrs` / `sim_ubs` / `sim_ns`
#'   (`n_sims` x subgroups, columns named by subgroup `name`),
#'   `subgroups`, `sim_config`, `cens_adjust`, `k_treat`, `hr_true`,
#'   `created`, `r_version`, and `forestsearch_version`. For a
#'   `"glm_dgm"`: `design = "glm"`, `cens_adjust = NULL`, and an
#'   additional `effect` field holding the fitter's scale metadata (see
#'   [subgroup_glm()]); survival results are constructed exactly as
#'   before, with no new fields.
#' @seealso [summary.subgroup_sims()], [benchmark_spec()],
#'   [subgroup_cox()], [subgroup_glm()], [validate_subgroups()]
#' @export
#' @importFrom foreach foreach
#' @importFrom doFuture %dofuture%
#' @examples
#' \dontrun{
#' sims <- run_subgroup_sims(
#'   dgm = dgm_uniform, subgroups = subgroups, n_sims = 1000,
#'   fit = subgroup_cox(survival::Surv(y_sim, event_sim) ~ treat_sim +
#'                        survival::strata(grade)),
#'   baseline = "fixed",
#'   analysis_time = 84, max_entry = 24,
#'   cens_adjust = cal_uniform$cens_adjust,
#'   cutpoints = list(age_med = median(gbsg$age)),
#'   hr_true = 0.70, k_treat = k_treat_uniform
#' )
#' saveRDS(sims, "results/extreme_sims_fixed_1000_payload.rds")
#' }
run_subgroup_sims <- function(dgm, subgroups, n_sims,
                              fit = subgroup_cox(
                                survival::Surv(y_sim, event_sim) ~ treat_sim),
                              baseline = c("resample", "fixed"),
                              n = NULL,
                              analysis_time, max_entry, cens_adjust,
                              cutpoints = list(),
                              benchmarks = benchmark_spec(),
                              min_n = 5L,
                              workers = NULL,
                              seed_base = 0L,
                              rand_seed_offset = 1e6L,
                              hr_true = NULL,
                              k_treat = NULL,
                              future_packages = c("forestsearch", "survival"),
                              validate = TRUE,
                              verbose = FALSE) {

  baseline <- match.arg(baseline)
  n_sims   <- as.integer(n_sims)
  stopifnot(length(n_sims) == 1L, !is.na(n_sims), n_sims >= 1L)
  is_glm <- inherits(dgm, "glm_dgm")
  # Resolve the default fit BEFORE it is forced: for a glm_dgm the
  # survival default (Surv(y_sim, event_sim) ~ treat_sim) would be
  # wrong, and because default arguments are promises, assigning here
  # means that default is never evaluated on the GLM path.
  # DGM-aware construction (Phase 4.4): carry the DGM's own
  # outcome_type / effect_measure so a binary DGM analyzes on its own
  # scale instead of the continuous/MD defaults. A continuous DGM
  # resolves to exactly the pre-4.4 values ("continuous", "MD",
  # adverse_outcome = TRUE), so the default path is behavior-identical
  # there. adverse_outcome is read from the DGM only for binary -- the
  # one type where generate_glm_dgm() defines it.
  if (is_glm && missing(fit)) fit <- subgroup_glm(
    outcome_type    = dgm$outcome_type,
    effect_measure  = dgm$effect_measure,
    # NOTE the nesting asymmetry (do not "normalize" these two):
    # adverse_outcome is a TOP-LEVEL glm_dgm field; offset_var lives in
    # model_params. offset_var is NULL for continuous/binary DGMs --
    # identical to the pre-4.5 default argument, so those default paths
    # are byte-identical. A count DGM built without offset_var hits the
    # count offset requirement at construction (unit-exposure
    # convention; see subgroup_glm()'s Rd).
    offset.name     = dgm$model_params$offset_var,
    adverse_outcome = if (identical(dgm$outcome_type, "binary"))
                        isTRUE(dgm$adverse_outcome) else TRUE
  )
  if (!is.function(fit)) stop("`fit` must be a function of a data frame.")
  if (!is.null(benchmarks) && !inherits(benchmarks, "benchmark_spec")) {
    stop("`benchmarks` must be a benchmark_spec() or NULL.")
  }

  # Design-specific simulate function and arguments.  The worker composes
  # do.call(sim_fun, c(list(dgm = dgm), sim_args, list(seed = seed_base + ss))),
  # so sim_args carries everything but dgm and seed.  The survival
  # branches mirror each vignette's simulate_from_dgm() call exactly
  # (resample passes `n`, fixed passes `baseline`; cens_adjust now rides
  # inside sim_args -- named arguments, so the composed call is identical
  # to the pre-seam form).  The GLM branch routes to
  # simulate_from_glm_dgm(n).
  if (is_glm) {
    if (!missing(analysis_time) || !missing(max_entry) ||
        !missing(cens_adjust)) {
      stop("`analysis_time`, `max_entry`, and `cens_adjust` are ",
           'survival-only arguments: a "glm_dgm" has no censoring or ',
           "entry process.")
    }
    if (baseline == "fixed") {
      stop('baseline = "fixed" is survival-only: a "glm_dgm" carries no ',
           'df_source panel. Omit `baseline` (it defaults to "resample").')
    }
    if (is.null(n)) {
      stop('`n` (patients per trial) is required for a "glm_dgm". ',
           "(simulate_from_glm_dgm(n = NULL) means the full ",
           "super-population -- an evaluation-frame convenience, not a ",
           "trial size.)")
    }
    sim_fun     <- simulate_from_glm_dgm
    sim_args    <- list(n = n)
    n_per_trial <- as.integer(n)
    eff <- attr(fit, "effect")
    if (!is.null(eff) && !is.null(dgm$effect_measure) &&
        !identical(eff$measure, dgm$effect_measure)) {
      warning(sprintf(
        paste0("`fit` measures '%s' but the DGM was built for '%s'; ",
               "check the fitter/DGM pairing."),
        eff$measure, dgm$effect_measure), call. = FALSE)
    }
  } else if (baseline == "resample") {
    if (is.null(n)) {
      stop('`n` (patients per trial) is required when baseline = "resample".')
    }
    sim_fun     <- simulate_from_dgm
    sim_args    <- list(n = n,
                        analysis_time = analysis_time, max_entry = max_entry,
                        cens_adjust = cens_adjust)
    n_per_trial <- as.integer(n)
  } else {
    if (!is.null(n)) {
      stop('`n` must be NULL when baseline = "fixed": per-trial N is ',
           "nrow(dgm$df_source) by construction.")
    }
    if (is.null(dgm$df_source)) {
      stop("dgm$df_source is missing: this DGM does not store the ",
           'fixed-baseline frame required for baseline = "fixed".')
    }
    sim_fun     <- simulate_from_dgm
    sim_args    <- list(baseline = baseline,
                        analysis_time = analysis_time, max_entry = max_entry,
                        cens_adjust = cens_adjust)
    n_per_trial <- nrow(dgm$df_source)
  }

  if (isTRUE(validate)) {
    probe <- if (!is.null(dgm$df_source)) dgm$df_source else dgm$df_super
    if (is.null(probe)) {
      warning("no df_source/df_super on the DGM; skipping subgroup validation.")
    } else {
      validate_subgroups(subgroups, probe, cutpoints, benchmarks)
    }
  }

  # Parallel plan.  workers = NULL: reuse a caller-established
  # non-sequential plan when one exists (no per-call worker-pool spawn/
  # teardown -- the multisession setup cost dominated the Phase 1
  # acceptance timings); otherwise, and for any numeric `workers`, set a
  # plan for this call and restore the previous one on exit.
  t0_plan     <- Sys.time()
  plan_reused <- FALSE
  if (is.null(workers)) {
    if (future::nbrOfWorkers() > 1L) {
      # A caller-established parallel plan is in place -- reuse it.
      # (nbrOfWorkers() rather than class-sniffing: future versions
      # differ in how plan() objects are classed, e.g. "FutureStrategy".)
      plan_reused <- TRUE
    } else {
      old_plan <- future::plan()
      on.exit(future::plan(old_plan), add = TRUE)
      future::plan(future::multisession,
                   workers = max(1L, as.integer(
                     ceiling(0.90 * future::availableCores()))))
    }
  } else {
    n_set <- if (workers < 1) {
      max(1L, as.integer(ceiling(workers * future::availableCores())))
    } else {
      as.integer(workers)
    }
    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)
    if (n_set <= 1L) {
      future::plan(future::sequential)
    } else {
      future::plan(future::multisession, workers = n_set)
    }
  }
  n_workers <- as.integer(future::nbrOfWorkers())
  t_plan    <- as.numeric(difftime(Sys.time(), t0_plan, units = "secs"))

  # Guard against the formula-environment trap diagnosed in Phase 2:
  # a fit whose formula is bound to the (populated) global environment,
  # or that serializes large, imposes a multi-second per-future cost.
  # subgroup_cox(lean = TRUE) -- the default -- avoids both conditions.
  if (n_workers > 1L) {
    fml <- attr(fit, "formula")
    fit_bytes <- tryCatch(length(serialize(fit, NULL)),
                          error = function(e) NA_integer_)
    if ((!is.null(fml) && identical(environment(fml), globalenv())) ||
        isTRUE(fit_bytes > 512 * 1024)) {
      warning(sprintf(
        paste0("`fit` will be shipped to every parallel worker and %s; ",
               "this slowed the Phase 2 renders ~8x. Build it with ",
               "subgroup_cox() (lean = TRUE) or give the formula a ",
               "minimal environment -- see ?subgroup_cox."),
        if (!is.null(fml) && identical(environment(fml), globalenv())) {
          "its formula is bound to the global environment"
        } else {
          sprintf("serializes to %.1f MB", fit_bytes / 1024^2)
        }), call. = FALSE)
    }
  }
  if (verbose) {
    message(sprintf(
      "run_subgroup_sims: %s design, %d trials, %d worker(s) [%s, %.1f s]",
      baseline, n_sims, n_workers,
      if (plan_reused) "plan reused" else "plan created", t_plan))
  }

  # Scoped RNG-misuse silencing.  Seeding is explicit per iteration
  # (seed_base + ss for the DGM draw; + rand_seed_offset for benchmark
  # membership), so future's misuse check is a false alarm here -- but
  # the fix must NOT be seed = TRUE: future would install an
  # L'Ecuyer-CMRG .Random.seed per future, switching the active RNG
  # kind, and set.seed(x, kind = NULL) seeds the CURRENT kind.
  # Demonstrated: set.seed(1000001); sample.int(686, 5) gives
  # 196 191 606 419 507 under Mersenne-Twister but 417 682 50 132 4
  # under L'Ecuyer -- every benchmark draw and DGM simulation inside
  # workers would silently change.  Ignoring the check changes no RNG
  # state at all; the option is saved and restored on exit.
  old_opt <- options(future.rng.onMisuse = "ignore")
  on.exit(options(old_opt), add = TRUE)

  # Simulation loop -- structure mirrors the vignettes' uniform-sims
  # chunk; per-iteration seeding is explicit, so seed = NULL here.
  t0 <- Sys.time()
  ss <- NULL  # foreach iteration variable (R CMD check appeasement)
  sim_results <- foreach::foreach(
    ss = seq_len(n_sims),
    .options.future = list(
      seed     = NULL,
      packages = future_packages
    )
  ) %dofuture% {
    .run_subgroup_sims_one(ss, dgm, sim_fun, sim_args,
                           seed_base, rand_seed_offset,
                           subgroups, cutpoints, benchmarks,
                           min_n, fit)
  }
  t_sims <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  sg_names <- vapply(subgroups, `[[`, character(1L), "name")
  sim_hrs  <- do.call(rbind, lapply(sim_results, `[[`, "hr"))
  sim_ubs  <- do.call(rbind, lapply(sim_results, `[[`, "ub"))
  sim_ns   <- do.call(rbind, lapply(sim_results, `[[`, "n"))
  colnames(sim_hrs) <- sg_names
  colnames(sim_ubs) <- sg_names
  colnames(sim_ns)  <- sg_names

  if (verbose) {
    message(sprintf("run_subgroup_sims: completed in %.1f s", t_sims))
  }

  fit_formula <- attr(fit, "formula")
  fit_desc <- if (!is.null(fit_formula)) {
    paste(deparse(fit_formula), collapse = " ")
  } else "custom function"
  if (is_glm) {
    out <- list(
      design      = "glm",
      n_sims      = n_sims,
      sim_hrs     = sim_hrs,
      sim_ubs     = sim_ubs,
      sim_ns      = sim_ns,
      subgroups   = subgroups,
      sim_config  = list(
        baseline         = "glm",
        n_per_trial      = n_per_trial,
        seed_base        = seed_base,
        rand_seed_offset = rand_seed_offset,
        min_n            = min_n,
        benchmarks       = benchmarks,
        workers          = n_workers,
        plan_reused      = plan_reused,
        t_plan_secs      = t_plan,
        t_sims_secs      = t_sims,
        fit              = fit_desc
      ),
      cens_adjust = NULL,
      k_treat     = k_treat,
      hr_true     = hr_true,
      effect      = attr(fit, "effect"),
      created     = Sys.time(),
      r_version   = R.version.string,
      forestsearch_version = tryCatch(
        as.character(utils::packageVersion("forestsearch")),
        error = function(e) NA_character_)
    )
  } else {
    out <- list(
      design      = baseline,
      n_sims      = n_sims,
      sim_hrs     = sim_hrs,
      sim_ubs     = sim_ubs,
      sim_ns      = sim_ns,
      subgroups   = subgroups,
      sim_config  = list(
        baseline         = baseline,
        n_per_trial      = n_per_trial,
        analysis_time    = analysis_time,
        max_entry        = max_entry,
        seed_base        = seed_base,
        rand_seed_offset = rand_seed_offset,
        min_n            = min_n,
        benchmarks       = benchmarks,
        workers          = n_workers,
        plan_reused      = plan_reused,
        t_plan_secs      = t_plan,
        t_sims_secs      = t_sims,
        fit              = fit_desc
      ),
      cens_adjust = cens_adjust,
      k_treat     = k_treat,
      hr_true     = hr_true,
      created     = Sys.time(),
      r_version   = R.version.string,
      forestsearch_version = tryCatch(
        as.character(utils::packageVersion("forestsearch")),
        error = function(e) NA_character_)
    )
  }
  class(out) <- c("subgroup_sims", "list")
  out
}

#' @export
print.subgroup_sims <- function(x, ...) {
  cat("<subgroup_sims>\n")
  cat(sprintf("  design      : %s\n", x$design))
  cat(sprintf("  n_sims      : %s\n", format(x$n_sims, big.mark = ",")))
  cat(sprintf("  subgroups   : %d (per-trial N = %s)\n",
              ncol(x$sim_hrs), x$sim_config$n_per_trial))
  cat(sprintf("  fit         : %s\n", x$sim_config$fit))
  cat(sprintf("  created     : %s | forestsearch %s\n",
              format(x$created, "%Y-%m-%d %H:%M"), x$forestsearch_version))
  invisible(x)
}
