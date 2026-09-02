# run_subgroup_sims.R -----------------------------------------------------
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

# Internal: one simulated trial -> per-subgroup (HR, UB, N) rows.
# Body order and semantics mirror the vignettes' run_one_sim() exactly:
# simulate -> flag_itt -> cut-point columns -> benchmark seed/draw ->
# subgroup subset()/fit loop with the < min_n skip rule.
.run_subgroup_sims_one <- function(ss, dgm, sim_args, cens_adjust,
                                   seed_base, rand_seed_offset,
                                   subgroups, cut_points, benchmarks,
                                   min_n, fit) {

  # 1. Draw one synthetic trial from the calibrated DGM.  `sim_args`
  #    carries exactly the design-specific argument set (resample: `n`;
  #    fixed: `baseline`), mirroring each vignette's call verbatim.
  df_s <- do.call(simulate_from_dgm, c(
    list(dgm = dgm),
    sim_args,
    list(cens_adjust = cens_adjust, seed = seed_base + ss)
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
#' @param dgm A DGM from [generate_aft_dgm_flex()]. For
#'   `baseline = "fixed"` it must carry `df_source`.
#' @param subgroups List of subgroup definitions, each a list with
#'   character scalars `id` (an expression over trial columns, exposed
#'   cut-points, `flag_itt`, and benchmark columns), `name`, and `grp`.
#' @param n_sims Number of simulated trials.
#' @param fit Per-subgroup analysis function `data -> c(estimate, upper)`;
#'   see [subgroup_cox()]. Default is an unstratified Cox fit on the
#'   simulated columns; the vignettes pass a grade-stratified formula.
#' @param baseline `"resample"` (random-X: draw `n` patients from the
#'   super-population per trial) or `"fixed"` (conditional-on-X: every
#'   trial is the `df_source` panel).
#' @param n Patients per trial. Required for `"resample"`; must be
#'   `NULL` for `"fixed"`, where N is `nrow(dgm$df_source)` by
#'   construction.
#' @param analysis_time,max_entry,cens_adjust Passed to
#'   [simulate_from_dgm()] (administrative censoring time, recruitment
#'   window, calibrated censoring shift).
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
#'   `k_treat`).
#' @param future_packages Packages loaded on workers, default
#'   `c("forestsearch", "survival")` as in the vignettes.
#' @param validate Run [validate_subgroups()] first (default `TRUE`).
#' @param verbose Print progress lines.
#'
#' @return An object of class `"subgroup_sims"`: a named list with
#'   `design`, `n_sims`, matrices `sim_hrs` / `sim_ubs` / `sim_ns`
#'   (`n_sims` x subgroups, columns named by subgroup `name`),
#'   `subgroups`, `sim_config`, `cens_adjust`, `k_treat`, `hr_true`,
#'   `created`, `r_version`, and `forestsearch_version`.
#' @seealso [summary.subgroup_sims()], [benchmark_spec()],
#'   [subgroup_cox()], [validate_subgroups()]
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
  if (!is.function(fit)) stop("`fit` must be a function of a data frame.")
  if (!is.null(benchmarks) && !inherits(benchmarks, "benchmark_spec")) {
    stop("`benchmarks` must be a benchmark_spec() or NULL.")
  }

  # Design-specific simulate_from_dgm() arguments, mirroring each
  # vignette's call exactly: resample passes `n` (no `baseline`);
  # fixed passes `baseline` (no `n` -- N is the df_source panel).
  if (baseline == "resample") {
    if (is.null(n)) {
      stop('`n` (patients per trial) is required when baseline = "resample".')
    }
    sim_args    <- list(n = n,
                        analysis_time = analysis_time, max_entry = max_entry)
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
    sim_args    <- list(baseline = baseline,
                        analysis_time = analysis_time, max_entry = max_entry)
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
    .run_subgroup_sims_one(ss, dgm, sim_args, cens_adjust,
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
      fit              = if (!is.null(fit_formula)) {
        paste(deparse(fit_formula), collapse = " ")
      } else "custom function"
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
