#' Collect Bootstrap Foreach Results
#'
#' Internal helper that converts the list returned by a foreach loop with
#' \code{.errorhandling = "pass"} into a single data.table, separating
#' successful iterations (data.table rows) from failed iterations
#' (condition objects).
#'
#' Implements the convention used in the package's Quarto calibration
#' suite: \code{vapply(..., inherits, ..., "error")} for error
#' separation, \code{data.table::rbindlist(..., fill = TRUE)} for
#' combine, and \code{stop()} with the first error message when every
#' iteration has failed. Partial failures generate a warning with the
#' failure count and first error message; the successful rows are
#' returned.
#'
#' @param results List of length \code{nb_boots}, as returned by a
#'   foreach loop with \code{.errorhandling = "pass"} and no
#'   \code{.combine} (so each element is either a data.table or a
#'   condition object).
#' @param nb_boots Integer. Number of bootstrap iterations attempted.
#'
#' @return A data.table containing the rbinded successful iterations.
#'   When all iterations have failed, this function does not return -
#'   it raises an error.
#'
#' @keywords internal
#' @noRd
.collect_bootstrap_results <- function(results, nb_boots) {
  is_err <- vapply(results, inherits, logical(1), what = "error")
  n_err  <- sum(is_err)

  if (n_err == nb_boots) {
    first_err <- results[[which(is_err)[1L]]]
    stop(
      "All ", nb_boots, " bootstrap iterations failed. First error:\n  ",
      conditionMessage(first_err),
      call. = FALSE
    )
  }

  if (n_err > 0L) {
    first_err <- results[[which(is_err)[1L]]]
    warning(
      sprintf(
        paste0("%d of %d bootstrap iterations failed (%.1f%%); these ",
               "are dropped from the result.\n  First error: %s"),
        n_err, nb_boots, 100 * n_err / nb_boots,
        conditionMessage(first_err)
      ),
      call. = FALSE
    )
  }

  data.table::rbindlist(results[!is_err], fill = TRUE)
}


#' Bootstrap Results for ForestSearch with Bias Correction
#'
#' Runs bootstrap analysis for ForestSearch, fitting Cox models and computing
#' bias-corrected estimates and valid CIs (see vignette for references)
#'
#' @param fs.est List. ForestSearch results object from \code{\link{forestsearch}}.
#'   Must contain:
#'   \itemize{
#'     \item \code{df.est}: Data frame with analysis data including \code{treat.recommend}
#'     \item \code{confounders.candidate}: Character vector of confounder names
#'     \item \code{args_call_all}: List of original forestsearch call arguments
#'   }
#' @param df_boot_analysis Data frame. Bootstrap analysis data with same structure
#'   as \code{fs.est$df.est}. Must contain columns for outcome, event, treatment,
#'   and the \code{treat.recommend} flag.
#' @param cox.formula.boot Formula. Cox model formula for bootstrap, typically
#'   created by \code{\link{build_cox_formula}}. Should be of form
#'   \code{Surv(outcome, event) ~ treatment}.
#' @param nb_boots Integer. Number of bootstrap samples to generate (e.g., 500-1000).
#'   More iterations provide better bias correction but increase computation time.
#' @param show_three Logical. If \code{TRUE}, prints detailed progress for the
#'   first three bootstrap iterations for debugging purposes. Default: \code{FALSE}.
#' @param H_obs Numeric. Observed log hazard ratio for subgroup H (harm/questionable group,
#'   \code{treat.recommend == 0}) from original sample. Used as reference for
#'   bias correction.
#' @param Hc_obs Numeric. Observed log hazard ratio for subgroup H^c (complement/recommend,
#'   \code{treat.recommend == 1}) from original sample. Used as reference for
#'   bias correction.
#' @param seed Integer. Random seed for reproducibility. Default 8316951L.
#'   Only used when \code{boot_index_mat} is \code{NULL}; when an
#'   explicit index matrix is supplied, it is the source of bootstrap
#'   row selection and \code{seed} has no effect on bootstrap indices.
#' @param estimator_fn Closure or \code{NULL}. Effect-estimator closure from
#'   \code{\link{make_effect_estimator}} for GLM outcomes.  When \code{NULL}
#'   (default), the Cox model via \code{cox.formula.boot} is used.
#' @param boot_index_mat Integer matrix or \code{NULL}. Optional
#'   pre-generated \eqn{B \times N} matrix of bootstrap row indices. When
#'   supplied, row \eqn{b} is used as the resampling vector for iteration
#'   \eqn{b} (\code{df_boot_analysis[boot_index_mat[b, ], ]}); the
#'   per-iteration call to \code{sample.int()} is bypassed. This is the
#'   path taken by \code{\link{forestsearch_bootstrap_dofuture}}, which
#'   generates the matrix once on the main process so that the same
#'   indices drive both the worker's bootstrap data construction and the
#'   \code{Ystar_mat} count matrix. When \code{NULL} (default), each
#'   iteration draws its own indices via \code{sample.int()} - retained
#'   for backward compatibility with direct callers.
#'
#' @return Data.table with one row per bootstrap iteration and columns:
#'   \describe{
#'     \item{boot_id}{Integer. Bootstrap iteration number (1 to \code{nb_boots})}
#'     \item{H_biasadj_1}{Bias-corrected estimate for H using method 1:
#'       \code{H_obs - (Hstar_star - Hstar_obs)}}
#'     \item{H_biasadj_2}{Bias-corrected estimate for H using method 2:
#'       \code{2*H_obs - (H_star + Hstar_star - Hstar_obs)}}
#'     \item{Hc_biasadj_1}{Bias-corrected estimate for H^c using method 1}
#'     \item{Hc_biasadj_2}{Bias-corrected estimate for H^c using method 2}
#'     \item{max_sg_est}{Numeric. Maximum subgroup hazard ratio found}
#'     \item{L}{Integer. Number of candidate factors evaluated}
#'     \item{max_count}{Integer. Maximum number of factor combinations}
#'     \item{events_H_0}{Integer. Number of events in control arm of original subgroup H on bootstrap sample}
#'     \item{events_H_1}{Integer. Number of events in treatment arm of original subgroup H on bootstrap sample}
#'     \item{events_Hc_0}{Integer. Number of events in control arm of original subgroup H^c on bootstrap sample}
#'     \item{events_Hc_1}{Integer. Number of events in treatment arm of original subgroup H^c on bootstrap sample}
#'     \item{events_Hstar_0}{Integer. Number of events in control arm of new subgroup H* on original data}
#'     \item{events_Hstar_1}{Integer. Number of events in treatment arm of new subgroup H* on original data}
#'     \item{events_Hcstar_0}{Integer. Number of events in control arm of new subgroup H^c* on original data}
#'     \item{events_Hcstar_1}{Integer. Number of events in treatment arm of new subgroup H^c* on original data}
#'     \item{tmins_search}{Numeric. Minutes spent on subgroup search in this iteration}
#'     \item{tmins_iteration}{Numeric. Total minutes for this bootstrap iteration}
#'     \item{Pcons}{Numeric. Consistency p-value for top subgroup}
#'     \item{any_found}{Integer 0/1. Identification flag; equivalent to
#'       \code{as.integer(!is.na(Pcons))}.  Added in v0.2.0 to mirror the
#'       \code{fold_summary$any_found} column produced by
#'       \code{\link{forestsearch_tenfold}}, so cross-API diagnostics
#'       (\code{sum(results$any_found == 1L)}) work identically for
#'       bootstrap and CV outputs.}
#'     \item{hr_sg}{Numeric. Hazard ratio for top subgroup}
#'     \item{N_sg}{Integer. Sample size of top subgroup}
#'     \item{E_sg}{Integer. Number of events in top subgroup}
#'     \item{K_sg}{Integer. Number of factors defining top subgroup}
#'     \item{g_sg}{Numeric. Subgroup group ID}
#'     \item{m_sg}{Numeric. Subgroup index}
#'     \item{M.1}{Character. First factor label}
#'     \item{M.2}{Character. Second factor label}
#'     \item{M.3}{Character. Third factor label}
#'     \item{M.4}{Character. Fourth factor label}
#'     \item{M.5}{Character. Fifth factor label}
#'     \item{M.6}{Character. Sixth factor label}
#'     \item{M.7}{Character. Seventh factor label}
#'     \item{grf_cuts_b}{Character. GRF policy-tree cut expressions
#'       returned on this bootstrap sample, collapsed with " | " when
#'       multiple cuts are produced.  \code{NA_character_} when GRF was
#'       not used, the bootstrap forestsearch call errored, or GRF
#'       returned no cuts.  Populated independently of whether a
#'       subgroup was identified: GRF may surface a candidate cut on
#'       bootstraps where the downstream consistency stage rejects all
#'       candidates.  Mirrors the \code{grf_cuts} column that
#'       \code{\link{forestsearch_tenfold}} captures in its
#'       \code{fold_summary} return slot.}
#'   }
#'   Rows where no valid subgroup was found will have \code{NA} for bias corrections.
#'   The returned object has a "timing" attribute with summary statistics.
#'
##' @section Bias Correction Methods:
#' Two bias correction approaches are implemented:
#' \enumerate{
#'   \item \strong{Method 1 (Simple Optimism)}:
#'     \deqn{H_{adj1} = H_{obs} - (H^*_{*} - H^*_{obs})}
#'     where \eqn{H^*_{*}} is the new subgroup HR on bootstrap data and
#'     \eqn{H^*_{obs}} is the new subgroup HR on original data.
#'
#'   \item \strong{Method 2 (Double Bootstrap)}:
#'     \deqn{H_{adj2} = 2 \times H_{obs} - (H_{*} + H^*_{*} - H^*_{obs})}
#'     where \eqn{H_{*}} is the original subgroup HR on bootstrap data.
#' }
#' where:
#' \itemize{
#'   \item \code{H_obs}: Original subgroup HR on original data
#'   \item \code{H_star}: Original subgroup HR on bootstrap data
#'   \item \code{Hstar_obs}: New subgroup (found in bootstrap) HR on original data
#'   \item \code{Hstar_star}: New subgroup (found in bootstrap) HR on bootstrap data
#' }
#'
#' @section Computational Details:
#' \itemize{
#'   \item Uses \code{doFuture} backend for parallel execution (configured externally)
#'   \item Sets reproducible seeds: \code{8316951 + boot * 100} for each iteration
#'   \item Each bootstrap iteration runs full ForestSearch pipeline including
#'     variable selection, subgroup search, and consistency evaluation
#'   \item Sequential execution within each bootstrap prevents nested parallelization
#'   \item Failed bootstrap iterations generate warnings but don't stop execution
#'   \item Confounders are removed from bootstrap data to force fresh variable selection
#' }
#'
#' @section Bootstrap Configuration:
#' Each bootstrap iteration modifies ForestSearch arguments to:
#' \itemize{
#'   \item \strong{Suppress output}: \code{details}, \code{show_candidate_summary},
#'     \code{plot.sg}, \code{plot.grf} all set to \code{FALSE}
#'   \item \strong{Force re-selection}: \code{grf_res} and \code{grf_cuts} set to \code{NULL}
#'   \item \strong{Prevent nested parallel}: \code{parallel_args$plan = "sequential"},
#'     \code{workers = 1}
#' }
#'
#' @section Performance Considerations:
#' \itemize{
#'   \item Typical runtime: 1-5 seconds per bootstrap iteration
#'   \item For 1000 bootstraps with 6 workers: ~3-10 minutes total
#'   \item Memory usage scales with dataset size and number of workers
#'   \item Consider reducing \code{nb_boots} for initial testing (e.g., 100)
#' }
#'
#' @section Error Handling:
#' The function gracefully handles three failure modes:
#' \enumerate{
#'   \item Bootstrap sample creation fails: Returns row with all \code{NA}
#'   \item ForestSearch fails to run: Warns and returns row with all \code{NA}
#'   \item ForestSearch runs but finds no subgroup: Returns row with all \code{NA}
#' }
#' All three cases ensure the foreach loop can still combine results via \code{rbind}.
#'
#' @note This function is designed to be called within a \code{foreach} loop
#'   with \code{\%dofuture\%} operator. It requires:
#'   \itemize{
#'     \item All functions in \code{\link{get_bootstrap_exports}} to be available
#'       in the parallel workers
#'     \item Packages listed in \code{BOOTSTRAP_REQUIRED_PACKAGES} to be installed
#'     \item Proper parallel backend setup via \code{\link{setup_parallel_SGcons}}
#'   }
#'
#' @seealso
#' \code{\link{forestsearch_bootstrap_dofuture}} for the wrapper function that
#'   sets up parallelization and calls this function
#' \code{\link{build_cox_formula}} for creating the Cox formula
#' \code{\link{fit_cox_models}} for initial Cox model fitting
#' \code{\link{get_Cox_sg}} for Cox model fitting on subgroups
#' \code{\link{get_dfRes}} for processing bootstrap results into confidence intervals
#' \code{\link{bootstrap_ystar}} for generating the Ystar matrix
#'
#' @examples
#' \dontrun{
#' # Typically called via forestsearch_bootstrap_dofuture()
#' # Manual usage for debugging:
#'
#' # 1. Fit initial ForestSearch model
#' fs_result <- forestsearch(
#'   df.analysis = mydata,
#'   outcome.name = "time",
#'   event.name = "status",
#'   treat.name = "treatment",
#'   confounders.name = c("age", "sex", "stage")
#' )
#'
#' # 2. Build Cox formula
#' cox_formula <- build_cox_formula("time", "status", "treatment")
#'
#' # 3. Get observed estimates
#' cox_fits <- fit_cox_models(fs_result$df.est, cox_formula)
#'
#' # 4. Set up parallel backend
#' library(doFuture)
#' plan(multisession, workers = 6)
#'
#' # 5. Run bootstrap (note: this is already parallelized internally)
#' boot_results <- bootstrap_results(
#'   fs.est = fs_result,
#'   df_boot_analysis = fs_result$df.est,
#'   cox.formula.boot = cox_formula,
#'   nb_boots = 100,
#'   show_three = TRUE,
#'   H_obs = cox_fits$H_obs,
#'   Hc_obs = cox_fits$Hc_obs
#' )
#'
#' # 6. Check results
#' summary(boot_results)
#'
#' # Proportion of bootstraps that found a subgroup
#' mean(!is.na(boot_results$H_biasadj_2))
#' }
#'
#' @family bootstrap functions
#' @importFrom foreach foreach
#' @importFrom data.table data.table rbindlist
#' @importFrom doFuture %dofuture%
#' @export
bootstrap_results <- function(fs.est, df_boot_analysis, cox.formula.boot,
                              nb_boots, show_three, H_obs, Hc_obs,
                              seed = 8316951L,
                              estimator_fn = NULL,
                              boot_index_mat = NULL) {
  # =========================================================================
  # SECTION: INITIALIZE TIMING
  # =========================================================================
  t_start_bootstrap <- proc.time()[3]

  set.seed(seed)

  NN <- nrow(df_boot_analysis)
  id0 <- seq_len(NN)

  # =========================================================================
  # SECTION: VALIDATE OPTIONAL boot_index_mat
  # =========================================================================
  # When supplied, the index matrix is the source of truth for bootstrap
  # row selection; the per-iteration sample.int() call is bypassed and
  # in_boot <- boot_index_mat[boot, ] is used instead.  This removes the
  # implicit doFuture-stream alignment with bootstrap_ystar() that the
  # NULL path relies on.
  if (!is.null(boot_index_mat)) {
    if (!is.matrix(boot_index_mat) || !is.numeric(boot_index_mat)) {
      stop("'boot_index_mat' must be a numeric matrix.", call. = FALSE)
    }
    if (nrow(boot_index_mat) != nb_boots) {
      stop("'boot_index_mat' must have nb_boots rows; got ",
           nrow(boot_index_mat), " rows for nb_boots = ", nb_boots, ".",
           call. = FALSE)
    }
    if (ncol(boot_index_mat) != NN) {
      stop("'boot_index_mat' must have nrow(df_boot_analysis) columns; got ",
           ncol(boot_index_mat), " columns for n = ", NN, ".",
           call. = FALSE)
    }
  }

  # =========================================================================
  # SECTION: HOIST ITERATION-INVARIANT BOOTSTRAP-FRAME PREPARATION
  # =========================================================================
  # The set of columns to drop (`fs.est$confounders.candidate` plus
  # "treat.recommend") and the analysis-frame projection `dfnew` do not
  # vary across bootstrap iterations and were previously recomputed B
  # times inside the foreach worker.  Compute them once here.
  #
  # `keep_cols` is a logical vector of column positions, used inside the
  # worker to subset `df_boot` (which is `df_boot_analysis[in_boot, ]`,
  # i.e. row-resampled but column-identical).  This replaces a per-
  # iteration `!(names(df_boot) %in% drop.vars)` lookup with a
  # precomputed positional index.
  drop_vars_boot <- c(fs.est$confounders.candidate, "treat.recommend")
  keep_cols_boot <- !(names(df_boot_analysis) %in% drop_vars_boot)
  dfnew_predict  <- df_boot_analysis[, keep_cols_boot, drop = FALSE]

  # Pre-extract the forestsearch argument list that each bootstrap
  # iteration clones and customises.  Holding this in a small named
  # binding (rather than reading `fs.est$args_call_all` from inside the
  # worker) prevents doFuture's globals-detection from serialising the
  # entire `fs.est` object - including `df.est`, `find.grps`, GRF
  # outputs, and other slots that the worker never reads - to each
  # parallel worker.  `fs.est` itself is no longer referenced inside
  # the foreach body after this change.
  args_FS_template <- fs.est$args_call_all

  foreach_results <- suppressWarnings({foreach::foreach(
    boot = seq_len(nb_boots),
    .options.future = list(
      seed = TRUE,
      globals = TRUE
    ),
    .errorhandling = "pass"
  ) %dofuture% {

    # =========================================================================
    # SECTION: BOOTSTRAP ITERATION TIMING (PER ITERATION)
    # =========================================================================
    t_iter_start <- proc.time()[3]

    # Simple console feedback (only for first few and milestones)
    if (boot %in% c(1, 10, 50, 100, 250, 500, 750, 1000) || boot == nb_boots) {
      message(sprintf("Bootstrap iteration %d/%d", boot, nb_boots))
    }

    show3 <- FALSE
    if (show_three) show3 <- (boot <= 3)

    # Source bootstrap row indices: explicit pre-generated matrix when
    # provided, otherwise fall back to per-iteration sample.int() (legacy
    # path retained for direct callers of bootstrap_results()).
    if (!is.null(boot_index_mat)) {
      in_boot <- boot_index_mat[boot, ]
    } else {
      in_boot <- sample.int(NN, size = NN, replace = TRUE)
    }
    df_boot <- df_boot_analysis[in_boot, ]
    df_boot$id_boot <- seq_len(nrow(df_boot))

    # =================================================================
    # Extract variable names from formula (survival) or args (GLM)
    # =================================================================
    is_glm_boot <- !is.null(estimator_fn)

    if (!is_glm_boot) {
      outcome_var <- all.vars(cox.formula.boot[[2]])[1]
      event_var <- all.vars(cox.formula.boot[[2]])[2]
      treat_var <- all.vars(cox.formula.boot[[3]])[1]
    } else {
      outcome_var <- args_FS_template$outcome.name
      event_var   <- args_FS_template$event.name
      treat_var   <- args_FS_template$treat.name
    }

    # =================================================================
    # Bootstrap data evaluated at ORIGINAL subgroup H
    # =================================================================

    # Check events/arm-size in subgroup H (treat.recommend == 0)
    df_H <- subset(df_boot, treat.recommend == 0)

    if (!is_glm_boot) {
      events_H_0 <- sum(df_H[df_H[[treat_var]] == 0, event_var], na.rm = TRUE)
      events_H_1 <- sum(df_H[df_H[[treat_var]] == 1, event_var], na.rm = TRUE)

      if (events_H_0 < 5 || events_H_1 < 5) {
        message(sprintf("  Bootstrap %d: Low events in subgroup H (Control=%d, Treat=%d)",
                        boot, events_H_0, events_H_1))
      }
    } else {
      events_H_0 <- sum(df_H[[treat_var]] == 0, na.rm = TRUE)
      events_H_1 <- sum(df_H[[treat_var]] == 1, na.rm = TRUE)
    }

    # Fit model on H subgroup of bootstrap data
    if (!is_glm_boot) {
      fitH_star <- get_Cox_sg(
        df_sg = df_H,
        cox.formula = cox.formula.boot,
        est.loghr = TRUE
      )
    } else {
      fitH_star <- fit_subgroup_effect(df_H, estimator_fn)
    }

    H_star <- fitH_star$est_obs

    # Check events/arm-size in subgroup Hc (treat.recommend == 1)
    df_Hc <- subset(df_boot, treat.recommend == 1)

    if (!is_glm_boot) {
      events_Hc_0 <- sum(df_Hc[df_Hc[[treat_var]] == 0, event_var], na.rm = TRUE)
      events_Hc_1 <- sum(df_Hc[df_Hc[[treat_var]] == 1, event_var], na.rm = TRUE)

      if (events_Hc_0 < 5 || events_Hc_1 < 5) {
        message(sprintf("  Bootstrap %d: Low events in subgroup Hc (Control=%d, Treat=%d)",
                        boot, events_Hc_0, events_Hc_1))
      }
    } else {
      events_Hc_0 <- sum(df_Hc[[treat_var]] == 0, na.rm = TRUE)
      events_Hc_1 <- sum(df_Hc[[treat_var]] == 1, na.rm = TRUE)
    }

    if (!is_glm_boot) {
      fitHc_star <- get_Cox_sg(
        df_sg = df_Hc,
        cox.formula = cox.formula.boot,
        est.loghr = TRUE
      )
    } else {
      fitHc_star <- fit_subgroup_effect(df_Hc, estimator_fn)
    }

    Hc_star <- fitHc_star$est_obs

    # =================================================================
    # Initialize bias corrections and other metrics as NA
    # =================================================================
    H_biasadj_1 <- H_biasadj_2 <- NA
    Hc_biasadj_1 <- Hc_biasadj_2 <- NA
    max_sg_est <- NA
    L <- NA
    max_count <- NA

    # Initialize event counts for NEW subgroups (found in bootstrap)
    events_Hstar_0 <- NA
    events_Hstar_1 <- NA
    events_Hcstar_0 <- NA
    events_Hcstar_1 <- NA

    events_H_0 <- NA
    events_H_1 <- NA
    events_Hc_0 <- NA
    events_Hc_1 <- NA

    # Initialize timing for forestsearch within this iteration
    tmins_search <- NA
    tmins_iteration <- NA

    # =================================================================
    # Prepare bootstrap dataframes - drop confounders and treat.recommend
    # =================================================================
    # `dfnew_predict` and the column-position vector `keep_cols_boot`
    # are precomputed before the foreach loop (invariant across
    # iterations).  Only `dfnew_boot` needs to be built per iteration,
    # since `df_boot` differs each time.
    dfnew_boot <- df_boot[, keep_cols_boot, drop = FALSE]

    # =================================================================
    # Configure forestsearch arguments for bootstrap
    # =================================================================
    # Per-iteration clone of the pre-extracted template; R's copy-on-
    # modify semantics make the subsequent `args_FS_boot$X <- ...`
    # mutations local to this worker iteration.
    args_FS_boot <- args_FS_template
    args_FS_boot$df.analysis <- dfnew_boot
    args_FS_boot$df.predict <- dfnew_predict

    # CATEGORY 1: OUTPUT SUPPRESSION
    args_FS_boot$details <- show3
    args_FS_boot$show_candidate_summary <- FALSE
    args_FS_boot$plot.sg <- FALSE
    args_FS_boot$plot.grf <- FALSE
    args_FS_boot$quiet <- TRUE

    # CATEGORY 2: VARIABLE RE-SELECTION
    args_FS_boot$grf_res <- NULL
    args_FS_boot$grf_cuts <- NULL
    # DINA inputs must also be nulled so each replicate re-fits DINA and
    # re-runs its selection from scratch -- otherwise a supplied/cached DINA
    # fit (use_dina screening or subgroup_method = "dina") would be reused on
    # every resample, defeating the selection-adjustment the bootstrap exists
    # to provide.
    args_FS_boot$dina_res <- NULL
    args_FS_boot$dina_cuts <- NULL

    # CATEGORY 3: SEQUENTIAL EXECUTION
    args_FS_boot$parallel_args$plan <- "sequential"
    args_FS_boot$parallel_args$workers <- 1L
    args_FS_boot$parallel_args$show_message <- FALSE

    # =================================================================
    # Run forestsearch on bootstrap sample (WITH TIMING)
    # =================================================================
    t_fs_start <- proc.time()[3]

    run_bootstrap <- try(do.call(forestsearch, args_FS_boot), TRUE)

    t_fs_end <- proc.time()[3]
    tmins_search <- (t_fs_end - t_fs_start) / 60

    if (inherits(run_bootstrap, "try-error")) {
      warning("Bootstrap ", boot, " failed: ", as.character(run_bootstrap))
    }

    # =================================================================
    # Compute bias corrections if bootstrap succeeded AND found subgroup
    # =================================================================
    if (!inherits(run_bootstrap, "try-error") && !is.null(run_bootstrap$sg.harm)) {

      # Extract prediction datasets from bootstrap ForestSearch run
      df_PredBoot <- run_bootstrap$df.predict
      dfboot_PredBoot <- run_bootstrap$df.est

      # Extract search metrics
      max_sg_est <- as.numeric(run_bootstrap$find.grps$max_sg_est)
      max_count <- run_bootstrap$find.grps$max_count
      L <- run_bootstrap$find.grps$L

      # Extract the identified subgroup

      # ==============================================================
      # Check events/arm-sizes in NEW subgroups found by bootstrap
      # ==============================================================

      # NEW subgroup H* on ORIGINAL data
      df_Hstar <- subset(df_PredBoot, treat.recommend == 0)

      if (!is_glm_boot) {
        events_Hstar_0 <- sum(df_Hstar[df_Hstar[[treat_var]] == 0, event_var], na.rm = TRUE)
        events_Hstar_1 <- sum(df_Hstar[df_Hstar[[treat_var]] == 1, event_var], na.rm = TRUE)

        if (events_Hstar_0 < 5 || events_Hstar_1 < 5) {
          message(sprintf("  Bootstrap %d: Low events in NEW subgroup H* (Control=%d, Treat=%d)",
                          boot, events_Hstar_0, events_Hstar_1))
        }
      } else {
        events_Hstar_0 <- sum(df_Hstar[[treat_var]] == 0, na.rm = TRUE)
        events_Hstar_1 <- sum(df_Hstar[[treat_var]] == 1, na.rm = TRUE)
      }

      # NEW subgroup Hc* on ORIGINAL data
      df_Hcstar <- subset(df_PredBoot, treat.recommend == 1)

      if (!is_glm_boot) {
        events_Hcstar_0 <- sum(df_Hcstar[df_Hcstar[[treat_var]] == 0, event_var], na.rm = TRUE)
        events_Hcstar_1 <- sum(df_Hcstar[df_Hcstar[[treat_var]] == 1, event_var], na.rm = TRUE)

        if (events_Hcstar_0 < 5 || events_Hcstar_1 < 5) {
          message(sprintf("  Bootstrap %d: Low events in NEW subgroup Hc* (Control=%d, Treat=%d)",
                          boot, events_Hcstar_0, events_Hcstar_1))
        }
      } else {
        events_Hcstar_0 <- sum(df_Hcstar[[treat_var]] == 0, na.rm = TRUE)
        events_Hcstar_1 <- sum(df_Hcstar[[treat_var]] == 1, na.rm = TRUE)
      }

      # ==============================================================
      # Compute bias corrections for subgroup H (harm/questionable)
      # ==============================================================

      # Hstar_obs: New subgroup (from bootstrap) evaluated on ORIGINAL data
      if (!is_glm_boot) {
        fitHstar_obs <- get_Cox_sg(
          df_sg = df_Hstar,
          cox.formula = cox.formula.boot,
          est.loghr = TRUE
        )
      } else {
        fitHstar_obs <- fit_subgroup_effect(df_Hstar, estimator_fn)
      }
      Hstar_obs <- fitHstar_obs$est_obs

      # Hstar_star: New subgroup (from bootstrap) evaluated on BOOTSTRAP data
      if (!is_glm_boot) {
        fitHstar_star <- get_Cox_sg(
          df_sg = subset(dfboot_PredBoot, treat.recommend == 0),
          cox.formula = cox.formula.boot,
          est.loghr = TRUE
        )
      } else {
        fitHstar_star <- fit_subgroup_effect(
          subset(dfboot_PredBoot, treat.recommend == 0), estimator_fn
        )
      }
      Hstar_star <- fitHstar_star$est_obs

      # Bias correction method 1: Simple optimism correction
      H_biasadj_1 <- H_obs - (Hstar_star - Hstar_obs)

      # Bias correction method 2: Double correction
      H_biasadj_2 <- 2 * H_obs - (H_star + Hstar_star - Hstar_obs)

      # ==============================================================
      # Compute bias corrections for subgroup H^c (complement/recommend)
      # ==============================================================

      # Hcstar_obs: New subgroup complement evaluated on ORIGINAL data
      if (!is_glm_boot) {
        fitHcstar_obs <- get_Cox_sg(
          df_sg = df_Hcstar,
          cox.formula = cox.formula.boot,
          est.loghr = TRUE
        )
      } else {
        fitHcstar_obs <- fit_subgroup_effect(df_Hcstar, estimator_fn)
      }
      Hcstar_obs <- fitHcstar_obs$est_obs

      # Hcstar_star: New subgroup complement evaluated on BOOTSTRAP data
      if (!is_glm_boot) {
        fitHcstar_star <- get_Cox_sg(
          df_sg = subset(dfboot_PredBoot, treat.recommend == 1),
          cox.formula = cox.formula.boot,
          est.loghr = TRUE
        )
      } else {
        fitHcstar_star <- fit_subgroup_effect(
          subset(dfboot_PredBoot, treat.recommend == 1), estimator_fn
        )
      }
      Hcstar_star <- fitHcstar_star$est_obs

      # Apply same correction methods for H^c
      Hc_biasadj_1 <- Hc_obs - (Hcstar_star - Hcstar_obs)
      Hc_biasadj_2 <- 2 * Hc_obs - (Hc_star + Hcstar_star - Hcstar_obs)
    }

    # =================================================================
    # CALCULATE ITERATION TIMING
    # =================================================================
    t_iter_end <- proc.time()[3]
    tmins_iteration <- (t_iter_end - t_iter_start) / 60

    # =================================================================
    # RETURN: data.table with event counts, TIMING, AND CONSISTENCY RESULTS
    # =================================================================

    # Initialize consistency result columns as NA (FAST: pre-allocation)
    Pcons <- hr_sg <- N_sg <- E_sg <- K_sg <- g_sg <- m_sg <- NA_real_
    M.1 <- M.2 <- M.3 <- M.4 <- M.5 <- M.6 <- M.7 <- NA_character_

    # Initialize Phase C GRF-cut capture.  Populated whenever the
    # bootstrap forestsearch call itself succeeds; independent of whether
    # a subgroup was identified (GRF may return a cut even when the
    # downstream consistency stage rejects every candidate).
    grf_cuts_b <- NA_character_

    # Extract consistency results if available
    if (!inherits(run_bootstrap, "try-error") &&
        !is.null(run_bootstrap$grp.consistency) &&
        !is.null(run_bootstrap$grp.consistency$out_sg) &&
        !is.null(run_bootstrap$grp.consistency$out_sg$result)) {

      sg_result <- run_bootstrap$grp.consistency$out_sg$result

      if (nrow(sg_result) > 0) {
        # Extract first row (FAST: single row access)
        first_row <- sg_result[1, ]

        # Extract values directly (FAST: direct column access)
        # Use data.table's [[ operator which is faster than $
        Pcons <- first_row[["Pcons"]]
        hr_sg <- first_row[["hr"]]
        N_sg <- first_row[["N"]]
        E_sg <- first_row[["E"]]
        K_sg <- first_row[["K"]]
        g_sg <- first_row[["g"]]
        m_sg <- first_row[["m"]]

        # Extract factor labels only up to maxk (FAST: direct access, no loops)
        # Get maxk from args or use K_sg as indicator
        actual_k <- if (!is.na(K_sg)) K_sg else 0

        if (actual_k >= 1) M.1 <- first_row[["M.1"]]
        if (actual_k >= 2) M.2 <- first_row[["M.2"]]
        if (actual_k >= 3) M.3 <- first_row[["M.3"]]
        if (actual_k >= 4) M.4 <- first_row[["M.4"]]
        if (actual_k >= 5) M.5 <- first_row[["M.5"]]
        if (actual_k >= 6) M.6 <- first_row[["M.6"]]
        if (actual_k >= 7) M.7 <- first_row[["M.7"]]
      }
    }

    # Extract GRF cuts independently of subgroup identification.
    # Mirrors the CV-side capture in forestsearch_tenfold(); collapses
    # multi-cut output with " | " to match the CV format exactly.
    if (!inherits(run_bootstrap, "try-error") &&
        !is.null(run_bootstrap$grf_cuts) &&
        length(run_bootstrap$grf_cuts) > 0L) {
      grf_cuts_b <- paste(run_bootstrap$grf_cuts, collapse = " | ")
    }

    # Single data.table creation (FAST: no cbind, no copies)
    dfres <- data.table::data.table(
      boot_id = boot,
      H_biasadj_1 = H_biasadj_1,
      H_biasadj_2 = H_biasadj_2,
      Hc_biasadj_1 = Hc_biasadj_1,
      Hc_biasadj_2 = Hc_biasadj_2,
      max_sg_est = max_sg_est,
      L = L,
      max_count = max_count,

      # Event counts for ORIGINAL subgroup evaluated on BOOTSTRAP sample
      events_H_0 = events_H_0,
      events_H_1 = events_H_1,
      events_Hc_0 = events_Hc_0,
      events_Hc_1 = events_Hc_1,

      # Event counts for NEW subgroup (found in bootstrap) evaluated on ORIGINAL data
      events_Hstar_0 = events_Hstar_0,
      events_Hstar_1 = events_Hstar_1,
      events_Hcstar_0 = events_Hcstar_0,
      events_Hcstar_1 = events_Hcstar_1,

      # TIMING COLUMNS
      tmins_search = tmins_search,
      tmins_iteration = tmins_iteration,

      # CONSISTENCY RESULTS FROM TOP SUBGROUP
      Pcons = Pcons,
      # any_found: integer 0/1, TRUE iff a subgroup was identified in
      # this iteration (equivalent to !is.na(Pcons)).  Mirrors the CV
      # fold_summary$any_found convention so cross-API diagnostics can
      # use the same column name.  Added in v0.2.0.
      any_found = as.integer(!is.na(Pcons)),
      hr_sg = hr_sg,
      N_sg = N_sg,
      E_sg = E_sg,
      K_sg = K_sg,
      g_sg = g_sg,              # Subgroup group ID
      m_sg = m_sg,              # Subgroup index
      M.1 = M.1,
      M.2 = M.2,
      M.3 = M.3,
      M.4 = M.4,
      M.5 = M.5,
      M.6 = M.6,
      M.7 = M.7,

      # GRF POLICY-TREE OUTPUT (Phase C)
      grf_cuts_b = grf_cuts_b
    )

    return(dfres)
  }
  })

  # =========================================================================
  # SECTION: COLLECT FOREACH RESULTS
  # =========================================================================
  # `.errorhandling = "pass"` (combined with no `.combine`) causes the
  # foreach to return a list of length nb_boots in which each element is
  # either a data.table (success) or a condition object (error).  The
  # helper separates the two, surfaces the failure count and first error
  # message via warning() (or stop() on 100% failure), and rbinds the
  # successes via data.table::rbindlist(..., fill = TRUE).  This mirrors
  # the convention established in the package's Quarto calibration suite.
  foreach_results <- .collect_bootstrap_results(foreach_results, nb_boots)

  # =========================================================================
  # SECTION: CALCULATE TOTAL BOOTSTRAP TIMING
  # =========================================================================
  t_end_bootstrap <- proc.time()[3]
  tmins_total_bootstrap <- (t_end_bootstrap - t_start_bootstrap) / 60

  # Add timing summary as attributes
  attr(foreach_results, "timing") <- list(
    total_minutes = tmins_total_bootstrap,
    total_hours = tmins_total_bootstrap / 60,
    n_boots = nb_boots,
    avg_minutes_per_boot = tmins_total_bootstrap / nb_boots,
    avg_seconds_per_boot = (tmins_total_bootstrap * 60) / nb_boots
  )

  # =========================================================================
  # SECTION: RETURN-VALUE SCHEMA VALIDATION
  # =========================================================================
  # Catch accidental column drops/renames before downstream code
  # (summarize_bootstrap_subgroups, print.fs_bootstrap, smoke tests)
  # silently misbehaves.  Covers the columns those consumers actually
  # read, not an exhaustive enumeration of every output column.
  required_cols <- c(
    "boot_id",
    "Pcons", "any_found",
    "hr_sg", "N_sg", "K_sg",
    "M.1", "M.2", "M.3", "M.4", "M.5", "M.6", "M.7",
    "H_biasadj_1", "H_biasadj_2", "Hc_biasadj_1", "Hc_biasadj_2",
    "grf_cuts_b"
  )
  missing_cols <- setdiff(required_cols, names(foreach_results))
  if (length(missing_cols) > 0L) {
    stop(sprintf(
      "bootstrap_results() return value is missing required column(s): %s.\n  This indicates a return-contract regression; downstream consumers (summarize_bootstrap_subgroups, print.fs_bootstrap) will fail or silently misbehave.",
      paste(shQuote(missing_cols), collapse = ", ")
    ), call. = FALSE)
  }

  return(foreach_results)
}
