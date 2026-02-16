#' @title ForestSearch: Exploratory Subgroup Identification
#'
#' @description
#' Identifies subgroups with differential treatment effects in clinical trials
#' using a combination of Generalized Random Forests (GRF), LASSO variable
#' selection, and exhaustive combinatorial search with split-sample validation.
#'
#' @param df.analysis Data frame. Analysis dataset with required columns.
#' @param outcome.name Character. Name of time-to-event outcome variable. Default "tte".
#' @param event.name Character. Name of event indicator (1=event, 0=censored). Default "event".
#' @param treat.name Character. Name of treatment variable (1=treatment, 0=control). Default "treat".
#' @param id.name Character. Name of subject ID variable. Default "id".
#' @param potentialOutcome.name Character. Name of potential outcome variable (optional).
#' @param flag_harm.name Character. Name of true harm flag for simulations (optional).
#' @param confounders.name Character vector. Names of candidate subgroup-defining variables.
#' @param parallel_args List. Parallel processing configuration:
#'   \describe{
#'     \item{plan}{Character. One of "multisession", "multicore", "callr", "sequential"}
#'     \item{workers}{Integer. Number of parallel workers}
#'     \item{show_message}{Logical. Show parallel setup messages}
#'   }
#' @param df.predict Data frame. Prediction dataset (optional).
#' @param df.test Data frame. Test dataset (optional).
#' @param is.RCT Logical. Is this a randomized controlled trial? Default TRUE.
#' @param seedit Integer. Random seed. Default 8316951.
#' @param est.scale Character. Estimation scale ("hr" or "rmst"). Default "hr".
#' @param use_lasso Logical. Use LASSO for variable selection. Default TRUE.
#' @param use_grf Logical. Use GRF for variable importance. Default TRUE.
#' @param grf_res GRF results object (optional, for reuse).
#' @param grf_cuts List. Custom GRF cut points (optional).
#' @param max_n_confounders Integer. Maximum confounders to consider. Default 1000.
#' @param grf_depth Integer. GRF tree depth. Default 2.
#' @param dmin.grf Integer. Minimum events for GRF. Default 12.
#' @param frac.tau Numeric. Fraction of tau for RMST. Default 0.6.
#' @param return_selected_cuts_only Logical. If TRUE (default), GRF returns only cuts from the
#'   tree depth that identified the selected subgroup meeting `dmin.grf`. If FALSE
#'   returns all cuts from all fitted trees (depths 1 through `grf_depth`).
#'   See \code{\link{grf.subg.harm.survival}} for details.
#' @param conf_force Character vector. Variables to force include (optional).
#' @param defaultcut_names Character vector. Default cut variable names (optional).
#' @param cut_type Character. Cut type ("default" or "custom"). Default "default".
#' @param exclude_cuts Character vector. Variables to exclude from cutting (optional).
#' @param replace_med_grf Logical. Replace median with GRF cuts. Default FALSE.
#' @param cont.cutoff Integer. Cutoff for continuous vs categorical. Default 4.
#' @param conf.cont_medians Named numeric vector. Median values for continuous variables (optional).
#' @param conf.cont_medians_force Named numeric vector. Forced median values (optional).
#' @param n.min Integer. Minimum subgroup size. Default 60.
#' @param hr.threshold Numeric. Minimum HR for candidate subgroups. Default 1.25.
#' @param hr.consistency Numeric. Minimum HR for consistency validation. Default 1.0.
#' @param sg_focus Character. Subgroup selection focus. One of "hr", "hrMaxSG", "maxSG",
#'   "hrMinSG", "minSG". Default "hr".
#' @param stop_threshold Numeric. Early stopping threshold for consistency
#'   evaluation. When a candidate subgroup's estimated consistency probability
#'   exceeds this threshold, evaluation stops early. Default 0.95.
#'   \strong{Note:} Automatically reset to NULL when \code{sg_focus} is
#'   "hrMaxSG" or "hrMinSG", as these criteria prioritize hazard ratio in
#'   selection and require full evaluation of all candidates.
#' @param fs.splits Integer. Number of splits for consistency evaluation (or maximum
#'   splits when \code{use_twostage = TRUE}). Default 1000.
#' @param m1.threshold Numeric. Maximum median survival threshold. Default Inf.
#' @param pconsistency.threshold Numeric. Minimum consistency proportion. Default 0.90.
#' @param showten_subgroups Logical. Show top 10 subgroups. Default FALSE.
#' @param d0.min Integer. Minimum control arm events. Default 12.
#' @param d1.min Integer. Minimum treatment arm events. Default 12.
#' @param max.minutes Numeric. Maximum search time in minutes. Default 3.
#' @param minp Numeric. Minimum prevalence threshold. Default 0.025.
#' @param details Logical. Print progress details. Default FALSE.
#' @param maxk Integer. Maximum number of factors per subgroup. Default 2.
#' @param by.risk Integer. Risk table interval. Default 12.
#' @param plot.sg Logical. Plot subgroup survival curves. Default FALSE.
#' @param plot.grf Logical. Plot GRF results. Default FALSE.
#' @param max_subgroups_search Integer. Maximum subgroups to evaluate. Default 10.
#' @param vi.grf.min Numeric. Minimum GRF variable importance. Default -0.2.
#' @param use_twostage Logical. Use two-stage sequential consistency algorithm for
#'   improved performance. Default FALSE for backward compatibility. When TRUE,
#'   \code{fs.splits} becomes the maximum number of splits, and early stopping
#'   is enabled. See Details.
#' @param twostage_args List. Parameters for two-stage algorithm (only used when
#'   \code{use_twostage = TRUE}):
#'   \describe{
#'     \item{n.splits.screen}{Integer. Splits for Stage 1 screening. Default 30.}
#'     \item{screen.threshold}{Numeric. Consistency threshold for Stage 1. Default
#'       is automatically calculated to provide ~2.5 SE margin.}
#'     \item{batch.size}{Integer. Splits per batch in Stage 2. Default 20.}
#'     \item{conf.level}{Numeric. Confidence level for early stopping. Default 0.95.}
#'     \item{min.valid.screen}{Integer. Minimum valid Stage 1 splits. Default 10.}
#'   }
#'
#' @return A list of class "forestsearch" containing:
#'   \describe{
#'     \item{grp.consistency}{Consistency evaluation results including:
#'       \itemize{
#'         \item out_sg: Selected subgroup based on sg_focus
#'         \item sg_focus: Focus criterion used
#'         \item df_flag: Treatment recommendations
#'         \item algorithm: "twostage" or "fixed"
#'         \item n_candidates_evaluated: Number evaluated
#'         \item n_passed: Number passing threshold
#'       }}
#'     \item{find.grps}{Subgroup search results}
#'     \item{confounders.candidate}{Candidate confounders considered}
#'     \item{confounders.evaluated}{Confounders after variable selection}
#'     \item{df.est}{Analysis data with treatment recommendations}
#'     \item{df.predict}{Prediction data with recommendations (if provided)}
#'     \item{df.test}{Test data with recommendations (if provided)}
#'     \item{minutes_all}{Total computation time}
#'     \item{grf_res}{GRF results object}
#'     \item{sg_focus}{Subgroup focus criterion used}
#'     \item{sg.harm}{Selected subgroup definition}
#'     \item{grf_cuts}{GRF cut points used}
#'     \item{prop_maxk}{Proportion of max combinations searched}
#'     \item{max_sg_est}{Maximum subgroup HR estimate}
#'     \item{grf_plot}{GRF plot object (if plot.grf = TRUE)}
#'     \item{args_call_all}{All arguments for reproducibility}
#'   }
#'
#' @details
#' \strong{Algorithm Overview:}
#' \enumerate{
#'   \item \strong{Variable Selection}: GRF identifies variables with treatment
#'     effect heterogeneity; LASSO selects most predictive
#'   \item \strong{Subgroup Discovery}: Exhaustive search over factor combinations
#'     up to \code{maxk}
#'   \item \strong{Consistency Validation}: Split-sample validation ensures
#'     reproducibility
#'   \item \strong{Selection}: Choose subgroup based on \code{sg_focus} criterion
#' }
#'
#' \strong{Two-Stage Consistency Algorithm:}
#' When \code{use_twostage = TRUE}, the consistency evaluation uses an optimized
#' algorithm that can provide 3-10x speedup:
#' \itemize{
#'   \item \strong{Stage 1}: Quick screening with \code{n.splits.screen} splits
#'     eliminates clearly non-viable candidates
#'   \item \strong{Stage 2}: Sequential batched evaluation with early stopping
#'     for candidates passing Stage 1
#' }
#'
#' The two-stage algorithm is recommended for:
#' \itemize{
#'   \item Exploratory analyses with many candidate subgroups
#'   \item Large \code{fs.splits} values (>200)
#'   \item Iterative model development
#' }
#'
#' For final regulatory submissions, \code{use_twostage = FALSE} may be preferred
#' for exact reproducibility.
#'
#' @examples
#' \dontrun{
#' # Example 1: Standard analysis (backward compatible)
#' result <- forestsearch(
#'   df.analysis = trial_data,
#'   sg_focus = "hr",
#'   hr.threshold = 1.25,
#'   pconsistency.threshold = 0.90,
#'   fs.splits = 400,
#'   details = TRUE
#' )
#'
#' # Example 2: Fast exploratory analysis with two-stage
#' result_fast <- forestsearch(
#'   df.analysis = trial_data,
#'   sg_focus = "maxSG",
#'   hr.threshold = 1.15,
#'   pconsistency.threshold = 0.85,
#'   fs.splits = 500,
#'   use_twostage = TRUE,
#'   details = TRUE
#' )
#'
#' # Example 3: Two-stage with custom parameters
#' result_custom <- forestsearch(
#'   df.analysis = trial_data,
#'   sg_focus = "hr",
#'   hr.threshold = 1.3,
#'   pconsistency.threshold = 0.95,
#'   fs.splits = 600,
#'   use_twostage = TRUE,
#'   twostage_args = list(
#'     n.splits.screen = 50,
#'     batch.size = 25,
#'     conf.level = 0.99
#'   ),
#'   parallel_args = list(plan = "multisession", workers = 4),
#'   details = TRUE
#' )
#' }
#'
#' @references
#' \itemize{
#'   \item FDA Guidance for Industry: Enrichment Strategies for Clinical Trials
#'   \item Athey & Imbens (2016). Recursive partitioning for heterogeneous
#'     causal effects. PNAS.
#'   \item Wager & Athey (2018). Estimation and inference of heterogeneous
#'     treatment effects using random forests. JASA.
#' }
#'
#' @seealso
#' \code{\link{subgroup.consistency}} for consistency evaluation details
#' \code{\link{forestsearch_bootstrap_dofuture}} for bootstrap inference
#' \code{\link{forestsearch_Kfold}} for cross-validation
#'
#' @importFrom survival coxph Surv
#' @importFrom grf causal_survival_forest variable_importance
#' @importFrom glmnet cv.glmnet
#' @importFrom data.table data.table setorder
#' @importFrom stats quantile sd median
#' @importFrom stats complete.cases
#' @importFrom future.apply future_lapply
#' @importFrom randomForest randomForest
#' @importFrom weightedsurv df_counting
#' @export

forestsearch <- function(df.analysis,
                         outcome.name = "tte",
                         event.name = "event",
                         treat.name = "treat",
                         id.name = "id",
                         potentialOutcome.name = NULL,
                         flag_harm.name = NULL,
                         confounders.name = NULL,
                         parallel_args = list(plan = "callr", workers = 6, show_message = TRUE),
                         df.predict = NULL,
                         df.test = NULL,
                         is.RCT = TRUE,
                         seedit = 8316951,
                         est.scale = "hr",
                         use_lasso = TRUE,
                         use_grf = TRUE,
                         grf_res = NULL,
                         grf_cuts = NULL,
                         max_n_confounders = 1000,
                         grf_depth = 2,
                         dmin.grf = 12,
                         frac.tau = 0.6,
                         return_selected_cuts_only = TRUE,
                         conf_force = NULL,
                         defaultcut_names = NULL,
                         cut_type = "default",
                         exclude_cuts = NULL,
                         replace_med_grf = FALSE,
                         cont.cutoff = 4,
                         conf.cont_medians = NULL,
                         conf.cont_medians_force = NULL,
                         n.min = 60,
                         hr.threshold = 1.25,
                         hr.consistency = 1.0,
                         sg_focus = "hr",
                         fs.splits = 1000,
                         m1.threshold = Inf,
                         pconsistency.threshold = 0.90,
                         stop_threshold = 0.95,
                         showten_subgroups = FALSE,
                         d0.min = 12,
                         d1.min = 12,
                         max.minutes = 3,
                         minp = 0.025,
                         details = FALSE,
                         maxk = 2,
                         by.risk = 12,
                         plot.sg = FALSE,
                         plot.grf = FALSE,
                         max_subgroups_search = 10,
                         vi.grf.min = -0.2,
                         # NEW: Two-stage consistency parameters
                         use_twostage = TRUE,
                         twostage_args = list()) {

  # ===========================================================================
  # SECTION 1: CAPTURE ALL ARGUMENTS FOR REPRODUCIBILITY
  # ===========================================================================

  args_names <- names(formals())
  args_call_all <- mget(args_names, envir = environment())

  # ===========================================================================
  # SECTION 2: VALIDATE INPUTS
  # ===========================================================================

  # Validate parallel arguments
  if (length(parallel_args) > 0) {
    allowed_plans <- c("multisession", "multicore", "callr", "sequential")
    plan_type <- parallel_args$plan
    n_workers <- parallel_args$workers
    max_cores <- parallel::detectCores()

    if (is.null(plan_type)) {
      stop("parallel_args$plan must be specified.")
    }

    if (!plan_type %in% allowed_plans) {
      stop("plan_type must be one of: ", paste(allowed_plans, collapse = ", "))
    }

    if (is.null(n_workers) || !is.numeric(n_workers) || n_workers < 1) {
      n_workers <- 1
    } else {
      n_workers <- min(n_workers, max_cores)
    }
  }

  # Reset stop_threshold for HR-prioritized sg_focus options
  if (!is.null(stop_threshold) && sg_focus %in% c("hrMaxSG", "hrMinSG")) {
    if (details) {
      cat(
        "Note: stop_threshold reset to NULL for sg_focus = '", sg_focus, "'.\n",
        "Early stopping disabled when HR is prioritized in selection;\n",
        "Also consider increasing max_subgroups_search if search does not appear exhaustive.\n\n",
        sep = ""
      )
    }
    stop_threshold <- NULL
    args_call_all$stop_threshold <- NULL
  }

  # Validate two-stage parameters
 if (!is.logical(use_twostage) || length(use_twostage) != 1) {
    stop("'use_twostage' must be a single logical value (TRUE or FALSE)")
  }

  if (use_twostage && length(twostage_args) > 0) {
    valid_ts_args <- c("n.splits.screen", "screen.threshold", "batch.size",
                       "conf.level", "min.valid.screen")
    invalid_args <- setdiff(names(twostage_args), valid_ts_args)
    if (length(invalid_args) > 0) {
      warning("Unknown twostage_args parameters ignored: ",
              paste(invalid_args, collapse = ", "))
    }
  }

  if (details && use_twostage) {
    cat("\n=== Two-Stage Consistency Evaluation Enabled ===\n")
    cat("Stage 1 screening splits:",
        ifelse(is.null(twostage_args$n.splits.screen), 30, twostage_args$n.splits.screen), "\n")
    cat("Maximum total splits:", fs.splits, "\n")
    cat("Batch size:",
        ifelse(is.null(twostage_args$batch.size), 20, twostage_args$batch.size), "\n")
    cat("================================================\n\n")
  }

  # ===========================================================================
  # SECTION 3: INITIALIZE TIMING AND DATA
  # ===========================================================================

  t.start_all <- proc.time()[3]

  df <- df.analysis
  grf_plot <- NULL

  # ===========================================================================
  # SECTION 3A: GRF CUT GENERATION (if use_grf = TRUE)
  # ===========================================================================

  # If using grf and cuts not already populated, run grf.subg.harm.survival
  if (use_grf && (is.null(grf_res) || is.null(grf_res$tree.cuts))) {
    if (details) {
      cat("GRF stage for cut selection with dmin, tau =", c(dmin.grf, frac.tau), "\n")
      if (return_selected_cuts_only) {
        cat("  return_selected_cuts_only = TRUE: using cuts from selected tree only\n")
      }
    }

    grf_res <- tryCatch(
      grf.subg.harm.survival(
        data = df.analysis,
        confounders.name = confounders.name,
        outcome.name = outcome.name,
        event.name = event.name,
        id.name = id.name,
        treat.name = treat.name,
        frac.tau = frac.tau,
        n.min = n.min,
        dmin.grf = dmin.grf,
        RCT = is.RCT,
        details = details,
        maxdepth = grf_depth,
        seedit = seedit,
        return_selected_cuts_only = return_selected_cuts_only
      ),
      error = function(e) {
        warning("GRF analysis failed: ", e$message)
        return(NULL)
      }
    )

    # Extract GRF cuts if successful
    if (!is.null(grf_res) && !inherits(grf_res, "try-error")) {
      grf_cuts <- grf_res$tree.cuts

      if (details) {
        cat("GRF cuts identified:", length(grf_cuts), "\n")
        if (length(grf_cuts) > 0) {
          cat("  Cuts:", paste(grf_cuts, collapse = ", "), "\n")
        }
        if (!is.null(grf_res$selected_depth)) {
          cat("  Selected tree depth:", grf_res$selected_depth, "\n")
        }
      }

      # Generate GRF plot if requested
      if (plot.grf && !is.null(grf_res$tree)) {
        grf_plot <- tryCatch({
          plot(grf_res$tree)
        }, error = function(e) {
          warning("GRF plot failed: ", e$message)
          NULL
        })
      }
    } else {
      grf_cuts <- NULL
      if (details) {
        cat("GRF analysis did not produce cuts\n")
      }
    }
  }

  # ===========================================================================
  # SECTION 4: DATA PREPARATION (get_FSdata)
  # ===========================================================================

  FSdata <- tryCatch(
    do.call(
      get_FSdata,
      filter_call_args(
        args_call_all,
        get_FSdata,
        list(df.analysis = df.analysis, grf_cuts = grf_cuts)
      )
    ),
    error = function(e) {
      warning("Error in get_FSdata: ", e$message)
      return(NULL)
    }
  )

  if (inherits(FSdata, "try-error") || is.null(FSdata)) {
    warning("FSdata failure - returning NULL result")
    return(list(sg.harm = NULL))
  }

  # Extract FSdata components
  lassoomit <- FSdata$lassoomit
  lassokeep <- FSdata$lassokeep
  df <- FSdata$df

  Y <- df[, outcome.name]
  Event <- df[, event.name]
  Treat <- df[, treat.name]

  FSconfounders.name <- FSdata$confs_names
  confs_labels <- FSdata$confs

  if (is.null(df.predict)) df.predict <- df

  # ===========================================================================
  # SECTION 5: GRF VARIABLE IMPORTANCE SCREENING
  # ===========================================================================

  if (!is.null(vi.grf.min)) {
    X <- as.matrix(df[, FSconfounders.name])
    X <- apply(X, 2, as.numeric)

    tau.rmst <- min(c(max(Y[Treat == 1 & Event == 1]), max(Y[Treat == 0 & Event == 1])))

    if (!is.RCT) {
      cs.forest <- try(suppressWarnings(
        grf::causal_survival_forest(X, Y, Treat, Event,
                                     horizon = 0.9 * tau.rmst, seed = 8316951)
      ), TRUE)
    } else {
      cs.forest <- try(suppressWarnings(
        grf::causal_survival_forest(X, Y, Treat, Event, W.hat = 0.5,
                                     horizon = 0.9 * tau.rmst, seed = 8316951)
      ), TRUE)
    }

    vi.cs <- round(grf::variable_importance(cs.forest), 4)
    vi.cs2 <- data.frame(confs_labels, FSconfounders.name, vi.cs)
    vi.order <- order(vi.cs, decreasing = TRUE)
    vi.cs2 <- vi.cs2[vi.order, ]

    conf.screen <- vi.cs2[, 2]
    vi_ratio <- vi.cs2[, 3] / max(vi.cs2[, 3])
    selected.vars <- which(vi_ratio > vi.grf.min)
    conf.screen <- conf.screen[selected.vars]
    conf.screen <- conf.screen[seq_len(min(length(conf.screen), max_n_confounders))]
  } else {
    conf.screen <- FSconfounders.name
  }

  # ===========================================================================
  # SECTION 6: SUBGROUP SEARCH
  # ===========================================================================

  # Prepare data for subgroup search
  # NOTE: df[, conf.screen] columns are factors with levels {0, 1}.
  # dummy() expands each 2-level factor into two indicator columns
  # (e.g., q1 -> q1.0, q1.1), so the search explores BOTH directions
  # of each cut (subgroup and complement).

  df.confounders <- df[, conf.screen]
  df.confounders <- dummy(df.confounders)


  id <- df[, c(id.name)]
  df.fs <- data.frame(Y, Event, Treat, id, df.confounders)
  Z <- as.matrix(df.confounders)
  colnames(Z) <- names(df.confounders)


  search_overrides <- list(
        Y = Y,
        Event = Event,
        Treat = Treat,
        Z = Z,
        d0.min = d0.min,
        d1.min = d1.min,
        n.min = n.min,
        hr.threshold = hr.threshold,
        max.minutes = max.minutes,
        minp = minp,
        maxk = maxk,
        parallel_workers = 1,  # Ensure sequential execution if inside parallel context
        details = details        # Optionally suppress details
  )

  # Merge and filter arguments
  search_args <- modifyList(args_call_all, search_overrides)
  # Optionally, filter only valid arguments for subgroup.search
  valid_args <- names(formals(subgroup.search))
  search_args <- search_args[names(search_args) %in% valid_args]

  find.grps <- tryCatch(
    do.call(subgroup.search, search_args),
    error = function(e) {
      warning("Error in subgroup.search: ", e$message)
      return(NULL)
    }
  )
  if (is.null(find.grps) || inherits(find.grps, "try-error")) {
    warning("Subgroup search failed")
    return(list(sg.harm = NULL, args_call_all = args_call_all))
  }


  # ===========================================================================
  # SECTION 7: INITIALIZE OUTPUT VARIABLES
  # ===========================================================================

  sg.harm <- NULL
  df.est_out <- NULL
  df.predict_out <- NULL
  df.test_out <- NULL
  grp.consistency <- NULL

  max_sg_est <- find.grps$max_sg_est
  prop_maxk <- find.grps$prop_max_count

  t.end_all <- proc.time()[3]
  t.min_all <- (t.end_all - t.start_all) / 60

  # ===========================================================================
  # SECTION 8: CHECK FOR VALID SUBGROUPS
  # ===========================================================================

  has_subgroups <- FALSE

  if (!is.null(find.grps) &&
      !inherits(find.grps, "try-error") &&
      !is.null(find.grps$out.found) &&
      !is.null(find.grps$out.found$hr.subgroups)) {

    hr_values <- find.grps$out.found$hr.subgroups$HR
    has_subgroups <- any(hr_values > hr.consistency, na.rm = TRUE)
  }

  # ===========================================================================
  # SECTION 9: CONSISTENCY EVALUATION (WITH TWO-STAGE OPTION)
  # ===========================================================================

  if (has_subgroups) {
    # Set plotting parameter if needed
    if (plot.sg && is.null(by.risk)) {
      by.risk <- round(max(Y) / 12, 0)
    }

    if (details) {
      n_candidates <- nrow(find.grps$out.found$hr.subgroups)
      cat("# of candidate subgroups (meeting all criteria) =", n_candidates, "\n")
    }

    # -------------------------------------------------------------------------
    # Build arguments for subgroup.consistency
    # -------------------------------------------------------------------------

    # Base override arguments
    consistency_overrides <- list(
      df = df.fs,
      hr.subgroups = find.grps$out.found$hr.subgroups,
      Lsg = find.grps$L,
      confs_labels = confs_labels,
      n.splits = fs.splits,
      stop_Kgroups = max_subgroups_search,
      # NEW: Pass two-stage parameters
      use_twostage = use_twostage,
      twostage_args = twostage_args
    )

    # Run subgroup consistency analysis with error handling
    grp.consistency <- tryCatch({
      do.call(
        subgroup.consistency,
        filter_call_args(
          args_call_all,
          subgroup.consistency,
          consistency_overrides
        )
      )
    }, error = function(e) {
      warning("Error in subgroup.consistency: ", e$message)
      return(NULL)
    })

    # Handle errors gracefully
    if (is.null(grp.consistency) || inherits(grp.consistency, "try-error")) {
      if (details) {
        cat("Consistency analysis failed - proceeding without results\n")
      }
      grp.consistency <- NULL
    }

    # Update timing
    t.end_all <- proc.time()[3]
    t.min_all <- (t.end_all - t.start_all) / 60

    if (details) {
      cat("Seconds and minutes forestsearch overall =", round(c(t.min_all*60,t.min_all), 4), "\n")
      if (!is.null(grp.consistency$algorithm)) {
        cat("Consistency algorithm used:", grp.consistency$algorithm, "\n")
      }
    }

    # -------------------------------------------------------------------------
    # Process results if consistency analysis succeeded
    # -------------------------------------------------------------------------

    if (!is.null(grp.consistency) && !is.null(grp.consistency$sg.harm)) {
      sg.harm <- grp.consistency$sg.harm

      # Extract prediction datasets
      temp <- grp.consistency$df_flag

      # Merge to analysis data
      df.est_out <- merge(df, temp, by = "id", all.x = TRUE)

      # Return df.predict
      if (!is.null(df.predict)) {
        df.predict_out <- merge(df.predict, temp, by = "id", all.x = TRUE)
      }

      # Return df.test
      if (!is.null(df.test)) {
        df.test_out <- get_dfpred(
          df.predict = df.test,
          sg.harm = grp.consistency$sg.harm,
          version = 2
        )
      }
    }
  } # End has_subgroups

  # ===========================================================================
  # SECTION 10: COMPILE AND RETURN OUTPUT
  # ===========================================================================

  out <- list(
    grp.consistency = grp.consistency,
    find.grps = find.grps,
    confounders.candidate = FSconfounders.name,
    confounders.evaluated = confs_labels,
    df.est = df.est_out,
    df.predict = df.predict_out,
    df.test = df.test_out,
    minutes_all = t.min_all,
    grf_res = grf_res,
    sg_focus = sg_focus,
    sg.harm = sg.harm,
    grf_cuts = grf_cuts,
    prop_maxk = prop_maxk,
    max_sg_est = max_sg_est,
    grf_plot = grf_plot,
    args_call_all = args_call_all,
    # NEW: Include algorithm information
    consistency_algorithm = if (!is.null(grp.consistency)) grp.consistency$algorithm else NA
  )

  # Return early if FSdata or find.grps failed
  if (inherits(FSdata, "try-error") || inherits(find.grps, "try-error")) {
    out <- list(sg.harm = NULL, args_call_all = args_call_all)
  }

  class(out) <- c("forestsearch", "list")
  return(out)
}


