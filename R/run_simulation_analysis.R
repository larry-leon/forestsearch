# =============================================================================
# run_simulation_analysis() — General simulation wrapper
# =============================================================================
#
# Replaces oc_analyses_gbsg.R::run_simulation_analysis(), which was coupled
# to simulate_from_gbsg_dgm() and hardcoded column names (y.sim, event.sim,
# flag.harm, treat).
#
# Design:
#   - Calls simulate_from_dgm() (the general simulator).
#   - Column names are explicit parameters; no GBSG-specific names anywhere.
#   - Helpers (extract_fs_estimates_gen, extract_grf_estimates_gen,
#     run_forestsearch_analysis_gen, run_grf_analysis_gen) are contained here
#     and receive column names from the top-level call.
#   - default_sim_params() / default_grf_params_gen() use general names.
#   - GBSG usage is just: run_simulation_analysis(dgm, ...) with defaults.
# =============================================================================


# =============================================================================
# Default parameter lists
# =============================================================================

# Null-coalescing helper (local to this file).
# Defined here rather than importing rlang::%||% because this file
# runs inside parallel workers where rlang may not be loaded.
if (!exists("%||%", inherits = FALSE))
  `%||%` <- function(a, b) if (!is.null(a)) a else b


#' Default ForestSearch parameters (general)
#'
#' Returns the default \code{fs_params} list used by
#' \code{run_simulation_analysis()}.  The \code{parallel_args} entry
#' defaults to \code{list(plan = "sequential")} because
#' \code{run_simulation_analysis()} is designed to be called inside a
#' \code{foreach() \%dofuture\%} loop (one replicate per worker).  Running
#' the inner \code{forestsearch()} multisession in that context produces
#' nested parallelism: each outer worker tries to spawn its own pool of
#' inner workers, which \code{parallelly} rejects with a 300\%-load hard
#' limit error.  Users who call \code{run_simulation_analysis()} once at
#' the top level (i.e., not inside \code{\%dofuture\%}) can opt back into
#' multisession by passing
#' \code{fs_params = list(parallel_args = list(plan = "multisession", workers = N))}.
#'
#' @keywords internal
default_sim_params <- function() {
  list(
    outcome.name             = "y_sim",
    event.name               = "event_sim",
    treat.name               = "treat_sim",
    id.name                  = "id",
    use_lasso                = TRUE,
    use_grf                  = FALSE,
    hr.threshold             = 1.25,
    hr.consistency           = 1.0,
    pconsistency.threshold   = 0.90,
    fs.splits                = 400,
    n.min                    = 60,
    d0.min                   = 12,
    d1.min                   = 12,
    maxk                     = 2,
    max.minutes              = 5,
    by.risk                  = 12,
    vi.grf.min               = -0.2,
    use_twostage             = FALSE,
    twostage_args            = list(),
    # Run inner forestsearch() sequentially by default.  See function
    # documentation above for rationale (nested-parallelism avoidance).
    parallel_args            = list(plan = "sequential")
  )
}

#' Default GRF parameters (general)
#' @keywords internal
default_grf_params_gen <- function() {
  list(
    outcome.name = "y_sim",
    event.name   = "event_sim",
    treat.name   = "treat_sim",
    id.name      = "id",
    n.min        = 60,
    dmin.grf     = 12,
    frac.tau     = 0.60,
    maxdepth     = 2,
    RCT          = TRUE,
    sg.criterion = "mDiff",
    seedit       = 8316951L
  )
}


# =============================================================================
# Internal: generalized-method dispatch helpers
# =============================================================================

#' NULL-preserving list merge
#'
#' Like \code{utils::modifyList()}, but an entry in \code{over} whose value is
#' \code{NULL} \emph{sets} the corresponding element of \code{base} to an
#' explicit \code{NULL} (kept in the list) rather than deleting it.  Required
#' so that an explicit \code{n.min = NULL} survives the parameter merges and
#' reaches \code{forestsearch()}, where it triggers the sample-size-adaptive
#' minimum-size rule; ordinary \code{modifyList()} would strip it and silently
#' fall back to the formal default.
#' @keywords internal
#' @noRd
.modify_keep_null <- function(base, over) {
  for (nm in names(over)) base[nm] <- over[nm]
  base
}

#' Expand a canonical method label into forestsearch() overrides
#'
#' Parses a \code{"<engine>-<focus>"} label (e.g. \code{"fs-effMaxSG"},
#' \code{"dina-eff"}, \code{"grf-tree"}) into a list of \code{forestsearch()}
#' argument overrides.  The engine token names the \emph{selector}; every arm
#' is ultimately run through \code{forestsearch()}, which supplies the shared
#' post-selection inference, so no "FS" component appears in the label.
#' @keywords internal
#' @noRd
.expand_method_label <- function(label) {
  parts <- strsplit(label, "-", fixed = TRUE)[[1L]]
  if (length(parts) < 2L) {
    stop(sprintf(
      "Cannot parse method label '%s'; expected '<engine>-<focus>' (e.g. 'fs-effMaxSG').",
      label), call. = FALSE)
  }
  engine <- parts[1L]
  focus  <- paste(parts[-1L], collapse = "-")
  over   <- list(selection_rule = "neighborhood")
  if (identical(engine, "fs")) {
    over$subgroup_method <- "consistency"
    over$sg_focus        <- focus
    over$use_grf         <- TRUE
    over$use_dina        <- TRUE
    over$use_lasso       <- FALSE
    over["n.min"]        <- list(NULL)   # explicit NULL -> adaptive rule
  } else if (identical(engine, "dina")) {
    over$subgroup_method <- "dina"
    over$sg_focus        <- focus
  } else if (identical(engine, "grf")) {
    over$subgroup_method <- "grf"
    if (identical(focus, "tree")) {
      over$grf_selection <- "tree"
    } else {
      over$grf_selection <- "frontier"
      over$sg_focus      <- focus
    }
  } else {
    stop(sprintf(
      "Unknown engine '%s' in method label '%s'; expected 'fs', 'dina', or 'grf'.",
      engine, label), call. = FALSE)
  }
  over
}

#' Normalize the `methods` argument to a named list of override lists
#'
#' Accepts either a character vector of canonical labels (expanded via
#' \code{.expand_method_label()}) or a pre-built named list of
#' \code{forestsearch()} override lists, returning a named list keyed by arm
#' label.  Returns \code{NULL} unchanged.
#' @keywords internal
#' @noRd
.normalize_methods <- function(methods) {
  if (is.null(methods)) return(NULL)
  if (is.character(methods)) {
    if (anyNA(methods) || any(!nzchar(methods))) {
      stop("All entries of a character 'methods' vector must be non-empty labels.",
           call. = FALSE)
    }
    specs <- lapply(methods, .expand_method_label)
    names(specs) <- methods
    return(specs)
  }
  if (is.list(methods)) {
    nm <- names(methods)
    if (is.null(nm) || any(!nzchar(nm))) {
      stop("When 'methods' is a list, every element must be named with its arm label.",
           call. = FALSE)
    }
    if (anyDuplicated(nm)) {
      stop("Duplicate arm labels in 'methods': ",
           paste(unique(nm[duplicated(nm)]), collapse = ", "), ".", call. = FALSE)
    }
    return(methods)
  }
  stop("'methods' must be NULL, a character vector of labels, or a named list of override lists.",
       call. = FALSE)
}

# =============================================================================
# Main function
# =============================================================================

#' Run One Simulation Replicate
#'
#' General replacement for the legacy \code{run_simulation_analysis()} that
#' was coupled to \code{simulate_from_gbsg_dgm()} and GBSG-specific column
#' names.  This version calls \code{\link{simulate_from_dgm}} and accepts
#' explicit column-name parameters, making it applicable to any DGM built
#' with \code{\link{generate_aft_dgm_flex}}.
#'
#' @param sim_id Integer. Simulation replicate index (used as seed offset).
#' @param dgm An \code{"aft_dgm_flex"} object from
#'   \code{\link{generate_aft_dgm_flex}} or \code{\link{setup_gbsg_dgm}},
#'   or a \code{"glm_dgm"} object from \code{\link{generate_glm_dgm}}.
#' @param n_sample Integer. Per-replicate sample size.
#' @param analysis_time Numeric. Calendar time of analysis on the DGM time
#'   scale.  Use \code{Inf} (default) for no administrative censoring —
#'   equivalent to the legacy \code{max_follow = Inf}.
#' @param cens_adjust Numeric. Log-scale shift to censoring times passed to
#'   \code{simulate_from_dgm(cens_adjust = ...)}. Replaces legacy
#'   \code{muC_adj}. Default \code{0}.
#' @param max_follow \strong{Deprecated.} Use \code{analysis_time} instead.
#'   If supplied, its value is forwarded to \code{analysis_time} with a
#'   warning. Retained for backward compatibility with legacy scripts.
#' @param muC_adj \strong{Deprecated.} Use \code{cens_adjust} instead.
#'   If supplied, its value is forwarded to \code{cens_adjust} with a
#'   warning. Retained for backward compatibility with legacy scripts.
#' @param confounders_base Character vector of base confounder names.
#' @param n_add_noise Integer. Number of independent N(0,1) noise variables
#'   to append. Default \code{0L}.
#' @param outcome_name Name of the observed time column in simulated data.
#'   Default \code{"y_sim"}.
#' @param event_name Name of the event indicator column. Default
#'   \code{"event_sim"}.
#' @param treat_name Name of the treatment column. Default \code{"treat_sim"}.
#' @param harm_col Name of the true-subgroup indicator column. Default
#'   \code{"flag_harm"}.
#' @param run_fs Logical. Run ForestSearch (LASSO). Default \code{TRUE}.
#' @param run_fs_grf Logical. Run ForestSearch (LASSO + GRF). Default
#'   \code{TRUE}.
#' @param run_grf Logical. Run standalone GRF. Default \code{TRUE}.
#' @param methods Optional specification of analysis arms that generalizes the
#'   \code{run_fs} / \code{run_fs_grf} / \code{run_grf} booleans.  When
#'   non-\code{NULL} the booleans are ignored and \strong{every} arm is run via
#'   \code{forestsearch()} (through \code{.run_fs_analysis_gen()}), so all arms
#'   --- including \code{dina-*} and \code{grf-*} --- inherit the bias-corrected
#'   post-selection inference.  Two forms are accepted:
#'   \itemize{
#'     \item A character vector of canonical \code{"<engine>-<focus>"} labels,
#'       e.g. \code{c("fs-effMaxSG", "dina-effMaxSG", "grf-tree")}, expanded
#'       internally to per-arm \code{forestsearch()} overrides.  Engines are
#'       \code{"fs"} (\code{subgroup_method = "consistency"}), \code{"dina"},
#'       and \code{"grf"} (\code{grf_selection = "frontier"}, or \code{"tree"}
#'       for the \code{grf-tree} label).  Defaults applied during expansion:
#'       \code{selection_rule = "neighborhood"}; for \code{fs-*} arms
#'       \code{use_grf = TRUE}, \code{use_dina = TRUE}, \code{use_lasso = FALSE},
#'       and \code{n.min = NULL} (the sample-size-adaptive rule).
#'     \item A named list of explicit \code{forestsearch()} argument-override
#'       lists, label = name, e.g.
#'       \code{list("fs-effMaxSG" = list(subgroup_method = "consistency", sg_focus = "effMaxSG", use_grf = TRUE, use_dina = TRUE, use_lasso = FALSE), "grf-tree" = list(subgroup_method = "grf", grf_selection = "tree"))}.
#'   }
#'   Each arm's overrides are merged on top of \code{default_sim_params()}
#'   overlaid with \code{fs_params}; per-arm overrides take precedence.  The
#'   merge preserves an explicit \code{n.min = NULL}.  Default \code{NULL}
#'   preserves the legacy boolean behavior.
#' @param fs_params Named list of ForestSearch parameter overrides.  Any
#'   element of \code{default_sim_params()} can be overridden, including
#'   \code{parallel_args} (see \emph{Parallel Processing} section below).
#'   Passing \code{NULL} or an empty list uses defaults throughout.
#' @param grf_params Named list of GRF parameter overrides.
#' @param cox_formula Optional Cox formula for unadjusted ITT.
#' @param cox_formula_adj Optional adjusted Cox formula.
#' @param n_sims_total Integer. Total simulations (for progress messages).
#' @param seed_base Integer. Base seed; replicate seed = \code{seed_base +
#'   sim_id}. Default \code{8316951L}.
#' @param verbose Logical. Print progress messages. Default \code{FALSE}.
#' @param verbose_n Integer. If set, only print for \code{sim_id <=
#'   verbose_n}. Default \code{NULL}.
#' @param debug Logical. Print detailed debug output. Default \code{FALSE}.
#' @param keep Character vector. Optional names of heavy diagnostic
#'   objects to attach as a list in \code{attr(result, "diagnostics")}.
#'   Allowed values (any combination):
#'   \describe{
#'     \item{\code{"fs_full"}}{The full \code{forestsearch()} result
#'       object (minus the \code{df.est} analysis data frame, which is
#'       stripped to control memory).  Includes \code{sg.harm},
#'       \code{grp.consistency}, \code{confounders.candidate},
#'       \code{confounders.evaluated}, GRF cuts, and timing.}
#'     \item{\code{"df.est"}}{The analysis data frame with
#'       \code{treat.recommend} column attached.  Heavy; use sparingly.}
#'     \item{\code{"candidate_table"}}{The data table of candidate
#'       subgroups evaluated
#'       (\code{fs_full$grp.consistency$out_sg$result}).  Medium-weight;
#'       useful for diagnosing which candidates were near-threshold.}
#'     \item{\code{"grf_full"}}{The full
#'       \code{grf.subg.harm.*()} result object.}
#'   }
#'   Default \code{character(0)} (attach nothing; current behavior).
#' @param keep_first_n Integer.  When \code{keep} is non-empty, attach
#'   the requested diagnostics only if \code{sim_id <= keep_first_n}.
#'   Default \code{Inf} (attach for all replicates).  Use this to
#'   capture full details for, e.g., the first 20 simulations without
#'   bloating the result object for long production runs.
#'
#' @return A \code{data.table} with one row per analysis method.
#'   Always-present scalar columns cover: true-subgroup size/effect
#'   from the DGM, per-method detection flag (\code{any.H}), estimated
#'   subgroup size (\code{size.H}, \code{size.Hc}), effect estimates
#'   (\code{hr.H.hat}, \code{hr.Hc.hat}, \code{hr.itt}, and their true
#'   counterparts), classification metrics against the DGM-stored
#'   subgroup (\code{sens}, \code{spec}, \code{ppv}, \code{npv}), and
#'   pairwise concordance between methods (\code{agree.*.*},
#'   \code{kappa.*.*}).  As of v0.2.0 the following ForestSearch
#'   diagnostic columns are also always populated (\code{NA} for GRF
#'   rows):
#'   \describe{
#'     \item{\code{sg.def}}{Character.  Cut-expression string of the
#'       identified subgroup, formed as
#'       \code{paste(fs_full$sg.harm, collapse = " & ")}.  Empty
#'       string when no subgroup was identified.  Read across all
#'       replicates to diagnose which cuts the method is selecting
#'       (e.g., \code{table(results$sg.def)}).}
#'     \item{\code{sg.n_factors}}{Integer.  Number of factors in the
#'       identified subgroup definition.}
#'     \item{\code{n_candidates_evaluated}}{Integer.  Number of
#'       candidate subgroups actually evaluated for consistency
#'       (before or at the early-stop candidate).}
#'     \item{\code{n_candidates_total}}{Integer.  Total number of
#'       candidate subgroups available.}
#'     \item{\code{n_passed}}{Integer.  Number of candidates meeting
#'       the consistency threshold.}
#'     \item{\code{consistency_algorithm}}{Character.  \code{"fixed"}
#'       or \code{"twostage"}.}
#'     \item{\code{early_stop_triggered}}{Logical.  Whether
#'       evaluation stopped early due to \code{stop_threshold}.}
#'     \item{\code{fs_minutes}}{Numeric.  Wall time for the FS run
#'       in minutes.}
#'   }
#'   The analogous GRF row carries \code{sg.def} (from
#'   \code{grf_result$sg.harm.id}) and \code{grf_selected_depth}.
#'   When \code{keep} is non-empty and \code{sim_id <= keep_first_n},
#'   the heavy objects are attached as
#'   \code{attr(result, "diagnostics")} — a named list with one
#'   element per keeper (\code{fs_full}, \code{grf_full}, etc.).
#'
#' @seealso \code{\link{simulate_from_dgm}},
#'   \code{\link{generate_aft_dgm_flex}}, \code{\link{setup_gbsg_dgm}}
#'
#' @section GLM Parameters:
#' GLM-specific parameters (\code{outcome_type}, \code{effect_measure},
#' \code{offset.name}) must be passed inside \code{fs_params}, not as
#' top-level arguments.  They route only to the estimation step
#' (\code{.extract_fs_estimates_gen}, \code{.extract_grf_estimates_gen}),
#' not to \code{forestsearch()} itself, which uses Cox PH for subgroup
#' identification in v0.1.x.  Passing these as top-level arguments will
#' result in them being silently ignored.
#'
#' @section Parallel Processing:
#' \code{run_simulation_analysis()} is designed to be called inside a
#' \code{foreach() \%dofuture\%} loop (one replicate per worker).  In
#' that idiom, \emph{outer} parallelism is provided by \code{\%dofuture\%}
#' (one replicate per worker) and the \emph{inner} \code{forestsearch()}
#' call should run sequentially within each worker.  Running both layers
#' multisession produces nested parallelism: each outer worker tries to
#' spawn its own pool of inner workers, which \code{parallelly} rejects
#' with a 300\%-load hard-limit error.
#'
#' To prevent this, \code{default_sim_params()} sets
#' \code{parallel_args = list(plan = "sequential")} as the default for
#' the inner \code{forestsearch()} call.  Users who call
#' \code{run_simulation_analysis()} once at the top level (e.g., for
#' interactive debugging) and want the inner pipeline to run in parallel
#' can opt back in by passing
#' \code{fs_params = list(parallel_args = list(plan = "multisession", workers = N))}.
#'
#' @importFrom data.table data.table rbindlist
#' @importFrom survival coxph Surv
#' @importFrom stats rnorm
#' @examples
#' \dontrun{
#' dgm <- setup_gbsg_dgm(model = "null", verbose = FALSE)
#'
#' # Inner forestsearch() runs sequentially by default; safe for both
#' # interactive use and embedding inside foreach() %dofuture% loops.
#' result <- run_simulation_analysis(sim_id = 1, dgm = dgm, n_sample = 500)
#'
#' # Opt back into inner multisession for a one-off top-level call:
#' result <- run_simulation_analysis(
#'   sim_id    = 1,
#'   dgm       = dgm,
#'   n_sample  = 500,
#'   fs_params = list(parallel_args = list(plan = "multisession",
#'                                         workers = 4))
#' )
#' }
#' @export
run_simulation_analysis <- function(
    sim_id,
    dgm,
    n_sample,
    analysis_time    = Inf,
    cens_adjust      = 0,
    # Deprecated legacy parameter names — emit a warning if used
    max_follow       = NULL,
    muC_adj          = NULL,
    confounders_base = c("v1", "v2", "v3", "v4", "v5", "v6", "v7"),
    n_add_noise      = 0L,
    outcome_name     = "y_sim",
    event_name       = "event_sim",
    treat_name       = "treat_sim",
    harm_col         = "flag_harm",
    run_fs           = TRUE,
    run_fs_grf       = TRUE,
    run_grf          = TRUE,
    methods          = NULL,
    fs_params        = list(),
    grf_params       = list(),
    cox_formula      = NULL,
    cox_formula_adj  = NULL,
    n_sims_total     = NULL,
    seed_base        = 8316951L,
    verbose          = FALSE,
    verbose_n        = NULL,
    debug            = FALSE,
    keep             = character(0),
    keep_first_n     = Inf
) {

  show_verbose <- verbose && (is.null(verbose_n) || sim_id <= verbose_n)

  # ── Validate keep / keep_first_n ─────────────────────────────────────────
  valid_keep <- c("fs_full", "df.est", "candidate_table", "grf_full")
  keep_invalid <- setdiff(keep, valid_keep)
  if (length(keep_invalid) > 0L) {
    stop(sprintf(
      "Invalid 'keep' values: %s. Allowed: %s.",
      paste(shQuote(keep_invalid), collapse = ", "),
      paste(shQuote(valid_keep), collapse = ", ")
    ), call. = FALSE)
  }
  keep_this_sim <- length(keep) > 0L && sim_id <= keep_first_n

  # ── Deprecated parameter shims ─────────────────────────────────────────────
  if (!is.null(max_follow)) {
    warning(
      "'max_follow' is deprecated. Use 'analysis_time' instead.\n",
      "  Mapping: analysis_time = max_follow (", max_follow, ")",
      call. = FALSE
    )
    analysis_time <- max_follow
  }
  if (!is.null(muC_adj)) {
    warning(
      "'muC_adj' is deprecated. Use 'cens_adjust' instead.\n",
      "  Mapping: cens_adjust = muC_adj (", muC_adj, ")",
      call. = FALSE
    )
    cens_adjust <- muC_adj
  }

  if (show_verbose) {
    message("\n", paste(rep("=", 60), collapse = ""))
    message(sprintf("Simulation %d", sim_id))
    if (!is.null(n_sims_total))
      message(sprintf("  Progress: %d / %d (%.1f%%)",
                      sim_id, n_sims_total, 100 * sim_id / n_sims_total))
    message(paste(rep("=", 60), collapse = ""))
  }

  # ── Simulate data ──────────────────────────────────────────────────────────
  if (show_verbose)
    message(sprintf("\n[1] Simulating data (n=%d, analysis_time=%s, cens_adjust=%g)...",
                    n_sample, analysis_time, cens_adjust))

  # Class-based dispatch: glm_dgm uses simulate_from_glm_dgm(),
  # aft_dgm_flex uses simulate_from_dgm().
  is_glm_dgm <- inherits(dgm, "glm_dgm")

  sim_data <- tryCatch(
    if (is_glm_dgm) {
      simulate_from_glm_dgm(
        dgm  = dgm,
        n    = n_sample,
        seed = seed_base + sim_id
      )
    } else {
      simulate_from_dgm(
        dgm           = dgm,
        n             = n_sample,
        seed          = seed_base + sim_id,
        analysis_time = analysis_time,
        cens_adjust   = cens_adjust
      )
    },
    error = function(e) stop("Simulation failed: ", e$message)
  )

  # GLM sims don't produce event_sim; create it for downstream compatibility.
  # For binary: event_sim = y_sim (the outcome IS the event).
  # For continuous/count: event_sim = 1 (all observed, no censoring).
  if (is_glm_dgm && !event_name %in% names(sim_data)) {
    if (!is.null(dgm$outcome_type) && dgm$outcome_type == "binary") {
      sim_data[[event_name]] <- sim_data[[outcome_name]]
    } else {
      sim_data[[event_name]] <- 1L
    }
  }

  if (show_verbose)
    message(sprintf("    n=%d  events=%d (%.1f%%)",
                    nrow(sim_data), sum(sim_data[[event_name]]),
                    100 * mean(sim_data[[event_name]])))

  # ── Optional noise variables ───────────────────────────────────────────────
  confounders_name <- confounders_base
  if (n_add_noise > 0L) {
    set.seed(seed_base + 1000L * sim_id)
    noise_nms <- paste0("noise", seq_len(n_add_noise))
    for (nm in noise_nms) sim_data[[nm]] <- stats::rnorm(nrow(sim_data))
    confounders_name <- c(confounders_base, noise_nms)
    if (show_verbose)
      message(sprintf("    Added %d noise variable(s)", n_add_noise))
  }

  # ── True subgroup properties ───────────────────────────────────────────────
  has_harm <- harm_col %in% names(sim_data)
  size_H_true  <- if (has_harm) sum(sim_data[[harm_col]])        else NA_integer_
  prop_H_true  <- if (has_harm) mean(sim_data[[harm_col]])       else NA_real_
  size_Hc_true <- if (has_harm) nrow(sim_data) - size_H_true    else NA_integer_
  prop_Hc_true <- if (has_harm) 1 - prop_H_true                 else NA_real_

  hr_H_dgm   <- get_dgm_hr(dgm, "hr_H")
  hr_Hc_dgm  <- get_dgm_hr(dgm, "hr_Hc")
  ahr_H_dgm  <- get_dgm_hr(dgm, "ahr_H")
  ahr_Hc_dgm <- get_dgm_hr(dgm, "ahr_Hc")

  if (show_verbose) {
    message(sprintf("\n[2] True subgroup: H n=%s (%.1f%%), Hc n=%s (%.1f%%)",
                    ifelse(is.na(size_H_true), "NA", size_H_true),
                    100 * prop_H_true,
                    ifelse(is.na(size_Hc_true), "NA", size_Hc_true),
                    100 * prop_Hc_true))
    message(sprintf("    DGM HR_H=%.3f  HR_Hc=%.3f  AHR_H=%.3f  AHR_Hc=%.3f",
                    hr_H_dgm, hr_Hc_dgm,
                    ifelse(is.na(ahr_H_dgm), NA, ahr_H_dgm),
                    ifelse(is.na(ahr_Hc_dgm), NA, ahr_Hc_dgm)))
  }

  df_pop <- data.table::data.table(
    sim          = sim_id,
    sizeH_true   = size_H_true,
    propH_true   = prop_H_true,
    sizeHc_true  = size_Hc_true,
    propHc_true  = prop_Hc_true
  )

  # ── Merge parameter defaults ───────────────────────────────────────────────
  fs_defaults  <- default_sim_params()
  grf_defaults <- default_grf_params_gen()
  grf_merged   <- modifyList(grf_defaults, grf_params)

  # Ensure outcome/event/treat names in FS defaults match what simulate produced
  fs_defaults$outcome.name <- outcome_name
  fs_defaults$event.name   <- event_name
  fs_defaults$treat.name   <- treat_name
  grf_merged$outcome.name  <- outcome_name
  grf_merged$event.name    <- event_name
  grf_merged$treat.name    <- treat_name

  # ── Run analyses ───────────────────────────────────────────────────────────
  results_list <- list()
  sg_hat_list  <- list()   # per-subject subgroup assignments for concordance

  # Generalized method dispatch.  When `methods` is supplied it takes over
  # from the run_fs / run_fs_grf / run_grf booleans, and every arm -- fs-*,
  # dina-*, grf-* -- is executed through forestsearch() so all share the
  # bias-corrected post-selection inference.
  method_specs <- .normalize_methods(methods)

  if (show_verbose) {
    to_run <- if (is.null(method_specs)) {
      c(if (run_fs) "FS", if (run_fs_grf) "FSlg", if (run_grf) "GRF")
    } else {
      names(method_specs)
    }
    message(sprintf("\n[3] Running: %s", paste(to_run, collapse = ", ")))
  }

  # Accumulator for heavy diagnostic objects (populated only when
  # keep_this_sim = TRUE).  Each element is itself a named list carrying
  # whichever keepers were requested for that method.
  diagnostics_list <- list()

  if (is.null(method_specs) && run_fs) {
    fs_p <- modifyList(
      modifyList(fs_defaults, list(use_lasso = TRUE, use_grf = FALSE)),
      fs_params
    )
    fs_out <- .run_fs_analysis_gen(
        data = sim_data, confounders_name = confounders_name,
        params = fs_p, dgm = dgm,
        cox_formula = cox_formula, cox_formula_adj = cox_formula_adj,
        outcome_name = outcome_name, event_name = event_name,
        treat_name = treat_name, harm_col = harm_col,
        label = "FS", verbose = show_verbose,
        keep = if (keep_this_sim) keep else character(0)
      )
    sg_hat_list[["FS"]] <- attr(fs_out, "sg_hat")
    if (keep_this_sim) {
      keepers <- attr(fs_out, "keepers")
      if (!is.null(keepers)) diagnostics_list[["FS"]] <- keepers
    }
    # Strip attributes before cbind so data.table semantics stay clean
    attr(fs_out, "sg_hat")  <- NULL
    attr(fs_out, "keepers") <- NULL
    results_list[["FS"]] <- cbind(df_pop, fs_out)
  }

  if (is.null(method_specs) && run_fs_grf) {
    fs_p <- modifyList(
      modifyList(fs_defaults, list(use_lasso = TRUE, use_grf = TRUE)),
      fs_params
    )
    fslg_out <- .run_fs_analysis_gen(
        data = sim_data, confounders_name = confounders_name,
        params = fs_p, dgm = dgm,
        cox_formula = cox_formula, cox_formula_adj = cox_formula_adj,
        outcome_name = outcome_name, event_name = event_name,
        treat_name = treat_name, harm_col = harm_col,
        label = "FSlg", verbose = show_verbose,
        keep = if (keep_this_sim) keep else character(0)
      )
    sg_hat_list[["FSlg"]] <- attr(fslg_out, "sg_hat")
    if (keep_this_sim) {
      keepers <- attr(fslg_out, "keepers")
      if (!is.null(keepers)) diagnostics_list[["FSlg"]] <- keepers
    }
    attr(fslg_out, "sg_hat")  <- NULL
    attr(fslg_out, "keepers") <- NULL
    results_list[["FSlg"]] <- cbind(df_pop, fslg_out)
  }

  if (is.null(method_specs) && run_grf) {
    # Resolve id_name from merged GRF params (falls back to "id")
    id_name_resolved <- grf_merged$id.name %||% "id"
    # Effective names for estimation (see FS sync comment above)
    grf_est_outcome <- fs_params[["outcome.name"]] %||% outcome_name
    grf_est_event   <- fs_params[["event.name"]]   %||% event_name
    grf_out <- .run_grf_analysis_gen(
        data = sim_data, confounders_name = confounders_name,
        params = grf_merged, dgm = dgm,
        cox_formula = cox_formula, cox_formula_adj = cox_formula_adj,
        outcome_name = grf_est_outcome, event_name = grf_est_event,
        treat_name = treat_name, harm_col = harm_col,
        label = "GRF", verbose = show_verbose, debug = debug,
        outcome_type   = fs_params$outcome_type,
        effect_measure = fs_params$effect_measure,
        offset_name    = fs_params$offset.name,
        id_name        = id_name_resolved,
        keep = if (keep_this_sim) keep else character(0)
      )
    sg_hat_list[["GRF"]] <- attr(grf_out, "sg_hat")
    if (keep_this_sim) {
      keepers <- attr(grf_out, "keepers")
      if (!is.null(keepers)) diagnostics_list[["GRF"]] <- keepers
    }
    attr(grf_out, "sg_hat")  <- NULL
    attr(grf_out, "keepers") <- NULL
    results_list[["GRF"]] <- cbind(df_pop, grf_out)
  }

  # ── Generalized method dispatch (methods != NULL) ────────────────────────
  # Each arm is run via forestsearch() through .run_fs_analysis_gen(), with
  # subgroup_method / grf_selection / sg_focus / screening knobs set per arm.
  # The arm label names the selector ("dina-effMaxSG", "grf-tree", ...);
  # forestsearch() is the shared inference backbone and is absent from labels.
  if (!is.null(method_specs)) {
    for (arm_label in names(method_specs)) {
      # Precedence: default_sim_params() -> fs_params -> per-arm overrides.
      # .modify_keep_null() preserves an explicit n.min = NULL so the
      # sample-size-adaptive rule reaches forestsearch().
      fs_p <- .modify_keep_null(
        .modify_keep_null(fs_defaults, fs_params),
        method_specs[[arm_label]]
      )
      arm_out <- .run_fs_analysis_gen(
        data = sim_data, confounders_name = confounders_name,
        params = fs_p, dgm = dgm,
        cox_formula = cox_formula, cox_formula_adj = cox_formula_adj,
        outcome_name = outcome_name, event_name = event_name,
        treat_name = treat_name, harm_col = harm_col,
        label = arm_label, verbose = show_verbose,
        keep = if (keep_this_sim) keep else character(0)
      )
      sg_hat_list[[arm_label]] <- attr(arm_out, "sg_hat")
      if (keep_this_sim) {
        keepers <- attr(arm_out, "keepers")
        if (!is.null(keepers)) diagnostics_list[[arm_label]] <- keepers
      }
      attr(arm_out, "sg_hat")  <- NULL
      attr(arm_out, "keepers") <- NULL
      results_list[[arm_label]] <- cbind(df_pop, arm_out)
    }
  }

  if (length(results_list) == 0) {
    warning(
      "No analyses were run. ",
      if (is.null(method_specs))
        "Check run_fs / run_fs_grf / run_grf."
      else
        "Check the 'methods' specification."
    )
    return(NULL)
  }

  # ── Pairwise classification concordance ──────────────────────────────────
  # When multiple methods run on the same sim_data, compute per-subject

  # agreement between each pair.  Stored as columns agree.M1.M2 and
  # kappa.M1.M2 on every row (values are per-sim, identical across methods).
  valid_sg <- names(sg_hat_list)[
    vapply(sg_hat_list, function(x) !is.null(x) && length(x) > 0, logical(1))
  ]

  if (length(valid_sg) >= 2) {
    pairs <- utils::combn(valid_sg, 2, simplify = FALSE)
    for (pair in pairs) {
      sg1 <- sg_hat_list[[pair[1]]]
      sg2 <- sg_hat_list[[pair[2]]]

      if (length(sg1) == length(sg2)) {
        agree <- mean(sg1 == sg2)

        # Cohen's kappa
        ct <- table(factor(sg1, levels = 0:1),
                    factor(sg2, levels = 0:1))
        p_obs <- sum(diag(ct)) / sum(ct)
        p_exp <- sum(rowSums(ct) * colSums(ct)) / sum(ct)^2
        kappa <- if (abs(1 - p_exp) < 1e-10) 1.0 else
          (p_obs - p_exp) / (1 - p_exp)

        col_agree <- paste0("agree.", pair[1], ".", pair[2])
        col_kappa <- paste0("kappa.", pair[1], ".", pair[2])

        for (nm in names(results_list)) {
          results_list[[nm]][[col_agree]] <- agree
          results_list[[nm]][[col_kappa]] <- kappa
        }
      }
    }
  }

  result <- data.table::rbindlist(results_list, fill = TRUE)

  # Attach heavy diagnostics (only populated when keep_this_sim = TRUE).
  # Consumers can retrieve via attr(result, "diagnostics").  Note:
  # data.table::rbindlist() may drop this attribute in downstream
  # collect_results() unless the caller preserves per-sim structure;
  # see collect_results(..., keep_diagnostics = TRUE).
  if (length(diagnostics_list) > 0L) {
    # Tag with sim_id for provenance when collected across replicates
    attr(diagnostics_list, "sim_id") <- sim_id
    attr(result, "diagnostics") <- diagnostics_list
  }

  if (show_verbose)
    message(sprintf("\n[4] Done: %d rows x %d cols\n%s",
                    nrow(result), ncol(result),
                    paste(rep("=", 60), collapse = "")))

  result
}


# =============================================================================
# Internal: ForestSearch analysis
# =============================================================================

#' @keywords internal
.run_fs_analysis_gen <- function(
    data, confounders_name, params, dgm,
    cox_formula, cox_formula_adj,
    outcome_name, event_name, treat_name, harm_col,
    label, verbose,
    keep = character(0)
) {
  if (verbose) {
    message(sprintf("\n  [%s] ForestSearch: n=%d events=%d (%.1f%%)",
                    label, nrow(data),
                    sum(data[[event_name]]), 100 * mean(data[[event_name]])))
  }

  fs_args <- list(
    df.analysis      = data,
    confounders.name = confounders_name,
    details          = verbose,
    plot.sg          = FALSE
  )

  # Parameters we ALWAYS manage internally — users cannot override
  # these via fs_params even if they try.  The rest flows from
  # formals(forestsearch) via filter_valid_args(), so any new
  # forestsearch() argument is automatically reachable.
  fs_internal_args <- c(
    "df.analysis",      # we pass `data` explicitly
    "confounders.name", # we pass `confounders_name` explicitly
    "details",          # driven by verbose argument above
    "plot.sg",          # always FALSE in simulation loops
    "quiet"             # set below to TRUE
  )

  user_fs_args <- filter_valid_args(
    params        = params,
    target_fun    = forestsearch,
    exclude       = fs_internal_args,
    warn_unknown  = FALSE  # silent drop matches legacy behavior
  )
  # NULL-preserving merge (not modifyList): an explicit n.min = NULL in
  # `params` must survive to forestsearch() so the sample-size-adaptive rule
  # fires; modifyList() would strip it and silently revert to the formal
  # default.  Affects both the legacy boolean path and the methods path.
  fs_args <- .modify_keep_null(fs_args, user_fs_args)

  fs_args$quiet <- TRUE
  fs_result <- tryCatch(
    do.call(forestsearch, fs_args),
    error = function(e) { warning(label, " failed: ", e$message); NULL }
  )

  has_result <- !is.null(fs_result) &&
    !is.null(fs_result$grp.consistency$out_sg$result) &&
    nrow(fs_result$grp.consistency$out_sg$result) > 0

  # Effective names for estimation: honour fs_params overrides.
  # Critical for Poisson+offset (IRR) where fs_params sets
  # outcome.name = "event_sim" (binary count) but the function-level
  # outcome_name still holds "y_sim" (survival time).  Without this
  # sync, the Poisson GLM receives continuous times as the response,
  # producing dpois() non-integer warnings and biased estimates.
  est_outcome <- params[["outcome.name"]] %||% outcome_name
  est_event   <- params[["event.name"]]   %||% event_name

  out <- .extract_fs_estimates_gen(
    df           = data,
    fs_res       = if (has_result) fs_result$grp.consistency$out_sg$result else NULL,
    fs_full      = if (has_result) fs_result else NULL,
    dgm          = dgm,
    cox_formula  = cox_formula,
    cox_formula_adj = cox_formula_adj,
    outcome_name = est_outcome,
    event_name   = est_event,
    treat_name   = treat_name,
    harm_col     = harm_col,
    analysis     = label,
    verbose      = verbose,
    outcome_type   = params$outcome_type,
    effect_measure = params$effect_measure,
    offset_name    = params$offset.name
  )

  # ── Populate FS-specific scalar diagnostics ──────────────────────────────
  # These are always attached; callers that don't care can ignore them.
  # Read safely — any path may be NULL when fs_result is NULL or partial.
  .sg_harm <- if (has_result) fs_result$sg.harm else NULL
  .gc      <- if (!is.null(fs_result)) fs_result$grp.consistency else NULL

  out$sg.def <- if (!is.null(.sg_harm) && length(.sg_harm) > 0L &&
                    !all(is.na(.sg_harm))) {
    paste(.sg_harm, collapse = " & ")
  } else {
    ""
  }
  out$sg.n_factors <- if (!is.null(.sg_harm))
    sum(!is.na(.sg_harm)) else 0L
  out$n_candidates_evaluated <- if (!is.null(.gc$n_candidates_evaluated))
    as.integer(.gc$n_candidates_evaluated) else NA_integer_
  out$n_candidates_total     <- if (!is.null(.gc$n_candidates_total))
    as.integer(.gc$n_candidates_total) else NA_integer_
  out$n_passed               <- if (!is.null(.gc$n_passed))
    as.integer(.gc$n_passed) else NA_integer_
  out$consistency_algorithm  <- if (!is.null(.gc$algorithm))
    as.character(.gc$algorithm) else NA_character_
  out$early_stop_triggered   <- if (!is.null(.gc$early_stop_triggered))
    as.logical(.gc$early_stop_triggered) else NA
  out$fs_minutes             <- if (!is.null(fs_result$minutes_all))
    as.numeric(fs_result$minutes_all) else NA_real_

  # ── Assemble keepers (heavy objects) when requested ─────────────────────
  keepers <- NULL
  if (length(keep) > 0L && !is.null(fs_result)) {
    keepers <- list()
    if ("fs_full" %in% keep) {
      # Strip df.est (potentially large analysis data frame) to control
      # memory; users can request it explicitly via keep = "df.est".
      fs_keep <- fs_result
      fs_keep$df.est     <- NULL
      fs_keep$df.predict <- NULL
      fs_keep$df.test    <- NULL
      keepers$fs_full <- fs_keep
    }
    if ("df.est" %in% keep && !is.null(fs_result$df.est)) {
      keepers$df.est <- fs_result$df.est
    }
    if ("candidate_table" %in% keep && has_result) {
      keepers$candidate_table <- fs_result$grp.consistency$out_sg$result
    }
  }

  # sg_hat is carried by the caller for concordance calculation;
  # keepers is populated (possibly NULL) for the diagnostics accumulator.
  attr(out, "keepers") <- keepers

  out
}


# =============================================================================
# Internal: extract FS estimates  (column-name-general)
# =============================================================================

#' @keywords internal
.extract_fs_estimates_gen <- function(
    df, fs_res, fs_full, dgm,
    cox_formula, cox_formula_adj,
    outcome_name, event_name, treat_name, harm_col,
    analysis, verbose,
    outcome_type = NULL, effect_measure = NULL, offset_name = NULL
) {

  # ── Detect GLM mode ──────────────────────────────────────────────────────
  is_glm <- !is.null(outcome_type) && outcome_type != "survival"

  out <- data.table::data.table(
    analysis    = analysis,
    any.H       = 0L,
    size.H      = NA_integer_,
    size.Hc     = nrow(df),
    hr.H.true   = NA_real_, hr.H.hat   = NA_real_,
    hr.Hc.true  = NA_real_, hr.Hc.hat  = NA_real_,
    ahr.H.true  = NA_real_, ahr.Hc.true = NA_real_,
    ahr.H.hat   = NA_real_, ahr.Hc.hat  = NA_real_,
    cde.H.true  = NA_real_, cde.Hc.true = NA_real_,
    cde.H.hat   = NA_real_, cde.Hc.hat  = NA_real_,
    hr.itt      = NA_real_,
    hr.adj.itt  = NA_real_,
    ppv = NA_real_, npv = NA_real_, sens = NA_real_, spec = NA_real_,
    p.cens = 1 - mean(df[[event_name]]),
    taumax = max(df[[outcome_name]])
  )

  # ── Build estimation closures ─────────────────────────────────────────────
  if (is_glm) {
    # Resolve GLM family from effect_measure
    glm_family <- switch(
      effect_measure %||% "OR",
      "OR"  = , "RD"  = stats::binomial(),
      "IRR" = , "IRD" = stats::poisson(),
      "MD"  =          stats::gaussian(),
      stats::binomial()
    )
    # Whether the effect is on log scale (exponentiate) or identity
    log_scale <- effect_measure %in% c("OR", "IRR", "RR")

    .glm_effect <- function(sub_df) {
      tryCatch({
        # GLM response is outcome_name — the column that matches the GLM
        # family.  For native GLM sims this is "y_sim"; for Poisson+offset
        # on survival data the caller passes "event_sim" (binary indicator)
        # so the response is a 0/1 count with log(time) as offset.
        if (!is.null(offset_name) && offset_name %in% names(sub_df)) {
          off_vec <- log(pmax(sub_df[[offset_name]], 1e-6))
          fit <- stats::glm(sub_df[[outcome_name]] ~ sub_df[[treat_name]],
                            family = glm_family, offset = off_vec)
        } else {
          fit <- stats::glm(sub_df[[outcome_name]] ~ sub_df[[treat_name]],
                            family = glm_family)
        }
        b <- coef(fit)[[2]]  # treatment coefficient
        if (log_scale) exp(b) else b
      }, error = function(e) NA_real_)
    }

    # ITT
    out$hr.itt <- .glm_effect(df)
  } else {
    # ── Survival (Cox) ITT ─────────────────────────────────────────────────
    out$hr.itt <- tryCatch(
      exp(survival::coxph(
        survival::Surv(df[[outcome_name]], df[[event_name]]) ~ df[[treat_name]]
      )$coefficients),
      error = function(e) NA_real_
    )
  }

  # Adjusted ITT (Cox only; skip for GLM)
  if (!is_glm && !is.null(cox_formula_adj)) {
    out$hr.adj.itt <- tryCatch(
      exp(survival::coxph(cox_formula_adj, data = df)$coefficients[treat_name]),
      error = function(e) NA_real_
    )
  }

  if (is.null(fs_res) || nrow(fs_res) == 0) {
    out$hr.Hc.hat <- out$hr.itt
    out$size.H  <- 0L
    out$size.Hc <- nrow(df)

    # Classification when nothing found: everyone is predicted Hc.
    # Matches fs.estimates.out() from the paper's original code.
    if (harm_col %in% names(df)) {
      true_H <- df[[harm_col]] == 1L
      n_true_H <- sum(true_H)

      if (n_true_H > 0) {
        # Alt: H exists but was not found.
        # sens = 0 (none of H captured), ppv = 0 (no Ĥ predicted).
        out$sens <- 0
        out$ppv  <- 0
        out$spec <- sum(!true_H) / nrow(df)  # = |Hc|/N
        out$npv  <- 1                          # all Hc correctly in Ĥc
      }
      # Null (n_true_H == 0): sens/ppv remain NA (undefined).
      # spec = 1, npv = 1 (everyone is Hc and predicted Hc).
      if (n_true_H == 0) {
        out$spec <- 1
        out$npv  <- 1
      }
    }

    attr(out, "sg_hat") <- rep(0L, nrow(df))
    return(out)
  }

  out$any.H <- 1L

  # Build subgroup indicator
  df_sg <- NULL
  if (!is.null(fs_full) && !is.null(fs_full$df.est) &&
      "treat.recommend" %in% names(fs_full$df.est)) {
    df$sg_hat <- as.integer(fs_full$df.est$treat.recommend == 0)
  } else {
    sg_factors <- character(0)
    for (i in seq_len(min(7, ncol(fs_res)))) {
      cn <- paste0("M.", i)
      if (cn %in% names(fs_res) && !is.na(fs_res[[cn]][1]))
        sg_factors <- c(sg_factors, fs_res[[cn]][1])
    }
    if (length(sg_factors) > 0)
      df$sg_hat <- create_subgroup_indicator(df, sg_factors)
    else {
      out$hr.Hc.hat <- out$hr.itt
      out$size.H  <- 0L
      out$size.Hc <- nrow(df)

      if (harm_col %in% names(df)) {
        true_H <- df[[harm_col]] == 1L
        n_true_H <- sum(true_H)
        if (n_true_H > 0) {
          out$sens <- 0;  out$ppv <- 0
          out$spec <- sum(!true_H) / nrow(df);  out$npv <- 1
        } else {
          out$spec <- 1;  out$npv <- 1
        }
      }
      attr(out, "sg_hat") <- rep(0L, nrow(df))
      return(out)
    }
  }

  out$size.H  <- sum(df$sg_hat, na.rm = TRUE)
  out$size.Hc <- sum(df$sg_hat == 0L, na.rm = TRUE)

  # ── Subgroup effect estimates ────────────────────────────────────────────
  if (is_glm) {
    if (out$size.H  > 10) out$hr.H.hat  <- .glm_effect(df[df$sg_hat == 1, ])
    if (out$size.Hc > 10) out$hr.Hc.hat <- .glm_effect(df[df$sg_hat == 0, ])
  } else {
    .cox_hr <- function(sub_df) tryCatch(
      exp(survival::coxph(
        survival::Surv(sub_df[[outcome_name]], sub_df[[event_name]]) ~
          sub_df[[treat_name]]
      )$coefficients),
      error = function(e) NA_real_
    )

    if (out$size.H  > 10) out$hr.H.hat  <- .cox_hr(df[df$sg_hat == 1, ])
    if (out$size.Hc > 10) out$hr.Hc.hat <- .cox_hr(df[df$sg_hat == 0, ])
  }

  # ── AHR and CDE (survival) / AOR (GLM) ───────────────────────────────────
  if (!is_glm) {
    if ("loghr_po" %in% names(df)) {
      out$ahr.H.hat  <- compute_ahr(df, df$sg_hat)
      out$ahr.Hc.hat <- compute_ahr(df, 1L - df$sg_hat)
    }
    if (all(c("theta_0", "theta_1") %in% names(df))) {
      out$cde.H.hat  <- compute_cde(df, df$sg_hat)
      out$cde.Hc.hat <- compute_cde(df, 1L - df$sg_hat)
    }
  } else {
    # GLM: marginal causal effect from potential outcomes (p0/p1 or mu0/mu1)
    em <- effect_measure %||% "OR"
    if (any(c("p0", "mu0") %in% names(df))) {
      out$ahr.H.hat  <- compute_aor(df, df$sg_hat, em)
      out$ahr.Hc.hat <- compute_aor(df, 1L - df$sg_hat, em)
      # CDE analogue (binary OR only)
      out$cde.H.hat  <- compute_cde_glm(df, df$sg_hat, em)
      out$cde.Hc.hat <- compute_cde_glm(df, 1L - df$sg_hat, em)
    }
  }

  # ── True-subgroup classification metrics ─────────────────────────────────
  if (harm_col %in% names(df)) {
    true_H <- df[[harm_col]] == 1L
    hat_H  <- df$sg_hat == 1L

    tp <- sum(true_H & hat_H);   fp <- sum(!true_H & hat_H)
    tn <- sum(!true_H & !hat_H); fn <- sum(true_H & !hat_H)

    out$sens <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
    out$spec <- if ((tn + fp) > 0) tn / (tn + fp) else NA_real_
    out$ppv  <- if ((tp + fp) > 0) tp / (tp + fp) else NA_real_
    out$npv  <- if ((tn + fn) > 0) tn / (tn + fn) else NA_real_

    # True-subgroup effect estimates
    if (is_glm) {
      if (sum(true_H) > 10)  out$hr.H.true  <- .glm_effect(df[true_H, ])
      if (sum(!true_H) > 10) out$hr.Hc.true <- .glm_effect(df[!true_H, ])

      # Marginal causal effect in true subgroup (GLM analogue of AHR)
      em <- effect_measure %||% "OR"
      if (any(c("p0", "mu0") %in% names(df))) {
        out$ahr.H.true  <- compute_aor(df, as.integer(true_H), em)
        out$ahr.Hc.true <- compute_aor(df, as.integer(!true_H), em)
        # CDE analogue (binary OR only)
        out$cde.H.true  <- compute_cde_glm(df, as.integer(true_H), em)
        out$cde.Hc.true <- compute_cde_glm(df, as.integer(!true_H), em)
      }
    } else {
      if (sum(true_H) > 10)  out$hr.H.true  <- .cox_hr(df[true_H, ])
      if (sum(!true_H) > 10) out$hr.Hc.true <- .cox_hr(df[!true_H, ])

      if ("loghr_po" %in% names(df)) {
        out$ahr.H.true  <- compute_ahr(df, df[[harm_col]])
        out$ahr.Hc.true <- compute_ahr(df, 1L - df[[harm_col]])
      }
      if (all(c("theta_0", "theta_1") %in% names(df))) {
        out$cde.H.true  <- compute_cde(df, df[[harm_col]])
        out$cde.Hc.true <- compute_cde(df, 1L - df[[harm_col]])
      }
    }
  }

  attr(out, "sg_hat") <- if ("sg_hat" %in% names(df)) df$sg_hat else
    rep(0L, nrow(df))
  out
}


# =============================================================================
# Internal: GRF analysis
# =============================================================================

#' @keywords internal
.run_grf_analysis_gen <- function(
    data, confounders_name, params, dgm,
    cox_formula, cox_formula_adj,
    outcome_name, event_name, treat_name, harm_col,
    label, verbose, debug,
    outcome_type = NULL, effect_measure = NULL, offset_name = NULL,
    id_name = "id",
    keep = character(0)
) {
  # ── GRF function dispatch ────────────────────────────────────────────────

  # GLM outcomes route to grf.subg.harm.glm(); survival (or NULL) routes
  # to the original grf.subg.harm.survival().  When outcome_type is NULL
  # (survival scenarios that don't set it), the survival path runs exactly
  # as before — no behavior change.
  is_glm <- !is.null(outcome_type) && outcome_type != "survival"
  grf_fn_name <- if (is_glm) "grf.subg.harm.glm" else "grf.subg.harm.survival"

  grf_fun <- tryCatch(
    get(grf_fn_name, mode = "function",
        envir = asNamespace("forestsearch")),
    error = function(e) NULL
  )

  if (is.null(grf_fun)) {
    warning(grf_fn_name, " not found. Skipping GRF analysis.")
    out <- .extract_grf_estimates_gen(
      df = data, grf_est = NULL, dgm = dgm,
      cox_formula = cox_formula, cox_formula_adj = cox_formula_adj,
      outcome_name = outcome_name, event_name = event_name,
      treat_name = treat_name, harm_col = harm_col,
      analysis = label, verbose = verbose, debug = debug,
      outcome_type = outcome_type, effect_measure = effect_measure,
      offset_name = offset_name, id_name = id_name
    )
    out$sg.def <- ""
    out$sg.n_factors <- 0L
    out$grf_selected_depth <- NA_integer_
    attr(out, "keepers") <- NULL
    return(out)
  }

  grf_args <- list(data = data, confounders.name = confounders_name,
                   details = verbose)

  # Parameters always managed internally — users cannot override these
  # via grf_params.  The rest flows from formals(grf_fun), so any new
  # grf.subg.harm.*() argument is automatically reachable.  grf_fun
  # was dispatched above to either grf.subg.harm.survival (is_glm FALSE)
  # or grf.subg.harm.glm (is_glm TRUE), so the allowlist tracks the
  # correct target signature automatically.
  grf_internal_args <- c(
    "data",
    "confounders.name",
    "details"
  )

  user_grf_args <- filter_valid_args(
    params        = params,
    target_fun    = grf_fun,
    exclude       = grf_internal_args,
    warn_unknown  = FALSE
  )
  grf_args <- utils::modifyList(grf_args, user_grf_args)

  # GLM-specific parameters (ignored by grf.subg.harm.survival)
  if (is_glm) {
    grf_args$outcome_type  <- outcome_type
    if (!is.null(effect_measure)) grf_args$effect_measure <- effect_measure
    if (!is.null(offset_name))    grf_args$offset.name    <- offset_name
    # grf.subg.harm.glm() does NOT accept event.name — remove it
    # to prevent "unused argument" error from do.call()
    grf_args$event.name <- NULL

    # Match forestsearch()'s adverse_outcome default: TRUE for binary/count.
    # This ensures standalone GRF uses the same Y-flip as FS's internal
    # GRF pre-screening, producing consistent subgroup identification.
    ao <- params$adverse_outcome
    if (is.null(ao)) ao <- outcome_type %in% c("binary", "count")
    grf_args$adverse_outcome <- ao

    # Pass tune_grf if set (default FALSE preserves existing behavior)
    tg <- params$tune_grf
    if (!is.null(tg)) grf_args$tune_grf <- tg
  }

  grf_result <- tryCatch(
    do.call(grf_fun, grf_args),
    error = function(e) { warning(label, " failed: ", e$message); NULL }
  )

  out <- .extract_grf_estimates_gen(
    df = data, grf_est = grf_result, dgm = dgm,
    cox_formula = cox_formula, cox_formula_adj = cox_formula_adj,
    outcome_name = outcome_name, event_name = event_name,
    treat_name = treat_name, harm_col = harm_col,
    analysis = label, verbose = verbose, debug = debug,
    outcome_type = outcome_type, effect_measure = effect_measure,
    offset_name = offset_name, id_name = id_name
  )

  # ── Populate GRF-specific scalar diagnostics ─────────────────────────────
  # For GRF results, sg.harm.id carries the cut-expression vector
  # (character).  See docs in grf_subg_harm_glm.R / grf_main.R.
  .sg_cuts <- if (!is.null(grf_result))
    grf_result$sg.harm.id %||% grf_result$sg_harm_id else NULL

  out$sg.def <- if (!is.null(.sg_cuts) && length(.sg_cuts) > 0L &&
                    !all(is.na(.sg_cuts)) && !all(.sg_cuts == "")) {
    paste(.sg_cuts, collapse = " & ")
  } else {
    ""
  }
  out$sg.n_factors <- if (!is.null(.sg_cuts))
    sum(!is.na(.sg_cuts) & .sg_cuts != "") else 0L
  out$grf_selected_depth <- if (!is.null(grf_result$selected_depth))
    as.integer(grf_result$selected_depth) else NA_integer_

  # ── Assemble GRF keepers when requested ─────────────────────────────────
  keepers <- NULL
  if (length(keep) > 0L && !is.null(grf_result) && "grf_full" %in% keep) {
    keepers <- list(grf_full = grf_result)
  }
  attr(out, "keepers") <- keepers

  out
}


# =============================================================================
# Internal: extract GRF estimates  (column-name-general)
# =============================================================================

#' @keywords internal
.extract_grf_estimates_gen <- function(
    df, grf_est, dgm,
    cox_formula, cox_formula_adj,
    outcome_name, event_name, treat_name, harm_col,
    analysis, verbose, debug,
    outcome_type = NULL, effect_measure = NULL, offset_name = NULL,
    id_name = "id"
) {

  # ── Detect GLM mode ──────────────────────────────────────────────────────
  is_glm <- !is.null(outcome_type) && outcome_type != "survival"

  out <- data.table::data.table(
    analysis    = analysis,
    any.H       = 0L,
    size.H      = NA_integer_,
    size.Hc     = nrow(df),
    hr.H.true   = NA_real_, hr.H.hat   = NA_real_,
    hr.Hc.true  = NA_real_, hr.Hc.hat  = NA_real_,
    ahr.H.true  = NA_real_, ahr.Hc.true = NA_real_,
    ahr.H.hat   = NA_real_, ahr.Hc.hat  = NA_real_,
    cde.H.true  = NA_real_, cde.Hc.true = NA_real_,
    cde.H.hat   = NA_real_, cde.Hc.hat  = NA_real_,
    hr.itt      = NA_real_,
    hr.adj.itt  = NA_real_,
    ppv = NA_real_, npv = NA_real_, sens = NA_real_, spec = NA_real_,
    p.cens = 1 - mean(df[[event_name]]),
    taumax = max(df[[outcome_name]])
  )

  # ── Build estimation closures ─────────────────────────────────────────────
  if (is_glm) {
    glm_family <- switch(
      effect_measure %||% "OR",
      "OR"  = , "RD"  = stats::binomial(),
      "IRR" = , "IRD" = stats::poisson(),
      "MD"  =          stats::gaussian(),
      stats::binomial()
    )
    log_scale <- effect_measure %in% c("OR", "IRR", "RR")

    .glm_effect_grf <- function(sub_df, oc, ec, tc) {
      tryCatch({
        # GLM response is oc (outcome_name).  For Poisson+offset on
        # survival data, the caller sets oc = "event_sim" (binary count)
        # with offset = log(time); for native GLM sims, oc = "y_sim".
        if (!is.null(offset_name) && offset_name %in% names(sub_df)) {
          off_vec <- log(pmax(sub_df[[offset_name]], 1e-6))
          fit <- stats::glm(sub_df[[oc]] ~ sub_df[[tc]],
                            family = glm_family, offset = off_vec)
        } else {
          fit <- stats::glm(sub_df[[oc]] ~ sub_df[[tc]],
                            family = glm_family)
        }
        b <- coef(fit)[[2]]
        if (log_scale) exp(b) else b
      }, error = function(e) NA_real_)
    }

    out$hr.itt <- .glm_effect_grf(df, outcome_name, event_name, treat_name)
  } else {
    out$hr.itt <- tryCatch(
      exp(survival::coxph(
        survival::Surv(df[[outcome_name]], df[[event_name]]) ~ df[[treat_name]]
      )$coefficients),
      error = function(e) NA_real_
    )
  }

  if (!is_glm && !is.null(cox_formula_adj)) {
    out$hr.adj.itt <- tryCatch(
      exp(survival::coxph(cox_formula_adj, data = df)$coefficients[treat_name]),
      error = function(e) NA_real_
    )
  }

  if (is.null(grf_est)) {
    out$hr.Hc.hat <- out$hr.itt
    out$size.H  <- 0L
    out$size.Hc <- nrow(df)

    # Classification when GRF fails: everyone is predicted Hc.
    if (harm_col %in% names(df)) {
      true_H <- df[[harm_col]] == 1L
      n_true_H <- sum(true_H)
      if (n_true_H > 0) {
        out$sens <- 0;  out$ppv <- 0
        out$spec <- sum(!true_H) / nrow(df);  out$npv <- 1
      } else {
        out$spec <- 1;  out$npv <- 1
      }
    }
    attr(out, "sg_hat") <- rep(0L, nrow(df))
    return(out)
  }

  # GRF subgroup detection
  sg_harm_id <- grf_est$sg.harm.id %||% grf_est$sg_harm_id
  has_sg     <- !is.null(sg_harm_id) && length(sg_harm_id) > 0 &&
                !all(is.na(sg_harm_id))
  has_tr     <- !is.null(grf_est$data) &&
                "treat.recommend" %in% names(grf_est$data)

  if (!has_sg || !has_tr) {
    out$hr.Hc.hat <- out$hr.itt
    out$size.H  <- 0L
    out$size.Hc <- nrow(df)

    if (harm_col %in% names(df)) {
      true_H <- df[[harm_col]] == 1L
      n_true_H <- sum(true_H)
      if (n_true_H > 0) {
        out$sens <- 0;  out$ppv <- 0
        out$spec <- sum(!true_H) / nrow(df);  out$npv <- 1
      } else {
        out$spec <- 1;  out$npv <- 1
      }
    }
    attr(out, "sg_hat") <- rep(0L, nrow(df))
    return(out)
  }

  out$any.H <- 1L
  grf_data  <- grf_est$data
  harm_ind  <- as.integer(grf_data$treat.recommend == 0)
  grf_data$sg_hat <- harm_ind

  out$size.H  <- sum(harm_ind)
  out$size.Hc <- sum(harm_ind == 0L)

  # Detect column names in grf_data (may differ from sim_data names)
  o_col <- if (outcome_name %in% names(grf_data)) outcome_name else
           if ("y.sim" %in% names(grf_data)) "y.sim" else outcome_name
  e_col <- if (event_name %in% names(grf_data)) event_name else
           if ("event.sim" %in% names(grf_data)) "event.sim" else event_name
  t_col <- if (treat_name %in% names(grf_data)) treat_name else
           if ("treat" %in% names(grf_data)) "treat" else treat_name

  if (is_glm) {
    if (out$size.H  > 10) out$hr.H.hat  <-
      .glm_effect_grf(grf_data[grf_data$sg_hat == 1, ], o_col, e_col, t_col)
    if (out$size.Hc > 10) out$hr.Hc.hat <-
      .glm_effect_grf(grf_data[grf_data$sg_hat == 0, ], o_col, e_col, t_col)
  } else {
    .cox_hr <- function(sub_df, oc, ec, tc) tryCatch(
      exp(survival::coxph(
        survival::Surv(sub_df[[oc]], sub_df[[ec]]) ~ sub_df[[tc]]
      )$coefficients),
      error = function(e) NA_real_
    )

    if (out$size.H  > 10) out$hr.H.hat  <-
      .cox_hr(grf_data[grf_data$sg_hat == 1, ], o_col, e_col, t_col)
    if (out$size.Hc > 10) out$hr.Hc.hat <-
      .cox_hr(grf_data[grf_data$sg_hat == 0, ], o_col, e_col, t_col)
  }

  # AHR / CDE (survival) / AOR (GLM)
  if (!is_glm) {
    .merge_sg <- function(base_df, sg_df, cols) {
      if (id_name %in% names(base_df) && id_name %in% names(sg_df))
        merge(base_df[, c(id_name, cols), drop = FALSE],
              sg_df[, c(id_name, "sg_hat")], by = id_name, all.y = TRUE)
      else if (nrow(base_df) == nrow(sg_df)) {
        base_df$sg_hat <- sg_df$sg_hat; base_df
      } else NULL
    }

    if ("loghr_po" %in% names(df)) {
      m <- .merge_sg(df, grf_data, "loghr_po")
      if (!is.null(m)) {
        out$ahr.H.hat  <- compute_ahr(m, m$sg_hat)
        out$ahr.Hc.hat <- compute_ahr(m, 1L - m$sg_hat)
      }
    }
    if (all(c("theta_0", "theta_1") %in% names(df))) {
      m <- .merge_sg(df, grf_data, c("theta_0", "theta_1"))
      if (!is.null(m)) {
        out$cde.H.hat  <- compute_cde(m, m$sg_hat)
        out$cde.Hc.hat <- compute_cde(m, 1L - m$sg_hat)
      }
    }
  } else {
    # GLM: marginal causal effect from potential outcomes
    po_cols <- intersect(c("p0", "p1", "mu0", "mu1"), names(df))
    if (length(po_cols) >= 2) {
      .merge_sg_po <- function(base_df, sg_df, cols) {
        if (id_name %in% names(base_df) && id_name %in% names(sg_df))
          merge(base_df[, c(id_name, cols), drop = FALSE],
                sg_df[, c(id_name, "sg_hat")], by = id_name, all.y = TRUE)
        else if (nrow(base_df) == nrow(sg_df)) {
          base_df$sg_hat <- sg_df$sg_hat; base_df
        } else NULL
      }
      m <- .merge_sg_po(df, grf_data, po_cols)
      if (!is.null(m)) {
        em <- effect_measure %||% "OR"
        out$ahr.H.hat  <- compute_aor(m, m$sg_hat, em)
        out$ahr.Hc.hat <- compute_aor(m, 1L - m$sg_hat, em)
        # CDE analogue (binary OR only)
        out$cde.H.hat  <- compute_cde_glm(m, m$sg_hat, em)
        out$cde.Hc.hat <- compute_cde_glm(m, 1L - m$sg_hat, em)
      }
    }
  }

  # True-subgroup metrics
  if (harm_col %in% names(df)) {
    .merge_sg_class <- function(base_df, sg_df, hcol) {
      if (id_name %in% names(base_df) && id_name %in% names(sg_df))
        merge(base_df[, c(id_name, hcol), drop = FALSE],
              sg_df[, c(id_name, "sg_hat")], by = id_name, all.y = TRUE)
      else if (nrow(base_df) == nrow(sg_df)) {
        base_df$sg_hat <- sg_df$sg_hat; base_df
      } else NULL
    }

    m_class <- .merge_sg_class(df, grf_data, harm_col)
    if (!is.null(m_class)) {
      true_H <- m_class[[harm_col]] == 1L
      hat_H  <- m_class$sg_hat == 1L

      tp <- sum(true_H & hat_H);   fp <- sum(!true_H & hat_H)
      tn <- sum(!true_H & !hat_H); fn <- sum(true_H & !hat_H)

      out$sens <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
      out$spec <- if ((tn + fp) > 0) tn / (tn + fp) else NA_real_
      out$ppv  <- if ((tp + fp) > 0) tp / (tp + fp) else NA_real_
      out$npv  <- if ((tn + fn) > 0) tn / (tn + fn) else NA_real_
    }

    true_H_df <- df[[harm_col]] == 1L
    if (is_glm) {
      if (sum(true_H_df) > 10)  out$hr.H.true  <-
        .glm_effect_grf(df[true_H_df, ], outcome_name, event_name, treat_name)
      if (sum(!true_H_df) > 10) out$hr.Hc.true <-
        .glm_effect_grf(df[!true_H_df, ], outcome_name, event_name, treat_name)

      # Marginal causal effect in true subgroup (GLM analogue of AHR)
      em <- effect_measure %||% "OR"
      if (any(c("p0", "mu0") %in% names(df))) {
        out$ahr.H.true  <- compute_aor(df, as.integer(true_H_df), em)
        out$ahr.Hc.true <- compute_aor(df, as.integer(!true_H_df), em)
        # CDE analogue (binary OR only)
        out$cde.H.true  <- compute_cde_glm(df, as.integer(true_H_df), em)
        out$cde.Hc.true <- compute_cde_glm(df, as.integer(!true_H_df), em)
      }
    } else {
      if (sum(true_H_df) > 10)  out$hr.H.true  <-
        .cox_hr(df[true_H_df, ], outcome_name, event_name, treat_name)
      if (sum(!true_H_df) > 10) out$hr.Hc.true <-
        .cox_hr(df[!true_H_df, ], outcome_name, event_name, treat_name)

      if ("loghr_po" %in% names(df)) {
        out$ahr.H.true  <- compute_ahr(df, df[[harm_col]])
        out$ahr.Hc.true <- compute_ahr(df, 1L - df[[harm_col]])
      }
      if (all(c("theta_0", "theta_1") %in% names(df))) {
        out$cde.H.true  <- compute_cde(df, df[[harm_col]])
        out$cde.Hc.true <- compute_cde(df, 1L - df[[harm_col]])
      }
    }
  }

  attr(out, "sg_hat") <- grf_data$sg_hat
  out
}