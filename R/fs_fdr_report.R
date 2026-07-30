# =============================================================================
# fs_fdr_report.R
#
# Operating false-discovery characterization for a *fitted* ForestSearch
# analysis, with optional harm confirmation on the de-biased effect.
#
# This is NOT threshold calibration.  It does not solve for (c1, c2) that
# control an FDR target.  For a GIVEN analysis run at a FIXED (c1, c2), it
# answers two questions by simulation under a null:
#
#   (A) What IS the false-discovery rate?
#         A = P(declare any harm subgroup | c1, c2) under H0.
#   (B) Does requiring harm confirmation on the de-biased HR reduce it?
#         B(g) = P(declare AND de-biased HR >= g).   Always B(g) <= A.
#       The gap A - B(g) is what the de-biasing (the submitted-paper method)
#       buys at confirmation threshold g.
#
# The (c1, c2) and every identifier knob are INHERITED from the fitted object
# (fs_analysis$args_call_all), so the reported FDR is the FDR of the analysis
# actually run -- not a re-specified variant.  The data-binding arguments (the
# outcome/event/treat/id/flag_harm names and the confounder set) are rebound to
# the DGM, because simulate_from_dgm() emits fixed column names.
#
# Scope: survival / Cox analyses under the GBSG-AFT null only.  The GLM
# (binary / count / continuous) null is a separate tool.
#
# Reuses package internals from fpr_calibration.R: `%||%` and `.wilson_ci()`.
#
# Author: Larry F. Leon  |  forestsearch >= 0.2.0
# =============================================================================

# GBSG covariate set emitted by setup_gbsg_dgm(); fallback candidate set when
# the fitted analysis carries no confounders.name.
.FS_GBSG_CONFOUNDERS <- c("er", "age", "meno", "pgr", "nodes", "size", "grade")


#' Operating False-Discovery Rate of a Fitted ForestSearch Analysis
#'
#' Characterizes the operating false-discovery rate (FDR) of a **fitted**
#' \code{\link{forestsearch}} analysis and tests whether requiring harm
#' confirmation on the de-biased effect reduces it.  This is a **measurement**
#' tool, not a calibration tool:
#' it does not solve for thresholds that control an FDR target.  For the
#' \eqn{(c_1, c_2)} pair the analysis actually used, it reports two quantities
#' under a null:
#'
#' - **A** (raw FDR): \eqn{A = P(\text{declare any harm subgroup} \mid c_1, c_2)}.
#' - **B(g)** (confirmed FDR): \eqn{B(g) = P(\text{declare AND de-biased } HR \ge g)}.
#'   Because harm confirmation can only remove declarations, \eqn{B(g) \le A}
#'   always; the gap \eqn{A - B(g)} is what the de-biasing step buys at
#'   confirmation threshold \eqn{g}.
#'
#' All thresholds and identifier knobs are inherited from
#' \code{fs_analysis$args_call_all}, so the FDR reported is that of the analysis
#' as run.  Multiplier resampling is forced on internally
#' (\code{mr_inference = TRUE}) regardless of the fitted setting, because B
#' requires the per-replicate de-biased HR.
#'
#' @param fs_analysis A fitted \code{forestsearch()} object.  Must carry
#'   \code{args_call_all} (the stored call arguments), including
#'   \code{hr.threshold} and \code{hr.consistency}.
#' @param df_analysis The analyst's data frame.  Used only to default the
#'   simulated sample size (\code{n = nrow(df_analysis)}); the null data are
#'   generated from the GBSG-AFT DGM, not resampled from this frame.
#' @param null Character vector; which null(s) to simulate.  \code{"hom"} is the
#'   homogeneous null (\code{k_inter = 0}: a real, uniform treatment effect,
#'   no modifier); \code{"harm"} is the null-harm-region null
#'   (\code{target_hr_harm = 1.0}: a real modifier tuned so the flagged region
#'   has HR = 1.0 while treatment stays active elsewhere).  Default both.
#' @param c_confirm Numeric vector of de-biased-HR harm-confirmation
#'   thresholds at which to report B, i.e. the values \code{t_confirm} is swept
#'   over.  Default \code{c(0.9, 1.0, 1.25)}.  On real data harm confirmation
#'   is inert; here it annotates the null's operating characteristics only.
#' @param nsims Integer; null replicates per null.  Default 1000.
#' @param n Integer or \code{NULL}; simulated trial size.  Default
#'   \code{nrow(df_analysis)}.  The published reference baseline used
#'   \code{n = 700}; pass that to reproduce it exactly.
#' @param confounders Character vector or \code{NULL}; candidate covariates.
#'   Default \code{NULL} inherits \code{confounders.name} from the fitted
#'   object (intersected per replicate with the DGM columns), falling back to
#'   the GBSG covariate set when the analysis carried none.
#' @param n_super Integer; super-population size for the DGM.  Default 25000.
#' @param analysis_time Numeric; administrative censoring time (months).
#'   Default 84.
#' @param cens_adjust Numeric; censoring shift passed to
#'   \code{simulate_from_dgm()}.  Default \code{log(1.5)}.
#' @param event_min Integer; minimum events for a replicate to be scored.
#'   Replicates below this (or with a single treatment arm) are counted as
#'   visible failures, never silently dropped.  Default 5.
#' @param n_workers Integer or \code{NULL}; parallel workers for the outer
#'   replicate loop.  Default \code{NULL} uses ~90% of physical cores.  The
#'   inner \code{forestsearch()} call always runs sequentially
#'   (\code{parallel_args = list(plan = "sequential")}) to avoid nested
#'   parallelism, which breaks forestsearch.
#' @param seed_base Integer; base seed.  Each null uses a disjoint seed block
#'   and each replicate an explicit seed, so results are plan-invariant.
#'   Default 8316951.
#' @param quiet Logical; suppress progress and per-null messages.  Default
#'   \code{FALSE}.
#'
#' @return An object of class \code{"fs_fdr_report"}: a list with
#' \describe{
#'   \item{\code{fdr}}{Data frame, one row per null x \code{c_confirm}, with
#'     \code{A}, \code{A_lo}, \code{A_hi} (raw FDR and 95% Wilson CI),
#'     \code{B}, \code{B_lo}, \code{B_hi} (confirmed FDR and CI), \code{reduction}
#'     (\code{A - B}), and the accounting columns \code{n_valid},
#'     \code{n_declared}, \code{n_deb_na}, \code{fs_errors}, \code{unhandled}.}
#'   \item{\code{hr_structure}}{Data frame of the realized hazard ratios each
#'     DGM encodes -- the structural check that the nulls are different worlds.}
#'   \item{\code{verdict}}{\code{"CLEAN"} if no replicate failed across any
#'     null, else \code{"WARNING"}.}
#'   \item{\code{n_errors}}{Total failed replicates (forestsearch errors +
#'     unhandled).}
#'   \item{\code{c1}, \code{c2}, \code{n}, \code{nsims}, \code{nulls}}{The
#'     inherited thresholds and run settings.}
#'   \item{\code{call_args}}{Reproducibility metadata.}
#' }
#'
#' @details
#' ## What "false discovery" means here
#' Under \code{"hom"} the treatment effect is real but uniform, so any declared
#' subgroup is a manufactured heterogeneity -- the winner's-curse-of-search
#' rate.  Under \code{"harm"} a real modifier exists and the flagged region
#' sits exactly at HR = 1.0, so a small upward fluctuation reads as harm; this
#' is the more adversarial null.  The declared cut need not overlap
#' \code{flag_harm}: under either null any declaration is a false discovery.
#'
#' ## The confirmation accounting (why n_deb_na matters)
#' B(g) counts a declaration only when its de-biased HR is finite and clears
#' \code{g}; a declaration whose de-biased HR could not be computed counts as
#' unconfirmed.  Such cases are counted in \code{n_deb_na} and reported, so
#' that \code{A - B} is never mistaken for pure de-biasing benefit when part
#' of it is de-biased-HR unavailability.
#'
#' ## Inheritance and rebinding
#' Identifier knobs (thresholds, \code{subgroup_method}, \code{consistency_method},
#' \code{conf_force}, \code{conf.cont_jcuts}, \code{maxk}, MR arguments, ...)
#' transfer verbatim from \code{args_call_all}.  The data-binding arguments --
#' the five column names and the confounder set -- are rebound to the DGM,
#' because \code{simulate_from_dgm()} emits fixed names
#' (\code{y_sim / event_sim / treat_sim / id / flag_harm}).  \code{mr_inference}
#' is forced \code{TRUE} so B is computable.
#'
#' @examples
#' \dontrun{
#' # A fitted survival/Cox analysis (thresholds and knobs live on the object)
#' fit <- forestsearch(
#'   df.analysis = gbsg_df,
#'   outcome.name = "time", event.name = "status",
#'   treat.name = "treat", id.name = "id",
#'   confounders.name = c("er", "age", "meno", "pgr", "nodes", "size", "grade"),
#'   is.RCT = TRUE, est.scale = "hr",
#'   subgroup_method = "consistency", consistency_method = "resample",
#'   hr.threshold = 1.25, hr.consistency = 1.0,
#'   pconsistency.threshold = 0.90, maxk = 2,
#'   mr_inference = TRUE)
#'
#' # Characterize its operating FDR and whether harm confirmation helps
#' rep <- fs_fdr_report(fit, gbsg_df, nsims = 1000, n_workers = 100)
#' rep                       # prints the realized-HR table and the A/B table
#' rep$fdr                   # one row per null x c_confirm
#'
#' # Reproduce the published reference (n = 700): A ~ 0.138 (hom), 0.278 (harm)
#' fs_fdr_report(fit, gbsg_df, n = 700, nsims = 1000, n_workers = 100)
#' }
#'
#' @seealso \code{\link{forestsearch}} (and its vocabulary section),
#'   \code{\link{fs_mr_inference}} for the MR step this forces on,
#'   \code{\link{mr_estimates_table}},
#'   \code{\link{forestsearch_bootstrap_dofuture}} for the full bootstrap,
#'   \code{\link{fpr_calibration}},
#'   \code{\link{setup_gbsg_dgm}}, \code{\link{calibrate_k_inter}},
#'   \code{\link{simulate_from_dgm}}
#'
#' @importFrom foreach foreach
#' @importFrom future plan
#' @importFrom doFuture %dofuture%
#' @importFrom parallel detectCores
#' @importFrom utils modifyList head
#' @export
fs_fdr_report <- function(fs_analysis,
                          df_analysis,
                          null          = c("hom", "harm"),
                          c_confirm        = c(0.9, 1.0, 1.25),
                          nsims         = 1000L,
                          n             = NULL,
                          confounders   = NULL,
                          n_super       = 25000L,
                          analysis_time = 84,
                          cens_adjust   = log(1.5),
                          event_min     = 5L,
                          n_workers     = NULL,
                          seed_base     = 8316951L,
                          quiet         = FALSE) {

  # ---- validate the fitted object and extract the inherited configuration ----
  if (is.null(fs_analysis$args_call_all))
    stop("fs_analysis$args_call_all is missing; pass a fitted forestsearch() ",
         "object.", call. = FALSE)
  params <- fs_analysis$args_call_all
  c1 <- params$hr.threshold
  c2 <- params$hr.consistency
  if (is.null(c1) || is.null(c2))
    stop("hr.threshold / hr.consistency not found in fs_analysis$args_call_all.",
         call. = FALSE)

  # ---- scope guard: survival / Cox only (the GBSG-AFT null is survival) ----
  otype <- params$outcome_type
  if (!is.null(otype) && otype %in% c("binary", "count", "continuous"))
    stop(sprintf(paste0("fs_fdr_report() supports survival/Cox analyses only; ",
                        "the GBSG-AFT null is not valid for outcome_type = ",
                        "'%s'. The GLM null is a separate tool."), otype),
         call. = FALSE)

  null <- match.arg(null, c("hom", "harm"), several.ok = TRUE)
  if (is.null(n)) n <- nrow(df_analysis)
  if (is.null(n_workers))
    n_workers <- max(1L, ceiling(0.90 *
                                   (parallel::detectCores(logical = FALSE) - 1L)))

  # ---- confounders: analyst's set (intersected per rep with the DGM columns),
  #      else the GBSG covariate set ----
  conf_inherit <- params$confounders.name
  conf_use <- if (!is.null(confounders)) {
    confounders
  } else if (!is.null(conf_inherit) && length(conf_inherit) > 0) {
    conf_inherit
  } else {
    .FS_GBSG_CONFOUNDERS
  }

  # ---- strip calibrator-managed + data-binding args; force MR ON ----
  # Identifier knobs transfer verbatim; the five names + confounders are rebound
  # to the DGM per replicate (simulate_from_dgm() emits fixed names).
  drop <- c("df.analysis", "hr.threshold", "hr.consistency", "seedit",
            "parallel_args", "outcome.name", "event.name", "treat.name",
            "id.name", "flag_harm.name", "confounders.name")
  fs_params <- params[setdiff(names(params), drop)]
  fs_params$mr_inference <- TRUE                   # B needs the de-biased HR
  if (is.null(fs_params$mr_inference_args))
    fs_params$mr_inference_args <- list(ci_method = "ij", draws = 1000L,
                                       include_complement = TRUE)
  fs_params$quiet   <- TRUE
  fs_params$details <- FALSE
  fs_params$plot.sg <- FALSE

  # names emitted by simulate_from_dgm() (the single place these appear)
  dgm_names <- list(outcome = "y_sim", event = "event_sim",
                    treat = "treat_sim", id = "id", harm = "flag_harm")

  # ---- build the requested DGM null(s): canonical GBSG-AFT structure ----
  if (!quiet) message("Building GBSG-AFT null(s) ...")
  dgms <- list()
  if ("hom" %in% null)
    dgms[["H0-hom  (k_inter = 0)"]] <-
      setup_gbsg_dgm(model = "alt", k_inter = 0, n_super = n_super)
  if ("harm" %in% null) {
    k_harm <- calibrate_k_inter(target_hr_harm = 1.0, model = "alt",
                                use_ahr = FALSE)
    if (!quiet) message(sprintf("  H0-harm k_inter = %.4f", k_harm))
    dgms[["H0-harm (target HR_H = 1.0)"]] <-
      setup_gbsg_dgm(model = "alt", k_inter = k_harm, n_super = n_super)
  }

  # per-DGM dataset generator (fixed n / censoring; explicit per-rep seed)
  gen_for <- function(dgm) {
    force(dgm)
    function(seed) {
      d <- simulate_from_dgm(dgm, n = n, analysis_time = analysis_time,
                             cens_adjust = cens_adjust, seed = seed)
      d[[dgm_names$id]] <- seq_len(nrow(d))
      d
    }
  }

  # ---- run each null through the visible-error harness ----
  if (!quiet)
    message(sprintf("Running %d reps/null x %d null(s) at n = %d on %d workers ...",
                    nsims, length(dgms), n, n_workers))
  runs <- lapply(seq_along(dgms), function(j) {
    .fs_fdr_run_null(
      gen_df = gen_for(dgms[[j]]), label = names(dgms)[j],
      fs_params = fs_params, c1 = c1, c2 = c2, dgm_names = dgm_names,
      conf_use = conf_use, nsims = nsims, n_workers = n_workers,
      seed_base = seed_base + 100000L * j,   # disjoint seed block per null
      c_confirm = c_confirm, n = n, event_min = event_min, quiet = quiet)
  })

  fdr <- do.call(rbind, lapply(runs, `[[`, "rows"))
  rownames(fdr) <- NULL

  # ---- realized-HR structural table (proves the nulls are different worlds) --
  hr_structure <- do.call(rbind, lapply(seq_along(dgms), function(j)
    .fs_fdr_hr_row(names(dgms)[j], dgms[[j]])))
  rownames(hr_structure) <- NULL

  # ---- CLEAN / WARNING verdict (visible-error accounting) ----
  tot_err <- sum(vapply(runs, function(r) r$fs_errors + r$unhandled, integer(1)))
  verdict <- if (tot_err == 0L) "CLEAN" else "WARNING"
  if (!quiet) {
    cat(sprintf("\nfailure accounting: forestsearch errors + unhandled = %d\n",
                tot_err))
    cat(if (verdict == "CLEAN")
          "  -> CLEAN: no failures across any null; FDR estimates trustworthy.\n"
        else
          "  -> WARNING: failures present; inspect messages before trusting.\n")
  }

  structure(
    list(fdr          = fdr,
         hr_structure = hr_structure,
         verdict      = verdict,
         n_errors     = tot_err,
         c1           = c1,
         c2           = c2,
         n            = n,
         nsims        = nsims,
         nulls        = null,
         call_args    = list(c_confirm = c_confirm, n_super = n_super,
                             analysis_time = analysis_time,
                             cens_adjust = cens_adjust, seed_base = seed_base,
                             n_workers = n_workers,
                             mr_inference_args = fs_params$mr_inference_args)),
    class = "fs_fdr_report")
}


# ---------------------------------------------------------------------------
# Internal: one replicate.  Builds the forestsearch call by inheriting the
# stripped identifier knobs and rebinding the DGM names/confounders, runs it
# with BOTH errors and consistency warnings captured (never masked), and
# returns the declaration flag plus the de-biased HR.
#   flag: 1 declared, 0 clean no-subgroup, -1 failure/degenerate
#   deb_hr: HR-scale de-biased estimate of the selected subgroup (NA if none)
#   msg: failure text (NA otherwise)
# ---------------------------------------------------------------------------
.fs_fdr_one_rep <- function(df, seed, fs_params, c1, c2, dgm_names,
                            conf_use, event_min = 5L) {
  # degenerate-simulation guard, counted as a visible failure
  if (sum(df[[dgm_names$event]]) < event_min ||
      length(unique(df[[dgm_names$treat]])) < 2L)
    return(list(flag = -1L, deb_hr = NA_real_,
                msg = "degenerate simulation (< event_min events or single arm)"))

  confs <- intersect(conf_use, names(df))
  call_args <- utils::modifyList(fs_params, list(
    df.analysis      = df,
    outcome.name     = dgm_names$outcome,
    event.name       = dgm_names$event,
    treat.name       = dgm_names$treat,
    id.name          = dgm_names$id,
    flag_harm.name   = dgm_names$harm,
    confounders.name = confs,
    hr.threshold     = c1,
    hr.consistency   = c2,
    seedit           = seed,
    parallel_args    = list(plan = "sequential")))   # NO nested parallelism

  fs_err <- FALSE
  err_msg <- NA_character_
  fs <- withCallingHandlers(
    tryCatch(do.call(forestsearch, call_args),
             error = function(e) {
               fs_err <<- TRUE
               err_msg <<- conditionMessage(e)
               NULL
             }),
    warning = function(w) {
      # A consistency/subscript warning is a hidden failure of the search --
      # promote it so it is counted, not muffled into a fake non-declaration.
      if (grepl("subgroup.consistency|subscript", conditionMessage(w))) {
        fs_err <<- TRUE
        err_msg <<- conditionMessage(w)
      }
      invokeRestart("muffleWarning")
    })
  if (fs_err || is.null(fs))
    return(list(flag = -1L, deb_hr = NA_real_, msg = err_msg))

  declared <- !is.null(fs$sg.harm) && length(fs$sg.harm) > 0
  deb_hr <- NA_real_
  if (declared && !is.null(fs$mr_inference))
    # MR returns debiased$est on the effect scale (HR here, via to_eff())
    deb_hr <- tryCatch(fs$mr_inference$debiased$est, error = function(e) NA_real_)
  list(flag = as.integer(declared), deb_hr = deb_hr, msg = NA_character_)
}


# ---------------------------------------------------------------------------
# Internal: run one null.  Parallel over replicates (outer multisession, inner
# forestsearch sequential); collect flags, de-biased HRs, and failure messages;
# print a diagnostic line with distinct failure messages; score A and B(g) with
# Wilson CIs.  Returns the per-c_confirm rows plus the failure accounting.
# ---------------------------------------------------------------------------
.fs_fdr_run_null <- function(gen_df, label, fs_params, c1, c2, dgm_names,
                             conf_use, nsims, n_workers, seed_base,
                             c_confirm, n, event_min, quiet) {
  future::plan("sequential"); gc()
  future::plan("multisession", workers = n_workers)

  i <- NULL   # suppress R CMD check NOTE for the foreach iterator variable
  outs <- foreach::foreach(
    i = seq_len(nsims), .errorhandling = "pass",
    .options.future = list(packages = c("forestsearch", "survival"),
                           seed = TRUE)) %dofuture% {
    future::plan("sequential")
    d <- gen_df(seed_base + i)
    .fs_fdr_one_rep(d, seed_base + i, fs_params, c1, c2, dgm_names,
                    conf_use, event_min)
  }
  future::plan("sequential"); gc()

  get_flag <- function(o) if (is.list(o) && !is.null(o$flag)) o$flag else NA_integer_
  get_deb  <- function(o) if (is.list(o) && !is.null(o$deb_hr)) o$deb_hr else NA_real_
  get_msg  <- function(o)
    if (is.list(o) && !is.null(o$msg)) o$msg
    else if (inherits(o, "condition")) conditionMessage(o) else NA_character_

  flags <- vapply(outs, get_flag, integer(1))
  debs  <- vapply(outs, get_deb,  numeric(1))
  msgs  <- vapply(outs, get_msg,  character(1))

  n_unhdl <- sum(is.na(flags))
  n_fserr <- sum(flags == -1L, na.rm = TRUE)
  ok      <- flags %in% c(0L, 1L)
  n_valid <- sum(ok)
  n_decl  <- sum(flags[ok] == 1L)

  if (!quiet) {
    cat(sprintf("[%s] valid=%d  declared=%d  fs_errors=%d  unhandled=%d  (nsims=%d)\n",
                label, n_valid, n_decl, n_fserr, n_unhdl, nsims))
    if (n_fserr > 0 || n_unhdl > 0) {
      bad <- unique(msgs[!is.na(msgs)])
      cat("   failure messages (distinct):\n")
      for (m in utils::head(bad, 5L)) cat("     -", m, "\n")
    }
  }

  # A: declare-any, with a 95% Wilson CI
  A    <- if (n_valid > 0) n_decl / n_valid else NA_real_
  a_ci <- .wilson_ci(n_decl, n_valid)

  # B(g): declared AND de-biased HR >= g (NA de-biased HR counts as unconfirmed)
  decl_ok <- ok & flags == 1L
  rows <- lapply(c_confirm, function(g) {
    n_confirm <- sum(decl_ok & is.finite(debs) & debs >= g)
    n_dbna <- sum(decl_ok & !is.finite(debs))
    B      <- if (n_valid > 0) n_confirm / n_valid else NA_real_
    b_ci   <- .wilson_ci(n_confirm, n_valid)
    data.frame(
      null       = label,
      n          = n,
      nsims      = nsims,
      n_valid    = n_valid,
      fs_errors  = n_fserr,
      unhandled  = n_unhdl,
      n_declared = n_decl,
      A          = round(A, 4),
      A_lo       = round(a_ci[1], 4),
      A_hi       = round(a_ci[2], 4),
      c_confirm     = g,
      n_deb_na   = n_dbna,
      B          = round(B, 4),
      B_lo       = round(b_ci[1], 4),
      B_hi       = round(b_ci[2], 4),
      reduction  = round(A - B, 4),
      stringsAsFactors = FALSE)
  })

  list(rows = do.call(rbind, rows),
       n_valid = n_valid, n_declared = n_decl,
       fs_errors = n_fserr, unhandled = n_unhdl,
       msgs = unique(msgs[!is.na(msgs)]))
}


# ---------------------------------------------------------------------------
# Internal: one row of the realized-HR structural table for a DGM.
# ---------------------------------------------------------------------------
.fs_fdr_hr_row <- function(tag, dgm) {
  hr <- dgm$hazard_ratios
  # k_inter is stored under effect_modifiers by setup_gbsg_dgm() /
  # .create_gbsg_dgm_() (sim_aft_gbsg.R), not model_params.
  ki <- dgm$effect_modifiers$k_inter %||% dgm$model_params$k_inter %||%
    dgm$k_inter %||% NA_real_
  data.frame(
    null    = tag,
    k_inter = round(ki, 3),
    HR_H    = round(dgm$hr_H_true  %||% hr$HR_H  %||% NA_real_, 3),
    HR_Hc   = round(dgm$hr_Hc_true %||% hr$HR_Hc %||% NA_real_, 3),
    AHR_H   = round(dgm$AHR_H_true  %||% NA_real_, 3),
    AHR_Hc  = round(dgm$AHR_Hc_true %||% NA_real_, 3),
    stringsAsFactors = FALSE)
}


#' Print method for \code{fs_fdr_report} objects
#'
#' @param x An \code{fs_fdr_report} object.
#' @param ... Unused; present for S3 compatibility.
#' @return The input \code{x}, invisibly.
#' @export
print.fs_fdr_report <- function(x, ...) {
  cat("\n=== fs_fdr_report ===\n")
  cat(sprintf("  Inherited thresholds : c1 = %.3f, c2 = %.3f\n", x$c1, x$c2))
  cat(sprintf("  n per replicate      : %d\n", x$n))
  cat(sprintf("  Replicates per null  : %d\n", x$nsims))
  cat(sprintf("  Verdict              : %s (failed replicates: %d)\n",
              x$verdict, x$n_errors))
  cat("\n  Realized null structure (HR by region):\n")
  print(x$hr_structure, row.names = FALSE)
  cat("\n  Operating FDR (A) and harm-confirmed FDR (B), 95% Wilson CIs:\n")
  print(x$fdr, row.names = FALSE)
  cat("\n  A = P(declare any harm subgroup) under H0.\n")
  cat("  B(g) = P(declare AND de-biased HR >= g).  reduction = A - B(g).\n")
  if (x$n_errors > 0L)
    cat("  WARNING: failed replicates present; inspect run messages.\n")
  invisible(x)
}
