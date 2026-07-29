# =============================================================================
# STAGE E -- END-TO-END, MEASURED (not modelled)
# =============================================================================
# Replaces the analytic cost model with an Rprof profile of a real
# forestsearch() call, so the share of runtime attributable to each stage is
# observed rather than inferred.

.e_row <- function(id, what, value, unit, note = "") {
  data.frame(stringsAsFactors = FALSE,
             id = id, what = what, value = value, unit = unit, note = note)
}


#' Attribute a real forestsearch() run to stages via Rprof
#'
#' Sequential by design: profiling across future workers does not aggregate,
#' so a sequential profile is the only honest attribution. The wall-clock
#' number that matters is measured separately in stage P.
endtoend_profile <- function(df, fs_args, interval = 0.02) {
  pf <- tempfile(fileext = ".out")
  args <- c(list(df.analysis = df), fs_args)
  args$parallel_args <- list(plan = "sequential", workers = 1L)

  future::plan(future::sequential)
  gc()
  Rprof(pf, interval = interval, line.profiling = FALSE)
  elapsed <- tryCatch(system.time(
    res <- do.call(forestsearch::forestsearch, args))[["elapsed"]],
    error = function(e) NA_real_)
  Rprof(NULL)

  s <- tryCatch(summaryRprof(pf), error = function(e) NULL)
  if (is.null(s)) return(.e_row("E1", "profile unavailable", NA, "s"))

  tot <- s$by.total
  pick <- function(pat) {
    hits <- grep(pat, rownames(tot))
    if (!length(hits)) return(0)
    max(tot$total.pct[hits])
  }

  # Patterns corrected on-machine against a real profile (the sandbox names
  # did not exist):
  #   * consistency runs as evaluate_consistency_twostage / .run_single_-
  #     consistency_split (NOT "subgroup.consistency", which is only the thin
  #     dispatcher frame);
  #   * candidate screening runs as search_combinations_parallel /
  #     evaluate_combination_with_status / fit_cox_for_subgroup;
  #   * GRF and LASSO candidate CONSTRUCTION are separate, dominant stages that
  #     the original two patterns attributed to neither -- they are added
  #     explicitly (E3b, E3c) so the shares sum to a plausible total instead of
  #     hiding ~50% of runtime in "other".
  rows <- list(
    .e_row("E1", "forestsearch() wall-clock (sequential)", elapsed, "s",
           "profiled run; see stage P for parallel wall-clock"),
    .e_row("E2", "share in consistency evaluation",
           pick("evaluate_consistency_twostage|run_single_consistency_split|subgroup.consistency|evaluate_subgroup_consistency|consistency_via_splits"),
           "% of total", "measured; tiny here because early stop fires at candidate 1"),
    .e_row("E3", "share in subgroup search (candidate screening)",
           pick("search_combinations_parallel|search_combinations|evaluate_combination_with_status|evaluate_combination|subgroup.search"),
           "% of total", "per-candidate Cox+survfit screening during enumeration"),
    .e_row("E3b", "share in GRF candidate construction",
           pick("causal_survival_forest|causal_forest|fit_causal_forest|regression_forest|grf::"),
           "% of total", "candidate construction -- NOT touched by any candidate change"),
    .e_row("E3c", "share in LASSO candidate construction",
           pick("cv.glmnet|glmnet"), "% of total",
           "candidate construction -- NOT touched by any candidate change"),
    .e_row("E4", "share in coxph() specifically", pick("coxph"),
           "% of total", "upper bound on what the fast engine can remove"),
    .e_row("E5", "share in survfit/summary.survfit",
           pick("survfit"), "% of total",
           "median columns; unused for selection when m1.threshold = Inf"),
    .e_row("E6", "share in model.frame/model.matrix/terms",
           pick("model.frame|model.matrix|terms"), "% of total",
           "pure formula dispatch -- the overhead coxph.fit() removes")
  )
  unlink(pf)
  do.call(rbind, rows)
}


#' Project the saving implied by the measured shares
#'
#' Deliberately conservative: applies each candidate's measured speedup only
#' to the share of runtime the profile attributes to it (Amdahl), so the
#' projection cannot exceed what the profile supports.
endtoend_projection <- function(prof, micro) {
  share <- function(id) {
    v <- prof$value[prof$id == id]
    if (!length(v) || is.na(v[1])) NA_real_ else v[1] / 100
  }
  spd <- function(id) {
    v <- micro$speedup[micro$id == id]
    if (!length(v)) NA_real_ else stats::median(v, na.rm = TRUE)
  }

  s_cox <- share("E4"); x_cox <- spd("M1")
  amdahl <- function(p, s) if (is.na(p) || is.na(s)) NA_real_ else
    1 / ((1 - p) + p / s)

  rbind(
    .e_row("E7", "projected speedup, fast Cox engine only",
           amdahl(s_cox, x_cox), "x",
           sprintf("Amdahl on measured coxph share %.1f%% at %.1fx",
                   100 * s_cox, x_cox)),
    .e_row("E8", "NOTE", NA, "",
           paste("Futility stopping reduces the NUMBER of fits and so",
                 "multiplies on top of E7; its factor is data-dependent and",
                 "is measured directly in stage P2, not projected here."))
  )
}


run_endtoend <- function(cfg, df, micro = NULL) {
  message("  E: profiling a sequential forestsearch() run")
  prof <- endtoend_profile(df, cfg$fs_args)
  if (!is.null(micro)) prof <- rbind(prof, endtoend_projection(prof, micro))
  prof
}
