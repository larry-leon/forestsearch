# stage0_common.R -- shared setup for the bootstrap-efficiency Stage-0
# measurement (task cc_task_bootstrap_stage0_2026-10-26.md).
#
# MEASUREMENT ONLY.  Nothing here edits any file of the forestsearch package,
# on disk or in the installed library.  s2 uses base::trace() on two functions
# IN THE RUNNING SESSION (forestsearch:::fit_cox_for_subgroup and
# survival::coxph.fit) to record their arguments, and untrace()s them before
# any timing; the installed package files are untouched.
#
# Audit basis: CRAN forestsearch 0.4.0 (Date/Publication 2026-10-10).

suppressPackageStartupMessages({
  library(forestsearch)
  library(survival)
  library(future)
  library(doFuture)
  library(data.table)
})

OUT <- Sys.getenv("STAGE0_OUT", "results")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

out_path <- function(...) file.path(OUT, ...)

save_result <- function(x, name) {
  saveRDS(x, out_path(paste0(name, ".rds")))
  invisible(x)
}

# key = value lines, appended to results/keyvalues.txt (the report quotes these)
kv <- function(...) {
  vals <- list(...)
  stopifnot(!is.null(names(vals)), all(nzchar(names(vals))))
  lines <- vapply(names(vals), function(k) {
    v <- vals[[k]]
    v <- if (is.numeric(v)) paste(format(v, digits = 6, trim = TRUE), collapse = ", ")
         else paste(as.character(v), collapse = ", ")
    sprintf("%s = %s", k, v)
  }, character(1))
  cat(lines, sep = "\n")
  cat(lines, file = out_path("keyvalues.txt"), sep = "\n", append = TRUE)
  invisible(lines)
}

# --- Fixture V: the gbsg configuration of the package's own vignette ---------
# inst/doc/survival-post-selection.qmd (CRAN 0.4.0), chunks `data` and `fit`,
# with two changes that give the per-replicate shape of the bootstrap:
#   * mr_inference = FALSE (bootstrap_results() strips mr_inference from every
#     replicate under the default mr_in_replicates = FALSE), and
#     mr_inference_args dropped with it;
#   * parallel_args = list(plan = "sequential") -- as the vignette already has,
#     and as every bootstrap replicate forces (CATEGORY 3 in bootstrap_results()).
# Everything else is the vignette's call verbatim.
stage0_fixture <- function() {
  df.analysis <- survival::gbsg
  df.analysis <- within(df.analysis, {
    id          <- seq_len(nrow(df.analysis))
    time_months <- rfstime / 30.4375
    grade3      <- ifelse(grade == "3", 1, 0)
  })
  confounders.name <- c("age", "meno", "size", "grade3", "nodes", "pgr", "er")
  args <- list(
    df.analysis = df.analysis,
    outcome.name = "time_months", event.name = "status",
    treat.name = "hormon", id.name = "id",
    confounders.name = confounders.name,
    is.RCT = TRUE, seedit = 8316951, est.scale = "hr",
    use_lasso = FALSE, use_grf = FALSE, use_dina = FALSE,
    conf.cont_jcuts = list(er = 10, pgr = 10), collapse_cuts = TRUE,
    conf_force = c("er <= 0", "pgr <= 0"),
    subgroup_method = "consistency", sg_focus = "effMaxSG",
    selection_rule = "neighborhood", effect_neighborhood = 0.20,
    mr_inference = FALSE,
    n.min = 60, d0.min = 10, d1.min = 10, maxk = 2,
    hr.threshold = 1.0, hr.consistency = 1.0, pconsistency.threshold = 0.80,
    fs.splits = 200, max.minutes = 3,
    details = FALSE, quiet = TRUE, plot.sg = FALSE,
    parallel_args = list(plan = "sequential")
  )
  list(label = "V: survival-post-selection vignette gbsg config, mr_inference = FALSE",
       df = df.analysis, args = args)
}

run_fixture_call <- function(fx) {
  suppressWarnings(suppressMessages(do.call(forestsearch::forestsearch, fx$args)))
}

# Candidate-flow counters of one forestsearch() result (survival consistency path)
candidate_counts <- function(fs) {
  fc <- fs$find.grps$filter_counts
  gc <- fs$grp.consistency
  list(
    n_evaluated_search   = fc$n_evaluated,
    n_reached_cox_fit    = fc$n_passed_sample_size,   # status >= 5: fit attempted
    n_cox_fit_ok         = fc$n_passed_cox,           # status >= 6
    n_passed_effect      = fc$n_passed_hr,            # status 7 (medians computed here)
    n_consistency_total  = gc$n_candidates_total,
    n_consistency_evaluated = gc$n_candidates_evaluated,
    n_consistency_passed = gc$n_passed,
    early_stop_triggered = isTRUE(gc$early_stop_triggered),
    sg_harm = paste(fs$sg.harm, collapse = " & ")
  )
}
