# stage0b_common.R -- shared setup for the bootstrap-efficiency Stage-0b
# measurement (task cc_task_bootstrap_stage0b_2026-10-26.md): the GLM paths
# and the consistency stage, all three outcome types.
#
# MEASUREMENT ONLY.  Nothing here edits any file of the forestsearch package,
# on disk or in the installed library.  b2 uses base::trace() on package
# functions IN THE RUNNING SESSION to record their arguments, and untrace()s
# them before any timing.  Every "formula path" timed below is the installed
# package's own function, called unmodified on the recorded arguments.
#
# Audit basis: CRAN forestsearch 0.4.0 (Date/Publication 2026-10-10).

suppressPackageStartupMessages({
  library(forestsearch)
  library(survival)
  library(future)
  library(doFuture)
  library(data.table)
})

OUT <- Sys.getenv("STAGE0B_OUT", "results")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
out_path <- function(...) file.path(OUT, ...)
save_result <- function(x, name) { saveRDS(x, out_path(paste0(name, ".rds"))); invisible(x) }

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

# ---- fixtures ----------------------------------------------------------------
# Survival: fixture V of Stage 0 (the package's own vignette gbsg configuration,
# mr_inference = FALSE, sequential).  Continuous and binary: the committed
# ACTG175 application configurations, transplanted by CC into
# fixture_continuous.R / fixture_binary.R (see the task document, section 2).
# Each fixture function returns list(label, df, args, source).
fixture_survival <- function() {
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
  list(label = "survival: vignette gbsg config (Stage-0 fixture V)", df = df.analysis, args = args,
       source = "inst/doc/survival-post-selection.qmd (CRAN 0.4.0), mr_inference = FALSE")
}

source("fixture_continuous.R")   # defines fixture_continuous()
source("fixture_binary.R")       # defines fixture_binary()

FIXTURES <- list(survival = fixture_survival, continuous = fixture_continuous, binary = fixture_binary)

run_fixture_call <- function(fx) {
  suppressWarnings(suppressMessages(do.call(forestsearch::forestsearch, fx$args)))
}

candidate_counts <- function(fs) {
  fc <- fs$find.grps$filter_counts
  gc <- fs$grp.consistency
  list(
    n_evaluated_search      = fc$n_evaluated,
    n_reached_fit           = fc$n_passed_sample_size,
    n_fit_ok                = fc$n_passed_cox,
    n_passed_effect         = fc$n_passed_hr,
    n_consistency_total     = gc$n_candidates_total,
    n_consistency_evaluated = gc$n_candidates_evaluated,
    n_consistency_passed    = gc$n_passed,
    early_stop_triggered    = isTRUE(gc$early_stop_triggered),
    sg_harm                 = paste(fs$sg.harm, collapse = " & ")
  )
}

# capture environment used by the tracers (lives in .GlobalEnv so a tracer
# evaluated inside a package function's frame can reach it)
new_cap <- function() {
  e <- new.env()
  e$n <- 0L
  e$items <- vector("list", 65536L)
  e$push <- function(item) {
    e$n <- e$n + 1L
    if (e$n > length(e$items)) length(e$items) <- 2L * length(e$items)
    e$items[[e$n]] <- item
  }
  e
}
cap_items <- function(e) e$items[seq_len(e$n)]

# median of timing reps, alternating order to avoid order bias
time_pair <- function(f_a, f_b, reps = 3L) {
  ta <- tb <- numeric(reps)
  for (r in seq_len(reps)) {
    if (r %% 2L == 1L) { ta[r] <- system.time(f_a())[["elapsed"]]; tb[r] <- system.time(f_b())[["elapsed"]] }
    else               { tb[r] <- system.time(f_b())[["elapsed"]]; ta[r] <- system.time(f_a())[["elapsed"]] }
  }
  list(a = ta, b = tb, a_med = median(ta), b_med = median(tb))
}
