# =============================================================================
# assembly_counts_2026-09-02.R -- section 5 substitution proof for the
# assembly-skip task: count the assembly's constructor
# (data.table::data.table) and the medians-site constructor (data.frame)
# during one F2 and one F5 call at the INSTALLED version, attributing each
# dispatch to the innermost forestsearch function on the stack.
# At 0.3.5 expected: F2 (unadjusted) 0 data.table constructions from
# fit_cox_for_subgroup; 121 data.frame constructions from
# evaluate_combination_with_status (the med_df medians site, a separate
# constructor); F5 (adjusted) still constructs one frame per fitted
# candidate (1,410).
# =============================================================================
suppressPackageStartupMessages(library(forestsearch))
suppressPackageStartupMessages(library(survival))
D <- "dev/glm-continuous-sims"
ns <- asNamespace("forestsearch")
cat(sprintf("installed forestsearch %s | %s\n",
            as.character(utils::packageVersion("forestsearch")), R.version.string))
cat(sprintf("git HEAD: %s\n", system("git rev-parse --short HEAD", intern = TRUE)))

# ---- F2/F5 gbsg fixture (verbatim from medians_baseline_2026-09-02.R) -------
df_surv <- survival::gbsg
df_surv <- within(df_surv, { id <- seq_len(nrow(df_surv)); time_months <- rfstime / 30.4375; grade3 <- ifelse(grade == "3", 1, 0) })
args_surv <- function(extra = list()) {
  base <- list(
    df.analysis = df_surv, outcome.name = "time_months", event.name = "status",
    treat.name = "hormon", id.name = "id",
    confounders.name = c("age", "meno", "size", "grade3", "nodes", "pgr", "er"),
    is.RCT = TRUE, seedit = 8316951, est.scale = "hr",
    use_lasso = FALSE, use_grf = FALSE, use_dina = FALSE,
    conf.cont_jcuts = list(er = 10, pgr = 10), collapse_cuts = TRUE,
    conf_force = c("er <= 0", "pgr <= 0"),
    subgroup_method = "consistency", sg_focus = "effMaxSG",
    selection_rule = "neighborhood", effect_neighborhood = 0.10,
    mr_inference = FALSE,
    n.min = 60, d0.min = 10, d1.min = 10, maxk = 2,
    hr.threshold = 1.0, hr.consistency = 1.0, pconsistency.threshold = 0.90,
    minp = 0.025, use_twostage = TRUE, stop_threshold = NULL,
    consistency_method = "resample", fs.splits = 1000, max.minutes = 3,
    vi.grf.min = NULL,
    show_candidate_summary = FALSE, details = FALSE, quiet = TRUE, plot.sg = FALSE,
    parallel_args = list(plan = "sequential", workers = 1L))
  utils::modifyList(base, extra)
}

# ---- constructor tracer -----------------------------------------------------
.ce <- new.env()
.ce$ns_fns <- ls(ns, all.names = TRUE)
attribute_stack <- quote({
  fnames <- sub("^.*::", "", vapply(sys.calls(), function(cl) {
    f <- cl[[1]]
    if (is.name(f)) as.character(f) else paste(deparse(f, nlines = 1L), collapse = "")
  }, character(1)))
  hit <- rev(fnames[fnames %in% .ce$ns_fns])
  if (length(hit)) hit[1] else "<non-forestsearch>"
})
count_case <- function(name, args) {
  cat(sprintf("---- %s\n", name))
  .ce$dt_by <- character(0); .ce$df_by <- character(0)
  suppressMessages(trace("data.table", where = asNamespace("data.table"), print = FALSE,
    tracer = bquote(.ce$dt_by <- c(.ce$dt_by, eval(.(attribute_stack))))))
  suppressMessages(trace("data.frame", where = baseenv(), print = FALSE,
    tracer = bquote(.ce$df_by <- c(.ce$df_by, eval(.(attribute_stack))))))
  set.seed(102L)
  fs <- suppressWarnings(do.call(forestsearch, args))
  suppressMessages(untrace("data.table", where = asNamespace("data.table")))
  suppressMessages(untrace("data.frame", where = baseenv()))
  cat(sprintf("  selected: %s\n", paste(fs$sg.harm, collapse = " & ")))
  cat("  data.table() dispatches by innermost forestsearch caller:\n")
  print(table(.ce$dt_by))
  cat("  data.frame() dispatches by innermost forestsearch caller:\n")
  print(table(.ce$df_by))
  list(name = name, dt_by = table(.ce$dt_by), df_by = table(.ce$df_by))
}
res_F2 <- count_case("F2_surv_gbsg (unadjusted)", args_surv())
res_F5 <- count_case("F5_surv_adjusted", args_surv(list(adjust_covariates = "age")))
saveRDS(list(F2 = res_F2, F5 = res_F5,
             pkg_version = as.character(utils::packageVersion("forestsearch")),
             built_at = Sys.time()),
        file.path(D, "assembly_counts_2026-09-02.rds"))
cat("\nwritten:", file.path(D, "assembly_counts_2026-09-02.rds"), "\n")
