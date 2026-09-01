# =============================================================================
# assembly_settleB_2026-09-02.R -- section 7 read-only settlement: is
# first-passer stop_threshold sound for maxeffCons?
# (dev/tasks/cc_task_assembly_skip_2026-09-02.md)
# -----------------------------------------------------------------------------
# The F3 anchor fixture (block verbatim from medians_baseline_2026-09-02.R /
# assembly_battery_2026-09-02.R) run twice at the installed version:
#   (a) as-is: stop_threshold defaults to pconsistency.threshold
#       (forestsearch() formal, R/forestsearch_main.R L1235) -> early stop at
#       the first passer;
#   (b) stop_threshold = NULL -> all candidates evaluated.
# Compared: the selected subgroup, and the selected row's values in the
# consistency output; the extra rows the full run carries are reported.
# NO code change either way -- measurement only.
# =============================================================================
suppressPackageStartupMessages(library(forestsearch))
stopifnot(requireNamespace("speff2trial", quietly = TRUE))
D <- "dev/glm-continuous-sims"
cat(sprintf("installed forestsearch %s | %s\n",
            as.character(utils::packageVersion("forestsearch")), R.version.string))
cat(sprintf("git HEAD: %s\n", system("git rev-parse --short HEAD", intern = TRUE)))

# ---- F3 anchor fixture (verbatim from medians_baseline_2026-09-02.R) --------
cont_vars <- c("age", "preanti", "wtkg", "karnof", "cd40", "cd80")
bin_vars  <- c("hemo", "homo", "drugs", "race", "gender", "symptom")
actg_oc <- subset(speff2trial::ACTG175, arms %in% c(1L, 3L))
actg_oc$id <- seq_len(nrow(actg_oc))
actg_oc$treat <- ifelse(actg_oc$arms == 1L, 1L, 0L)
actg_oc <- actg_oc[!is.na(actg_oc$cd420), ]
actg_oc$y_decline <- actg_oc$cd40 - actg_oc$cd420
for (v in bin_vars) actg_oc[[v]] <- as.factor(actg_oc[[v]])
stopifnot(nrow(actg_oc) == 1083L)
args_anchor <- list(
  df.analysis      = actg_oc,
  confounders.name = c(cont_vars, bin_vars),
  outcome.name     = "y_decline",
  treat.name       = "treat",
  id.name          = "id",
  outcome_type     = "continuous",
  effect_measure   = "MD",
  adverse_outcome  = TRUE,
  seedit           = 8316951,
  sg_focus         = "maxeffCons",
  selection_rule   = "neighborhood",
  consistency_method = "resample",
  effect.threshold       = 10,
  consistency.threshold  = 5,
  pconsistency.threshold = 0.90,
  use_twostage     = TRUE,
  conf.cont_jcuts  = list(age = 10, preanti = 10, wtkg = 10,
                          karnof = 10, cd40 = 10, cd80 = 10),
  cut_type         = "default",
  maxk             = 2L,
  n.min            = 60L,
  d0.min           = 10L,
  d1.min           = 10L,
  fs.splits        = 500L,
  use_lasso        = FALSE,
  use_grf          = FALSE,
  use_dina         = FALSE,
  is.RCT           = TRUE,
  parallel_args    = list(plan = "sequential"),
  details          = FALSE,
  quiet            = TRUE)

run_anchor <- function(label, extra) {
  cat(sprintf("---- anchor: %s\n", label))
  # NOTE: modifyList() DELETES entries whose value is NULL, which would turn an
  # explicit stop_threshold = NULL back into the formal default -- append with
  # c() instead so the NULL is actually passed to forestsearch().
  args <- c(args_anchor, extra)
  set.seed(103L)
  t0 <- proc.time()[[3]]
  fs <- suppressWarnings(do.call(forestsearch, args))
  wall <- proc.time()[[3]] - t0
  gc_ <- fs$grp.consistency
  outsg <- tryCatch(as.data.frame(gc_$out_sg$result), error = function(e) NULL)
  cat(sprintf("  wall %.2f s; selected: %s; n_evaluated %s of %s; early_stop %s (cand %s); consistency rows %s\n",
              wall, if (is.null(fs$sg.harm)) "<none>" else paste(fs$sg.harm, collapse = " & "),
              gc_$n_candidates_evaluated %||% NA, gc_$n_candidates_total %||% NA,
              gc_$early_stop_triggered %||% NA, gc_$early_stop_candidate %||% NA,
              if (!is.null(outsg)) nrow(outsg) else NA))
  list(label = label, wall = wall, sg.harm = fs$sg.harm, fs = fs, out_sg = outsg)
}
`%||%` <- function(a, b) if (is.null(a)) b else a

a <- run_anchor("default_stop_threshold (= pconsistency.threshold 0.90)", list())
b <- run_anchor("stop_threshold_NULL (all evaluated)", list(stop_threshold = NULL))

cat("\n==== settlement ====\n")
cat("selection (default stop): ", paste(a$sg.harm, collapse = " & "), "\n")
cat("selection (full eval):    ", paste(b$sg.harm, collapse = " & "), "\n")
cat("selected subgroup identical:", identical(a$sg.harm, b$sg.harm), "\n")

# consistency-output comparison: the selected row's values, and what extra the
# full run carries.  Locate the consistency result table on each fit.
find_cons_table <- function(fs) {
  x2 <- tryCatch(as.data.frame(fs$grp.consistency$out_sg$result), error = function(e) NULL)
  if (!is.null(x2) && nrow(x2) > 0) return(x2)
  NULL
}
ta <- find_cons_table(a$fs); tb <- find_cons_table(b$fs)
cat("\nconsistency table (default stop):", if (is.null(ta)) "<not found>" else sprintf("%d rows", nrow(ta)), "\n")
if (!is.null(ta)) print(ta, digits = 10)
cat("\nconsistency table (full eval):", if (is.null(tb)) "<not found>" else sprintf("%d rows", nrow(tb)), "\n")
if (!is.null(tb)) print(utils::head(tb, 12), digits = 10)
if (!is.null(ta) && !is.null(tb) && nrow(ta) >= 1) {
  common <- intersect(names(ta), names(tb))
  sel_a <- ta[1, common, drop = FALSE]
  sel_b <- tb[1, common, drop = FALSE]
  rownames(sel_a) <- rownames(sel_b) <- NULL
  cat("\nselected row identical across runs (common columns):",
      isTRUE(all.equal(sel_a, sel_b, tolerance = 0)) && identical(unname(as.list(sel_a)), unname(as.list(sel_b))), "\n")
  cat("identical():", identical(sel_a, sel_b), "\n")
}
saveRDS(list(default_stop = list(sg.harm = a$sg.harm, wall = a$wall, table = ta,
                                 n_evaluated = a$fs$n_candidates_evaluated),
             full_eval    = list(sg.harm = b$sg.harm, wall = b$wall, table = tb,
                                 n_evaluated = b$fs$n_candidates_evaluated),
             selection_identical = identical(a$sg.harm, b$sg.harm),
             pkg_version = as.character(utils::packageVersion("forestsearch")),
             built_at = Sys.time()),
        file.path(D, "assembly_settleB_2026-09-02.rds"))
cat("\nwritten:", file.path(D, "assembly_settleB_2026-09-02.rds"), "\n")
