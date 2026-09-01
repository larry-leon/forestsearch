# =============================================================================
# medians_profile_2026-09-02.R -- section 7 profiler for the medians-on-
# survivors task (dev/tasks/cc_task_medians_on_survivors_2026-09-02.md)
# -----------------------------------------------------------------------------
# One F2 gbsg single call under Rprof(interval = 0.010, line.profiling = TRUE),
# the established procedure and buckets (bootstrap_reprofile_2026-09-01.R,
# BUCKETS verbatim).  Env-switched like the battery:
#   Stage A (0.3.3):  Rscript ... .R          -> medians_profile_F2_before_2026-09-02.Rprof
#   Stage C (0.3.4):  FS_MEDIANS_OUT=postchange Rscript ... .R
#                                             -> medians_profile_F2_after_2026-09-02.Rprof
# =============================================================================
suppressPackageStartupMessages(library(forestsearch))
suppressPackageStartupMessages(library(survival))
D <- "dev/glm-continuous-sims"
OUT <- Sys.getenv("FS_MEDIANS_OUT", "baseline")
TAG <- if (OUT == "baseline") "before" else "after"
cat(sprintf("stage: %s (%s) | installed forestsearch %s | %s\n",
            OUT, TAG, as.character(utils::packageVersion("forestsearch")), R.version.string))
cat(sprintf("git HEAD: %s\n", system("git rev-parse --short HEAD", intern = TRUE)))

# ---- F2 gbsg fixture (verbatim from medians_baseline_2026-09-02.R) ----------
df_surv <- survival::gbsg
df_surv <- within(df_surv, { id <- seq_len(nrow(df_surv)); time_months <- rfstime / 30.4375; grade3 <- ifelse(grade == "3", 1, 0) })
args_surv <- list(
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

# ---- Rprof bucket parser: BUCKETS is the 08-30 list VERBATIM ----------------
BUCKETS <- list(
  medians      = c("survfit", "survfit.formula", "survfitKM", "summary.survfit", "survfit.Surv"),
  grf          = c("causal_forest", "causal_survival_forest", "fit_causal_forest_glm", "variable_importance",
                   "regression_forest", "fit_grf_variable_importance"),
  fit          = c("fit_glm_for_subgroup", "fit_cox_for_subgroup"),
  cuts         = c("get_FSdata"),
  consistency  = c("evaluate_subgroup_consistency", "evaluate_consistency_twostage", "consistency_resample",
                   ".consistency_via_splits"),
  future_overhead = c("future_lapply", "future_xapply", "getGlobalsAndPackages", "findGlobals", "globalsOf",
                      "makeFuture", "run.Future", "resolve", "value", "result", "Future", "future", "sequential",
                      "SequentialFuture", "getGlobalsAndPackagesXApply", "packagesOf", "cleanup"),
  dedup_select = c("remove_near_duplicate_subgroups", ".maxeff_membership_dedup", "sort_subgroups", "sort_subgroups_preview"),
  enumeration  = c("subgroup.search", "generate_combination_indices", "evaluate_combination_with_status",
                   "meets_prevalence_threshold", "extract_idx_flagredundancy", "get_covs_in", "has_positive_variance",
                   "format_search_results", "get_subgroup_membership", "calculate_max_combinations"),
  consistency_other = c("subgroup.consistency"))
parse_rprof <- function(file) {
  ln <- readLines(file)
  hdr <- grepl("sample.interval=", ln, fixed = TRUE) | grepl("^#File ", ln)
  interval <- as.numeric(sub(".*sample.interval=", "", ln[grepl("sample.interval=", ln, fixed = TRUE)][1])) / 1e6
  samples <- ln[!hdr]
  stacks <- lapply(strsplit(samples, " ", fixed = TRUE), function(tok) {
    tok <- tok[nzchar(tok)]; tok <- tok[!grepl("^[0-9]+#[0-9]+$", tok)]; gsub('"', "", tok) })
  bucket_of <- function(st) { for (b in names(BUCKETS)) if (any(st %in% BUCKETS[[b]])) return(b); "everything_else" }
  bk <- vapply(stacks, bucket_of, character(1))
  tab <- table(factor(bk, levels = c(names(BUCKETS), "everything_else")))
  data.frame(bucket = names(tab), samples = as.integer(tab), secs = as.integer(tab) * interval,
             share = as.numeric(tab) / length(bk), stringsAsFactors = FALSE) -> out
  list(table = out, n_samples = length(bk), interval = interval, total_secs = length(bk) * interval,
       stacks = stacks, bucket_by_sample = bk)
}
top_n <- function(sr, which = c("by.self", "by.total"), n = 20) {
  which <- match.arg(which); x <- sr[[which]]; x <- x[order(-x[[if (which == "by.self") "self.time" else "total.time"]]), ]
  head(x, n)
}

# ---- warm-up, then the profiled call ----------------------------------------
set.seed(102L)
gc(); t0 <- proc.time()[[3]]
fs_w <- do.call(forestsearch, args_surv)
wall_warm <- proc.time()[[3]] - t0
cat(sprintf("warm-up call: %.2f s wall; selected: %s\n", wall_warm, paste(fs_w$sg.harm, collapse = " & ")))

pf <- file.path(D, sprintf("medians_profile_F2_%s_2026-09-02.Rprof", TAG))
set.seed(102L)
gc(); Rprof(pf, interval = 0.010, line.profiling = TRUE, memory.profiling = FALSE)
t0 <- proc.time()[[3]]
fs_p <- do.call(forestsearch, args_surv)
wall_p <- proc.time()[[3]] - t0
Rprof(NULL)
stopifnot(identical(fs_p$sg.harm, fs_w$sg.harm))
sr <- summaryRprof(pf, lines = "hide")
pb <- parse_rprof(pf)
cat(sprintf("profiled call: %.2f s wall; %d samples x %.3f s = %.2f s sampled (%.1f%% of wall)\n",
            wall_p, pb$n_samples, pb$interval, pb$total_secs, 100 * pb$total_secs / wall_p))
cat("bucket table (share of sampled time):\n"); print(pb$table, digits = 4, row.names = FALSE)
cat("top 20 by SELF time:\n");  print(top_n(sr, "by.self"))
cat("top 20 by TOTAL time:\n"); print(top_n(sr, "by.total"))
saveRDS(list(stage = OUT, tag = TAG, wall_warm = wall_warm, wall_profiled = wall_p,
             buckets = pb$table, n_samples = pb$n_samples, interval = pb$interval,
             pkg_version = as.character(utils::packageVersion("forestsearch")),
             built_at = Sys.time()),
        file.path(D, sprintf("medians_profile_F2_%s_2026-09-02.rds", TAG)))
cat("\nwritten:", pf, "\n")
