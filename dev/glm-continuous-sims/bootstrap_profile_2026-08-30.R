# =============================================================================
# bootstrap_profile_2026-08-30.R -- where does one forestsearch() call spend
# its time?  (task: dev/tasks/cc_task_bootstrap_profile_2026-08-30.md)
# -----------------------------------------------------------------------------
# Measurement only.  No R/ change, no optimisation.  Two fixtures:
#   cont -- the MD40 DGM the OC-wrapper scripts rebuild (ACTG175 CD4 change,
#           calibrated harm), ONE trial at n = 500 (sim_id 1: seed 8316951 + 1),
#           the drivers' forestsearch() arguments (13 confounders incl. str2).
#   surv -- survival::gbsg (686 rows) with the effMaxSG application's arguments
#           (quarto/applications/gbsg/analysis_gbsg_survival_effMaxSG.qmd).
# Both are run in the configuration a bootstrap REPLICATE runs them in
# (bootstrap_analysis_dofuture.R L573-585): parallel_args plan = "sequential",
# mr_inference = FALSE, no plots, quiet.  The MR share of the drivers' own
# single call is measured separately (cont, vi.grf.min = NULL, mr on) as a
# side line, because the replicates do not carry it.
#
# Per configuration, two calls:
#   (i)  an INSTRUMENTED call -- trace() on the namespace functions, counting
#        candidates reaching a fit (with per-fit elapsed time), candidates
#        reaching the consistency screen, closed-form NA rates, literal-split
#        fallbacks, and the dedup's rows before/after.  trace() is an in-memory
#        modification of the loaded namespace; nothing on disk changes.
#   (ii) a PROFILED call -- Rprof(interval = 0.010, line.profiling = TRUE,
#        memory.profiling = FALSE) with every trace removed; buckets are
#        assigned per sample from the whole stack, first match in priority
#        order (medians > GRF > fit > cuts > consistency > dedup/selection >
#        enumeration > else), so a survfit() under fit_cox_for_subgroup() is
#        counted as "medians", not "fit".
# Then the bootstrap check: forestsearch_bootstrap_dofuture() at nb_boots = 20,
# plan = "sequential", on the continuous fixture at the installed default.
#
# Usage:  Rscript dev/glm-continuous-sims/bootstrap_profile_2026-08-30.R
# Writes (dev/glm-continuous-sims/): bootstrap_profile_2026-08-30.rds,
#   bootstrap_profile_<config>_2026-08-30.Rprof (raw sample files), and the log.
# =============================================================================
suppressPackageStartupMessages(devtools::load_all(".", quiet = TRUE))
stopifnot(requireNamespace("speff2trial", quietly = TRUE))
D <- "dev/glm-continuous-sims"
ns <- asNamespace("forestsearch")
vi_default <- formals(forestsearch)$vi.grf.min
cat(sprintf("installed forestsearch %s; vi.grf.min default = %s\n",
            as.character(utils::packageVersion("forestsearch")), deparse(vi_default)))
NB_BOOTS <- 20L

# ---- fixtures -----------------------------------------------------------------
actg_frame <- function() {
  actg_df <- subset(speff2trial::ACTG175, arms %in% c(1L, 3L))
  actg_df$id <- seq_len(nrow(actg_df))
  actg_df$treat <- 1L - ifelse(actg_df$arms == 1L, 1L, 0L)
  actg_df$cd4_change <- actg_df$cd420 - actg_df$cd40
  actg_df <- actg_df[!is.na(actg_df$cd420), ]
  actg_df$z1  <- as.factor(ifelse(actg_df$age > 34, 1L, 0L))
  actg_df$z2  <- as.factor(ifelse(actg_df$preanti <= 744.5, 1L, 0L))
  actg_df$z3  <- as.factor(ifelse(actg_df$wtkg <= 75, 1L, 0L))
  actg_df$z4  <- as.factor(ifelse(actg_df$karnof <= median(actg_df$karnof), 1L, 0L))
  actg_df$z5  <- as.factor(ifelse(actg_df$cd40 <= median(actg_df$cd40), 1L, 0L))
  actg_df$z6  <- as.factor(ifelse(actg_df$cd80 <= median(actg_df$cd80), 1L, 0L))
  actg_df$z7  <- as.factor(actg_df$hemo);  actg_df$z8  <- as.factor(actg_df$homo)
  actg_df$z9  <- as.factor(actg_df$drugs); actg_df$z10 <- as.factor(actg_df$race)
  actg_df$z11 <- as.factor(actg_df$gender);actg_df$z12 <- as.factor(actg_df$symptom)
  # the driver converts the analysis binaries to factors before the DGM (L301)
  for (v in c("hemo", "homo", "drugs", "race", "gender", "symptom", "str2")) actg_df[[v]] <- as.factor(actg_df[[v]])
  actg_df
}
PAY <- "quarto/simulations/actg175/continuous/mr_md_harm/fs_maxeffCons_mr_md40_knoise0_n500_s1000_d5000/fs_maxeffCons_mr_md40_knoise0_n500_res_1_1000.rds"
truth <- readRDS(PAY)$truth
dgm <- generate_glm_dgm(
  data = actg_frame(), factor_vars = paste0("z", 1:12), outcome_var = "cd4_change",
  treatment_var = "treat", outcome_type = "continuous", effect_measure = "MD",
  subgroup_vars = c("z1", "z2"), subgroup_cuts = list(z1 = 1L, z2 = 1L),
  model = "alt", k_treat = 1, k_inter = truth$beta_inter, n_super = 5000L,
  seed = 8316951L, verbose = FALSE)
stopifnot(abs(dgm$hazard_ratios$harm_subgroup - truth$effect_Q) < 1e-9,
          abs(dgm$subgroup_info$proportion - truth$prevalence_Q) < 1e-9)
SEED_BASE <- 8316951L; SIM_ID <- 1L
df_cont <- simulate_from_glm_dgm(dgm, n = 500L, seed = SEED_BASE + SIM_ID)
cat(sprintf("cont fixture: n = %d, columns %s\n", nrow(df_cont), paste(names(df_cont), collapse = " ")))
confs_cont <- c("age", "preanti", "wtkg", "karnof", "cd40", "cd80",
                "hemo", "homo", "drugs", "race", "gender", "symptom", "str2")
stopifnot(all(confs_cont %in% names(df_cont)), all(c("y_sim", "treat_sim", "id") %in% names(df_cont)))
# the driver's call (sim_fs_maxeffCons_mr_md40_knoise0_n500_batch_1_1000.qmd L561-590), replicate configuration
args_cont <- function(vi, mr = FALSE, pa = list(plan = "sequential", workers = 1L)) list(
  df.analysis = df_cont, confounders.name = confs_cont,
  outcome.name = "y_sim", treat.name = "treat_sim", id.name = "id",
  outcome_type = "continuous", effect_measure = "MD",
  effect.threshold = 30, consistency.threshold = 10, pconsistency.threshold = 0.90,
  fs.splits = 400L, n.min = 60L, d0.min = 12L, d1.min = 12L, maxk = 2L,
  vi.grf.min = vi, sg_focus = "maxeffCons", selection_rule = "neighborhood",
  effect_neighborhood = 0.10, stop_threshold = NULL, consistency_method = "resample",
  conf.cont_jcuts = list(age = 10, preanti = 10),
  use_lasso = FALSE, use_dina = FALSE, use_grf = FALSE, use_twostage = TRUE,
  is.RCT = TRUE, adverse_outcome = FALSE, details = FALSE, quiet = TRUE,
  seedit = SEED_BASE + SIM_ID, parallel_args = pa,
  mr_inference = mr,
  mr_inference_args = if (mr) list(ci_method = "ij", draws = 5000L, include_complement = TRUE) else NULL)

df_surv <- survival::gbsg
df_surv <- within(df_surv, { id <- seq_len(nrow(df_surv)); time_months <- rfstime / 30.4375; grade3 <- ifelse(grade == "3", 1, 0) })
cat(sprintf("surv fixture: gbsg n = %d, events %d\n", nrow(df_surv), sum(df_surv$status)))
# analysis_gbsg_survival_effMaxSG.qmd L539-585 with the effMaxSG focus profile resolved
# (hrMaxSG: thresholds live, selection_rule neighborhood, stop_threshold NULL, minp 0.025, twostage TRUE)
args_surv <- function(vi) list(
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
  vi.grf.min = vi,
  show_candidate_summary = FALSE, details = FALSE, quiet = TRUE, plot.sg = FALSE,
  parallel_args = list(plan = "sequential", workers = 1L))

# ---- instrumentation via trace() -----------------------------------------------
.prof_env <- new.env()
reset_counters <- function() {
  .prof_env$fit_calls <- 0L; .prof_env$fit_times <- numeric(0); .prof_env$fit_null <- 0L
  .prof_env$cr_calls <- 0L; .prof_env$cr_na <- 0L
  .prof_env$split_calls <- 0L
  .prof_env$cons_calls <- 0L; .prof_env$cons_pass <- 0L
  .prof_env$dedup_before <- NA_integer_; .prof_env$dedup_after <- NA_integer_
  .prof_env$med_calls <- 0L; .prof_env$med_times <- numeric(0)
  .prof_env$grf_calls <- 0L
}
TRACED <- character(0)
tr <- function(what, tracer = NULL, exit = NULL) {
  if (!exists(what, envir = ns, inherits = FALSE)) return(invisible(NULL))
  suppressMessages(trace(what, tracer = tracer, exit = exit, where = ns, print = FALSE))
  TRACED <<- c(TRACED, what); invisible(NULL)
}
set_traces <- function() {
  for (f in c("fit_glm_for_subgroup", "fit_cox_for_subgroup")) tr(f,
    tracer = quote(.prof_env$t0_fit <- proc.time()[[3]]),
    exit = quote({ .prof_env$fit_calls <- .prof_env$fit_calls + 1L
                   .prof_env$fit_times <- c(.prof_env$fit_times, proc.time()[[3]] - .prof_env$t0_fit)
                   if (is.null(returnValue())) .prof_env$fit_null <- .prof_env$fit_null + 1L }))
  tr("consistency_resample",
     exit = quote({ .prof_env$cr_calls <- .prof_env$cr_calls + 1L
                    if (is.na(returnValue()$rate_closed)) .prof_env$cr_na <- .prof_env$cr_na + 1L }))
  tr(".consistency_via_splits", tracer = quote(.prof_env$split_calls <- .prof_env$split_calls + 1L))
  for (f in c("evaluate_subgroup_consistency", "evaluate_consistency_twostage")) tr(f,
    tracer = quote(.prof_env$cons_calls <- .prof_env$cons_calls + 1L),
    exit = quote(if (!is.null(returnValue())) .prof_env$cons_pass <- .prof_env$cons_pass + 1L))
  tr("remove_near_duplicate_subgroups",
     tracer = quote(.prof_env$dedup_before <- nrow(hr_subgroups)),
     exit = quote(.prof_env$dedup_after <- nrow(returnValue())))
  tr(".maxeff_membership_dedup",
     tracer = quote(.prof_env$dedup_before <- nrow(found.hrs)),
     exit = quote(.prof_env$dedup_after <- nrow(returnValue())))
  tr("fit_causal_forest_glm", tracer = quote(.prof_env$grf_calls <- .prof_env$grf_calls + 1L))
}
clear_traces <- function() { for (f in unique(TRACED)) suppressMessages(untrace(f, where = ns)); TRACED <<- character(0) }

# ---- Rprof bucket parser --------------------------------------------------------
BUCKETS <- list(
  medians      = c("survfit", "survfit.formula", "survfitKM", "summary.survfit", "survfit.Surv"),
  grf          = c("causal_forest", "causal_survival_forest", "fit_causal_forest_glm", "variable_importance",
                   "regression_forest", "fit_grf_variable_importance"),
  fit          = c("fit_glm_for_subgroup", "fit_cox_for_subgroup"),
  cuts         = c("get_FSdata"),
  consistency  = c("evaluate_subgroup_consistency", "evaluate_consistency_twostage", "consistency_resample",
                   ".consistency_via_splits"),
  # the future machinery around the consistency screen (future_lapply batches: globals/package
  # resolution, future creation and resolution) -- outside the evaluator itself
  future_overhead = c("future_lapply", "future_xapply", "getGlobalsAndPackages", "findGlobals", "globalsOf",
                      "makeFuture", "run.Future", "resolve", "value", "result", "Future", "future", "sequential",
                      "SequentialFuture", "getGlobalsAndPackagesXApply", "packagesOf", "cleanup"),
  dedup_select = c("remove_near_duplicate_subgroups", ".maxeff_membership_dedup", "sort_subgroups", "sort_subgroups_preview"),
  enumeration  = c("subgroup.search", "generate_combination_indices", "evaluate_combination_with_status",
                   "meets_prevalence_threshold", "extract_idx_flagredundancy", "get_covs_in", "has_positive_variance",
                   "format_search_results", "get_subgroup_membership", "calculate_max_combinations"),
  consistency_other = c("subgroup.consistency"))    # remainder of subgroup.consistency() outside the above
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
  list(table = out, n_samples = length(bk), interval = interval, total_secs = length(bk) * interval, stacks = stacks)
}
top_n <- function(sr, which = c("by.self", "by.total"), n = 20) {
  which <- match.arg(which); x <- sr[[which]]; x <- x[order(-x[[if (which == "by.self") "self.time" else "total.time"]]), ]
  head(x, n)
}
cnt_snapshot <- function() list(
  fit_calls = .prof_env$fit_calls, fit_null = .prof_env$fit_null,
  fit_mean_secs = if (length(.prof_env$fit_times)) mean(.prof_env$fit_times) else NA_real_,
  fit_median_secs = if (length(.prof_env$fit_times)) median(.prof_env$fit_times) else NA_real_,
  fit_total_secs = sum(.prof_env$fit_times),
  cons_calls = .prof_env$cons_calls, cons_pass = .prof_env$cons_pass,
  cr_calls = .prof_env$cr_calls, cr_na = .prof_env$cr_na, split_calls = .prof_env$split_calls,
  dedup_before = .prof_env$dedup_before, dedup_after = .prof_env$dedup_after, grf_calls = .prof_env$grf_calls)

# ---- one configuration ------------------------------------------------------------
run_config <- function(name, args) {
  cat(sprintf("\n################ %s\n", name))
  # (i) instrumented
  reset_counters(); set_traces()
  gc(); t0 <- proc.time()[[3]]
  fs_i <- do.call(forestsearch, args)
  wall_i <- proc.time()[[3]] - t0
  clear_traces()
  cnt <- cnt_snapshot()
  fg <- fs_i$find.grps
  cat(sprintf("instrumented call: %.2f s wall; selected: %s\n", wall_i, paste(fs_i$sg.harm, collapse = " & ")))
  cat(sprintf("  candidates reaching a fit: %d (fit returned NULL: %d); mean %.4f s, median %.4f s per fit; fits total %.2f s (%.1f%% of wall)\n",
              cnt$fit_calls, cnt$fit_null, cnt$fit_mean_secs, cnt$fit_median_secs, cnt$fit_total_secs, 100 * cnt$fit_total_secs / wall_i))
  cat(sprintf("  hr.subgroups rows (passed the effect screen): %s; dedup rows before -> after: %s -> %s\n",
              if (!is.null(fg$out.found$hr.subgroups)) nrow(fg$out.found$hr.subgroups) else NA, cnt$dedup_before, cnt$dedup_after))
  cat(sprintf("  candidates reaching the consistency screen: %d (passed: %d); consistency_resample() calls %d, closed form NA %d, literal-split fallbacks %d\n",
              cnt$cons_calls, cnt$cons_pass, cnt$cr_calls, cnt$cr_na, cnt$split_calls))
  cat(sprintf("  GRF variable-importance forest fits: %d\n", cnt$grf_calls))
  if (!is.null(fg$n_candidates_evaluated)) cat(sprintf("  find.grps: n_candidates_total %s, n_candidates_evaluated %s, n_passed %s, early_stop %s\n",
              fg$n_candidates_total, fg$n_candidates_evaluated, fg$n_passed, fg$early_stop_triggered))
  if (!is.null(fg$out.found$prop_max_count)) cat(sprintf("  search coverage prop_max_count = %.3f (max.minutes cut-off if < 1)\n", fg$out.found$prop_max_count))
  # warm-up: untrace() restores the original closures, which the JIT then recompiles on
  # first use; a bootstrap worker pays that once per worker, not per replicate, so it is
  # excluded from the profiled call by one un-timed call here.
  invisible(do.call(forestsearch, args))
  # (ii) profiled -- 0.010 s is this platform's minimum sampling interval (Rprof warns below it)
  pf <- file.path(D, sprintf("bootstrap_profile_%s_2026-08-30.Rprof", name))
  gc(); Rprof(pf, interval = 0.010, line.profiling = TRUE, memory.profiling = FALSE)
  t0 <- proc.time()[[3]]
  fs_p <- do.call(forestsearch, args)
  wall_p <- proc.time()[[3]] - t0
  Rprof(NULL)
  stopifnot(identical(fs_p$sg.harm, fs_i$sg.harm))
  sr <- summaryRprof(pf, lines = "hide")
  pb <- parse_rprof(pf)
  cat(sprintf("profiled call: %.2f s wall; %d samples x %.3f s = %.2f s sampled\n", wall_p, pb$n_samples, pb$interval, pb$total_secs))
  cat("bucket table (share of sampled time):\n"); print(pb$table, digits = 4, row.names = FALSE)
  cat("top 20 by SELF time:\n");  print(top_n(sr, "by.self"))
  cat("top 20 by TOTAL time:\n"); print(top_n(sr, "by.total"))
  list(name = name, wall_instrumented = wall_i, wall_profiled = wall_p, counters = cnt, buckets = pb$table,
       n_samples = pb$n_samples, interval = pb$interval, by_self = top_n(sr, "by.self", 30), by_total = top_n(sr, "by.total", 40),
       sg.harm = fs_p$sg.harm, fs = fs_p, rprof_file = pf,
       vi.grf.min = args$vi.grf.min, mr_inference = isTRUE(args$mr_inference))
}

res <- list()
res$cont_viNULL  <- run_config("cont_viNULL",  args_cont(NULL))                 # bootstrap-replicate parallel_args: plan sequential, workers 1
res$cont_viNULL_seqloop <- run_config("cont_viNULL_seqloop", args_cont(NULL, pa = list(plan = "sequential")))  # the DRIVER's parallel_args: fails the plan+workers check -> plain loop
res$cont_vim02   <- run_config("cont_vim02",   args_cont(-0.2))
res$surv_viNULL  <- run_config("surv_viNULL",  args_surv(NULL))
res$cont_viNULL_mr <- run_config("cont_viNULL_mr", args_cont(NULL, mr = TRUE))   # side line: the driver's own call carries MR

# ---- the bootstrap check ---------------------------------------------------------
cat(sprintf("\n################ bootstrap check: nb_boots = %d, plan = sequential, cont fixture at vi.grf.min = %s\n", NB_BOOTS, deparse(vi_default)))
fs_est <- res$cont_viNULL$fs
stopifnot(!is.null(fs_est$sg.harm), "treat.recommend" %in% names(fs_est$df.est))
gc(); t0 <- proc.time()[[3]]
fb <- suppressWarnings(forestsearch_bootstrap_dofuture(fs_est, nb_boots = NB_BOOTS, seed = SEED_BASE + SIM_ID,
                                                       details = FALSE, parallel_args = list(plan = "sequential", workers = 1L)))
wall_fb <- proc.time()[[3]] - t0
per_rep <- wall_fb / NB_BOOTS
single <- res$cont_viNULL$wall_profiled
cat(sprintf("bootstrap wall %.1f s; per replicate %.2f s; profiled single call %.2f s; ratio per-rep / single = %.3f\n",
            wall_fb, per_rep, single, per_rep / single))
# per-replicate search time recorded by the bootstrap itself, when present
bt <- NULL
for (nm in c("boot_results", "results", "bootstrap_results")) if (!is.null(fb[[nm]]) && is.data.frame(fb[[nm]])) { bt <- fb[[nm]]; break }
if (is.null(bt)) { cand <- Filter(function(x) is.data.frame(x) && any(grepl("tmins", names(x))), fb); if (length(cand)) bt <- cand[[1]] }
if (!is.null(bt) && any(grepl("tmins_search", names(bt)))) {
  cat(sprintf("  bootstrap's own per-replicate search time: mean %.2f s, median %.2f s (tmins_search * 60); iteration mean %.2f s\n",
              60 * mean(bt$tmins_search, na.rm = TRUE), 60 * median(bt$tmins_search, na.rm = TRUE),
              if ("tmins_iteration" %in% names(bt)) 60 * mean(bt$tmins_iteration, na.rm = TRUE) else NA_real_))
  cat(sprintf("  replicates with an identified subgroup: %d of %d\n", sum(!is.na(bt$tmins_search)), nrow(bt)))
} else cat("  (no per-replicate timing table found in the bootstrap object; names: ", paste(names(fb), collapse = " "), ")\n")
cat(sprintf("projection, nb_boots = 1000, sequential: %.0f s = %.1f min from the profiled single call; %.0f s = %.1f min from the measured per-replicate mean\n",
            1000 * single, 1000 * single / 60, 1000 * per_rep, 1000 * per_rep / 60))

out <- list(results = lapply(res, function(r) r[setdiff(names(r), "fs")]),
            bootstrap = list(nb_boots = NB_BOOTS, wall = wall_fb, per_rep = per_rep, single = single, ratio = per_rep / single,
                             table = bt),
            vi_default = deparse(vi_default), fixtures = list(cont = list(n = 500L, seed = SEED_BASE + SIM_ID, confounders = confs_cont, payload_truth = truth),
                                                             surv = list(n = nrow(df_surv), data = "survival::gbsg", app = "quarto/applications/gbsg/analysis_gbsg_survival_effMaxSG.qmd")),
            built_at = Sys.time(), pkg_version = as.character(utils::packageVersion("forestsearch")),
            R_version = R.version.string, cores = parallel::detectCores())
saveRDS(out, file.path(D, "bootstrap_profile_2026-08-30.rds"))
cat("\nwritten:", file.path(D, "bootstrap_profile_2026-08-30.rds"), "\n")
