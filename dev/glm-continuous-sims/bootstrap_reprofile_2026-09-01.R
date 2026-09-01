# =============================================================================
# bootstrap_reprofile_2026-09-01.R -- re-profile the search at 0.3.2: where the
# bootstrap's cost sits now that the dispatch overhead is gone
# (task: dev/tasks/cc_task_bootstrap_reprofile_2026-09-01.md)
# -----------------------------------------------------------------------------
# Measurement only.  No R/ change, no optimisation.  Three fixtures:
#   cont   -- the MD40 DGM at n = 500, the drivers' arguments, replicate
#             configuration (plan = "sequential", workers = 1L -- at 0.3.2 this
#             takes the PLAIN LOOP: the 0.3.1 fix short-circuits on the plan
#             alone, subgroup_consistency_main.R Section 6).  Same fixture and
#             configuration as the 08-30 profile's cont_viNULL, so shares and
#             absolute seconds are like-for-like.
#   surv   -- survival::gbsg, the effMaxSG application's arguments, same
#             configuration (08-30's surv_viNULL).
#   anchor -- NEW, for scale only: the ACTG175 applied anchor, the fixed-family
#             maxeffCons call of quarto/applications/actg175/
#             analysis_actg175_continuous_oc.qmd section 2, verbatim
#             (parallel_args = list(plan = "sequential"), seed 8316951).
#             One profiled call; no bootstrap on it.
# Bucket definitions are the 08-30 script's BUCKETS list verbatim.  Added:
#   - a fit-bucket sub-split (solve = lm.fit/glm.fit/coxph.fit/agreg.fit vs
#     machinery = model.frame/model.matrix, vcov/summary, formula construction,
#     subset/assembly, other) for task section 4(a);
#   - medians-computed vs candidates-reaching-dedup counts for section 4(b);
#   - the bootstrap check at nb_boots = 20 (section 5);
#   - two timed multisession (workers = 6) calls on the continuous fixture,
#     default batch size vs one-batch-per-worker (section 6).
#
# Usage:  Rscript dev/glm-continuous-sims/bootstrap_reprofile_2026-09-01.R
# Writes (dev/glm-continuous-sims/): bootstrap_reprofile_2026-09-01.rds,
#   bootstrap_reprofile_<config>_2026-09-01.Rprof, and the log.
# =============================================================================
suppressPackageStartupMessages(devtools::load_all(".", quiet = TRUE))
stopifnot(requireNamespace("speff2trial", quietly = TRUE))
D <- "dev/glm-continuous-sims"
ns <- asNamespace("forestsearch")
vi_default <- formals(forestsearch)$vi.grf.min
cat(sprintf("installed forestsearch %s; vi.grf.min default = %s\n",
            as.character(utils::packageVersion("forestsearch")), deparse(vi_default)))
cat(sprintf("git HEAD: %s\n", system("git rev-parse --short HEAD", intern = TRUE)))
NB_BOOTS <- 20L

# ---- fixtures (cont + surv verbatim from bootstrap_profile_2026-08-30.R) ----
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
cat(sprintf("cont fixture: n = %d\n", nrow(df_cont)))
confs_cont <- c("age", "preanti", "wtkg", "karnof", "cd40", "cd80",
                "hemo", "homo", "drugs", "race", "gender", "symptom", "str2")
stopifnot(all(confs_cont %in% names(df_cont)), all(c("y_sim", "treat_sim", "id") %in% names(df_cont)))
args_cont <- function(vi, mr = FALSE, pa = list(plan = "sequential", workers = 1L), extra_pa = NULL) {
  if (!is.null(extra_pa)) pa <- c(pa, extra_pa)
  list(
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
}

df_surv <- survival::gbsg
df_surv <- within(df_surv, { id <- seq_len(nrow(df_surv)); time_months <- rfstime / 30.4375; grade3 <- ifelse(grade == "3", 1, 0) })
cat(sprintf("surv fixture: gbsg n = %d, events %d\n", nrow(df_surv), sum(df_surv$status)))
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

# anchor fixture: quarto/applications/actg175/analysis_actg175_continuous_oc.qmd
# data-prep (L54-64) and the anchor chunk (L88-120) verbatim; params$seed = 8316951
cont_vars <- c("age", "preanti", "wtkg", "karnof", "cd40", "cd80")
bin_vars  <- c("hemo", "homo", "drugs", "race", "gender", "symptom")
actg_oc <- subset(speff2trial::ACTG175, arms %in% c(1L, 3L))
actg_oc$id <- seq_len(nrow(actg_oc))
actg_oc$treat <- ifelse(actg_oc$arms == 1L, 1L, 0L)
actg_oc <- actg_oc[!is.na(actg_oc$cd420), ]
actg_oc$y_decline <- actg_oc$cd40 - actg_oc$cd420
for (v in bin_vars) actg_oc[[v]] <- as.factor(actg_oc[[v]])
stopifnot(nrow(actg_oc) == 1083L)
cat(sprintf("anchor fixture: ACTG175 arms 1/3, N = %d\n", nrow(actg_oc)))
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

# ---- instrumentation via trace() (08-30 verbatim) ---------------------------
.prof_env <- new.env()
reset_counters <- function() {
  .prof_env$fit_calls <- 0L; .prof_env$fit_times <- numeric(0); .prof_env$fit_null <- 0L
  .prof_env$cr_calls <- 0L; .prof_env$cr_na <- 0L
  .prof_env$split_calls <- 0L
  .prof_env$cons_calls <- 0L; .prof_env$cons_pass <- 0L
  .prof_env$dedup_before <- NA_integer_; .prof_env$dedup_after <- NA_integer_
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
# ---- section 4(a): sub-split of the FIT bucket's samples --------------------
# solve = the actual numerical fit; the rest is machinery, broken into the
# task's named parts.  Priority: solve first (a sample inside lm.fit is the
# solve even though model.frame is also on the stack), then the machinery
# parts in order.
FIT_SPLIT <- list(
  solve            = c("lm.fit", "glm.fit", "coxph.fit", "agreg.fit"),
  model_frame_matrix = c("model.frame", "model.frame.default", "model.matrix", "model.matrix.default",
                         ".getXlevels", "model.response", "model.extract", "na.omit", "makepredictcall"),
  vcov_summary     = c("vcov", "vcov.lm", "summary.lm", "summary.coxph", "summary", "chol2inv", "qr.coef",
                       "qr.lm"),
  formula_construction = c("as.formula", "formula", "formula.character", "terms", "terms.formula",
                           "str2lang", ".build_adjusted_formula", "Surv"),
  subset_assembly  = c("[.data.frame", "[[.data.frame", "[.data.table", "data.table", "cbind",
                       "as.data.frame", "data.frame"),
  wrapper_other    = character(0))   # remainder: lm/glm/coxph wrappers, tryCatch, coef, sqrt/diag
fit_split <- function(pb) {
  fit_stacks <- pb$stacks[pb$bucket_by_sample == "fit"]
  cls <- vapply(fit_stacks, function(st) {
    for (b in names(FIT_SPLIT)) if (length(FIT_SPLIT[[b]]) && any(st %in% FIT_SPLIT[[b]])) return(b)
    "wrapper_other" }, character(1))
  tab <- table(factor(cls, levels = names(FIT_SPLIT)))
  data.frame(part = names(tab), samples = as.integer(tab), secs = as.integer(tab) * pb$interval,
             share_of_fit = if (length(cls)) as.numeric(tab) / length(cls) else NA_real_,
             stringsAsFactors = FALSE)
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

# ---- one configuration ------------------------------------------------------
run_config <- function(name, args) {
  cat(sprintf("\n################ %s\n", name))
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
  cat(sprintf("  medians computed (survival fits returning non-NULL): %d; candidates reaching dedup: %s; ratio %.1f\n",
              cnt$fit_calls - cnt$fit_null, cnt$dedup_before,
              if (!is.na(cnt$dedup_before) && cnt$dedup_before > 0) (cnt$fit_calls - cnt$fit_null) / cnt$dedup_before else NA_real_))
  # warm-up (untrace() restores closures; JIT recompiles once)
  invisible(do.call(forestsearch, args))
  # profiled call -- 0.010 s is this platform's minimum sampling interval
  pf <- file.path(D, sprintf("bootstrap_reprofile_%s_2026-09-01.Rprof", name))
  gc(); Rprof(pf, interval = 0.010, line.profiling = TRUE, memory.profiling = FALSE)
  t0 <- proc.time()[[3]]
  fs_p <- do.call(forestsearch, args)
  wall_p <- proc.time()[[3]] - t0
  Rprof(NULL)
  stopifnot(identical(fs_p$sg.harm, fs_i$sg.harm))
  sr <- summaryRprof(pf, lines = "hide")
  pb <- parse_rprof(pf)
  fsplit <- fit_split(pb)
  cat(sprintf("profiled call: %.2f s wall; %d samples x %.3f s = %.2f s sampled (%.1f%% of wall)\n",
              wall_p, pb$n_samples, pb$interval, pb$total_secs, 100 * pb$total_secs / wall_p))
  cat("bucket table (share of sampled time):\n"); print(pb$table, digits = 4, row.names = FALSE)
  cat("fit-bucket sub-split (section 4a):\n"); print(fsplit, digits = 4, row.names = FALSE)
  cat("top 20 by SELF time:\n");  print(top_n(sr, "by.self"))
  cat("top 20 by TOTAL time:\n"); print(top_n(sr, "by.total"))
  list(name = name, wall_instrumented = wall_i, wall_profiled = wall_p, counters = cnt, buckets = pb$table,
       fit_split = fsplit,
       n_samples = pb$n_samples, interval = pb$interval, total_sampled = pb$total_secs,
       by_self = top_n(sr, "by.self", 30), by_total = top_n(sr, "by.total", 40),
       sg.harm = fs_p$sg.harm, fs = fs_p, rprof_file = pf, vi.grf.min = args$vi.grf.min)
}

res <- list()
res$cont   <- run_config("cont",   args_cont(NULL))   # replicate configuration; plain loop at 0.3.2
res$surv   <- run_config("surv",   args_surv(NULL))
res$anchor <- run_config("anchor", args_anchor)
cat(sprintf("\nanchor selected subgroup: %s (qmd gate expects {age <= 37} & !{cd40 <= 507}): %s\n",
            paste(res$anchor$sg.harm, collapse = " & "),
            setequal(res$anchor$sg.harm, c("{age <= 37}", "!{cd40 <= 507}"))))

# ---- then/now comparison against the 08-30 rds ------------------------------
old <- readRDS(file.path(D, "bootstrap_profile_2026-08-30.rds"))
cmp <- function(new, old_r, label) {
  ob <- old_r$buckets; nb <- new$buckets
  m <- merge(nb, ob, by = "bucket", suffixes = c("_now", "_0830"), sort = FALSE)
  m <- m[match(nb$bucket, m$bucket), c("bucket", "secs_0830", "secs_now", "share_0830", "share_now")]
  cat(sprintf("\n---- then/now buckets: %s (0.3.0 wall %.2f s -> 0.3.2 wall %.2f s)\n",
              label, old_r$wall_profiled, new$wall_profiled))
  print(m, digits = 4, row.names = FALSE); m
}
cmp_cont <- cmp(res$cont, old$results$cont_viNULL, "cont (replicate configuration)")
cmp_surv <- cmp(res$surv, old$results$surv_viNULL, "surv gbsg")

# ---- section 5: the bootstrap check ----------------------------------------
cat(sprintf("\n################ bootstrap check: nb_boots = %d, plan = sequential, cont fixture, vi.grf.min = NULL\n", NB_BOOTS))
fs_est <- res$cont$fs
stopifnot(!is.null(fs_est$sg.harm), "treat.recommend" %in% names(fs_est$df.est))
gc(); t0 <- proc.time()[[3]]
fb <- suppressWarnings(forestsearch_bootstrap_dofuture(fs_est, nb_boots = NB_BOOTS, seed = SEED_BASE + SIM_ID,
                                                       details = FALSE, parallel_args = list(plan = "sequential", workers = 1L)))
wall_fb <- proc.time()[[3]] - t0
per_rep <- wall_fb / NB_BOOTS
single <- res$cont$wall_profiled
cat(sprintf("bootstrap wall %.1f s; per replicate %.2f s; profiled single call %.2f s; ratio per-rep / single = %.3f\n",
            wall_fb, per_rep, single, per_rep / single))
bt <- NULL
for (nm in c("boot_results", "results", "bootstrap_results")) if (!is.null(fb[[nm]]) && is.data.frame(fb[[nm]])) { bt <- fb[[nm]]; break }
if (is.null(bt)) { cand <- Filter(function(x) is.data.frame(x) && any(grepl("tmins", names(x))), fb); if (length(cand)) bt <- cand[[1]] }
if (!is.null(bt) && any(grepl("tmins_search", names(bt)))) {
  cat(sprintf("  bootstrap's own per-replicate search time: mean %.2f s, median %.2f s; iteration mean %.2f s\n",
              60 * mean(bt$tmins_search, na.rm = TRUE), 60 * median(bt$tmins_search, na.rm = TRUE),
              if ("tmins_iteration" %in% names(bt)) 60 * mean(bt$tmins_iteration, na.rm = TRUE) else NA_real_))
  cat(sprintf("  replicates with an identified subgroup: %d of %d\n", sum(!is.na(bt$tmins_search)), nrow(bt)))
}
cat(sprintf("vs the narrow-dispatch task's 4.95 s per replicate: ratio %.3f\n", per_rep / 4.95))
cat(sprintf("projection, nb_boots = 1000, sequential: %.0f s = %.1f min (narrow-dispatch projected ~82 min)\n",
            1000 * per_rep, 1000 * per_rep / 60))

# ---- section 6: batch size under a real parallel plan ----------------------
n_screen <- res$cont$counters$cons_calls
batch_one_per_worker <- as.integer(ceiling(n_screen / 6))
cat(sprintf("\n################ batch-size measurement: multisession, workers = 6; n candidates at the screen = %d\n", n_screen))
time_call <- function(args) { gc(); t0 <- proc.time()[[3]]; fs <- do.call(forestsearch, args); c(wall = proc.time()[[3]] - t0, sg = paste(fs$sg.harm, collapse = " & ")) }
a_def <- args_cont(NULL, pa = list(plan = "multisession", workers = 6L))
a_big <- args_cont(NULL, pa = list(plan = "multisession", workers = 6L, batch_size = batch_one_per_worker))
r_def <- time_call(a_def)
r_big <- time_call(a_big)
rounds_def <- ceiling(n_screen / min(6L * 2L, n_screen))
rounds_big <- ceiling(n_screen / batch_one_per_worker)
cat(sprintf("default batch size (min(workers*2, n) = %d): wall %.2f s, %d future_lapply rounds; selected %s\n",
            min(12L, n_screen), as.numeric(r_def["wall"]), rounds_def, r_def["sg"]))
cat(sprintf("batch_size = %d (one batch per worker): wall %.2f s, %d future_lapply rounds; selected %s\n",
            batch_one_per_worker, as.numeric(r_big["wall"]), rounds_big, r_big["sg"]))

out <- list(results = lapply(res, function(r) r[setdiff(names(r), "fs")]),
            cmp = list(cont = cmp_cont, surv = cmp_surv),
            bootstrap = list(nb_boots = NB_BOOTS, wall = wall_fb, per_rep = per_rep, single = single,
                             ratio = per_rep / single, table = bt),
            batch = list(workers = 6L, n_screen = n_screen,
                         default = list(wall = as.numeric(r_def["wall"]), batch = min(12L, n_screen), rounds = rounds_def),
                         big = list(wall = as.numeric(r_big["wall"]), batch = batch_one_per_worker, rounds = rounds_big)),
            vi_default = deparse(vi_default),
            built_at = Sys.time(), pkg_version = as.character(utils::packageVersion("forestsearch")),
            git_head = system("git rev-parse --short HEAD", intern = TRUE),
            R_version = R.version.string, cores = parallel::detectCores())
saveRDS(out, file.path(D, "bootstrap_reprofile_2026-09-01.rds"))
cat("\nwritten:", file.path(D, "bootstrap_reprofile_2026-09-01.rds"), "\n")
