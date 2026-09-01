# =============================================================================
# fastpath_baseline_2026-09-02.R -- Stage A/C driver for the unadjusted
# fast-path task (dev/tasks/cc_task_fastpath_unadjusted_2026-09-02.md)
# -----------------------------------------------------------------------------
# Runs the full fixture battery against the INSTALLED forestsearch and saves
# pruned comparison sets.  The same script is both stages:
#   Stage A (0.3.2):  Rscript ... .R                    -> fastpath_baseline_2026-09-02.rds
#   Stage C (0.3.3):  FS_FASTPATH_OUT=postchange Rscript ... .R
#                                                       -> fastpath_postchange_2026-09-02.rds
# Fixtures:
#   F1 cont MD40 (replicate configuration) . F2 surv gbsg effMaxSG . F3 anchor
#     -- the three re-profile fixtures, blocks verbatim from
#        bootstrap_reprofile_2026-09-01.R (predicate TRUE expected)
#   F4 cont + adjust_covariates . F5 surv + adjust_covariates .
#   F6 cont + ps_method = "logistic", ps_adjust_method = "iptw"
#     -- routing fixtures (predicate FALSE expected)
#   E1-E6 -- synthetic edge fixtures (see task section 3)
# Plus a 20-replicate sequential bootstrap on F1 (seed 8316952).
#
# Volatile-field rule: the ONLY fields excluded from the saved comparison sets
# are wall-clock timings -- list fields `time_search` (find.grps) and
# `minutes_all` (top level), and bootstrap frame columns `tmins_search`,
# `tmins_iteration`, `tmins_all`.  data.tables are normalised to data.frame
# (their .internal.selfref external pointer differs across sessions; contents
# are compared in full).  Every dropped path is recorded in the rds.
# Self-consistency: F1, E1 and the F1 bootstrap run TWICE within this session
# and the two pruned sets must be identical() -- proof the exclusion list is
# complete, run in both stages so the call sequence (and any global RNG use)
# is identical across stages.
# =============================================================================
suppressPackageStartupMessages(library(forestsearch))
suppressPackageStartupMessages(library(survival))
stopifnot(requireNamespace("speff2trial", quietly = TRUE))
D <- "dev/glm-continuous-sims"
OUT <- Sys.getenv("FS_FASTPATH_OUT", "baseline")
ns <- asNamespace("forestsearch")
cat(sprintf("stage output: %s | installed forestsearch %s | vi.grf.min default %s | %s\n",
            OUT, as.character(utils::packageVersion("forestsearch")),
            deparse(formals(forestsearch)$vi.grf.min), R.version.string))
cat(sprintf("git HEAD: %s\n", system("git rev-parse --short HEAD", intern = TRUE)))

# ---- fixtures F1-F3: blocks verbatim from bootstrap_reprofile_2026-09-01.R --
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
confs_cont <- c("age", "preanti", "wtkg", "karnof", "cd40", "cd80",
                "hemo", "homo", "drugs", "race", "gender", "symptom", "str2")
stopifnot(all(confs_cont %in% names(df_cont)))
args_cont <- function(extra = list()) {
  base <- list(
    df.analysis = df_cont, confounders.name = confs_cont,
    outcome.name = "y_sim", treat.name = "treat_sim", id.name = "id",
    outcome_type = "continuous", effect_measure = "MD",
    effect.threshold = 30, consistency.threshold = 10, pconsistency.threshold = 0.90,
    fs.splits = 400L, n.min = 60L, d0.min = 12L, d1.min = 12L, maxk = 2L,
    vi.grf.min = NULL, sg_focus = "maxeffCons", selection_rule = "neighborhood",
    effect_neighborhood = 0.10, stop_threshold = NULL, consistency_method = "resample",
    conf.cont_jcuts = list(age = 10, preanti = 10),
    use_lasso = FALSE, use_dina = FALSE, use_grf = FALSE, use_twostage = TRUE,
    is.RCT = TRUE, adverse_outcome = FALSE, details = FALSE, quiet = TRUE,
    seedit = SEED_BASE + SIM_ID, parallel_args = list(plan = "sequential", workers = 1L),
    mr_inference = FALSE)
  utils::modifyList(base, extra)
}

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

# ---- edge fixtures E1-E6 ----------------------------------------------------
# E1: floor-boundary continuous.  The size gate is STRICT -- subgroup_search.R
# L613 `if (nx <= n.min) return(list(status = 4L, ...))` -- so a 75-slice at
# n.min = 75 would be excluded, not fitted.  n.min = 74 places the z1 == 1
# slice (n = 75) exactly at the smallest admissible size, which is the task's
# stated purpose ("fits at minimum size, both arms present").
set.seed(20260902)
e1_df <- data.frame(id = 1:150,
                    z1 = rep(0:1, each = 75),
                    z2 = sample(0:1, 150, TRUE))
e1_df$treat <- as.integer(unlist(lapply(split(seq_len(150), e1_df$z1),
                                        function(ix) rep_len(0:1, length(ix)))))
e1_df$y <- rnorm(150) + 0.8 * e1_df$z1 * e1_df$treat
args_e_cont <- function(df, extra = list()) {
  base <- list(
    df.analysis = df, confounders.name = c("z1", "z2"),
    outcome.name = "y", treat.name = "treat", id.name = "id",
    outcome_type = "continuous", effect_measure = "MD",
    effect.threshold = 0.1, consistency.threshold = 0.05, pconsistency.threshold = 0.5,
    fs.splits = 100L, n.min = 74L, maxk = 2L, vi.grf.min = NULL,
    sg_focus = "maxeffCons", selection_rule = "neighborhood",
    stop_threshold = NULL, consistency_method = "resample",
    use_lasso = FALSE, use_dina = FALSE, use_grf = FALSE, use_twostage = TRUE,
    is.RCT = TRUE, adverse_outcome = TRUE, details = FALSE, quiet = TRUE,
    seedit = 20260902, parallel_args = list(plan = "sequential", workers = 1L))
  utils::modifyList(base, extra)
}

# E2: zero-variance slice -- y constant within z1 == 1 (both arms), so the
# slice fit has rss = 0; summary.lm's "essentially perfect fit" warning path.
set.seed(20260903)
e2_df <- e1_df
e2_df$y <- ifelse(e2_df$z1 == 1, 5, rnorm(150))

# E3: rank-deficient slice -- z1 IS the treatment indicator, so the z1 == 1
# slice contains a single arm; lm(y ~ treat) gives an NA treat coefficient.
# Reachable on the continuous path: status 3 (d0/d1 floors) is skipped
# entirely for continuous outcomes (subgroup_search.R L601-609), and the only
# gate before the fit is `nx <= n.min` (L613).  n.min = 60 < 75 admits it.
set.seed(20260904)
e3_df <- data.frame(id = 1:150, treat = rep(0:1, each = 75),
                    z2 = sample(0:1, 150, TRUE))
e3_df$z1 <- e3_df$treat
e3_df$y <- rnorm(150)

# E4: heavy ties (survival): integer event times 1:6.
set.seed(20260905)
e4_df <- data.frame(id = 1:120, Y = sample(1:6, 120, TRUE),
                    E = rbinom(120, 1, 0.8), treat = rep_len(0:1, 120),
                    z1 = rep(0:1, each = 60), z2 = sample(0:1, 120, TRUE))
args_e_surv <- function(df, extra = list()) {
  base <- list(
    df.analysis = df, outcome.name = "Y", event.name = "E",
    treat.name = "treat", id.name = "id", confounders.name = c("z1", "z2"),
    is.RCT = TRUE, seedit = 20260902, est.scale = "hr",
    use_lasso = FALSE, use_grf = FALSE, use_dina = FALSE,
    subgroup_method = "consistency", sg_focus = "effMaxSG",
    selection_rule = "neighborhood",
    n.min = 40, d0.min = 5, d1.min = 5, maxk = 2,
    hr.threshold = 1.0, hr.consistency = 1.0, pconsistency.threshold = 0.5,
    minp = 0.025, use_twostage = TRUE, stop_threshold = NULL,
    consistency_method = "resample", fs.splits = 100, vi.grf.min = NULL,
    details = FALSE, quiet = TRUE, plot.sg = FALSE,
    parallel_args = list(plan = "sequential", workers = 1L))
  utils::modifyList(base, extra)
}

# E5: degenerate Cox fit that PASSES the floors.  True perfect separation (all
# events in one arm) is floor-unreachable: it forces d0 = 0 (or d1 = 0), and
# meets_event_criteria() (subgroup_search.R L707-709) requires
# `d0 >= d0.min && d1 >= d1.min` with d0.min, d1.min >= 1 -- quoted in the
# report.  The floor-passing degeneracy is monotone likelihood: in the z1 == 1
# slice every treated event precedes every control event, so the treatment
# coefficient diverges (coxph warns "Loglik converged before variable").
set.seed(20260906)
e5_z1 <- data.frame(id = 1:60, z1 = 1L, treat = rep(0:1, each = 30))
e5_z1$Y <- c(runif(20, 100, 110), rep(120, 10),    # controls: late events / censored
             runif(30, 1, 10))                     # treated: all early events
e5_z1$E <- c(rep(1L, 20), rep(0L, 10), rep(1L, 30))
e5_z0 <- data.frame(id = 61:130, z1 = 0L, treat = rep_len(0:1, 70),
                    Y = rexp(70) * 50, E = rbinom(70, 1, 0.7))
e5_df <- rbind(e5_z1, e5_z0)
e5_df$z2 <- sample(0:1, 130, TRUE)

# E6: boundary event count -- treated events in the z1 == 1 slice exactly at
# d1.min = 5; the >= gate admits it and the fit runs at the boundary.
set.seed(20260907)
e6_z1 <- data.frame(id = 1:60, z1 = 1L, treat = rep(0:1, each = 30))
e6_z1$E <- c(rbinom(30, 1, 0.6), rep(0L, 25), rep(1L, 5))  # control ~18 events; treated exactly 5
e6_z1$Y <- runif(60, 1, 100)
e6_z0 <- data.frame(id = 61:130, z1 = 0L, treat = rep_len(0:1, 70),
                    E = rbinom(70, 1, 0.7), Y = runif(70, 1, 100))
e6_df <- rbind(e6_z1, e6_z0)
e6_df$z2 <- sample(0:1, 130, TRUE)
stopifnot(sum(e6_df$E[e6_df$z1 == 1 & e6_df$treat == 1]) == 5L)

# ---- pruner -----------------------------------------------------------------
VOLATILE_FIELDS  <- c("time_search", "minutes_all")
VOLATILE_COLUMNS <- c("tmins_search", "tmins_iteration", "tmins_all")
.pe <- new.env(); .pe$dropped <- character(0)
prune <- function(x, path = "") {
  if (is.environment(x)) { .pe$dropped <- c(.pe$dropped, paste0(path, " <environment>")); return("<environment>") }
  if (is.function(x))    { .pe$dropped <- c(.pe$dropped, paste0(path, " <function>"));    return("<function>") }
  if (is.call(x))        { .pe$dropped <- c(.pe$dropped, paste0(path, " <call>")); return(paste(deparse(x), collapse = " ")) }
  if (inherits(x, "formula")) { .pe$dropped <- c(.pe$dropped, paste0(path, " <formula>")); return(paste(deparse(x), collapse = " ")) }
  if (inherits(x, "data.table")) x <- as.data.frame(x)
  if (is.data.frame(x)) {
    drop <- intersect(names(x), VOLATILE_COLUMNS)
    if (length(drop)) { .pe$dropped <- c(.pe$dropped, paste0(path, "$", drop, " [column]"))
                        x <- x[, setdiff(names(x), drop), drop = FALSE] }
    return(x)
  }
  if (is.list(x)) {
    nms <- names(x)
    if (!is.null(nms)) {
      kill <- nms %in% VOLATILE_FIELDS
      if (any(kill)) { .pe$dropped <- c(.pe$dropped, paste0(path, "$", nms[kill], " [field]"))
                       x <- x[!kill]; nms <- names(x) }
    }
    for (i in seq_along(x)) {
      lab <- if (!is.null(nms) && nzchar(nms[i])) nms[i] else as.character(i)
      x[i] <- list(prune(x[[i]], paste0(path, "$", lab)))
    }
    return(x)
  }
  x
}
first_diff <- function(a, b, path = "") {
  if (identical(a, b)) return(NULL)
  if (is.list(a) && is.list(b) && !is.data.frame(a) && !is.data.frame(b) &&
      identical(names(a), names(b)) && length(a) == length(b)) {
    for (i in seq_along(a)) {
      d <- first_diff(a[[i]], b[[i]], paste0(path, "$", if (!is.null(names(a))) names(a)[i] else i))
      if (!is.null(d)) return(d)
    }
    return(list(path = paste0(path, " [attributes differ]"), a = attributes(a), b = attributes(b)))
  }
  if (is.data.frame(a) && is.data.frame(b) && identical(names(a), names(b))) {
    for (cn in names(a)) if (!identical(a[[cn]], b[[cn]])) {
      j <- which(a[[cn]] != b[[cn]] | xor(is.na(a[[cn]]), is.na(b[[cn]])))[1]
      return(list(path = paste0(path, "$", cn, if (!is.na(j)) paste0("[", j, "]")),
                  a = if (!is.na(j)) a[[cn]][j] else a[[cn]], b = if (!is.na(j)) b[[cn]][j] else b[[cn]]))
    }
    return(list(path = paste0(path, " [frame attributes differ]"), a = attributes(a), b = attributes(b)))
  }
  list(path = path, a = a, b = b)
}

# ---- instrumented runner ----------------------------------------------------
reset_counters <- function() { .pe$fit_calls <- 0L; .pe$fit_null <- 0L; .pe$cons_calls <- 0L
                               .pe$dedup_before <- NA_integer_; .pe$dedup_after <- NA_integer_ }
TR <- character(0)
tr <- function(what, tracer = NULL, exit = NULL) {
  if (!exists(what, envir = ns, inherits = FALSE)) return(invisible(NULL))
  suppressMessages(trace(what, tracer = tracer, exit = exit, where = ns, print = FALSE))
  TR <<- c(TR, what); invisible(NULL)
}
set_traces <- function() {
  for (f in c("fit_glm_for_subgroup", "fit_cox_for_subgroup", "fit_glm_for_subgroup_fast")) tr(f,
    exit = quote({ .pe$fit_calls <- .pe$fit_calls + 1L
                   if (is.null(returnValue())) .pe$fit_null <- .pe$fit_null + 1L }))
  for (f in c("evaluate_subgroup_consistency", "evaluate_consistency_twostage")) tr(f,
    tracer = quote(.pe$cons_calls <- .pe$cons_calls + 1L))
  tr("remove_near_duplicate_subgroups",
     tracer = quote(.pe$dedup_before <- nrow(hr_subgroups)),
     exit = quote(.pe$dedup_after <- nrow(returnValue())))
}
clear_traces <- function() { for (f in unique(TR)) suppressMessages(untrace(f, where = ns)); TR <<- character(0) }

run_case <- function(name, args, seed_before) {
  cat(sprintf("---- %s\n", name))
  reset_counters(); set_traces()
  warns <- character(0)
  set.seed(seed_before)
  t0 <- proc.time()[[3]]
  fs <- withCallingHandlers(
    tryCatch(do.call(forestsearch, args), error = function(e) structure(list(error = conditionMessage(e)), class = "fs_error")),
    warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") })
  wall <- proc.time()[[3]] - t0
  clear_traces()
  counters <- list(fit_calls = .pe$fit_calls, fit_null = .pe$fit_null, cons_calls = .pe$cons_calls,
                   dedup_before = .pe$dedup_before, dedup_after = .pe$dedup_after)
  if (inherits(fs, "fs_error")) {
    cat(sprintf("  CALL ERRORED: %s\n", fs$error))
    return(list(name = name, error = fs$error, warnings = warns, counters = counters, wall = wall, pruned = NULL))
  }
  cat(sprintf("  wall %.2f s; selected: %s; fits %d (NULL %d); screen %d; dedup %s -> %s; warnings %d\n",
              wall, if (is.null(fs$sg.harm)) "<none>" else paste(fs$sg.harm, collapse = " & "),
              counters$fit_calls, counters$fit_null, counters$cons_calls,
              counters$dedup_before, counters$dedup_after, length(warns)))
  list(name = name, sg.harm = fs$sg.harm, warnings = warns, counters = counters, wall = wall,
       pruned = prune(fs, name))
}

res <- list()
res$F1 <- run_case("F1_cont",   args_cont(),   101L)
res$F2 <- run_case("F2_surv",   args_surv(),   102L)
res$F3 <- run_case("F3_anchor", args_anchor,   103L)
res$F4 <- run_case("F4_cont_adjusted", args_cont(list(adjust_covariates = c("age", "wtkg"))), 104L)
res$F5 <- run_case("F5_surv_adjusted", args_surv(list(adjust_covariates = "age")), 105L)
res$F6 <- run_case("F6_cont_iptw", args_cont(list(ps_method = "logistic", ps_adjust_method = "iptw")), 106L)
res$E1 <- run_case("E1_floor_boundary", args_e_cont(e1_df), 111L)
res$E2 <- run_case("E2_zero_variance",  args_e_cont(e2_df), 112L)
res$E3 <- run_case("E3_rank_deficient", args_e_cont(e3_df, list(n.min = 60L)), 113L)
res$E4 <- run_case("E4_heavy_ties",     args_e_surv(e4_df), 114L)
res$E5 <- run_case("E5_cox_degenerate", args_e_surv(e5_df), 115L)
res$E6 <- run_case("E6_floor_events",   args_e_surv(e6_df), 116L)

# ---- self-consistency: rerun F1 and E1, expect identical pruned sets --------
res_F1b <- run_case("F1_cont_selfcheck", args_cont(), 101L)
res_E1b <- run_case("E1_selfcheck", args_e_cont(e1_df), 111L)
for (pair in list(list("F1", res$F1, res_F1b), list("E1", res$E1, res_E1b))) {
  ok <- identical(pair[[2]]$pruned, pair[[3]]$pruned)
  cat(sprintf("self-consistency %s: %s\n", pair[[1]], ok))
  if (!ok) { d <- first_diff(pair[[2]]$pruned, pair[[3]]$pruned); cat("  FIRST DIFF at ", d$path, "\n"); print(d$a); print(d$b) }
}

# ---- F1 bootstrap (20 replicates, sequential), twice ------------------------
run_boot <- function(label) {
  cat(sprintf("---- bootstrap %s\n", label))
  set.seed(999L)
  # rebuild the F1 fit exactly (the bootstrap consumes the full fs object)
  fs_est <- do.call(forestsearch, args_cont())
  t0 <- proc.time()[[3]]
  fb <- suppressWarnings(forestsearch_bootstrap_dofuture(fs_est, nb_boots = 20L, seed = SEED_BASE + SIM_ID,
                                                         details = FALSE,
                                                         parallel_args = list(plan = "sequential", workers = 1L)))
  wall <- proc.time()[[3]] - t0
  cat(sprintf("  wall %.1f s; per replicate %.2f s\n", wall, wall / 20))
  list(wall = wall, per_rep = wall / 20, pruned = prune(fb, paste0("boot_", label)))
}
boot_a <- run_boot("first")
boot_b <- run_boot("selfcheck")
ok <- identical(boot_a$pruned, boot_b$pruned)
cat(sprintf("self-consistency bootstrap: %s\n", ok))
if (!ok) { d <- first_diff(boot_a$pruned, boot_b$pruned); cat("  FIRST DIFF at ", d$path, "\n"); print(d$a); print(d$b) }

out <- list(results = res, boot = boot_a, boot_selfcheck_ok = ok,
            selfcheck = list(F1 = identical(res$F1$pruned, res_F1b$pruned),
                             E1 = identical(res$E1$pruned, res_E1b$pruned)),
            dropped_paths = sort(unique(.pe$dropped)),
            pkg_version = as.character(utils::packageVersion("forestsearch")),
            git_head = system("git rev-parse --short HEAD", intern = TRUE),
            R_version = R.version.string, built_at = Sys.time())
saveRDS(out, file.path(D, sprintf("fastpath_%s_2026-09-02.rds", OUT)))
cat(sprintf("\ndropped paths (%d unique):\n", length(out$dropped_paths)))
cat(head(out$dropped_paths, 40), sep = "\n")
cat("\nwritten:", file.path(D, sprintf("fastpath_%s_2026-09-02.rds", OUT)), "\n")
