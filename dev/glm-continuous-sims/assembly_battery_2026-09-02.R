# =============================================================================
# medians_baseline_2026-09-02.R -- Stage A/C driver for the medians-on-survivors
# task (dev/tasks/cc_task_medians_on_survivors_2026-09-02.md)
# -----------------------------------------------------------------------------
# Runs the fixture battery against the INSTALLED forestsearch and saves pruned
# comparison sets.  The same script is both stages:
#   Stage A (0.3.3):  Rscript ... .R                    -> medians_baseline_2026-09-02.rds
#   Stage C (0.3.4):  FS_ASSEMBLY_OUT=postchange Rscript ... .R
#                                                       -> medians_postchange_2026-09-02.rds
# Fixtures (task section 3):
#   F2 surv gbsg effMaxSG (primary) . F5 surv gbsg + adjust_covariates = "age"
#   F1 cont MD40 . F3 anchor          -- untouched-path guards (m1 = m0 = NA)
#   E-ties   -- the fast-path E4 heavy-ties construction, verbatim
#   E-named  -- a PASSING subgroup whose control arm never reaches median
#               survival (m0 NA today; the moved call must reproduce it)
#   E-finite -- gbsg with m1.threshold = 60 months (the finite-filter consumer)
#   E-zero   -- E-ties data with hr.threshold = 100: effect screen passes nobody
#   E-degen  -- the fast-path E5 monotone-likelihood construction, verbatim
# Plus a 20-replicate sequential gbsg bootstrap (seed 8316952), run twice.
#
# Volatile-field rule (REUSED VERBATIM from fastpath_baseline_2026-09-02.R):
# the ONLY fields excluded from the saved comparison sets are wall-clock
# timings -- list fields `time_search` (find.grps) and `minutes_all` (top
# level), and bootstrap frame columns `tmins_search`, `tmins_iteration`,
# `tmins_all`.  data.tables are normalised to data.frame.  Every dropped path
# is recorded in the rds.
# Self-consistency: F2, E-ties and the gbsg bootstrap run TWICE within this
# session and the two pruned sets must be identical().
# Elision instrumentation: survival::survfit dispatches are counted per
# fixture via trace() on the survival namespace, with the innermost
# forestsearch caller recorded for attribution.  These counts are REPORTING
# ONLY -- they differ across stages by design and are not part of the
# equality comparison sets.
# =============================================================================
suppressPackageStartupMessages(library(forestsearch))
suppressPackageStartupMessages(library(survival))
stopifnot(requireNamespace("speff2trial", quietly = TRUE))
D <- "dev/glm-continuous-sims"
OUT <- Sys.getenv("FS_ASSEMBLY_OUT", "baseline")
ns <- asNamespace("forestsearch")
cat(sprintf("stage output: %s | installed forestsearch %s | vi.grf.min default %s | %s\n",
            OUT, as.character(utils::packageVersion("forestsearch")),
            deparse(formals(forestsearch)$vi.grf.min), R.version.string))
cat(sprintf("git HEAD: %s\n", system("git rev-parse --short HEAD", intern = TRUE)))

# ---- fixtures F1/F3: blocks verbatim from fastpath_baseline_2026-09-02.R ----
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

# ---- edge fixtures ----------------------------------------------------------
# E-ties: the fast-path E4 heavy-ties construction, verbatim (integer event
# times 1:6, n = 120).
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

# E-degen: the fast-path E5 monotone-likelihood construction, verbatim.  In the
# z1 == 1 slice every treated event precedes every control event, so the
# treatment coefficient diverges (coxph warns "Loglik converged before
# variable") -- HR is huge, so the candidate PASSES the hr.threshold = 1.0
# screen (verified from the Stage A hr.subgroups, see checks below).
set.seed(20260906)
e5_z1 <- data.frame(id = 1:60, z1 = 1L, treat = rep(0:1, each = 30))
e5_z1$Y <- c(runif(20, 100, 110), rep(120, 10),    # controls: late events / censored
             runif(30, 1, 10))                     # treated: all early events
e5_z1$E <- c(rep(1L, 20), rep(0L, 10), rep(1L, 30))
e5_z0 <- data.frame(id = 61:130, z1 = 0L, treat = rep_len(0:1, 70),
                    Y = rexp(70) * 50, E = rbinom(70, 1, 0.7))
e5_df <- rbind(e5_z1, e5_z0)
e5_df$z2 <- sample(0:1, 130, TRUE)

# E-named: a set where a PASSING subgroup's control arm never reaches median
# survival.  z1 == 1 slice (n = 80): controls have 6 events among 40 (KM never
# drops below 0.5 -> m0 is NA), treated have 36 events among 40, spread over
# times overlapping the control events (no monotone likelihood).  HR >> 1
# clears hr.threshold = 1.0; d0 = 6 >= 5, d1 = 36 >= 5, n = 80 > 40 clear the
# floors.  The moved computation must reproduce the identical NA and the dedup
# key must paste the identical key string.
set.seed(20260908)
en_z1 <- data.frame(id = 1:80, z1 = 1L, treat = rep(0:1, each = 40))
en_z1$Y <- c(runif(6, 30, 100), rep(100, 34),        # controls: 6 events, 34 censored at 100
             runif(36, 5, 80), runif(4, 80, 100))    # treated: 36 events, 4 censored
en_z1$E <- c(rep(1L, 6), rep(0L, 34), rep(1L, 36), rep(0L, 4))
en_z0 <- data.frame(id = 81:150, z1 = 0L, treat = rep_len(0:1, 70),
                    Y = rexp(70) * 40, E = rbinom(70, 1, 0.7))
en_df <- rbind(en_z1, en_z0)
en_df$z2 <- sample(0:1, 150, TRUE)
stopifnot(sum(en_df$E[en_df$z1 == 1 & en_df$treat == 0]) == 6L)

# E-finite: the gbsg configuration with m1.threshold set finite (60 months) --
# no live driver exercises this consumer; the Stage A checks below verify the
# value actually excludes at least one row at 0.3.3.
M1_THRESHOLD_FINITE <- 60

# ---- pruner (verbatim from fastpath_baseline_2026-09-02.R) ------------------
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
                               .pe$dedup_before <- NA_integer_; .pe$dedup_after <- NA_integer_
                               .pe$survfit_calls <- 0L; .pe$survfit_by <- character(0) }
.pe$ns_fns <- ls(ns, all.names = TRUE)
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
  # Elision instrumentation: count survival::survfit dispatches and record the
  # innermost forestsearch function on the call stack for attribution.
  suppressMessages(trace("survfit", where = asNamespace("survival"), print = FALSE,
    tracer = quote({
      .pe$survfit_calls <- .pe$survfit_calls + 1L
      fnames <- sub("^.*::", "", vapply(sys.calls(), function(cl) {
        f <- cl[[1]]
        if (is.name(f)) as.character(f) else paste(deparse(f, nlines = 1L), collapse = "")
      }, character(1)))
      hit <- rev(fnames[fnames %in% .pe$ns_fns])
      .pe$survfit_by <- c(.pe$survfit_by, if (length(hit)) hit[1] else "<non-forestsearch>")
    })))
}
clear_traces <- function() {
  for (f in unique(TR)) suppressMessages(untrace(f, where = ns)); TR <<- character(0)
  suppressMessages(try(untrace("survfit", where = asNamespace("survival")), silent = TRUE))
}

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
  survfit_count <- list(n = .pe$survfit_calls, by = table(.pe$survfit_by))
  if (inherits(fs, "fs_error")) {
    cat(sprintf("  CALL ERRORED: %s\n", fs$error))
    return(list(name = name, error = fs$error, warnings = warns, counters = counters,
                survfit = survfit_count, wall = wall, pruned = NULL, hr_subgroups = NULL))
  }
  cat(sprintf("  wall %.2f s; selected: %s; fits %d (NULL %d); screen %d; dedup %s -> %s; warnings %d; survfit %d\n",
              wall, if (is.null(fs$sg.harm)) "<none>" else paste(fs$sg.harm, collapse = " & "),
              counters$fit_calls, counters$fit_null, counters$cons_calls,
              counters$dedup_before, counters$dedup_after, length(warns), survfit_count$n))
  hr_sg <- tryCatch(as.data.frame(fs$find.grps$out.found$hr.subgroups), error = function(e) NULL)
  list(name = name, sg.harm = fs$sg.harm, warnings = warns, counters = counters,
       survfit = survfit_count, wall = wall, pruned = prune(fs, name), hr_subgroups = hr_sg)
}

res <- list()
res$F2      <- run_case("F2_surv_gbsg",     args_surv(),   102L)
res$F5      <- run_case("F5_surv_adjusted", args_surv(list(adjust_covariates = "age")), 105L)
res$F1      <- run_case("F1_cont",          args_cont(),   101L)
res$F3      <- run_case("F3_anchor",        args_anchor,   103L)
res$Eties   <- run_case("Eties_heavy_ties", args_e_surv(e4_df), 114L)
res$Enamed  <- run_case("Enamed_na_median", args_e_surv(en_df), 121L)
res$Efinite <- run_case("Efinite_m1_60",    args_surv(list(m1.threshold = M1_THRESHOLD_FINITE)), 122L)
res$Ezero   <- run_case("Ezero_no_passers", args_e_surv(e4_df, list(hr.threshold = 100)), 123L)
res$Edegen  <- run_case("Edegen_monotone",  args_e_surv(e5_df), 115L)

# ---- fixture-validity checks (Stage A design premises, printed both stages) -
chk <- function(label, ok) cat(sprintf("check %-28s %s\n", label, if (isTRUE(ok)) "OK" else "** FAILED **"))
sg_en <- res$Enamed$hr_subgroups
chk("Enamed: NA median in a passer",
    !is.null(sg_en) && nrow(sg_en) > 0 && any(is.na(sg_en$m1) | is.na(sg_en$m0)))
sg_fin <- res$Efinite$hr_subgroups
n_excl <- if (!is.null(sg_fin)) sum(sg_fin$HR >= 1.0 & (is.na(sg_fin$m1) | sg_fin$m1 > M1_THRESHOLD_FINITE)) else NA_integer_
n_kept <- if (!is.null(sg_fin)) sum(sg_fin$HR >= 1.0 & !is.na(sg_fin$m1) & sg_fin$m1 <= M1_THRESHOLD_FINITE) else NA_integer_
cat(sprintf("Efinite: m1.threshold = %s excludes %s screen-passing rows, keeps %s\n",
            M1_THRESHOLD_FINITE, n_excl, n_kept))
chk("Efinite: filter excludes >= 1", isTRUE(n_excl >= 1L))
chk("Efinite: filter keeps >= 1",    isTRUE(n_kept >= 1L))
sg_z <- res$Ezero$hr_subgroups
chk("Ezero: screen passes nobody", is.null(sg_z) || nrow(sg_z) == 0)
sg_d <- res$Edegen$hr_subgroups
chk("Edegen: degenerate passer present",
    !is.null(sg_d) && nrow(sg_d) > 0 && any(sg_d$HR > 100))
chk("F1: continuous medians all NA",
    is.null(res$F1$hr_subgroups) || all(is.na(res$F1$hr_subgroups$m1) & is.na(res$F1$hr_subgroups$m0)))
chk("F1/F3: zero survfit calls",
    res$F1$survfit$n == 0L && res$F3$survfit$n == 0L)

# ---- self-consistency: rerun F2 and E-ties, expect identical pruned sets ----
res_F2b <- run_case("F2_selfcheck",    args_surv(),        102L)
res_Etb <- run_case("Eties_selfcheck", args_e_surv(e4_df), 114L)
for (pair in list(list("F2", res$F2, res_F2b), list("Eties", res$Eties, res_Etb))) {
  ok <- identical(pair[[2]]$pruned, pair[[3]]$pruned)
  cat(sprintf("self-consistency %s: %s\n", pair[[1]], ok))
  if (!ok) { d <- first_diff(pair[[2]]$pruned, pair[[3]]$pruned); cat("  FIRST DIFF at ", d$path, "\n"); print(d$a); print(d$b) }
}

# ---- gbsg bootstrap (20 replicates, sequential), twice ----------------------
run_boot <- function(label) {
  cat(sprintf("---- bootstrap %s\n", label))
  set.seed(999L)
  # rebuild the F2 fit exactly (the bootstrap consumes the full fs object)
  fs_est <- do.call(forestsearch, args_surv())
  t0 <- proc.time()[[3]]
  fb <- suppressWarnings(forestsearch_bootstrap_dofuture(fs_est, nb_boots = 20L, seed = 8316952L,
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

# hr_subgroups frames are extracted for the checks above but are already part
# of the pruned sets; keep the standalone copies out of the gated object so the
# comparison surface stays exactly the pruned fs object + warnings + counters.
res_gated <- lapply(res, function(r) r[setdiff(names(r), c("hr_subgroups", "survfit", "wall"))])
res_report <- lapply(res, function(r) r[c("name", "wall", "survfit")])

out <- list(results = res_gated, report = res_report,
            boot = boot_a, boot_selfcheck_ok = ok,
            selfcheck = list(F2 = identical(res$F2$pruned, res_F2b$pruned),
                             Eties = identical(res$Eties$pruned, res_Etb$pruned)),
            efinite = list(threshold = M1_THRESHOLD_FINITE, n_excluded = n_excl, n_kept = n_kept),
            dropped_paths = sort(unique(.pe$dropped)),
            pkg_version = as.character(utils::packageVersion("forestsearch")),
            git_head = system("git rev-parse --short HEAD", intern = TRUE),
            R_version = R.version.string, built_at = Sys.time())
saveRDS(out, file.path(D, sprintf("assembly_%s_2026-09-02.rds", OUT)))
cat(sprintf("\ndropped paths (%d unique):\n", length(out$dropped_paths)))
cat(head(out$dropped_paths, 40), sep = "\n")
cat("\nwritten:", file.path(D, sprintf("assembly_%s_2026-09-02.rds", OUT)), "\n")
