# =============================================================================
# closeout_planrepro_2026-09-02.R -- section 3 verification for the bootstrap
# close-out task (dev/tasks/cc_task_bootstrap_closeout_2026-09-02.md):
# direct-call plan reproducibility.
# -----------------------------------------------------------------------------
# Does a DIRECT forestsearch() call produce identical results under
# parallel_args = sequential/1 vs multisession/6 for the same seed?
# (Replicates always take the plain sequential loop, so this is direct-call
# reproducibility hygiene -- the recorded seeding-unification item.)
# Fixture blocks VERBATIM from the committed battery
# (medians_baseline_2026-09-02.R / assembly_battery_2026-09-02.R): F1
# continuous (seedit = 8316952 as the drivers set) and one F2 gbsg repeat so
# the verdict is stated per outcome type.  Pruner verbatim from the battery.
# NO fix, NO unification, NO default change here.
# =============================================================================
suppressPackageStartupMessages(library(forestsearch))
suppressPackageStartupMessages(library(survival))
stopifnot(requireNamespace("speff2trial", quietly = TRUE))
D <- "dev/glm-continuous-sims"
cat(sprintf("installed forestsearch %s | %s\n",
            as.character(utils::packageVersion("forestsearch")), R.version.string))
cat(sprintf("git HEAD: %s\n", system("git rev-parse --short HEAD", intern = TRUE)))

# ---- F1 fixture (verbatim from medians_baseline_2026-09-02.R) ---------------
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

# ---- pruner (verbatim from medians_baseline_2026-09-02.R) -------------------
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

# ---- run both plans per outcome type ----------------------------------------
run_plan <- function(name, args_fn, pargs, seed_before) {
  cat(sprintf("---- %s\n", name))
  set.seed(seed_before)
  t0 <- proc.time()[[3]]
  fs <- suppressWarnings(do.call(forestsearch, args_fn(list(parallel_args = pargs))))
  wall <- proc.time()[[3]] - t0
  cat(sprintf("  wall %.2f s; selected: %s\n", wall,
              if (is.null(fs$sg.harm)) "<none>" else paste(fs$sg.harm, collapse = " & ")))
  list(name = name, wall = wall, sg.harm = fs$sg.harm, pruned = prune(fs, name))
}
verdict <- function(label, a, b) {
  # The two calls differ in ONE argument by construction -- parallel_args --
  # and forestsearch() echoes its arguments verbatim in args_call_all, so that
  # echo is excluded before comparison (it is the experimental manipulation,
  # not a result).  Everything else, results included, must be identical().
  pa <- a$pruned; pb <- b$pruned
  pa$args_call_all$parallel_args <- NULL
  pb$args_call_all$parallel_args <- NULL
  ok <- identical(pa, pb)
  cat(sprintf("\n%s: sequential vs multisession/6 identical (parallel_args argument echo excluded): %s\n", label, ok))
  if (!ok) { d <- first_diff(pa, pb); cat("  FIRST DIFF at ", d$path, "\n")
             print(d$a, digits = 22); print(d$b, digits = 22) }
  ok
}

f1_seq <- run_plan("F1_sequential",    args_cont, list(plan = "sequential",   workers = 1L), 101L)
f1_par <- run_plan("F1_multisession6", args_cont, list(plan = "multisession", workers = 6L,  show_message = FALSE), 101L)
f2_seq <- run_plan("F2_sequential",    args_surv, list(plan = "sequential",   workers = 1L), 102L)
f2_par <- run_plan("F2_multisession6", args_surv, list(plan = "multisession", workers = 6L,  show_message = FALSE), 102L)

ok_f1 <- verdict("F1 continuous", f1_seq, f1_par)
ok_f2 <- verdict("F2 gbsg survival", f2_seq, f2_par)
saveRDS(list(F1 = list(seq = f1_seq, par = f1_par, identical = ok_f1),
             F2 = list(seq = f2_seq, par = f2_par, identical = ok_f2),
             pkg_version = as.character(utils::packageVersion("forestsearch")),
             built_at = Sys.time()),
        file.path(D, "closeout_planrepro_2026-09-02.rds"))
cat("\nwritten:", file.path(D, "closeout_planrepro_2026-09-02.rds"), "\n")
