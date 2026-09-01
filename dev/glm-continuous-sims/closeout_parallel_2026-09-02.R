# =============================================================================
# assembly_settleA_2026-09-02.R -- section 6 read-only settlement: can
# forestsearch_bootstrap_dofuture() parallelize across replicates
# reproducibly?  (dev/tasks/cc_task_assembly_skip_2026-09-02.md)
# -----------------------------------------------------------------------------
# The continuous F1 bootstrap configuration (fixture block verbatim from
# medians_baseline_2026-09-02.R / assembly_battery_2026-09-02.R), nb_boots =
# 40, seed 8316952: once sequential, once under the wrapper's own outer plan
# (parallel_args = list(plan = "multisession", workers = 20) ->
# setup_parallel_SGcons() -> future::plan()).  Payloads pruned with the
# battery's own exclusion and compared with identical().  NO code change
# either way -- measurement only.
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

# ---- the measurement --------------------------------------------------------
run_boot40 <- function(label, pargs) {
  cat(sprintf("---- bootstrap 100 reps: %s\n", label))
  set.seed(999L)
  fs_est <- do.call(forestsearch, args_cont())
  t0 <- proc.time()[[3]]
  fb <- suppressWarnings(forestsearch_bootstrap_dofuture(fs_est, nb_boots = 100L, seed = 8316952L,
                                                         details = FALSE,
                                                         parallel_args = pargs))
  wall <- proc.time()[[3]] - t0
  cat(sprintf("  wall %.1f s; per replicate %.2f s\n", wall, wall / 100))
  list(label = label, wall = wall, per_rep = wall / 100, pruned = prune(fb, paste0("boot40_", label)))
}
seq40 <- run_boot40("sequential",   list(plan = "sequential",   workers = 1L,  show_message = FALSE))
par40 <- run_boot40("multisession48", list(plan = "multisession", workers = 48L, show_message = FALSE))
ok <- identical(seq40$pruned, par40$pruned)
cat(sprintf("\npayloads identical (sequential vs multisession/48): %s\n", ok))
if (!ok) { d <- first_diff(seq40$pruned, par40$pruned); cat("  FIRST DIFF at ", d$path, "\n"); print(d$a, digits = 22); print(d$b, digits = 22) }
cat(sprintf("speedup: %.1fx; implied B = 1000 wall at 48 workers: %.1f min (sequential: %.1f min)\n",
            seq40$wall / par40$wall, par40$wall / 100 * 1000 / 60,
            seq40$per_rep * 1000 / 60))
saveRDS(list(seq = seq40, par = par40, identical = ok,
             pkg_version = as.character(utils::packageVersion("forestsearch")),
             built_at = Sys.time()),
        file.path(D, "closeout_parallel_2026-09-02.rds"))
cat("written:", file.path(D, "closeout_parallel_2026-09-02.rds"), "\n")
