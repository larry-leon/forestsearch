# =============================================================================
# run_gbsg_app_null.R -- realized declaration rate under a uniform marginal
# Cox HR of 0.75, GBSG application family
# (dev/tasks/TASK_gbsg_app_null_declaration_2026-09-23.md)
#
# One cell: n = 686, FS only, no MR.  Replicate trials drawn from a null
# (uniform-effect) AFT DGM built on survival::gbsg and calibrated to a
# super-population marginal Cox HR of 0.75, each run through forestsearch() at
# the application's own settings, counting the fraction that declare.
#
# Usage (from the package root):
#   FS_APPNULL_MODE=timing Rscript quarto/simulations/gbsg_app_null/run_gbsg_app_null.R
#   FS_APPNULL_MODE=full   Rscript quarto/simulations/gbsg_app_null/run_gbsg_app_null.R
# Optional: FS_APPNULL_WORKERS (default 14), FS_APPNULL_NSIMS (default 10 / 1000).
# =============================================================================

# forestsearch must be INSTALLED at HEAD's R/ (not load_all()): the doFuture
# multisession workers are separate R processes that see only the installed
# package.
suppressPackageStartupMessages({
  library(forestsearch)
  library(survival)
  library(foreach)
  library(doFuture)
})

mode      <- Sys.getenv("FS_APPNULL_MODE", "timing")
stopifnot(mode %in% c("timing", "full"))
n_workers <- as.integer(Sys.getenv("FS_APPNULL_WORKERS", "14"))
n_sims    <- as.integer(Sys.getenv("FS_APPNULL_NSIMS",
                                   if (mode == "timing") "10" else "1000"))
out_dir   <- "quarto/simulations/gbsg_app_null/results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- the application's settings -------------------------------------------
# Read from fs-glms-interpretable/quarto/gbsg/analysis_gbsg_mr.qmd (read-only):
# the frame at L151-166, the forestsearch() call at L219-233.  MR arguments are
# out of scope (TASK 7).
doc_seed  <- 8316951L          # seedit = 8316951 (L223)
n_sample  <- 686L              # stopifnot(nrow(df) == 686) (L164)

df <- survival::gbsg
df <- within(df, {
  id <- seq_len(nrow(df))
  time_months <- rfstime / 30.4375
  grade3 <- ifelse(grade == "3", 1, 0)
  treat <- hormon })
stopifnot(nrow(df) == 686, !anyNA(df), sum(df$treat == 1) == 246,
          sum(df$status) == 299)

confounders <- c("age", "meno", "size", "grade3", "nodes", "pgr", "er")  # L163

# Everything the application passes to forestsearch() except the frame, the
# column names (the simulated frame carries y_sim / event_sim / treat_sim / id)
# and MR.  Settings the call leaves implicit take forestsearch()'s formal
# default and are written out here, because default_sim_params() would
# otherwise substitute its own (use_twostage = FALSE, use_lasso = TRUE,
# d0.min = d1.min = 12, fs.splits = 400, hr.threshold = 1.25).
fs_params_app <- list(
  is.RCT = TRUE, est.scale = "hr",
  use_lasso = FALSE, use_grf = FALSE, use_dina = FALSE,
  max_n_confounders = 1000, conf_force = c("er <= 0", "pgr <= 0"),
  cut_type = "default", cont.cutoff = 4, conf.cont_jcuts = list(er = 10, pgr = 10),
  collapse_cuts = TRUE,
  subgroup_method = "consistency", sg_focus = "effMaxSG",
  selection_rule = "neighborhood", effect_neighborhood = 0.20,
  n.min = 60, d0.min = 10, d1.min = 10, maxk = 2,
  hr.threshold = 1.00, hr.consistency = 1.00, pconsistency.threshold = 0.90,
  consistency_method = "resample", m1.threshold = Inf, minp = 0.025,
  fs.splits = 1000,
  # implicit in the application's call -> forestsearch() formal defaults
  use_twostage = TRUE, outcome_type = "survival",
  mr_inference = FALSE,
  parallel_args = list(plan = "sequential")
)
# stop_threshold = NULL (L233) must reach forestsearch() as an explicit NULL:
# the run_fs path merges fs_params with modifyList(), which drops it (the
# formal default is pconsistency.threshold = 0.90); the `methods` path merges
# with .modify_keep_null(), which keeps it.  Hence methods = list(FS = list()).
fs_params_app["stop_threshold"] <- list(NULL)

# ---- the null DGM (TASK 2) -------------------------------------------------
base_args <- list(
  data            = df,
  continuous_vars = c("age", "size", "nodes", "pgr", "er"),
  factor_vars     = c("meno", "grade3"),
  outcome_var     = "time_months", event_var = "status", treatment_var = "hormon",
  model           = "null",
  n_super         = 5000L,
  seed            = doc_seed,
  verbose         = FALSE
)
k_treat <- calibrate_k_treat(target_hr_overall = 0.75, base_args = base_args,
                             use_ahr = FALSE, tol_rel = 1, verbose = FALSE)
dgm <- do.call(generate_aft_dgm_flex, c(base_args, list(k_treat = k_treat)))
hr_overall <- unname(dgm$hazard_ratios$overall)

# Gate 2 (STOP on failure)
gate2 <- c(hr_within_1pct = abs(hr_overall / 0.75 - 1) < 0.01,
           flag_harm_zero = all(dgm$df_super$flag_harm == 0))
cat(sprintf("k_treat = %.6f | overall marginal Cox HR = %.6f | AHR = %.6f\n",
            k_treat, hr_overall, dgm$hazard_ratios$AHR))
print(gate2)
if (!all(gate2)) stop("Gate 2 failed")

# ---- per-candidate null (finding, full mode only) --------------------------
# The family is enumerated on df_super by the steps of fs_oc_family_enumerate()
# sections 2-4 (R/fs_oc_family.R), which accepts glm_dgm objects only: cuts
# from get_FSdata() on df_super, both directions via dummy(), all <= maxk
# combinations, the empty / minp / rmin / size floors as population
# proportions at n = 686, identical memberships collapsed.  The event floors
# d0.min / d1.min are sample-event counts and have no population analogue;
# they are not applied.
enumerate_family <- function(dgm, n) {
  ds <- dgm$df_super
  ds$.y <- 0; ds$.e <- 1
  FSdata <- get_FSdata(
    df.analysis = ds, use_lasso = FALSE, use_grf = FALSE, grf_cuts = NULL,
    confounders.name = confounders, cont.cutoff = 4,
    conf_force = fs_params_app$conf_force, conf.cont_medians = NULL,
    conf.cont_medians_force = NULL, conf.cont_jcuts = fs_params_app$conf.cont_jcuts,
    dina_cuts = NULL, collapse_cuts = TRUE,
    collapse_cuts_args = eval(formals(forestsearch)$collapse_cuts_args),
    defaultcut_names = eval(formals(forestsearch)$defaultcut_names),
    cut_type = "default", exclude_cuts = eval(formals(forestsearch)$exclude_cuts),
    outcome.name = ".y", event.name = ".e", details = FALSE,
    outcome_type = "survival")
  Zdf <- dummy(FSdata$df[, FSdata$confs_names, drop = FALSE])
  Z <- as.matrix(Zdf); storage.mode(Z) <- "numeric"
  L <- ncol(Z)
  col_lab <- .fs_oc_column_labels(names(Zdf), FSdata$confs_names, FSdata$confs)
  maxk <- 2L; minp <- 0.025; rmin <- eval(formals(subgroup.search)[["rmin"]])
  combo <- generate_combination_indices(L, maxk)
  tot <- calculate_max_combinations(L, maxk)
  memb <- list(); lab <- character(0)
  cnt <- c(cut_columns = L, enumerated = 0L, empty = 0L, minp = 0L, rmin = 0L,
           size = 0L, kept = 0L)
  for (kk in seq_len(tot)) {
    ci <- get_covs_in(kk, maxk, L, combo$counts_1, combo$indices_1,
                      combo$counts_2, combo$indices_2,
                      combo$counts_3, combo$indices_3)
    if (sum(ci) < 1L || sum(ci) > maxk) next
    cnt["enumerated"] <- cnt["enumerated"] + 1L
    sel <- which(ci == 1); x <- Z[, sel, drop = FALSE]
    if (!has_positive_variance(x)) { cnt["empty"] <- cnt["empty"] + 1L; next }
    if (!meets_prevalence_threshold(x, minp)) { cnt["minp"] <- cnt["minp"] + 1L; next }
    if (.fs_oc_redundant(x, rmin / n)) { cnt["rmin"] <- cnt["rmin"] + 1L; next }
    m <- get_subgroup_membership(Z, ci)
    if (mean(m) < 60 / n) { cnt["size"] <- cnt["size"] + 1L; next }
    cnt["kept"] <- cnt["kept"] + 1L
    memb[[length(memb) + 1L]] <- m
    lab <- c(lab, paste(col_lab[sel], collapse = " & "))
  }
  mm <- do.call(cbind, memb)
  key <- apply(mm, 2L, function(v) paste(which(v), collapse = ","))
  first <- !duplicated(key)
  list(lab = lab[first], memb = mm[, first, drop = FALSE],
       counts = c(cnt, duplicate = sum(!first), M = sum(first)),
       cuts = FSdata$confs)
}

# The enumeration uses package internals; resolve them from the namespace.
environment(enumerate_family) <- list2env(
  list(confounders = confounders, fs_params_app = fs_params_app),
  parent = asNamespace("forestsearch"))

# True marginal Cox HR of each candidate, on the estimand the calibration uses
# (calculate_hazard_ratios(), R/generate_aft_dgm_helpers.R): every subject's
# two potential outcomes under a common extreme-value error, uncensored.  The
# error is drawn n_eps = 20 times per subject and stacked, so the per-candidate
# value carries little Monte Carlo noise from the draw.
candidate_truth <- function(dgm, fam, n_eps = 20L, seed = 20260923L) {
  ds <- dgm$df_super; mp <- dgm$model_params
  set.seed(seed)
  N <- nrow(ds)
  eps <- log(stats::rexp(N * n_eps))
  lp1 <- rep(ds$lin_pred_1, n_eps); lp0 <- rep(ds$lin_pred_0, n_eps)
  T1 <- exp(mp$mu + mp$tau * eps + lp1); T0 <- exp(mp$mu + mp$tau * eps + lp0)
  fit_hr <- function(idx) {
    ii <- rep(idx, n_eps)
    d <- data.frame(time = c(T1[ii], T0[ii]), treat = rep(1:0, each = sum(ii)))
    exp(unname(coef(coxph(Surv(time, rep(1, nrow(d))) ~ treat, data = d,
                          ties = "breslow"))))
  }
  list(overall = fit_hr(rep(TRUE, N)),
       hr_g = vapply(seq_len(ncol(fam$memb)),
                     function(j) fit_hr(fam$memb[, j]), numeric(1)))
}

# ---- one replicate (TASK 4) ------------------------------------------------
seed_base <- doc_seed
one_rep <- function(b) {
  t0 <- proc.time()[["elapsed"]]
  warns <- character(0)
  res <- tryCatch(
    withCallingHandlers(
      run_simulation_analysis(
        sim_id = b, dgm = dgm, n_sample = n_sample,
        confounders_base = confounders,
        methods = list(FS = list()),
        fs_params = utils::modifyList(fs_params_app, list(seedit = seed_base + b)),
        run_fs = TRUE, run_grf = FALSE, run_fs_grf = FALSE,
        seed_base = seed_base, verbose = FALSE),
      warning = function(w) {
        warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning")
      }),
    error = function(e) e)
  sec <- proc.time()[["elapsed"]] - t0
  fs_failed <- any(grepl("^FS failed:", warns))
  if (inherits(res, "error")) {
    return(data.frame(sim = b, seed = seed_base + b, errored = TRUE,
                      error_msg = conditionMessage(res), any.H = NA_integer_,
                      size.H = NA_integer_, sg.def = NA_character_,
                      p.cens = NA_real_, n_candidates_total = NA_integer_,
                      n_passed = NA_integer_, hr.itt = NA_real_,
                      warnings = paste(warns, collapse = " | "), seconds = sec))
  }
  r <- as.data.frame(res)
  data.frame(sim = b, seed = seed_base + b, errored = fs_failed,
             error_msg = if (fs_failed) paste(grep("^FS failed:", warns, value = TRUE),
                                              collapse = " | ") else NA_character_,
             any.H = as.integer(r$any.H), size.H = as.integer(r$size.H),
             sg.def = as.character(r$sg.def), p.cens = r$p.cens,
             n_candidates_total = as.integer(r$n_candidates_total),
             n_passed = as.integer(r$n_passed), hr.itt = r$hr.itt,
             warnings = paste(warns, collapse = " | "), seconds = sec)
}

# ---- run -------------------------------------------------------------------
cat(sprintf("mode = %s | n_sims = %d | workers = %d | n = %d\n",
            mode, n_sims, n_workers, n_sample))
plan("multisession", workers = n_workers)
t_run <- proc.time()[["elapsed"]]
rows <- foreach(b = seq_len(n_sims), .options.future = list(seed = TRUE,
                packages = c("forestsearch", "survival"))) %dofuture% one_rep(b)
wall <- proc.time()[["elapsed"]] - t_run
plan("sequential")
results <- do.call(rbind, rows)

cat(sprintf("wall-clock %.1f s | mean per-replicate %.1f s | errored %d\n",
            wall, mean(results$seconds), sum(results$errored)))
cat(sprintf("declared %d of %d\n", sum(results$any.H == 1L, na.rm = TRUE), n_sims))
proj_h <- wall / n_sims * 1000 / 3600
if (mode == "timing") cat(sprintf("projected 1,000-replicate wall-clock: %.2f h\n", proj_h))

payload <- list(
  task = "TASK_gbsg_app_null_declaration_2026-09-23", mode = mode,
  results = results, wall_seconds = wall, n_workers = n_workers,
  n_sims = n_sims, n_sample = n_sample, seed_base = seed_base,
  seeds = "simulate_from_dgm seed = seedit = seed_base + sim; DGM seed = seed_base",
  k_treat = k_treat, hazard_ratios = dgm$hazard_ratios, gate2 = gate2,
  fs_params = fs_params_app, confounders = confounders,
  git_head = system("git rev-parse HEAD", intern = TRUE),
  forestsearch_version = as.character(utils::packageVersion("forestsearch")),
  forestsearch_built = utils::packageDescription("forestsearch")$Built,
  R = R.version.string, platform = R.version$platform,
  gbsg = c(n = nrow(df), events = sum(df$status), event_rate = mean(df$status)))

if (mode == "full") {
  fam <- enumerate_family(dgm, n_sample)
  tru <- candidate_truth(dgm, fam)
  payload$family <- list(lab = fam$lab, counts = fam$counts, cuts = fam$cuts,
                         Pg = colMeans(fam$memb), hr_true = tru$hr_g,
                         hr_overall_check = tru$overall)
  cat("family counts:\n"); print(fam$counts)
  cat(sprintf("overall check %.4f | HR_true(g) min %.4f median %.4f max %.4f | > 0.75: %d of %d\n",
              tru$overall, min(tru$hr_g), median(tru$hr_g), max(tru$hr_g),
              sum(tru$hr_g > 0.75), length(tru$hr_g)))
}

saveRDS(payload, file.path(out_dir, sprintf("gbsg_app_null_%s.rds", mode)))
utils::write.csv(results, file.path(out_dir, sprintf("gbsg_app_null_%s.csv", mode)),
                 row.names = FALSE)
