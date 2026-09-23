# =============================================================================
# run_gbsg_pstar_grid.R -- p* calibration grid under a uniform marginal Cox HR
# of 0.75, on the GBSG application's own covariates and allocation
# (dev/tasks/TASK_gbsg_pstar_grid_2026-09-23.md)
#
# FS only, no MR.  The simulation calls simulate_from_dgm(baseline = "fixed")
# and forestsearch() directly (run_simulation_analysis() cannot carry the
# allocation arguments).  Each replicate is run at a floor p*, and the
# per-replicate max Pcons over the c1-surviving candidates gives the
# declaration rate at every p* >= the floor: rate(p*) = mean(max Pcons >= p*).
#
# Modes (FS_PSTAR_MODE):
#   gateA  -- DGM, calibration check on df_source, candidate family and
#             per-candidate truth (uncensored and under the DGM's censoring)
#   gateB  -- Cell 1, 200 replicates at p* = 0.90 and at p* = 0.50, same seeds;
#             exact-agreement check of the floor equivalence
#   gateC  -- Cell 1, 10 replicates at the floor on FS_PSTAR_WORKERS; wall
#             clock and 5,000-replicate projection; one replicate re-run
#             sequentially under L'Ecuyer-CMRG to check reproducibility
#   full   -- FS_PSTAR_CELL = 1 or 2, 5,000 replicates at the floor.
#             Requires FS_PSTAR_GO = 1 (Larry's compute go).
#
# Usage (from the package root):
#   FS_PSTAR_MODE=gateA Rscript quarto/simulations/gbsg_app_null/run_gbsg_pstar_grid.R
#   FS_PSTAR_MODE=gateB FS_PSTAR_WORKERS=48 Rscript ...
#   FS_PSTAR_MODE=gateC FS_PSTAR_WORKERS=48 Rscript ...
#   FS_PSTAR_MODE=full FS_PSTAR_CELL=1 FS_PSTAR_GO=1 FS_PSTAR_WORKERS=48 Rscript ...
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

mode <- Sys.getenv("FS_PSTAR_MODE", "gateA")
stopifnot(mode %in% c("gateA", "gateB", "gateC", "full"))
n_workers <- as.integer(Sys.getenv("FS_PSTAR_WORKERS", "48"))
out_dir   <- "quarto/simulations/gbsg_app_null/results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

pstar_floor <- 0.50
abort_sec   <- 2 * 3600     # hard abort for any full cell

# ---- the application's settings (verbatim from run_gbsg_app_null.R, e6477a22)
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
fs_params_app["stop_threshold"] <- list(NULL)
# ---- end of the inherited settings block ------------------------------------

# The 2026-09-23 run reached forestsearch() through run_simulation_analysis()'s
# methods path: default_sim_params() (column names y_sim / event_sim /
# treat_sim / id, and by.risk, vi.grf.min, max.minutes, twostage_args) merged
# NULL-preservingly with fs_params_app, filtered to forestsearch()'s formals,
# plus details = FALSE, plot.sg = FALSE, quiet = TRUE.  Rebuild the same
# argument list here so the only change is p* and the simulation call.
ns <- asNamespace("forestsearch")
fs_args_base <- ns$.modify_keep_null(ns$default_sim_params(), fs_params_app)
fs_args_base <- fs_args_base[names(fs_args_base) %in% names(formals(forestsearch))]
fs_args_base <- c(list(confounders.name = confounders, details = FALSE,
                       plot.sg = FALSE, quiet = TRUE), fs_args_base)
stopifnot("stop_threshold" %in% names(fs_args_base),
          is.null(fs_args_base$stop_threshold))

# Settings of TASK section 1 that must be visible, as received, in every gate
# replicate's args_call_all.  effect_measure = NULL is checked separately (it
# is not passed; forestsearch() resolves it from outcome_type).  sg_focus is
# passed as "effMaxSG" and stored in args_call_all in its canonical form
# "hrMaxSG" (forestsearch_main.R, the sg_focus normalisation sync).
expected_settings <- list(
  hr.threshold = 1.00, hr.consistency = 1.00, sg_focus = "hrMaxSG",
  selection_rule = "neighborhood", effect_neighborhood = 0.20,
  consistency_method = "resample", use_twostage = TRUE,
  outcome_type = "survival", maxk = 2, n.min = 60, d0.min = 10, d1.min = 10,
  fs.splits = 1000, conf.cont_jcuts = list(er = 10, pgr = 10),
  cont.cutoff = 4, conf_force = c("er <= 0", "pgr <= 0"),
  collapse_cuts = TRUE, m1.threshold = Inf, minp = 0.025, is.RCT = TRUE,
  est.scale = "hr", use_lasso = FALSE, use_grf = FALSE, use_dina = FALSE)

check_settings <- function(aca, pstar) {
  ok <- vapply(names(expected_settings), function(nm)
    isTRUE(all.equal(aca[[nm]], expected_settings[[nm]])), logical(1))
  c(ok,
    stop_threshold_null = "stop_threshold" %in% names(aca) &&
      is.null(aca$stop_threshold),
    pstar = isTRUE(all.equal(aca$pconsistency.threshold, pstar)))
}

# ---- the null DGM (as run_gbsg_app_null.R builds it) -----------------------
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
gateA_file <- file.path(out_dir, "gbsg_pstar_gateA.rds")
if (mode == "gateA") {
  k_treat_orig <- calibrate_k_treat(target_hr_overall = 0.75, base_args = base_args,
                                    use_ahr = FALSE, tol_rel = 1, verbose = FALSE)
  k_treat <- k_treat_orig
} else {
  gA <- readRDS(gateA_file)
  stopifnot(isTRUE(gA$gateA_pass))
  k_treat <- gA$k_treat_used
}
dgm <- do.call(generate_aft_dgm_flex, c(base_args, list(k_treat = k_treat)))
stopifnot(nrow(dgm$df_source) == 686L, sum(dgm$df_source$treat) == 246L)

# Both potential outcomes per subject under a common extreme-value error,
# drawn n_eps times per subject and stacked (the construction of
# calculate_hazard_ratios() and of the 2026-09-23 truth).  With cens = TRUE
# each stacked subject also carries censoring times drawn from the DGM's own
# censoring model (common censoring error across arms; analysis_time = Inf,
# so no administrative censoring), and the Cox fit uses the observed
# min(T, C) and T <= C.
po_frame <- function(ds, mp, n_eps, seed) {
  set.seed(seed)
  N <- nrow(ds)
  eps <- log(stats::rexp(N * n_eps))
  T1 <- exp(mp$mu + mp$tau * eps + rep(ds$lin_pred_1, n_eps))
  T0 <- exp(mp$mu + mp$tau * eps + rep(ds$lin_pred_0, n_eps))
  cm <- mp$censoring
  stopifnot(cm$type == "weibull")
  ec <- log(stats::rexp(N * n_eps))
  C1 <- exp(cm$mu + cm$tau * ec + rep(ds$lin_pred_cens_1, n_eps))
  C0 <- exp(cm$mu + cm$tau * ec + rep(ds$lin_pred_cens_0, n_eps))
  list(N = N, n_eps = n_eps, T1 = T1, T0 = T0, C1 = C1, C0 = C0)
}
po_hr <- function(po, idx = rep(TRUE, po$N), cens = FALSE) {
  ii <- rep(idx, po$n_eps)
  m <- sum(ii)
  if (cens) {
    tt <- c(pmin(po$T1[ii], po$C1[ii]), pmin(po$T0[ii], po$C0[ii]))
    ee <- c(po$T1[ii] <= po$C1[ii], po$T0[ii] <= po$C0[ii])
  } else {
    tt <- c(po$T1[ii], po$T0[ii]); ee <- rep(TRUE, 2L * m)
  }
  exp(unname(coef(coxph(Surv(tt, as.numeric(ee)) ~ rep(1:0, each = m),
                        ties = "breslow"))))
}
hr_source <- function(dgm, n_eps = 20L, seed = 20260923L)
  po_hr(po_frame(dgm$df_source, dgm$model_params, n_eps, seed))

# ---- candidate family on df_source (gateA) ---------------------------------
# Same route as run_gbsg_app_null.R's enumerate_family(), now on the observed
# 686 rows, so the floors are exact sample counts at n = 686.  Also returns
# the MR route's count (every <= maxk combination with membership >= n.min,
# keyed by label, no other floor and no membership collapse), which is how
# the application's family_size_prereduction = 1,744 was formed
# (forestsearch_main.R, MR block) -- for the diagnostic comparison.
enumerate_family <- function(ds, n) {
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
  n_mr_route <- 0L
  cnt <- c(cut_columns = L, enumerated = 0L, empty = 0L, minp = 0L, rmin = 0L,
           size = 0L, kept = 0L)
  for (kk in seq_len(tot)) {
    ci <- get_covs_in(kk, maxk, L, combo$counts_1, combo$indices_1,
                      combo$counts_2, combo$indices_2,
                      combo$counts_3, combo$indices_3)
    if (sum(ci) < 1L || sum(ci) > maxk) next
    cnt["enumerated"] <- cnt["enumerated"] + 1L
    if (sum(get_subgroup_membership(Z, ci)) >= 60) n_mr_route <- n_mr_route + 1L
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
       counts = c(cnt, duplicate = sum(!first), M = sum(first),
                  mr_route = n_mr_route),
       cuts = FSdata$confs)
}
environment(enumerate_family) <- list2env(
  list(confounders = confounders, fs_params_app = fs_params_app),
  parent = asNamespace("forestsearch"))

run_gateA <- function() {
  ds <- dgm$df_source
  hr_src <- hr_source(dgm)
  hr_src_seeds <- vapply(0:4, function(s) hr_source(dgm, seed = 20260923L + s),
                         numeric(1))
  k_used <- k_treat
  recal <- NULL
  if (abs(hr_src / 0.75 - 1) >= 0.01) {
    # Recalibrate against df_source: calibrate_k_treat()'s mechanism
    # (.calibrate_by_root on a rebuilt DGM per k), with the objective read on
    # df_source by the same potential-outcome construction.
    f <- function(k) hr_source(do.call(generate_aft_dgm_flex,
                                       c(base_args, list(k_treat = k)))) - 0.75
    recal <- stats::uniroot(f, interval = c(-5, 5), extendInt = "yes",
                            tol = 1e-6)
    k_used <- recal$root
    dgm <<- do.call(generate_aft_dgm_flex, c(base_args, list(k_treat = k_used)))
    ds <- dgm$df_source
    hr_src <- hr_source(dgm)
  }
  mp <- dgm$model_params
  d_lp <- ds$lin_pred_1 - ds$lin_pred_0
  hr_patient <- exp(-d_lp / mp$tau)
  gate <- c(hr_within_1pct = abs(hr_src / 0.75 - 1) < 0.01,
            flag_harm_zero_source = all(ds$flag_harm == 0),
            flag_harm_zero_super = all(dgm$df_super$flag_harm == 0))
  cat(sprintf("k_treat original = %.6f | used = %.6f | recalibrated = %s\n",
              k_treat_orig, k_used, !is.null(recal)))
  cat(sprintf("df_super overall marginal Cox HR = %.6f (DGM) | AHR = %.6f\n",
              dgm$hazard_ratios$overall, dgm$hazard_ratios$AHR))
  cat(sprintf("df_source marginal Cox HR (20 eps, seed 20260923) = %.6f\n", hr_src))
  cat(sprintf("  MC spread over seeds 20260923+0:4: %s\n",
              paste(sprintf("%.4f", hr_src_seeds), collapse = " ")))
  cat(sprintf("patient-level conditional HR: min %.6f max %.6f\n",
              min(hr_patient), max(hr_patient)))
  print(gate)
  if (!all(gate)) stop("Gate A failed")

  t0 <- proc.time()[["elapsed"]]
  fam <- enumerate_family(ds, n_sample)
  po <- po_frame(ds, mp, 20L, 20260923L)
  hr_unc <- vapply(seq_len(ncol(fam$memb)),
                   function(j) po_hr(po, fam$memb[, j], cens = FALSE), numeric(1))
  hr_cen <- vapply(seq_len(ncol(fam$memb)),
                   function(j) po_hr(po, fam$memb[, j], cens = TRUE), numeric(1))
  ov <- c(uncensored = po_hr(po), censored = po_hr(po, cens = TRUE))
  cens_rate_po <- 1 - mean(c(po$T1 <= po$C1, po$T0 <= po$C0))
  cat(sprintf("family + truth: %.1f s\n", proc.time()[["elapsed"]] - t0))
  cat("family counts:\n"); print(fam$counts)
  diff_pct <- 100 * (fam$counts[["mr_route"]] / 1744 - 1)
  cat(sprintf("MR-route count %d vs application 1,744: %+.1f%% | task-route M = %d (%+.1f%%)\n",
              fam$counts[["mr_route"]], diff_pct, fam$counts[["M"]],
              100 * (fam$counts[["M"]] / 1744 - 1)))
  summ <- function(h) c(min = min(h), median = stats::median(h), max = max(h),
                        n_above_075 = sum(h > 0.75), share_above_075 = mean(h > 0.75))
  truth_tab <- rbind(uncensored = summ(hr_unc), censored = summ(hr_cen))
  cat(sprintf("overall on df_source: uncensored %.4f | censored %.4f | PO censoring rate %.3f\n",
              ov[1], ov[2], cens_rate_po))
  print(round(truth_tab, 4))
  top <- order(-hr_cen)[1:5]
  print(data.frame(lab = fam$lab[top], prev = round(colMeans(fam$memb)[top], 3),
                   hr_unc = round(hr_unc[top], 4), hr_cen = round(hr_cen[top], 4)))

  saveRDS(list(
    task = "TASK_gbsg_pstar_grid_2026-09-23", mode = "gateA",
    gateA_pass = all(gate), gate = gate,
    k_treat_orig = k_treat_orig, k_treat_used = k_used, recalibration = recal,
    hr_super = dgm$hazard_ratios, hr_source = hr_src,
    hr_source_seeds = hr_src_seeds, hr_patient = unique(round(hr_patient, 10)),
    family = list(lab = fam$lab, counts = fam$counts, cuts = fam$cuts,
                  Pg = colMeans(fam$memb), hr_true_uncensored = hr_unc,
                  hr_true_censored = hr_cen, overall = ov,
                  po_censoring_rate = cens_rate_po),
    family_vs_app = c(app_prereduction = 1744, mr_route = fam$counts[["mr_route"]],
                      diff_pct = diff_pct),
    truth_table = truth_tab,
    fs_args_base = fs_args_base, confounders = confounders,
    git_head = system("git rev-parse HEAD", intern = TRUE),
    forestsearch_version = as.character(utils::packageVersion("forestsearch")),
    forestsearch_built = utils::packageDescription("forestsearch")$Built,
    R = R.version.string, platform = R.version$platform), gateA_file)
  if (abs(diff_pct) > 10) stop("family count differs from 1,744 by more than 10%")
  invisible(NULL)
}

# ---- one replicate ----------------------------------------------------------
cell_args <- function(cell) switch(as.character(cell),
  "1" = list(draw_treatment = TRUE, rand_ratio = 246 / 440),
  "2" = list(draw_treatment = FALSE, rand_ratio = 1),
  stop("unknown cell"))

one_rep <- function(b, cell, pstar, check = TRUE) {
  t0 <- proc.time()[["elapsed"]]
  ca <- cell_args(cell)
  df_b <- simulate_from_dgm(
    dgm = dgm, n = NULL, baseline = "fixed",
    draw_treatment = ca$draw_treatment, rand_ratio = ca$rand_ratio,
    analysis_time = Inf, cens_adjust = 0, seed = doc_seed + b)
  args <- fs_args_base
  args$df.analysis <- df_b
  args$pconsistency.threshold <- pstar
  args$seedit <- doc_seed + b
  warns <- character(0)
  fit <- tryCatch(
    withCallingHandlers(do.call(forestsearch, args),
      warning = function(w) {
        warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning")
      }),
    error = function(e) e)
  sec <- proc.time()[["elapsed"]] - t0
  err_warn <- grep("^Error|failed", warns, value = TRUE)
  errored <- inherits(fit, "error") || !is.null(fit$error_log) ||
    length(err_warn) > 0L
  err_msg <- if (inherits(fit, "error")) conditionMessage(fit) else
    if (!is.null(fit$error_log)) paste(unlist(fit$error_log), collapse = ": ") else
    if (length(err_warn)) paste(err_warn, collapse = " | ") else NA_character_
  gc_ <- if (!inherits(fit, "error")) fit$grp.consistency else NULL
  res <- gc_$out_sg$result
  declared <- !is.null(res) && nrow(res) > 0L
  set_ok <- if (check && !inherits(fit, "error"))
    check_settings(fit$args_call_all, pstar) else NA
  data.frame(
    sim = b, seed = doc_seed + b, cell = cell, pstar = pstar,
    errored = errored, error_msg = err_msg,
    declared = as.integer(declared),
    max_pcons = if (declared) max(as.numeric(res$Pcons)) else NA_real_,
    n_rows_out_sg = if (declared) nrow(res) else 0L,
    n_candidates_total = if (!is.null(gc_$n_candidates_total))
      as.integer(gc_$n_candidates_total) else NA_integer_,
    n_passed = if (!is.null(gc_$n_passed)) as.integer(gc_$n_passed) else NA_integer_,
    sg.def = if (declared && !is.null(fit$sg.harm))
      paste(fit$sg.harm, collapse = " & ") else "",
    n = nrow(df_b), n_treated = sum(df_b$treat_sim),
    events = sum(df_b$event_sim), event_rate = mean(df_b$event_sim),
    settings_ok = if (isTRUE(is.na(set_ok[1]))) NA else all(set_ok),
    settings_bad = if (isTRUE(is.na(set_ok[1]))) NA_character_ else
      paste(names(set_ok)[!set_ok], collapse = ","),
    effect_measure = if (!inherits(fit, "error"))
      paste(format(fit$args_call_all$effect_measure), collapse = "") else NA_character_,
    warnings = paste(warns, collapse = " | "), seconds = sec)
}

run_reps <- function(B, cell, pstar, workers) {
  plan("multisession", workers = workers)
  on.exit(plan("sequential"), add = TRUE)
  t_run <- proc.time()[["elapsed"]]
  rows <- foreach(b = seq_len(B), .options.future = list(seed = TRUE,
                  packages = c("forestsearch", "survival"))) %dofuture%
    one_rep(b, cell, pstar)
  wall <- proc.time()[["elapsed"]] - t_run
  list(results = do.call(rbind, rows), wall = wall)
}

meta <- function() list(
  git_head = system("git rev-parse HEAD", intern = TRUE),
  forestsearch_version = as.character(utils::packageVersion("forestsearch")),
  forestsearch_built = utils::packageDescription("forestsearch")$Built,
  R = R.version.string, platform = R.version$platform,
  n_workers = n_workers, k_treat = k_treat, fs_args_base = fs_args_base,
  seeds = "simulate_from_dgm seed = seedit = 8316951 + b; DGM seed = 8316951")

# ---- modes ------------------------------------------------------------------
if (mode == "gateA") run_gateA()

if (mode == "gateB") {
  B <- as.integer(Sys.getenv("FS_PSTAR_NSIMS", "200"))
  cat(sprintf("gateB | cell 1 | B = %d | workers = %d\n", B, n_workers))
  r90 <- run_reps(B, 1, 0.90, n_workers)
  cat(sprintf("p* = 0.90: wall %.1f s\n", r90$wall))
  r50 <- run_reps(B, 1, pstar_floor, n_workers)
  cat(sprintf("p* = 0.50: wall %.1f s\n", r50$wall))
  a <- r90$results; f <- r50$results
  stopifnot(identical(a$sim, f$sim))
  ind_floor <- as.integer(!is.na(f$max_pcons) & f$max_pcons >= 0.90)
  disagree <- which(ind_floor != a$declared)
  ncand_mismatch <- which(a$n_candidates_total != f$n_candidates_total |
                          is.na(a$n_candidates_total) != is.na(f$n_candidates_total))
  checks <- c(
    exact_agreement = length(disagree) == 0L,
    n_candidates_match = length(ncand_mismatch) == 0L,
    no_errors = !any(a$errored) && !any(f$errored),
    settings_ok = all(a$settings_ok) && all(f$settings_ok),
    same_data = identical(a$n_treated, f$n_treated) && identical(a$events, f$events))
  cat(sprintf("declared at 0.90: %d / %d | I(max Pcons >= 0.90) at floor: %d / %d\n",
              sum(a$declared), B, sum(ind_floor), B))
  cat(sprintf("floor: declared (any Pcons >= 0.50) %d / %d\n", sum(f$declared), B))
  cat("disagreeing replicates:", if (length(disagree)) disagree else "none", "\n")
  cat("n_candidates_total mismatches:",
      if (length(ncand_mismatch)) ncand_mismatch else "none", "\n")
  cat("settings_bad (0.90):", unique(a$settings_bad), "| (floor):", unique(f$settings_bad), "\n")
  cat("effect_measure in args_call_all:", unique(c(a$effect_measure, f$effect_measure)), "\n")
  print(checks)
  saveRDS(c(list(task = "TASK_gbsg_pstar_grid_2026-09-23", mode = "gateB",
                 results_090 = a, results_floor = f, wall_090 = r90$wall,
                 wall_floor = r50$wall, disagree = disagree,
                 ncand_mismatch = ncand_mismatch, checks = checks), meta()),
          file.path(out_dir, "gbsg_pstar_gateB.rds"))
  utils::write.csv(rbind(a, f), file.path(out_dir, "gbsg_pstar_gateB.csv"),
                   row.names = FALSE)
  if (!all(checks)) stop("Gate B failed")
}

if (mode == "gateC") {
  B <- as.integer(Sys.getenv("FS_PSTAR_NSIMS", "10"))
  cat(sprintf("gateC | cell 1 | floor p* = %.2f | B = %d | workers = %d\n",
              pstar_floor, B, n_workers))
  r <- run_reps(B, 1, pstar_floor, n_workers)
  res <- r$results
  proj <- r$wall / B * 5000
  cat(sprintf("wall-clock %.1f s | mean per-replicate %.1f s | errored %d\n",
              r$wall, mean(res$seconds), sum(res$errored)))
  cat(sprintf("projected 5,000-replicate wall-clock (x %d / %d): %.0f s = %.1f min\n",
              5000, B, proj, proj / 60))
  cat(sprintf("projection by per-replicate mean / workers: %.0f s\n",
              mean(res$seconds) * 5000 / n_workers))
  # Post-condition 7: replicate 1 re-run sequentially under L'Ecuyer-CMRG
  RNGkind("L'Ecuyer-CMRG")
  seq1 <- one_rep(1L, 1, pstar_floor)
  keep <- c("declared", "max_pcons", "n_candidates_total", "n_passed", "sg.def",
            "n_treated", "events")
  repro <- isTRUE(all.equal(seq1[keep], res[1, keep], check.attributes = FALSE))
  cat("sequential L'Ecuyer re-run of replicate 1 reproduces parallel:", repro, "\n")
  print(rbind(parallel = res[1, keep], sequential = seq1[keep]))
  saveRDS(c(list(task = "TASK_gbsg_pstar_grid_2026-09-23", mode = "gateC",
                 results = res, wall_seconds = r$wall, projected_5000_seconds = proj,
                 repro_seq_lecuyer = repro, seq_rerun = seq1), meta()),
          file.path(out_dir, "gbsg_pstar_gateC.rds"))
  if (proj > abort_sec) cat("WARNING: projection exceeds the 2 h abort\n")
}

if (mode == "full") {
  stopifnot(Sys.getenv("FS_PSTAR_GO") == "1")
  cell <- as.integer(Sys.getenv("FS_PSTAR_CELL"))
  stopifnot(cell %in% 1:2)
  B <- as.integer(Sys.getenv("FS_PSTAR_NSIMS", "5000"))
  cat(sprintf("full | cell %d | floor p* = %.2f | B = %d | workers = %d\n",
              cell, pstar_floor, B, n_workers))
  plan("multisession", workers = n_workers)
  t_run <- proc.time()[["elapsed"]]
  chunk <- max(n_workers * 4L, 1L)
  rows <- list()
  for (s in seq(1L, B, by = chunk)) {
    bb <- s:min(s + chunk - 1L, B)
    rows <- c(rows, foreach(b = bb, .options.future = list(seed = TRUE,
                    packages = c("forestsearch", "survival"))) %dofuture%
                one_rep(b, cell, pstar_floor))
    el <- proc.time()[["elapsed"]] - t_run
    cat(sprintf("  %d / %d done, %.0f s\n", max(bb), B, el))
    if (el > abort_sec) { plan("sequential"); stop("2 h hard abort") }
  }
  wall <- proc.time()[["elapsed"]] - t_run
  plan("sequential")
  res <- do.call(rbind, rows)
  cat(sprintf("wall-clock %.1f s | errored %d of %d\n", wall, sum(res$errored), B))
  saveRDS(c(list(task = "TASK_gbsg_pstar_grid_2026-09-23", mode = "full",
                 cell = cell, results = res, wall_seconds = wall), meta()),
          file.path(out_dir, sprintf("gbsg_pstar_cell%d.rds", cell)))
  utils::write.csv(res, file.path(out_dir, sprintf("gbsg_pstar_cell%d.csv", cell)),
                   row.names = FALSE)
}
