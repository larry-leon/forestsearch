#!/usr/bin/env Rscript
# =========================================================================
# Phase 1 acceptance gate: run_subgroup_sims() must be BIT-IDENTICAL to
# the committed extreme-subgroups vignettes (commit 002f4f37) for both
# baseline designs, and summary() must reproduce the Section 6.6.1
# results table to the character.
#
# All blocks marked "verbatim" below are machine-extracted from the two
# vignette .qmd files -- do not edit them here; edit the vignettes and
# regenerate if they ever change.
#
# Usage:   Rscript dev/accept_phase1_bitident.R [n_sims]     # default 50
# Requires: forestsearch installed WITH the Phase 1 wrapper files.
# =========================================================================
suppressMessages({ library(forestsearch); library(survival) })

args  <- commandArgs(trailingOnly = TRUE)
n_acc <- if (length(args) >= 1) as.integer(args[1]) else 50L

.fails <- 0L
ok <- function(cond, msg) {
  cat(sprintf("[%s] %s\n", if (isTRUE(cond)) "PASS" else "FAIL", msg))
  if (!isTRUE(cond)) .fails <<- .fails + 1L
  invisible(cond)
}

cat("forestsearch version:", as.character(packageVersion("forestsearch")), "\n")
if (!exists("run_subgroup_sims")) {
  stop("run_subgroup_sims() not found -- install the package with the ",
       "Phase 1 files (Gate 0) before running this script.")
}
cat("Acceptance n_sims  :", n_acc, "per design\n\n")
t_all <- Sys.time()

# ---- data-prep (verbatim, shared by both vignettes) ---------------------
# Load GBSG data and convert time to months (DGM time scale)
data(gbsg, package = "survival")

gbsg$time_months <- gbsg$rfstime / 30.4375
gbsg$treat        <- gbsg$hormon

# Kaplan–Meier median event time.  The naive median over events only —
# median(time[event == 1]) — is biased downward in the presence of
# censoring, since the KM curve accounts for the at-risk set at each
# event time.  If the KM curve does not cross 0.5 the median is
# returned as NA (event time "not reached").
km_overall <- survfit(Surv(time_months, status) ~ 1, data = gbsg)
km_median  <- summary(km_overall)$table[["median"]]

cat("N =", nrow(gbsg), "  Events:", sum(gbsg$status),
    sprintf("(%.0f%%)", 100 * mean(gbsg$status)), "\n")
cat("Median follow-up, (months):",
    round(median(gbsg$time_months), 1), "\n")
cat("Median event time, KM (months)  :",
    if (is.na(km_median)) "not reached" else round(km_median, 1), "\n")

# ---- per-design sim_config lists (verbatim values from each vignette) ---
cfg_res <- list(
  n_per_trial      = nrow(gbsg),
  analysis_time    = 84,
  max_entry        = 24,
  seed_base        = 0L,
  rand_seed_offset = 1e6L
)
cfg_fix <- list(
  baseline         = "fixed",
  analysis_time    = 84,
  max_entry        = 24,
  seed_base        = 0L,
  rand_seed_offset = 1e6L
)

# ---- build-dgm (verbatim, fixed-doc version; identical construction) ----
.t_dgm_start <- Sys.time()

# ── All arguments to generate_aft_dgm_flex() except k_treat ───────────────
# k_treat is the calibration parameter; everything else is fixed up front.
base_args_uniform <- list(
  # ── Required: source data & outcome ─────────────────────────────────────
  data                = gbsg,
  continuous_vars     = c("age", "size", "nodes", "pgr", "er"),
  factor_vars         = c("meno", "grade"),
  outcome_var         = "time_months",
  event_var           = "status",
  # ── Treatment ───────────────────────────────────────────────────────────
  treatment_var       = "hormon",
  draw_treatment      = FALSE,                                  # default
  # ── HTE specification — uniform treatment effect, no embedded subgroup ──
  model               = "null",                                 # uniform HR throughout
  subgroup_vars       = NULL,                                   # default — no embedded HTE
  subgroup_cuts       = NULL,                                   # default
  k_inter             = 1,                                      # default (irrelevant under "null")
  spline_spec         = NULL,                                   # default
  set_beta_spec       = list(set_var = NULL, beta_var = NULL),  # default
  # ── Censoring model (independent of treatment effect) ───────────────────
  continuous_vars_cens = NULL,                                  # default → inherit
  factor_vars_cens     = NULL,                                  # default → inherit
  select_censoring     = TRUE,                                  # default — AIC select
  cens_type            = "weibull",                             # default
  cens_params          = list(),                                # default
  cens_intercept_only  = FALSE,                                 # default
  # ── Super population & misc ─────────────────────────────────────────────
  n_super             = 5000,
  standardize         = FALSE,                                  # default
  seed                = 99,
  verbose             = FALSE
)

# ── Step 1: calibrate k_treat so that the AHR is exactly 0.70 ─────────────
# We calibrate to the AHR (average hazard ratio over the super-population),
# not the marginal Cox HR, because the simulation analyses below fit a Cox
# model stratified by `grade` — and under Cox non-collapsibility the
# marginal Cox HR on fully-observed potential outcomes is attenuated toward
# 1 relative to both the AHR and the stratified Cox HR.  Calibrating to the
# AHR makes the simulation median HR track 0.70 closely; calibrating to the
# marginal Cox HR (`use_ahr = FALSE`, the default) would pull the simulation
# median below 0.70.  See `?calibrate_k_treat` for both options.
k_treat_uniform <- calibrate_k_treat(
  target_hr_overall = 0.70,
  base_args         = base_args_uniform,
  k_treat_range     = c(-5, 5),
  tol               = 1e-6,
  use_ahr           = TRUE,
  verbose           = TRUE
)

# ── Step 2: build the DGM with the calibrated k_treat ─────────────────────
dgm_uniform <- do.call(
  generate_aft_dgm_flex,
  c(base_args_uniform, list(k_treat = k_treat_uniform))
)

t_dgm <- as.numeric(difftime(Sys.time(), .t_dgm_start, units = "secs"))

cat("Uniform DGM (AHR-calibrated)\n",
    "  k_treat calibrated to       : ", round(k_treat_uniform, 4), "\n",
    "  AHR (target 0.70)           : ",
    round(dgm_uniform$hazard_ratios$AHR, 4), "\n",
    "  marginal Cox HR (attenuated): ",
    round(dgm_uniform$hazard_ratios$overall, 4), "\n",
    "  implied median event T      : ",
    round(exp(dgm_uniform$model_params$mu), 1), " months\n", sep = "")

# ── Hard stop: fixed-baseline mode requires dgm$df_source ─────────────────
# The whole vignette simulates via baseline = "fixed", which reads the
# source frame stored by the constructor (forestsearch versions that
# include Step 7b).  Fail here, at build time, rather than 1,000 lines
# later inside the parallel loop.
if (is.null(dgm_uniform$df_source)) {
  stop("dgm_uniform$df_source is missing: this forestsearch version does ",
       "not store the fixed-baseline frame. Update the package ",
       "(baseline = \"fixed\" support) and rebuild the DGM.")
}
cat("Fixed-baseline panel  : nrow(df_source) =",
    nrow(dgm_uniform$df_source), "patients (the GBSG trial data)\n")

# ---- calibrate cens_adjust, random-X design (verbatim, random doc) ------
sim_config <- cfg_res
.t_cal_start <- Sys.time()

cal_uniform <- calibrate_cens_adjust(
  dgm           = dgm_uniform,
  target        = "rate",                       # "rate" | "km_median"
  n             = 1000,
  rand_ratio    = 1,                            # default
  analysis_time = sim_config$analysis_time,
  max_entry     = sim_config$max_entry,
  seed          = 42,
  interval      = c(-3, 3),                     # default — uniroot search interval
  tol           = 1e-4,                         # default — uniroot tolerance
  n_eval        = 2000,
  verbose       = TRUE
)

t_cal <- as.numeric(difftime(Sys.time(), .t_cal_start, units = "secs"))
cal_res <- cal_uniform

# ---- calibrate cens_adjust, fixed-X design (verbatim, fixed doc) --------
sim_config <- cfg_fix
.t_cal_start <- Sys.time()

# Calibrate under the same fixed-baseline design as the main loop.
# calibrate_cens_adjust() forwards `...` to simulate_from_dgm(), so
# baseline = "fixed" applies inside the objective — but that mode requires
# n == nrow(df_source), hence n and n_eval are pinned to the panel size
# (they are NOT free Monte Carlo sizes here, unlike the random-X vignette).
n_panel <- nrow(dgm_uniform$df_source)

cal_uniform <- calibrate_cens_adjust(
  dgm           = dgm_uniform,
  target        = "rate",                       # "rate" | "km_median"
  n             = n_panel,                      # pinned: fixed panel size
  rand_ratio    = 1,                            # default
  analysis_time = sim_config$analysis_time,
  max_entry     = sim_config$max_entry,
  seed          = 42,
  interval      = c(-3, 3),                     # default — uniroot search interval
  tol           = 1e-4,                         # default — uniroot tolerance
  n_eval        = n_panel,                      # pinned: must equal n under "fixed"
  verbose       = TRUE,
  baseline      = sim_config$baseline           # forwarded to simulate_from_dgm()
)

t_cal <- as.numeric(difftime(Sys.time(), .t_cal_start, units = "secs"))
cal_fix <- cal_uniform

cat(sprintf("\ncens_adjust: resample = %.6f | fixed = %.6f\n\n",
            cal_res$cens_adjust, cal_fix$cens_adjust))

# ---- define-uniform-subgroups (verbatim) --------------------------------
set.seed(8316951)
gbsg$flag_itt <- 1L

# ── Scalar cut-points for continuous variables ────────────────────────────
# Computed once from the observed data; re-exposed as columns inside each
# simulated df_s so that subset expressions like "age <= age_med" resolve.
age_med  <- median(gbsg$age)          # ~53 yrs
size_med <- median(gbsg$size)         # ~25 mm
nodes_q3 <- quantile(gbsg$nodes, 0.75)  # top-quartile node count (~5)
pgr_med  <- median(gbsg$pgr)          # median PGR
er_cut   <- 20                         # ER-low clinical threshold (fmol / mg)

subgroups <- list(

  # ── A. Full trial (ITT) ─────────────────────────────────────────────────
  list(id = "flag_itt == 1",                       name = "All Patients",              grp = "ITT"),

  # ── B. Single clinical / factor subgroups ───────────────────────────────
  list(id = "meno == 1",                           name = "Post-menopausal",           grp = "Clinical"),
  list(id = "meno == 0",                           name = "Pre-menopausal",            grp = "Clinical"),
  list(id = "grade == 3",                          name = "Grade 3",                   grp = "Clinical"),
  list(id = "grade != 3",                          name = "Grade 1/2",                 grp = "Clinical"),

  # ── C. Single continuous-variable cuts ──────────────────────────────────
  list(id = "age <= age_med",                      name = "Age (young)",               grp = "Continuous"),
  list(id = "age > age_med",                       name = "Age (older)",               grp = "Continuous"),
  list(id = "age <= 50",                           name = "Age <= 50",                 grp = "Continuous"),
  list(id = "age > 50",                            name = "Age > 50",                  grp = "Continuous"),
  list(id = "size <= size_med",                    name = "Tumour size (small)",       grp = "Continuous"),
  list(id = "size > size_med",                     name = "Tumour size (large)",       grp = "Continuous"),
  list(id = "nodes == 0",                          name = "Node-negative",             grp = "Continuous"),
  list(id = "nodes > 0 & nodes <= 3",              name = "Nodes 1-3",                grp = "Continuous"),
  list(id = "nodes > nodes_q3",                    name = "High nodes (>Q3)",          grp = "Continuous"),
  list(id = "pgr <= pgr_med",                      name = "PGR low",                  grp = "Continuous"),
  list(id = "pgr > pgr_med",                       name = "PGR high",                 grp = "Continuous"),
  list(id = "er <= er_cut",                        name = "ER-low (<20)",              grp = "Continuous"),
  list(id = "er > er_cut",                         name = "ER-high (>=20)",            grp = "Continuous"),

  # ── D. Menopausal status x Grade (2x2) ──────────────────────────────────
  list(id = "meno == 0 & grade == 3",              name = "Pre-meno / Grade 3",        grp = "Interaction: meno x grade"),
  list(id = "meno == 1 & grade == 3",              name = "Post-meno / Grade 3",       grp = "Interaction: meno x grade"),
  list(id = "meno == 0 & grade != 3",              name = "Pre-meno / Grade 1-2",      grp = "Interaction: meno x grade"),
  list(id = "meno == 1 & grade != 3",              name = "Post-meno / Grade 1-2",     grp = "Interaction: meno x grade"),

  # ── E. Menopausal status x Age ───────────────────────────────────────────
  list(id = "meno == 0 & age <= 50",               name = "Pre-meno / Age<=50",        grp = "Interaction: meno x age"),
  list(id = "meno == 0 & age > 50",                name = "Pre-meno / Age>50",         grp = "Interaction: meno x age"),
  list(id = "meno == 1 & age <= age_med",          name = "Post-meno / Age (young)",   grp = "Interaction: meno x age"),
  list(id = "meno == 1 & age > age_med",           name = "Post-meno / Age (older)",   grp = "Interaction: meno x age"),

  # ── F. Menopausal status x ER status ────────────────────────────────────
  list(id = "meno == 0 & er <= er_cut",            name = "Pre-meno / ER-low",         grp = "Interaction: meno x ER"),
  list(id = "meno == 0 & er > er_cut",             name = "Pre-meno / ER-high",        grp = "Interaction: meno x ER"),
  list(id = "meno == 1 & er <= er_cut",            name = "Post-meno / ER-low",        grp = "Interaction: meno x ER"),
  list(id = "meno == 1 & er > er_cut",             name = "Post-meno / ER-high",       grp = "Interaction: meno x ER"),

  # ── G. Grade x Node status ───────────────────────────────────────────────
  list(id = "grade == 3 & nodes == 0",             name = "Grade 3 / Node-neg",        grp = "Interaction: grade x nodes"),
  list(id = "grade == 3 & nodes > 0",              name = "Grade 3 / Node-pos",        grp = "Interaction: grade x nodes"),
  list(id = "grade != 3 & nodes == 0",             name = "Grade 1-2 / Node-neg",      grp = "Interaction: grade x nodes"),
  list(id = "grade != 3 & nodes > 0",              name = "Grade 1-2 / Node-pos",      grp = "Interaction: grade x nodes"),
  list(id = "grade == 3 & nodes > nodes_q3",       name = "Grade 3 / High nodes",      grp = "Interaction: grade x nodes"),

  # ── H. Grade x ER status ─────────────────────────────────────────────────
  list(id = "grade == 3 & er <= er_cut",           name = "Grade 3 / ER-low",          grp = "Interaction: grade x ER"),
  list(id = "grade == 3 & er > er_cut",            name = "Grade 3 / ER-high",         grp = "Interaction: grade x ER"),
  list(id = "grade != 3 & er <= er_cut",           name = "Grade 1-2 / ER-low",        grp = "Interaction: grade x ER"),
  list(id = "grade != 3 & er > er_cut",            name = "Grade 1-2 / ER-high",       grp = "Interaction: grade x ER"),

  # ── I. Grade x PGR ──────────────────────────────────────────────────────
  list(id = "grade == 3 & pgr <= pgr_med",         name = "Grade 3 / PGR low",         grp = "Interaction: grade x PGR"),
  list(id = "grade == 3 & pgr > pgr_med",          name = "Grade 3 / PGR high",        grp = "Interaction: grade x PGR"),
  list(id = "grade != 3 & pgr <= pgr_med",         name = "Grade 1-2 / PGR low",       grp = "Interaction: grade x PGR"),
  list(id = "grade != 3 & pgr > pgr_med",          name = "Grade 1-2 / PGR high",      grp = "Interaction: grade x PGR"),

  # ── J. Age x ER x Menopausal (3-way) ───────────────────────────────────
  list(id = "meno == 0 & age <= 50 & er <= er_cut",  name = "Pre-meno/Yng/ER-low",    grp = "3-way"),
  list(id = "meno == 0 & age <= 50 & er > er_cut",   name = "Pre-meno/Yng/ER-high",   grp = "3-way"),
  list(id = "meno == 1 & grade == 3 & er <= er_cut", name = "Post-meno/G3/ER-low",    grp = "3-way"),
  list(id = "meno == 1 & grade == 3 & nodes > 0",    name = "Post-meno/G3/Node-pos",  grp = "3-way"),
  list(id = "grade == 3 & nodes > 0 & er <= er_cut", name = "G3/Node-pos/ER-low",     grp = "3-way"),

  # ── K. Size x Node combinations ─────────────────────────────────────────
  list(id = "size > size_med & nodes > 0",         name = "Large/Node-pos",            grp = "Interaction: size x nodes"),
  list(id = "size > size_med & nodes == 0",        name = "Large/Node-neg",            grp = "Interaction: size x nodes"),
  list(id = "size <= size_med & nodes > 0",        name = "Small/Node-pos",            grp = "Interaction: size x nodes"),
  list(id = "size <= size_med & nodes == 0",       name = "Small/Node-neg",            grp = "Interaction: size x nodes"),

  # ── L. Random benchmark subgroups ───────────────────────────────────────
  # No clinical meaning; regenerated each simulation. Calibrate variability
  # by size against the clinical and interaction subgroups above.
  list(id = "random60 == 1",                       name = "random60",                  grp = "Random (N~60)"),
  list(id = "random40 == 1",                       name = "random40",                  grp = "Random (N~40)"),
  list(id = "random20 == 1",                       name = "random20",                  grp = "Random (N~20)"),
  list(id = "random15 == 1",                       name = "random15",                  grp = "Random (N~15)")
)

# ── Show observed sizes from GBSG ────────────────────────────────────────
# Scalar cut-points exposed as columns so eval() resolves them correctly.
gbsg_eval           <- gbsg
gbsg_eval$age_med   <- age_med
gbsg_eval$size_med  <- size_med
gbsg_eval$nodes_q3  <- nodes_q3
gbsg_eval$pgr_med   <- pgr_med
gbsg_eval$er_cut    <- er_cut
gbsg_eval$random60  <- 0L
gbsg_eval$random40  <- 0L
gbsg_eval$random20  <- 0L
gbsg_eval$random15  <- 0L

sg_sizes <- sapply(subgroups, function(s) {
  tryCatch(
    sum(eval(parse(text = s$id), envir = gbsg_eval)),
    error = function(e) NA_integer_
  )
})
is_random_sg <- grepl("^random", sapply(subgroups, `[[`, "name"))
rand_sizes   <- c(60L, 40L, 20L, 15L)

sg_overview <- data.frame(
  Subgroup = sapply(subgroups, `[[`, "name"),
  Group    = sapply(subgroups, `[[`, "grp"),
  N_source = ifelse(is_random_sg, rand_sizes[cumsum(is_random_sg)], sg_sizes),
  stringsAsFactors = FALSE
)
print(sg_overview, row.names = FALSE)
cat("
Total subgroups:", nrow(sg_overview), "
")

# ---- cut_points bundle (verbatim lines from uniform-sims) ---------------
cut_points <- list(age_med = age_med, size_med = size_med,
                    nodes_q3 = nodes_q3, pgr_med = pgr_med, er_cut = er_cut)

# ---- transcribed per-trial helpers (verbatim; renamed per design) -------
cox_sR <- function(data) {
  fit <- tryCatch(
    coxph(Surv(y_sim, event_sim) ~ treat_sim + strata(grade), data = data),
    error = function(e) NULL
  )
  if (is.null(fit) || nrow(summary(fit)$conf.int) == 0) {
    return(c(NA_real_, NA_real_))
  }
  ci <- summary(fit)$conf.int[1, ]
  c(ci["exp(coef)"], ci["upper .95"])
}

# Run one simulated trial and return per-subgroup (HR, UB, N) row vectors.
# `cut_points`   list of scalars exposed as columns of df_s so subgroup
#                expressions like "age <= age_med" resolve correctly.
# `cox_fn`       per-subgroup analysis function — swap for sensitivity
#                analyses (unstratified, X+sR, etc.) without touching
#                the loop.
run_one_sim <- function(ss, dgm, cens_adjust, sim_config,
                         subgroups, cut_points, cox_fn) {

  # 1. Draw one synthetic trial from the calibrated DGM.
  #    baseline = "fixed": the covariate panel is dgm$df_source (every
  #    source patient exactly once), identical for every ss — only
  #    treatment, entry, event times, and censoring are re-drawn.  n is
  #    implied by the panel and deliberately not passed.
  df_s <- simulate_from_dgm(
    dgm           = dgm,
    baseline      = sim_config$baseline,
    analysis_time = sim_config$analysis_time,
    max_entry     = sim_config$max_entry,
    cens_adjust   = cens_adjust,
    seed          = sim_config$seed_base + ss
  )

  # 2. Expose flag_itt and scalar cut-points as columns
  df_s$flag_itt <- 1L
  for (nm in names(cut_points)) df_s[[nm]] <- cut_points[[nm]]

  # 3. Regenerate random-benchmark subgroups for this trial.
  #    Membership is re-drawn each replicate, but over the FIXED panel:
  #    sizes are exactly 60/40/20/15 every trial (n_s is constant), so the
  #    random benchmarks isolate pure membership noise at fixed N.
  set.seed(sim_config$seed_base + ss + sim_config$rand_seed_offset)
  n_s   <- nrow(df_s)
  r_idx <- sample.int(n_s, min(60L, n_s), replace = FALSE)
  df_s$random60 <- as.integer(seq_len(n_s) %in% r_idx)
  df_s$random40 <- as.integer(seq_len(n_s) %in% r_idx[seq_len(min(40L, length(r_idx)))])
  df_s$random20 <- as.integer(seq_len(n_s) %in% r_idx[seq_len(min(20L, length(r_idx)))])
  df_s$random15 <- as.integer(seq_len(n_s) %in% r_idx[seq_len(min(15L, length(r_idx)))])

  # 4. Fit cox_fn() in every pre-defined subgroup
  sg_ids <- vapply(subgroups, `[[`, character(1L), "id")
  n_sg   <- length(subgroups)
  hr_row <- rep(NA_real_, n_sg)
  ub_row <- rep(NA_real_, n_sg)
  n_row  <- rep(NA_real_, n_sg)

  for (gg in seq_len(n_sg)) {
    df_sg <- tryCatch(
      subset(df_s, eval(parse(text = sg_ids[[gg]]))),
      error = function(e) NULL
    )
    if (is.null(df_sg) || nrow(df_sg) < 5) next
    n_row[gg]  <- nrow(df_sg)
    r          <- cox_fn(df_sg)
    hr_row[gg] <- r[1]
    ub_row[gg] <- r[2]
  }

  list(hr = hr_row, ub = ub_row, n = n_row)
}

run_one_sim_res <- function(ss, dgm, cens_adjust, sim_config,
                         subgroups, cut_points, cox_fn) {

  # 1. Draw one synthetic trial from the calibrated DGM
  df_s <- simulate_from_dgm(
    dgm           = dgm,
    n             = sim_config$n_per_trial,
    analysis_time = sim_config$analysis_time,
    max_entry     = sim_config$max_entry,
    cens_adjust   = cens_adjust,
    seed          = sim_config$seed_base + ss
  )

  # 2. Expose flag_itt and scalar cut-points as columns
  df_s$flag_itt <- 1L
  for (nm in names(cut_points)) df_s[[nm]] <- cut_points[[nm]]

  # 3. Regenerate random-benchmark subgroups for this trial
  set.seed(sim_config$seed_base + ss + sim_config$rand_seed_offset)
  n_s   <- nrow(df_s)
  r_idx <- sample.int(n_s, min(60L, n_s), replace = FALSE)
  df_s$random60 <- as.integer(seq_len(n_s) %in% r_idx)
  df_s$random40 <- as.integer(seq_len(n_s) %in% r_idx[seq_len(min(40L, length(r_idx)))])
  df_s$random20 <- as.integer(seq_len(n_s) %in% r_idx[seq_len(min(20L, length(r_idx)))])
  df_s$random15 <- as.integer(seq_len(n_s) %in% r_idx[seq_len(min(15L, length(r_idx)))])

  # 4. Fit cox_fn() in every pre-defined subgroup
  sg_ids <- vapply(subgroups, `[[`, character(1L), "id")
  n_sg   <- length(subgroups)
  hr_row <- rep(NA_real_, n_sg)
  ub_row <- rep(NA_real_, n_sg)
  n_row  <- rep(NA_real_, n_sg)

  for (gg in seq_len(n_sg)) {
    df_sg <- tryCatch(
      subset(df_s, eval(parse(text = sg_ids[[gg]]))),
      error = function(e) NULL
    )
    if (is.null(df_sg) || nrow(df_sg) < 5) next
    n_row[gg]  <- nrow(df_sg)
    r          <- cox_fn(df_sg)
    hr_row[gg] <- r[1]
    ub_row[gg] <- r[2]
  }

  list(hr = hr_row, ub = ub_row, n = n_row)
}

run_one_sim_fix <- function(ss, dgm, cens_adjust, sim_config,
                         subgroups, cut_points, cox_fn) {

  # 1. Draw one synthetic trial from the calibrated DGM.
  #    baseline = "fixed": the covariate panel is dgm$df_source (every
  #    source patient exactly once), identical for every ss — only
  #    treatment, entry, event times, and censoring are re-drawn.  n is
  #    implied by the panel and deliberately not passed.
  df_s <- simulate_from_dgm(
    dgm           = dgm,
    baseline      = sim_config$baseline,
    analysis_time = sim_config$analysis_time,
    max_entry     = sim_config$max_entry,
    cens_adjust   = cens_adjust,
    seed          = sim_config$seed_base + ss
  )

  # 2. Expose flag_itt and scalar cut-points as columns
  df_s$flag_itt <- 1L
  for (nm in names(cut_points)) df_s[[nm]] <- cut_points[[nm]]

  # 3. Regenerate random-benchmark subgroups for this trial.
  #    Membership is re-drawn each replicate, but over the FIXED panel:
  #    sizes are exactly 60/40/20/15 every trial (n_s is constant), so the
  #    random benchmarks isolate pure membership noise at fixed N.
  set.seed(sim_config$seed_base + ss + sim_config$rand_seed_offset)
  n_s   <- nrow(df_s)
  r_idx <- sample.int(n_s, min(60L, n_s), replace = FALSE)
  df_s$random60 <- as.integer(seq_len(n_s) %in% r_idx)
  df_s$random40 <- as.integer(seq_len(n_s) %in% r_idx[seq_len(min(40L, length(r_idx)))])
  df_s$random20 <- as.integer(seq_len(n_s) %in% r_idx[seq_len(min(20L, length(r_idx)))])
  df_s$random15 <- as.integer(seq_len(n_s) %in% r_idx[seq_len(min(15L, length(r_idx)))])

  # 4. Fit cox_fn() in every pre-defined subgroup
  sg_ids <- vapply(subgroups, `[[`, character(1L), "id")
  n_sg   <- length(subgroups)
  hr_row <- rep(NA_real_, n_sg)
  ub_row <- rep(NA_real_, n_sg)
  n_row  <- rep(NA_real_, n_sg)

  for (gg in seq_len(n_sg)) {
    df_sg <- tryCatch(
      subset(df_s, eval(parse(text = sg_ids[[gg]]))),
      error = function(e) NULL
    )
    if (is.null(df_sg) || nrow(df_sg) < 5) next
    n_row[gg]  <- nrow(df_sg)
    r          <- cox_fn(df_sg)
    hr_row[gg] <- r[1]
    ub_row[gg] <- r[2]
  }

  list(hr = hr_row, ub = ub_row, n = n_row)
}

# ---- transcribed sequential loops (recombine as in uniform-sims) --------
loop_mats <- function(one_fun, cfg, cens_adjust) {
  res <- lapply(seq_len(n_acc), one_fun, dgm = dgm_uniform,
                cens_adjust = cens_adjust, sim_config = cfg,
                subgroups = subgroups, cut_points = cut_points,
                cox_fn = cox_sR)
  sg <- vapply(subgroups, `[[`, character(1L), "name")
  g  <- function(f) { x <- do.call(rbind, lapply(res, `[[`, f))
                      colnames(x) <- sg; x }
  list(hr = g("hr"), ub = g("ub"), n = g("n"))
}

t0 <- Sys.time()
V_res <- loop_mats(run_one_sim_res, cfg_res, cal_res$cens_adjust)
t_loop_res <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
t0 <- Sys.time()
V_fix <- loop_mats(run_one_sim_fix, cfg_fix, cal_fix$cens_adjust)
t_loop_fix <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat(sprintf("Transcribed loops: resample %.1f s | fixed %.1f s\n",
            t_loop_res, t_loop_fix))

# ---- package runner (default parallel workers) --------------------------
fit_v <- subgroup_cox(Surv(y_sim, event_sim) ~ treat_sim + strata(grade))

t0 <- Sys.time()
R_res <- run_subgroup_sims(
  dgm = dgm_uniform, subgroups = subgroups, n_sims = n_acc, fit = fit_v,
  baseline = "resample", n = nrow(gbsg),
  analysis_time = 84, max_entry = 24, cens_adjust = cal_res$cens_adjust,
  cutpoints = cut_points, hr_true = 0.70, k_treat = k_treat_uniform)
R_fix <- run_subgroup_sims(
  dgm = dgm_uniform, subgroups = subgroups, n_sims = n_acc, fit = fit_v,
  baseline = "fixed",
  analysis_time = 84, max_entry = 24, cens_adjust = cal_fix$cens_adjust,
  cutpoints = cut_points, hr_true = 0.70, k_treat = k_treat_uniform)
R_fix_seq <- run_subgroup_sims(
  dgm = dgm_uniform, subgroups = subgroups, n_sims = n_acc, fit = fit_v,
  baseline = "fixed", workers = 1,
  analysis_time = 84, max_entry = 24, cens_adjust = cal_fix$cens_adjust,
  cutpoints = cut_points, hr_true = 0.70, k_treat = k_treat_uniform)
t_runner <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat(sprintf("Runner (2 parallel studies + 1 sequential): %.1f s\n\n",
            t_runner))

# ---- GATES: bit-identity ------------------------------------------------
ok(identical(V_res$hr, R_res$sim_hrs), "resample: sim_hrs bit-identical")
ok(identical(V_res$ub, R_res$sim_ubs), "resample: sim_ubs bit-identical")
ok(identical(V_res$n,  R_res$sim_ns),  "resample: sim_ns  bit-identical")
ok(identical(V_fix$hr, R_fix$sim_hrs), "fixed   : sim_hrs bit-identical")
ok(identical(V_fix$ub, R_fix$sim_ubs), "fixed   : sim_ubs bit-identical")
ok(identical(V_fix$n,  R_fix$sim_ns),  "fixed   : sim_ns  bit-identical")
ok(identical(R_fix_seq$sim_hrs, R_fix$sim_hrs) &&
     identical(R_fix_seq$sim_ubs, R_fix$sim_ubs) &&
     identical(R_fix_seq$sim_ns,  R_fix$sim_ns),
   "fixed   : parallel == sequential (real pipeline)")

# ---- GATE: summary() vs transcribed results table (fixed design) --------
n_sims_uniform <- n_acc
sim_hrs <- V_fix$hr; sim_ubs <- V_fix$ub; sim_ns <- V_fix$n
sg_names <- colnames(sim_hrs)

# (verbatim uniform-results-table chunk, final print() dropped)
# Median HR and 1-99th percentile ECI across simulations
hr_q  <- t(apply(sim_hrs, 2, quantile, probs = c(0.01, 0.50, 0.99),
                 na.rm = TRUE))
ub_q  <- t(apply(sim_ubs, 2, quantile, probs = c(0.01, 0.50, 0.99),
                 na.rm = TRUE))

# n_valid: number of simulations where Cox estimation succeeded
n_valid     <- apply(sim_hrs, 2, function(x) sum(!is.na(x)))
pct_valid   <- round(100 * n_valid / n_sims_uniform, 1)

# Under baseline = "fixed" subgroup sizes are constant across trials, so
# this "mean" is the exact per-trial N for every covariate-defined subgroup
# (and exactly 60/40/20/15 for the random benchmarks, whose membership —
# not size — is re-drawn each replicate).
mean_n      <- round(apply(sim_ns,  2, mean, na.rm = TRUE))

# Tail probabilities and unconditional medians used across the results
# table and every forest panel.  Computed here once so every downstream
# chunk references the same vector — denominated by n_valid (not
# n_sims_uniform) so NAs from estimation failures are excluded.
pr_hr_lt050 <- apply(sim_hrs, 2, function(x) mean(x < 0.5, na.rm = TRUE))
pr_hr_gt1   <- apply(sim_hrs, 2, function(x) mean(x > 1.0, na.rm = TRUE))
pr_ub_ge2   <- apply(sim_ubs, 2, function(x) mean(x >= 2.0, na.rm = TRUE))
pr_ub_ge3   <- apply(sim_ubs, 2, function(x) mean(x >= 3.0, na.rm = TRUE))

# Unconditional median displays (formatted as character for annotation columns)
mhr_uncond <- ifelse(is.na(hr_q[, 2]), "-",
                     as.character(round(hr_q[, 2], 2)))
mub_uncond <- ifelse(is.na(ub_q[, 2]), "-",
                     as.character(round(ub_q[, 2], 2)))

fmt_pct <- function(x) ifelse(is.nan(x), "-", paste0(round(100 * x, 1), "%"))

results_tbl <- data.frame(
  Subgroup    = sg_names,
  Mean_N      = ifelse(is.nan(mean_n), "-", as.character(mean_n)),
  N_valid     = paste0(n_valid, " (", pct_valid, "%)"),
  # ── HR-panel columns ────────────────────────────────────────────────────
  Pr_HR_lt050 = fmt_pct(pr_hr_lt050),
  Pr_HR_gt1   = fmt_pct(pr_hr_gt1),
  mHR         = mhr_uncond,
  # ── UB-panel columns ────────────────────────────────────────────────────
  Pr_UB_ge2   = fmt_pct(pr_ub_ge2),
  Pr_UB_ge3   = fmt_pct(pr_ub_ge3),
  mUB         = mub_uncond,
  stringsAsFactors = FALSE
)
colnames(results_tbl) <- c(
  "Subgroup", "N", "#(% converged)",
  "Pr(HR<0.5)", "Pr(HR>1.0)", "mHR",
  "Pr(UB>=2)", "Pr(UB>=3)", "mUB"
)

S <- summary(R_fix, hr_true = 0.70)
ok(identical(names(results_tbl), names(S$results_tbl)),
   "summary: results table column names identical")
for (j in seq_along(results_tbl)) {
  ok(identical(results_tbl[[j]], S$results_tbl[[j]]),
     sprintf("summary: column '%s' identical", names(results_tbl)[j]))
}
ok(S$n_single + S$n_combo == sum(S$ok),
   "summary: single/combo panels partition the ok subgroups")
hi <- S$sg_names[S$highrisk$include][S$highrisk$ord]
ok(length(hi) >= 1 && hi[1] == "All Patients",
   "summary: high-risk panel anchors ITT first")

# ---- verdict ------------------------------------------------------------
cat(sprintf("\nTotal acceptance runtime: %.1f min\n",
            as.numeric(difftime(Sys.time(), t_all, units = "mins"))))
if (.fails == 0L) {
  cat("\n==== PHASE 1 ACCEPTANCE: ALL GATES PASSED ====\n")
} else {
  cat(sprintf("\n==== PHASE 1 ACCEPTANCE: %d GATE(S) FAILED ====\n", .fails))
  quit(status = 1L)
}
