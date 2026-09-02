#!/usr/bin/env Rscript
# =========================================================================
# PHASE 2 TIMING DIAGNOSTIC (read-only; writes nothing, commits nothing)
#
# Question 1 -- what actually happened in the ported renders?
#   Reads workers / plan_reused / t_plan_secs / t_sims_secs from the two
#   *_ported.rds payloads the Gate 3 renders wrote.
#
# Question 2 -- if the pool engaged, where is the loop overhead?
#   Rebuilds the fixed-design DGM (the slow part, ~4-6 min), then times
#   run_subgroup_sims() at a small n under an outer multisession plan
#   with three fits that differ ONLY in what serializes to workers:
#     A. fit whose formula environment is a HEAVY non-global env
#        (reproduces the suspected knitr/quarto chunk-env condition)
#     B. fit built at top level (formula env = globalenv, by-reference)
#     C. fit with a lean baseenv closure (nothing to drag)
#   plus a sequential baseline.  Identical results across A/B/C are
#   asserted; the three wall times localize the cost.
#
# Usage:  Rscript dev/diag_phase2_timing.R [n_diag]      # default 400
# =========================================================================
suppressMessages({ library(forestsearch); library(survival) })
args   <- commandArgs(trailingOnly = TRUE)
n_diag <- if (length(args) >= 1) as.integer(args[1]) else 400L

res_dir <- "quarto/extreme_subgroups/fixed_random/results"

cat("==== Q1: sim_config recorded by the Gate 3 renders ====\n")
for (d in c("resample", "fixed")) {
  f <- file.path(res_dir, sprintf("extreme_sims_%s_10000_ported.rds", d))
  if (!file.exists(f)) { cat("  MISSING:", f, "\n"); next }
  sc <- readRDS(f)$sim_config
  cat(sprintf(
    "  %-8s workers = %s | plan_reused = %s | t_plan = %.2f s | t_sims = %.1f s\n",
    d, sc$workers, sc$plan_reused, sc$t_plan_secs, sc$t_sims_secs))
}
cat("\nInterpretation: plan_reused = FALSE or workers = 1  -> reuse never\n")
cat("engaged in-render (context problem).  plan_reused = TRUE with large\n")
cat("workers and large t_sims -> loop-internal overhead; Q2 localizes it.\n\n")

# ---- rebuild fixed-design DGM (verbatim vignette chunks) ----------------
cat("==== Rebuilding DGM + fixed-design calibration (~4-6 min) ====\n")
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

sim_config <- list(
  baseline         = "fixed",
  analysis_time    = 84,
  max_entry        = 24,
  seed_base        = 0L,
  rand_seed_offset = 1e6L
)

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

cut_points <- list(age_med = age_med, size_med = size_med,
                    nodes_q3 = nodes_q3, pgr_med = pgr_med, er_cut = er_cut)

# ---- three fits differing only in serialization footprint ---------------
# A: formula created inside a heavy non-global environment (as a knitr
#    chunk env would be), so serializing the fit drags that env by value.
heavy <- new.env(parent = globalenv())
assign("dgm_uniform_copy", dgm_uniform, envir = heavy)
assign("gbsg_copy",        gbsg,        envir = heavy)
assign("base_args_copy",   base_args_uniform, envir = heavy)
fml_heavy <- eval(quote(Surv(y_sim, event_sim) ~ treat_sim + strata(grade)),
                  envir = heavy)
stopifnot(identical(environment(fml_heavy), heavy))
fit_A <- subgroup_cox(fml_heavy)

# B: formula at top level of this script (env = globalenv, by reference).
fit_B <- subgroup_cox(Surv(y_sim, event_sim) ~ treat_sim + strata(grade))

# C: lean closure, no formula, environment = baseenv().
fit_C <- local(
  function(data) {
    fitc <- tryCatch(
      survival::coxph(
        survival::Surv(y_sim, event_sim) ~ treat_sim +
          survival::strata(grade), data = data),
      error = function(e) NULL)
    if (is.null(fitc) || nrow(summary(fitc)$conf.int) == 0) {
      return(c(NA_real_, NA_real_))
    }
    ci <- summary(fitc)$conf.int[1, ]
    c(ci["exp(coef)"], ci["upper .95"])
  }, envir = baseenv())

cat(sprintf("\nserialized fit sizes:  A (heavy env) = %.4f MB | B = %.4f MB | C = %.4f MB\n",
  length(serialize(fit_A, NULL)) / 1024^2,
  length(serialize(fit_B, NULL)) / 1024^2,
  length(serialize(fit_C, NULL)) / 1024^2))

# ---- timed runs ----------------------------------------------------------
future::plan(future::multisession,
             workers = max(1L, as.integer(ceiling(
               0.90 * future::availableCores()))))
# multisession spawns its cluster lazily on FIRST USE, not at plan();
# without this warm-up the first timed variant absorbs the full
# worker-pool startup (~145 s observed at 116 workers) and the
# per-variant comparison is confounded.  The warm-up reports spawn as
# its own measured line -- the ported renders pay this once inside
# t_sims, which is what Gate 4's wall-time budget covers.
cat("warming worker pool (spawns the full cluster) ...\n")
t0 <- Sys.time()
invisible(future::value(future::future(Sys.getpid())))
cat(sprintf("pool warm-up (%d workers): %.1f s\n",
            future::nbrOfWorkers(),
            as.numeric(difftime(Sys.time(), t0, units = "secs"))))
# A trivial future does NOT trigger the ~140 s one-time cost the first
# real runner-shaped workload pays on a fresh pool (probed: not spawn --
# 0.4 s for 116 PIDs -- and not loadNamespace fan-out -- 4 s; consistent
# with first-use lazy-load forcing plus first-touch memory across
# workers).  Warm with a tiny real run so A/B/C time the fit mechanism,
# and report the initialization as its own measured line.  The ported
# renders pay this once inside t_sims; Gate 4's wall budget covers it.
cat("first-workload warm-up (run_subgroup_sims) ...\n")
# The initialization is PER-WORKER (~1.27 s/worker, dispatch-serialized;
# rate confirmed across three independent runs) and triggers only on a
# real runner-shaped shipment.  n_sims = 2 x nbrOfWorkers() covers every
# worker even under scheduling variation -- n_sims = 8 warmed exactly 8
# of 116 and left A holding the other 108.  RNG-misuse notes formerly
# suppressed locally here are now silenced package-side: the runner
# scopes future.rng.onMisuse around its loop, so this diagnostic's full
# log doubles as the gate (zero UNRELIABLE lines expected).
n_warm <- max(8L, 2L * as.integer(future::nbrOfWorkers()))
t0 <- Sys.time()
invisible(run_subgroup_sims(
  dgm = dgm_uniform, subgroups = subgroups, n_sims = n_warm, fit = fit_C,
  baseline = "fixed", analysis_time = 84, max_entry = 24,
  cens_adjust = cal_uniform$cens_adjust, cutpoints = cut_points,
  workers = NULL, hr_true = 0.70))
cat(sprintf(
  "first-workload warm-up (n_sims = %d, covers every worker): %.1f s  (one-time per fresh pool)\n",
  n_warm, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
run1 <- function(fit, label, workers = NULL) {
  t0 <- Sys.time()
  r <- run_subgroup_sims(
    dgm = dgm_uniform, subgroups = subgroups, n_sims = n_diag, fit = fit,
    baseline = "fixed", analysis_time = 84, max_entry = 24,
    cens_adjust = cal_uniform$cens_adjust, cutpoints = cut_points,
    workers = workers, hr_true = 0.70, verbose = TRUE)
  wt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cat(sprintf("  %-28s wall = %6.1f s | t_plan = %5.2f | t_sims = %6.1f | reused = %s | workers = %d\n",
    label, wt, r$sim_config$t_plan_secs, r$sim_config$t_sims_secs,
    r$sim_config$plan_reused, r$sim_config$workers))
  r
}
cat(sprintf("\n==== Q2: timed runs, n_diag = %d (outer plan: multisession) ====\n", n_diag))
R_A <- run1(fit_A, "A: heavy-env formula")
R_B <- run1(fit_B, "B: globalenv formula")
R_C <- run1(fit_C, "C: lean baseenv closure")
R_S <- run1(fit_C, "S: sequential baseline", workers = 1)
future::plan(future::sequential)

stopifnot(identical(R_A$sim_hrs, R_B$sim_hrs),
          identical(R_B$sim_hrs, R_C$sim_hrs),
          identical(R_C$sim_hrs, R_S$sim_hrs))
cat("\n[PASS] A/B/C/S matrices identical -- fits are numerically equivalent\n")
cat("\nReading the three parallel times:\n")
cat("  A >> B ~ C  -> formula-environment export is the render cost;\n")
cat("                 fix: environment hygiene in subgroup_cox().\n")
cat("  A ~ B ~ C, all >> compute floor -> per-future overhead unrelated\n")
cat("                 to fit; next probe: chunking / globals of dgm.\n")
cat("  all fast     -> in-script fine; problem is quarto-context only\n")
cat("                 (see the probe qmd in the kickoff).\n")
