# accept_phase41_survival_bitident.R --------------------------------------
# SENTINEL: P41-BITIDENT-v1-2026-09-02
#
# Phase 4.1 survival byte-identity gate. Two modes across the install
# boundary:
#   Rscript dev/accept_phase41_survival_bitident.R capture
#     -- run with the OLD installed forestsearch (pre-Phase-4.1 build);
#        asserts subgroup_glm is NOT exported, runs the reference study
#        for both baselines, saves /tmp/p41_baseline.rds.
#   Rscript dev/accept_phase41_survival_bitident.R compare
#     -- run after devtools::document() + install of the Phase-4.1 code;
#        asserts subgroup_glm IS exported, re-runs the identical study,
#        and requires identical() on the FULL result objects after
#        nulling only the volatile fields (created, t_plan_secs,
#        t_sims_secs). Field names must match exactly and no `effect`
#        field may appear on survival results.
#
# Uses library(forestsearch) ONLY (never load_all: the point is the
# installed builds). cens_adjust is a fixed constant -- this is an
# identity gate, not a calibration.

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args) >= 1L) args[[1L]] else ""
if (!mode %in% c("capture", "compare")) {
  stop("usage: Rscript dev/accept_phase41_survival_bitident.R capture|compare")
}
suppressPackageStartupMessages({ library(forestsearch); library(survival) })
baseline_rds <- "/tmp/p41_baseline.rds"

has_glm <- "subgroup_glm" %in% getNamespaceExports("forestsearch")
cat(sprintf("installed forestsearch %s | subgroup_glm exported: %s\n",
            as.character(utils::packageVersion("forestsearch")), has_glm))
if (mode == "capture" && has_glm) {
  stop("capture mode requires the PRE-Phase-4.1 installed build ",
       "(subgroup_glm must not be exported yet).")
}
if (mode == "compare" && !has_glm) {
  stop("compare mode requires the Phase-4.1 build to be installed ",
       "(subgroup_glm not found in the namespace).")
}

run_study <- function() {
  data(gbsg, package = "survival")
  gbsg$time_months <- gbsg$rfstime / 30.4375
  set.seed(99)
  dgm <- generate_aft_dgm_flex(
    data = gbsg,
    continuous_vars = c("age", "size", "nodes", "pgr", "er"),
    factor_vars = c("meno", "grade"),
    outcome_var = "time_months", event_var = "status",
    treatment_var = "hormon", subgroup_vars = NULL,
    model = "null", n_super = 1500, seed = 99, verbose = FALSE)

  age_med <- stats::median(gbsg$age)
  sgs <- list(
    list(id = "flag_itt == 1",      name = "All Patients", grp = "ITT"),
    list(id = "meno == 1",          name = "Post-meno",    grp = "Clinical"),
    list(id = "meno == 0",          name = "Pre-meno",     grp = "Clinical"),
    list(id = "grade == 3",         name = "Grade 3",      grp = "Clinical"),
    list(id = "age <= age_med",     name = "Age young",    grp = "Continuous"),
    list(id = "age > age_med",      name = "Age older",    grp = "Continuous"),
    list(id = "meno==1 & grade==3", name = "Post/G3",
         grp = "Interaction: mg"),
    list(id = "random60 == 1",      name = "random60",  grp = "Random (N~60)"),
    list(id = "random40 == 1",      name = "random40",  grp = "Random (N~40)"),
    list(id = "random20 == 1",      name = "random20",  grp = "Random (N~20)"),
    list(id = "random15 == 1",      name = "random15",  grp = "Random (N~15)"))

  fit <- subgroup_cox(
    survival::Surv(y_sim, event_sim) ~ treat_sim + survival::strata(grade))

  common <- list(dgm = dgm, subgroups = sgs, n_sims = 25L, fit = fit,
                 analysis_time = 84, max_entry = 24, cens_adjust = 0.45,
                 cutpoints = list(age_med = age_med),
                 benchmarks = benchmark_spec(), min_n = 5L,
                 workers = 1, hr_true = 0.70, k_treat = 0.95,
                 validate = TRUE, verbose = FALSE)
  res <- do.call(run_subgroup_sims,
                 c(common, list(baseline = "resample", n = 200L)))
  fix <- do.call(run_subgroup_sims,
                 c(common, list(baseline = "fixed")))
  list(resample = res, fixed = fix)
}

null_volatile <- function(x) {
  x$created <- NULL
  x$sim_config$t_plan_secs <- NULL
  x$sim_config$t_sims_secs <- NULL
  x
}

t0 <- Sys.time()
cur <- run_study()
cat(sprintf("study ran in %.1f s\n",
            as.numeric(difftime(Sys.time(), t0, units = "secs"))))

if (mode == "capture") {
  saveRDS(cur, baseline_rds)
  cat("BASELINE CAPTURED ->", baseline_rds, "\n")
  cat(sprintf("resample sim_hrs[1:2, 1:3]:\n"))
  print(cur$resample$sim_hrs[1:2, 1:3])
  quit(status = 0)
}

# --- compare -------------------------------------------------------------
if (!file.exists(baseline_rds)) {
  stop("baseline missing: run the capture mode first (", baseline_rds, ")")
}
ref <- readRDS(baseline_rds)
for (bl in c("resample", "fixed")) {
  a <- null_volatile(ref[[bl]])
  b <- null_volatile(cur[[bl]])
  if (!identical(names(a), names(b))) {
    stop("FIELD-NAME DRIFT on the survival ", bl, " path:\n  old: ",
         paste(names(a), collapse = ", "), "\n  new: ",
         paste(names(b), collapse = ", "))
  }
  if (!is.null(cur[[bl]]$effect)) {
    stop("survival ", bl, " result unexpectedly carries an `effect` field.")
  }
  if (!identical(a, b)) {
    cat("MISMATCH DIAGNOSIS (", bl, "): field-by-field identical()?\n")
    for (nm in names(a)) cat(sprintf("  %-22s %s\n", nm,
                                     identical(a[[nm]], b[[nm]])))
    stop("survival ", bl, ": old vs new objects are NOT identical.")
  }
  cat("PASS:", bl, "-- full objects identical (volatile fields nulled)\n")
}
cat("ACCEPT: Phase 4.1 survival byte-identity holds on both baselines.\n")
