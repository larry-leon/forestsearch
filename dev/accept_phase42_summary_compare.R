# accept_phase42_summary_compare.R -----------------------------------------
# SENTINEL: P42-ACC-v2-2026-09-02  (r2: guard-order fixtures)
#
# Phase 4.2 acceptance gate. Two modes across the install boundary:
#   Rscript dev/accept_phase42_summary_compare.R capture
#     -- run with the Phase-4.1 installed build (summary() has no
#        est_thresholds formal yet); computes summaries and the compare
#        frame for a reference survival study on both baselines,
#        including the printed tables, and saves /tmp/p42_baseline.rds.
#   Rscript dev/accept_phase42_summary_compare.R compare
#     -- run after installing Phase 4.2; requires full identical() on
#        every captured object and every captured print, then runs the
#        GLM behavior sections (headers, hand-checked tails, "-"
#        degradation, high-risk panel, overrides, metadata guard,
#        effect attribute) against a real generate_glm_dgm() study.
#
# library(forestsearch) only (installed builds are the test subjects).

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args) >= 1L) args[[1L]] else ""
if (!mode %in% c("capture", "compare")) {
  stop("usage: Rscript dev/accept_phase42_summary_compare.R capture|compare")
}
suppressPackageStartupMessages({ library(forestsearch); library(survival) })
baseline_rds <- "/tmp/p42_baseline.rds"
fail <- function(...) stop(sprintf(...), call. = FALSE)
ok   <- function(...) cat("PASS:", sprintf(...), "\n")

has_42 <- "est_thresholds" %in%
  names(formals(utils::getS3method("summary", "subgroup_sims")))
cat(sprintf("installed forestsearch %s | effect-aware summary(): %s\n",
            as.character(utils::packageVersion("forestsearch")), has_42))
if (mode == "capture" && has_42)
  fail("capture mode requires the pre-4.2 installed build.")
if (mode == "compare" && !has_42)
  fail("compare mode requires the Phase-4.2 build to be installed.")

run_survival_study <- function() {
  gbsg <- survival::gbsg
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
    list(id = "age <= age_med",     name = "Age young",  grp = "Continuous"),
    list(id = "age > age_med",      name = "Age older",  grp = "Continuous"),
    list(id = "meno==1 & grade==3", name = "Post/G3",
         grp = "Interaction: mg"),
    list(id = "random60 == 1", name = "random60", grp = "Random (N~60)"),
    list(id = "random40 == 1", name = "random40", grp = "Random (N~40)"),
    list(id = "random20 == 1", name = "random20", grp = "Random (N~20)"),
    list(id = "random15 == 1", name = "random15", grp = "Random (N~15)"))
  fit <- subgroup_cox(
    survival::Surv(y_sim, event_sim) ~ treat_sim + survival::strata(grade))
  common <- list(dgm = dgm, subgroups = sgs, n_sims = 25L, fit = fit,
                 analysis_time = 84, max_entry = 24, cens_adjust = 0.45,
                 cutpoints = list(age_med = age_med),
                 benchmarks = benchmark_spec(), min_n = 5L,
                 workers = 1, hr_true = 0.70, k_treat = 0.95,
                 validate = TRUE, verbose = FALSE)
  res_r <- do.call(run_subgroup_sims,
                   c(common, list(baseline = "resample", n = 200L)))
  res_f <- do.call(run_subgroup_sims, c(common, list(baseline = "fixed")))
  list(res_r = res_r, res_f = res_f)
}

pack <- function(st) {
  s_r <- summary(st$res_r); s_f <- summary(st$res_f)
  list(
    s_r = s_r, s_f = s_f,
    p_r = utils::capture.output(print(s_r)),
    p_f = utils::capture.output(print(s_f)),
    cmp = compare_subgroup_sims(st$res_r, st$res_f,
                                expect_designs = c("resample", "fixed"))
  )
}

st <- run_survival_study()
cur <- pack(st)

if (mode == "capture") {
  saveRDS(cur, baseline_rds)
  cat("BASELINE CAPTURED ->", baseline_rds, "\n")
  quit(status = 0)
}

# --- compare: legacy byte-identity ---------------------------------------
if (!file.exists(baseline_rds))
  fail("baseline missing: run capture mode first (%s)", baseline_rds)
ref <- readRDS(baseline_rds)
for (nm in names(ref)) {
  if (!identical(ref[[nm]], cur[[nm]])) {
    if (is.list(ref[[nm]])) {
      for (f in union(names(ref[[nm]]), names(cur[[nm]])))
        cat(sprintf("  %s$%-16s identical: %s\n", nm, f,
                    identical(ref[[nm]][[f]], cur[[nm]][[f]])))
    }
    fail("legacy byte-identity broken at '%s'", nm)
  }
}
for (s in list(cur$s_r, cur$s_f)) {
  if (!is.null(s$effect) || !is.null(s$thresholds) || !is.null(s$labels))
    fail("legacy summary grew effect-aware fields")
}
if (!is.null(attr(cur$cmp, "effect")))
  fail("legacy compare grew an effect attribute")
ok("legacy summaries, prints, and compare frame byte-identical; no new fields")

# --- GLM behavior sections ------------------------------------------------
gbsg <- survival::gbsg
gbsg$y_cont <- log(gbsg$rfstime + 1)
dgm <- generate_glm_dgm(
  data = gbsg, factor_vars = c("meno", "grade"),
  continuous_vars = c("age", "size"), outcome_var = "y_cont",
  treatment_var = "hormon", outcome_type = "continuous",
  model = "null", subgroup_vars = NULL, n_super = 800L, seed = 99L,
  verbose = FALSE)
age_med <- stats::median(dgm$df_super$age)
sgs <- list(
  list(id = "flag_itt == 1",  name = "All Patients", grp = "ITT"),
  list(id = "meno == 1",      name = "Post-meno",    grp = "Clinical"),
  list(id = "meno == 0",      name = "Pre-meno",     grp = "Clinical"),
  list(id = "age <= age_med", name = "Age young",    grp = "Continuous"),
  list(id = "age > age_med",  name = "Age older",    grp = "Continuous"),
  list(id = "random60 == 1",  name = "random60", grp = "Random (N~60)"),
  list(id = "random12 == 1",  name = "random12", grp = "Random (N~12)"),
  list(id = "random3 == 1",   name = "random3",  grp = "Random (N~3)"))
res_g <- run_subgroup_sims(
  dgm = dgm, subgroups = sgs, n_sims = 40L, n = 180L,
  cutpoints = list(age_med = age_med),
  benchmarks = benchmark_spec(sizes = c(60L, 12L, 3L)), min_n = 5L,
  workers = 1, hr_true = dgm$hazard_ratios$overall, verbose = FALSE)

s <- summary(res_g)
want <- c("Subgroup", "N", "#(% converged)",
          "Pr(MD<lo)", "Pr(MD>0)", "mMD",
          "Pr(UB(MD)>=t1)", "Pr(UB(MD)>=t2)", "mUB(MD)")
if (!identical(colnames(s$results_tbl), want))
  fail("GLM headers wrong: [%s]",
       paste(colnames(s$results_tbl), collapse = " | "))
hand_gt0 <- apply(res_g$sim_hrs, 2, function(x) mean(x > 0, na.rm = TRUE))
if (!identical(s$pr_hr_gt1, hand_gt0)) fail("Pr(MD>0) != hand")
if (!all(is.na(s$pr_hr_lt050)) || !all(is.na(s$pr_ub_ge2)))
  fail("disabled tails not all-NA")
if (!all(s$results_tbl[["Pr(MD<lo)"]] == "-"))
  fail("disabled tail column not '-'")
if (!(s$highrisk$n == 1L &&
      identical(s$sg_names[s$highrisk$include], "All Patients")))
  fail("high-risk panel should be ITT-only under NA UB tails")
if (is.null(s$effect) || !identical(s$labels$est, "MD"))
  fail("effect-aware fields missing on the GLM summary")
po <- utils::capture.output(print(s))
if (!grepl("true MD:", po[1], fixed = TRUE))
  fail("GLM print header lacks the MD label: %s", po[1])
ok("GLM summary: headers, hand tails, '-' columns, highrisk, extras, print")

s2 <- summary(res_g, ub_thresholds = c(0.3, 0.6))
hand_ub <- apply(res_g$sim_ubs, 2, function(x) mean(x >= 0.3, na.rm = TRUE))
if (!identical(s2$pr_ub_ge2, hand_ub) ||
    !identical(colnames(s2$results_tbl)[7], "Pr(UB(MD)>=0.3)"))
  fail("explicit ub override wrong")
s3 <- summary(st$res_r, ub_thresholds = c(1.5, 2.5))
if (!identical(colnames(s3$results_tbl)[4:9],
               c("Pr(HR<0.5)", "Pr(HR>1.0)", "mHR",
                 "Pr(UB>=1.5)", "Pr(UB>=2.5)", "mUB")))
  fail("survival override headers drifted: [%s]",
       paste(colnames(s3$results_tbl)[4:9], collapse = " | "))
ok("overrides: GLM thresholds hand-verified; survival keeps legacy literals")

nan2na <- function(v) { v[is.nan(v)] <- NA_real_; v }
cmp_g <- compare_subgroup_sims(res_g, res_g)
if (!identical(unname(cmp_g$hr1_r), unname(nan2na(100 * hand_gt0))))
  fail("GLM compare hr1 != hand")
if (!all(is.na(cmp_g$ub2_r)))
  fail("GLM compare disabled UB tail not all-NA")
if (!identical(attr(cmp_g, "effect")$measure, "MD"))
  fail("GLM compare missing effect attribute")
# Guard order: the fixtures have DIFFERENT panels, so this also proves
# the metadata guard runs before the panel-alignment guard (4.2r2).
e <- tryCatch({ compare_subgroup_sims(res_g, st$res_r); NULL },
              error = function(e) conditionMessage(e))
if (is.null(e) || !grepl("incompatible effect", e))
  fail("MD-vs-HR compare should hard-error first, got: %s",
       if (is.null(e)) "<no error>" else e)
# Same metadata (legacy pair), renamed panel: the panel guard still fires.
res_rn <- st$res_r
colnames(res_rn$sim_hrs)[2] <- "RENAMED"
e2 <- tryCatch({ compare_subgroup_sims(res_rn, st$res_f); NULL },
               error = function(e) conditionMessage(e))
if (is.null(e2) || !grepl("Subgroup names differ", e2))
  fail("panel-alignment guard lost: %s",
       if (is.null(e2)) "<no error>" else e2)
ok("GLM compare: hand-checked, NA tails, effect attr, guard order (metadata first, panel second)")

cat("\nACCEPT: Phase 4.2 summary/compare gates all pass.\n")
