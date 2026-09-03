# accept_phase43_plot.R -----------------------------------------------------
# SENTINEL: P43-ACC-v2-2026-09-03 [delivery sentinel: p44r2-51b6d40b]
#
# Phase 4.3 acceptance gate (v2: standing at any HEAD >= 4.3). Sections:
#   A. Legacy byte-identity: the pre-4.3 plot code is PINNED to the 4.2
#      commit -- `git show 310039f6:R/plot_subgroup_sims.R` -- with the
#      extracted bytes asserted against a recorded sha256, then sourced
#      beside the repo's plot file, each with
#      gg_forest replaced by a capture stub (function(...) list(...)).
#      For a reference survival study, all twelve (baseline x metric x
#      panel) argument lists and forest_height() must be identical().
#   B. GLM behavior via the stub: identity scale, null at 0, engine-
#      delegated axes, trimmed annotation columns/widths, exact
#      footnotes, override paths.
#   C. Real-engine smoke: the INSTALLED plot() on the GLM summary must
#      build and render (ggsave to a temp PNG) for hr and ub panels --
#      the identity-scale end-to-end proof the stub cannot give.

suppressPackageStartupMessages({ library(forestsearch); library(survival) })
fail <- function(...) stop(sprintf(...), call. = FALSE)
ok   <- function(...) cat("PASS:", sprintf(...), "\n")

new_path <- "R/plot_subgroup_sims.R"
if (!file.exists(new_path)) fail("run from the package root")
if (!any(grepl("Phase 4.3 (scale-aware plots)", readLines(new_path),
               fixed = TRUE)))
  fail("repo plot file is not the placed 4.3 version")
old_path <- tempfile(fileext = ".R")
# v2 (p44r2): the legacy reference is pinned to the Phase 4.2 commit
# instead of asserting HEAD == 4.2 (structurally unrunnable at any HEAD
# past the 4.3 commit). The extraction writes raw bytes (stdout = file,
# no readLines round trip) and is self-verifying: the recorded sha256
# arbitrates even under an abbreviated-hash collision or a history
# rewrite, and strictly subsumes the old "no Phase 4.3 marker" check.
legacy_commit <- "310039f6"
legacy_sha256 <-
  "4b472d2f198d0c0a929e8a1bf3dc0616f4c973f306005ce6d8da2e4f857b3406"
full_hash <- suppressWarnings(system2(
  "git", c("rev-parse", "--verify", paste0(legacy_commit, "^{commit}")),
  stdout = TRUE, stderr = TRUE))
if (!is.null(attr(full_hash, "status")) || length(full_hash) != 1L)
  fail("legacy commit %s not found in this repository", legacy_commit)
cat("legacy reference commit:", full_hash, "\n")
st <- system2("git",
              c("show", paste0(legacy_commit, ":R/plot_subgroup_sims.R")),
              stdout = old_path)
if (!identical(st, 0L)) fail("git show failed for the legacy reference")
if (!requireNamespace("digest", quietly = TRUE))
  fail("the 'digest' package is required for the legacy content check")
got_sha <- digest::digest(file = old_path, algo = "sha256")
if (!identical(got_sha, legacy_sha256))
  fail("legacy reference content hash mismatch: got %s", got_sha)
ok("legacy reference pinned to %s (sha256 verified)", full_hash)

mk_env <- function(path) {
  e <- new.env(parent = globalenv())
  assign("gg_forest", function(...) list(...), envir = e)
  sys.source(path, envir = e)
  e
}
envO <- mk_env(old_path)
envN <- mk_env(new_path)

# --- reference survival study (installed runner + summary) ----------------
gbsg <- survival::gbsg
gbsg$time_months <- gbsg$rfstime / 30.4375
set.seed(99)
dgm <- generate_aft_dgm_flex(
  data = gbsg, continuous_vars = c("age", "size", "nodes", "pgr", "er"),
  factor_vars = c("meno", "grade"), outcome_var = "time_months",
  event_var = "status", treatment_var = "hormon", subgroup_vars = NULL,
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
               benchmarks = benchmark_spec(), min_n = 5L, workers = 1,
               hr_true = 0.70, k_treat = 0.95, verbose = FALSE)
res_r <- do.call(run_subgroup_sims,
                 c(common, list(baseline = "resample", n = 200L)))
res_f <- do.call(run_subgroup_sims, c(common, list(baseline = "fixed")))

for (nm in c("resample", "fixed")) {
  res <- if (nm == "resample") res_r else res_f
  S <- summary(res)
  for (m in c("hr", "ub")) for (p in c("single", "combo", "highrisk")) {
    a <- envO$plot.subgroup_sims_summary(S, metric = m, panel = p)
    b <- envN$plot.subgroup_sims_summary(S, metric = m, panel = p)
    if (!identical(a, b)) {
      for (f in union(names(a), names(b)))
        cat(sprintf("  %-12s identical: %s\n", f,
                    identical(a[[f]], b[[f]])))
      fail("A. legacy args differ: %s/%s/%s", nm, m, p)
    }
  }
  if (!identical(envO$forest_height(S, "combo"),
                 envN$forest_height(S, "combo")))
    fail("A. forest_height drifted (%s)", nm)
}
ok("A. legacy: all 12 gg_forest argument lists + forest_height identical")

# --- B. GLM behavior via the stub -----------------------------------------
gbsg$y_cont <- log(gbsg$rfstime + 1)
dgm_g <- generate_glm_dgm(
  data = gbsg, factor_vars = c("meno", "grade"),
  continuous_vars = c("age", "size"), outcome_var = "y_cont",
  treatment_var = "hormon", outcome_type = "continuous", model = "null",
  subgroup_vars = NULL, n_super = 800L, seed = 99L, verbose = FALSE)
age_med <- stats::median(dgm_g$df_super$age)
sgs_g <- sgs[1:7]
sgs_g <- c(sgs_g, list(
  list(id = "random60 == 1", name = "random60", grp = "Random (N~60)"),
  list(id = "random12 == 1", name = "random12", grp = "Random (N~12)"),
  list(id = "random3 == 1",  name = "random3",  grp = "Random (N~3)")))
res_g <- run_subgroup_sims(
  dgm = dgm_g, subgroups = sgs_g, n_sims = 40L, n = 180L,
  cutpoints = list(age_med = age_med),
  benchmarks = benchmark_spec(sizes = c(60L, 12L, 3L)), min_n = 5L,
  workers = 1, hr_true = dgm_g$hazard_ratios$overall, verbose = FALSE)
Sg <- summary(res_g)
ht <- sprintf("%.2f", res_g$hr_true)

a <- envN$plot.subgroup_sims_summary(Sg, metric = "hr", panel = "single")
if (!identical(a$xlog, FALSE) || !identical(a$vert_lines, 0) ||
    !identical(a$xlab, "MD"))
  fail("B. GLM hr: xlog/vert/xlab wrong")
if ("xlim" %in% names(a) || "ticks_at" %in% names(a))
  fail("B. GLM hr: xlim/ticks should be engine-delegated (absent)")
if (!identical(names(a$annot), c("N", "MD>0", "mMD")) ||
    !identical(a$widths, c(3.5, 5, 1.1, 1.2, 1.1)))
  fail("B. GLM hr: annot/widths wrong")
if (!identical(a$footnote,
      paste0("Median MD (point) + 1st-99th percentile ECI  |  ",
             "40 simulated trials, N = 180 per trial  |  ",
             "True MD = ", ht, " (red dashed)")))
  fail("B. GLM hr footnote: %s", a$footnote)
b <- envN$plot.subgroup_sims_summary(Sg, metric = "ub", panel = "highrisk")
if (!grepl("^High-risk filter disabled", b$footnote) ||
    !identical(b$subgroups, "All Patients"))
  fail("B. GLM highrisk: %s", b$footnote)
S3 <- summary(res_r, ub_thresholds = c(1.5, 2.5))
d <- envN$plot.subgroup_sims_summary(S3, metric = "ub", panel = "single")
if (!identical(d$xlog, TRUE) || !identical(d$xlim, c(0.30, 9.0)) ||
    !identical(names(d$annot), c("N", "UB>=1.5", "UB>=2.5", "mUB")))
  fail("B. survival override: ratio constants / annot wrong")
ok("B. GLM + override behavior via the stub: all assertions hold")

# --- C. real-engine render smoke ------------------------------------------
for (m in c("hr", "ub")) {
  p <- tryCatch(plot(Sg, metric = m, panel = "single"),
                error = function(e) e)
  if (inherits(p, "error"))
    fail("C. installed plot(%s) errored: %s", m, conditionMessage(p))
  f <- tempfile(fileext = ".png")
  r <- tryCatch({ ggplot2::ggsave(f, p, width = 8, height = 5, dpi = 72)
                  TRUE },
                error = function(e) conditionMessage(e))
  if (!isTRUE(r)) fail("C. render (%s) failed: %s", m, r)
  if (!file.exists(f) || file.size(f) < 1000)
    fail("C. render (%s) produced no usable file", m)
}
ok("C. real gg_forest renders identity-scale hr and ub panels")

cat("\nACCEPT: Phase 4.3 plot gates all pass.\n")
