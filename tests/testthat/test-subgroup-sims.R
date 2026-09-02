# Unit tests for run_subgroup_sims() companions: benchmark_spec(),
# subgroup_cox(), validate_subgroups(), argument guards that fire before
# any simulation, and summary.subgroup_sims() on a hand-built object.
#
# The end-to-end bit-identity of run_subgroup_sims() against the
# extreme-subgroups vignettes (real GBSG DGM, both baselines, parallel ==
# sequential) is covered by dev/accept_phase1_bitident.R, which is run as
# a release gate rather than on every test invocation.

test_that("benchmark_spec stores and validates its fields", {
  b <- benchmark_spec()
  expect_s3_class(b, "benchmark_spec")
  expect_identical(b$sizes, c(60L, 40L, 20L, 15L))
  expect_true(b$nested)
  expect_identical(b$prefix, "random")
  b2 <- benchmark_spec(sizes = c(30, 10), nested = FALSE, prefix = "rnd")
  expect_identical(b2$sizes, c(30L, 10L))
  expect_false(b2$nested)
  expect_error(benchmark_spec(sizes = c(10, 10)), "distinct")
  expect_error(benchmark_spec(sizes = 0), "positive")
})

test_that("subgroup_cox reproduces coxph conf.int extraction and NA fallback", {
  set.seed(11)
  m <- 120L
  d <- data.frame(
    treat_sim = rbinom(m, 1, 0.5),
    grade     = sample(1:3, m, replace = TRUE)
  )
  tt <- rexp(m, exp(-0.4 * d$treat_sim - 1))
  cc <- rexp(m, 1 / 30)
  d$y_sim     <- pmin(tt, cc)
  d$event_sim <- as.integer(tt <= cc)

  f <- subgroup_cox(survival::Surv(y_sim, event_sim) ~ treat_sim +
                      survival::strata(grade))
  r <- f(d)
  fit <- survival::coxph(
    survival::Surv(y_sim, event_sim) ~ treat_sim + survival::strata(grade),
    data = d)
  ci <- summary(fit)$conf.int[1, ]
  expect_equal(unname(r), unname(c(ci["exp(coef)"], ci["upper .95"])))
  expect_s3_class(attr(f, "formula"), "formula")

  # Estimation failure path -> c(NA, NA), no error
  expect_identical(unname(f(d[0, ])), c(NA_real_, NA_real_))
})

test_that("validate_subgroups names offenders and passes clean panels", {
  d <- data.frame(age = c(40, 60, 70), meno = c(0L, 1L, 1L))
  cp <- list(age_med = 55)
  good <- list(
    list(id = "flag_itt == 1",  name = "All Patients", grp = "ITT"),
    list(id = "age <= age_med", name = "Young",        grp = "Continuous"),
    list(id = "random60 == 1",  name = "random60",     grp = "Random (N~60)"))
  expect_true(validate_subgroups(good, d, cp))

  bad <- c(good, list(list(id = "no_col == 1", name = "Bad", grp = "x")))
  expect_error(validate_subgroups(bad, d, cp), "Bad.*no_col")

  wrongtype <- list(list(id = "age", name = "NotLogical", grp = "x"))
  expect_error(validate_subgroups(wrongtype, d, cp), "logical")

  dup <- list(list(id = "meno == 1", name = "Same", grp = "a"),
              list(id = "meno == 0", name = "Same", grp = "b"))
  expect_error(validate_subgroups(dup, d, cp), "duplicate")

  nofield <- list(list(id = "meno == 1", name = "OnlyName"))
  expect_error(validate_subgroups(nofield, d, cp), "grp")
})

test_that("run_subgroup_sims argument guards fire before any simulation", {
  stub <- list(df_source = data.frame(age = 1:10),
               df_super  = data.frame(age = 1:20))
  sg <- list(list(id = "age > 0", name = "All", grp = "ITT"))

  expect_error(
    run_subgroup_sims(stub, sg, 2, baseline = "resample",
                      analysis_time = 84, max_entry = 24, cens_adjust = 0),
    "`n`.*required")
  expect_error(
    run_subgroup_sims(stub, sg, 2, baseline = "fixed", n = 5,
                      analysis_time = 84, max_entry = 24, cens_adjust = 0),
    "must be NULL")
  expect_error(
    run_subgroup_sims(list(df_super = data.frame(age = 1:5)), sg, 2,
                      baseline = "fixed",
                      analysis_time = 84, max_entry = 24, cens_adjust = 0),
    "df_source is missing")
  expect_error(
    run_subgroup_sims(stub, sg, 2, baseline = "resample", n = 5,
                      analysis_time = 84, max_entry = 24, cens_adjust = 0,
                      benchmarks = list(sizes = 60)),
    "benchmark_spec")
  expect_error(
    run_subgroup_sims(stub, sg, 2, fit = "not a function",
                      baseline = "resample", n = 5,
                      analysis_time = 84, max_entry = 24, cens_adjust = 0),
    "function")
  # Automatic validation stops on a bad id before simulate_from_dgm runs
  expect_error(
    run_subgroup_sims(stub, list(list(id = "nope == 1", name = "B",
                                      grp = "x")),
                      2, baseline = "resample", n = 5,
                      analysis_time = 84, max_entry = 24, cens_adjust = 0),
    "invalid subgroup definition")
})

test_that("summary.subgroup_sims reproduces the vignette statistics", {
  sg <- list(
    list(id = "a", name = "All Patients", grp = "ITT"),
    list(id = "b", name = "Cl",           grp = "Clinical"),
    list(id = "c", name = "Ix",           grp = "Interaction: meno x grade"),
    list(id = "d", name = "3w",           grp = "3-way"),
    list(id = "e", name = "random60",     grp = "Random (N~60)"),
    list(id = "f", name = "EmptySG",      grp = "Continuous"))
  hr <- matrix(c(0.7, 0.8, 0.6, 0.9,  0.4, 0.6, 1.2, 0.5,
                 0.3, 1.1, 0.45, 0.9, 0.55, 0.65, 0.75, 0.85,
                 1.4, 0.2, 0.9, NA,   NA, NA, NA, NA), 4, 6,
               dimnames = list(NULL, vapply(sg, `[[`, "", "name")))
  ub <- hr * 2.5
  ns <- matrix(rep(c(90, 45, 20, 12, 60, NA), each = 4), 4, 6,
               dimnames = dimnames(hr))
  obj <- structure(
    list(design = "fixed", n_sims = 4L, sim_hrs = hr, sim_ubs = ub,
         sim_ns = ns, subgroups = sg, hr_true = 0.7),
    class = c("subgroup_sims", "list"))

  S <- summary(obj)
  expect_s3_class(S, "subgroup_sims_summary")
  expect_equal(unname(S$pr_ub_ge2["Cl"]), mean(ub[, "Cl"] >= 2))
  expect_equal(unname(S$pr_hr_lt050["Ix"]),
               mean(hr[, "Ix"] < 0.5, na.rm = TRUE))
  expect_equal(unname(S$hr_q["All Patients", 2]),
               unname(quantile(hr[, 1], 0.5)))
  # Structurally empty column: NA median, excluded from panels, "-" in table
  expect_true(is.na(S$hr_q["EmptySG", 2]))
  expect_false(S$ok[6])
  expect_identical(S$results_tbl$N[6], "-")
  expect_identical(S$results_tbl$mHR[6], "-")
  # Convergence formatting
  expect_identical(S$results_tbl[["#(% converged)"]][5], "3 (75%)")
  # Category mapping incl. the Interaction prefix rule
  expect_identical(S$cat, c("ITT", "Clinical", "Interaction", "3-way",
                            "Random", "Continuous"))
  # Panel partition
  expect_setequal(S$sg_names[S$ok][S$idx_hr_single], c("All Patients", "Cl"))
  expect_setequal(S$sg_names[S$ok][S$idx_hr_combo], c("Ix", "3w", "random60"))
  expect_identical(S$n_single + S$n_combo, sum(S$ok))
  # High-risk panel: ITT anchored first, then descending Pr(UB>=2)
  hi_names <- S$sg_names[S$highrisk$include][S$highrisk$ord]
  expect_identical(hi_names[1], "All Patients")
  expect_true(all(diff(
    S$pr_ub_ge2[S$highrisk$include][S$highrisk$ord][-1]) <= 0))
  # Print methods run
  expect_output(print(S), "subgroup_sims_summary")
  expect_output(print(obj), "subgroup_sims")
})

# ── Phase 2: plot method and forest_height ───────────────────────────────
# Exact argument assembly per (metric, panel) is verified against the
# vignette chunk formulas in the development harness; here we cover the
# arithmetic, the summary field, real gg_forest() execution, and guards.

.make_S <- function() {
  sg <- list(
    list(id = "a", name = "All Patients", grp = "ITT"),
    list(id = "b", name = "Cl",           grp = "Clinical"),
    list(id = "c", name = "Ix",           grp = "Interaction: meno x grade"),
    list(id = "d", name = "3w",           grp = "3-way"),
    list(id = "e", name = "random60",     grp = "Random (N~60)"),
    list(id = "f", name = "EmptySG",      grp = "Continuous"))
  hr <- matrix(c(0.7, 0.8, 0.6, 0.9,  0.4, 0.6, 1.2, 0.5,
                 0.3, 1.1, 0.45, 0.9, 0.55, 0.65, 0.75, 0.85,
                 1.4, 0.2, 0.9, NA,   NA, NA, NA, NA), 4, 6,
               dimnames = list(NULL, vapply(sg, `[[`, "", "name")))
  obj <- structure(
    list(design = "fixed", n_sims = 4L, sim_hrs = hr, sim_ubs = hr * 2.5,
         sim_ns = matrix(rep(c(90, 45, 20, 12, 60, NA), each = 4), 4, 6,
                         dimnames = dimnames(hr)),
         subgroups = sg, hr_true = 0.7),
    class = c("subgroup_sims", "list"))
  summary(obj)
}

test_that("summary exposes n_per_trial (sim_config or ITT fallback)", {
  S <- .make_S()
  expect_identical(S$n_per_trial, 90L)   # legacy object -> ITT mean N
})

test_that("forest_height converts panel rows to inches", {
  S <- .make_S()
  expect_identical(forest_height(S, "single"),
                   round(S$n_single * 0.45 + 1.5, 1))
  expect_identical(forest_height(S, "combo", 0.32, 1.2),
                   round(S$n_combo * 0.32 + 1.2, 1))
  expect_identical(forest_height(S, "highrisk"),
                   round(S$highrisk$n * 0.45 + 1.5, 1))
  expect_error(forest_height(S, "nope"))
})

test_that("plot() renders every (metric, panel) with real gg_forest", {
  S <- .make_S()
  for (mm in c("hr", "ub")) for (pp in c("single", "combo", "highrisk")) {
    p <- plot(S, metric = mm, panel = pp)
    expect_true(inherits(p, "gg") || inherits(p, "patchwork"),
                info = paste(mm, pp))
  }
  # Overrides merge over defaults without duplicate-argument errors
  expect_no_error(plot(S, "ub", "combo", ref_col = "black",
                       footnote = "custom", base_size = 9))
})

test_that("plot() requires hr_true", {
  S <- .make_S(); S$hr_true <- NULL
  expect_error(plot(S, "hr", "single"), "hr_true")
})

# ── Phase 2.1: lean formula environment in subgroup_cox() ────────────────

test_that("subgroup_cox lean re-homes the formula off the calling env", {
  f <- subgroup_cox(survival::Surv(y_sim, event_sim) ~ treat_sim)
  e <- environment(attr(f, "formula"))
  expect_false(identical(e, globalenv()))
  expect_identical(parent.env(e), asNamespace("survival"))
  # The closure itself is srcref-stripped with a minimal environment,
  # so its serialized size is install-flag independent (getSrcref would
  # otherwise be non-NULL under devtools' keep-source installs).
  expect_null(utils::getSrcref(f))
  expect_identical(parent.env(environment(f)), asNamespace("survival"))
  expect_identical(ls(environment(f)), "formula")
  # lean = FALSE preserves capture
  f2 <- subgroup_cox(survival::Surv(y_sim, event_sim) ~ treat_sim,
                     lean = FALSE)
  expect_identical(environment(attr(f2, "formula")), environment())
})

test_that("lean and non-lean fits are numerically identical (incl. strata)", {
  set.seed(21)
  m <- 150L
  d <- data.frame(treat_sim = rbinom(m, 1, 0.5),
                  grade = sample(1:3, m, TRUE),
                  age = rnorm(m, 55, 8))
  tt <- rexp(m, exp(-0.5 * d$treat_sim - 1)); cc <- rexp(m, 1 / 25)
  d$y_sim <- pmin(tt, cc); d$event_sim <- as.integer(tt <= cc)
  fml <- survival::Surv(y_sim, event_sim) ~ treat_sim +
    survival::strata(grade)
  expect_identical(subgroup_cox(fml)(d),
                   subgroup_cox(fml, lean = FALSE)(d))
})

test_that("lean trades caller-symbol resolution, gracefully", {
  set.seed(22)
  m <- 150L
  d <- data.frame(treat_sim = rbinom(m, 1, 0.5), age = rnorm(m, 55, 8))
  tt <- rexp(m, exp(-0.5 * d$treat_sim - 1)); cc <- rexp(m, 1 / 25)
  d$y_sim <- pmin(tt, cc); d$event_sim <- as.integer(tt <= cc)
  cutoff <- 50
  fml <- survival::Surv(y_sim, event_sim) ~ treat_sim + I(age > cutoff)
  # lean = FALSE resolves `cutoff` from here; lean = TRUE cannot, and the
  # tryCatch fallback returns the NA pair rather than erroring.
  expect_true(all(is.finite(subgroup_cox(fml, lean = FALSE)(d))))
  expect_identical(unname(subgroup_cox(fml)(d)),
                   c(NA_real_, NA_real_))
})

test_that("lean fit serializes small even when the formula was born heavy", {
  heavy <- new.env(parent = globalenv())
  assign("ballast", matrix(rnorm(4e5), 500), envir = heavy)
  fml <- eval(quote(survival::Surv(y_sim, event_sim) ~ treat_sim),
              envir = heavy)
  lean_b  <- length(serialize(subgroup_cox(fml), NULL))
  heavy_b <- length(serialize(subgroup_cox(fml, lean = FALSE), NULL))
  expect_lt(lean_b, 100 * 1024)
  expect_gt(heavy_b, 1024 * 1024)
})

# ── Phase 3: compare_subgroup_sims() ─────────────────────────────────────

.make_pair <- function() {
  mk <- function(design, shift) {
    sg <- c("All Patients", "Cl", "EmptySG")
    hr <- matrix(c(0.7, 0.8, 0.6, 0.9,
                   0.4 + shift, 0.6, 1.2, 0.5,
                   NA, NA, NA, NA), 4, 3, dimnames = list(NULL, sg))
    list(design = design, n_sims = 4L, sim_hrs = hr, sim_ubs = hr * 2.5,
         sim_ns = matrix(rep(c(90, 45, NA), each = 4), 4, 3,
                         dimnames = dimnames(hr)))
  }
  list(r = mk("resample", 0), f = mk("fixed", 0.05))
}

test_that("compare_subgroup_sims reproduces the memo statistics", {
  p <- .make_pair()
  cmp <- compare_subgroup_sims(p$r, p$f,
                               expect_designs = c("resample", "fixed"))
  expect_identical(
    names(cmp),
    c("subgroup", "N_r", "N_f", "ub2_r", "ub2_f", "ub3_r", "ub3_f",
      "mUB_r", "mUB_f", "hr05_r", "hr05_f", "hr1_r", "hr1_f",
      "mHR_r", "mHR_f"))
  expect_identical(cmp$ub2_f[2],
                   100 * mean(p$f$sim_ubs[, "Cl"] >= 2, na.rm = TRUE))
  expect_identical(cmp$mHR_r[1], median(p$r$sim_hrs[, 1]))
  # Structurally empty: NaN -> NA in N and tail columns, NA medians
  expect_true(is.na(cmp$N_r[3]) && is.na(cmp$ub2_f[3]) &&
                is.na(cmp$mUB_r[3]))
  expect_identical(attr(cmp, "designs"),
                   c(x = "resample", y = "fixed"))
  expect_identical(attr(cmp, "n_sims"), c(x = 4L, y = 4L))
  expect_s3_class(cmp, "data.frame")   # plain frame, static-mode type
})

test_that("compare_subgroup_sims agrees with summary() fields", {
  p <- .make_pair()
  obj <- structure(c(p$f, list(subgroups = list(
    list(id = "a", name = "All Patients", grp = "ITT"),
    list(id = "b", name = "Cl", grp = "Clinical"),
    list(id = "c", name = "EmptySG", grp = "Continuous")))),
    class = c("subgroup_sims", "list"))
  S <- summary(obj, hr_true = 0.7)
  cmp <- compare_subgroup_sims(p$r, p$f)
  nn <- function(v) { v[is.nan(v)] <- NA_real_; unname(v) }
  expect_equal(cmp$ub2_f, nn(100 * S$pr_ub_ge2))
  expect_equal(cmp$hr05_f, nn(100 * S$pr_hr_lt050))
  expect_equal(cmp$mUB_f, unname(S$ub_q[, 2]))   # median == type-7 Q(.5)
})

test_that("compare_subgroup_sims guards alignment and designs", {
  p <- .make_pair()
  bad <- p$f
  colnames(bad$sim_hrs)[2] <- "Renamed"
  expect_error(compare_subgroup_sims(p$r, bad), "Subgroup names differ")
  expect_error(
    compare_subgroup_sims(p$f, p$r,
                          expect_designs = c("resample", "fixed")),
    "Design labels")
  expect_error(compare_subgroup_sims(list(sim_hrs = 1), p$f),
               "sim_hrs / sim_ubs / sim_ns")
  cmp2 <- compare_subgroup_sims(p$r, p$f, suffixes = c(".x", ".y"))
  expect_true(all(c("N.x", "mHR.y") %in% names(cmp2)))
})
