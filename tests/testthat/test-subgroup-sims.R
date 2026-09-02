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
