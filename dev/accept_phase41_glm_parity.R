# accept_phase41_glm_parity.R ---------------------------------------------
# SENTINEL: P41-GLMPAR-v3-2026-09-03 [delivery sentinel: c1r1-483d82b8]
#
# Phase 4.1 GLM acceptance gate, run against the INSTALLED Phase-4.1
# build (library(forestsearch); never load_all -- doFuture workers need
# the installed package). Contents:
#   A. Runner vs transcribed hand loop: a real glm_dgm from
#      generate_glm_dgm() (gbsg-derived continuous outcome, null model),
#      run_subgroup_sims() with the default subgroup_glm() fit versus an
#      INDEPENDENTLY written per-subgroup slow chain
#      (lm + coef + vcov, z upper, NA guards) over the same simulate/
#      benchmark/seed transcript -- identical() on all three matrices.
#      This simultaneously certifies the estimator layer's fast-path
#      bit-identity contract at scale.
#   B. Reused-plan parallel run identical to sequential; plan_reused
#      recorded.
#   C. Degenerate behavior in-loop: the size-3 benchmark is all-NA via
#      min_n; single-arm draws in the size-6 benchmark yield NA fits.
#   D. Result contract: design = "glm", effect metadata, NULL
#      cens_adjust, no analysis_time in sim_config, fit provenance.
#   E. Policing: five specced errors + the DGM/fitter mismatch warning.
#   F. Shipment sanity: serialized fitter < 512 KB; lean formula env.

suppressPackageStartupMessages({
  library(forestsearch); library(survival); library(future)
})
stopifnot("subgroup_glm" %in% getNamespaceExports("forestsearch"))
cat(sprintf("installed forestsearch %s\n",
            as.character(utils::packageVersion("forestsearch"))))
fail <- function(...) stop(sprintf(...), call. = FALSE)
ok   <- function(...) cat("PASS:", sprintf(...), "\n")

# --- A. real glm_dgm + runner --------------------------------------------
data(gbsg, package = "survival")
gbsg$y_cont <- log(gbsg$rfstime + 1)          # continuous stand-in outcome
dgm <- generate_glm_dgm(
  data = gbsg,
  factor_vars     = c("meno", "grade"),
  continuous_vars = c("age", "size"),
  outcome_var     = "y_cont",
  treatment_var   = "hormon",
  outcome_type    = "continuous",
  model           = "null",
  subgroup_vars   = NULL,
  n_super         = 800L,
  seed            = 99L,
  verbose         = FALSE)
stopifnot(inherits(dgm, "glm_dgm"),
          identical(dgm$effect_measure, "MD"))

age_med <- stats::median(dgm$df_super$age)
sgs <- list(
  list(id = "flag_itt == 1",      name = "All Patients", grp = "ITT"),
  list(id = "meno == 1",          name = "Post-meno",    grp = "Clinical"),
  list(id = "meno == 0",          name = "Pre-meno",     grp = "Clinical"),
  list(id = "grade == 3",         name = "Grade 3",      grp = "Clinical"),
  list(id = "age <= age_med",     name = "Age young",    grp = "Continuous"),
  list(id = "age > age_med",      name = "Age older",    grp = "Continuous"),
  list(id = "meno==1 & grade==3", name = "Post/G3", grp = "Interaction: mg"),
  list(id = "random60 == 1",      name = "random60", grp = "Random (N~60)"),
  list(id = "random12 == 1",      name = "random12", grp = "Random (N~12)"),
  list(id = "random6 == 1",       name = "random6",  grp = "Random (N~6)"),
  list(id = "random3 == 1",       name = "random3",  grp = "Random (N~3)"))
cuts  <- list(age_med = age_med)
bench <- benchmark_spec(sizes = c(60L, 12L, 6L, 3L))
n_sims <- 40L; n_trial <- 180L; seed_base <- 0L; off <- 1e6L
hr_true <- dgm$hazard_ratios$overall

res <- run_subgroup_sims(
  dgm = dgm, subgroups = sgs, n_sims = n_sims, n = n_trial,
  cutpoints = cuts, benchmarks = bench, min_n = 5L,
  workers = 1, hr_true = hr_true, verbose = FALSE)

# Transcribed hand loop with an independent slow-chain fitter.
hfit <- function(d) {
  fit <- tryCatch(stats::lm(y_sim ~ treat_sim, data = d),
                  error = function(e) NULL)
  if (is.null(fit)) return(c(NA_real_, NA_real_))
  cf <- stats::coef(fit)
  if (!("treat_sim" %in% names(cf)) || is.na(cf[["treat_sim"]]))
    return(c(NA_real_, NA_real_))
  se <- tryCatch(sqrt(diag(stats::vcov(fit)))[["treat_sim"]],
                 error = function(e) NA_real_)
  if (is.na(se)) return(c(NA_real_, NA_real_))
  out <- c(cf[["treat_sim"]],
           cf[["treat_sim"]] + stats::qnorm(0.975) * se)
  if (!all(is.finite(out))) return(c(NA_real_, NA_real_))
  unname(out)
}
sg_ids <- vapply(sgs, `[[`, character(1L), "id")
sg_nm  <- vapply(sgs, `[[`, character(1L), "name")
n_sg <- length(sgs)
H <- U <- Nn <- matrix(NA_real_, n_sims, n_sg, dimnames = list(NULL, sg_nm))
for (ss in seq_len(n_sims)) {
  df_s <- simulate_from_glm_dgm(dgm, n = n_trial, seed = seed_base + ss)
  df_s$flag_itt <- 1L
  for (nm in names(cuts)) df_s[[nm]] <- cuts[[nm]]
  set.seed(seed_base + ss + off)
  m <- nrow(df_s)
  r_idx <- sample.int(m, min(max(bench$sizes), m), replace = FALSE)
  for (s in bench$sizes) {
    df_s[[paste0("random", s)]] <-
      as.integer(seq_len(m) %in% r_idx[seq_len(min(s, length(r_idx)))])
  }
  for (gg in seq_len(n_sg)) {
    df_sg <- tryCatch(subset(df_s, eval(parse(text = sg_ids[[gg]]))),
                      error = function(e) NULL)
    if (is.null(df_sg) || nrow(df_sg) < 5L) next
    Nn[ss, gg] <- nrow(df_sg)
    r <- hfit(df_sg)
    H[ss, gg] <- r[1]; U[ss, gg] <- r[2]
  }
}
if (!identical(res$sim_hrs, H)) fail("GLM parity: sim_hrs differ")
if (!identical(res$sim_ubs, U)) fail("GLM parity: sim_ubs differ")
if (!identical(res$sim_ns,  Nn)) fail("GLM parity: sim_ns differ")
ok("A. runner vs hand loop: all three matrices identical (n_sims=%d)",
   n_sims)

# --- B. reused-plan parallel ---------------------------------------------
old_plan <- future::plan()
future::plan(future::multisession,
             workers = max(2L, min(4L, future::availableCores())))
res_par <- run_subgroup_sims(
  dgm = dgm, subgroups = sgs, n_sims = n_sims, n = n_trial,
  cutpoints = cuts, benchmarks = bench, min_n = 5L,
  workers = NULL, hr_true = hr_true, verbose = FALSE)
future::plan(old_plan)
if (!identical(res_par$sim_hrs, H) || !identical(res_par$sim_ubs, U) ||
    !identical(res_par$sim_ns, Nn))
  fail("parallel reused-plan run differs from sequential")
if (!isTRUE(res_par$sim_config$plan_reused))
  fail("plan_reused not recorded on the reuse path")
ok("B. reused-plan parallel identical; plan_reused = TRUE")

# --- C. degenerates in-loop ----------------------------------------------
if (!all(is.na(res$sim_hrs[, "random3"])))
  fail("min_n skip: random3 should be all-NA (size 3 < min_n 5)")
n_na6 <- sum(is.na(res$sim_hrs[, "random6"]) & !is.na(res$sim_ns[, "random6"]))
ok("C. random3 all-NA by min_n; random6 single-arm NA fits observed = %d",
   n_na6)

# --- D. result contract ---------------------------------------------------
stopifnot(identical(res$design, "glm"),
          is.null(res$cens_adjust),
          is.null(res$sim_config$analysis_time),
          is.null(res$sim_config$max_entry),
          identical(res$sim_config$baseline, "glm"),
          identical(res$sim_config$n_per_trial, n_trial),
          identical(res$effect$measure, "MD"),
          identical(res$effect$outcome_type, "continuous"),
          identical(res$effect$log_scale, FALSE),
          identical(res$effect$null_value, 0),
          identical(res$effect$ub_label, "UB(MD)"),
          identical(res$sim_config$fit, "y_sim ~ treat_sim"),
          identical(res$hr_true, hr_true))
ok("D. GLM result contract holds (design/effect/cens_adjust/sim_config)")

# --- E. policing ----------------------------------------------------------
`%||%` <- function(a, b) if (is.null(a)) b else a
expect_err <- function(expr, pat, what) {
  e <- tryCatch({ expr; NULL }, error = function(e) conditionMessage(e))
  if (is.null(e) || !grepl(pat, e))
    fail("%s: expected error /%s/, got: %s", what, pat, e %||% "<none>")
  ok("E. %s errors as specced", what)
}
# v3 (arc C' step C1): the "glm + baseline=fixed is survival-only"
# assertion is retired -- C1 implements fixed for glm_dgm, and the
# fixed-baseline contract (strict wrapper n-rule, df_source guards,
# panel invariants) is owned by dev/accept_glm_fixed_baseline.R.
expect_err(run_subgroup_sims(dgm, sgs, 2L, n = 50L, analysis_time = 84,
                             workers = 1),
           "survival-only arguments", "glm + analysis_time")
expect_err(run_subgroup_sims(dgm, sgs, 2L, workers = 1),
           "required for a \"glm_dgm\"", "glm + n=NULL")
# v2 (p44r2): Phase 4.4 lifts the binary gate -- the former error
# expectation becomes a constructs-assertion (measure/log_scale only;
# threshold stamps are accept_phase44 G6's territory).
fb <- subgroup_glm(outcome_type = "binary")
if (!is.function(fb) ||
    !identical(attr(fb, "effect")$measure, "OR") ||
    !isTRUE(attr(fb, "effect")$log_scale))
  fail("subgroup_glm(binary) must construct an OR/log-scale fitter")
ok("E. subgroup_glm binary constructs (OR, log-scale) per Phase 4.4")
expect_err(subgroup_glm(effect_measure = "OR"),
           "not available", "subgroup_glm OR-on-continuous")
# E-v2 (p44r2): under the Phase 4.4 DGM-aware default fit, a DGM/fitter
# measure mismatch on the DEFAULT path can no longer reach the runner's
# warning -- construction stops at subgroup_glm() validation. Both
# halves are asserted: the default path stops, and the runner's
# "fitter/DGM pairing" warning still fires for an EXPLICITLY passed
# mismatched fitter.
expect_err({
  dgm2 <- dgm; dgm2$effect_measure <- "OR"
  run_subgroup_sims(dgm2, sgs, 2L, n = 50L, workers = 1,
                    cutpoints = cuts, benchmarks = bench)
}, "not available for outcome_type", "default-fit measure mismatch stops")
w <- tryCatch({
  dgm2 <- dgm; dgm2$effect_measure <- "OR"
  run_subgroup_sims(dgm2, sgs, 2L, n = 50L, fit = subgroup_glm(),
                    workers = 1, cutpoints = cuts, benchmarks = bench)
  NULL
}, warning = function(w) conditionMessage(w))
if (is.null(w) || !grepl("fitter/DGM pairing", w))
  fail("DGM/fitter measure-mismatch warning missing")
ok("E. default-fit mismatch stops; explicit-fitter mismatch warns")

# --- F. shipment sanity ---------------------------------------------------
fg <- subgroup_glm()
sz <- length(serialize(fg, NULL))
if (sz > 512 * 1024) fail("subgroup_glm serializes %.1f KB", sz / 1024)
if (identical(environment(attr(fg, "formula")), globalenv()))
  fail("provenance formula bound to globalenv")
ok("F. subgroup_glm serializes %.1f KB; formula environment lean",
   sz / 1024)

cat("\nACCEPT: Phase 4.1 GLM parity, contract, policing, and shipment ",
    "gates all pass.\n", sep = "")
