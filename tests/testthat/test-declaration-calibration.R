# =============================================================================
# fs_declaration_calibration(): calibrated declaration threshold (kappa_hat)
# and family-wise size of the p-star screen (fw_size).
# TASK_declaration_calibration_2026-09-22_v2, Step 3 (tests 1-8) and the
# field-level post-conditions of Step 4 (PC4-PC8).  PC1-PC3 (baseline digest
# identity, formals, diff scope) are checked against the pre-change pin by a
# separate script; see the report.
# =============================================================================

.md5 <- function(x) {
  f <- tempfile()
  on.exit(unlink(f))
  writeBin(serialize(x, NULL, version = 3), f)
  unname(tools::md5sum(f))
}

.strip_timing <- function(x) {
  if (is.list(x) && !is.data.frame(x)) {
    drop <- names(x) %in% c("timing_seconds", "minutes_all", "time_search")
    if (any(drop)) x <- x[!drop]
    for (i in seq_along(x)) if (!is.null(x[[i]])) x[[i]] <- .strip_timing(x[[i]])
  }
  x
}

# A bare MR-result shaped list around a synthetic field, so the exported
# function can be exercised on fields built by the internal helper alone.
.decl_fit <- function(fld, bh, p_star = 0.90, c_cons = 0, c_screen = NULL) {
  nm <- paste0("g", seq_along(bh))
  list(declaration_field = list(
    Mstar = fld$Mstar, Zstar = fld$Zstar,
    beta_hat = stats::setNames(bh, nm),
    sigma_D = stats::setNames(fld$sigma_D, nm),
    family_id = nm,
    meta = list(multiplier = "gaussian", B = fld$B, c_cons = c_cons,
                p_star = p_star, c_screen = c_screen,
                column_sd = fld$column_sd, zstar_mean = fld$zstar_mean,
                field_cor = fld$field_cor,
                shared_multipliers = fld$shared_multipliers)))
}

# Closed-form bivariate normal (Section 6.1): P(X <= a, Y <= a; rho).
.phi2 <- function(a, rho) {
  if (rho == 1) return(stats::pnorm(a))
  if (rho == 0) return(stats::pnorm(a)^2)
  stats::integrate(function(z) stats::dnorm(z) *
                     stats::pnorm((a - rho * z) / sqrt(1 - rho^2)),
                   lower = -Inf, upper = a, rel.tol = 1e-12)$value
}
.fw_target <- function(rho) 1 - .phi2(stats::qnorm(0.95), rho)
# fw_size is taken at the rounded screen: p* = 0.90 at digits = 2 admits on
# Pcons >= 0.895 (TASK_declcal_rounding_alignment_2026-09-23)
.z_eff <- stats::qnorm((1 + 0.895) / 2)
.fw_target_eff <- function(rho) 1 - .phi2(.z_eff, rho)
.kappa_target <- function(rho) {
  stats::uniroot(function(k) .phi2(k, rho) - 0.95, c(1.0, 3.0),
                 tol = 1e-10)$root
}

# ---------------------------------------------------------------------------
# Fixtures: GBSG (survival, the tested path) with the field retained, and the
# same fit with the field left unset; S1.7 configuration B (OLS, nested).
# ---------------------------------------------------------------------------

.gb <- local({
  d <- survival::gbsg
  d$id <- seq_len(nrow(d))
  d$time_months <- d$rfstime / 30.4375
  d$grade3 <- ifelse(d$grade == "3", 1, 0)
  d
})

.gb_fit <- function(mr_args) {
  suppressWarnings(forestsearch(
    .gb, outcome.name = "time_months", event.name = "status",
    treat.name = "hormon", id.name = "id",
    confounders.name = c("age", "meno", "size", "grade3", "nodes", "pgr", "er"),
    use_lasso = FALSE, use_grf = FALSE, sg_focus = "hr", maxk = 2,
    hr.threshold = 1.1, hr.consistency = 1.0, pconsistency.threshold = 0.90,
    n.min = 60, d0.min = 12, d1.min = 12, use_twostage = FALSE,
    seedit = 8316951, quiet = TRUE,
    parallel_args = list(plan = "sequential", workers = 1L),
    mr_inference = TRUE, mr_inference_args = mr_args))
}

.mr_base <- list(draws = 2000L, ci_method = "ij")
.fit_on  <- .gb_fit(utils::modifyList(.mr_base, list(
  keep_declaration_field = TRUE, keep_field_matrix = TRUE)))
.fit_off <- .gb_fit(.mr_base)

.cfgB <- local({
  set.seed(20260922L)
  n <- 400L
  A <- stats::rbinom(n, 1L, 0.5)
  X <- sample.int(4L, n, replace = TRUE)
  Y <- c(0, 0, -1, -1)[X] * A + stats::rnorm(n)      # k = 0: a tie in B
  list(df = data.frame(id = seq_len(n), Y = Y, A = A, X = X),
       cands = list(g1 = which(X == 1L), g2 = which(X <= 2L)))
})
.spec_md <- list(outcome_type = "continuous", effect_measure = "MD",
                 treat.name = "A", outcome.name = "Y", event.name = NULL,
                 offset.name = NULL, adjust_covariates = NULL,
                 adverse_outcome = TRUE)
.adm_B <- list(effect_floor = NULL,
               consistency = list(c_cons = 0, p_star = 0.90))
.B8 <- 200000L

# ---------------------------------------------------------------------------
# Test 1 -- relabelling exactness
# ---------------------------------------------------------------------------
test_that("1: with kappa_hat -> z_pstar the rule reproduces the executed screen", {
  dc <- fs_declaration_calibration(.fit_on)
  expect_true(dc$reduction$replay_check)          # replay == executed candidate set
  expect_identical(dc$reduction$n_unmatched, 0L)
  s <- dc$screened
  expect_gt(length(s), 1L)
  thr <- pmax(dc$c_screen, dc$c_cons + dc$z_pstar * dc$sigma_D[s])
  relabelled <- s[dc$beta_hat[s] >= thr]
  expect_gt(length(dc$admitted_current), 0L)
  expect_setequal(relabelled, dc$admitted_current)
  # the closed form, written out: admission <=> T >= z_{(1+p*)/2}
  rate <- pmax(0, 2 * stats::pnorm(dc$T_hat[s]) - 1)
  expect_setequal(s[rate >= dc$p_star & dc$beta_hat[s] >= dc$c_screen],
                  relabelled)
})

# ---------------------------------------------------------------------------
# Test 2 -- single candidate: kappa_hat -> qnorm(1 - alpha)
# ---------------------------------------------------------------------------
test_that("2: a single-candidate field calibrates to qnorm(1 - alpha)", {
  set.seed(101L)
  n <- 300L; B <- 20000L
  db <- matrix(stats::rnorm(n) * stats::rexp(n), n, 1L)
  xi <- matrix(stats::rnorm(n * B), n, B)
  fld <- forestsearch:::.fs_decl_field(db, xi)
  for (a in c(0.05, 0.10)) {
    dc <- fs_declaration_calibration(.decl_fit(fld, bh = 0.1), alpha = a)
    q <- stats::qnorm(1 - a)
    mcsd <- sqrt(a * (1 - a) / B) / stats::dnorm(q)
    expect_lt(abs(dc$kappa_hat - q), 4 * mcsd)
  }
})

# ---------------------------------------------------------------------------
# Test 3 -- independent columns: fw_size -> 1 - (1 - alpha1)^G
# ---------------------------------------------------------------------------
test_that("3: independent columns give fw_size = 1 - (1 - alpha1)^G", {
  set.seed(202L)
  G <- 5L; m <- 80L; n <- G * m; B <- 20000L
  db <- matrix(0, n, G)                        # disjoint supports -> independent
  for (g in seq_len(G)) db[(g - 1L) * m + seq_len(m), g] <- stats::rnorm(m)
  xi <- matrix(stats::rnorm(n * B), n, B)
  fld <- forestsearch:::.fs_decl_field(db, xi)
  # the screen as implemented: round(Pcons, 2) >= p*, i.e. Pcons >= p* - 0.005
  pcons_eff <- c("0.8" = 0.795, "0.9" = 0.895)
  for (ps in c(0.80, 0.90)) {
    dc <- fs_declaration_calibration(.decl_fit(fld, bh = rep(0, G), p_star = ps))
    alpha1 <- 1 - stats::pnorm(stats::qnorm((1 + pcons_eff[[as.character(ps)]]) / 2))
    tgt <- 1 - (1 - alpha1)^G
    expect_lt(abs(dc$fw_size - tgt), 4 * sqrt(tgt * (1 - tgt) / B))
  }
})

# ---------------------------------------------------------------------------
# Test 4 -- monotonicity on a nested family
# ---------------------------------------------------------------------------
test_that("4: fw_size >= alpha1 and kappa_hat is non-decreasing in family size", {
  set.seed(303L)
  n <- 400L; B <- 20000L; Gmax <- 6L
  # nested, overlapping supports: candidate g is the first 100 + 50 g subjects
  base <- stats::rnorm(n)
  db <- vapply(seq_len(Gmax), function(g) {
    v <- numeric(n); idx <- seq_len(100L + 50L * g); v[idx] <- base[idx]; v
  }, numeric(n))
  xi <- matrix(stats::rnorm(n * B), n, B)
  fld_all <- forestsearch:::.fs_decl_field(db, xi, keep_matrix = TRUE)
  alpha1 <- 1 - stats::pnorm(stats::qnorm(0.95))
  kap <- fw <- numeric(Gmax)
  for (G in seq_len(Gmax)) {
    fld <- forestsearch:::.fs_decl_field(db[, seq_len(G), drop = FALSE], xi)
    dc <- fs_declaration_calibration(.decl_fit(fld, bh = rep(0, G)))
    kap[G] <- dc$kappa_hat; fw[G] <- dc$fw_size
  }
  expect_true(all(diff(kap) >= 0))
  expect_true(all(diff(fw) >= 0))
  expect_true(all(fw >= alpha1 - 4 * sqrt(alpha1 * (1 - alpha1) / B)))
  expect_true(all(fw[-1] >= alpha1))
})

# ---------------------------------------------------------------------------
# Test 5 -- the pre-reduction family is conservative
# ---------------------------------------------------------------------------
test_that("5: kappa_hat(prereduction) >= kappa_hat(reduced) when the reduction bit", {
  pre <- fs_declaration_calibration(.fit_on, family = "prereduction")
  red <- fs_declaration_calibration(.fit_on, family = "reduced")
  expect_gte(length(pre$reduction$removed), 1L)
  expect_gt(pre$n_family_prereduction, pre$n_family_reduced)
  for (a in c(0.05, 0.10)) {
    expect_gte(fs_declaration_calibration(.fit_on, alpha = a)$kappa_hat,
               fs_declaration_calibration(.fit_on, alpha = a,
                                          family = "reduced")$kappa_hat)
  }
  expect_gte(pre$fw_size, red$fw_size)
  expect_match(red$family_label, "CONDITIONAL ON THE REALIZED FAMILY")
  expect_output(print(red), "conditional on the realized")
  expect_output(print(pre), "prereduction")
})

# ---------------------------------------------------------------------------
# Test 6 -- purity
# ---------------------------------------------------------------------------
test_that("6: the fit is unchanged by the call", {
  before <- .md5(.fit_on)
  invisible(fs_declaration_calibration(.fit_on))
  invisible(fs_declaration_calibration(.fit_on, alpha = 0.10, family = "reduced"))
  expect_identical(.md5(.fit_on), before)
})

# ---------------------------------------------------------------------------
# Test 7 -- default-off
# ---------------------------------------------------------------------------
test_that("7: unset, nothing is stored and the call refuses by name", {
  f <- formals(forestsearch:::fs_mr_inference)
  expect_identical(f$keep_declaration_field, FALSE)
  expect_identical(f$keep_field_matrix, FALSE)
  expect_null(.fit_off$mr_inference$declaration_field)
  # the capture perturbs nothing else: with it on, the MR result minus the
  # new element is the MR result with it off
  on <- .fit_on$mr_inference
  on$declaration_field <- NULL
  expect_identical(.strip_timing(on), .strip_timing(.fit_off$mr_inference))
  expect_error(fs_declaration_calibration(.fit_off), "keep_declaration_field")
  # the reduced diagnostic needs the matrix and says so
  no_mat <- .gb_fit(utils::modifyList(.mr_base,
                                      list(keep_declaration_field = TRUE)))
  expect_null(no_mat$mr_inference$declaration_field$Zstar)
  expect_error(fs_declaration_calibration(no_mat, family = "reduced"),
               "keep_field_matrix")
  expect_equal(fs_declaration_calibration(no_mat)$kappa_hat,
               fs_declaration_calibration(.fit_on)$kappa_hat)
})

# ---------------------------------------------------------------------------
# Test 8 -- closed-form two-candidate acceptance check (S1.7, configuration B)
# ---------------------------------------------------------------------------
.mr8 <- fs_mr_inference(
  .cfgB$df, .cfgB$cands, .spec_md, selected_members = .cfgB$cands$g1,
  admission = .adm_B, reselection = "maxeff", draws = .B8,
  multiplier = "gaussian", ci_method = "ij", seed = 7L,
  keep_declaration_field = TRUE)
.db8 <- forestsearch:::.fs_mr_assemble(.cfgB$df, .cfgB$cands, .spec_md)$B
.rho8 <- local({
  s <- sqrt(colSums(.db8^2))
  sum(.db8[, 1] * .db8[, 2]) / (s[1] * s[2])
})

.check8 <- function(dc, rho) {
  expect_lt(abs(dc$field_cor[1, 2] - rho), 4 / sqrt(.B8))              # 8a
  expect_lt(abs(dc$fw_size - .fw_target_eff(rho)), 0.0027)              # 8b
  expect_lt(abs(dc$kappa_hat - .kappa_target(rho)), 0.025)              # 8c
  expect_gte(dc$fw_size, 1 - stats::pnorm(.z_eff))                      # 8d
  expect_lte(dc$fw_size, 1 - stats::pnorm(.z_eff)^2)
  expect_gte(dc$kappa_hat, 1.644854); expect_lte(dc$kappa_hat, 1.954508)
}

test_that("8: the exported path matches the closed form (Gaussian multipliers)", {
  expect_gte(.rho8, 0.05)                           # non-degenerate: nested B
  dc <- fs_declaration_calibration(.mr8)
  expect_identical(dc$multiplier_law, "gaussian")
  expect_identical(dc$B, .B8)
  .check8(dc, .rho8)
})

test_that("8: the internal helper on a hand-supplied db matches the closed form", {
  set.seed(808L)
  xi <- matrix(stats::rnorm(nrow(.db8) * .B8), nrow(.db8), .B8)
  fld <- forestsearch:::.fs_decl_field(.db8, xi, keep_matrix = FALSE)
  dc <- fs_declaration_calibration(.decl_fit(fld, bh = c(0, 0)))
  .check8(dc, .rho8)
})

test_that("8d: the closed-form targets bracket and decrease in rho", {
  grid <- c(0, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99, 1)
  fw <- vapply(grid, .fw_target, numeric(1))
  kp <- vapply(grid, .kappa_target, numeric(1))
  expect_equal(fw[c(1, 8)], c(0.0975, 0.05), tolerance = 1e-6)
  expect_equal(kp[c(1, 8)], c(1.954508, 1.644854), tolerance = 1e-6)
  expect_true(all(diff(fw) < 0))
  expect_true(all(diff(kp) < 0))
  # the task's reference table (Section 6.1)
  expect_equal(fw, c(0.097500, 0.096287, 0.092865, 0.087811, 0.080401,
                     0.068132, 0.055811, 0.050000), tolerance = 1e-5)
  expect_equal(kp, c(1.954508, 1.950821, 1.938467, 1.916332, 1.877299,
                     1.797586, 1.698675, 1.644854), tolerance = 1e-6)
})

# ---------------------------------------------------------------------------
# Post-conditions 4-8 (Step 4)
# ---------------------------------------------------------------------------
test_that("PC4: the field is centred with unit column scale", {
  for (dc in list(fs_declaration_calibration(.fit_on),
                  fs_declaration_calibration(.mr8))) {
    expect_lt(abs(dc$zstar_mean), 4 / sqrt(dc$B))
    expect_true(all(abs(dc$column_sd - 1) < 0.10))
  }
})

test_that("PC5: one shared multiplier vector per draw; ncol(Zstar) = family size", {
  fld <- .fit_on$mr_inference$declaration_field
  expect_true(fld$meta$shared_multipliers)
  expect_identical(ncol(fld$Zstar), length(fld$family_id))
  expect_identical(ncol(fld$Zstar),
                   fs_declaration_calibration(.fit_on)$n_family_prereduction)
  red <- fs_declaration_calibration(.fit_on, family = "reduced")
  expect_identical(length(red$beta_hat), red$n_family_reduced)
  # regenerate the main stream from the recorded seed: the field is exactly
  # the shared-multiplier assembly of it
  set.seed(7L)
  xi <- forestsearch:::.fs_mr_multipliers(nrow(.db8), .B8, "gaussian")
  fld8 <- forestsearch:::.fs_decl_field(.db8, xi, keep_matrix = FALSE)
  expect_identical(fld8$Mstar, .mr8$declaration_field$Mstar)
  # element-wise: Zstar[b, g] = sum_i xi[i, b] db[i, g] / sigma_D(g)
  small <- forestsearch:::.fs_decl_field(.db8, xi[, 1:50], keep_matrix = TRUE)
  for (b in c(1L, 17L, 50L)) for (g in 1:2) {
    expect_equal(small$Zstar[b, g],
                 sum(xi[, b] * .db8[, g]) / sqrt(sum(.db8[, g]^2)))
  }
})

test_that("PC6: sigma_D is the value the closed-form screen used", {
  # GBSG: the selected candidate, through consistency_resample() with the
  # screen's own Cox warm start
  dc <- fs_declaration_calibration(.fit_on)
  ids <- .fit_on$grp.consistency$df_flag$id[.fit_on$grp.consistency$sg.harm.id == 1]
  sub <- .gb[.gb$id %in% ids, ]
  init <- survival::coxph(survival::Surv(time_months, status) ~ hormon,
                          data = sub)$coefficients[1]
  rr <- consistency_resample(sub, method = "closed", tte.name = "time_months",
                             event.name = "status", treat.name = "hormon",
                             cox_init = init)
  sel <- dc$admitted_current
  expect_equal(unname(dc$sigma_D[sel]), rr$sigma_D, tolerance = 1e-6)
  expect_equal(unname(dc$beta_hat[sel]), rr$beta_hat, tolerance = 1e-6)
  fld <- .fit_on$mr_inference$declaration_field
  expect_equal(unname(fld$sigma_D), unname(fld$meta$sigma_D_field))
  # OLS: both candidates
  d8 <- fs_declaration_calibration(.mr8)
  for (g in 1:2) {
    r8 <- consistency_resample(.cfgB$df[.cfgB$cands[[g]], ], method = "closed",
                               outcome_type = "continuous",
                               effect_measure = "MD", treat.name = "A",
                               outcome.name = "Y", consistency_threshold = 0)
    expect_equal(unname(d8$sigma_D[g]), r8$sigma_D, tolerance = 1e-10)
  }
})

test_that("PC7: n_family_prereduction >= n_family_reduced", {
  dc <- fs_declaration_calibration(.fit_on)
  expect_gte(dc$n_family_prereduction, dc$n_family_reduced)
})
