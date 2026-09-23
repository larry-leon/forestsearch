# =============================================================================
# Protected null level c0 for the calibrated declaration threshold.
# TASK_declcal_c0_rchange_2026-09-22, Step 3 (tests 1-7) and the reduction
# post-condition PC5 (GBSG and continuous GLM).  PC1-PC3 (baseline digests,
# formals, diff scope) are checked by dev/verification/postcond_fits.R against
# the pre-change pin; see the report.
# =============================================================================

.md5c <- function(x) {
  f <- tempfile()
  on.exit(unlink(f))
  writeBin(serialize(x, NULL, version = 3), f)
  unname(tools::md5sum(f))
}

.strip_t <- function(x) {
  if (is.list(x) && !is.data.frame(x)) {
    drop <- names(x) %in% c("timing_seconds", "minutes_all", "time_search")
    if (any(drop)) x <- x[!drop]
    for (i in seq_along(x)) if (!is.null(x[[i]])) x[[i]] <- .strip_t(x[[i]])
  }
  x
}

# Base-R bivariate normal CDF with unequal margins, upper limit a:
# Phi2(a, b; rho) = int_{-Inf}^{a} dnorm(x) pnorm((b - rho x) / sqrt(1 - rho^2)) dx
.phi2ab <- function(a, b, rho) {
  stats::integrate(function(x) stats::dnorm(x) *
                     stats::pnorm((b - rho * x) / sqrt(1 - rho^2)),
                   lower = -Inf, upper = a, rel.tol = 1e-12)$value
}

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------
.c0_grid <- c(0.70, 0.75, 0.80, 0.85, 1.0)      # c2 = hr.consistency = 1.0

.gbc <- local({
  d <- survival::gbsg
  d$id <- seq_len(nrow(d))
  d$time_months <- d$rfstime / 30.4375
  d$grade3 <- ifelse(d$grade == "3", 1, 0)
  d
})

.gbc_fit <- function(mr_args) {
  suppressWarnings(forestsearch(
    .gbc, outcome.name = "time_months", event.name = "status",
    treat.name = "hormon", id.name = "id",
    confounders.name = c("age", "meno", "size", "grade3", "nodes", "pgr", "er"),
    use_lasso = FALSE, use_grf = FALSE, sg_focus = "hr", maxk = 2,
    hr.threshold = 1.1, hr.consistency = 1.0, pconsistency.threshold = 0.90,
    n.min = 60, d0.min = 12, d1.min = 12, use_twostage = FALSE,
    seedit = 8316951, quiet = TRUE,
    parallel_args = list(plan = "sequential", workers = 1L),
    mr_inference = TRUE, mr_inference_args = mr_args))
}

.mr_b <- list(draws = 2000L, ci_method = "ij", keep_declaration_field = TRUE)
.fit_c0   <- .gbc_fit(utils::modifyList(.mr_b, list(keep_field_matrix = TRUE,
                                                    declaration_c0 = .c0_grid)))
.fit_noc0 <- .gbc_fit(utils::modifyList(.mr_b, list(keep_field_matrix = TRUE)))
.fit_bare <- .gbc_fit(.mr_b)                   # field, no matrix, no c0

# Continuous GLM (identity scale, MD): c2 = hr.consistency = 0.25
.cont_fit <- local({
  cd <- .make_continuous_data(N = 400L, MD_harm = 2)
  cargs <- .fs_args_for("continuous",
    confounders = c("age", "biomarker", "biomarker_hi", "sex"),
    extra = list(use_grf = FALSE, use_lasso = FALSE, quiet = TRUE,
                 hr.threshold = 0.5, hr.consistency = 0.25,
                 mr_inference = TRUE,
                 mr_inference_args = list(draws = 300L,
                                          keep_declaration_field = TRUE,
                                          keep_field_matrix = TRUE,
                                          declaration_c0 = c(0, 0.1, 0.25))))
  suppressWarnings(do.call(forestsearch, c(list(df.analysis = cd), cargs)))
})

# S1.7 configuration B (OLS mean difference, nested g1 in g2), c2 = 0
.cfgB0 <- local({
  set.seed(20260922L)
  n <- 400L
  A <- stats::rbinom(n, 1L, 0.5)
  X <- sample.int(4L, n, replace = TRUE)
  Y <- c(0, 0, -1, -1)[X] * A + stats::rnorm(n)
  list(df = data.frame(id = seq_len(n), Y = Y, A = A, X = X),
       cands = list(g1 = which(X == 1L), g2 = which(X <= 2L)))
})
.spec_md0 <- list(outcome_type = "continuous", effect_measure = "MD",
                  treat.name = "A", outcome.name = "Y", event.name = NULL,
                  offset.name = NULL, adjust_covariates = NULL,
                  adverse_outcome = TRUE)
.adm_B0 <- list(effect_floor = NULL,
                consistency = list(c_cons = 0, p_star = 0.90))

# ---------------------------------------------------------------------------
# Test 1 -- reduction to the current construction at c0 = c2 (and PC5)
# ---------------------------------------------------------------------------
.check_reduction <- function(fit, c0_all) {
  fld <- fit$mr_inference$declaration_field
  meta <- fld$meta
  c2 <- if (isTRUE(meta$log_scale)) exp(meta$c_cons) else meta$c_cons
  key <- as.character(c2)
  expect_true(key %in% colnames(fld$Mstar_c0))
  expect_identical(unname(fld$Mstar_c0[, key]), fld$Mstar)
  d0 <- fs_declaration_calibration(fit)
  d1 <- fs_declaration_calibration(fit, c0 = c0_all)
  row <- match(key, as.character(d1$c0$table$c0))
  expect_identical(d1$c0$source, "capture")
  expect_true(d1$c0$table$is_c2[row])
  expect_identical(d1$c0$table$kappa_hat[row], d0$kappa_hat)
  expect_identical(d1$c0$table$fw_size[row], d0$fw_size)
  expect_identical(d1$c0$admitted_calibrated[[key]], d0$admitted_calibrated)
  # the c0 block is purely additive
  d1$c0 <- NULL
  expect_identical(d1, d0)
}

test_that("1: at c0 = c2 every quantity is the unshifted one (GBSG)", {
  .check_reduction(.fit_c0, .c0_grid)
})

test_that("1 / PC5: at c0 = c2 every quantity is the unshifted one (continuous GLM)", {
  meta <- .cont_fit$mr_inference$declaration_field$meta
  expect_false(isTRUE(meta$log_scale))
  expect_identical(meta$c_cons, 0.25)
  .check_reduction(.cont_fit, c(0, 0.1, 0.25))
})

test_that("1: the field-matrix route reproduces the capture exactly", {
  cap <- fs_declaration_calibration(.fit_c0, c0 = .c0_grid)
  mat <- fs_declaration_calibration(.fit_noc0, c0 = .c0_grid)
  expect_identical(mat$c0$source, "field_matrix")
  expect_identical(unname(mat$c0$Mstar_c0), unname(cap$c0$Mstar_c0))
  expect_identical(mat$c0$table, cap$c0$table)
})

# ---------------------------------------------------------------------------
# Test 2 -- monotonicity in c0
# ---------------------------------------------------------------------------
test_that("2: kappa_hat and fw_size are non-decreasing in c0", {
  for (a in c(0.05, 0.10)) {
    tb <- fs_declaration_calibration(.fit_c0, alpha = a, c0 = .c0_grid)$c0$table
    expect_identical(tb$c0, .c0_grid)
    expect_true(all(diff(tb$kappa_hat) >= 0))
    expect_true(all(diff(tb$fw_size) >= 0))
    expect_true(all(diff(tb$n_admitted_calibrated) <= 0))
  }
})

# ---------------------------------------------------------------------------
# Test 3 -- shift correctness on a hand-supplied db
# ---------------------------------------------------------------------------
test_that("3: Mstar_c0 is the row maximum of Zstar - delta, element-wise", {
  set.seed(4040L)
  n <- 200L; B <- 3000L; G <- 4L
  db <- matrix(0, n, G)
  for (g in seq_len(G)) {
    idx <- seq_len(60L + 35L * g)                   # nested supports
    db[idx, g] <- stats::rnorm(length(idx))
  }
  xi <- matrix(stats::rnorm(n * B), n, B)
  sig <- sqrt(colSums(db^2))
  c_cons <- log(1.0)
  c0 <- c(0.70, 0.85, 1.0)
  c0_cmp <- forestsearch:::.fs_decl_c0_cmp(c0, c_cons, log_scale = TRUE)
  delta <- forestsearch:::.fs_decl_c0_shift(c0_cmp, c_cons, sig)
  fld <- forestsearch:::.fs_decl_field(db, xi, keep_matrix = TRUE,
                                       shift = delta)
  expect_identical(colnames(fld$Mstar_shift), as.character(c0))
  for (k in seq_along(c0)) {
    dk <- (c_cons - log(c0[k])) / sig                # independent construction
    expect_equal(delta[, k], dk)
    ref <- apply(fld$Zstar - rep(dk, each = B), 1L, max)
    expect_equal(fld$Mstar_shift[, k], ref, tolerance = 0)
  }
  expect_identical(fld$Mstar_shift[, "1"], fld$Mstar)
  # a vector shift is the one-column case; no shift leaves the return as before
  v <- forestsearch:::.fs_decl_field(db, xi, shift = delta[, 1])
  expect_identical(unname(v$Mstar_shift[, 1]), unname(fld$Mstar_shift[, 1]))
  expect_null(forestsearch:::.fs_decl_field(db, xi)$Mstar_shift)
  expect_false("Mstar_shift" %in% names(forestsearch:::.fs_decl_field(db, xi)))
})

# ---------------------------------------------------------------------------
# Test 4 -- scale guard
# ---------------------------------------------------------------------------
test_that("4: c0 > c2 errors naming both; a log on a ratio path errors", {
  expect_error(fs_declaration_calibration(.fit_c0, c0 = 1.1),
               "c0 = 1.1 exceeds c2 = 1")
  expect_error(fs_declaration_calibration(.fit_c0, c0 = c(0.8, 1.2)),
               "c0 <= c2")
  expect_error(fs_declaration_calibration(.fit_c0, c0 = log(0.75)),
               "natural ratio scale")
  expect_error(fs_declaration_calibration(.fit_c0, c0 = 0),
               "natural ratio scale")
  # identity path (continuous GLM, c2 = 0.25): only c0 <= c2 is checkable
  expect_error(fs_declaration_calibration(.cont_fit, c0 = 0.3),
               "exceeds c2 = 0.25")
  # at capture time, through fs_mr_inference() directly
  expect_error(
    fs_mr_inference(.cfgB0$df, .cfgB0$cands, .spec_md0,
                    selected_members = .cfgB0$cands$g1, admission = .adm_B0,
                    reselection = "maxeff", draws = 50L,
                    multiplier = "gaussian", ci_method = "ij", seed = 7L,
                    keep_declaration_field = TRUE, declaration_c0 = 0.5),
    "c0 = 0.5 exceeds c2 = 0")
  # no shifted maxima and no matrix: refuse by name, never fall back
  expect_error(fs_declaration_calibration(.fit_bare, c0 = 0.75),
               "keep_field_matrix")
  expect_error(fs_declaration_calibration(.fit_bare, c0 = 0.75),
               "declaration_c0")
  # a level with no capture column and no matrix is refused, not substituted
  expect_error(fs_declaration_calibration(.fit_bare, c0 = 0.9),
               "never substituted")
})

# ---------------------------------------------------------------------------
# Test 5 -- closed-form two-candidate check with a shift (S1.7, config. B)
# ---------------------------------------------------------------------------
.B5 <- 200000L
.c0_5 <- -0.05                                   # c2 = 0 on the MD scale
.db5 <- forestsearch:::.fs_mr_assemble(.cfgB0$df, .cfgB0$cands, .spec_md0)$B
.sig5 <- sqrt(colSums(.db5^2))
.rho5 <- sum(.db5[, 1] * .db5[, 2]) / (.sig5[1] * .sig5[2])
.dl5 <- (0 - .c0_5) / .sig5                      # delta_1, delta_2
# fw_size is taken at the rounded screen: p* = 0.90 at digits = 2 admits on
# Pcons >= 0.895 (TASK_declcal_rounding_alignment_2026-09-23)
.z_eff5 <- stats::qnorm((1 + 0.895) / 2)
.fw5_target <- 1 - .phi2ab(.z_eff5 + .dl5[1], .z_eff5 + .dl5[2], .rho5)
.kappa5_target <- stats::uniroot(
  function(k) .phi2ab(k + .dl5[1], k + .dl5[2], .rho5) - 0.95,
  c(0, 3), tol = 1e-10)$root

.check5 <- function(dc) {
  tb <- dc$c0$table
  expect_lt(abs(tb$fw_size - .fw5_target), 0.0027)
  expect_lt(abs(tb$kappa_hat - .kappa5_target), 0.025)
}

test_that("5: the shifted closed form is non-degenerate and below the unshifted", {
  expect_gte(.rho5, 0.05)
  expect_true(all(.dl5 > 0))
  expect_lt(.fw5_target, 1 - .phi2ab(.z_eff5, .z_eff5, .rho5))  # unshifted target
  expect_lt(.kappa5_target, 1.867648)
})

test_that("5: the exported path matches the shifted closed form (Gaussian)", {
  mr <- fs_mr_inference(
    .cfgB0$df, .cfgB0$cands, .spec_md0, selected_members = .cfgB0$cands$g1,
    admission = .adm_B0, reselection = "maxeff", draws = .B5,
    multiplier = "gaussian", ci_method = "ij", seed = 7L,
    keep_declaration_field = TRUE, declaration_c0 = .c0_5)
  dc <- fs_declaration_calibration(mr, c0 = .c0_5)
  expect_identical(dc$c0$source, "capture")
  expect_identical(dc$B, .B5)
  .check5(dc)
})

test_that("5: the internal helper on a hand-supplied db matches the shifted closed form", {
  set.seed(808L)
  xi <- matrix(stats::rnorm(nrow(.db5) * .B5), nrow(.db5), .B5)
  fld <- forestsearch:::.fs_decl_field(.db5, xi, keep_matrix = FALSE,
                                       shift = cbind(.dl5))
  m <- fld$Mstar_shift[, 1]
  expect_lt(abs(mean(m > .z_eff5) - .fw5_target), 0.0027)
  expect_lt(abs(stats::quantile(m, 0.95, type = 1, names = FALSE) -
                  .kappa5_target), 0.025)
})

# ---------------------------------------------------------------------------
# Test 6 -- default-off
# ---------------------------------------------------------------------------
test_that("6: unset, the fit is the baseline fit and carries no Mstar_c0", {
  expect_null(formals(forestsearch:::fs_mr_inference)$declaration_c0)
  expect_true("declaration_c0" %in% names(formals(forestsearch:::fs_mr_inference)))
  expect_null(formals(fs_declaration_calibration)$c0)
  fld0 <- .fit_noc0$mr_inference$declaration_field
  expect_false("Mstar_c0" %in% names(fld0))
  expect_false(any(c("c0", "c0_cmp") %in% names(fld0$meta)))
  expect_false("c0" %in% names(fs_declaration_calibration(.fit_noc0)))
  # the c0 capture perturbs nothing else
  on <- .fit_c0
  on$mr_inference$declaration_field$Mstar_c0 <- NULL
  on$mr_inference$declaration_field$meta$c0 <- NULL
  on$mr_inference$declaration_field$meta$c0_cmp <- NULL
  on$args_call_all$mr_inference_args$declaration_c0 <- NULL
  expect_identical(.strip_t(on), .strip_t(.fit_noc0))
})

# ---------------------------------------------------------------------------
# Test 7 -- purity
# ---------------------------------------------------------------------------
test_that("7: the fit is unchanged by fs_declaration_calibration(fit, c0 = )", {
  b1 <- .md5c(.fit_c0); b2 <- .md5c(.fit_noc0)
  invisible(fs_declaration_calibration(.fit_c0, c0 = .c0_grid))
  invisible(fs_declaration_calibration(.fit_c0, alpha = 0.10, c0 = 0.8,
                                       family = "reduced"))
  invisible(fs_declaration_calibration(.fit_noc0, c0 = .c0_grid))
  expect_identical(.md5c(.fit_c0), b1)
  expect_identical(.md5c(.fit_noc0), b2)
})

test_that("print shows one row per c0 and labels the c2 row", {
  dc <- fs_declaration_calibration(.fit_c0, c0 = .c0_grid)
  expect_output(print(dc), "Protected null level c0")
  expect_output(print(dc), "= c2, unshifted")
  expect_identical(nrow(dc$c0$table), length(.c0_grid))
})
