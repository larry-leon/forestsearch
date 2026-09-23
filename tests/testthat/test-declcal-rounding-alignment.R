# =============================================================================
# fs_declaration_calibration() on the screen's rounded admission rule.
# TASK_declcal_rounding_alignment_2026-09-23: the effective threshold helper
# (.fs_pcons_eff), the settable p-star pair, digits / consistency_method
# resolution, and Gate 3 (helper vs round() on a dense grid).
# =============================================================================

.ra_field <- function(G = 5L, n = 200L, B = 4000L, seed = 20260923L) {
  set.seed(seed)
  db <- matrix(stats::rnorm(n * G) * stats::rexp(n * G), n, G)
  xi <- matrix(stats::rnorm(n * B), n, B)
  forestsearch:::.fs_decl_field(db, xi)
}

.ra_fit <- function(fld, p_star = 0.90, aca = NULL, fs_class = FALSE) {
  nm <- paste0("g", seq_len(fld$G))
  bh <- stats::setNames(seq(-0.1, 0.3, length.out = fld$G), nm)
  mr <- list(declaration_field = list(
    Mstar = fld$Mstar, Zstar = fld$Zstar, beta_hat = bh,
    sigma_D = stats::setNames(fld$sigma_D, nm), family_id = nm,
    meta = list(multiplier = "gaussian", B = fld$B, c_cons = 0,
                p_star = p_star, c_screen = NULL, log_scale = TRUE,
                column_sd = fld$column_sd, zstar_mean = fld$zstar_mean,
                field_cor = fld$field_cor,
                shared_multipliers = fld$shared_multipliers)))
  if (!fs_class) return(mr)
  structure(list(mr_inference = mr, args_call_all = aca),
            class = "forestsearch")
}

test_that("helper: effective thresholds, including a representation-error p*", {
  pe <- forestsearch:::.fs_pcons_eff
  expect_equal(pe(0.90, 2), 0.895)
  expect_equal(pe(0.90, 3), 0.8995)
  expect_equal(pe(0.9936, 2), 0.995)
  expect_equal(pe(0.9936, 3), 0.9935)
  expect_equal(pe(0.9936, 4), 0.99355)
  expect_equal(pe(0.07, 2), 0.065)          # 0.07 * 100 = 7.000000000000001
  expect_equal(pe(0.905, 2), 0.905)         # off-grid p* rounds up to 0.91
})

test_that("Gate 3: helper agrees with round() except at the exact boundary double", {
  pe <- forestsearch:::.fs_pcons_eff
  for (d in 2:4) for (p in c(0.90, 0.95, 0.99, 0.9936)) {
    thr <- pe(p, d)
    h <- 10^(-d)
    grid <- sort(unique(c(seq(thr - 2 * h, thr + 2 * h, length.out = 20001),
                          thr, thr * (1 + c(-1, 1) * 2^-52))))
    grid <- grid[grid >= 0 & grid <= 1]
    bad <- grid[(round(grid, d) >= p) != (grid >= thr)]
    # every disagreement is the boundary value itself, where round() goes
    # down on a stored double just below the half; reported, not special-cased
    expect_true(all(bad == thr), info = sprintf("d = %d, p* = %s", d, p))
    expect_lte(length(bad), 1L)
  }
})

test_that("settable p*: smallest grid value reaching kappa, and NA when none does", {
  st <- forestsearch:::.fs_decl_settable
  pe <- forestsearch:::.fs_pcons_eff
  for (d in 2:4) for (k in c(1.2, 1.645, 2.0, 2.7262897414, 2.9)) {
    s <- st(k, d)
    if (!is.finite(s$p_star)) next
    z_at <- function(p) stats::qnorm((1 + pe(p, d)) / 2)
    expect_gte(s$z_eff, k)
    expect_equal(s$z_gap, s$z_eff - k)
    expect_equal(s$pcons_eff, pe(s$p_star, d))
    expect_lt(z_at(s$p_star - 10^(-d)), k)             # minimal on the grid
  }
  expect_equal(st(2.7262897414, 2)$p_star, 1)
  expect_equal(st(2.7262897414, 4)$p_star, 0.9937)
  none <- st(5, 2)                                      # beyond p* = 1.00
  expect_true(is.na(none$p_star) && is.na(none$z_eff))
  tb <- forestsearch:::.fs_decl_settable_table(c(2.7262897414, 5), 2L)
  expect_identical(tb$pstar_achievable, c(TRUE, FALSE))
  expect_identical(tb$digits_fine[1], 4L)
  expect_lt(tb$z_gap_fine[1], 0.01)
})

test_that("fw_size is taken at the rounded threshold; kappa_hat is unchanged", {
  fld <- .ra_field()
  dc <- fs_declaration_calibration(.ra_fit(fld), alpha = 0.10, c0 = 0.8)
  z_eff <- stats::qnorm((1 + 0.895) / 2)
  expect_equal(dc$pcons_eff, 0.895)
  expect_equal(dc$z_pstar, z_eff)
  expect_identical(dc$fw_size, mean(fld$Mstar > z_eff))
  expect_gte(dc$fw_size, mean(fld$Mstar > stats::qnorm(0.95)))
  expect_identical(dc$kappa_hat,
                   stats::quantile(fld$Mstar, 0.90, type = 1, names = FALSE))
  expect_false("pstar_implied" %in% names(dc$c0$table))
  expect_true(all(c("pstar_settable", "z_eff_settable", "z_gap",
                    "digits_fine", "pstar_fine") %in% names(dc$c0$table)))
  expect_gte(dc$c0$table$z_gap, 0)
})

test_that("digits and consistency_method are read from args_call_all, with recorded fallbacks", {
  fld <- .ra_field()
  bare <- fs_declaration_calibration(.ra_fit(fld))
  expect_identical(bare$digits, 2L)
  expect_match(bare$digits_source, "default")
  expect_identical(bare$consistency_method, "resample")
  expect_match(bare$consistency_method_source, "assumed")

  d3 <- fs_declaration_calibration(.ra_fit(fld, fs_class = TRUE, aca = list(
    pconsistency.digits = 3L, consistency_method = "resample")))
  expect_identical(d3$digits, 3L)
  expect_identical(d3$digits_source, "fit$args_call_all$pconsistency.digits")
  expect_equal(d3$pcons_eff, 0.8995)

  sp <- fs_declaration_calibration(.ra_fit(fld, fs_class = TRUE, aca = list(
    pconsistency.digits = 2L, consistency_method = "split")), c0 = 0.8)
  expect_identical(sp$consistency_method, "split")
  expect_true(is.na(sp$fw_size) && is.na(sp$z_pstar))
  expect_null(sp$admitted_pstar)
  expect_true(all(is.na(sp$c0$table$fw_size)))
  expect_true(all(is.na(sp$c0$table$pstar_settable)))
  expect_identical(sp$kappa_hat, bare$kappa_hat)        # kappa is screen-free
  expect_output(print(sp), "split")
})
