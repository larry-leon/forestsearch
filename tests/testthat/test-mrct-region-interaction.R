# Region-dependent treatment effect (TASK_fs-dgm-region-interaction_2026-11-04):
# inject_mrct_structure(region_loghr_shift, region_harm_loghr_shift) and the
# four-cell (b, h, b_R, h_R) form of ahr_limit() / ahr_design(). The S2 builders
# are in helper-mrct-analytic.R.
#
#   psi(X) = log b + log(h/b) 1{H*} + gamma_R 1{R} + gamma_HR 1{H*} 1{R},
#   gamma_R = log(b_R / b),  gamma_HR = log(h_R / h) - log(b_R / b).

skip_on_cran()

# Fingerprint of a data.frame or list: md5 of every element's name and values
# written as hexadecimal doubles (%a), so equal fingerprints mean bit-equal values.
mrct_fingerprint <- function(x) {
  f <- tempfile()
  on.exit(unlink(f))
  writeLines(unlist(lapply(names(x), function(n)
    c(n, if (is.numeric(x[[n]])) sprintf("%a", x[[n]]) else as.character(x[[n]])))), f)
  unname(tools::md5sum(f))
}

# The S2 shifts that give region hazard ratios b_R, h_R from b = 0.55, h = 2.
s2_shifts <- function(b_R, h_R, b = 0.55, h = 2)
  list(region_loghr_shift = log(b_R / b),
       region_harm_loghr_shift = log(h_R / h) - log(b_R / b))

mrct_build_S2_shifted <- function(b_R, h_R) {
  key <- sprintf("S2_shift_%.12g_%.12g", b_R, h_R)
  if (is.null(mrct_memo[[key]]))
    mrct_memo[[key]] <- do.call(mrct_build_S2_base, s2_shifts(b_R, h_R))
  mrct_memo[[key]]
}

psi_theta_cols <- c("theta_0", "theta_1", "loghr_po", "lin_pred_0", "lin_pred_1",
                    "lin_pred_obs", "flag_harm", "z_region")

# ---- T1 backward compatibility ---------------------------------------------

test_that("T1: with the default shifts the builder output is the pre-change output", {
  mrct_skip_if_no_build()
  dgm <- mrct_build_S2_base()
  # Explicit zeros are the default call, bit for bit (every element except the
  # region model's closures, whose environments differ between calls).
  dgm0 <- mrct_build_S2_base(region_loghr_shift = 0, region_harm_loghr_shift = 0)
  expect_identical(dgm0$df_super, dgm$df_super)
  expect_identical(dgm0[setdiff(names(dgm), "mrct")], dgm[setdiff(names(dgm), "mrct")])
  expect_identical(dgm0$mrct[setdiff(names(dgm$mrct), "region_model")],
                   dgm$mrct[setdiff(names(dgm$mrct), "region_model")])
  expect_null(dgm$mrct$region_effect)
  expect_identical(names(dgm$mrct), c("seed_n", "expand", "region", "classes", "spline",
                                      "region_model", "latent_R"))

  # Against the builder before the change: fingerprints recorded at commit
  # 9eef6a64 (R 4.6.2, x86_64-pc-linux-gnu, reference LAPACK 3.12.0). The
  # build's floating point depends on the platform's latentcor fit and linear
  # algebra, so this part is not run on CRAN.
  testthat::skip_on_cran()
  ds <- dgm$df_super
  expect_identical(mrct_fingerprint(ds[psi_theta_cols]), "02532cd4c105b319b2a7d77b447671b5")
  expect_identical(mrct_fingerprint(ds), "44b789c74f98180fa3d8dd0547d798e7")
  expect_identical(mrct_fingerprint(dgm$hazard_ratios), "576e2ffd92da479aa50513d249079151")
})

# ---- T2 regional null --------------------------------------------------------

test_that("T2: b_R = h_R = 1 gives psi = 0 in the region and a null region limit", {
  mrct_skip_if_no_build()
  dgm <- mrct_build_S2_shifted(1, 1)
  base <- mrct_build_S2_base()
  ds <- dgm$df_super; reg <- ds$z_region == 1
  expect_true(any(reg & ds$flag_harm == 1) && any(reg & ds$flag_harm == 0))
  expect_lt(max(abs(ds$loghr_po[reg])), 1e-12)
  # outside the region nothing moves, bit for bit
  for (v in psi_theta_cols) expect_identical(ds[[v]][!reg], base$df_super[[v]][!reg], label = v)
  # the control arm and the prognostic structure are untouched everywhere
  expect_identical(ds$theta_0, base$df_super$theta_0)
  expect_identical(ds$lin_pred_0, base$df_super$lin_pred_0)
  expect_identical(dgm$mrct$region_effect, s2_shifts(1, 1))

  # the DGM's own effects (per-row route) on the region set
  pop <- ahr_population(scale_dgm_hazard(dgm, 3), analysis_time = 40, max_entry = 11)
  e <- ahr_limit(pop, "region", variance = FALSE)
  expect_lt(abs(e$beta), 1e-6)
  expect_equal(e$ahr, 1, tolerance = 1e-12)
  # the four-cell route on the unshifted population agrees
  e4 <- ahr_limit(mrct_pop_S2(), "region", b = 0.55, h = 2, b_R = 1, h_R = 1, variance = FALSE)
  expect_lt(abs(e4$beta), 1e-6)
})

# ---- T3 four-cell AHR identity -------------------------------------------------

test_that("T3: exp(mean psi) on a set straddling H* is b^p00 h^p10 b_R^p01 h_R^p11", {
  mrct_skip_if_no_build()
  pop <- mrct_pop_S2()
  cd4 <- mrct_build_S2()$df_super$cd40
  S <- cd4 >= 150 & cd4 <= 400            # straddles the H* boundary cd40 <= 230
  H <- pop$band; R <- pop$region
  p <- c(p00 = mean(!H[S] & !R[S]), p10 = mean(H[S] & !R[S]),
         p01 = mean(!H[S] & R[S]), p11 = mean(H[S] & R[S]))
  expect_true(all(p > 0))
  for (m in list(c(0.55, 2, 0.8, 1.5), c(0.55, 2, 1, 1), c(0.7, 1.3, 2.5, 0.4))) {
    b <- m[1]; h <- m[2]; b_R <- m[3]; h_R <- m[4]
    e <- ahr_limit(pop, S, b = b, h = h, b_R = b_R, h_R = h_R, variance = FALSE)
    expect_equal(e$ahr, b^p[["p00"]] * h^p[["p10"]] * b_R^p[["p01"]] * h_R^p[["p11"]],
                 tolerance = 1e-12)
  }
  # b_R = b, h_R = h is the two-cell model, bit for bit
  expect_identical(ahr_limit(pop, S, b = 0.55, h = 2, b_R = 0.55, h_R = 2, variance = FALSE),
                   ahr_limit(pop, S, b = 0.55, h = 2, variance = FALSE))
  expect_identical(ahr_design(pop, 0.55, 2, 1000)[c("sets", "corr", "prob")],
                   ahr_design(pop, 0.55, 2, 1000, b_R = 0.55, h_R = 2)[c("sets", "corr", "prob")])
})

# ---- T4 interaction algebra ------------------------------------------------------

test_that("T4: per-row psi of the gamma form equals the (b_R, h_R) form on a grid", {
  mrct_skip_if_no_build()
  pop <- mrct_pop_S2()
  H <- pop$band; R <- pop$region
  grid <- expand.grid(b = c(0.4, 0.55, 1), h = c(0.8, 2), b_R = c(0.3, 0.55, 1.6),
                      h_R = c(0.5, 2, 3))
  for (i in seq_len(nrow(grid))) {
    g <- grid[i, ]
    gam <- s2_shifts(g$b_R, g$h_R, b = g$b, h = g$h)
    psi_gamma <- log(g$b) + log(g$h / g$b) * H +
      .region_effect_shift(H, R, gam$region_loghr_shift, gam$region_harm_loghr_shift)
    psi_cell <- .ahr_psi(pop, g$b, g$h, NULL, g$b_R, g$h_R)
    expect_equal(psi_cell, psi_gamma, tolerance = 1e-12)
    expect_identical(unique(psi_cell[!H & !R]), log(g$b))
    expect_identical(unique(psi_cell[R & H]), if (g$b_R == g$b && g$h_R == g$h) log(g$b) + log(g$h / g$b) else log(g$h_R))
  }
  # the builder applies the gamma form to its own effects
  base <- mrct_build_S2_base()$df_super
  ds <- mrct_build_S2_shifted(0.8, 1.5)$df_super
  gam <- s2_shifts(0.8, 1.5)
  shift <- gam$region_loghr_shift * (ds$z_region == 1) +
    gam$region_harm_loghr_shift * (ds$flag_harm == 1) * (ds$z_region == 1)
  expect_equal(ds$loghr_po, base$loghr_po + shift, tolerance = 1e-12)
  tau <- mrct_build_S2_base()$model_params$tau
  expect_equal(ds$lin_pred_1, base$lin_pred_1 - tau * shift, tolerance = 1e-12)
  expect_equal(ds$lin_pred_obs, base$lin_pred_obs - tau * shift * (ds$treat == 1),
               tolerance = 1e-12)
  # region hazard ratios: b_R outside H*, h_R inside it
  reg <- ds$z_region == 1
  expect_equal(exp(unique(round(ds$loghr_po[reg & ds$flag_harm == 0], 10))), 0.8, tolerance = 1e-9)
  expect_equal(exp(unique(round(ds$loghr_po[reg & ds$flag_harm == 1], 10))), 1.5, tolerance = 1e-9)
})

test_that("the region shifts are validated", {
  mrct_skip_if_no_build()
  expect_error(mrct_build_S2_base(region_loghr_shift = NA), "one finite number")
  expect_error(mrct_build_S2_base(region_harm_loghr_shift = c(0, 1)), "one finite number")
  small <- function(...) inject_mrct_structure(
    mrct_seed_data(), outcome_var = "time_months", event_var = "cens",
    treatment_var = "treat_mrct", continuous_vars = "cd40", factor_vars = "symptom",
    x_pred = "cd40", spline_spec = list(knot = 300, zeta = 800, log_hrs = rep(log(0.55), 3)),
    n_super = 500, ...)
  expect_error(small(region_treat_loghr = 0.2, region_loghr_shift = 0.1), "mutually exclusive")
  expect_error(small(region_harm_loghr_shift = 0.1), "needs x3")
  expect_error(ahr_limit(mrct_pop_S2(), "region", b_R = 1), "need b and h")
  expect_error(ahr_limit(mrct_pop_S2(), "region", b = 0.55, h = 2, b_R = -1), "positive")
})
