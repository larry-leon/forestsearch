# =============================================================================
# Signed orientation of fs_oc_family_enumerate(): opposite-sign region effects
# supported; same-sign families identical to the former oriented-abs reading;
# a zero planted effect rejected; the null path untouched in structure.
# =============================================================================

# Hand-built continuous DGMs in the shape of the fs_oc_family_enumerate()
# example: covariates + potential-outcome columns + flag_harm.  sign_Qc
# controls the complement's effect direction relative to the (negative)
# planted Q effect.
make_signed_dgm <- function(sign_Qc = c("same", "opposite", "zero_Q")) {
  sign_Qc <- match.arg(sign_Qc)
  set.seed(41)
  N <- 800
  age <- round(rnorm(N, 35, 9)); pre <- round(rexp(N, 1 / 500))
  V   <- factor(rbinom(N, 1, 0.42), levels = 0:1)
  inQ <- as.integer(age > 34 & pre <= 745)
  mu0 <- 40 + 0.2 * age
  mu1 <- switch(sign_Qc,
    same     = mu0 - 26 - 40 * inQ,   # m_tau[Q], m_tau[Qc] both negative
    opposite = mu0 + 9  - 40 * inQ,   # m_tau[Q] < 0 < m_tau[Qc]
    zero_Q   = mu0 - 26 + 26 * inQ    # m_tau[Q] = 0, m_tau[Qc] = -26
  )
  structure(list(
    df_super = data.frame(age = age, preanti = pre, V = V,
                          mu0 = mu0, mu1 = mu1, flag_harm = inQ),
    outcome_type = "continuous", effect_measure = "MD",
    model_params = list(sigma = 127.5)), class = c("glm_dgm", "list"))
}

make_null_dgm <- function() {
  d <- make_signed_dgm("same")
  d$df_super$mu1 <- d$df_super$mu0 - 26      # homogeneous effect
  d$df_super$flag_harm <- 0L                 # Q empty
  d
}

signed_args <- list(confounders.name = c("age", "preanti", "V"),
                    conf.cont_jcuts = list(age = 4, preanti = 4), n.min = 60)

region_mtau <- function(dgm) {
  r <- fs_dgm_scale(dgm)$regions
  stats::setNames(r$m_tau, r$region)
}

test_that("same-sign families reproduce the oriented-abs reading exactly", {
  dgm <- make_signed_dgm("same")
  mt  <- region_mtau(dgm)
  expect_true(sign(mt[["Q"]]) == sign(mt[["Qc"]]))  # the abs-domain case

  fam <- fs_oc_family_enumerate(dgm, signed_args, n = 500)

  # the pre-edit reference, constructed from the abs formula itself
  tauQc_abs <- abs(mt[["Qc"]])
  bint_abs  <- abs(mt[["Q"]]) - abs(mt[["Qc"]])
  expect_identical(fam$beta_g, unname(tauQc_abs + bint_abs * fam$PQg))

  # orientation provenance agrees with the abs reading on this domain
  expect_identical(fam$orientation$s, sign(mt[["Q"]]))
  expect_identical(fam$orientation$tauQc, unname(tauQc_abs))
  expect_identical(fam$orientation$bint, unname(bint_abs))
  expect_identical(fam$orientation$m_tau_Q, unname(mt[["Q"]]))
  expect_identical(fam$orientation$m_tau_Qc, unname(mt[["Qc"]]))
})

test_that("opposite-sign families enumerate with the signed mixture", {
  dgm <- make_signed_dgm("opposite")
  mt  <- region_mtau(dgm)
  expect_true(sign(mt[["Q"]]) != sign(mt[["Qc"]]))

  # the former same-sign stop() must not fire
  fam <- fs_oc_family_enumerate(dgm, signed_args, n = 500)

  # beta_g is the signed mixture, computed directly from the scale table
  s <- sign(mt[["Q"]])
  beta_direct <- s * (mt[["Qc"]] + (mt[["Q"]] - mt[["Qc"]]) * fam$PQg)
  expect_equal(fam$beta_g, unname(beta_direct), tolerance = 1e-12)

  # monotone in purity (bint > 0 once oriented), and benefit-direction
  # candidates carry oriented-negative means
  expect_true(fam$orientation$bint > 0)
  expect_true(all(diff(fam$beta_g[order(fam$PQg)]) >= 0))
  expect_true(all(fam$beta_g[fam$PQg == 0] < 0))
  expect_true(fam$orientation$tauQc < 0)   # may be negative on this domain
})

test_that("a planted effect of exactly zero is rejected with guidance", {
  dgm <- make_signed_dgm("zero_Q")
  mt  <- region_mtau(dgm)
  expect_identical(unname(mt[["Q"]]), 0)
  expect_error(fs_oc_family_enumerate(dgm, signed_args, n = 500),
               "no harm direction to orient by")
})

test_that("the null path is untouched in structure", {
  fam <- fs_oc_family_enumerate(make_null_dgm(), signed_args, n = 500)
  expect_true(fam$null)
  expect_false("orientation" %in% names(fam))
  expect_identical(names(fam),
                   c("lab", "Pg", "PQg", "beta_g", "se_g", "sens_g", "spec_g",
                     "ovl", "M", "PQ", "memb", "null", "scale", "n",
                     "args_used", "cuts", "counts"))
  expect_identical(fam$PQg, rep(0, fam$M))
  expect_identical(fam$beta_g, rep(26, fam$M))  # |common effect|, oriented
})
