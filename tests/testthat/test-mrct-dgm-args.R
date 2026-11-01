# Parts B and D of TASK_fs-mrct-analytic-layer_2026-11-01 (section 5, tests 4
# and 5): scale_dgm_hazard(), simulate_from_dgm(entry_delay =), and the
# latent correlation of expand_covariates() / inject_mrct_structure(). The
# builders are in helper-mrct-analytic.R.

# ---- scale_dgm_hazard() ---------------------------------------------------

test_that("scale_dgm_hazard() sets mu from the fitted value and records m", {
  skip_if_not_installed("speff2trial")
  base <- mrct_build_S2_base()
  mu0 <- base$model_params$mu; tau <- base$model_params$tau
  d3 <- scale_dgm_hazard(base, 3)
  expect_identical(d3$model_params$mu, mu0 - tau * log(3))
  expect_identical(d3$mrct$mu_fitted, mu0)
  expect_identical(d3$mrct$hazard_multiplier, 3)
  expect_identical(attr(d3, "m"), 3)
  # not cumulative
  expect_identical(scale_dgm_hazard(scale_dgm_hazard(base, 2), 3)$model_params$mu, d3$model_params$mu)
  # m = 1: unchanged except for the record
  d1 <- scale_dgm_hazard(base, 1)
  expect_identical(d1$model_params, base$model_params)
  expect_identical(d1$df_super, base$df_super)
  d1$mrct$mu_fitted <- NULL; d1$mrct$hazard_multiplier <- NULL; attr(d1, "m") <- NULL
  expect_identical(d1, base)
  # the attribute route of the mrct-consistency builders
  b2 <- base; attr(b2, "mu_fitted") <- mu0; attr(b2, "m") <- 4
  b2$model_params$mu <- mu0 - tau * log(4)
  expect_identical(scale_dgm_hazard(b2, 3)$model_params$mu, d3$model_params$mu)
  # validation
  expect_error(scale_dgm_hazard(base, 0), "positive")
  expect_error(scale_dgm_hazard(base, c(1, 2)), "positive")
  expect_error(scale_dgm_hazard(list(), 2), "generate_aft_dgm_flex")
  b3 <- base; attr(b3, "m") <- 3
  expect_error(scale_dgm_hazard(b3, 2), "fitted intercept")
})

# ---- 4. simulate_from_dgm(entry_delay = ) ---------------------------------

test_that("entry_delay = NULL is the current call; a delay shifts the delayed patients only", {
  skip_if_not_installed("speff2trial")
  dgm <- mrct_build_S2()
  args <- list(dgm = dgm, n = 1000, max_entry = 11, analysis_time = 40, cens_adjust = 50, seed = 123)
  s0 <- do.call(simulate_from_dgm, args)
  sN <- do.call(simulate_from_dgm, c(args, list(entry_delay = NULL)))
  expect_identical(sN, s0)
  a <- 6
  sD <- do.call(simulate_from_dgm, c(args, list(entry_delay = list(var = "z_region", delay = a))))
  # with cens_adjust = 50 the only censoring is administrative: c_time = 40 - entry
  entry0 <- 40 - s0$c_time; entryD <- 40 - sD$c_time
  reg <- sD$z_region == 1
  expect_true(any(reg) && any(!reg))
  expect_identical(sD$id, s0$id)
  expect_identical(sD$treat_sim, s0$treat_sim)
  expect_identical(sD$t_true, s0$t_true)
  expect_true(all(entryD[reg] >= a - 1e-9 & entryD[reg] <= a + 11 + 1e-9))
  expect_true(all(entryD[!reg] >= -1e-9 & entryD[!reg] <= 11 + 1e-9))
  expect_equal(entryD[reg], entry0[reg] + a, tolerance = 1e-9)
  expect_identical(sD$c_time[!reg], s0$c_time[!reg])
  # a numeric vector of per-patient delays, in sampling order
  sV <- do.call(simulate_from_dgm, c(args, list(entry_delay = a * as.numeric(sD$z_region == 1))))
  expect_identical(sV, sD)
  # entry_var keeps precedence
  dg2 <- dgm; dg2$df_super$entrytime <- 5
  sE <- do.call(simulate_from_dgm, c(utils::modifyList(args, list(dgm = dg2)),
                                     list(entry_var = "entrytime", entry_delay = list(var = "z_region", delay = a))))
  expect_true(all(abs(sE$c_time - 35) < 1e-9))
  # validation
  expect_error(do.call(simulate_from_dgm, c(args, list(entry_delay = list(var = "nope", delay = 1)))), "entry_delay")
  expect_error(do.call(simulate_from_dgm, c(args, list(entry_delay = 1:3))), "length n")
  expect_error(do.call(simulate_from_dgm, c(args, list(entry_delay = list(var = "z_region", delay = -1)))), "non-negative")
})

# ---- 5. the latent correlation (Part D) -----------------------------------

test_that("supplying the stored latent_R reproduces the build without a fit", {
  mrct_skip_if_no_build()
  base <- mrct_build_S2_base()
  R <- base$mrct$latent_R
  expect_true(is.matrix(R))
  local_mocked_bindings(.latentcor_available = function() stop("latentcor must not be consulted"))
  again <- mrct_build_S2_base(latent_R = R)
  expect_identical(again$df_super, base$df_super)
  expect_identical(again$model_params, base$model_params)
  expect_identical(again$mrct$latent_R, R)
  expect_identical(again$hazard_ratios, base$hazard_ratios)
  expect_length(attr(again, "build_warnings"), 0L)
  # the order of the supplied names does not matter
  p <- rev(seq_len(nrow(R)))
  again2 <- mrct_build_S2_base(latent_R = R[p, p])
  expect_identical(again2$df_super, base$df_super)
})

test_that("a supplied latent_R is validated", {
  set.seed(1)
  X <- data.frame(a = rnorm(50), b = rnorm(50), c = rbinom(50, 1, 0.4))
  R <- diag(3); dimnames(R) <- list(names(X), names(X))
  expect_silent(out <- expand_covariates(X, names(X), 20, latent_R = R))
  expect_identical(attr(out, "latent_R"), R)
  bad <- R; bad[1, 2] <- 0.5
  expect_error(expand_covariates(X, names(X), 20, latent_R = bad), "symmetric")
  bad <- R; dimnames(bad) <- list(c("a", "b", "x"), c("a", "b", "x"))
  expect_error(expand_covariates(X, names(X), 20, latent_R = bad), "names must match")
  bad <- R; bad[1, 2] <- bad[2, 1] <- 1; bad[1, 3] <- bad[3, 1] <- -1; bad[2, 3] <- bad[3, 2] <- 1
  expect_error(expand_covariates(X, names(X), 20, latent_R = bad), "positive definite")
  bad <- R; diag(bad) <- 2
  expect_error(expand_covariates(X, names(X), 20, latent_R = bad), "unit diagonal")
  expect_error(expand_covariates(X, names(X), 20, latent_R = unname(R)), "names")
})

test_that("the normal-score fallback warns once when latentcor is unavailable", {
  set.seed(1)
  X <- data.frame(a = rnorm(50), b = rnorm(50), c = rbinom(50, 1, 0.4))
  local_mocked_bindings(.latentcor_available = function() FALSE)
  msgs <- character(0)
  out <- withCallingHandlers(expand_covariates(X, names(X), 20, seed = 3),
                             warning = function(w) { msgs <<- c(msgs, conditionMessage(w)); invokeRestart("muffleWarning") })
  expect_length(msgs, 1L)
  expect_match(msgs, "latentcor package is not installed")
  expect_match(msgs, "normal-score")
  # the computation is the fallback as before: the normal-score correlation
  Zs <- sapply(X, function(v) qnorm((rank(v, ties.method = "average") - 0.5) / length(v)))
  R <- as.matrix(Matrix::nearPD(cor(Zs), corr = TRUE)$mat); dimnames(R) <- list(names(X), names(X))
  expect_identical(attr(out, "latent_R"), R)
  # and through inject_mrct_structure(): exactly one fallback warning per build
  skip_if_not_installed("speff2trial")
  b <- mrct_build_S2_base(verbose = FALSE)
  w <- attr(b, "build_warnings")
  expect_identical(sum(grepl("latent correlation not fitted", w)), 1L)
})

test_that("bootstrap expansion neither fits nor warns", {
  X <- data.frame(a = rnorm(30), b = rnorm(30))
  local_mocked_bindings(.latentcor_available = function() FALSE)
  expect_silent(expand_covariates(X, names(X), 10, method = "bootstrap"))
})
