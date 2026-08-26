# Tests for fs_dgm_scale() and fs_scale_se().
#
# All deterministic and fast: every expectation is a closed-form algebraic
# identity or an error path.  No Monte Carlo, so nothing here is flaky on CRAN.
#
# The generalization across outcome families was validated separately by
# Monte Carlo against 40,000 fixed-count trials (binary/RD and count/IRD agreed
# to 0.12% and 0.10%); those runs are too slow for the test suite and the
# algebraic identities below pin the same quantities exactly.

# ---- fixtures ---------------------------------------------------------------

make_continuous_dgm <- function(N = 2000, nQ = 700,
                                tau_Qc = -26.255235876,
                                b_int = -13.744764124,
                                sigma = 127.500125111) {
  in_q <- c(rep(1L, nQ), rep(0L, N - nQ))
  mu0  <- seq(40, 80, length.out = N)          # deterministic, non-constant
  tau  <- tau_Qc + b_int * in_q
  structure(
    list(
      df_super = data.frame(mu0 = mu0, mu1 = mu0 + tau, flag_harm = in_q),
      outcome_type = "continuous", effect_measure = "MD",
      model_params = list(sigma = sigma)
    ),
    class = c("glm_dgm", "list")
  )
}

make_binary_dgm <- function(N = 1000) {
  in_q <- rep(c(1L, 0L), each = N / 2)
  p0   <- seq(0.2, 0.6, length.out = N)
  p1   <- p0 + 0.05 + 0.10 * in_q
  structure(
    list(df_super = data.frame(p0 = p0, p1 = p1, flag_harm = in_q),
         outcome_type = "binary", effect_measure = "RD",
         model_params = list()),
    class = c("glm_dgm", "list")
  )
}

make_count_dgm <- function(N = 1000) {
  in_q <- rep(c(1L, 0L), each = N / 2)
  mu0  <- seq(2, 8, length.out = N)
  mu1  <- mu0 * 1.2 + 0.5 * in_q
  structure(
    list(df_super = data.frame(mu0 = mu0, mu1 = mu1, flag_harm = in_q),
         outcome_type = "count", effect_measure = "IRD",
         model_params = list()),
    class = c("glm_dgm", "list")
  )
}

var_fp <- function(a) mean((a - mean(a))^2)
cov_fp <- function(a, b) mean((a - mean(a)) * (b - mean(b)))


# ---- structure --------------------------------------------------------------

test_that("fs_dgm_scale returns the three default regions", {
  sc <- fs_dgm_scale(make_continuous_dgm())
  expect_s3_class(sc, "fs_dgm_scale")
  expect_identical(sc$regions$region, c("Q", "Qc", "S"))
  expect_identical(sc$regions$n_g, c(700L, 1300L, 2000L))
  expect_equal(sc$regions$P_g, c(0.35, 0.65, 1))
  expect_equal(sc$p_treat, 0.5)
})

test_that("region labels can be renamed", {
  sc <- fs_dgm_scale(make_continuous_dgm(), labels = c("H", "Hc", "all"))
  expect_identical(sc$regions$region, c("H", "Hc", "all"))
})


# ---- closed-form identities -------------------------------------------------

test_that("tau is constant on Q and on Qc, so V_g[tau] and C_g[mu0,tau] vanish", {
  sc <- fs_dgm_scale(make_continuous_dgm())
  q  <- sc$regions[sc$regions$region %in% c("Q", "Qc"), ]
  expect_equal(q$V_tau, c(0, 0))
  expect_equal(q$C_mu0_tau, c(0, 0))
})

test_that("V_S[tau] equals beta_int^2 * pi * (1 - pi)", {
  b_int <- -13.744764124
  dgm   <- make_continuous_dgm(b_int = b_int)
  pi_q  <- mean(dgm$df_super$flag_harm == 1)
  sc    <- fs_dgm_scale(dgm)
  expect_equal(sc$regions$V_tau[sc$regions$region == "S"],
               b_int^2 * pi_q * (1 - pi_q))
})

test_that("C_S[mu0,tau] equals beta_int * pi * (m_Q[mu0] - m_S[mu0])", {
  b_int <- -13.744764124
  dgm   <- make_continuous_dgm(b_int = b_int)
  d     <- dgm$df_super
  pi_q  <- mean(d$flag_harm == 1)
  sc    <- fs_dgm_scale(dgm)
  expect_equal(
    sc$regions$C_mu0_tau[sc$regions$region == "S"],
    b_int * pi_q * (mean(d$mu0[d$flag_harm == 1]) - mean(d$mu0))
  )
})

test_that("the continuous bracket decomposes as sigma^2 + V[mu0] + C + V[tau]/2", {
  sigma <- 127.500125111
  sc    <- fs_dgm_scale(make_continuous_dgm(sigma = sigma))
  expect_equal(
    sc$regions$bracket,
    sigma^2 + sc$regions$V_mu0 + sc$regions$C_mu0_tau + 0.5 * sc$regions$V_tau
  )
})

test_that("bracket(Q) is never below the noise floor sigma^2", {
  sigma <- 127.500125111
  sc    <- fs_dgm_scale(make_continuous_dgm(sigma = sigma))
  expect_gte(sc$regions$bracket[sc$regions$region == "Q"], sigma^2)
})

test_that("V_eff equals 2 * (V_arm0 + V_arm1) at 1:1 allocation", {
  sc <- fs_dgm_scale(make_continuous_dgm())
  expect_equal(sc$regions$V_eff, 2 * (sc$regions$V_arm0 + sc$regions$V_arm1))
})


# ---- outcome families -------------------------------------------------------

test_that("binary arm variance is m_g[p(1-p)] + V_g[p]", {
  dgm <- make_binary_dgm()
  d   <- dgm$df_super
  g   <- d$flag_harm == 1
  sc  <- fs_dgm_scale(dgm)
  r   <- sc$regions[sc$regions$region == "Q", ]
  expect_equal(r$V_arm0, mean(d$p0[g] * (1 - d$p0[g])) + var_fp(d$p0[g]))
  expect_equal(r$V_arm1, mean(d$p1[g] * (1 - d$p1[g])) + var_fp(d$p1[g]))
  expect_true(is.na(sc$sigma))
})

test_that("count arm variance is m_g[mu] + V_g[mu] (Poisson)", {
  dgm <- make_count_dgm()
  d   <- dgm$df_super
  g   <- d$flag_harm == 1
  r   <- fs_dgm_scale(dgm)$regions
  r   <- r[r$region == "Q", ]
  expect_equal(r$V_arm0, mean(d$mu0[g]) + var_fp(d$mu0[g]))
  expect_equal(r$V_arm1, mean(d$mu1[g]) + var_fp(d$mu1[g]))
})

test_that("ratio-scale effect measures are rejected, not approximated", {
  dgm <- make_binary_dgm()
  for (em in c("OR", "RR", "IRR")) {
    dgm$effect_measure <- em
    expect_error(fs_dgm_scale(dgm), "identity-scale")
  }
})


# ---- flexible regions -------------------------------------------------------

test_that("arbitrary regions accept logical vectors and integer indices alike", {
  dgm  <- make_continuous_dgm()
  N    <- nrow(dgm$df_super)
  keep <- seq_len(500)
  lg   <- logical(N); lg[keep] <- TRUE
  a <- fs_dgm_scale(dgm, regions = list(r = lg))
  b <- fs_dgm_scale(dgm, regions = list(r = keep))
  expect_equal(a$regions, b$regions)
  expect_identical(a$regions$n_g, 500L)
})

test_that("a custom region reproduces the default when it is Q", {
  dgm <- make_continuous_dgm()
  d   <- dgm$df_super
  def <- fs_dgm_scale(dgm)
  cus <- fs_dgm_scale(dgm, regions = list(Q = d$flag_harm == 1))
  expect_equal(cus$regions$bracket, def$regions$bracket[1])
})

test_that("malformed region specifications are rejected", {
  dgm <- make_continuous_dgm()
  expect_error(fs_dgm_scale(dgm, regions = list(rep(TRUE, 2000))), "named")
  expect_error(fs_dgm_scale(dgm, regions = list(z = c(TRUE, FALSE))), "length")
  expect_error(fs_dgm_scale(dgm, regions = list(z = c(1, 99999))), "indices")
  expect_error(fs_dgm_scale(dgm, regions = list(z = "a")), "logical vector")
  expect_error(fs_dgm_scale(dgm, regions = list(z = rep(FALSE, 2000))), "empty")
})

test_that("non-glm_dgm input and bad rand_ratio are rejected", {
  expect_error(fs_dgm_scale(list()), "glm_dgm")
  expect_error(fs_dgm_scale(make_continuous_dgm(), rand_ratio = 0), "positive")
  expect_error(fs_dgm_scale(make_continuous_dgm(), rand_ratio = c(1, 2)),
               "single positive")
})


# ---- allocation -------------------------------------------------------------

test_that("unequal allocation changes V_eff as V_1/p + V_0/(1-p)", {
  dgm <- make_continuous_dgm()
  sc2 <- fs_dgm_scale(dgm, rand_ratio = 2)
  expect_equal(sc2$p_treat, 2 / 3)
  expect_equal(sc2$regions$V_eff,
               sc2$regions$V_arm1 / (2 / 3) + sc2$regions$V_arm0 / (1 / 3))
})

test_that("1:1 allocation is the minimum-variance allocation here (V0 == V1 on Q)", {
  dgm <- make_continuous_dgm()      # tau constant on Q => V_arm0 == V_arm1
  v1  <- fs_dgm_scale(dgm, rand_ratio = 1)$regions$V_eff[1]
  v2  <- fs_dgm_scale(dgm, rand_ratio = 2)$regions$V_eff[1]
  expect_lt(v1, v2)
})


# ---- fs_scale_se ------------------------------------------------------------

test_that("fs_scale_se is sqrt(V_eff / (n P_g))", {
  sc <- fs_dgm_scale(make_continuous_dgm())
  n  <- 500
  expect_equal(unname(fs_scale_se(sc, n)),
               sqrt(sc$regions$V_eff / (n * sc$regions$P_g)))
  expect_named(fs_scale_se(sc, n), c("Q", "Qc", "S"))
  expect_equal(unname(fs_scale_se(sc, n, "Q")),
               sqrt(sc$regions$V_eff[1] / (n * sc$regions$P_g[1])))
})

test_that("the Jensen inflation is positive and decays with n", {
  sc <- fs_dgm_scale(make_continuous_dgm())
  r  <- vapply(c(250, 500, 1000, 2000), function(n) {
    unname(fs_scale_se(sc, n, "Q", jensen = TRUE) /
             fs_scale_se(sc, n, "Q"))
  }, numeric(1))
  expect_true(all(r > 1))
  expect_true(all(diff(r) < 0))
  expect_lt(r[4], 1.005)
})

test_that("fs_scale_se rejects unknown regions", {
  sc <- fs_dgm_scale(make_continuous_dgm())
  expect_error(fs_scale_se(sc, 500, "nope"), "unknown region")
  expect_error(fs_scale_se(list(), 500), "fs_dgm_scale")
})


# ---- print ------------------------------------------------------------------

test_that("print returns its argument invisibly", {
  sc <- fs_dgm_scale(make_continuous_dgm())
  expect_output(out <- print(sc), "sampling scale")
  expect_identical(out, sc)
})
