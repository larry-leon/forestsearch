# ============================================================================
# test-synthetic-data-documented-behaviour.R
#
# Pins what the help pages of generate_gbsg_bootstrap_general(),
# generate_bootstrap_synthetic() and generate_bootstrap_with_noise() say
# (TASK_cran_prep_docs_legacy_2026-11-12, Part 3). Expected values come from
# the bootstrap resample, not from the functions: after set.seed(seed), the
# first random call in generate_bootstrap_synthetic() and in
# generate_bootstrap_with_noise() is sample(nrow(data), n, replace = TRUE),
# so each test rebuilds those row indices itself.
# ============================================================================

resample_idx <- function(n_data, n, seed) {
  set.seed(seed)
  sample(n_data, n, replace = TRUE)
}

# The call in generate_bootstrap_synthetic()'s first example, with n, seed and
# noise_level as arguments.
gbsg_example1_call <- function(n, seed, noise_level) {
  generate_bootstrap_synthetic(
    data = survival::gbsg,
    continuous_vars = c("age", "size", "nodes", "pgr", "er", "rfstime"),
    cat_vars = c("meno", "hormon", "status"),
    ordinal_vars = c("grade"),
    id_var = "pid",
    n = n,
    seed = seed,
    noise_level = noise_level
  )
}

# A made-up data frame: an id (named pid, the default id_var), a continuous
# column, a numeric 0/1 column and a three-level factor.
make_mixed_data <- function() {
  set.seed(2026)
  data.frame(
    pid = seq_len(2000),
    x   = rnorm(2000, mean = 50, sd = 10),
    b   = rbinom(2000, 1, 0.4),
    f   = factor(sample(c("a", "b", "c"), 2000, replace = TRUE))
  )
}

test_that("T1: generate_gbsg_bootstrap_general() keeps grade as resampled", {
  out <- suppressMessages(generate_gbsg_bootstrap_general(n = 200, seed = 11))
  idx <- resample_idx(nrow(survival::gbsg), 200, 11)
  expect_identical(out$grade, survival::gbsg$grade[idx])
})

test_that("T2: the wrapper returns the data of generate_bootstrap_synthetic()'s first example call", {
  wrapper <- suppressMessages(
    generate_gbsg_bootstrap_general(n = 200, seed = 11, noise_level = 0.15)
  )
  direct <- suppressMessages(gbsg_example1_call(n = 200, seed = 11, noise_level = 0.15))
  expect_identical(wrapper, direct)
})

test_that("T3: under auto-detection a numeric 0/1 column is continuous and stays as resampled; a factor is perturbed", {
  dat <- make_mixed_data()
  out <- generate_bootstrap_with_noise(dat, seed = 7, noise_level = 0.1)
  idx <- resample_idx(nrow(dat), nrow(dat), 7)
  # Value equality: the continuous step returns the integer 0/1 column as double.
  expect_equal(out$b, dat$b[idx])
  p_f <- mean(as.character(out$f) != as.character(dat$f[idx]))
  expect_gt(p_f, 0)
  expect_lt(p_f, 0.15)
})

test_that("T4: a numeric 0/1 column named in cat_vars, with continuous_vars given, is perturbed", {
  dat <- make_mixed_data()
  out <- generate_bootstrap_with_noise(dat, continuous_vars = "x", cat_vars = "b",
                                       seed = 7, noise_level = 0.1)
  idx <- resample_idx(nrow(dat), nrow(dat), 7)
  p_b <- mean(out$b != dat$b[idx])
  expect_gt(p_b, 0)
  expect_lt(p_b, 0.15)
})
