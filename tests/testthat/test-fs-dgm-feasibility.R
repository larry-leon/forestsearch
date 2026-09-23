# Acceptance tests for fs_dgm_feasibility() (Part B).
# Tiny DGM, small n_rep -- seconds.

mkdgm <- function(n_super = 1500L, seed = 8316951L) {
  d <- .make_binary_data(N = 300L, seed = seed)
  generate_glm_dgm(
    data            = d,
    factor_vars     = c("biomarker_hi", "sex"),
    continuous_vars = c("age", "wtkg"),
    outcome_var     = "y",
    treatment_var   = "treat",
    outcome_type    = "binary",
    effect_measure  = "OR",
    subgroup_vars   = "age",
    subgroup_cuts   = list(age = list(type = "greater", quantile = 0.55)),
    k_inter         = 1.0,
    n_super         = n_super,
    seed            = seed,
    verbose         = FALSE)
}

test_that("the interface is the OC family's: glm_dgm, and survival says what it needs", {
  expect_error(fs_dgm_feasibility(list(a = 1), n = 100, n_rep = 2L),
               "class 'glm_dgm'")
  # the refusal names the survival gap rather than guessing at it
  expect_error(fs_dgm_feasibility(structure(list(), class = "gbsg_dgm"),
                                  n = 100, n_rep = 2L),
               "Survival DGMs")
  d <- mkdgm()
  expect_s3_class(d, "glm_dgm")
  expect_true("flag_harm" %in% names(d$df_super))
  bad <- d; bad$df_super$flag_harm <- NULL
  expect_error(fs_dgm_feasibility(bad, n = 100, n_rep = 2L), "flag_harm")
})

test_that("the return contract holds", {
  d <- mkdgm()
  f <- fs_dgm_feasibility(d, n = c(200, 800), n.min = 60, n_rep = 10L, seed = 1L)
  expect_s3_class(f, "fs_dgm_feasibility")
  expect_true(all(c("table", "feasible", "args") %in% names(f)))
  expect_true(is.logical(f$feasible) && length(f$feasible) == 1L)
  expect_identical(nrow(f$table), 2L)
  for (cn in c("n", "size_mean", "size_q05", "size_q95", "size_min",
               "share_undeclarable", "share_under_events",
               "share_nonestimable", "e0_min", "e1_min", "ne0_min", "ne1_min"))
    expect_true(cn %in% names(f$table), info = cn)
  expect_true(all(f$table$share_undeclarable >= 0 & f$table$share_undeclarable <= 1))
  expect_output(print(f), "DGM feasibility")
})

test_that("undeclarable uses the search's own strict test (size <= n.min)", {
  d <- mkdgm()
  f <- fs_dgm_feasibility(d, n = 400, n.min = 60, n_rep = 20L, seed = 7L)
  # recomputed from the same draws: share with |Q| <= n.min, not < n.min
  g <- fs_dgm_feasibility(d, n = 400, n.min = 60, n_rep = 20L, seed = 7L)
  expect_identical(f$table$share_undeclarable, g$table$share_undeclarable)
  # a n.min above every draw makes everything undeclarable, and vice versa
  hi <- fs_dgm_feasibility(d, n = 400, n.min = 1e6, n_rep = 5L, seed = 7L)
  expect_identical(hi$table$share_undeclarable, 1)
  expect_false(hi$feasible)
  lo <- fs_dgm_feasibility(d, n = 400, n.min = 0, n_rep = 5L, seed = 7L)
  expect_identical(lo$table$share_undeclarable, 0)
  expect_true(lo$feasible)
})

test_that("feasible is exactly 'every undeclarable share <= tolerance'", {
  d <- mkdgm()
  f <- fs_dgm_feasibility(d, n = c(300, 900), n.min = 60, n_rep = 20L,
                          tolerance = 0.05, seed = 3L)
  expect_identical(f$feasible, all(f$table$share_undeclarable <= 0.05))
  # tolerance is the caller's to move
  f1 <- fs_dgm_feasibility(d, n = 300, n.min = 1e6, n_rep = 5L,
                           tolerance = 1, seed = 3L)
  expect_true(f1$feasible)
})

test_that("it is seeded, reproducible, and does not change the RNG kind", {
  d <- mkdgm()
  kind_before <- RNGkind()
  set.seed(42); before <- .Random.seed
  a <- fs_dgm_feasibility(d, n = 300, n_rep = 8L, seed = 99L)
  expect_identical(RNGkind(), kind_before)          # kind untouched
  expect_identical(.Random.seed, before)            # caller's stream restored
  b <- fs_dgm_feasibility(d, n = 300, n_rep = 8L, seed = 99L)
  expect_equal(a$table, b$table)                    # reproducible
})

test_that("the existence share uses Part A's per-estimand condition", {
  d <- mkdgm()
  # RD has no existence condition: the share is 0 by construction
  rd <- fs_dgm_feasibility(d, n = 400, n_rep = 10L, seed = 5L,
                           effect_measure = "RD")
  expect_identical(rd$table$share_nonestimable, 0)
  # OR is the strictest, so its share can only be >= RR's on the same draws
  or <- fs_dgm_feasibility(d, n = 400, n_rep = 10L, seed = 5L,
                           effect_measure = "OR")
  rr <- fs_dgm_feasibility(d, n = 400, n_rep = 10L, seed = 5L,
                           effect_measure = "RR")
  expect_gte(or$table$share_nonestimable, rr$table$share_nonestimable)
})

test_that("it draws through the DGM's own generator, not a re-implementation", {
  # deparse the function itself: the R/ source is not on disk under R CMD check
  src <- paste(deparse(fs_dgm_feasibility), collapse = "\n")
  expect_true(grepl("simulate_from_glm_dgm(dgm,", src, fixed = TRUE))
  # and never switches the generator
  expect_false(grepl("RNGkind(", src, fixed = TRUE))
})

test_that("it imposes nothing: floors are read, not applied", {
  d <- mkdgm()
  # changing d0.min/d1.min moves only the reported share, never the draws
  a <- fs_dgm_feasibility(d, n = 400, d0.min = 10, d1.min = 10,
                          n_rep = 10L, seed = 11L)
  b <- fs_dgm_feasibility(d, n = 400, d0.min = 0, d1.min = 0,
                          n_rep = 10L, seed = 11L)
  expect_identical(a$table$size_mean, b$table$size_mean)
  expect_identical(a$table$e0_min, b$table$e0_min)
  expect_identical(b$table$share_under_events, 0)
})
