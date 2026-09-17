# GRF membership evaluation codes covariates as the forest matrix does
# (TASK_grf_dina_fixes_2026-09-16, P1).  No forest is fitted: the structured
# subgroup definitions are constructed by hand.

.grf_cut_def <- function(...) {
  cuts <- list(...)
  cj <- data.frame(
    variable = vapply(cuts, `[[`, character(1), "variable"),
    op = vapply(cuts, `[[`, character(1), "op"),
    value = vapply(cuts, `[[`, numeric(1), "value"),
    stringsAsFactors = FALSE
  )
  list(
    conjunctions = list(cj), labels = NULL,
    definition = NA_character_, is_disjunction = FALSE
  )
}

.grf_coding_df <- function() {
  data.frame(
    z = factor(c(0, 1, 1, 0, 1, 0), levels = c("0", "1")),
    x = c(10, 20, 30, 40, 50, 60)
  )
}

test_that("a cut on a 0/1 factor covariate evaluates without NA or warning", {
  df <- .grf_coding_df()
  def <- .grf_cut_def(list(variable = "z", op = "<=", value = 0))
  expect_no_warning(tr <- .grf_evaluate_subgroup(def, df))
  expect_false(anyNA(tr))
  expect_identical(tr == 0L, c(TRUE, FALSE, FALSE, TRUE, FALSE, TRUE))

  def_gt <- .grf_cut_def(list(variable = "z", op = ">", value = 0))
  expect_no_warning(tr_gt <- .grf_evaluate_subgroup(def_gt, df))
  expect_identical(tr_gt == 0L, c(FALSE, TRUE, TRUE, FALSE, TRUE, FALSE))
})

test_that("a cut on a numeric covariate is unchanged", {
  df <- .grf_coding_df()
  def <- .grf_cut_def(list(variable = "x", op = "<=", value = 30))
  expect_no_warning(tr <- .grf_evaluate_subgroup(def, df))
  expect_identical(tr, c(0L, 0L, 0L, 1L, 1L, 1L))
})

test_that("a conjunction of a factor cut and a numeric cut is evaluated jointly", {
  df <- .grf_coding_df()
  def <- .grf_cut_def(
    list(variable = "z", op = ">", value = 0),
    list(variable = "x", op = "<=", value = 40)
  )
  expect_no_warning(tr <- .grf_evaluate_subgroup(def, df))
  expect_identical(tr == 0L, c(FALSE, TRUE, TRUE, FALSE, FALSE, FALSE))
})

test_that("the evaluator's coding is the forest matrix's coding", {
  df <- data.frame(
    z = factor(c(0, 1, 1, 0), levels = c("0", "1")),
    s = factor(c("b", "a", "b", "a")),
    ch = c("1", "0", "0", "1"),
    x = c(1.5, 2.5, 3.5, 4.5),
    stringsAsFactors = FALSE
  )
  X <- .build_grf_X(df, names(df))
  # as.matrix() stores the whole frame as double, so compare on that mode.
  for (v in names(df)) {
    expect_identical(as.double(.grf_code_column(df[[v]])), unname(X[, v]))
  }
  expect_identical(.grf_code_column(df$z), c(0, 1, 1, 0))
  expect_identical(.grf_code_column(df$s), c(2L, 1L, 2L, 1L))
  expect_identical(.grf_code_column(df$x), df$x)
})

test_that("a cut on a non-numeric-level factor uses the integer codes", {
  df <- data.frame(s = factor(c("b", "a", "b", "a")))
  def <- .grf_cut_def(list(variable = "s", op = "<=", value = 1))
  expect_no_warning(tr <- .grf_evaluate_subgroup(def, df))
  expect_identical(tr == 0L, c(FALSE, TRUE, FALSE, TRUE))
})
