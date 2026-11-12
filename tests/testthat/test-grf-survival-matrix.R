# Survival GRF covariate matrix (audit finding A-059).
#
# grf.subg.harm.survival() built its matrix with
# apply(data[, confounders.name, drop = FALSE], 2, as.numeric).  When any
# confounder is a factor or character column, apply() first turns the frame
# into a character matrix: numeric columns are formatted to 7 significant
# digits and non-numeric labels become NA.  The selected cut is evaluated on
# the coded raw columns (.grf_evaluate_subgroup() / .grf_code_column()), so
# the two disagreed.  The matrix is now built through .build_grf_X(), the GLM
# path's coding.
#
# The first three blocks build no forest: fit_causal_forest() is mocked to
# capture the matrix it is handed and stop.

.grf_capture_X <- function(data, confounders.name) {
  cap <- new.env()
  testthat::local_mocked_bindings(
    fit_causal_forest = function(X, ...) {
      cap$X <- X
      stop("matrix captured", call. = FALSE)
    }
  )
  expect_error(
    grf.subg.harm.survival(
      data = data, confounders.name = confounders.name,
      outcome.name = "time", event.name = "event", id.name = "id",
      treat.name = "treat", n.min = 5
    ),
    "matrix captured"
  )
  cap$X
}

.grf_matrix_df <- function(n = 40L) {
  set.seed(11)
  data.frame(
    id    = seq_len(n),
    time  = rexp(n, 0.1),
    event = rep(c(1L, 0L, 1L, 1L), length.out = n),
    treat = rep(0:1, length.out = n),
    x     = rnorm(n) * 123.456789123,      # needs more than 7 digits
    u     = runif(n),
    k     = sample(1:9, n, replace = TRUE),
    l     = rep(c(TRUE, FALSE, FALSE), length.out = n),
    z     = factor(rep(c("0", "1"), length.out = n), levels = c("0", "1")),
    s     = factor(rep(c("m", "f", "f", "m", "f"), length.out = n)),
    ch    = rep(c("b", "a", "c"), length.out = n),
    stringsAsFactors = FALSE
  )
}

test_that("all-numeric confounders: matrix identical to the start-SHA build", {
  df <- .grf_matrix_df()
  # Start-SHA expression, verbatim (grf_main.R:226 before the fix).
  old_build <- function(data, nm) apply(data[, nm, drop = FALSE], 2, as.numeric)
  shapes <- list(
    list(df, c("x", "u")),
    list(df, c("k")),
    list(df, c("x", "k", "l")),
    list(df, c("l", "u")),
    list(df[-c(2, 5), ], c("u", "x", "k")),             # non-automatic row names
    list(`rownames<-`(df, paste0("r", seq_len(nrow(df)))), c("x", "u"))
  )
  for (sh in shapes) {
    X <- .grf_capture_X(sh[[1]], sh[[2]])
    expect_identical(X, old_build(sh[[1]], sh[[2]]),
                     info = paste(sh[[2]], collapse = ","))
  }
})

test_that("with a factor present, numeric columns keep their exact values", {
  df <- .grf_matrix_df()
  nm <- c("x", "z", "u", "s", "k", "ch")
  X  <- .grf_capture_X(df, nm)
  expect_identical(colnames(X), nm)
  expect_identical(storage.mode(X), "double")
  expect_false(anyNA(X))
  expect_identical(unname(X[, "x"]), df$x)
  expect_identical(unname(X[, "u"]), df$u)
  expect_identical(unname(X[, "k"]), as.double(df$k))
  # The start-SHA build lost digits here (as.matrix() formats the numbers).
  old_x <- suppressWarnings(apply(df[, nm, drop = FALSE], 2, as.numeric))[, "x"]
  expect_false(identical(unname(old_x), df$x))
})

test_that("factor and character columns are coded as the GLM path and the evaluator code them", {
  df <- .grf_matrix_df()
  nm <- c("x", "z", "s", "ch")
  X  <- .grf_capture_X(df, nm)
  expect_identical(X, apply(.build_grf_X(df, nm), 2, as.numeric))
  for (v in nm) {
    expect_identical(unname(X[, v]), as.double(.grf_code_column(df[[v]])),
                     info = v)
  }
  expect_identical(unname(X[, "z"]), as.double(as.character(df$z)))
  expect_identical(unname(X[, "s"]), as.double(as.integer(df$s)))
})

test_that("fitted membership agrees with the evaluated cut when factors are present", {
  skip_on_cran()
  old_cores <- Sys.getenv("_R_CHECK_LIMIT_CORES_", unset = NA)
  Sys.setenv("_R_CHECK_LIMIT_CORES_" = "TRUE")   # grf num.threads = 2
  on.exit(if (is.na(old_cores)) Sys.unsetenv("_R_CHECK_LIMIT_CORES_")
          else Sys.setenv("_R_CHECK_LIMIT_CORES_" = old_cores), add = TRUE)

  d <- .make_survival_data(N = 600L, HR_harm = 2.5, seed = 42L)
  set.seed(43)
  d$region <- factor(sample(c("EU", "US", "Asia"), nrow(d), TRUE))
  d$bm     <- rnorm(nrow(d), 50, 15)
  confs <- list(
    numeric_level_factors = c("age", "stage", "sex", "noise", "bm"),
    with_label_factor     = c("age", "stage", "sex", "noise", "region", "bm")
  )
  for (cn in names(confs)) {
    for (sel in c("tree", "frontier")) {
      r <- suppressWarnings(grf.subg.harm.survival(
        data = d, confounders.name = confs[[cn]], outcome.name = "time",
        event.name = "event", id.name = "id", treat.name = "treat",
        n.min = 60, dmin.grf = 0, maxdepth = 2, grf_selection = sel))
      info <- paste(cn, sel)
      expect_true(length(r$sg_def$conjunctions) > 0L, info = info)
      ev <- .grf_evaluate_subgroup(r$sg_def, d)
      # Tree: membership is predict(tree, X); frontier: size is counted on X.
      expect_identical(as.integer(r$data$treat.recommend), ev, info = info)
      expect_identical(as.integer(r$grf.gsub$Nsg), sum(ev == 0L), info = info)
    }
  }
})
