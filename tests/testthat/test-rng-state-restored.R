# =============================================================================
# The caller's random-number state is restored by exported forestsearch calls.
#
# subgroup.consistency() calls set.seed(seed) (seed = 8316951 unless supplied;
# forestsearch() does not supply it).  The outermost exported call that can
# reach it saves the caller's .Random.seed on entry and restores it on exit
# (.fs_preserve_rng()); nested exported calls leave the state alone, so the
# draws inside one call come from the same stream as before.
# Compute: small sequential fits (N = 400).
# =============================================================================

skip_on_cran()

# Built once, here: the helper calls set.seed() itself.
.rs_data <- .make_survival_data(N = 400L, HR_harm = 2.5, seed = 42L)

.rs_args <- function(extra = list()) {
  a <- .fs_args_for("survival", confounders = c("age", "stage", "sex", "noise"),
                    extra = c(list(use_grf = FALSE), extra))
  a[intersect(names(a), names(formals(forestsearch)))]
}

.rs_fit <- function(extra = list()) {
  suppressWarnings(suppressMessages(
    do.call(forestsearch, c(list(df.analysis = .rs_data), .rs_args(extra)))))
}

# Drop timing fields (names mentioning time/minutes/seconds) at any depth.
.rs_strip <- function(x) {
  if (is.list(x) && !is.data.frame(x)) {
    nm <- names(x)
    if (!is.null(nm)) {
      x <- x[!grepl("time|minutes|seconds|elapsed|^t\\.", nm, ignore.case = TRUE)]
    }
    for (i in seq_along(x)) {
      xi <- x[[i]]
      if (!is.null(xi) && !is.function(xi) && !is.environment(xi)) {
        x[[i]] <- .rs_strip(xi)
      }
    }
  }
  x
}

.rs_state <- function() {
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    get(".Random.seed", envir = globalenv(), inherits = FALSE)
  } else {
    NULL
  }
}

test_that("forestsearch() leaves the caller's .Random.seed as it found it", {
  set.seed(123)
  before <- .rs_state()
  fit <- .rs_fit()
  expect_false(is.null(fit$sg.harm))          # the consistency stage ran
  expect_identical(.rs_state(), before)
  # the caller's next draw is the one after set.seed(123), not one from the
  # stream subgroup.consistency() seeded
  set.seed(123); expected <- runif(1)
  set.seed(123); invisible(.rs_fit())
  expect_identical(runif(1), expected)
})

test_that("forestsearch(mr_inference = TRUE) leaves the caller's .Random.seed as it found it", {
  set.seed(2024)
  before <- .rs_state()
  fit <- .rs_fit(list(mr_inference = TRUE,
                      mr_inference_args = list(draws = 200L, ci_method = "ij")))
  expect_false(is.null(fit$mr_inference))
  expect_identical(.rs_state(), before)
})

test_that("an absent .Random.seed is absent again after forestsearch()", {
  old <- .rs_state()
  on.exit(if (!is.null(old)) assign(".Random.seed", old, envir = globalenv()),
          add = TRUE)
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    rm(".Random.seed", envir = globalenv())
  }
  invisible(.rs_fit())
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
})

test_that("a direct subgroup.consistency() call leaves the caller's .Random.seed as it found it", {
  # Capture the arguments forestsearch() passes, then call it directly.
  captured <- NULL
  orig <- get("subgroup.consistency", envir = asNamespace("forestsearch"))
  local({
    testthat::local_mocked_bindings(
      subgroup.consistency = function(...) {
        captured <<- list(...)
        orig(...)
      },
      .package = "forestsearch")
    invisible(.rs_fit())
  })
  expect_false(is.null(captured))
  set.seed(77)
  before <- .rs_state()
  res <- suppressWarnings(suppressMessages(do.call(subgroup.consistency, captured)))
  expect_identical(.rs_state(), before)
  # and it is unchanged by the caller's state (it seeds itself)
  set.seed(1)
  res2 <- suppressWarnings(suppressMessages(do.call(subgroup.consistency, captured)))
  expect_identical(.rs_strip(res2), .rs_strip(res))
})

test_that("two consecutive identical calls give identical results", {
  set.seed(9)
  a <- .rs_fit()
  b <- .rs_fit()
  expect_identical(.rs_strip(a), .rs_strip(b))
})

test_that("the restore also runs on an error exit, and later calls restore again", {
  # The mechanism itself: a function that starts like the exported ones,
  # reseeds, draws and then fails.
  f <- function() {
    forestsearch:::.fs_preserve_rng()
    set.seed(1)
    stats::runif(3)
    stop("boom")
  }
  set.seed(5)
  before <- .rs_state()
  expect_error(f(), "boom")
  expect_identical(.rs_state(), before)
  invisible(.rs_fit())
  expect_identical(.rs_state(), before)
})

test_that(".fs_rng_hold() keeps later exported calls in the same evaluation from restoring", {
  # What a future worker sees: no enclosing exported call (flag FALSE), then
  # .fs_rng_hold() as forestsearch_tenfold()'s iteration body calls it before
  # its fold loop.  The fits then continue one stream, as in the main process.
  ns <- asNamespace("forestsearch")
  st <- get(".fs_rng", envir = ns)
  old_flag <- st$active
  on.exit(st$active <- old_flag, add = TRUE)
  st$active <- FALSE
  set.seed(11)
  invisible(.rs_fit())
  restored <- .rs_state()
  set.seed(11); expect_identical(restored, .rs_state())   # outermost: restored
  forestsearch:::.fs_rng_hold()
  expect_true(st$active)
  set.seed(11); before <- .rs_state()
  invisible(.rs_fit())
  expect_false(identical(.rs_state(), before))             # held: stream continues
})
