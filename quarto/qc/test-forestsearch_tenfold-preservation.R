# tests/testthat/test-forestsearch_tenfold-preservation.R
#
# Regression tests for the fold_summary / resCV_all preservation feature
# introduced in v0.3.0.  These tests lock in:
#   (1) fold_summary is always produced, with correct shape and columns.
#   (2) resCV_all is NULL by default (keep_resCV = FALSE).
#   (3) resCV_all is a list-of-data.frames of length `sims` when TRUE.
#   (4) Backward compatibility: every legacy return slot still populated.
#   (5) print.fs_tenfold surfaces the new fold-level summary.

# A small, fast survival example (GBSG-inspired but trimmed so each
# forestsearch() call takes seconds, not minutes).
.make_small_fs <- function(n = 200L, seed = 123L) {
  set.seed(seed)
  df <- data.frame(
    id       = seq_len(n),
    y        = rexp(n, rate = 0.02),
    status   = rbinom(n, 1L, 0.6),
    treat    = rbinom(n, 1L, 0.5),
    age      = rnorm(n, mean = 55, sd = 10),
    marker   = rbinom(n, 1L, 0.3),
    grade    = rbinom(n, 1L, 0.25)
  )

  forestsearch::forestsearch(
    df.analysis             = df,
    confounders.name        = c("age", "marker", "grade"),
    outcome.name            = "y",
    event.name              = "status",
    treat.name              = "treat",
    id.name                 = "id",
    is.RCT                  = TRUE,
    use_lasso               = FALSE,
    use_grf                 = FALSE,
    use_twostage            = FALSE,
    sg_focus                = "maxSG",
    fs.splits               = 50L,
    pconsistency.threshold  = 0.80,
    hr.threshold            = 1.25,
    hr.consistency          = 1.00,
    n.min                   = 30,
    d0.min                  = 5,
    d1.min                  = 5,
    maxk                    = 2L,
    max_subgroups_search    = 5L,
    parallel_args           = list(plan = "sequential",
                                   workers = 1L,
                                   show_message = FALSE),
    quiet                   = TRUE,
    details                 = FALSE
  )
}

test_that("forestsearch_tenfold returns fold_summary with correct shape", {
  skip_on_cran()                          # CV is slow
  fs <- .make_small_fs()
  skip_if(is.null(fs$sg.harm),
          "Test dataset did not identify any subgroup; rerun with new seed.")

  out <- forestsearch::forestsearch_tenfold(
    fs.est        = fs,
    sims          = 2L,
    Kfolds        = 3L,
    details       = FALSE,
    parallel_args = list(plan = "sequential",
                         workers = 1L,
                         show_message = FALSE)
  )

  # Shape
  expect_s3_class(out, "fs_tenfold")
  expect_true("fold_summary" %in% names(out))
  expect_s3_class(out$fold_summary, "data.frame")
  expect_equal(nrow(out$fold_summary), out$sims * out$Kfolds)

  # Columns
  expect_named(
    out$fold_summary,
    c("sim", "fold", "n_test", "sg1", "sg2", "any_found"),
    ignore.order = FALSE
  )

  # Types (matters for rbind-ability and downstream table ops)
  expect_type(out$fold_summary$sim,       "integer")
  expect_type(out$fold_summary$fold,      "integer")
  expect_type(out$fold_summary$n_test,    "integer")
  expect_type(out$fold_summary$sg1,       "character")
  expect_type(out$fold_summary$sg2,       "character")
  expect_type(out$fold_summary$any_found, "integer")

  # any_found is 0/1, consistent with sg1/sg2 presence
  expect_true(all(out$fold_summary$any_found %in% c(0L, 1L)))
  expect_equal(
    out$fold_summary$any_found,
    as.integer(
      !is.na(out$fold_summary$sg1) | !is.na(out$fold_summary$sg2)
    )
  )

  # n_test is positive and sums to approximately nrow(df.est) per sim
  per_sim_totals <- tapply(out$fold_summary$n_test,
                           out$fold_summary$sim, sum)
  expect_equal(as.numeric(per_sim_totals),
               rep(nrow(fs$df.est), out$sims))
})

test_that("keep_resCV = FALSE returns NULL resCV_all (default)", {
  skip_on_cran()
  fs <- .make_small_fs()
  skip_if(is.null(fs$sg.harm))

  out <- forestsearch::forestsearch_tenfold(
    fs.est        = fs,
    sims          = 2L,
    Kfolds        = 3L,
    details       = FALSE,
    parallel_args = list(plan = "sequential",
                         workers = 1L,
                         show_message = FALSE)
  )

  expect_true("resCV_all" %in% names(out))
  expect_null(out$resCV_all)
})

test_that("keep_resCV = TRUE returns list-of-data.frames with length sims", {
  skip_on_cran()
  fs <- .make_small_fs()
  skip_if(is.null(fs$sg.harm))

  out <- forestsearch::forestsearch_tenfold(
    fs.est        = fs,
    sims          = 2L,
    Kfolds        = 3L,
    details       = FALSE,
    keep_resCV    = TRUE,
    parallel_args = list(plan = "sequential",
                         workers = 1L,
                         show_message = FALSE)
  )

  expect_type(out$resCV_all, "list")
  expect_length(out$resCV_all, out$sims)

  for (resCV_i in out$resCV_all) {
    expect_s3_class(resCV_i, "data.frame")
    # Each sim's resCV should have all original subjects (one row each)
    expect_equal(nrow(resCV_i), nrow(fs$df.est))
    # Required columns present
    expect_true(all(
      c("cvindex", "sg1", "sg2",
        "treat.recommend", "treat.recommend.original") %in% names(resCV_i)
    ))
    # cvindex spans 1..Kfolds
    expect_setequal(unique(resCV_i$cvindex), seq_len(out$Kfolds))
  }
})

test_that("legacy return fields remain populated (backward compat)", {
  skip_on_cran()
  fs <- .make_small_fs()
  skip_if(is.null(fs$sg.harm))

  out <- forestsearch::forestsearch_tenfold(
    fs.est        = fs,
    sims          = 2L,
    Kfolds        = 3L,
    details       = FALSE,
    parallel_args = list(plan = "sequential",
                         workers = 1L,
                         show_message = FALSE)
  )

  # Every v0.1/v0.2 return field still present
  legacy_fields <- c("sens_summary", "find_summary",
                     "sens_out", "find_out",
                     "timing_minutes", "sims", "Kfolds")
  for (nm in legacy_fields) {
    expect_true(nm %in% names(out),
                info = paste("Missing legacy field:", nm))
  }
  expect_type(out$sens_summary, "double")
  expect_type(out$find_summary, "double")
  expect_true(is.matrix(out$sens_out))
  expect_true(is.matrix(out$find_out))
  expect_equal(nrow(out$sens_out), out$sims)
  expect_equal(nrow(out$find_out), out$sims)
})

test_that("fold_summary agreement matches find_summary Any metric", {
  skip_on_cran()
  fs <- .make_small_fs()
  skip_if(is.null(fs$sg.harm))

  out <- forestsearch::forestsearch_tenfold(
    fs.est        = fs,
    sims          = 3L,
    Kfolds        = 3L,
    details       = FALSE,
    parallel_args = list(plan = "sequential",
                         workers = 1L,
                         show_message = FALSE)
  )

  # Per-sim "any" rate computed from fold_summary
  per_sim_any <- tapply(
    out$fold_summary$any_found,
    out$fold_summary$sim,
    mean, na.rm = TRUE
  )

  # Should match find_out[, "Any"] (column naming varies -- may be
  # lower-case "any" or "Any Found"; take first column by convention
  # since find_out's first column is the "any" metric).
  find_any_col <- grep("^[Aa]ny", colnames(out$find_out), value = TRUE)[1]
  skip_if(is.na(find_any_col),
          "Could not locate Any column in find_out for cross-check.")

  expect_equal(
    as.numeric(per_sim_any),
    as.numeric(out$find_out[, find_any_col]),
    tolerance = 1e-8
  )
})

test_that("print.fs_tenfold surfaces fold-level summary lines", {
  skip_on_cran()
  fs <- .make_small_fs()
  skip_if(is.null(fs$sg.harm))

  out <- forestsearch::forestsearch_tenfold(
    fs.est        = fs,
    sims          = 2L,
    Kfolds        = 3L,
    details       = FALSE,
    parallel_args = list(plan = "sequential",
                         workers = 1L,
                         show_message = FALSE)
  )

  printed <- capture.output(print(out))
  expect_true(any(grepl("Folds with subgroup identified", printed)))
  expect_true(any(grepl("Distinct leading subgroups across folds", printed)))
  expect_true(any(grepl("Per-subject resCV:", printed)))
  # Default call should report resCV not preserved
  expect_true(any(grepl("not preserved", printed)))
})

test_that("print.fs_tenfold reports preserved resCV when keep_resCV = TRUE", {
  skip_on_cran()
  fs <- .make_small_fs()
  skip_if(is.null(fs$sg.harm))

  out <- forestsearch::forestsearch_tenfold(
    fs.est        = fs,
    sims          = 2L,
    Kfolds        = 3L,
    details       = FALSE,
    keep_resCV    = TRUE,
    parallel_args = list(plan = "sequential",
                         workers = 1L,
                         show_message = FALSE)
  )

  printed <- capture.output(print(out))
  expect_true(any(grepl("preserved for 2 simulation", printed)))
})
