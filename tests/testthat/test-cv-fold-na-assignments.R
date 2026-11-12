# ============================================================================
# test-cv-fold-na-assignments.R
#
# forestsearch_KfoldOut(): held-out fold assignments {1, NA} or {0, NA}.
#
# A held-out subject gets an NA fold assignment when a covariate that
# defines the fold's subgroup is missing for that subject.  If every
# non-NA held-out subject is on one side, treat.recommend takes two
# distinct values (one level and NA), so n_sgfound is 2, but table()
# drops NA and the fold-by-original table has a single row.  The fix takes
# the two-level branch only when that table has two rows and two columns,
# so these shapes use the existing single-level rule instead of stopping
# at tabit[2, 2].
#
# The single-level rule (forestsearch_KfoldOut(), case "nrow(tabit) == 1,
# both full-data levels present") applies the two-level definitions
#   sens_H  = n(fold 0, orig 0) / n(orig 0)
#   sens_Hc = n(fold 1, orig 1) / n(orig 1)
#   ppv_H   = n(fold 0, orig 0) / n(fold 0)
#   ppv_Hc  = n(fold 1, orig 1) / n(fold 1)
# with the missing fold level counted as zero subjects and 0/0 = 0.  All
# counts are over subjects whose fold assignment is not NA (table() leaves
# NA out).  Hand-built inputs; no search.
# ============================================================================

.sens_names_f2 <- c("sens_H", "sens_Hc", "ppv_H", "ppv_Hc")

# 100 subjects in 4 folds of 25; 0 = in H, 1 = in the complement Hc.
#   positions  1-8  : full-data H  (orig 0), fold assignment NA
#   positions  9-20 : full-data Hc (orig 1), fold assignment NA
#   positions 21-47 : full-data H  (orig 0), fold assignment `lev` (27 subjects)
#   positions 48-100: full-data Hc (orig 1), fold assignment `lev` (53 subjects)
# Full data: 35 in H, 65 in Hc.  Non-NA fold assignments: 80 subjects,
# 27 of them in full-data H and 53 in full-data Hc.
.f2_orig <- c(rep(0, 8), rep(1, 12), rep(0, 27), rep(1, 53))
.f2_rec  <- function(lev) c(rep(NA_real_, 20), rep(lev, 80))

.f2_kfoldout <- function(rec) {
  cvidx <- rep(1:4, each = 25L)
  fold_sg1 <- rep("{age <= 50}", 4L)
  res <- list(
    cv_args = list(confounders.name = c("age", "sex"), outcome.name = "time",
                   event.name = "event", treat.name = "treat"),
    sg0.name = "Not recommend", sg1.name = "Recommend", Kfolds = 4L,
    sg_analysis = "{age <= 50}",
    resCV = data.frame(
      cvindex = cvidx, sg1 = fold_sg1[cvidx], sg2 = NA_character_,
      treat.recommend = rec, treat.recommend.original = .f2_orig,
      stringsAsFactors = FALSE))
  forestsearch_KfoldOut(res, outall = FALSE)$sens_metrics_original
}

test_that("fold assignments {1, NA}: single-level agreement values, no stop", {
  # Premise: two distinct values with NA, one level after table() drops NA
  expect_identical(length(unique(.f2_rec(1))), 2L)

  s <- .f2_kfoldout(.f2_rec(1))
  expect_identical(names(s), .sens_names_f2)
  # Non-NA table: fold 1 row only, (orig 0, orig 1) = (27, 53); fold 0 row = 0.
  #   sens_H  = 0 / 27  = 0
  #   sens_Hc = 53 / 53 = 1
  #   ppv_H   = 0 / 0   -> 0 (empty fold level, 0/0 = 0)
  #   ppv_Hc  = 53 / 80 (full-data Hc share among non-NA fold assignments)
  expect_equal(unname(s), c(0, 1, 0, 53 / 80))
})

test_that("fold assignments {0, NA}: single-level agreement values, no stop", {
  expect_identical(length(unique(.f2_rec(0))), 2L)

  s <- .f2_kfoldout(.f2_rec(0))
  expect_identical(names(s), .sens_names_f2)
  # Non-NA table: fold 0 row only, (orig 0, orig 1) = (27, 53); fold 1 row = 0.
  #   sens_H  = 27 / 27 = 1
  #   sens_Hc = 0 / 53  = 0
  #   ppv_H   = 27 / 80 (full-data H share among non-NA fold assignments)
  #   ppv_Hc  = 0 / 0   -> 0 (empty fold level, 0/0 = 0)
  expect_equal(unname(s), c(1, 0, 27 / 80, 0))
})
