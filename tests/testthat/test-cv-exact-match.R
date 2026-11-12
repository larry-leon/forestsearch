# ============================================================================
# test-cv-exact-match.R
#
# A-103: in CV_sgs(), when the full-data (analysis) subgroup has ONE factor,
# "Exact" must mean the fold selected that same subgroup -- the factor alone,
# with the fold's other slot empty (NA).  A fold whose subgroup merely
# CONTAINS the factor (a two-factor subgroup {A} & {B}) is not the same
# subgroup and scores 0 on exact_match.  The containment indicator is kept
# for one_match ("At least 1") and cov1_exact ("Cov1 exact").
#
# Folds, in order (sg1, sg2):
#   1. {A}          (A,  NA)  -- same subgroup as the analysis subgroup {A}
#   2. {A} & {B}    (A,  B )  -- contains A, but a different subgroup
#   3. {B} & {A}    (B,  A )  -- contains A (second slot), different subgroup
#   4. {B}          (B,  NA)  -- does not contain A
#   5. none         (NA, NA)  -- no subgroup found, scored 0
# ============================================================================

.A_lab <- "{age <= 50}"
.B_lab <- "{sex >= 1}"
.sg1_folds <- c(.A_lab, .A_lab, .B_lab, .B_lab, NA)
.sg2_folds <- c(NA, .B_lab, .A_lab, NA, NA)

test_that("CV_sgs one-factor subgroup: exact_match requires the same subgroup", {
  out <- CV_sgs(sg1 = .sg1_folds, sg2 = .sg2_folds,
                confs = c("age", "sex"), sg_analysis = .A_lab)

  # Expected by hand from the definition (task spec F4): exact = 1 only when
  # the fold's subgroup is {A} alone (other slot NA).  Fold 1 is {A} alone;
  # folds 2 and 3 have a second factor B; fold 4 has no A; fold 5 is
  # "none", which is scored 0.
  expect_equal(as.numeric(out$exact_match), c(1, 0, 0, 0, 0))

  # Expected by hand: "factor appears in the fold's subgroup" (either slot).
  # Folds 1, 2, 3 contain A; fold 4 does not; fold 5 is none -> 0.
  expect_equal(as.numeric(out$one_match),  c(1, 1, 1, 0, 0))
  expect_equal(as.numeric(out$cov1_exact), c(1, 1, 1, 0, 0))
})

test_that("forestsearch_KfoldOut find_metrics: Exact is the share of same-subgroup folds", {
  # Synthetic 5-fold CV result, one fold per row pattern above; 12 rows per
  # fold.  Only sg1 / sg2 per fold drive find_metrics.
  n_per <- 12L
  K <- 5L
  df <- data.frame(
    id = seq_len(K * n_per),
    tte = rep(seq(1, 2, length.out = n_per), K),
    event = rep(c(1L, 0L, 1L), length.out = K * n_per),
    treat = rep(0:1, length.out = K * n_per),
    cvindex = rep(seq_len(K), each = n_per)
  )
  df$treat.recommend.original <- rep(c(0L, 0L, 1L, 1L), length.out = nrow(df))
  df$treat.recommend <- df$treat.recommend.original
  df$sg1 <- .sg1_folds[df$cvindex]
  df$sg2 <- .sg2_folds[df$cvindex]

  res <- list(
    resCV = df, Kfolds = K, sg_analysis = .A_lab,
    sg0.name = "Not recommend", sg1.name = "Recommend",
    cv_args = list(confounders.name = c("age", "sex"), outcome.name = "tte",
                   event.name = "event", treat.name = "treat")
  )
  ko <- forestsearch_KfoldOut(res, outall = FALSE)

  # Expected by hand: find_metrics are means over folds of the per-fold
  # indicators above.  Exact = mean(1, 0, 0, 0, 0) = 1/5;
  # At least 1 = Cov1 exact = mean(1, 1, 1, 0, 0) = 3/5.
  expect_equal(unname(ko$find_metrics["Exact"]), 1 / 5)
  expect_equal(unname(ko$find_metrics["At least 1"]), 3 / 5)
  expect_equal(unname(ko$find_metrics["Cov1 exact"]), 3 / 5)
})
