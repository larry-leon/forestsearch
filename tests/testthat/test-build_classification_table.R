# Tests for build_classification_table(): NPV for a replicate with no
# identified subgroup.  Leon et al. (2024), Section 3: with no subgroup the
# estimated complement is everyone, so sens(Hc) = spec = 1 and
# ppv(Hc) = npv = |Hc|/N, the share of that replicate's subjects truly in the
# complement (1 under the null).  Hand-built results; no search, no
# simulation.

.bct_dgm <- function() {
  list(hazard_ratios = list(overall = 0.8, harm_subgroup = 1.5,
                            no_harm_subgroup = 0.7),
       df_super = data.frame(flag.harm = c(1L, 0L, 0L, 0L)))
}

# Two replicates with a subgroup, two without.  The two without carry NA
# classification values, as output from an older package version does, and
# their own |Hc|/N in propHc_true, as every run_simulation_analysis() row does.
.bct_results <- function(null = FALSE) {
  data.frame(
    analysis    = "FS",
    any.H       = c(1L, 1L, 0L, 0L),
    size.H      = c(60L, 70L, 0L, 0L),
    propHc_true = if (null) rep(1, 4) else c(0.80, 0.775, 0.75, 0.70),
    sens        = if (null) rep(NA_real_, 4) else c(0.6, 0.5, NA, NA),
    spec        = c(0.85, 0.95, NA, NA),
    ppv         = if (null) c(0, 0, NA, NA) else c(0.5, 0.7, NA, NA),
    npv         = if (null) c(1, 1, NA, NA) else c(0.9, 0.8, NA, NA),
    stringsAsFactors = FALSE)
}

.bct_value <- function(results, metric, hyp = "alt") {
  sc <- list(s = list(results = results, label = "M1", n_sample = 200L,
                      dgm = .bct_dgm(), hypothesis = hyp))
  tab <- build_classification_table(sc, digits = 8)
  d <- as.data.frame(tab[["_data"]])
  d$FS[d$Metric == metric]
}

test_that("no-find, alternative: a non-detected replicate contributes npv = |Hc|/N and spec = 1", {
  r <- .bct_results()
  expect_equal(.bct_value(r, "ppv(Hc)"),  (0.9 + 0.8 + 0.75 + 0.70) / 4)
  expect_equal(.bct_value(r, "sens(Hc)"), (0.85 + 0.95 + 1 + 1) / 4)
  expect_equal(.bct_value(r, "sens(H)"),  (0.6 + 0.5 + 0 + 0) / 4)
  expect_equal(.bct_value(r, "ppv(H)"),   (0.5 + 0.7 + 0 + 0) / 4)

  one <- r[3, , drop = FALSE]
  expect_equal(.bct_value(one, "ppv(Hc)"),  0.75)
  expect_equal(.bct_value(one, "sens(Hc)"), 1)

  # values already stored (current extractors) are used as they are
  stored <- r; stored$spec[3:4] <- 1; stored$npv[3:4] <- c(0.75, 0.70)
  stored$sens[3:4] <- 0; stored$ppv[3:4] <- 0
  expect_equal(.bct_value(stored, "ppv(Hc)"), .bct_value(r, "ppv(Hc)"))
})

test_that("no-find, null: the same input gives npv = 1 and spec = 1 for a non-detected replicate", {
  # Under the null the fills do not run; the extractors store spec = 1 and
  # npv = 1 for a replicate with no subgroup.
  r <- .bct_results(null = TRUE)
  r$spec[3:4] <- 1; r$npv[3:4] <- 1
  expect_equal(.bct_value(r, "ppv(Hc)", hyp = "null"),  1)
  expect_equal(.bct_value(r, "sens(Hc)", hyp = "null"), (0.85 + 0.95 + 1 + 1) / 4)

  one <- r[3:4, , drop = FALSE]
  expect_equal(.bct_value(one, "ppv(Hc)",  hyp = "null"), 1)
  expect_equal(.bct_value(one, "sens(Hc)", hyp = "null"), 1)
})

test_that("no-find, alternative: a row without |Hc|/N keeps the fill of 1", {
  # No truth column: run_simulation_analysis() leaves propHc_true NA, so
  # |Hc|/N is not available and the row is filled with 1 (see the comment
  # in build_classification_table()).
  r <- .bct_results(); r$propHc_true[4] <- NA_real_
  expect_equal(.bct_value(r, "ppv(Hc)"), (0.9 + 0.8 + 0.75 + 1) / 4)
  r2 <- .bct_results(); r2$propHc_true <- NULL
  expect_equal(.bct_value(r2, "ppv(Hc)"), (0.9 + 0.8 + 1 + 1) / 4)
})
