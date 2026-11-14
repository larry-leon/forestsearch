# =============================================================================
# compute_frontier_cis(): orientation when the fit stores no adverse_outcome
# (A3, cran-prep part 5).
#
# compute_frontier_cis() rebuilds its effect estimator from the fit's
# args_call_all.  When the fit stores adverse_outcome it is used as is.  When
# none is stored, the estimator fell back to TRUE ("higher values indicate
# harm") whatever the outcome type, while forestsearch() resolves an unset
# adverse_outcome by its own rule, outcome_type %in% c("binary", "count")
# (R/forestsearch_main.R: in Section 1A, and again in the block under the
# "Resolve adverse_outcome default" comment).  For a continuous fit at the
# default orientation the rule gives FALSE, so the CI
# table reported MD on y while the fit reports MD on -y: same magnitude,
# opposite sign.
#
# How the fallback is reached: forestsearch() resolves adverse_outcome to a
# logical before it captures its formals into args_call_all
# (mget(names(formals())), Section 1B), so a fit returned by forestsearch()
# always stores it.  The fallback is reached only by an object without the
# entry (an object saved by an older version, or one built by hand).  The
# test therefore removes the entry explicitly from a fresh fit.
#
# Expected values are computed by hand from the data (no package code):
#   - unadjusted OLS of y on a 0/1 treatment gives the treatment coefficient
#     mean(y | T = 1) - mean(y | T = 0) and the pooled-variance standard
#     error sqrt(s2 * (1/n0 + 1/n1)), s2 = (SS0 + SS1) / (n0 + n1 - 2);
#   - orientation (?forestsearch, @param adverse_outcome): FALSE, the
#     continuous default, means higher values indicate benefit; the package's
#     convention is effect > 0 = treatment harm, so under FALSE the fit's MD
#     is the coefficient for -y, i.e. -(mean(y | T = 1) - mean(y | T = 0)),
#     with the same standard error;
#   - naive CI (?compute_frontier_cis, "Naive CI": full-sample Wald CI;
#     difference measures are on the natural scale): estimate -/+ z * se,
#     z = qnorm(0.975) for ci_level = 0.95;
#   - Split CI (?compute_frontier_cis, "Split CI": a 95% interval centered at
#     the full-sample estimate): its midpoint is the estimate.
# The hand value for the selected subgroup is also tied to the fit's own
# reported estimate (max_sg_est), which compute_frontier_cis() does not touch.
#
# Tolerance 1e-10, set before the first run: the package's OLS (QR) and the
# hand difference of means / pooled SE agree to floating-point rounding
# (order 1e-15 here); an orientation error is 2 * |MD|, of order 1.
# =============================================================================

skip_on_cran()

.fco_pa  <- list(plan = "sequential", workers = 1L, show_message = FALSE)
.fco_tol <- 1e-10

# Planted region where treatment RAISES y by 1.5 (age <= 50 & biomarker_hi);
# y is then negated, so under the continuous default (higher y = better) the
# same region is where treatment harms and the fit's MD there is positive.
.fco_data <- function() {
  df <- .make_continuous_data(N = 300L, MD_harm = 1.5, seed = 42L)
  df$biomarker_hi <- as.numeric(as.character(df$biomarker_hi))
  df$sex <- as.numeric(as.character(df$sex))
  df$y <- -df$y
  df
}

# One continuous fit at the default orientation (adverse_outcome unset),
# built once for the file.
.fco_fixture <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      cache <<- suppressWarnings(suppressMessages(forestsearch(
        df.analysis = .fco_data(),
        confounders.name = c("age", "biomarker", "biomarker_hi", "sex"),
        outcome.name = "y", treat.name = "treat", id.name = "id",
        outcome_type = "continuous", effect_measure = "MD",
        effect.threshold = 0.5, consistency.threshold = 0,
        pconsistency.threshold = 0.5, fs.splits = 20, n.min = 30, maxk = 2,
        use_lasso = FALSE, use_grf = FALSE, sg_focus = "hr", is.RCT = TRUE,
        details = FALSE, quiet = TRUE, seedit = 1L,
        parallel_args = .fco_pa)))
    }
    cache
  }
})

# The fit with its stored adverse_outcome removed (see header).
.fco_unstored <- function(fit) {
  fit$args_call_all$adverse_outcome <- NULL     # deletes the list element
  fit
}

# Hand MD (see header): OLS coefficient and pooled SE, oriented.
.fco_md <- function(d, adverse) {
  y0 <- d$y[d$treat == 0]; y1 <- d$y[d$treat == 1]
  n0 <- length(y0); n1 <- length(y1)
  s2 <- (sum((y0 - mean(y0))^2) + sum((y1 - mean(y1))^2)) / (n0 + n1 - 2)
  raw <- mean(y1) - mean(y0)
  list(est = if (adverse) raw else -raw, se = sqrt(s2 * (1 / n0 + 1 / n1)))
}

# Membership of one frontier cut label, by hand.  Only the two label forms
# this fixture produces are handled: "{name <= number}" (numeric cut) and
# "{name}" (a 0/1 indicator, member when == 1), each optionally negated by a
# leading "!".  Any other form fails the premise below.
.fco_label_re <- "^!?\\{[A-Za-z_.][A-Za-z0-9_.]*( <= -?[0-9.]+)?\\}$"
.fco_mask <- function(df, label) {
  neg  <- startsWith(label, "!")
  body <- sub("^!?\\{(.*)\\}$", "\\1", label)
  m <- if (grepl(" <= ", body, fixed = TRUE)) {
    p <- strsplit(body, " <= ", fixed = TRUE)[[1]]
    df[[p[1]]] <= as.numeric(p[2])
  } else {
    df[[body]] == 1
  }
  if (neg) !m else m
}

test_that("premise: forestsearch() stores the resolved orientation", {
  fit <- .fco_fixture()
  # adverse_outcome was not passed; forestsearch() resolves it to
  # "continuous" %in% c("binary", "count") = FALSE and stores that value.
  expect_identical(fit$args_call_all$adverse_outcome, FALSE)
  expect_identical(fit$outcome_type, "continuous")
  expect_identical(fit$effect_measure, "MD")
  # The frontier is not empty, so the CI table has rows to check.
  expect_gt(nrow(fit$grp.consistency$out_sg$pareto_frontier), 0L)
  # Removing the entry leaves no element of that name, and `$` finds no
  # partial match either, so compute_frontier_cis() reaches its fallback.
  u <- .fco_unstored(fit)
  expect_false("adverse_outcome" %in% names(u$args_call_all))
  expect_null(u$args_call_all$adverse_outcome)
})

test_that("no stored orientation: selected subgroup = hand MD = max_sg_est", {
  fit <- .fco_fixture()
  ci  <- compute_frontier_cis(.fco_unstored(fit), n_splits = 20L, seed = 1L)

  # The selected subgroup is the first row of the sorted candidate table
  # (top_result = result_new[1, ]); its members are treat.recommend == 0 in
  # the fit's own df.est.
  m_sel <- as.integer(fit$grp.consistency$out_sg$result$m[1])
  expect_true(m_sel %in% ci$m)
  row <- ci[ci$m == m_sel, ]

  d   <- fit$df.est
  ref <- .fco_md(d[d$treat.recommend == 0, ], adverse = FALSE)
  z   <- stats::qnorm(0.975)

  # Hand value against the fit's own estimate (independent of the CI table).
  expect_equal(ref$est, fit$max_sg_est, tolerance = .fco_tol)
  expect_gt(ref$est, 0)                  # planted harm: MD > 0 = harm

  expect_equal(row$estimate, ref$est, tolerance = .fco_tol)
  expect_identical(sign(row$estimate), sign(fit$max_sg_est))
  expect_equal(row$naive_lcl, ref$est - z * ref$se, tolerance = .fco_tol)
  expect_equal(row$naive_ucl, ref$est + z * ref$se, tolerance = .fco_tol)
  expect_equal((row$split_lcl + row$split_ucl) / 2, ref$est,
               tolerance = .fco_tol)
})

test_that("no stored orientation: every frontier row = hand MD, oriented", {
  fit <- .fco_fixture()
  ci  <- compute_frontier_cis(.fco_unstored(fit), n_splits = 20L, seed = 1L)
  pf  <- fit$grp.consistency$out_sg$pareto_frontier
  d   <- fit$df.est
  z   <- stats::qnorm(0.975)
  m_cols <- grep("^M\\.", names(pf), value = TRUE)

  expect_identical(nrow(ci), nrow(pf))
  for (k in seq_len(nrow(pf))) {
    cuts <- unlist(pf[k, m_cols, with = FALSE], use.names = FALSE)
    cuts <- cuts[!is.na(cuts) & nzchar(cuts)]
    expect_true(all(grepl(.fco_label_re, cuts)),
                label = paste("label forms, row", k))
    keep <- Reduce(`&`, lapply(cuts, .fco_mask, df = d))
    ref  <- .fco_md(d[keep, ], adverse = FALSE)
    row  <- ci[ci$m == as.integer(pf$m[k]), ]
    lab  <- paste0("frontier m = ", pf$m[k])
    expect_equal(row$estimate, ref$est, tolerance = .fco_tol,
                 label = paste(lab, "estimate"))
    expect_identical(sign(row$estimate), sign(ref$est),
                     label = paste(lab, "sign"))
    expect_equal(row$naive_lcl, ref$est - z * ref$se, tolerance = .fco_tol,
                 label = paste(lab, "naive_lcl"))
    expect_equal(row$naive_ucl, ref$est + z * ref$se, tolerance = .fco_tol,
                 label = paste(lab, "naive_ucl"))
  }
})

test_that("no stored orientation gives the table of the stored default", {
  # forestsearch()'s rule for a continuous fit is FALSE, the value this fit
  # stores; with the same seed the two tables are the same object.
  fit <- .fco_fixture()
  # Relational check: it compares two compute_frontier_cis() outputs, not a
  # hand value.  The hand-value checks in the two tests above carry the gate
  # (independent expected values).
  expect_identical(
    compute_frontier_cis(.fco_unstored(fit), n_splits = 20L, seed = 1L),
    compute_frontier_cis(fit, n_splits = 20L, seed = 1L))
})
