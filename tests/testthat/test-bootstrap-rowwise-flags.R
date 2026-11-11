# A-002: on the consistency path, forestsearch() attaches treat.recommend to
# df.est by row position and evaluates df.predict by the selected rule, so a
# bootstrap resample (ids repeated) gives Eq. 7's two terms on the right data.
# Transplanted from reproducer R01 (REPORT_code_audit_2026-11-09.md appendix).

.a002_data <- function() {
  set.seed(1); n <- 600
  d <- data.frame(id = seq_len(n), x1 = rbinom(n, 1, 0.5), x2 = rnorm(n),
                  x3 = rnorm(n), treat = rbinom(n, 1, 0.5))
  t0 <- rexp(n, 0.1 * exp(ifelse(d$x1 == 1, log(3), log(0.7)) * d$treat))
  cc <- runif(n, 5, 30)
  d$tte <- pmin(t0, cc); d$event <- as.integer(t0 <= cc)
  idx <- sample.int(n, n, replace = TRUE)
  list(d = d, boot = d[idx, ])
}

.a002_fit <- function(df.analysis, df.predict = NULL) {
  old_plan <- future::plan(future::sequential)
  on.exit(future::plan(old_plan), add = TRUE)
  suppressWarnings(suppressMessages(forestsearch(
    df.analysis = df.analysis, df.predict = df.predict,
    confounders.name = c("x1", "x2", "x3"), outcome.name = "tte",
    event.name = "event", treat.name = "treat", id.name = "id",
    use_lasso = FALSE, use_grf = FALSE, maxk = 1,
    parallel_args = list(plan = "sequential", workers = 1, show_message = FALSE),
    details = FALSE)))
}

test_that("bootstrap resample: df.est keeps resample rows, df.predict is rule-evaluated", {
  skip_on_cran()
  dat <- .a002_data(); d <- dat$d; boot <- dat$boot
  fs <- .a002_fit(boot, df.predict = d)
  expect_false(is.null(fs$sg.harm))

  # 1. resample size, not sum(k^2)
  expect_equal(nrow(fs$df.est), nrow(boot))
  # 2. every original row, flagged
  expect_equal(nrow(fs$df.predict), nrow(d))
  expect_false(anyNA(fs$df.predict$treat.recommend))
  # 3. df.predict flags are the selected rule evaluated on d
  expect_equal(fs$df.predict$treat.recommend,
               get_dfpred(d, fs$sg.harm, version = 2)$treat.recommend)
})

test_that("original data (unique ids): row-wise flags equal the former merge on id", {
  skip_on_cran()
  d <- .a002_data()$d
  fs0 <- .a002_fit(d)
  expect_false(is.null(fs0$sg.harm))

  # 4. the merge the code used before the change, computed inline
  flag <- fs0$grp.consistency$df_flag
  merged <- merge(d, flag, by.x = "id", by.y = "id", all.x = TRUE)
  expect_equal(fs0$df.est$treat.recommend[order(fs0$df.est$id)],
               merged$treat.recommend[order(merged$id)])
})
