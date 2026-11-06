# test-mr-method-tables.R
#
# mr_method_comparison() and mr_cost_table() (TASK_wrapper_batch_2026-11-06, B3): on a
# standalone MR result with a stand-in bootstrap, equal to the survival-post-selection
# article's code, transplanted below. The gate on a real fit and bootstrap ran once.

.mmt_mr <- function(seed = 4242L, ...) {
  set.seed(11L); n <- 300L
  df <- data.frame(time = stats::rexp(n, 0.2), event = stats::rbinom(n, 1L, 0.8),
                   treat = stats::rbinom(n, 1L, 0.5), x1 = stats::rnorm(n), x2 = stats::rnorm(n))
  cands <- list(g1 = which(df$x1 > 0), g2 = which(df$x2 > 0), g3 = which(df$x1 + df$x2 > 0.3))
  spec <- list(outcome_type = "survival", effect_measure = "HR", treat.name = "treat",
               outcome.name = "time", event.name = "event", offset.name = NULL,
               adjust_covariates = NULL, adverse_outcome = TRUE)
  fs_mr_inference(df = df, candidates = cands, spec = spec, selected_members = cands$g1,
                  admission = list(effect_floor = NULL, consistency = NULL), reselection = "maxeff",
                  draws = 200L, seed = seed, include_complement = TRUE, ci_method = "ij", ...)
}
.mmt_fit <- function(g = .mmt_mr(), ...) structure(list(mr_inference = g, ...), class = c("forestsearch", "list"))
.mmt_boot <- function(h = 1.234, hc = 0.876, B = 25L, minutes = 0.5) {
  list(H_estimates = list(H2 = h), Hc_estimates = list(H2 = hc), nb_boots = B,
       results = structure(data.frame(x = 1), timing = list(total_minutes = minutes)))
}
# ---- transplanted: vignettes/articles/survival-post-selection.qmd, three-columns
three_columns <- function(method, fit, boot) {
  g <- fit$mr_inference; cc <- g$complement
  n_H <- g$n_selected; n_Hc <- cc$n
  f2 <- function(e, l, u) sprintf("%.2f (%.2f, %.2f)", e, l, u)
  data.frame(
    method = method, region = c("H", "Hc"),
    n = c(sprintf("%d (%d%%)", n_H, round(100 * n_H / (n_H + n_Hc))),
          sprintf("%d (%d%%)", n_Hc, round(100 * n_Hc / (n_H + n_Hc)))),
    naive = c(f2(g$naive$est, g$naive$lower, g$naive$upper),
              f2(cc$naive$est, cc$naive$lower, cc$naive$upper)),
    bootstrap = sprintf("%.2f", c(boot$H_estimates$H2, boot$Hc_estimates$H2)),
    multiplier = c(f2(g$debiased$est, g$debiased$lower, g$debiased$upper),
                   f2(cc$debiased$est, cc$debiased$lower, cc$debiased$upper)))
}

test_that("mr_method_comparison() equals the article's three_columns(), cell for cell", {
  fits <- list(FS = .mmt_fit(), DINA = .mmt_fit(.mmt_mr(seed = 99L)))
  boots <- list(DINA = .mmt_boot(0.9, 0.7), FS = .mmt_boot())      # matched by name
  out <- mr_method_comparison(fits, boots)
  ref <- rbind(three_columns("FS", fits$FS, boots$FS), three_columns("DINA", fits$DINA, boots$DINA))
  expect_identical(class(out), "data.frame")
  expect_identical(out, ref)
})

test_that("est.scale and digits", {
  fit <- .mmt_fit(); g <- fit$mr_inference; boot <- .mmt_boot()
  out <- mr_method_comparison(list(FS = fit), list(FS = boot), digits = 3)
  expect_identical(out$bootstrap, c("1.234", "0.876"))
  expect_identical(out$naive[1], sprintf("%.3f (%.3f, %.3f)", g$naive$est, g$naive$lower, g$naive$upper))
  inv <- mr_method_comparison(list(FS = fit), list(FS = boot), est.scale = "1/hr", digits = 3)
  expect_identical(inv$naive[1], sprintf("%.3f (%.3f, %.3f)", 1 / g$naive$est, 1 / g$naive$upper, 1 / g$naive$lower))
  expect_identical(inv$bootstrap, sprintf("%.3f", 1 / c(1.234, 0.876)))
  # H2 stored on a fit's 1/hr scale is shown as asked
  fit_inv <- .mmt_fit(args_call_all = list(est.scale = "1/hr"))
  expect_identical(mr_method_comparison(list(FS = fit_inv), list(FS = boot), digits = 3)$bootstrap,
                   sprintf("%.3f", 1 / c(1.234, 0.876)))
  expect_identical(mr_method_comparison(list(FS = fit_inv), list(FS = boot), est.scale = "1/hr", digits = 3)$bootstrap,
                   c("1.234", "0.876"))
  expect_error(mr_method_comparison(list(FS = fit), list(FS = boot), est.scale = "rmst"))
  md <- fit; md$mr_inference$log_scale <- FALSE
  expect_error(mr_method_comparison(list(FS = md), list(FS = boot), est.scale = "1/hr"), "ratio measures only")
})

test_that("only the two-term IJ interval reaches the multiplier column", {
  fit <- .mmt_fit(.mmt_mr(ij_residual = "winner"))
  out <- mr_method_comparison(list(FS = fit), list(FS = .mmt_boot()))
  expect_match(out$multiplier[1], "\\(NA, NA\\)$")
  expect_false(grepl(sprintf("%.2f", fit$mr_inference$debiased$lower_w), out$multiplier[1], fixed = TRUE))
})

test_that("mr_cost_table() equals the article's cost chunk", {
  fits <- list(FS = .mmt_fit(), GRF = .mmt_fit(.mmt_mr(seed = 7L)))
  boots <- list(FS = .mmt_boot(B = 25L), GRF = .mmt_boot(B = 10L))
  B <- c(FS = 25L, GRF = 10L); boot_secs <- c(FS = 123.4, GRF = 56.7)
  # ---- transplanted: vignettes/articles/survival-post-selection.qmd, cost
  cost <- data.frame(
    method = names(fits),
    mr_seconds = sapply(fits, function(f) f$mr_inference$timing_seconds),
    bootstrap_B = B[names(fits)],
    bootstrap_seconds = boot_secs[names(fits)])
  cost$seconds_per_replicate <- cost$bootstrap_seconds / cost$bootstrap_B
  cost$projected_B_1000 <- round(cost$seconds_per_replicate * 1000)
  # ---- the function
  out <- mr_cost_table(fits, boots, boot_seconds = boot_secs)
  expect_identical(class(out), "data.frame")
  expect_identical(names(out), names(cost))
  for (cn in names(cost)) expect_equal(out[[cn]], unname(cost[[cn]]), label = cn, tolerance = 0)
})

test_that("mr_cost_table(): stored wall time, project_B and the checks", {
  fits <- list(FS = .mmt_fit()); boots <- list(FS = .mmt_boot(B = 20L, minutes = 0.5))
  out <- mr_cost_table(fits, boots, project_B = 500L)
  expect_identical(out$bootstrap_seconds, 30)
  expect_identical(out$seconds_per_replicate, 1.5)
  expect_identical(out$projected_B_500, 750)
  no_time <- list(FS = list(nb_boots = 2L, results = data.frame(x = 1)))
  expect_error(mr_cost_table(fits, no_time), "stores no wall time")
  expect_error(mr_cost_table(fits, boots, boot_seconds = 3), "named by method")
  expect_error(mr_cost_table(fits, boots, project_B = 0), "positive")
  expect_error(mr_method_comparison(list(.mmt_fit()), boots), "named list")
  expect_error(mr_method_comparison(fits, list(Z = .mmt_boot())), "names of fits")
  expect_error(mr_method_comparison(list(FS = list(sg.harm = "x")), boots), "no mr_inference")
  no_cc <- .mmt_fit(.mmt_mr()); no_cc$mr_inference$complement <- NULL
  expect_error(mr_method_comparison(list(FS = no_cc), boots), "no de-biased complement")
  expect_error(mr_method_comparison(fits, list(FS = list())), "H_estimates\\$H2")
})

test_that("a real fit and a two-replicate bootstrap: equal to the article's code", {
  skip_on_cran()                                     # about 10 s
  set.seed(1); n <- 500
  d <- data.frame(id = seq_len(n), treat = stats::rbinom(n, 1, 0.5), x1 = stats::rnorm(n), x2 = stats::rnorm(n))
  d$time <- stats::rexp(n, 0.1 * ifelse(d$treat == 1, ifelse(d$x1 < -0.6, 3, 0.45), 1))
  d$event <- as.integer(d$time < 30); d$time <- pmin(d$time, 30)
  fit <- forestsearch(
    df.analysis = d, outcome.name = "time", event.name = "event", treat.name = "treat", id.name = "id",
    confounders.name = c("x1", "x2"), use_grf = FALSE, use_lasso = FALSE, maxk = 1, fs.splits = 200,
    sg_focus = "hrMaxSG", mr_inference = TRUE, mr_inference_args = list(draws = 200, field_R_out = 50, field_R_in = 25),
    parallel_args = list(plan = "sequential"), quiet = TRUE, details = FALSE, plot.sg = FALSE,
    show_candidate_summary = FALSE)
  boot <- suppressWarnings(forestsearch_bootstrap_dofuture(fs.est = fit, nb_boots = 2L, seed = 8316951L,
                                                           parallel_args = list(plan = "sequential")))
  expect_identical(mr_method_comparison(list(FS = fit), list(FS = boot)), three_columns("FS", fit, boot))
  cost <- mr_cost_table(list(FS = fit), list(FS = boot), boot_seconds = c(FS = 10))
  expect_identical(cost$method, "FS")
  expect_identical(cost$bootstrap_B, 2)
  expect_identical(cost$mr_seconds, fit$mr_inference$timing_seconds)
  expect_gt(mr_cost_table(list(FS = fit), list(FS = boot))$bootstrap_seconds, 0)
})
