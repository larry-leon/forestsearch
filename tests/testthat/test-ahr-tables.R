# test-ahr-tables.R
#
# ahr_truths_table() and ahr_design_table() (TASK_wrapper_batch_2026-11-06, B2): on a
# small build the rows equal the direct ahr_limit() and ahr_design() calls. The
# comparison against the committed analytic records ran once, as a gate.

# The fixture's censoring fit (survreg) stops on Windows (win-builder, 2026-11-13).
skip_on_cran()

.at_dgm <- function(h = 3) {
  set.seed(101L); n <- 500L
  d <- data.frame(age = stats::rnorm(n, 50, 10), lab = stats::rnorm(n, 300, 100),
                  f1 = stats::rbinom(n, 1L, 0.4), treat = stats::rbinom(n, 1L, 0.5))
  d$tte <- stats::rexp(n, 0.05 * exp(-0.5 * d$treat))
  d$event <- as.integer(d$tte < 40); d$tte <- pmin(d$tte, 40)
  suppressWarnings(inject_mrct_structure(
    d, outcome_var = "tte", event_var = "event", treatment_var = "treat",
    continuous_vars = c("age", "lab"), factor_vars = "f1", x_pred = "lab",
    spline_spec = list(knot = 300, zeta = 800, log_hrs = rep(log(0.55), 3)),
    region = list(prevalence = 0.25, or_pred = 1, x1_vars = "f1", or_x1 = 3, loghr = -log(3),
                  band = list(cut = 230, or = 1)),
    x3 = list(vars = "lab", cuts = list(lab = 230), loghr = log(h / 0.55)),
    n_super = 2000, expand = "bootstrap", seed = 7L,
    continuous_vars_cens = c("age", "lab"), factor_vars_cens = "f1"))
}
.at_memo <- new.env()
.at_pops <- function() {
  if (is.null(.at_memo$p))
    .at_memo$p <- list(A = ahr_population(.at_dgm(3), 30, 6), B = ahr_population(.at_dgm(2), 30, 6))
  .at_memo$p
}

test_that("ahr_truths_table() equals the direct ahr_limit() calls", {
  p <- .at_pops()
  tt <- ahr_truths_table(p, sets = c("all", "region", "region_Gstar"), b = 0.55, h = c(A = 3, B = 2))
  expect_identical(class(tt), "data.frame")
  expect_identical(names(tt), c("cell", "set", "horizon", "marghr", "marghr_inf", "ahr", "cde"))
  expect_identical(tt$cell, rep(c("A", "B"), each = 3))
  expect_identical(tt$set, rep(c("all", "region", "region_Gstar"), 2))
  for (i in seq_len(nrow(tt))) {
    h <- if (tt$cell[i] == "A") 3 else 2
    e <- ahr_limit(p[[tt$cell[i]]], tt$set[i], b = 0.55, h = h, variance = FALSE)
    expect_identical(c(tt$marghr[i], tt$marghr_inf[i], tt$ahr[i], tt$cde[i]),
                     c(e$hr, e$marghr_inf, e$ahr, e$cde))
    expect_identical(tt$horizon[i], 30)
  }
})

test_that("sets: logical vectors, labels, and the DGM's own effects", {
  p <- .at_pops()$A
  mask <- p$region & p$band
  tt <- ahr_truths_table(p, sets = list("all", H_in_region = mask, reg = "region"), b = 0.55, h = 3)
  expect_identical(tt$set, c("all", "H_in_region", "reg"))
  expect_identical(tt$cell, rep("cell", 3))
  e <- ahr_limit(p, mask, b = 0.55, h = 3, variance = FALSE)
  expect_identical(tt$marghr[2], e$hr)
  own <- ahr_truths_table(p, sets = "region")
  expect_identical(own$marghr, ahr_limit(p, "region", variance = FALSE)$hr)
  dgm <- .at_dgm(3)
  via_dgm <- ahr_truths_table(dgm, sets = "region", analysis_time = 30, max_entry = 6)
  expect_identical(via_dgm$marghr, own$marghr)
  expect_error(ahr_truths_table(dgm), "analysis_time and max_entry")
  expect_error(ahr_truths_table(p, sets = list(mask)), "needs a name")
  expect_error(ahr_truths_table(list(p), sets = "all"), "unique names")
  expect_error(ahr_truths_table(1), "must be an ahr_population")
  expect_error(ahr_truths_table(list(A = p), sets = "all", h = c(Z = 1, Y = 2), b = 0.5), "'A' is not in it")
})

test_that("ahr_design_table() equals ahr_design() and the seven rows are the article's", {
  near <- function(a, b, tol = 1e-9) all(abs(a - b) <= tol * pmax(1, abs(b)))
  p <- .at_pops()
  dt <- ahr_design_table(p, b = 0.55, h = c(A = 3, B = 2), n = 800)
  expect_identical(class(dt), "data.frame")
  expect_identical(dt$cell, c("A", "B"))
  for (i in 1:2) {
    h <- c(3, 2)[i]
    prob <- ahr_design(p[[i]], b = 0.55, h = h, n = 800)$prob
    for (q in names(prob)) expect_identical(dt[[q]][i], prob[[q]], label = q)
    expect_identical(dt$p_Gstar_lt1_given_trigger[i], prob[["p_trigger_and_Gstar_lt1"]] / prob[["p_trigger"]])
    # transplanted: the application article's pre-table chunk
    pA <- prob[["p_itt_ub_lt1"]]; pB <- prob[["p_region_itt_gt1"]]
    pAB <- prob[["p_trigger"]];   pABC <- prob[["p_trigger_and_Gstar_lt1"]]
    rows_pre <- c(pA, pB, pA - pAB, pAB, pABC, pAB - pABC, (pA - pAB) + pABC)
    stopifnot(near(rows_pre[3] + rows_pre[4], rows_pre[1], 1e-12),
              near(rows_pre[5] + rows_pre[6], rows_pre[4], 1e-12),
              near(rows_pre[3] + rows_pre[5], rows_pre[7], 1e-12))
    expect_identical(unlist(dt[i, grep("^r[1-7]_", names(dt))], use.names = FALSE), rows_pre)
    expect_identical(c(dt$b[i], dt$h[i], dt$n[i]), c(0.55, h, 800))
  }
  expect_error(ahr_design_table(p, b = 0.55, n = 800), "b and h are required")
  expect_error(ahr_design_table(p, b = 0.55, h = 2, n = 0), "one positive number")
})
