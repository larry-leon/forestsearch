# =============================================================================
# vi.grf.min / max_n_confounders deprecation (0.4.0.9000): both formals are
# accepted, warned about once per call and ignored; no GRF forest is fitted in
# the forest-search path; fs_family_report() echoes the passed values on one
# inert row.
# =============================================================================

# Counts grf forest fits made anywhere inside expr, via base::trace() on the
# grf namespace (forestsearch calls grf::causal_forest / grf::causal_survival_forest
# through `::`, which resolves in that namespace).  Session-only; untraced on exit.
count_grf_fits <- function(expr) {
  cnt <- new.env(); cnt$n <- 0L
  # The tracer is evaluated in the traced function's frame, so it must carry
  # the counter environment inside the call (a closure named here is not found).
  tracer <- bquote(assign("n", get("n", envir = .(cnt)) + 1L, envir = .(cnt)))
  fns <- c("causal_forest", "causal_survival_forest")
  suppressMessages(for (f in fns) trace(f, tracer = tracer, where = asNamespace("grf"), print = FALSE))
  on.exit(suppressMessages(for (f in fns) untrace(f, where = asNamespace("grf"))), add = TRUE)
  value <- force(expr)
  list(value = value, n_fits = cnt$n)
}

.vi_df <- .make_continuous_data(N = 300L, MD_harm = 1.5, seed = 42L)

.run <- function(...) {
  suppressMessages(forestsearch(
    df.analysis = .vi_df, confounders.name = c("age", "biomarker", "biomarker_hi", "sex"),
    outcome.name = "y", event.name = "y", treat.name = "treat", id.name = "id",
    outcome_type = "continuous", effect_measure = "MD", adverse_outcome = TRUE,
    effect.threshold = 0.5, consistency.threshold = 0.0, pconsistency.threshold = 0.5,
    fs.splits = 20, n.min = 30, maxk = 2, use_lasso = FALSE, use_grf = FALSE,
    sg_focus = "eff", is.RCT = TRUE, details = FALSE,
    quiet = TRUE, seedit = 1L, parallel_args = list(plan = "sequential"), ...))
}

.membership <- function(fit) {
  if (is.null(fit$df.est)) return(NULL)
  o <- order(fit$df.est$id)
  data.frame(id = fit$df.est$id[o], treat.recommend = fit$df.est$treat.recommend[o])
}

.expect_same_fit <- function(a, b) {
  expect_identical(a$sg.harm, b$sg.harm)
  expect_identical(.membership(a), .membership(b))
  expect_identical(a$grp.consistency$out_sg$result, b$grp.consistency$out_sg$result)
  expect_identical(a$find.grps$out.found$hr.subgroups, b$find.grps$out.found$hr.subgroups)
}

test_that("deprecation warnings: non-default values warn, a default call does not", {
  expect_warning(.run(vi.grf.min = -0.2), "deprecated")
  expect_warning(.run(max_n_confounders = 5), "deprecated")
  w <- character(0)
  withCallingHandlers(.run(), warning = function(cond) {
    w <<- c(w, conditionMessage(cond)); invokeRestart("muffleWarning")
  })
  expect_false(any(grepl("deprecated", w)))
})

test_that("inert: deprecated values leave the fit unchanged", {
  base <- suppressWarnings(.run())
  .expect_same_fit(suppressWarnings(.run(vi.grf.min = -0.2)), base)
  .expect_same_fit(suppressWarnings(.run(max_n_confounders = 5)), base)
})

test_that("no forest: vi.grf.min = -0.2 fits no grf forest", {
  r <- count_grf_fits(suppressWarnings(.run(vi.grf.min = -0.2)))
  expect_identical(r$n_fits, 0L)
})

test_that("echo: fs_family_report() shows the passed value on the inert row", {
  r <- fs_family_report(suppressWarnings(.run(vi.grf.min = -0.2)))
  row <- r[r$stage == "GRF variable-importance screening (deprecated)", ]
  expect_identical(nrow(row), 1L)
  expect_identical(row$status, "inert")
  expect_match(row$values, "vi.grf.min = -0.2", fixed = TRUE)
})
