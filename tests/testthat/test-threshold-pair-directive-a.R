# ============================================================================
# test-threshold-pair-directive-a.R
#
# Acceptance tests for TASK_directive_A_2026-09-18.md.
#
# c1 (effect.threshold / hr.threshold) admits a candidate subgroup to the
# family; c2 (consistency.threshold / hr.consistency) is the per-split floor it
# must then clear.  Directive A makes them a pair on the forest-search
# consistency path with a ratio estimand:
#
#   * c2 > c1 is a stop(), naming the values and the spellings the caller used;
#   * c1 supplied without c2 derives c2 = 0.80 * c1 on the RATIO scale,
#     announced once per fit and written back through the SECTION 2B-ii sync so
#     every bootstrap replicate and CV fold resolves the same value;
#   * both spellings of one threshold at disagreeing values is an error.
#
# Scope is narrow on purpose: subgroup_method = "consistency" with survival
# (HR) or binary "OR".  dina, grf, RD, IRD, MD and IRR are inert.
#
# Most of this file executes no fit: it drives the source-lifted resolver in
# helper-threshold-sync.R over the argument-list matrix, and calls
# .fs_resolve_threshold_pair() directly.  The two micro-fits at the end are the
# only compute, and they are the smallest that exercise the real entry point.
# ============================================================================

.DA_IN_SCOPE <- c("survival", "binary-unset", "binary-OR")

.da_probe <- function(validate) {
  probe_threshold_sync(sync = TRUE, validate = validate, extended = TRUE)
}

.da_scope_rows <- function(p) {
  p$estimand %in% .DA_IN_SCOPE & p$method == "consistency"
}


# ---------------------------------------------------------------------------
# The matrix
# ---------------------------------------------------------------------------

test_that("the probe matrix covers the Directive A cell families", {
  p <- .da_probe(TRUE)

  expect_identical(nrow(p), 175L)
  # The 63 sync cells are still there, unchanged in identity.
  expect_identical(
    sum(p$method == "consistency" & p$values %in% c("default", "custom", "-") &
          p$subset != "c2"),
    63L)
  expect_setequal(unique(p$values),
                  c("-", "default", "custom", "c1derive", "c2high", "c2gt",
                    "disagree"))
  expect_setequal(unique(p$method), c("consistency", "dina", "grf"))
})


test_that("nothing outside the consistency path on a ratio estimand moves", {
  # The comparison is re-measured, not read from a frozen file: the same probe
  # with the SECTION 1A3 statements excluded IS the pre-Directive-A tree, and
  # it reproduces dev/reports/baseline_directive_A_2026-09-18.csv (asserted
  # below).  Every cell outside the scope must be byte-identical between the
  # two runs.
  b <- .da_probe(FALSE)
  p <- .da_probe(TRUE)

  out <- !.da_scope_rows(p)
  expect_identical(sum(out), 112L)      # 175 cells, 63 of them in scope
  expect_identical(b[out, ], p[out, ])

  # Named explicitly, because these are the dispositions that govern:
  for (est in c("binary-RD", "continuous-MD", "count-IRD", "count-IRR")) {
    k <- p$estimand == est
    expect_identical(b[k, ], p[k, ], info = est)
  }
  for (m in c("dina", "grf")) {
    k <- p$method == m
    expect_identical(sum(k), 14L)
    expect_identical(b[k, ], p[k, ], info = m)
  }
})


test_that("in scope, exactly the derivation and error cells move", {
  b <- .da_probe(FALSE)
  p <- .da_probe(TRUE)
  k <- .da_scope_rows(p)
  moved <- vapply(seq_len(nrow(p)),
                  function(i) !identical(b[i, ], p[i, ]), logical(1))

  expect_identical(sum(k), 63L)
  # Of the 63 in-scope cells, 36 move and 27 are untouched (the pair rule is
  # silent wherever both thresholds are supplied and c2 <= c1).
  expect_identical(sum(moved), 36L)
  expect_true(all(k[moved]))

  # 18 derive, 18 error -- and nothing else in the whole matrix errors or
  # announces anything.
  expect_identical(sum(!is.na(p$parent_error)), 18L)
  expect_identical(sum(p$n_parent_messages > 0L), 18L)
  expect_true(all(k[!is.na(p$parent_error)]))
  expect_true(all(k[p$n_parent_messages > 0L]))

  # Errors: c2 above c1 (supplied, or the default c1), and the two spellings
  # disagreeing.  Derivations: every cell that supplies c1 alone.
  err <- table(p$values[!is.na(p$parent_error)])
  expect_identical(as.integer(err[c("c2gt", "c2high", "disagree")]),
                   c(6L, 6L, 6L))
  expect_identical(length(err), 3L)
  der <- table(p$values[p$n_parent_messages > 0L])
  expect_identical(as.integer(der[c("default", "custom", "c1derive")]),
                   c(6L, 6L, 6L))
  expect_identical(length(der), 3L)

  # 6 of the 18 derivations sit at the rule's fixed point (c1 = 1.25 -> c2 =
  # 1.0): the resolution is unmoved there and only the announcement is new.
  fixed <- p$n_parent_messages > 0L & p$values == "default"
  cols <- c("parent_screening", "parent_consistency", "parent_scr_natural",
            "parent_con_natural", "parent_scale")
  expect_identical(b[fixed, cols], p[fixed, cols])
})


test_that("a replicate resolves the derived c2, and never re-announces it", {
  p <- .da_probe(TRUE)

  # Gate 6c: replicate == parent in every cell of the matrix, derived ones
  # included.
  expect_identical(sum(p$violation), 0L)

  der <- p$n_parent_messages > 0L
  expect_identical(p$parent_con_natural[der], p$rep_con_natural[der])
  expect_identical(p$acall_consistency.threshold[der],
                   p$parent_con_natural[der])

  # The announcement fires once, in the parent frame, and cannot fire per
  # replicate: the replay supplies consistency.threshold explicitly, so
  # user_set_consistency is TRUE there and the derivation is not re-entered.
  expect_true(all(p$n_parent_messages[der] == 1L))
  expect_identical(sum(p$n_rep_messages), 0L)
})


test_that("the probe's validate = FALSE run is the committed baseline", {
  f <- "../../dev/reports/baseline_directive_A_2026-09-18.csv"
  skip_if_not(file.exists(f), "baseline CSV not in this tree")
  b <- utils::read.csv(f, colClasses = "character")
  got <- .da_probe(FALSE)
  got[] <- lapply(got, as.character)
  expect_identical(got, b)
})


# ---------------------------------------------------------------------------
# The rule itself
# ---------------------------------------------------------------------------

.da_pair <- function(c1, c2, c1_new = NULL, c2_new = NULL,
                     c1_legacy = NULL, c2_legacy = NULL,
                     outcome_type = "survival", effect_measure = NULL,
                     subgroup_method = "consistency",
                     user_set_threshold = TRUE, user_set_consistency = TRUE,
                     quiet = FALSE) {
  forestsearch:::.fs_resolve_threshold_pair(
    outcome_type = outcome_type, effect_measure = effect_measure,
    subgroup_method = subgroup_method, c1 = c1, c2 = c2,
    c1_new = c1_new, c2_new = c2_new,
    c1_legacy = c1_legacy, c2_legacy = c2_legacy,
    user_set_threshold = user_set_threshold,
    user_set_consistency = user_set_consistency, quiet = quiet)
}


test_that("c2 > c1 stops, naming the values and the spellings used", {
  # New spellings, the task document's worked example verbatim.
  expect_error(
    .da_pair(c1 = 0.90, c2 = 1.00, c1_new = 0.90, c2_new = 1.00),
    "c2 > c1 not allowed for FS: consistency.threshold = 1.00 exceeds effect.threshold = 0.90",
    fixed = TRUE)

  # Legacy spellings: the message follows the caller.
  expect_error(
    .da_pair(c1 = 0.90, c2 = 1.00, c1_legacy = 0.90, c2_legacy = 1.00),
    "c2 > c1 not allowed for FS: hr.consistency = 1.00 exceeds hr.threshold = 0.90",
    fixed = TRUE)

  # Mixed, and against a c1 left at its default (the caller wrote no c1, so
  # the preferred spelling is reported).
  expect_error(
    .da_pair(c1 = 1.25, c2 = 1.50, c2_legacy = 1.50,
             user_set_threshold = FALSE),
    "c2 > c1 not allowed for FS: hr.consistency = 1.50 exceeds effect.threshold = 1.25",
    fixed = TRUE)

  # Binary OR is in scope on the same terms.
  expect_error(
    .da_pair(c1 = 0.90, c2 = 1.00, c1_new = 0.90, c2_new = 1.00,
             outcome_type = "binary", effect_measure = "OR"),
    "c2 > c1 not allowed for FS", fixed = TRUE)
})


test_that("c2 = c1 passes, and c2 < c1 passes", {
  expect_silent(expect_identical(.da_pair(c1 = 1.25, c2 = 1.25), 1.25))
  expect_silent(expect_identical(.da_pair(c1 = 1.25, c2 = 1.00), 1.00))
})


test_that("a silent c2 derives 0.80 * c1 on the ratio scale", {
  # The task's worked values.  Multiplicative on the natural scale, which is
  # an ADDITIVE shift of log(0.80) once forestsearch() takes logs -- never
  # 0.80 * log(c1).
  for (v in list(c(1.25, 1.00), c(1.00, 0.80), c(0.90, 0.72))) {
    got <- suppressMessages(
      .da_pair(c1 = v[1], c2 = 1.0, c1_new = v[1],
               user_set_consistency = FALSE))
    expect_equal(got, v[2], tolerance = 1e-12)
    # The log-scale statement of the same rule.
    expect_equal(log(got), log(v[1]) + log(0.80), tolerance = 1e-12)
  }
})


test_that("the derivation is announced once, naming the value and its c1", {
  expect_message(
    .da_pair(c1 = 0.90, c2 = 1.0, c1_new = 0.90, user_set_consistency = FALSE),
    "derived as 0.80 * effect.threshold = 0.80 * 0.9 = 0.72", fixed = TRUE)
  expect_message(
    .da_pair(c1 = 0.90, c2 = 1.0, c1_legacy = 0.90,
             user_set_consistency = FALSE),
    "derived as 0.80 * hr.threshold = 0.80 * 0.9 = 0.72", fixed = TRUE)

  # Exactly one message, and quiet = TRUE suppresses it without changing the
  # value.
  m <- testthat::capture_messages(
    .da_pair(c1 = 1.50, c2 = 1.0, c1_new = 1.50, user_set_consistency = FALSE))
  expect_identical(length(m), 1L)
  expect_silent(
    got <- .da_pair(c1 = 1.50, c2 = 1.0, c1_new = 1.50,
                    user_set_consistency = FALSE, quiet = TRUE))
  expect_equal(got, 1.20, tolerance = 1e-12)
})


test_that("an explicitly supplied c2 is never overridden", {
  for (sp in c("new", "legacy")) {
    got <- expect_silent(.da_pair(
      c1 = 1.50, c2 = 1.10,
      c1_new = if (sp == "new") 1.50, c2_new = if (sp == "new") 1.10,
      c1_legacy = if (sp == "legacy") 1.50,
      c2_legacy = if (sp == "legacy") 1.10))
    expect_identical(got, 1.10)
  }
  # Supplied at the value the derivation would have produced: still not a
  # derivation, so nothing is announced.
  expect_silent(.da_pair(c1 = 1.50, c2 = 1.20, c1_new = 1.50, c2_new = 1.20))
})


test_that("both spellings of one threshold, disagreeing, is an error", {
  expect_error(
    .da_pair(c1 = 1.50, c2 = 1.00, c1_new = 1.50, c1_legacy = 1.25),
    paste0("effect.threshold = 1.5 and hr.threshold = 1.25 are two spellings ",
           "of the same threshold (c1) and disagree; supply one."),
    fixed = TRUE)
  expect_error(
    .da_pair(c1 = 1.25, c2 = 1.20, c2_new = 1.20, c2_legacy = 1.00),
    paste0("consistency.threshold = 1.2 and hr.consistency = 1 are two ",
           "spellings of the same threshold (c2) and disagree; supply one."),
    fixed = TRUE)
  # Agreeing is not an error.
  expect_silent(.da_pair(c1 = 1.25, c2 = 1.00, c1_new = 1.25, c1_legacy = 1.25,
                         c2_new = 1.00, c2_legacy = 1.00))
})


test_that("dina, grf and the non-OR estimands are inert", {
  degenerate <- list(c1 = 0.90, c2 = 1.00, c1_new = 0.90, c2_new = 1.00)

  # No error, no derivation, c2 returned untouched.
  for (m in c("dina", "grf")) {
    expect_identical(
      do.call(.da_pair, c(degenerate, list(subgroup_method = m))), 1.00)
    expect_identical(
      do.call(.da_pair, c(degenerate,
                          list(subgroup_method = m,
                               user_set_consistency = FALSE))), 1.00)
  }
  for (em in c("RD", "RR", "IRD")) {
    expect_identical(
      do.call(.da_pair, c(degenerate,
                          list(outcome_type = "binary", effect_measure = em))),
      1.00)
  }
  for (ot in c("continuous", "count")) {
    em <- if (ot == "continuous") "MD" else "IRR"
    expect_identical(
      do.call(.da_pair, c(degenerate,
                          list(outcome_type = ot, effect_measure = em,
                               user_set_consistency = FALSE))), 1.00)
  }
})


test_that("a non-scalar or NA threshold reaches the branches that handle it", {
  # isTRUE() on the comparison, not a bare `if`: these must pass through here
  # and fail (or not) exactly where they failed before.
  expect_identical(.da_pair(c1 = NA_real_, c2 = 1.00), 1.00)
  expect_identical(.da_pair(c1 = 1.25, c2 = NA_real_), NA_real_)
})


# ---------------------------------------------------------------------------
# Step 4 -- the deletion gate, as source facts
# ---------------------------------------------------------------------------

test_that("the consistency-stage skip branch stays reachable by legal calls", {
  # Route 1: an empty candidate family.  format_search_results() returns
  # out.found = NULL when nothing passed the search, so has_subgroups is FALSE
  # whatever the thresholds are.
  empty <- forestsearch:::format_search_results(
    results_list = list(), Z = matrix(0, 0, 0), details = FALSE,
    t.sofar = 0, L = 0L, max_count = 0L)
  expect_null(empty$out.found)

  # Route 2: sg_focus = "maxeff" on the consistency engine disables the search
  # effect floor, so candidate effects are NOT bounded below by c1 and
  # any(hr_values > c2) can be FALSE with c2 <= c1.
  ap <- forestsearch:::.fs_admission_applies("maxeff", "consistency")
  expect_false(unname(ap[["effect"]]))
  expect_true(unname(
    forestsearch:::.fs_admission_applies("hr", "consistency")[["effect"]]))

  # So the branch is not unreachable and was not deleted.
  expect_true(any(grepl("has_subgroups <- any(hr_values > check_threshold",
                        deparse(body(forestsearch)), fixed = TRUE)))
})


# ---------------------------------------------------------------------------
# The two micro-fits
# ---------------------------------------------------------------------------

test_that("a binary-OR fit with c1 = 0.90 and no c2 runs the stage at 0.72", {
  df <- .make_binary_data(N = 150L, seed = 7L)
  args <- .fs_args_for("binary", confounders = c("age", "biomarker_hi", "sex"),
                       extra = list(hr.threshold = 0.90, quiet = FALSE,
                                    fs.splits = 20L, use_grf = FALSE,
                                    use_lasso = FALSE))
  args$hr.consistency <- NULL           # genuinely unsupplied: missing() TRUE

  msgs <- character(0)
  fit <- withCallingHandlers(
    do.call(forestsearch, c(list(df.analysis = df), args)),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage")
    },
    warning = function(w) invokeRestart("muffleWarning"))

  # The consistency stage received 0.72, on the natural scale and on the log
  # scale the comparison actually uses.
  expect_equal(fit$threshold_config$consistency_natural, 0.72,
               tolerance = 1e-12)
  expect_equal(fit$threshold_config$consistency, log(0.72), tolerance = 1e-12)
  expect_equal(fit$threshold_config$screening_natural, 0.90, tolerance = 1e-12)

  # Announced exactly once for the whole fit.
  expect_identical(sum(grepl("derived as 0.80 *", msgs, fixed = TRUE)), 1L)

  # And carried into args_call_all, so a replicate or a fold resolves 0.72.
  expect_equal(fit$args_call_all$consistency.threshold, 0.72,
               tolerance = 1e-12)
  expect_equal(fit$args_call_all$hr.consistency, 0.72, tolerance = 1e-12)
})


test_that("c2 > c1 errors before any model is fit", {
  # df.analysis is too small for any model to be fitted, and the columns the
  # search would need are absent: if the pair check ran after any fitting, the
  # error raised would be a data error, not this one.
  df <- data.frame(id = 1:3, treat = c(0L, 1L, 0L), y = c(0L, 1L, 0L))
  args <- .fs_args_for("binary", extra = list(hr.threshold = 0.90,
                                              hr.consistency = 1.00))

  expect_error(
    do.call(forestsearch, c(list(df.analysis = df), args)),
    "c2 > c1 not allowed for FS: hr.consistency = 1.00 exceeds hr.threshold = 0.90",
    fixed = TRUE)

  # The same call with a legal pair gets past the check and fails (or returns)
  # on the data instead -- proof that the check is what fired above.
  args$hr.consistency <- 0.80
  err <- tryCatch({
    suppressWarnings(suppressMessages(
      do.call(forestsearch, c(list(df.analysis = df), args))))
    NA_character_
  }, error = function(e) conditionMessage(e))
  expect_false(isTRUE(grepl("c2 > c1 not allowed", err, fixed = TRUE)))
})


# ---------------------------------------------------------------------------
# The rider (separable): the latent match.arg default
# ---------------------------------------------------------------------------

test_that("make_effect_estimator()'s binary choice order defaults to OR", {
  # match.arg() returns the FIRST choice for a NULL argument, so the order of
  # the choices vector is a default in its own right.  It is unreachable today
  # -- the resolution above it always yields a length-1 effect_measure -- but
  # it must not contradict that resolution.
  b <- gsub("[[:space:]]+", " ",
            paste(deparse(body(make_effect_estimator)), collapse = " "))
  expect_identical(
    length(gregexpr('c("OR", "RD", "RR", "IRR", "IRD")', b,
                    fixed = TRUE)[[1L]]), 1L)
  expect_identical(
    gregexpr('c("RD", "OR", "RR", "IRR", "IRD")', b, fixed = TRUE)[[1L]][1L],
    -1L)

  # The reachable default is unchanged, and an explicit measure still wins.
  df <- .make_binary_data(N = 60L, seed = 3L)
  f_default <- make_effect_estimator(outcome_type = "binary",
                                     treat.name = "treat", outcome.name = "y")
  f_or <- make_effect_estimator(outcome_type = "binary", effect_measure = "OR",
                                treat.name = "treat", outcome.name = "y")
  f_rd <- make_effect_estimator(outcome_type = "binary", effect_measure = "RD",
                                treat.name = "treat", outcome.name = "y")
  expect_equal(f_default(df)$estimate, f_or(df)$estimate)
  expect_false(isTRUE(all.equal(f_default(df)$estimate, f_rd(df)$estimate)))
})
