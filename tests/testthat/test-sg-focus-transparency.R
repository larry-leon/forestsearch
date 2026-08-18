# =============================================================================
# sg_focus transparency: alias announcement + the selection invariant
# =============================================================================
# The announcement change set is required to be SELECTION-NEUTRAL.  These
# tests pin both halves: that the messages say the right thing, and that
# turning them on changed no subgroup and no estimate.

# ---------------------------------------------------------------------------
# The alias table and its two readers
# ---------------------------------------------------------------------------

test_that("the alias table and .normalize_sg_focus() agree in both directions", {
  al <- forestsearch:::.FS_SG_FOCUS_ALIASES
  expect_identical(al,
                   c(effMaxSG = "hrMaxSG", effMinSG = "hrMinSG",
                     eff = "hr", maxcons = "hr"))
  # Forwards: every alias normalizes to its recorded canonical form.
  for (a in names(al)) {
    expect_identical(forestsearch:::.normalize_sg_focus(a), unname(al[[a]]),
                     info = a)
  }
  # Every canonical focus is a fixed point.
  for (cf in forestsearch:::.FS_SG_FOCUS_CANONICAL) {
    expect_identical(forestsearch:::.normalize_sg_focus(cf), cf, info = cf)
  }
  # Backwards: the announcement's alias class.
  expect_identical(forestsearch:::.sg_focus_aliases_of("hr"),
                   c("eff", "maxcons"))
  expect_identical(forestsearch:::.sg_focus_aliases_of("hrMaxSG"), "effMaxSG")
  expect_identical(forestsearch:::.sg_focus_aliases_of("hrMinSG"), "effMinSG")
  expect_identical(forestsearch:::.sg_focus_aliases_of("maxeff"), character(0))
})

test_that("non-character / non-scalar input passes through the normalizer", {
  expect_identical(forestsearch:::.normalize_sg_focus(42L), 42L)
  expect_identical(forestsearch:::.normalize_sg_focus(c("eff", "hr")),
                   c("eff", "hr"))
})

# ---------------------------------------------------------------------------
# .announce_sg_focus()
# ---------------------------------------------------------------------------

.ann <- function(...) {
  msgs <- character(0)
  withCallingHandlers(
    forestsearch:::.announce_sg_focus(...),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    })
  msgs
}

test_that("an alias emits exactly one message naming the full alias class", {
  m <- .ann("eff", "hr", "consistency", quiet = FALSE)
  expect_length(m, 1L)
  expect_identical(
    trimws(m),
    "sg_focus 'eff' resolves to canonical rule 'hr' (aliases: eff, maxcons).")

  m2 <- .ann("maxcons", "hr", "consistency", quiet = FALSE)
  expect_length(m2, 1L)
  expect_match(m2, "resolves to canonical rule 'hr'", fixed = TRUE)

  m3 <- .ann("effMaxSG", "hrMaxSG", "consistency", quiet = FALSE)
  expect_length(m3, 1L)
  expect_match(m3, "(aliases: effMaxSG).", fixed = TRUE)
})

test_that("a canonical focus emits nothing on the consistency path", {
  expect_length(.ann("hr", "hr", "consistency", quiet = FALSE), 0L)
  expect_length(.ann("maxSG", "maxSG", "consistency", quiet = FALSE), 0L)
  expect_length(.ann("maxeffCons", "maxeffCons", "consistency",
                     quiet = FALSE), 0L)
})

test_that("maxeff / maxeffCons announce the effect-argmax collapse on dina and grf", {
  for (m in c("dina", "grf")) {
    for (f in c("maxeff", "maxeffCons")) {
      out <- .ann(f, f, m, quiet = FALSE)
      expect_length(out, 1L)
      expect_identical(
        trimws(out),
        sprintf(paste0("sg_focus '%s' ranks as the effect argmax on the %s ",
                       "path (no Pcons is computed)."), f, m),
        info = paste(m, f))
    }
  }
})

test_that("at most one message fires: the two conditions are disjoint", {
  # No alias resolves to maxeff/maxeffCons, so no input can trigger both.
  al <- forestsearch:::.FS_SG_FOCUS_ALIASES
  expect_length(intersect(unname(al), c("maxeff", "maxeffCons")), 0L)
  for (m in c("consistency", "dina", "grf")) {
    for (a in names(al)) {
      expect_lte(length(.ann(a, unname(al[[a]]), m, quiet = FALSE)), 1L)
    }
  }
})

test_that("quiet = TRUE silences the announcement completely", {
  expect_length(.ann("eff", "hr", "dina", quiet = TRUE), 0L)
  expect_length(.ann("maxeffCons", "maxeffCons", "dina", quiet = TRUE), 0L)
})

test_that("the announcement uses message(), never cat()", {
  # cat() goes to stdout, which suppressMessages() cannot silence and which
  # the sim cells must not see.
  out <- capture.output(
    suppressMessages(
      forestsearch:::.announce_sg_focus("eff", "hr", "dina", quiet = FALSE)))
  expect_identical(out, character(0))
})

# ---------------------------------------------------------------------------
# The invariant: announcing changes no selection
# ---------------------------------------------------------------------------

test_that("dina selects identically under eff / hr / maxcons / maxeffCons", {
  skip_on_cran()
  df <- .make_survival_data(N = 200L, HR_harm = 2.0, tau = 60, seed = 42L)
  args <- .fs_args_for("survival",
                       confounders = c("age", "stage", "sex", "noise"),
                       extra = list(fs.splits = 100L, maxk = 2L, n.min = 40L,
                                    pconsistency.threshold = 0.60,
                                    hr.threshold = 1.10, use_grf = FALSE,
                                    quiet = TRUE, seedit = 42L,
                                    subgroup_method = "dina"))

  sels <- lapply(c("eff", "hr", "maxcons", "maxeffCons"), function(f) {
    fs <- suppressMessages(suppressWarnings(do.call(
      forestsearch, c(list(df.analysis = df),
                      modifyList(args, list(sg_focus = f))))))
    fs$sg.harm
  })

  # The documented fact: on DINA these four are ONE rule.  Asserting it here
  # makes the collapse a tested property rather than a code comment.
  expect_identical(sels[[2]], sels[[1]])
  expect_identical(sels[[3]], sels[[1]])
  expect_identical(sels[[4]], sels[[1]])
  expect_false(is.null(sels[[1]]))
})

test_that("consistency selects identically under the eff / hr / maxcons aliases", {
  skip_on_cran()
  df <- .make_survival_data(N = 200L, HR_harm = 2.0, tau = 60, seed = 42L)
  args <- .fs_args_for("survival",
                       confounders = c("age", "stage", "sex", "noise"),
                       extra = list(fs.splits = 100L, maxk = 2L, n.min = 40L,
                                    pconsistency.threshold = 0.60,
                                    hr.threshold = 1.10, use_grf = FALSE,
                                    quiet = TRUE, seedit = 42L,
                                    subgroup_method = "consistency"))

  sels <- lapply(c("eff", "hr", "maxcons"), function(f) {
    fs <- suppressMessages(suppressWarnings(do.call(
      forestsearch, c(list(df.analysis = df),
                      modifyList(args, list(sg_focus = f))))))
    fs$sg.harm
  })
  expect_identical(sels[[2]], sels[[1]])
  expect_identical(sels[[3]], sels[[1]])
  expect_false(is.null(sels[[1]]))
})

test_that("forestsearch() announces the alias at run time and quiet suppresses it", {
  skip_on_cran()
  df <- .make_survival_data(N = 200L, HR_harm = 2.0, tau = 60, seed = 42L)
  args <- .fs_args_for("survival",
                       confounders = c("age", "stage", "sex", "noise"),
                       extra = list(fs.splits = 50L, maxk = 1L, n.min = 40L,
                                    pconsistency.threshold = 0.60,
                                    hr.threshold = 1.10, use_grf = FALSE,
                                    seedit = 42L,
                                    subgroup_method = "dina",
                                    sg_focus = "eff"))

  msgs <- character(0)
  suppressWarnings(withCallingHandlers(
    do.call(forestsearch, c(list(df.analysis = df),
                            modifyList(args, list(quiet = FALSE)))),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    }))
  expect_true(any(grepl("resolves to canonical rule 'hr'", msgs, fixed = TRUE)))

  # quiet = TRUE: the announcement must not appear.
  msgs_q <- character(0)
  suppressWarnings(withCallingHandlers(
    do.call(forestsearch, c(list(df.analysis = df),
                            modifyList(args, list(quiet = TRUE)))),
    message = function(m) {
      msgs_q <<- c(msgs_q, conditionMessage(m))
      invokeRestart("muffleMessage")
    }))
  expect_false(any(grepl("resolves to canonical rule", msgs_q, fixed = TRUE)))
})

# ---------------------------------------------------------------------------
# The GRF frontier rule map is exhaustive (Item 7)
# ---------------------------------------------------------------------------

test_that("every canonical focus has an explicit GRF frontier branch", {
  # .assert_sg_focus_dispatch_complete() reads the dispatch sites themselves.
  # If the frontier switch loses a branch, this is where it surfaces.
  expect_true(forestsearch:::.assert_sg_focus_dispatch_complete())

  branches <- forestsearch:::.find_switch_branches(forestsearch, "frontier_rule")
  expect_setequal(branches, forestsearch:::.FS_SG_FOCUS_CANONICAL)
})

test_that("the GRF frontier switch has no silent default arm", {
  # Previously an unrecognised focus fell through to "effMaxSG" with no
  # condition raised -- a run ranked by a rule the caller never named.
  src <- paste(deparse(body(forestsearch)), collapse = "\n")
  expect_match(src, "GRF frontier rule: no branch for", fixed = TRUE)
  # And the error-swallowing re-normalization is gone.
  expect_false(grepl(".sgf <- tryCatch", src, fixed = TRUE))
})
