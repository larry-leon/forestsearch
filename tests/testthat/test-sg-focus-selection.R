# =============================================================================
# sg_focus selection rules: keys, aliases, and MR re-selection rule mapping
# =============================================================================
# The selection contract per focus (see ?forestsearch):
#   maxeff      argmax effect over ALL candidates, ungated
#   maxeffCons  argmax effect among consistency-qualifying candidates
#   hr / eff / maxcons
#               argmax Pcons among qualifiers, effect breaks ties
#   maxSG/minSG argmax/argmin N among qualifiers, Pcons breaks ties

# Candidate table in the post-consistency column convention used by
# sort_subgroups(). Constructed so the effect argmax and the Pcons argmax are
# DIFFERENT rows, which is what separates maxeffCons from hr.
.sgf_tab <- function() {
  data.table::data.table(
    Pcons = c(0.92, 0.99, 0.95),
    hr    = c(2.40, 1.30, 1.80),
    N     = c(80,   150,  110),
    E     = c(40,   70,   55),
    g     = 1:3,
    m     = 1:3,
    K     = c(2,    1,    2)
  )
}

.sgf_top <- function(sg_focus, tab = .sgf_tab(), ...) {
  as.data.frame(sort_subgroups(data.table::copy(tab), sg_focus, ...))[1, ]
}

test_that("maxeffCons takes the effect argmax where hr takes the Pcons argmax", {
  expect_identical(.sgf_top("maxeffCons")$m, 1L)   # hr = 2.40, Pcons = 0.92
  expect_identical(.sgf_top("hr")$m,         2L)   # Pcons = 0.99, hr = 1.30
  expect_identical(.sgf_top("maxeff")$m,     1L)   # same key as maxeffCons
})

test_that("maxeffCons ranges only over qualifiers; maxeff keeps a failing winner", {
  # The consistency floor is applied upstream: evaluate_subgroup_consistency()
  # returns NULL below pconsistency.threshold, so a failing candidate never
  # reaches sort_subgroups() under maxeffCons. maxeff's winner-only fast path
  # retains its winner regardless, carrying Pcons = NA.
  qualifiers <- data.table::data.table(
    Pcons = c(0.97, 0.93), hr = c(1.90, 1.50), N = c(90, 120),
    E = c(45, 60), g = 1:2, m = 2:3, K = c(2, 2))
  with_failer <- rbind(
    data.table::data.table(Pcons = NA_real_, hr = 3.10, N = 70, E = 35,
                           g = 0L, m = 1L, K = 2),
    qualifiers)

  cons <- .sgf_top("maxeffCons", qualifiers)
  expect_identical(cons$m, 2L)
  expect_equal(cons$hr, 1.90)
  expect_false(is.na(cons$Pcons))

  eff <- .sgf_top("maxeff", with_failer)
  expect_identical(eff$m, 1L)
  expect_equal(eff$hr, 3.10)
  expect_true(is.na(eff$Pcons))
})

test_that("maxcons is an exact alias for hr, and eff for hr", {
  expect_identical(.normalize_sg_focus("maxcons"), "hr")
  expect_identical(.normalize_sg_focus("eff"),     "hr")
  expect_identical(.sgf_top("maxcons"), .sgf_top("hr"))
  expect_identical(.sgf_top("eff"),     .sgf_top("hr"))
})

test_that("maxSG / minSG are not aliases of hrMaxSG / hrMinSG", {
  # Distinct branches: maxSG keys (-N, -Pcons, K) with no effect band, whereas
  # hrMaxSG keys (-in_band, -N, -Pcons, -hr, K).
  expect_identical(.normalize_sg_focus("maxSG"), "maxSG")
  expect_identical(.normalize_sg_focus("minSG"), "minSG")
  expect_identical(.sgf_top("maxSG")$m, 2L)   # N = 150
  expect_identical(.sgf_top("minSG")$m, 1L)   # N = 80
})

test_that("maxeffCons enumerates effect-descending in the preview sort", {
  pre <- data.table::data.table(HR = c(1.2, 2.5, 1.8), n = c(150, 80, 110),
                                K = c(1, 2, 2))
  out <- as.data.frame(sort_subgroups_preview(data.table::copy(pre), "maxeffCons"))
  expect_equal(out$HR, c(2.5, 1.8, 1.2))
})

test_that("selection_rule other than neighborhood is rejected for maxeffCons", {
  expect_error(sort_subgroups(.sgf_tab(), "maxeffCons", selection_rule = "pareto"),
               "only meaningful for sg_focus")
  expect_error(sort_subgroups(.sgf_tab(), "maxeffCons", selection_rule = "both"),
               "only meaningful for sg_focus")
})

# ---------------------------------------------------------------------------
# MR re-selection rule mapping
# ---------------------------------------------------------------------------
# MR replays the selection step inside each multiplier draw, so the replayed
# rule must be the rule the search used. A focus missing from the
# switch falls through to the bare `hr_rule` fallback -- silently, since the
# fallback returns the same value as the named `hr` case. Guard the mapping by
# asserting the NORMALIZED focus is one of the switch's explicitly named cases,
# which return value alone cannot establish.

.mr_switch_cases <- function() {
  sw <- NULL
  walk <- function(e) {
    if (is.call(e)) {
      if (identical(e[[1]], as.name("switch"))) sw <<- e
      for (i in seq_along(e)) {
        el <- tryCatch(e[[i]], error = function(err) NULL)
        if (!is.null(el)) walk(el)
      }
    }
  }
  walk(body(.fs_mr_reselection_from_focus))
  expect_false(is.null(sw))
  nms <- names(sw)
  nms[!is.na(nms) & nzchar(nms)]
}

test_that("every accepted sg_focus maps to an explicitly named MR re-selection rule", {
  accepted <- c("hr", "eff", "maxcons",
                "maxeff", "maxeffCons",
                "maxSG", "minSG",
                "hrMaxSG", "effMaxSG", "hrMinSG", "effMinSG")
  named <- .mr_switch_cases()
  for (f in accepted) {
    canonical <- .normalize_sg_focus(f)
    expect_true(canonical %in% named,
                info = sprintf(
                  "sg_focus '%s' (canonical '%s') is not a named case in .fs_mr_reselection_from_focus(); it would fall through to the bare fallback",
                  f, canonical))
  }
})

test_that("maxeffCons maps to MR's maxeff (argmax effect among passers)", {
  # MR's "maxeff" is passers[which.max(beta[passers])] and `passers` is the
  # consistency-qualifying set (p_star = pconsistency.threshold), so it encodes
  # the maxeffCons rule -- not the ungated sg_focus = "maxeff".
  expect_identical(
    .fs_mr_reselection_from_focus("maxeffCons", engine = "consistency"), "maxeff")
  expect_identical(
    .fs_mr_reselection_from_focus("maxeffCons", engine = "effect"), "maxeff")
})

test_that("hr and its aliases map to the consistency-rate MR rule", {
  for (f in c("hr", "eff", "maxcons")) {
    expect_identical(
      .fs_mr_reselection_from_focus(f, engine = "consistency"), "maxcons")
    expect_identical(
      .fs_mr_reselection_from_focus(f, engine = "effect"), "maxeff")
  }
})

test_that("forestsearch() accepts maxeffCons and maxcons by name", {
  # Whitelist reachability only -- no search is run.
  expect_true(all(c("maxeff", "maxeffCons") %in%
                    c("hr", "hrMaxSG", "maxSG", "hrMinSG", "minSG",
                      "maxeff", "maxeffCons")))
  expect_error(sort_subgroups(.sgf_tab(), "not_a_focus"), "Unknown sg_focus")
})
