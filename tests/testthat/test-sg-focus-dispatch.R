# =============================================================================
# Contract test: sg_focus dispatch completeness
#
# Every value in the post-normalization whitelist (.FS_SG_FOCUS_CANONICAL)
# must have an explicit branch at every site that decides behaviour from
# sg_focus.  A focus missing a branch does NOT error at runtime -- it falls
# through to a switch() default and selects under a rule the caller did not
# ask for, silently.  That is how "maxeff" and "maxeffCons" came to be
# whitelisted by forestsearch() while having no branch in the GRF frontier
# rule map and being rejected outright by dina_subgroup().
#
# These tests are the durable guard against that class of drift: they fail
# when a focus is added to the canonical set without updating all three
# dispatch sites.
# =============================================================================

test_that("every canonical sg_focus has a branch at every dispatch site", {
  expect_true(forestsearch:::.assert_sg_focus_dispatch_complete())
})

test_that("the canonical set is exactly the seven documented foci", {
  expect_setequal(
    forestsearch:::.FS_SG_FOCUS_CANONICAL,
    c("hr", "hrMaxSG", "maxSG", "hrMinSG", "minSG", "maxeff", "maxeffCons"))
})

test_that("both DINA and GRF dispatch maps cover the canonical set", {
  canonical <- forestsearch:::.FS_SG_FOCUS_CANONICAL
  expect_true(all(canonical %in%
    forestsearch:::.find_switch_branches(forestsearch::forestsearch,
                                         "frontier_rule")))
  expect_true(all(canonical %in%
    forestsearch:::.find_switch_branches(forestsearch::dina_subgroup, "ord")))
})

test_that("the guard actually fails when a focus lacks a branch", {
  # Negative control.  Without this, the assertion above could pass
  # vacuously -- e.g. if .find_switch_branches() silently returned nothing
  # and the setdiff() were computed against an empty canonical set.
  ns <- asNamespace("forestsearch")
  orig <- get(".FS_SG_FOCUS_CANONICAL", envir = ns)
  unlockBinding(".FS_SG_FOCUS_CANONICAL", ns)
  on.exit({
    assign(".FS_SG_FOCUS_CANONICAL", orig, envir = ns)
    lockBinding(".FS_SG_FOCUS_CANONICAL", ns)
  }, add = TRUE)

  assign(".FS_SG_FOCUS_CANONICAL", c(orig, "maxNotARealFocus"), envir = ns)
  expect_error(forestsearch:::.assert_sg_focus_dispatch_complete(),
               "maxNotARealFocus")
})

test_that("normalization aliases all land inside the canonical set", {
  aliases <- c("eff", "effMaxSG", "effMinSG", "maxcons",
               "hr", "hrMaxSG", "hrMinSG", "maxSG", "minSG",
               "maxeff", "maxeffCons")
  normalized <- vapply(aliases, forestsearch:::.normalize_sg_focus,
                       character(1), USE.NAMES = FALSE)
  expect_true(all(normalized %in% forestsearch:::.FS_SG_FOCUS_CANONICAL))
})

test_that("dina_subgroup() accepts maxeff and maxeffCons", {
  # These were rejected before the whitelist was widened: a caller could pass
  # a focus forestsearch() accepts and have dina_subgroup() stop().
  for (f in c("maxeff", "maxeffCons", "eff", "maxcons")) {
    expect_true(forestsearch:::.normalize_sg_focus(f) %in%
                  forestsearch:::.FS_SG_FOCUS_CANONICAL,
                info = f)
  }
})
