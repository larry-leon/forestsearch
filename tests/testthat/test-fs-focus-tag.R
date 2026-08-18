# =============================================================================
# fs_focus_tag(): the stem-tag matrix, asserted exhaustively
# =============================================================================
# The matrix below is the contract, not an illustration.  It is the thing the
# simulation drivers used to each re-derive by hand, and the DINA maxeffCons
# cell is where they got it wrong: the run ranks by the effect argmax there,
# so the stem must say "eff", not "maxeffCons".

.FT_MATRIX <- read.table(
  header = TRUE, stringsAsFactors = FALSE, text = "
sg_focus     consistency  dina      grf
eff          maxcons      eff       eff
hr           maxcons      eff       eff
maxcons      maxcons      eff       eff
maxeff       maxeff       eff       eff
maxeffCons   maxeffCons   eff       eff
effMaxSG     effMaxSG     effMaxSG  effMaxSG
hrMaxSG      effMaxSG     effMaxSG  effMaxSG
effMinSG     effMinSG     effMinSG  effMinSG
hrMinSG      effMinSG     effMinSG  effMinSG
maxSG        maxSG        maxSG     maxSG
minSG        minSG        minSG     minSG
")

test_that("fs_focus_tag() reproduces the documented 11 x 3 matrix exactly", {
  for (i in seq_len(nrow(.FT_MATRIX))) {
    f <- .FT_MATRIX$sg_focus[i]
    for (m in c("consistency", "dina", "grf")) {
      expect_identical(
        fs_focus_tag(m, f), .FT_MATRIX[[m]][i],
        info = sprintf("fs_focus_tag(%s, %s)", m, f))
    }
  }
  # Guard the matrix itself: 11 spellings, no duplicates.
  expect_identical(nrow(.FT_MATRIX), 11L)
  expect_identical(anyDuplicated(.FT_MATRIX$sg_focus), 0L)
})

test_that("the matrix covers every accepted spelling, canonical and alias", {
  # If a canonical focus or an alias is ever added, this fails until the
  # matrix above grows a row for it -- which is the point.
  accepted <- union(forestsearch:::.FS_SG_FOCUS_CANONICAL,
                    names(forestsearch:::.FS_SG_FOCUS_ALIASES))
  expect_identical(sort(accepted), sort(.FT_MATRIX$sg_focus))
})

test_that("dina and grf collapse all five effect-family spellings onto 'eff'", {
  five <- c("eff", "hr", "maxcons", "maxeff", "maxeffCons")
  for (m in c("dina", "grf")) {
    expect_identical(unique(vapply(five, fs_focus_tag,
                                   character(1), subgroup_method = m)),
                     "eff",
                     info = m)
  }
  # ...and that consistency does NOT, which is the distinction the tag exists
  # to preserve.
  expect_identical(fs_focus_tag("consistency", "maxeffCons"), "maxeffCons")
  expect_identical(fs_focus_tag("consistency", "maxeff"),     "maxeff")
  expect_identical(fs_focus_tag("consistency", "eff"),        "maxcons")
})

test_that("fs_focus_tag() agrees with the tag the dina path actually uses", {
  # dina_subgroup()'s ordering switch and the GRF frontier_rule switch both
  # key the effect argmax as "eff".  If either is renamed, the stem tag and
  # the rule that ran part company silently -- so pin the string.
  expect_identical(fs_focus_tag("dina", "maxeffCons"), "eff")
  expect_identical(fs_focus_tag("grf",  "maxeff"),     "eff")
})

test_that("fs_focus_tag() rejects non-scalar / NA input", {
  expect_error(fs_focus_tag("dina", c("eff", "hr")), "single non-NA")
  expect_error(fs_focus_tag("dina", NA_character_),  "single non-NA")
  expect_error(fs_focus_tag(NA_character_, "eff"),   "single non-NA")
  expect_error(fs_focus_tag(1L, "eff"),              "single non-NA")
})

test_that("an unrecognised focus passes through unchanged", {
  # The default arm is load-bearing: maxSG/minSG rely on it.  Documented so a
  # future 'validate here' change does not break them.
  expect_identical(fs_focus_tag("dina", "maxSG"), "maxSG")
  expect_identical(fs_focus_tag("consistency", "not_a_focus"), "not_a_focus")
})
