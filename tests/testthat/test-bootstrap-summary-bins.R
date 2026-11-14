# A-147: summarize_bootstrap_subgroups() bins the bootstrap consistency
# (Pcons) and subgroup size (N_sg) into ranges whose labels describe
# left-closed bins: "<0.5", "0.5-0.7", ..., "0.85-0.95", ">=0.95" and
# "<50", "50-99", ..., "200-299", ">=300".  A value exactly on a break
# must therefore be counted in the bin whose label starts at that value
# (0.95 is ">=0.95", 50 is "50-99"), not in the bin below it.
#
# Every expected value below is read off the bin labels alone:
#   - a label "a-b" holds values v with a <= v < b (the next label
#     starts at b, so b itself belongs to the next label);
#   - "<a" holds v < a; ">=a" holds v >= a.
# The labels are the ones documented in the code and returned in the
# `Consistency Range` / `Size Range` columns.  No expected value was
# taken from the output of the fixed code.

skip_on_cran()

# One synthetic "found" bootstrap row per call, so exactly one bin of each
# distribution has Count == 1 and that bin's label is the bin of the value.
.a147_bins_of <- function(pcons, n_sg) {
  res <- data.table::data.table(
    any_found = 1L,
    Pcons     = pcons,
    N_sg      = n_sg,
    M.1       = "{age>=50}",
    M.2       = NA_character_
  )
  out <- summarize_bootstrap_subgroups(res, nb_boots = 1L, maxk = 2)
  c(
    consistency = out$consistency_dist$`Consistency Range`[
      out$consistency_dist$Count == 1L],
    size = out$size_dist$`Size Range`[out$size_dist$Count == 1L]
  )
}

test_that("A-147: each consistency value falls in the bin its label names", {
  pcons <- c(0.95, 0.90, 0.85, 1.00, 0.949)
  # By hand, from the labels (a <= v < b):
  #   0.95  : on the 0.95 break; ">=0.95" starts at 0.95        -> ">=0.95"
  #   0.90  : 0.85 <= 0.90 < 0.95                                -> "0.85-0.95"
  #   0.85  : on the 0.85 break; "0.85-0.95" starts at 0.85      -> "0.85-0.95"
  #   1.00  : 1.00 >= 0.95 (the top bin keeps its upper end, so a
  #           consistency of 1 is still counted)                 -> ">=0.95"
  #   0.949 : 0.85 <= 0.949 < 0.95                               -> "0.85-0.95"
  expected <- c(">=0.95", "0.85-0.95", "0.85-0.95", ">=0.95", "0.85-0.95")
  # Sizes held at 120 (an interior point of "100-149") so only Pcons varies.
  got <- vapply(pcons, function(p) .a147_bins_of(p, 120L)[["consistency"]],
                character(1))
  expect_identical(got, expected)
})

test_that("A-147: each subgroup size falls in the bin its label names", {
  n_sg <- c(50L, 49L, 100L, 300L, 299L)
  # By hand, from the labels (a <= v < b):
  #   50  : on the 50 break; "50-99" starts at 50                 -> "50-99"
  #   49  : 49 < 50                                               -> "<50"
  #   100 : on the 100 break; "100-149" starts at 100             -> "100-149"
  #   300 : on the 300 break; ">=300" starts at 300               -> ">=300"
  #   299 : 200 <= 299 < 300                                      -> "200-299"
  expected <- c("50-99", "<50", "100-149", ">=300", "200-299")
  # Pcons held at 0.92 (an interior point of "0.85-0.95") so only N_sg varies.
  got <- vapply(n_sg, function(n) .a147_bins_of(0.92, n)[["size"]],
                character(1))
  expect_identical(got, expected)
})

test_that("A-147: counts over all five rows at once match the labels", {
  res <- data.table::data.table(
    any_found = 1L,
    Pcons     = c(0.95, 0.90, 0.85, 1.00, 0.949),
    N_sg      = c(50L, 49L, 100L, 300L, 299L),
    M.1       = "{age>=50}",
    M.2       = NA_character_
  )
  out <- summarize_bootstrap_subgroups(res, nb_boots = 5L, maxk = 2)

  # Labels in the documented order (the rows of each distribution table).
  # These two label checks pass before and after the A-147 fix (on the
  # pre-fix code they pass while the per-value, Count and Percent checks
  # fail): the labels are the `labels` of the two cut() calls, and table()
  # of a factor keeps every level, zero counts included, in level order; the
  # fix changed only right = FALSE.  They anchor the row order that the
  # positional Count and Percent checks below rely on.
  expect_identical(
    out$consistency_dist$`Consistency Range`,
    c("<0.5", "0.5-0.7", "0.7-0.8", "0.8-0.85", "0.85-0.95", ">=0.95")
  )
  expect_identical(
    out$size_dist$`Size Range`,
    c("<50", "50-99", "100-149", "150-199", "200-299", ">=300")
  )

  # Counts tallied by hand from the per-value bins in the two tests above:
  #   ">=0.95"    : 0.95, 1.00          -> 2
  #   "0.85-0.95" : 0.90, 0.85, 0.949   -> 3
  #   all lower consistency bins        -> 0
  expect_identical(out$consistency_dist$Count, c(0L, 0L, 0L, 0L, 3L, 2L))
  #   "<50": 49 -> 1; "50-99": 50 -> 1; "100-149": 100 -> 1;
  #   "150-199": none -> 0; "200-299": 299 -> 1; ">=300": 300 -> 1
  expect_identical(out$size_dist$Count, c(1L, 1L, 1L, 0L, 1L, 1L))

  # Percent is the percentage of the successful iterations (the
  # 100 * Count / n_found definition in the function body, which the fix
  # does not touch); n_found = 5 here, so by hand 3/5 = 60, 2/5 = 40,
  # 1/5 = 20 and 0/5 = 0.
  expect_equal(out$consistency_dist$Percent, c(0, 0, 0, 0, 60, 40))
  expect_equal(out$size_dist$Percent, c(20, 20, 20, 0, 20, 20))
})
