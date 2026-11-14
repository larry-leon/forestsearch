# C1: format_subgroup_summary_tables() colour-codes the rows of the
# consistency (Pcons) distribution produced by
# summarize_bootstrap_subgroups().  The highlight sets must use the labels
# the summary produces; a label written differently (">= 0.95" with a
# space, where the summary writes ">=0.95") silently matches no row, so
# that row is never highlighted.
#
# The checks below look at behaviour, not at the source text: which rows
# of the formatted gt table receive a style, read from the table's
# `_styles` record (one row per styled cell, with `rownum` and the style
# list).
#
# Expected values, by hand:
#   - The summary's labels and their row order are the documented bins
#     "<0.5", "0.5-0.7", "0.7-0.8", "0.8-0.85", "0.85-0.95", ">=0.95"
#     (left-closed: a label "a-b" holds a <= v < b; see
#     test-bootstrap-summary-bins.R).  table() of the binned factor keeps
#     every level, so the distribution always has these six rows.
#   - The colour each row should carry is read off the table's own
#     footnote, "Color coding: Green (>=0.85), Yellow (0.7-0.85),
#     Red (<0.7)":
#       "<0.5"      : every value < 0.5  < 0.7         -> red
#       "0.5-0.7"   : 0.5 <= v < 0.7, so v < 0.7       -> red
#       "0.7-0.8"   : 0.7 <= v < 0.8, inside 0.7-0.85  -> yellow
#       "0.8-0.85"  : 0.8 <= v < 0.85, inside 0.7-0.85 -> yellow
#       "0.85-0.95" : 0.85 <= v < 0.95, so v >= 0.85   -> green (bold)
#       ">=0.95"    : v >= 0.95 >= 0.85                -> green (bold)
#   - The fills are the three cell_fill() colours of the function's
#     tab_style() calls for green, yellow and red ("#d4edda", "#fff3cd",
#     "#f8d7da", the success / warning / danger pastel shades); the fix
#     does not touch them.  gt stores colours upper-cased, so they are
#     compared upper-cased.  Only the green style adds bold text.
#   - So every one of the six rows carries exactly one fill.
# What the tests check: which rows of the formatted table get each fill
# and which get bold text.  In the first test every row's fills must equal
# the one above (the second checks the ">=0.95" row only), so a row with
# no fill or with a second fill fails; a summary label missing from its
# highlight set (as ">=0.95" was) shows as a row without its fill.
# What they do not check: a label in a highlight set that matches no
# summary row (say an extra ">= 0.95" next to ">=0.95") styles no row, so
# it changes nothing these tests can see, and they still pass.
# No expected value was taken from the output of the fixed code.

skip_on_cran()

.c1_green  <- toupper("#d4edda")
.c1_yellow <- toupper("#fff3cd")
.c1_red    <- toupper("#f8d7da")

# One synthetic "found" bootstrap row per Pcons value; the other columns
# are the minimum summarize_bootstrap_subgroups() needs (as in
# test-bootstrap-summary-bins.R).
.c1_summary <- function(pcons) {
  res <- data.table::data.table(
    any_found = 1L,
    Pcons     = pcons,
    N_sg      = 120L,
    M.1       = "{age>=50}",
    M.2       = NA_character_
  )
  summarize_bootstrap_subgroups(res, nb_boots = length(pcons), maxk = 2)
}

# For each row of the formatted consistency table: its label, the fills
# applied to it (sorted, collapsed with "+", "" when none) and whether any
# style applied to it sets bold text.
.c1_row_styles <- function(tbl) {
  labels <- as.character(tbl[["_data"]][["Consistency Range"]])
  st <- tbl[["_styles"]]
  st <- st[st$locname == "data", , drop = FALSE]
  fills <- character(length(labels))
  bold <- logical(length(labels))
  for (i in seq_along(labels)) {
    s_i <- st$styles[st$rownum == i]
    f_i <- unique(unlist(lapply(s_i, function(s) s$cell_fill$color)))
    fills[i] <- paste(sort(toupper(f_i)), collapse = "+")
    bold[i] <- any(vapply(s_i, function(s) {
      identical(s$cell_text$weight, "bold")
    }, logical(1)))
  }
  data.frame(label = labels, fill = fills, bold = bold,
             stringsAsFactors = FALSE)
}

test_that("C1: every consistency row the summary produces gets its colour", {
  # One value inside each bin, by hand (a <= v < b):
  #   0.30 -> "<0.5"; 0.60 -> "0.5-0.7"; 0.75 -> "0.7-0.8";
  #   0.82 -> "0.8-0.85"; 0.90 -> "0.85-0.95"; 0.97 -> ">=0.95"
  s <- .c1_summary(c(0.30, 0.60, 0.75, 0.82, 0.90, 0.97))
  tabs <- format_subgroup_summary_tables(s, nb_boots = 6L)
  rs <- .c1_row_styles(tabs$consistency_dist)

  # Premise: the rows are the six labels the summary documents.
  expect_identical(
    rs$label,
    c("<0.5", "0.5-0.7", "0.7-0.8", "0.8-0.85", "0.85-0.95", ">=0.95")
  )

  # Each row carries exactly one fill, the one the footnote names.
  expect_identical(
    rs$fill,
    c(.c1_red, .c1_red, .c1_yellow, .c1_yellow, .c1_green, .c1_green)
  )
  # Bold goes with green only.
  expect_identical(rs$bold, c(FALSE, FALSE, FALSE, FALSE, TRUE, TRUE))

  # The row the task names: the top consistency bin is highlighted green.
  expect_identical(rs$fill[rs$label == ">=0.95"], .c1_green)
})

test_that("C1: a run where every bootstrap has Pcons >= 0.95 highlights that row", {
  # All five values are >= 0.95, so by hand every one falls in ">=0.95"
  # (1.00 included: the top bin keeps its upper end), Count = 5 there and
  # 0 in the five other rows.
  s <- .c1_summary(c(0.95, 0.96, 0.98, 0.99, 1.00))
  expect_identical(s$consistency_dist$Count, c(0L, 0L, 0L, 0L, 0L, 5L))

  tabs <- format_subgroup_summary_tables(s, nb_boots = 5L)
  rs <- .c1_row_styles(tabs$consistency_dist)
  # The one row that holds every bootstrap is the green, bold row.
  expect_identical(rs$fill[rs$label == ">=0.95"], .c1_green)
  expect_true(rs$bold[rs$label == ">=0.95"])
})
