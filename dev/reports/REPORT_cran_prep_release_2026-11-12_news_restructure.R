# news_restructure.R -- TASK_cran_prep_release_2026-11-12, Part 3 step 2.
# Usage: Rscript --vanilla news_restructure.R <NEWS.md in> <summary.md> <NEWS.md out>
#
# Moves every bullet above the "# forestsearch 0.4.0" heading, verbatim, under
# one "# forestsearch 0.5.0" heading, in five groups (the task's order):
# Breaking changes; Bug fixes; New features; Superseded and deprecated;
# Documentation.  A bullet is a top-level "* " line with every following line
# up to the next blank line, top-level "* " line or heading (no bullet in this
# region has an internal blank line; some sections are tight lists, with no
# blank line between bullets).  In the output every bullet is followed by one
# blank line.  Within a group the bullets keep their order in the input.  The
# opening summary is read from <summary.md>.  The old development-version
# paragraph (the only prose above the first section heading) is replaced by
# the summary.  Everything from "# forestsearch 0.4.0" on is copied unchanged.
#
# Group rules (input section -> group), with named exceptions:
#   0.4.1.9011 "New features"      -> New features
#   0.4.1.9011 "Behaviour changes" -> Bug fixes, except
#       bootstrap orientation (scripts must drop a negation step),
#       caller's random-number state (a later unseeded draw changes),
#       fpr_calibration() default and new error             -> Breaking changes
#       mr_estimates_table() without a bootstrap template      -> New features
#   0.4.1.9011 "Removed"           -> Breaking changes
#   0.4.1.9011 "Superseded"        -> Superseded and deprecated
#   0.4.1.9011 "Documentation"     -> Documentation
#   0.4.1.9011 "Internal"          -> Documentation (test files; no behaviour change)
#   0.4.1 "Deprecated"             -> Superseded and deprecated
#   0.4.1 "Internal": fit_causal_forest_glm() removed -> Superseded and deprecated;
#       survival (>= 3.5-8) required -> Breaking changes;
#       lean Cox fit helper and probe -> New features (with Performance)
#   0.4.1 "Performance"            -> New features
#   0.4.1 "Bug fixes"              -> Bug fixes
#   0.4.1 "Documentation"          -> Documentation

a <- commandArgs(TRUE)
stopifnot(length(a) == 3)
x <- readLines(a[1], encoding = "UTF-8")
summ <- readLines(a[2], encoding = "UTF-8")

cut <- which(x == "# forestsearch 0.4.0")
stopifnot(length(cut) == 1)
head_part <- x[seq_len(cut - 1)]
tail_part <- x[cut:length(x)]

groups <- c("Breaking changes", "Bug fixes", "New features",
            "Superseded and deprecated", "Documentation")

exceptions <- list(
  "* `forestsearch_bootstrap_dofuture()` on GLM outcomes now builds the estimator" = "Breaking changes",
  "* `forestsearch()` and the other exported functions that can reach the fixed"   = "Breaking changes",
  "* `fpr_calibration()` with neither `outcome_type_approx` nor"                   = "Breaking changes",
  "* `mr_estimates_table()` no longer needs a bootstrap template:"                  = "New features",
  "* The internal helper `fit_causal_forest_glm()` is removed (unused since"        = "Superseded and deprecated",
  "* `survival (>= 3.5-8)` is now required: the lean unadjusted search-stage"       = "Breaking changes",
  "* The lean unadjusted Cox fit lives in the internal helper"                       = "New features"
)
section_rule <- c(
  "New features" = "New features", "Behaviour changes" = "Bug fixes",
  "Removed" = "Breaking changes", "Superseded" = "Superseded and deprecated",
  "Documentation" = "Documentation", "Internal" = "Documentation",
  "Deprecated" = "Superseded and deprecated", "Performance" = "New features",
  "Bug fixes" = "Bug fixes")

version <- NA_character_; section <- NA_character_
bullets <- list(); i <- 1L; n <- length(head_part)
while (i <= n) {
  l <- head_part[i]
  if (startsWith(l, "# ")) { version <- sub("^# ", "", l); section <- NA_character_; i <- i + 1L; next }
  if (startsWith(l, "### ")) { section <- sub("^### ", "", l); i <- i + 1L; next }
  if (startsWith(l, "* ")) {
    j <- i
    while (j + 1L <= n && nzchar(head_part[j + 1L]) &&
           !startsWith(head_part[j + 1L], "* ") &&
           !startsWith(head_part[j + 1L], "#")) j <- j + 1L
    blk <- head_part[i:j]
    stopifnot(!is.na(section))
    grp <- exceptions[[blk[1]]]
    if (is.null(grp)) grp <- section_rule[[section]]
    bullets[[length(bullets) + 1L]] <- list(text = blk, version = version,
                                            section = section, group = grp, line = i)
    i <- j + 1L; next
  }
  # blank line, or the development-version paragraph (replaced by the summary)
  i <- i + 1L
}
used <- vapply(bullets, function(b) b$text[1], "")
stopifnot(all(names(exceptions) %in% used))   # every exception names a real bullet

out <- c("# forestsearch 0.5.0", "", summ, "")
for (g in groups) {
  bs <- Filter(function(b) b$group == g, bullets)
  if (!length(bs)) next
  out <- c(out, paste("###", g), "")
  for (b in bs) out <- c(out, b$text, "")
}
out <- c(out, tail_part)
con <- file(a[3], open = "wb")
writeLines(out, con, useBytes = TRUE)   # LF endings; input ends with a newline
close(con)

map <- data.frame(
  n = seq_along(bullets),
  input_line = vapply(bullets, `[[`, 0L, "line"),
  from = paste(vapply(bullets, `[[`, "", "version"), "/", vapply(bullets, `[[`, "", "section")),
  to = vapply(bullets, `[[`, "", "group"),
  first_line = substr(used, 1, 70))
print(map, right = FALSE)
cat("\nbullets:", length(bullets), "\n")
print(table(factor(map$to, levels = groups)))
