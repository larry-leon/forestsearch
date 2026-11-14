DONE

# REPORT — MRCT articles: general installation note (2026-11-14)

Session `fs-mrct-install-note`, host `pop-os`, `~/Documents/GitHub/forestsearch`, branch `master`.
Task: `dev/tasks/TASK_mrct_articles_install_note_2026-11-14.md`. Start SHA
`f675541dd05966567d2ca565d41b19123f5e9784` (= `origin/master`). Text only. No render, no `R CMD check`,
no test suite, no agents.

## Step 0 — gate (passed)

Branch `master`; `git status --porcelain` empty after the task commit. Each gate string occurred
exactly once in each file (counted as occurrences, `grep -oF | wc -l`):

| string | mrct-consistency.qmd | mrct-consistency-simulations.qmd |
|---|---|---|
| `The CRAN release of forestsearch (0.4.0) does not yet contain the functions this article uses` | 1 | 1 |
| `This article needs the development version of forestsearch (see 'Installing forestsearch'); ` | 1 | 1 |

## Steps 1 and 2 — edits

- Edit 1: the callout (`::: {.callout-note}` … `:::`), lines 62–88 of `mrct-consistency.qmd` and
  59–85 of `mrct-consistency-simulations.qmd`, replaced with the text between the task's marker
  lines, taken from the committed task file by a script (not retyped). The new callout occurs
  verbatim exactly once in each file.
- Edit 2: the `stop(` line of the `install-check` chunk, `the development version of forestsearch`
  → `forestsearch 0.5.0 or later`. Nothing else in the chunk changed.

## Step 3 — checks

1. **Diff.** `git diff --stat` listed only the two files (37 lines each; 48 insertions, 26
   deletions in all). The changed lines are the same in both files. Hunk headers:
   - `vignettes/articles/mrct-consistency.qmd`
     - `@@ -62,29 +62,40 @@ companion article [mrct-consistency-simulations](mrct-consistency-simulations.ht`
     - `@@ -94,7 +105,7 @@ needed <- c("inject_mrct_structure", "ahr_population", "ahr_design", "get_dfpred`
   - `vignettes/articles/mrct-consistency-simulations.qmd`
     - `@@ -59,29 +59,40 @@ t_article <- proc.time()[["elapsed"]]`
     - `@@ -92,7 +103,7 @@ needed <- c("inject_mrct_structure", "scale_dgm_hazard", "ahr_population", "ahr_`

   In each file the first hunk is the callout, the second the one `stop(` line.

2. **Strings.** Both gate strings now occur 0 times in each file. **The count of
   `forestsearch 0.5.0 or later` does not match the task's figure:** it occurs **2** times in each
   file, not 3. The task's own callout text writes its second mention as
   `packageVersion("forestsearch")   # this article needs 0.5.0 or later`, so the literal string
   `forestsearch 0.5.0 or later` appears only on the callout's first line and on the `stop(` line.
   The shorter string `0.5.0 or later` occurs 3 times in each file: twice in the callout, once in the
   chunk. That is the split the task describes in parentheses. The text was kept exactly as given.

   | file | `forestsearch 0.5.0 or later` | `0.5.0 or later` (lines) |
   |---|---|---|
   | mrct-consistency.qmd | 2 | 3 (65, 78, 108) |
   | mrct-consistency-simulations.qmd | 2 | 3 (62, 75, 106) |

3. **Parse.** `knitr::purl()` (knitr 1.52) into the session scratch directory under `/tmp`, then
   `parse()`:
   - `mrct-consistency.qmd`: 399 purled lines, 150 top-level expressions, parsed without error.
   - `mrct-consistency-simulations.qmd`: 690 purled lines, 247 top-level expressions, parsed without
     error.

   Both purled scripts contain the new `stop(` message. The callout's two ```` ```r ```` blocks are
   plain Markdown code, not R chunks, so `purl()` leaves them out. They were also parsed on their own
   (4 expressions each, both files, no error). Nothing was rendered.

4. **CRAN tarball.** `.Rbuildignore` line 20: `^vignettes/articles$`. Neither article goes into the
   package tarball.

5. **Other files under `vignettes/` (listed only, nothing changed).** No other file mentions the
   development version or a GitHub install. Four places state a CRAN release number, or "the current
   CRAN release", in an installation note or install check:
   - `vignettes/articles/extreme_subgroups.qmd:210–211` (section "Running this article yourself"):
     "This article is written for the current CRAN release of **forestsearch** (0.4.0) and runs as a
     single file."
   - `vignettes/articles/extreme_subgroups.qmd:26` (setup chunk comment, above a
     `packageVersion("forestsearch") < "0.4.0"` check whose `stop()` says "This article needs
     forestsearch >= 0.4.0. Update it with install.packages(\"forestsearch\")."): "# This article is
     written for the current CRAN release of forestsearch."
   - `vignettes/articles/guohe-nested-continuum.qmd:37–38` (chunk comment, above the same `>= 0.4.0`
     check and `stop()` at 39–42): "# This article is written for the current CRAN release of
     forestsearch; the # two functions it uses are exported from 0.4.0 on."
   - `vignettes/articles/actg175-glm-outcomes.qmd:25` (setup chunk, above the same `>= 0.4.0` check
     and `stop()` at 26–29): "# Written for forestsearch >= 0.4.0, the release the package website
     renders with."

   These match the search but are not installation notes, so they are left out of the list above:
   `survival-post-selection.qmd:20`, a version-history note ("In forestsearch 0.4.0 this was a CRAN
   vignette; from 0.4.1 it is an article on the package website …"); `methodology.qmd:83`, a link to
   the source code on GitHub; `extreme_subgroups.qmd:224–229`, a `download.file()` of the `.qmd` from
   raw.githubusercontent.com, which fetches the article file and does not install the package;
   `actg175-glm-outcomes.qmd:569, 606` and `survival-post-selection.qmd:386`, timings recorded under
   forestsearch 0.2.0. `vignettes/forestsearch.qmd` has no installation note.

## Step 4 — this session's commits

| SHA | message | files |
|---|---|---|
| `a3ec983f` | dev: task MRCT articles install note (as received) | `dev/tasks/TASK_mrct_articles_install_note_2026-11-14.md` |
| `9f91302b` | docs: MRCT articles, general installation note (Steps 1 and 2) | `vignettes/articles/mrct-consistency.qmd`, `vignettes/articles/mrct-consistency-simulations.qmd` |
| (the commit that adds this report) | dev: MRCT articles install note report (Step 4) | `dev/reports/REPORT_mrct_articles_install_note_2026-11-14.md` |

## Final checks

These were run after the two article files were committed, before this report's own commit. The
closing message re-checks them after that commit.

- `git diff --name-only f675541d..HEAD`: the task file and the two article files only.
- `git status --porcelain`: empty.
- Nothing pushed: `origin/master` is still `f675541d` (the start SHA), and `master` is 2 commits
  ahead. No pull, merge or tag: the reflog shows only this session's two commits since the start,
  and no tag contains the start SHA.

A copy of this report is in `~/Downloads`.
