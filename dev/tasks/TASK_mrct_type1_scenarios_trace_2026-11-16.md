# TASK — MRCT articles: where are the type-1 error scenarios of the extended model? (2026-11-16)

Session name: `fs-mrct-type1-trace`

Host `pop-os`, repository `~/Documents/GitHub/forestsearch`, branch **`master`**.
Larry asked on 2026-11-16: the two MRCT articles have "completely lost the type-1 error
simulation scenarios based on the subgroup effect extended model". This task finds out what
happened. Read this whole file first.

## What is known (from the two files Larry sent to chat on 2026-11-16)

| File | Lines | md5 |
|---|---|---|
| `mrct-consistency.qmd` | 1,001 | `71aa9b658466713efcd58c370427eeed` |
| `mrct-consistency-simulations.qmd` | 1,428 | `2447dd482e9f0ccc5320bb5f6a0c9925` |

- The application article still **describes** the extended model: treatment-by-region
  interactions `\gamma_R` and `\gamma_{H^*\!R}`, effects `(b_R, h_R)` inside the region
  (about lines 205-212 and 243-252), and says its configurations "with no regional treatment
  benefit" quantify "the false-positive rate of the subgroup pathway" (about lines 153-157).
  It then says "Every cell of this article is the base model."
- The simulations article has eight cells, all of the base model: N0, N1, S1, S2, S3, S4,
  S2_D1, S2_D2. Both builds pass `region_treat_loghr = NULL`. It holds no cell, table or
  record for a regional treatment effect.

So the text promises extended-model scenarios that neither article reports. What is not known:
whether those scenarios were ever in these articles in this repository, were on a branch or in
a file that never reached `master`, or exist only outside this repository.

## Rules

- **Read-only.** No file in the repository changes except this task file and the report. No
  article, no `R/`, no test, no record. No render, no simulation, no package install, no
  `R CMD check`, no test suite. No agents. One short report.
- In any other repository: read only. No `git fetch`, no checkout, no change of any kind.
- Work on `master`. **No push, pull, merge or tag.**
- No stop gates. Record what you find and go on.
- Evidence comes from the files and the git history, quoted with commit SHAs, paths and line
  numbers. No guesses about intent.

## Search terms

Use these throughout (case-insensitive where it applies):
`region_treat_loghr`, `gamma_R`, `b_R`, `h_R`, `type-1`, `type 1`, `type I`, `type-I`,
`false-positive`, `false positive`, `no regional treatment`, `extension`, `extended model`,
`treatment-by-region`.

## Step 0 — start

1. Copy this file into `dev/tasks/` and commit it (path-limited). Record the SHA of `master`.
2. Record `git status --short`. For `vignettes/articles/mrct-consistency.qmd` and
   `vignettes/articles/mrct-consistency-simulations.qmd`: line count and md5 of the working
   copy and of the copy at `HEAD`. Say whether each equals the md5 in the table above.

## Step 1 — the history of the two articles

1. For each article (and any earlier path it had: use `--follow`), list every commit that
   touched it on any branch (`--all`): SHA, date, subject, lines added and removed.
2. For each search term, `git log --all -S` (and `-G` where a pattern is needed) limited to
   `vignettes/`, `inst/extdata/` and `R/`. For every commit that **removed** matching lines
   from either article, quote the removed lines (a few lines of context) and the commit
   subject.
3. Name the commits with the largest deletions from each article and say, in one line each,
   what they removed (section headings are enough).
4. State plainly: was there ever, at any commit on any branch, a version of either article
   that ran, tabulated or read records for cells with a regional treatment effect
   (`region_treat_loghr` not `NULL`, or `b_R` / `h_R` given values)? If yes: the last commit
   that had it, the commit that removed it, and the section headings and cell names involved.

## Step 2 — outside `master`

1. Branches, local and remote, not merged into `master`: list them; for each, whether it
   touches either article or `inst/extdata/mrct-consistency/`, and whether the search terms
   appear in its version of the articles.
2. `git stash list`, with one line on any stash touching the articles.
3. Other copies on disk: files named `mrct-consistency*.qmd` (and `*.html`) anywhere under
   the repository (including `dev/`, `quarto/`, `legacy_vignette/`), under `~/Downloads`, and
   under the other repositories in `~/Documents/GitHub/`. For each: path, size, modification
   time, line count, and the number of matches of the search terms.

## Step 3 — the earlier renders

For `dev/article_renders/pre_0.5.0/mrct-consistency*.html` (and any other earlier render
found in Step 2.3): list the section headings, and say whether any heading or table concerns
extended-model or type-1 error scenarios. Compare the headings with those of the current
sources and list headings present in a render but absent from the current source.

## Step 4 — the shipped records

List the files in `inst/extdata/mrct-consistency/`. For each CSV with a `cell` or `setting`
column, list its distinct values. Say whether any record concerns cells other than the eight
named above. Quote what `PROVENANCE.md` says about any such cells.

## Step 5 — the package

Where is `region_treat_loghr` used in `R/`? From the source and the help text: what values
does it accept, and does `inject_mrct_structure()` (or another exported function) build a
mechanism with `b_R` / `h_R` as the application article describes? Does any test exercise it?

## Step 6 — the mrct-consistency repository (read-only)

If `~/Documents/GitHub/mrct-consistency` exists: record its branch and last commit; search
its tracked files for the search terms (file names and match counts, the ten files with most
matches); list the directories and result files that hold simulations of cells with a
regional treatment effect (paths, dates, cell names, replicates, trial size); and name the
document sections (technical report, manuscript, supplement) that report them. If the
repository is not there, say so.

## Step 7 — report and close

1. Write `dev/reports/REPORT_mrct_type1_scenarios_trace_2026-11-16.md`. First line, one of:
   - `REMOVED FROM THE ARTICLES in <sha>`
   - `NEVER ON MASTER: found in <branch or file>`
   - `NEVER IN THE ARTICLES: results exist in <where>`
   - `NOT FOUND ANYWHERE SEARCHED`
   Then one short section per step. Facts only: no recommendation and no fix.
2. Commit the report (path-limited). Copy it to `~/Downloads`.
3. Final checks (report the result): this session's commits touch only this task file and the
   report; nothing else in this repository or any other changed; nothing was pushed, pulled,
   merged or tagged.
4. Print the first line and the answer to Step 1.4 as the closing message.
