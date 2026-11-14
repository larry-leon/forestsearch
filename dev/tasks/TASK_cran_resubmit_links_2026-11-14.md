# TASK — CRAN resubmission: links that resolve today (2026-11-14)

Session name: `fs-cran-resubmit-links`

Host `pop-os`, repository `~/Documents/GitHub/forestsearch`, branch **`master`**, version 0.5.0.
Read this whole file first.

## Why

The last check (`dev/reports/REPORT_cran_resubmit_vignette_2026-11-14.md`) has one NOTE: the URL
`https://larry-leon.github.io/forestsearch/articles/gbsg-table1.html` returns 404. It is linked
from the package vignette and from `README.md`, both of which are in the CRAN tarball. The page
does not exist until the website is rebuilt, and the website build on GitHub has failed on every
push since 2026-11-13. Larry's reading is that the site renders the articles with the CRAN
release of forestsearch (0.4.0), which lacks functions the MRCT articles need, so the site
cannot rebuild until 0.5.0 is on CRAN. CRAN would flag a dead link, so the two files in the
tarball must link to a page that exists today: the articles index,
`https://larry-leon.github.io/forestsearch/articles/`.

## Rules

- **Text only,** in `vignettes/forestsearch.qmd`, `README.md` and `cran-comments.md`. Nothing
  under `R/`, `man/`, `tests/` or `inst/` changes. The articles under `vignettes/articles/` keep
  their links to the new article (they are not in the tarball, and the links resolve once the
  site is rebuilt).
- **Lean run.** No agents. One short report, plus the raw check log.
- **No test suite.** Step 4 runs the one release check Larry authorized for this task.
- Work on `master`. Commits are path-limited. **No push, pull, merge or tag.**
- No install into the user or site library. Scratch: `/tmp/cran12`, scratch build first on
  `R_LIBS=/tmp/cran12/lib`.
- Untracked `AGENTS.md` and `CRAN-SUBMISSION` are Larry's: leave them untouched. Read every
  clean-tree condition as "no modified tracked files".

## Step 0 — start

1. Copy this file into `dev/tasks/` and commit it. Record the start SHA.
2. **Gate (the only stop):** branch `master`; `DESCRIPTION` says `Version: 0.5.0`; no modified
   tracked files; the string `articles/gbsg-table1.html` occurs exactly twice in
   `vignettes/forestsearch.qmd`.

## Step 1 — how the website build installs the package (read-only)

Read `.github/workflows/pkgdown.yaml`. Report, quoting the lines: whether the site build installs
forestsearch from the repository or from CRAN, and which step is step 6 of the `pkgdown` job (the
step the failed runs stop at). If `gh` is installed and signed in, also quote the last 30 lines of
the failed step from the most recent failed run (`gh run view --log-failed`); if it is not, say
so and go on. Change nothing.

## Step 2 — the links

**`vignettes/forestsearch.qmd`, place 1.** Replace

~~~~text
to end — is the article
[Interpretable subgroup identification with forestsearch](https://larry-leon.github.io/forestsearch/articles/gbsg-table1.html)
on the package website.
~~~~

with

~~~~text
to end — is the article *Interpretable subgroup identification with
forestsearch*, listed under Articles on the
[package website](https://larry-leon.github.io/forestsearch/articles/).
~~~~

**`vignettes/forestsearch.qmd`, place 2.** Replace

~~~~text
- The article
  [Interpretable subgroup identification with forestsearch](https://larry-leon.github.io/forestsearch/articles/gbsg-table1.html)
  continues this analysis: the selection-adjusted bounds for this region and
~~~~

with

~~~~text
- The article *Interpretable subgroup identification with forestsearch*,
  listed under Articles on the
  [package website](https://larry-leon.github.io/forestsearch/articles/),
  continues this analysis: the selection-adjusted bounds for this region and
~~~~

**`README.md`.** In each place that links to `articles/gbsg-table1.html`, name the article in
italics and link the words "package website" to
`https://larry-leon.github.io/forestsearch/articles/`, in the same form as above ("the article
*Interpretable subgroup identification with forestsearch*, listed under Articles on the
[package website](…)"). Keep the rest of each sentence. Quote each sentence before and after in
the report.

**Check:** `gbsg-table1.html` no longer occurs in `vignettes/forestsearch.qmd` or `README.md`.
`curl` the articles index URL and record its status (200 expected). One commit.

## Step 3 — render

`git archive` HEAD's package paths, `R CMD INSTALL` into `/tmp/cran12/lib`, and render the
vignette once in scratch. Record that it renders and its wall clock.

## Step 4 — the release check (once)

1. From `/tmp/cran12/build`: `R CMD build` on the repository, then `R CMD check --as-cran` on
   the tarball. Do not set `NOT_CRAN` or `OMP_THREAD_LIMIT`. Nothing is written into the
   repository.
2. Record errors, warnings and every NOTE verbatim, the total check time, and that neither
   `AGENTS.md` nor `CRAN-SUBMISSION` is in the tarball.
3. Copy `00check.log` to `dev/reports/REPORT_cran_resubmit_links_2026-11-14_check_log.txt`.

## Step 5 — `cran-comments.md`

Under "R CMD check results": set the counts line from Step 4, and make the sentence that
follows it state Step 4's status exactly (for example "The local check above ended with
`Status: OK`."). Nothing else in the file changes. One commit.

## Step 6 — report and close

1. Write `dev/reports/REPORT_cran_resubmit_links_2026-11-14.md`: first line READY FOR
   WIN-BUILDER (the check has 0 errors, 0 warnings and no URL NOTE) or NOT READY; Step 1's
   findings; Step 2's before and after; the check result with every NOTE verbatim; the
   `cran-comments.md` text as committed; this session's commits.
2. Commit the report and the check log. Copy the report to `~/Downloads`.
3. **Final checks** (report the result): this session's commits touch only
   `vignettes/forestsearch.qmd`, `README.md`, `cran-comments.md`, this task file, the report and
   the check log; no modified tracked files; nothing was pushed, merged or tagged.
4. Delete `/tmp/cran12/lib`. Keep `/tmp/cran12/build`.
5. Print the first line and Step 1's findings as the closing message.
