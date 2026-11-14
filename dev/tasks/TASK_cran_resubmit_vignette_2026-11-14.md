# TASK — CRAN resubmission: short vignette, the long one becomes an article, MKL safeguard (2026-11-14)

Session name: `fs-cran-resubmit-vignette`

Host `pop-os`, repository `~/Documents/GitHub/forestsearch`, branch **`master`**, version 0.5.0.
Larry approved this on 2026-11-14 ("proceed with shortening the vignette as you see fit"). Read
this whole file first.

## Background

CRAN's incoming pre-test of 0.5.0 (2026-11-14) passed on Windows and Debian (`Status: OK`) and was
flagged for one thing: `Overall checktime 12 min > 10 min` on Windows. In that Windows log the
vignette re-build took 217 s, the R code checks 146 s, the examples 61 s, the manuals 90 s and
the tests 14 s. The vignette is the only large step the package controls.

The same email lists an additional issue for the released 0.4.0 under MKL: ten test failures,
all in tests that require bit-identical floating-point results (including a repeated evaluation
in one session that was not bit-identical). None of those test files runs on CRAN in 0.5.0. The
one CRAN-run test file, `test-core-forestsearch.R`, compares two same-seed calls with
`expect_identical`, which could fail the same way.

## What this task does

1. The current vignette `vignettes/forestsearch.qmd` (Table 1 of León and Anderson 2026b, run end
   to end) becomes the website article `vignettes/articles/gbsg-table1.qmd`.
2. A short vignette takes its place. Its text is the file
   `~/Downloads/forestsearch_vignette_short_2026-11-14.qmd`, written by chat and approved by
   Larry. It runs the forest-search identification only and shows the inference call unevaluated.
3. The core test's repeat comparison gets a tolerance.
4. One release check, and `cran-comments.md` for the resubmission.

## Rules

- **Lean run.** The main session does everything. No agents. One short report, plus the raw
  check log as its only companion file.
- **No full test suite with `NOT_CRAN=true`.** Step 7 runs the one release check Larry
  authorized for this task.
- **`R/`: documentation text only, and only under Step 3.3.** No function body changes. No
  expected value in any test is edited.
- The short vignette's text is not rewritten. Step 5 says what may be fixed if it fails to run.
- Work on `master`. Commits are path-limited. **No push, pull, merge or tag.**
- No install into the user or site library. Scratch: `/tmp/cran11`; test runs and renders use a
  scratch build first on `R_LIBS=/tmp/cran11/lib`.
- Only Step 0 has a stop gate. After it, record a failure and go on.

## Step 0 — start

1. Copy this file into `dev/tasks/` and commit it. Record the start SHA.
2. **Gate (the only stop):** branch `master`; `DESCRIPTION` says `Version: 0.5.0`;
   `git status --porcelain` is empty; `~/Downloads/forestsearch_vignette_short_2026-11-14.qmd`
   exists (record its md5).

## Step 1 — the long vignette becomes an article

1. `git mv vignettes/forestsearch.qmd vignettes/articles/gbsg-table1.qmd`.
2. In the moved file, in the prose only, replace the word "vignette" with "article" (for example
   "**This vignette** runs", "the vignette's numbers", "the vignette defines one"). The YAML
   header and its `%\Vignette…` lines stay as they are. List each replacement. Nothing else in
   the file changes.

## Step 2 — the short vignette

Copy `~/Downloads/forestsearch_vignette_short_2026-11-14.qmd` to `vignettes/forestsearch.qmd`,
byte for byte (confirm with md5).

## Step 3 — the site, the links, NEWS

1. **`_pkgdown.yml`:** list the new article `articles/gbsg-table1` in the articles index and in
   the Articles menu, first among the analysis articles, in the same way
   `survival-post-selection` is listed. Its menu text: "The GBSG analysis (Table 1), end to end".
2. **References to the vignette.** Search `README.md`, `NEWS.md`, `vignettes/`, `R/` and
   `_pkgdown.yml` for links to `articles/forestsearch.html`, for `vignette("forestsearch")`, and
   for text that says the vignette reproduces Table 1, runs the three identifiers, or runs the
   declaration calibration. Where the text describes that content, point it to
   `articles/gbsg-table1.html` (full URL `https://larry-leon.github.io/forestsearch/articles/gbsg-table1.html`
   where a full URL is used). A plain pointer to the vignette as the starting point stays. List
   every place found and what was done.
3. **Help pages:** if a roxygen comment under `R/` describes the vignette's old content, change
   those words only, run `devtools::document()`, and confirm that only the `.Rd` files of the
   edited help pages changed. List each change. If there is none, say so.
4. **`NEWS.md`:** add one bullet under `# forestsearch 0.5.0`, in the Documentation group:

   > The package vignette is now a short introduction: identification by forest search on the
   > GBSG trial, with the inference call shown. The complete analysis it used to run (Table 1 of
   > León and Anderson 2026b: the three identifiers, the selection-adjusted bounds and the
   > declaration calibration) is the website article *Interpretable subgroup identification with
   > forestsearch* (`articles/gbsg-table1`).

5. One commit for Steps 1 to 3.

## Step 4 — the core test: a tolerance on the repeat comparison

In `tests/testthat/test-core-forestsearch.R`, second block only:

- `expect_identical(roundtrip(.core_strip_timing(b)), roundtrip(.core_strip_timing(a)))` becomes
  `expect_equal(roundtrip(.core_strip_timing(b)), roundtrip(.core_strip_timing(a)), tolerance = 1e-8)`;
- the block's name and the header comment say "the same result" where they say "an identical
  result" or "an identical object", and the header comment gains one line: the comparison uses a
  tolerance because some linear-algebra libraries (MKL) do not reproduce floating-point results
  bit for bit, even within a session.

Nothing else in the file changes. Checks: the file passes in CRAN mode and with `NOT_CRAN=true`;
as a negative control, a second call with another `seedit` still fails the block. One commit.

## Step 5 — run what changed

1. `git archive` HEAD's package paths and `R CMD INSTALL` them into `/tmp/cran11/lib`.
2. Render the short vignette in scratch. Record the render wall clock and the fit time it
   prints. If it does not run: fix only an evident error in the code of a chunk, list the fix
   exactly, and change no prose. If the fix is not evident, record the error and go on.
3. Report these from the rendered page, verbatim: the identified subgroup line, the funnel
   table, and the unadjusted table. (The prose states 75, 82 and 61 patients and a progesterone
   cut of 26; say whether the rendered numbers agree.)
4. Render the moved article once by the site-render route (recorded time for the old vignette:
   88 s) and confirm it builds.

## Step 6 — `cran-comments.md`

1. Insert this section after "## Submission":

   > ## Resubmission
   >
   > This is a resubmission of 0.5.0. The first submission passed the incoming pre-tests on
   > Windows and Debian (Status: OK) and was flagged for the overall check time on Windows
   > (12 min > 10 min). In this version:
   >
   > - the vignette is now a short introduction; the long analysis it ran is an article on the
   >   package website and is not part of the package;
   > - on CRAN a small core set of tests of the main function is run; the full test suite is run
   >   locally before each release.
   >
   > The additional issue reported for 0.4.0 under MKL came from tests that required
   > bit-identical floating-point results. Those tests are no longer run on CRAN, and the
   > remaining CRAN-run comparison of repeated results uses a numerical tolerance.

2. Under "Test environments", make the win-builder entry read exactly
   `- win-builder: TO BE FILLED BY MAINTAINER`.
3. The counts line under "R CMD check results" is set from Step 7. Nothing else changes.
4. One commit, after Step 7.

## Step 7 — the release check (once)

1. From `/tmp/cran11/build`, run `R CMD build` on the repository and
   `R CMD check --as-cran` on the tarball. Do not set `NOT_CRAN` or `OMP_THREAD_LIMIT`.
   Nothing is written into the repository.
2. Record: errors, warnings, every NOTE verbatim; the timing lines (R code, examples, tests,
   vignette re-build, manuals) and the total; the tarball size. Beside each timing, give the
   figure from the last check (`REPORT_cran_core_tests_2026-11-13.md`).
3. Copy `00check.log` to `dev/reports/REPORT_cran_resubmit_vignette_2026-11-14_check_log.txt`.

## Step 8 — report and close

1. Write `dev/reports/REPORT_cran_resubmit_vignette_2026-11-14.md`. Short:
   - first line: READY FOR WIN-BUILDER or NOT READY (ready means: the short vignette renders,
     the core test passes in CRAN mode, and the release check has 0 errors and 0 warnings);
   - top table: one row per step — result, commit, key numbers;
   - Step 1.2's and Step 3's lists;
   - Step 5's times and the three items quoted from the rendered vignette;
   - the check result, every NOTE verbatim, the timing table (now beside before);
   - the full text of `cran-comments.md` as committed;
   - this session's commits. Facts only.
2. Commit the report and the check log (path-limited). Copy the report to `~/Downloads`.
3. **Final checks** (report the result): this session's commits touch only
   `vignettes/forestsearch.qmd`, `vignettes/articles/gbsg-table1.qmd`, `_pkgdown.yml`,
   `NEWS.md`, `README.md` or other files changed under Step 3.2 and 3.3 (listed),
   `tests/testthat/test-core-forestsearch.R`, `cran-comments.md`, this task file, the report and
   the check log; `DESCRIPTION` still says `Version: 0.5.0`; `git status --porcelain` is empty;
   nothing was pushed, merged or tagged.
4. Delete `/tmp/cran11/lib`. Keep `/tmp/cran11/build`.
5. Print the first line, the top table and the timing table as the closing message.

## Amendment 1 (Larry, in chat, 2026-11-14, at the Step 0 gate)

Context: the Step 0 gate failed as written because `git status --porcelain` listed two untracked
files, `AGENTS.md` and `CRAN-SUBMISSION`.

> Proceed, with one addition. Leave both files untracked and untouched. First add the line
> ^AGENTS\.md$ to .Rbuildignore, as its own commit, so that any build from this checkout leaves
> AGENTS.md out. Then continue the task as written, reading the Step 0 and final clean-tree
> conditions as "no modified tracked files"; list the two untracked files in the final status
> check. Build the repository as written in Step 7 and confirm AGENTS.md is not in the tarball.
