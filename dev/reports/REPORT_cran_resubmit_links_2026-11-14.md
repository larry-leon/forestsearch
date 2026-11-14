READY FOR WIN-BUILDER

# REPORT — CRAN resubmission: links that resolve today (2026-11-14)

Task: `dev/tasks/TASK_cran_resubmit_links_2026-11-14.md`. Session `fs-cran-resubmit-links`, host
`pop-os`, branch `master`, version 0.5.0. Start SHA `670a4664`; task commit `91c6cc29`. Main session
only: no agents, no test suite. Every run used the scratch build `R_LIBS=/tmp/cran12/lib` (`git archive`
of HEAD's package paths at `1731ddb4`, `R CMD INSTALL`). R 4.6.2, Quarto 1.10.18.

READY because the release check has 0 errors, 0 warnings and no NOTE at all (`Status: OK`). The URL
NOTE of the previous check (`articles/gbsg-table1.html`, 404) is gone: the vignette and `README.md` now
link the articles index, which returns HTTP 200 today.

**Step 0 gate held:** branch `master`; `DESCRIPTION` says `Version: 0.5.0`; 0 modified tracked files
(`git status --porcelain` empty); `articles/gbsg-table1.html` occurred exactly twice in
`vignettes/forestsearch.qmd` (lines 79 and 304). `AGENTS.md` and `CRAN-SUBMISSION` were not at the
repository root this session (`ls`: no such file), so there was nothing of Larry's to leave untouched.

| Step | Result | Commit | Key numbers |
|---|---|---|---|
| 0 | Gate held | `91c6cc29` task | start `670a4664` |
| 1 | Read only | — | site installs forestsearch from CRAN; step 6 of the run is "Build site"; `gh` not installed |
| 2 | Done | `1731ddb4` | 4 places (2 vignette, 2 README); `gbsg-table1.html` now 0 in both files; articles index HTTP 200 |
| 3 | Renders | — | install 17 s; vignette render 7.06 s wall, exit 0 |
| 4 | 0 errors, 0 warnings, 0 NOTEs, `Status: OK` | log in the report commit | check 268 s (4 min 28 s; previous check 4 min 29.1 s); build 36 s; tarball 1,638,479 bytes; no `AGENTS.md` or `CRAN-SUBMISSION` in it |
| 5 | Done | `9c8305d3` | counts line now `0 errors \| 0 warnings \| 0 notes`; status sentence already right, unchanged |
| 6 | This report | report commit | final checks below |

## Step 1 — how the website build installs the package (read only)

`.github/workflows/pkgdown.yaml`, unchanged:

```text
34      - uses: r-lib/actions/setup-r-dependencies@v2
35        with:
36          extra-packages: any::pkgdown, any::katex, cran::forestsearch
37          needs: website
38
39      - name: Build site
40        run: pkgdown::build_site_github_pages(new_process = FALSE, install = FALSE)
41        shell: Rscript {0}
```

- **forestsearch comes from CRAN, not from the repository.** Line 36 asks for `cran::forestsearch`, and
  line 40 passes `install = FALSE`, so pkgdown does not install the checked-out source. The articles are
  therefore rendered with the CRAN release (0.4.0 today), which agrees with Larry's reading. The
  `setup-r-dependencies` step (line 34) sets no `packages:` input, so it installs the action's default,
  which per the r-lib/actions documentation is the checkout's dependencies, not the checkout itself (not
  verified here).
- **Step 6 is "Build site".** The `pkgdown` job lists six steps: `actions/checkout@v4` (line 26),
  `setup-pandoc` (28), `setup-r` (30), `setup-r-dependencies` (34), `Build site` (39–41) and
  `Deploy to GitHub pages 🚀` (43–49). GitHub's run view counts "Set up job" as step 1, so its step 6 is
  `Build site`, line 40 quoted above. (Counting the YAML entries alone, the sixth would be the deploy,
  which runs only after a successful build; a run stopping there would already have a built site.)
- **Failed-run log: not read.** `gh` is not installed on this host (`gh: command not found`), so the last
  30 lines of the failed step could not be quoted.
- Seen while checking Step 2's URL: the live articles index (navbar "Released version" 0.4.1) lists
  *Interpretable subgroup identification with forestsearch* today, linking to `forestsearch.html`, the
  0.4.1-era render of the old vignette (the same GBSG analysis before it moved). The moved article
  `vignettes/articles/gbsg-table1.qmd` keeps that title in its YAML, so after the site is rebuilt the same
  entry will link to `gbsg-table1.html`. The italic title in the new text names an index entry both today
  and after the rebuild.

## Step 2 — the links

Line numbers are those at the start SHA `670a4664`.

**`vignettes/forestsearch.qmd`, place 1 (lines 77–80).** Before:

```text
calibration, which together are Table 1 of León and Anderson (2026b) run end
to end — is the article
[Interpretable subgroup identification with forestsearch](https://larry-leon.github.io/forestsearch/articles/gbsg-table1.html)
on the package website.
```

After:

```text
calibration, which together are Table 1 of León and Anderson (2026b) run end
to end — is the article *Interpretable subgroup identification with
forestsearch*, listed under Articles on the
[package website](https://larry-leon.github.io/forestsearch/articles/).
```

**`vignettes/forestsearch.qmd`, place 2 (lines 303–305).** Before:

```text
- The article
  [Interpretable subgroup identification with forestsearch](https://larry-leon.github.io/forestsearch/articles/gbsg-table1.html)
  continues this analysis: the selection-adjusted bounds for this region and
```

After:

```text
- The article *Interpretable subgroup identification with forestsearch*,
  listed under Articles on the
  [package website](https://larry-leon.github.io/forestsearch/articles/),
  continues this analysis: the selection-adjusted bounds for this region and
```

Both vignette replacements are the task's text, verbatim.

**`README.md`, place 1 (Quick start, lines 76–79).** Before:

```text
This is the forest-search analysis of the German Breast Cancer Study Group
trial in Table 1 of León and Anderson (2026b). The article
[Interpretable subgroup identification with forestsearch](https://larry-leon.github.io/forestsearch/articles/gbsg-table1.html)
runs it in full, with the two model-based identifiers, and reads every number.
```

After:

```text
This is the forest-search analysis of the German Breast Cancer Study Group
trial in Table 1 of León and Anderson (2026b). The article
*Interpretable subgroup identification with forestsearch*, listed under
Articles on the [package website](https://larry-leon.github.io/forestsearch/articles/),
runs it in full, with the two model-based identifiers, and reads every number.
```

**`README.md`, place 2 (Documentation, first bullet, lines 137–138).** Before:

```text
- The article [Interpretable subgroup identification with forestsearch](https://larry-leon.github.io/forestsearch/articles/gbsg-table1.html):
  the analysis above, end to end.
```

After:

```text
- The article *Interpretable subgroup identification with forestsearch*,
  listed under Articles on the
  [package website](https://larry-leon.github.io/forestsearch/articles/):
  the analysis above, end to end.
```

`README.md` is the source (there is no `README.Rmd`). The rest of each sentence is unchanged.

**Check.** `gbsg-table1.html` occurs 0 times in `vignettes/forestsearch.qmd` and 0 times in `README.md`.
`curl https://larry-leon.github.io/forestsearch/articles/` returned **HTTP 200** on 2026-11-14 11:46:28
PST (the old page, `articles/gbsg-table1.html`, still returns 404). The articles under
`vignettes/articles/` keep their links to the new article, as the task says. One commit, `1731ddb4`
(`README.md` +5/−2, `vignettes/forestsearch.qmd` +6/−5).

## Step 3 — render

`git archive 1731ddb4` of the package paths (`.Rbuildignore DESCRIPTION NAMESPACE LICENSE NEWS.md
README.md R inst man tests vignettes/forestsearch.qmd`, 7.5 MB) into `/tmp/cran12/src`; `R CMD INSTALL
--library=/tmp/cran12/lib` exit 0 in 17 s; under `R_LIBS=/tmp/cran12/lib` the first library path is
`/tmp/cran12/lib` and `packageVersion("forestsearch")` is 0.5.0. `quarto render forestsearch.qmd` on a
copy in `/tmp/cran12/render`: **renders**, exit 0, **7.06 s** wall. The HTML has no `gbsg-table1`, and
both new "package website" links point to `https://larry-leon.github.io/forestsearch/articles/`.

## Step 4 — the release check (once)

From `/tmp/cran12/build`, with `R_LIBS=/tmp/cran12/lib` and neither `NOT_CRAN` nor `OMP_THREAD_LIMIT`
set (both unset in the script; neither was in the environment, and there is no `~/.R/check.Renviron`):

- `R CMD build ~/Documents/GitHub/forestsearch` (HEAD `1731ddb4`, tree clean): exit 0, **36 s**,
  `forestsearch_0.5.0.tar.gz`, 1,638,479 bytes, 749 entries.
- `R CMD check --as-cran forestsearch_0.5.0.tar.gz`: exit 0, **268 s** (4 min 28 s, 11:48:20 to
  11:52:48 PST).

**Result: 0 errors, 0 warnings, 0 NOTEs. `Status: OK`.** There is no NOTE to quote: `00check.log` has
no line containing `NOTE`, `WARNING` or `ERROR`. Its URL-bearing check passed:

```text
* checking CRAN incoming feasibility ... [7s/15s] OK
```

Timed lines from the log: install 20 s, R code for possible problems 45 s, examples 35 s, examples with
`--run-donttest` 62 s, PDF manual 13 s, HTML manual 18 s; tests and the vignette re-build are under the
10 s print threshold. The built `README.md` and `inst/doc/forestsearch.html` contain no `gbsg-table1`.
Neither `AGENTS.md` nor `CRAN-SUBMISSION` is in the tarball (`tar -tzf` lists neither; neither file
exists in the checkout today). Nothing was written into the repository: after the check, `git status
--porcelain` listed only this report, untracked at the time.

`00check.log` is copied verbatim (71 lines, `cmp` identical) to
`dev/reports/REPORT_cran_resubmit_links_2026-11-14_check_log.txt`.

## Step 5 — `cran-comments.md`

One line changed under "R CMD check results": `0 errors | 0 warnings | 1 note` →
`0 errors | 0 warnings | 0 notes`. The sentence after it already read "The local check above ended with
`Status: OK`.", which is Step 4's status exactly, so it is unchanged. Nothing else in the file changed
(`git show --stat 9c8305d3`: 1 file, +1/−1). Commit `9c8305d3`. The file as committed (44 lines):

```markdown
## Submission

This is an update from forestsearch 0.4.0 (the current CRAN version) to
0.5.0. NEWS.md lists the changes under "forestsearch 0.5.0".

This release removes `create_gbsg_dgm()` and `simulate_from_gbsg_dgm()`,
which were deprecated in 0.4.0 (the first CRAN release to export them,
where they already warned that they were deprecated).

## Resubmission

This is a resubmission of 0.5.0. The first submission passed the incoming pre-tests on
Windows and Debian (Status: OK) and was flagged for the overall check time on Windows
(12 min > 10 min). In this version:

- the vignette is now a short introduction; the long analysis it ran is an article on the
  package website and is not part of the package;
- on CRAN a small core set of tests of the main function is run; the full test suite is run
  locally before each release.

The additional issue reported for 0.4.0 under MKL came from tests that required
bit-identical floating-point results. Those tests are no longer run on CRAN, and the
remaining CRAN-run comparison of repeated results uses a numerical tolerance.

## Test environments

- Local: Pop!_OS 24.04 LTS, x86_64-pc-linux-gnu, R 4.6.2 (2026-10-29);
  `R CMD check --as-cran` on the tarball built by `R CMD build`.
- win-builder: TO BE FILLED BY MAINTAINER

## R CMD check results

0 errors | 0 warnings | 0 notes

The local check above ended with `Status: OK`.

On CRAN a small core set of tests of the main function `forestsearch()` is run. The full
test suite is run locally before each release.

## Reverse dependencies

There are no reverse dependencies: CRAN's package index (read on
2026-11-12) lists no package that depends on, imports, links to or suggests
forestsearch.
```

The win-builder line is still "TO BE FILLED BY MAINTAINER", as before this task.

## This session's commits

All on `master`, on top of the start SHA `670a4664`, each path-limited:

| Commit | Message | Paths |
|---|---|---|
| `91c6cc29` | task: cran-resubmit-links 2026-11-14 (as received) | `dev/tasks/TASK_cran_resubmit_links_2026-11-14.md` |
| `1731ddb4` | docs: link the articles index instead of the unbuilt gbsg-table1 page (Step 2) | `README.md`, `vignettes/forestsearch.qmd` |
| `9c8305d3` | docs: cran-comments.md check counts from the links check (Step 5) | `cran-comments.md` |
| (this commit) | dev: CRAN resubmission links report and check log (Step 6) | this report, `dev/reports/REPORT_cran_resubmit_links_2026-11-14_check_log.txt` |

## Final checks

Run with this report and the check log staged, then again after the commit, with the same outcome
(`/tmp/cran12/final_checks.sh`, read only):

- **Paths: PASS.** This session's commits (`670a4664..HEAD`, `--no-renames`) touch exactly
  `vignettes/forestsearch.qmd`, `README.md`, `cran-comments.md`,
  `dev/tasks/TASK_cran_resubmit_links_2026-11-14.md`, this report and
  `dev/reports/REPORT_cran_resubmit_links_2026-11-14_check_log.txt`. Nothing under `R/`, `man/`, `tests/`
  or `inst/` changed (0 paths).
- **Tree: PASS.** 0 modified tracked files; no untracked files.
- **No push, merge or tag: PASS.** `origin/master` is still `670a4664`; its last reflog entry is
  "update by push" at 11:28:49 PST, before this session's first commit (11:45:36). `HEAD` is 4 ahead of
  it and 0 behind. No merge commits in the range. No tag points into the range; the newest tag is still
  `v0.4.0` (2026-10-09). No `git push`, `pull`, `fetch`, `merge` or `tag` command was run.

Clean-up: `/tmp/cran12/lib` deleted after these checks; `/tmp/cran12/build` kept (tarball, `build.log`,
`check.log`, `forestsearch.Rcheck/`). `/tmp/cran12/src` and `/tmp/cran12/render` (the export and the
scratch render) are also left in place. The report is copied to `~/Downloads` under the same name.
