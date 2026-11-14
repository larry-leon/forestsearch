READY FOR WIN-BUILDER

# REPORT — CRAN resubmission: short vignette, the long one becomes an article, MKL safeguard (2026-11-14)

Task: `dev/tasks/TASK_cran_resubmit_vignette_2026-11-14.md` (with Amendment 1). Session
`fs-cran-resubmit-vignette`, host `pop-os`, branch `master`, version 0.5.0. Start SHA `e9efc7ff`;
task commit `ca53c0e5`. Main session only, no agents. Every test run and render used the scratch build
`R_LIBS=/tmp/cran11/lib` (`git archive` of HEAD's package paths, `R CMD INSTALL`). R 4.6.2, testthat 3.3.2,
Quarto 1.10.18, pkgdown 2.2.1.

READY because the short vignette renders, the core test passes in CRAN mode (2 blocks, 42 expectations,
0 failures), and the release check has 0 errors and 0 warnings. The one NOTE is a 404 for the new
article's URL, which is not yet on the website.

**Step 0 gate, Amendment 1.** The gate failed as written: `git status --porcelain` listed two untracked
files, `AGENTS.md` (2026-11-14 10:14) and `CRAN-SUBMISSION` (09:53, SHA `e9efc7ff`). Larry answered in chat
(appended verbatim to the task file, `4997efe2`): leave both untouched and untracked; add `^AGENTS\.md$` to
`.Rbuildignore` as its own commit (`e2d447c5`); read the Step 0 and final clean-tree conditions as "no
modified tracked files"; confirm `AGENTS.md` is not in the Step 7 tarball. The gate then held: branch
`master`, `Version: 0.5.0`, 0 modified tracked files, short vignette present, md5
`93acaa320ac102a4825e45e0c692ba65`.

| Step | Result | Commit | Key numbers |
|---|---|---|---|
| 0 | Gate held under Amendment 1 (it failed as written) | `ca53c0e5` task, `4997efe2` Amendment 1, `e2d447c5` `.Rbuildignore` | start `e9efc7ff`; 0 modified tracked files; 2 untracked files (Larry's) |
| 1 | Done | `12332599` | `git mv` to `vignettes/articles/gbsg-table1.qmd`; 3 word replacements; YAML unchanged |
| 2 | Done | `12332599` | md5 `93acaa320ac102a4825e45e0c692ba65` for source and copy |
| 3 | Done | `12332599` | `_pkgdown.yml` 2 entries; 26 places found: 20 edits (16 pointers, 4 wordings), 6 kept; `R/` unchanged; 1 NEWS bullet |
| 4 | PASS | `a4fff743` | CRAN mode and `NOT_CRAN=true`: 2 blocks, 42 expectations, 0 failures, 1.4 s; negative control fails the block |
| 5 | Renders | — | install 19.3 s; vignette render 7.1 s wall, fit 1.7 s; numbers agree with prose; article by site route 106 s (old vignette's record: 88 s) |
| 6 | Done | `46c1fe50` | Resubmission section; win-builder line; counts `0 errors \| 0 warnings \| 1 note` |
| 7 | 0 errors, 0 warnings, 1 NOTE | log in the report commit | check 4 min 29.1 s (was 5 min 39.0 s); build 36.2 s; vignette re-build under 10 s (was 78 s); tarball 1,638,428 bytes; no `AGENTS.md` in it |
| 8 | This report | report commit | final checks below |

## Step 1.2 — "vignette" → "article" in the moved file's prose (3 replacements)

Line numbers are the same before and after the move (`vignettes/forestsearch.qmd` at `e9efc7ff` →
`vignettes/articles/gbsg-table1.qmd`).

| Line | Before | After |
|---|---|---|
| 119 | `**This vignette** runs the analysis reported as Table 1 …` | `**This article** runs the analysis reported as Table 1 …` |
| 216 | `… so the vignette's numbers are reproducible.` | `… so the article's numbers are reproducible.` |
| 323 | `… so the vignette defines one; …` | `… so the article defines one; …` |

The YAML header, including the `vignette: >` key and its three `%\Vignette…` lines (14–17), is unchanged.
`diff` of the old vignette against the article has 6 changed lines (3 removed, 3 added), and nothing else.

## Step 3 — the site, the links, NEWS

**3.1 `_pkgdown.yml`.** Articles menu: after "Methodology", before "Post-Selection Inference", added
`- text: The GBSG analysis (Table 1), end to end` / `href: articles/gbsg-table1.html`. Articles index,
group "Post-Selection Inference": `- articles/gbsg-table1` added first, unquoted, as
`articles/survival-post-selection` is. `yaml::read_yaml()` parses the file; menu item 4 and the group's
contents read back as written.

**3.2 References.** Searched `README.md`, `NEWS.md`, `vignettes/` and `R/` (with `command grep`, not the
shell's ignore-aware `grep` wrapper) and `_pkgdown.yml` for `articles/forestsearch.html`,
`vignette("forestsearch")`, and text saying the vignette reproduces Table 1, runs the three identifiers or
runs the declaration calibration. The new short vignette was not searched for edits (its text is not
rewritten). Line numbers are those at `e9efc7ff`. "Article link" below is
`[Interpretable subgroup identification with forestsearch](https://larry-leon.github.io/forestsearch/articles/gbsg-table1.html)`,
the title and link the short vignette itself uses.

| # | Place | What it said | Done |
|---|---|---|---|
| 1 | `_pkgdown.yml:32` | navbar "Getting Started" → `articles/forestsearch.html` | kept (plain pointer to the vignette) |
| 2 | `_pkgdown.yml:40` | Articles menu "Introduction to forestsearch" → `articles/forestsearch.html` | kept (plain pointer) |
| 3 | `README.md:77` | `vignette("forestsearch")` runs it [Table 1's forest search] in full, with the two model-based identifiers | pointer → "The article" + article link |
| 4 | `README.md:136` | `vignette("forestsearch")`: the analysis above, end to end | pointer → "The article" + article link |
| 5 | `NEWS.md:722` (0.5.0, Documentation) | New getting-started vignette, `vignette("forestsearch")`: … Table 1, run end to end — the three identifiers … declaration calibration | lead → "New article on the package website, *Interpretable subgroup identification with forestsearch* (`articles/gbsg-table1`):"; rest of the bullet unchanged |
| 6 | `NEWS.md:739` (0.5.0, Documentation) | The README … shows the analysis of `vignette("forestsearch")` | → "of the website article *Interpretable subgroup identification with forestsearch* (`articles/gbsg-table1`)" |
| 7 | `vignettes/articles/survival-post-selection.qmd:79` | `vignette("forestsearch")` runs the companion analysis … the search, the one-sided bounds and the declaration calibration | pointer → "The article" + article link |
| 8 | same file `:116` | The companion vignette uses ε = 0.20 … (N 75) | "vignette" → "article" (refers back to #7) |
| 9 | same file `:120` | … region of that vignette's candidate table | "vignette's" → "article's" (refers back to #7) |
| 10 | same file `:132` | … the one-sided bounds of the companion vignette | "vignette" → "article" (the bounds are Table 1's) |
| 11 | same file `:449` | … the declaration calibration … are the subject of … `vignette("forestsearch")` | pointer → "the article" + article link |
| 12 | same file `:453` | `vignette("forestsearch")`: … the three identifiers … the declaration calibration | pointer → "The article" + article link |
| 13 | `vignettes/articles/guohe-nested-continuum.qmd:875` | Forest search's screened identifier is evaluated … in `vignette("forestsearch")` | kept: the short vignette still runs forest search's screened identifier |
| 14 | `vignettes/articles/forestsearch-2024.Rmd:29` | current methodology … with its declaration calibration … documented in `vignette("forestsearch")` and the post-selection article | pointer → "the article" + article link |
| 15 | `vignettes/articles/extreme_subgroups.qmd:257–260` | The Getting Started vignette is also installed with the package: `vignette("forestsearch")` | kept (plain pointer) |
| 16 | `vignettes/articles/methodology.qmd:64` | current methodology [incl. bounds and declaration calibration, lines 60–62] is documented in `vignette("forestsearch")` | pointer → "the article" + article link |
| 17 | same file `:163` | `fs_declaration_calibration()` … `vignette("forestsearch")` runs it on the GBSG analysis | pointer → "the article" + article link |
| 18 | same file `:989` | Table 1 of León and Anderson (2026b), run end to end in `vignette("forestsearch")` | pointer → "the article" + article link |
| 19 | same file `:1326` | the GBSG forest search of `vignette("forestsearch")` — … inference and the field bounds included | pointer → "the article" + article link |
| 20 | same file `:1328` | … and the vignette prints its own timing | "vignette" → "article" |
| 21 | same file `:1340` | the current products [corrected estimate, field bounds, Bonferroni pair] (… 6.1; `vignette("forestsearch")`) | pointer → "the article" + article link |
| 22 | `vignettes/articles/paper_simulations.qmd:59` | … bounds, two-subgroup claim and declaration calibration … documented in `vignette("forestsearch")` | pointer → "the article" + article link |
| 23 | `vignettes/articles/actg175-glm-outcomes.qmd:123` | The binary search uses the settings of the breast cancer analysis in `vignette("forestsearch")` | kept: the short vignette's "The specification" section carries those settings |
| 24 | same file `:342` | `fs_declaration_calibration()` … (see `vignette("forestsearch")` for the construction) | pointer → "the article" + article link |
| 25 | same file `:620` | `vignette("forestsearch")` runs the same constructions … with the declaration calibration | pointer → "The article" + article link |
| 26 | `R/bootstrap_analysis_dofuture.R:61` | roxygen "(see vignette for references)" on the bootstrap | kept: generic, not the Table 1 content |

Other `vignette` hits under `R/` (`run_subgroup_sims.R`, `plot_subgroup_sims.R`, `compare_*`,
`simulation_tables.R`, `gg_forest.R`, `bootstrap_summaries_helpers.R`) are generic or name the
extreme-subgroups and simulation Quarto documents, not the package vignette's content. After the edits, the only `vignette("forestsearch")`
and `articles/forestsearch.html` pointers left are #1, #2, #13, #15 and #23.

Totals: 20 edits (16 pointers retargeted, 4 wordings), 5 kept, 1 roxygen place found and kept.
README.md no longer names `vignette("forestsearch")`. Both of its places described the full analysis.

**3.3 Help pages.** No roxygen comment under `R/` describes the vignette's old content (the one roxygen
mention, #26, is generic). No `R/` change, so `devtools::document()` was not run and no `.Rd` changed.

**3.4 NEWS.** The bullet was added as given, in `# forestsearch 0.5.0` → `### Documentation`, directly
after the bullet in #5.

## Step 4 — the core test

`tests/testthat/test-core-forestsearch.R`, second block: `expect_identical(roundtrip(.core_strip_timing(b)),
roundtrip(.core_strip_timing(a)))` became `expect_equal(…, tolerance = 1e-8)`. The block name is now
"a second call with the same seed returns the same result". The header says "returns the same result"
(was "an identical object") and gains one line, given as one physical line of 153 characters:
`# The comparison uses a tolerance because some linear-algebra libraries (MKL) do not reproduce
floating-point results bit for bit, even within a session.` Nothing else in the file changed (ASCII only,
final newline kept).

Checks (installed scratch build, `test_file()` on an export of `tests/`):

- CRAN mode (`NOT_CRAN` unset, `on_cran()` TRUE): 2 blocks, 41 + 1 expectations, 0 failures, 1.4 s.
- `NOT_CRAN=true` (`on_cran()` FALSE): the same, 1.4 s.
- Negative control: the second call with `seedit = 20261114` fails the block (1 failure). waldo reports one
  difference: `` `actual$args_call_all$seedit`: 20261114.0 `` against `8316951.0`. A side check found that
  the default call's result does not depend on the seed. Both seeds select `!{size <= 29} {er <= 8}`, and
  with the recorded `seedit` aligned the two stripped objects are `identical()`. The default
  `consistency_method = "resample"` is evaluated in closed form, so only the recorded argument tells the
  two calls apart.

## Step 5 — run what changed

1. `git archive a4fff743` (package paths) → `R CMD INSTALL --library=/tmp/cran11/lib`: 19.3 s.
2. Short vignette, `quarto render` in `/tmp/cran11/vig` (copy md5 `93acaa32…`): exit 0, **7.1 s wall
   clock**. The fit time it prints is **1.7 s**. No chunk failed, so no fix was needed.
3. Quoted verbatim from the rendered page (`pandoc -t plain`):

   The identified subgroup line:

   ```
   #> Identified subgroup: {er <= 0} & {pgr <= 26}
   #> Patients in it: 75 of 686. Fit: 1.7 s
   ```

   The funnel table:

   ```
     stage                                                   candidates
     ----------------------------------------------------- ------------
     enumerated (at most two cuts, size-eligible)                 2,211
     fittable (ten events in each arm, Cox fit obtained)          1,410
     clearing the screening threshold, HR > 1.0                     121
     passing the consistency screen, p_cons >= 0.90                   9
   ```

   The unadjusted table:

   ```
                             N   events unadjusted HR (95% CI)
     ------------------- ----- -------- ------------------------
     identified region      75       41 2.222 (1.176, 4.199)
     complement            611      258 0.605 (0.463, 0.791)
   ```

   **The numbers agree with the prose.** The candidate table and the printed lines give 75 for
   max N(ε) (`#> max N(eps) takes the largest candidate within the band: N 75, HR 2.22.`), 82 for max N
   and 61 for max β (`#> max N alone would take N 82 (HR 1.95); max beta the N 61.`). The cut is
   `{pgr <= 26}`. The 61 is `{er <= 0} & {size <= 35}`, HR 2.537, the largest in the family.
4. Moved article by the site-render route: package-path export of HEAD plus `_pkgdown.yml` and `.Rprofile`
   in `/tmp/cran11/site`, `R_LIBS=/tmp/cran11/lib Rscript -e 'pkgdown::init_site();
   pkgdown::build_article("articles/gbsg-table1")'`. Exit 0, **106 s wall clock** including `init_site()`
   (the task's record for the old vignette: 88 s). The log has no `Skipping`. The page
   `articles/gbsg-table1.html` has the pkgdown marker on line 2 and the navbar `Released version` 0.5.0,
   and it says "This article runs…". Its printed fit times: forest search with inference 39.1 s, DINA
   8.1 s, GRF 25.4 s.

## Step 7 — the release check

From `/tmp/cran11/build`, `R_LIBS=/tmp/cran11/lib` (`NOT_CRAN` and `OMP_THREAD_LIMIT` unset):
`R CMD build ~/Documents/GitHub/forestsearch`, then `R CMD check --as-cran forestsearch_0.5.0.tar.gz`.

- **Status: 1 NOTE.** 0 errors, 0 warnings. The NOTE, verbatim:

  ```
  * checking CRAN incoming feasibility ... [7s/16s] NOTE
  Maintainer: ‘Larry Leon <larry.leon.05@post.harvard.edu>’

  Found the following (possibly) invalid URLs:
    URL: https://larry-leon.github.io/forestsearch/articles/gbsg-table1.html
      From: inst/doc/forestsearch.html
            README.md
      Status: 404
      Message: Not Found
  ```

  The URL is the new article's page, which is not on the website yet. On 2026-11-14, `curl` gave 404 for
  it, and 200 for `articles/forestsearch.html`, which serves the Table 1 analysis with navbar 0.4.1. The
  site workflow `.github/workflows/pkgdown.yaml` runs on push to `master`. The URL appears in the short
  vignette's own text (twice, unchanged) and in `README.md` (Step 3.2 #3, #4).
- Tarball `forestsearch_0.5.0.tar.gz`: **1,638,428 bytes**, 749 entries. No `AGENTS.md`, no
  `CRAN-SUBMISSION`, no `vignettes/articles/`. Its `vignettes/forestsearch.qmd` md5 is `93acaa32…`, and its
  `tests/` is identical to HEAD `a4fff743`'s. `R CMD build` removed no directory this time.
- Nothing was written into the repository. No file outside `.git/` is newer than the check's start
  marker, and `git status` is unchanged.
- Check log: `dev/reports/REPORT_cran_resubmit_vignette_2026-11-14_check_log.txt` (copy of `00check.log`).

| Timing line | This check | Last check (`REPORT_cran_core_tests_2026-11-13.md`) |
|---|---|---|
| CRAN incoming feasibility | [7s/16s] NOTE | [6s/16s] OK |
| package can be installed | [20s/20s] | [20s/20s] |
| R code for possible problems | [46s/46s] | [47s/47s] |
| examples | [35s/35s] | [34s/34s] |
| examples with --run-donttest | [61s/61s] | [62s/62s] |
| tests | no time printed (under 10 s); `testthat.Rout` elapsed 4.486 s, `[ FAIL 0 \| WARN 0 \| SKIP 163 \| PASS 42 ]` | no time printed; elapsed 4.741 s, same counts |
| re-building of vignette outputs | **no time printed (under 10 s)**; the Step 5.2 render took 7.1 s | [78s/77s] |
| PDF version of manual | [13s/13s] | [13s/13s] |
| HTML version of manual | [18s/18s] | [17s/17s] |
| **check total** (shell `time`) | **4 min 29.1 s** | 5 min 39.0 s |
| build (shell `time`) | 36.2 s | 1 min 44.9 s |
| tarball | 1,638,428 bytes | 1,662,444 bytes |

## Step 6 — `cran-comments.md` as committed (`46c1fe50`)

The Resubmission section was inserted after the Submission section, before "## Test environments", with
the task's text and line breaks. Under the task's "nothing else changes", the sentence "The local check
above ended with `Status: OK`." stays as it was. It now disagrees with this check's `Status: 1 NOTE` and
with the counts line.

````markdown
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

0 errors | 0 warnings | 1 note

The local check above ended with `Status: OK`.

On CRAN a small core set of tests of the main function `forestsearch()` is run. The full
test suite is run locally before each release.

## Reverse dependencies

There are no reverse dependencies: CRAN's package index (read on
2026-11-12) lists no package that depends on, imports, links to or suggests
forestsearch.
````

## This session's commits

| Commit | Subject |
|---|---|
| `ca53c0e5` | dev: task CRAN resubmission, short vignette (as received) |
| `4997efe2` | dev: task CRAN resubmission, Amendment 1 (as received in chat) |
| `e2d447c5` | build: ignore AGENTS.md (Amendment 1) |
| `12332599` | docs: short vignette; the GBSG Table 1 analysis becomes an article (Steps 1 to 3) |
| `a4fff743` | test: core test compares the repeated call with a tolerance (Step 4) |
| `46c1fe50` | docs: cran-comments.md for the 0.5.0 resubmission (Step 6) |
| (this commit) | dev: this report and the check log (Step 8) |

Nothing was pushed, pulled, merged or tagged.

## Final checks (run with this report and the check log staged)

1. The session's commits, `e9efc7ff..HEAD` plus the staged report and check log (`--no-renames`), touch
   16 paths, all in the allowed list: `vignettes/forestsearch.qmd`, `vignettes/articles/gbsg-table1.qmd`,
   `_pkgdown.yml`, `NEWS.md`, `README.md`, the five articles changed under Step 3.2
   (`survival-post-selection.qmd`, `forestsearch-2024.Rmd`, `methodology.qmd`, `paper_simulations.qmd`,
   `actg175-glm-outcomes.qmd`), `tests/testthat/test-core-forestsearch.R`, `cran-comments.md`,
   `.Rbuildignore` (Amendment 1), the task file, this report and the check log. Paths outside the list:
   0. Paths under `R/`, `man/`, `inst/`, `DESCRIPTION` or `NAMESPACE`: 0. No file changed under Step 3.3.
2. `DESCRIPTION` says `Version: 0.5.0`: yes.
3. Modified tracked files: 0 (read under Amendment 1). Staged: only this report and the check log.
   Untracked: `AGENTS.md`, `CRAN-SUBMISSION` (Larry's, left untouched).
4. Nothing pushed, pulled, merged or tagged: `origin/master` is still `e9efc7ff`; no remote branch contains
   `ca53c0e5`; no tag points at a session commit; the reflog since `e9efc7ff` has only this session's 6
   commits (no merge, pull or push).

Step 8.4 follows this commit: delete `/tmp/cran11/lib`, keep `/tmp/cran11/build`.
