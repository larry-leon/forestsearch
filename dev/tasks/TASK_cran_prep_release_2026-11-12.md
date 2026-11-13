# TASK — CRAN preparation, part 6: simulation effect scales, release 0.5.0 (2026-11-12)

Complete on its own. Run from the forestsearch repo root (`~/Documents/GitHub/forestsearch`).
Authorized by Larry on 2026-11-12, after the part 5 report. This task runs unattended and
never stops to ask.

It continues on branch `cran-prep`. Larry has merged `cran-prep` into `master` through part 5
(version 0.4.1.9011).

## What this task does

| Part | Item | Kind |
|---|---|---|
| 1 | `run_simulation_analysis()` reports risk difference, rate difference and risk ratio on the scale each label says | Changes behaviour, for those three measures only |
| 2 | Which committed simulations used those three measures | Read-only |
| 3 | Release preparation for version 0.5.0 | Version, NEWS, `cran-comments.md`; checks; no code change |

**Left to Larry, not done here:** uploading to win-builder or R-hub, the CRAN submission,
tagging, rebuilding the pkgdown site, re-rendering any article.

## Rules

- No `git fetch`, `pull`, `push`, `merge` or `tag`. Every `git add` names explicit paths.
- Nothing is uploaded anywhere. Reading CRAN's package index and package sources is allowed.
- No install into the user or site library. Temporary libraries and scratch go under `/tmp`.
- fs-post-selection and fs-glms-interpretable: read-only, no branch change.
- `vignettes/` and `quarto/`: read-only. Nothing is re-rendered. No simulation campaign is run.
- No new package dependency. Do not install a tool to run a check; if a tool is missing,
  report that check as NOT RUN.
- Evidence is the current source. Quote file:line in the report.
- Gates are stop-on-failure for their own part, never stop-to-ask. FAILED: the part's files
  are restored and the output reported. BLOCKED: a pre-set condition was not met and nothing
  was changed. The report is always written.
- Do not adjust a seed, a tolerance or an expected value to make a test pass.
- Reproducers: synthetic data, sequential plan, `pkgload::load_all()` from the tree, cap
  2 minutes each.
- No timing estimates. Report measured wall clocks.

## Step 0 — branch and record

1. `git status --porcelain` is empty. **Gate.**
2. If branch `cran-prep` exists, check it out. Otherwise create it from `master`.
3. **Gate:** `git diff --quiet master cran-prep -- R man tests DESCRIPTION NAMESPACE NEWS.md vignettes _pkgdown.yml`
   exits 0, and `dev/reports/REPORT_cran_prep_displays_se_docs_2026-11-12.md` exists. If not,
   stop: BLOCKED.
4. Copy `~/Downloads/TASK_cran_prep_release_2026-11-12.md` to `dev/tasks/` and commit it by
   explicit path (`dev: add CRAN preparation task, part 6`). This commit is the **start SHA**.
   Export it with `git archive` (package paths only) to `/tmp`.
5. Record: start SHA, `DESCRIPTION` version (expected 0.4.1.9011), each sibling repository's
   branch, SHA and `git status --porcelain`.

## Part 1 — simulation effect scales (`R/run_simulation_analysis.R`)

**Found by part 5's sweep** (its report §4, A4): in the forest search extractor's
`.glm_effect()` (near 995–1022) and its twin in the GRF extractor (near 1392–1419), the value
returned for three measures is on another scale:

| Label | What is returned today | What the truth columns (`compute_aor()`) use |
|---|---|---|
| `"RD"` | the logistic coefficient (log odds ratio) | `m1 - m0` |
| `"IRD"` | the Poisson coefficient (log rate ratio) | `m1 - m0` |
| `"RR"` | the exponentiated logistic coefficient (odds ratio) | `m1 / m0` |

1. **Confirm at source** and quote the lines: what each label returns today, what
   `compute_aor()` returns for it, and how the rate ratio handles exposure (offset or not).
2. **List every consumer** of the affected columns (`hr.itt`, `hr.H.hat`, `hr.Hc.hat`,
   `hr.H.true`, `hr.Hc.true` and any others the two functions fill) under `R/`: the summary,
   table and plot functions, and what each does with the value for these three measures.
   **If any consumer applies a transformation that would be wrong for the corrected value,
   stop this part: BLOCKED, with the list.** Do not change a consumer.
3. **Change, in both extractors:** from the same two-arm fit the function already makes,
   return
   - `"RD"`: `p1 - p0`, the two fitted arm proportions;
   - `"RR"`: `p1 / p0`;
   - `"IRD"`: `r1 - r0`, the two fitted arm rates, on the exposure convention the rate ratio
     already uses, so that `r1 / r0` equals the value returned for `"IRR"`.
   A value that is not finite is returned as the function's existing failure value.
4. **Not changed:** odds ratio, rate ratio, mean difference and survival; the direction
   (treatment minus control, the simulation layer's convention); the truth columns.
5. **Test** (new file; written first; fails before, passes after): on hand-built two-arm data
   with known counts, each of the three measures equals the value worked out by hand from the
   arm proportions or rates, in both extractors' functions. State the working in comments.
6. **Unaffected case:** for odds ratio, rate ratio, mean difference and survival, the function
   values, and a `run_simulation_analysis()` run of a few replicates per outcome type with a
   fixed seed, are identical to the start SHA.
7. **Scope:** the commit touches `R/run_simulation_analysis.R`, its `man/` pages if roxygen
   changes, and the new test file. Correct any roxygen that describes these columns.
8. Commit by explicit path
   (`fix: run_simulation_analysis() reports RD, IRD and RR on their own scales`).

## Part 2 — which simulations used the three measures (read-only)

Search forestsearch `quarto/` (mark archive and legacy directories), `vignettes/`, `tests/`,
and, read-only, both sibling repositories, for every `run_simulation_analysis()` run (direct,
or through a wrapper that calls it) on a non-survival mechanism whose effect measure is a
risk difference, a rate difference or a risk ratio. Say how the measure is set at each site.

| Repo | File:line | Measure | Current or archive | Outputs that show the affected columns (table, figure) | Tracked payloads that store them (count) |
|---|---|---|---|---|---|

State the totals, and state plainly if there are none. Change nothing.

## Part 3 — release preparation for 0.5.0

Part 3 runs whatever Part 1's status is. The report's first lines state Part 1's status so
that Larry can decide whether to submit.

1. **Version.** `DESCRIPTION`: `Version: 0.5.0`, and the `Date` field if there is one. List,
   without editing, every other place a version string appears (`inst/CITATION`, `README`,
   `_pkgdown.yml`, `cran-comments.md`).
2. **NEWS.** Put every bullet that sits above the `# forestsearch 0.4.0` heading under one
   `# forestsearch 0.5.0` heading, grouped under: Breaking changes; Bug fixes; New features;
   Superseded and deprecated; Documentation. Move bullets verbatim; merge two only if they are
   identical. Add a short opening summary (at most eight lines) written from the bullets,
   naming the two removed functions and the main changes in results a user could see.
   **Gate:** a script shows the multiset of bullets before equals the multiset after (identical
   duplicates aside), with the list of any merged. Save it beside the report.
3. **`cran-comments.md`** (create or update), truthful and short:
   - the check result from step 5, with R version and platform;
   - test environments: list only those actually run in this task; say that win-builder and
     R-hub have not been run;
   - the reverse-dependency result from step 4;
   - one sentence: this release removes `create_gbsg_dgm()` and `simulate_from_gbsg_dgm()`,
     which were deprecated in 0.4.0;
   - each NOTE from step 5, with its explanation.
   Confirm that the file is excluded from the build (`.Rbuildignore`).
4. **Reverse dependencies.** From CRAN's index, list the packages that depend on, import,
   link to or suggest forestsearch. For each, download its source to `/tmp` and search it for
   `create_gbsg_dgm`, `simulate_from_gbsg_dgm` and `get_dgm_with_output`. Report the list and
   any use. Do not run their checks. If the index cannot be reached, report this step BLOCKED.
5. **Release check.** Commit steps 1 to 3 first (`chore: release 0.5.0`), then export the
   head, build the tarball with `R CMD build` in `/tmp`, and run `R CMD check --as-cran` on it
   (cap 60 minutes). **Gate:** 0 errors, 0 warnings.
   - Quote every NOTE in full and class it: CRAN incoming feasibility, or other.
   - A NOTE of another kind may be repaired only if its cause is in `DESCRIPTION`, `NEWS.md`
     or `cran-comments.md`; otherwise list it as NEEDS ATTENTION.
   - Report: tarball size; the wall clock of the whole check and of its tests, examples and
     vignette steps; whether the tests ran in CRAN mode (`NOT_CRAN` unset).
6. **Full suite** (`devtools::test()`, `NOT_CRAN=true`), cap 30 minutes. **Gate:** 0 failures,
   0 errors. Report blocks and expectations beside part 5's (736 blocks, 9,242 expectations).
7. **Default-path gate.** Re-run the six equivalence cases with part 4's scripts at the head
   against the start-SHA export. All six identical.
8. **CRAN policy audit (read-only; list, do not edit).** One table each:
   - **Cores.** Every example, every test that runs on CRAN, and the vignette CRAN builds,
     that can start more than two workers or cores: file:line and how many.
   - **Time.** The ten slowest test files in CRAN mode and the five slowest examples.
   - **Side effects.** Any example, test or vignette that writes outside `tempdir()`, changes
     `options()` or `par()` without restoring them, or installs a package.
   - **Spelling and links.** Run the spelling and link checkers if they are installed; list
     the findings. Fix a misspelling only in a file this task already edits.
9. Commit any repair from step 5 by explicit path (`chore: release 0.5.0, check repairs`).

## Report and close

1. Write `dev/reports/REPORT_cran_prep_release_2026-11-12.md`, in bullets and tables:
   - first lines: the status of Parts 1, 2 and 3, the version, and whether the package is
     ready to submit as far as this task can tell;
   - provenance and commits;
   - Part 1: sites before and after, the consumer list, the red and green runs, the hand
     working, the unaffected case, one before/after example; or the reason if BLOCKED;
   - Part 2: the table and totals;
   - Part 3: the NEWS gate, the full text of `cran-comments.md`, the reverse-dependency list,
     every NOTE, the timings, the suite, the default-path gate, the audit tables;
   - **Larry's remaining steps**, in order: win-builder and R-hub, submission, tag, merge,
     pkgdown rebuild;
   - OPEN ITEMS, limited to what this task leaves open.
2. Commit the report and its scripts by explicit path (`dev: CRAN preparation report, part 6`)
   and copy the report to `~/Downloads`.
3. **Final gates:** `git status --porcelain` is empty; both sibling repositories are
   unchanged; `git diff <start> --stat` lists only files the report names, and nothing under
   `vignettes/` or `quarto/`. Nothing is pushed, merged, tagged or uploaded.
