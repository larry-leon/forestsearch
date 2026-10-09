# TASK — CRAN release, Stage 3: win-builder follow-up (forestsearch 0.4.0)

- Date: 2026-10-09
- Repo: `larry-leon/forestsearch`, branch `feature/glm-extension`, at or after the Stage 2 report commit.
- Gates are stop-on-failure.

## Why

win-builder (R-release 4.6.1 and R-devel r90650) both returned `Status: 1 NOTE`, and their tests took **11 minutes** on Windows (`testthat.Rout`: elapsed 674 s; whole check ~20 min). The same suite takes 133 s on the Mac.

1. **The NOTE (both platforms), verbatim:**
   ```
   * checking HTML version of manual ... NOTE
   Found the following HTML validation problems:
   forestsearch.html:28329:58 (inject_mrct_structure.Rd:69): Error: <region> is not recognized!
   forestsearch.html:28329:58 (inject_mrct_structure.Rd:69): Warning: discarding unexpected <region>
   ```
   The Mac check missed it because the local HTML Tidy is too old.
2. **Check time:** CRAN asks that a package check take well under 10 minutes. An 11-minute test run on Windows is a likely reviewer rejection.

## R/ callout

- `R/mrct_inject_dgm.R`: roxygen text only (the `<region>` string) — documentation, no code line changes.
- `tests/testthat/*`: adding `skip_on_cran()` only — no test logic changed, no test removed. Tests still run in full locally (`devtools::test()` sets `NOT_CRAN=true`).
- Nothing changes behaviour or the method.

## Authorizations (Larry's kickoff of this file)

- One local `R CMD build` + `R CMD check --as-cran` (Step 4), with the temp-library setup used in Stage 1 (CRAN `weightedsurv` and `latentcor` into a scratch library only; binaries allowed for latentcor's Fortran dependencies).
- One upload: `devtools::check_win_devel(pkg = ".")` (Step 5). No other uploads; no CRAN submission.
- Not authorized: `git fetch`/`pull`/`push`, tags, installs into the user library, any other R/ edit.

## First action

1. Copy this file to `dev/tasks/TASK_cran_release_wb_followup_2026-10-09.md` and commit it alone.

## Preconditions (stop on failure)

- P1. `git status --porcelain` shows no tracked-file modifications.
- P2. `Version: 0.4.0`; HEAD contains `dev/reports/REPORT_cran_release_remote_checks_2026-10-09.md`. Record HEAD SHA.
- P3. `pandoc`, `quarto`, `pdflatex` on `PATH`.

## Steps

### Step 1 — HTML NOTE (≤ 15 min)
- Read `man/inject_mrct_structure.Rd` line 69 and the roxygen source in `R/mrct_inject_dgm.R` that generates it. Record both lines verbatim.
- Make the literal `<region>` render as text: wrap it in backticks (`` `<region>` ``) in the roxygen source; if it is already inside a code span or a different construct is responsible, use the minimal edit that makes Rd emit an escaped angle bracket, and record why.
- `roxygen2::roxygenise()`.
- Search all of `man/*.Rd` for any other raw `<word>` that would reach HTML unescaped (outside `\code{}`, `\preformatted{}`, `\verb{}`, examples, and `\url`/`\href`); list hits and fix them the same way.
- Gate 1: `tools::Rd2HTML("man/inject_mrct_structure.Rd", out = <scratch file>)` output contains `&lt;region&gt;` and no raw `<region>`; NAMESPACE byte-identical; in `R/`, every diff line matches `^[+-]\s*#'`.

### Step 2 — Time the test files (≤ 20 min, no edits)
- For each `tests/testthat/test-*.R`, time it in a fresh R process with `NOT_CRAN` unset (so existing `skip_on_cran()` calls apply): `Rscript -e 'devtools::load_all(quiet=TRUE); t <- system.time(testthat::test_file("<f>", reporter = "silent")); cat(t[["elapsed"]])'`. Record elapsed per file, sorted descending, and the total.

### Step 3 — Skip the heavy files on CRAN (≤ 15 min)
- Target: the CRAN-run test time on the Mac ≤ **30 s** (Windows ran ~5× slower: 133 s → 674 s, so this keeps Windows near 2–3 min).
- Going down the Step 2 ranking, add a file-level `skip_on_cran()` as the first executable line of each file until the remaining (non-skipped) total is ≤ 30 s. Never edit test bodies.
- Do not skip the two Stage 1 files `test-save-leff-calibration.R` and `test-check-thread-cap.R` (both fast).
- Record the skipped files and the projected remaining total.
- Gate 3: re-run Step 2's timing over all files: non-skipped total ≤ 30 s; and `devtools::test()` (full suite, `NOT_CRAN=true`) gives `FAIL 0` — this is the only full-suite run, authorized here.

### Step 4 — Local release check (build ≤ 45 min, check ≤ 75 min)
- Same commands and environment as Stage 1 Step 9 (temp library; `R_LIBS`, `_R_CHECK_CRAN_INCOMING_REMOTE_=FALSE`, `_R_CHECK_LIMIT_CORES_=TRUE`, `OMP_THREAD_LIMIT=2`, both timing-ratio thresholds 2.5).
- Gate 4: 0 ERRORs, 0 WARNINGs; record the status line, the tests elapsed time, and NOTEs verbatim. Copy `00check.log` to `dev/reports/cran_wbfollowup_00check_2026-10-09.log`. Delete temp library and scratch dirs.

### Step 5 — Commit, then one win-builder R-devel upload (≤ 10 min)
- Commit Steps 1 and 3 (explicit paths; two commits: docs; test skips) and the check log **before** uploading.
- `devtools::check_win_devel(pkg = ".")` once. Record console output verbatim. Do not wait for the email.

### Step 6 — cran-comments.md
Rewrite the results sections:
- `## Test environments`: local macOS aarch64 (R 4.5.2), win-builder R-release 4.6.1 and R-devel (r90650), first round; win-builder R-devel second round pending (after the fixes); mac-builder unavailable (HTTP 502 on two attempts, 2026-10-09).
- `## R CMD check results`: `0 errors | 0 warnings | 0 notes` expected after the fixes; state that the first win-builder round's single NOTE (HTML validation, `inject_mrct_structure.Rd`) is fixed in this version; local check status from Step 4 verbatim.
- One line: long-running tests are skipped on CRAN with `skip_on_cran()`; the full suite runs locally.
- Keep `## Submission` and `## Reverse dependencies` as they are.
- Commit.

## Deliverable

`dev/reports/REPORT_cran_release_wb_followup_2026-10-09.md` (bullets): the two Rd lines and fix; Step 2 timing table; skipped files and remaining total; Step 4 status and tests time; Step 5 console output; cran-comments summary. Commit it, then copy it to `~/Downloads`.

## Post-conditions

- PC1. Changed files since the task-doc commit: this task doc, the report, the check log, `cran-comments.md`, `R/mrct_inject_dgm.R` (and any other R/ file whose roxygen produced a Step 1 hit), the regenerated `man/*.Rd`, and `tests/testthat/test-*.R` files whose only change is the added `skip_on_cran()` line.
- PC2. `NAMESPACE` byte-identical to the Stage 2 HEAD; `Version: 0.4.0` unchanged.
- PC3. Every changed test file's diff is exactly one added line matching `^\+\s*(testthat::)?skip_on_cran\(\)\s*$`.
- PC4. `git status --porcelain` clean; no tarball, `.Rcheck` or temp library in the repo.
