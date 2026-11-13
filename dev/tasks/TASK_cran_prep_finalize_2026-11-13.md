# TASK — CRAN preparation, finalize: last fix, test timing, full suite, release check (2026-11-13)

Session name: `fs-cran-prep-finalize`

Host `pop-os`, repository `~/Documents/GitHub/forestsearch`, branch `cran-prep`, version 0.5.0.
Larry said "finalize" on 2026-11-13. This is the one task in which the full test suite and
`R CMD check --as-cran` run. It replaces steps 5, 6 and the step 8 timing table of
`TASK_cran_prep_release_2026-11-12.md`, which were skipped. Read this whole file first.

## Rules

- **Lean run.** The main session does everything. No verifier, polish or report-drafting
  agents. One short report, plus the raw check log as its only companion file.
- **One pass that reports everything.** Only Step 0 has a stop gate. After it, a failure is
  recorded and the task goes on, so that one run shows every problem. The report's first line
  says READY or NOT READY.
- **Packages.** No install into the user or site library. The user library holds forestsearch
  0.4.1.9002, and test workers load whatever is installed, so every test run in this task puts
  a scratch build first: `R_LIBS=/tmp/cran8/lib`. Scratch: `/tmp/cran8`.
- **Commits** are path-limited; stage nothing else. One change, one commit.
- **No `R/` change other than Step 1.** No expected test value is edited anywhere.
- **No push, pull, merge or tag.** The version stays 0.5.0.
- fs-post-selection and fs-glms-interpretable are not touched. Keep `~/fs050_check/payloads/`.
  You may delete the leftover `/tmp/chk050/lib` and `/tmp/chk050/src`.

## Step 0 — start

1. Copy this file into `dev/tasks/` and commit it. Record the start SHA.
2. **Gate (the only stop):** branch `cran-prep`; `DESCRIPTION` says `Version: 0.5.0`;
   `git status --porcelain` is empty.
3. `git archive` the start SHA into `/tmp/cran8/start` (the "start export").
4. Run the six default-path cases on the start export, as parts 6 and 7 ran them
   (`REPORT_cran_prep_code_fixes_2026-11-12_gate_run.R` and `_gate_compare.R`, unchanged).

## Step 1 — `R/`, change of behaviour: the GRF row of count simulations (option B)

Background: `dev/reports/REPORT_cran_prep_policy_fixes_2026-11-12.md`, "Fix 3". The forest fit
in `grf.subg.harm.glm()` does not use `effect_measure`; for count outcomes the callee only
validates it and uses it for its own summary row. `.run_grf_analysis_gen()`
(`R/run_simulation_analysis.R:1340`) passes the simulation's measure unchanged, so `"IRR"`
and `"IRD"` are rejected and every count GRF row fails.

- **Larry's decision:** for a count outcome, `.run_grf_analysis_gen()` passes `"log_IRR"` to
  `grf.subg.harm.glm()`, whether the simulation's measure is `"IRR"` or `"IRD"`. Every other
  outcome type passes what it passes now. Nothing inside `grf.subg.harm.glm()` changes.
- **Check before editing:** confirm by reading that no value the simulation reports for a GRF
  row (`hr.*`, `ahr.*`, classification columns) is taken from the callee's own summary row;
  the simulation's extractor must compute them. If one is taken from that row, **do not make
  this fix**: record the line and go to Step 2.
- Add a short source comment at the call: the forest fit does not use the value, and the
  simulation's extractor computes the reported effect on the simulation's own scale.
- Do not edit roxygen comments. No `man/` change.
- **Test** (new file `tests/testthat/test-sim-grf-count-row.R`; fails before, passes after):
  - a small count DGM with an exposure offset, one replicate each for `"IRR"` and `"IRD"`: the
    GRF analysis completes, with no error or warning matching `effect_measure must be`
    (before the fix, record the message);
  - 2 threads; `skip_if_not_installed("grf")`; `skip_on_cran()` if the file takes over
    5 seconds (measure it).
  - Also run one count replicate **without** an offset column and record what happens
    (completes, or the message it gives). This is a record, not part of the test.
- **Unaffected case:** `run_simulation_analysis()` with 3 replicates each for binary (OR),
  continuous (MD) and count (IRR), `seed_base` 8316951. Every row is `identical()` between
  the start export and the fixed tree (timing fields dropped) except the count GRF rows.
- **Also:** the five existing extractor test files and `test-sim-effect-scales.R` pass; the six
  default-path cases are identical to Step 0.
- **NEWS:** one bullet under `# forestsearch 0.5.0`, "Bug fixes". No other NEWS line changes.
- One commit (fix, test, NEWS bullet). If any of the above fails: revert the uncommitted
  edits, record why, and go on.

## Step 2 — tests only: time the suite as CRAN runs it, and skip slow files on CRAN

1. `git archive` HEAD into `/tmp/cran8/src` and `R CMD INSTALL` it into `/tmp/cran8/lib`.
   Confirm that the main process and one worker load forestsearch 0.5.0 from there.
2. Run the whole suite once **in CRAN mode**: `testthat::on_cran()` must be `TRUE` inside the
   run (`devtools::test()` sets `NOT_CRAN=true`, so do not use it here), with
   `_R_CHECK_LIMIT_CORES_=TRUE`. Record the evidence that the run was in CRAN mode, the total
   wall clock, and the wall clock per test file.
3. **Larry's rule:** every test file that takes more than 30 seconds in CRAN mode is skipped
   on CRAN. Add `skip_on_cran()` in the form that skips the whole file under testthat 3.3.2.
   Nothing else in those files changes.
4. Re-run only the changed files in CRAN mode and confirm each is skipped and no longer takes
   the time. Record the new CRAN-mode total (measured or by subtraction; say which).
5. One commit for all the skips. If no file is over 30 seconds, make no commit and say so.
6. In the report: the ten slowest files before the rule, the files skipped, and the total
   before and after.

## Step 3 — the full suite

- Run the full suite once with `NOT_CRAN=true` and `R_LIBS=/tmp/cran8/lib` (rebuild the
  scratch library first only if `R/` changed after Step 2.1).
- Record blocks, expectations, failures, errors, skips, warnings and the wall clock.
  Baseline: part 5 had 736 blocks, 9,242 expectations, 0 failures; parts 6 and 7 and Step 1
  added test files since.
- For each failure or error: the file, the test name and the message. Do not edit code or
  expected values. Go on to Step 4 either way.

## Step 4 — the release check

1. From `/tmp/cran8/build`, run `R CMD build` on the repository, then
   `R CMD check --as-cran` on the tarball. Do not set `NOT_CRAN` or `OMP_THREAD_LIMIT`.
   Nothing is written into the repository (no tarball, no `.Rcheck` directory).
2. Record: errors, warnings, and **every NOTE verbatim**; the tarball size and the installed
   size; the total check time.
3. **Timing table** (from the check's own output): examples total and every example over
   5 seconds; tests total; vignette re-build time; total.
4. **If a NOTE reports CPU time well above elapsed time** (threads): run the check once more
   with `OMP_THREAD_LIMIT=2` and report both results. Change no file for this.
5. **Repairs allowed without asking,** each as its own commit, followed by one re-run of the
   check:
   - an `.Rbuildignore` line for a stray top-level file or directory that the check names;
   - a URL the check reports as permanently moved, in a file that is part of the build (text
     only).
   Anything else is recorded, not repaired.
6. Copy `00check.log` to
   `dev/reports/REPORT_cran_prep_finalize_2026-11-13_check_log.txt`.

## Step 5 — `cran-comments.md`

Rewrite only the sections "Test environments" and "R CMD check results":

- **Test environments:** the local line (OS, R version, `R CMD check --as-cran` on the built
  tarball). Then two lines left for Larry, exactly: `- win-builder (devel and release): TO BE
  FILLED BY MAINTAINER` and `- R-hub: TO BE FILLED BY MAINTAINER`.
- **R CMD check results:** the counts of errors, warnings and NOTEs from the final check run,
  and each NOTE verbatim with one factual line under it. No claim that is not in the check
  output.
- The "Submission" and "Reverse dependencies" sections stay as they are.
- One commit.

## Step 6 — report and close

1. Write `dev/reports/REPORT_cran_prep_finalize_2026-11-13.md`. Short:
   - **first line:** READY or NOT READY. READY means: the final check has 0 errors and
     0 warnings; the full suite has 0 failures and 0 errors; no CRAN-run test file is over
     30 seconds. If NOT READY, the reasons in bullets;
   - **top table:** one row per step — result, commit SHA (if any), the key numbers;
   - every NOTE verbatim;
   - the timing table;
   - the Step 2 file table;
   - Step 1: the red and green counts, the unaffected case, and the no-offset record;
   - the start SHA and this session's commits;
   - the full text of `cran-comments.md` as committed;
   - anything left open, in bullets. Facts only.
2. Commit the report and the check log (path-limited). Copy the report to `~/Downloads`.
3. **Final checks** (report the result; these do not undo anything):
   - this session's commits (list the SHAs) touch only: `R/run_simulation_analysis.R`,
     `tests/testthat/` files named in Steps 1 and 2, `NEWS.md`, `cran-comments.md`,
     `.Rbuildignore` or a URL file if Step 4.5 was used, this task file, the report and the
     check log;
   - `DESCRIPTION` still says `Version: 0.5.0`; no file under `man/` changed;
   - `git status --porcelain` is empty;
   - nothing was pushed, merged or tagged.
4. Delete `/tmp/cran8/lib` and `/tmp/cran8/src`. Keep `/tmp/cran8/build`.
5. Print the first line, the top table and every NOTE as the closing message.

## Amendment 1 (Larry, 2026-11-13, in chat during Step 4)

Amendment 1 (Larry, 2026-11-13): in the Step 4.5 repair, also add an .Rbuildignore line that excludes tests/testthat/Rplots.pdf, in the same commit as the legacy_vignette line. Then re-run the build and check once as Step 4.5 says, and report the tarball size before and after and confirm that neither stray is in the new tarball. Append this amendment to the task file as its own commit. Everything else in the task stands.
