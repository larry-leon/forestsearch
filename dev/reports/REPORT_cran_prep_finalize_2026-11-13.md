READY

# REPORT — CRAN preparation, finalize: last fix, test timing, full suite, release check (2026-11-13)

Task: `dev/tasks/TASK_cran_prep_finalize_2026-11-13.md`, with Amendment 1 (Larry, in chat during
Step 4). Session `fs-cran-prep-finalize`, host `pop-os`, branch `cran-prep`, version 0.5.0.
R 4.6.2 (2026-10-29), testthat 3.3.2, Pop!_OS 24.04 LTS. No agents were used.

READY because: the final check (run 2) has 0 errors, 0 warnings and 0 NOTEs; the full suite has
0 failures and 0 errors; no CRAN-run test file is over 30 seconds (the slowest left is 21.4 s).

| Step | Result | Commit | Key numbers |
|---|---|---|---|
| 0 Start | gate PASS | `d57378f4` (task file) | branch `cran-prep`, `Version: 0.5.0`, `git status --porcelain` empty; six default-path cases on the start export, 17.2 s |
| 1 Count GRF row | PASS, with one reading of the unaffected case (§1) | `87688d39` | red 1 test / 8 expectations, 4 failed; green 8 / 8, 1.57 s; extractor files 274 / 274; `test-sim-effect-scales.R` 54 / 54; six cases identical to Step 0 |
| 2 CRAN-mode timing | 2 files over 30 s, now skipped on CRAN | `b1394deb` | CRAN-mode total 208.4 s before, 103.1 s after (by subtraction) |
| 3 Full suite, `NOT_CRAN=true` | PASS | none | 112 files, 743 blocks, 9,310 expectations, 0 failures, 0 errors, 3 skips, 38 warnings; 689.8 s |
| 4 Check, run 1 | 0 errors, 0 warnings, 1 NOTE | none | tarball 8,531,649 bytes; installed 4,760 KB; check 7 min 10.1 s |
| 4.5 Repair (+ Amendment 1) | two `.Rbuildignore` lines | `844926db` (amendment), `f48c8f5f` (repair) | `^legacy_vignette$`, `^tests/testthat/Rplots\.pdf$` |
| 4 Check, run 2 (final) | **Status: OK**: 0 errors, 0 warnings, 0 NOTEs | none | tarball 1,660,903 bytes; installed 4,760 KB; check 6 min 56.4 s |
| 5 `cran-comments.md` | rewritten (two sections) | `9e9a7a73` | 0 errors, 0 warnings, 0 notes |
| 6 Report | this file + check log | (this commit) | |

Step 4.4 did not apply: no NOTE reports CPU time above elapsed time, so the `OMP_THREAD_LIMIT=2`
run was not made.

## NOTEs (verbatim)

Run 1:

```
* checking top-level files ... NOTE
Non-standard file/directory found at top level:
  ‘legacy_vignette’
```

Run 2 (final): none. Run 1's and run 2's `00check.log` differ only in the current time, the
timings and this NOTE. The companion file
`REPORT_cran_prep_finalize_2026-11-13_check_log.txt` is run 2's `00check.log`.

## Timing table (release check, run 2; run 1 in brackets)

| Item | CPU / elapsed |
|---|---|
| Examples (default pass) | 33 s / 33 s (34 / 34) |
| Examples over 5 s, default pass | none; largest `fs_to_guohe` 0.318 s; the 222 examples sum to 2.03 s |
| Examples with `--run-donttest` | 58 s / 59 s (59 / 60) |
| Examples over 5 s, `--run-donttest` pass | `mr_method_tables`: 9.373 user + 0.051 sys, 9.418 elapsed (its `\donttest{}` block) |
| Tests (`testthat.R`) | 93 s / 103 s (97 / 108) |
| Vignette re-build | 74 s / 73 s (74 / 73) |
| Total `R CMD check` (shell `time`) | 6 min 56.4 s (7 min 10.1 s) |
| `R CMD build` (shell `time`) | 1 min 41.0 s (1 min 47.0 s) |

- The check's example-timing NOTE (threshold 5 s) is computed on the default pass only.
- The `--run-donttest` pass overwrites `forestsearch-Ex.timings`, so the default pass's copy was
  saved during run 2 by a file watcher.
- Other timed lines (run 2): incoming feasibility 6 s / 14 s, install 19 s / 19 s, R code
  42 s / 42 s, PDF manual 12 s / 12 s, HTML manual 17 s / 17 s.

## Step 2 — CRAN-mode timing

- Scratch library: `git archive` of `87688d39` into `/tmp/cran8/src`, then `R CMD INSTALL` into
  `/tmp/cran8/lib`. With `R_LIBS=/tmp/cran8/lib`, the main process and one `multisession` worker
  both loaded forestsearch 0.5.0 from `/tmp/cran8/lib/forestsearch`. `R/` did not change after
  this, so Step 3 used the same library.
- Run: a driver that does what `tests/testthat.R` → `test_check()` does: `library(forestsearch)`
  from the scratch library, `TESTTHAT_IS_CHECKING=true`, `test_dir(load_package = "installed")`.
  It adds a per-file wall-clock reporter. `NOT_CRAN` was unset, `_R_CHECK_LIMIT_CORES_=TRUE`, and
  the run was from `/tmp/cran8/src/tests`.
- CRAN-mode evidence, recorded by the reporter at the run's start and end: `on_cran()` TRUE,
  `NOT_CRAN` unset, `interactive()` FALSE. Also, 115 skips gave the reason "On CRAN"; 12 of them
  were file-level.
- testthat 3.3.2 does not export `on_cran()`, so it was called as `testthat:::on_cran()`.
- Before the rule: 112 files, 617 blocks, 4,957 expectations, 0 failed, 0 errors, 118 skipped
  blocks, 18 warnings.
- Total wall clock: 208.4 s from the driver, 210.0 s from the shell.

Ten slowest files before the rule (CRAN mode):

| File | Seconds | Blocks | Expectations | Rule |
|---|---|---|---|---|
| `test-mrct-case-study.R` | 74.8 | 7 | 64 | skipped on CRAN |
| `test-mrct-region-sims-mr.R` | 30.5 | 4 | 71 | skipped on CRAN |
| `test-mrct-region-sims-record.R` | 21.4 | 6 | 172 | |
| `test-mr-estimates-table-nullboot.R` | 9.8 | 4 | 33 | |
| `test-betaHhat-exact.R` | 8.7 | 6 | 38 | |
| `test-ahr-analytic.R` | 8.3 | 8 | 208 | |
| `test-dfbeta-lm-influence-path.R` | 4.3 | 5 | 31 | |
| `test-cross-outcome-parity.R` | 3.3 | 9 | 38 | |
| `test-threshold-remap-messages.R` | 3.3 | 3 | 27 | |
| `test-missing-membership-messages.R` | 3.0 | 10 | 60 | |

- Form of the skip: a top-level `skip_on_cran()` after each file's header comment. Nothing else in
  either file changed (+2 lines each). A scratch file confirmed this form under testthat 3.3.2: in
  CRAN mode no test block runs, and with `NOT_CRAN=true` everything runs.
- Re-run of the two files in CRAN mode: both were skipped at file level ("On CRAN"), in 0.007 s
  and 0.015 s.
- New CRAN-mode total, by subtraction: 208.4 − 74.8 − 30.5 + 0.02 = 103.1 s. The check's own
  test step in run 2 measured 93 s CPU / 103 s elapsed.
- `test-mrct-region-sims-mr.R` was 0.5 s over the 30-second line.

## Step 1 — count GRF row

**Check before editing (passed).**
- `.extract_grf_estimates_gen()` computes every value the simulation reports for a GRF row:
  `any.H`, `size.*`, `hr.*`, `ahr.*`, `cde.*`, `sens`, `spec`, `ppv`, `npv`.
- Its inputs are `grf_est$sg.harm.id`, `grf_est$data` (`treat.recommend` and the outcome columns),
  the simulated data, and the simulation's own `effect_measure`.
- `sg.def` and `sg.n_factors` come from `sg.harm.id`.
- No value is read from the callee's `effect_measure` or from a callee summary row.
- `grf.subg.harm.glm()` uses `effect_measure` only to validate it, in the no-offset warning and a
  verbose message, and to store it in its result (`R/grf_subg_harm_glm.R:325-352, 385`).

**The fix.**
- Edit at `R/run_simulation_analysis.R:1340`: for `outcome_type = "count"`, pass `"log_IRR"`.
  Otherwise the call is as before.
- A count with no measure already reached the callee's default `"log_IRR"`, so that case is
  unchanged.
- The only callee object the simulation keeps is the `grf_full` keeper (when `keep` asks for it).
  For a count, its `$effect_measure` is now `"log_IRR"`.

**Red and green.**
- Red: the start code with the new test gave 1 test, 8 expectations, 4 failed. The message before
  the fix: `GRF failed: For outcome_type = 'count', effect_measure must be 'log_IRR' or 'log_RR'.
  Got: 'IRR'.` (and `'IRD'`).
- Green: 8 / 8, 1.57 s (0.65 s in the CRAN-mode suite), so the file has no `skip_on_cran()`.

**No-offset record (count IRR, no `offset.name`).**
- Before the fix: `GRF failed: … Got: 'IRR'.`
- After the fix: the run completes. GRF identified 574 of 600 patients,
  `(age <= 40 & biomarker_hi > 0) | (age > 40 & age <= 84)`.
- It also gives the callee's warning: `outcome_type = 'count' with effect_measure = 'log_IRR' but
  offset.name = NULL. Fitting Poisson GLM without offset -- effects will be on the log-rate (not
  log-IRR) scale. Supply offset.name = '<person_time_column>' for incidence rate ratios.`

**Unaffected case.**
- Setup: `run_simulation_analysis()` with FS, FSlg and GRF, 3 replicates each, `seed_base`
  8316951, `n_sample` 600, on the start export and the fixed tree, with `fs_minutes` dropped.
- Binary (OR) and continuous (MD): 18 / 18 rows `identical()`.
- Count (IRR): the 3 GRF rows differ, as expected. Before the fix they were failed rows; after, all
  three identify a subgroup.
- The 6 count FS and FSlg rows differ only in `agree.FS.GRF`, `kappa.FS.GRF`, `agree.FSlg.GRF` and
  `kappa.FSlg.GRF`. These are per-replicate columns written on every row
  (`R/run_simulation_analysis.R:732-733`) that pair a method with GRF. Before the fix
  `kappa.FS.GRF` was 0 in all three replicates; after, it is 0.589, 0.568 and 0.027.
- Every other column of those rows is identical, including `agree.FS.FSlg` and `kappa.FS.FSlg`.
- As written ("every row identical except the count GRF rows"), the gate is not met for those 6
  rows. I kept the fix because the differing columns are GRF output (open item).

**Other checks.**
- The five existing extractor test files are the files that name the extractors at part 6's start
  (`9657afd4`): `test-build_classification_table.R` 13, `test-fs_mr_oc_summary.R` 57,
  `test-missing-membership-messages.R` 60, `test-nofind-classification.R` 96,
  `test-npv-empty-complement.R` 48. That is 274 / 274, part 6's figure.
  `test-sim-effect-scales.R` passed 54 / 54.
- Six default-path cases on the fixed tree: identical to Step 0 after dropping timing fields
  (PASS, 17.1 s).
- NEWS: one bullet, first under `# forestsearch 0.5.0` → `### Bug fixes`. The diff is 11 added
  lines and 0 deleted.

## Step 3 — full suite

- Route as in part 5: `devtools::test()` in the repository at `b1394deb`, with `NOT_CRAN=true` and
  `R_LIBS=/tmp/cran8/lib`. `on_cran()` was FALSE.
- Result: 112 files, 743 blocks, 9,310 expectations, 0 failures, 0 errors, 3 skips, 38 warnings.
  Wall clock 689.8 s (11 min 31.7 s from the shell).
- Against part 5 (109 files, 736 blocks, 9,242 expectations): +3 files, +7 blocks,
  +68 expectations. These are `test-sim-effect-scales.R` (4 / 54, part 6),
  `test-mrct-handlers-restored.R` (2 / 6, part 7) and `test-sim-grf-count-row.R` (1 / 8, Step 1).
- The 3 skipped blocks are part 5's three. The warnings are 38, as in part 5.
- Slowest file: `test-mr-inference.R`, 164.7 s.

## Step 4 — release check

- `R CMD build ~/Documents/GitHub/forestsearch` from `/tmp/cran8/build`, then
  `R CMD check --as-cran forestsearch_0.5.0.tar.gz`.
- `NOT_CRAN` and `OMP_THREAD_LIMIT` were unset. `R_LIBS=/tmp/cran8/lib` was set (the task's
  library rule).
- Nothing was written into the repository: `R CMD build` copies the whole checkout (14 GB,
  `quarto/` 5.9 GB) to a temporary directory before it drops the build-ignored files.
- Run 1's tarball contained `legacy_vignette/` (27 entries; untracked, excluded through
  `.git/info/exclude`) and `tests/testthat/Rplots.pdf` (git-ignored, left by test runs).
- Repair `f48c8f5f` added two lines: `^legacy_vignette$`, named by the NOTE, and
  `^tests/testthat/Rplots\.pdf$`, per Amendment 1.
- Run 2's tarball: 1,660,903 bytes (run 1: 8,531,649). `tar -tzf … | grep -c -E
  "legacy_vignette|Rplots\.pdf"` gives 0, so neither stray is in it.
- Installed size (`du -sk` of the installed package in `.Rcheck`): 4,760 KB in both runs.

## Commits

- Start SHA: `d57378f4cc712b3be42262333c2125df2096cf9d`, the task-file commit. Its parent is
  `d68aa78d`, HEAD at the start.

| SHA | Commit |
|---|---|
| `d57378f4` | dev: task CRAN preparation, finalize (as received) |
| `87688d39` | fix: run_simulation_analysis() passes log_IRR to the count GRF fit (finalize, Step 1) |
| `b1394deb` | test: skip the two CRAN-mode files over 30 s on CRAN (finalize, Step 2) |
| `844926db` | dev: task CRAN preparation finalize, amendment 1 (Rplots.pdf build-ignore line) |
| `f48c8f5f` | chore: build-ignore legacy_vignette/ and tests/testthat/Rplots.pdf (finalize, Step 4.5 + amendment 1) |
| `9e9a7a73` | docs: cran-comments.md test environment and check results (finalize, Step 5) |
| (this commit) | dev: CRAN preparation report, finalize |

## `cran-comments.md` as committed (`9e9a7a73`)

```
## Submission

This is an update from forestsearch 0.4.0 (the current CRAN version) to
0.5.0. NEWS.md lists the changes under "forestsearch 0.5.0".

This release removes `create_gbsg_dgm()` and `simulate_from_gbsg_dgm()`,
which were deprecated in 0.4.0 (the first CRAN release to export them,
where they already warned that they were deprecated).

## Test environments

- Local: Pop!_OS 24.04 LTS, x86_64-pc-linux-gnu, R 4.6.2 (2026-10-29);
  `R CMD check --as-cran` on the tarball built by `R CMD build`.
- win-builder (devel and release): TO BE FILLED BY MAINTAINER
- R-hub: TO BE FILLED BY MAINTAINER

## R CMD check results

0 errors | 0 warnings | 0 notes

The local check above ended with `Status: OK`.

## Reverse dependencies

There are no reverse dependencies: CRAN's package index (read on
2026-11-12) lists no package that depends on, imports, links to or suggests
forestsearch.
```

## Final checks

Run with this report and the check log staged. All four were run again after this commit; that
result is in the session's closing message.

1. This session's commits (`d57378f4` `87688d39` `b1394deb` `844926db` `f48c8f5f` `9e9a7a73` and
   this one) touch only these files, all from the allowed list (PASS):
   - `R/run_simulation_analysis.R`, `NEWS.md`, `cran-comments.md`, `.Rbuildignore`;
   - `tests/testthat/test-sim-grf-count-row.R`, `tests/testthat/test-mrct-case-study.R`,
     `tests/testthat/test-mrct-region-sims-mr.R`;
   - the task file, this report and the check log.
2. `DESCRIPTION` says `Version: 0.5.0`, and no file under `man/` changed (PASS).
3. `git status --porcelain` showed only the two staged files before this commit.
4. Nothing was pushed, merged or tagged (PASS):
   - `origin/cran-prep` is `d68aa78d`, pushed at 07:31:59, before this session's first commit
     at 07:35:08;
   - no session commit is on it;
   - there are no merge commits;
   - no tag contains a session commit (tags: `v0.1.0`, `v0.4.0`).

## Open

- **Step 1 unaffected case, as written:** the 6 count FS and FSlg rows are not `identical()`. They
  differ only in the four per-replicate GRF agreement columns (`agree.FS.GRF`, `kappa.FS.GRF`,
  `agree.FSlg.GRF`, `kappa.FSlg.GRF`). The fix is committed. If that reading is not accepted,
  `git revert 87688d39` undoes the fix, its test and the NEWS bullet.
- testthat 3.3.2 does not export `on_cran()`. CRAN mode was evidenced with `testthat:::on_cran()`.
- `test-mrct-region-sims-mr.R` was 30.5 s in CRAN mode, 0.5 s over the line, and is now skipped
  on CRAN.
- `mr_method_tables` takes 9.4 s in the `--run-donttest` pass (its `\donttest{}` block). The
  check gave no NOTE for it.
- win-builder and R-hub have not been run. `cran-comments.md` leaves their lines for the
  maintainer.
- After a count GRF fit, the `grf_full` keeper stores `effect_measure = "log_IRR"`.
- Scratch: `/tmp/cran8/lib` and `/tmp/cran8/src` were deleted. Kept: `/tmp/cran8/build` (run 2's
  tarball and `.Rcheck`), and `/tmp/cran8/start` and `/tmp/cran8/out` (scripts, logs, run 1's
  tarball). `/tmp/chk050/lib` and `/tmp/chk050/src` were already absent.
  `~/fs050_check/payloads/` was not touched.
