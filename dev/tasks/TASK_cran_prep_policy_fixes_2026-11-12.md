# TASK — CRAN preparation, part 7: two CRAN-policy fixes and the count GRF row (2026-11-12)

Session name: `fs-cran-prep-policy-fixes`

Host `pop-os`, repository `~/Documents/GitHub/forestsearch`, branch `cran-prep`. Three fixes,
all approved by Larry on 2026-11-12. Read this whole file before starting.

## Rules

- **Lean run.** The main session implements and gates each fix itself. No verifier, polish or
  report-drafting agents. One short report, no detail file.
- **No `R CMD check`, no `devtools::check()`, no vignette build, no full test suite.** Each
  fix is gated by its own test file, its named unaffected case and, for an `R/` change, the
  six default-path cases. The last green baseline is the part 6 report
  (`dev/reports/REPORT_cran_prep_release_2026-11-12.md`).
- **One fix, one commit** (the fix, its test and its NEWS bullet together). Path-limited
  commits; stage nothing else.
- **Another CC session is running simulations on this machine** (scratch `/tmp/chk050`, up to
  100 workers, payloads in `~/fs050_check`). Do not touch `/tmp/chk050`, `~/fs050_check`,
  fs-post-selection or fs-glms-interpretable. Leave that session's commits and files in
  forestsearch alone. If a commit fails on `index.lock`, wait 10 seconds and retry, up to 5
  times.
- **Scratch:** `/tmp/cran7`. No install into the user or site library. Use
  `pkgload::load_all()`. Tests and reproducers use at most 2 workers and 2 threads.
- **No push, pull, merge or tag.** The version stays 0.5.0.
- **Do not edit roxygen comments;** no `man/` file changes in this task.
- **Gates** stop on failure; they never stop to ask. A gate that fails for one fix stops that
  fix only: revert its uncommitted edits, record why, and go on to the next fix.

## Step 0 — start

1. Copy this file into `dev/tasks/` and commit it (path-limited). Record the start SHA.
2. **Gate:** branch is `cran-prep`; `DESCRIPTION` says `Version: 0.5.0`;
   `git status --porcelain` shows no modified tracked file under `R/`, `tests/`, `man/` or
   `NEWS.md`.
3. `git archive` the start SHA into `/tmp/cran7/start` (the "start export").
4. Run the six default-path cases on the start export, as part 6 ran them
   (`REPORT_cran_prep_code_fixes_2026-11-12_gate_run.R` and `_gate_compare.R`, unchanged).

## Fix 1 — tests only: a CRAN-run test starts up to 6 workers

- **Where:** `tests/testthat/test-mrct-case-study.R:42-44`, in `cs_runs()`:
  `min(6, max(2, detectCores()))`.
- **Change:** that test file never starts more than 2 workers. No `R/` change. No expected
  value in the file changes.
- **Gate:** the file's tests pass with the worker line as the only edit. Record whether the
  tests ran or were skipped (they need `speff2trial` and `latentcor`), and the file's wall
  clock before and after.
- **If an expectation fails with 2 workers:** do not edit expected values. Revert, and report
  the failing expectation.
- No NEWS bullet.

## Fix 2 — `R/`, change of behaviour: `mrct_region_sims()` leaves an option changed

- **Where:** `R/mrct_simulation.R:414`, `progressr::handlers("progress")`. The part 6 audit
  found that it sets the option `progressr.handlers` and does not restore it.
- **Change:** after `mrct_region_sims()` returns or stops, the caller's
  `progressr.handlers` setting is what it was before the call. The progress display during
  the call is unchanged. Suggested form: save `getOption("progressr.handlers")` before that
  line and restore it with `on.exit(options(progressr.handlers = <saved>), add = TRUE)`. If
  reading the source shows that the call changes anything else that persists, restore that
  too and say so.
- **Returned values do not change.**
- **Test** (new file `tests/testthat/test-mrct-handlers-restored.R`; fails before, passes
  after):
  - with the option unset before the call, it is unset after;
  - with a handler set by the caller before the call (for example `"void"`), the same setting
    is there after;
  - the test restores the session's own setting when it ends.
  - Use the smallest configuration the existing helper
    (`helper-mrct-region-sims-record.R`) provides, with 2 workers. Measure the test's wall
    clock; if it is over 5 seconds, add `skip_on_cran()`.
- **Unaffected case:** the value returned by `mrct_region_sims()` for that helper
  configuration is `identical()` between the start export and the fixed tree (same seed,
  timing fields dropped).
- **Gates:** red before and green after; the unaffected case; the existing MRCT test files
  that call `mrct_region_sims()` pass; the six default-path cases identical to Step 0.
- **NEWS:** one bullet under `# forestsearch 0.5.0`, "Bug fixes". No other NEWS line changes.

## Fix 3 — `R/`, change of behaviour: the GRF row fails in count simulations

- **Where:** `R/run_simulation_analysis.R:1340` (head; `:1294` at the part 6 start).
  `.run_grf_analysis_gen()` passes `effect_measure = "IRR"` to `grf.subg.harm.glm()`, which
  accepts only `"log_IRR"` and `"log_RR"`, so the GRF row of every count simulation fails.
- **Read first, and put this table in the report:** for every outcome type and effect
  measure, the value `.run_grf_analysis_gen()` passes, and whether `grf.subg.harm.glm()`
  accepts it. For the count values the callee accepts, say in one line each what they do.
- **Change:** `.run_grf_analysis_gen()` passes a value the callee accepts.
  - `"IRR"` maps to the callee value that estimates the same quantity the simulation
    extractor reports for IRR: the rate ratio, with the exposure offset when the offset column
    is in the data (the extractor's own rule).
  - If the table shows another value the callee rejects, fix it the same way and list it.
  - **Stop this fix** (revert, report the table, change nothing) if the source does not
    determine the mapping, or if the fix would need a change inside `grf.subg.harm.glm()`.
- **Nothing else changes:** the extractor, the truth, the FS rows, and the GRF rows for
  binary, continuous and survival outcomes.
- **Test** (new file `tests/testthat/test-sim-grf-count-row.R`; fails before, passes after):
  - on a small count DGM, the GRF analysis of one replicate completes without the
    argument-matching error (before the fix, record the error message);
  - the value passed is one of the callee's accepted choices, read from the callee itself and
    not typed into the test;
  - 2 threads; `skip_if_not_installed("grf")`.
- **Unaffected case:** `run_simulation_analysis()` with 3 replicates each for binary (OR),
  continuous (MD) and count (IRR), `seed_base` 8316951. Every row is `identical()` between
  the start export and the fixed tree (timing fields dropped) **except** the count GRF rows.
  Reuse part 6's unaffected-case script if it is committed.
- **Gates:** red before and green after; the unaffected case; the five existing extractor
  test files and `test-sim-effect-scales.R` pass; the six default-path cases identical to
  Step 0.
- **NEWS:** one bullet under `# forestsearch 0.5.0`, "Bug fixes". No other NEWS line changes.

## Step 4 — report and close

1. Write `dev/reports/REPORT_cran_prep_policy_fixes_2026-11-12.md`. Short:
   - **top table:** one row per fix — PASS / STOPPED, commit SHA, test file and counts (red
     before, green after), unaffected case result, default-path result;
   - the start SHA and this session's commits;
   - the Fix 3 table;
   - for Fix 1, whether the tests ran or were skipped;
   - measured wall clocks of the new and changed test files;
   - anything left open, in bullets. No recommendation beyond the facts.
2. Commit the report (path-limited) and copy it to `~/Downloads`.
3. **Final gates:**
   - this session's commits (list the SHAs) touch only: `tests/testthat/test-mrct-case-study.R`,
     `tests/testthat/test-mrct-handlers-restored.R`, `tests/testthat/test-sim-grf-count-row.R`,
     `R/mrct_simulation.R`, `R/run_simulation_analysis.R`, `NEWS.md`, this task file and the
     report;
   - `DESCRIPTION` still says `Version: 0.5.0`;
   - no file under `man/` changed;
   - nothing was pushed, merged or tagged.
4. Print the top table and the final-gate result as the closing message.
