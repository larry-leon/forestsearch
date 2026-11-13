# REPORT — CRAN preparation, part 7: two CRAN-policy fixes and the count GRF row (2026-11-12)

Task: `dev/tasks/TASK_cran_prep_policy_fixes_2026-11-12.md`. Session `fs-cran-prep-policy-fixes`,
host `pop-os`, branch `cran-prep`, version 0.5.0.

| Fix | Result | Commit | Test file and counts | Unaffected case | Default path (6 cases) |
|---|---|---|---|---|---|
| 1 `test-mrct-case-study.R` capped at 2 workers | PASS | `4482a0d9` | `test-mrct-case-study.R`: 7 tests / 64 expectations, 0 failed, 0 skipped (no red step: tests-only fix) | n/a (no `R/` change) | n/a (no `R/` change) |
| 2 `mrct_region_sims()` restores `progressr.handlers` | PASS | `f6bb74a7` | `test-mrct-handlers-restored.R` (new): red on start 2 tests / 6 expectations, **2 failed**; green on fix 2 tests / 6 expectations, 0 failed | `identical()` TRUE (no timing columns in the returned table) | PASS, all six identical to Step 0 |
| 3 count GRF row | **STOPPED** | none | none written | not run | not run |

- Start SHA (after the task-file commit): `9c94727bc91d5eb951fb2dcac55cc4ff9eb66886`.
- This session's commits: `9c94727b` (task file), `4482a0d9` (fix 1), `f6bb74a7` (fix 2), plus the
  commit of this report.
- Step 0 gate: branch `cran-prep`, `Version: 0.5.0`, no modified tracked file under `R/`, `tests/`,
  `man/`, `NEWS.md`. PASS.

## Fix 1

- Edit: `workers = min(6L, max(2L, parallel::detectCores()))` → `workers = 2L` (the only change).
- The tests **ran** (`speff2trial` and `latentcor` are installed); none skipped.
- Workers load the installed build. The installed copy in the user library is 0.4.1.9002, so
  the start export was installed into a scratch library (`/tmp/cran7/lib`, `R_LIBS`) and the
  MRCT test files ran with it. Nothing was installed into the user or site library.
- Wall clock after: 149.9 s (first run), 90.8 s (second run, in the fix 2 gate). The machine was
  also running the other session's simulations.
- Wall clock before: **not re-measured**. The original line starts 6 workers on this host, which
  the task's 2-worker limit rules out. The last recorded figure is 84.4 s
  (`REPORT_cran_prep_displays_se_docs_2026-11-12_suite_diff.txt`, with up to 6 workers).

## Fix 2

- `R/mrct_simulation.R`: save `getOption("progressr.handlers")` before
  `progressr::handlers("progress")`, restore with `on.exit(..., add = TRUE)`. `options(x = NULL)`
  leaves an unset option unset.
- Other persistent state in `mrct_region_sims()`: the future plan, already restored by the existing
  `on.exit(future::plan(old_plan))`. No other option, environment variable or RNG kind is set in the
  function body.
- Test wall clock 15.3 s (> 5 s), so both tests have `skip_on_cran()`. The green run used
  `NOT_CRAN=true`: 14.6 s.
- Existing files that call `mrct_region_sims()`, all green on the fixed tree (`NOT_CRAN=true`,
  scratch library):
  `test-mrct-core-summary.R` 4/32, `test-mrct-region-sims-mr.R` 4/71,
  `test-mrct-region-sims-record.R` 6/172, `test-plan-restored.R` 8/20,
  `test-mrct-case-study.R` 7/64, `test-mrct-default-confounders.R` 1/8 (tests/expectations, 0
  failed, 0 skipped).
- NEWS: one bullet, first under `# forestsearch 0.5.0` → `### Bug fixes`.

## Fix 3 — table, and why it stopped

`.run_grf_analysis_gen()` (`R/run_simulation_analysis.R:1340`) passes `fs_params$effect_measure`
unchanged. Survival (or `outcome_type` NULL) goes to `grf.subg.harm.survival()`, which gets no
`effect_measure`. In `grf.subg.harm.glm()` the value is checked only for `outcome_type = "count"`
(`R/grf_subg_harm_glm.R:346`). Apart from that check and the no-offset warning, the value is
only stored in the result. The causal-forest fit does not use it. Acceptance was confirmed with
a direct call on the start export (`.make_count_data()` / `.make_binary_data()`, offset `ftime`).

| Outcome type | Effect measure (simulation) | Value passed | Accepted by `grf.subg.harm.glm()` |
|---|---|---|---|
| survival / NULL | — | none (survival GRF) | n/a |
| binary | OR | `"OR"` | yes (not validated; stored only) |
| binary | RR | `"RR"` | yes (not validated; stored only) |
| binary | RD | `"RD"` | yes (not validated; stored only) |
| binary | NULL | none → callee default `"log_OR"` | yes |
| continuous | MD | `"MD"` | yes (not validated; stored only) |
| continuous | NULL | none → callee default `"MD"` | yes |
| count | IRR | `"IRR"` | **no**: `For outcome_type = 'count', effect_measure must be 'log_IRR' or 'log_RR'. Got: 'IRR'.` |
| count | IRD | `"IRD"` | **no**: `... Got: 'IRD'.` |
| count | NULL | none → callee default `"log_IRR"` | yes |

The count values the callee accepts:
- `"log_IRR"`: log incidence rate ratio. It warns when `offset.name` is NULL and otherwise
  continues; in `create_glm_row()` it means a Poisson (or quasi/negbin) GLM with the log-exposure
  offset.
- `"log_RR"`: log relative risk. In `create_glm_row()` it means a Poisson GLM with a sandwich
  variance (modified Poisson). It gives no no-offset warning.

The call is wrapped in `tryCatch(..., warning(label, " failed: ", ...))`. Each count GRF row gets
the extractor's no-result values plus that warning.

**Why it stopped.** `IRR` → `"log_IRR"` follows from the source. The extractor's IRR is the
exponentiated Poisson coefficient with the log-exposure offset when the offset column is in the
data, and that is the quantity `"log_IRR"` names. The callee also rejects `IRD`, which the
simulation supports (`run_simulation_analysis()` docs; extractor family switch at
`:1441`/`:1450`). The extractor reports a rate difference for `IRD`, and the callee has no value
that estimates a rate difference. So the source does not determine a mapping for `IRD`. Under the
task's stop rule, no change was made: no edit, no test, no NEWS line. The working tree had no
Fix 3 edits to revert.

## Left open

- Fix 3 is not applied. Count GRF rows still fail for both `IRR` and `IRD`. The `IRR` → `"log_IRR"`
  mapping is determined by the source. `IRD` needs a decision; no callee value matches it.
- Fix 1 "before" wall clock was not measured (the 2-worker limit). The recorded earlier figure is
  given above.
- `test-mrct-case-study.R:7` still says "install before running"; that comment was not touched.
- The installed `forestsearch` in the user library is 0.4.1.9002. MRCT test files that run
  workers use whatever is installed unless a scratch library is put first on `R_LIBS`.
