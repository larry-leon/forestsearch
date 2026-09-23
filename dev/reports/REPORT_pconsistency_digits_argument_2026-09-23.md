# REPORT -- `pconsistency.digits` exposed as a `forestsearch()` argument (2026-09-23)

Task: `dev/tasks/TASK_pconsistency_digits_argument_2026-09-23.md` (committed `de2dd97f`).
Outcome: **all gates pass** (Gates 4 and 5 under the amended criteria below).

## Environment

| item | value |
|---|---|
| HEAD at start | `de2dd97f` (task doc commit; parent `4856db0b`) |
| installed forestsearch | 0.3.5.9000, built 2026-09-23 02:32 UTC |
| R | 4.6.1 (2026-06-24), x86_64-pc-linux-gnu |
| machine | pop-os, Linux 7.1.5 |
| pre-existing untracked (never staged) | `quarto/simulations/actg175/binary_020/mr_or_harm/fs_effMaxSG_mr_field_or075_n500_nb20_{redes,relaunch}_d5000/`, `quarto/simulations/actg175/binary_020/smoke_{redes,relaunch}.html`, `quarto/simulations/gbsg_020/scripts_dinamr/logs/nullmr_findings.err` |

All GBSG fits used `devtools::load_all()` with `parallel_args = list(plan = "sequential")`,
so they exercise the source tree, not the installed build.

## Amendments to the committed task doc (from Larry, 2026-09-23, mid-task)

None of these are in the committed task doc.

1. **30-minute hard abort lifted.** It was set for simulation compute and did not
   account for two `R CMD check` runs plus the full suite. Gates 4 and 5 ran to
   completion.
2. **Option 1 for the Gate 4 coverage-guard failure:** classify
   `pconsistency.digits` in `R/fs_family_report.R` alongside
   `pconsistency.threshold`. This changes `fs_family_report()`'s output, a
   behaviour change confined to a reporting function.
3. **Do not fix the pre-existing `effect_measure` failure**
   (`test-fs-family-report.R:115`). Record it as a finding.
4. **Gate 4 pass criterion amended:** the suite passes if its only remaining
   failure is the pre-existing `effect_measure` one, confirmed as failing on the
   pre-change baseline.
5. **Post-condition 7 amended:** `R/fs_family_report.R` is a permitted file.
6. **Gate 5 criterion amended:** the baseline is 2 ERRORs, 1 WARNING, 3 NOTEs,
   not clean. The gate passes if the modified tree produces no new items against
   that baseline.
7. **Roxygen de-duplication (documentation only), done after Gates 4 and 5 and
   before the commits.** `evaluate_subgroup_consistency()` and
   `evaluate_consistency_twostage()` drop their local `@param pconsistency.digits`
   and take `@inheritParams subgroup.consistency`, whose block is the canonical
   text. `forestsearch()` keeps its own block because it carries the "Passed to
   `subgroup.consistency()`" note. Then re-run `devtools::document()` and
   `R CMD check` against the baseline.

## Change

- `R/forestsearch_main.R`: new formal `pconsistency.digits = 2`, placed after
  `stop_threshold` (the threshold group), and passed explicitly in
  `consistency_overrides` to `subgroup.consistency()`. `args_call_all` is
  `mget(names(formals()))`, so the new formal is recorded without further code.
- `.make_eval_*` factories: **not touched.** They already forward
  `pconsistency.digits` (`R/subgroup_consistency_helpers.R`, `.make_eval_subgroup_consistency`
  and `.make_eval_consistency_twostage`).
- No change to any rounding, comparison (`p.consistency < pconsistency.threshold`)
  or selection-sort code.
- Roxygen rewritten at the four sites named in the task: `subgroup.consistency()`
  (canonical), the `@noRd` factory block, and (before amendment 7) the two
  exported evaluators. The new `forestsearch()` entry uses the same text plus
  the pass-through sentence.
- `R/fs_family_report.R` (amendment 2): `pconsistency.digits` added to the
  consistency-screen row's `arguments` and `values`
  (`... pconsistency.threshold = 0.9; pconsistency.digits = 2; ...`).
- `NEWS.md` entry under the development version.

## Step 1 baseline

The GBSG application at the settings `run_gbsg_app_null.R` (`e6477a22`) uses:
`survival::gbsg`, `treat.name = "hormon"`, `seedit = 8316951`, `p* = 0.90`,
`consistency_method = "resample"`, `fs.splits = 1000`, `stop_threshold = NULL`,
MR off. This was run on the unmodified source, and `grp.consistency$out_sg` was
saved to `~/Downloads/pcons_digits_baseline.rds` (not committed). Declared
subgroup: `{er <= 0} & {pgr <= 26}`, N 75, HR 2.221839, Pcons 0.99. 9 candidate
rows.

## Gate 1: behaviour unchanged at the default. PASS

| run | `identical(out_sg, baseline)` | declared | N | HR |
|---|---|---|---|---|
| `pconsistency.digits` omitted | TRUE (all 6 components: `result`, `pareto_frontier`, `sg.harm`, `sg.harm_label`, `df_flag`, `sg.harm.id`) | `{er <= 0} & {pgr <= 26}` | 75 | 2.221839 |
| `pconsistency.digits = 2` | TRUE (all 6 components) | `{er <= 0} & {pgr <= 26}` | 75 | 2.221839 |

## Gate 2: the argument arrives. PASS

`args_call_all$pconsistency.digits` is 2 on the omitted run, 2 on the explicit-2
run and 6 on the digits = 6 run. The baseline fit had no such entry.

## Gate 3: the pass-through is real. Recorded (no expected value)

`pconsistency.digits = 6` ran without error.

| candidate | Pcons (6 digits) | on 0.01 grid? | Pcons (2 digits) |
|---|---|---|---|
| declared `{er <= 0} & {pgr <= 26}` (N 75, HR 2.221839) | 0.985229 | no | 0.99 |
| max-effect `{er <= 0} & {size <= 35}` (N 61, HR 2.536918) | 0.989668 | no | 0.99 |

- The declared subgroup is unchanged: `{er <= 0} & {pgr <= 26}`, N 75, HR 2.221839.
- At 2 digits the declared and max-effect candidates tie at 0.99. At 6 digits
  the tie is gone and the same subgroup is still selected.
- The admitted table drops from 9 rows to 8. `{pgr <= 0} & {er <= 9}` (Pcons 0.90
  at 2 digits, HR 1.684) is no longer admitted, so its unrounded proportion lies
  in [0.895, 0.90). This is a live example of the rounding deciding admission at
  `p* = 0.90`.

## Gate 4: test suite. PASS (amended criterion)

| tree | `devtools::test()` line | failures |
|---|---|---|
| pre-change (clean export of `de2dd97f`) | `[ FAIL 1 \| WARN 32 \| SKIP 3 \| PASS 5710 ]` | `test-fs-family-report.R:115` (`effect_measure`) |
| modified, before amendment 2 | `[ FAIL 2 \| WARN 32 \| SKIP 3 \| PASS 5709 ]` | `:114` (`pconsistency.digits` unclassified, **new**) + `:115` |
| modified, after amendment 2 | `[ FAIL 1 \| WARN 32 \| SKIP 3 \| PASS 5710 ]` | `:115` only |

The final line is identical to the pre-change line. `test-directive-c.R`, which
treats `fs_family_report` as a p* echo site, passes after the change (182/182).

**`test-search-reproducibility.R`:** passes at the default (6/6). Run
additionally from a scratch copy with `pconsistency.digits = 6L` added to its
`forestsearch()` call, plus a probe test: passes (7/7). The probe confirmed the
6-digit Pcons values are off the 0.01 grid (0.966667, 1.000000, 0.733333,
0.566667, 0.733333) and that `args_call_all` carries 6L. The worker-count
invariance therefore does **not** depend on rounding masking differences.

## Gate 5: `R CMD check --as-cran`. PASS (amended criterion)

`rcmdcheck::rcmdcheck(args = "--as-cran")`, with `RSTUDIO_PANDOC` set, on
tracked-files exports: `git archive` of `de2dd97f` for the baseline and
`git ls-files` of the working tree for the modified runs.

| tree | status | tests line inside check |
|---|---|---|
| baseline `de2dd97f` | 2 ERRORs, 1 WARNING, 3 NOTEs | `[ FAIL 2 \| WARN 22 \| SKIP 77 \| PASS 5299 ]` |
| modified (Gates 1-4 state) | 2 ERRORs, 1 WARNING, 3 NOTEs | `[ FAIL 2 \| WARN 22 \| SKIP 77 \| PASS 5299 ]` |
| modified + amendment 7 | 2 ERRORs, 1 WARNING, 3 NOTEs | `[ FAIL 2 \| WARN 22 \| SKIP 77 \| PASS 5299 ]` |

With the `[Ns/Ms]` timing figures masked, every error, warning and note in both
modified runs appears word-for-word in the baseline. There are no new items.
The baseline items, none of them this task's:

- ERROR: `fs_dgm_feasibility` example (`--run-donttest`).
- ERROR: tests, from `test-fs-dgm-feasibility.R:103` (`readLines` of package
  source, unavailable in the check tree) and `test-fs-family-report.R:115`
  (`effect_measure`).
- WARNING: non-ASCII characters in `R/fs_bias_coverage.R`.
- NOTE: CRAN incoming (maintainer; large version component).
- NOTE: `fs_plot_bias_coverage` undefined globals (`b cell cov cov1 cov2 estimator obs r ref`).
- NOTE: HTML manual, `tidy` not installed.

## Amendment 7: Rd diff check

After `devtools::document()`, only `man/evaluate_subgroup_consistency.Rd` and
`man/evaluate_consistency_twostage.Rd` changed against the pre-amendment Rd set.
No other `.Rd` file changed. In both files the only differing line is the end of
the `pconsistency.digits` item:

```
< deliberate coarsening rather than a formatting step. Default 2.}
---
> deliberate coarsening rather than a formatting step.
> Default: 2}
```

The explanatory text is identical. The inherited canonical text (from
`subgroup.consistency()`) ends "Default: 2", where the removed local copies
ended " Default 2.". The canonical block was left unchanged, as instructed. No
other parameter was inherited or altered.

## Findings

1. **Pre-existing `effect_measure` failure** (`test-fs-family-report.R:115`):
   `effect_measure` is both classified by `fs_family_report()` and on the test's
   `.FR_OUT_OF_SCOPE` list. It fails on the pre-change tree and was not touched
   (amendment 3).
2. **The coverage guard did its job.** Adding the formal failed `:114` until the
   argument was classified. That was the guard's stated purpose ("forcing a
   decision").
3. **`devtools::document()` reports a pre-existing error** in
   `R/fs_bias_coverage.R:17` (inline code `` `r = se_mean / sd_emp` `` parsed as
   R). It is unrelated and was not touched.
4. **At `p* = 0.90` on the GBSG application, rounding admits one candidate that
   would otherwise fail** (`{pgr <= 0} & {er <= 9}`, Gate 3). The declared
   subgroup does not depend on it.

## Files changed

`R/forestsearch_main.R`, `R/subgroup_consistency_main.R` (roxygen),
`R/subgroup_consistency_helpers.R` (roxygen), `R/fs_family_report.R`,
`man/forestsearch.Rd`, `man/subgroup.consistency.Rd`,
`man/evaluate_subgroup_consistency.Rd`, `man/evaluate_consistency_twostage.Rd`,
`NEWS.md`, `dev/tasks/TASK_pconsistency_digits_argument_2026-09-23.md`, this
report. The Step 1 baseline `.rds` and the Gate 1-3 `.rds` files stay in
`~/Downloads`, uncommitted.

## Wall clock

Task start 12:32:37 (task-doc commit) to report commit (see the commit
timestamp). Main components: two full `devtools::test()` runs on the modified
tree (493 s, 487 s) and one on the baseline tree, three `R CMD check --as-cran`
runs (about 10 min each), and GBSG fits of about 4 s each.
