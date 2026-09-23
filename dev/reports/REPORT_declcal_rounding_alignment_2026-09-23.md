# REPORT -- kappa-hat / FW-hat aligned with the screen's rounded admission rule (2026-09-23)

Task: `dev/tasks/TASK_declcal_rounding_alignment_2026-09-23.md` (committed `97b36dd6`).
Outcome: **all five gates pass** (Gate 5 under the amended criterion below; Gate 4
passes on an empty admitted set, see its caveat).

## Environment

| item | value |
|---|---|
| HEAD at start | `97b36dd6` (task doc commit; parent `0b9d1bdb`) |
| installed forestsearch | 0.3.5.9000, built 2026-09-23 02:32 UTC |
| R | 4.6.1 (2026-06-24), x86_64-pc-linux-gnu |
| machine | pop-os, Linux 7.1.5 |
| pre-existing untracked (never staged) | `quarto/simulations/actg175/binary_020/mr_or_harm/fs_effMaxSG_mr_field_or075_n500_nb20_{redes,relaunch}_d5000/`, `quarto/simulations/actg175/binary_020/smoke_{redes,relaunch}.html`, `quarto/simulations/gbsg_020/scripts_dinamr/logs/nullmr_findings.err` |

Every GBSG fit used `devtools::load_all()` with `parallel_args = list(plan = "sequential")`,
so it exercised the source tree, not the installed build. Each fit took about 42 s, well
inside the task's 30-minute compute abort. Gate 5's full suite and `R CMD check` are
verification the task requires, not task compute.

## Amendments to the committed task doc (from Larry, 2026-09-23)

None of these are in the committed task doc.

1. **Dependency satisfied.** `pconsistency.digits` is a `forestsearch()` argument
   (committed `06ac5391`, an ancestor of HEAD), so `args_call_all` carries it. On the
   GBSG fit, `digits_source = "fit$args_call_all$pconsistency.digits"` and the value is 2.
2. **Gate 5 criterion.** The `R CMD check` baseline is now 1 NOTE, "Version contains
   large components (0.3.5.9000)", which is deliberately left open as the dev-version
   marker. The gate passes if the modified tree gives exactly that NOTE and nothing
   else. The version was not bumped.
3. **Post-condition 8 amended (mid-task)** to add `tests/testthat/` to the permitted
   files, including the new `test-declcal-rounding-alignment.R`. Gate 5 requires
   updating tests that encode the old exact-threshold FW-hat, so the original list
   contradicted it.

## Step 1 -- baseline

Unmodified tree (HEAD `97b36dd6`). The GBSG application ran at its own settings as
`quarto/simulations/gbsg_app_null/run_gbsg_app_null.R` (`e6477a22`) sets it up: the real
`survival::gbsg` frame, `p* = 0.90`, `consistency_method = "resample"`, and `effMaxSG`.
It used the application's MR arguments from `fs-glms-interpretable/quarto/gbsg/analysis_gbsg_mr.qmd`
(5000 Poisson draws, `keep_declaration_field = TRUE`, and `declaration_c0 = c(0.70, 0.75, 0.80, 0.85)`).
Then `fs_declaration_calibration(fit, alpha = 0.10, c0 = 0.75)` was run on it. The
payload is saved at `~/Downloads/declcal_rounding_baseline_2026-09-23.rds` and is not
committed.

| quantity | baseline |
|---|---|
| kappa-hat(c0 = 0.75, alpha = 0.10) | 2.7262897414 (Mstar_c0 from the capture) |
| FW-hat(c0 = 0.75) | 0.6508 |
| pstar_implied | 0.9935949254 |
| admitted by kappa-hat (pre-reduction family, 1744) | none (0) |
| admitted by the executed screen | 9; `sg.harm` = `{er <= 0} & {pgr <= 26}` |
| unshifted kappa-hat(0.10) / FW-hat | 3.5029697507 / 0.9700 |

## Section 2 -- the helper

`.fs_pcons_eff(p_star, digits)` in `R/fs_declaration_calibration.R` is the only place the
effective threshold `g - 0.5 * 10^-digits` is written, with `g` equal to p* rounded up to
the grid (post-condition 7). `grep` over `R/` finds the expression only in that helper
and its roxygen. The settable-p* search does not repeat it: it steps along the grid and
evaluates every candidate through the helper.

**One deviation from the literal expression.** `ceiling(p_star * 10^digits)` is applied
to `round(p_star * 10^digits, 6)`, not to the raw product. The raw product carries
representation error that `ceiling()` pushes up a whole grid step. For example,
`0.07 * 100` is `7.000000000000001`, so the literal form gives `g = 0.08`. None of the
task's Gate 3 values triggers this: 0.90, 0.95 and 0.99 scale to exact integers, and
0.9936 scales to 99.36, 993.6 or 9936. The acceptance tests pin `.fs_pcons_eff(0.07, 2) == 0.065`.

**`consistency_method = "split"`, verified from source.** The Cox MR path runs under a
split-consistency search (`R/forestsearch_main.R:3696-3703`: `.mr_cox_ok` does not
require resample), so a split fit can carry a declaration field. On that path, `Pcons`
is `round(mean(flag), digits)` over literal splits (`R/subgroup_consistency_helpers.R:365`),
and admission is not a threshold on `T`. The calibration therefore handles split
explicitly:

- `fw_size`, `z_pstar`, the c0 table's `fw_size` and the settable columns are `NA`.
- `admitted_pstar` is `NULL`.
- `kappa_hat` and `admitted_calibrated` are still returned, since neither depends on
  the screen's rule.
- The print method says so, and the roxygen documents it.

The method is read from `args_call_all$consistency_method`. A bare `fs_mr_inference()`
result has no `args_call_all`. It is treated as resample, because its field is the
closed-form statistic, and `consistency_method_source` records that assumption.

## Section 3 -- changes to `R/fs_declaration_calibration.R`

- **FW-hat** is taken at `z_pstar = qnorm((1 + pcons_eff) / 2)`, with `pcons_eff = .fs_pcons_eff(p_star, digits)`,
  in both the unshifted and c0 paths. The exact-threshold value is neither computed nor
  stored (post-condition 6). `z_pstar` keeps its name but now holds the effective z
  cutoff, and `pcons_eff` is returned beside it.
- **`admitted_pstar`**, the relabelled current rule, uses the same effective threshold.
  This is required for it to *be* the current rule (see the Gate 4 supplement).
- **kappa-hat** is untouched.
- **`pstar_implied` is removed** from the c0 table. In its place are:
  - `pstar_settable`, `pstar_achievable`, `pcons_eff_settable`, `z_eff_settable` and
    `z_gap` (positive = conservative), all at the fit's digits;
  - `digits_fine`, `pstar_fine` and `z_gap_fine`: the smallest digits in 1..12 where the
    z gap is below 0.01.

  When no p* <= 1 reaches kappa-hat at the fit's digits, `pstar_achievable = FALSE`, the
  columns are `NA`, and the print method says so plainly.
- **`digits`** comes from `args_call_all$pconsistency.digits`. When that is absent, the
  fallback is `formals(subgroup.consistency)$pconsistency.digits` (2), following the
  precedent in `declcalc0_run.R:316-320`. `digits_source` records which route was used.
- New top-level elements: `digits`, `digits_source`, `pcons_eff`, `consistency_method`
  and `consistency_method_source`.
- `fw_size` keeps the strict `>` comparison: `mean(Mstar > z_pstar)`. The screen admits
  at `>=`, but `Mstar` is continuous, so the difference has measure zero. It was left
  unchanged to keep the diff to what the task asks for.

## Gates

### Gate 1 -- kappa-hat untouched: PASS

| | baseline | modified | `identical()` |
|---|---|---|---|
| kappa-hat(c0 = 0.75, alpha = 0.10) | 2.7262897414 | 2.7262897414 | TRUE |
| kappa-hat(0.10), unshifted | 3.5029697507 | 3.5029697507 | TRUE |

`admitted_calibrated` (c0 = 0.75) and `admitted_current` are also `identical()` to the baseline.

### Gate 2 -- FW-hat rises: PASS

| | baseline (exact z = 1.644854) | modified (effective z = 1.621082) |
|---|---|---|
| FW-hat(c0 = 0.75, alpha = 0.10) | 0.6508 | **0.6672** |
| FW-hat unshifted | 0.9700 | 0.9744 |

### Gate 3 -- helper vs `round()`: PASS, two boundary disagreements enumerated

The check covered digits 2, 3 and 4 against p* in {0.90, 0.95, 0.99, 0.9936}. Each grid
held 200,001 evenly spaced points over `pcons_eff +/- 2 * 10^-digits`, plus
`pcons_eff` itself and its two neighbouring doubles, clipped to [0, 1]. That is
2,300,038 points in all. `round(Pcons, d) >= p*` and `Pcons >= pcons_eff` disagree at
exactly two points, both of them the exact boundary double:

| digits | p* | Pcons (stored double) | `round(Pcons, 2)` | round rule | helper rule |
|---|---|---|---|---|---|
| 2 | 0.95 | 0.94499999999999995 | 0.94 | FALSE | TRUE |
| 2 | 0.99 | 0.98499999999999999 | 0.98 | FALSE | TRUE |

The stored doubles for 0.945 and 0.985 sit just below the half. R's `round()` goes down
there, while the helper's threshold is the same double, so `>=` admits. The other ten
(digits, p*) cells have no disagreement, including 0.895, where the stored double is
just above the half. These points are reported, not special-cased.
`test-declcal-rounding-alignment.R` asserts that every disagreement equals the
boundary value and that each cell has at most one.

### Gate 4 -- the settable pair round-trips: PASS (on an empty set; see caveat)

The pair reported for the GBSG fit, c0 = 0.75, alpha = 0.10:

| | value |
|---|---|
| `pstar_settable` at digits = 2 | **1.00** (achievable) |
| `pcons_eff_settable` / `z_eff_settable` | 0.995 / 2.807034 |
| `z_gap` | +0.080744 (conservative) |
| `digits_fine` / `pstar_fine` / `z_gap_fine` | **4** / **0.9937** / +0.002849 |

`forestsearch()` was re-run at the application's settings with
`pconsistency.threshold = 1.00, pconsistency.digits = 2`, and again at
`pconsistency.threshold = 0.9937, pconsistency.digits = 4`. Both admit **0** candidates
and declare nothing (`n_passed = 0`, `sg.harm` empty; MR was skipped because there was
no subgroup). kappa-hat also admits 0: none of the 1744 pre-reduction candidates and
none of the 120 screened ones, where the maximum T is 2.564499, below 2.726290. The
sets match.

**Caveat.** The match is on the empty set, so it confirms the guidance is usable here
but does not stress the case where candidates sit between kappa-hat and `z_eff_settable`.
At digits = 2 the settable screen is 0.081 z above kappa-hat. It can therefore admit a
strict subset of what kappa-hat admits, and on another fit with candidates in
[2.726, 2.807) it would. The roxygen says so: the pair is guidance on what to set, not
an identity.

**Supplement: the rounded relabelling reproduces the executed screen; the exact one did
not.** At the fit's own p* = 0.90, digits = 2, the executed screen admitted 9 of the 120
screened candidates. Relabelling at the effective z (1.621082) reproduces those 9
exactly. Relabelling at the exact z (1.644854) misses one: `q2.1 & q4.1`, with
T = 1.638366 and unrounded Pcons = 0.898655, which rounds to 0.90. That candidate is
`{pgr <= 0} & {er <= 9}` (N 71, E 42, HR 1.683928). It is the only one of the 9 admitted
rows without `q1.1` = `{er <= 0}` and the only one at a rounded Pcons of 0.90. It is the
same candidate that `REPORT_pconsistency_digits_argument_2026-09-23.md` (lines 101-103)
bounded to [0.895, 0.90) by re-running at 3 digits and watching it drop out. **Two
independent routes found the same candidate:**

1. the admission-table route, from the digits task;
2. the closed-form T of the declaration field, here.

The value 0.898655 lies inside that bound. The baseline's `admitted_pstar` omitted it;
the modified one includes it, so `admitted_pstar` gains exactly `q2.1 & q4.1` (22 -> 23
on the pre-reduction family).

### Gate 5 -- tests and `R CMD check`: PASS

- **Acceptance files** (`load_all()` + `test_file()`):

  | file | tests | expectations | failures |
  |---|---|---|---|
  | `test-declcal-rounding-alignment.R` | 5 | 115 | 0 |
  | `test-declaration-calibration.R` | 14 | 75 | 0 |
  | `test-declaration-c0.R` | 12 | 71 | 0 |

- **Full suite** (`devtools::test()`): 53 files, 473 tests, 5,861 expectations. 0 failed,
  0 errors, 3 skipped, 32 warnings.
- **`rcmdcheck::rcmdcheck(args = "--as-cran")`**, with the PDF manual built and
  `RSTUDIO_PANDOC` set: 0 errors, 0 warnings, 1 NOTE. The NOTE is
  `checking CRAN incoming feasibility ... NOTE` / "Version contains large components
  (0.3.5.9000)", which is exactly the amended baseline and nothing else.

## Test changes, individually

1. **New:** `tests/testthat/test-declcal-rounding-alignment.R`, 5 tests:
   - helper values, including the representation-error case `0.07`;
   - Gate 3 on a 20,001-point grid per cell;
   - minimality of the settable p*, the unachievable case (kappa 5 at digits 2), and
     `digits_fine = 4` for the GBSG kappa;
   - `fw_size` at the rounded threshold, kappa-hat unchanged, and `pstar_implied` gone;
   - digits and consistency_method resolution, including the split path returning `NA` / `NULL`.
2. **`test-declaration-calibration.R`, test 3** (independent columns). `alpha1` was
   `1 - pnorm(qnorm((1 + ps) / 2))`, the exact threshold; it is now taken at
   `pcons_eff` = 0.795 / 0.895. The old form passed only inside its 4-MC-SD tolerance.
3. **`test-declaration-calibration.R`, `.fw_target` / `.check8`** (test 8, both the
   exported path and the helper path). They failed against the new definition, as
   expected.
   - Added `.z_eff = qnorm((1 + 0.895) / 2)` and `.fw_target_eff(rho)`.
   - 8b compares `fw_size` with `.fw_target_eff`.
   - 8d's FW bounds `[0.05, 0.0975]` become `[1 - pnorm(.z_eff), 1 - pnorm(.z_eff)^2]` = [0.0525, 0.1022].
   - `.fw_target` itself is kept: test 8d's closed-form reference table ("targets
     bracket and decrease in rho") is a statement about the bivariate normal at
     `qnorm(0.95)`, not about the package's FW-hat.
   - The kappa checks (8c) are unchanged.
4. **`test-declaration-c0.R`, test 5.** Three occurrences, all updated:
   - `.fw5_target` used `qnorm(0.95)` and now uses `.z_eff5 = qnorm((1 + 0.895) / 2)`.
     The target moves from 0.045178 to 0.047572.
   - The hard-coded unshifted prerequisite `0.078751` is replaced by the same closed
     form computed at `.z_eff5`, `1 - .phi2ab(.z_eff5, .z_eff5, .rho5)` = **0.082443**.
     With rho5 = 0.734879, the expression gives 0.078751 at `qnorm(0.95)`, so the
     replacement is the committed constant's own formula moved to the rounded threshold.
   - The helper-path line `mean(m > qnorm(0.95))` becomes `mean(m > .z_eff5)`.

   Both shifted checks already passed within the 0.0027 tolerance, but they asserted the
   old definition.

**On the diagnostic line "old unshifted c0 constant check: 0.078751 new: NA".** That
line was not a computation, in the test's scope or outside it. It was a throwaway
`cat()` in which I typed the literal `NA` as a placeholder and never filled it in. No
expression evaluated to `NA`, and nothing reached an expectation. The replacement value,
computed afterwards, is 0.082443, as above. The suite is green because the updated test
computes that value inline.

Test 1 of `test-declaration-calibration.R` was **not** changed. It still asserts
`rate >= p_star` (unrounded) against the relabelled set. It passes on its fixture
(`sg_focus = "hr"`, `hr.threshold = 1.1`) because no screened candidate there falls in
[0.895, 0.90), and it asserts no FW-hat or implied-p* value. It is flagged below.

## Post-conditions

| # | check | status |
|---|---|---|
| 1 | kappa-hat identical to baseline, >= 4 decimals | PASS (`identical()`; 2.7262897414) |
| 2 | new FW-hat > baseline, both reported | PASS (0.6508 -> 0.6672) |
| 3 | helper vs `round()`, disagreements enumerated | PASS (2 of 2,300,038, both boundary doubles) |
| 4 | settable pair reproduces kappa-hat's admitted set | PASS (both empty; caveat above) |
| 5 | tests pass; check clean; updated tests listed | PASS |
| 6 | no exact-threshold FW-hat retained | PASS: not computed, stored or printed |
| 7 | effective threshold in one place | PASS: `.fs_pcons_eff()` only |
| 8 (amended) | files confined to the permitted list | PASS: `R/fs_declaration_calibration.R`, `man/fs_declaration_calibration.Rd`, `NEWS.md`, `tests/testthat/` (3 files), `dev/tasks/`, this record |

## Out of scope, recorded

- **MR admission** (`R/fs_mr_inference.R:660-661`) still uses
  `qnorm((1 + p_star) / 2)` and carries the same misalignment. Not changed, by decision.
- **`declcalc0_run.R` / `declcal_c0approx_run.R`** compute their own
  `pstar_implied = 2 * pnorm(kappa) - 1` and an inline `round(rate, digits) >= p_star`.
  They are not changed, but they no longer match the package's c0 table columns. They
  could call `forestsearch:::.fs_pcons_eff()` when next touched.
- **`dev/verification/report_values*.R`** read `pstar_implied` and `z_pstar`. They are
  historical scripts and were not changed; `pstar_implied` no longer exists, and
  `z_pstar` now means the effective cutoff.

## Side issues (flagged, not fixed)

1. `test-declaration-calibration.R` test 1: the second `expect_setequal` and its comment,
   "admission <=> T >= z_{(1+p*)/2}", encode the exact threshold. The test passes only
   because its fixture has no candidate in the rounding band. It should compare
   `round(rate, 2) >= p_star`.
2. The GBSG application's analysis document (`fs-glms-interpretable`) reports FW-hat and
   `pstar_implied` from the old definition. Re-rendering belongs to the applications chat.
