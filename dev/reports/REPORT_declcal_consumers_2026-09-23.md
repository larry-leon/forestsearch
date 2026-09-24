# REPORT — Declaration-calibration consumers: inventory, stopped before fixes

- **Task:** `dev/tasks/TASK_declcal_consumers_2026-09-23.md` (`1b11967c`), with Larry's amendment (verbatim below).
- **Outcome: STOPPED at §1, by rule.** Two stop conditions hold, independently:
  1. **`R/` callout.** The widened search finds a second derivation of the screen's threshold **inside `R/`**:
     the OC predictor's resample gate uses the exact cutoff `qnorm((1 + p) / 2)`, not the rounded screen
     (`R/fs_oc_predict.R:279`, `R/fs_oc_grid.R:582`). The task says: "If a fix turns out to require an `R/`
     change, stop and report rather than making it." Removing every second derivation, as the amendment
     requires, needs this `R/` change.
  2. **Amendment's size rule.** The widened inventory is about 30 sites in about 20 files across six directories
     (below), against the original three starting points. Several fixes need a decision on meaning rather than
     an edit (§3 below). Gates 2 and 3 (full `devtools::test()`, `R CMD check --as-cran`) come on top. That is
     materially larger than the 1 h abort, so per the amendment: "stop and report the inventory rather than
     starting the fixes."
- **Nothing was fixed.** No script, test or `R/` file is modified. The stale test (§3 of the task) is also
  untouched, because it is one of the fixes. Gates 1–3 did not run; Gate 4 (no `R/` change) holds.

## 0. Record

| item | value |
|---|---|
| HEAD at start | `1b11967c` (task doc); `R/` last changed at `7713942e` |
| R, platform | R 4.6.1, pop-os, Linux 7.1.5 |
| Untracked before the task (never staged) | the two `actg175/.../_d5000/` directories, `actg175/binary_020/smoke_redes.html`, `smoke_relaunch.html`, `gbsg_020/scripts_dinamr/logs/nullmr_findings.err` |
| Clock | 18:46 start; stopped at the inventory |

## Amendment (Larry, 2026-09-23, received during the Section 5 sweep; recorded verbatim)

> 1. §1's search must also find locally-computed equivalents, not only references to the removed field: any
>    site computing 2*pnorm(kappa)-1, qnorm((1+p)/2), or its own version of the implied p-star or the screen's
>    threshold. Report each with file, line, and what it computes.
> 2. §2: a site that computes its own value is a worse defect than one that reads a removed field. A vanished
>    line is at least visible; a locally-computed pstar_implied_05 still prints a plausible number, now at the
>    exact cutoff and inconsistent with the aligned package, with nothing signalling it. Fix these by reading
>    the value from the calibration output, or by calling .fs_pcons_eff() where the screen's threshold is what
>    is wanted. Do not leave a second derivation anywhere.
> 3. Gate 1 is extended: for each fixed consumer, assert both that no formatted field is zero-length AND that
>    any implied-p-star or threshold value it reports matches the value the calibration returns. A number that
>    is merely present is not the test.
>
> If the widened scope makes the task materially larger than its 1 h abort allows, stop and report the
> inventory rather than starting the fixes.

## 1. Inventory (from source, this checkout)

Search: `pstar_implied`, `fw_size`, `z_pstar`, `pconsistency_digits`, `0.651`, `0.6508`, `0.9936`, and the
locally-computed forms `2 * pnorm(.) - 1`, `qnorm((1 + .) / 2)`, `round(rate, digits) >= p_star`, over `*.R`,
`*.qmd`, `*.Rmd`, `*.sh`. Excluded: `man/`, and the two sibling worktrees under `.claude/worktrees/`
(`mr-terminology`, `replication-check`), which are other checkouts, not this tree.

**Kinds:** **RF** reads a removed or changed field. **IP** computes its own implied p-star. **TH** computes its
own screen threshold (exact cutoff unless noted). **HC** hard-codes a value. **OK** is not a screen quantity
(a CI z, or the definition of Pcons itself) and is listed only because the pattern matched it.

### 1a. `R/` (the callout: not touched)

| file:line | kind | what it computes |
|---|---|---|
| `R/fs_oc_predict.R:279` | TH | `.fs_oc_gate()`, resample gate: `z_p <- qnorm((1 + pconsistency) / 2)`, the admission `Bhat - c2 >= z_p * se_g`. **Exact cutoff**: the OC prediction admits on a different threshold from the executed (rounded) screen. |
| `R/fs_oc_grid.R:582` | TH | `.fs_oc_reduce()`, the same exact-cutoff eligibility for the OC grid. |
| `R/fs_oc_grid.R:383` | TH | root-bracket padding `lo - zp * max(se_g)`; exact `zp`, harmless (a bracket), but a third copy. |
| `R/consistency_resample.R:453` | OK | `rate_closed = 2 * pnorm(delta / sigma_D) - 1`: the definition of Pcons, not a threshold. |
| `R/dina.R:1005`, `R/subgroup_search.R:838` | OK | CI z-values (`level`, 0.95). |
| `R/fs_declaration_calibration.R`, `R/fs_mr_inference.R:678` | — | the canonical sites (`.fs_pcons_eff()` and its callers). |

### 1b. Declaration-calibration scripts, `quarto/simulations/gbsg_020/scripts_dinamr/` (the task's own scope)

| file:line | kind | what it computes / does |
|---|---|---|
| `declcal_run.R:214-215`, `declcalc0_run.R:230-231`, `declcal_c0approx_run.R:237-238` | HC/TH | `z_exact <- qnorm((1 + 0.90) / 2)`; `z_round <- qnorm((1 + 0.895) / 2)`: the rounded threshold hard-coded as 0.895 (correct only at p* 0.90, digits 2) |
| `declcal_run.R:346-347`, `declcalc0_run.R:370-371`, `declcal_c0approx_run.R:369-370` | IP | `pstar_implied_05/10 <- 2 * pnorm(kappa_hat) - 1`, the **exact-cutoff** inverse, **stored in the payload schema** (`declcal_run.R:207`, `declcalc0_run.R:223`, `declcal_c0approx_run.R:224`) |
| `declcalc0_run.R:390`, `declcal_c0approx_run.R:389` | IP | per-c0 `pstar_implied_05_<c0> <- 2 * pnorm(k05) - 1`, stored |
| `declcal_run.R:363-364`, `declcalc0_run.R:413-414`, `declcal_c0approx_run.R:403-404` | TH | inline rounded screen `round(pmax(0, 2 * pnorm(Tp) - 1), digits) >= p_star`. This is the rounded rule, so it agrees in value, but it is a second copy (the task names `declcalc0_run.R:414`) |
| `declcal_run.R:298`, `declcalc0_run.R:322`, `declcal_c0approx_run.R:321` | RF (schema) | `r$pconsistency_digits`, a payload column the script fills from `fit$args_call_all`; not a removed field. Listed because the task searches for it. |
| `declcal_findings.R:56, 79-85` | IP (read) | reads the payload's `pstar_implied_05` and formats its median/IQR and quantile table ("how strict the calibration is"), so it **prints the exact-cutoff inverse** |
| `declcal_findings.R:111` | RF (schema) | prints `unique(r$pconsistency_digits)` (a payload column; fine) |
| `declcalc0_findings.R:29-30` | IP | `ptxt(k) = sprintf("%.5f", 2 * pnorm(k) - 1)`, formatted into the c0 tables, with a footnote stating the exact-cutoff definition |
| `declcalc0_findings.R:103, 111, 207-214` | IP (read) + HC/TH | reads `pstar_implied_05[_c0]`; `qnorm((1 + 0.895) / 2)` hard-coded as "the as-executed cutoff" |
| `declcal_c0approx_findings.R:111` | IP | `2 * pnorm(median kappa) - 1` formatted into the approximate table |

### 1c. `dev/verification/` (the task's own starting point)

| file:line | kind | what it does |
|---|---|---|
| `report_values_c0.R:12` | **RF: silent vanish** | `f(t05$pstar_implied[i])`: `pstar_implied` is gone from the c0 table (`test-declcal-rounding-alignment.R:90` asserts its absence), so `f(NULL)` (`formatC`, `:6`) is `character(0)` and `cat()` prints nothing for that field: it **disappears from each printed row** (checked: `cat("A", f(NULL), "B")` prints `A  B`). This is the failure mode the task was written for. |
| `report_values_c0.R:13, 22` | RF (changed) | `fw_size` read and printed; it still exists but now means the rounded rule. A `.fw5_target` comparison at `:22` is only as current as that constant. |
| `report_values.R:12-13, 29-38, 49` | RF (changed) | `fw_size` read and printed, and compared against targets (`fwt`) |
| `report_values.R:15, 50` | RF (changed) + TH | `z_pstar` read (now the rounded z); `:50` rebuilds the admission `beta_hat >= pmax(c_screen, c_cons + z_pstar * sigma_D)` inline, a second copy of the screen from the returned `z_pstar` |

### 1d. Other directories (in the amendment's scope: "anywhere")

| file:line | kind | what it computes |
|---|---|---|
| `quarto/simulations/gbsg_app_null/pstar_grid_findings.R:103` | HC | prints "FW_0.10(0.75) **0.651**": the pre-alignment FW-hat; `ac860c6d` moved it to 0.6672 |
| `quarto/resampling/fdr_family_multiplier.R:169` | TH | `z <- qnorm((1 + p_star) / 2)`; `t_g <- pmax(c_screen, c_consistency + z * sigma_D)`, exact |
| `quarto/resampling/consistency_resampling_theory.qmd:904, 1136` | TH | the same exact flag threshold, in the theory document's worked code (rendered) |
| `quarto/resampling/consistency_resample.R:164`, `validate_consistency_adjusted.R:92`, `gbsg_consistency_demo_standalone.R:33` | OK | the closed-form Pcons definition |
| `quarto/simulations/actg175/binary_methods/fdr_mr_inference.R:112, 199, 280` | TH | exact flag threshold in a script-level MR/FDR implementation (the pre-package binary methods) |
| `quarto/simulations/actg175/binary_methods/fdr_family_multiplier.R:169` | TH | the same, a copy of the `quarto/resampling` file |
| `quarto/simulations/actg175/continuous/scripts_mdf1/reselection_check.R:54` | TH | rebuilds MR's admission from `adm$consistency$p_star` at the **exact** cutoff; a check of MR that now disagrees with `7713942e` |
| `quarto/simulations/actg175/continuous/oc_wrapper_verification.qmd:45` | TH | `z_p <- qnorm((1 + pcons) / 2)`, the OC wrapper's gate re-derived, exact (mirrors `R/fs_oc_predict.R`) |
| `dev/glm-continuous-sims/verification/mr_mechanism_A1.R:203, 267`; `mr_mechanism_probe_superset.R:26`; `mr_mechanism_probe_restrict.R:26` | TH | MR admission rebuilt at the exact cutoff (verification of the pre-alignment MR) |
| `dev/glm-continuous-sims/sigma_d_diagnostic_2026-08-29.R:99` | OK | prints `qnorm(0.95)` as a labelled constant |
| `dev/identifier-alignment/code_theory_audit.qmd:457` | TH (prose) | quotes `z = qnorm((1 + p_star)/2)` as "§2.4 verbatim" |

### 1e. Tests

| file:line | kind | note |
|---|---|---|
| `tests/testthat/test-declaration-calibration.R:111-124` (test 1) | TH | the stale test named in task §3: `:122-124` rebuild the screen as `rate <- 2 * pnorm(T) - 1; rate >= dc$p_star`, the **unrounded** rule |
| `test-declaration-calibration.R:53, 155-158`; `test-declaration-c0.R:230` | HC | `.z_eff <- qnorm((1 + 0.895) / 2)`, independent oracles, updated by `96f84ad8`; a pinned constant in a test is deliberate |
| `test-mr-admission-rounded.R:52-53, 83-84`; `test-declcal-rounding-alignment.R:32-108` | — | the alignment's own tests; they compare exact against `.fs_pcons_eff()` on purpose |
| `test-fs-oc-predict.R:214` | TH | `thr <- max(30, 10 + qnorm((1 + pc) / 2) * 14)`: pins the OC predictor's **exact** gate, so it would change with any `R/fs_oc_predict.R` fix |

## 2. What the fixes would be (not started)

- **Mechanical, script-only.** `report_values_c0.R:12` (drop `pstar_implied`; print `pstar_settable` /
  `pstar_achievable` from the c0 table); `report_values.R:50` (take the admitted set from the calibration output
  instead of rebuilding it); the three inline `round(rate, digits) >= p_star` screens (call `.fs_pcons_eff()`);
  the `z_round <- qnorm((1 + 0.895) / 2)` constants (`qnorm((1 + .fs_pcons_eff(p_star, digits)) / 2)`);
  `pstar_grid_findings.R:103` (0.651 → read from the payload or update and pin); the stale test 1.
- **Needs an `R/` change.** `R/fs_oc_predict.R:279` and `R/fs_oc_grid.R:582` (and `:383`) should gate on
  `.fs_pcons_eff(pconsistency, digits)`. That also changes `test-fs-oc-predict.R:214`, and it needs a `digits`
  argument threaded into the OC functions.

## 3. Decisions for Larry before the fixes

1. **The committed `pstar_implied_05` payload columns.** The declcal / declcalc0 / c0approx payloads already
   carry `pstar_implied_05[_c0]` computed at the exact inverse, and the findings scripts read those columns.
   "Read the value from the calibration output" has no target there: the payloads predate the settable
   columns, and the task forbids re-running the campaigns. The options:
   - (a) derive `pstar_settable` from each row's stored `kappa_hat` with the package's settable helper
     `.fs_decl_settable()` (`R/fs_declaration_calibration.R:205`; one derivation, the package's);
   - (b) keep the column, relabel it plainly as the exact-cutoff inverse, and add the settable value beside it;
   - (c) leave the committed findings as historical records and fix only the run scripts for future runs.

   **Recommendation: (a) for the findings scripts plus the run-script fix.** The package's helper becomes the
   only derivation, and no payload changes.
2. **Historical verification and theory code** (`dev/glm-continuous-sims/verification/*`,
   `actg175/binary_methods/fdr_*`, `actg175/continuous/scripts_mdf1/reselection_check.R`,
   `quarto/resampling/*`). These reproduce the pre-alignment MR and the theory document's exact-cutoff
   exposition, and their recorded outputs are tied to committed reports. Rewriting them changes what those
   records mean. **Recommendation:** leave them, and add a one-line header comment to each saying it encodes
   the exact cutoff that `7713942e` / `96f84ad8` replaced. The theory `.qmd` is a manuscript-adjacent
   document, which is your call.
3. **The OC predictor `R/` change** (§2): a separate task with its own test update.

## Gates

- Gate 1 (no silent emptiness; value matches calibration): **not run**, no consumer was fixed.
  `report_values_c0.R:12` is identified from source as a live silent-vanish site.
- Gate 2 (`devtools::test()`), Gate 3 (`R CMD check --as-cran`): **not run**, no file under `R/` or `tests/`
  changed.
- Gate 4 (no `R/` change): **holds**; `git status --short -- R/` is empty.

## Post-conditions

1. Every §1 hit listed with file, line and disposition: yes. Every disposition is "not fixed: stopped at
   inventory", with the reason above.
2–5. Not met; they require the fixes.
6. No `R/` file modified: met.
7. Nothing written to `fs-glms-interpretable`: met.
