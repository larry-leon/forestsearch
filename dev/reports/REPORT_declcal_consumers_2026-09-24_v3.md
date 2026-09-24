# REPORT — Declaration-calibration consumers v3: complete after the dispositions; all gates pass

- **Task:** `dev/tasks/TASK_declcal_consumers_2026-09-24_v3.md` (`dd3b8331`), resumed under
  `dev/tasks/DISPOSITIONS_declcal_consumers_2026-09-24.md` (`23c6af1b`).
- **Outcome: COMPLETE. Gates 1, 2 (read under classes (a)–(d)), 3 and 4 pass.** Outside `R/`, no run script
  derives the threshold. The first pass stopped at Gate 2 by the `R/` rule; that pass is kept below as history
  (§§1–8). The dispositions and the final gates are in §9.

## 0. Record

| item | value |
|---|---|
| HEAD | `dd3b8331` (task doc v3); before it `fdd2f500` |
| R, platform | R 4.6.1 (2026-06-24), x86_64-pc-linux-gnu, pop-os |
| Installed forestsearch | 0.3.5.9000, built 2026-09-24 01:50:39 UTC. `R/` is unchanged since that build (`7713942e` is the last `R/` change). |
| `.fs_decl_settable` in namespace | TRUE |
| Untracked before the task (never staged) | `actg175/binary_020/mr_or_harm/…_redes_d5000/`, `…_relaunch_d5000/`, `actg175/binary_020/smoke_redes.html`, `smoke_relaunch.html`, `gbsg_020/scripts_dinamr/logs/nullmr_findings.err` |
| Clock | 04:34 UTC start, stopped about 05:00 UTC (well inside the 2 h abort) |

## 1. Inventory dispositions at the stop (one row per site in `REPORT_declcal_consumers_2026-09-23.md`; history, §9.1 gives the final class of each site)

Kinds: RF = reads a removed or changed field; IP = own implied p*; TH = own threshold; HC = hard-coded value;
OK = not a screen quantity.

| file | line (at 4ef26719) | kind | class | what changed | Gate 1 |
|---|---|---|---|---|---|
| `gbsg_020/scripts_dinamr/declcal_run.R` | 214-215 | HC/TH | run script | `z_exact`/`z_round` constants removed. `z_round <- fs_declaration_calibration(mr)$z_pstar` per replicate (now :348-351), with an assert that the calibration's digits equal the fit's. | PASS |
| same | 346-347, schema 207 | IP | run script | `pstar_implied_05/10` replaced by `pstar_settable_05/10 <- .fs_decl_settable(kappa_hat, digits)$p_star` (:353-354). Schema adds `z_pstar`. | PASS |
| same | 348 | TH (exact) | run script | `alpha_FW_hat_1645` dropped (the exact-cutoff column, per §2a). `alpha_FW_hat_1621 <- dcal$fw_size`. | PASS |
| same | 363-369 | TH | run script | Inline `round(rate, digits) >= p_star` replaced by `post %in% dcal$admitted_pstar` (:369). `declared_conv_exact` and `n_band` (exact-cutoff comparators) dropped. Meta `z_exact`/`z_round` replaced by `threshold_source`. The RATES line no longer prints `exact`. | PASS |
| same | 298 | RF (schema) | run script | none needed (the payload's own `pconsistency_digits`) | n/a |
| `declcalc0_run.R` | 230-231, 370-373, 390, 394-395, 413-419, 223 | HC/TH/IP | run script | Same as `declcal_run.R`. Per c0: `pstar_settable_05_<c0>`; `fw_1645_<c0>` dropped; `fw_1621_<c0>` at `dcal$z_pstar`. `declared_conv_exact` removed from the identity-gate columns (it would have compared `NULL` and silently passed). | PASS |
| `declcal_c0approx_run.R` | 237-238, 369-372, 389, 393, 403-409, 224 | HC/TH/IP | run script | Same as `declcalc0_run.R` (this driver has no `fw_1621_<c0>`) | PASS |
| `declcal_findings.R` | 56, 79-85 | IP (read) | run script | Reads of the stored exact-cutoff inverse replaced by `.fs_decl_settable_table(kappa_hat_05, pconsistency_digits)`: settable `p*` at the fit's digits, the share with none ≤ 1, and the package's finer pair `pstar_fine [digits_fine]`. Table 7.3(b) now gives `pstar_fine` quantiles. | PASS |
| same | 111 | RF (schema) | run script | none needed | n/a |
| `declcalc0_findings.R` | 29-34 | IP | run script | `ptxt()` and the footnote now use the settable pair. The footnote no longer states the exact inverse. | PASS |
| same | 103, 111 | IP (read) | run script | Reads replaced by the settable pair derived from `kappa_hat_05[_c0]` | PASS |
| same | 207-214 | IP + HC | run script | Settable/`pstar_fine` distribution. `qnorm((1 + 0.895) / 2)` replaced by `z_exec(cl)`, which reads the payload's recorded `meta$z_round` (else `results$z_pstar`). That value equals the package's `z_pstar` (1.62108225085241; `all.equal` TRUE). | PASS |
| `declcal_c0approx_findings.R` | 111 | IP | run script | `2 * pnorm(median kappa) - 1` replaced by `pset1()`: settable `p*`, or "none <= 1 (fine p [digits d])" | PASS |
| `dev/verification/report_values_c0.R` | 12 | RF (silent vanish) | run script | `f(t05$pstar_implied[i])` replaced by `f(t05$pstar_settable[i], 2)`; the header column is renamed. | PASS |
| same | 13, 22 | RF (changed) | run script | Checked; no change needed. `fw_size` is the package's rounded value, and `.fw5_target` is already at the rounded cutoff (test fixture). | PASS |
| same | 25 (not in the inventory) | TH (exact) | run script | `mean(m > qnorm(0.95))`, an exact-cutoff fw compared against the rounded `.fw5_target`, replaced by `mean(m > d5$z_pstar)`. The discrepancy is now 0.000233. | PASS |
| `dev/verification/report_values.R` | 12-13, 29-38 | RF (changed) | run script | `fwt <- .fw_target(.rho8)` (the exact z = `qnorm(0.95)` target) replaced by `1 - .phi2(d8$z_pstar, .rho8)`. T8 fw disc: gaussian 0.000422, helper 0.000602, poisson −0.000208. | PASS |
| same | 15 | RF (changed) | run script | Label only: "z_pstar (rounded rule)" | PASS |
| same | 49-50 | RF + TH | run script | The inline admission rebuild is replaced by `intersect(s, dg$admitted_pstar)`; "GLM relabel reproduces: TRUE". | PASS |
| `gbsg_app_null/pstar_grid_findings.R` | 103 | HC | run script | `0.651` replaced by `fw_app_075 <- 0.6672`, with provenance. The value was re-checked as `mean(Mstar_c0 > z_pstar)` on the alignment report's capture (`~/Downloads/declcal_rounding_baseline_2026-09-23.rds`): 0.6672. | PASS |
| `quarto/resampling/fdr_family_multiplier.R` | 169 | TH | archived (theory companion, "PROTOTYPE / SKETCH") | **header not added (stopped)** | — |
| `quarto/resampling/consistency_resampling_theory.qmd` | 904, 1136 | TH | archived (theory note) | **header not added (stopped)** | — |
| `quarto/resampling/consistency_resample.R:164`, `validate_consistency_adjusted.R:92`, `gbsg_consistency_demo_standalone.R:33` | — | OK (Pcons definition) | see §5.2 | not touched. These are Gate 2 pattern hits that fit none of (a)–(c). | — |
| `actg175/binary_methods/fdr_mr_inference.R` | 112, 199, 280 | TH | **run script** (sourced by the committed `actg175_binary_m1*_adjusted*.qmd` simulation documents) | **not fixed: needs the §5.1 decision** | — |
| `actg175/binary_methods/fdr_family_multiplier.R` | 169 | TH | **run script** (same) | **not fixed: §5.1** | — |
| `actg175/continuous/scripts_mdf1/reselection_check.R` | 54 | TH | **run script** (verification) | **not fixed: §5.1** | — |
| `actg175/continuous/oc_wrapper_verification.qmd` | 45 | TH (mirrors the OC gate) | archived (rendered verification note). Its exact cutoff mirrors `R/fs_oc_predict.R:279`, which is out of scope. | **header not added (stopped)** | — |
| `dev/glm-continuous-sims/verification/mr_mechanism_A1.R` | 203, 267 | TH | **run script** (verification; `…_earlystop.R` and `…_b1_price.R` source it) | **not fixed: §5.1** | — |
| `…/mr_mechanism_probe_superset.R`, `…_probe_restrict.R` | 26 | TH | **run script** (verification) | **not fixed: §5.1** | — |
| `dev/glm-continuous-sims/sigma_d_diagnostic_2026-08-29.R` | 99 | OK (printed constant) | archived (dated diagnostic) | **header not added (stopped)** | — |
| `dev/identifier-alignment/code_theory_audit.qmd` | 457 | TH (prose) | archived | **header not added (stopped)** | — |
| `R/fs_oc_predict.R:279`, `R/fs_oc_grid.R:582, 383` | — | TH | out of scope (separate decision) | untouched. **Record entry (dispositions §1):** the forward map from (p\*, digits) to the rounded z cutoff (the suggested `.fs_z_eff(p, digits)`) is the map the six archived scripts would need and the same map the OC gate lacks. `R/fs_oc_grid.R:383` (the root-bracket padding `lo - zp * max(se_g)`, exact `zp`) is the third OC site. This is a record, not a task. | — |
| `tests/testthat/test-declaration-calibration.R` | 111-124 (test 1) | TH | test (§3) | fixed (§3 below) | Gate 3 form PASS |
| `test-declaration-calibration.R:53, 155-158`, `test-declaration-c0.R:230`, `test-mr-admission-rounded.R`, `test-declcal-rounding-alignment.R` | — | HC | test oracle (b) | untouched | — |
| `test-fs-oc-predict.R:214` | — | TH | out of scope (c) | untouched | — |

**How Gate 1 was checked.**

- **Drivers.** For each driver, the HEAD copy and the fixed copy replayed the same replicates at the committed
  cells' knobs, in a scratch session with installed forestsearch. That was reps 2 and 32 of C3 at `B_cal` 500
  for `declcal`/`declcalc0` (rep 2 is one with `n_band > 0` in the committed payload), and reps 1 and 2 of B5 at
  `B_cal` 2000 for `c0approx`. **Every shared column is identical**, including `declared_conv`,
  `n_admitted_conv` and `alpha_FW_hat_1621`: 29, 65 and 49 shared columns. Every new field has length 1.
  `pstar_settable_*` is `identical()` to `.fs_decl_settable(kappa_hat, digits)$p_star`, and `z_pstar` equals
  the calibration's own `z_pstar`.
- **Findings scripts.** The three findings scripts ran on the committed payloads in a scratch mirror, so no
  committed output was overwritten. `declcal_c0approx_findings.R` reads per-cell capture files that are not
  committed; I rebuilt them from the committed `results/declcal_c0approx_res.rds` in scratch. All exited 0,
  and no markdown row has an empty cell. The diff against the committed logs is confined to the p* columns.
  Spot checks match the package: A1 median `pstar_fine` 0.99961; B5 c0 0.70 settable median 0.97, none ≤ 1
  in 0 reps.
- **The two `dev/verification` scripts** ran to exit 0 with no empty field.

**What the settable column now says.** For every unshifted declcal cell, κ̂₀.₀₅ is about 3.5–3.6. No p* ≤ 1
reaches that at `pconsistency.digits = 2`: the largest cutoff there is z = 2.807, at p* = 1.00. So the settable
p* is "none" in 100% of replicates. The package's finer pair puts it at about 0.9996 with 3–6 digits. The old
column printed 0.9996 to five decimals as though it were settable at the executed digits; it was not.

## 2. Gate 2 at the stop — FAILED (history; the final Gate 2 is §9.5)

The inventory's patterns, re-run over `*.R *.qmd *.Rmd *.sh` outside `R/` and `man/`, give 43 hits. The
patterns were `pstar_implied`, `0.651`, `0.6508`, `0.9936`, `2 * pnorm(.) - 1`, `qnorm((1 + .`,
`round(rate, digits) >= p_star` and `0.895`. Reads of the calibration's own `fw_size` / `z_pstar` /
`pconsistency_digits` are the fix itself, so I counted them as reads, not derivations.

- **(b) test oracles: 22.** `test-declaration-calibration.R` (6, including test 1's rounded oracle and band
  assertion), `test-declaration-c0.R` (2), `test-declcal-rounding-alignment.R` (9), `test-mr-admission-rounded.R`
  (5); `test-fs-oc-predict.R` is counted under (c).
- **(c) OC pin: 1.** `test-fs-oc-predict.R`.
- **(a) archived with a header line present: 0.** No header was added because the task stopped. Six hits would
  become (a) once headed: `consistency_resampling_theory.qmd` (2), `quarto/resampling/fdr_family_multiplier.R`,
  `oc_wrapper_verification.qmd`, `code_theory_audit.qmd` and `sigma_d_diagnostic_2026-08-29.R`.
- **In no class: 12 hits.** Nine are in six run scripts:
  - `actg175/binary_methods/fdr_mr_inference.R` (3) and `fdr_family_multiplier.R` (1)
  - `reselection_check.R` (1)
  - `mr_mechanism_A1.R` (2), `mr_mechanism_probe_superset.R` (1) and `mr_mechanism_probe_restrict.R` (1)
  - and three more are the Pcons-definition lines `quarto/resampling/consistency_resample.R:164`,
    `validate_consistency_adjusted.R:92` and `gbsg_consistency_demo_standalone.R:33`
- **Text-only hits: 2.** `declcal_findings.R:71` is a label describing the stored `n_band` column as "rate in
  [0.895, 0.900)". `quarto/GuoHe/guohe_sec52_sim.R:102` is a data vector that contains 0.895, a false positive
  of the `0.895` pattern.

Sentence of record: **Outside `R/`, run scripts still derive the threshold: 9 hits in 6 verification or
prototype run scripts, plus 3 Pcons-definition hits in no class. Remaining hits: 0 (a) archived
(6 pending headers), 22 (b) test oracles, 1 (c) OC pin.** None of the 9 fixed consumers has a derivation
hit (`declcal_findings.R:71` is a text label), and `pstar_implied` appears in no consumer (post-condition 2 holds).

## 3. Test 1

- **Old expectation.** On `.fit_on` (p* 0.90), `s[rate >= p_star & beta_hat >= c_screen]` equals the relabelled
  set: the unrounded rule. It passed on current code only because `.fit_on` has no screened candidate with
  Pcons in [0.895, 0.90). The top screened rates are 0.92806, 0.88180, 0.86597 and 0.86457, so the test could
  not tell the two rules apart.
- **New expectation.** `s[round(rate, digits) >= p_star & beta_hat >= c_screen]` equals the relabelled set, which
  equals `admitted_current`. `digits` is the fixture's `pconsistency.digits`, else 2, and the oracle is written
  out without `.fs_pcons_eff()`. The test runs on both `.fit_on` and a new band fixture. On the band fixture it
  asserts, with a message saying the test loses its power without one, that at least one screened candidate is
  admitted by the rounded rule and not by the unrounded one.
- **Fixture change (recorded).** `.gb_fit()` gains `pstar = 0.90` (the default, so every other fixture is
  unchanged). The new fixture is `.fit_band <- .gb_fit(<field on>, pstar = 0.87)`. There, `q5.0 & q18.1`
  (Pcons 0.86597) rounds to 0.87 and is admitted, while the unrounded rule rejects it. It adds about 8 s.
- **Run.** `devtools::test(filter = "declaration-calibration")`: `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 82 ]`,
  under 10 min. This ran in §3 to validate the test. Gate 3 is formally after the Gate 2 stop, so it is
  reported here rather than as a passed gate.
- **Old-rule failure.** Scratch copy outside the repo, with the oracle swapped to `rate >= dc$p_star`:

```
── 1. Failure ('test-declaration-calibration-OLDRULE.R:136:5'): 1: with kappa_ha
Expected `rounded` to have the same values as `relabelled`.
Actual: "q4.0 & q18.1", "q18.1 & q20.0"
Expected: "q4.0 & q18.1", "q5.0 & q18.1", "q18.1 & q20.0"
Absent: "q5.0 & q18.1"

── 2. Failure ('test-declaration-calibration-OLDRULE.R:139:7'): 1: with kappa_ha
the band fixture has no screened candidate on which round(Pcons, digits) >= p* and Pcons >= p* disagree; without one, test 1 cannot tell the rounded rule from the unrounded one and loses its power
```

## 4. Gates 3–4 and post-conditions

- **Gate 3:** the test file passes (§3). This is not a gate result, because Gate 2 failed first.
- **Gate 4:** holds. `git status --short -- R/` is empty.
- **Post-conditions.**
  - PC1 (every site dispositioned): no; 6 run scripts are pending the §5 decision.
  - PC2 (no consumer references `pstar_implied`): met.
  - PC3 (Gate 1): met for every fixed consumer.
  - PC4 (Gate 2): **not met**.
  - PC5 (test 1): the band fixture and the old-rule failure are shown; Gate 3 was not formally reached.
  - PC6 (no `R/` change): met.
  - PC7 (catalogue files unmodified): met.
  - PC8 (nothing written to `fs-glms-interpretable`): met.

## 5. Decisions for Larry (answered by the dispositions; see §9)

1. **The six remaining run scripts** (`actg175/binary_methods/fdr_mr_inference.R` and `fdr_family_multiplier.R`,
   `reselection_check.R`, `mr_mechanism_A1.R`, and the two probes). Options:
   - (a) An `R/` helper, `.fs_z_eff(p_star, digits)`, meaning `qnorm((1 + .fs_pcons_eff(p, d)) / 2)`. It would
     be used at `fs_mr_inference.R:678`, `fs_declaration_calibration.R:579`, and inside `.fs_decl_settable()`.
     Scripts would call it with no local formula. This is the natural companion to the OC-gate fix and belongs
     in that task.
   - (b) Add a Gate 2 class for "maps `.fs_pcons_eff()` to z", then fix the scripts with the local
     `qnorm((1 + .fs_pcons_eff(…)) / 2)`. Most need a `digits` argument threaded through the prototype
     signatures.
   - (c) Reclassify them as archived and add header lines only. They verify pre-alignment behaviour, and their
     outputs are tied to committed reports: `REPORT_stopB_md_harm_grid`, the actg175 binary-methods documents
     and the MDF1 Stage 1b check. "Fixing" them changes what those records mean. This was v1's recommendation
     §3.2.

   **Recommendation: (c) now, and (a) in the OC-gate task.** No future run in the current plan invokes these
   scripts. A future re-use would pick up (a).
2. **The Pcons-definition hits** (`2 * pnorm(delta / sigma_D) - 1` in three `quarto/resampling` scripts) are not
   a threshold. The inventory marked them OK, but Gate 2 has no class for them. **Recommendation:** add a class
   "(d) Pcons definition, not a threshold", or treat them as archived theory companions with headers.
3. **My judgment call under §2a (please confirm).** The drivers no longer record the exact-cutoff comparators
   `alpha_FW_hat_1645`, `declared_conv_exact`, `n_band` and `fw_1645_<c0>`, only `pstar_implied_*`. Recording
   them needs a local exact threshold, which is a Gate 2 hit. The findings scripts still read those columns from
   the committed payloads, and they are tied to those payloads. A future re-run of the findings on a new payload
   would find those columns absent. The drivers' own RATES and identity-gate lines were adjusted so they don't
   silently read `NULL`.
4. **Two exact comparators left in the findings scripts** (`declcal_findings.R:107`
   `kappa_hat_10 < qnorm(0.95)`; `declcalc0_findings.R:213` `pre_ex <- max_T_pre >= qnorm(0.95)`). They are
   labelled "1.6449" and pair with the stored exact columns. They are not an inventory pattern, and I did not
   change them (scope). Tell me if they should move to the executed cutoff.

## 6. Working-tree changes at the stop (history; §9.7 has the commits)

Modified, not staged:

- **Run-script fixes:** `dev/verification/report_values.R`, `report_values_c0.R`; `gbsg_020/scripts_dinamr/`
  `declcal_run.R`, `declcalc0_run.R`, `declcal_c0approx_run.R`, `declcal_findings.R`, `declcalc0_findings.R`,
  `declcal_c0approx_findings.R`; `gbsg_app_null/pstar_grid_findings.R`.
- **Test fix:** `tests/testthat/test-declaration-calibration.R`.
- **This report.**

Proposed commits once §5 is settled, in the task's order: run-script fixes → archived headers → test fix →
report.

## 7. Catalogue entries owed to the next gbsg_020 closeout (not edited here)

- **Fixed drivers:** `scripts_dinamr/declcal_run.R`, `declcalc0_run.R`, `declcal_c0approx_run.R`. Schema:
  `pstar_settable_*` and `z_pstar` are added; the exact-cutoff columns are dropped.
- **Fixed read-outs:** `scripts_dinamr/declcal_findings.R`, `declcalc0_findings.R`,
  `declcal_c0approx_findings.R`. Their p* columns now report the settable pair.
- **Archived outputs, header line added (§9.1):** `scripts_dinamr/logs/declcal_findings.txt`,
  `logs/declcalc0_findings.txt`, `logs/declcal_c0approx.txt`. The `.rds` payloads
  (`results/declcal_*`, `declcal_c0approx_res.rds`, `declcal_findings.rds`, `declcalc0_findings.rds`) cannot
  carry a header line, so they are annotated through the catalogue entry only.
- **Outside gbsg_020 (their own closeout rule):** `gbsg_app_null/logs/pstar_grid_findings.txt`, header line
  added.

## 8. Side issues (not acted on)

- Committed `.md` reports quote the exact-cutoff implied p* or 0.651. Examples:
  `REPORT_gbsg_app_null_declaration_2026-09-23.md`, `REPORT_gbsg_pstar_grid_2026-09-23.md`, and the declcal
  campaign reports. The inventory's search covered `*.R/*.qmd/*.Rmd/*.sh` only, so they are outside this task.
- In `report_values_c0.R`, an unreachable settable p* prints as `NA`. The package's own print says `none`.

## 9. Dispositions (resumption, 2026-09-24)

Larry's answers to §5 are in `dev/tasks/DISPOSITIONS_declcal_consumers_2026-09-24.md`, committed as `23c6af1b`.
The resumption started at 05:07 UTC and finished at about 05:15 UTC, inside the 1 h abort. `R/` is not touched.

**Background shell.** The one still running at the stop was this task's read-only
`find / -name declcal_c0approx_B1_res_1_20.rds`, which looked for the uncommitted per-cell capture files. It
was stopped. It wrote nothing into the repo, and nothing from it is staged.

### 9.1 Final class of each site

- **Fixed run scripts (9 files, unchanged from §1):**
  - the three drivers `declcal_run.R`, `declcalc0_run.R` and `declcal_c0approx_run.R`;
  - the three findings scripts `declcal_findings.R`, `declcalc0_findings.R` and `declcal_c0approx_findings.R`;
  - `dev/verification/report_values.R` and `report_values_c0.R`;
  - `gbsg_app_null/pstar_grid_findings.R`.

  The findings scripts carry the §9.3 and §9.4 additions.
- **Class (a), archived, with a header line and no code change (dispositions §1).** These eight files carry the
  line "Verifies pre-alignment (exact-cutoff) behaviour, before the 2026-09-23 alignment (96f84ad8, 7713942e);
  not valid against the aligned package. See dev/tasks/TASK_declcal_consumers_2026-09-24_v3.md.":
  - `actg175/binary_methods/fdr_mr_inference.R` and `fdr_family_multiplier.R`
  - `actg175/continuous/scripts_mdf1/reselection_check.R`
  - `dev/glm-continuous-sims/verification/mr_mechanism_A1.R`, `_probe_superset.R` and `_probe_restrict.R`
  - `_probe_earlystop.R` and `_b1_price.R`, which read A1 up to its `## --- run ---` marker, so the added line
    does not shift what they evaluate

  `git diff dd3b8331 --numstat` shows `1 0` (one line added) on each of the eight.
- **Class (a), archived, with the v3 §2b line "Values computed at the exact cutoff, before the 2026-09-23
  alignment (96f84ad8, 7713942e). See dev/tasks/TASK_declcal_consumers_2026-09-24_v3.md."** In `.qmd` files it
  sits as an HTML comment immediately after the YAML front matter.
  - Theory and diagnostic files: `quarto/resampling/fdr_family_multiplier.R`,
    `consistency_resampling_theory.qmd`, `oc_wrapper_verification.qmd`, `code_theory_audit.qmd` and
    `sigma_d_diagnostic_2026-08-29.R`.
  - Committed outputs of the fixed findings scripts: `gbsg_020/scripts_dinamr/logs/declcal_findings.txt`,
    `declcalc0_findings.txt`, `declcal_c0approx.txt` and `gbsg_app_null/logs/pstar_grid_findings.txt`. The
    `.rds` payloads cannot carry a line; §7 lists them.
- **Class (d), Pcons definitions, with no change and no header.** Each of these three lines computes the
  consistency proportion from T and compares it with no threshold:
  - `quarto/resampling/consistency_resample.R:164`
  - `quarto/resampling/validate_consistency_adjusted.R:92`
  - `quarto/resampling/gbsg_consistency_demo_standalone.R:33`
- **Out of scope:** the OC gate. The record entry is in the §1 table (`R/fs_oc_grid.R:383` is the third site).

### 9.2 Dropped exact-cutoff columns: confirmed (dispositions §3)

The drivers stay as fixed. The committed payloads keep `alpha_FW_hat_1645`, `declared_conv_exact`, `n_band` and
`fw_1645_<c0>` as records.

### 9.3 Column checks in the findings scripts

Each findings script now defines `need_cols()` and the full list of columns it reads, at the top. It checks
every payload as soon as it is read, and stops with a message naming any missing column:

- `declcal_findings.R` checks `$results` and `$aux`.
- `declcalc0_findings.R` checks `$results` (including every per-c0 column), `$aux`, and the committed declcal
  reference payloads.
- `declcal_c0approx_findings.R` checks the per-cell captures and the references.

**Demonstration.** I made a scratch copy of `declcal_bnull_A1_res_1_2000.rds` with `declared_conv_exact`
removed, in a scratch results directory whose other payloads are symlinks to the committed ones. Running
`declcal_findings.R` there exits with status 1 and writes no output:

```
Error: ../results/declcal_bnull_A1_res_1_2000.rds $results lacks column(s) read by this script: declared_conv_exact
Execution halted
```

### 9.4 `qnorm(0.95)` / `1.6449` / `1.644854` in the nine fixed consumers (dispositions §4)

| file:line | use | action |
|---|---|---|
| `declcal_findings.R:107` (now :128) | **threshold**: `kappa_hat_10 < qnorm(0.95)` | Now `kappa_hat_10 < z_exec(cl)`, the payload's recorded executed cutoff. The column is relabelled "reps with kappa_hat_10 < executed cutoff z_pstar". The counts are unchanged in every cell. |
| same, :102 | column label "< 1.6449" | relabelled (above) |
| `declcalc0_findings.R:213` (now :233) | **threshold**: `pre_ex <- max_T_pre >= qnorm(0.95)` | Now `>= z_exec(cl)`. Its paired estimator in the same column block changes from `fw_1645_c0` to `fw_1621_c0`, because an fw at the exact cutoff would not measure a rate at the executed cutoff. The header and the free-check (a) note are relabelled. Example: B1 at c0 0.70 changes from 0.3329 vs 0.1220 to 0.3465 vs 0.1295. |
| `declcalc0_findings.R:209` | column label "(max_T_pre >= 1.6449)" | relabelled "(max_T_pre >= z_pstar)" |
| `declcalc0_findings.R:253` | label describing the stored `declared_conv_exact` ("max_T_post >= 1.6449") | Left as is: it describes a committed record column, not a threshold use. |
| `declcal_c0approx_findings.R:124` | label defining the stored `fw_1645` column ("1{M*(c0) > 1.6449}") | Left as is, for the same reason. |
| `fw_1645` / `alpha_FW_hat_1645` column names in the findings scripts | reads of stored record columns | Left as is (dispositions §3: the payloads keep them as records). |
| driver comments `:216/:232/:239` naming the dropped `alpha_FW_hat_1645` | comment | left as is |
| Wilson intervals (`qnorm(0.975)`) and `z975` in the drivers | interval critical values | untouched |

No threshold use of `qnorm(0.95)`, `1.6449` or `1.644854` remains in the nine fixed consumers.

### 9.5 Gates

- **Gate 1 (changed files only).** The three findings scripts ran again on the committed payloads in the scratch
  mirror. All exited 0, with no zero-length or empty field. Against the pre-dispositions run, only the §9.4
  rows and labels differ.
  - The recorded cutoff `meta$z_round` is `identical()` to the package's `z_pstar` (1.62108225085241) in the
    A1, C3 and declcalc0 B1 payloads.
  - An independent recompute of the B1 pre-family rate at that cutoff gives 0.1295, as printed.
  - The drivers, `dev/verification` scripts and `pstar_grid_findings.R` did not change, so they were not re-run.
- **Gate 2 (final), under the inventory's own patterns.** The patterns are `pstar_implied`, `0.651`, `0.6508`,
  `0.9936`, `2 * pnorm(.) - 1`, `qnorm((1 + .` and `round(rate, digits) >= p_star`. §2 had added `0.895`, which is
  not an inventory pattern; its only non-test hits were a label and a data vector, so it is dropped here. The
  search covers `*.R *.qmd *.Rmd *.sh` outside `R/` and `man/`: **36 hits, every one classified.**
  - (a) 15: `fdr_mr_inference.R` 3, `mr_mechanism_A1.R` 2, `consistency_resampling_theory.qmd` 2, and 1 each in
    `actg175/.../fdr_family_multiplier.R`, `quarto/resampling/fdr_family_multiplier.R`, `reselection_check.R`,
    the two probes, `oc_wrapper_verification.qmd`, `code_theory_audit.qmd` and `sigma_d_diagnostic`. Each file
    has its header line.
  - (b) 17: `test-declcal-rounding-alignment.R` 7, `test-mr-admission-rounded.R` 5,
    `test-declaration-calibration.R` 4 (including test 1's rounded oracle and band assertion) and
    `test-declaration-c0.R` 1.
  - (c) 1: `test-fs-oc-predict.R`.
  - (d) 3.

  Reads of the calibration's own outputs (`fw_size`, `z_pstar`, `pconsistency_digits`) appear only in the nine
  fixed consumers, as the fix itself.

  **Sentence of record: Outside `R/`, no run script derives the threshold; remaining hits: 15 (a) archived,
  17 (b) test oracles, 1 (c) OC pin, 3 (d) Pcons definitions.**
- **Gate 3.** The test file did not change after §3's run (`[ FAIL 0 | WARN 0 | SKIP 0 | PASS 82 ]`), so it was
  not re-run, per the dispositions.
- **Gate 4.** `git status --short -- R/` is empty.

### 9.6 Test 1, recorded only

- `test-declaration-calibration.R` has **no** `skip_on_cran()` or other skip. Its top-level fixtures (the GBSG
  fits, including `.fit_band`, and the `B = 200000` field) run on CRAN.
- The band fit (`.fit_band`, p* 0.87, 2000 draws, field on) takes **5.8 s** wall clock on pop-os; a comparable
  field-off fit takes 5.1 s.

### 9.7 Post-conditions and commits

**Post-conditions.**

- **v3 PC1–PC8 hold**, with Gate 2 read under classes (a)–(d):
  - every site is dispositioned;
  - no consumer references `pstar_implied`;
  - Gate 1 holds;
  - Gate 2 passes;
  - Gate 3 and test 1's band and old-rule failure are shown in §3;
  - no `R/` change;
  - `gbsg_020/status_curated.md` and `current_status.md` are unmodified;
  - nothing was written to `fs-glms-interpretable`.
- **Each of the eight §1 files differs from `dd3b8331` by the one added header line.**
- **Every findings script has the column check,** and the check stops on the scratch payload (§9.3).
- **No threshold use of the exact cutoff** remains in the nine fixed consumers (§9.4).
- **Only named paths are staged,** and the pre-existing untracked files stay unstaged.

**Commits,** by explicit paths, in this order: run-script fixes; archived-file header lines; test fix; this
report.
