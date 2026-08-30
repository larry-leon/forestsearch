# CC task — breadth, stage 1: the interaction ladder on the MD40 DGM, and a locked forecast

**File:** `dev/tasks/cc_task_oc_breadth_ladder_2026-08-30.md` · **Issued:** 2026-08-30 by chat, commissioned by Larry
**Repo:** `~/Documents/GitHub/forestsearch` (pop-os) · **Branch:** `feature/glm-extension`
**Transport:** `~/Downloads` → `dev/tasks/`, committed first (§1).
**Reads:** `REPORT_oc_wrapper_confs_2026-08-30.md` + `oc_wrapper_grid_corrected_2026-08-30.rds` and its script `oc_wrapper_grid_corrected_2026-08-30.R` (the corrected MD40 run: family, settings, seed, the null inversions); `REPORT_residual_2026-08-30.md` §2 (the `se_scaled / se_direct` diagnostic, reused here); the MD40 driver(s) under `quarto/simulations/actg175/continuous/` for how the DGM was built.

**Larry's decision, 2026-08-30 — breadth by effect regime, in two stages.** The wrapper has been validated on one DGM in a regime where the screen saturates: the planted interaction (13.7) sits under a winner's-curse inflation of ≈ 75, so no `c1` separates power from type-I error, and every exact agreement was of a saturated quantity. The regime where the method would be recommended — a harm large enough that some `c1` gives 80% power at the null's 5% type-I — is untested and is where the wrapper's prediction is qualitatively different: a crossover, not a ceiling.

**Stage 1 (this task) is analytic only and ends in a forecast.** Build a ladder of DGMs identical to MD40 except for `k_inter`, run the wrapper on each, and from the ladder name the smallest harm at which the search has 80% power at the null's 5% threshold, and the `c1` that delivers it — together with every operating characteristic the wrapper predicts at that point. **The committed report is the forecast, locked before any replicate is drawn.** Stage 2 — one measured cell at that harm — is a separate go/no-go and is **not** part of this task.

**Why this is cheap — verified from source, 2026-08-30.** `k_inter` enters `generate_glm_dgm()` at one line, `mu1[in_Q] <- mu1[in_Q] + beta_inter` (L462). `flag_harm` is defined on the source data before sampling (L303–320), `df_super` is drawn under `set.seed(seed)` (L395–397), `mu0` is the baseline prediction (L454), and the residual SD is `sigma(fit_base)` (L489–491) — none depends on `k_inter`. `fs_oc_family_enumerate()` builds membership from `df_super`'s covariates and `flag_harm` alone (`fs_oc_family.R` L278–375), and the scale enters through `bint = m_Q − tauQc` (L253–256) with `seQ1000` from `V_eff[Q]`, and Q's treatment effect is constant within Q. **So across the ladder `lab`, `Pg`, `ovl`, `M`, `sens_g`, `spec_g` and `se_g` are bit-identical and only `beta_g = tauQc + bint·PQg` moves.** That is both the reason a rung costs minutes and the gate in §3.

---

## ⚠ CATEGORY

**No `R/` change of any kind.** No new export, no default change, no edit to any package file, driver or document. Writes: scratch scripts, their `.rds`/`.log`, and the report, under `dev/glm-continuous-sims/`. Plus this task document.

**Compute:** §2's reproduction gate ≈ 10 min (one draw set). §3's ladder: six rungs, one draw set each under the production gate, ≈ 10–15 min per rung. §5's confirmatory rung under both gates ≈ 25 min. §4's diagnostic: `fs_dgm_scale()` over M regions per rung, seconds to a minute each. §6's sensitivity: one draw set. **Estimate 2–3 hours.** No replicates, no simulation study, no renders.

**Unattended.** Gates stop the task; never ask, never work around.

---

## 1. Provenance and first commit — GATE

Standard block (`hostname; pwd; git rev-parse --abbrev-ref HEAD; git rev-parse --short HEAD; git status --porcelain; git log --oneline -6`) plus the installed version (expect **0.3.1**; the dispatch change touches nothing the wrapper calls, and §2 proves it). *GATE:* branch `feature/glm-extension`, clean tree, `2276ee10` in the log. Copy this document into `dev/tasks/` and commit it alone; the `~/Downloads` stem arrives with hyphens stripped, so report both names.

## 2. Reproduce MD40 exactly — GATE

Everything in this task is built on the MD40 DGM being rebuilt bit-for-bit. Establish from source, quoting it:

1. **How the driver built MD40** — the data preparation from `speff2trial::ACTG175` (arms, outcome, treatment coding, every derived column), and the exact `generate_glm_dgm()` or `calibrate_glm_interaction()` call: `factor_vars`, `continuous_vars`, `subgroup_vars`, `subgroup_cuts`, `k_treat`, `n_super`, `seed`, `k_inter` or `target_effect`, `adverse_outcome`. Say which driver you read and whether the three drivers (alt n = 500, alt n = 700, null n = 500) agree on the build.
2. **The `forestsearch_args`** the corrected run used, from `oc_wrapper_grid_corrected_2026-08-30.R`, and the corrected `.rds`'s `settings` — `c1`, `c2`, `pconsistency`, `draws`, `block`, `seed` (expect 20260825). Use these, unchanged, for every rung.

Then, in `dev/glm-continuous-sims/oc_breadth_ladder_2026-08-30.R` (parts as you see fit; log each):

- Rebuild MD40 **exactly as the driver did.** Record `dgm$model_params$k_inter` — call it `k40` — and `fs_dgm_scale(dgm)`'s `m_tau` for Q and Qc.
- If the driver used `calibrate_glm_interaction()`, also rebuild it directly with `generate_glm_dgm(..., k_inter = k40)` and the same other arguments. *GATE:* `identical()` on `df_super` (every column, `mu1` included). This establishes that the direct route, which the ladder uses, reproduces the calibrated object.
- Enumerate the family at n = 500 with `fs_oc_family_enumerate()`. *GATE:* `identical()` against the corrected run's stored n = 500 family on **all** of `lab`, `Pg`, `PQg`, `beta_g`, `se_g`, `sens_g`, `spec_g`, `ovl`, `M`, `memb`. Expect M = 1696.
- One draw-level check: re-run **the exact call the corrected script made** for the n = 500 alternative point under `"resample"` at the driver's `(c1, c2)` — same function (`fs_oc_predict()` or `fs_oc_grid()`), same `draws`, same `block`, seed 20260825 — on the rebuilt DGM. *GATE:* `identical()` to the corrected run's stored object at that point. Blocked and one-block draws are not bit-identical to each other (the layout differs), so the call must match the stored one's form; say which it was. This is the standing guard that 0.3.1 changed nothing the wrapper touches.

## 3. The ladder — GATE per rung

Rungs are specified by the **oriented true effect in Q**, `abs(m_tau[Q])` from `fs_dgm_scale()`, which for MD40 is 40. The relation is exactly linear (`generate_glm_dgm.R` L462), so for a target `q`:

    k_inter(q) = k40 + s * (q - 40)      with  s = sign(m_tau[Q]) of the MD40 scale table

**Rungs: q ∈ {60, 80, 100, 120, 140, 160}.** For each, build the DGM with `generate_glm_dgm()` and every argument identical to §2's direct rebuild except `k_inter`, then:

- *GATE:* `abs(m_tau[Q])` equals `q` to 1e-9, and `abs(m_tau[Qc])` equals MD40's.
- Enumerate the family at n = 500. *GATE:* `lab`, `Pg`, `PQg`, `se_g`, `sens_g`, `spec_g`, `ovl`, `M`, `memb` are `identical()` to MD40's, and `beta_g − beta_g_MD40` equals `(q − 40) · PQg` to 1e-9. Any other difference means something in the DGM moved with `k_inter` that source said would not — stop.
- `fs_oc_grid()` at n = 500, `"resample"`, the driver's `c2` and `pconsistency`, seed 20260825, the stored `draws` and `block`, over **`c1` from 0 to 300 by 1**. Save per rung the grid's `table` and the full `fs_oc_predict` objects at the named points only — the driver's `c1`, and §4's `c1_05` and `c1_10` once known; **do not commit the per-threshold `results` list** (hundreds of MB across the ladder). If any single `.rds` exceeds 20 MB, leave it out of the commit and say so.

Report per rung, in one table: `q`, `k_inter`, `bint`, and at the driver's own `c1` the full `fs_oc_predict` row (`det_rate ± se`, `EnH`, `Esens`, `Espec`, `Eppv`, `Enpv`, `EbetaH`, `Enaive_bias`, `mass_below`); and the inverted `c1` for declaration 80 / 90 / 95% via `fs_oc_invert(solve_for = "c1")` at the same seed, the way the corrected script did it — its one-block draw set agrees with the blocked grid to Monte-Carlo precision, not bit-for-bit; show the grid's `det_rate` at the nearest integer `c1` beside each inversion as the cross-check.

## 4. The null side, and the `se_g` diagnostic — no new DGM

**Null thresholds.** From the corrected run's stored null inversions at n = 500 (`"resample"`), read the `c1` at which the null false-declaration rate is **5%** — call it `c1_05` — and the stored `ceiling`. If a **10%** inversion (`c1_10`) is not stored, obtain it with one `fs_oc_invert(solve_for = "c1", target = c(0.05, 0.10))` call on the **stored** null family at the stored settings and seed, and *GATE* that call's 0.05 result `identical()` to the stored 0.05 inversion (same function, same seed, same one-block draws) before using its 0.10 value. If the null family is re-enumerated rather than read, gate it `identical()` to the stored one first. **Do not re-derive anything else about the null.**

**Power at the null thresholds.** For every rung of §3, read from its grid `det_rate` at `c1_05` and at `c1_10` (nearest grid point; state the rounding). This is the ladder's central column: **power at fixed type-I**, per harm.

**The `se_g` diagnostic — the §8 condition, measured on the ladder.** `se_g` scales Q's `V_eff` by prevalence. For a mixed-purity candidate the bracket gains `½·k_inter²·PQg(1−PQg)` plus the baseline-mean covariance term (`fs_dgm_scale.R` roxygen, the balanced-arm reading), so the approximation that measured within [0.992, 1.015] at 13.7 is expected to widen with `k_inter`. For **MD40 and every rung**: `fs_dgm_scale(dgm, regions = <each column of family$memb>)`, `se_direct = sqrt(V_eff[g] / (n·Pg))` at n = 500, `ratio = se_g / se_direct`. Report per rung: the ratio's range, median, and its distribution in three purity bands (`PQg ≥ 0.95`; `0.25 ≤ PQg < 0.95`; `PQg < 0.25`), plus the correlation with `Pg`. State plainly, per rung, whether the band has moved outside the 2% measured at MD40. **This is reported, not acted on; §6 is the only place it is used, and only as a sensitivity.**

## 5. The forecast rung — GATE

From §4's power-at-`c1_05` column, find the smallest `q` on the ladder with power ≥ 0.80. Then:

- If two rungs bracket 0.80, **interpolate `q` linearly** between them, **round up to the nearest multiple of 5**, and build that DGM as a seventh rung under §3's gates — the forecast is stated at a built rung, never at an interpolation.
- If no rung reaches 0.80 at `c1_05`, **extend the ladder upward by 20 at a time** (180, 200, 220) under §3's gates until one does, to a maximum of `q = 220`. If 220 does not reach it, the forecast rung is the smallest `q` reaching 0.80 at **`c1_10`** instead, found the same way — and say so at the top of the report.
- If the lowest rung (60) already exceeds 0.80 at `c1_05`, extend **downward** by 10 (50) and take the smallest rung at or above 0.80.

At the forecast rung run **both gates** — `"resample"` and `"split"` — on the same seed and settings, and report:

1. `c1*`: the largest `c1` with predicted power ≥ 0.80 under `"resample"` (`fs_oc_invert`, target 0.80 — the exact order statistic), with the achieved rate and MC SE, and the grid's `det_rate` at the nearest integer `c1` as the cross-check; and the null false-declaration rate at that same `c1*`, obtained as `fs_oc_predict(family = <the stored null family>, n = 500, c1 = c1*, c2 = <driver c2>, "resample", seed 20260825, stored draws)` — the exact rate at that threshold, with its MC SE.
2. The complete `fs_oc_predict` row at `(n = 500, c1*, driver c2)` under both gates: `det_rate`, `EnH`, `Esens`, `Espec`, `Eppv`, `Enpv`, `EbetaH`, `Enaive_bias`, `mass_below`, each with its MC SE where it has one.
3. The selection distribution given declaration: `sel_c` mass on Q itself, on candidates with `PQg ≥ 0.95`, and on the top 15 candidates by `sel_c` with their `Pg`, `PQg` and `beta_g`.
4. The same row at the **driver's own `c1`**, so the saturated screen's picture at this harm is on the record beside the crossover.

*GATE:* §3's identity gates at the forecast rung; `c1*` attainable (not `NA`).

## 6. Sensitivity — the direct-`V_eff` family, at the forecast rung only

Take the forecast rung's family object and replace `se_g` with §4's `se_direct` for every candidate (a modified `fs_oc_family` passed via `family =`; nothing in `R/` is touched). Re-run `fs_oc_predict()` at `(n = 500, c1*, driver c2)`, `"resample"`, same seed, and report beside §5's row: power at `c1*`, `EnH`, `Eppv`, `Esens`, `EbetaH`, and the `sel_c` mass on Q.

**State the status in these words:** this is a sensitivity showing what the prevalence-scaling approximation costs at this harm. It is **not** a proposed constructor, it is **not** adopted, and the forecast of §5 stands as issued with the wrapper as it is. Whether `se_g` should ever come from the direct per-candidate `V_eff` is Larry's decision, taken on these numbers.

## 7. Report — the forecast is the deliverable

`dev/glm-continuous-sims/REPORT_oc_breadth_ladder_2026-08-30.md`:

0. Provenance (raw), and §2's build quotations.
1. **The forecast, first, in one table** — the rung (`q`, `k_inter`), `c1*`, predicted power ± SE, predicted null type-I at `c1*` ± SE, and every §5 quantity under both gates. Then, in prose, the three predictions Stage 2 will be scored against: (a) the declaration rate at `c1*`; (b) the null false-declaration rate at `c1*`, which Stage 2 checks against the **existing** null payload without a new null run; (c) the selected rule's `E[β(Ĥ)]`, specificity, PPV, sensitivity and `E|Ĥ|` — **carrying the documented limitation**: `E|Ĥ|`, PPV, sensitivity and the naive bias inherit a between-rule gap of +2.11 / +0.61 / +1.65 subjects at MD40 that five mechanisms failed to explain, and the population-versus-sample offset on PPV and sensitivity is expected; whether that gap moves with the harm is a finding of Stage 2, not a prediction of Stage 1.
2. §3's ladder table and §4's power-at-fixed-type-I column — this is the figure the paper wants: power against `c1` per harm, type-I from the null.
3. §4's `se_g` diagnostic per rung, with the plain statement of where the band left 2%.
4. §6's sensitivity beside §5, with the status paragraph verbatim.
5. Every gate's arithmetic shown.
6. Ten-line summary.

Commit scripts, `.rds` outputs, logs and the report by explicit path; **the commit hash of the report is the lock — quote it in the ten-line summary.** No push. No `devtools::install()`.

## 8. Out of scope

- No `R/` change, no new constructor, no edit to any driver, application or document. `oc_wrapper_verification.qmd` is not touched.
- **No Stage 2** — no replicate draws, no `forestsearch()` runs, no driver execution. The measured cell is a separate go/no-go after Larry reads the forecast.
- No `n = 700` rungs, no second DGM population, no noise covariates, no binary/OR work.
- No re-derivation of the corrected MD40 run beyond §2's gates; no residual mechanism of any kind — the residual is closed as a documented limitation.
- No push, no install, no renders.
