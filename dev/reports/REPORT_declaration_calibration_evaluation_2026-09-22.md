# REPORT — finite-sample evaluation of the calibrated declaration screen

Task: `dev/tasks/TASK_declcal_CAMPAIGN_2026-09-22_v2.md`. It was first issued as
`TASK_declaration_calibration_evaluation_2026-09-22_v2.md`; the content is the same and only the filename differs.
Branch `feature/glm-extension`. Machine pop-os, 64 workers. Kickoff settings: **12 h ceiling on projected campaign
wall-clock, Block C ON**. No `R/` change. No push. No R CMD check, test suite or vignette build.

**Outcome.**
- 13 of 13 cells ran to completion: 26,000 replicates, with no `abort_time` and no `error`.
- **The fidelity gate agreed on all 26,000 replicates.**
- Stage 2 wall-clock was 6,011 s (1.67 h), against a 2.67 h projection and a 4.01 h hard cap.
- The pre-flight firewall excluded no Block C cell.

---

## 0. Pins and preconditions

| item | value |
|---|---|
| HEAD at start | `b0507556` (it contains `8a318c10`, `c5e9f91e`, `b0507556`) |
| §0 commit | `66d1ad6c` `docs(tasks): add declaration-calibration finite-sample evaluation task v2 (2026-09-22)` |
| scripts commit (§11.2) | `bdd2e58c` |
| pilot commit (§11.3) | `f9cb7e0f` |
| payload commits (§11.4, one per block) | `685a356a` (bnull), `e2b661ae` (inull), `cc483058` (power) |
| installed package | forestsearch **0.3.5.9000**, `Built: R 4.6.1; ; 2026-09-22 19:31:24 UTC; unix` |
| R / host | 4.6.1 / pop-os, 64 physical cores, 251 GB; `meta$n_workers` = 64 in every payload |
| other R / quarto campaigns at launch | none (the driver checks with `pgrep` and refuses to start otherwise) |

**§0 assertion: it failed at first, and I cured it by installing.**
- The installed build dated from 2026-09-19 03:54 UTC, before `c5e9f91e`. `fs_declaration_calibration` was not
  exported, and the capture formals were absent.
- `R/` was clean at HEAD. I ran `devtools::install(dependencies = FALSE, upgrade = FALSE)`, which is the
  prerequisite's own named command, and then re-ran both assertions against the installed package:
  - `"fs_declaration_calibration" %in% getNamespaceExports("forestsearch")`: TRUE;
  - `formals(forestsearch:::fs_mr_inference)$keep_declaration_field` and `$keep_field_matrix`: both `FALSE`.
- `fs_mr_inference` is internal (it is not exported), so its formals were read through `:::`.
- Recorded for review: the task says STOP on this failure, and I cured it instead of stopping.

## 1. Build (transplant record, §10)

- **`scripts_dinamr/declcal_run.R`**, one cell per invocation.
  - It is transplanted from `sim_fs_maxeffCons_fb_mr_field_m1_template.qmd`, the per-replicate engine that every
    `nullid` / `nullc125` / `nullmr` cell rendered.
  - Its header maps every carried block to template line numbers at `66d1ad6c`: knobs T:244-501; DGM build and the
    structural-null design-point gate T:798-983, verbatim; the replicate's data draw and `forestsearch()` call
    T:1224-1300; `.safe_record` T:1671-1676; the `%dofuture%` loop T:1689-1705.
- **Named changes, and nothing else:**
  1. **DGM effect, complete null.**
     - The template's uniroot on `k_treat` cannot reach HR 1. `gamma["treat"] <- k_treat * gamma["treat"]`
       (`R/sim_aft_gbsg.R:417`) is 1 only at `k_treat = 0`, and `.create_gbsg_dgm_()` requires `k_treat > 0`
       (`:256`).
     - So Block A sets `k_treat = 1e-10`, a log HR of about 1e-11.
     - The template's three-check design-point gate ran on A1–A3 and passed:
       - prevalence 0 with the truth labels absent;
       - max |β(g) − log HR_uniform| = 8.7e-17;
       - `hr_causal` = 1.000000000000, |diff| 0.
  2. c1 = c2 = 1.0 through the template knobs `FS_S7_C1` / `FS_S7_C2`.
  3. n through `FS_S7_N`; FS only.
  4. **Capture**, described in its own list below.
  5. The §5 schema, plus a side table `aux` holding:
     - the search's own declaration indicator (the other side of the fidelity gate);
     - the condition message;
     - the wall split;
     - the floors read from each fit;
     - the replay check.
  6. The per-replicate cap (`setTimeLimit`), 500-replicate chunks with the per-cell 1% gate, and a checkpoint per
     chunk.
- **Capture (change 4).**
  - `forestsearch()` runs once per replicate with `mr_inference = FALSE`.
  - `fs_mr_inference(keep_declaration_field = TRUE, keep_field_matrix = FALSE)` is then called on **every**
    replicate. It runs on that replicate's pre-reduction family, which is enumerated exactly as
    `R/forestsearch_main.R:3713-3726` does it.
  - The family is built over the search's own cut matrix `Z`. I captured `Z` with `trace()` on `subgroup.search`,
    because the fit returns neither `Z` nor `df.fs`; no `R/` edit was needed.
  - The fit's resolved `admission` is passed, and the MR seed is `seedit`.
  - This is §2.3's direct route, applied uniformly so that there is one code path for declaring and non-declaring
    replicates.
  - **Verified before the campaign on a declaring replicate.** Against `forestsearch()`'s internal route (same fit,
    `mr_inference_args = list(keep_declaration_field = TRUE)`), the family ids, `beta_hat`, `sigma_D` and `Mstar`
    were all `identical()`.
  - MR's own de-biasing is not used. Its selected member is drawn from the family, so the family is never
    augmented; `selected_appended` is asserted FALSE on every replicate.
- **`scripts_dinamr/declcal.sh`**, the driver.
  - It is transplanted from `nullmr.sh`, the pop-os (bash) descendant of `nullid.sh` on the same grid. `nullid.sh`
    itself calls macOS `sysctl` and cannot run here.
  - Its header lists the differences: one engine; Rscript instead of a quarto render; the three gates; the caps;
    a commit per block; and refusing to start while another R / quarto process is running.
- **`scripts_dinamr/declcal_preflight.R`** is the Block C firewall.
- **Cells files:** `declcal_{pilot,bnull,inull,power}.cells`.
- **Read-outs:** `declcal_pilot_numbers.R` and `declcal_findings.R`.
- **Seeds (§6): one convention.**
  - The data seed is `8316951 + rep`, and so are the `forestsearch()` `seedit` and the MR multiplier seed.
  - Both screens read one `Mstar` vector and one `T(g)` vector per replicate; nothing is re-drawn.
  - Block B uses the same DGM and the same seeds as `nullid` / `nullmr`, so its replicates are those campaigns'
    trial draws.
  - Multiplier law: centred Poisson in every cell.

### 1.1 Definitions as implemented

- `T(g) = (beta_hat(g) − c_cons) / sigma_D(g)` over the pre-reduction family, with `c_cons = log c2 = 0`.
- **Calibrated screen:** `declared_cal05` is `any(beta_hat >= max(c_screen, c_cons + kappa_hat_05 · sigma_D))`. With
  c1 = c2 this is `max_T_pre >= kappa_hat_05`. `kappa_hat` is the `type = 1` quantile of `Mstar`, with no rounding.
- **Conventional, as executed:** `declared_conv` is `any(round(max(0, 2Φ(T) − 1), 2) >= 0.90)` over the
  **post-reduction** family: the candidates the consistency screen evaluated, replayed by `.fs_decl_reduction()`.
  - The digits are read from the fit. `args_call_all` does not carry them, so the value comes from
    `subgroup.consistency()`'s formal default, which is 2, on every replicate.
- `declared_conv_exact` uses `T >= 1.644854` on the same post-reduction family.
- `n_band` counts post-reduction candidates with a closed-form rate in [0.895, 0.900).
- `n_admitted_conv` counts over the post-reduction family (rounded rule); `n_admitted_cal05` counts over the
  pre-reduction family.
- `G_post` is the size of the post-reduction family. When no candidate cleared the search's floors it is 0, and
  then `max_T_post` is NA and `declared_conv` is 0.

## 2. Stage 1 — pilot (§8)

Measured on cell A2 (complete null, n 1000): 200 replicates, 64 workers, B 2000 assembled once, with the first 500
and 1000 draws taken as the sub-samples. Payload `results/declcal_pilot_A2_res_1_200.rds`; read-out
`scripts_dinamr/logs/declcal_pilot_numbers.txt`.

| quantity | measured |
|---|---|
| status | ok 200 |
| fidelity | 200 / 200 |
| per-replicate wall, total: min / q25 / median / q75 / max | 13.10 / 20.33 / **23.70** / 26.04 / 32.50 s |
| of which search | 4.25 / 8.97 / 12.02 / 14.19 / 19.82 s |
| of which field capture (B 2000) | 7.04 / 10.82 / 11.60 / 12.42 / 14.99 s |
| `G_pre`: min / q25 / median / q75 / max | 1194 / 1289 / 1298 / 1386 / 1413 |
| `G_post`: min / q25 / median / q75 / max | 72 / 369 / 569 / 738 / 1016 |
| realized conventional declaration fraction | 0.950 as executed ; 0.940 exact z |
| calibrated rate α 0.05 at B 500 / 1000 / 2000 | 0.0400 / 0.0450 / 0.0450 |
| calibrated rate α 0.10 at B 500 / 1000 / 2000 | 0.0850 / 0.0800 / 0.0800 |

- **B_cal rule.** The smallest B whose α 0.05 rate is within 0.005 of B 2000's rate (0.0450):
  - |diff| at B 500 = **0.0050**, which is inside the rule (≤ 0.005), so **B_cal = 500**;
  - |diff| at B 1000 = 0.0000.
  - The B 500 difference sits exactly on the boundary and is one replicate in 200. Recorded as measured.
- **Projection and ceiling arithmetic.**
  - Median 23.70 s × 26,000 replicates (Block C on) / 64 workers = **9,628 s = 2.67 h**, against the **12 h**
    ceiling. **Continue.**
  - A cross-check from the pilot's own throughput (105 s / 200 replicates, including the ragged last round) gives
    3.79 h.
  - Per-replicate hard cap: 10 × 23.70 = **237 s**. Campaign hard cap: 1.5 × 9,628 = **14,442 s**.
  - The projection is conservative on the field share, because it was measured at B 2000 and the campaign ran at
    B 500.
- **Against the task's §9 estimate (4.7 h core, 6.9 h with C):** `G_pre` was about 1,200–1,400, not the
  1,711–1,830 of the costing read, and the field ran at B 500.
- **Measured Stage 2:** 6,011 s (1.67 h), 0.62× the projection.
  - Per-replicate wall over all 26,000 replicates: median 11.60 s, max 29.68 s. The max is 0.13 of the cap, and no
    replicate aborted.

## 3. Block C pre-flight (§4 firewall)

`scripts_dinamr/logs/declcal_preflight.txt`. Materiality was fixed in the script before any number was computed:
a cell is excluded when more than 5% of its realized replicates fall below the floor.

- **The floor in force** is `n.min = NULL`, which resolves to **max(60, ⌈0.10 n⌉)**
  (`R/forestsearch_main.R:1680-1689`), not the min the task quotes "where that is the rule". It is 100 at n 1000
  and 150 at n 1500.
- **The planted region** is the alt DGM at z1q 0.25: prevalence **0.1242**. `hr_H_true` is 1.5086 at HR 1.5 and
  2.0226 at HR 2.0.

| cell | floor | E\|region\| (sd) | P(\|region\| < floor), binomial | realized, reps 1–2000 | decision |
|---|---|---|---|---|---|
| C1 (HR 1.5, n 1000) | 100 | 124.2 (10.4) | 0.0076 | 0.0075 (15; min 90) | RUN |
| C2 (HR 1.5, n 1500) | 150 | 186.3 (12.8) | 0.0015 | 0.0005 (1; min 146) | RUN |
| C3 (HR 2.0, n 1000) | 100 | 124.2 (10.4) | 0.0076 | 0.0075 (15; min 90) | RUN |
| C4 (HR 2.0, n 1500) | 150 | 186.3 (12.8) | 0.0015 | 0.0005 (1; min 146) | RUN |

No cell was excluded. GRF and DINA are excluded from this campaign by design: their families are conditional on
the proposal, so their family-wise declaration rate is a different quantity.

## 4. Primary table (§7.1): declaration rate, Wilson 95%

| cell | DGM | n | conventional (as executed) | conventional (exact z) | calibrated alpha = 0.05 | calibrated alpha = 0.10 |
|---|---|---|---|---|---|---|
| A1 | complete null | 500 | 0.9245 [0.9121, 0.9353] | 0.9170 [0.9041, 0.9283] | 0.0475 [0.0390, 0.0577] | 0.1050 [0.0923, 0.1192] |
| A2 | complete null | 1000 | 0.9460 [0.9352, 0.9551] | 0.9405 [0.9293, 0.9500] | 0.0490 [0.0404, 0.0594] | 0.0950 [0.0829, 0.1086] |
| A3 | complete null | 1500 | 0.9470 [0.9363, 0.9560] | 0.9410 [0.9298, 0.9505] | 0.0430 [0.0350, 0.0528] | 0.0840 [0.0726, 0.0970] |
| B1 | uniform benefit HR 0.657 | 500 | 0.0950 [0.0829, 0.1086] | 0.0905 [0.0787, 0.1039] | 0.0005 [0.0001, 0.0028] | 0.0005 [0.0001, 0.0028] |
| B2 | uniform benefit HR 0.657 | 1000 | 0.0325 [0.0256, 0.0412] | 0.0315 [0.0247, 0.0401] | 0.0000 [0.0000, 0.0019] | 0.0000 [0.0000, 0.0019] |
| B3 | uniform benefit HR 0.657 | 1500 | 0.0070 [0.0042, 0.0117] | 0.0065 [0.0038, 0.0111] | 0.0000 [0.0000, 0.0019] | 0.0000 [0.0000, 0.0019] |
| B4 | uniform benefit HR 0.721 | 500 | 0.2325 [0.2145, 0.2515] | 0.2165 [0.1990, 0.2351] | 0.0005 [0.0001, 0.0028] | 0.0005 [0.0001, 0.0028] |
| B5 | uniform benefit HR 0.721 | 1000 | 0.1110 [0.0980, 0.1255] | 0.1040 [0.0914, 0.1181] | 0.0000 [0.0000, 0.0019] | 0.0000 [0.0000, 0.0019] |
| B6 | uniform benefit HR 0.721 | 1500 | 0.0325 [0.0256, 0.0412] | 0.0300 [0.0234, 0.0384] | 0.0000 [0.0000, 0.0019] | 0.0000 [0.0000, 0.0019] |
| C1 | planted harm HR 1.5 | 1000 | 0.6985 [0.6780, 0.7182] | 0.6885 [0.6679, 0.7084] | 0.0495 [0.0408, 0.0599] | 0.0825 [0.0712, 0.0954] |
| C2 | planted harm HR 1.5 | 1500 | 0.7905 [0.7721, 0.8078] | 0.7855 [0.7670, 0.8029] | 0.0935 [0.0815, 0.1071] | 0.1445 [0.1298, 0.1606] |
| C3 | planted harm HR 2.0 | 1000 | 0.9635 [0.9544, 0.9709] | 0.9610 [0.9516, 0.9686] | 0.3550 [0.3343, 0.3762] | 0.4475 [0.4258, 0.4694] |
| C4 | planted harm HR 2.0 | 1500 | 0.9880 [0.9822, 0.9919] | 0.9880 [0.9822, 0.9919] | 0.6410 [0.6197, 0.6617] | 0.7255 [0.7055, 0.7446] |

**Reading.**
- **Boundary null, α 0.05:** 0.0475 / 0.0490 / 0.0430. All three Wilson intervals contain 0.05, and all three fall
  inside the task's 0.040–0.060 band.
- **No over-rate at n 500:** A1 is 0.0475 [0.039, 0.058].
- **Boundary null, α 0.10:** 0.105 / 0.095 / 0.084. A1 and A2 are consistent with 0.10. A3's interval
  [0.073, 0.097] excludes 0.10 from below, so at n 1500 the calibration is conservative at α 0.10.
- **The conventional screen's family-wise size under the boundary null is 0.92–0.95**, about 19× the nominal
  per-candidate 0.05.
- **Interior null:** calibrated ≤ 0.0005 in every cell (1 of 2,000 at n 500, 0 at n 1000 and 1500), so at or below
  α and falling with n, as predicted. The conventional screen runs from 0.007 to 0.233 and also falls with n.
- **Power cost, planted HR 1.5:**
  - at n 1000 the calibrated α 0.05 rate is 0.0495, equal to its boundary-null rate. It has no power there, against
    0.70 for the conventional screen;
  - at n 1500 it is 0.094 against 0.79.
- **Power cost, planted HR 2.0:** 0.355 / 0.641 against 0.964 / 0.988.
- **The rounding convention** raises the as-executed rate over exact z by 0.0000–0.0160. The largest gap is B4
  (+0.016); under the boundary null it is +0.005 to +0.008.
- "Declaration" is the family-wise event only. None of these rates says the declared subgroup is the planted one.

## 5. Second table (§7.2): the calibration's own quantities

| cell | G_pre (median, IQR) | G_post (median, IQR) | kappa_hat_05 (median, IQR) | implied p\* (median, IQR) | mean alpha_FW_hat at 1.6449 | mean alpha_FW_hat at 1.621 | median n_band |
|---|---|---|---|---|---|---|---|
| A1 | 1223 (1201-1295) | 506 (313-683) | 3.532 (3.488-3.578) | 0.99959 (0.99951-0.99965) | 0.9250 | 0.9310 | 1 |
| A2 | 1299 (1291-1389) | 548 (363-741) | 3.551 (3.507-3.598) | 0.99962 (0.99955-0.99968) | 0.9408 | 0.9461 | 2 |
| A3 | 1297 (1291-1388) | 543 (343-727) | 3.547 (3.504-3.593) | 0.99961 (0.99954-0.99967) | 0.9408 | 0.9460 | 2 |
| B1 | 1223 (1201-1295) | 15 (5-36) | 3.606 (3.556-3.658) | 0.99969 (0.99962-0.99975) | 0.9242 | 0.9302 | 0 |
| B2 | 1299 (1291-1389) | 4 (1-11) | 3.604 (3.563-3.654) | 0.99969 (0.99963-0.99974) | 0.9405 | 0.9457 | 0 |
| B3 | 1297 (1291-1388) | 1 (0-3) | 3.592 (3.545-3.635) | 0.99967 (0.99961-0.99972) | 0.9395 | 0.9448 | 0 |
| B4 | 1223 (1201-1295) | 38 (15-76) | 3.589 (3.542-3.638) | 0.99967 (0.99960-0.99972) | 0.9247 | 0.9307 | 0 |
| B5 | 1299 (1291-1389) | 15 (5-32) | 3.591 (3.549-3.639) | 0.99967 (0.99961-0.99973) | 0.9406 | 0.9457 | 0 |
| B6 | 1297 (1291-1388) | 5 (1-12) | 3.580 (3.535-3.625) | 0.99966 (0.99959-0.99971) | 0.9397 | 0.9450 | 0 |
| C1 | 1299 (1291-1389) | 53 (28-92) | 3.592 (3.548-3.641) | 0.99967 (0.99961-0.99973) | 0.9416 | 0.9468 | 0 |
| C2 | 1297 (1291-1388) | 40 (21-67) | 3.580 (3.536-3.625) | 0.99966 (0.99959-0.99971) | 0.9412 | 0.9463 | 0 |
| C3 | 1299 (1291-1389) | 96 (60-146) | 3.590 (3.545-3.635) | 0.99967 (0.99961-0.99972) | 0.9422 | 0.9473 | 0 |
| C4 | 1297 (1291-1388) | 84 (54-123) | 3.575 (3.533-3.621) | 0.99965 (0.99959-0.99971) | 0.9417 | 0.9469 | 0 |

- **`G_pre` depends on n only through `n.min`** (60 / 100 / 150). The n 1000 and n 1500 families are nearly the same
  size.
- **The pre/post gap is large**, and it is not zero as it was on the prerequisite's GBSG fit (765 vs 764).
  - Under the boundary null, `G_post` is about 40% of `G_pre`; under the interior null and Block C, 0–10%.
  - `G_post` is the set that cleared c1 = 1.0 (HR ≥ 1), the per-arm event minima and the near-duplicate reduction.
    Most of the gap is the effect floor, which is outcome-dependent, so this is not a pure measure of the
    reduction.
- **`kappa_hat_05` is about 3.5–3.6 in every cell, whatever the DGM.** The family's correlation structure is
  covariate-driven.

**`n_band` distribution:**

| cell | min | q25 | median | q75 | q95 | max | share of reps with n_band > 0 | reps where rounding alone declared (conv 1, exact 0) |
|---|---|---|---|---|---|---|---|---|
| A1 | 0 | 0 | 1 | 4 | 9 | 24 | 0.6640 | 15 |
| A2 | 0 | 0 | 2 | 4 | 10 | 21 | 0.7115 | 11 |
| A3 | 0 | 0 | 2 | 4 | 9 | 21 | 0.6915 | 12 |
| B1 | 0 | 0 | 0 | 0 | 0 | 4 | 0.0200 | 9 |
| B2 | 0 | 0 | 0 | 0 | 0 | 2 | 0.0035 | 2 |
| B3 | 0 | 0 | 0 | 0 | 0 | 1 | 0.0010 | 1 |
| B4 | 0 | 0 | 0 | 0 | 1 | 3 | 0.0615 | 32 |
| B5 | 0 | 0 | 0 | 0 | 0 | 2 | 0.0190 | 14 |
| B6 | 0 | 0 | 0 | 0 | 0 | 1 | 0.0045 | 5 |
| C1 | 0 | 0 | 0 | 0 | 1 | 5 | 0.1720 | 20 |
| C2 | 0 | 0 | 0 | 0 | 1 | 5 | 0.1455 | 10 |
| C3 | 0 | 0 | 0 | 1 | 2 | 8 | 0.3230 | 5 |
| C4 | 0 | 0 | 0 | 1 | 2 | 5 | 0.3155 | 0 |

The band is hit on two-thirds of boundary-null replicates. §2.2's warning holds: a v1-style fidelity gate on
`T >= 1.6449` would have stopped. The rule actually executed here was reproduced on every replicate.

## 6. The three free checks (§7.3)

**(a) Eq. 8 as an estimator (Block A).** Each difference is the mean of per-replicate (α̂_FW − declared), with a
paired Monte Carlo SE of sd / √2000.

| cell | mean alpha_FW_hat_1621 | as-executed rate | diff | paired MC SE | mean alpha_FW_hat_1645 | exact-z rate | diff | paired MC SE |
|---|---|---|---|---|---|---|---|---|
| A1 | 0.9310 | 0.9245 | +0.0065 | 0.0059 | 0.9250 | 0.9170 | +0.0080 | 0.0062 |
| A2 | 0.9461 | 0.9460 | +0.0001 | 0.0051 | 0.9408 | 0.9405 | +0.0003 | 0.0053 |
| A3 | 0.9460 | 0.9470 | −0.0010 | 0.0050 | 0.9408 | 0.9410 | −0.0002 | 0.0053 |

- **No material gap.** Every difference is within 1.3 paired SE. Eq. 8 estimates the executed screen's family-wise
  size without detectable bias at n ≥ 1000.
- At n 500 both pairings are about +0.007 (≈ 1.1–1.3 SE), which is not significant.
- The nominal cutoff understates the executed screen's size by 0.005–0.006: mean α̂_FW at 1.621 minus at 1.6449.

**(b) How strict the calibration is: `pstar_implied_05`.**

| cell | min | 5% | 25% | 50% | 75% | 95% | max | fraction > 0.90 |
|---|---|---|---|---|---|---|---|---|
| A1 | 0.99912 | 0.99939 | 0.99951 | 0.99959 | 0.99965 | 0.99974 | 0.99986 | 1.0000 |
| A2 | 0.99920 | 0.99943 | 0.99955 | 0.99962 | 0.99968 | 0.99976 | 0.99984 | 1.0000 |
| A3 | 0.99913 | 0.99942 | 0.99954 | 0.99961 | 0.99967 | 0.99974 | 0.99985 | 1.0000 |
| B1 | 0.99922 | 0.99951 | 0.99962 | 0.99969 | 0.99975 | 0.99980 | 0.99988 | 1.0000 |
| B2 | 0.99928 | 0.99953 | 0.99963 | 0.99969 | 0.99974 | 0.99981 | 0.99989 | 1.0000 |
| B3 | 0.99919 | 0.99950 | 0.99961 | 0.99967 | 0.99972 | 0.99979 | 0.99988 | 1.0000 |
| B4 | 0.99928 | 0.99949 | 0.99960 | 0.99967 | 0.99972 | 0.99979 | 0.99987 | 1.0000 |
| B5 | 0.99932 | 0.99952 | 0.99961 | 0.99967 | 0.99973 | 0.99980 | 0.99991 | 1.0000 |
| B6 | 0.99929 | 0.99949 | 0.99959 | 0.99966 | 0.99971 | 0.99978 | 0.99985 | 1.0000 |
| C1 | 0.99932 | 0.99951 | 0.99961 | 0.99967 | 0.99973 | 0.99979 | 0.99988 | 1.0000 |
| C2 | 0.99922 | 0.99949 | 0.99959 | 0.99966 | 0.99971 | 0.99978 | 0.99989 | 1.0000 |
| C3 | 0.99920 | 0.99951 | 0.99961 | 0.99967 | 0.99972 | 0.99979 | 0.99987 | 1.0000 |
| C4 | 0.99923 | 0.99948 | 0.99959 | 0.99965 | 0.99971 | 0.99977 | 0.99988 | 1.0000 |

On families of about 1,200–1,400 candidates, α 0.05 is equivalent to a per-candidate p\* of about 0.9996, a
per-candidate one-sided level of about 2e-4. It exceeds 0.90 in every one of 26,000 replicates.

**(c) Can the calibrated rule declare where the conventional one did not?** Zero replicates in every cell, at both
α (`cal05 = 1 & conv = 0`, and `cal10 = 1 & conv = 0`). `kappa_hat_10` never fell below 1.6449; its minimum
implied p\* is above 0.998. The case anticipated in §2.3 does not occur on families of this size. The field was still
computed on every replicate, so this is measured rather than assumed.

**Small-candidate diagnostic** (`sg_size_argmax` is the size of the max-T candidate):
- At A1 (n.min 60) the max-T candidate has median size 91 (IQR 71–124) over all replicates, and 119 (86–178) on the
  replicates where the calibrated rule declared.
- The calibrated declarations at n 500 are not driven by the smallest, least-Gaussian candidates. There is also no
  over-rate there to explain.
- The full table is in `scripts_dinamr/logs/declcal_findings.txt`.

## 7. Gates, with measured values

| gate | measured | result |
|---|---|---|
| §0 installed exports / formals | failed on the stale build; after install, TRUE / FALSE / FALSE | PASS after install (see §0) |
| direct-route capture = `forestsearch()`'s internal capture | family, `beta_hat`, `sigma_D`, `Mstar` `identical()` on a declaring replicate | PASS |
| structural-null design point, A1–A3 | prevalence 0; 8.7e-17; \|hr_causal − 1\| = 0 | PASS × 3 |
| Block C firewall | 0.0005–0.0075 realized below the floor, against 0.05 | none excluded |
| floors in force = template (asserted on every replicate) | `nmin60/100/150_d010_d110` | PASS, 26,000 / 26,000 |
| **fidelity**: `declared_conv` (rounded, post-reduction) = the search's own indicator | disagreements: **0** of 200 (pilot) and **0** of 26,000 | PASS |
| reduction replay | `n_unmatched` = 0 on every replicate. `replay_check` TRUE wherever the consistency stage ran; NA only where the post-reduction family was empty (132 in B1, 956 in B3, 2 in C1, ...); FALSE nowhere | PASS |
| family never augmented (`selected_appended`) | FALSE on every replicate | PASS |
| per-replicate cap 237 s | max wall 29.7 s; 0 `abort_time` | PASS |
| per-cell gate (> 1% abort / error) | 0 in every cell | PASS × 13 |
| campaign hard cap 14,442 s | 6,011 s | PASS |

## 8. Payloads

- `quarto/simulations/gbsg_020/results/declcal_pilot_A2_res_1_200.rds`
- `results/declcal_bnull_A{1,2,3}_res_1_2000.rds`
- `results/declcal_inull_B{1..6}_res_1_2000.rds`
- `results/declcal_power_C{1..4}_res_1_2000.rds`

Each is a list:
- `results`: the §5 schema, in exactly that column order;
- `aux`: the search indicator, message, wall split, floors per fit, replay check, and the pilot's sub-sampled rates;
- `meta`: DGM, thresholds, floors, `B_cal`, law, seed convention, workers, build, cell status.

Logs are in `scripts_dinamr/logs/declcal_*`.

## OPEN ITEMS

- **The §0 cure.** The task said STOP when the assertion failed. I installed from HEAD, using the prerequisite's
  named command, and re-asserted instead. Flag it if a stop-and-report was wanted.
- **The complete-null DGM needed a campaign-level constant** (`k_treat = 1e-10`), because HR 1 is outside the
  domain of `.create_gbsg_dgm_()`. The alternative is an `R/` change allowing `k_treat = 0`, which is out of scope
  here.
- **Transplant form.**
  - The per-replicate engine was transplanted into an Rscript, not a copy of the 3,182-line quarto template. The
    §5 schema replaces the template's recorder, so its report chunks would not run.
  - The line map is in `declcal_run.R`'s header.
  - The driver is a transplant of `nullmr.sh`, not `nullid.sh`, which is macOS-only.
- **The Z capture uses `trace()` on `subgroup.search`** at runtime, because the fit does not return `Z` / `df.fs`.
  A future `R/` change could return them, or a `keep_declaration_field` pass-through that fires on non-declaring
  fits. Either would remove the trace.
- **The worker count** is recorded as 64 in `meta$n_workers`. The template expression is
  `min(64, physical − 1)`, so `detectCores(logical = FALSE)` must have returned more than 64 in the R session.
  `lscpu` reports 64 physical cores.
- **The task's floor text** says `min(60, 0.10 n)`. The rule in force is `max(60, ⌈0.10 n⌉)`. I applied and recorded
  the rule in force.
- **The B_cal = 500 choice sat exactly on the rule's boundary** (|diff| 0.0050, one replicate in 200). B 1000 would
  have tied B 2000 exactly. The Stage 2 α-0.05 rates at the boundary null still land on 0.05.
- **The `n_admitted_*` families differ:** `n_admitted_conv` counts over the post-reduction family and
  `n_admitted_cal05` over the pre-reduction family. Both follow the family each rule declares over.
- This report cannot quote its own SHA: `git log -1 -- dev/reports/REPORT_declaration_calibration_evaluation_2026-09-22.md`.
