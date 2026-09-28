# REPORT — `gbsg_020` DGM: the trial-wide (ITT) hazard ratio per cell (read-only), 2026-09-28

**Task:** `dev/tasks/TASK_gbsg020_itt_estimand_2026-09-28.md` (committed `3d796ab1`, as received; supersedes `TASK_gbsg020_dgm_estimands_2026-09-24.md`, never run).
**Repo / branch:** `forestsearch`, `feature/glm-extension`. HEAD at the start of the reading: `3d796ab1` (parent `bbec91e4`).
**Kind:** read-only. Nothing re-run, nothing computed. R was used once, in a session that wrote nothing, to read `names()`, `class()` and stored constants (`truth`, `meta`) from committed bundles. md5 checksums of every bundle read were taken before and after and are identical (26 files).

## Gates

- G1 (no tracked file modified before any step): PASS. `git status --porcelain` at start showed only untracked entries, listed here and left alone:
  `quarto/simulations/actg175/binary_020/mr_or_harm/fs_effMaxSG_mr_field_or075_n500_nb20_redes_d5000/`,
  `quarto/simulations/actg175/binary_020/mr_or_harm/fs_effMaxSG_mr_field_or075_n500_nb20_relaunch_d5000/`,
  `quarto/simulations/actg175/binary_020/smoke_redes.html`, `quarto/simulations/actg175/binary_020/smoke_relaunch.html`,
  `quarto/simulations/gbsg_020/scripts_dinamr/logs/nullmr_findings.err`.
- G2 (repo and branch): PASS.
- G3 (`gbsg_020` directory and at least one committed results bundle): PASS. 948 tracked files under `quarto/simulations/gbsg_020/results/`.

## Files read

- `quarto/simulations/gbsg_020/current_status.md` (§1, §2.1, §2.8).
- `quarto/simulations/gbsg_020/sim_fs_maxeffCons_fb_mr_field_m1_template.qmd` (lines 357–397, 796–870, 925–962, 1855–1910).
- `R/sim_aft_gbsg.R` (`.create_gbsg_dgm_()` lines 233–250, 280–330, 405–425, 470–565, 666–724; `calibrate_k_inter()` lines 1004–1043).
- `R/setup_gbsg_dgm.R` (lines 85–135).
- `R/oc_analyses.R` (`compute_dgm_cde()` lines 131–176).
- Committed bundles, all under `quarto/simulations/gbsg_020/results/` (the designated FS comparator per §2.1 at every grid cell, all three n; the FS `nullmr` bundles at the two null design points, all three n; and two FS `nullid` bundles as a cross-check):
  - 12.4 %: `fs_maxeffCons_fb_mr_field_m1_h150_knoise0_n{500,1000,1500}_p12ext_combined_1_2000.rds`; `..._h175_..._n{500,1000,1500}_tier2_combined_1_2000.rds`; `..._h100_..._n500_tier2_combined_1_2000.rds`, `..._h100_..._n{1000,1500}_p12ext_combined_1_2000.rds`.
  - 31 %: `fs_effMaxSG_fb_mr_field_m1_h{150,175}_knoise0_n500_z1q60_nb20_e1stud_combined_1_2000.rds`; `..._h{150,175,100}_..._n{1000,1500}_z1q60_nb20_cert20_combined_1_2000.rds`; `..._h100_..._n500_z1q60_nb20_cert20_combined_1_2000.rds`.
  - Null: `fs_effMaxSG_fb_mr_field_m1_h066_knoise0_n{500,1000,1500}_null657_nb20_nullmr_res_1_2000.rds`; `..._h072_..._n{500,1000,1500}_null721_nb20_nullmr_res_1_2000.rds`; `..._h066_..._n500_null657_nb20_nomr_nullid_res_1_2000.rds`; `..._h072_..._n500_null721_nb20_nomr_nullid_res_1_2000.rds`.

Every bundle is a `list` with top-level names `results, truth, meta`; every `truth` has exactly the fields `hr_causal, marg_H, marg_Hc, cde_H, cde_Hc`.

## Step 1 — Q1: the trial-wide marginal Cox hazard ratio per cell

**Where it lives.** The trial-wide marginal Cox HR is `dgm$hr_causal`, stored in every bundle as `truth$hr_causal` (template line 844, comment "overall causal HR"). It is built in `.create_gbsg_dgm_()` as `exp(coef)` of `survival::coxph(Surv(time, event) ~ treat)` fitted on the stacked potential outcomes (`T_1`, `T_0`) of every subject of the super-population (`R/sim_aft_gbsg.R:505–520`; the stacked frame at 511–515). Treatment is the only covariate, so it is the marginal (population-averaged) Cox HR over the whole trial population, not a conditional or patient-level quantity.

**Design constant, not realized.** It is computed once at DGM build, before any replicate, from the super-population drawn under `set.seed(seed)` (`R/sim_aft_gbsg.R:478`) with `seed = seed_base` and fixed `n_super` (template 829–831). It does not depend on `n`. Read from the bundles, the stored value is identical to all 15 significant digits across the three n bundles of every row below, and identical between the `nullmr` and `nullid` bundles at each null point.

**Targeted or derived.**
- `alt` (harm and attenuated-benefit cells): the calibration target is the HR **inside the planted region**, `dgm$hr_H_true` (`calibrate_k_inter()`, `R/sim_aft_gbsg.R:1028`, `use_ahr = FALSE` at template 817). `hr_causal` is recorded but not targeted: it is what the planted region at `FS_S7_HR` and the unmodified complement produce together.
- `null` (uniform-benefit cells): `hr_causal` **is** the target. `k_treat` is solved by `uniroot` so that `setup_gbsg_dgm(model = "null", ...)$hr_causal` equals `FS_S7_HR` (template 819–826), and the design-point gate asserts it before any replicate (template 929–956).

| Cell type | Prevalence | Trial-wide marginal Cox HR (`truth$hr_causal`) | Constant or realized | Source (field or line) |
|---|---|---|---|---|
| harm, planted HR 1.50 | 12.4 % (`z1q` 0.25; `harm_prevalence_super` 0.12418) | 0.704145 | design constant, derived (not targeted) | `R/sim_aft_gbsg.R:517–520`; template :844; `..._h150_..._p12ext_combined_1_2000.rds` `truth$hr_causal` |
| harm, planted HR 1.75 | 12.4 % | 0.710198 | design constant, derived | same lines; `..._h175_..._tier2_combined_1_2000.rds` |
| harm, planted HR 1.50 | 31 % (`z1q` 0.60; `harm_prevalence_super` 0.30655) | 0.869414 | design constant, derived | same lines; `..._h150_..._z1q60_nb20_{e1stud,cert20}_combined_1_2000.rds` |
| harm, planted HR 1.75 | 31 % | 0.894386 | design constant, derived | same lines; `..._h175_..._z1q60_nb20_{e1stud,cert20}_combined_1_2000.rds` |
| attenuated-benefit, planted HR 1.00 | 12.4 % | 0.684734 | design constant, derived | same lines; `..._h100_..._{tier2,p12ext}_combined_1_2000.rds` |
| attenuated-benefit, planted HR 1.00 | 31 % | 0.792179 | design constant, derived | same lines; `..._h100_..._z1q60_nb20_cert20_combined_1_2000.rds` |
| uniform-benefit (`null`), marginal HR 0.657 | no prevalence dimension (no region; `harm_prevalence_super` 0) | 0.657 (stored 0.657000012) | design constant, **targeted** (`k_treat` 1.27232145415221) | template :819–826 (calibration), :929–956 (gate); `R/sim_aft_gbsg.R:517–520`; `..._h066_..._null657_nb20_nullmr_res_1_2000.rds` `truth$hr_causal`, `meta$k_treat` |
| uniform-benefit (`null`), marginal HR 0.721 | no prevalence dimension | 0.721 (stored 0.721000000) | design constant, **targeted** (`k_treat` 0.99182531577674) | same lines; `..._h072_..._null721_nb20_nullmr_res_1_2000.rds` |

Values are the stored doubles rounded to six decimals. Full stored values, for the record: 0.704144840402614, 0.710198276661431, 0.869414265761502, 0.894385838687509, 0.684733834817233, 0.792178852282265, 0.657000011553082, 0.720999999972782.

No cell type is `not established from source`: every row carries a stored field and the line that constructs it.

Two facts about the record, stated without interpretation:
- The `alt` bundles carry `meta$target_hr_harm` and `meta$harm_z1_quantile` but **no** `meta$dgm_model`, `meta$k_inter` or `meta$k_treat` (read as NA): those fields postdate the bundles (template 1898–1903, "NA there, and read with `%||%` 'alt'"). The `null` bundles carry all three (`dgm_model = "null"`, `k_inter = 1`, `k_treat` as above).
- No bundle stores the trial-wide HR anywhere other than `truth$hr_causal`; `meta` holds the knobs, not the estimand.

## Step 2 — Q2 and Q3, confirmation against source

**Q2 — complement marginal HR 0.657 at 12.4 % and 0.721 at 31 %: confirmed.**
- `truth$marg_Hc` is `dgm$hr_Hc_true` (template :846), the same Cox-on-stacked-potential-outcomes construction restricted to `flag.harm == 0` (`R/sim_aft_gbsg.R:524, 537–540`).
- Stored: `marg_Hc = 0.656891415008813` in all nine 12.4 % bundles (`meta$harm_z1_quantile 0.25`, `meta$harm_prevalence_super 0.12418`); `marg_Hc = 0.720557356454637` in all nine 31 % bundles (`harm_z1_quantile 0.60`, `harm_prevalence_super 0.30655`). `marg_Hc` is identical across HR 1.50, 1.75 and 1.00 within a prevalence.
- Precision note, recorded as a fact: the catalog's 0.657 and 0.721 are the three-decimal roundings of 0.656891 and 0.720557. The `null` design targets 0.657 and 0.721 exactly (stored `hr_causal` 0.657000012 and 0.721000000), so the null uniform effect and the alt complement effect agree at three decimals and differ beyond that.

**Q3 — `k_inter` calibrated to the target HR inside the region, `k_treat = 1`, complement carries the base effect unmodified: confirmed.**
- `k_inter` target: template :816–818 calls `calibrate_k_inter(target_hr_harm, model = "alt", use_ahr = FALSE, z1_quantile)`; the objective reads `dgm$hr_H_true` (`R/sim_aft_gbsg.R:1028`), the marginal Cox HR inside `flag.harm == 1` (:523, 527–530). Stored `marg_H`: 1.50861 / 1.76908 / 1.00045 at 12.4 %, 1.49903 / 1.74616 / 0.99987 at 31 %.
- `k_treat = 1`: template :814 sets `k_treat <- 1` and only the `null` branch (:819–826) changes it; `calibrate_k_inter()`'s own default is `k_treat = 1` (`R/sim_aft_gbsg.R:1007`) and the template call passes none.
- Base effect unmodified on the complement: `gamma["treat"] <- k_treat * gamma["treat"]` (`R/sim_aft_gbsg.R:417`) is the identity at `k_treat = 1`; `k_inter` scales only `gamma["zh"]` (:420); `zh = treat * z1 * z3` (:287) is zero on every complement subject. The complement's treatment coefficient is therefore the GBSG-fitted `gamma["treat"]`, and `marg_Hc` differs between prevalences only because the complement is a different subpopulation.
- The `k_treat = 1` confirmation is from source, not from the stored record: the `alt` bundles do not carry `meta$k_treat` (see Step 1).

No disagreement between the catalog's §1 and source was found.

## Step 3 — Q4: fields of the stored `truth` object

Listed to the extent needed to place Q1 on one scale. Q1's values (`hr_causal`) and Q2's (`marg_Hc`, `marg_H`) are all on the same scale, the marginal Cox HR from treatment-only Cox fits on stacked potential outcomes; only `cde_*` is on a different scale.

| Field | What it holds | Scale | Source |
|---|---|---|---|
| `hr_causal` | Trial-wide HR, Cox(treat) on stacked potential outcomes of the full super-population | marginal Cox | `R/sim_aft_gbsg.R:517–520`; template :844 |
| `marg_H` | Same construction restricted to the planted region H (θ†(H)); `NA` under `null` | marginal Cox | `R/sim_aft_gbsg.R:527–530`, :547; template :845 |
| `marg_Hc` | Same construction restricted to the complement Hc (θ†(Hc)); under `null` set equal to `hr_causal` | marginal Cox | `R/sim_aft_gbsg.R:537–540`, :548; template :846 |
| `cde_H` | Controlled direct effect θ‡(H): `mean(exp(theta_1)) / mean(exp(theta_0))` over H; `NA` under `null` | patient-level potential-outcome hazard ratio, not marginal Cox | `R/oc_analyses.R:163`; template :847 |
| `cde_Hc` | θ‡(Hc), same over Hc; under `null` it is the uniform patient-level HR (stored 0.582908 at the 0.657 point, 0.656562 at the 0.721 point, matching catalog §1) | patient-level, not marginal Cox | `R/oc_analyses.R:164`; template :848 |

## Post-conditions

- P1: `git status` shows exactly the two added files (task in `dev/tasks/`, this report in `dev/reports/`) and no modified tracked file: checked at commit time.
- P2: nothing re-run, nothing computed, no payload written; md5 of all 26 bundles read identical before and after.
- P3: every Q1 row carries a file:line.
- Copy placed at `~/Downloads/gbsg020_itt_estimand_2026-09-28/`.
