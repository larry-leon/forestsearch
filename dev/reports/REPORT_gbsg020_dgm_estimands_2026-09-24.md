# REPORT — `gbsg_020` DGM estimands: the ITT effect and the prevalence pairing (read-only)

Task: `dev/tasks/TASK_gbsg020_dgm_estimands_2026-09-24.md` (committed as received, `93f0be27`). Branch `feature/glm-extension`.

**No simulation was re-run and no estimate or rate was computed.** R was used only in `Rscript --vanilla` sessions. They read the top-level `names()` and classes of committed bundles and the stored values of their `truth` and `meta` fields, and wrote nothing. Every bundle read was checksummed before and after, and all 20 were unchanged (P2).

## 0. Gates and paths read

- **G0.1.** No tracked file was modified. These untracked files were already present. They were left alone and are not part of this task:
  - `quarto/simulations/actg175/binary_020/mr_or_harm/fs_effMaxSG_mr_field_or075_n500_nb20_redes_d5000/`
  - `quarto/simulations/actg175/binary_020/mr_or_harm/fs_effMaxSG_mr_field_or075_n500_nb20_relaunch_d5000/`
  - `quarto/simulations/actg175/binary_020/smoke_redes.html`
  - `quarto/simulations/actg175/binary_020/smoke_relaunch.html`
  - `quarto/simulations/gbsg_020/scripts_dinamr/logs/nullmr_findings.err`
- **G0.2.** Repository `forestsearch`, branch `feature/glm-extension`: pass.
- **G0.3.** `quarto/simulations/gbsg_020/` exists, with committed payloads under `results/`: pass.

**Source files read.** "Template" below means the first entry. Paths are relative to `quarto/simulations/gbsg_020/` unless they start with `R/` or `~`.

| path | used for |
|---|---|
| `sim_fs_maxeffCons_fb_mr_field_m1_template.qmd` (last commit `6292e8c2`) | the campaign template; `current_status.md` line 12 calls it "One template". Relevant lines: 357–397 (knobs), 796–880 (DGM build and `truth`) |
| `R/sim_aft_gbsg.R` | `.create_gbsg_dgm_()`, from line 236; `calibrate_k_inter()`, from line 1004 |
| `R/setup_gbsg_dgm.R` | lines 88–120, the wrapper that calls `.create_gbsg_dgm_()` |
| `R/oc_analyses.R` | `compute_dgm_cde()`, lines 131–176 |
| `current_status.md` | §1, lines 11–22; §2.1, lines 26–30 |
| `~/Downloads/HANDOFF_declaration_calibration_c0_2026-09-22.md` | line 160, the "`k_treat = 1`" statement |

**Bundles read (`results/`, structure and stored constants only):**

- all 18 designated Section 5 comparator bundles, named per `current_status.md` §2.1:
  - `fs_maxeffCons_fb_mr_field_m1_h{100,150,175}_knoise0_n{500,1000,1500}_{p12ext,tier2}_combined_1_2000.rds`, the 12.4% cells;
  - `fs_effMaxSG_fb_mr_field_m1_h{100,150,175}_knoise0_n{500,1000,1500}_z1q60_nb20_{cert20,e1stud}_combined_1_2000.rds`, the 31% cells;
- uniform benefit: `fs_effMaxSG_fb_mr_field_m1_h066_knoise0_n1000_null657_nb20_nomr_nullid_res_1_2000.rds` and `..._h072_..._null721_nb20_nomr_nullid_res_1_2000.rds`;
- for the handoff check in Q3: `declcal_power_C1_res_1_2000.rds` and `declcal_power_C3_res_1_2000.rds`.

**Manuscript not read.** The `fs-glms-interpretable` manuscript is in another repository and is not a Step 1 source. Its wording is therefore not checked here. Disagreements are recorded against the handoff and against this repository's record only.

---

## Q1 — The ITT estimand

**Field.** The trial-wide marginal Cox HR is `truth$hr_causal`. Template lines 843–844 set `hr_causal = dgm$hr_causal`, commented "# overall causal HR". `R/sim_aft_gbsg.R` lines 517–520 define `dgm$hr_causal`: `exp(coxph(Surv(time, event) ~ treat, data = df_po)$coefficients)`. `df_po` is the stacked treated and control potential outcomes of the whole super-population, with every row an event (lines 509–514). This is the marginal Cox scale.

**A design constant, not a per-replicate value:**

- Template line 842: "Fixed (population) truth targets, read defensively from the DGM object."
- `truth` is a single top-level list in each bundle, next to `results` and `meta`.
- Within a cell type, the stored value is identical at n 500, 1000 and 1500, and across the `p12ext`/`tier2` and `cert20`/`e1stud` campaigns (bundle reads below).

**Stored `truth$hr_causal`:**

| cell type | planted region (target) | prevalence (`meta$harm_prevalence_super`) | `truth$hr_causal` | how it is set |
|---|---|---|---|---|
| uniform benefit | none (structural null) | 0 | **0.6570000116** (`null657`); **0.721** (`null721`) | **target**: `FS_S7_HR` solved for through `k_treat` (template lines 819–827) |
| attenuated benefit | HR 1.00 | 0.12418 | **0.6847338348** | **derived**, not targeted |
| attenuated benefit | HR 1.00 | 0.30655 | **0.7921788523** | derived |
| harm | HR 1.50 | 0.12418 | **0.7041448404** | derived |
| harm | HR 1.50 | 0.30655 | **0.8694142658** | derived |
| harm | HR 1.75 | 0.12418 | **0.7101982767** | derived |
| harm | HR 1.75 | 0.30655 | **0.8943858387** | derived |

- **Uniform-benefit cells.** `hr_causal` is the design's target. Template lines 384–387 say "FS_S7_HR is then the target SUPER-POPULATION MARGINAL Cox HR (dgm$hr_causal …)", and lines 824–826 run `uniroot` on `.hr_at(kt) - target_hr_harm`.
- **Attenuated-benefit and harm cells.** `hr_causal` is recorded but not targeted. Template lines 815–818 calibrate `k_inter` only. `k_treat` stays at 1 (line 814).
- **Per-replicate realized values.** Some campaigns record a per-replicate whole-trial fit, `itt_est` (`current_status.md` line 170). That is an estimate, not the DGM's estimand. This report reads no per-replicate column.

**Disagreement with the record.** The design-level ITT HR of the planted-region cells, `truth$hr_causal`, is **not quoted anywhere in `current_status.md`**. A search for `hr_causal`, "overall causal" and the stored values finds `hr_causal` only as the null design's target (line 18). §1 (lines 13–15) names only the prevalences, the planted-region HR and the complement HRs 0.657 and 0.721. The value exists only in the bundles' `truth` field, and the template prints it at line 852 ("overall causal HR").

---

## Q2 — The prevalence pairing

**12.4% gives `marg_Hc` 0.657 and 31% gives 0.721.**

| prevalence | `harm_z1_quantile` | `truth$marg_Hc` (stored) | cells where it is stored |
|---|---|---|---|
| 12.4% (`meta$harm_prevalence_super` 0.12418) | 0.25 | **0.6568914150** | all nine 12.4% cells: HR 1.00 / 1.50 / 1.75 × n 500 / 1000 / 1500 |
| 31% (0.30655) | 0.60 | **0.7205573565** | all nine 31% cells: HR 1.00 / 1.50 / 1.75 × n 500 / 1000 / 1500 |

This matches `current_status.md` line 15: "HR 0.657 at 12.4%, 0.721 at 31%".

**Where it is defined.** Template line 846 has `marg_Hc = dgm$hr_Hc_true`, commented "# theta-dagger (Hc)". `R/sim_aft_gbsg.R` lines 536–540 define it for `model == "alt"` as the marginal Cox HR on the stacked potential outcomes restricted to `flag.harm == 0`.

**Which cell types have it:**

- **Both harm and attenuated-benefit cells.** Both run `model = "alt"`, so it is the same field computed the same way.
- **Within a prevalence, it is one number.** It is bit-identical across the three planted-region HRs and the three sample sizes, so the harm and attenuated-benefit cells carry the same value.
- **Uniform-benefit cells.** `marg_Hc` is also stored, with a different definition. Under `model = "null"`, `R/sim_aft_gbsg.R` line 548 sets `hr_Hc_true <- hr_causal`, so there `marg_Hc` equals the whole-trial HR: `0.6570000116` and `0.721`.

**Disagreement to record.** The alt design's complement values and the null design's targets are different objects, and they differ past the third decimal:

- `marg_Hc` at 12.4% is 0.6568914150; the `null657` target is 0.657 (stored 0.6570000116);
- `marg_Hc` at 31% is 0.7205573565; the `null721` target is 0.721.

`current_status.md` line 19 says the null points are "fixed by their marginal Cox HR (0.657, 0.721), because that is what the alt design's complement effects are". The null designs were solved to the 3-decimal labels, not to the stored `marg_Hc` values.

---

## Q3 — What the DGM holds fixed

**What is held.** The **planted region's marginal Cox HR**. `calibrate_k_inter()` (`R/sim_aft_gbsg.R` lines 1004–1045) solves `k_inter` so that `dgm$hr_H_true` equals `target_hr_harm`. The template calls it with `use_ahr = FALSE` (template lines 816–818), so the target is the marginal Cox HR inside the region, not the AHR. `k_inter` scales only `gamma["zh"]` (line 420). `zh` = treat × z1 × z3 is zero for every complement subject (lines 287, 359, 364).

**What is derived:**

- **The trial-wide effect** (`hr_causal`). It is not pinned. At a fixed prevalence it moves with the planted HR: at 12.4%, 0.6847 / 0.7041 / 0.7102 at HR 1.00 / 1.50 / 1.75 (Q1 table).
- **The complement's effect** (`marg_Hc`). It is not pinned or solved for either. Nothing in the template or in `.create_gbsg_dgm_()` targets it.

**Is the complement effect "the base treatment effect left unmodified"?**

- **Yes, in the sense of the multiplier.** `k_treat = 1` (template line 814), and `gamma["treat"] <- k_treat * gamma["treat"]` (line 417) leaves it unscaled.
- **But the base coefficient itself depends on the prevalence.** `gamma` comes from a Weibull AFT fit to the GBSG data (lines 375–385), with design matrix `covs_true = c("treat", "z1", ..., "zh")` (line 335). `z1` is defined from `z1_quantile` (lines 279–280), which is also what sets the prevalence. The fit is therefore re-run on a different `z1` and `zh` at each prevalence, so `gamma["treat"]`, and with it `b0["treat"] = -gamma["treat"]/sigma` (line 441), can differ between 12.4% and 31%.
- **The stored fields show that it does differ.** For a complement subject, `theta_1 - theta_0 = b0["treat"]` (lines 467–469, with `zh` = 0). So `cde_Hc = mean(exp(theta_1[Hc])) / mean(exp(theta_0[Hc]))` (`R/oc_analyses.R` line 164) is the complement's patient-level HR, exp(b0["treat"]). Its stored value is **0.5847778746 at 12.4%** and **0.6563949895 at 31%**.
- **So `marg_Hc` does not differ between prevalences "only because the complement is a different subpopulation".** The patient-level complement effect itself also differs between the two prevalences, through the refit of the base AFT model.
- **Neither option in the question holds as stated.** The trial-wide effect is not pinned, and the complement is not solved for.

**Disagreement with the handoff.** HANDOFF line 160 says: "C's complement runs at `k_treat = 1`; B5 matches on the patient-level scale, B2 on the marginal Cox scale".

- **"`k_treat = 1`" is confirmed.** `meta$k_treat` = 1 in `declcal_power_C1` and `C3`.
- **"B2 on the marginal Cox scale" is confirmed.** Block C's stored `marg_Hc` is 0.6568914150 in both C1 (HR 1.5) and C3 (HR 2.0). That is the 12.4% complement value (their `meta$harm_prevalence_super` is 0.12418), and it matches B2's marginal target of 0.657.
- **"B5 matches on the patient-level scale" does not match the source.**
  - B5's patient-level HR is 0.656562, as `cde_Hc` of the `null721` bundle (0.6565617728).
  - The 12.4% complement's patient-level HR (`cde_Hc`) is 0.5847778746 in the Section 5 12.4% bundles. The `declcal` Block C `truth` carries no `cde_*` fields, so it cannot be read there.
  - The number B5's patient-level HR is close to is the 12.4% complement's **marginal** HR (0.656891), a comparison across two different scales.

---

## Q4 — The truth object

The same five-field list is stored as the bundles' top-level `truth`. It is built at template lines 843–849.

| field | holds | scale | source |
|---|---|---|---|
| `hr_causal` | Cox HR of treatment on the stacked potential outcomes of the whole super-population: the trial-wide (ITT) estimand | **marginal Cox** | template line 844; `R/sim_aft_gbsg.R` lines 517–520 |
| `marg_H` | the same Cox HR restricted to the planted region (`flag.harm == 1`); the calibration target of `k_inter`. `NA` under `null` | **marginal Cox** (θ†, H) | template line 845; `R/sim_aft_gbsg.R` lines 526–533 and 547 |
| `marg_Hc` | the same Cox HR restricted to the complement (`flag.harm == 0`). Under `null` it is set equal to `hr_causal` | **marginal Cox** (θ†, Hᶜ) | template line 846; `R/sim_aft_gbsg.R` lines 536–543 and 548 |
| `cde_H` | mean(exp θ₁)/mean(exp θ₀) over the planted region: a ratio of average potential-outcome hazard multipliers. `NA` under `null` | **other: controlled direct effect** (θ‡, H) | template line 847; `R/oc_analyses.R` lines 163, 175 |
| `cde_Hc` | the same over the complement. θ₁ − θ₀ is constant there (b0["treat"]), so it equals the complement's uniform **patient-level** HR | **other: CDE** (θ‡, Hᶜ). Numerically the patient-level HR, by construction | template line 848; `R/oc_analyses.R` lines 164, 176; `R/sim_aft_gbsg.R` lines 467–469 |

**Notes:**

- **`truth` does not carry the AHR** (`dgm$AHR`, exp of the mean `loghr_po`; `R/sim_aft_gbsg.R` line 556), nor `AHR_H_true` / `AHR_Hc_true`. `hazard_ratios$AHR*` exists on the DGM object (lines 665–672), but the template's `truth` list does not copy it.
- **The prevalence is not in `truth`.** It is in `meta$harm_prevalence_super` (0.12418 / 0.30655), next to `meta$harm_z1_quantile` (0.25 / 0.60) and `meta$target_hr_harm`.
- **The `declcal` payloads carry a three-field `truth` under `meta$truth`**: `hr_causal`, `marg_H` and `marg_Hc`, with no `cde_*`.
