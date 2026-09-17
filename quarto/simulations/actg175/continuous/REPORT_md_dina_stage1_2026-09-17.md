# REPORT — DINA on the ACTG175 continuous (MD) design: campaign `mddina`, Stage 1 (edits, smoke, calibration → Gate 1)

Date: 2026-09-17. Machine: `pop-os` (64 physical cores, 251 GB; R 4.6.1, reference BLAS). Branch `feature/glm-extension`. Task: `dev/tasks/TASK_md_dina_campaign_2026-09-17.md` (committed as received, `189d4ec2`), which transplants the mechanics of `dev/tasks/TASK_md_grf_2026-09-16.md` ("the GRF task") and `TASK_md_grf_resume_2026-09-16.md`. References: `REPORT_grf_dina_fixes_2026-09-16.md` ("the fix record"; P2 landed, `064fce91`), `REPORT_md_dina_grf_stage0_2026-09-16.md` ("S0"), the `mdgrf` records and the `mdsgnb20` campaign. Installed forestsearch 0.3.5, `Built: R 4.6.1; ; 2026-09-17 04:47:31 UTC; unix`; no `R/` change, no install. Every render and R process ran with `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1`. Line numbers: `R/` at HEAD (unchanged since `064fce91`); the MD template at `0927669c` unless stated; the survival template `m1` = `quarto/simulations/gbsg_020/sim_fs_maxeffCons_fb_mr_field_m1_template.qmd` at HEAD.

**Governing constraint.** DINA's candidate family is generated from a fitted surface, so the fixed-family condition does not hold. Every coverage figure this campaign produces is coverage of the estimand conditional on the proposed family; comparisons with FS and GRF are descriptive, not a contest, and state the confound (identifier, family construction, detected set).

**Outcome.** Every Stage 1 gate passes (§1.1, §1.2, §1.3, §1.6(a), (b), (c)). The §1.7 projection for Stage 2 at W = 63 is **6.10 h**, under 8 h, so **the advance go applies**: Stage 2 launches after this record's commit with `MDSG_WORKERS=63 MDSG_TIMEOUT=6864 MDSG_CEILING=32939`.

## 1.1 Provenance and first commit — GATE PASS

```
pop-os
feature/glm-extension
853c43c6
853c43c6 actg175/continuous closeout: generate current_status.md at 688d4b68
688d4b68 actg175/continuous status_curated.md: campaign mdgrf (payloads, records, the conditional-family reading convention, where to start); ...
94d9cf53 mdgrf record (TASK_md_grf_2026-09-16 §3.3): GRF on the MD design, four cells x 2,000 replicates, ...
[tracked modifications: none]
mdgrf closeout in HEAD
2086198 11-11:27:12 /bin/bash -c source .../snapshot-bash-...sh ... eval 'cd ~/.../gbsg_020 && export PATH="/usr/lib/rstudio/resources/app/bin/quarto/bin:$PATH" && ... ( while true; do free -m | awk ...; sleep 3; done > "$SCR/mem_uburst.txt" ) & MP=$!; ... quarto render ...'
```
- The one line the `ps | grep` pattern returns is a bash subshell started 11 days earlier by another session: a memory-sampling loop (`free -m` / `sleep 3`) whose command text contains the word `quarto`. Its only child is `sleep 3`; `pgrep` finds no `R`, `Rscript` or `quarto` process. No R, Rscript or quarto process was running (finding F1). The loop was left alone.
- First commit: `189d4ec2 dev/tasks: TASK_md_dina_campaign_2026-09-17.md as received (campaign mddina, Stages 1-3)`.

## 1.2 Package check — GATE PASS

- `git diff --quiet 019be60f..HEAD -- R/ DESCRIPTION NAMESPACE` succeeds (the last `R/` commit is `064fce91`, an ancestor of `019be60f`).
- `packageDescription("forestsearch")$Built` = `R 4.6.1; ; 2026-09-17 04:47:31 UTC; unix` in the main session and in two doFuture multisession workers (pids 4167955, 4167956).

## 1.3 Gate 0 from the records — PASS

1. **Identifier path.**
   - Template knob lines (`0927669c`): `:139` `subgroup_method <- .env_chr("FS_MD_METHOD", "consistency")  # "consistency" | "dina" | "grf"`; `:140` `stopifnot(subgroup_method %in% c("consistency", "dina", "grf"))`; `:698` `subgroup_method = subgroup_method,` (in `base_args`); `:720–721` `dina = list(dina_select_statistic = dina_select_statistic, dina_args = dina_args),` (the `method_args` switch); `:282–283` `dina_args <- list()`, `dina_select_statistic <- "effect"`. The DINA value is `FS_MD_METHOD=dina`. Under it `use_dina` is not passed: it belongs to the `consistency` arm only (`:714–719`).
   - F5 in the fix record (`:68`): *"DINA on the MD design, md40 n500, sim_id 1, `adverse_outcome = FALSE` (the template's), S0.3's settings (`FS_MD_METHOD=dina`, `dina_args = list()`, `dina_select_statistic = "effect"`, `effMaxSG`, ε 0.20, `neighborhood`), `details = TRUE`; MR off."*, dispatched as *"`do.call(forestsearch, c(base_args, method_args))`"* from the template's own block (`:62`).
   - From source: `R/forestsearch_main.R:2238–2239` `if (subgroup_method == "dina") { dsel <- .forestsearch_dina_select(` … and the section ends with `return(out)` (`:2355`), before the consistency search. `.forestsearch_dina_select()` computes `tau_sign <- .dina_tau_sign(outcome_type, adverse_outcome)` (`R/forestsearch_helpers.R:1441`) and passes `tau_sign = tau_sign` to `dina_subgroup()` (`:1496`). These are the lines P2 changed. `.dina_tau_sign()` (`:1376`) returns −1 for `continuous` with `adverse_outcome = FALSE`. The `use_dina` screening path ("SECTION 3B: DINA CUT GENERATION (if use_dina = TRUE)", `main:2649`, `if (use_dina && is.null(dina_cuts))` `:2657`, `m_diff_sel` `:2697`, the `dina_subgroup()` call `:2698–2709` without `tau_sign`) lies in the consistency section, which the DINA branch never reaches. **The DINA value runs DINA's identifier path.** Not a stop.
2. **Floors, F5 after P2** (fix record `:461–467` and its table `:498–504`): proposal floor *"`m_diff = 30` on −tau-hat (harm)"* (details line *"Harm floor: m_diff = 30.0000  on -tau-hat (harm-oriented; adverse_outcome = FALSE)"*); admission floor *"30 on the harm-oriented MD"*; searched / proposed 8,324 / 2,690; admitted 1,768; *"tau_hat (as ranked) min 30.0064 max 111.3391 | all >= m_diff 30: TRUE"*.
3. **Forwarding (S0.4).** The DINA branch calls `.fs_apply_mr(..., mr_inference_args = mr_inference_args, seedit = seedit)` (`main:2338–2349`). In `.fs_apply_mr()` (`R/fs_mr_inference_methods.R:141`) each argument in the dispositions is forwarded from `mr_inference_args`: `ci_method` (`:171`, fallback `"ij"`), `draws` (`:168`, fallback 2000L), `include_complement` (`:170`), `confirm_rule` (`:164`), `t_confirm` (`:163`, passed through), `return_reselection` (`:176`), `field_uniform` (`:180`), `field_complement` (`:182`), `field_scale_complement` (`:186`), `ij_residual` (`:188`). None is hard-coded. The fallbacks differ from the FS branch's (`ci_method` `"ij"` against `"field"`, `main:3421`), but the template passes all of these explicitly (`:342–348`: `ci_method`, `draws = 5000L`, `include_complement = TRUE`, `confirm_rule`, `field_uniform`, `field_complement`, `ij_residual`, `field_scale_complement`, `return_reselection`; `t_confirm` NULL, so not passed), so no fallback applies. Not a stop.
4. **Re-selection (S0.4).** DINA side: `.dina_reselect_on_effect()` orders the admitted set under `hrMaxSG` (the canonical form of `effMaxSG`) by `.compute_inclusion_band(hr_vec = eff[ok], n_vec = sz[ok], selection_rule, effect_neighborhood)` and then size (`R/forestsearch_helpers.R:1251–1256`), on the effect estimator's harm-oriented MD. MR side: `reselection_default = .fs_mr_reselection_from_focus(sg_focus, engine = "effect")` (`main:2347`), which maps `effMaxSG = "effMaxSG"` (`R/fs_mr_inference_methods.R:104`), with `effect_neighborhood` and `selection_rule_default = selection_rule` (`main:2346–2348`), banded by the same `.compute_inclusion_band()` in `.fs_mr_select()`. S0 §5.4: *"The engines' effect re-selection (5.1) uses the same focus, the same band helper and the same harm-oriented MD scale: **aligned**"*. **Aligned** for DINA on this outcome. Not a stop.
5. **Recorder fields.** S0.5 (§6.1) lists *"`n_family`, `n_cons_qual`, `band_n` …, `admitted_n` …, `p_hat_sum`, `p_hat_top1..3` …, the nine `fld_recov_*`"*. After the GRF task's E3 (`894da993`) the template carries all of them (`p_hat_top1..3` as `p_top1..3` / `p_lab1..3`). **The template lacks none of them.** In `m1` the `admitted_n` fill is GRF-only (*"NA_integer_ on every non-GRF path, so consistency and DINA records are unchanged"*, m1 `:997`), and the committed `dinamr` bundles carry `n_cons_qual` and `band_n` all-NA on DINA (0 of 1,049 detected rows finite in `dina_effMaxSG_fb_mr_field_m1_h100_knoise0_n1000_nb20_dinamr_combined_1_2000.rds`). E2 (§1.5) therefore adds DINA's proposal and admission counts, which §1.6(c) needs per replicate and which have no survival precedent (finding F3).

## 1.4 The survival reference (read only)

Sources: **MD** = the MD template at `189d4ec2` (= `894da993`); **m1** at HEAD; **runner** = `quarto/simulations/gbsg_020/scripts_dinamr/campaign.sh:9–11` (`KN=(FS_S7_METHOD=dina FS_S7_FOCUS=effMaxSG FS_S7_NBHD=0.20 FS_S7_FIELD_COMPLEMENT=TRUE FS_S7_FIELD_SCALEC=selected FS_S7_FIELD_DECOMP=TRUE FS_S7_FIELD_RECOV=TRUE FS_S7_IJ_RESIDUAL=two_term FS_S7_FB=none FS_S7_CAMPAIGN=dinamr FS_S7_WORKERS=12)`, with `FS_S7_ER_JCUTS` unset by `env -u`; `README.md`: *"`FS_S7_ER_JCUTS` is deliberately unset — inert on DINA"*); **meta** = `results/dina_effMaxSG_fb_mr_field_m1_h100_knoise0_n1000_nb20_dinamr_combined_1_2000.rds` (`subgroup_method dina`, `sg_focus effMaxSG`, `focus_tag effMaxSG`, `effect_neighborhood 0.2`, `ci_method field`, `field_complement TRUE`, `field_decompose TRUE`, `field_scale_complement selected`, `field_recovery TRUE`, `ij_residual two_term`, `mr_draws 5000`, `seed_base 8316951`, `consistency_method resample`, `stop_threshold NULL`, `forestsearch_version 0.3.5`; no `dina_args`, `dina_select_statistic` or `selection_rule` key).

| Setting | MD template (`189d4ec2`) | Survival DINA campaign `dinamr` (m1 / runner / meta) | Campaign takes |
|---|---|---|---|
| DINA's selection through the knob | `FS_MD_METHOD` (`:139–140`), value `dina` → `subgroup_method = "dina"` in `base_args` (`:698`) → `main:2238` | `FS_S7_METHOD` (m1 `:310–311`), runner `FS_S7_METHOD=dina`; meta `subgroup_method dina` | `FS_MD_METHOD=dina` |
| `dina_args` | `list()` (`:282`), passed in the `dina` arm (`:721`) | `list()` (m1 `:537`, arm `:1066`); resolves to family `gaussian` on MD / `cox` on survival (`.map_dina_family()`, `helpers:1011–1019`), `seed = seedit`, `selected_only` TRUE, `max_depth` 2, `grid_probs` 0.1–0.9 (`.resolve_dina_args()`, `helpers:1037–`) | `list()` |
| `dina_select_statistic` | `"effect"` (`:283`, arm `:720`) | `"effect"` (m1 `:538`, arm `:1065`) | `"effect"` |
| survival-only DINA arguments | — | none passed: `dina_args` is empty, so `cens_type` / `cens_params` (survival-only fit keys, `helpers:1043`, `:1121`) take their defaults; `FS_S7_ER_JCUTS` (a survival cut-grid knob) is unset by the runner | none to omit |
| `use_dina` | `FALSE` (`:269`), consistency arm only (`:717`) | `FALSE` (m1 `:585`), consistency arm only (`:1062`) | not passed on the DINA path |
| proposal floor | `effect.threshold = 30` (`:257`, `:694`), aliased into `hr.threshold` (`main:1352`); `m_diff = hr.threshold` for gaussian (`helpers:1436`) on `tau_sign × tau_hat` with `tau_sign = −1` (`:1441`) | `hr.threshold = 0.90` (m1 `:568`); `m_diff = log(hr.threshold)` for cox, `tau_sign = 1` | the MD template's (no knob; follows the effect threshold) |
| admission floor | `effect_floor = hr.threshold` = 30 on the harm-oriented MD (`helpers:1222–1224`, `admission$effect_floor`) | `effect_floor` on log HR | the MD template's |
| how DINA uses the band | a sort key: `hrMaxSG` orders the admitted set by `.compute_inclusion_band()`, then size (`helpers:1251–1256`); natively `dina_subgroup()` mirrors it on oriented `tau_hat` | same code | package |
| `n.min` | `60L` (`:262`) | `NULL` (m1 `:573`; the survival default) | MD's |
| focus / band / rule | `FS_MD_FOCUS` / `FS_MD_NBHD` / `"neighborhood"` (`:245`) | `FS_S7_FOCUS` / `FS_S7_NBHD` / `"neighborhood"` literal (m1 `:547`); runner `effMaxSG`, `0.20`; meta `effMaxSG`, `0.2` | `effMaxSG`, `0.20`, `neighborhood` |
| `vi.grf.min` | `-0.2` (`:270`, `:697`); inert on the DINA path (the variable-importance screen is in the consistency section, `main:2848–2916`, after DINA returns at `:2355`) | not passed | MD's (inert) |
| MR arguments on the DINA path | `mr_inference_args` (`:342–348`): `ci_method = "field"` (`FS_MD_CI`), `draws = 5000L`, `include_complement = TRUE`, `confirm_rule = "point"`, `field_uniform = FALSE`, `field_complement = TRUE`, `ij_residual = "two_term"`, `field_scale_complement = "selected"`, `return_reselection = TRUE`; `t_confirm` NULL | m1 `:683–691`: the same set plus `field_decompose = mr_field_decompose` (`:648`) and `field_recovery = mr_field_recovery` (`:668`); runner sets both TRUE; meta `field_decompose TRUE`, `field_recovery TRUE` | MD's (finding F4) |
| DINA seeding | `seedit = sd_i` (`:704`) → `.resolve_dina_args(seed_default = seedit)` | `seedit = seed_base + sim_id` (m1 `:1044`) | same |
| workers | `FS_MD_WORKERS` (cap physical − 1) | `FS_S7_WORKERS=12` (Mac) | smoke 20; campaign §1.7 |

**Differences, as findings:** `n.min` (60 against the survival default); `field_decompose` and `field_recovery` (TRUE in `dinamr`, not passed here, as in `mdsgnb20` and `mdgrf`); `vi.grf.min` passed and inert; the floors sit on different scales by outcome (log HR against harm-oriented MD); the `dinamr` meta records no DINA argument or `selection_rule` key. The DINA arguments proper (`dina_args`, `dina_select_statistic`) are identical.

## 1.5 Template and generator edits — no `R/`

**Template, commit `0927669c`** (bottom-up; additions only, `git diff --stat` 30 insertions, 1 deletion — the deletion is the closing `))` of the poolability key vector, re-added after the new keys; every R chunk parses):
- **E1 — DINA argument block: no edit.** The §1.4 DINA arguments are already in the template from `894da993`: `:282–283` (m1 `:537–538` verbatim: `dina_args <- list()          # extra args forwarded to the DINA selector (dina_args =)` / `dina_select_statistic <- "effect"  # only if subgroup_method == "dina": "dina" | "effect"`), applied only under the DINA value by the `method_args` switch `:720–721` (m1 `:1065–1066` verbatim). No floor knob was added.
- **E2 — recorder** (DINA-guarded; no m1 precedent, finding F3): `.na_record()` gains `dina_searched_n`, `dina_proposed_n`, `dina_tau_min` after `id_secs` (`:630–634`). The fill (`:751–767`, after the GRF `admitted_n` fill and before the no-detection return) reads `fs.est$grp.consistency$out_sg` when `subgroup_method == "dina"`: `admitted_n <- out_sg$admitted_n` (set by `.dina_reselect_on_effect()`, `R/forestsearch_helpers.R:1238`), `dina_searched_n <- n_candidates_searched`, `dina_proposed_n <- n_candidates_qualifying` (`R/dina_subgroup.R:627`), `dina_tau_min <- min(candidates$tau_hat)`, i.e. the smallest oriented tau-hat among the proposed. The selection object exists only when DINA selects (the no-subgroup contract returns `grp.consistency = NULL`, `helpers:1546–1548`), so these fields are NA on DINA non-detections by construction.
- **E3 — `meta` and poolability keys:** `meta` gains `dina_select_statistic` and `dina_args` (deparsed to the string `"list()"`) after the GRF keys (`:1103–1107`); the combine-mode key vector gains the same two (`:1163–1164`). `subgroup_method` was already both.
- **Default-path check:** `git diff` shows the three new recorder columns (NA on the default path), a fill guarded by `identical(subgroup_method, "dina")`, and two `meta` / poolability keys that every path records. The `meta` additions are unguarded, as `894da993`'s GRF keys were (finding F5). §1.6(a) and (b) show the FS and GRF paths unchanged to the last bit.

**Catalog generator, commit `7e26ad42`.** The generator is `scripts_mdsgnb20/current_status_regen.R`, with its inventory table in `scripts_mdsgnb20/status_inventory.R`; the task's `<dir>/current_status_regen.R` does not exist (finding F2). `current_status_regen.R`'s campaign list gains `mdgrf` and `mddina`, copied from `mdsgnb20`'s row with the tag, the stem prefix (`grf_` / `dina_`), the scripts directory and the document pattern changed. The git-log spec reads the stem prefix (it was a literal `fs_`). Two sentences that named "both" / "the two" campaigns now read "every" / "the". `status_inventory.R` gains the `mdgrf` and `mddina` rules: bundles, scripts, logs, combine renders, extracts, records, heartbeat and halt file. The `mddina` record patterns exclude `REPORT_md_dina_grf_stage0_*`. A trial run lists the `mdgrf` rows (4/4 combined, 8/8 batch, 4/4 renders, the extract).

## 1.6 Scripts — commit `d65f3b9b`, under `scripts_mddina/` (transplants of `scripts_mdgrf/`)

`mem_sampler.sh` (unchanged but for its header); `run_mddina.sh` (the GRF runner with `FS_MD_METHOD=dina` in place of `FS_MD_METHOD=grf FS_MD_DMIN_GRF=30`, the `dina_` stem, `TAG=mddina`, `DATE=2026-09-17`, `REPORT_md_dina_gate2_2026-09-17.md`, the conditional-family sentence naming DINA, and **no `TRAILER`**: `commit()` uses `-m "$msg"`, the carried fix); `gate2.R` (the GRF Gate 2 pointed at the DINA bundles; `meta` checks `subgroup_method dina`, `dina_select_statistic effect`, `dina_args list()` in place of the four GRF keys; the E2 checks — fields filled on every declared replicate, `dina_tau_min ≥ 30`, `1 ≤ admitted_n ≤ dina_proposed_n` — with the E2 fields on non-detections reported as structural NA; `n_family` relabelled as MR's family); `smoke_identity.R` (modes `identity`, `grfreg` — the identity comparison against `mdgrf`'s combined bundle — and `dina`). The smoke and calibration drivers and the floor check live untracked in `logs_mddina/`.

## 1.6 Smoke — (a) PASS, (b) PASS, (c) PASS

Renders (`logs_mddina/smoke.sh`; md40 n500, sim_id 1–20, 20 workers, `FS_MD_FOCUS=effMaxSG FS_MD_NBHD=0.20 FS_MD_FIELD_SCALEC=selected FS_MD_CI=field FS_MD_FB=none`; template `0927669c`):
```
SMOKE fs: rc=0 wall_s=92 peak_mb=18326 2026-09-17T08:14:07Z
SMOKE grf: rc=0 wall_s=81 peak_mb=15028 2026-09-17T08:15:28Z
SMOKE dina: rc=0 wall_s=167 peak_mb=21994 2026-09-17T08:18:15Z
SMOKE DONE
```

**(a) FS regression** (`scripts_mddina/smoke_identity.R 40 500 mddinasmokefs 20 identity`, tag `mddinasmokefs`, knob at its default):
```
== SMOKE identity: md=40 n=500 tag=mddinasmokefs | fs_effMaxSG_mr_field_md40_knoise0_n500_nb20_mddinasmokefs_res_1_20.rds vs mdsgnb20 sim_id 1-20 ==
  [PASS] smoke bundle exists: mr_md_harm/fs_effMaxSG_mr_field_md40_knoise0_n500_nb20_mddinasmokefs_d5000/fs_effMaxSG_mr_field_md40_knoise0_n500_nb20_mddinasmokefs_res_1_20.rds
  [PASS] mdsgnb20 bundle exists: mr_md_harm/fs_effMaxSG_mr_field_md40_knoise0_n500_nb20_mdsgnb20_d5000/fs_effMaxSG_mr_field_md40_knoise0_n500_nb20_mdsgnb20_combined_1_2000.rds
  [PASS] both carry sim_id 1-20
  [PASS] no CONFIG-ERROR rows (status: DETECTED 20)
  new columns in this template, reported apart (19): n_family, n_cons_qual, band_n, admitted_n, p_hat_sum, fld_recov_sens_H, fld_recov_ppv_H, fld_recov_sens_Hc, fld_recov_npv_Hc, fld_recov_q10, fld_recov_q50, fld_recov_q90, fld_recov_share1, fld_recov_n_used, warn_msg, id_secs, dina_searched_n, dina_proposed_n, dina_tau_min
  [PASS] all 137 paired mdsgnb20 columns present in the smoke bundle
  pairing: 20 of 20 rows identical on 126 numeric + 11 character columns (max rel diff among identical rows 0.00e+00)
  classification: no enumerated rows
  [PASS] zero selection flips (0 enumerated rows, all in a mdsgnb20 Gate 2 class)
  reported: mr_msg differs on 0 of 20 rows
  reported: err_msg differs on 0 of 20 rows
  reported: fb_err differs on 0 of 20 rows
  reported: FB columns (13) -- finite fb_H_est: mdsgnb20 0 rows, smoke 0
  reported (new columns on this path): n_family non-NA 20 | n_cons_qual non-NA 20 | band_n non-NA 20 | admitted_n non-NA 0 | p_hat_sum non-NA 20 | fld_recov_sens_H non-NA 0 | fld_recov_ppv_H non-NA 0 | fld_recov_sens_Hc non-NA 0 | fld_recov_npv_Hc non-NA 0 | fld_recov_q10 non-NA 0 | fld_recov_q50 non-NA 0 | fld_recov_q90 non-NA 0 | fld_recov_share1 non-NA 0 | fld_recov_n_used non-NA 0 | warn_msg non-NA 0 | id_secs non-NA 20 | dina_searched_n non-NA 0 | dina_proposed_n non-NA 0 | dina_tau_min non-NA 0
  [PASS] stem and meta carry the consistency engine (knob at its default)
  [PASS] nine fld_Hc_*_s and nine fld_joint_s_* columns present
  [PASS] nine fld_Hc_*_s finite on all 20 filled replicates
  [PASS] nine fld_joint_s_* finite on all 20 filled replicates
  [PASS] fld_Hc_lo1s_s <= fld_Hc_up1s_s and fld_Hc_lo2s_s <= fld_Hc_hi2s_s
  [PASS] Bonferroni harm bound identical between joint and joint_s where draw counts agree (20 of 20 rows agree; max |diff| 0.00e+00)
  [PASS] field-s inverted around the same beta-tilde^c (max |diff| 1.78e-15)
  [PASS] meta field_scale_complement = selected
SMOKE identity md=40 n=500: PASS
```
20 of 20 rows identical to `mdsgnb20` on 126 numeric + 11 character columns (max relative difference 0), no enumerated row, zero selection flips. **PASS.**

**(b) GRF regression** (`… mddinasmokegrf 20 grfreg`, tag `mddinasmokegrf`, `FS_MD_METHOD=grf FS_MD_DMIN_GRF=30` with `mdgrf`'s knobs):
```
== SMOKE grfreg: md=40 n=500 tag=mddinasmokegrf | grf_effMaxSG_mr_field_md40_knoise0_n500_nb20_mddinasmokegrf_res_1_20.rds vs mdgrf sim_id 1-20 ==
  [PASS] smoke bundle exists: mr_md_harm/grf_effMaxSG_mr_field_md40_knoise0_n500_nb20_mddinasmokegrf_d5000/grf_effMaxSG_mr_field_md40_knoise0_n500_nb20_mddinasmokegrf_res_1_20.rds
  [PASS] mdgrf bundle exists: mr_md_harm/grf_effMaxSG_mr_field_md40_knoise0_n500_nb20_mdgrf_d5000/grf_effMaxSG_mr_field_md40_knoise0_n500_nb20_mdgrf_combined_1_2000.rds
  [PASS] both carry sim_id 1-20
  [PASS] no CONFIG-ERROR rows (status: DETECTED 20)
  new columns in this template, reported apart (3): dina_searched_n, dina_proposed_n, dina_tau_min
  [PASS] all 152 paired mdgrf columns present in the smoke bundle
  pairing: 20 of 20 rows identical on 140 numeric + 12 character columns (max rel diff among identical rows 0.00e+00)
  classification: no enumerated rows
  [PASS] zero selection flips (0 enumerated rows, all in a mdgrf Gate 2 class)
  reported: mr_msg differs on 0 of 20 rows
  reported: err_msg differs on 0 of 20 rows
  reported: fb_err differs on 0 of 20 rows
  reported: FB columns (13) -- finite fb_H_est: mdgrf 0 rows, smoke 0
  reported (new columns on this path): dina_searched_n non-NA 0 | dina_proposed_n non-NA 0 | dina_tau_min non-NA 0
  [PASS] stem and meta carry the GRF engine with mdgrf's arguments (dmin_grf 30, frontier, depth 2, effect)
  [PASS] nine fld_Hc_*_s and nine fld_joint_s_* columns present
  [PASS] nine fld_Hc_*_s finite on all 20 filled replicates
  [PASS] nine fld_joint_s_* finite on all 20 filled replicates
  [PASS] fld_Hc_lo1s_s <= fld_Hc_up1s_s and fld_Hc_lo2s_s <= fld_Hc_hi2s_s
  [PASS] Bonferroni harm bound identical between joint and joint_s where draw counts agree (20 of 20 rows agree; max |diff| 0.00e+00)
  [PASS] field-s inverted around the same beta-tilde^c (max |diff| 7.11e-15)
  [PASS] meta field_scale_complement = selected
SMOKE grfreg md=40 n=500: PASS
```
20 of 20 rows identical to `mdgrf` on 140 numeric + 12 character columns (max relative difference 0), no enumerated row, zero selection flips. **PASS.**

**(c) DINA** (`… mddinasmoke 20 dina`, tag `mddinasmoke`, `FS_MD_METHOD=dina`, 20 workers):
```
== SMOKE dina: md=40 n=500 tag=mddinasmoke | dina_effMaxSG_mr_field_md40_knoise0_n500_nb20_mddinasmoke_res_1_20.rds vs mdsgnb20 sim_id 1-20 ==
  [PASS] smoke bundle exists: mr_md_harm/dina_effMaxSG_mr_field_md40_knoise0_n500_nb20_mddinasmoke_d5000/dina_effMaxSG_mr_field_md40_knoise0_n500_nb20_mddinasmoke_res_1_20.rds
  [PASS] mdsgnb20 bundle exists: mr_md_harm/fs_effMaxSG_mr_field_md40_knoise0_n500_nb20_mdsgnb20_d5000/fs_effMaxSG_mr_field_md40_knoise0_n500_nb20_mdsgnb20_combined_1_2000.rds
  [PASS] both carry sim_id 1-20
  [PASS] no CONFIG-ERROR rows (status: DETECTED 20)
  new columns in this template, reported apart (19): n_family, n_cons_qual, band_n, admitted_n, p_hat_sum, fld_recov_sens_H, fld_recov_ppv_H, fld_recov_sens_Hc, fld_recov_npv_Hc, fld_recov_q10, fld_recov_q50, fld_recov_q90, fld_recov_share1, fld_recov_n_used, warn_msg, id_secs, dina_searched_n, dina_proposed_n, dina_tau_min
  [PASS] n_true identical on every row
  [PASS] oracle columns (or_H_est,or_H_lo,or_H_hi,or_H_se,or_Hc_est,or_Hc_lo,or_Hc_hi,or_Hc_se) within 1e-8 relative on every row (max 0.00e+00)
  meta: subgroup_method dina | dina_select_statistic effect | dina_args list() | focus effMaxSG | nbhd 0.2 | rule neighborhood | ci field | scalec selected | pkg 0.3.5 | host pop-os | workers 20 | R 4.6.1 | effect_threshold 30
  [PASS] meta carries identifier dina with its arguments (effect, list()), focus effMaxSG, band 0.20, the rule, ci_method field, field_scale_complement selected, pkg 0.3.5, host pop-os; stem dina_effMaxSG_
  FACT declared: 20 of 20 | NO-DETECTION 0
  FACT sim_id 1: sg_def [{cd40 >= 400} & {cd80 >= 1040}] | n_sel 78 | n_harm 78 | nv_H_est 109.2217 | searched 8324 | proposed 2690 | admitted_n 1768 | n_family 2690 | tau_min 30.0064 | sens 0.132 ppv 0.308
  FACT fix record F5 after P2 (quoted): sg [{cd40 >= 400} & {cd80 >= 1040}] n 78 | searched 8324 | proposed 2690 | admitted_n 1768 | tau_hat min 30.0064
  [PASS] sim_id 1 selects the fix record's F5 after-P2 subgroup ({cd40 >= 400} & {cd80 >= 1040}, n 78): the campaign's DINA arguments equal F5's
  [PASS] E2 fields filled on every declared replicate (dina_searched_n 20/20, dina_proposed_n 20/20, dina_tau_min 20/20, admitted_n 20/20)
  STRUCTURAL E2 fields on non-detections all NA: TRUE
  [PASS] every proposed candidate at oriented tau-hat >= 30 on declared replicates (min over replicates 30.0001)
  [PASS] 1 <= admitted_n <= dina_proposed_n on every declared replicate
  FACT family sizes (declared): searched 8314-8768 | proposed min 308 median 2174 mean 3056.3 max 6474 | admitted min 247 median 1586 mean 2468.2 max 6287 | n_family (MR) min 308 median 2174 max 6474
  [PASS] n_family filled on every declared replicate with a gate (values 308-6474)
  [PASS] p_hat_H and p_hat_sum recorded (p_hat_H mean 0.077)
  STRUCTURAL n_cons_qual all-NA: TRUE | band_n all-NA: TRUE (DINA has no consistency screen)
  WARNINGS (verbatim, distinct, with per-row counts on 0 of 20 rows):
  [PASS] zero factor-comparison warnings (0 rows carry one)
  FACT timing: fit_mr_secs mean 63.1 median 50.1 max 123.2 | id_secs (DINA fit, incl. re-selection) mean 5.69 median 4.11 max 13.41 | fld_H_secs mean 34.8 | fld_Hc_secs mean 2.73
  FACT vs mdsgnb20 (sim 1-20): mean |Hhat| (n_harm) 104.9 vs 120.3 | sens 0.230 vs 0.309 | ppv 0.364 vs 0.423 | declared 20 vs 20
  [PASS] nine fld_Hc_*_s and nine fld_joint_s_* columns present
  [PASS] nine fld_Hc_*_s finite on all 20 filled replicates
  [PASS] nine fld_joint_s_* finite on all 20 filled replicates
  [PASS] fld_Hc_lo1s_s <= fld_Hc_up1s_s and fld_Hc_lo2s_s <= fld_Hc_hi2s_s
  [PASS] Bonferroni harm bound identical between joint and joint_s where draw counts agree (20 of 20 rows agree; max |diff| 0.00e+00)
  [PASS] field-s inverted around the same beta-tilde^c (max |diff| 7.11e-15)
  [PASS] meta field_scale_complement = selected
SMOKE dina md=40 n=500: PASS
```

**Floors as applied, and the proposed family, on every smoke replicate.** `logs_mddina/floor_check.R` (untracked, read-only) purls the template's `setup-knobs`, `build-dgm` and `machinery` chunks with `FS_MD_METHOD=dina`. For each smoke replicate it regenerates the data as `record_replicate()` does and re-identifies with the template's argument block, MR off and `details = TRUE`. It reads the "Harm floor" details line, `fs.est$admission$effect_floor`, and `grp.consistency$out_sg` (proposed candidates with their oriented `tau_hat` and `sel_effect`). Each result is checked against the smoke bundle.
```
forestsearch 0.3.5 Built R 4.6.1; ; 2026-09-17 04:47:31 UTC; unix | subgroup_method dina | dina_args list() | dina_select_statistic effect | effect.threshold 30 | adverse_outcome FALSE
sim 1 details (verbatim):
   sg_focus 'effMaxSG' resolves to canonical rule 'hrMaxSG' (aliases: effMaxSG).
   [forestsearch] DINA selection (subgroup_method = "dina")
     Family:              gaussian
     sg_focus:            hrMaxSG
     selection_rule:      neighborhood
     effect_neighborhood: 0.2
     Harm floor:          m_diff = 30.0000  on -tau-hat (harm-oriented; adverse_outcome = FALSE)
     n.min:               60
     DINA frontier candidates (per-covariate non-dominated):
     Candidates searched:  8324
     Candidates qualifying (>= floor, >= n.min): 2690
     SELECTED: {cd40 >=  400} & {cd80 >= 1040}  (n = 78, mean tau-hat = 74.2003)
sim 1 admission floor as applied: fs.est$admission$effect_floor = 30 (effect scale: harm-oriented MD, the effect estimator on -Y)
sim 1 tau_sign recorded in out_sg$call: tau_sign | .dina_tau_sign("continuous", FALSE) = -1
sim  1: Candidates searched:  8324 | Harm floor:          m_diff = 30.0000  on -tau-hat (harm-oriented; adverse_outcome = FALSE) | proposed 2690 | tau_hat min 30.0064 (< 30: 0) | sel_effect finite 2690, >= admission floor 30: 1768 | admitted_n 1768 | harm-oriented MD > 0 among proposed 2675 | {cd40 >= 400} & {cd80 >= 1040} | factor warnings 0 | matches bundle (sg, admitted, proposed): TRUE
sim  2: Candidates searched:  8768 | Harm floor:          m_diff = 30.0000  on -tau-hat (harm-oriented; adverse_outcome = FALSE) | proposed 1513 | tau_hat min 30.0001 (< 30: 0) | sel_effect finite 1513, >= admission floor 30: 1293 | admitted_n 1293 | harm-oriented MD > 0 among proposed 1504 | {wtkg >= 77.111999999999995} & {cd40 >= 359} | factor warnings 0 | matches bundle (sg, admitted, proposed): TRUE
sim  3: Candidates searched:  8330 | Harm floor:          m_diff = 30.0000  on -tau-hat (harm-oriented; adverse_outcome = FALSE) | proposed 5738 | tau_hat min 30.1641 (< 30: 0) | sel_effect finite 5738, >= admission floor 30: 4709 | admitted_n 4709 | harm-oriented MD > 0 among proposed 5722 | {preanti <= 284.59999999999991} & {wtkg <= 67.998000000000005} | factor warnings 0 | matches bundle (sg, admitted, proposed): TRUE
sim  4: Candidates searched:  8600 | Harm floor:          m_diff = 30.0000  on -tau-hat (harm-oriented; adverse_outcome = FALSE) | proposed 2840 | tau_hat min 30.0115 (< 30: 0) | sel_effect finite 2840, >= admission floor 30: 988 | admitted_n 988 | harm-oriented MD > 0 among proposed 2632 | {cd80 <= 1057.1999999999998} & {race >= 1} | factor warnings 0 | matches bundle (sg, admitted, proposed): TRUE
sim  5: Candidates searched:  8570 | Harm floor:          m_diff = 30.0000  on -tau-hat (harm-oriented; adverse_outcome = FALSE) | proposed 1739 | tau_hat min 30.0006 (< 30: 0) | sel_effect finite 1739, >= admission floor 30: 1157 | admitted_n 1157 | harm-oriented MD > 0 among proposed 1730 | {preanti <= 0} & {cd40 <= 300} | factor warnings 0 | matches bundle (sg, admitted, proposed): TRUE
sim  6: Candidates searched:  8582 | Harm floor:          m_diff = 30.0000  on -tau-hat (harm-oriented; adverse_outcome = FALSE) | proposed 1427 | tau_hat min 30.0002 (< 30: 0) | sel_effect finite 1427, >= admission floor 30: 1237 | admitted_n 1237 | harm-oriented MD > 0 among proposed 1422 | {age >= 42} & {cd40 <= 406} | factor warnings 0 | matches bundle (sg, admitted, proposed): TRUE
sim  7: Candidates searched:  8580 | Harm floor:          m_diff = 30.0000  on -tau-hat (harm-oriented; adverse_outcome = FALSE) | proposed 308 | tau_hat min 30.0116 (< 30: 0) | sel_effect finite 308, >= admission floor 30: 247 | admitted_n 247 | harm-oriented MD > 0 among proposed 307 | {preanti >= 833.00000000000023} & {cd80 >= 688.10000000000014} | factor warnings 0 | matches bundle (sg, admitted, proposed): TRUE
sim  8: Candidates searched:  8330 | Harm floor:          m_diff = 30.0000  on -tau-hat (harm-oriented; adverse_outcome = FALSE) | proposed 6337 | tau_hat min 30.0328 (< 30: 0) | sel_effect finite 6337, >= admission floor 30: 4432 | admitted_n 4432 | harm-oriented MD > 0 among proposed 6279 | {cd80 >= 1040.4000000000001} & {str2 >= 1} | factor warnings 0 | matches bundle (sg, admitted, proposed): TRUE
sim  9: Candidates searched:  8570 | Harm floor:          m_diff = 30.0000  on -tau-hat (harm-oriented; adverse_outcome = FALSE) | proposed 6004 | tau_hat min 30.0118 (< 30: 0) | sel_effect finite 6004, >= admission floor 30: 5507 | admitted_n 5507 | harm-oriented MD > 0 among proposed 6001 | {preanti <= 334.19999999999993} & {symptom <= 0} | factor warnings 0 | matches bundle (sg, admitted, proposed): TRUE
sim 10: Candidates searched:  8388 | Harm floor:          m_diff = 30.0000  on -tau-hat (harm-oriented; adverse_outcome = FALSE) | proposed 4964 | tau_hat min 30.0135 (< 30: 0) | sel_effect finite 4964, >= admission floor 30: 4467 | admitted_n 4467 | harm-oriented MD > 0 among proposed 4950 | {age >= 42} & {cd80 >= 681} | factor warnings 0 | matches bundle (sg, admitted, proposed): TRUE
sim 11: Candidates searched:  8602 | Harm floor:          m_diff = 30.0000  on -tau-hat (harm-oriented; adverse_outcome = FALSE) | proposed 6474 | tau_hat min 30.1229 (< 30: 0) | sel_effect finite 6474, >= admission floor 30: 6287 | admitted_n 6287 | harm-oriented MD > 0 among proposed 6467 | {cd80 <= 919} & {str2 >= 1} | factor warnings 0 | matches bundle (sg, admitted, proposed): TRUE
sim 12: Candidates searched:  8342 | Harm floor:          m_diff = 30.0000  on -tau-hat (harm-oriented; adverse_outcome = FALSE) | proposed 2047 | tau_hat min 30.0051 (< 30: 0) | sel_effect finite 2047, >= admission floor 30: 1632 | admitted_n 1632 | harm-oriented MD > 0 among proposed 2037 | {cd80 <= 605} | factor warnings 0 | matches bundle (sg, admitted, proposed): TRUE
sim 13: Candidates searched:  8562 | Harm floor:          m_diff = 30.0000  on -tau-hat (harm-oriented; adverse_outcome = FALSE) | proposed 742 | tau_hat min 30.0138 (< 30: 0) | sel_effect finite 742, >= admission floor 30: 353 | admitted_n 353 | harm-oriented MD > 0 among proposed 730 | {age <= 39} & {cd80 >= 1091.3000000000002} | factor warnings 0 | matches bundle (sg, admitted, proposed): TRUE
sim 14: Candidates searched:  8574 | Harm floor:          m_diff = 30.0000  on -tau-hat (harm-oriented; adverse_outcome = FALSE) | proposed 316 | tau_hat min 30.0691 (< 30: 0) | sel_effect finite 316, >= admission floor 30: 300 | admitted_n 300 | harm-oriented MD > 0 among proposed 316 | {preanti <= 30.600000000000023} & {homo >= 1} | factor warnings 0 | matches bundle (sg, admitted, proposed): TRUE
sim 15: Candidates searched:  8344 | Harm floor:          m_diff = 30.0000  on -tau-hat (harm-oriented; adverse_outcome = FALSE) | proposed 2302 | tau_hat min 30.0104 (< 30: 0) | sel_effect finite 2302, >= admission floor 30: 1636 | admitted_n 1636 | harm-oriented MD > 0 among proposed 2291 | {age <= 30} & {wtkg <= 77} | factor warnings 0 | matches bundle (sg, admitted, proposed): TRUE
sim 16: Candidates searched:  8624 | Harm floor:          m_diff = 30.0000  on -tau-hat (harm-oriented; adverse_outcome = FALSE) | proposed 1732 | tau_hat min 30.0068 (< 30: 0) | sel_effect finite 1732, >= admission floor 30: 893 | admitted_n 893 | harm-oriented MD > 0 among proposed 1704 | {age <= 30.700000000000017} & {homo >= 1} | factor warnings 0 | matches bundle (sg, admitted, proposed): TRUE
sim 17: Candidates searched:  8350 | Harm floor:          m_diff = 30.0000  on -tau-hat (harm-oriented; adverse_outcome = FALSE) | proposed 514 | tau_hat min 30.0022 (< 30: 0) | sel_effect finite 514, >= admission floor 30: 458 | admitted_n 458 | harm-oriented MD > 0 among proposed 512 | {wtkg <= 69.627600000000001} & {cd40 >= 339.5} | factor warnings 0 | matches bundle (sg, admitted, proposed): TRUE
sim 18: Candidates searched:  8384 | Harm floor:          m_diff = 30.0000  on -tau-hat (harm-oriented; adverse_outcome = FALSE) | proposed 6244 | tau_hat min 30.0070 (< 30: 0) | sel_effect finite 6244, >= admission floor 30: 5846 | admitted_n 5846 | harm-oriented MD > 0 among proposed 6231 | {wtkg <= 83.469920000000002} & {cd40 <= 251} | factor warnings 0 | matches bundle (sg, admitted, proposed): TRUE
sim 19: Candidates searched:  8632 | Harm floor:          m_diff = 30.0000  on -tau-hat (harm-oriented; adverse_outcome = FALSE) | proposed 5437 | tau_hat min 30.0327 (< 30: 0) | sel_effect finite 5437, >= admission floor 30: 4615 | admitted_n 4615 | harm-oriented MD > 0 among proposed 5409 | {age >= 34} & {preanti <= 52.5} | factor warnings 0 | matches bundle (sg, admitted, proposed): TRUE
sim 20: Candidates searched:  8314 | Harm floor:          m_diff = 30.0000  on -tau-hat (harm-oriented; adverse_outcome = FALSE) | proposed 1759 | tau_hat min 30.0154 (< 30: 0) | sel_effect finite 1759, >= admission floor 30: 1540 | admitted_n 1540 | harm-oriented MD > 0 among proposed 1753 | {preanti <= 335} & {cd80 >= 990} | factor warnings 0 | matches bundle (sg, admitted, proposed): TRUE
TOTAL: proposed candidates below oriented tau-hat 30: 0 | factor-comparison warnings: 0 | replicates matching the bundle: 20 of 20
```
(The purl step prints `object 'fb_in_tables' not found` from an inline expression outside the three chunks; it does not affect them.)

- **Proposal floor as applied:** `m_diff = 30.0000 on -tau-hat (harm-oriented; adverse_outcome = FALSE)` on all 20 replicates; `.dina_tau_sign("continuous", FALSE) = -1`, recorded in `out_sg$call` as the `tau_sign` argument. **Admission floor as applied:** `fs.est$admission$effect_floor = 30` on the harm-oriented MD. Both sit at 30 on the harm-oriented scale.
- **Every proposed candidate has oriented tau-hat ≥ 30:** 0 below 30 over all 20 replicates' proposed families; the per-replicate minimum ranges 30.0001–30.1641.
- **Zero factor-comparison warnings** (captured `warn_msg` NA on 20 of 20 rows; 0 in the re-identification). Factor covariates enter DINA as numeric (`.coerce_covariates_numeric()`, `helpers:1421–1423`); selections such as `{race >= 1}`, `{homo >= 1}` and `{str2 >= 1}` evaluate without warning.

**Facts (§1.6(c)):**
- Declared: 20 of 20.
- **sim_id 1:** `{cd40 >= 400} & {cd80 >= 1040}`, n 78; 8,324 searched, 2,690 proposed, 1,768 admitted, `tau_hat` min 30.0064. This is the fix record's F5 after P2 exactly (*"sg [{cd40 >= 400} & {cd80 >= 1040}] n 78"*; *"Candidates searched: 8324"*; *"Candidates qualifying (>= floor, >= n.min): 2690"*; *"admitted_n 1768"*; *"tau_hat (as ranked) min 30.0064"*). The campaign's DINA arguments equal F5's (`dina_args = list()`, `dina_select_statistic = "effect"`, `effMaxSG`, ε 0.20, `neighborhood`, the template's `adverse_outcome = FALSE`; F5 ran MR off, which does not enter the selection), so this agreement is a gate: **PASS.**
- `fit_mr_secs` mean 63.1 / median 50.1 / max 123.2 s. DINA fit time (`id_secs`, identification including the effect re-selection) mean 5.69 / median 4.11 / max 13.41 s. Field pass 34.8 s, complement field 2.73 s.
- Family sizes over declared replicates: searched 8,314–8,768; proposed min 308 / median 2,174 / mean 3,056.3 / max 6,474; admitted min 247 / median 1,586 / mean 2,468.2 / max 6,287; MR's `n_family` equals the proposed count on every replicate (308–6,474).
- Beside FS on the same 20 draws (descriptive): mean |Ĥ| 104.9 against 120.3, sensitivity 0.230 against 0.309, PPV 0.364 against 0.423.

*GATE §1.6:* (a), (b) and (c) hold.

## 1.7 Calibration

md40 n700, `FS_MD_METHOD=dina`, the campaign knobs, `FS_MD_WORKERS` = 16, 32, 63 with 3 × W replicates (tags `mddinacal16/32/63`), memory sampled every 5 s (`logs_mddina/calib.sh`). The summary, `logs_mddina/calib_summary.R`, applies the GRF task's §1.7 method as `mdsgnb20` and `mdgrf` did: fixed overhead = the 16-worker wall minus 3 rounds × the mean; loop cost = (wall − overhead) / replicates; per-cell cost scaled by `mdsgnb20`'s `fit_mr_secs` ratios; 90 s per combine render.
```
CALIB W=16 reps=48 rc=0 wall_s=450 peak_mb=27723 2026-09-17T08:27:55Z
CALIB W=32 reps=96 rc=0 wall_s=543 peak_mb=51535 2026-09-17T08:36:58Z
CALIB W=63 reps=189 rc=0 wall_s=795 peak_mb=101665 2026-09-17T08:50:13Z
CALIB DONE
total wall 1788 s
```
| W | replicates | render wall (s) | fit_mr_secs mean / median / p90 / max (s) | ratio to 16-worker mean | DINA fit (id_secs) mean / median / p90 (s) | fld_H_secs / fld_Hc_secs mean (s) | peak summed RSS (MB) | replicates per minute | declared | rows with warnings |
|---|---|---|---|---|---|---|---|---|---|---|
| 16 | 48 | 450 | 89.80 / 95.43 / 152.14 / 171.41 | 1.000 | 7.55 / 8.46 / 13.03 | 44.44 / 4.69 | 27723 | 6.40 | 48 | 0 |
| 32 | 96 | 543 | 96.66 / 95.58 / 176.84 / 202.26 | 1.076 | 8.03 / 7.76 / 15.13 | 46.86 / 5.01 | 51535 | 10.61 | 96 | 0 |
| 63 | 189 | 795 | 163.06 / 147.65 / 321.43 / 409.40 | 1.816 | 12.74 / 11.23 / 25.06 | 73.48 / 9.24 | 101665 | 14.26 | 189 | 0 |

Fixed render overhead (16-worker render, wall minus 3 rounds x mean): 180.6 s. Wall-based loop cost per replicate: W 16: 5.612 s; W 32: 3.775 s; W 63: 3.251 s.

Projection at W = 63 (loop cost 3.251 s per replicate at md40 n700, scaled per cell by mdsgnb20 fit_mr_secs ratios md40_n500 0.6626, md120_n500 0.8017, null_n500 0.6352, md40_n700 1.0000; overhead 181 s per batch render; combine render 90 s):

| cell | mdsgnb20 fit_mr_secs mean (s) | ratio | batch of 1,000 (s) | batch (min) | cell: 2 batches + combine (min) |
|---|---|---|---|---|---|
| md40_n500 | 56.24 | 0.6626 | 2335 | 38.9 | 79.3 |
| md120_n500 | 68.04 | 0.8017 | 2787 | 46.4 | 94.4 |
| null_n500 | 53.91 | 0.6352 | 2246 | 37.4 | 76.4 |
| md40_n700 | 84.88 | 1.0000 | 3432 | 57.2 | 115.9 |

**Projection: 21959 s = 366.0 min = 6.10 h.  Ceiling (1.5x): 32939 s = 549.0 min = 9.15 h.  Per-render timeout (2 x the longest projected batch, at least 20 min): 6864 s = 114.4 min.**
MDSG_WORKERS=63 MDSG_TIMEOUT=6864 MDSG_CEILING=32939
replicates per minute: W 16 6.40; W 32 10.61; W 63 14.26 | peak summed RSS (MB): W 16 27723; W 32 51535; W 63 101665 of 257519
ADVANCE_GO_CONDITION projection_h=6.100 under_8h=TRUE

- **W = 63.** It minimizes projected wall: 14.26 replicates per minute, against 10.61 at 32 and 6.40 at 16. **What limits it:** the template's worker cap of physical cores − 1 = 63 (`sim_fs_maxeffCons_mr_field_md_template.qmd:115–116`), i.e. cores. The scaling falloff is steeper than GRF's (per-replicate `fit_mr_secs` 1.82× the 16-worker mean at 63, p90 321 s), but throughput still rises to the cap. **Memory does not bind:** peak summed RSS 101,665 MB of 251 GB at W = 63, about 1.6 GB per worker (finding F7).
- **Projection for Stage 2 at W = 63: 21,959 s (366.0 min, 6.10 h)**, for four cells × (two 1,000-replicate batches + combine). **Ceiling** 1.5 × projection = **32,939 s (549.0 min, 9.15 h)**. **Per-render timeout** 2 × the longest projected batch (57.2 min) = **6,864 s (114.4 min)**, above the 20-min floor. Runner environment: `MDSG_WORKERS=63 MDSG_TIMEOUT=6864 MDSG_CEILING=32939`.
- Every calibration replicate declared (48, 96, 189); none carries a warning.
- Stage 1 compute: smoke renders 340 s, the floor check (20 fits, MR off), calibration 1,788 s — under the 2-h ceiling.

## 1.8 Gate 1 — PASS; the advance go applies

Every Stage 1 gate is green (§1.1, §1.2, §1.3, §1.6(a), (b), (c)) and the §1.7 projection (6.10 h) is under 8 hours. **The advance go, as the task states it (Dispositions):** "If every Stage 1 gate is green and the §1.7 projection for Stage 2 is under 8 hours, do not stop at Gate 1. Note the advance go and its condition in the Gate 1 record, run Stage 2 at the worker count chosen in §1.7 with the record's ceiling and per-render timeout, and run Stage 3 once Stage 2 is green." **Its condition is met.** Stage 2 launches after this record's commit: `scripts_mddina/run_mddina.sh` with `MDSG_WORKERS=63 MDSG_TIMEOUT=6864 MDSG_CEILING=32939`, knob `FS_MD_METHOD=dina`, `FS_MD_FOCUS=effMaxSG FS_MD_NBHD=0.20 FS_MD_FIELD_SCALEC=selected FS_MD_CI=field FS_MD_CAMPAIGN=mddina FS_MD_FB=none`, cells in `mdsgnb20`'s order (md40 n500, md120 n500, null n500, md40 n700). No cell, replicate, gate or knob is reduced.

## Findings

- **F1 (§1.1).** The task's `ps | grep -E '[e]xec/R|[R]script|[q]uarto'` returns an 11-day-old bash memory-sampling loop from another session, whose command text contains `quarto`. It runs no R, Rscript or quarto process (its only child is `sleep 3`). It samples `free -m` every 3 s into a scratchpad file of that session; it was not touched.
- **F2 (§1.5).** The catalog generator is `scripts_mdsgnb20/current_status_regen.R` (with `scripts_mdsgnb20/status_inventory.R` and `check_current_status.sh`), not `<dir>/current_status_regen.R`. Both generator files were edited. `current_status_regen.R`'s git-log spec had a literal `fs_` stem prefix, now read per campaign, without which the `mdgrf` / `mddina` renders would not be counted.
- **F3 (§1.3.5, E2).** S0.5's field list is complete in the template after `894da993`. The survival DINA recorder has no DINA-specific count (`admitted_n` is GRF-only in m1). The three `dina_*` fields and the DINA fill of `admitted_n` are authored, not transplanted: they carry the §1.6(c) quantities per replicate. They exist on declared replicates only, because DINA's selection object is dropped on a non-detection (`helpers:1546–1548`). On a DINA non-detection the bundle therefore cannot say whether the proposal or the admission came up empty.
- **F4 (§1.4).** `dinamr` passed `field_decompose = TRUE` and `field_recovery = TRUE`; this campaign takes the MD template's (neither passed), by the task's rule, as `mdsgnb20` and `mdgrf` did. The nine `fld_recov_*` columns stay NA.
- **F5 (§1.5).** The `meta` keys `dina_select_statistic` / `dina_args` are recorded on every path (as the GRF keys are), so FS and GRF bundles written by this template carry them as inert metadata.
- **F6 (§1.6(c)).** MR's family on the DINA path is the whole proposed family (`n_family` = `dina_proposed_n` on 20 of 20), up to 6,474 candidates. That is five times GRF's (1,185–1,372) and the reason DINA with MR costs more per replicate than GRF on this design (63.1 s against GRF's 32.5 s at 20 workers), whereas on survival it cost less.
- **F7 (§1.7).** DINA's per-replicate cost is heavy-tailed at W = 63 (median 147.7 s, p90 321.4 s, max 409.4 s), so the fixed-overhead estimate (180.6 s, from the 16-worker render's wall minus 3 × the mean) absorbs straggler time; the projection is correspondingly conservative. Memory reached 101.7 GB at W = 63 (the FS campaign's calibration reached 71.8 GB, GRF's 56.6 GB).
- **F8 (§1.4).** The committed `dinamr` bundle `meta` records no `dina_args`, `dina_select_statistic` or `selection_rule` key; the §1.4 values for them come from m1 and the runner.
- **F9 (Stage 3 preparation).** `summary_continuous_field_mddina.qmd`, the §3.1 transplant of `summary_continuous_field_mdgrf.qmd` (untracked until Stage 3), was dry-rendered against the md40 n500 DINA smoke bundle (`MDSG_SUMMARY_TAG=mddinasmoke MDSG_SUMMARY_GLOB=res_1_20 MDSG_SUMMARY_OUT=logs_mddina/dryrun`), exit 0: 281 DINA rows, 7 FS and 9 GRF comparator rows for the one cell, copied with their commits.

## Untracked outputs of this stage (for §3.5 deletion)

`mr_md_harm/fs_effMaxSG_mr_field_md40_knoise0_n500_nb20_mddinasmokefs_d5000/`, `mr_md_harm/grf_effMaxSG_mr_field_md40_knoise0_n500_nb20_mddinasmokegrf_d5000/`, `mr_md_harm/dina_effMaxSG_mr_field_md40_knoise0_n500_nb20_mddinasmoke_d5000/`, `mr_md_harm/dina_effMaxSG_mr_field_md40_knoise0_n700_nb20_mddinacal{16,32,63}_d5000/`, and in `logs_mddina/`: `smoke.sh`, `smoke_driver.log`, `{fs,grf,dina}.{log,peak_mb}`, `smoke_{fs,grf,dina}.html`, `smoke_identity_{fs,grf,dina}.txt`, `floor_check.R`, `floor_check.txt`, `calib.sh`, `calib_driver.log`, `calib_summary.R`, `calib_summary.txt`, `mddinacal{16,32,63}.{log,peak_mb,html}`, `dryrun/`.

## git log --oneline for this stage

```
189d4ec2 dev/tasks: TASK_md_dina_campaign_2026-09-17.md as received (campaign mddina, Stages 1-3)
0927669c MD template E1-E3 (TASK_md_dina_campaign_2026-09-17 §1.5): E1 no edit -- the DINA argument block (dina_args = list(), dina_select_statistic = "effect"
7e26ad42 actg175/continuous catalog generator (TASK_md_dina_campaign_2026-09-17 §1.5, carried fix): current_status_regen.R gains campaign rows mdgrf and mddina
d65f3b9b scripts_mddina (TASK_md_dina_campaign_2026-09-17 §1.6, transplants of scripts_mdgrf): mem_sampler.sh unchanged; run_mddina.sh with FS_MD_METHOD=dina (
<this record: the next commit>
```
