# REPORT — GRF on the ACTG175 continuous (MD) design: campaign `mdgrf`, Stage 1 (edits, smoke → Gate 1: **STOP at §1.6(c)**)

Date: 2026-09-16. Machine: `pop-os` (AMD Ryzen Threadripper PRO 5995WX, 64 physical cores, 251 GB; reference BLAS; R 4.6.1). Branch `feature/glm-extension`. Task: `dev/tasks/TASK_md_grf_2026-09-16.md` (committed as received, `e5ee1008`), on Larry's "proceed as you recommend" after the DINA/GRF Stage 0 (`S0` = `REPORT_md_dina_grf_stage0_2026-09-16.md`, `ee47ee61`). Installed forestsearch 0.3.5, `Built: R 4.6.1; ; 2026-09-16 05:57:14 UTC; unix`; no `R/` change, no install. Every render ran with `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1`. Line numbers: the MD template at `e5ee1008` (before the §1.5 edits), the survival template `m1` = `quarto/simulations/gbsg_020/sim_fs_maxeffCons_fb_mr_field_m1_template.qmd` at HEAD, `R/` at HEAD.

**Governing constraint.** GRF's candidate family is generated from a fitted surface, so the fixed-family condition does not hold. Every coverage figure this campaign would produce is coverage of the estimand conditional on the proposed family; comparisons with FS are descriptive, not a contest, and state the confound (identifier, family construction, detection set).

**Outcome.** Stage 1 gates §1.1, §1.2, §1.3, §1.6(a) and §1.6(b) pass. **§1.6(c) resolves as a stop:** the two factor-comparison warnings come from `.grf_evaluate_subgroup()` comparing a factor column of the analysis frame with `<=` / `>`, which yields NA membership; on this design that drops every enumerated candidate that uses one of the seven binary covariates at the effect re-selection (645 of 1,257 on sim_id 1, 237 of them above the DR floor), so the candidate evaluation is not correct and the task's rule ("if the comparison could mis-evaluate or drop a candidate on this design, stop") applies. The template edits and scripts are committed (green); the calibration (§1.7) was not run (it follows the gate); Stage 2 was not launched. Fixing the evaluation is an `R/` change (finding F2), out of scope; none is proposed here.

## 1.1 Provenance and first commit — GATE PASS

```
pop-os
feature/glm-extension
f3f188f9
f3f188f9 ACTG175 continuous intervals under effMaxSG, eps = 0.20: closeout (TASK_actg175_continuous_intervals_2026-09-16) -- ...
da25abcc ACTG175 continuous: intervals document under effMaxSG, eps = 0.20 (TASK_actg175_continuous_intervals_2026-09-16) -- ...
e1c3c7a6 Add TASK_actg175_continuous_intervals_2026-09-16 as received
[tracked modifications: none]
DINA/GRF Stage 0 in HEAD
[R / Rscript / quarto / deno processes, by process name: none]
```
- The task's `ps | grep -E '[e]xec/R|[R]script|[q]uarto'` pattern matches the text of the shell wrapper running it; by process name (`ps -eo comm`) no R, Rscript, quarto or deno process existed (finding F12).
- First commit: `e5ee1008 Add TASK_md_grf_2026-09-16 as received`.

## 1.2 Package check — GATE PASS

- `git diff --quiet 0071c17e..HEAD -- R/ DESCRIPTION NAMESPACE`: no change since the install.
- `packageDescription("forestsearch")$Built` = `R 4.6.1; ; 2026-09-16 05:57:14 UTC; unix` in the main session and in two `doFuture` multisession workers (pids 3786571, 3786572): agree.

## 1.3 Gate 0 from the Stage 0 record — PASS (quotations from `S0`)

- **Alignment (S0.2, §3 "Alignment").** *"GRF: its proposal surface is already harm-oriented under `adverse_outcome = FALSE`, `dmin.grf = 0` sits on it, and the admission floor is FS's 30 on the harm-oriented MD: aligned as the package stands, with `dmin.grf` on a different (DR-score) scale from FS's floor."* and the table row: *"`dmin.grf` (default 0; `main:1232`) on the DR harm-effect scale in RMST/outcome units: `elig <- cand[cand$effect >= dmin, ]` (`R/grf_subgroup_labels.R:358`) … re-ranking admission `effect_floor` 30 on harm-oriented MD (`helpers:1632–1635`, `admitted_n` `:1654`)"*. S0.2 does not state that no argument places GRF's floor on the harm-oriented scale: the admission floor 30 on the harm-oriented MD is `effect.threshold` (aliased into `hr.threshold`, `main:1352`; `effect_floor <- as.numeric(hr.threshold)`, `helpers:2345`), and `dmin.grf` is the frontier pre-filter's floor on the DR-score harm effect, which for a continuous outcome is in outcome (MD) units. The disposition `dmin.grf = 30` is applied as the knob `FS_MD_DMIN_GRF` (§1.5 E2). Not a stop; finding F1 records that the pre-filter is inert under the effect re-selection.
- **Forwarding (S0.4, §5.3).** *"forwarded from `mr_inference_args` with a fallback — `draws` (2000L), `include_complement` (TRUE), `confirm_rule` ("point"), `t_confirm` (passed through, NULL → near-null), `return_reselection`, `field_R_out`, `field_R_in`, `field_uniform`, `field_M_cap`, `field_complement`, `field_decompose`, `field_scale_complement`, `ij_residual`, `field_recovery` … `ci_method = .g(mr_inference_args$ci_method, "ij")` (`:171`) — the fallback is `"ij"`, not `fs_mr_inference()`'s `"field"`, so DINA and GRF get the field block only when `ci_method = "field"` is passed. Of the MD template's arguments, `ci_method`, `draws`, `include_complement`, `confirm_rule`, `field_uniform`, `field_complement`, `ij_residual`, `field_scale_complement` and `return_reselection` are all passed explicitly (template `:320–326`) and are forwarded; `t_confirm` only when set."* Every MR argument in the dispositions is forwarded; none is hard-coded on the GRF path to a value different from the FS branch's. Not a stop.
- **Re-selection (S0.4, §5.4).** *"The engines' effect re-selection (5.1) uses the same focus, the same band helper and the same harm-oriented MD scale: **aligned**"* (`.fs_mr_reselection_from_focus(sg_focus, engine = "effect")` → `effMaxSG`, `R/fs_mr_inference_methods.R:104`; `.compute_inclusion_band()` in `.fs_mr_select()`, `R/fs_mr_inference.R:135–149`). Not a stop.
- **Recorder fields (S0.5, §6.1).** *"recorder fields the survival DINA/GRF campaigns record and this template lacks: `n_family`, `n_cons_qual`, `band_n` (m1 `:988`), `admitted_n` (`:997`, fill `:1092–1094`), `p_hat_sum`, `p_hat_top1..3` (`:998–999`), the nine `fld_recov_*` (`:954–955`); the MD template records `p_hat_H`, `p_top1..3`, `p_lab1..3` (`:531–538`) but not `p_hat_sum`."*
- **Warnings (S0.3, §4).** *"`'<=' not meaningful for factors`, `'>' not meaningful for factors`"* — "ran with two warnings … the candidate-cut evaluation compares the analysis frame's factor covariates; a finding, F6".

## 1.4 The survival reference (read only)

Sources: **MD** = `sim_fs_maxeffCons_mr_field_md_template.qmd` at `e5ee1008`; **m1** at HEAD; **runner** = `quarto/simulations/gbsg_020/scripts_dinamr/grfmr.sh:19–21` (`KN=(FS_S7_METHOD=grf FS_S7_FOCUS=effMaxSG FS_S7_NBHD=0.20 FS_S7_FIELD_COMPLEMENT=TRUE FS_S7_FIELD_SCALEC=selected FS_S7_FIELD_DECOMP=TRUE FS_S7_FIELD_RECOV=TRUE FS_S7_IJ_RESIDUAL=two_term FS_S7_FB=none FS_S7_CAMPAIGN=grfmr FS_S7_WORKERS=12)`; its header `:12–15`: "`dmin.grf = 0.0`, `grf_selection = "frontier"` and `grf_select_statistic = "effect"` are TEMPLATE LITERALS (template lines 503-506), not `FS_S7_*` knobs"); **meta** = `quarto/simulations/gbsg_020/results/grf_effMaxSG_fb_mr_field_m1_h100_knoise0_n1000_nb20_grfmr_combined_1_2000.rds` (`subgroup_method grf`, `sg_focus effMaxSG`, `effect_neighborhood 0.2`, `ci_method field`, `field_complement TRUE`, `field_decompose TRUE`, `field_scale_complement selected`, `field_recovery TRUE`, `ij_residual two_term`, `mr_draws 5000`, `seed_base 8316951`, `consistency_method resample`, `stop_threshold NULL`, `forestsearch_version 0.3.5`; no `dmin.grf`, `grf_depth` or `selection_rule` key).

| Setting | MD template (`e5ee1008`) | Survival GRF campaign `grfmr` (m1 / runner / meta) | Campaign takes |
|---|---|---|---|
| identifier knob | literal `subgroup_method <- "consistency"` (`:134`) | `subgroup_method <- .env_chr("FS_S7_METHOD", "consistency")`; `stopifnot(subgroup_method %in% c("consistency", "dina", "grf"))` (m1 `:310–311`); runner `FS_S7_METHOD=grf` | E1: `FS_MD_METHOD`, default `consistency` |
| how each identifier is selected | `forestsearch(subgroup_method = )`: `"dina"` → `.forestsearch_dina_select()` (`main:2238`), `"grf"` → `.forestsearch_grf_select()` (`:2369`, `:2422`), else the consistency search; stem via `method_tag` (`:176`) and `fs_focus_tag()` (`:179`) | the same package dispatch; m1 `method_tag` `:395`, `focus_tag` `:416`; `method_args` switch `:1057–1069` | package |
| `grf_selection` | not passed | `"frontier"` (m1 `:533`) | `"frontier"` |
| `grf_select_statistic` | not passed | `"effect"` (m1 `:534`) | `"effect"` |
| `grf_depth` | not passed | `2L` (m1 `:535`) | `2L` |
| `dmin.grf` | not passed → `0.0` on GLM outcomes when missing (`main:1973–1975`) | `0.0` (m1 `:536`) | **30** (disposition), knob `FS_MD_DMIN_GRF` — the one departure from the survival arguments |
| `frac.tau` | not passed (0.8 default, `main:1233`) | not passed | omitted: survival-only (the GRF time horizon, `main:407–411`); `.build_grf_glm_args()` (`R/forestsearch_helpers.R`) does not forward it to `grf.subg.harm.glm()` |
| forest and honesty settings, tree count | package-level, both: `grf::causal_forest(X, Y_grf, W, W.hat = rep(0.5, n) [RCT], seed = seedit, tune.parameters = "none" [tune_grf = FALSE])` (`R/grf_subg_harm_glm.R:457–467`); grf 2.6.1 defaults `num.trees = 2000`, `honesty = TRUE` | same (no template control) | package |
| GRF seeding | `seedit = sd_i` with `sd_i <- seed_base + sim_id` (`:621`, `:662`); inside `set.seed(seedit)` and `causal_forest(seed = seedit)` (`glm:460`, `:463`) | `seedit = seed_base + sim_id` (m1 `:1044`) | same |
| frontier filter and band | package: `.grf_frontier_select()` `elig <- cand[cand$effect >= dmin, ]`, band `.compute_inclusion_band()`, largest in band (`R/grf_subgroup_labels.R:358–388`); effect re-selection with the admission floor `dmin_eff = effect_floor` (`helpers:1632–1641`), `admitted_n` (`:1651`, `:1654`) | same | package |
| `vi.grf.min` | `-0.2` (`:264`, `:653`); inert under `subgroup_method = "grf"` (S0 F4) | not passed (m1 has no `vi.grf` line) | MD's (inert) — finding F9 |
| `n.min` | `60L` (`:256`) | `NULL` (m1 `:573`; the survival GRF default applies) | MD's |
| focus / band / rule | `FS_MD_FOCUS` (`:139`) / `FS_MD_NBHD` (`:148`) / `"neighborhood"` (`:239`) | `FS_S7_FOCUS` / `FS_S7_NBHD` / `"neighborhood"` literal (m1 `:547`); runner `effMaxSG`, `0.20`; meta `effMaxSG`, `0.2` | `effMaxSG`, `0.20`, `neighborhood` |
| MR arguments on the GRF path | `mr_inference_args` (`:323–329`): `ci_method = "field"` (`FS_MD_CI`), `draws = 5000L`, `include_complement = TRUE`, `confirm_rule = "point"`, `field_uniform = FALSE`, `field_complement = TRUE`, `ij_residual = "two_term"`, `field_scale_complement = "selected"`, `return_reselection = TRUE`; `t_confirm` NULL (near-null, 0 for MD) | m1 `:683–691`: the same set plus `field_decompose = mr_field_decompose` and `field_recovery = mr_field_recovery`; runner sets both TRUE; meta `field_decompose TRUE`, `field_recovery TRUE` | MD's — findings F7, F8 (`field_decompose`, `field_recovery` not passed, as in `mdsgnb20`) |
| workers | `FS_MD_WORKERS` (cap physical − 1) | `FS_S7_WORKERS=12` (Mac) | smoke 20 |

## 1.5 Template edits — commit `894da993` (template only; `R/` untouched)

Applied bottom-up (Python, exact-line anchors) so the `e5ee1008` line numbers held; every R chunk parses; `git diff --stat` 124 insertions, 25 deletions (plus the `id_secs` lines).
- **E1 — identifier knob.** `:134` `subgroup_method <- "consistency"` replaced by m1 `:310–311` under the `FS_MD_` prefix: `subgroup_method <- .env_chr("FS_MD_METHOD", "consistency")` + `stopifnot(subgroup_method %in% c("consistency", "dina", "grf"))`. `method_tag` (`:176`) and `focus_tag` (`:179`) already key on it, so the stem becomes `grf_effMaxSG_mr_field_md40_knoise0_n500_nb20_<tag>` (observed). `use_dina` / `use_grf` are consistency-only `method_args` (m1 `:1061–1062`), so under `grf` they are not passed, as in m1.
- **E2 — GRF argument block.** After `:264`, m1 `:533–538` verbatim (`grf_selection <- "frontier"`, `grf_select_statistic <- "effect"`, `grf_depth <- 2L`, `dina_args <- list()`, `dina_select_statistic <- "effect"`) with `dmin.grf <- .env_num("FS_MD_DMIN_GRF", 30)` in place of m1's literal `0.0`, guarded `stopifnot(is.finite(dmin.grf), dmin.grf >= 0)`. The call (`:646–669`) is rebuilt as m1 `:1039–1073`: `base_args` (every former argument except the consistency-only set, plus `subgroup_method = subgroup_method`) + `method_args <- switch(subgroup_method, consistency = list(consistency_method, use_lasso, use_grf, use_twostage, use_dina, conf.cont_jcuts, fs.splits, maxk, d0.min, d1.min), dina = …, grf = list(grf_selection, grf_depth, dmin.grf, grf_select_statistic))` + `do.call(forestsearch, c(base_args, method_args))`. GRF seeding stays `seedit = sd_i` (m1 `:1044`). Warnings are captured with `withCallingHandlers()` into the new `warn_msg` column (distinct messages with counts) instead of `suppressWarnings()` (finding F5).
- **E3 — recorder.** After `:592`: `n_family`, `n_cons_qual`, `band_n`, `admitted_n`, `p_hat_sum` (m1 `:988–999`), the nine `fld_recov_*` (m1 `:952–955`), `warn_msg`, `id_secs` (= `60 * fs.est$minutes_all`, the engine's own clock before the MR gate; finding F6). Fills: `admitted_n` before the no-detection return (m1 `:1087–1094`); `n_cons_qual` / `band_n` from `grp.consistency$out_sg$result` after `DETECTED` (m1 `:1164–1171`); `n_family <- g$n_family` (m1 `:1148`); `p_hat_sum <- sum(ph)` in the re-selection block (m1 `:1151`); `fld_recov_*` from `f$recovery` (m1 `:1279–1291`; structurally NA here, finding F7). `p_hat_top1..3` are not duplicated: the MD template's `p_top1..3` / `p_lab1..3` (`:590–592`) carry the same quantities.
- **E4 — `meta` and poolability keys.** `meta` gains `dmin_grf`, `grf_selection`, `grf_depth`, `grf_select_statistic` after `subgroup_method` (`:979`); the combine-mode key vector (`:1028–1032`) gains the same four (`subgroup_method` was already a key).
- **Default-path check.** The argument set passed on the default path is the former literal call's, plus the explicit `subgroup_method = "consistency"` (the package default); the only removed lines are the literal call and the literal `subgroup_method`. The call form (`do.call`) and the warning capture are not "additions guarded by the knob" (finding F4); §1.6(a) verifies the rendered default path is unchanged to the last bit.

## 1.6 Scripts — commit `f0b9c844`, under `scripts_mdgrf/` (transplants of `scripts_mdsgnb20/`)

`mem_sampler.sh` (unchanged); `run_mdgrf.sh` (the FS runner with `FS_MD_METHOD=grf FS_MD_DMIN_GRF=30` added to its knob line, the `grf_effMaxSG_` stem, `TAG=mdgrf`, `REPORT_md_grf_gate2_2026-09-16.md`, the conditional-family sentence in the Gate 2 report header); `gate2.R` (the FS Gate 2 pointed at the GRF bundles: comparator `mdsgnb20` in both directions, `meta` with the four GRF keys, `admitted_n` recorded prominently and required ≥ 1 on declared replicates, `n_family` stated as MR's kept family, `n_cons_qual` / `band_n` reported as structural NA, `p_hat_sum`, MR failures ≤ 40, `warn_msg` summary, `id_secs` timing); `smoke_identity.R` (mode `identity` for §1.6(a) against `mdsgnb20` sim_id 1–20, every column except `*_secs`, new columns reported apart; mode `grf` for §1.6(b)). The smoke driver `logs_mdgrf/smoke.sh` is untracked.

## 1.6 Smoke — (a) PASS, (b) PASS, (c) **STOP**

Renders (`logs_mdgrf/smoke.sh`; md40 n500, sim_id 1–20, 20 workers, `FS_MD_FOCUS=effMaxSG FS_MD_NBHD=0.20 FS_MD_FIELD_SCALEC=selected FS_MD_CI=field FS_MD_FB=none`; template `894da993`):

| render | knob | tag | exit | wall (s) | peak summed RSS (MB) |
|---|---|---|---|---|---|
| (a) FS regression | default | `mdgrfsmokefs` | 0 | 91 | 18,358 |
| (b) GRF | `FS_MD_METHOD=grf FS_MD_DMIN_GRF=30` | `mdgrfsmoke` | 0 | 71 | 12,985 |

**(a) FS regression** (`scripts_mdgrf/smoke_identity.R 40 500 mdgrfsmokefs 20 identity`, `logs_mdgrf/smoke_identity_fs.txt`): 137 paired `mdsgnb20` recorder columns (126 numeric + 11 character; `*_secs`, messages and FB columns apart) on sim_id 1–20: **20 of 20 rows identical, max relative difference 0.00e+00, zero enumerated rows, zero selection flips**; `mr_msg` / `err_msg` / `fb_err` identical on every row; the 16 new columns reported apart (on the FS path `n_cons_qual`, `band_n`, `n_family`, `p_hat_sum`, `id_secs` finite on 20/20; `admitted_n` NA on 20/20; `warn_msg` NA on 20/20 — the FS path raised no warning); stem and `meta` carry the consistency engine; the field-s checks pass (nine + nine finite on 20/20; `lo1s_s ≤ up1s_s`; Bonferroni harm bound identical between `joint` and `joint_s` on 20/20; `est2_s + lam_mean_s` vs `est2 + lam_mean` max |diff| 1.78e-15). **PASS.**

**(b) GRF** (`… mdgrfsmoke 20 grf`, `logs_mdgrf/smoke_identity_grf.txt`):
- ran without error; 20 of 20 `DETECTED`, no `CONFIG-ERROR`;
- `n_true` identical and the eight oracle columns within 1e-8 (max 0.00e+00) of `mdsgnb20` on all 20 sim_ids — **same draws**;
- `meta`: `subgroup_method grf`, `dmin_grf 30`, `grf_selection frontier`, `grf_depth 2`, `grf_select_statistic effect`, `sg_focus effMaxSG`, `effect_neighborhood 0.2`, `selection_rule neighborhood`, `ci_method field`, `field_scale_complement selected`, `pkg_version 0.3.5`, `hostname pop-os`, `n_workers 20`, `r_version 4.6.1`; stem `grf_effMaxSG_mr_field_md40_knoise0_n500_nb20_mdgrfsmoke`;
- the floor as applied: `meta$dmin_grf = 30` on the DR-score harm effect (the frontier pre-filter); `meta$effect_threshold = 30` on the harm-oriented MD (the admission floor; `effect_floor`, `helpers:2345`);
- E3 fields: `admitted_n` recorded on 20/20 (137–568; sim_id 1: 234), `n_family` 604–680 (sim_id 1: 612), `p_hat_H` and `p_hat_sum` on 20/20 (p̂(Ĥ) mean 0.101); `n_cons_qual` and `band_n` all-NA (structural: no consistency screen on GRF); the nine `fld_Hc_*_s` and nine `fld_joint_s_*` finite on 20/20 filled replicates; `fld_Hc_lo1s_s ≤ fld_Hc_up1s_s`; Bonferroni harm bound identical between `joint` and `joint_s` (20/20 draw counts agree, max |diff| 0); field-s inverted around the same β̃ᶜ (max |diff| 1.78e-15);
- **facts:** declared 20 of 20; **sim_id 1 selects `{preanti <= 792.8} & {cd40 > 364}`, n = 156, naive oriented MD 72.0396, `admitted_n` 234, sensitivity 0.346, PPV 0.404 — S0.3's replicate exactly** (`{preanti <= 792.8} & {cd40 > 364}`, n 156 / 344, MD 72.04, admitted 234, sens 0.346, PPV 0.404); `fit_mr_secs` mean 22.8 / median 22.7 / max 24.3 s; GRF fit time (`id_secs`, the identifier incl. its effect re-selection) mean 3.14 / median 3.13 / max 3.37 s; field pass 13.8 s, complement field 0.74 s; against FS on the same 20 draws: mean |Ĥ| 123.5 vs 120.3, |Ĥ| larger on 10 and smaller on 10, sensitivity 0.311 vs 0.309, PPV 0.435 vs 0.423 (descriptive; different identifiers and families);
- **warnings, verbatim (captured per replicate; on 20 of 20 rows):** `‘<=’ not meaningful for factors` (986–1166 times per replicate) and `‘>’ not meaningful for factors` (305–412 times per replicate).
**PASS** on its checks.

**(c) The warnings — STOP.**
- *Source.* `R/grf_subgroup_labels.R:224–228`, in `.grf_evaluate_subgroup(def, df)`: `x <- df[[v]]; member <- switch(op, "<=" = x <= val, ">" = x > val, …)` — the split variable's column of the **data frame** is compared with the numeric cut. The forest's covariate matrix is coerced first (`R/grf_subg_harm_glm.R:884–895`: a factor with all-numeric levels becomes `as.numeric(as.character(x))`), so candidates on the binary covariates are enumerated on `{0, 1}` (`.grf_dr_candidates()`, `labels:255–277`) — but `.grf_evaluate_subgroup()` applies no such coercion, and `factor <= 0` in R is NA with the warning `'<=' not meaningful for factors`.
- *What is compared.* On this design the analysis frame carries `hemo`, `homo`, `drugs`, `race`, `gender`, `symptom`, `str2` as factors (`sapply(df[confounders_analysis], class)` on sim_id 1: the six continuous covariates integer/numeric, the seven binary ones `factor`; the DGM keeps the frame's factor columns, `R/generate_glm_dgm.R:278–280`, `:404`; `simulate_from_glm_dgm()` samples its rows, `R/simulate_from_dgm.R:231–232`). Every candidate whose `v1` or `v2` is one of them is evaluated to NA membership.
- *Consequence, measured* (read-only probe on sim_id 1, md40 n500, the template's `setup-knobs` / `build-dgm` / `machinery` chunks purled into the scratchpad outside the repo, `forestsearch(subgroup_method = "grf", …, mr_inference = FALSE)` with the campaign's arguments; `n_true` 182 as in the bundles): 1,257 candidates enumerated; **612 use continuous covariates only and all 612 are scorable (`sel_effect` finite); 645 use a binary covariate and none is scorable** — at the effect re-selection, `mem <- which(.grf_evaluate_subgroup(sgd_i, df) == 0L)` is empty for them and `if (length(mem) < 6L) next` skips them (`R/forestsearch_helpers.R:1612–1614`), so they never enter `cand_hr`, the admission or the band; 237 of the 645 clear the DR floor 30 (604 at floor 0). `admitted_n` = 234 and the selection are identical at `dmin.grf` 30 and 0. Direct check: `{gender <= 0}` on the factor column → NA on 500 of 500 rows (with the warning); on the same column coerced to numeric → 74 members. In the 20 smoke replicates no selected rule uses a binary covariate; MR's `n_family` (604–680) coincides with the scorable continuous-only count (612 on sim_id 1), so MR's family excludes them too (finding F3).
- *Decision.* The comparison drops every candidate on the seven binary covariates on this design, `str2` (the binary proxy for `preanti`, a true-region covariate) among them. The task's rule: "If the comparison could mis-evaluate or drop a candidate on this design, stop." **Stop.** The correction — evaluating membership on the same coerced scale the X builder uses — is an `R/` change (finding F2); none is proposed here.

*GATE §1.6:* (a) holds, (b) holds, **(c) is resolved as a stop.**

## 1.7 Calibration — not run

It follows the §1.6 gate; with the gate stopped, no calibration render was authorized to any purpose. (For reference, the survival `grfmr` cost 15–17 s per replicate at 12 Mac workers, and here the GRF smoke ran 22.8 s per replicate at 20 workers including the field pass.)

## 1.8 Gate 1 — STOP; the advance go does not apply

The advance-go condition ("if every Stage 1 gate is green and the §1.7 projection for Stage 2 is under 8 hours") is not met: §1.6(c) is not green. Stage 2 was not launched; `run_mdgrf.sh` was not executed; `LOG_mdgrf_progress.txt` and `HALT_mdgrf.md` do not exist; the catalog (`status_curated.md` / `current_status.md`) is unchanged.

## Findings

- **F1 (§1.3, `dmin.grf`).** `dmin.grf` is the floor of the frontier pre-filter on the DR-score harm effect (`labels:358`), not a floor on the harm-oriented MD; the MD admission floor 30 is `effect.threshold`'s. Under `grf_select_statistic = "effect"` the pre-filter only picks the native DR winner, which the effect re-selection over the full scorable set overrides: on sim_id 1, `dmin.grf` 30 and 0 give the same selection and `admitted_n` 234. The disposition's `dmin.grf = 30` is applied and recorded, and is inert on this path whenever any candidate is scorable.
- **F2 (§1.6(c), the stop).** `.grf_evaluate_subgroup()` evaluates cuts on the raw data-frame column without the factor coercion the forest's X builder applies (`glm:884–895`); on this design's seven binary factor covariates every candidate's membership is NA, and the effect re-selection drops them (645 of 1,257 on sim_id 1; 237 above the DR floor). This is the source of S0's F6 warnings. An `R/` change would be needed (coerce the frame's factor columns as the X builder does before comparing); none is proposed.
- **F3 (§1.6(b)).** MR's `n_family` on the GRF path (604–680) equals the count of scorable candidates (612 on sim_id 1), not the enumerated pool (1,257) that `gate2G.R` names `n_family` on survival: the factor-covariate candidates are absent from MR's family as well. Observed, not traced in source.
- **F4 (§1.5).** The call is now `do.call(forestsearch, c(base_args, method_args))` (m1's form) rather than the literal call; the argument set on the default path is unchanged plus the explicit `subgroup_method = "consistency"`; §1.6(a) shows bit-identical output (max relative difference 0 on 126 numeric columns, 20/20 rows).
- **F5 (§1.5).** Warnings inside `forestsearch()` are captured into `warn_msg` (distinct messages with counts) instead of being suppressed; on the FS path the 20 smoke replicates raised none.
- **F6 (§1.5).** `id_secs` (60 × `minutes_all`, the identifier's clock before the MR gate) was added to report the GRF fit time the task asks for; it is not among S0.5's fields.
- **F7 (§1.5).** The nine `fld_recov_*` columns are structurally NA: `field_recovery` is not passed by the MD template (nor was it by `mdsgnb20`); they are carried for recorder parity with the survival DINA/GRF campaigns.
- **F8 (§1.4).** The survival `grfmr` campaign passed `field_decompose = TRUE` and `field_recovery = TRUE`; this campaign takes the MD template's (neither passed), by the task's rule (as `mdsgnb20`'s F2).
- **F9 (§1.4).** `vi.grf.min = -0.2` is passed on the GRF path and inert there (S0 F4); m1 does not pass it.
- **F10 (§1.6(b), cost).** GRF with MR costs 22.8 s per replicate at 20 workers on this design (identifier 3.1 s, field 13.8 s) against FS's 38.7 s at 20 workers under the same rule (`mdsgnb20` Stage 1 §1.6); the smoke reached 12,985 MB peak at 20 workers.
- **F11 (§1.6(b), fact).** On the same 20 draws GRF's |Ĥ| is larger than FS's on 10 replicates and smaller on 10 (means 123.5 vs 120.3), with sensitivity 0.311 vs 0.309 and PPV 0.435 vs 0.423 — descriptive only, and conditional on families that exclude the binary covariates (F2).
- **F12 (§1.1).** The task's process check matches its own shell wrapper's text; by process name nothing was running.

## Untracked outputs of this stage (listed; not deleted, since §3.5 was not reached)

`mr_md_harm/fs_effMaxSG_mr_field_md40_knoise0_n500_nb20_mdgrfsmokefs_d5000/` (1 batch bundle), `mr_md_harm/grf_effMaxSG_mr_field_md40_knoise0_n500_nb20_mdgrfsmoke_d5000/` (1 batch bundle), `logs_mdgrf/` (`smoke.sh`, `smoke_driver.log`, `fs.log`, `grf.log`, `fs.peak_mb`, `grf.peak_mb`, `smoke_fs.html`, `smoke_grf.html`, `smoke_identity_fs.txt`, `smoke_identity_grf.txt`). The probe lived in the session scratchpad outside the repo. `logs_mdsgnb20/` (untracked, pre-existing from the FS campaign) was not touched.

## git log --oneline for this stage

```
f0b9c844 scripts_mdgrf (TASK_md_grf_2026-09-16 §1.6, transplants of scripts_mdsgnb20): ...
894da993 MD template E1-E4 (TASK_md_grf_2026-09-16 §1.5, transplanted from the survival m1 template): ...
e5ee1008 Add TASK_md_grf_2026-09-16 as received
<this record: the next commit>
```
