# REPORT — Continuous/MD field port, Stage 0 (discovery)

Date: 2026-09-07. Machine: Mac Studio (Apple M4 Max). Branch: `feature/glm-extension-mac` from `22c5f854`. Package installed from this tree: forestsearch 0.3.5 (= DESCRIPTION). Task: `dev/tasks/TASK_continuous_field_mac_2026-09-07.md` (M-1–M-4 at defaults). No compute beyond one timing replicate (0d); no R/ edits in this stage.

## 0a. The continuous twin and its committed bundles

**Template.** There is no file named "template" on the continuous side. The twin is the batch document `quarto/simulations/actg175/continuous/sim_fs_maxeffCons_mr_md40_knoise0_n500_batch_1_1000.qmd` (1,587 lines; chunks `setup-knobs`:66, `build-dgm`:281, `machinery`:448, `run-batch`:799, then the summary layer). The n = 700 document differs from it in exactly one line (`n_sample <- 700L`, line 133; verified by `diff`); the md120 document adds the Stage-2 breadth block (direct `generate_glm_dgm()` build at the locked `k_inter`, optional `_c1star` stem under `FS_STAGE2_RUN=direct`). The null cell is the same document with `null_cell <- TRUE` (line 139; stem token `_mdnull_`, line 158).

**Stem construction** (lines 151–160, quoted):
```r
rds_stem <- sprintf("%s_%s_mr_md%02d_knoise%d_n%d",
                    method_tag, focus_tag, abs(target_md_harm),
                    k_random_noise, n_sample)
quickrun <- FALSE
if (isTRUE(null_cell)) rds_stem <- sub("_md[0-9]+_", "_mdnull_", rds_stem)
if (isTRUE(quickrun)) rds_stem <- paste0(rds_stem, "_quickrun")
results_dir <- file.path("mr_md_harm", sprintf("%s_s%d_d%d", rds_stem, n_sims, mr_draws))
```
No campaign tag, no env knobs, no save guard (`saveRDS(...)` at line 876 writes into `results_dir` unconditionally).

**Committed cells (M-1 candidates).** All under `quarto/simulations/actg175/continuous/mr_md_harm/`, `n_sims = 1000`, `sim_id` 1–1000, `seed_base = 8316951` (replicate seed `seed_base + sim_id`, `RNGkind("L'Ecuyer-CMRG")` at line 539), `mr_draws = 5000`, `sg_focus = "maxeffCons"`, `consistency_method = "resample"`, `effect.threshold = 30`, `consistency.threshold = 10`, `pconsistency = 0.90`, `maxk = 2`, `n.min = 60`, `conf.cont_jcuts = list(age = 10, preanti = 10)`, `include_str2 = TRUE`, `nb_boots = NULL` (no FB):

| cell | bundle | rows | detected | pkg | mean fit+MR s (Linux, 115 workers) |
|---|---|---|---|---|---|
| md40, n = 500 | `fs_maxeffCons_mr_md40_knoise0_n500_s1000_d5000/…_res_1_1000.rds` | 1000 | 1000 | 0.2.2 | 26.3 |
| md40, n = 700 | `fs_maxeffCons_mr_md40_knoise0_n700_s1000_d5000/…_res_1_1000.rds` | 1000 | 999 | 0.2.2 | 63.9 |
| md120, n = 500 | `fs_maxeffCons_mr_md120_knoise0_n500_s1000_d5000/…_res_1_1000.rds` | 1000 | 1000 | 0.3.1 | 30.0 |
| md120, n = 500, c1* | `…_md120_knoise0_n500_c1star_s1000_d5000/…_c1star_res_1_1000.rds` (effect.threshold = 135.741) | 1000 | 786 | 0.3.1 | 16.5 |
| null, n = 500 | `fs_maxeffCons_mr_mdnull_knoise0_n500_s1000_d5000/…_res_1_1000.rds` | 1000 | 998 | 0.2.2 | 25.5 |

Also committed: `…_md40_knoise0_n500_s100_d5000` (100 reps, 0.2.0) and the FB bundle `fb_mr_md_harm/fs_maxeffCons_fb_mr_md40_knoise0_n500_s100_d5000/…_res_1_100.rds` (100 reps, `nb_boots = 300`, pkg 0.2.0; FB ≈ 5,040 s/replicate) plus three `_quickrun` FB bundles (s20/s100/s1000). For M-3 the only joinable FB is md40 n = 500, `sim_id` 1–100.

**M-1 proposal for confirmation.** The four design cells md40 n500, md40 n700, md120 n500, null n500 (the c1* variant is a breadth-forecast scoring run at a non-standard threshold and is left out).

**DGM and harm orientation.** ACTG175 arms 1/3 (`speff2trial::ACTG175`), treatment coded ddI = 1 (`treat <- 1L - (arms == 1)`), outcome `cd4_change = cd420 − cd40`, `adverse_outcome = FALSE` (line 216: higher CD4 change is better, harm = negative MD). True region `z1 = age > 34 & z2 = preanti <= 744.5`, prevalence 0.3446 on the 5,000-row super-population. `calibrate_glm_interaction()` to `target_md_harm = −40` (structural `effect_Q = −40.0000000000`, `effect_Qc = −26.2552358760`, ITT −30.99, `beta_inter = −13.7`); md120 is built directly at the locked `k_inter = −93.7447641240`. Two scales coexist in the bundle: the MR gate's columns (`nv_/mr_*`), oracle and FB are **oriented** (positive = harm; the gate's working scale is `−cd4_change` via `consistency_resample.R:255–257`); `betaHhat_H/Hc` and `truth$effect_*` are **raw** (negative = harm). The summary layer bridges with `.orient = −1` (line 999).

**Truth attachment.** `fs_build_eval_frame(dgm, outcome_type = "continuous")` (exact finite-population frame; zero Monte Carlo error) and the post-loop `fs_attach_betaHhat(results, eval_df, focus = "harm", outcome_type = "continuous", effect_measure = "MD")` (line 849) giving `betaHhat_H`, `betaHhat_Hc`, `betaHhat_status`, `nH_eval`, `nHc_eval`. Reference set per block (the twin's θ†/θ‡ analogues, lines 1040–1050): `betaHhat` (per-replicate, exact), `struct` (`effect_Q`/`effect_Qc`), `marg` (fitted `lm` on one realized super-population draw, `marg_H = −35.83 (SE 6.24)`, `marg_Hc = −22.54 (SE 4.57)`; collapsibility z = +0.67 / +0.81), and the per-replicate oracle.

**Gate call** (lines 561–591, quoted in part):
```r
fs.est <- suppressWarnings(forestsearch(
  df.analysis = df, confounders.name = confs,
  outcome.name = outcome_name, treat.name = treat_name, id.name = id_name,
  outcome_type = "continuous", effect_measure = "MD",
  effect.threshold = md_threshold, consistency.threshold = md_consistency,
  pconsistency.threshold = pconsistency, fs.splits = fs_splits,
  n.min = n_min, d0.min = d0_min, d1.min = d1_min, maxk = maxk,
  vi.grf.min = vi_grf_min, sg_focus = sg_focus, selection_rule = selection_rule,
  effect_neighborhood = effect_neighborhood, stop_threshold = stop_threshold,
  consistency_method = consistency_method, conf.cont_jcuts = fs_conf.cont_jcuts,
  use_lasso = use_lasso, use_dina = use_dina, use_grf = use_grf,
  use_twostage = use_twostage, is.RCT = is_rct, adverse_outcome = adverse_outcome,
  details = FALSE, quiet = FALSE, seedit = sd_i, parallel_args = inner_parallel,
  mr_inference = TRUE,
  mr_inference_args = list(ci_method = "ij", draws = mr_draws, include_complement = TRUE)))
```

**Bundle columns (59; the identity anchors).** `sim_id status detected mr_ok err_msg mr_msg n_sel n_harm n_true sg_def covs betaHhat_H betaHhat_Hc fb_secs fit_mr_secs fb_err fb_src1 fb_src2 fb_nres`; oracle `or_H_{est,lo,hi,se}`, `or_Hc_*`; naive `nv_H_*`, `nv_Hc_*`; FB `fb_H_*`, `fb_Hc_*` (all NA in the MR-only bundles); MR (IJ) `mr_H_{est,lo,hi,se_ij}`, `mr_Hc_*`; `ij_source sens spec ppv npv betaHhat_status nH_eval nHc_eval`. Coverage indicators are not stored per replicate; they are computed in the summary layer (`.cover()`) and in the bundle's `oc$estimation` (`fs_mr_oc_summary()`: `cov_beta/cov_oracle/cov_struct/cov_marg` per block × estimator). No `fld_*`, no `p_hat`, no `mr_harm_flag`, no `label` column. Bundle elements: `results truth scale oc meta` (`meta` has no `ci_method`, `campaign_tag`, `field_complement`, `pkg_version` is `0.2.2`/`0.3.1`).

## 0b. The gate's field path on the continuous specification

- **Dispatch.** The consistency engine's MR call for every outcome type is the single site `R/forestsearch_main.R:3388–3421` (`.mr_glm_ok` / `.mr_cox_ok`), which forwards `ci_method`, `field_uniform`, `return_reselection`, `field_M_cap`, `field_complement`, `ij_residual`, `seed`; the GLM `gspec` carries `adverse_outcome = adverse_outcome`. `field_R_out`/`field_R_in` are not forwarded (package defaults 1000/500), exactly as on the survival template. No change needed.
- **Scale.** `.consistency_glm_pieces()` returns `log_scale = effect_measure %in% c("OR","RR","IRR")` (`R/consistency_resample.R:323`), so MD has `log_scale = FALSE`; `fs_mr_inference()` sets `to_eff <- function(x) if (log_scale) exp(x) else x` (`R/fs_mr_inference.R:543`) — the identity. Every field quantity is `to_eff()`'d from the working scale: `est2`, `lower_1s = to_eff(beta_deb − q95)` (:851–852), the two-sided quantile pair, the SE-type pair; the complement's `upper_1s = to_eff(bdc − q05)` (:1059); the joint pair (:1094–1118). `lambda_mean`, `lambda_sd = se_field` and the seven quantiles are working-scale (= MD-scale) summaries. Nothing in the field block (:789–925), `.fs_mr_field_complement()` (:989) or `.fs_mr_field_joint()` (:1094) assumes positivity or a log link.
- **Harm orientation.** The working scale is oriented (`−cd4_change`, positive = harm) because `adverse_outcome = FALSE` flips the outcome at `consistency_resample.R:255–257`. Hence `fld_H_lo1s` is the harm-direction lower bound on the oriented MD scale ("harm ≥ L"), and `fld_Hc_up1s` the complement's oriented upper bound ("harm ≤ U", i.e. benefit ≥ −U on the raw scale): the same one-sided-upper convention as the survival complement, applied on the oriented scale that every `*_est/lo/hi` column already uses. Targets for coverage are the oriented `betaHhat_H/Hc` (`.orient * betaHhat_*`).
- **`return_reselection = TRUE`** (`:957–961`) attaches `reselection$p_hat` (named by family label) and `winner`; `p_hat[selected_index]` is p̂(Ĥ). The recording convention exists upstream in `quarto/GuoHe/mr_vs_guohe_sim.R:210–217` (`p_hat_H`, `p_hat_top3_labels`, `p_hat_top3`); the survival template at 22c5f854 does not record it (the concurrent Linux nb20 task adds it there).
- **`fs_sim_bias_coverage()` scale assumptions** (`R/fs_bias_coverage.R`): `lt <- log(tgt)` (:95), `le <- log(cc$e)` (:135), `bias_log <- mean(le − lt)` (:137), and the normal-based one-sided bound `exp(log(e) ± z1·se)` guarded by `e > 0` (:128–129). On the oriented MD scale the complement's estimates can be near or below zero, so the log path would drop rows and mis-state bias. **Required change:** `scale = c("log", "identity")`, add-only, default `"log"` byte-identical; under `"identity"`: bias = mean(e − target), SD/SE on the natural scale, one-sided bound `e ± z1·se`, no positivity guard. Output column names unchanged (`bias_log` then holds the identity-scale bias; documented). `fs_plot_bias_coverage()` reads only `b`, `r`, coverages — scale-free, no change.
- **Fixture baseline (pre-change).** `fs_sim_bias_coverage()` on the seven committed `s7`/`map1` bundles vs `dev/tasks/bias_coverage_points.csv`: all 56 observed quantities within 0.005 (max 0.0045), 27 of 28 references within 0.001, the known `harm 1.75, n=500, field` c2_pred row at 0.00104 — identical to the adjudicated state in `REPORT_bias_coverage_display_2026-09-06.md` (fixture rounding). The 28 block × side tables from these bundles are saved for a byte-identity comparison after the change.

## 0c. Port map (survival template → continuous twin)

| survival template addition (`gbsg_020/sim_fs_maxeffCons_fb_mr_field_m1_template.qmd`) | where it goes in the continuous document |
|---|---|
| `.env_chr/.env_int/.env_num` and the `FS_S7_*` knobs (:227–420) | new `FS_MD_*` knobs in `setup-knobs`: `NSIMS`, `START`, `MODE`, `N`, `MD` (target, `null` for the null cell), `CAMPAIGN`, `FIELD_COMPLEMENT`, `IJ_RESIDUAL`, `FB` (`none`/`join`), `FB_PATH`, `WORKERS`, `QUICKRUN`, `SAVE_COMBINED`, `WINNER_ROWS` |
| `campaign_tag` in the stem (:371–379), `_quickrun` guard | stem `fs_maxeffCons_mr_field_md%02d_knoise%d_n%d_<campaign>` (token `mr_field`, so no glob can pool with the committed `mr_md40` bundles); `results_dir = mr_md_harm/<stem>_s<n>_d<draws>` as the twin |
| `.refuse_if_tracked()` (:423–435) before both `saveRDS` sites (:1300, :1417) | before the batch save (twin :876) and the combine save (:968) |
| MR knobs → `mr_inference_args` (:531–560): `ci_method = "field"`, `field_complement`, `ij_residual`, `return_reselection` | replaces the twin's inline list (:590–591); `ci_method` is a knob defaulting to `"field"` with `"ij"` available for the identity run |
| recorder: `fld_H_*` (23 cols), `fld_H_kappa…` (uniform, 8), `fld_Hc_*` (24), `mr_*_se_w/_lo_w/_hi_w/_se_wf/…` (12), `fld_joint_*` (9), `mr_harm_flag`, `label` (:750–810) | `.na_record()` (twin :463–490), appended after `ij_source`; plus `p_hat_H`, `p_top1..3`, `p_lab1..3` (GuoHe convention) |
| field extraction from `g$field` / `$complement` / `$joint` (:905–975) | `record_replicate()` after the complement block (twin :640) |
| `.ci_check()` interval invariant incl. field pairs (:1160–1195) | after the run loop (twin :846), before `fs_attach_betaHhat()` |
| `mr-settings-readout` chunk (:1457) | new chunk before `counts-header` |
| Table-2-layout rows: `est_keys_for()` adds `MR (field)` (est2, `lo2s/hi2s`, `se = fld_*_se`); winner rows hidden unless `WINNER_ROWS` | twin's `est_keys` (:1012) becomes per-block `est_keys_for(sfx)`; `.se_col_blk()` gains `fld`; `.build_block()` / `_med()` / boxplots / callout take the field row; bias stays absolute (MD units) with an added SD-units column |
| `field-coverage-wilson` (one-sided on the exposed side: H lower, Hc upper) (:2138) | ported with `e ± z95·se` in place of `exp(log e ± …)` |
| `field-bias-coverage` / `field-bias-coverage-complement` (display) (:2212, :2348) | ported with `scale = "identity"`; cell label `md=%d, n=%d` |
| `field-retained-bias`, `field-complement-diagnostics` (:2232, :2282) | ported on the MD scale (bias in units and in SD units; λ-SD/SD; IJ-SE/SD; one-sided margins and half-widths in MD units) plus p̂(Ĥ) summary |
| `field-joint-pair` (:2369) | ported (margins in MD units) |
| combine-mode poolability keys (:1364–1372): `ci_method`, `campaign_tag`, `target_hr_harm`, … | twin's key vector (:893) gains `ci_method`, `campaign_tag`, `target_md_harm`, `null_cell`, `field_complement`; meta gains the template provenance block |
| FB `"join"` mode with the naive-identity licence (:1245–1275) | ported; joins `fb_*` by `sim_id` from the committed FB bundle (md40 n500, 1–100) under `FS_MD_FB=join` |

The port is a **new document** beside the twin (the committed batch documents and bundles are not edited): `sim_fs_maxeffCons_mr_field_md_template.qmd`, driven by env knobs like the survival template.

## 0d. Mac cost anchors

- Hardware: Apple M4 Max, 14 physical cores (10 performance + 4 efficiency), 14 logical, 36 GB RAM. R 4.5.2, Quarto 1.10.18; speff2trial / patchwork / gt / doFuture present.
- The twin's own worker formula (`ceiling(0.90 × (14 − 1))`) gives 12; M-4's default is 13 (physical − 1). Both are measured under load in Stage 1c.
- **Single-replicate timing** (twin's `setup-knobs`, `build-dgm`, `machinery` chunks evaluated verbatim, `record_replicate(1)` in the main process, unloaded): DGM build 0.6 s; replicate 1 **2.4 s** wall (`fit_mr_secs` 2.4; `ci_method = "ij"`, 5,000 draws). Same replicate on the Linux bundle: 15.5 s under 115 concurrent workers.
- **Identity peek:** replicate 1 reproduced the committed md40 n = 500 row on every numeric pre-existing column (excluding timings) to a maximum relative difference of **3.6e-14**; `sg_def` identical (`!{cd40 <= 415} & !{cd80 <= 1022}`). Package version drift (bundles 0.2.2/0.3.1 vs 0.3.5) is the residual risk for the 5-per-cell identities in Stage 1b.
- Survival anchors for the field: ≈ 15 s/replicate at n = 500 for the field pass (R_out = 1000, R_in = 500) and ≈ 3 s for the complement, on the Linux box; the Mac numbers are measured in 1c.

## Gate 0

- 0a–0d resolved. The gate's field path runs on the continuous specification without change (identity `to_eff`, oriented working scale, all pass-throughs at the consistency call site).
- The only R/ change is the `scale` argument on `fs_sim_bias_coverage()` in `R/fs_bias_coverage.R` (add-only, default preserving).
- The committed bundles carry naive, oracle, MR (IJ) and complement columns for all four cells: the identity anchor exists.

**Gate 0: PASS.** Proceeding to Stage 1.
