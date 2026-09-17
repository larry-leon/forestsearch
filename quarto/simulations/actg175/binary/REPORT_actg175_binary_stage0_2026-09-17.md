# REPORT — ACTG175 binary (OR) simulation study under the current constructions: Stage 0 (read-only)

Date: 2026-09-17 (UTC). Machine: `pop-os` (64 physical cores, 251 GB; R 4.6.1). Branch `feature/glm-extension`.

**Tasks:**
- `dev/tasks/TASK_actg175_binary_stage0_2026-09-17_v2.md` ("v2", `d526f7ba`), which stopped at its §1 on the manuscript path (the earlier version of this record, `399e4b42`).
- `dev/tasks/TASK_actg175_binary_stage0_resume_2026-09-17.md` ("the resume task", `e064f017`), which executes v2 §2–§8 with the manuscript path restated (R1–R7).

**Read-only:** no `R/`, template, script, document or payload edit; no install; no render; no campaign. Nothing was written to fs-glms-interpretable or fs-post-selection. The only computation is §7, run in a temporary directory that has since been removed.

**Manuscript sources (R1):**
- `<ms>` = `~/Documents/GitHub/fs-glms-interpretable/dev/reference/post_selection/manuscript_jrssb_initialsubmit/`, written **[ms]** below.
- The secondary, `~/Documents/GitHub/fs-post-selection/jrssb_submission/`, written **[jrssb]**, was consulted only for a README. It has none (its contents: cover letter, proof, main and supplement PDFs, a `tex_zipped` archive). Its supplement PDF is byte-identical to [ms]'s (`cmp`).
- forestsearch paths are written **[fs]**, at HEAD.

## 1. Provenance — R2 GATE PASS

```
pop-os
feature/glm-extension
399e4b42                          (R2 HEAD; contains 399e4b42)
[forestsearch tracked modifications: 0]
[R / Rscript / quarto processes, by process name: none]
<ms> exists (fs-glms-interpretable; 232 tracked files under it)
fs-glms-interpretable  HEAD 74d3b45  status --porcelain lines 0
fs-post-selection      HEAD bc82bc3  status --porcelain lines 0
Built: R 4.6.1; ; 2026-09-17 04:47:31 UTC; unix
```
First commit: `e064f017 dev/tasks: TASK_actg175_binary_stage0_resume_2026-09-17.md as received`.

**How the DINA task ended.** `TASK_md_dina_campaign_2026-09-17` closed out completely:
- Gate 1 passed and the advance go applied (`REPORT_md_dina_stage1_2026-09-17.md`, `37395393`).
- Gate 2 passed 65 of 65 checks in each of the four cells, with no halt file (`REPORT_md_dina_gate2_2026-09-17.md`; campaign complete `30b282eb`, 27,727 s).
- Stage 3 is recorded in `REPORT_md_dina_2026-09-17.md` (`58be2fda`). The catalog closeout is `1d122b21`, and `check_current_status.sh --commit` passes.

**v2's first §1** (`399e4b42`) stopped because `dev/reference/post_selection/manuscript_jrssb_initialsubmit/` does not exist in forestsearch; R1 restates it.

## 2. S0.1 — Where the binary study lives

**2.1 In the manuscript** ([ms]):
- **Supplement includes.** `postselection_subgroup_inference_supplement.qmd:1130` `{{< include _sim_additional_results_intro.qmd >}}` carries the section "Additional simulation results" and its paragraph "**Binary outcomes on the odds-ratio scale.**" (`_sim_additional_results_intro.qmd:78–84`). `:1160` `{{< include _sim_mr_coverage_or075_actg175_conditional_target.qmd >}}` and `:1162` `{{< include _sim_mr_coverage_or075_actg175_marginal_target.qmd >}}` carry the two figures. `:1175` `{{{< include _sim_mr_coverage_or15_actg175.qmd >}}}` is escaped, so it is inactive.
- **Chunks.** `fig-mr-coverage-or075c` (conditional fragment `:273`) and `fig-mr-coverage-or075m` (marginal fragment `:273`). Both fragments read one grid:
  - `:41` `run_tag <- "maxeffCons_actg175_or075_seedtab_s1000"`
  - `:43` `rds_file <- sprintf("mr_coverage_grid_%s.rds", run_tag)`
  - `:47` `payload_dir <- "_payloads/"`
  - `:50` `target_pick <- "C_betaHhat"` (conditional) or `"C_dagger"` (marginal)
  - `:51` `estimator_pick <- c("MR", "naive")`
  - `:54–56` `bias_scale <- "rel"`, `bias_agg <- "median"`, `layout <- "grid2x2"`
- **Figure numbers.** The supplement renames figures with `\renewcommand{\thefigure}{S\arabic{figure}}` (`.tex:167`). In the rendered PDF (`pdftotext`), **Figure S9** is `fig-mr-coverage-or075c` (the conditional estimand θ(Ĥ)) and **Figure S10** is `fig-mr-coverage-or075m` (the marginal θ†). Both captions read *"(1,000 simulations per cell … NULL MR draws; detectors: consistency, dina, grf)"* (`.tex:3041–3056`, `:3068–3082`).
- **Payload.** `_payloads/mr_coverage_grid_maxeffCons_actg175_or075_seedtab_s1000.rds`. It is byte-identical (`cmp`) to [fs] `quarto/simulations/actg175/binary_020/mr_sweep/maxeffCons_actg175_or075_seedtab_s1000/mr_coverage_grid_maxeffCons_actg175_or075_seedtab_s1000.rds`. [ms] `_payloads/` also holds the older `mr_coverage_grid_actg175_or075_s{100,500,1000}.rds`.
- **README entry.** [ms] has no README; `payload_manifest.qmd` is its manifest. Its entry for this study (`:122–124`) names **`mr_coverage_grid_actg175_or075_s1000.rds`**, consumed by *"_sim_mr_coverage_or075_actg175_{marg,cond} (S8.7); both intros"*, built by *"forestsearch quarto/simulations/actg175/binary/mr_coverage_sweep_or075_s1k.qmd (run_tag actg175_or075_s1000)"*. That is not the grid the fragments read (finding F2). `:188–189` lists the inactive OR-1.5 grid as missing. [jrssb] has no README.

**2.2 In forestsearch** (`git log --all` over `*or075*`: the branch history is linear on `feature/glm-extension`, and no other ref adds a file):

| repo | ref | path | last commit (date) | tracked / present | role |
|---|---|---|---|---|---|
| fs-glms-interpretable | `74d3b45` | [ms] `postselection_subgroup_inference_supplement.qmd` (+ `.tex`, `.pdf`) | (dir) `47880a0` (2026-09-11) | yes / yes | supplement source, includes `:1130`, `:1160`, `:1162` |
| fs-glms-interpretable | `74d3b45` | [ms] `_sim_additional_results_intro.qmd` | (dir) `47880a0` | yes / yes | fragment: the binary paragraph (`:78–84`) and its live numbers (`:46–60`) |
| fs-glms-interpretable | `74d3b45` | [ms] `_sim_mr_coverage_or075_actg175_conditional_target.qmd`, `…_marginal_target.qmd` | (dir) `47880a0` | yes / yes | fragments: Figures S9, S10 |
| fs-glms-interpretable | `74d3b45` | [ms] `_supp_intro_sim_summary.qmd` | (dir) `47880a0` | yes / yes | fragment: front-matter summary reading the same grid (`:38`, `:67–73`) |
| fs-glms-interpretable | `74d3b45` | [ms] `_payloads/mr_coverage_grid_maxeffCons_actg175_or075_seedtab_s1000.rds` | (dir) `47880a0` | yes / yes | payload read by S9, S10 |
| fs-glms-interpretable | `74d3b45` | [ms] `_payloads/mr_coverage_grid_actg175_or075_s{100,500,1000}.rds` | (dir) `47880a0` | yes / yes | older payloads (read by no active fragment) |
| fs-glms-interpretable | `74d3b45` | [ms] `payload_manifest.qmd` | (dir) `47880a0` | yes / yes | README-equivalent (names the older grid) |
| fs-post-selection | `bc82bc3` | [jrssb] `postselection_subgroup_inference_supplement.pdf` | `5652143` (2026-08-20) | yes / yes | submitted PDF, identical to [ms]'s |
| forestsearch | HEAD | `quarto/simulations/actg175/binary_020/maxeffCons_mr_coverage_sweep_or075.qmd` | `1d42f6da` (2026-08-17) | yes / yes | **driver that produced S9/S10's grid**; DGM built inline |
| forestsearch | HEAD | `quarto/simulations/actg175/binary_020/maxeffCons_mr_coverage_sweep_or075.html` | `ba90bb1f` (2026-08-18) | yes / yes | its render (DGM truth printout, "Grid built 2026-08-18 02:00:21 from 21 cell(s)") |
| forestsearch | HEAD | `quarto/simulations/actg175/binary_020/mr_sweep/maxeffCons_actg175_or075_seedtab_s1000/` (21 `*_mr_n*_res.rds` + the grid) | `ba90bb1f` (2026-08-18) | yes (22) / yes | **payload** |
| forestsearch | HEAD | `quarto/simulations/actg175/binary_020/_sim_mr_coverage_or075.qmd` | `4b6690b3` (2026-08-18) | yes / yes | fragment copy |
| forestsearch | HEAD | `quarto/simulations/actg175/binary/maxeffCons_mr_coverage_sweep_or075.qmd` | `ac1cc50b` (2026-08-17) | yes / yes | earlier driver version (`n_sims 100`, `n_grid` by 1500, run_tag `maxeffCons_actg175_or075_s100`) |
| forestsearch | HEAD | `quarto/simulations/actg175/binary/mr_coverage_sweep_or075_s1k.qmd` (and `mr_coverage_sweep_or075{,_s100,_s500}.qmd`) | `f8bdb80f` (2026-08-07) | yes / yes | earlier drivers (run_tag `actg175_or075_s1000`, the manifest's) |
| forestsearch | HEAD | `quarto/simulations/actg175/binary/{effMaxSG,maxcons}_mr_coverage_sweep_or075.qmd` | `f8bdb80f` | yes / yes | rule variants (effMaxSG at ε 0.10: n 500 and 2000 only) |
| forestsearch | HEAD | `quarto/simulations/actg175/binary/mr_sweep/{actg175_or075_s100,actg175_or075_s500,legacy-actg175_or075_s1000,effMaxSG_actg175_or075_s1000,maxcons_actg175_or075_s1000,maxeffCons_actg175_or075_s1000}/` | `4cc13ec2`–`47eb3f6a` (2026-07-02 to 08-04) | yes / yes | older payloads |
| forestsearch | HEAD | `quarto/simulations/actg175/binary/mr_sweep_legacy/*` (or075, or10, or15, or15_effMaxSG) | — | yes / yes | legacy payloads (harm branch OR 1.5 present) |
| forestsearch | HEAD | `quarto/simulations/actg175/binary/build_actg175_glm_dgm.R` | `d135997a` (2026-06-24) | yes / yes | standalone DGM builder (OR 1.5, n_super 25,000); **not** used by the S9/S10 driver |
| forestsearch | HEAD | `quarto/simulations/actg175/binary/_sim_mr_coverage_or075.html` | `96a677bb` (2026-06-25) | yes / yes | earlier fragment render |

No record (`REPORT_*`) for this study existed before this one.

**2.3 Which files produced Figures S9 and S10.** The driver [fs] `quarto/simulations/actg175/binary_020/maxeffCons_mr_coverage_sweep_or075.qmd` (`1d42f6da`) wrote the 21 cell bundles and the grid under `binary_020/mr_sweep/maxeffCons_actg175_or075_seedtab_s1000/` (committed `ba90bb1f`). The grid was copied into [ms] `_payloads/`, where the two fragments render S9 (conditional) and S10 (marginal).

**2.4 Record location** (R3): `quarto/simulations/actg175/binary/`. The producing driver and payload are in the sibling directory `binary_020/` (finding F1).

## 3. S0.2 — The design as committed

Source: [fs] `binary_020/maxeffCons_mr_coverage_sweep_or075.qmd` at HEAD (`:line`), its render, and the bundles' `meta`.

**Data and subgroup.**
- **Trial and arms:** `speff2trial::ACTG175`, arms 1 vs 3, `treat = (arms == 1)` (ZDV+ddI = 1, ddI = 0); missing `cd420` dropped (`:235–238`).
- **Outcome:** `y_neg <- 1L - as.integer(cd420 > cd40)`, "no CD4 improvement" (`:239`).
- **Covariates:** `cont_vars` age, preanti, wtkg, karnof, cd40, cd80; `bin_vars` hemo, homo, drugs, race, gender, symptom (`:158–159`; `str2` is not in the pool).
  - The DGM builds `bin_vars` as factors.
  - The recorder coerces them to 0/1 numeric before every `forestsearch()` call: `for (v in bin_vars) if (… is.factor(df[[v]])) df[[v]] <- as.numeric(as.character(df[[v]]))` (`:465–467`).
  - The evaluation frame is coerced the same way (`:335–337`).
- **H:** `subgroup_cuts <- list(wtkg = list(type = "greater", quantile = 0.70), cd40 = list(type = "greater", quantile = 0.70))` (`:164–168`).
- **Build:** `calibrate_glm_interaction(…, outcome_var = "y_neg", target_effect = 0.75, outcome_type = "binary", effect_measure = "OR", k_treat = 1, adverse_outcome = FALSE, k_inter_range = c(0.3, 1.5), grid_step = 0.025, n_super = 100000L, seed = 8316951)` (`:249–266`).
  - `adverse_outcome = FALSE` there is the construction direction; the analysis uses `adverse_outcome <- TRUE`, with OR > 1 = harm (`:176`).
  - Calibrated `k_inter` (`dgm$model_params$beta_inter`, from §7's rebuild): 0.1480184.
- **Prevalence and super-population:** prevalence(H) 9.63% (render); super-population 100,000 (`:152`).

**Truths** (grid `G$truth`; render; `:268–274`):

| | θ† (marginal OR) | θ‡ (CDE) |
|---|---|---|
| H | 0.7499999955 | 0.7321189340 |
| Hᶜ | 0.6560116720 | 0.6313905111 |
| overall marginal OR | 0.6665537395 | |

- **Conditional target column:** `C_betaHhat` in the grid, from the per-replicate `betaHhat_H` / `betaHhat_Hc` attached by `fs_attach_betaHhat(results, eval_df, focus = "harm", outcome_type = "binary", effect_measure = "OR", …)` (`:705–709`). The evaluation frame is the full 100,000-subject super-population, `eval_seed = 20260628` (`:333–334`).
- **Evaluation-frame check:** the render prints *"theta-dagger check : 0.741 / 0.650 (eval frame at true H/Hc; matches marginal OR above)"* against the truths 0.750 / 0.656 (finding F6).

**Designs in the payloads:**
- **OR-0.75 protective grid:** present (this study).
- **Null branch:** present in the driver (`dgm_model <- "null"`, `generate_glm_dgm(model = "null")`, `:276–310`), but no null payload for this run tag.
- **OR-1.5 harm branch:** present in other committed payloads (`binary/mr_sweep_legacy/actg175_or15_s500`, `actg175_or15_effMaxSG_s{100,1000}`); absent from this run tag.

**The n grid:** `seq(500L, 2000L, by = 250L)`, i.e. 500, 750, …, 2000 (`:71`); 3 identifiers × 7 sizes = 21 cells (grid `manifest`, 21 rows).

**Replicates per cell actually run: 1,000.** Every one of the 21 bundles has 1,000 rows (sim_id 1–1000) and `meta$n_sims = 1000`; the driver sets `n_sims <- 1000L` (`:83`); the render prints *"Counts: n_sims=1000, nb_boots=NULL (MR-only=yes), mr_draws=5000."*
- The supplement text's "500 simulations per cell" (`_sim_additional_results_intro.qmd:80`) does not match.
- The S9/S10 captions' "1,000" does.

**MR draws: 5,000**, multiplier Poisson (the `fs_mr_inference()` default; the driver passes none):
- `mr_draws <- 5000L` (`:85`), passed as `mr_inference_args = list(ci_method = "ij", draws = mr_draws, include_complement = TRUE)` (`:510–511`).
- Every bundle's `meta$mr_draws` and the grid's `coverage$mr_draws` are 5000.
- **The caption's "NULL"** comes from the fragment: `gd_txt <- .rng(sl$gate_draws)` (fragment `:110`) reads a column `gate_draws` that the grid does not have (its column is `mr_draws`), and the caption prints it (finding F3).

**Seeds and RNG kind:**
- A pre-generated table: `seed_base <- 8316951L; set.seed(seed_base); SEED_TABLE <- sample.int(.Machine$integer.max - 1L, 5000L)`, looked up by global `sim_id` (`:124–133`), under R's default generator (Mersenne-Twister / Inversion / Rejection).
- Per replicate: `RNGkind("L'Ecuyer-CMRG"); set.seed(sd_i)`, then `simulate_from_glm_dgm(dgm, n, seed = sd_i)` (`:446–450`), and `seedit = sd_i` for the fit (`:489`).
- The DGM is built with `seed = seed_base` before any kind switch.
- `meta$seed_scheme` = "pre-generated table indexed by global sim_id".
- sim_id 1's seed is 1530735852.

**Settings per detector** (`:174–207`, `:484–522`; bundle `meta`):
- **Common:** `outcome_type = "binary"`, `effect_measure = "OR"`, `adverse_outcome = TRUE`; `sg_focus = "maxeffCons"`; `selection_rule = "neighborhood"`; `effect_neighborhood = 0.10`; `stop_threshold = NULL`.
- **Thresholds (OR scale):** `effect.threshold = 0.90`, `consistency.threshold = 0.80`, `pconsistency.threshold = 0.90`.
- **Search:** `n.min = 60`, `d0.min = d1.min = 10`, `maxk = 2`, `conf.cont_jcuts = list(cd40 = 10, wtkg = 10)`, `max_subgroups_search = Inf`.
- **Consistency:** `consistency_method = "resample"`, `fs.splits = 500`, `use_twostage = TRUE`, `use_lasso/use_grf/use_dina = FALSE`, `vi.grf.min = -0.2`.
- **GRF:** `grf_selection = "frontier"`, `grf_select_statistic = "effect"`, `grf_depth = 2`, `dmin.grf = 0`.
- **DINA:** `dina_args = list()`, `dina_select_statistic = "effect"`.

**Estimators recorded** (`:405–429`, `:577–631`):
- naive `nv_*`, and MR (IJ) `t2_*` with `t2_*_se_ij`, on H and Hᶜ.
- oracle `ora_*`: a logistic refit on the true H / Hᶜ, OR-scale Wald CI, the log-OR SE (`:363–380`, `:628–631`).
- FB `t1_*`: wired but dormant (`nb_boots <- NULL`, `:84`; all NA).
- No field, field-s, joint, or p̂ columns.

**Detection and targets:**
- **Detection:** `detected = 1` when `fs.est$sg.harm` is non-empty (`:568–571`), with the `DETECTED` / `NO-DETECTION` / `CONFIG-ERROR` status. All 21 bundles have zero `CONFIG-ERROR` rows.
- **Summaries:** coverage and bias are computed over detected replicates (`rH <- res[res$detected %in% 1L, ]`, `:805`).
- **Targets:** `C_dagger` = `truth$marg_<block>`, `C_ddagger` = `truth$cde_<block>`, `C_betaHhat` = `betaHhat_<block>`, `C_oracle` = the per-replicate `ora_<block>_est` (`:809–813`); coverage is the two-sided interval `[lo, hi]` (`:768`).
- **Detection rates** (grid manifest), n = 500 → 2000:
  - consistency 0.768, 0.829, 0.894, 0.891, 0.896, 0.893, 0.903;
  - DINA 0.935, 0.936, 0.928, 0.919, 0.893, 0.870, 0.848;
  - GRF 0.998, 0.999, 1.000, 1.000, 0.999, 1.000, 0.999.

**Per-replicate cost as recorded.** Machine: `meta`: host `Mac-Studio-M1-Ultra.local`, 12 workers, R 4.5.2, forestsearch 0.2.0 at `1d42f6da`, `parallel_mode = "sims"`. MR ran with `ci_method = "ij"` (no field pass).

| detector | `t2_secs` mean by n = 500 / 750 / 1000 / 1250 / 1500 / 1750 / 2000 (s) | cell `elapsed_sec` (s) |
|---|---|---|
| consistency | 5.8 / 8.1 / 10.5 / 12.4 / 14.1 / 15.6 / 17.3 | 513 / 700 / 907 / 1059 / 1208 / 1342 / 1479 |
| DINA | 2.1 / 1.8 / 1.6 / 1.5 / 1.5 / 1.4 / 1.3 | 209 / 199 / 167 / 157 / 151 / 143 / 131 |
| GRF | 4.3 / 5.6 / 6.8 / 7.8 / 8.5 / 9.6 / 10.8 | 374 / 483 / 580 / 663 / 726 / 818 / 919 |

## 4. S0.3 — The current package on the binary path (source at HEAD; `R/` last changed `064fce91`)

**The field for a binomial GLM:**
- **Per-candidate pieces for logistic outcomes.** `.fs_mr_pieces()` dispatches non-survival outcomes to `.consistency_glm_pieces()` (`R/fs_mr_inference.R:61–73`).
  - That function fits `OR = stats::glm(…, family = stats::binomial("logit"))` (`R/consistency_resample.R:271`).
  - It returns `beta_hat` (log-OR), the `.dfbeta_glm()` influence and `log_scale = effect_measure %in% c("OR", "RR", "IRR")` (`:323`).
- **Scale.** The field is computed on the log-OR working scale and reported as OR through `to_eff <- function(x) if (log_scale) exp(x) else x` (`R/fs_mr_inference.R:607`), e.g. `lower_1s = to_eff(beta_deb - qs[5])` (`:932`).
- **Outcome-type guard:** none. The field block is gated only by `if (ci_method == "field")` (`:862`), and `fs_mr_inference()` reads `outcome_type` only in the pieces dispatch (`:62`).
- **MR's `consistency_method` requirement.** On a GLM outcome the FS branch runs MR only under `consistency_method = "resample"`: `.mr_glm_ok <- consistency_method == "resample" && !is.null(estimator_fn)` (`R/forestsearch_main.R:3318`), with the skip message *"MR on a GLM outcome requires consistency_method = "resample"…"* (`:3333–3336`). DINA and GRF have no such condition.

**Field-s and `joint_s`.** `.fs_mr_field_complement()` returns the `_s` companions and `joint_s` whenever `field_scale_complement = "selected"` (`:1172–1174`, `:1291–1292`, `upper_1s_s = to_eff(bdc - qss[1])` `:1312`), with no outcome-type branch. They carry the same fields as on continuous, on the OR scale. Observed: all six §7 fits return `field$lower_1s`, `field$complement$upper_1s_s` and `field$joint_s$bonf_lower_H` / `bonf_upper_Hc`, all finite.

**`adverse_outcome` for binary:**
- **Default:** `adverse_outcome = NULL` resolves to TRUE for binary (`R/forestsearch_main.R:1287`, `:1343`).
- **FS:** the effect estimator flips `Y -> 1 - Y` when `adverse_outcome = FALSE` (`R/glm_effect_estimators.R:303–305`), as do the resample pieces (`R/consistency_resample.R:255–261`), so OR > 1 is harm under either value.
- **GRF:** `grf.subg.harm.glm()` uses `Y_grf <- 1L - Y_grf` when `adverse_outcome = TRUE` (`R/grf_subg_harm_glm.R:443–446`). Its frontier harm effect is `mean(Γ₀) − mean(Γ₁)` (`:510–514`), harm-oriented under either value.
- **DINA (after P2):** `.dina_tau_sign()` (`R/forestsearch_helpers.R:1376`) returns −1 for binary with `adverse_outcome = FALSE` and 1 otherwise; it is applied to tau-hat before the `m_diff` floor (`:1441`, `:1496`). Under this study's `adverse_outcome = TRUE` the sign is 1 and the log-OR floor `log(hr.threshold)` applies to the raw fit (`:1436–1437`).

**Factor covariates:**
- **GRF (after P1):** `.grf_code_column()` (`R/grf_helpers.R:719`) codes each column as the forest matrix does, for both the forest and the evaluator.
- **DINA:** `.coerce_covariates_numeric()` (`R/forestsearch_helpers.R:1327–1345`, applied at `:1421–1423`) turns all-numeric-level factors into numeric.
- The study's driver also coerces `bin_vars` itself (§3).

**The `effMaxSG` band for binary:**
- **Expression:** `.compute_inclusion_band()` (`R/subgroup_consistency_helpers.R:778–799`) keeps `hr_vec >= (1 - effect_neighborhood) * max(hr_vec)` on the natural effect.
- **Identifier side:**
  - FS: `if (isTRUE(effect_log_scale)) hr_vec <- exp(hr_vec)` (`:590`, `:679`, `:863`).
  - DINA: `effect_log_scale <- fit$family %in% c("binomial", …)` then `exp(cand_tau)` (`R/dina_subgroup.R:406`, `:518`), and `.dina_reselect_on_effect()` uses `eff <- if (log_scale) exp(eff_link)` (`R/forestsearch_helpers.R:1245`).
  - GRF: `.grf_reselect_on_effect()` uses `cand_hr$effect <- if (log_scale) exp(eff_link)` (`:1658`).
- **MR side:** `.fs_mr_select()` uses `eff <- if (log_scale) exp(beta[passers])` in the same helper (`R/fs_mr_inference.R:136–142`).
- **Aligned:** the same helper, on the OR scale, on both sides.

**`fs_sim_bias_coverage()` on the OR scale:**
- **Scale:** `scale = c("log", "identity")` (`R/fs_bias_coverage.R:79`); `"log"` is the ratio-measure scale.
- **Columns:** it reads `nv_<b>_*`, `mr_<b>_*` and `fld_<b>_*` directly (`:18–20`). Field-s needs the renamed copy used on the continuous path, because it reads `fld_Hc_*`, not `fld_Hc_*_s`.
- **Oracle target:** `target = "oracle"` reads `or_<block>_est` (`:111`), so bundles with the `ora_` prefix need a rename.

## 5. S0.4 — What a re-run needs (listed, not made)

Mapping: the committed MD machinery — `quarto/simulations/actg175/continuous/sim_fs_maxeffCons_mr_field_md_template.qmd` (identifier knob `FS_MD_METHOD`, field and field-s recorder, DINA and GRF fields), `scripts_mdsgnb20/`, `scripts_mdgrf/`, `scripts_mddina/`, `summary_continuous_field_md*.qmd` and the `md_*_metrics.csv` extracts — set against the binary driver above.

| gap | MD machinery | binary driver |
|---|---|---|
| outcome family and DGM builder | `calibrate_glm_interaction()` / `generate_glm_dgm()` on continuous `cd4_change`, `factor_vars = z1..z12`, H from `z1`/`z2` (age > 34, preanti ≤ 744.5) | `calibrate_glm_interaction()` on binary `y_neg` with `continuous_vars` and upper-quantile cuts on wtkg and cd40; construction `adverse_outcome = FALSE`, analysis TRUE |
| targets | `betaHhat_H/Hc` (raw) and `truth$effect_Q/Qc` (structural), `marg_*` fitted | θ† (`marg_*`) and θ‡ (`cde_*`) beside `betaHhat_*`; no structural/fitted pair |
| oracle prefix | `or_*` | `ora_*` (`fs_sim_bias_coverage(target = "oracle")` reads `or_*`) |
| MR and FB prefixes | `mr_*`, `fb_*` | `t2_*`, `t1_*` |
| recorder columns | `fld_H_*`, `fld_Hc_*`, nine `fld_Hc_*_s`, `fld_joint_*`, nine `fld_joint_s_*`, p̂ and labels, `n_family`, `admitted_n`, `dina_*`, `warn_msg`, `id_secs` | none of these (MR ran with `ci_method = "ij"`) |
| seeds | `seed_base + sim_id` | a pre-generated table indexed by `sim_id` |
| rule | `effMaxSG`, ε 0.20, `neighborhood` | `maxeffCons`, ε 0.10, `neighborhood` |
| thresholds | MD 30 / 10, harm-oriented | OR 0.90 / 0.80 (sub-null) |
| bound-location threshold scale | the MD ladder τ = 0 … 100 on the harm-oriented MD | needs an OR-scale ladder (the field bounds are ORs) |
| summary scale and orientation | `scale = "identity"`, orientation −1 on `betaHhat_*` | `scale = "log"`, no flip (`adverse_outcome = TRUE`, OR > 1 = harm) |
| parallel and campaign mechanics | batch / combine renders, runner, Gate 2, extract | one document renders all 21 cells, resumable per cell |

**Harm and null designs, as facts:**
- **Harm:** the builder takes any `target_effect` within the reach of `k_inter_range`. Committed binary drivers set `target_or_h` to 1.0 (`mr_coverage_sweep_or10.qmd:125`), 1.5 (`mr_coverage_sweep_or15.qmd:125`, `…_or15_effMaxSG.qmd:125`), 2.0 (`…_or20.qmd:131`) and 3.5 (`…_or35_effMaxSG.qmd:125`), all with `k_inter_range = c(0.3, 1.5)`. Payloads exist for OR 1.0, 1.5 and 3.5 (`binary/mr_sweep_legacy/`, `binary/mr_sweep/`); `actg175_or20_s1000` holds 3 files and no grid.
- **Null:** the driver's `dgm_model = "null"` branch (`generate_glm_dgm(…, model = "null", subgroup_vars = NULL)`, `:279–294`) sets every block's truth to the population OR / CDE.
- **Standalone builder:** `build_actg175_glm_dgm.R` (OR 1.5, n_super 25,000; its header: prevalence 9.88%, θ†(H) 1.505, θ†(Hᶜ) 0.656).

## 6. S0.5 — Cost references (quoted; nothing run for this section)

- **This study** (§3): MR with IJ only, on a Mac Studio M1 Ultra with 12 workers. `t2_secs` mean 5.8–17.3 s for consistency, 1.3–2.1 s for DINA and 4.3–10.8 s for GRF, over n = 500–2000.
- **`mdsgnb20`** (FS on the MD design, field on; pop-os, 63 workers): `fit_mr_secs` mean 56.24 / 68.04 / 53.91 / 84.88 s (md40 n500 / md120 / null / md40 n700). Source: `REPORT_md_grf_stage1_resume_2026-09-16.md`'s projection table, read from the `mdsgnb20` bundles; the campaign's own record is `REPORT_md_field_rerun_2026-09-15.md`.
- **`mdgrf`** (`REPORT_md_grf_2026-09-16.md`): 7,423 s for 8,000 replicates; `fit_mr_secs` mean 45.1 / 51.2 / 44.1 / 53.5 s, of which field 27.2 / 32.3 / 26.3 / 29.3 s.
- **`mddina`** (`REPORT_md_dina_2026-09-17.md`): 27,727 s for 8,000 replicates; `fit_mr_secs` mean 134.2 / 304.9 / 100.3 / 199.7 s, of which field 69.8 / 149.2 / 54.4 / 88.0 s.
- **§7's six fits** (pop-os, one worker each, field on, 5,000 draws): below.

## 7. S0.6 — Data check and one replicate per detector

**Condition — holds.**
- **Build date:** the installed `Built` (`2026-09-17 04:47:31 UTC`) postdates the last `R/` commit (`064fce91`, 2026-09-17 04:47:07 UTC).
- **Namespace comparison:** HEAD's source was extracted with `git archive` into `$(mktemp -d)` and loaded with `pkgload::load_all(export_all = TRUE)`. Functions were compared after `utils::removeSource()` with environments neutralized. Result: *"installed objects: 678 | HEAD source objects: 678 … common: 678 | differing: 0"*.
- **Inputs:** §2 located the recipe and the payload.

**7.1 Data check — GATE PASS.** The OR-0.75 design was rebuilt with the committed recipe (§3; seed table under R's default generator, the DGM before any kind switch, L'Ecuyer-CMRG per replicate). It was compared with the committed `fs_mr_n500_res.rds` and `fs_mr_n2000_res.rds`, sim_id 1 (the three identifiers share the draws).
```
forestsearch 0.3.5 | Built R 4.6.1; ; 2026-09-17 04:47:31 UTC; unix | R 4.6.1 | RNGkind at start: Mersenne-Twister/Inversion/Rejection
DGM rebuilt in 7.6 s | prevalence(H) 0.0963 | k_inter (beta_inter) 0.1480184 | N_super 100000
  truth or_causal rebuilt 0.6665537395 | committed 0.6665537395 | rel diff 8.33e-16
  truth marg_H    rebuilt 0.7499999955 | committed 0.7499999955 | rel diff 4.44e-16
  truth marg_Hc   rebuilt 0.6560116720 | committed 0.6560116720 | rel diff 1.18e-15
  truth cde_H     rebuilt 0.7321189340 | committed 0.7321189340 | rel diff 3.03e-16
  truth cde_Hc    rebuilt 0.6313905111 | committed 0.6313905111 | rel diff 1.41e-15
DATA n=500 sim_id 1: seed 1530735852 (committed 1530735852) | n_true 47 (committed 47) | ora_H_est 1.6153846154 (committed 1.6153846154) | ora_Hc_est 0.7397959184 (committed 0.7397959184) | max rel diff 5.10e-15
DATA n=2000 sim_id 1: seed 1530735852 (committed 1530735852) | n_true 189 (committed 189) | ora_H_est 0.6467351431 (committed 0.6467351431) | ora_Hc_est 0.6459288688 (committed 0.6459288688) | max rel diff 1.55e-14
DATA CHECK GATE: data-level columns within 1e-8: TRUE | truth targets within 1e-8: TRUE
```

**7.2 Fits.**
- **What ran:** `forestsearch()` on those two replicates for FS, GRF and DINA, one worker each, with a 30-min limit each (no fit reached it). The committed thresholds, cut grid and detector arguments (§3) were used, with:
  - `sg_focus = "effMaxSG"`, `effect_neighborhood = 0.20`;
  - `selection_rule = "neighborhood"` (the `mdsgnb20` bundles' `meta`);
  - MR on with `ci_method = "field"`, `draws = 5000L`, `include_complement = TRUE`, `field_complement = TRUE`, `field_scale_complement = "selected"`, `return_reselection = TRUE`.
- **Columns reported:**
  - "MR family" is `mr_inference$n_family`.
  - "admitted/qualifying" is the consistency-qualifying table size for FS, `grf_res$admitted_n` for GRF, and `out_sg$admitted_n` for DINA (with `proposed` = DINA's proposal count).
  - The oriented log-OR is `log(naive OR)` (`adverse_outcome = TRUE`, so OR > 1 is harm).
  - Field bounds are on the OR scale.
  - The committed selection beside each fit was made under the study's own rule (`maxeffCons`, ε 0.10, IJ), so it is a reference, not a check.
```
FIT consistency n= 500 | wall 39.6 s | error: none | warnings: none
    selection [{race} & !{karnof <= 95}] n 73 | oriented log-OR (naive) 0.8303 | MR family 2238 | admitted/qualifying 37 | field$lower_1s TRUE (0.2937) | complement$upper_1s_s TRUE (1.1332) | joint_s bonf lower_H TRUE upper_Hc TRUE
    committed payload (maxeffCons, eps 0.10, ci_method ij), sim_id 1: status DETECTED | sg_def [{symptom} & {wtkg <= 79}] | n_sel 61 | naive OR 2.377622 | t2_secs 6.0
FIT grf         n= 500 | wall 26.6 s | error: none | warnings: none
    selection [{age <= 28} & {preanti <= 777.20000000000005}] n 89 | oriented log-OR (naive) 0.6484 | MR family 1051 | admitted/qualifying 330 | field$lower_1s TRUE (0.2920) | complement$upper_1s_s TRUE (1.1375) | joint_s bonf lower_H TRUE upper_Hc TRUE
    committed payload (maxeffCons, eps 0.10, ci_method ij), sim_id 1: status DETECTED | sg_def [{karnof > 90} & {race > 0}] | n_sel 73 | naive OR 2.294118 | t2_secs 3.3
FIT dina        n= 500 | wall 28.1 s | error: none | warnings: none
    selection [{age <= 28} & {preanti <= 777.20000000000005}] n 89 | oriented log-OR (naive) 0.6484 | MR family 1157 | admitted/qualifying 842 | proposed 1157 | field$lower_1s TRUE (0.3028) | complement$upper_1s_s TRUE (1.1313) | joint_s bonf lower_H TRUE upper_Hc TRUE
    committed payload (maxeffCons, eps 0.10, ci_method ij), sim_id 1: status DETECTED | sg_def [{karnof >= 100} & {race >= 1}] | n_sel 73 | naive OR 2.294118 | t2_secs 3.0
FIT consistency n=2000 | wall 90.5 s | error: none | warnings: none
    selection [!{preanti <= 741} & !{wtkg <= 78}] n 145 | oriented log-OR (naive) 0.3844 | MR family 2978 | admitted/qualifying 3 | field$lower_1s TRUE (0.2824) | complement$upper_1s_s TRUE (0.7565) | joint_s bonf lower_H TRUE upper_Hc TRUE
    committed payload (maxeffCons, eps 0.10, ci_method ij), sim_id 1: status DETECTED | sg_def [{cd40 <= 210} & !{wtkg <= 70}] | n_sel 107 | naive OR 1.813268 | t2_secs 14.6
FIT grf         n=2000 | wall 48.3 s | error: none | warnings: none
    selection [{preanti > 0} & {wtkg <= 64.5}] n 276 | oriented log-OR (naive) 0.1354 | MR family 1345 | admitted/qualifying 49 | field$lower_1s TRUE (0.2127) | complement$upper_1s_s TRUE (0.7459) | joint_s bonf lower_H TRUE upper_Hc TRUE
    committed payload (maxeffCons, eps 0.10, ci_method ij), sim_id 1: status DETECTED | sg_def [{preanti > 842} & {race > 0}] | n_sel 74 | naive OR 1.4 | t2_secs 7.7
FIT dina        n=2000 | wall 11.1 s | error: none | warnings: none
    selection [{wtkg >= 68.040000000000006} & {cd40 <= 214}] n 129 | oriented log-OR (naive) 0.3589 | MR family 70 | admitted/qualifying 19 | proposed 70 | field$lower_1s TRUE (0.3015) | complement$upper_1s_s TRUE (0.7501) | joint_s bonf lower_H TRUE upper_Hc TRUE
    committed payload (maxeffCons, eps 0.10, ci_method ij), sim_id 1: status DETECTED | sg_def [{wtkg >= 71} & {cd40 <= 214}] | n_sel 111 | naive OR 1.561966 | t2_secs 1.2
```

| fit | wall (s, 1 worker) | warnings | selection (n) | MR family | admitted | field lower_1s / field-s upper_1s_s / joint_s Bonferroni |
|---|---|---|---|---|---|---|
| FS n 500 | 39.6 | none | `{race} & !{karnof <= 95}` (73) | 2,238 | 37 qualifying | present, finite |
| GRF n 500 | 26.6 | none | `{age <= 28} & {preanti <= 777.2}` (89) | 1,051 | 330 | present, finite |
| DINA n 500 | 28.1 | none | `{age <= 28} & {preanti <= 777.2}` (89) | 1,157 | 842 of 1,157 proposed | present, finite |
| FS n 2000 | 90.5 | none | `!{preanti <= 741} & !{wtkg <= 78}` (145) | 2,978 | 3 qualifying | present, finite |
| GRF n 2000 | 48.3 | none | `{preanti > 0} & {wtkg <= 64.5}` (276) | 1,345 | 49 | present, finite |
| DINA n 2000 | 11.1 | none | `{wtkg >= 68.04} & {cd40 <= 214}` (129) | 70 | 19 of 70 proposed | present, finite |

7.3 The temporary directory (`/tmp/tmp.rYm8ZDL99s`, holding the extracted source, the two scripts and their output) was removed after this record's numbers were taken.

## 8. Facts for Larry's decisions (facts only)

- **The supplement's design as committed:**
  - **Cells:** the ACTG175 binary OR-0.75 protective design, H = {wtkg > q70} ∩ {cd40 > q70}, prevalence 9.63%.
  - **Grid and replicates:** FS, DINA and GRF at n = 500, 750, …, 2000 (21 cells), with **1,000 replicates per cell**. The text's 500 is wrong; the captions' 1,000 is right.
  - **Draws:** **5,000 MR draws, Poisson multiplier**. The captions' "NULL" is a fragment column-name slip.
  - **Truths:** θ†(H) 0.750, θ†(Hᶜ) 0.656, θ‡(H) 0.732, θ‡(Hᶜ) 0.631.
  - **Rule and thresholds:** `maxeffCons` with ε 0.10, `neighborhood`, OR thresholds 0.90 / 0.80, `adverse_outcome = TRUE`.
  - **Estimators:** MR with IJ intervals only, naive and oracle; no FB.
  - **Targets:** θ†, θ‡, θ(Ĥ) and the oracle, over detected replicates.
  - **Provenance:** built on a Mac M1 Ultra (12 workers) with forestsearch 0.2.0 at `1d42f6da`; payload `binary_020/mr_sweep/maxeffCons_actg175_or075_seedtab_s1000/`.
- **What the package supports on binary today:**
  - The field, field-s and the `joint_s` pair run on the binary OR path with no outcome guard, on the log-OR scale reported as OR, and are finite on all six §7 fits.
  - MR on FS needs `consistency_method = "resample"`.
  - The `effMaxSG` band is the same helper on the natural OR scale on both the identifier and MR sides.
  - DINA's orientation under binary is 1 at `adverse_outcome = TRUE` and −1 at FALSE (P2).
  - GRF codes factor covariates consistently (P1).
- **Rule and settings previously used vs the current campaign standard:** `maxeffCons` at ε 0.10 with IJ-only MR and no field, against `effMaxSG` at ε 0.20 with `neighborhood`, the field on Ĥ, field-s on Ĥᶜ and the Bonferroni pair. At the current standard, the six fits select differently from the committed rows (§7.2).
- **Template gaps:** the §5 table.
- **Cost per replicate:**
  - Recorded (Mac, 12 workers, IJ only): consistency 5.8–17.3 s, DINA 1.3–2.1 s, GRF 4.3–10.8 s.
  - Measured today (pop-os, 1 worker, field on, sim_id 1): FS 39.6 / 90.5 s, GRF 26.6 / 48.3 s, DINA 28.1 / 11.1 s at n = 500 / 2000.
  - MD campaigns at 63 workers: FS 54–85 s, GRF 44–53 s, DINA 100–305 s.
- **Harm and null designs available:** OR 1.0, 1.5, 2.0 and 3.5 drivers are committed (payloads for 1.0, 1.5 and 3.5), and the driver has a null branch (`dgm_model = "null"`). All use the same builder, with `k_inter_range = c(0.3, 1.5)`.

## 9. Record location and catalog

This record is at `quarto/simulations/actg175/binary/` (R3). The directory has no `current_status.md` generator, so there is no catalog step (R4).

## 10. Findings

- **F1.** The driver and payload behind Figures S9 and S10 are in `quarto/simulations/actg175/binary_020/`, not `binary/`. `binary/maxeffCons_mr_coverage_sweep_or075.qmd` is an earlier 100-replicate, two-size version of the same driver.
- **F2.** [ms] `payload_manifest.qmd:122–124` names the grid `mr_coverage_grid_actg175_or075_s1000.rds` and the driver `binary/mr_coverage_sweep_or075_s1k.qmd`. The active fragments read `mr_coverage_grid_maxeffCons_actg175_or075_seedtab_s1000.rds`, from `binary_020/maxeffCons_mr_coverage_sweep_or075.qmd`. Both grids are in [ms] `_payloads/`.
- **F3.** The S9/S10 captions' "NULL MR draws" comes from `gd_txt <- .rng(sl$gate_draws)` (fragment `:110`): the grid's column is `mr_draws` (5000). The paragraph's "500 simulations per cell" (`_sim_additional_results_intro.qmd:80`) disagrees with the 1,000 rows per cell in every bundle.
- **F4.** The study ran MR with `ci_method = "ij"`: no field, field-s or joint columns exist in its payload, so a re-run is needed to report the current constructions.
- **F5.** The study's thresholds are sub-null on the OR scale (0.90 / 0.80) and its rule is `maxeffCons` at ε 0.10. The OR-0.75 H is a protective region, so "harm" selections here are false in direction by design.
- **F6.** The render's evaluation-frame check prints θ† 0.741 / 0.650 at the true H / Hᶜ against the calibrated 0.750 / 0.656; the driver labels it a match.
- **F7.** At n = 500, GRF and DINA select the same subgroup on sim_id 1 (`{age <= 28} & {preanti <= 777.2}`, n 89) under the current standard. Their committed rows also agree with each other under the old rule (`{karnof > 90} & {race > 0}` / `{karnof >= 100} & {race >= 1}`, n 73).
- **F8.** At n = 2000, DINA's proposed family on sim_id 1 is 70 candidates (19 admitted), against 1,157 at n = 500 — the family depends on the fitted surface.
- **F9.** [ms] and [jrssb] have no README; `payload_manifest.qmd` is the only manifest.
- **F10.** v2's `ps | grep` pattern matches this session's shell wrappers and an unrelated 11-day-old sampling loop; by process name nothing was running.

## Commits

```
d526f7ba dev/tasks: TASK_actg175_binary_stage0_2026-09-17_v2.md as received
399e4b42 ACTG175 binary Stage 0 v2: STOP at §1 (the first version of this record)
e064f017 dev/tasks: TASK_actg175_binary_stage0_resume_2026-09-17.md as received
<this record, completed>
```
