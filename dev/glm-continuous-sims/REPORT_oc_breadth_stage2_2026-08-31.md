# REPORT — breadth, stage 2: the measured cell at the forecast harm, scored against the locked forecast

**Task:** `dev/tasks/cc_task_oc_breadth_stage2_2026-08-31.md` (arrived in `~/Downloads` as `cc_task_oc_breadth_stage2_20260831.md`, copied as `cc_task_oc_breadth_stage2_2026-08-31.md`, committed alone as `b2d2a2de`).
**Forecast scored:** `REPORT_oc_breadth_ladder_2026-08-30.md`, commit **`c9cb0ca2`**, read from `oc_breadth_ladder_2026-08-30_forecast120.rds` / `_gate.rds` and `oc_wrapper_grid_corrected_2026-08-30.rds`.
**Writes:** one new driver `quarto/simulations/actg175/continuous/sim_fs_maxeffCons_mr_md120_knoise0_n500_batch_1_1000.qmd` and its two payloads under `mr_md_harm/`; under `dev/glm-continuous-sims/`: `oc_breadth_stage2_gate_2026-08-31.R` (+ `_gate.log`, `_gate.rds`), `oc_breadth_stage2_score_2026-08-31.R` (+ `_score.log`, `_score.rds`), `oc_breadth_stage2_2026-08-31_driver.diff`, `oc_breadth_stage2_2026-08-31_driver_purl.R`, `_run1.log`, `_run2.log`, this report.
**No `R/` change. No null run (§2(d)). No push, no install, no render.**

## 0. The forecast as read (6 decimals), and provenance

Read from `oc_breadth_ladder_2026-08-30_forecast120.rds` (commit `c9cb0ca2`), `gate.log` lines `[2a]`:

- forecast rung `q = 120`; **`k_inter = −93.744764`** (= `k40` −13.744764 + s·(120 − 40), s = −1; stored to full precision as −93.7447641240); **`c1* = 135.741162`** (135.7411624608); `c1_05 = 133.234686`; `c1_10 = 125.700212`.

| quantity at (n = 500, c1*, c2 = 10, pconsistency 0.9) | resample | split | null at c1* (resample) |
|---|---:|---:|---:|
| det_rate | **0.800000** | 0.799205 | **0.038675** |
| det_rate_se (MC) | 0.000894 | 0.000896 | 0.000431 |
| EnH | 72.835590 | 72.814450 | 66.554890 |
| Esens | 0.332644 | 0.332651 | — |
| Espec | 0.952637 | 0.952705 | 0.866890 |
| Eppv | 0.775219 | 0.775531 | 0.000000 |
| Enpv | 0.731992 | 0.732007 | 1.000000 |
| EbetaH | 98.927990 | 98.957240 | 26.255240 |
| Enaive_bias | 62.263700 | 62.241820 | 118.410300 |
| mass_below | 1.000000 | 1.000000 | 1.000000 |
| sel_c mass on PQg ≥ 0.95 | 0.344819 | 0.344861 | — |

At the driver's `c1 = 30` (resample): det 1.000000, EnH 73.577948, Esens 0.327692, Espec 0.947767, Eppv 0.755575, EbetaH 97.086444, naive 56.918245. The forecast rung's resample grid (`$grid$resample$table`, `c1 = 0..300`) is the predicted alternative curve (det_rate 1.000000 at 30, 0.837190 at 133, 0.796300 at 136). The null curve artefacts in the locked record are the corrected null sweep (resample, `c1 = 20..120` by 5) and the null inversions at targets 0.05/0.10/0.20/0.50/0.80/0.90/0.95 → `c1` 133.234686 / 125.700212 / 117.113793 / 101.640168 / 87.292147 / 80.167269 / 74.433141, plus the forecast's null point at `c1*`. No grid was re-run; no null prediction exists in the record above `c1 = 120` other than those points.

**§1 provenance.** `pop-os`, `/home/larryleon/Documents/GitHub/forestsearch`, `feature/glm-extension`, `8ff85f9c`, `git status --porcelain` empty, log `8ff85f9c c9cb0ca2 f5388486 2276ee10 fbd564de 911d73ba`; `git show --stat c9cb0ca2` names `REPORT_oc_breadth_ladder_2026-08-30.md`; `packageVersion("forestsearch") = 0.3.1`. PASS. Task doc committed alone as `b2d2a2de`.

## 2. Source gates

**(b) The MD40 alternative n = 500 driver** (`sim_fs_maxeffCons_mr_md40_knoise0_n500_batch_1_1000.qmd`):

- DGM (L320–336): `calibrate_glm_interaction(data = actg_df, factor_vars = z1..z12, outcome_var = "cd4_change", treatment_var = "treat", target_effect = target_md_harm (−40), outcome_type = "continuous", effect_measure = "MD", subgroup_vars = c("z1","z2"), subgroup_cuts = list(z1 = 1L, z2 = 1L), k_inter_range = c(0, 120), grid_step = 2, n_super = 5000L, seed = seed_base, verbose = FALSE)` — the call Stage 1 §2 reproduced; the null branch is `generate_glm_dgm(<same>, model = "null")`.
- Seeds (L119–123, L536–539): `seed_base <- 8316951L`; per replicate `sd_i <- seed_base + sim_id`; `RNGkind("L'Ecuyer-CMRG"); set.seed(sd_i)`; `simulate_from_glm_dgm(dgm, n = n_sample, seed = sd_i)`; `seedit = sd_i` to the search. Replicate count `n_sims <- 1000L`, `sim_id_start <- 1L`.
- `forestsearch_args` (L560–590): `effect.threshold = md_threshold` (**30**), `consistency.threshold = md_consistency` (**10**), `pconsistency.threshold = 0.90`, `fs.splits = 400L`, `n.min = 60L`, `d0.min = d1.min = 12L`, `maxk = 2L`, `vi.grf.min = -0.2`, `sg_focus = "maxeffCons"`, `selection_rule = "neighborhood"`, `effect_neighborhood = 0.10`, `stop_threshold = NULL`, `consistency_method = "resample"`, `conf.cont_jcuts = list(age = 10, preanti = 10)`, `use_lasso/dina/grf = FALSE`, `use_twostage = TRUE`, `is.RCT = TRUE`, `adverse_outcome = FALSE`, `seedit = sd_i`, **`parallel_args = inner_parallel` = `list(plan = "sequential")`** under `parallel_mode = "sims"`, `mr_inference = TRUE`, `mr_inference_args = list(ci_method = "ij", draws = 5000L, include_complement = TRUE)`; `confounders.name` = the 13-vector `age preanti wtkg karnof cd40 cd80 hemo homo drugs race gender symptom str2`.
- Outer plan (L107, L822–827): `n_workers <- ceiling(0.90 * (detectCores(logical = FALSE) − 1))` = **115** on this 128-core host; `registerDoFuture(); plan(multisession, workers = n_workers)`; `foreach(...) %dofuture% record_replicate(i)`.
- **The recorder** (`.na_record`, L454–487; `record_replicate`, L534–720) stores per replicate: `status`, `detected`, `n_harm` (= `sum(grp.consistency$sg.harm.id == 1)`), `n_true`, **`sg_def`** (= `paste(fs.est$sg.harm, collapse = " & ")`), `sens/spec/ppv/npv` on the analysis sample, **`nv_H_est`** (= `fs.est$mr_inference$naive$est`, oriented), `nv_H_lo/hi/se`, `mr_H_*`, `nv_Hc_*`, `mr_Hc_*`, oracle columns, `fit_mr_secs`, and after the loop `fs_attach_betaHhat()` adds the exact `betaHhat_H/Hc` on `df_super`; `.payload$oc <- fs_mr_oc_summary(.payload)` is the OC block. Membership is not stored as a vector; it is reconstructible from `sg_def` and the seed.
- *GATE:* the winner's fitted effect and its definition are stored. `nv_H_est` **is the search's own fitted effect for the winner**: on replicate 1 the consistency table's row-1 `hr` = 144.373062662 and `nv_H_est` = 144.3730626617 (difference 2.8e-14; replicate 2: 3.1e-13 — the MR naive estimate refits the same unadjusted `lm(y ~ treat)` on the same subgroup, `fs_mr_inference.R:484`, `consistency_resample.R:240ff`). **PASS.**

**(c) Where `effect.threshold` enters the search** (every non-comment site on this path, `grep -rn` over `R/`):

1. `forestsearch_main.R:1314–1316` — alias: `if (!is.null(effect.threshold)) hr.threshold <- effect.threshold`; `:1736` `effect_threshold <- hr.threshold` (identity scale for MD, no remap when user-set, `:1790–1800`).
2. `forestsearch_main.R:2949` — `search_overrides$hr.threshold = effect_threshold`; `:2970` `search_overrides$effect_threshold <- effect_threshold`; `:2958` `disable_effect_floor = !.admit_applies[["effect"]]` (FALSE for `maxeffCons`; the floor is disabled only for `sg_focus = "maxeff"`, `:1487–1530`).
3. `subgroup_search.R:126–129` — `screen_threshold <- effect_threshold`; `:156` passed as `hr.threshold` to `search_combinations_parallel`; **`:623–628` status 6: `if (!disable_effect_floor && glm_result$hr <= hr.threshold) return(list(status = 6L, result = NULL))`** — the candidate table is filtered, after the size floors (status 3/4: `d0.min/d1.min` skipped for continuous, `nx <= n.min`) and the fit (status 5), on the candidate's fitted effect, strict `>` to pass.
4. `forestsearch_main.R:3098` → `subgroup_consistency_main.R:545–548` — `found.hrs <- hr.subgroups[hr.subgroups$HR >= hr.threshold, ]` on entry to `subgroup.consistency()`: a re-filter of the same table on the same statistic (redundant with 3 up to the `>` / `>=` boundary), **before** the consistency loop; `hr.threshold` appears nowhere else in that file.
5. `forestsearch_main.R:1969–1973` — `.fs_resolve_admission(...)` carries `effect_floor = effect_threshold` into `admission_resolved`, consumed by MR (`fs_mr_inference.R:413–425`, `t_g <- pmax(effect_floor, c_cons + z·sdv)`; the per-draw re-selection domain) and by the DINA/GRF selectors (`forestsearch_helpers.R:1222`, `:1632` — not on this path). MR runs **after** the winner is fixed ("post-selection inference on the completed search; it cannot change the identified subgroup"); it affects `mr_H_*` (debiased) columns, not declaration, the winner, or `nv_H_est`.

Not touched by the threshold: the consistency evaluation itself (`hr.consistency`/`consistency_threshold` = 10 and `pconsistency` = 0.9), the early-stopping rule (`stop_threshold = NULL`, `stop_Kgroups`), the near-duplicate removal (statistics-keyed inside `subgroup.consistency`, no threshold argument), the size floors (`n.min`, `d0.min`, `d1.min`, `minp`, `rmin` — evaluated before status 6), the GRF screen (`vi.grf.min`), the two-stage pruning. **The expectation holds** for everything that determines declaration and the winner: raising `effect.threshold` from 30 to `c1*` removes candidates with fitted effect ≤ `c1*` from the table before consistency; the consistency loop, ranking (`maxeffCons`: highest effect meeting `pconsistency`) and stopping act on the reduced table, so the winner at `c1*` is the winner at 30 whenever the latter's effect is ≥ `c1*`, and no declaration at `c1*` occurs otherwise. The one refinement is item 5: the same value also enters MR's admission set, post hoc — irrelevant to §4's scoring, which uses `nv_H_est` and `sg_def`. Empirically (gate log): replicate 1 had 1395 candidates entering consistency at `c1 = 30` and 2 at `c1*`; replicate 2, 1248 and (declared) the same winner; both re-runs at `c1*` declared iff `nv_H_est ≥ c1*`, with identical `sg_def` and identical `sg.harm.id` vectors.

**(d) The null payload** (`fs_maxeffCons_mr_mdnull_knoise0_n500_res_1_1000.rds`, `seed_base 8316951`, `effect_threshold 30`, `n_workers 115`, pkg 0.2.2): 998 DETECTED / 2 NO-DETECTION; `nv_H_est` non-NA on 998 / 998 detected, `sg_def` and `n_harm` likewise. **The null is scored post hoc; no null run was made.**

## 3. The measured cell

**The driver.** `sim_fs_maxeffCons_mr_md120_knoise0_n500_batch_1_1000.qmd`, a copy of the MD40 file; `diff -u` (`oc_breadth_stage2_2026-08-31_driver.diff`, 4 hunks):

```
-target_md_harm <- -40                        +target_md_harm <- -120   (the scenario tag -> stem md120)
+stage2_forecast <- readRDS("../../../../dev/glm-continuous-sims/oc_breadth_ladder_2026-08-30_forecast120.rds")
+stopifnot(stage2_forecast$q == 120, abs(stage2_forecast$k_inter + 93.7447641240) < 1e-9)
+stage2_direct <- identical(Sys.getenv("FS_STAGE2_RUN"), "direct")
+if (isTRUE(stage2_direct)) rds_stem <- paste0(rds_stem, "_c1star")
-md_threshold   <- 30                         +md_threshold <- if (isTRUE(stage2_direct)) stage2_forecast$c1_star else 30
-  calibrate_glm_interaction(... target_effect = target_md_harm, k_inter_range = c(0,120), grid_step = 2 ...)
+  generate_glm_dgm(data = actg_df, factor_vars = dgm_factor_vars, outcome_var = "cd4_change", treatment_var = "treat",
+    outcome_type = "continuous", effect_measure = "MD", subgroup_vars = dgm_subgroup_vars, subgroup_cuts = dgm_subgroup_cuts,
+    model = "alt", k_treat = 1, k_inter = stage2_forecast$k_inter, adverse_outcome = FALSE, n_super = dgm_n_super, seed = seed_base, verbose = FALSE)
```
Everything else — the seed table, `n_sims = 1000`, the `forestsearch_args`, the recorder, `parallel_mode = "sims"` × 115 workers — is byte-identical (the diff shows nothing else). `k_inter` and `c1*` are read from the locked `.rds`, not typed; the driver's own `stopifnot(abs(truth$effect_Q − target_md_harm) < 1e-8)` held (`effect_Q = −120.0000000000`, log line 2 of both runs). The driver was executed **without a render**: `knitr::purl()` of the `.qmd`, truncated before the `summary-prep` chunk (`oc_breadth_stage2_2026-08-31_driver_purl.R`: the `setup-knobs`, `build-dgm`, `machinery`, `run-batch` chunks exactly as purled), run from the driver's directory with `Rscript` under the installed 0.3.1 (`meta$pkg_version = "0.3.1"` in both payloads), `FS_STAGE2_RUN=comparability` / `=direct`.

**DGM gates (gate log `[3]`):** `fs_dgm_scale`: `|m_tau[Q]| = 120.000000000000` (gate 1e-9 TRUE), `|m_tau[Qc]| = 26.255235876036` = MD40's (TRUE); family enumerated at n = 500, M = 1696, `identical()` to the stored forecast-rung family on `lab, Pg, PQg, se_g, sens_g, spec_g, M` (gate.rds record) and on `beta_g` (forecast120.rds): all TRUE. `df_super` differs from MD40's in exactly one column, `mu1`, by exactly −80 on Q and 0 off Q (`mu0` identical). PASS.

**Compute projection (stated before launch).** MD40's own record: 1000 reps in **455 s** on 115 workers (fit+MR 26.3 s/rep at 0.2.2, overhead k = 1.92); the gate script measured 16.3 s (c1 = 30) / 13.6 s (c1*) per call at 0.3.1. Projection: ceil(1000/115) × 16.3 × 1.92 ≈ 4.7 min per run; two runs ≤ 15 min even at MD40's 455 s each. Below 6 h → launched. **Measured:** run 1 (`effect.threshold = 30`) 22:49:25 → 22:57:14, **469 s** (fit+MR mean 30.0 s/rep, max 49 — the gate script was re-run concurrently, single-threaded); run 2 (`effect.threshold = c1*`) 22:57:40 → 23:03:54, **374 s** (16.5 s/rep, max 27.8). Sequential, both detached (`nohup`), no worker error, 1000 rows each.

**Paired replicates.** `simulate_from_glm_dgm()` (`R/simulate_from_glm_dgm.R:62ff`) does `set.seed(seed)`, then `sample()` the rows, then `rbinom()` the treatment, then `y = mu + rnorm(n, sd = sigma)`; nothing between depends on the DGM's `mu1`. So the RNG stream between sampling and outcome generation is consumed identically and **the residual draws are aligned**: on replicates 1–5 the sampled rows and `treat_sim` are `identical()` between the MD40 and MD120 DGMs and `y_sim` differs by exactly −80 on Q-treated subjects and exactly 0 everywhere else (gate log `[3]`). The two cells differ per replicate only by the shift of Q's treated outcomes. (The null payload shares the same rows and treatment too — its replicate 1–3 winners in the payload are the same strings as MD40's.)

## 4. Scoring (`oc_breadth_stage2_score_2026-08-31.R`, `_score.log`)

Run 1: 1000 / 1000 DETECTED at `c1 = 30`; `T̂ = nv_H_est` ranges 84.10–226.54; `declared_c1* = declared_30 & (T̂ ≥ c1*)`; no `T̂` equals `c1*` exactly. Population functionals of each realized rule are evaluated on `df_super` with `.fs_resolve_membership()` (sgdef §3 machinery; all rules resolved `ok`); the payload's exact `betaHhat_H` equals `tauQc + bint·PQg` on every declared replicate (< 1e-8).

**4.1 Power at `c1*`.** **786 / 1000 = 0.7860 (SE 0.0130)** against the forecast **0.8000 (MC SE 0.0009)**; gap −0.0140 = −1.08 combined SE.

**4.2 The measured declaration curve** (`rate(declared_30 & T̂ ≥ c1)`, `c1 = 0..300`, beside the rung's predicted `det_rate`; the figure lives in `score.rds$curve`):

| c1 | measured | predicted | gap | c1 | measured | predicted | gap |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0–60 | 1.000 | 1.0000 | 0 | 140 | 0.720 | 0.7337 | −0.0137 |
| 100 | 0.996 | 0.9977 | −0.0017 | 145 | 0.626 | 0.6479 | −0.0219 |
| 110 | 0.982 | 0.9874 | −0.0054 | 150 | 0.536 | 0.5554 | −0.0194 |
| 120 | 0.951 | 0.9535 | −0.0025 | 160 | 0.375 | 0.3708 | +0.0042 |
| 125 | 0.916 | 0.9201 | −0.0041 | 170 | 0.220 | 0.2206 | −0.0006 |
| 130 | 0.853 | 0.8723 | −0.0193 | 180 | 0.117 | 0.1168 | +0.0002 |
| 133 | 0.818 | 0.8372 | −0.0192 | 200 | 0.023 | 0.0233 | −0.0003 |
| 135 | 0.798 | 0.8104 | −0.0124 | 220 | 0.001 | 0.0032 | −0.0022 |
| 136 | 0.781 | 0.7963 | −0.0153 | 250–300 | 0.000 | ≤ 0.0001 | 0 |

Largest absolute gap **0.0246 at `c1 = 146`** (measured 0.605, predicted 0.630, −1.59 combined SE); |gap| > 0.02 only on `c1 ∈ [145, 151]`. The measured curve sits slightly below the predicted one over 130–155 (−0.012 to −0.025, each within 1.7 SE) and on it elsewhere.

**4.3 The direct run as the gate.** Run 2 declared **786 / 1000 = 0.7860 (SE 0.0130)** — the same number as run 1's post-hoc rate. Per replicate: run 2 declares iff `declared_c1*` from run 1 on **1000 / 1000** (0 disagreements); among the 786 joint declarations `sg_def` identical 786 / 786, reconstructed membership identical 786 / 786, `n_harm` identical 786, `T̂` identical to 1e-9 786. *GATE (§2(c) said the expectation holds): agreement is complete.* **PASS.** Run 1's post-hoc scoring is the measured cell; run 2 confirms it exactly.

**4.4 Type-I at `c1*`.** Null payload: **40 / 1000 = 0.0400 (SE 0.0062)** against **0.03868 (MC SE 0.00043)**; gap +0.0013 = 0.21 combined SE. The null's measured curve beside the record's null predictions:

| c1 (null inversion / c1*) | predicted | measured | SE | gap (SE) |
|---:|---:|---:|---:|---:|
| 133.235 (target 0.05) | 0.0500 | 0.049 | 0.0068 | −0.001 (−0.15) |
| 125.700 (0.10) | 0.1000 | 0.095 | 0.0093 | −0.005 (−0.54) |
| 117.114 (0.20) | 0.2000 | 0.182 | 0.0122 | −0.018 (−1.47) |
| 101.640 (0.50) | 0.5000 | 0.462 | 0.0158 | −0.038 (−2.40) |
| 87.292 (0.80) | 0.8000 | 0.773 | 0.0132 | −0.027 (−2.03) |
| 80.167 (0.90) | 0.9000 | 0.882 | 0.0102 | −0.018 (−1.76) |
| 74.433 (0.95) | 0.9500 | 0.939 | 0.0076 | −0.011 (−1.45) |
| **135.741 (c1\*)** | **0.0387** | **0.040** | 0.0062 | +0.001 (+0.21) |

Against the corrected null sweep (`c1 = 20..120` by 5): measured − predicted is +0.001 to +0.003 on 20–70 (measured 0.998 vs 0.9969), then −0.016 (75), −0.019 (80), −0.032 (85), −0.033 (90), **−0.045 (95, the largest)**, −0.038 (100), −0.040 (105), −0.028 (110), −0.019 (115), −0.010 (120). The null's measured curve runs 2–4.5 points *below* the prediction in the middle of the `c1` range (2–3 SE at 95–105) and rejoins it above 115; at the thresholds that matter here (`c1_05`, `c1*`) it is on the prediction.

**4.5 The declared rule, conditional on `declared_c1*` (n = 786):**

| quantity | measured (sample) | SE | measured (population of realized rules) | SE | forecast resample | forecast split |
|---|---:|---:|---:|---:|---:|---:|
| E\|Ĥ\| | **74.963** | 0.560 | **75.294** | 0.610 | 72.836 | 72.814 |
| PPV | 0.7870 | 0.0076 | 0.7740 | 0.0078 | 0.7752 | 0.7755 |
| sensitivity | 0.3473 | 0.0049 | 0.3429 | 0.0050 | 0.3326 | 0.3327 |
| specificity | 0.9535 | 0.0017 | 0.9505 | 0.0018 | 0.9526 | 0.9527 |
| E[β(Ĥ)] (population true effect; truth in Q 120) | — | — | **98.815** | 0.730 | 98.928 | 98.957 |
| naive bias `T̂ − β(Ĥ)` | **62.269** | 0.881 | — | — | 62.264 | 62.242 |

Population-versus-sample offsets (sample − population, paired): size −0.33 (SE 0.29), PPV +0.0129 (0.0015), sensitivity +0.0043 (0.0012), specificity +0.0029 (0.0004) — the expected direction (sample proportions at a selected rule biased upward), the same size as at MD40 (+0.005 / +0.001 / +0.002 there; PPV's offset is larger here).

**Limitation, carried verbatim from the forecast:** *the population-versus-sample offset on PPV and sensitivity is expected (the driver's sample proportions at a selected rule are biased upward), and at MD40 `E|Ĥ|`, PPV, sensitivity and the naive bias inherit a between-rule gap of +2.11 subjects that five mechanisms failed to explain.*

**The between-rule gap at this harm** (confs-compare §3 definition: measured-frequency-weighted population size of the realized rules minus the analytic `sel_c`-weighted `E|Ĥ|`):

| gate | analytic E\|Ĥ\| | pop. size of realized rules | **between-rule** | within-rule (sample − pop) | measured E\|Ĥ\| | MD40 between-rule | `sel_c`-reweighted check |
|---|---:|---:|---:|---:|---:|---:|---:|
| resample | 72.836 | 75.294 | **+2.459** (SE 0.61) | −0.331 | 74.963 | +2.111 | 75.258 |
| split | 72.814 | 75.294 | **+2.480** | −0.331 | 74.963 | +2.116 | 75.259 |

The gap at q = 120 is **+2.46 subjects against +2.11 at MD40** — larger by 0.35 (0.6 SE of the q = 120 estimate): within noise of "held", with the point estimate on the "grew" side. As at MD40, it is between-rule (the realized rules are larger in the population than the analytic weights predict; the reweighted check — analytic weights on the measured signatures' sizes — reproduces the measured 75.26, so it is the weights over signatures, not the rules within them). No mechanism is proposed.

**4.6 Composition.** Measured selection mass on rules with `PQg ≥ 0.95`: **294 / 786 = 0.3740 (SE 0.0173)** against the forecast **0.3448** (+0.029, 1.7 SE); mass on `PQg = 1`: 0.360. Q itself (`age > 34 & preanti ≤ 744.5`) was selected 0 times — it is not a family member (Stage 1), and no realized rule reproduces it. Top 15 realized rules (verbatim; population `Pg`, `PQg`, `β_g`; closest analytic label by Jaccard on `df_super`):

| realized rule | n | freq | Pg | PQg | β_g | closest analytic (Jaccard) |
|---|---:|---:|---:|---:|---:|---|
| !{age ≤ 37} & {preanti ≤ 0} | 13 | 0.0165 | 0.138 | 1.000 | 120.0 | age > 37 & preanti ≤ 0 (1.000) |
| !{age ≤ 37} & !{str2} | 12 | 0.0153 | 0.139 | 1.000 | 120.0 | age > 37 & str2 != 1 (1.000) |
| !{age ≤ 35} & !{str2} | 8 | 0.0102 | 0.171 | 1.000 | 120.0 | age > 35 & str2 != 1 (1.000) |
| !{age ≤ 35} & {age ≤ 39} | 8 | 0.0102 | 0.171 | 0.744 | 96.0 | age > 35 & age ≤ 39 (1.000) |
| !{age ≤ 33} & {preanti ≤ 0} | 7 | 0.0089 | 0.203 | 0.929 | 113.4 | age > 33 & preanti ≤ 0 (1.000) |
| !{age ≤ 34} & {preanti ≤ 0} | 7 | 0.0089 | 0.189 | 1.000 | 120.0 | age > 33 & preanti ≤ 0 (0.929) |
| !{age ≤ 38} & !{str2} | 7 | 0.0089 | 0.121 | 1.000 | 120.0 | age > 37 & str2 != 1 (0.869) |
| !{age ≤ 33} & !{homo} | 5 | 0.0064 | 0.151 | 0.684 | 90.4 | age > 33 & homo != 1 (1.000) |
| !{age ≤ 35} & !{homo} | 5 | 0.0064 | 0.126 | 0.738 | 95.5 | age > 35 & homo != 1 (1.000) |
| !{age ≤ 36} & {preanti ≤ 0} | 5 | 0.0064 | 0.158 | 1.000 | 120.0 | age > 35 & preanti ≤ 0 (0.932) |
| !{age ≤ 38} & {karnof ≤ 90} | 5 | 0.0064 | 0.147 | 0.729 | 94.6 | age > 37 & karnof ≤ 95 (0.910) |
| {preanti ≤ 0} & !{age ≤ 37} | 5 | 0.0064 | 0.138 | 1.000 | 120.0 | age > 37 & preanti ≤ 0 (1.000) |
| !{age ≤ 34} & !{str2} | 4 | 0.0051 | 0.191 | 1.000 | 120.0 | age > 33 & str2 != 1 (0.930) |
| !{age ≤ 35} & {preanti ≤ 0} | 4 | 0.0051 | 0.170 | 1.000 | 120.0 | age > 35 & preanti ≤ 0 (1.000) |
| !{age ≤ 36} & !{homo} | 4 | 0.0051 | 0.116 | 0.734 | 95.0 | age > 35 & homo != 1 (0.917) |

Realized age cuts at 34, 36, 38 have no analytic twin (the population grid has 33/35/37/39) — the sample-quantile-cut variation the sgdef/residual reports recorded. Top 15 analytic by `sel_c` (resample), with the measured frequency of the candidate's signature:

| analytic candidate | sel_c | Pg | PQg | β_g | measured freq of its signature |
|---|---:|---:|---:|---:|---:|
| age > 39 & preanti ≤ 188 | 0.0795 | 0.125 | 1.000 | 120.0 | 0.384 (age> & preanti≤, all) |
| age > 42 & preanti ≤ 875 | 0.0596 | 0.127 | 0.926 | 113.1 | " |
| age > 39 & preanti ≤ 402 | 0.0382 | 0.150 | 1.000 | 120.0 | " |
| age > 37 & preanti ≤ 0 | 0.0363 | 0.138 | 1.000 | 120.0 | " |
| age > 37 & preanti ≤ 39 | 0.0320 | 0.146 | 1.000 | 120.0 | " |
| age > 35 & homo != 1 | 0.0287 | 0.126 | 0.738 | 95.5 | 0.027 |
| age > 37 & str2 != 1 | 0.0277 | 0.139 | 1.000 | 120.0 | 0.050 |
| age > 33 & cd80 ≤ 645 | 0.0270 | 0.124 | 0.715 | 93.3 | 0.051 |
| age > 35 & cd40 ≤ 260 | 0.0269 | 0.124 | 0.723 | 94.1 | 0.045 |
| age > 35 & wtkg > 83 | 0.0260 | 0.120 | 0.724 | 94.1 | 0.041 |
| age > 39 & preanti ≤ 672 | 0.0233 | 0.175 | 1.000 | 120.0 | 0.384 |
| age > 33 & cd40 > 420 | 0.0210 | 0.122 | 0.691 | 91.1 | 0.041 |
| age > 35 & preanti ≤ 0 | 0.0208 | 0.170 | 1.000 | 120.0 | 0.384 |
| age > 35 & age ≤ 39 | 0.0194 | 0.171 | 0.744 | 96.0 | 0.038 |
| age > 33 & race | 0.0187 | 0.133 | 0.691 | 91.1 | 0.023 |

Over signatures: `age> & preanti≤` analytic 0.378 / measured 0.384; `age> & cd80≤` 0.058 / 0.051; `age> & cd40>` 0.055 / 0.041; `age> & cd40≤` 0.053 / 0.045; `age> & str2=0` 0.049 / 0.050; `age> & wtkg>` 0.045 / 0.041; `age> & cd80>` 0.037 / **0.060**; `age<= & age>` 0.034 / 0.038; `age> & homo=0` 0.034 / 0.027; `age> & wtkg≤` 0.024 / 0.033. **TVD over signatures 0.1216 (split 0.1220)**; multinomial floor at 786 declarations 0.0946 (2.5–97.5 %: 0.076–0.116); **excess 0.027** (MD40 corrected: 0.045 at n = 500). Analytic mass on never-selected signatures 0.047 (173 analytic signatures, 79 measured); measured mass absent from the family 0.0013 (1 replicate).

**4.7 At the driver's `c1 = 30`** (saturated; measured columns by the corrected report's `measured_from_payload`, MD40 columns from `oc_wrapper_grid_corrected_2026-08-30.rds`):

| quantity | MD120 measured (SE) | MD120 forecast resample / split | MD40 measured (SE) | MD40 analytic |
|---|---:|---:|---:|---:|
| det_rate | 1.0000 | 1.0000 / 1.0000 | 1.0000 | 0.9991 |
| E\|Ĥ\| | **75.571 (0.515)** | 73.578 / 73.560 | 72.335 (0.415) | 70.786 |
| PPV (sample) | 0.7616 (0.0072) | 0.7556 / 0.7558 | 0.4055 (0.0080) | 0.4000 |
| sensitivity (sample) | 0.3390 (0.0045) | 0.3277 / 0.3277 | 0.1706 (0.0035) | 0.1638 |
| specificity (sample) | 0.9473 (0.0016) | 0.9478 / 0.9478 | 0.8686 (0.0019) | 0.8701 |
| E[β(Ĥ)] | 96.443 (0.684) | 97.086 / 97.108 | 31.765 (0.108) | 31.754 |
| naive bias | 56.902 (0.846) | 56.918 / 56.864 | 74.273 (0.577) | 75.529 |

Population rows at `c1 = 30` (n = 1000): pop. size 76.18 (between-rule **+2.60**, within −0.61), PPV_pop 0.7487, sens_pop 0.3357, spec_pop 0.9441. The saturated picture at this harm has the same shape as MD40's: declaration, specificity, E[β(Ĥ)] and the naive bias on the prediction; E|Ĥ| +2.0 (sample) / +2.6 (population) above it; sensitivity +0.011 (2.5 SE).

## 5. Verdict

Combined SE = √(measured SE² + forecast MC SE²); "within noise" = |measured − predicted| ≤ 2 combined SE.

1. **Power at `c1*`: 0.786 ± 0.013 vs 0.800 — within noise** (−0.014, 1.1 SE). Run 2 measured the same 0.786 directly.
2. **Type-I at `c1*`: 0.040 ± 0.006 vs 0.0387 — within noise** (+0.001, 0.2 SE).
3. **E[β(Ĥ)]: 98.81 ± 0.73 vs 98.93 — within noise** (−0.11). **Specificity: 0.9505 (population) / 0.9535 (sample) vs 0.9526 — within noise** (−0.002 / +0.001).
4. **PPV: 0.7740 (population) vs 0.7752 — within noise** (−0.001); sample 0.7870, +0.012 — **population-versus-sample** (the paired offset is +0.013 ± 0.002). **Sensitivity: 0.3429 (population) vs 0.3326 — a gap of +0.010** (2.1 SE, at the boundary); sample 0.3473, +0.015 — population-versus-sample on top of it. **E|Ĥ|: 74.96 (sample) / 75.29 (population) vs 72.84 — a gap of +2.1 / +2.5 subjects** (3.8 / 4.0 SE). **Naive bias: 62.27 ± 0.88 vs 62.26 — within noise** (+0.006). **Between-rule gap at this harm: +2.46 ± 0.61 against +2.11 — within noise of MD40's value** (+0.35, 0.6 SE); it held, with the point estimate on the larger side.
5. **The declaration curves.** Alternative: on the prediction over 0–300, largest gap −0.025 at `c1 = 146`, nowhere beyond 1.7 SE — within noise everywhere. Null: on the prediction at `c1_05` and `c1*` and above 115; **a gap of −0.03 to −0.045 over `c1 = 85–105`** (2–3 SE at 95–105: the null's measured curve falls faster than predicted through the middle of its range, then rejoins).

**Does the wrapper predict the crossover regime as well as it predicted the saturated one?** Yes, on the two quantities that define it: at a threshold never run at the search — `c1* = 135.74`, chosen analytically from a 2×10⁵-draw prediction — the measured power is 0.786 against a forecast of 0.800 and the measured null false-declaration rate is 0.040 against 0.0387, both within one to one-and-a-quarter binomial SE, and the direct run at that threshold reproduces the post-hoc scoring replicate for replicate. The conditional operating characteristics behave exactly as they did in the saturated regime: E[β(Ĥ)], specificity, the naive bias and PPV (population) are on the forecast; E|Ĥ| carries the same between-rule gap (+2.5 against +2.1 subjects, within noise of each other), and sensitivity inherits a tenth of it. The whole measured alternative curve tracks the predicted one to within 0.025 across `c1 = 0..300`, and the composition over signatures is closer to the analytic selection distribution than at MD40 (TVD excess 0.027 against 0.045). The one feature the saturated cell could not show is the shape of the null curve between `c1 = 85` and 105, where the measurement runs 3–4.5 points under the prediction before meeting it again at the 5 % threshold. Recorded; no recommendation, no mechanism, no next task.

## 6. Ten-line summary

1. Task committed alone as `b2d2a2de` (arrived hyphen-stripped as `cc_task_oc_breadth_stage2_20260831.md`); §1 gate passed on 0.3.1 with `c9cb0ca2` naming the ladder report.
2. Forecast read from the locked `.rds`, not typed: `k_inter = −93.7447641240`, `c1* = 135.7411624608`; power 0.800000 ± 0.000894, null 0.038675 ± 0.000431, E|Ĥ| 72.836, PPV 0.7752, sens 0.3326, spec 0.9526, E[β(Ĥ)] 98.928, naive 62.264.
3. Source: `effect.threshold` filters the candidate table (status 6, `subgroup_search.R:626`) and is re-applied on entry to `subgroup.consistency()` (`:545–548`); otherwise only MR's post-selection admission set — not the consistency loop, stopping, dedup, floors or GRF. Expectation holds; recorder stores `nv_H_est` (= the winner's `hr` to 3e-13) and `sg_def`; null payload has both → no null run.
4. Driver `sim_fs_maxeffCons_mr_md120_…qmd`: three-hunk diff (direct `generate_glm_dgm` at the read `k_inter`, stem `md120`/`_c1star`, run-2 threshold); DGM gates exact (|m_tau[Q]| = 120.000000000000, family identical, `df_super` differs only in `mu1` by −80 on Q); replicates paired with MD40's down to the residual draws.
5. Projection ≈ 5 min/run against a 6 h gate; measured 469 s + 374 s, two runs of 1000, no failures, executed by purl + `Rscript` (no render).
6. **Power at c1\*: 786/1000 = 0.786 ± 0.013 vs 0.800 — within noise**; run 2 declared 786 directly, agreeing with post-hoc scoring on 1000/1000 replicates with identical winners on all 786.
7. **Type-I at c1\*: 40/1000 = 0.040 ± 0.006 vs 0.0387 — within noise.**
8. At c1\*: E[β(Ĥ)] 98.8 vs 98.9, spec 0.951 vs 0.953, PPV (pop) 0.774 vs 0.775, naive 62.27 vs 62.26 — within noise; sensitivity +0.010 (2.1 SE); **E|Ĥ| +2.1 (sample) / +2.5 (pop)** — the between-rule gap is **+2.46 against +2.11** (held, 0.6 SE); PQg ≥ 0.95 mass 0.374 vs 0.345; TVD excess over the multinomial floor 0.027.
9. Curves: alternative within 0.025 of the prediction on 0..300; null on the prediction at `c1_05`/c1\* but 0.03–0.045 under it over `c1 = 85–105`.
10. Committed as `6504e0ea` with the driver, both payloads (232 K / 204 K, tracked as the MD40 payloads are), scripts, logs and `.rds`; this line was added in the immediately following commit.
