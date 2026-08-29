# REPORT — the confounder-list mismatch settled, the analytic family corrected, and whether the correction closes the residuals

**Task:** `dev/tasks/cc_task_oc_wrapper_confs_2026-08-30.md` (commit `7e4a5dfd`)
**Repo / branch:** `~/Documents/GitHub/forestsearch` (pop-os), `feature/glm-extension`
**Follows and builds on:** `REPORT_oc_wrapper_sgdef_2026-08-29.md` (commits `69ccc84d`, `b8bd0f0a`, `f2ca69e6`) — that report **stands**; its tabulation, population evaluation and within-rule comparison are read from `sgdef_selection_2026-08-29.rds`, not redone.
**Category:** no `R/` change, no export, no edit to any package file, driver or document. Artefacts, all under `dev/glm-continuous-sims/`: `oc_wrapper_grid_corrected_2026-08-30.R` (+ four `.log`: `_alt500`, `_alt700`, `_null500`, `_merge`), `oc_wrapper_grid_corrected_2026-08-30.rds`, `oc_wrapper_confs_compare_2026-08-30.R` (+ `.log`, `.rds`), this report. The 08-29 scripts and `.rds` are untouched. No simulations, no replicates, no renders, no push.

---

## 0. Provenance (§1 — GATE passed)

```
pop-os · /home/larryleon/Documents/GitHub/forestsearch · feature/glm-extension
HEAD at start f2ca69e6 (git merge-base --is-ancestor f2ca69e6 HEAD: yes)
git status --porcelain: empty
git log -6: f2ca69e6 b8bd0f0a 69ccc84d 648c2af8 e3725912 6574beb5
packageVersion forestsearch 0.2.6
```

`ls ~/Downloads/cc_task_oc_wrapper_confs_2026*`: **exactly one match**, `cc_task_oc_wrapper_confs_20260830.md` (hyphens stripped) → copied to `dev/tasks/cc_task_oc_wrapper_confs_2026-08-30.md`, committed **`7e4a5dfd`**.

---

## 1. §2 — which covariate set did the search actually run?

Established from the driver source at the payloads' own commits and from the tracked payloads; every piece quoted.

### 1.1 The drivers' argument

The alt n = 500 and n = 700 drivers are `quarto/simulations/actg175/continuous/sim_fs_maxeffCons_mr_md40_knoise0_n{500,700}_batch_1_1000.qmd`; `diff` between them is **one line** (`n_sample <- 500L` vs `700L`, L133). Both payloads (`fs_maxeffCons_mr_md40_knoise0_n{500,700}_res_1_1000.rds`) were committed in `bb75cca6`; at that commit and at HEAD the two files read identically on every line below.

```r
k_random_noise <- 0                                                            # L116
analysis_continuous_vars <- c("age","preanti","wtkg","karnof","cd40","cd80")   # L242
analysis_binary_vars     <- c("hemo","homo","drugs","race","gender","symptom") # L243
include_str2 <- TRUE                                                           # L267
if (isTRUE(include_str2))                                                      # L268
  analysis_binary_vars <- c(analysis_binary_vars, "str2")                      # L269
for (v in analysis_binary_vars) actg_df[[v]] <- as.factor(actg_df[[v]])        # L301
confounders_analysis <- c(analysis_continuous_vars, analysis_binary_vars)      # L302
...
confs <- c(confounders_analysis,                                               # L550
           if (k_random_noise > 0) paste0("noise", seq_len(k_random_noise)))   # L551
fs.est <- suppressWarnings(forestsearch(                                       # L561
  df.analysis = df, confounders.name = confs,                                  # L562
```

So `confounders.name` = **`c("age","preanti","wtkg","karnof","cd40","cd80","hemo","homo","drugs","race","gender","symptom","str2")`, length 13**, in that order, for both alt drivers.

**The null cell.** There is no separate null driver: the null payload (`fs_maxeffCons_mr_mdnull_knoise0_n500_res_1_1000.rds`, committed `2b180813`) comes from the **same n = 500 file** with `null_cell` flipped (L139 `null_cell <- FALSE` is the committed text; L158 `if (isTRUE(null_cell)) rds_stem <- sub("_md[0-9]+_", "_mdnull_", rds_stem)` is what produces the `mdnull` stem; L320–326 select `generate_glm_dgm(model = "null")`). The flip was **not committed** — at `2b180813` the driver still reads `null_cell <- FALSE` — but at that commit `include_str2 <- TRUE` (L267) and `k_random_noise <- 0` (L116), and `null_cell` has no path into `analysis_binary_vars` / `confs`. The null run's `confounders.name` is therefore the same 13-vector *by the source*, and §1.4 confirms it *by the realized rules*. **The three drivers agree: 13 variables, `str2` last.**

### 1.2 Candidate versus evaluated — can `forestsearch()` reduce 13 → 12?

`R/forestsearch_main.R` L3351–3352:

```r
confounders.candidate = FSconfounders.name,   # = FSdata$confs_names  (L2738)
confounders.evaluated = confs_labels,         # = FSdata$confs        (L2739)
```

Both come out of `get_FSdata()`; `confounders.evaluated` is the vector of **cut labels** (`"age <= 35"`, …, `"str2"`), not a screened variable list. What could remove a variable between the 13 passed and the family searched, and whether it does here:

| mechanism | where | applies here? |
|---|---|---|
| GRF cut generation (`use_grf`) | `forestsearch_main.R` §3A L2462–2492 | `use_grf <- FALSE` (driver L214) — off |
| LASSO screen (`use_lasso`) | `get_fsdata.R` L362–395 | `use_lasso <- FALSE` (L214) — off |
| DINA (`use_dina`) | `forestsearch_main.R` L2590–2620 | `use_dina <- FALSE` (L214) — off |
| GRF variable-importance screen (`vi.grf.min`, `max_n_confounders`) | `forestsearch_main.R` §5 L2752–2831 | **runs** — the driver passes `vi.grf.min = -0.2` (L215/L568), not `NULL`, so `fit_causal_forest_glm()` is fitted per replicate; but the keep rule is `vi_ratio > vi.grf.min` with `vi_ratio = vi / max(vi) ∈ [0, 1]` (L2817–2819), which is **true for every column at −0.2**, and the truncation `seq_len(min(length, max_n_confounders))` is at the formal default **`max_n_confounders = 1000`** (L1167) against 74 cut columns. On GRF failure the fall-through is the full set (L2826–2831). **A no-op.** |
| continuous/categorical classification (`cont.cutoff`) | `get_fsdata.R` L273–277, `is.continuous(x, cutoff = 4)` = `length(unique(x)) >= 4` | `str2` has 2 unique values → categorical → enters as the bare cut `"str2"` (L486 `confs <- c(conf.categorical, conf.cont_Medcuts)`), exactly as `hemo` … `symptom` do |
| `exclude_cuts`, `conf_force`, `defaultcut_names`, `conf.cont_medians*` | `get_fsdata.R` | all `NULL` (not passed; formals default `NULL`) |
| `collapse_cuts` (near-redundant *continuous* cut collapse) | `get_fsdata.R` L13–27 | `TRUE` by default on both sides; acts on continuous cut grids, "categorical, indicator, bare-name" cuts are exempt (roxygen L18) — cannot drop `str2` |
| degenerate cuts | `fs_oc_family_enumerate()` `empty` stage; the search's own empty-membership guard | drops *rules* with empty membership, not variables |
| statistics-keyed near-duplicate removal, `rmin`, `minp`, `n.min` | the search / the enumeration | drop *rules*, not variables |

**A 13 → 12 drop is not legitimate here; nothing in the pipeline reduces the list.** `forestsearch()` evaluated all 13 on every replicate.

### 1.3 What the payloads say

`payload$meta` (all three payloads, same 18 fields): `n_sample, n_sims, nb_boots, mr_draws, sg_focus, selection_rule, consistency_method, subgroup_method, target_md_harm, effect_threshold, consistency_threshold, adverse_outcome, seed_base, sim_id_start, parallel_mode, n_workers, pkg_version, built_at`. **Neither `confounders.candidate` nor `confounders.evaluated` is recorded** — a recursive `str()` of each payload minus `$results` contains no `str2` and no `confound`. The only trace of the list is `results$sg_def` (and `results$covs`).

### 1.4 Evidence from the realized rules (from `sgdef_selection_2026-08-29.rds`, not recomputed)

`grepl("str2", reps$sg_def)` over detected replicates: **alt500 44 / 1000 (0.044); alt700 31 / 999 (0.031); null500 50 / 998 (0.050)** — the sgdef report's 0.044 / 0.031 / 0.050 `covariate_use` row exactly. A variable that appears in a realized winner was unambiguously in the evaluated family, in all three cells including the null.

### 1.5 What the analytic side used

`oc_wrapper_grid_2026-08-29.R` L60–66 (and, reading `G$fs_args` from it, `oc_wrapper_grid_sweep_2026-08-29.R` L48, `oc_wrapper_null_2026-08-29.R` L63, `sgdef_selection_2026-08-29.R`; and `fixture_run_fs_oc_predict_2026-08-28.R` per the fixture report's own description "`confounders.name` = the 12 analysis variables"):

```r
fs_args <- list(
  confounders.name = c("age", "preanti", "wtkg", "karnof", "cd40", "cd80",
                       "hemo", "homo", "drugs", "race", "gender", "symptom"),   # length 12
  conf.cont_jcuts = list(age = 10, preanti = 10),
  n.min = 60L, maxk = 2L, sg_focus = "maxeffCons",
  effect.threshold = 30, consistency.threshold = 10,
  pconsistency.threshold = 0.90)
```

The stored families' `cuts` (`G$families[["500"]]$cuts`, identical at 700 and at the null): 36 cut expressions — `age <= {25,28,30,31,33,35,37,39,42,47}`, `preanti <= {0,39,188,402,672,875,1055}`, `wtkg <= {75,67,83}`, `karnof <= {95,90}`, `cd40 <= {349,339,260,420}`, `cd80 <= {989,895,645,1210}`, `hemo`, `homo`, `drugs`, `race`, `gender`, `symptom` — **no `str2`**; `counts[["cut_columns"]] = 72` = 2 × 36. `str2` is on `dgm$df_super` (integer 0/1, `P(str2 = 0) = 0.4284`, `cor(str2, preanti) = 0.68`), so nothing prevented its enumeration; it was simply not in the list.

### 1.6 Every other enumeration-relevant argument, driver vs analytic reconstruction

| argument | driver (`forestsearch()` call, L561–581, or formal default) | analytic (`fs_args` → `fs_oc_family_enumerate()`, which reads `formals(forestsearch)` for anything not supplied, `R/fs_oc_family.R` L199–231) | |
|---|---|---|---|
| `confounders.name` | 13 (§1.1) | 12 | **DIFFER** |
| `conf.cont_jcuts` | `list(age = 10, preanti = 10)` (L213, L573) | same | agree |
| `cut_type` | not passed → `"default"` | not supplied → `"default"` | agree |
| `cont.cutoff` | not passed → `4` | not supplied → `4` | agree |
| `maxk` | `2L` (L176, L567) | `2L` | agree |
| `n.min` | `60L` (L176, L567) | `60L` (`args_used$n.min = 60`, floor 0.12 / 0.0857) | agree |
| `n.min.frac` | default `0.1`, inert (`n.min` supplied) | default, inert | agree |
| `minp` | not passed → `0.025`; kept (only `sg_focus = "maxeff"` zeroes it, L1469–1475) | `.arg("minp")` = `0.025`; same `maxeff`-only override (L220) | agree |
| `rmin` | not a formal; `subgroup.search` default `5` (L2922–2925: only `maxeff` sets 0) | `formals(subgroup.search)$rmin` = `5` (L216–217) | agree |
| `exclude_cuts` | not passed → `NULL` | `NULL` | agree |
| `conf_force` | not passed → `NULL` | `NULL` | agree |
| `collapse_cuts` / `collapse_cuts_args` | default `TRUE` / `list()` | default `TRUE` / `list()` | agree |
| `sg_focus` | `"maxeffCons"` (L127, L568) | `"maxeffCons"` | agree |
| `effect.threshold` | `30` (L171, L565); `meta$effect_threshold = 30` | `30` | agree |
| `consistency.threshold` | `10` (L173, L565); `meta$consistency_threshold = 10` | `10` | agree |
| `pconsistency.threshold` | `0.90` (L175, L566) | `0.90` | agree |
| `consistency_method` | `"resample"` (L165; `meta`) | both gates evaluated | (by design) |

Driver arguments with **no analytic counterpart** (recorded, not mismatches): `d0.min = d1.min = 12L`, `fs.splits = 400L` (unused on the `"resample"` path), `vi.grf.min = -0.2` (no-op, §1.2), `use_twostage = TRUE`, `selection_rule = "neighborhood"` / `effect_neighborhood = 0.10` (inert under `maxeffCons`), `stop_threshold = NULL`, `seedit = sd_i`, `is.RCT = TRUE`, `adverse_outcome = FALSE`. The only mismatch found is the one the sgdef report flagged; **no second mismatch.**

### 1.7 Verdict

**(A) The analytic list was wrong.** The three drivers passed 13 confounders including `str2`; nothing in `forestsearch()` or `get_FSdata()` removes it; the payloads do not record the list; the realized winners contain `str2` in 4.4 / 3.1 / 5.0 % of replicates; every stored family was enumerated on 12 without it. → §3.

---

## 2. §3 — the corrected family, beside the superseded one

`oc_wrapper_grid_corrected_2026-08-30.R` is the 08-29 grid + sweep + null scripts with **one input changed**: `fs_args$confounders.name` is the driver's 13-vector verbatim (`stopifnot(identical(fs_args[-1], G_OLD$fs_args[-1]))` guards that nothing else moved). Seed `20260825`, draws `2e5`, `block = Inf`, same `c1` grid `seq(20, 120, 5)`, same targets. The fixture rebuild gates passed in every cell (alt: all three `|diff| < 1e-9`, `fs_dgm_scale` regions identical to the payload *and* to the 08-29 `.rds`; null: all six gate fields `TRUE`); the measured columns recomputed from `payload$oc` / `$results` are `all.equal()` to those stored 08-29 (same payloads, `stopifnot`). The three cells ran as concurrent processes (`_alt500`, `_alt700`, `_null500` logs) and were merged into one `.rds` (`_merge` log). The sweep point at `c1 = 30` is `identical()` to the single-point `fs_oc_predict()` in all six (cell × gate) cases. Inversions use the order-statistic path (vector targets, one draw set per gate).

**Guards, re-run with the input unchanged and the code unchanged (GATE):** `fidelity_fs_oc_predict_2026-08-28.R` → `FIDELITY GATE: PASS (bit-identical)` (exit 0); `prerefactor_reference_2026-08-29.R check` → `REFACTOR GUARD: PASS (identical to the 0.2.4 reference)` (exit 0). Pass.

### 2.1 Families — stage counts

| cell | cut cols | enumerated | empty | minp | rmin | size | kept | duplicate | **M** | floor `n.min/n` | secs |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 500, 08-29 (12 vars) | 72 | 2628 | 119 | 0 | 99 | 631 | 1779 | 178 | **1601** | 0.120 | 8.5 |
| **500, corrected (13)** | **74** | **2775** | 127 | 0 | 107 | 656 | 1885 | 189 | **1696** | 0.120 | 9.5 |
| 700, 08-29 | 72 | 2628 | 119 | 0 | 94 | 435 | 1980 | 196 | **1784** | 0.0857 | 9.4 |
| **700, corrected** | **74** | **2775** | 127 | 0 | 102 | 448 | 2098 | 208 | **1890** | 0.0857 | 11.1 |
| null 500, 08-29 | 72 | 2628 | 119 | 0 | 99 | 631 | 1779 | 178 | **1601** | 0.120 | 8.6 |
| **null 500, corrected** | **74** | **2775** | 127 | 0 | 107 | 656 | 1885 | 189 | **1696** | 0.120 | 9.4 |

The 37th cut is `str2` (bare binary, both directions); it adds 147 combinations and **95 / 106 / 95 members** to the families. Every other cut expression is unchanged (population quantiles are independent of the list).

### 2.2 Single-point evaluations at the driver's `(c1, c2) = (30, 10)` — old · corrected · Δ · measured

**alt500** (M 1601 → 1696; predict 408 s / 785 s, peak 12.0 / 17.9 GB):

| quantity | resample old | **resample corr.** | Δ | split old | **split corr.** | Δ | measured (SE) |
|---|---:|---:|---:|---:|---:|---:|---:|
| det_rate | 0.99906 | 0.99906 | 0 | 0.99999 | 1.00000 | +0.00001 | 1.000 (0) |
| E\|Ĥ\| | 70.88 | **70.79** | −0.10 | 70.83 | **70.78** | −0.04 | 72.34 (0.41) |
| E[sens] | 0.1638 | 0.1638 | +0.00001 | 0.1637 | 0.1635 | −0.0001 | 0.1706 (0.0035) |
| E[spec] | 0.8698 | 0.8701 | +0.0003 | 0.8699 | 0.8700 | +0.0001 | 0.8686 (0.0019) |
| E[PPV] | 0.3994 | 0.4000 | +0.0006 | 0.3996 | 0.3994 | −0.0002 | 0.4055 (0.0080) |
| E[NPV] | 0.6643 | 0.6643 | +0.0001 | 0.6643 | 0.6643 | 0 | 0.6650 (0.0016) |
| E[β(Ĥ)] | 31.745 | 31.754 | +0.009 | 31.748 | 31.745 | −0.003 | 31.765 (0.108) |
| naive bias | 75.44 | **75.53** | +0.09 | 75.36 | **75.38** | +0.03 | 74.27 (0.58) |
| mass_below | 0.2653 | 0.2647 | −0.0007 | 0.2650 | 0.2659 | +0.0009 | — |
| M | 1601 | **1696** | +95 | 1601 | **1696** | +95 | — |

**alt700** (M 1784 → 1890; 503 s / 973 s, peak 13.4 / 17.0 GB):

| quantity | resample old | **resample corr.** | Δ | split old | **split corr.** | Δ | measured (SE) |
|---|---:|---:|---:|---:|---:|---:|---:|
| det_rate | 0.99981 | 0.99976 | −0.00005 | 1.00000 | 1.00000 | 0 | 0.999 (0.001) |
| E\|Ĥ\| | 73.03 | **73.04** | +0.01 | 73.07 | **73.08** | +0.02 | 73.63 (0.42) |
| E[sens] | 0.1217 | 0.1219 | +0.0003 | 0.1218 | 0.1222 | +0.0003 | 0.1229 (0.0028) |
| E[spec] | 0.9048 | 0.9049 | +0.0001 | 0.9048 | 0.9049 | +0.0001 | 0.9042 (0.0015) |
| E[PPV] | 0.4027 | 0.4034 | +0.0007 | 0.4029 | 0.4038 | +0.0009 | 0.4016 (0.0087) |
| E[NPV] | 0.6621 | 0.6622 | +0.0001 | 0.6622 | 0.6623 | +0.0001 | 0.6618 (0.0012) |
| E[β(Ĥ)] | 31.790 | 31.800 | +0.010 | 31.793 | 31.806 | +0.012 | 31.763 (0.117) |
| naive bias | 75.53 | **75.57** | +0.04 | 75.49 | **75.49** | +0.01 | 77.09 (0.56) |
| mass_below | 0.2613 | 0.2605 | −0.0008 | 0.2622 | 0.2599 | −0.0023 | — |
| M | 1784 | **1890** | +106 | 1784 | **1890** | +106 | — |

**null500** (M 1601 → 1696; 412 s / 789 s):

| quantity | resample old | **resample corr.** | Δ | split old | **split corr.** | Δ | measured (SE) |
|---|---:|---:|---:|---:|---:|---:|---:|
| false declaration | 0.99700 (0.00012) | **0.99694** (0.00012) | −0.00006 | 0.99999 | 0.99999 | 0 | 0.998 (0.0014) |
| P1 range | 0.1211–0.3784 | 0.1215–0.3781 | | 0.3460–0.3778 | 0.3457–0.3774 | | |
| L_eff | 12.22 | 12.19 | −0.03 | 25.72 | 25.76 | +0.04 | |
| E\|Ĥ\| | 71.00 | **70.90** | −0.10 | 70.90 | **70.83** | −0.07 | 72.05 (0.40) |
| E[spec] | 0.8580 | 0.8582 | +0.0002 | 0.8582 | 0.8583 | +0.0001 | 0.8559 (0.0008) |
| E[β(Ĥ)] | 26.2552 | 26.2552 | 0 | 26.2552 | 26.2552 | 0 | 26.2552 |
| naive bias | 76.22 | **76.32** | +0.10 | 76.05 | **76.07** | +0.02 | 74.85 (0.57) |
| mass_below | 1 | 1 | 0 | 1 | 1 | 0 | — |
| M | 1601 | **1696** | +95 | 1601 | **1696** | +95 | — |

### 2.3 The `c1` sweep and the inversions

Sweep (`fs_oc_grid`, 21 points, both gates; draw 381 / 766 s at n = 500, 475 / 945 s at n = 700, 386 / 770 s at the null; sweep 74–86 s): along the whole curve the corrected declaration rate differs from the 08-29 one by **|Δ| ≤ 0.0030 at n = 500, ≤ 0.0021 at n = 700, ≤ 0.0034 at the null** — one to three MC SEs (≈ 0.001 in the falling part) — with no systematic sign; at the ceiling (`c1 ≤ 40`) the differences are ≤ 0.00006. The full 42-row tables per cell are in `oc_wrapper_confs_compare_2026-08-30.log` / `.rds`.

Inversions (order-statistic path; 405–843 s per gate):

| cell | gate | target | `c1` 08-29 | **`c1` corrected** | Δ | | cell | gate | target | `c1` 08-29 | **corrected** | Δ |
|---|---|---:|---:|---:|---:|---|---|---|---:|---:|---:|---:|
| 500 | resample | 0.80 | 91.90 | 92.04 | +0.14 | | 700 | resample | 0.80 | 93.11 | 93.19 | +0.07 |
| 500 | resample | 0.90 | 84.72 | 84.95 | +0.23 | | 700 | resample | 0.90 | 86.52 | 86.64 | +0.12 |
| 500 | resample | 0.95 | 78.99 | 79.14 | +0.16 | | 700 | resample | 0.95 | 81.28 | 81.28 | −0.00 |
| 500 | split | 0.80 | 91.85 | 91.90 | +0.05 | | 700 | split | 0.80 | 93.14 | 93.10 | −0.04 |
| 500 | split | 0.90 | 84.57 | 84.68 | +0.11 | | 700 | split | 0.90 | 86.52 | 86.60 | +0.08 |
| 500 | split | 0.95 | 78.83 | 78.81 | −0.02 | | 700 | split | 0.95 | 81.31 | 81.29 | −0.01 |

Null (type-I targets): resample 0.05 / 0.10 / 0.20 / 0.50 / 0.80 / 0.90 / 0.95 → `c1` 133.24 / 125.70 / 117.11 / 101.64 / 87.29 / 80.17 / 74.43 (08-29: 133.11 / 125.76 / 117.11 / 101.57 / 87.15 / 79.96 / 74.21; Δ −0.06 … +0.22); split 133.04 / 125.69 / 117.05 / 101.60 / 87.12 / 79.85 / 74.04 (08-29: 133.24 / 125.71 / 117.01 / 101.54 / 87.07 / 79.81 / 74.04; Δ −0.20 … +0.06). Ceilings 0.99694 (resample) / 1.0000 (split) vs 0.99700 / 1.0000. All attainable. **Every inverted `c1` moves by less than 0.25** — below the 0.2 gate-to-gate spread the grid report already quoted.

---

## 3. §4 — does the correction close the residuals?

All from `oc_wrapper_confs_compare_2026-08-30.R` (signature parser copied verbatim from `sgdef_selection_2026-08-29.R`; the measured side and the per-rule population values read from `sgdef_selection_2026-08-29.rds`). "Corrected" is the `"resample"` gate unless stated; `"split"` is within 0.001 on every signature quantity.

### 3.1 `str2` mass — **accounted for**

| cell | corrected analytic on `str2` signatures (resample / split) | measured (n / N; binomial SE) | difference |
|---|---:|---:|---:|
| alt500 | 0.0382 / 0.0378 | 0.0440 (44 / 1000; 0.0065) | −0.006 (0.9 SE) |
| alt700 | 0.0288 / 0.0284 | 0.0310 (31 / 999; 0.0055) | −0.002 (0.4 SE) |
| null500 | 0.0357 / 0.0352 | 0.0501 (50 / 998; 0.0069) | −0.014 (2.1 SE) |

The corrected family puts 3–4 % of its selection mass on `str2` rules, against the measured 3–5 %: agreement within noise at both alternatives, 2 SE low at the null. Within the `str2` block the analytic prefers `str2=0 & wtkg>` and `homo=0 & str2=0` (0.006–0.007 each) where the search picked them once or not at all, and under-weights `cd80<= & str2=1` (0.0016–0.0017 vs 0.005–0.006) — the same kind of within-block tilt as elsewhere, on 1–6 replicates.

### 3.2 The between-rule size gap — **not closed; unchanged (slightly wider)**

The sgdef report's decisive comparison, on the corrected family (population size of the rules the search picked, `mean(nPg_pop)`, minus the analytic `sel_c`-weighted expectation `EnH`; and the re-weighting check: measured signatures' population sizes weighted by analytic `sel_c`, which reproduced the 08-29 analytic figures to 71.97 / 72.82 / 71.73 before being applied here):

| cell | pop. size of realized rules | analytic E\|Ĥ\| 08-29 → **corrected** | between-rule 08-29 → **corrected** | re-weighted 08-29 → **corrected** | within-rule (unchanged) |
|---|---:|---:|---:|---:|---:|
| alt500 | 72.90 | 70.88 → **70.79** | **+2.02 → +2.11** (split +2.07 → +2.12) | 71.97 → 71.88 | −0.56 |
| alt700 | 73.65 | 73.03 → **73.04** | **+0.62 → +0.61** (split +0.59 → +0.57) | 72.82 → 72.68 | −0.02 |
| null500 | 72.55 | 71.00 → **70.90** | **+1.55 → +1.65** (split +1.65 → +1.72) | 71.73 → 71.60 | −0.50 |

Adding `str2` **lowers** the analytic E|Ĥ| at n = 500 and at the null: the `str2` candidates the analytic argmax selects are slightly *smaller* than the family average (`sel_c`-weighted `n·Pg` 70.4 vs 70.8 at n = 500; 70.1 vs 70.9 at the null; at n = 700 they are larger, 77.9 vs 73.0, and E|Ĥ| rises by 0.01). The between-rule gap therefore moves by +0.1 / −0.01 / +0.1 — the wrong direction at two of three cells and within a tenth of a subject everywhere. **The confounder list does not account for the between-rule size gap.**

### 3.3 E|Ĥ| against measurement — **unchanged**

| cell | measured (SE) | analytic 08-29 (gap) | **corrected (gap)** |
|---|---:|---:|---:|
| alt500 | 72.34 (0.41) | 70.88 (+1.45, 3.5 SE) | **70.79 (+1.55, 3.7 SE)** |
| alt700 | 73.63 (0.42) | 73.03 (+0.60, 1.4 SE) | **73.04 (+0.58, 1.4 SE)** |
| null500 | 72.05 (0.40) | 71.00 (+1.05, 2.6 SE) | **70.90 (+1.15, 2.8 SE)** |

Not closed, not partly closed: the residual moves by ≤ 0.1 subject, a quarter of its MC SE.

### 3.4 Naive bias — **unchanged; the sign flip survives**

| cell | measured (SE) | analytic 08-29 (measured − analytic) | **corrected (measured − corrected)** | split corrected |
|---|---:|---:|---:|---:|
| alt500 | 74.27 (0.58) | 75.44 (**−1.17**) | **75.53 (−1.26)** | 75.38 (−1.11) |
| alt700 | 77.09 (0.56) | 75.53 (**+1.55**) | **75.57 (+1.52)** | 75.49 (+1.59) |
| null500 | 74.85 (0.57) | 76.22 (**−1.37**) | **76.32 (−1.47)** | 76.07 (−1.22) |

The corrected mix moves the analytic naive bias by +0.03 to +0.10 — a tenth of the discrepancy, and *away* from the measurement at n = 500 and the null. The n = 500 / n = 700 sign flip is untouched. The sgdef report's mechanism (aggregate = mixture over signatures whose components differ by ±5–9) is consistent with this: the `str2` block carries ~4 % of the mass with purity 0.45–0.49 (`PQg` of the selected `str2` candidates) — close to the family's 0.40 — so it shifts the mixture by little.

### 3.5 PPV and sensitivity — **unchanged, as expected**

| cell | quantity | measured (SE) | analytic 08-29 (gap) | **corrected (gap)** |
|---|---|---:|---:|---:|
| alt500 | E[PPV] | 0.4055 (0.0080) | 0.3994 (+0.0061) | **0.4000 (+0.0054)** |
| alt500 | E[sens] | 0.1706 (0.0035) | 0.1638 (+0.0068) | **0.1638 (+0.0068)** |
| alt500 | E[spec] | 0.8686 (0.0019) | 0.8698 (−0.0012) | **0.8701 (−0.0015)** |
| alt700 | E[PPV] | 0.4016 (0.0087) | 0.4027 (−0.0011) | **0.4034 (−0.0019)** |
| alt700 | E[sens] | 0.1229 (0.0028) | 0.1217 (+0.0012) | **0.1219 (+0.0009)** |
| null500 | E[spec] | 0.8559 (0.0008) | 0.8580 (−0.0021) | **0.8582 (−0.0023)** |

PPV at n = 500 moves +0.0006 (the within-rule component of +0.0047 ± 0.0016 stays), sensitivity by 10⁻⁵. Nothing here is attributable to the list.

### 3.6 The distribution comparison on the corrected family

**Top 15 corrected analytic signatures beside 08-29 analytic and measured** — alt500 (`"split"` within 0.001 on every row):

| signature | corrected | 08-29 | measured (n) | | signature | corrected | 08-29 | measured (n) |
|---|---:|---:|---:|---|---|---:|---:|---:|
| age> & preanti<= | 0.035 | 0.039 | 0.036 (36) | | age> & wtkg<= | 0.020 | 0.020 | 0.024 (24) |
| age> & cd40> | 0.032 | 0.032 | 0.023 (23) | | age<= & preanti> | 0.018 | 0.019 | 0.019 (19) |
| age<= & age> | 0.032 | 0.031 | 0.036 (36) | | cd80<= & preanti> | 0.017 | 0.018 | 0.015 (15) |
| age> & cd80<= | 0.031 | 0.032 | 0.025 (25) | | cd80> & preanti> | 0.015 | 0.015 | 0.009 (9) |
| age> & cd40<= | 0.028 | 0.029 | 0.033 (33) | | cd40> & wtkg> | 0.015 | 0.015 | 0.008 (8) |
| age> & cd80> | 0.023 | 0.023 | 0.028 (28) | | cd40> & preanti> | 0.015 | 0.018 | 0.018 (18) |
| age> & wtkg> | 0.023 | 0.023 | 0.019 (19) | | cd80> & wtkg> | 0.014 | 0.015 | 0.015 (15) |
| preanti<= & preanti> | 0.020 | 0.021 | 0.012 (12) | | | | | |

alt700: `age<= & age>` 0.074 (08-29 0.074; measured 0.065), `preanti<= & preanti>` 0.054 (0.057; 0.057), `age> & preanti<=` 0.028 (0.033; 0.035) — the top three unchanged in order. null500: `age<= & age>` 0.033 (0.032; 0.027), `age<= & preanti>` 0.025 (0.026; 0.027), `preanti<= & preanti>` 0.019 (0.019; **0.008**) still the visible outlier.

Where the new mass comes from: the eight largest *decreases* are all `preanti<=` / `preanti>` signatures (`preanti<= & wtkg>` −0.004, `homo=0 & preanti<=` −0.004, `age> & preanti<=` −0.004, `cd40> & preanti>` −0.003, …) — `str2` is a coarsening of `preanti` (`cor = 0.68`; every `str2 = 0` patient has `preanti ≤ 744.5`), so the `str2` candidates take their mass from their `preanti` near-twins, not from the `age` rules.

| cell | analytic mass on never-selected signatures 08-29 → **corrected** | measured mass absent from the family 08-29 → **corrected** (signatures) | still-absent signatures, corrected | TVD 08-29 → **corrected** (resample / split) | TVD excl. `str2`, renormalised | multinomial floor 08-29 → corrected (mean; 2.5–97.5 %) | excess over floor 08-29 → **corrected** |
|---|---:|---:|---|---:|---:|---:|---:|
| alt500 | 0.0139 → 0.0168 | **0.0460 (21) → 0.0030 (3)** | `drugs=1 & symptom=0`, `karnof<= & race=1`, `str2=1 & symptom=1` (1 each) | 0.189 / 0.191 → **0.189 / 0.189** | 0.172 → 0.176 | 0.136 (0.118–0.155) → 0.144 (0.126–0.164) | 0.053 → **0.045** |
| alt700 | 0.0147 → 0.0199 | **0.0430 (21) → 0.0140 (8)** | `age>` (3), `drugs=1 & homo=0` (3), `gender=0 & wtkg>` (2), `karnof<= & symptom=1` (2), `cd80> & str2=1`, `gender=0 & race=0`, `hemo=1 & race=0`, `str2=0 & symptom=1` (1 each) | 0.174 / 0.173 → **0.171 / 0.169** | 0.160 → 0.164 | 0.132 (0.114–0.150) → 0.139 (0.120–0.157) | 0.042 → **0.033** |
| null500 | 0.0190 → 0.0190 | **0.0521 (23) → 0.0040 (4)** | `drugs=1 & symptom=0`, `karnof<= & race=1`, `str2=0 & symptom=0`, `str2=1 & symptom=1` (1 each) | 0.198 / 0.198 → **0.197 / 0.196** | 0.178 → 0.185 | 0.139 (0.121–0.157) → 0.145 (0.127–0.164) | 0.059 → **0.052** |

**Why the remaining absentees are absent.** "Absent from the analytic family" is, as in the sgdef report, *zero analytic selection mass on the signature*; checked against the corrected family's labels it splits into two cases. **(a) Not representable — below the size floor.** The realized rules evaluated on `df_super` sit **below `n.min = 60`** in the population: `drugs=1 & symptom=0` 54.1, `karnof<= & race=1` 53.3, `str2=1 & symptom=1` 50.5 subjects at n = 500 (and the last is not in the null family either); `karnof<= & symptom=1` 58.7, `drugs=1 & homo=0` 55.0, `gender=0 & wtkg>` 51.8, `gender=0 & race=0` 54.0, `hemo=1 & race=0` 46.5, `str2=0 & symptom=1` 56.0 at n = 700 — floor-boundary cases the population enumeration drops and a sample occasionally admits. **(b) Representable, never selected in 2×10⁵ draws.** `age>` alone at n = 700 — the family holds nine one-factor `age > v` members (v = 25 … 42), population size 66–72 for the realized `age > 45/46`, and the analytic argmax gives *all nine* `sel_c = 0` while the search picked one-factor `age>` in 3 of 999 replicates; `str2=0 & symptom=0` (in both families as `symptom != 1 & str2 != 1`, 174 subjects, `sel_c = 0`; measured once at the null); `cd80> & str2=1` (in the family at three `cd80` cuts, `sel_c` 0 at n = 700 and 2.5×10⁻⁵ at the null; measured once at n = 700 at the sample cut `cd80 > 1216`). Case (b) is a genuine, small between-rule difference of the same kind as §3.6's tilt — the search occasionally lands on large one-factor or coarse rules the analytic never selects — not a representability failure. Total still-absent mass 0.3 / 1.4 / 0.4 %: **the representability hole is closed** to the floor-boundary residue the sgdef report already characterised, plus ~0.5 % of case (b).

**What the TVD says.** Total variation is essentially unchanged (0.189 → 0.189; 0.174 → 0.171; 0.198 → 0.197): the 3–5 % of `str2` mass that was "absent" is now present but distributed over `str2` signatures differently from the search's picks, and the mass removed from the `preanti` signatures lands on rules the search did not choose, so the TVD *excluding* `str2` rises slightly (0.172 → 0.176 etc.). The multinomial floor rises with the larger support (183 → 208 signatures at n = 500), so the **excess over the floor falls from 0.053 / 0.042 / 0.059 to 0.045 / 0.033 / 0.052** — about a fifth of it was the missing variable; four fifths remain. The sgdef report's tilt groups are unmoved: `age<= & age>` 0.031 → 0.032 (measured 0.036), `age> & (cd40|cd80|wtkg)` 0.159 → 0.157 (0.152), `preanti<= & preanti>` 0.021 → 0.020 (0.012), `(cd40|cd80) × (wtkg|preanti)` 0.188 → 0.175 (0.176) at n = 500, and likewise at n = 700 and the null. The over-selection of the age-band rule and the under-selection of `preanti<= & preanti>` (0.020 vs 0.012 at n = 500; 0.019 vs 0.008 at the null) are exactly as before.

### 3.7 Per-residual verdict

| residual (sgdef report) | 08-29 | corrected | **does the confounder correction account for it?** |
|---|---:|---:|---|
| measured mass on signatures the family cannot contain (`str2`) | 4.6 / 4.3 / 5.2 % | 0.3 / 1.4 / 0.4 % | **Yes** — the `str2` mass is now in the family at 3.8 / 2.9 / 3.6 % against 4.4 / 3.1 / 5.0 % measured (0.9 / 0.4 / 2.1 SE) |
| between-rule E\|Ĥ\| gap | +2.02 / +0.62 / +1.55 | +2.11 / +0.61 / +1.65 | **No** — moves ≤ 0.1, wrong direction at two cells |
| E\|Ĥ\| vs measured | +1.45 / +0.60 / +1.05 | +1.55 / +0.58 / +1.15 | **No** — unchanged |
| naive bias vs measured, with sign | −1.17 / +1.55 / −1.37 | −1.26 / +1.52 / −1.47 | **No** — moves ≤ 0.1; the sign flip survives |
| PPV (within-rule, n = 500) | +0.0061 (within +0.0047) | +0.0054 | **No** (and none was expected) |
| sensitivity | +0.0068 / +0.0012 | +0.0068 / +0.0009 | **No** |
| TVD excess over the multinomial floor | 0.053 / 0.042 / 0.059 | 0.045 / 0.033 / 0.052 | **Partly** — ≈ one fifth; the age-band / `preanti`-band tilt is untouched |

**Mixed, and mostly negative.** The list correction fixes the one thing it could fix — the family now contains what the search selected — and leaves every functional residual where it was: the analytic argmax still weights the same signatures the same way relative to the search, and that weighting, not the candidate set, is where the +1.5 / +0.6 / +1.1 in E|Ĥ| and the ±1.2–1.6 in naive bias live. The M = 1601 / 1784 numbers were, in every functional, within 0.1 of the M = 1696 / 1890 ones.

---

## 4. What is superseded, and what stands

Everything computed on the 12-variable family is **numerically superseded** by the corrected `.rds`; every *conclusion* drawn from those numbers survives, because no functional moved by more than 0.1 (subjects, bias units) or 0.001 (rates). Concretely:

**`REPORT_oc_wrapper_fixture_run_2026-08-28.md`** — superseded: the `"split"` n = 500 column of its comparison table (M = 1601: E|Ĥ| 70.8, sens 0.164, spec 0.870, PPV 0.400, NPV 0.664, E[β] 31.75, naive 75.4; the 634 s runtime; "333 of 1601 candidates above 0.5"; "`confounders.name` = the 12 analysis variables" in its category line is now known to be the wrong list). Stands: everything qualitative — family size (not the gate) is what moved the document's M = 16 figures; the wrapper column sits next to the measured column on every row.

**`REPORT_oc_wrapper_gate_and_n700_2026-08-29.md`** — superseded: the stage-count table (72 / 2628 / … / 1601, 1784 → §2.1 here), the timing table (328 / 633 / 403 / 783 s, 11.4–16.2 GB → 408 / 785 / 503 / 973 s, 12.0–17.9 GB), both comparison tables (n = 500 and n = 700 → §2.2 here), and verdict line 5 ("M = 1601 … M = 1784"). Stands: the two gates agree to three decimals at the alternative; `"resample"` runs in about half the wall-clock of `"split"`; the memory notes.

**`REPORT_oc_wrapper_grid_2026-08-29.md`** — superseded: the sweep timing table (M and seconds), the twelve-row inversion table (`c1` 91.90 / 84.72 / 78.99 / 91.85 / 84.57 / 78.83 / 93.11 / 86.52 / 81.28 / 93.14 / 86.52 / 81.31 → §2.3 here, all within 0.25), the ceiling 0.9991 / 0.9998 (→ 0.99906 / 0.99976). Stands: the refactor guard and fidelity results (re-passed here); "`c1` for 80 % is 92–93, for 95 % 79–81"; "the two gates give the same `c1` to within 0.2".

**`REPORT_oc_wrapper_null_2026-08-29.md`** — superseded: "M = **1601**, the same 1601 candidates as the alternative" and its stage counts; the null comparison table (0.9970 / 1.0000 / 0.998; E|Ĥ| 71.0; naive 76.22 / 76.05; `L_eff` 12.2 / 25.7 → §2.2 here); the fourteen-row inversion table (133.11 … 74.21 / 133.24 … 74.04 → §2.3); the reduction-guard timings quoted for M = 1601. Stands: the gates differ only at the ceiling and by 0.003 (0.99694 vs 1.0000 now); the type-I curve agreement along `c1`; the measured 0.998 sits between the gates; the M = 16 `worked-null` was a family artefact; `L_eff` 12–26.

**`REPORT_oc_wrapper_sgdef_2026-08-29.md`** — **stands** in §2 (tabulation), §3 (population evaluation), §4 (within-rule comparison, the decomposition table's *measured* and *population-of-realized* columns, and the between-rule verdict). Superseded: the decomposition table's *analytic* column and the between-rule figures +2.02 / +0.62 / +1.55 (→ +2.11 / +0.61 / +1.65; the conclusion is unchanged), the naive-bias analytic figures 75.44 / 75.53 / 76.22 (→ 75.53 / 75.57 / 76.32); §5's whole comparison table — the "analytic mass never selected", "measured mass absent … of which `str2`" columns, the TVD 0.17–0.20 / 0.13–0.14 floor, and explanation (i) ("because the wrapper's family was enumerated without `str2`") — replaced by §3.6 here; the §6 "open" item on whether the between-rule difference is a family effect is now **partly settled**: the missing variable is not it. The report's §1 fact (13 vs 12) is confirmed as stated.

**`oc_wrapper_verification.qmd`** — reads `oc_wrapper_grid_2026-08-29.rds` and `oc_wrapper_null_2026-08-29.rds`, so **every inline number in it** (M = 1601 / 1784, the stage counts, the comparison tables, the sweep curve, the inversion tables, the null section including `L_eff` and the M = 1601 sentences) is from the superseded family. **Not re-pointed in this session** (out of scope); it needs to read `oc_wrapper_grid_corrected_2026-08-30.rds` (`$alt`, `$null` mirror the two old objects' fields) and its prose "M = 1601" sentences need the new counts. The document's conclusions do not change.

**`HANDOFF_oc_wrapper_2026-08-29_v1.md`** — superseded numerically: "M = 1601 at n = 500, M = 1784 at n = 700" and every table and inverted-`c1` value quoted from the 08-29 `.rds` (the alternative tables, the null table with 0.9970, the 91.9 / 84.7 / 79.0 and 93.1 / 86.5 / 81.3 inversions, the null 133 / 101.5 / 74). Stands: §5's categorisation (a) exact / (b) signed gaps / (c) naive bias — with (b) already re-attributed by the sgdef report (E|Ĥ| is between-rule, PPV within-rule) and now known **not** to be a confounder-list effect either; the "no threshold delivers both 80 % power and small type-I" reading; the closing observation that family *size* was what moved the document's figures — with the addendum that the last 6 % of family size (95 members) moves nothing.

---

## 5. Close-out

`git status --porcelain` before staging: the new files only (`oc_wrapper_grid_corrected_2026-08-30.R`, its four `.log`, `.rds`; `oc_wrapper_confs_compare_2026-08-30.R`, `.log`, `.rds`; this report). Staged by explicit path; commit hash recorded in the follow-up commit. **No push. No install** (nothing in `R/` changed). Wall-clock: the three cells ran concurrently, 66 min for the longest (n = 700); guards 3 min.

**Concurrency check.** The merge and comparison steps were triggered twice, concurrently, by two waiters (an operator error, not a script property). Both steps are deterministic from stored inputs (`merge` from the three part files; `compare` from the three `.rds`), so the comparison was re-run alone afterwards: its log is **byte-identical** (`cmp`) to the one in the repository, its `.rds` is `identical()` to the committed one with the `built_at` timestamp stripped, and a fresh `merge` into a scratch directory is likewise `identical()` to the committed `oc_wrapper_grid_corrected_2026-08-30.rds` minus `built_at`. Nothing was replaced.

Commits: `7e4a5dfd` task doc; `58024dfe` script, logs, rds and this report; the hash line in the follow-up commit.

---

## 6. Verdict (ten lines)

1. **(A):** the drivers passed 13 confounders (`str2` last, `include_str2 <- TRUE`, L267–269 / L302 / L550 / L562); nothing in `forestsearch()` / `get_FSdata()` removes a variable (the VI screen at `−0.2` with `max_n_confounders = 1000` is a no-op); the payloads record no list; the realized winners use `str2` in 44 / 31 / 50 replicates; every stored family was built on 12. No second mismatch in any other enumeration argument.
2. Corrected families: **M = 1696 / 1890 / 1696** (was 1601 / 1784 / 1601), 74 cut columns, +95 / +106 / +95 `str2` members; same seed and draws; fidelity and refactor guards pass.
3. Every functional moves by ≤ 0.1 subject / bias unit and ≤ 0.001 in rate; every inverted `c1` by < 0.25; the sweep curves by ≤ 0.003.
4. `str2` mass: **accounted for** — 3.8 / 2.9 / 3.6 % analytic vs 4.4 / 3.1 / 5.0 % measured (≤ 2.1 SE).
5. Between-rule E|Ĥ| gap: **not** accounted for — +2.02 / +0.62 / +1.55 → +2.11 / +0.61 / +1.65.
6. E|Ĥ| and naive bias against measurement: **unchanged**; the n = 500 / n = 700 sign flip of the naive bias survives intact.
7. PPV and sensitivity: unchanged, as expected; the within-rule PPV component stays.
8. Distribution: the representability hole closes (absent measured mass 4.6 / 4.3 / 5.2 % → 0.3 / 1.4 / 0.4 %, the residue below the size floor), the TVD does not (0.189 / 0.174 / 0.198 → 0.189 / 0.171 / 0.197), the excess over the multinomial floor falls by about a fifth; the age-band / `preanti`-band tilt is exactly as before — the new mass is taken from the `preanti` near-twins of `str2`.
9. Hence the between-rule weighting difference is **not a candidate-set effect**; what remains is on the argmax side (sample-quantile cuts, the statistics-keyed near-duplicate removal, the `rmin` translation) or a different weighting altogether — not separated here, no recommendation made.
10. Superseded: every number computed on the 12-variable family (§4 lists them per document); every conclusion stands; `oc_wrapper_verification.qmd` still reads the 08-29 `.rds` and needs re-pointing later.
