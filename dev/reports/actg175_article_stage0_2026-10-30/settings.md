# Settings — the ACTG 175 calls against the GBSG Table 1 calls and `?forestsearch` Examples 2–3

Sources. GBSG: `dev/reports/vignettes_stage0_2026-10-27/table1_calls/{fs,dina,grf}_call.R`.
Examples: `R/forestsearch_main.R` lines 1476–1517 (Example 2, binary OR) and 1519–1560 (Example 3, continuous MD).
ACTG 175 calls: `calls_2026b/{binary,continuous}_fit.R`, `calls_2026a/identifier_{fs,dina,grf}.R` (verbatim).
Only arguments that differ from at least one comparator are listed; "—" = not passed (package default).
Data roles (`df.analysis`, `outcome.name`, `event.name`, `treat.name`, `id.name`, `confounders.name`) differ by
construction and are listed once in the first table. `seedit = 8316951` and `is.RCT = TRUE` are common to all
five ACTG 175 calls and the GBSG calls (Examples 2 and 3 omit `is.RCT`). The ladder adds `quiet = TRUE`,
`details = FALSE`, `plot.sg = FALSE`, `plot.grf = FALSE`, `parallel_args = list(plan = "sequential")` to every call.

## 1. 2026b binary (Table S16 OR column) — `analysis_actg175_binary_mr.qmd`

| argument | 2026b binary | GBSG FS (Table 1) | Example 2 |
|---|---|---|---|
| data / outcome | ACTG 175 arms 1 vs 3, `y_binary_neg` (= event) | GBSG, `time_months` / `status` | as 2026b |
| covariates | 6 continuous + 6 binary (numeric 0/1) | 7 (`age … er`) | as 2026b |
| `outcome_type`, `effect_measure`, `adverse_outcome` | `"binary"`, `"OR"`, `TRUE` | — (survival), `est.scale = "hr"` | as 2026b |
| `conf_force` | — | `c("er <= 0", "pgr <= 0")` | — |
| `conf.cont_jcuts` | `list(cd40 = 10, wtkg = 10)` | `list(er = 10, pgr = 10)` | as 2026b |
| `max_n_confounders`, `cut_type`, `cont.cutoff`, `collapse_cuts` | 1000, `"default"`, 4, `TRUE` | same | — |
| `effect_neighborhood` | 0.20 | 0.20 | **0.10** |
| thresholds | `effect.threshold = 1.00`, `consistency.threshold = 1.00` | `hr.threshold = 1.00`, `hr.consistency = 1.00` | `effect.threshold = 1.0`, `consistency.threshold = 1.0` |
| `max_subgroups_search` | `Inf` | — | — |
| `m1.threshold`, `minp`, `stop_threshold` | `Inf`, 0.025, `NULL` | same | — |
| `fs.splits` | 1000 | 1000 | **500** |
| `mr_inference_args` | field list (as GBSG), `declaration_c0 = c(0.70, 0.75, 0.80, 0.85)` | same list and c0 | `list(draws = 5000L)` only |
| `parallel_args` | — (default) | — | `multisession`, 4 workers |

## 2. 2026b continuous (Table S16 MD column) — `analysis_actg175_continuous_mr.qmd`

| argument | 2026b continuous | GBSG FS (Table 1) | Example 3 |
|---|---|---|---|
| data / outcome | ACTG 175 arms 1 vs 3, `y_decline = cd40 - cd420`; no `event.name` | GBSG | as 2026b |
| binary covariates | **factors** (`as.factor`) | — | factors |
| `outcome_type`, `effect_measure`, `adverse_outcome` | `"continuous"`, `"MD"`, `TRUE` | survival, `est.scale = "hr"` | as 2026b |
| `conf.cont_jcuts` | all six continuous at 10 | `er`, `pgr` at 10 | as 2026b |
| `conf_force`, `max_n_confounders`, `cont.cutoff`, `collapse_cuts` | — | set | — |
| `cut_type` | `"default"` | `"default"` | — |
| `effect.threshold`, `consistency.threshold` | **10**, **5** (CD4 cells) | `hr.threshold = 1`, `hr.consistency = 1` | 10, 5 |
| `use_twostage` | **`TRUE`** | — | — |
| `fs.splits` | 500 | 1000 | 500 |
| `m1.threshold`, `minp`, `stop_threshold` | — | `Inf`, 0.025, `NULL` | — |
| `mr_inference_args` | `ci_method = "field"`, draws 5000, `include_complement`, `field_complement`, `field_scale_complement = "selected"`, `return_reselection`, `keep_declaration_field = TRUE`, `keep_field_matrix = FALSE`, `declaration_c0 = c(-5, -10, -15)`; **no** `confirm_rule`, `field_uniform`, `ij_residual`, `field_decompose`, `field_recovery` | full list incl. those five | as 2026b without `keep_*` and `declaration_c0` |
| `parallel_args` | — (default) | — | `sequential` |

## 3. 2026a FS (Table S3 FS rows) — `analysis_actg175_binary_multimethod_fixed_family.qmd` @ 43b051b6

| argument | 2026a FS | GBSG FS (Table 1) | Example 2 |
|---|---|---|---|
| outcome | binary OR, adverse | survival | binary OR, adverse |
| `conf.cont_jcuts` | `cd40`, `wtkg` at 10 | `er`, `pgr` at 10 | as 2026a |
| GRF/DINA knobs passed but inert (`use_grf = FALSE`, `subgroup_method = "consistency"`) | `grf_depth = 2`, `dmin.grf = 0`, `frac.tau = 0.80`, `vi.grf.min = -0.2`, `return_selected_cuts_only = TRUE`, `tune_grf = FALSE`, `dina_args = list()`, `grf_selection = "frontier"`, `grf_select_statistic = "effect"`, `dina_select_statistic = "effect"` | — | — |
| `cut_type`, `cont.cutoff` | `"default"`, 4 | same | — |
| `max_n_confounders`, `collapse_cuts`, `conf_force` | — | 1000, `TRUE`, forced | — |
| `effect_neighborhood` | **0.10** | 0.20 | 0.10 |
| thresholds | `effect.threshold = 1.0`, `consistency.threshold = 1.0`, `pconsistency.threshold = 0.90` | HR 1.00 / 1.00 / 0.90 | same as 2026a |
| `consistency_method` | `"resample"` | `"resample"` | `"resample"` |
| `fs.splits`, `use_twostage` | 500, **`TRUE`** | 1000, — | 500, — |
| `m1.threshold`, `minp`, `stop_threshold` | — | `Inf`, 0.025, `NULL` | — |
| `mr_inference_args` | `list(draws = 5000L, family_native_neighborhood = NULL)` (package defaults otherwise) | full field list | `list(draws = 5000L)` |
| `show_candidate_summary` | `TRUE` | — | — |
| `parallel_args` | multisession, `n_cores` | — | multisession, 4 |
| bootstrap | `forestsearch_bootstrap_dofuture(nb_boots = 1000, show_three = FALSE, details = FALSE, digits = 2)`, default seed | — | — |

## 4. 2026a DINA (Table S3 DINA rows)

| argument | 2026a DINA | GBSG DINA (Table 1) |
|---|---|---|
| `subgroup_method`, `dina_select_statistic` | `"dina"`, `"effect"` | same |
| `dina_args` | `list(seed = 8316951L, n_folds = 2L)` | `list()` |
| `sg_focus`, `selection_rule` | `"effMaxSG"`, `"neighborhood"` | same |
| `effect_neighborhood` | **0.10** | 0.20 |
| `effect.threshold` | 1.0 (OR; → log-OR floor) | `hr.threshold = 1.00` |
| `n.min` | 60 | 60 |
| `d0.min`, `d1.min`, `maxk`, `fs.splits`, `m1.threshold`, `minp`, `stop_threshold` | — | 10, 10, 2, 1000, `Inf`, 0.025, `NULL` |
| `use_lasso/grf/dina`, `max_n_confounders`, `conf_force`, `cut_type`, `cont.cutoff`, `conf.cont_jcuts`, `collapse_cuts` | — | set |
| `mr_inference_args` | `list(draws = 5000L, family_native_neighborhood = NULL)` | field list |
| `parallel_args` | — | — |

## 5. 2026a GRF (Table S3 GRF rows)

| argument | 2026a GRF | GBSG GRF (Table 1) |
|---|---|---|
| `subgroup_method`, `grf_selection`, `grf_select_statistic` | `"grf"`, `"frontier"`, `"effect"` | same |
| `grf_depth`, `dmin.grf` | 2, 0.0 | 2, 0.0 |
| `frac.tau` | 0.8 | — |
| `sg_focus`, `selection_rule` | `"effMaxSG"`, `"neighborhood"` | same |
| `effect_neighborhood` | **0.10** | 0.20 |
| thresholds | none passed (package default `effect.threshold`) | `hr.threshold = 1.00` |
| `n.min` | 60 | 60 |
| `d0.min`, `d1.min`, `maxk`, `fs.splits`, `m1.threshold`, `minp`, `stop_threshold` | — | set |
| `use_lasso/grf/dina`, `max_n_confounders`, `conf_force`, `cut_type`, `cont.cutoff`, `conf.cont_jcuts`, `collapse_cuts` | — | set |
| `mr_inference_args` | `list(draws = 5000L, family_native_neighborhood = NULL)` | field list |
