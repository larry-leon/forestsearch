# REPORT — applications inventory (read-only) — 2026-08-31

All gates passed. Read-only throughout: no renders, no cache files created/moved/removed, nothing under `R/` touched. The only files created are the task document (`dev/tasks/cc_task_applications_inventory_2026-08-31.md`, committed `dfc8f3c9`) and this report. Every claim below traces to a quoted line (`path:line`) or a command output captured in this task.

**Headline findings (details in §9):**
1. `HANDOFF_continuous_2026-08-27_v5.md` does not exist anywhere in the tree or in git history (`git log --all --diff-filter=A -- '*HANDOFF_continuous*'` returns nothing). Gate 3 passed on the repro report alone.
2. The 08-27 record was rendered under **0.2.2** (not 0.2.0, which is the version of the *reference* renders).
3. The `vi.grf.min` reliance premise is only partly true: fits 1–2 of every multimethod document and the `compare_selection_rules()` run pin `-0.2`; the gbsg `fs_dina`/`fs_grf`, the actg175 `fs_anchor`, and every call in the four gbsg sgfocus-family documents (`_effMaxSG`, `_sgfocus`, `_frozen_family`, `_maxeff_mrconfirm`) and the actg175 template do **not** pass it.
4. **No gbsg payload exists on disk** (`quarto/applications/gbsg/_payloads/` is absent). The only current-form payload is `actg175/_payloads/template_actg175_continuous/…` (version 0.2.0, `est_scale = "md"`). No binary (`"or"`) payload exists outside `actg175/_archive/`. The deep-dive pair therefore had to be: the template (md) and an archive binary payload (or).
5. Exactly one LOO cache exists: `quarto/GuoHe/cv_out_analysis_gbsg_survival_effMaxSG_hrMaxSG_neighborhood.rds` (tracked, committed `43b051b6` 2026-08-16, alongside a 0.2.0 render). It carries no version field. On an as-is re-render `analysis_gbsg_survival_effMaxSG.qmd` **would load it silently** (`file.exists()` guard, no version check). The three sibling gbsg documents would recompute (no matching key on disk; `_maxeff_mrconfirm` has `loo_use_cache <- FALSE`) **and each would write a new cache file into `quarto/GuoHe/`** — a render side effect on a tracked directory.
6. `analysis_actg175_continuous_compare_all.qmd` does not write a `<stem>_payload.rds`; it writes `selected_subgroups_continuous.rds` + `comparison_continuous.rds`. The crossanalysis document's binary inputs (`selected_subgroups_binary.rds`, `comparison_binary.rds`) have **no active writer** — the only writer is archived.

---

## 0 Environment

| item | value |
|---|---|
| hostname | `pop-os` |
| repo | `~/Documents/GitHub/forestsearch` |
| branch (`git rev-parse --abbrev-ref HEAD`) | `feature/glm-extension` |
| `git log -1 --oneline` (before task-doc commit) | `dfc8f3c9 tasks: add applications inventory task (2026-08-31)` is the commit made on arrival; the pre-existing HEAD was its parent |
| `git status --porcelain \| wc -l` | `0` (clean tree after the task-doc commit) |
| `Rscript -e 'cat(as.character(packageVersion("forestsearch")), "\|", R.version.string)'` | `0.3.1 | R version 4.6.1 (2026-06-24)` |
| `DESCRIPTION` | `Version: 0.3.1` (line 3) |
| worktrees (`git worktree list`) | main; `.claude/worktrees/mr-terminology` (`feature/mr-terminology`); `.claude/worktrees/replication-check` (`wt/replication-check`) |

Gate results: **Gate 1** pass (branch correct, tree clean). **Gate 2** pass (`quarto/applications/` exists; 72 `.qmd`). **Gate 3** pass on `dev/glm-continuous-sims/REPORT_repro_check_vs_0.2.0_2026-08-27.md`; `HANDOFF_continuous_2026-08-27_v5.md` not found by `find dev -iname '*repro*' -o -iname '*2026-08-27*'`, by `find . -iname '*HANDOFF*'`, or in git history.

Note on mtimes: every tracked file under `quarto/applications/` and `quarto/GuoHe/` carries mtime `2026-08-17 21:06` (e.g. `stat` of `analysis_gbsg_survival_sgfocus.qmd`), while `R/forestsearch_main.R` is `2026-08-29 10:52`. The 08-17 21:06 stamp is a checkout time, not a render/write time; git commit dates are used for provenance below.

## 1 Document census

`find quarto/applications -name '*.qmd' | sort` → **68** files; `git ls-files 'quarto/applications/**' | grep '\.qmd$'` → **68** files, same set. **No untracked `.qmd`.**

Archive-like paths (flagged): `actg175/_archive/` (20 files), `gbsg/_archive/` (13), `gbsg/_broken/` (16) = 49. Active: **19** — 8 under `actg175/`, 5 under `gbsg/`, 6 at top level. Exact active list with last commit (`git log -1 --format='%h %ad %s' --date=short -- <path>`):

| # | path | last commit |
|---|---|---|
| A1 | `actg175/analysis_actg175_binary_multimethod_fixed_family.qmd` | `43b051b6 2026-08-16 update gbsg and actg175 analysis runs` |
| A2 | `actg175/analysis_actg175_binary_multimethod_frontend.qmd` | `cf4d6432 2026-08-15 chore(gbsg): export reproducibility payloads, drop max_subgroups_search from documents` |
| A3 | `actg175/analysis_actg175_binary_multimethod_psi_v2_2.qmd` | `cf4d6432 2026-08-15` (same) |
| A4 | `actg175/analysis_actg175_binary_multimethod_psi_v3a.qmd` | `cf4d6432 2026-08-15` (same) |
| A5 | `actg175/analysis_actg175_binary_sgfocus.qmd` | `889dc05a 2026-08-15 Guard payload n_* counts against absent data objects` |
| A6 | `actg175/analysis_actg175_continuous_compare_all.qmd` | `cf4d6432 2026-08-15` (same) |
| A7 | `actg175/analysis_actg175_crossanalysis_summary.qmd` | `db5be171 2026-08-13 fix(quarto): results_dir and dirout compose instead of overriding` |
| A8 | `actg175/template_actg175_continuous.qmd` | `61dd99df 2026-08-17 Create template_actg175_continuous.qmd` |
| G1 | `gbsg/analysis_gbsg_survival_effMaxSG.qmd` | `48c419eb 2026-08-15 Key the LOO cache on document, focus and selection_rule` |
| G2 | `gbsg/analysis_gbsg_survival_frozen_family.qmd` | `48c419eb 2026-08-15` (same) |
| G3 | `gbsg/analysis_gbsg_survival_maxeff_mrconfirm.qmd` | `48c419eb 2026-08-15` (same) |
| G4 | `gbsg/analysis_gbsg_survival_multimethod.qmd` | `cf4d6432 2026-08-15` (same as A2) |
| G5 | `gbsg/analysis_gbsg_survival_sgfocus.qmd` | `48c419eb 2026-08-15` (same) |
| T1 | `biomarker_comparison.qmd` | `b81e011d 2026-04-14 Reorganize quarto/ into topic directories` |
| T2 | `count_data_demo.qmd` | `b81e011d 2026-04-14` |
| T3 | `forestsearch_scenario_validation.qmd` | `b81e011d 2026-04-14` |
| T4 | `forestsearch_scenario_validation_100.qmd` | `b81e011d 2026-04-14` |
| T5 | `validation_glm_simulation_study.qmd` | `b81e011d 2026-04-14` |
| T6 | `validation_hte_tests_crump.qmd` | `ede6b683 2026-05-28 create quarto presentations directory and claude talk` |

(19 active + 49 archive/broken = 68; every non-active entry is under `_archive/` or `_broken/`.)

Archive/broken last commits: all `actg175/_archive/*` → `dd60b5fe 2026-08-13 chore(actg175): archive superseded application files`; all `gbsg/_archive/*` and `gbsg/_broken/*` → `9475181e 2026-08-13 chore(quarto): archive superseded and broken GBSG application .qmd files`.

Tracked `.html` without a matching `.qmd` (orphans): `actg175/analysis_actg175_binary_multimethod.html`, `actg175/analysis_actg175_binary_multimethod_psi_v2_2new.html`, and three under `gbsg/_archive/` (`2026-07-08_*`).

## 2 Per-document extraction

Line-numbered blocks are verbatim (`awk` extraction from `forestsearch(` to the balancing paren). Summary fields are read off the quotes.

### 2.1 A1 `actg175/analysis_actg175_binary_multimethod_fixed_family.qmd`

Calls `forestsearch()` **4 times**.

```
  557| fs <- forestsearch(
  558|   # --- Required data + roles ----------------------------------------------
  559|   df.analysis               = actg_df,
  560|   confounders.name          = confounders.name,
  561|   outcome.name              = outcome.name,
  562|   event.name                = outcome.name,    # binary: outcome IS the event
  563|   treat.name                = treat.name,
  564|   id.name                   = id.name,
  565|   is.RCT                    = TRUE,
  566|   seedit                    = analysis_seed,
  567|   # --- GLM parameters -----------------------------------------------------
  568|   outcome_type              = fs_outcome_type,
  569|   effect_measure            = fs_effect_measure,
  570|   adverse_outcome           = TRUE,
  571|   # --- Screening engines (candidate generation) ---------------------------
  572|   use_lasso                 = FALSE,
  573|   use_grf                   = FALSE,
  574|   grf_depth                 = fs_grf_depth,
  575|   dmin.grf                  = fs_dmin_grf,
  576|   frac.tau                  = 0.80,
  577|   vi.grf.min                = -0.2,
  578|   return_selected_cuts_only = TRUE,
  579|   tune_grf                  = FALSE,
  580|   use_dina                  = FALSE,
  581|   dina_args                 = fs_dina_args,
  582|   # --- Candidate-cut construction -----------------------------------------
  583|   cut_type                  = "default",
  584|   cont.cutoff               = 4L,
  585|   conf.cont_jcuts           = list(cd40 = 10, wtkg = 10),
  586|   # --- Selection engine ---------------------------------------------------
  587|   subgroup_method           = fs_subgroup_method,        # consistency | dina | grf
  588|   grf_selection             = fs_grf_selection,          # used when subgroup_method = "grf"
  589|   grf_select_statistic      = fs_grf_select_statistic,   # frontier only
  590|   dina_select_statistic     = fs_dina_select_statistic,  # used when subgroup_method = "dina"
  591|   sg_focus                  = fs_sg_focus,
  592|   selection_rule            = fs_selection_rule,
  593|   effect_neighborhood       = 0.10,
  594|   # --- Tier-2 de-biased gate ----------------------------------------------
  595|   mr_inference               = run_mr,
  596|   mr_inference_args          = list(draws = mr_draws,
  597|                                    family_native_neighborhood = fs_family_native_nbhd),
  598|   # --- Subgroup constraints -----------------------------------------------
  599|   n.min                     = 60L,
  600|   d0.min                    = 10L,
  601|   d1.min                    = 10L,
  602|   maxk                      = 2L,
  603|   # --- Thresholds (OR scale) ----------------------------------------------
  604|   effect.threshold          = fs_or_threshold,
  605|   consistency.threshold     = fs_or_consistency,
  606|   pconsistency.threshold    = fs_pconsistency,
  607|   consistency_method        = fs_consistency_method,   # "resample" (gate requires it)
  608|   # --- Consistency search -------------------------------------------------
  609|   fs.splits                 = 500L,
  610|   use_twostage              = TRUE,
  611|   # --- Output / diagnostics -----------------------------------------------
  612|   show_candidate_summary    = TRUE,
  613|   details                   = TRUE,
  614|   quiet                     = FALSE,
  615|   plot.sg                   = FALSE,   # KM not applicable; use plot_sg_glm_outcomes()
  616|   # --- Parallel -----------------------------------------------------------
  617|   parallel_args = list(
  618|     plan = "multisession",
  619|     workers = n_cores,
  620|     show_message = TRUE
  621|   )
  622| )
```

```
 1020| fs_dina_screen <- forestsearch(
 1021|   # --- Required data + roles ----------------------------------------------
 1022|   df.analysis               = actg_df,
 1023|   confounders.name          = confounders.name,
 1024|   outcome.name              = outcome.name,
 1025|   event.name                = outcome.name,
 1026|   treat.name                = treat.name,
 1027|   id.name                   = id.name,
 1028|   is.RCT                    = TRUE,
 1029|   seedit                    = analysis_seed,
 1030|   # --- GLM parameters -----------------------------------------------------
 1031|   outcome_type              = fs_outcome_type,
 1032|   effect_measure            = fs_effect_measure,
 1033|   adverse_outcome           = TRUE,
 1034|   # --- Screening engines: GRF + LASSO + DINA ------------------------------
 1035|   use_lasso                 = TRUE,
 1036|   use_grf                   = TRUE,
 1037|   grf_depth                 = 2L,
 1038|   dmin.grf                  = 0,
 1039|   frac.tau                  = 0.80,
 1040|   vi.grf.min                = -0.2,
 1041|   return_selected_cuts_only = TRUE,
 1042|   tune_grf                  = FALSE,
 1043|   use_dina                  = TRUE,
 1044|   dina_args = list(
 1045|     seed = dina_seed, n_folds = dina_nfolds,
 1046|     scope = "wide", max_per_covariate = 3L, max_subgroups = 10L,
 1047|     selected_only = dina_screen_selected_only
 1048|   ),
 1049|   # --- Candidate-cut construction -----------------------------------------
 1050|   cut_type                  = "default",
 1051|   cont.cutoff               = 4L,
 1052|   conf.cont_jcuts           = list(cd40 = 10, wtkg = 10),
 1053|   # --- Selection engine ---------------------------------------------------
 1054|   subgroup_method           = "consistency",
 1055|   sg_focus                  = dina_screen_sg_focus,
 1056|   selection_rule            = "neighborhood",
 1057|   effect_neighborhood       = 0.10,
 1058|   # --- Subgroup constraints -----------------------------------------------
 1059|   n.min                     = 60L,
 1060|   d0.min                    = 10L,
 1061|   d1.min                    = 10L,
 1062|   maxk                      = 2L,
 1063|   # --- Thresholds (OR scale) ----------------------------------------------
 1064|   effect.threshold          = dina_screen_or_threshold,
 1065|   consistency.threshold     = dina_screen_or_consistency,
 1066|   pconsistency.threshold    = fs_pconsistency,
 1067|   consistency_method        = fs_consistency_method,
 1068|   # --- Consistency search -------------------------------------------------
 1069|   fs.splits                 = 500L,
 1070|   use_twostage              = TRUE,
 1071|   # --- Output / diagnostics -----------------------------------------------
 1072|   show_candidate_summary    = TRUE,
 1073|   details                   = TRUE,
 1074|   quiet                     = FALSE,
 1075|   plot.sg                   = FALSE,
 1076|   # --- Parallel -----------------------------------------------------------
 1077|   parallel_args = list(
 1078|     plan = "multisession",
 1079|     workers = n_cores,
 1080|     show_message = TRUE
 1081|   )
 1082| )
```

```
 1221| fs_dina <- forestsearch(
 1222|   # --- Required data + roles ----------------------------------------------
 1223|   df.analysis         = actg_df,
 1224|   confounders.name    = confounders.name,
 1225|   outcome.name        = outcome.name,
 1226|   event.name          = outcome.name,
 1227|   treat.name          = treat.name,
 1228|   id.name             = id.name,
 1229|   is.RCT              = TRUE,
 1230|   seedit              = analysis_seed,
 1231|   # --- GLM parameters -----------------------------------------------------
 1232|   outcome_type        = fs_outcome_type,
 1233|   effect_measure      = fs_effect_measure,
 1234|   adverse_outcome     = TRUE,
 1235|   # --- DINA-selection mode ------------------------------------------------
 1236|   subgroup_method     = "dina",
 1237|   dina_args           = list(seed = dina_seed, n_folds = dina_nfolds),  # fit tuning only
 1238|   dina_select_statistic = fs_dina_select_statistic,
 1239|   # --- Selection criteria (passed through to dina_subgroup) ---------------
 1240|   sg_focus            = fs_dina_sg_focus,
 1241|   selection_rule      = "neighborhood",
 1242|   effect_neighborhood = 0.10,
 1243|   n.min               = 60L,
 1244|   effect.threshold    = fs_dina_or_threshold,  # -> log-OR harm floor
 1245|   # --- Tier-2 de-biased gate ----------------------------------------------
 1246|   mr_inference         = run_mr,
 1247|   mr_inference_args    = list(draws = mr_draws,
 1248|                              family_native_neighborhood = fs_dina_family_native_nbhd),
 1249|   # --- Output -------------------------------------------------------------
 1250|   details             = TRUE
 1251|   # NOTE: subgroup_method = "dina" bypasses GRF, LASSO and the consistency
 1252|   # search; use_grf / use_lasso / use_twostage / fs.splits / maxk and related
 1253|   # parameters do not apply in this mode.
 1254| )
```

```
 1443| fs_grf <- forestsearch(
 1444|   # --- Required data + roles ----------------------------------------------
 1445|   df.analysis         = actg_df,
 1446|   confounders.name    = confounders.name,
 1447|   outcome.name        = outcome.name,
 1448|   event.name          = outcome.name,
 1449|   treat.name          = treat.name,
 1450|   id.name             = id.name,
 1451|   is.RCT              = TRUE,
 1452|   seedit              = analysis_seed,   # also seeds the forest
 1453|   # --- GLM parameters -----------------------------------------------------
 1454|   outcome_type        = fs_outcome_type,
 1455|   effect_measure      = fs_effect_measure,
 1456|   adverse_outcome     = TRUE,
 1457|   # --- GRF-selection mode -------------------------------------------------
 1458|   subgroup_method     = "grf",
 1459|   grf_selection       = fs_grf_demo_selection,   # "tree" | "frontier"
 1460|   grf_select_statistic = fs_grf_select_statistic, # frontier only
 1461|   grf_depth           = fs_grf_demo_depth,
 1462|   dmin.grf            = fs_grf_demo_dmin,          # DR-score harm floor (RD scale)
 1463|   frac.tau            = 0.8,
 1464|   # --- Selection criteria -------------------------------------------------
 1465|   sg_focus            = fs_grf_sg_focus,          # -> GRF frontier/tree rule
 1466|   selection_rule      = "neighborhood",
 1467|   effect_neighborhood = 0.10,
 1468|   n.min               = 60L,
 1469|   # --- Tier-2 de-biased gate ----------------------------------------------
 1470|   mr_inference         = run_mr,
 1471|   mr_inference_args    = list(draws = mr_draws,
 1472|                              family_native_neighborhood = fs_grf_family_native_nbhd),
 1473|   # --- Output -------------------------------------------------------------
 1474|   details             = TRUE
 1475|   # NOTE: subgroup_method = "grf" bypasses LASSO, DINA and the consistency
 1476|   # search; use_lasso / use_dina / use_twostage / fs.splits / maxk and related
 1477|   # parameters do not apply in this mode.
 1478| )
```

Field summary (setup-chunk lines):
- `vi.grf.min`: **`-0.2` explicit** in `fs` (L577) and `fs_dina_screen` (L1040); **absent** in `fs_dina` (L1221–1254) and `fs_grf` (L1443–1478). Under 0.3.1 the latter two take `NULL`; note both run `subgroup_method = "dina"`/`"grf"`, which the in-call comments (L1251–1253, L1475–1477) say bypass the consistency search.
- `sg_focus`: L85 `fs_sg_focus <- "effMaxSG"   # GLM-natural alias of "hrMaxSG"`; L178 `dina_screen_sg_focus <- "effMaxSG"`; L188 `fs_dina_sg_focus <- "effMaxSG"`; L204 `fs_grf_sg_focus <- "effMaxSG"`.
- selection rule: L86 `fs_selection_rule <- "neighborhood"`; literal `"neighborhood"` in calls 2–4.
- `consistency_method`: L93 `fs_consistency_method <- "resample"`.
- `parallel_args`: `list(plan = "multisession", workers = n_cores, show_message = TRUE)`; L56 `n_cores <- max(1L, floor(0.80 * parallel::detectCores(logical = FALSE)))`.
- seeds: L344 `analysis_seed <- 8316951L`; L345 `set.seed(analysis_seed)`; `seedit = analysis_seed` in every call; L996 `dina_seed <- 8316951L`.
- `max_subgroups_search`: comments only — L124 `# 2026-08 no max_subgroups_search truncation.`; L133 `# that capped candidate evaluation at 50 (max_subgroups_search = 50L). That`; L139 `# 5000 MR draws, max_subgroups_search = 50L.`
- Bootstrap: present, **ungated** chunk `` ```{r bootstrap} `` (L770); L773–775 `fs_bc <- forestsearch_bootstrap_dofuture(fs.est = fs, nb_boots = NB,`; L61 `NB <- 1000L`. Two further bootstraps for `fs_dina` (L1284 `nb_boots = NB`) and `fs_grf` (L1507).
- MR: L99 `run_mr <- TRUE     # add mr_inference = TRUE to the forestsearch call`; L100 `mr_draws      <- 5000L    # multiplier-bootstrap draws for the gate`; passed as `mr_inference = run_mr` in calls 1, 3, 4.
- LOO: L69 `run_loo <- FALSE    # Leave-one-out / N-fold OOB (oob-cv)`; chunk L1734 `` ```{r oob-cv, eval = run_loo} ``; no cache path — computes inline (L1737–1756, `Kfolds = nrow(fs$df.est)`).
- Payload: L117 `results_dir <- NULL`; L118 `dirout      <- NULL`; L1937 `.dirout   <- if (is.null(dirout)) .qmd_stem else dirout`; L1938 `.base     <- if (is.null(results_dir)) file.path(.qmd_dir, "_payloads") else results_dir`; L1939 `.out_dir  <- file.path(.base, .dirout)`; L2121 `.payload_file <- file.path(.out_dir, paste0(.dirout, "_payload.rds"))`; L2117 `est_scale = "or"`.
- Data: L369 `actg_df <- subset(speff2trial::ACTG175, arms %in% c(1L, 3L))`.

### 2.2 A2 `actg175/analysis_actg175_binary_multimethod_frontend.qmd`

Calls `forestsearch()` 4 times at L550–615, L1013–1075, L1214–1247, L1436–1471. The four call bodies are **argument-for-argument identical to A1's** (same text, shifted by −7 lines) with one difference in the first call:

```
  565|   use_lasso                 = TRUE,
  566|   use_grf                   = TRUE,
  ...
  570|   vi.grf.min                = -0.2,
```

- `vi.grf.min`: `-0.2` at L570, L1033; absent in `fs_dina`, `fs_grf`.
- `sg_focus` L85 `"effMaxSG"`; rule L86 `fs_selection_rule        <- "both"`; consistency L93 `"resample"`; parallel same as A1; seeds L337–338 `analysis_seed <- 8316951L; set.seed(analysis_seed)`, L989 `dina_seed <- 8316951L`.
- `max_subgroups_search`: no hits.
- Bootstrap: ungated `{r bootstrap}`, L768 `nb_boots = NB`, L61 `NB <- 1000L`.
- MR: L99 `run_mr <- TRUE`; L100 `mr_draws      <- 2000L`.
- LOO: L69 `run_loo <- FALSE`; L1727 `{r oob-cv, eval = run_loo}`; no cache.
- Payload: L117–118 `results_dir <- NULL / dirout <- NULL`; L1930–1931 same construction as A1; L2053 `.payload_file <- file.path(.out_dir, paste0(.dirout, "_payload.rds"))`; L2049 `est_scale = "or"`.

### 2.3 A3 `actg175/analysis_actg175_binary_multimethod_psi_v2_2.qmd`

Calls at L535–600, L998–1060, L1199–1232, L1421–1456 — text identical to A1's four calls (`use_lasso = FALSE, use_grf = FALSE` in call 1).
- `vi.grf.min` `-0.2` at L555, L1018; absent in `fs_dina`, `fs_grf`.
- `sg_focus` L85 `"effMaxSG"`; rule L86 `"both"`; consistency L93 `"resample"`; seeds L322–323, L974; NB L61 `1000L`; MR L99 `TRUE`, L100 `mr_draws <- 5000L`; LOO L69 `FALSE`, L1712 `{r oob-cv, eval = run_loo}`; payload L117–118, L1915–1916, L2099; L2095 `est_scale = "or"`. `max_subgroups_search`: none.

### 2.4 A4 `actg175/analysis_actg175_binary_multimethod_psi_v3a.qmd`

Line-for-line identical extraction to A3 (calls L535–600, L998–1060, L1199–1232, L1421–1456; `vi.grf.min` L555, L1018; rule L86 `"both"`; `mr_draws` L100 `5000L`; NB L61; LOO L69 `FALSE`; payload L1915–1916, L2099; L2095 `est_scale = "or"`).

### 2.5 A5 `actg175/analysis_actg175_binary_sgfocus.qmd`

One call:

```
  649| fs <- forestsearch(
  650|   # --- Required data + roles ----------------------------------------------
  651|   df.analysis               = actg_df,
  652|   confounders.name          = confounders.name,
  653|   outcome.name              = outcome.name,
  654|   event.name                = outcome.name,    # binary: outcome IS the event
  655|   treat.name                = treat.name,
  656|   id.name                   = id.name,
  657|   is.RCT                    = TRUE,
  658|   seedit                    = analysis_seed,
  659|   # --- GLM parameters -----------------------------------------------------
  660|   outcome_type              = fs_outcome_type,
  661|   effect_measure            = fs_effect_measure,
  662|   adverse_outcome           = TRUE,
  663|   # --- Screening engines (candidate generation) ---------------------------
  664|   use_lasso                 = TRUE,
  665|   use_grf                   = TRUE,
  666|   grf_depth                 = 2L,
  667|   dmin.grf                  = 0,
  668|   frac.tau                  = 0.80,
  669|   vi.grf.min                = -0.2,
  670|   return_selected_cuts_only = TRUE,
  671|   tune_grf                  = FALSE,
  672|   use_dina                  = FALSE,
  673|   # --- Candidate-cut construction -----------------------------------------
  674|   cut_type                  = "default",
  675|   cont.cutoff               = 4L,
  676|   conf.cont_jcuts           = list(cd40 = 10, wtkg = 10),
  677|   # --- Selection engine: driven entirely by the focus profile -------------
  678|   subgroup_method           = "consistency",
  679|   sg_focus                  = focus$supplied,
  680|   selection_rule            = focus$selection_rule,
  681|   effect_neighborhood       = focus$effect_neighborhood,
  682|   # --- Tier-2 de-biased gate ----------------------------------------------
  683|   # MR re-selects under the SAME rule the search used --
  684|   # .fs_mr_reselection_from_focus() maps hrMaxSG -> effMaxSG,
  685|   # hrMinSG -> effMinSG, hr -> maxcons, and so on -- so no reselection
  686|   # argument is needed here.
  687|   mr_inference              = run_mr,
  688|   mr_inference_args         = list(draws = mr_draws),
  689|   # --- Subgroup constraints (retained under every focus) ------------------
  690|   n.min                     = 60L,
  691|   d0.min                    = 10L,
  692|   d1.min                    = 10L,
  693|   maxk                      = 2L,
  694|   # --- Thresholds + search settings, pre-resolved by the focus profile ----
  695|   effect.threshold          = focus$effect_threshold,
  696|   consistency.threshold     = focus$consistency_threshold,
  697|   pconsistency.threshold    = focus$pconsistency,
  698|   minp                      = focus$minp,
  699|   use_twostage              = focus$use_twostage,
  700|   stop_threshold            = focus$stop_threshold,
  701|   consistency_method        = fs_consistency_method,
  702|   fs.splits                 = 500L,
  703|   # --- Output / diagnostics -----------------------------------------------
  704|   show_candidate_summary    = TRUE,
  705|   details                   = TRUE,
  706|   quiet                     = FALSE,
  707|   plot.sg                   = FALSE,   # KM not applicable; see plot_sg_glm_outcomes()
  708|   # --- Parallel -----------------------------------------------------------
  709|   parallel_args = list(
  710|     plan = "multisession",
  711|     workers = n_cores,
  712|     show_message = TRUE
  713|   )
  714| )
```

- `vi.grf.min`: `-0.2` (L669).
- `sg_focus`: L129 `fs_sg_focus            <- "effMaxSG"` → `focus$supplied`; rule L137 `fs_selection_rule      <- "both"        # "neighborhood" | "pareto" | "both"`; consistency L150 `"resample"`; parallel L709–713 (L351 `n_cores <- max(1L, floor(0.80 * parallel::detectCores(logical = FALSE)))`); seeds L411–412 `analysis_seed <- 8316951L; set.seed(analysis_seed)`.
- `max_subgroups_search`: L586 `                   use_twostage = FALSE, max_subgroups_search = Inf, minp = 0)` (provenance-assertion block, not a call argument).
- Bootstrap: ungated `{r bootstrap}` L887; L890–892 `fs_bc <- forestsearch_bootstrap_dofuture(fs.est = fs, nb_boots = NB,`; L354 `NB       <- 1000L   # full-bootstrap resamples (production: 500-2000)`.
- MR: L355 `run_mr   <- TRUE`; L356 `mr_draws <- 2000L   # multiplier draws`.
- LOO: L359 `run_loo <- FALSE    # leave-one-out / N-fold OOB`; L1232 `{r oob-cv, eval = run_loo}`; no cache (`Kfolds = nrow(fs$df.est)` L1238).
- Payload: L82–83 `results_dir <- NULL / dirout <- NULL`; L1353–1357 (same construction; `.payload_file <- file.path(.out_dir, paste0(.dirout, "_payload.rds"))`); L1433 `est_scale = "or"`. Subgroup definition stored at L1381 `sg_harm  = if (is.null(.fs)) NULL else .fs$sg.harm`.
- Data L435 `actg_df <- subset(speff2trial::ACTG175, arms %in% c(1L, 3L))`.

### 2.6 A6 `actg175/analysis_actg175_continuous_compare_all.qmd`

Two `forestsearch` entry points: a direct call (`fs_anchor`) and the wrapper `compare_selection_rules()` which the document says runs `forestsearch()` per combination (L23, L286).

```
  404| fs_anchor <- forestsearch(
  405|   df.analysis      = actg_df,
  406|   confounders.name = confounders.name,
  407|   outcome.name     = adverse_outcome,
  408|   treat.name       = treat.name,
  409|   id.name          = id.name,
  410|   outcome_type     = fs_outcome_type,
  411|   effect_measure   = fs_effect_measure,
  412|   adverse_outcome  = TRUE,
  413|   seedit           = analysis_seed,
  414| 
  415|   sg_focus         = "maxeffCons",
  416|   selection_rule   = "neighborhood",
  417|   consistency_method = fs_consistency_method,
  418| 
  419|   effect.threshold       = fs_md_threshold,
  420|   consistency.threshold  = fs_md_consistency,
  421|   pconsistency.threshold = fs_pconsistency,
  422|   use_twostage     = TRUE,
  423|   conf.cont_jcuts  = fs_conf.cont_jcuts,
  424|   cut_type         = "default",
  425|   maxk             = fs_maxk,
  426|   n.min            = fs_n_min,
  427|   d0.min           = fs_d0_min,
  428|   d1.min           = fs_d1_min,
  429|   fs.splits        = fs_splits,
  430| 
  431|   # MR-compatible front end: required by .validate_mr_configuration()
  432|   use_lasso        = FALSE,
  433|   use_grf          = FALSE,
  434|   use_dina         = FALSE,
  435| 
  436|   is.RCT           = fs_is_RCT,
  437|   parallel_args    = list(plan = "multisession", workers = n_cores),
  438|   details          = FALSE,
  439|   quiet            = FALSE,
  440| 
  441|   mr_inference     = TRUE,
  442|   mr_inference_args = list(ci_method = "ij", draws = fs_mr_draws,
  443|                            include_complement = TRUE)
  444| )
```

Wrapper call (pass-through arguments), L297–350 (excerpt of the relevant lines):

```
  297| comparison <- compare_selection_rules(
  298|   df.analysis      = actg_df,
  299|   sg_focus         = sg_focus_vec,
  300|   selection_rule   = selection_rule_vec,
  ...
  321|   seedit           = analysis_seed,
  ...
  326|   consistency_method = fs_consistency_method,
  ...
  328|   use_grf          = fs_use_grf,
  ...
  330|   use_lasso        = fs_use_lasso,
  ...
  337|   parallel_args    = list(
  338|     plan         = "multisession",
  339|     workers      = n_cores,
  340|     show_message = TRUE
  341|   ),
  ...
  346|   vi.grf.min       = fs_vi_grf_min,
  347|   dmin.grf         = fs_dmin_grf,
```

- `vi.grf.min`: **absent** in `fs_anchor` (L404–444; but `use_grf = FALSE` there, L433); **explicit** in the wrapper run via L146 `fs_vi_grf_min      <- -0.2` (with L143–144 `fs_use_lasso <- TRUE`, `fs_use_grf <- TRUE`).
- `sg_focus` / rule: L202 `sg_focus_vec       <- c("effMaxSG", "effMaxSG", "effMinSG", "effMinSG", "eff", "maxSG","minSG", "maxeffCons")`; L203 `selection_rule_vec <- c("pareto",   "both", "pareto", "both", "neighborhood", "neighborhood", "neighborhood", "neighborhood")`; anchor L415–416 `"maxeffCons"` / `"neighborhood"`.
- `consistency_method`: L116 `fs_consistency_method <- "resample"`.
- `parallel_args`: L437 and L337–341; L80 `n_cores <- max(1L, floor(0.95 * parallel::detectCores(logical = FALSE)))` (0.95, not 0.80).
- seeds: L83–84 `analysis_seed <- 8316951L; set.seed(analysis_seed)`; `seedit = analysis_seed`; L153 `frontier_params$seed = 1L` (CI seed, L306 `ci_seed`).
- `max_subgroups_search`: no hits.
- Bootstrap: **none** (`grep -nE 'NB|nb_boots|boot'` returns nothing).
- MR: anchor only — L441 `mr_inference     = TRUE`, L133 `fs_mr_draws        <- 5000L`; wrapper call passes no `mr_inference`.
- LOO: **none**.
- Payload write — **not the `<stem>_payload.rds` form**: L725–728
  ```
  725| saveRDS(selected_tidy,
  726|         file = file.path(output_dir, "selected_subgroups_continuous.rds"))
  727| saveRDS(comparison,
  728|         file = file.path(output_dir, "comparison_continuous.rds"))
  ```
  with L173–174 `results_dir <- NULL / dirout <- NULL`, L179–181 `.dirout <- if (is.null(dirout)) qmd_stem else dirout; .base <- if (is.null(results_dir)) file.path(qmd_dir, "_payloads") else results_dir; output_dir <- file.path(.base, .dirout)`.
- Data: L225 `actg_df <- subset(speff2trial::ACTG175, arms %in% c(1L, 3L))`; L237 `actg_df <- actg_df[!is.na(actg_df$cd420), ]`.

### 2.7 A7 `actg175/analysis_actg175_crossanalysis_summary.qmd`

Does **not** call `forestsearch()`. Inputs:
```
   73| results_dir_base <- file.path(qmd_dir, "_payloads")
   74| dirout_stem      <- "analysis_actg175_continuous_compare_all"
   75| input_dir        <- file.path(results_dir_base, dirout_stem)
  ...
   86| load_if_exists <- function(path, label) {
   87|   if (file.exists(path)) {
   88|     obj <- readRDS(path)
  ...
   99| tidy_cont <- load_if_exists(
  100|   file.path(input_dir, "selected_subgroups_continuous.rds"),
  101|   "continuous")
  102| tidy_bin  <- load_if_exists(
  103|   file.path(input_dir, "selected_subgroups_binary.rds"),
  104|   "binary")
  ...
  251| cmp_cont <- if (file.exists(file.path(input_dir, "comparison_continuous.rds"))) {
  252|   readRDS(file.path(input_dir, "comparison_continuous.rds"))
  253| } else NULL
  254| cmp_bin  <- if (file.exists(file.path(input_dir, "comparison_binary.rds"))) {
  255|   readRDS(file.path(input_dir, "comparison_binary.rds"))
  256| } else NULL
```
(§7 has the rest.)

### 2.8 A8 `actg175/template_actg175_continuous.qmd`

```
  101| fs_result <- forestsearch(
  102|   df.analysis      = actg_df,
  103|   confounders.name = confounders,
  104|   outcome.name     = "y_decline",
  105|   treat.name       = "treat",
  106|   id.name          = "id",
  107|   outcome_type     = "continuous",
  108|   effect_measure   = "MD",
  109|   adverse_outcome  = TRUE,
  110|   seedit           = 8316951L,
  111|   n.min            = 60L,
  112|   d0.min           = 10L,
  113|   d1.min           = 10L,
  114|   hr.threshold     = c1,
  115|   hr.consistency   = c2,
  116|   pconsistency.threshold = 0.80,
  117|   fs.splits        = 500L,
  118|   maxk             = 2L,
  119|   sg_focus         = "minSG",
  120|   selection_rule   = "neighborhood",
  121|   use_lasso        = FALSE,
  122|   use_grf          = FALSE,
  123|   is.RCT           = TRUE,
  124|   cont.cutoff      = 4L,
  125|   details          = TRUE,
  126|   quiet            = FALSE,
  127|   parallel_args    = list(plan = "sequential", workers = 1,
  128|                           show_message = FALSE)
  129| )
```
- `vi.grf.min` absent (prose L89: "Because no screener runs, the GRF screening parameters (`vi.grf.min`, …"); `use_grf = FALSE`.
- `sg_focus = "minSG"`, `selection_rule = "neighborhood"`; `consistency_method` absent (package default); parallel sequential; seeds L11 `set.seed(8316951L)`, `seedit = 8316951L`.
- No bootstrap, no MR, no LOO; `max_subgroups_search` none.
- Payload: L28–29 `results_dir <- NULL / dirout <- NULL`; L215–219 construction identical to A1; L219 `.payload_file <- file.path(.out_dir, paste0(.dirout, "_payload.rds"))`; L270 `est_scale = "md"`; L230 `sg_harm  = if (is.null(.fs)) NULL else .fs$sg.harm`.

### 2.9 G1 `gbsg/analysis_gbsg_survival_effMaxSG.qmd`

```
  539| fs <- forestsearch(
  540|   df.analysis,
  541|   outcome.name          = outcome.name,
  542|   event.name            = event.name,
  543|   treat.name            = treat.name,
  544|   id.name               = id.name,
  545|   confounders.name      = confounders.name,
  546|   is.RCT                = TRUE,
  547|   seedit                = 8316951,
  548|   est.scale             = "hr",
  549|   # candidate generation: consistency search only (no LASSO/GRF/DINA screens)
  550|   use_lasso             = FALSE,
  551|   use_grf               = FALSE,
  552|   use_dina              = FALSE,
  553|   # continuous-cut construction + pre-specified receptor candidates
  554|   conf.cont_jcuts       = list(er = 10, pgr = 10),
  555|   collapse_cuts         = TRUE,
  556|   conf_force            = c("er <= 0", "pgr <= 0"),
  557|   # selection: driven entirely by the focus profile
  558|   subgroup_method       = "consistency",
  559|   sg_focus              = focus$supplied,
  560|   selection_rule        = focus$selection_rule,
  561|   effect_neighborhood   = focus$effect_neighborhood,
  562|   # Tier-2 de-biased gate (MR) computed inside the fit.  MR re-selects under
  563|   # the SAME rule the search used -- .fs_mr_reselection_from_focus() maps
  564|   # hrMaxSG -> effMaxSG, hrMinSG -> effMinSG, hr -> maxcons, and so on -- so
  565|   # no reselection argument is needed here.
  566|   mr_inference          = TRUE,
  567|   mr_inference_args     = list(draws = mr_draws),
  568|   # subgroup constraints (retained under every focus: they define which
  569|   # candidates are ESTIMABLE, not which one wins)
  570|   n.min                 = 60, d0.min = 10, d1.min = 10, maxk = 2,
  571|   # thresholds + search settings, pre-resolved by the focus profile
  572|   hr.threshold          = focus$hr_threshold,
  573|   hr.consistency        = focus$hr_consistency,
  574|   pconsistency.threshold = focus$pconsistency,
  575|   minp                  = focus$minp,
  576|   use_twostage          = focus$use_twostage,
  577|   stop_threshold        = focus$stop_threshold,
  578|   consistency_method    = fs_consistency_method,
  579|   fs.splits             = 1000,
  580|   max.minutes           = 3,
  581|   show_candidate_summary = TRUE, details = TRUE, quiet = FALSE,
  582|   plot.sg               = TRUE,
  583|   parallel_args = list(plan = "multisession", workers = n_cores,
  584|                        show_message = TRUE)
  585| )
```
- `vi.grf.min`: **absent** (`use_grf = FALSE`, L551).
- `sg_focus`: L110 `fs_sg_focus            <- "effMaxSG"`; rule L116 `fs_selection_rule      <- "neighborhood"        # "neighborhood" | "pareto" | "both"`; both flow through `.fs_focus_profile()` (L239–241). Consistency L126 `"resample"`. Parallel L583–584. Seed `seedit = 8316951` (L547); no `set.seed` in the field grep.
- `max_subgroups_search`: L484 `                   use_twostage = FALSE, max_subgroups_search = Inf, minp = 0)` (assertion block); L800 `  "(\`max_subgroups_search -> Inf\`).\n\n",` (prose).
- Bootstrap: ungated `{r bootstrap}` (L734 region; L738 `nb_boots = NB`); L346 `NB          <- 1000L    # full-bootstrap resamples (production: 500-2000)`.
- MR: in-fit, `mr_inference = TRUE` (L566); L347 `mr_draws    <- 5000L    # MR multiplier draws`. Also Guo–He: L348 `GH_B <- 2000L`, L354–355 `gh_requested <- TRUE; run_gh <- isTRUE(gh_requested) && focus$run_guohe` (skipped unless focus is `maxeff`, L277).
- LOO: L357 `run_loo      <- TRUE    # leave-one-out (N-fold) out-of-sample subgroup estimate`; L358 `run_cv_kfold <- TRUE`, L359 `Ksims <- 50L`; chunk L999 `{r loo, eval = run_loo}`; cache construction L66 `gh_dir <- Sys.getenv("GH_DIR", unset = "../../../quarto/GuoHe")` and:
  ```
 1017| .cache_doc   <- tryCatch(tools::file_path_sans_ext(basename(knitr::current_input())),
 1018|                          error = function(e) "analysis")
 1019| .cache_focus <- if (is.null(fs$sg_focus)) "unknownfocus" else as.character(fs$sg_focus)[1]
 1020| .cache_rule  <- {
 1021|   .r <- fs$args_call_all$selection_rule
 1022|   if (is.null(.r) || !nzchar(as.character(.r)[1])) "neighborhood" else as.character(.r)[1]
 1023| }
 1024| loo_cache <- Sys.getenv(
 1025|   "LOO_CACHE",
 1026|   unset = file.path(gh_dir, sprintf("cv_out_%s_%s_%s.rds",
 1027|                                     .cache_doc, .cache_focus, .cache_rule)))
 1028| if (file.exists(loo_cache)) {
 1029|   cv_out <- readRDS(loo_cache)
 1030|   timings$loo <- NA_real_
 1031|   cat("LOO loaded from cache:", loo_cache, "\n")
 1032| } else {
 1033|   t0 <- proc.time()
 1034|   fs_OOB <- forestsearch_Kfold(
 1035|     fs.est  = fs,
 1036|     details = FALSE,
 1037|     Kfolds  = nrow(df.analysis),   # N-fold = leave-one-out
 1038|     parallel_args = list(plan = "multisession", workers = loo_workers,
 1039|                          show_message = TRUE)
 1040|   )
 1041|   plan("sequential")
 1042|   timings$loo <- (proc.time() - t0)[["elapsed"]]
 1043|   cv_out <- forestsearch_KfoldOut(res = fs_OOB, details = FALSE, outall = TRUE)
 1044|   saveRDS(cv_out, loo_cache)
 1045| }
  ```
  (quoted from the line-identical G5; G1's L1020–1030 were verified identical.) The guard is `file.exists()` only — **no version check**.
- Payload: L83–84 `results_dir <- NULL / dirout <- NULL`; L1217–1221 construction identical to G5 quote (§2.13); L1301 `est_scale = "hr"`; L1303 `forestsearch_version = tryCatch(as.character(utils::packageVersion("forestsearch")), …`.
- Data: `df.analysis <- gbsg` (L391 in G5; G1 identical structure at L391–404 — `cat("Sample size:", nrow(df.analysis))` L404).

### 2.10 G2 `gbsg/analysis_gbsg_survival_frozen_family.qmd`

```
  256| fs <- forestsearch(
  257|   df.analysis,
  258|   outcome.name       = outcome.name,
  259|   event.name         = event.name,
  260|   treat.name         = treat.name,
  261|   id.name            = id.name,
  262|   confounders.name   = confounders.name,   # categorical only
  263|   is.RCT             = TRUE,
  264|   seedit             = 8316951,
  265|   est.scale          = "hr",
  266|   use_lasso          = FALSE,
  267|   use_grf            = FALSE,
  268|   use_dina           = FALSE,
  269|   conf_force         = frozen_cuts,        # all cuts fixed
  270|   # conf.cont_jcuts intentionally OMITTED -- no quantile cuts
  271|   collapse_cuts      = TRUE,
  272|   subgroup_method    = "consistency",
  273|   sg_focus           = "maxeff",
  274|   n.min              = 60, d0.min = 10, d1.min = 10, maxk = 2,
  275|   mr_inference        = TRUE,
  276|   mr_inference_args   = list(draws = mr_draws),
  277|   consistency_method = "resample",
  278|   # stop_threshold: INERT here.  forestsearch_main.R:1304 sets it to NULL under
  279|   # sg_focus = "maxeff", together with pconsistency.threshold -> 0,
  280|   # use_twostage -> FALSE, max_subgroups_search -> Inf and minp -> 0, and warns
  281|   # listing the overrides.  The engine never sees the value below.
  282|   fs.splits          = 1000, stop_threshold = 0.95,
  283|   use_twostage       = TRUE, max.minutes = 3,
  284|   show_candidate_summary = TRUE, details = TRUE, quiet = FALSE,
  285|   plot.sg            = TRUE,
  286|   parallel_args = list(plan = "multisession", workers = n_cores,
  287|                        show_message = TRUE)
  288| )
```
- `vi.grf.min` **absent**; `sg_focus = "maxeff"`; `selection_rule` **not passed** (cache-key fallback `"neighborhood"`, L560–563); `consistency_method = "resample"` literal; seed `seedit = 8316951`.
- `max_subgroups_search`: L213 `             "max_subgroups_search", "minp")` (name list in the provenance table), L280 comment.
- Bootstrap: ungated `{r bootstrap}` L351; L353–354 `fs_bc <- forestsearch_bootstrap_dofuture(fs.est = fs, nb_boots = NB, …`; L85 `NB         <- 1000L`.
- MR: `mr_inference = TRUE` in-fit; L86 `mr_draws <- 5000L`; L87 `GH_B <- 1000L`.
- LOO: L90 `run_loo    <- TRUE`; L537 `{r loo, eval = run_loo}`; L62 `gh_dir <- Sys.getenv("GH_DIR", unset = "../../../quarto/GuoHe")`; L564–567 same `sprintf("cv_out_%s_%s_%s.rds", .cache_doc, .cache_focus, .cache_rule)`; L568 `if (file.exists(loo_cache)) {`; L584 `saveRDS(cv_out, loo_cache)`.
- Payload: L79–80 `results_dir <- NULL / dirout <- NULL`; L743–747; L827 `est_scale = "hr"`.

### 2.11 G3 `gbsg/analysis_gbsg_survival_maxeff_mrconfirm.qmd`

```
  157| fs <- forestsearch(
  158|   df.analysis,
  159|   outcome.name          = outcome.name,
  160|   event.name            = event.name,
  161|   treat.name            = treat.name,
  162|   id.name               = id.name,
  163|   confounders.name      = confounders.name,
  164|   is.RCT                = TRUE,
  165|   seedit                = 8316951,
  166|   est.scale             = "hr",
  167|   # candidate generation: consistency search only (no LASSO/GRF/DINA screens)
  168|   use_lasso             = FALSE,
  169|   use_grf               = FALSE,
  170|   use_dina              = FALSE,
  171|   # continuous-cut construction + pre-specified receptor candidates
  172|   conf.cont_jcuts       = list(er = 10, pgr = 10),
  173|   collapse_cuts         = TRUE,
  174|   conf_force            = c("er <= 0", "pgr <= 0"),
  175|   # selection: the effect-maximiser argmax (Guo & He primitive)
  176|   subgroup_method       = "consistency",
  177|   sg_focus              = "maxeff",
  178|   # subgroup constraints
  179|   n.min                 = 60, d0.min = 10, d1.min = 10, maxk = 2,
  180|   # Tier-2 de-biased gate (MR) computed inside the fit; run_mr toggles it. The
  181|   # independent confirmation fit (see the MR-confirmation section) forces
  182|   # mr_inference = TRUE regardless, by replaying fs$args_call_all.
  183|   mr_inference           = run_mr,
  184|   mr_inference_args      = list(draws = mr_draws),
  185|   # consistency search
  186|   consistency_method    = "resample",
  187|   fs.splits             = 1000, stop_threshold = 0.95,
  188|   use_twostage          = TRUE, max.minutes = 3,
  189|   show_candidate_summary = TRUE, details = TRUE, quiet = FALSE,
  190|   plot.sg               = TRUE,
  191|   parallel_args = list(plan = "multisession", workers = n_cores,
  192|                        show_message = TRUE)
  193| )
```
Plus a second fit (L254 prose: "second `forestsearch()` fit with `mr_inference = TRUE` forced on. It replays the …" — chunk L273 `{r mr-confirm, eval = run_mr_confirm}`, L100 `run_mr_confirm     <- TRUE`).
- `vi.grf.min` absent; `sg_focus = "maxeff"`; no `selection_rule`; `consistency_method = "resample"`; seed 8316951.
- Bootstrap: ungated `{r bootstrap}` L204; L206–209 `fs_bc <- forestsearch_bootstrap_dofuture(fs.est = fs, nb_boots = NB, seed = 8316951, …`; L84 `NB <- 1000L`.
- MR: L93 `run_mr     <- TRUE`; L85 `mr_draws <- 5000L`; chunk L240 `{r mr, eval = run_mr}`.
- LOO: L103 `run_loo    <- TRUE`; L117 `loo_use_cache <- FALSE`; L483 `{r loo, eval = run_loo}`; L48 `gh_dir <- …"../../../quarto/GuoHe"`; L510–513 same key; **L514 `if (loo_use_cache && file.exists(loo_cache)) {`** → with `loo_use_cache <- FALSE` the cache is never read, the N-fold is always recomputed, and L530 `saveRDS(cv_out, loo_cache)` **writes** `quarto/GuoHe/cv_out_analysis_gbsg_survival_maxeff_mrconfirm_maxeff_neighborhood.rds` on every render.
- Payload: L65–66 NULL/NULL; L689–693; L777 `est_scale = "hr"`.

### 2.12 G4 `gbsg/analysis_gbsg_survival_multimethod.qmd`

Four calls. Call 1 (L693–789) excerpt of the discriminating lines (full text in the extraction file; every argument is named with its default in a trailing comment):
```
  693| fs <- forestsearch(
  ...
  706|   seedit                    = 8316951,        # default 8316951
  707|   est.scale                 = "hr",           # default "hr"
  709|   use_lasso                 = FALSE,           # default TRUE
  710|   use_grf                   = FALSE,           # default TRUE
  ...
  716|   vi.grf.min                = -0.2,           # default -0.2
  ...
  736|   subgroup_method           = fs_subgroup_method,  # set in setup chunk: consistency | dina | grf
  ...
  740|   sg_focus                  = fs_sg_focus,         # set in setup chunk
  741|   selection_rule            = "neighborhood",  # default "neighborhood"
  742|   effect_neighborhood       = 0.10,            # default 0.10
  744|   mr_inference               = run_mr,           # set in setup chunk
  745|   mr_inference_args          = list(draws = mr_draws,
  746|                                    family_native_neighborhood = fs_family_native_nbhd),
  ...
  758|   consistency_method        = fs_consistency_method,
  ...
  762|   fs.splits                 = 1000,           # default 1000
  763|   stop_threshold            = NULL,           # NULL is what this configuration RESOLVES to, not a
  ...
  784|   parallel_args = list(
  785|     plan = "multisession",
  786|     workers = n_cores,
  787|     show_message = TRUE
  788|   )
  789| )
```
Call 2 `fs_dina_screen` L1263–1355: `use_lasso = TRUE` (L1279), `use_grf = TRUE` (L1280), `dmin.grf = 4` (L1284), **`vi.grf.min = -0.2` (L1286)**, `use_dina = TRUE` (L1289), `sg_focus = dina_screen_sg_focus` (L1310), `selection_rule = "neighborhood"` (L1311), `consistency_method = fs_consistency_method` (L1324), parallel L1350–1354.
Call 3 `fs_dina` L1502–1534: `subgroup_method = "dina"` (L1515), `sg_focus = fs_dina_sg_focus` (L1520), `selection_rule = "neighborhood"` (L1521), `mr_inference = run_mr` (L1526) — **no `vi.grf.min`**.
Call 4 `fs_grf` L1743–1777: `subgroup_method = "grf"` (L1756), `sg_focus = fs_grf_sg_focus` (L1764), `selection_rule = "neighborhood"` (L1765), `mr_inference = run_mr` (L1769) — **no `vi.grf.min`**.
- Note the in-line comment at L716/L1286 says `# default -0.2` — stale under 0.3.x (NEWS.md L52–54: "The default of `forestsearch(vi.grf.min = )` changed from `-0.2` to `NULL`").
- `sg_focus`: L164 `fs_sg_focus <- "effMaxSG"`, L236, L246, L262 all `"effMaxSG"`; rule: literal `"neighborhood"` in every call (no `fs_selection_rule` variable); consistency L170 `"resample"`; seeds `seedit = 8316951`, L1239 `dina_seed <- 8316951`.
- `max_subgroups_search`: L412 `             "max_subgroups_search", "minp")`, L435 `  "maxk", "max_subgroups_search", "minp",` (provenance name lists).
- Bootstrap: ungated `{r bootstrap}`; L940–948 `fs_bc <- forestsearch_bootstrap_dofuture(mr_in_replicates = FALSE, fs.est = fs, nb_boots = NB, …`; L131 `NB <- 1000`; also L1582, L1817 for `fs_dina`/`fs_grf`.
- MR: L180 `run_mr <- TRUE`; L181 `mr_draws      <- 5000L`.
- LOO: L144–145 `run_cv  <- FALSE` / `run_loo <- FALSE`; L2135 `{r oob-cv, eval = run_loo}`; inline compute (L2139–2154, `Kfolds = nrow(df.analysis)`), **no cache path**.
- Payload: L198–199 NULL/NULL; L2442–2444; L2638 `.payload_file <- file.path(.out_dir, paste0(.dirout, "_payload.rds"))`; L2632 `est_scale = "hr"`; L2617 `labels = list(fs = .meta_fs$label, dina = .meta_dina$label, grf = .meta_grf$label, dina_standalone = dina_standalone)`.
- Data: L500 `df.analysis <- gbsg`; L518 `cat("Sample size:", nrow(df.analysis), "\n")`.

### 2.13 G5 `gbsg/analysis_gbsg_survival_sgfocus.qmd`

Call L539–585 is **textually identical** to G1's (§2.9). Everything else identical to G1 except L116 `fs_selection_rule      <- "both"        # "neighborhood" | "pareto" | "both"`. Payload construction:
```
 1213| .qmd_dir  <- tryCatch(dirname(knitr::current_input(dir = TRUE)),
 1214|                       error = function(e) getwd())
 1215| .qmd_stem <- tryCatch(tools::file_path_sans_ext(basename(knitr::current_input())),
 1216|                       error = function(e) "analysis")
 1217| .dirout   <- if (is.null(dirout)) .qmd_stem else dirout
 1218| .base     <- if (is.null(results_dir)) file.path(.qmd_dir, "_payloads") else results_dir
 1219| .out_dir  <- file.path(.base, .dirout)
 1220| dir.create(.out_dir, recursive = TRUE, showWarnings = FALSE)
 1221| .payload_file <- file.path(.out_dir, paste0(.dirout, "_payload.rds"))
 ...
 1243| payload <- list(
 1244|   table  = .table,
 1245|   labels = list(
 1246|     sg_harm  = if (is.null(.fs)) NULL else .fs$sg.harm,
 1247|     sg_focus = if (is.null(.fs)) NULL else .fs$sg_focus,
 1248|     focus    = .get("focus")),
 ...
 1288|     loo = if (isTRUE(.get("run_loo")) && !is.null(.loo)) list(
 1289|       row          = .loo,
 1290|       sens_metrics = if (is.null(.cvo)) NULL else .cvo$sens_metrics_original,
 1291|       find_metrics = if (is.null(.cvo)) NULL else .cvo$find_metrics,
 1292|       cache        = .get("loo_cache")) else NULL,
 ...
 1301|   est_scale = "hr",
 1302|   built_at  = Sys.time(),
 1303|   forestsearch_version = tryCatch(as.character(utils::packageVersion("forestsearch")),
 1304|                                   error = function(e) NA_character_))
```
Data L391 `df.analysis <- gbsg`.

### 2.14–2.19 Top-level documents

- **T1 `biomarker_comparison.qmd`** — no `forestsearch()` call. Input L140 `gbsg <- gbsg |>`; bootstrap of its own (L46 `sim_n_boot_mfpi <- 200L`, L339 `set.seed(2024L)`).
- **T2 `count_data_demo.qmd`** — one call on simulated data:
  ```
  653| fs <- forestsearch(
  654|   df.analysis      = df,
  ...
  660|   outcome_type     = "count",
  661|   effect_measure   = "IRR",
  ...
  671|   use_lasso        = TRUE,
  672|   use_grf          = TRUE,
  ...
  675|   seedit           = 8316951,
  ...
  678|   parallel_args    = sim_parallel
  679| )
  ```
  No `vi.grf.min` (→ `NULL` under 0.3.1 with `use_grf = TRUE`); no `sg_focus`/`selection_rule`/`consistency_method`; L35 `sim_parallel  <- list(plan = "multisession", workers = sim_n_workers,`; bootstrap L903 `nb = sim_n_boot` with L30 `sim_n_boot <- 50L`; no MR, no LOO cache, no payload.
- **T3/T4 `forestsearch_scenario_validation[_100].qmd`** — no direct `forestsearch(` call; parameters go through `base_fs_params` (T3 L313–320: `… maxk = maxk, vi.grf.min = vi_grf_min, seedit = sim_seed, …`), with L42 `vi_grf_min    <- -0.2` in both. `set.seed(42)` (T3 L367, T4 L468).
- **T5 `validation_glm_simulation_study.qmd`** — two calls on simulated data (L529–555, L621–646), `use_grf = TRUE`, **no `vi.grf.min`**, `plan = "sequential"`; no payload/LOO/bootstrap.
- **T6 `validation_hte_tests_crump.qmd`** — only illustrative snippets: L277 `  fs <- forestsearch(...)` and L1310 `fs <- forestsearch(df.analysis, ...)` (inside a non-executed fenced block, L1308 opens with bare ```` ``` ````). Data L947 `df_gbsg <- gbsg`.

### 2.20 Archive / broken (compact)

Per-file counts (`forestsearch(` assignments | `vi.grf.min` hits | `max_subgroups_search` hits | `cv_out` hits | `_payload` hits) were captured for all 49 files (command output). Notable: the four `actg175/_archive/20260730_*psi_v2_2*` files each have 4 calls, 2 `vi.grf.min`, **4 `max_subgroups_search`**, 3 `cv_out`, and a payload writer of the **old flat form** — `20260730_analysis_actg175_binary_multimethod_psi_v2_2A.qmd:2090: .payload_file <- file.path(.out_dir, "actg175_table2A_payload.rds")`. `gbsg/_archive/2026-07-30_analysis_gbsg_survival_maxeff.qmd` has 12 `cv_out` hits (the pre-48c419eb focus-only cache key). Every `gbsg/_archive` and `_broken` file with a `cv_out` hit (18 files) predates the doc/focus/rule key. `_archive/20260520_analysis_actg175_binary_compare_all.qmd:624–626` is the only writer of `selected_subgroups_binary.rds` / `comparison_binary.rds` (§7).

## 3 Summary table

| doc | `forestsearch()`? | `vi.grf.min` | `sg_focus` | rule | consistency | `parallel_args` | boot B | MR draws | LOO | `results_dir`/`dirout` override | payload path |
|---|---|---|---|---|---|---|---|---|---|---|---|
| A1 fixed_family | 4 | `-0.2` ×2; absent in `fs_dina`,`fs_grf` | effMaxSG | neighborhood | resample | multisession, 0.80·phys | 1000 ×3 | 5000 | `run_loo=FALSE`; no cache | both `NULL` (L117–118) | `_payloads/<stem>/<stem>_payload.rds`, `or` |
| A2 frontend | 4 | same pattern | effMaxSG | **both** | resample | same | 1000 ×3 | **2000** | FALSE; no cache | NULL/NULL | same, `or` |
| A3 psi_v2_2 | 4 | same | effMaxSG | both | resample | same | 1000 ×3 | 5000 | FALSE; no cache | NULL/NULL | same, `or` |
| A4 psi_v3a | 4 | same | effMaxSG | both | resample | same | 1000 ×3 | 5000 | FALSE; no cache | NULL/NULL | same, `or` |
| A5 binary_sgfocus | 1 | `-0.2` | effMaxSG | both | resample | same | 1000 | 2000 | FALSE; no cache | NULL/NULL (L82–83) | same, `or` |
| A6 continuous_compare_all | 1 direct + wrapper (8 combos) | anchor **absent** (`use_grf=FALSE`); wrapper `-0.2` | 8-vector (L202); anchor maxeffCons | 8-vector (L203); anchor neighborhood | resample | multisession, **0.95**·phys | **none** | anchor 5000; wrapper none | **none** | NULL/NULL (L173–174) | **`selected_subgroups_continuous.rds` + `comparison_continuous.rds`** in `_payloads/<stem>/` |
| A7 crossanalysis | no | — | — | — | — | — | — | — | — | hard-coded `results_dir_base`/`dirout_stem` (L73–74) | reads only |
| A8 template | 1 | absent (`use_grf=FALSE`) | minSG | neighborhood | default (not passed) | sequential | none | none | none | NULL/NULL (L28–29) | same form, `md` |
| G1 effMaxSG | 1 | absent (`use_grf=FALSE`) | effMaxSG | neighborhood | resample | multisession, 0.80·phys | 1000 (+GH_B 2000, gated off) | 5000 (in-fit `TRUE`) | **`run_loo=TRUE`; cache `gh_dir/cv_out_<doc>_<focus>_<rule>.rds`** | NULL/NULL (L83–84) | same form, `hr` |
| G2 frozen_family | 1 | absent | maxeff | not passed (→`neighborhood` in key) | resample | same | 1000 (+GH_B 1000) | 5000 (in-fit) | TRUE; cache | NULL/NULL (L79–80) | same, `hr` |
| G3 maxeff_mrconfirm | 1 (+confirm refit) | absent | maxeff | not passed | resample | same | 1000 (seed 8316951) | 5000 (`run_mr=TRUE`) | TRUE; **`loo_use_cache=FALSE`** → recompute + write | NULL/NULL (L65–66) | same, `hr` |
| G4 multimethod | 4 | `-0.2` ×2; absent in `fs_dina`,`fs_grf` | effMaxSG | neighborhood (literal) | resample | same | 1000 ×3 | 5000 | FALSE; no cache | NULL/NULL (L198–199) | same, `hr` |
| G5 sgfocus | 1 | absent | effMaxSG | **both** | resample | same | 1000 | 5000 (in-fit) | TRUE; cache | NULL/NULL | same, `hr` |
| T2 count_data_demo | 1 (sim) | absent, `use_grf=TRUE` | — | — | — | multisession | 50 | — | — | — | none |
| T3/T4 scenario_validation | via `base_fs_params` | `-0.2` (L42) | — | — | — | — | — | — | — | — | none |
| T5 validation_glm_sim | 2 (sim) | absent, `use_grf=TRUE` | — | — | — | sequential | — | — | — | — | none |
| T1, T6 | no | — | — | — | — | — | — | — | — | — | — |

## 4 Payloads on disk

`find quarto/applications -type d -name '_payloads*'` → **`quarto/applications/actg175/_payloads`** only. **No `gbsg/_payloads`** (`ls: cannot access 'quarto/applications/gbsg/_payloads': No such file or directory`). No version-tagged payload directories.

`find quarto/applications -name '*_payload.rds' -printf …`:

| path | mtime | size |
|---|---|---|
| `actg175/_payloads/template_actg175_continuous/template_actg175_continuous_payload.rds` | 2026-08-17 21:06 | 282154 bytes |
| `actg175/_archive/20260730_actg175_table2A_payload.rds` | 2026-08-17 21:06 | 1917 bytes |
| `actg175/_archive/20260730_actg175_table2Amac_payload.rds` | 2026-08-17 21:06 | 1935 bytes |

Two more `.rds` payloads don't match the `*_payload.rds` glob: `_archive/20260730_actg175_table2_payload_mac.rds` (1932 bytes), `_archive/20260730_actg175_table2_payload_mac_w2.rds` (1932 bytes). Also under `quarto/applications/_data/`: `selected_subgroups_binary.rds` (590 bytes, mtime 2026-05-14), `selected_subgroups_continuous.rds` (609 bytes).

`readRDS` summary (`version | built_at | est_scale | names | names(extras)`):

| payload | version | built_at | est_scale | names | extras |
|---|---|---|---|---|---|
| template_actg175_continuous | **0.2.0** | 2026-08-17 09:17:37 | **md** | table,labels,meta,extras,est_scale,built_at,forestsearch_version | fixed_family,fs_tables,grf |
| _archive table2A | 0.2.0 | 2026-07-08 07:47:26 | or | same 7 | concordance,harm_rates |
| _archive table2Amac | 0.2.0 | 2026-07-08 09:23:04 | or | same 7 | concordance,harm_rates |
| _archive table2_payload_mac | 0.2.0 | 2026-07-08 16:40:41 | or | same 7 | concordance,harm_rates |
| _archive table2_payload_mac_w2 | 0.2.0 | 2026-07-08 16:02:51 | or | same 7 | concordance,harm_rates |

**Deep dive 1 — `template_actg175_continuous_payload.rds`** (the only current-form payload; there is no gbsg payload to deep-dive, see §9):

```
List of 7
 $ table               :List of 20   (gt_tbl)
 $ labels              :List of 3
  ..$ sg_harm : Named chr [1:2] "!{cd40 <= 421}" "!{wtkg <= 82}"
  .. ..- attr(*, "names")= chr [1:2] "q18.0" "q12.0"
  ..$ sg_focus: chr "minSG"
  ..$ focus   : NULL
 $ meta                :List of 9
  ..$ n_total   : int 1083
  ..$ n_events  : NULL
  ..$ outcome   : chr "y_decline"
  ..$ event     : NULL
  ..$ treat     : chr "treat"
  ..$ covariates: chr [1:12] "age" "preanti" "wtkg" "karnof" ...
  ..$ c1        : num 10
  ..$ c2        : num 5
  ..$ args_call :List of 85
 $ extras              :List of 3
  ..$ fixed_family:List of 4
  ..$ fs_tables   :List of 2
  ..$ grf         :List of 4
 $ est_scale           : chr "md"
 $ built_at            : POSIXct[1:1], format: "2026-08-17 09:17:37"
 $ forestsearch_version: chr "0.2.0"
```
(`$ table` sub-fields elided; they are gt internals.) **Selected-subgroup definition lives in `labels$sg_harm`** — a named character vector of clause strings (names are the cut ids). `meta$args_call` (85 elements) carries the fitted call. This is the `sgfocus`-family shape (A5, A8, G1, G2, G3, G5 all write `labels$sg_harm = .fs$sg.harm`, e.g. G5 L1246).

**Deep dive 2 — `_archive/20260730_actg175_table2A_payload.rds`** (binary, `"or"`; the multimethod shape):

```
List of 7
 $ table               :'data.frame':	6 obs. of  13 variables:
  ..$ method   : chr [1:6] "FS" "FS" "DINA" "DINA" ...
  ..$ region   : chr [1:6] "H" "Hc" "H" "Hc" ...
  ..$ n        : int [1:6] 72 1011 98 985 95 988
  ..$ pct      : num [1:6] 7 93 9 91 9 91
  ..$ naive_est: num [1:6] 3.584 0.587 1.936 0.6 2.249 ...
  ..$ naive_lo : num [1:6] 1.315 0.454 0.837 0.463 0.972 ...
  ..$ naive_hi : num [1:6] 9.774 0.759 4.475 0.779 5.203 ...
  ..$ fb_est   : num [1:6] 2.431 0.606 0.862 0.636 1.091 ...
  ..$ fb_lo    : num [1:6] 0.804 0.411 0.232 0.396 0.396 ...
  ..$ fb_hi    : num [1:6] 7.355 0.893 3.203 1.022 3.01 ...
  ..$ mr_est   : num [1:6] 1.803 0.617 1.072 0.609 1.623 ...
  ..$ mr_lo    : num [1:6] 0.637 0.37 0.479 0.364 0.643 ...
  ..$ mr_hi    : num [1:6] 5.1 1.03 2.4 1.02 4.1 ...
 $ labels              :List of 4
  ..$ fs             : chr "!{wtkg <= 86} & !{cd40 <= 380}"
  ..$ dina           : chr "{preanti >= 849.40000000000009} & {cd40 >= 338}"
  ..$ grf            : chr "{wtkg > 84.369600000000005} & {cd40 > 368.20000000000005}"
  ..$ dina_standalone:List of 3
 $ meta                :List of 7
  ..$ B         : int 1000
  ..$ gate_draws: int 5000
  ..$ n_total   : int 1083
  ..$ cores     : int 128
  ..$ machine   :List of 18
  ..$ covariates: chr [1:12] "age" "preanti" "wtkg" "karnof" ...
  ..$ timings   :List of 3
 $ extras              :List of 2
  ..$ concordance:List of 4
  ..$ harm_rates :List of 2
 $ est_scale           : chr "or"
 $ built_at            : POSIXct[1:1], format: "2026-07-08 07:47:26"
 $ forestsearch_version: chr "0.2.0"
```
**Selected-subgroup definitions live in `labels$fs` / `labels$dina` / `labels$grf`** — single definition strings per method; `table$n` gives H/Hc sizes. No per-subject membership vector is stored in either shape, so a membership-based comparison must re-derive membership from the definition string against the data (or from `meta$args_call` + a refit). Note the current multimethod writers (A1 L2107–2108, G4 L2620–2621) name the draws field `mr_draws`, whereas this 07-08 payload has `gate_draws` — a shape drift between the archived payload and what a re-render would write.

The two deep-dive shapes are the only two top-level shapes present (all five payloads share the 7 top-level names).

## 5 LOO caches

`gh_dir` resolution: every gbsg sgfocus-family doc sets `gh_dir <- Sys.getenv("GH_DIR", unset = "../../../quarto/GuoHe")` (G1 L66, G2 L62, G3 L48, G5 L66). Relative to the qmd dir `quarto/applications/gbsg/`, `realpath -m` → **`/home/larryleon/Documents/GitHub/forestsearch/quarto/GuoHe`**. (The path is only correct if knitr's working directory is the qmd directory.)

Repo-wide `find . -name 'cv_out_*.rds' -not -path './.git/*'`:

| cache file | mtime | size | parsed `<doc>_<focus>_<rule>` | consumer on re-render |
|---|---|---|---|---|
| `quarto/GuoHe/cv_out_analysis_gbsg_survival_effMaxSG_hrMaxSG_neighborhood.rds` | 2026-08-17 21:06 (git: `43b051b6` 2026-08-16) | 1730 bytes | doc=`analysis_gbsg_survival_effMaxSG`, focus=`hrMaxSG`, rule=`neighborhood` | **G1** (`effMaxSG` → canonical `hrMaxSG`, rule `neighborhood`) — would be **loaded silently** (`file.exists()` only, L1028) |
| `.claude/worktrees/replication-check/quarto/GuoHe/cv_out_frozen_family.rds` | 2026-07-30 13:21 | 1189 bytes | old focus-less key (`frozen_family`) | none in the main tree (separate worktree, pre-`48c419eb` key form) |
| `.claude/worktrees/mr-terminology/quarto/GuoHe/cv_out_frozen_family.rds` | 2026-07-29 23:25 | 1189 bytes | same | none |

Contents of the one live cache (`str`, max.level 2): names `itt_tab, SG_tab_original, SG_tab_Kfold, CV_summary, sens_metrics_original, find_metrics, SGs_found, tab_all`; `CV_summary` vectors are length **686** and `SGs_found` is `chr [1:686, 1:2]` (N = 686 LOO folds); **no version, date, or commit field**. It is tracked in git and was committed with the 0.2.0 renders (`43b051b6`, whose `analysis_gbsg_survival_effMaxSG.html` reports `forestsearch 0.2.0`).

Keys that the other three consumers would look for (none on disk → recompute, then `saveRDS` into `quarto/GuoHe/`):
- G5 → `cv_out_analysis_gbsg_survival_sgfocus_hrMaxSG_both.rds`
- G2 → `cv_out_analysis_gbsg_survival_frozen_family_maxeff_neighborhood.rds` (rule falls back to `"neighborhood"` because no `selection_rule` is passed; the `.cache_focus` string is whatever `fs$sg_focus` holds for `"maxeff"` — not verified here)
- G3 → `cv_out_analysis_gbsg_survival_maxeff_mrconfirm_maxeff_neighborhood.rds` — written on every render regardless (L514 `loo_use_cache && …`, L117 `FALSE`).

A1–A5 and G4 have `run_loo <- FALSE` and no cache path; nothing is consumed or written.

## 6 The 08-27 record, quoted

Sources located: `dev/glm-continuous-sims/REPORT_repro_check_vs_0.2.0_2026-08-27.md` (113 lines). No task document for it exists under `dev/tasks/` (`ls dev/tasks` shows no `*repro*`; `git log --all -- 'dev/tasks/*repro*'` is empty). `HANDOFF_continuous_2026-08-27_v5.md`: **does not exist** (not in tree, not in git history). The 08-30 report `dev/glm-continuous-sims/REPORT_vi_grf_default_2026-08-30.md` is quoted where it names siblings.

### 6.1 Documents rendered
`REPORT_repro_check_vs_0.2.0_2026-08-27.md:13-16`:
```
| application | files | .qmd last commit | .html last commit | rendered |
|---|---|---|---|---|
| GBSG survival multi-method | `quarto/applications/gbsg/analysis_gbsg_survival_multimethod.{qmd,html}` | `cf4d6432` 2026-08-15 21:12 | `43b051b6` 2026-08-16 12:07 | Aug 15, 2026 |
| ACTG175 continuous compare-all | `quarto/applications/actg175/analysis_actg175_continuous_compare_all.{qmd,html}` | `cf4d6432` 2026-08-15 21:12 | `43b051b6` 2026-08-16 12:07 | Aug 16, 2026 |
```
Exactly **two** documents. "And their siblings" (from `dev/tasks/cc_task_vi_grf_default_2026-08-30.md:43`) is resolved by `REPORT_vi_grf_default_2026-08-30.md:81`: "Siblings: `analysis_gbsg_survival_effMaxSG.qmd` and `_frozen_family.qmd` (both with rendered `.html`) never mention `vi.grf.min`; the four `analysis_actg175_binary_multimethod_*.qmd` pass it in fits 1–2 …" — the siblings were **inspected, not rendered**, on 08-27.

### 6.2 Version + commit
`:1-5`: "reproduction check of the two application references against forestsearch 0.2.2 … the currently installed package (0.2.2, built from `ae6f4025`'s `R/`, HEAD `a4ec8c6d`) reproduces both pre-GLM-extension application references." `:108-110`: "`packageVersion("forestsearch")` = 0.2.2, Built 2026-08-27 18:03 UTC — from `ae6f4025`'s tree, the last commit touching `R/` (none since, through HEAD `a4ec8c6d`)." → rendered under **0.2.2 at HEAD `a4ec8c6d`** (R/ = `ae6f4025`). The "0.2.0" in the filename is the reference version (§6.5).

### 6.3 Comparison method
`:56-59`: "Tag-stripped whole-file diffs: 712 changed lines (gbsg), 122 (actg175). 640 lines were noise (per-render random `gt` CSS ids, render dates, wall-clock timing echoes and timing-table cells, worker/batch topology echoes, package attach messages, sessionInfo platform lines). Everything substantive:" followed by six attributed differences (`:61-87`). `:6-7`: "**No selection, estimate, CI, or bootstrap summary changed.**" → the method was a **diff of rendered `.html` output**. Payload fields, memberships, definition strings: **not stated in `REPORT_repro_check_vs_0.2.0_2026-08-27.md`** as compared objects (payloads are mentioned only as exports: `:37-38` "Payload exports resolved to scratch-stem `_payloads/` directories, so no tracked artifact was touched").

### 6.4 LOO cache handling
**Not stated in `REPORT_repro_check_vs_0.2.0_2026-08-27.md`.** (Consistent with §2: both rendered documents have `run_loo <- FALSE` — G4 L145, and A6 has no LOO — so no cache was consumed or written.)

### 6.5 Baseline provenance
`:27-30`: "The references were rendered on a **different machine, R, and BLAS**: R 4.5.2 on `aarch64-apple-darwin20` (macOS, Accelerate/vecLib), forestsearch **0.2.0**. The re-render used R 4.6.1 on x86_64 Linux (reference BLAS/LAPACK 3.12), forestsearch **0.2.2**, with 102/121 workers where the references used 11/13." `:18-23`: "Neither `.qmd` was modified after its `.html` was rendered, so the comparison is source-identical. … The files above are the pair committed together in `43b051b6` "update gbsg and actg175 analysis runs"". → Baseline = the **tracked `.html` renders** at `43b051b6`, produced under 0.2.0 on macOS/aarch64. Which 0.2.0 *payloads*: **not stated** (baseline was HTML, not payloads). Commit of the 0.2.0 package build: **not stated**.

### 6.6 Wall-clock
`:35-38`: "Re-render mechanics: scratch-named `.qmd` copies in their own directories, sequential renders with the installed package, exits 0/0, wall-clocks 12m53.5s (gbsg; reference ~15.3 min) and 6m00.2s (actg175)." → 12m53.5s + 6m00.2s = **18m53.7s ≈ 19 min** (the "≈19-minute total" is a sum computed here; the report itself never states a total). No per-section timings are stated.

## 7 Crossanalysis state

Reader (A7): L102–104 `tidy_bin  <- load_if_exists(file.path(input_dir, "selected_subgroups_binary.rds"), "binary")` guarded by L86–95 `load_if_exists` (`if (file.exists(path)) { obj <- readRDS(path) … } else { cat(sprintf("  [MISSING] …")); NULL }`); L254–256 `cmp_bin <- if (file.exists(file.path(input_dir, "comparison_binary.rds"))) { readRDS(…) } else NULL`; L259–263 prints "Skipping per-subject agreement: … comparison_binary.rds not available." when absent. `input_dir` = `<qmd_dir>/_payloads/analysis_actg175_continuous_compare_all` (L73–75). Its own prose table (L41) attributes `selected_subgroups_binary.rds` to `analysis_actg175_binary_compare_all.qmd`, which is **archived**.

Repo-wide `find`: `./quarto/applications/_data/selected_subgroups_binary.rds` (plus copies in the two `.claude/worktrees`). **No `comparison_binary.rds` anywhere.** The `_data/` copy is **not** in the directory A7 reads (`actg175/_payloads/analysis_actg175_continuous_compare_all/`, which does not exist on disk — only `_payloads/template_actg175_continuous/` does).

Writers (`git grep -n`): only `quarto/applications/actg175/_archive/20260520_analysis_actg175_binary_compare_all.qmd:624` `file = file.path(output_dir, "selected_subgroups_binary.rds"))` and `:626` `file = file.path(output_dir, "comparison_binary.rds"))`. **No active document writes either file.** → on an as-is render A7 would load `selected_subgroups_continuous.rds` only if A6 has been rendered first, report `[MISSING]` for the binary tidy file, and skip the per-subject agreement section.

## 8 Render-cost estimates — no timing runs

Basis codes: **(i)** per-document timing quoted in §6.6; **(ii)** apportionment of the 18m54s total by structure, labelled as such; **(iii)** N × single-search time; **U** = unknown, not measured in this task. All 08-27 timings are for 0.2.2 on this host with 102/121 workers and `run_loo = FALSE`.

| doc | structure (from §2) | estimate | basis |
|---|---|---|---|
| G4 multimethod | 3 fits + 3 FB (B=1000) + 3 MR (5000) + 1 DINA screen; LOO off | **12m53.5s** | (i) `:37` |
| A6 continuous_compare_all | 8 wrapper combos + 1 MR anchor (5000); no FB/LOO | **6m00.2s** | (i) `:37` |
| A1–A4 binary multimethod (each) | same skeleton as G4 (3 fits, 3×B=1000, MR 2000–5000), binary GLM, N = 1083 | **≈ 13 min each, ≈ 52 min for four** | (ii) G4's time carried over on identical structure; binary-vs-Cox engine cost not measured (U) |
| A5 binary_sgfocus | 1 fit + FB 1000 + MR 2000 + Pareto section; LOO off | **≈ 4–5 min** | (ii) ≈ ⅓ of G4 (one fit/FB/MR triple instead of three) |
| G1 effMaxSG | 1 fit + FB 1000 + MR 5000 + repeated 10-fold CV (Ksims 50) + LOO (**cached** → 0) | **≈ 4–5 min + 10-fold CV (U)** | (ii) as A5; K-fold CV cost U; LOO 0 s if the stale cache is left in place |
| G5 sgfocus | same as G1 but LOO **not cached** | **≈ 4–5 min + CV (U) + LOO** | (ii) + (iii) below |
| G2 frozen_family | 1 fit + FB 1000 + MR 5000 + GH (B 1000) + CV + LOO (not cached) | as G5 + GH (U) | (ii)/(iii) |
| G3 maxeff_mrconfirm | 2 fits (+confirm) + FB 1000 + MR 5000 + GH 2000 + CV + LOO (**always** recomputed) | as G2 + one extra fit | (ii)/(iii) |
| A8 template | 1 sequential fit, no FB/MR/LOO | small (U) | — |
| A7 crossanalysis | reads only | seconds | — |

LOO recompute cost (iii): **N = 686** for gbsg — confirmed three ways: G4 L500 / G5 L391 `df.analysis <- gbsg` with `Kfolds = nrow(df.analysis)` (G4 L2148, G5 L1037); `survival::gbsg` has 686 rows (`Rscript` output); the live cache's `CV_summary` vectors are length 686. Single-search time: **U — not stated in the 08-27 report and not measured in this task.** So LOO(gbsg) = 686 × t_search, with fold-level parallelism over `loo_workers = n_cores` (G5 L363; each fold's inner refit is pinned sequential, L361–362). For actg175, N = **1083** after the arms 1/3 subset (payload `meta$n_total = 1083`, both deep dives; A6 additionally drops `NA cd420` at L237, post-filter N not determined here) — irrelevant to LOO on an as-is render because every actg175 document has `run_loo <- FALSE`.

## 9 Deviations from the record's expectations

1. **`HANDOFF_continuous_2026-08-27_v5.md` does not exist** (tree or history). §6 is sourced from the repro report and the 08-30 vi.grf report only.
2. **08-27 render version was 0.2.2, not 0.2.0.** The reference renders were 0.2.0 (macOS/aarch64, R 4.5.2); the 08-27 check rendered 0.2.2 on this host and compared **HTML diffs**, not payloads or memberships (§6.3). The record's "0.2.0 payloads produced on which machine at which commit" has no answer in the source: the baseline was `.html` at `43b051b6`, and no 0.2.0 gbsg payload exists on disk.
3. **`vi.grf.min` reliance is partial, not general.** Explicit `-0.2`: A1–A5 fits 1–2, G4 fits 1–2, A6 wrapper run, T3/T4. Absent (→ `NULL` under 0.3.1): A1–A4 & G4 `fs_dina`/`fs_grf` (`subgroup_method = "dina"/"grf"`), A6 `fs_anchor` and A8 and G1/G2/G3/G5 (all `use_grf = FALSE`), T2 and T5 (`use_grf = TRUE`, simulated data). The in-call comments `# default -0.2` at G4 L716/L1286 are stale against NEWS.md L52–54.
4. **No gbsg payload on disk; only one current-form payload at all** (`template_actg175_continuous`, 0.2.0, `md`). The expected gbsg `"hr"` / actg175-binary `"or"` deep-dive pair could not be taken from current-form payloads; the `"or"` example is a 07-08 archive payload of the old flat-path form (`_archive/…_table2A_payload.rds`, writer `_archive/20260730_…psi_v2_2A.qmd:2090` → `"actg175_table2A_payload.rds"`, not `<stem>_payload.rds`).
5. **Payload path/shape deviations:** A6 writes `selected_subgroups_continuous.rds` + `comparison_continuous.rds`, not `<stem>_payload.rds`. Two shapes coexist among writers: multimethod (`labels$fs/dina/grf`, `meta$B`, `meta$mr_draws`) vs sgfocus-family (`labels$sg_harm`, `meta$NB`, `meta$loo`). The archived binary payload's `meta$gate_draws` differs from the current writers' `meta$mr_draws`. `results_dir`/`dirout` overrides default `NULL` in every writer (as expected); A7 uses hard-coded `results_dir_base`/`dirout_stem` instead.
6. **Cache key form confirmed** as `cv_out_<doc>_<focus>_<rule>.rds` with **no version** (G1/G2/G3/G5) — and the one live cache would be **silently loaded by G1** under 0.3.1 (guard is `file.exists()` only). G2, G3, G5 would each **write a new cache into the tracked `quarto/GuoHe/`** on render; G3 does so unconditionally (`loo_use_cache <- FALSE`). G4 and all actg175 docs use no cache at all. Two stale focus-less caches (`cv_out_frozen_family.rds`) sit only in `.claude/worktrees/`.
7. **`max_subgroups_search` occurrences** — all comments, provenance-assertion blocks, or name lists; none is a call argument in any active document (A1 L124/133/139; A5 L586; G1 L484/800; G5 L484/800; G2 L213/280; G4 L412/435). Archive files still pass it (4 hits each in the `20260730_*psi_v2_2*` files).
8. **`est_scale` values as expected** (`"hr"` ×5 gbsg, `"or"` ×5 actg175 binary, `"md"` template) — but the only on-disk current payload is `"md"`.
9. **Cross-document parameter drift** worth knowing before specifying a render: `selection_rule` is `"neighborhood"` in A1/G1/G4 but `"both"` in A2/A3/A4/A5/G5; `mr_draws` is 2000 in A2/A5 vs 5000 elsewhere; A1/A3/A4 run `use_lasso = use_grf = FALSE` in fit 1 whereas A2 and A5 run `TRUE`; A6 uses 0.95 of physical cores vs 0.80 elsewhere. `gh_dir` is a relative path that resolves correctly only when the working directory is the qmd's directory.
10. **Crossanalysis binary inputs have no active writer** (§7); `_data/selected_subgroups_binary.rds` (2026-05-14) is not on A7's read path.
11. Tracked orphan `.html` without `.qmd`: `actg175/analysis_actg175_binary_multimethod.html`, `actg175/analysis_actg175_binary_multimethod_psi_v2_2new.html`.
