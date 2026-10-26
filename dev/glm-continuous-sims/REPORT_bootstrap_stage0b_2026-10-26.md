# REPORT — bootstrap efficiency Stage 0b: the GLM paths and the consistency stage, measured, with the modification map (no package code modified)

**Audit basis.** CRAN `forestsearch` 0.4.0 (Date/Publication 2026-10-10; `Packaged: 2026-10-10 04:12:05 UTC; larryleon`). The installed build is the CRAN build (`Repository: CRAN`, `RemoteType: standard`, `Built: R 4.6.1; ; 2026-10-25 23:17:25 UTC; unix`). All line references in §4 are to the 0.4.0 sources; the nine `R/` files under study plus `NAMESPACE` are byte-identical to CRAN 0.4.0 in the working tree.

**Task:** `dev/tasks/cc_task_bootstrap_stage0b_2026-10-26.md` (committed alone as `9fc39fb1`; received in `~/Downloads` as `cc_task_bootstrap_stage0b_2026-10-26.md`, plus a byte-identical duplicate `cc_task_bootstrap_stage0b_2026-10-26 (1).md`). **Bundle:** `dev/glm-continuous-sims/stage0b_glm_2026-10-26/` (tarball route). **Run:** 2026-10-26, plain terminal (Claude Code shell, not RStudio) on `pop-os`, `bash run_stage0b.sh`; b0 → b3 all exit 0 (start ≈08:15:50, b3 done 08:19:18 PDT, ≈3.5 min wall). No gate fired.

## 1. Provenance

| item | value |
|---|---|
| HEAD at §1 | `d6613995c882d15a7168674bb7e9a7ba4ced3564` (b0 recorded `9fc39fb1…`, because it ran after the task-document commit, which touches only `dev/tasks/`) |
| branch | `master` |
| `git status --porcelain` at §1 (verbatim) | `?? quarto/extreme_subgroups/fixed_random/` (pre-existing dirt; never staged) |
| installed `forestsearch` | 0.4.0 · `Packaged: 2026-10-10 04:12:05 UTC; larryleon` · `Repository: CRAN` |
| survival · data.table | 3.8.12 · 1.18.6.1 |
| R | R version 4.6.1 (2026-06-24), Linux, 128 available cores |
| speff2trial | installed (b0 gate 2 passed) |
| `~/Downloads/stage0b_glm_2026-10-26.tgz` sha256 | `5f61d5182116b0e1668a9687025dcccb5bea5c6c2d0cea09009816478aab5a3c` (matches) |

**Bundle checksums.** After extraction, all eight matched the manifest. After the §2.4 transplant (`sha256sum` re-run):

```
98f5bc1af8bb91e556cc60477ed760a2649d062bd95aaa71c0d2f199b032903b  b0_provenance.R
08001177c58c43bb9c2223ccda64bfbcac2bf0f0acd29e791c9742c06212c8e6  b1_profiles.R
453573052f8e73d19902749da7999398d1539e8776d54f0317db662cb06ecfea  b2_lean.R
fcfa83c8775a4c74c9abc5010f7676fc36e9956015aa0c766f6404549696419d  b3_summary.R
024ae519c9d80081c7a21e0672b2d6521f594337e67b65246cd7e7996d4438b4  run_stage0b.sh
9429af0edb65b30c28d50991f43b0a497454dcd9e228ae42c14763fb91b9b101  stage0b_common.R
1374466802109dbdfcc613bfa5f2c3788204f634cab601ea78ffd28a60e5e219  fixture_binary.R       (transplanted; was 08daeb07…)
beda3110ddd2254b4da09235ae617cee8676d0a7ea92f861a1b0a71557f260df  fixture_continuous.R   (transplanted; was 16044dfb…)
```

The seven non-fixture files still match the manifest.

**Repo files under study vs CRAN 0.4.0 MD5** (`results/provenance_md5.csv`):

| file | CRAN 0.4.0 | repo | match |
|---|---|---|---|
| DESCRIPTION | 19b9459a03ed5c4a0c12943eb57d82b6 | 8c1cec4441bf8c2a8dff1ef720ec34b5 | DIFFERS |
| NAMESPACE | 00257d44d6f21d10f9af68b9f0dc227c | same | match |
| R/Cox_estimation_helpers.R | ffa9744bf4b3264041aa4c787cd5958f | same | match |
| R/consistency_resample.R | 861c6b4e7d00c8176bbab6b7f863f724 | same | match |
| R/forestsearch_main.R | 11d80e74c4fcdba43fe8bbf9a1287d3f | same | match |
| R/glm_effect_estimators.R | 8e8239f94e5f832af33e4d9cc0ff64d3 | same | match |
| R/subgroup_consistency_helpers.R | f32b55b62d683337f6e1b36b3aa3bc17 | same | match |
| R/subgroup_consistency_main.R | 1285150ed4891883468d01d5ee981b78 | same | match |
| R/subgroup_search.R | f06d6bba8135f41a57962a659ebbe601 | same | match |
| R/bootstrap_analysis_dofuture.R | a891d66e7337b9a1ebdf524463102978 | same | match |
| R/forestsearch_cross_validation.R | 359955f9ba9bd850c978230aee3b4429 | same | match |

`DESCRIPTION` differs only by the post-release development bump (b0's `git log` since 2026-10-09: `99d0f987 chore: begin 0.4.0.9000 development`, plus the citation commits `2405f11f`, `c373fd83`, and the release commit `48674761`). It is the same as at Stage 0. The timed code is the installed CRAN build regardless.

### 1.1 Fixture sources (§2.4 transplant)

Both GLM fixtures were located in this repo under `quarto/applications/actg175/`, the first route in the order of preference. Neither shipped default was kept. `~/GitHub/fs-glms-interpretable/quarto/actg175/` does not exist on this machine and was not needed.

- **Continuous (MD):** `quarto/applications/actg175/analysis_actg175_continuous_intervals.qmd`, chunks `data-prep` + `anchor`. This is the `y_decline` analysis with `adverse_outcome = TRUE`, `effect.threshold = 10`, `consistency.threshold = 5`, and twelve confounders. Of the two documents in the directory with that block, it is the one whose selection rule is effMaxSG at ε = 0.20 under `selection_rule = "neighborhood"`, which is the rule the task's default describes. The other document, `analysis_actg175_continuous_oc.qmd`, carries the same block with `sg_focus = "maxeffCons"`. The intervals document is also the latest one (commit of 2026-09-16), and its `anchor` chunk gates its own result.
- **Binary (OR):** `quarto/applications/actg175/analysis_actg175_binary_multimethod_fixed_family.qmd`. Chunks transplanted: `parallel-setup` (only the `fs_*` settings the call references), `analysis-roles`, `data-treat`, `data-outcomes`, `data-covariates`, and `forestsearch-main`. Five binary OR documents exist in the directory. This one is the most recently committed (`43b051b6`, 2026-08-16), its rendered HTML is committed, and its configuration note names it "the fixed-family configuration the method's Theorem 1 certifies". The other four are the `_frontend`, `_psi_v2_2` and `_psi_v3a` multimethod variants and `_sgfocus`.

Overrides applied, and nothing else: `mr_inference = FALSE` (with `mr_inference_args` dropped), `details = FALSE`, `quiet = TRUE`, `plot.sg = FALSE`, `plot.grf = FALSE`, and `parallel_args = list(plan = "sequential")`. In the continuous document, `details = FALSE`, `quiet = TRUE` and `parallel_args = sequential` are already its own values. `plot.grf` is absent from both documents and was added as specified. In the continuous fixture, `seedit = params$seed` resolves to the document's YAML `params: seed: 8316951` through a local `params <- list(seed = 8316951)`. Kept as in the documents: `use_grf = FALSE`, `use_lasso = FALSE`, `vi.grf.min` (binary: `-0.2`; continuous: not set, so the package default `NULL` applies), `use_twostage = TRUE`, `consistency_method = "resample"`, and every threshold.

**Transplanted continuous fixture (verbatim, `fixture_continuous.R` function body):**

```r
fixture_continuous <- function() {
  if (!requireNamespace("speff2trial", quietly = TRUE))
    stop("speff2trial is required for the ACTG175 fixtures")
  params <- list(seed = 8316951)

  # ---- data-prep (verbatim) ----
  cont_vars <- c("age", "preanti", "wtkg", "karnof", "cd40", "cd80")
  bin_vars  <- c("hemo", "homo", "drugs", "race", "gender", "symptom")
  confounders.name <- c(cont_vars, bin_vars)

  actg_df <- subset(speff2trial::ACTG175, arms %in% c(1L, 3L))
  actg_df$id <- seq_len(nrow(actg_df))
  actg_df$treat <- ifelse(actg_df$arms == 1L, 1L, 0L)
  actg_df <- actg_df[!is.na(actg_df$cd420), ]
  actg_df$y_decline <- actg_df$cd40 - actg_df$cd420
  for (v in bin_vars) actg_df[[v]] <- as.factor(actg_df[[v]])

  N <- nrow(actg_df)
  stopifnot(N == 1083L)

  # ---- anchor (verbatim parameter block, overrides marked) ----
  args <- list(
    df.analysis      = actg_df,
    confounders.name = confounders.name,
    outcome.name     = "y_decline",
    treat.name       = "treat",
    id.name          = "id",
    outcome_type     = "continuous",
    effect_measure   = "MD",
    adverse_outcome  = TRUE,
    seedit           = params$seed,
    sg_focus         = "effMaxSG",
    effect_neighborhood = 0.20,
    selection_rule   = "neighborhood",
    consistency_method = "resample",
    effect.threshold       = 10,
    consistency.threshold  = 5,
    pconsistency.threshold = 0.90,
    use_twostage     = TRUE,
    conf.cont_jcuts  = list(age = 10, preanti = 10, wtkg = 10,
                            karnof = 10, cd40 = 10, cd80 = 10),
    cut_type         = "default",
    maxk             = 2L,
    n.min            = 60L,
    d0.min           = 10L,
    d1.min           = 10L,
    fs.splits        = 500L,
    use_lasso        = FALSE,
    use_grf          = FALSE,
    use_dina         = FALSE,
    is.RCT           = TRUE,
    parallel_args    = list(plan = "sequential"),   # override (document: same)
    details          = FALSE,                       # override (document: same)
    quiet            = TRUE,                        # override (document: same)
    mr_inference     = FALSE,                       # override (document: TRUE; mr_inference_args dropped)
    plot.sg          = FALSE,                       # override (document: absent)
    plot.grf         = FALSE                        # override (document: absent)
  )
  list(label = "continuous: ACTG175 MD, committed application (effMaxSG, eps = 0.20)",
       df = actg_df, args = args,
       source = "quarto/applications/actg175/analysis_actg175_continuous_intervals.qmd, chunks data-prep + anchor")
}
```

**Transplanted binary fixture (verbatim, `fixture_binary.R` function body):**

```r
fixture_binary <- function() {
  if (!requireNamespace("speff2trial", quietly = TRUE))
    stop("speff2trial is required for the ACTG175 fixtures")

  # ---- parallel-setup: settings referenced by the call (verbatim) ----
  fs_sg_focus              <- "effMaxSG"   # GLM-natural alias of "hrMaxSG"
  fs_selection_rule        <- "neighborhood"       # neighborhood ∩ Pareto (effMaxSG/effMinSG)
  fs_or_threshold          <- 1.0
  fs_or_consistency        <- 1.0
  fs_pconsistency          <- 0.90
  fs_consistency_method    <- "resample"
  fs_subgroup_method <- "consistency"  # "consistency" | "dina" | "grf"
  fs_grf_selection   <- "frontier"     # only if subgroup_method == "grf": "tree" | "frontier"
  fs_grf_select_statistic <- "effect"  # only if subgroup_method=="grf" & grf_selection=="frontier": "dr" | "effect".
  fs_grf_depth       <- 2L             # GRF depth (screening + GRF selection)
  fs_dmin_grf        <- 0              # GRF DR-score screening floor (RD scale for binary; 0 = permissive)
  fs_dina_args       <- list()         # extra args forwarded to DINA selection
  fs_dina_select_statistic <- "effect" # only if subgroup_method=="dina": "dina" | "effect".

  # ---- analysis-roles (verbatim) ----
  analysis_seed <- 8316951L
  set.seed(analysis_seed)
  cont_vars <- c("age", "preanti", "wtkg", "karnof", "cd40", "cd80")
  bin_vars  <- c("hemo", "homo", "drugs", "race", "gender", "symptom")
  confounders.name <- c(cont_vars, bin_vars)
  id.name   <- "id"
  treat.name   <- "treat"
  outcome.name <- "y_binary_neg"   # FS adverse-scale outcome (1 = no improvement)
  fs_outcome_type   <- "binary"
  fs_effect_measure <- "OR"

  # ---- data-treat / data-outcomes / data-covariates (verbatim) ----
  actg_df <- subset(speff2trial::ACTG175, arms %in% c(1L, 3L))
  actg_df$id <- seq_len(nrow(actg_df))
  actg_df$treat <- ifelse(actg_df$arms == 1L, 1L, 0L)
  actg_df <- actg_df[!is.na(actg_df$cd420), ]
  actg_df$y_binary_pos <- as.integer(actg_df$cd420 > actg_df$cd40)
  actg_df$y_binary_neg <- 1L - actg_df$y_binary_pos
  for (v in bin_vars) actg_df[[v]] <- as.numeric(actg_df[[v]])

  # ---- forestsearch-main (verbatim parameter block, overrides marked) ----
  args <- list(
    df.analysis               = actg_df,
    confounders.name          = confounders.name,
    outcome.name              = outcome.name,
    event.name                = outcome.name,    # binary: outcome IS the event
    treat.name                = treat.name,
    id.name                   = id.name,
    is.RCT                    = TRUE,
    seedit                    = analysis_seed,
    outcome_type              = fs_outcome_type,
    effect_measure            = fs_effect_measure,
    adverse_outcome           = TRUE,
    use_lasso                 = FALSE,
    use_grf                   = FALSE,
    grf_depth                 = fs_grf_depth,
    dmin.grf                  = fs_dmin_grf,
    frac.tau                  = 0.80,
    vi.grf.min                = -0.2,
    return_selected_cuts_only = TRUE,
    tune_grf                  = FALSE,
    use_dina                  = FALSE,
    dina_args                 = fs_dina_args,
    cut_type                  = "default",
    cont.cutoff               = 4L,
    conf.cont_jcuts           = list(cd40 = 10, wtkg = 10),
    subgroup_method           = fs_subgroup_method,        # consistency | dina | grf
    grf_selection             = fs_grf_selection,          # used when subgroup_method = "grf"
    grf_select_statistic      = fs_grf_select_statistic,   # frontier only
    dina_select_statistic     = fs_dina_select_statistic,  # used when subgroup_method = "dina"
    sg_focus                  = fs_sg_focus,
    selection_rule            = fs_selection_rule,
    effect_neighborhood       = 0.10,
    mr_inference              = FALSE,                     # override (document: run_mr = TRUE; mr_inference_args dropped)
    n.min                     = 60L,
    d0.min                    = 10L,
    d1.min                    = 10L,
    maxk                      = 2L,
    effect.threshold          = fs_or_threshold,
    consistency.threshold     = fs_or_consistency,
    pconsistency.threshold    = fs_pconsistency,
    consistency_method        = fs_consistency_method,   # "resample" (gate requires it)
    fs.splits                 = 500L,
    use_twostage              = TRUE,
    show_candidate_summary    = TRUE,
    details                   = FALSE,                     # override (document: TRUE)
    quiet                     = TRUE,                      # override (document: FALSE)
    plot.sg                   = FALSE,   # KM not applicable; use plot_sg_glm_outcomes()
    plot.grf                  = FALSE,                     # override (document: absent)
    parallel_args             = list(plan = "sequential")  # override (document: multisession, n_cores)
  )
  list(label = "binary: ACTG175 OR, committed application (fixed family, effMaxSG, eps = 0.10)",
       df = actg_df, args = args,
       source = "quarto/applications/actg175/analysis_actg175_binary_multimethod_fixed_family.qmd, chunks parallel-setup + analysis-roles + data-* + forestsearch-main")
}
```

**Diff against the shipped defaults.** The full unified diffs are committed as `stage0b_glm_2026-10-26/logs/fixture_continuous.R.diff` and `logs/fixture_binary.R.diff`. Parameter by parameter:

| item | continuous: shipped default → transplanted | binary: shipped default → transplanted |
|---|---|---|
| treatment coding | `treat <- 1L - as.integer(arms == 1)` (arm 1 = 0) → `ifelse(arms == 1L, 1L, 0L)` (**arm 1 = ZDV+ddI = 1; reversed relative to the default**) | same reversal |
| rows | all arms 1/3 → `!is.na(cd420)` filter added (N = 1083 either way) | same |
| outcome | `y_decline = cd40 - cd420` (unchanged) | `y_bin = cd420 <= cd40` → `y_binary_neg = 1 - (cd420 > cd40)` (same values) |
| binary covariates | numeric → `as.factor()` | numeric → `as.numeric()` (unchanged in effect) |
| confounder order | age, wtkg, karnof, preanti, cd40, cd80, … → age, preanti, wtkg, karnof, cd40, cd80, … | same reorder |
| `event.name` | — | added (`= outcome.name`) |
| thresholds | `hr.threshold = 10`, `hr.consistency = 5` → `effect.threshold = 10`, `consistency.threshold = 5` (aliases, same values) | `hr.threshold = 1.0` → `effect.threshold = 1.0`; **`hr.consistency 0.8` → `consistency.threshold 1.0`** |
| `pconsistency.threshold` | **0.80 → 0.90** | **0.80 → 0.90** |
| `effect_neighborhood` | 0.20 (unchanged) | **0.20 → 0.10** |
| cuts | `conf.cont_jcuts` 10L on six vars + `collapse_cuts = TRUE` → 10 on the same six vars, `cut_type = "default"`, no `collapse_cuts` | 10L on six vars + `collapse_cuts = TRUE` → **`list(cd40 = 10, wtkg = 10)` only**, `cut_type = "default"`, `cont.cutoff = 4L` |
| `vi.grf.min` | not set (unchanged; package default NULL) | **not set → `-0.2`** (the GRF variable-importance forest runs; see §6) |
| GRF/DINA settings | — | `grf_depth`, `dmin.grf`, `frac.tau`, `return_selected_cuts_only`, `tune_grf`, `dina_args`, `grf_selection`, `grf_select_statistic`, `dina_select_statistic`, `subgroup_method = "consistency"` added (inert with `use_grf = use_dina = FALSE`, except through `vi.grf.min`) |
| `d0.min`/`d1.min` | — → 10L / 10L | 10 / 10 (unchanged) |
| `fs.splits` | — → 500L | — → 500L |
| `use_twostage` | — (package default) → `TRUE` | — → `TRUE` |
| `max.minutes` | 3 → removed (package default) | 3 → removed |
| `show_candidate_summary` | — | → `TRUE` (document's value; prints nothing under `quiet = TRUE` in the logs) |

**Identified subgroup vs the documents' committed records** (read-only):
- Continuous: `!{cd40 <= 507} & {gender}`. This **matches** the `anchor` chunk's own gate, `stopifnot(setequal(fs_anchor$sg.harm, c("!{cd40 <= 507}", "{gender}")), n_H == 79L, …)`, which is the result recorded in `REPORT_actg175_applied_effmaxsg_stage0_2026-09-16.md` §3. The committed payload `_payloads/analysis_actg175_continuous_oc_intervals/…_effmaxsg_payload.rds` was built by 0.3.5 on 2026-09-16.
- Binary: `!{wtkg <= 86} & !{cd40 <= 380}`. This **matches** the committed rendered output `analysis_actg175_binary_multimethod_fixed_family.html`, whose selected row 1 is OR 3.584, n = 72, `!{wtkg <= 86} & !{cd40 <= 380}`. No binary `.rds` payload is committed under `_payloads/`, so the HTML is the record.

**Against chat's dry run** (`~/Downloads/stage0b_dryrun_container_2026-10-26.md`, approximate defaults): the GLM candidate flows differ by more than an order of magnitude. Continuous passed effect: 146 here vs 4170 in the dry run. Continuous consistency candidates: 131 here vs 4032. Binary consistency candidates: 136 here vs 4085. This is consistent with the defaults' reversed treatment coding (the treatment-effect sign flips, so most cells pass a harm threshold), together with the tighter thresholds above. I did not test this beyond the diff. The dry run's GLM shares therefore describe a different, much heavier consistency stage, and none of its GLM figures carry over.

## 2. `results/stage0b_summary.md` (verbatim)

# Stage 0b -- GLM paths and the consistency stage (CRAN forestsearch 0.4.0)

Installed: forestsearch 0.4.0 (2026-10-10 04:12:05 UTC; larryleon; repository CRAN) · survival 3.8.12 · data.table 1.18.6.1 · R version 4.6.1 (2026-06-24) · cores 128
Repo: /home/larryleon/Documents/GitHub/forestsearch @ 9fc39fb17cc68ba99e1c2c7aacda0645b59fc4a0 (master); files matching CRAN 0.4.0: 10, differing: 1

## 1. The single sequential call and the replicate, per fixture

| fixture | source | n | call s (median of 3) | sequential replicate s (nb = 5) | search: evaluated / reached fit / passed effect | consistency: total / evaluated / passed | early stop | identified subgroup |
|---|---|---|---|---|---|---|---|---|
| survival | inst/doc/survival-post-selection.qmd (CRAN 0.4.0), mr_inference = FALSE | 686 | 3.154 | 3.139 | 2211 / 1410 / 121 | 120 / 120 / 16 | FALSE | {er <= 0} & {pgr <= 26} |
| continuous | quarto/applications/actg175/analysis_actg175_continuous_intervals.qmd, chunks data-prep + anchor | 1083 | 1.544 | 1.963 | 5886 / 4702 / 146 | 131 / 131 / 8 | FALSE | !{cd40 <= 507} & {gender} |
| binary | quarto/applications/actg175/analysis_actg175_binary_multimethod_fixed_family.qmd, chunks parallel-setup + analysis-roles + data-* + forestsearch-main | 1083 | 5.450 | 6.371 | 3240 / 2577 / 145 | 136 / 136 / 5 | FALSE | !{wtkg <= 86} & !{cd40 <= 380} |

Rprof buckets, median of 3 profiles, seconds of total time (a function reached through do.call() is not attributable by name; its children are):

| fixture | call s | search stage | of which per-candidate fit | consistency stage | of which pieces (fit + dfbeta) | consistency stage minus pieces (extraction, bookkeeping) | dedup / sort | GRF | LASSO | survfit medians |
|---|---|---|---|---|---|---|---|---|---|---|
| survival | 3.11 | 2.55 | 2.01 | 0.52 | 0.25 | 0.27 | 0.00 | 0.00 | 0.00 | 0.12 |
| continuous | 1.46 | 1.15 | 0.59 | 0.23 | 0.11 | 0.12 | 0.02 | 0.00 | 0.00 | 0.00 |
| binary | 6.33 | 5.14 | 4.63 | 0.28 | 0.15 | 0.13 | 0.01 | 1.63 | 0.00 | 0.00 |

## 2. Candidates: machinery removed per call, identity, and the 30% bar (per-replicate shape)

Machinery = (package function on the recorded inputs) - (lean equivalent on the same inputs), direct timers, median of 3 alternating reps. Share = machinery / measured sequential replicate. Identity = identical() on what the package consumes.

| candidate | fixture | what | n calls | package s | lean s | machinery s | share of replicate | identity | clears 30%? |
|---|---|---|---|---|---|---|---|---|---|
| **A** | survival | search-stage Cox fit (Stage 0, committed) | 1410 | 1.887 | 0.378 | **1.509** | **47.0%** | 1410/1410 | YES |
| **B** | survival | consistency pieces (fit + dfbeta), narrow frame | 120 | 0.261 | 0.079 | **0.182** | **5.8%** | rate identical 120/120; NULL sets agree TRUE | no |
| **H (+B)** | survival | whole consistency evaluator: lean extraction + lean pieces | 120 | 0.498 | 0.085 | **0.413** | **13.2%** | Pcons identical: TRUE (passed 16 vs 16; lean fallbacks 0) | no |
| **H alone** | survival | lean extraction + lean warm start, package pieces | 120 | 0.498 | 0.289 | **0.209** | **6.7%** | Pcons identical: TRUE | no |
| **E** | continuous | consistency pieces (fit + dfbeta), narrow frame | 131 | 0.096 | 0.013 | **0.083** | **4.2%** | rate identical 131/131; NULL sets agree TRUE | no |
| **H+E** | continuous | whole consistency evaluator: lean extraction + lean pieces | 131 | 0.215 | 0.015 | **0.200** | **10.2%** | Pcons identical: TRUE (passed 8 vs 8; lean fallbacks 0) | no |
| **H alone** | continuous | lean extraction, package pieces | 131 | 0.215 | 0.100 | **0.115** | **5.9%** | Pcons identical: TRUE | no |
| **G** | binary | consistency pieces (fit + dfbeta), narrow frame | 136 | 0.152 | 0.045 | **0.107** | **1.7%** | rate identical 136/136; NULL sets agree TRUE | no |
| **F** | binary | search-stage OR fit (glm.fit on vectors, SE from QR) | 2577 | 4.295 | 1.167 | **3.128** | **49.1%** | identical 2577/2577; kinds agree TRUE | YES |
| **H+G** | binary | whole consistency evaluator: lean extraction + lean pieces | 136 | 0.268 | 0.048 | **0.220** | **3.5%** | Pcons identical: TRUE (passed 5 vs 5; lean fallbacks 0) | no |
| **H alone** | binary | lean extraction, package pieces | 136 | 0.268 | 0.155 | **0.113** | **1.8%** | Pcons identical: TRUE | no |

Continuous search-stage fit (already lean since 0.3.3): 0.582 s per call over 4702 candidates -- no candidate.

## 3. Combined per fixture (per-replicate shape; machinery seconds are additive across stages)

| fixture | replicate s | search-stage machinery s | consistency-stage machinery s (whole evaluator) | total s | share | clears 30%? |
|---|---|---|---|---|---|---|
| survival | 3.139 | 1.509 | 0.413 | **1.922** | **61.2%** | YES |
| continuous | 1.963 | 0.000 | 0.200 | **0.200** | **10.2%** | no |
| binary | 6.371 | 3.128 | 0.220 | **3.348** | **52.5%** | YES |

The survival search-stage figure (A) is Stage 0's committed measurement on this machine (1.509 s of a 3.2112 s replicate); the survival replicate measured here is the denominator of the combined row.


## 4. Modification map (CRAN 0.4.0 line references; nothing here was changed)

| candidate | function(s) that would change | file:lines (0.4.0) | what replaces what | new helper | measured status |
|---|---|---|---|---|---|
| A | `fit_cox_for_subgroup()` unadjusted branch | R/subgroup_search.R:809-845 | formula `coxph()` -> `survival::coxph.fit()` with captured conventions; fallback to the formula block on error | `.fs_cox_lean_ci()` | Stage 0: 1410/1410 identical; 1.509 s/call |
| B | `.consistency_cox_pieces()` (and the warm-start fits in `evaluate_subgroup_consistency()` / `evaluate_consistency_twostage()`) | R/consistency_resample.R:55-91; R/subgroup_consistency_helpers.R:1457-1467, 1827-1837 | two formula `coxph()` + `residuals(type = 'dfbeta')` with model-frame re-evaluation -> `coxph.fit()` x2 + `residuals.coxph()` on a constructed object | shares `.fs_cox_lean_ci()`'s fit | pieces rate identical 120/120; machinery 0.182 s/call |
| E | `.consistency_glm_pieces()` MD branch + `.dfbeta_glm()` | R/consistency_resample.R:292, 168-176 | `lm()` + `lm.influence()` -> `lm.fit()` + hat values from the fit's QR | `.fs_glm_pieces_lean()` | pieces rate identical 131/131; machinery 0.083 s/call |
| F | routing predicate in `forestsearch()` + a fast OR closure + `fit_glm_for_subgroup_fast()` dispatch | R/forestsearch_main.R:2085-2104; R/glm_effect_estimators.R (new closure beside `.make_lm_estimator_fast`, :1058); R/subgroup_search.R:975-1003 | data-frame slice + `glm()` + `vcov()` -> vectors + `glm.fit()` + SE from the QR; existence check via the package's own `.fs_existence_reason()` | `.make_or_estimator_fast()`, `.fs_glm_two_group()` | search fits identical 2577/2577; machinery 3.128 s/call |
| G | `.consistency_glm_pieces()` OR branch + `.dfbeta_glm()` | R/consistency_resample.R:271-272, 168-176 | `glm()` + `lm.influence()` -> `glm.fit()` + hat values from the fit's QR | shares `.fs_glm_pieces_lean()` | pieces rate identical 136/136; machinery 0.107 s/call |
| H | member extraction in `evaluate_consistency_twostage()` and `evaluate_subgroup_consistency()` | R/subgroup_consistency_helpers.R:1803-1811 and 1431-1439 (+ index row extraction :1778, :1391) | `as.data.table(df)` + sequential row filters per candidate -> one precomputed logical column cache and an index matrix, subset only the columns the pieces need | none (hoisting into `subgroup.consistency()`) | Pcons identical, H alone: survival TRUE, continuous TRUE, binary TRUE; machinery (H alone) 0.209 / 0.115 / 0.113 s/call |
| C, D | none | -- | closed by Stage 0 measurement | -- | -- |

Count/IRR is not measured: its search fit and pieces are the binary code paths with `poisson(log)` and an offset; the binary figures are the proxy. RR/RD (multi-tier fallbacks), adjusted and propensity-weighted fits, DINA and GRF identifiers are outside every candidate above.

## 5. Reading

- survival: replicate 3.14 s; search-stage machinery 1.51 s, consistency-stage machinery 0.41 s; combined 61.2% (bar 30.0%: cleared). Consistency-stage identity: TRUE.
- continuous: replicate 1.96 s; search-stage machinery 0.00 s, consistency-stage machinery 0.20 s; combined 10.2% (bar 30.0%: not cleared). Consistency-stage identity: TRUE.
- binary: replicate 6.37 s; search-stage machinery 3.13 s, consistency-stage machinery 0.22 s; combined 52.5% (bar 30.0%: cleared). Consistency-stage identity: TRUE.
- Shares depend on each fixture's candidate flow (the consistency stage's cost is proportional to the number of candidates evaluated under full evaluation); the fixture sources above say which configuration each row describes.


## 3. Reading (numbers only)

All three fixtures are the committed configurations. Survival is Stage 0's fixture V (the vignette gbsg configuration). Continuous and binary are the ACTG175 application blocks transplanted in §1.1. **No shipped default was kept.**

- **survival** (replicate 3.139 s): **A alone clears** the bar. Against the Stage 0 replicate it is 47.0% (summary row A: 1.509 s of 3.2112 s). Against this run's replicate it is 1.509 / 3.139 = 48.1%. B (5.8%), H alone (6.7%) and H (+B) (13.2%) do not clear alone. Combined A + H (+B) is **61.2%**. Identity: A 1410/1410 (Stage 0); B rate identical 120/120 with NULL sets in agreement; Pcons identical TRUE for both H (+B) and H alone (16 vs 16 passed, 0 lean fallbacks). Consistency stage: 120 candidates in total, **120 evaluated under full evaluation, no early stop**. Rprof: 0.52 s of the 3.11 s call, of which 0.25 s is pieces and 0.27 s is the rest.
- **continuous** (replicate 1.963 s): **no candidate clears**, alone or in combination. E is 4.2%, H alone 5.9%, H+E 10.2%; the combined figure is **10.2%** (the search stage has no candidate, since the fit has been lean since 0.3.3: 0.582 s over 4702 candidates). Identity: E rate identical 131/131 with NULL sets in agreement; Pcons identical TRUE for both H+E and H alone (8 vs 8, 0 fallbacks). Consistency stage: 131 / **131 evaluated**, no early stop; Rprof puts it at 0.23 s of a 1.46 s call. Most of the call is the search stage (1.15 s). Of that, 0.59 s is the per-candidate fit, and about 0.56 s is enumeration and bookkeeping outside the per-candidate fit, which no candidate in the map addresses.
- **binary** (replicate 6.371 s): **F alone clears** at **49.1%** (3.128 s; identical 2577/2577; result kinds agree, with 2577 fits, 0 non-estimable and 0 NULL). G (1.7%), H alone (1.8%) and H+G (3.5%) do not clear. Combined F + H+G is **52.5%**. Identity: G rate identical 136/136; Pcons identical TRUE for both H+G and H alone (5 vs 5, 0 fallbacks). Consistency stage: 136 / **136 evaluated**, no early stop; Rprof puts it at 0.28 s of a 6.33 s call. The "GRF 1.63 s" bucket double-counts a single forest: `fit_causal_forest_glm` and its child `causal_forest` are each 0.82 s (12.9% of the profiled call). This is the variable-importance screening forest that `vi.grf.min = -0.2` runs (`R/forestsearch_main.R:3284`) even though `use_grf = FALSE`. It is outside every candidate.
- **Where the consistency stage's share comes from:** in all three fixtures every candidate that enters consistency is evaluated (total = evaluated, early stop FALSE), and there are only 120–136 of them, so the consistency stage is 0.23–0.52 s per call. The whole-evaluator machinery of 0.20–0.41 s is real but small next to the search stage. The large consistency-stage shares in chat's dry run came from 4032–4085 candidates under its approximate defaults (§1.1).

## 4. Modification map

The table is from the summary (§2, section 4), annotated. Every row is the **same computation with a bit-identity target**, and **every measured identity held**: A 1410/1410 (Stage 0), B 120/120, E 131/131, F 2577/2577, G 136/136, and H Pcons identical on all three fixtures, both alone and with lean pieces. No row changes behaviour, and no row moves code without also re-implementing it.

| candidate | category | functions touched: exported / internal | signature change? |
|---|---|---|---|
| A | same computation (bit-identity target; held) | `fit_cox_for_subgroup()`: internal | none |
| B | same computation (held) | `.consistency_cox_pieces()`: internal. The warm-start fits sit inside `evaluate_subgroup_consistency()` / `evaluate_consistency_twostage()`, which are **exported**, but their bodies change, not their signatures. | none |
| E | same computation (held) | `.consistency_glm_pieces()` and `.dfbeta_glm()`: internal | none |
| F | same computation (held) | the routing predicate in `forestsearch()` (**exported**; body only), a new internal closure `.make_or_estimator_fast()` and `.fs_glm_two_group()`, and `fit_glm_for_subgroup_fast()` dispatch (internal) | none |
| G | same computation (held) | `.consistency_glm_pieces()` and `.dfbeta_glm()`: internal | none |
| H | same computation (held); in the form the map describes, it also **moves code** (the cache is built once in `subgroup.consistency()`) | `evaluate_consistency_twostage()` and `evaluate_subgroup_consistency()` are **exported**, and so is `subgroup.consistency()` | **possibly.** If the cache is hoisted into `subgroup.consistency()`, it has to reach the exported evaluators, which needs a new argument; an optional argument defaulting to `NULL` would be additive. The alternative is to rebuild the cache inside each call. See OPEN ITEMS. |
| C, D | closed by Stage 0 | — | — |

## 5. Script deviations

None. All eight scripts ran as shipped. The only authored change is the §2.4 fixture transplant in `fixture_continuous.R` and `fixture_binary.R`. The driver's console output, which I redirected outside the bundle during the run, was moved into `logs/run_stage0b_console.log`. The two fixture diffs were saved to `logs/`.

## 6. OPEN ITEMS

1. **Binary GRF variable-importance forest.** About 0.82 s per call (12.9% of the profiled binary call) is `fit_causal_forest_glm()`, run because the document sets `vi.grf.min = -0.2` with `use_grf = FALSE`. It runs once per bootstrap replicate. No candidate in the map covers it, and it was not measured as one.
2. **Continuous search-stage remainder.** About 0.56 s of the 1.15 s continuous search stage lies outside the per-candidate fit: enumeration, `evaluate_combination_with_status` bookkeeping, and `future_lapply` overhead under the sequential plan. Nothing in this map measures it as a candidate, so continuous stays at 10.2%.
3. **H hoisting and exported signatures.** See §4. In addition, b2's lean evaluator builds its logical column cache once, outside the timed loop, so the H machinery figures exclude that one-time cost (one `== 1` per cut column per call).
4. **Mixed survival denominator.** Row A quotes Stage 0's 47.0% (denominator 3.2112 s). The combined row uses this run's replicate of 3.139 s, against which A alone is 48.1%.
5. **Replicate vs call.** For the GLM fixtures the sequential replicate is heavier than the single call: continuous 1.963 vs 1.544 s, binary 6.371 vs 5.450 s. Every share uses the replicate, as specified.
6. **Shipped defaults reversed the treatment coding** (arm 1 = 0). Chat's dry-run GLM figures describe a different candidate flow and should not be compared with these.
7. `dev/glm-continuous-sims/current_status.md` does not exist (as at Stage 0). I did not create it.
8. Out of scope by design: count/IRR, RR, RD, adjusted and propensity-weighted fits, DINA and GRF identifiers. The summary's note gives the binary paths as the proxy for count.

## 7. Verdict (ten lines)

1. Measurement only. No package file changed (machine check in the closing block).
2. All three fixtures are the committed configurations, and both ACTG175 fixtures reproduce their documents' recorded subgroups.
3. Survival: A alone removes 47–48% of a replicate; with B and H the total is 61.2%. Every identity held.
4. Binary: F alone removes 49.1% of a replicate (3.128 s, 2577/2577 identical); with G and H the total is 52.5%.
5. Continuous: nothing clears. All candidates together give 10.2%, because its search fit is already lean.
6. Under the application configurations the consistency stage is small: 120–136 candidates, 0.2–0.5 s per call.
7. The consistency-stage candidates (B, E, G, H) are each below 15% on every fixture, and every one is bit-identical.
8. H in hoisted form touches exported evaluators. Every other row is a body-only change to internal functions or to `forestsearch()`.
9. The largest unmapped cost found is the binary VI forest (≈0.82 s per call).
10. No recommendation is made. Whether to commission any change is Larry's decision.

## Closing block (post-conditions)


Run after commit `a41f43d7` (bundle + report); this block was then appended in a final report-only commit, after which the same three checks were re-run; their output is unchanged (only the report file changed, under dev/).

```
$ git diff --name-only d6613995c882d15a7168674bb7e9a7ba4ced3564 HEAD | grep -v '^dev/'
(no output)
$ git status --porcelain | grep -v '^?? ' | grep -v '^.. dev/'
(no output)
$ git status --porcelain | grep '^?? ' | grep -v ' dev/'
?? quarto/extreme_subgroups/fixed_random/
```

The third line lists only the pre-existing dirt recorded in §1. Commits, in order: `9fc39fb1` (task document alone), `a41f43d7` (bundle sources after §2.4, results/, logs/, this report), then the report-only commit adding this block.
