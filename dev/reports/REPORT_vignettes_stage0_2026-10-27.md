# REPORT — vignettes stage 0: Table 1 provenance, package facts, timing ladder (2026-10-27)

Task: `dev/tasks/cc_task_vignettes_stage0_2026-10-27.md` (downloaded as `~/Downloads/cc_task_vignettes_stage0_2026-10-27.md` — hyphenated name; no hyphen-stripped copy existed). Branch `feature/vignettes-0.4.1` from master `174778d0`; task commit `60d7ca33`. Record directory `dev/reports/vignettes_stage0_2026-10-27/`. All gates passed; nothing in `R/`, `vignettes/`, `man/`, `DESCRIPTION`, `NEWS.md`, `_pkgdown.yml`, `.Rbuildignore`, `tests/` was touched; fs-glms-interpretable was only read.

## Provenance

- Host `pop-os`, 128 cores (`nproc`), R 4.6.1 (2026-06-24), library `/home/larryleon/R/x86_64-pc-linux-gnu-library/4.6`. Full record in `provenance.txt`.
- Installed forestsearch `0.4.0.9000`, built `R 4.6.1; ; 2026-10-27 21:01:28 UTC; unix`. **GATE passed:** installed `NEWS.md` is byte-identical to the repo's (`cmp` silent).
- survival 3.8.12, grf 2.6.1, future 1.80.0, doFuture 1.4.0, quarto (R) 1.5.1, knitr 1.52, rmarkdown 2.32; quarto CLI 1.10.18 (on PATH in this shell and at the RStudio path).
- Pre-existing `git status --short`: clean. Branch/ancestry GATE passed (`174778d0` is an ancestor of HEAD).

## Inventory A — where Table 1 came from

- **Sibling repo** `~/Documents/GitHub/fs-glms-interpretable`, HEAD **`19c6ac41`** (expected `a772a2fe`; that commit is in its history two commits back — recorded, not a stop). Working tree clean.
- **No `_payloads/` subdirectories exist** under `quarto/gbsg/` or `quarto/actg175/`. Each document writes `<doc>_payload.rds` beside itself; those are what's listed below.

`quarto/gbsg/`:

| .qmd | title | html committed | payload (size; `built_at`; `forestsearch_version`) |
|---|---|---|---|
| `analysis_gbsg_mr.qmd` | GBSG: selection-adjusted bounds for the max N(ε) subgroup | yes | 674,893 B; 2026-09-23 23:53:10; 0.3.5.9000 |
| `analysis_gbsg_dina_mr.qmd` | GBSG, DINA: selection-adjusted bounds for the max N(ε) subgroup | yes | 3,540 B; 2026-09-23 22:55:20; 0.3.5.9000 |
| `analysis_gbsg_grf_mr.qmd` | GBSG, GRF: selection-adjusted bounds for the max N(ε) subgroup | yes | 33,587 B; 2026-09-23 22:55:55; 0.3.5.9000 |
| `analysis_gbsg_mr_pstar099.qmd` | GBSG: selection-adjusted bounds for the max N(ε) subgroup (p\* 0.99 variant) | yes | 674,727 B; 2026-09-24 00:04:27; 0.3.5.9000 |
| `analysis_gbsg_trio_summary.qmd` | GBSG: the three identifiers side by side | yes | 2,408 B; 2026-09-17 20:07:53; not carried |

`quarto/actg175/`:

| .qmd | title | html committed | payload |
|---|---|---|---|
| `analysis_actg175_binary_mr.qmd` | ACTG175: … odds-ratio path | yes | 813,266 B; 2026-09-23 23:54:32; 0.3.5.9000 |
| `analysis_actg175_binary_mr_pstar099.qmd` | same, p\* 0.99 | yes | 813,120 B; 2026-09-24 00:18:05; 0.3.5.9000 |
| `analysis_actg175_continuous_mr.qmd` | ACTG175 continuous: … mean-difference path | yes | 1,046,847 B; 2026-09-23 23:56:27; 0.3.5.9000 |
| `analysis_actg175_continuous_mr_pstar099.qmd` | same, p\* 0.99 | **no html, no payload** | — |

- **Which document produced which column** (located in each payload's `$table` / `$bounds`, matched to full precision):
  - **FS** ← `analysis_gbsg_mr.qmd` (chunk `fit`; payload rendered 2026-09-23 23:53, commit `a481b23`). Region `{er <= 0} & {pgr <= 26}`, N 75, corrected 1.27618, field 0.63356, field-s 0.80357, Bonferroni 0.50868 / 0.83959. This is the candidate the declaration-calibration task named (c1 = c2 = 1.00, p\* 0.90). **GATE passed.**
  - **DINA** ← `analysis_gbsg_dina_mr.qmd` (chunk `fit`; payload 2026-09-23 22:55, commit `bee0ab1`). Corrected 1.07609, field 0.63779, field-s 0.84206, Bonferroni 0.57247 / 0.88682.
  - **GRF** ← `analysis_gbsg_grf_mr.qmd` (chunk `fit`; payload 2026-09-23 22:55, commit `299fa11`). Corrected 1.11099, field 0.51047, field-s 0.84629, Bonferroni 0.40368 / 0.87440.
  - **No single document assembles all three.** `analysis_gbsg_trio_summary.qmd` (2026-09-17) does place them side by side, but its FS corrected values are 1.2806 / 0.6446, from an older render, so it is not Table 1's source. Its DINA/GRF rows equal the 09-23 payloads.
  - The DINA/GRF `.qmd` files were last committed 2026-09-17 (`8014cfa`), before their 09-23 re-renders. The FS `.qmd` was committed with its render. So all three committed sources are the code that produced the payloads.
- **Transplants** (`table1_calls/`): `fs_call.R`, `dina_call.R`, `grf_call.R`. Each contains chunks `label-map` + `setup` + `fit` verbatim, with a header giving source path, fence line ranges and sibling HEAD. `label-map` is included because the DINA/GRF `fit` chunks call `describe_subgroup()`. All three prepare data from `survival::gbsg` with no file read.
- **Settings table** (values as written in each `fit` chunk; "—" = not passed):

| setting | FS | DINA | GRF |
|---|---|---|---|
| identifier | `subgroup_method = "consistency"` | `"dina"`; `dina_select_statistic = "effect"`, `dina_args = list()` | `"grf"`; `grf_selection = "frontier"`, `grf_select_statistic = "effect"`, `grf_depth = 2L`, `dmin.grf = 0.0` |
| consistency method | `consistency_method = "resample"` (closed form, §2.1) | — (no consistency screen) | — |
| `fs.splits` | 1000 (inert, §2.1) | 1000 | 1000 |
| `pconsistency.threshold` | 0.90 | — (default 0.9; not applied) | — |
| `hr.threshold` | 1.00 | 1.00 | 1.00 |
| `hr.consistency` | 1.00 | — (default 1) | — |
| `sg_focus` | `"effMaxSG"` | `"effMaxSG"` | `"effMaxSG"` |
| `selection_rule` / `effect_neighborhood` | `"neighborhood"` / 0.20 | same | same |
| `n.min` / `d0.min` / `d1.min` / `maxk` | 60 / 10 / 10 / 2 | same | same |
| `conf.cont_jcuts` | `list(er = 10, pgr = 10)` | same | same |
| `conf_force` | `c("er <= 0", "pgr <= 0")` | same | same |
| `cut_type` / `cont.cutoff` / `collapse_cuts` | `"default"` / 4 / TRUE | same | same |
| `use_lasso` / `use_grf` / `use_dina` | FALSE / FALSE / FALSE | same | same |
| `m1.threshold` / `minp` / `stop_threshold` | Inf / 0.025 / NULL | same | same |
| `seedit` | 8316951 | 8316951 | 8316951 |
| `is.RCT` / `est.scale` | TRUE / `"hr"` | same | same |
| `mr_inference` | TRUE | TRUE | TRUE |
| `mr_inference_args` | `ci_method = "field"`, `draws = mr_draws` (`mr_draws <- 5000L`), `include_complement = TRUE`, `confirm_rule = "point"`, `field_uniform = FALSE`, `field_complement = TRUE`, `ij_residual = "two_term"`, `field_decompose = TRUE`, `field_scale_complement = "selected"`, `field_recovery = FALSE`, `return_reselection = TRUE`, `keep_declaration_field = TRUE`, `keep_field_matrix = FALSE`, `declaration_c0 = c(0.70, 0.75, 0.80, 0.85)` | same minus `keep_declaration_field`, `keep_field_matrix`, `declaration_c0` | same as DINA |
| field outer × inner | not passed → 1000 × 500 (not settable on the FS path, §2.2) | not passed → 1000 × 500 | not passed → 1000 × 500 |
| complement scaling / Bonferroni | `field_scale_complement = "selected"`; Bonferroni fixed in package (γ 0.025 per side) | same | same |
| `parallel_args` | — (package default `multisession`) | — | — |
| `max.minutes` / `vi.grf.min` / `max_n_confounders` | — / — / **1000** | — / — / 1000 | — / — / 1000 |
| `tune_grf` | — | — | — (default FALSE) |

- **Render wall-clock** (from `REPORT_applications_rerun_2026-09-24.md` §"Wall-clock per document": whole `quarto render`, fit in parentheses):
  - FS 54.0 s (44.9 s)
  - DINA 18.4 s (8.1 s)
  - GRF 32.9 s (25.2 s)

  The report states "Machine: Linux x86_64, R 4.6.1, reference BLAS 3.12.0, 1 core." The payloads' `meta$fit_sec` agree: FS 44.55 s (the 23:53 follow-up render), DINA 8.06 s, GRF 25.20 s. **Caveat:** no document passes `parallel_args`, so the package default `multisession` plan was in force at render. The payloads' `cores = 1L` is a literal written in the payload chunk, not a measurement. On the closed-form path the consistency screen consumes no RNG, and (per `R/subgroup_consistency_main.R:707-715`) its results are invariant to the plan, so this affects timing only, not the numbers.

## Inventory B — package facts (detail in `package_facts.md`)

- **§2.1 GATE passed.** `forestsearch(consistency_method = c("resample", "split"))`, default `"resample"` (`R/forestsearch_main.R:1686`, `:1718`). Despite its name, `"resample"` is the closed form `max(0, 2*pnorm((beta_hat - log c2)/sigma_D) - 1)` (`R/consistency_resample.R:455-457`; roxygen `R/forestsearch_main.R:1029-1045`: "Despite the name, "resample" performs no numerical resampling at evaluation time -- the rate is a closed-form expression").
  - Fallback to literal splitting happens per candidate only when the closed form is `NA` (`R/subgroup_consistency_helpers.R:1485-1500` fixed loop; `:1863-1907` two-stage, `resample_unsupported <- is.na(p.consistency)` at `:1890`). It happens for the whole run only if `consistency_resample()` is absent or a GLM path lacks a resample spec (`R/subgroup_consistency_main.R:427-438`).
  - `fs.splits` is used only on fallback.
  - The result records the method (`fit$args_call_all$consistency_method`) and algorithm (`fit$consistency_algorithm`) but **no per-candidate flag or fallback count**. The ladder measured it instead: **120 closed-form evaluations, 0 fallbacks** in every FS rung.
  - **Larry's clarification is confirmed from source**, with one wording point for the vignette. The Table 1 call selects the closed form explicitly by `consistency_method = "resample"`, so the vignette needs to explain that this name means the closed form.
- **§2.2 Formals.**
  - `formals_forestsearch.txt` has all 86 formals of `forestsearch()` with defaults, plus the formals of `fs_mr_inference()` (exported).
  - The documented `mr_inference_args` elements are listed in `package_facts.md`.
  - **`field_R_out` / `field_R_in` are not forwarded on the FS path** (`R/forestsearch_main.R:3828-3896`), but are forwarded on the DINA/GRF path (`R/fs_mr_inference_methods.R:178-179`).
- **§2.3 Element paths.** With `g <- fit$mr_inference`, `cc <- g$complement`, `fb <- g$field`, `fc <- fb$complement`, `js <- fb$joint_s`:

  | quantity | path |
  |---|---|
  | unadjusted HR and CI, Ĥ | `g$naive$est`, `$lower`, `$upper` |
  | unadjusted HR and CI, Ĥᶜ | `cc$naive$est`, `$lower`, `$upper` |
  | corrected, Ĥ | `g$debiased$est` |
  | corrected, Ĥᶜ | `cc$debiased$est` |
  | field lower bound | `fb$lower_1s` |
  | field-s upper bound | `fc$upper_1s_s` |
  | Bonferroni pair | `js$bonf_lower_H`, `js$bonf_upper_Hc` |
  | MR family size | `g$n_family` |
  | MR draws used | `g$settings$draws` |
  | field outer / inner | `fb$R_out`, `fb$R_in` |
  | N | `g$n_selected`, `cc$n` |

  - **Events are not stored.** They are computed from `fit$df.est$treat.recommend == 0` membership.
  - **The seed is not stored** in `g$settings`. It is `mr_inference_args$seed`, defaulting to `seedit`.
  - `print.forestsearch` / `summary.forestsearch` print the field and field-s bounds, the Bonferroni pair, the unadjusted Ĥ HR, N, Pcons and the IJ interval. **They do not print** the corrected estimates, the unadjusted CIs, the Ĥᶜ unadjusted HR, events, draw counts or the seed.
  - `fit$mr_inference` has no print method.
  - `mr_estimates_table()` exists but requires bootstrap results and has no field bounds.
- **`declaration_c0`** is reachable via `mr_inference_args = list(keep_declaration_field = TRUE, declaration_c0 = …)`, and it stores `g$declaration_field$Mstar_c0`. `fs_declaration_calibration(fit, c0 = …)` returns κ̂, the implied p\* screen and per-c0 `pstar_set` / `pstar_fine`. On S1 (500 draws): κ̂ 3.77, screen Pcons ≥ 0.895, 9 admitted; c0 0.70 → p\* 1 settable at 2 digits; 0.75/0.80/0.85 → none at 2 digits.
- **§2.4 Build facts.**
  - `DESCRIPTION` and `.Rbuildignore` are reproduced verbatim in `package_facts.md`; `_pkgdown.yml` (595 lines) is copied to Downloads.
  - `VignetteBuilder: knitr, quarto`.
  - Imports: `grf`, `policytree`. Suggests: `quarto`, `knitr`, `rmarkdown`, `speff2trial`.
  - CRAN vignettes are `forestsearch.Rmd` (knitr::rmarkdown) and `survival-post-selection.qmd` (quarto::html). `vignettes/articles` is build-ignored.
  - `articles/extreme_subgroups.qmd` has no `vignette:` block.
  - No test references `vignettes/` or `inst/doc`.

## The ladder

Single process, `plan = "sequential"`, installed build, transplanted calls verbatim plus the stated overrides (`parallel_args`, `quiet`, `details`, `plot.sg`, `plot.grf`, and `mr_inference` / `mr_inference_args$draws` per rung). Field 1000 × 500 in every MR rung (package default, and the transplanted calls' own). Wall = `/usr/bin/time -v` elapsed for the whole `Rscript` (R startup + data prep + fit); fit = `proc.time()` around the `forestsearch()` call. **No warnings in any rung.** This is expected, not a gap: `max_n_confounders = 1000` equals its default, and `max.minutes` / `vi.grf.min` are not passed.

| rung | id | MR draws | wall s | fit s | max RSS MB | region | N (ev) Ĥ / Ĥᶜ | unadj Ĥ (95% CI) | corr Ĥ | field lb | unadj Ĥᶜ (95% CI) | corr Ĥᶜ | field-s ub | Bonf |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| S0 | FS | — | 3.61 | 1.76 | 336 | er≤0 & pgr≤26 | 75 (41) / 611 (258) | — | — | — | — | — | — | — |
| S1 | FS | 500 | 27.83 | 25.94 | 565 | er≤0 & pgr≤26 | 75 (41) / 611 (258) | 2.222 (1.169, 4.221) | 1.277 | 0.634 | 0.605 (0.465, 0.788) | 0.652 | 0.812 | 0.509 / 0.849 |
| S1r | FS | 500 | 26.97 | 25.08 | 546 | same | same | 2.222 (1.169, 4.221) | 1.277 | 0.634 | 0.605 (0.465, 0.788) | 0.652 | 0.812 | 0.509 / 0.849 |
| S2 | FS | 2000 | 30.89 | 29.04 | 684 | same | same | 2.222 (1.169, 4.221) | 1.283 | 0.637 | 0.605 (0.465, 0.788) | 0.644 | 0.803 | 0.511 / 0.839 |
| S3 | FS | 5000 | 38.39 | 36.46 | 1,050 | same | same | 2.222 (1.169, 4.221) | **1.276** | **0.634** | 0.605 (0.465, 0.788) | **0.645** | **0.804** | **0.509 / 0.840** |
| D1 | DINA | 500 | 9.50 | 7.68 | 358 | grade3≥1 & pgr≤20 | 113 (61) / 573 (238) | 1.315 (0.746, 2.317) | 1.054 | 0.625 | 0.641 (0.489, 0.839) | 0.661 | 0.845 | 0.561 / 0.890 |
| G1 | GRF | 500 | 22.57 | 20.77 | 448 | er≤0 | 82 (45) / 604 (254) | 1.951 (1.037, 3.673) | 1.091 | 0.501 | 0.615 (0.472, 0.801) | 0.669 | 0.853 | 0.396 / 0.881 |
| D3 | DINA | 5000 | 9.73 | 7.93 | 437 | grade3≥1 & pgr≤20 | 113 (61) / 573 (238) | 1.315 (0.746, 2.317) | **1.076** | **0.638** | 0.641 (0.489, 0.839) | **0.658** | **0.842** | **0.572 / 0.887** |
| G3 | GRF | 5000 | 26.75 | 24.90 | 626 | er≤0 | 82 (45) / 604 (254) | 1.951 (1.037, 3.673) | **1.111** | **0.510** | 0.615 (0.472, 0.801) | **0.663** | **0.846** | **0.404 / 0.874** |

Cumulative ladder wall: 196.2 s. No rung near the 20-min cap. Machine-readable: `ladder.csv` (driver: wall, RSS, status) and `ladder_quantities.csv` (R: fit time, counts, quantities at 6 dp). The full comparison is in `comparison_table1.csv`.

**Search counters and p_cons (S0)**

- **GATE A passed:** region `{er <= 0} & {pgr <= 26}`, N 75, 41 events.
- Counters vs Table S13, side by side:

  | Table S13 | value | ladder field | value |
  |---|---|---|---|
  | enumerated | 2,211 | `filter_counts$n_evaluated` | 2,211 |
  | fittable | 1,410 | `n_passed_sample_size` = `n_passed_cox` | 1,410 |
  | clearing screen | 121 | `n_passed_hr` | 121 |
  | passing consistency | 9 | `grp.consistency$n_passed` | 9 |
  | frontier | 3 | **no counter in the result** | — |
  | re-selection family | 1,744 | `mr_inference$n_family` | 1,744 |

  - Intermediate counters: variance 2,069, prevalence 2,069, redundancy 1,943, events 1,454.
  - Definitional note: 121 cleared the HR screen but the consistency stage evaluated 120 (`n_candidates_evaluated` = `n_candidates_total` = 120). One screened candidate is removed before consistency. I did not trace which one or why (consistency-stage duplicate removal is the likely mechanism; not verified).
  - 8 of the 120 raw closed-form rates are ≥ 0.90. 9 pass because Pcons is rounded to `pconsistency.digits` = 2 first (0.8987 → 0.90).
  - The 9 passing candidates with Pcons are in `postprocess_log.txt`. The full 120-row p_cons table is in `logs/S0.log` and `rungs/rung_S0.rds$pcons_closed`.
- **Per-candidate fallback flag:** all 120 are closed form; 0 `NA`; 0 calls to `.consistency_via_splits`.

**Comparisons**

- **S1 vs S1r:** `all.equal` TRUE at tolerance 0 on the region, N/events and every Table 1 quantity. It is also TRUE on the whole `mr_inference` object with timings stripped. **Deterministic.**
- **S3 vs Table 1, FS column:** all nine quantities equal the published 3-dp values and **equal the committed payload to full precision** (abs diff 0). Log-scale differences against the 3-dp published values are rounding only (|Δ| ≤ 0.0007). So the installed 0.4.0.9000 reproduces the 0.3.5.9000 render exactly, seed and draw layout included.
- **D3 / G3 vs the DINA / GRF columns:** region, N, events and all nine quantities exact against the payloads (abs diff 0). D1 and G1 select the same regions; only the corrected estimates and bounds move with the draw count.
- **Reduced draws vs published** (S1 / D1 / G1 at 500, S2 at 2000), largest |log diff| on the MR quantities:
  - S1: 0.010 (Ĥᶜ corrected 0.652 vs 0.645; field-s 0.812 vs 0.804; Bonferroni Ĥᶜ 0.849 vs 0.840)
  - S2: 0.005
  - D1: 0.021 (corrected 1.054 vs 1.076; field 0.625 vs 0.638; Bonferroni Ĥ 0.561 vs 0.572)
  - G1: 0.020 (corrected 1.091 vs 1.111; Bonferroni Ĥ 0.396 vs 0.404)

  Unadjusted quantities, regions, N and events do not change with draws. At 500 draws the third decimal of most MR quantities differs from Table 1. At 2000 (FS), all six MR quantities (6 of 12) still differ at the third decimal.

## Budget arithmetic

| bundle | rungs | process wall | fit only |
|---|---|---|---|
| **full** (manuscript draws) | S3 + D3 + G3 | **74.9 s** | 69.3 s |
| **reduced** (500 draws) | S1 + D1 + G1 | **59.9 s** | 54.4 s |
| search alone (FS) | S0 | 3.6 s | 1.8 s |

- Reference points: both current vignettes together 45 s; `devtools::check()` 11 m 53 s, tests 403 s. The full bundle adds about 75 s of single-core compute (plus rendering) to a vignette build, under 11% of today's check time.
- **Draw count is not the cost driver.** Going from 500 to 5000 multiplier draws costs about 15 s across the three identifiers (MR stage: FS 0.56 → 4.56 s; GRF 0.21 → 2.39 s; DINA 0.02 → 0.36 s). The field's Monte Carlo is the fixed cost: FS 12.9 s Ĥ + 2.8 s Ĥᶜ, GRF 12.1 + 1.5 s, DINA 6.7 + 0.1 s. That is about 36 s of the full bundle, and it is unchanged by `draws`.
- On the FS path the field size cannot be reduced through `forestsearch()` (§2.2). A `params` switch on `draws` therefore saves 20% at most, and costs exact agreement with Table 1.

## Script deviations (verbatim record)

- **Ladder data rows.** The driver writes `ladder.csv` (rung, exit code, status, wall, `/usr/bin/time` elapsed, max RSS kB, cumulative). The per-rung quantities go to a separate `ladder_quantities.csv` written by `ladder.R`, rather than one merged CSV. Per-rung `.rds`: `rungs/args_<rung>.rds` (resolved arguments minus the data frame, the verbatim call and the evaluated call) and `rungs/rung_<rung>.rds` (the record). The full fit objects are in `fits/fit_<rung>.rds`.
- **Candidate labels missing from the traced p_cons table.** The two-stage evaluator does not have `this.m_label` / `m` under those names in the frame `consistency_resample()` is called from, so the `m` / `label` columns of `pcons_closed` are `NA`. n, d, β̂, σ_D and the rate are recorded. The labels of the 9 passing candidates come from `grp.consistency$out_sg$result`. Not re-run: the ladder's run budget was spent, and labels are not needed for the fallback count.
- **Fallback detection is by tracing.** `trace()` was put on `subgroup.consistency`, `consistency_resample` and `.consistency_via_splits` (exit tracers that only record). Overhead is negligible, but it is inside the timed region.
- **Post-processing script.** `postprocess.R` was added (reads saved fits; runs no fit) for §2.3 and the comparisons. Its log is `postprocess_log.txt`, about 350 kB, mostly the raw print of `fit$mr_inference`, kept as the evidence for "no print method".
- **`summary(fit)` printed twice in that log** (`summary.forestsearch` prints and then my `print(s)` printed again). This is a script artifact, not a package fact.
- **`/usr/bin/time` was present**, so the in-R `system.time` fallback was not needed. In-R `proc.time()` and `gc()` max-used are in each log too.

## OPEN ITEMS (findings only; nothing changed)

1. **`field_R_out` / `field_R_in` are not forwarded by `forestsearch()` on the consistency (FS) MR path** (`R/forestsearch_main.R:3828-3896`). They are forwarded on the DINA/GRF path (`R/fs_mr_inference_methods.R:178-179`). The `mr_inference_args` roxygen documents neither. The field is the dominant per-fit cost, so this decides whether a reduced-cost vignette mode is possible for FS. The FS document's own settings table states "package defaults, not settable through `forestsearch()`".
2. **Default `ci_method` differs by path:** `"field"` on the FS path (`R/forestsearch_main.R:3854`) vs `"ij"` in `.fs_apply_mr()` for DINA/GRF (`R/fs_mr_inference_methods.R:171`, `.g(mr_inference_args$ci_method, "ij")`). The `mr_inference_args` roxygen says `"field"` is the default. The Table 1 calls pass `"field"` explicitly, so they are unaffected.
3. **No per-candidate closed-form/fallback record in the result.** A vignette sentence like "no candidate fell back to splitting" cannot be read off the fit today. Tracing was needed.
4. **The MR seed is not stored** in `fit$mr_inference$settings`, and **events in Ĥ/Ĥᶜ are not stored**. A Table 1 assembled from the fit has to recompute events from `df.est` and take the seed from `args_call_all$seedit`.
5. **No exported helper assembles the Table 1 row.** `print` / `summary` omit the corrected estimates and unadjusted CIs; `mr_estimates_table()` requires bootstrap results. `fit$mr_inference` has no print method (raw-list print of about 3,000 lines).
6. **`print` / `summary` show `sg_focus: hrMaxSG`** for a call made with `effMaxSG` (normalized name), and `Admission set: effect floor 0; consistency floor 0 at p* = 0.9` (floors on the log scale). Both are worth checking for vignette readability.
7. **Table S13's "3 on the frontier" has no counter in the result object.**
8. **One screened candidate is dropped between the HR screen (121) and consistency evaluation (120)**, mechanism not traced. This bears on how the vignette states the funnel against Table S13.

## Commit and Downloads

- Committed (see close-out commit): `ladder.R`, `run_ladder.sh`, `postprocess.R`, `logs/` (9), `logs_driver.txt`, `ladder.csv`, `ladder_quantities.csv`, `comparison_table1.csv`, `postprocess_log.txt`, `provenance.txt`, `package_facts.md`, `formals_forestsearch.txt`, `table1_calls/` (3), `rungs/` (18 `.rds`, all < 5 kB), `fits/` (9 `.rds`, 45–350 kB, all under 1 MB and not gitignored), and this report. Nothing withheld.
- `ls -la ~/Downloads/vignettes_stage0_2026-10-27/` — see below.

`.Rbuildignore` is copied as `Rbuildignore.txt`, because a dotfile is hidden in a file picker. `methodology.*` and `paper_simulations.*` each match only the `.Rmd`. The three sibling `.qmd` files are copies; the originals are untouched.

```
total 484
drwxrwxr-x 2 larryleon larryleon   4096 Oct 27 20:32 .
drwxr-xr-x 4 larryleon larryleon  20480 Oct 27 20:32 ..
-rw-r--r-- 1 larryleon larryleon  34457 Oct 27 20:32 analysis_gbsg_dina_mr.qmd
-rw-r--r-- 1 larryleon larryleon  35307 Oct 27 20:32 analysis_gbsg_grf_mr.qmd
-rw-r--r-- 1 larryleon larryleon  49636 Oct 27 20:32 analysis_gbsg_mr.qmd
-rw-rw-r-- 1 larryleon larryleon   5009 Oct 27 20:32 comparison_table1.csv
-rw-r--r-- 1 larryleon larryleon   2577 Oct 27 20:32 DESCRIPTION
-rw-rw-r-- 1 larryleon larryleon   7372 Oct 27 20:32 dina_call.R
-rw-r--r-- 1 larryleon larryleon  41211 Oct 27 20:32 forestsearch.Rmd
-rw-rw-r-- 1 larryleon larryleon   3307 Oct 27 20:32 formals_forestsearch.txt
-rw-rw-r-- 1 larryleon larryleon   6905 Oct 27 20:32 fs_call.R
-rw-rw-r-- 1 larryleon larryleon   7444 Oct 27 20:32 grf_call.R
-rw-rw-r-- 1 larryleon larryleon    449 Oct 27 20:32 ladder.csv
-rw-rw-r-- 1 larryleon larryleon   2129 Oct 27 20:32 ladder_quantities.csv
-rw-r--r-- 1 larryleon larryleon  52722 Oct 27 20:32 methodology.Rmd
-rw-rw-r-- 1 larryleon larryleon  16892 Oct 27 20:32 package_facts.md
-rw-rw-r-- 1 larryleon larryleon 106952 Oct 27 20:32 paper_simulations.Rmd
-rw-r--r-- 1 larryleon larryleon  18403 Oct 27 20:32 _pkgdown.yml
-rw-r--r-- 1 larryleon larryleon    577 Oct 27 20:32 Rbuildignore.txt
-rw-rw-r-- 1 larryleon larryleon  23074 Oct 27 20:32 REPORT_vignettes_stage0_2026-10-27.md
-rw-r--r-- 1 larryleon larryleon   9218 Oct 27 20:32 survival-post-selection.qmd
```
