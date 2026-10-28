# REPORT — post-selection vignette stage 0: Table S2 provenance, package facts, MR and bootstrap timing (2026-10-28)

Task: `dev/tasks/cc_task_postsel_stage0_2026-10-28.md` (committed `51fc55cc`). Record: `dev/reports/postsel_stage0_2026-10-28/`.
**Outcome:** every gate passed (after the harness fix in Deviation 1, which was made before GATE A was ever evaluated). GATE A passes at M-FS. All seven rungs completed in 109 s total. Every naive and MR cell of Table S2 reproduces at 2 dp for all three identifiers.

## Provenance

- **Repo and branch:** host pop-os, branch `feature/vignettes-0.4.1`.
  - HEAD at task start: `d8756010`; `99e8a64a` is an ancestor (GATE 0.1 pass).
  - HEAD after the task-doc commit: `51fc55cc`.
- **Pre-existing untracked paths** (not staged): `vignettes/articles/methodology.html`, `vignettes/articles/paper_simulations.html`.
- **Install check (GATE 0.3 pass):**
  - Installed forestsearch is `0.4.0.9000`, built 2026-10-27 21:01:28 UTC.
  - `git diff 936b6734 HEAD -- R/ DESCRIPTION NAMESPACE inst/` is empty (0 lines).
- **Environment:** R 4.6.1, survival 3.8.12, grf 2.6.1, future 1.80.0, doFuture 1.4.0, gt 1.4.0, nproc 128. Full detail in `provenance.txt`.

## §1 Inventory A — where Table S2 came from

### Manuscript repo

- Found at `~/Documents/GitHub/fs-post-selection` (`~/GitHub/fs-post-selection` does not exist). It is a git repo: HEAD `bc82bc32c8a6ef66a6fc252992278f102b44b819` (2026-09-15, "Create ujag122_supplemental_files.zip"), and `git status --short` is clean (GATE pass).
- **Fragment:** `manuscript/_supp_gbsg_table.qmd`, chunk `tbl-gbsg-multimethod`.
  - It is included at `postselection_subgroup_inference_supplement.qmd:1126` (`{{< include _supp_gbsg_table.qmd >}}`).
  - The caption is built as `cap_tbl`: "De-biased subgroup hazard ratios for the GBSG trial; parenthetical values are 95% confidence intervals. Here $\widehat{H}$ is ...".
  - It reads `P$table` from `_payloads/analysis_gbsg_survival_multimethod_payload.rds` and renders it with kableExtra. Columns: naive_* → Naive, fb_* → Full bootstrap, mr_* → Multiplier.
  - Copied to `manuscript_source/_supp_gbsg_table.qmd`, together with the include excerpt and a copy of the payload.
- **Payload** `manuscript/_payloads/analysis_gbsg_survival_multimethod_payload.rds`:
  - md5 `98be34fc565aba9ca94e58d12f51726f`; committed in fs-post-selection `49290de` (2026-08-17).
  - `built_at` 2026-08-15 22:50:30 PDT; `forestsearch_version` "0.2.0".
  - Machine: Mac-Studio-3.local, Apple M4 Max, 14 physical cores, 11 workers used, R 4.5.2, Accelerate BLAS.
  - `meta$B` 1000, `meta$mr_draws` 5000.
- **Manifest entry** (`manuscript/payload_manifest.qmd:78–83`, file last changed `d7b4c2b` 2026-08-17; copied to `manuscript_source/manifest_excerpt.txt`):
  - Status: "keep - post-audit (built 08-15)".
  - Generator: "forestsearch quarto/applications/gbsg/analysis_gbsg_survival_multimethod.qmd, chunk export-section6-2-payload (writes <stem>_payload.rds; export machinery cf4d6432, latest run 43b051b6)".
  - The grep hits for `supp_gbsg` and `De-biased subgroup hazard ratios` outside these files are archives (`_archive/`, `manuscript_legacy/`), survey reports, and the rendered `.tex`.

### Producing document (GATE pass)

- **Document:** forestsearch `quarto/applications/gbsg/analysis_gbsg_survival_multimethod.qmd` as of **`cf4d6432`** (2026-08-15 21:12), the last commit to touch the .qmd before the 22:50 payload.
  - The rendered html was committed in **`43b051b6`** (2026-08-16 12:07, "update gbsg and actg175 analysis runs"). It shows "August 15, 2026" and "forestsearch 0.2.0, read at render time".
  - At HEAD the .qmd differs from `cf4d6432` only at lines 144–145 (`run_cv` / `run_loo` FALSE → TRUE, `6ad5a724`). In the 0.2.0 render K-fold and LOO were therefore off.
- **Payloads in forestsearch carrying the same table:**
  - The 0.2.0 payloads (`9b8d92ae`, built 2026-08-13 02:39, and `19f0a5bb`, built 2026-08-15 14:18, both Mac) give a `$table` **`identical()`** to the manuscript payload.
  - The `_payloads_2026-08-31` (0.3.1), `_payloads_2026-09-01` and `_payloads_2026-09-01_complete` (0.3.5) payloads (pop-os) agree under `all.equal()` but are not identical.
  - The manuscript's copy itself is the 22:50 render, and its forestsearch-side file was deleted from tracking in `6769cfd0`.
- **Bootstrap source:** the B = 1,000 bootstrap was run **in the same render**, not separately. Chunks `bootstrap`, `fs-dina-bootstrap`, `fs-grf-bootstrap`, with `forestsearch_bootstrap_dofuture(mr_in_replicates = FALSE, nb_boots = NB = 1000, multisession, 11 workers)`. The payload's `fb_*` columns are `boot$H_estimates$H2/H2_lower/H2_upper` and the `Hc_estimates` equivalents.
- **Row mapping in Table S2:**
  - FS = the main fit `fs` (`subgroup_method = "consistency"`, chunk `forestsearch-main`).
  - DINA = `fs_dina` (`subgroup_method = "dina"`, chunk `fs-dina-select`).
  - GRF = `fs_grf` (`subgroup_method = "grf"`, chunk `fs-grf-select`).
  - The DINA-screening fit (A) and the standalone DINA (B) are not in the table.

### Transplants (`table_s2_calls/`, verbatim, header with path, chunks and commit)

- **`data_prep.R`:** chunks `setup` (92–99), `parallel-setup` (106–270), `load-packages` (314–331), `data-setup` (498–522), `dina-config` (1233–1243). `setup` was added after attempt 1; see Deviation 1.
- **`fs_call.R`:** `forestsearch-main` (687–809).
- **`dina_call.R`:** `fs-dina-select` (1499–1544).
- **`grf_call.R`:** `fs-grf-select` (1740–1787).
- **`bootstrap_calls.R`:** `bootstrap` (931–963), `fs-dina-bootstrap` (1563–1597), `fs-grf-bootstrap` (1802–1832).
- **`table_assembly.R`:** `bootstrap-summary`, `debias-MR-table`, and their DINA/GRF twins, plus `export-section6-2-payload` (2434–2641; `.engine_rows()` builds the Table S2 rows).

### Settings table (Table S2 document @ `cf4d6432` vs Table 1 calls of stage 0 2026-10-27, `table1_calls/`)

| identifier | argument | Table S2 doc | Table 1 (2026b) |
|---|---|---|---|
| FS | `effect_neighborhood` | 0.10 | 0.20 |
| FS | `mr_inference_args` | `list(draws = 5000, family_native_neighborhood = NULL)` | `list(ci_method = "field", draws = 5000, include_complement = TRUE, confirm_rule = "point", field_uniform = FALSE, field_complement = TRUE, ij_residual = "two_term", field_decompose = TRUE, field_scale_complement = "selected", field_recovery = FALSE, return_reselection = TRUE, keep_declaration_field = TRUE, keep_field_matrix = FALSE, declaration_c0 = c(0.70, 0.75, 0.80, 0.85))` |
| FS | `ci_method` resolved | not passed → **"ij" at 0.2.0** (doc line 993) / **"field" at HEAD** | "field" |
| FS | other args passed only by the S2 doc | `use_twostage = TRUE`, `twostage_args = list()`, `max.minutes = 3`, `vi.grf.min = -0.2` (deprecated; warns at HEAD), `grf_depth = 2L`, `dmin.grf = 4`, `frac.tau = 0.80`, `return_selected_cuts_only = TRUE` (GRF args inert with `use_grf = FALSE`), `show_candidate_summary = TRUE`, `details = TRUE`, `plot.sg = TRUE`, multisession 11 workers | not passed (formal defaults) |
| FS | same in both | `conf_force = c("er <= 0", "pgr <= 0")`, `conf.cont_jcuts = list(er = 10, pgr = 10)`, `n.min 60 / d0.min 10 / d1.min 10 / maxk 2`, `hr.threshold 1.0`, `hr.consistency 1.0`, `pconsistency.threshold 0.90`, `consistency_method "resample"`, `fs.splits 1000`, `sg_focus "effMaxSG"`, `selection_rule "neighborhood"`, `stop_threshold NULL`, `seedit 8316951`, `use_lasso/use_grf/use_dina FALSE` | same |
| DINA | `effect_neighborhood` | 0.10 | 0.20 |
| DINA | `hr.threshold` | **0.90** (`fs_dina_hr_threshold`) | 1.00 |
| DINA | `dina_args` | `list(seed = 8316951, n_folds = 2)` | `list()` |
| DINA | `mr_inference_args` / resolved `ci_method` | `list(draws = 5000, family_native_neighborhood = NULL)` / **"ij"** (`.fs_apply_mr` default) | full list as FS / "field" |
| DINA | `conf_force`, `conf.cont_jcuts`, `max_n_confounders`, `cut_type`, `cont.cutoff`, `collapse_cuts`, `use_*`, `d0.min`, `d1.min`, `maxk`, `m1.threshold`, `minp`, `fs.splits`, `stop_threshold`, `est.scale` | **not passed** (HEAD defaults: `conf_force NULL`, `conf.cont_jcuts NULL`, others as Table 1) | passed (as FS) |
| DINA | `n_family` (MR) | 84 | 30 |
| GRF | `effect_neighborhood` | 0.10 | 0.20 |
| GRF | `mr_inference_args` / resolved `ci_method` | as DINA / **"ij"** | full list / "field" |
| GRF | `hr.threshold` | **not passed** (HEAD default 1.25; the doc's note says it does not apply in GRF mode) | 1.00 |
| GRF | `conf_force`, `conf.cont_jcuts` and the cut-grid args | not passed | passed |
| GRF | same in both | `grf_selection "frontier"`, `grf_select_statistic "effect"`, `grf_depth 2L`, `dmin.grf 0.0`, `frac.tau 0.8`, `sg_focus "effMaxSG"`, `n.min 60` | same |
| GRF | `n_family` (MR) | 858 | 858 |
| bootstrap (all three) | — | `forestsearch_bootstrap_dofuture(mr_in_replicates = FALSE, fs.est = <fit>, nb_boots = 1000, show_three = FALSE, details = FALSE (FS) / TRUE (DINA, GRF), parallel_args = list(plan = "multisession", workers = 11, show_message = TRUE))`; `seed` not passed → default `8316951L` | no bootstrap |

**Recorded document timings** (0.2.0, Mac Studio M4 Max; html `43b051b6` Computational Timing table and payload `meta$timings`):

| identifier | fit (s) | MR loop (s) | FB B=1000 bootstrap (s), 11 workers |
|---|---|---|---|
| FS | 7.4 | 0.38 | 605.9 |
| DINA | 0.5 | 0.11 | 71.8 |
| GRF | 3.1 | 0.13 | 225.7 |

Render total: 920.6 s. At 0.2.0 the FS fit did not run the field block, which is why its fit time (7.4 s) is far below the 31.8 s measured here.

## §2 Package facts (`package_facts_s2.md`)

- **Bootstrap call:** `forestsearch_bootstrap_dofuture(fs.est, nb_boots, seed = 8316951L, details = FALSE, show_three = FALSE, parallel_args = list(), digits = 4, mr_in_replicates = FALSE)`.
  - It needs `df.est`, `args_call_all` and an identified subgroup.
  - It accepts DINA and GRF fits: `bootstrap_results()` keeps the raw confounders for `"dina"` / `"grf"` (`bootstrap_analysis_dofuture.R:384–395`).
  - Each replicate re-runs `do.call(forestsearch, args_call_all + overrides)` from scratch, with the front-end caches nulled and MR stripped under the default (line 607).
  - Per-replicate correction: `H_biasadj_2 = 2*H_obs - (H_star + Hstar_star - Hstar_obs)`. Across replicates, `get_dfRes` takes the mean and computes the IJ SE with the finite-B correction. If `varhat <= termc` it falls back to the uncorrected variance.
  - Returns an `fs_bootstrap` object. The corrected estimate and CI are at `boot$H_estimates$H2/H2_lower/H2_upper` and the `Hc_estimates` equivalents; formatted strings are at `boot$SG_CIs$H_bc` / `$Hc_bc`.
- **`mr_estimates_table(fs, boot_results, est.scale = "hr", digits = 2)`** returns a **gt** table with the naive column (from the FB template) and the MR de-biased Ĥ/Ĥᶜ. It does **not** carry the FB corrected cells. Those are in `summarize_bootstrap_results(sgharm, boot_results, create_plots = FALSE, est.scale = "hr", digits = 2)$table`, a separate gt table with the same layout. No package function produces the three-column Table S2 layout; the document built it locally (`.engine_rows()`).
- **MR paths:**
  - Two-sided intervals: `g$debiased$est/lower/upper` and `g$complement$debiased$est/lower/upper`, built from the IJ SE under both `"ij"` and `"field"` (`fs_mr_inference.R:770, 854`). Under `"wald"` they use the robust SE.
  - Field-only extras live in `g$field` (`lower_1s`, `lower_2s/upper_2s`, `est2`, `complement$upper_1s_s`, `joint_s`).
- **Default change since 0.2.0:** the consistency branch now defaults `ci_method` to `"field"` (`forestsearch_main.R:3854`). The DINA and GRF branches keep `"ij"` by design (`fs_mr_inference_methods.R:137–139, 171`).
- **Rendering:**
  - gt tables render as self-contained inline HTML; the 0.2.0 render printed both.
  - They are `gt_tbl` objects, so they cannot be passed to `knitr::kable`. A kable version would be rebuilt from the numeric paths.
  - kableExtra (used by the manuscript fragment) is in neither Imports nor Suggests.
  - Nothing found that would break a self-contained HTML vignette.

## §3 The ladder

Single process, `plan = "sequential"`, installed build. Each rung ran as `timeout 1200 /usr/bin/time -v Rscript ladder.R <rung>`.

- **M rungs** use the transplanted calls verbatim plus these overrides: `parallel_args = list(plan = "sequential")`, `quiet = TRUE`, `details = FALSE`, `plot.sg = FALSE`, `plot.grf = FALSE`.
- **B rungs** use the document's bootstrap call verbatim with `nb_boots = 5`, `parallel_args = list(plan = "sequential")` and `details = FALSE`, on the saved M fit. The seed is the document's: the default `8316951L`.

`ladder.csv`:

| rung | status | wall (s) | peak RSS (kB) | cumulative (s) |
|---|---|---|---|---|
| M-FS | completed | 33.71 | 798,856 | 33.71 |
| M-DINA | completed | 3.55 | 424,728 | 37.26 |
| M-GRF | completed | 13.13 | 632,724 | 50.39 |
| M-FS-ij | completed | 17.80 | 770,700 | 68.19 |
| B-FS | completed | 12.73 | 403,264 | 80.92 |
| B-DINA | completed | 6.04 | 389,940 | 86.96 |
| B-GRF | completed | 21.75 | 471,540 | 108.71 |

### M rungs (`ladder_quantities_M.csv`; fit time = the `forestsearch()` call alone)

| rung | ci_method | fit (s) | region | N/E Ĥ | N/E Ĥᶜ | n_family | MR loop (s) | field Ĥ / Ĥᶜ (s) |
|---|---|---|---|---|---|---|---|---|
| M-FS | field (resolved default) | 31.83 | `{er <= 0} & {size <= 35}` | 61 / 34 | 625 / 265 | 1744 | 4.68 | 13.00 / 2.76 |
| M-DINA | ij (resolved default) | 1.73 | `{grade3 >= 1} & {pgr <= 10}` | 89 / 50 | 597 / 249 | 84 | 0.49 | — |
| M-GRF | ij (resolved default) | 11.28 | `{er <= 0}` | 82 / 45 | 604 / 254 | 858 | 2.41 | — |
| M-FS-ij | ij | 15.98 | `{er <= 0} & {size <= 35}` | 61 / 34 | 625 / 265 | 1744 | 4.62 | — |

**GATE A (M-FS): PASS.**
- Region `{er <= 0} & {size <= 35}`, N 61.
- Naive 2.54 (1.25, 5.17).
- Multiplier 1.44 (0.62, 3.36).

**M-FS field bounds** (Ĥ):
- `lower_1s` 0.699; two-sided `lower_2s` / `upper_2s` 0.573 / 3.096; `est2` 1.459.
- Ĥᶜ studentized `upper_1s_s` 0.805.
- Bonferroni joint pair 0.573 (Ĥ) / 0.835 (Ĥᶜ).

**M-FS-ij vs M-FS:**
- `all.equal(..., tolerance = 0)` on `g$debiased` (est, lower, upper, se_ij), `g$complement$debiased` (same fields) and `g$naive` is TRUE, and `identical()` is TRUE.
- Wall-clock: dropping the field block saves 15.85 s of fit time (15.91 s rung wall), almost all of the field blocks' 13.00 + 2.76 s.

### Table S2 comparison (`comparison_table_s2.csv`; payload = 0.2.0 render, now = HEAD)

| method | SG | n payload / now | Naive (now) | 2 dp | max abs diff | Multiplier (now) | 2 dp | max abs diff | Full bootstrap: Table S2 (B=1000) | B=5 here |
|---|---|---|---|---|---|---|---|---|---|---|
| FS | Ĥ | 61 / 61 | 2.536918 (1.245398, 5.167788) | ✓ | 9e-16 | 1.443947 (0.620609, 3.359577) | ✓ | **1.8e-3** | 1.94 (0.89, 4.21) | 2.68 (0.39, 18.31) |
| FS | Ĥᶜ | 625 / 625 | 0.607938 (0.468596, 0.788715) | ✓ | 1e-16 | 0.644104 (0.401602, 1.033039) | ✓ | **1.0e-3** | 0.63 (0.43, 0.93) | 0.69 (0.006, 78.90) |
| DINA | Ĥ | 89 / 89 | 1.599358 (0.879607, 2.908058) | ✓ | 1e-16 | 1.114573 (0.452056, 2.748052) | ✓ | 1e-15 | 0.89 (0.38, 2.07) | 1.10 (0.72, 1.69) |
| DINA | Ĥᶜ | 597 / 597 | 0.629072 (0.482280, 0.820542) | ✓ | 1e-16 | 0.659100 (0.399799, 1.086577) | ✓ | 2e-16 | 0.66 (0.40, 1.11) | 0.77 (0.009, 63.90) |
| GRF | Ĥ | 82 / 82 | 1.951393 (1.036610, 3.673451) | ✓ | 2e-15 | 1.053848 (0.474586, 2.340133) | ✓ | 4e-16 | 1.04 (0.49, 2.24) | 1.48 (0.78, 2.80) |
| GRF | Ĥᶜ | 604 / 604 | 0.614995 (0.472315, 0.800778) | ✓ | 0 | 0.660887 (0.410252, 1.064643) | ✓ | 2e-16 | 0.66 (0.39, 1.13) | 0.71 (0.009, 54.14) |

- **Recorded, not gated:**
  - DINA matches Table S2: `{grade3 >= 1} & {pgr <= 10}`, N 89, 1.11 (0.45, 2.75).
  - GRF matches Table S2: `{er <= 0}`, N 82, 1.05 (0.47, 2.34).
  - All naive cells, and the DINA and GRF MR cells, agree with the payload to machine precision.
  - The FS MR cells (both regions) agree at 2 dp but differ from the 0.2.0 payload at the third decimal. FS `n_family` (1744), selection and naive are unchanged.
- **The B = 5 full-bootstrap column is not a reproduction and is not comparable.** The package warns "nb_boots < 100 may produce unreliable confidence intervals". The Ĥᶜ intervals at B = 5 are degenerate: lower bounds near 0.01, upper bounds of 54–79. All 5 replicates identified a subgroup for each identifier, so the call runs for DINA and GRF.

### Bootstrap rungs (`ladder_quantities_B.csv`, `budget*.csv`)

| rung | call (s) | total / 5 (s) | replicate 1 (s) | per-replicate iteration times (s) | Σ iterations (s) | fixed overhead (s) |
|---|---|---|---|---|---|---|
| B-FS | 10.93 | 2.19 | 1.28 | 1.28 / 1.77 / 2.29 / 1.50 / 1.53 | 8.37 | 2.56 |
| B-DINA | 4.25 | 0.85 | 0.27 | 0.27 / 0.33 / 0.58 / 0.37 / 0.10 | 1.66 | 2.59 |
| B-GRF | 19.94 | 3.99 | 3.44 | 3.44 / 3.56 / 3.50 / 3.62 / 3.29 | 17.41 | 2.53 |

Replicate 1 is not the slowest in any rung. Each call carries about 2.55 s of fixed overhead that does not scale with B (index matrix, observed fits, `get_dfRes`, table build).

### Budget arithmetic (single core; no recommendation)

- **The three M fits as the document ran them:** FS 31.83 + DINA 1.73 + GRF 11.28 = **44.84 s** of fit time (rung walls 33.71 + 3.55 + 13.13 = 50.39 s, including R start-up and data prep).
- **FS-ij saving:** 15.85 s of fit time, bringing FS to 15.98 s and the three fits to 28.99 s.
- **B affordable per identifier**, bootstrap call alone (M fit excluded). Two estimates are given: (a) the task's total/5 rate, and (b) overhead-adjusted, `floor((T - overhead) / mean iteration)`.

| identifier | per-replicate (total/5) | (a) 30 s / 60 s / 120 s | mean iteration (s) | (b) 30 s / 60 s / 120 s |
|---|---|---|---|---|
| FS | 2.19 s | 13 / 27 / 54 | 1.673 | 16 / 34 / 70 |
| DINA | 0.85 s | 35 / 70 / 141 | 0.332 | 82 / 172 / 353 |
| GRF | 3.99 s | 7 / 15 / 30 | 3.482 | 7 / 16 / 33 |

These are 5-replicate means. The FS iteration times range from 1.28 to 2.29 s and DINA from 0.10 to 0.58 s, so the B figures carry that spread. The measurements were taken on pop-os under light load (`/proc/loadavg` read just after the ladder: 0.36 / 0.49 / 0.35 on 128 threads).

## Deviations (verbatim)

1. **Attempt 1 of the ladder crashed at M-FS before the fit ran:** `Error: object 'timings' not found` (exit 1, 1.39 s). The driver printed "STOP: M-FS (GATE A) did not pass; ladder halted."
   - **Cause:** my `data_prep.R` transplant omitted the document's `setup` chunk (lines 92–99: `timings <- list()`, `t_vignette_start <- proc.time()`). The `parallel-setup` chunk writes `timings$kfold <- NA_real_` and so depends on it.
   - **Not a gate result:** GATE A's checks (region, N, values) were never evaluated, and no gate condition was reinterpreted.
   - **Fix:** I added the `setup` chunk verbatim to `data_prep.R` (header updated to list it) and smoke-checked that the data prep evaluates.
   - **Rerun:** the ladder was rerun from a clean record. Attempt 1's log, `ladder.csv` and driver output are kept in `logs_attempt1/`.
   - Flagging for chat in case it prefers a crash at the gate rung to be treated as a stop.
2. **M-FS-ij direction.** The document passes no `ci_method`. At HEAD the FS call therefore resolves to `"field"` (at 0.2.0 it resolved to `"ij"`). "M-FS as the document ran it" was taken as the verbatim call, which resolves to `"field"`, so M-FS-ij sets `"ij"`.
3. **M-FS warning.** M-FS emitted one warning: `vi.grf.min` is deprecated and ignored (the transplanted call passes `vi.grf.min = -0.2`). It was left as is, since the call is verbatim.
4. **Two budget estimates.** The budget table gives both the task's total/5 rate and an overhead-adjusted figure. The second was added because the fixed ~2.55 s per call is a large share at B = 5.

## OPEN ITEMS (findings only)

1. **FS MR drift since 0.2.0.** The FS multiplier cells differ from the Table S2 payload at the third decimal (Ĥ 1.443947 vs 1.442118; upper 3.359577 vs 3.358359; Ĥᶜ upper 1.033039 vs 1.032018). They are equal at 2 dp. DINA and GRF MR, all naive cells, the selections and `n_family` are bit-identical. Something in the consistency-branch MR path changed between 0.2.0 and 0.4.0.9000. The cause was not investigated (out of scope for a read-only inventory).
2. **The FS call now runs the field block by default.** That takes the FS fit from ~16 s to ~32 s single-core here. The DINA and GRF calls stay on `"ij"`. The resulting asymmetry in the same document (FS carries `g$field`, DINA and GRF do not) is a property of the HEAD defaults.
3. **The task header's settings do not match the document for DINA and GRF.**
   - "DINA with a mean-effect floor of log 1.0" vs the document's `fs_dina_hr_threshold = 0.90`.
   - "the same candidate family as Table 1 ... `{er <= 0}` and `{pgr <= 0}` pre-specified" vs the DINA and GRF calls, which pass neither `conf_force` nor `conf.cont_jcuts`.
   - The DINA MR family is 84 here against 30 in the Table 1 DINA call.
4. **No package function renders Table S2's layout.** `mr_estimates_table()` shows naive + MR; `summarize_bootstrap_results()$table` shows naive + FB. A vignette reproducing Table S2 would assemble the three columns itself, as the document's `.engine_rows()` does.
5. **The `mr_estimates_table()` roxygen is stale.** It says the IJ interval comes "under the default `ci_method = "ij"`" (`R/mr_estimates_table.R:31`). The consistency default is `"field"`; the rendered interval is IJ either way.
6. **At B = 5 the bootstrap Ĥᶜ IJ interval is degenerate** (lower ≈ 0.01, upper 54–79 for all three identifiers), while Ĥ is finite. Any small-B vignette illustration will show this unless the interval is suppressed or B is raised.
7. **The FS fit in the 0.2.0 document took 7.4 s on 11 workers.** The 31.8 s here is single-core and includes the field block. The two are not comparable as a speed-up.

## Copied to `~/Downloads/postsel_stage0_2026-10-28/`

(listing appended below after the copy)

Two copies of the producing document are included: `_cf4d6432` (as rendered for Table S2) and `_HEAD` (differs at lines 144–145 only). `survival-post-selection.qmd` is `git show HEAD:` (identical to the working tree). The report's own size below is from before this listing was appended; the copy was refreshed afterwards.

```
/home/larryleon/Downloads/postsel_stage0_2026-10-28:
total 320
drwxrwxr-x 4 larryleon larryleon   4096 Oct 28 11:08 .
drwxr-xr-x 9 larryleon larryleon  20480 Oct 28 11:08 ..
-rw-rw-r-- 1 larryleon larryleon 115357 Oct 28 11:08 analysis_gbsg_survival_multimethod_cf4d6432.qmd
-rw-r--r-- 1 larryleon larryleon 115357 Oct 28 11:08 analysis_gbsg_survival_multimethod_HEAD.qmd
-rw-rw-r-- 1 larryleon larryleon    384 Oct 28 11:08 ladder.csv
drwxrwxr-x 2 larryleon larryleon   4096 Oct 28 11:08 manuscript_source
-rw-rw-r-- 1 larryleon larryleon   9673 Oct 28 11:08 package_facts_s2.md
-rw-rw-r-- 1 larryleon larryleon  22005 Oct 28 11:08 REPORT_postsel_stage0_2026-10-28.md
-rw-rw-r-- 1 larryleon larryleon  10100 Oct 28 11:08 survival-post-selection.qmd
drwxrwxr-x 2 larryleon larryleon   4096 Oct 28 11:08 table_s2_calls

/home/larryleon/Downloads/postsel_stage0_2026-10-28/manuscript_source:
total 32
drwxrwxr-x 2 larryleon larryleon 4096 Oct 28 11:08 .
drwxrwxr-x 4 larryleon larryleon 4096 Oct 28 11:08 ..
-rw-r--r-- 1 larryleon larryleon 1934 Oct 28 11:08 analysis_gbsg_survival_multimethod_payload.rds
-rw-rw-r-- 1 larryleon larryleon  945 Oct 28 11:08 manifest_excerpt.txt
-rw-r--r-- 1 larryleon larryleon 9263 Oct 28 11:08 _supp_gbsg_table.qmd
-rw-rw-r-- 1 larryleon larryleon  544 Oct 28 11:08 supplement_include_excerpt.txt

/home/larryleon/Downloads/postsel_stage0_2026-10-28/table_s2_calls:
total 68
drwxrwxr-x 2 larryleon larryleon  4096 Oct 28 11:08 .
drwxrwxr-x 4 larryleon larryleon  4096 Oct 28 11:08 ..
-rw-rw-r-- 1 larryleon larryleon  4551 Oct 28 11:08 bootstrap_calls.R
-rw-rw-r-- 1 larryleon larryleon 12741 Oct 28 11:08 data_prep.R
-rw-rw-r-- 1 larryleon larryleon  3374 Oct 28 11:08 dina_call.R
-rw-rw-r-- 1 larryleon larryleon  8328 Oct 28 11:08 fs_call.R
-rw-rw-r-- 1 larryleon larryleon  3589 Oct 28 11:08 grf_call.R
-rw-rw-r-- 1 larryleon larryleon 12756 Oct 28 11:08 table_assembly.R
```
