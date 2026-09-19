# REPORT — forestsearch 0.3.5.9000: the `R/` change review, and the fs-glms-interpretable revalidation questions Q1–Q3

**Date:** 2026-09-19 · **Task:** `dev/tasks/TASK_code_review_revalidation_2026-09-19.md` (committed as received, `86862f56`)
**Kind:** read-only. `R/` was not touched; no package was installed; nothing was re-run; no test file was executed;
no `R CMD check`; no data generator was called. Facts only — no recommendations, proposals or tasks.

**Pin.** `PIN = 64cc49813f097f410e7a0c44453c7a487cedcf78`, branch `feature/glm-extension`,
`2026-09-19 03:40:44 -0700`, *"report(binary): record the external mid-run push, and the tracked logs_or files"*.
Everything below is read at the pin (`git show $PIN:<path>`, `git grep $PIN`), never from the working tree,
because a binary stage-2 run may be live in another session in this tree.

---

## 1. Answers

### Q1 — the four MR default changes: was each Larry's decision, and where is it recorded?

| default | current value | verdict | source |
|---|---|---|---|
| `ci_method` `"ij"` → `"field"` | `"field"` (`R/fs_mr_inference.R:559`; `forestsearch()` fallback `R/forestsearch_main.R:3775`) | **RECORDED DECISION** | `dev/tasks/TASK_cimethod_note_2026-09-09.md:3` *"Approver: Larry (2026-09-09)"*; `:6` *"**Larry's constraint**, and the point of Gate D2c below: the IJ two-sided interval must not be lost"*; commit `beda000b` |
| `field_complement` → `TRUE` | `TRUE` (`R/fs_mr_inference.R:566`; fallback `:3790`) | **RECORDED DECISION** | `dev/tasks/TASK_cert20_2026-09-08.md:13` *"## Part D — Defaults become the recommendations (classified: changes behaviour; **decided by Larry F-2**)"*; `:3` *"Approver: Larry (F-1/F-2 decided 2026-09-08 …)"*; commit `fb62705c` |
| `field_scale_complement` → `"selected"` | `c("selected","none")` (`R/fs_mr_inference.R:568`; fallback `:3798`) | **RECORDED DECISION** | same two lines; plus `dev/notes/REVIEW_E1_fields_2026-09-08.md:99` (item **D-1**, recommendation (a), *"Larry's call."*) and `dev/notes/HANDOFF_mr_field_linux_2026-09-08.md:20` (*"decided by Larry 2026-09-08"*, and *"**Package default stays `"none"`** … a default flip is a separate decision"* — the flip taken hours later as F-2); commit `fb62705c` |
| `return_reselection` → `TRUE` | `TRUE` (`R/fs_mr_inference.R:561`; fallback `:3784`) | **RECORDED DECISION** | `TASK_cert20_2026-09-08.md:13`, `:3`, `:15`; commit `fb62705c` |

One fact the sites add to the question: **the DINA and GRF branches do not carry the new `ci_method` default.**
`.fs_apply_mr()` falls back to `"ij"` (`R/fs_mr_inference_methods.R:171`), a value that predates the pin's
range (introduced `b22c743c`, 2026-06-09; renamed `a909cece`, 2026-07-30) and was not flipped by `beda000b`.
`1d9401cb` states this deliberately: *"`ci_method` is unchanged and keeps the wrapper's own `"ij"` default"*
(`NEWS.md:210-211`). The other three defaults **are** inherited on those branches, because `.fs_apply_mr()`
reads them from `formals(fs_mr_inference)` at call time (`:176-187`).

### Q2 — was the HR/OR-only scope of threshold rule (A) Larry's decision?

**Yes, and the implementation matches the record exactly, with one date discrepancy.**
`dev/tasks/TASK_directive_A_2026-09-18.md:3` — *"**Authorized by:** Larry, 2026-09-18"*; `:14` — *"## Scope —
dispositions already taken (these govern)"*; `:16` — *"**HR and binary OR only.** RD, IRD, MD and IRR: no error,
no derivation, resolution unchanged"*; `:18-19` — *"The error and the derivation act only where the consistency
stage runs … Under `"dina"` and `"grf"` c2 is inert"*. The code agrees line for line
(`.fs_resolve_threshold_pair()`, `R/forestsearch_main.R:97-100` gate, `:135` derivation, `:140-143` stop).
**Discrepancy:** the chat record dates the decision 2026-09-17; every committed record dates it **2026-09-18**
(the task document's own header, and `dev/reports/STATUS_thresholds_workstream_2026-09-18.md:20`). No
2026-09-17 record of it exists in the repository.

### Q3 — did any committed pre-fix GRF run feeding the companion manuscript pass a factor covariate?

**No.** No run in the assessed set passed a factor-class covariate to the fixed evaluator before the fix
(`0cd33f7b`, 2026-09-16 21:34:01 -0700).
- The two **pre-fix** GRF campaigns (`grfmr`, and `idsweep`'s GRF rows, both on `gbsg_020`) search seven
  **integer** columns — `confounders_base <- c("er","age","meno","pgr","nodes","size","grade")`, identical at
  both producing commits and at the fix's parent — and nothing coerces them. `grfmr`'s exposure is already
  settled by a committed record: `quarto/simulations/gbsg_020/REPORT_grf_factor_exposure_2026-09-16.md:78`
  (**F1**), commit `35816257`.
- The **post-fix** GRF runs are `actg175` continuous `mdgrf` (its own catalog entry names *"the fix it waited
  for"*) and `actg175` binary `orgrf` (added 2026-09-19; its template coerces `bin_vars` to numeric first).
- The one GRF-identifier document in the manuscript set, `quarto/gbsg/analysis_gbsg_grf_mr.qmd`
  (fs-glms-interpretable), has an all-numeric covariate list and a **post-fix** payload (`8014cfa`,
  2026-09-17 17:40).
- The only manuscript document that passes factors (`quarto/actg175/analysis_actg175_continuous_mr.qmd:215`,
  six `as.factor` binaries) runs `use_grf = FALSE` and no GRF identifier at all, so it never reaches the
  evaluator.

### Part 1 summary lines — one per governing task document

- **(no task document; "Phase 1–4.5" arc, 2026-09-02 → 2026-09-03, `9bd8f96c..c1752b99`, 12 commits)** — built the
  subgroup-simulation wrapper family (`run_subgroup_sims()`, `summary`/`print`/`plot`/`compare_subgroup_sims()`,
  `subgroup_cox()`, `subgroup_glm()`, `benchmark_spec()`, `validate_subgroups()`, `forest_height()`), then
  lifted it to GLM outcomes, and added `df_source` / `baseline = "fixed"` to the GLM DGM pair. Class: **add-only**
  (ten new exports, five new S3 methods), with one **messages/display only** commit (`b54d82d7`).
- **`TASK_mr_vs_guohe_2026-09-04_addendum-A.md`** (`736088e3`) — `fs_mr_inference(return_reselection=)`, D6.
  Class: **add-only**.
- **`TASK_mr_field_vs_guohe_2026-09-05.md`** (`87880a24`) — `ci_method = "field"` added as a choice.
  Class: **add-only** (new construction on an opt-in path).
- **`TASK_mr_field_section7_2026-09-05.md`** (`a6702fd8`) — complement interval SE under `"field"`: Wald → IJ.
  Class: **changes the method** (on the then-opt-in field path).
- **`TASK_mr_field_uniform_2026-09-05.md`** (`b1316d22`, `3ad6ee13`) — `fs_mr_field_uniform()` κ(Σ̂) calibration
  and its `forestsearch()` pass-through; H2/H3 adjustment `M_cap` 12→40, `kappa` grid [1,2]→[1,3].
  Class: **add-only**, then **changes behaviour** (two knob defaults inside the new function).
- **`TASK_gbsg_frozen_intervals_2026-09-05.md`** (`c1b957b2`, `c1752b99` docs) — pass-throughs.
  Class: **add-only**.
- **`TASK_bias_coverage_display_2026-09-06.md`** (`72bd714a`) — `fs_sim_bias_coverage()` +
  `fs_plot_bias_coverage()`. Class: **add-only** (two new exports).
- **`TASK_mr_field_complement_2026-09-06.md`** (`ef3e609a`) — `field_complement=`, `side=`.
  Class: **add-only**.
- **`TASK_complement_refinements_2026-09-06.md`** (`8fbb3bdc`) — `ij_residual` variants and `field$joint`.
  Class: **add-only**.
- **`TASK_continuous_field_mac_2026-09-07.md` / `HANDOFF_…`** (`2f118042`, and `c0f48a7c` the Mac merge) —
  `scale = c("log","identity")` on `fs_sim_bias_coverage()`. Class: **add-only** (default byte-identical).
- **`PROPOSAL_complement_field_scale_2026-09-08_v2.md` + `TASK_field_studentize_stage1_e0_2026-09-08.md`**
  (`7245e898`) — `field_decompose`. Class: **add-only**.
- **`TASK_field_studentize_e1_2026-09-08.md`** (`3880022e`) — `field_scale_complement` and the `_s` companions.
  Class: **add-only** (default `"none"` at the time).
- **`TASK_cert20_2026-09-08.md` Part D** (`fb62705c`) — `field_complement`, `field_scale_complement`,
  `return_reselection` defaults become the recommendations. Class: **changes behaviour** (three defaults).
- **`TASK_cimethod_note_2026-09-09.md` Part D2** (`beda000b`) — `ci_method` default `"ij"` → `"field"`.
  Class: **changes behaviour** (a default).
- **`TASK_print_vignette_2026-09-09.md`** (`df251a96`, `105186cd`) — certified survival products reported from
  `print`/`summary.forestsearch()`; the low-p̂ caveat. Class: **messages/display only** (adds code; existing
  output byte-identical when no MR results are present).
- **`TASK_field_recovery_2026-09-09.md`** (`e819a58a`, `0f1eb18c`) — `field_recovery` diagnostics; `sens_H`
  reported. Class: **add-only**.
- **`TASK_o1_forwarding_2026-09-10.md`** (`1d9401cb`) — `.fs_apply_mr()` forwards the full 25-argument set.
  Class: **add-only** (no default changes on any branch; `ci_method` left at the wrapper's `"ij"`).
- **`TASK_grf_dina_fixes_2026-09-16.md`** (`0cd33f7b` P1, `064fce91` P2) — `.grf_code_column()` shared by
  `.build_grf_X()` and `.grf_evaluate_subgroup()`; `dina_subgroup(tau_sign=)` with `forestsearch()` passing
  `-1` for continuous/binary with `adverse_outcome = FALSE`. Class: **changes behaviour** (a selection path,
  on both).
- **`TASK_threshold_naming_docs_2026-09-18.md`** (`cc1714e0`, `0671525f`, `458cef48`) — `c1`/`c2`/`p*` named in
  every exported function; `threshold_config` documented; the frontier caps documented as a display trim.
  Class: **docs/roxygen only**.
- **`TASK_binary_default_or_2026-09-18.md`** (`5b1dac43`, `1ba22ded`, `fda9dcb0`, `7e4e4c87`) — binary
  `effect_measure` resolves to `"OR"`; the unreachable second site deleted. Class: **changes behaviour** (a
  default) + **dead-code removal**.
- **`TASK_threshold_sync_2026-09-18.md`** (`571b202c`, `85847aa0`) — SECTION 2B-ii: replicates and CV folds
  resolve what the parent resolved; the estimation entry points default binary to `"OR"`.
  Class: **changes behaviour** (a resolution path, on replay only; plus a default).
- **`TASK_directive_A_2026-09-18.md`** (`822cd69a`, `0037ef6d`, `70aab091`, `787d138e`) — `c2 > c1` is an error
  on the FS consistency path; a silent `c2` derives `0.80 * c1`; `"OR"` first in
  `make_effect_estimator()`'s binary choices. Class: **changes behaviour** ×2 (a new error, a new default
  derivation) + **add-only/alignment**.
- **`TASK_directive_C_2026-09-18.md`** (`c5bd3b5d`, `b8243432`, `31255426`, `c7e95160`, `c8b475b4`, `b9a32705`) —
  DINA refuses `RD`/`IRD`; inert frontier keys warn; `dina_frontier()` caps default `Inf` with a trim warning;
  the frontier print retitled; the `c2`/`p*` echoes annotated under dina/grf. Class: **changes behaviour** ×3
  (a new error, a new warning, two defaults) + **messages/display only** ×2.
- **`TASK_estimability_boundary_2026-09-18.md`** (`1719056f`, `f9b794f6`) — `NA`-with-reason where the estimand
  does not exist; `fs_dgm_feasibility()` added and exported. Class: **changes behaviour** (an estimate path) +
  **add-only** (one new export).

---

## 2. Part 0 — the pin

| | |
|---|---|
| `PIN` | `64cc49813f097f410e7a0c44453c7a487cedcf78` |
| branch | `feature/glm-extension` |
| `git log -1` | `2026-09-19 03:40:44 -0700` · *report(binary): record the external mid-run push, and the tracked logs_or files* |
| `DESCRIPTION` `Version` at `PIN` | **0.3.5.9000** |
| installed `Version` | 0.3.5.9000 |
| installed `Built` | `R 4.6.1; ; 2026-09-19 03:54:02 UTC; unix` |
| installed `Packaged` | absent from `packageDescription()` |
| library | `/home/larryleon/R/x86_64-pc-linux-gnu-library/4.6/forestsearch` |
| R | R version 4.6.1 (2026-06-24) |
| platform | `x86_64-pc-linux-gnu` |
| `pkgload` | installed, 1.5.3 |

**Installed namespace against `PIN`'s `R/`.** `PIN` has no `src/` directory and `pkgload` is installed, so the
comparison ran: `git archive $PIN DESCRIPTION NAMESPACE R` into a `mktemp` directory, every closure in
`asNamespace("forestsearch")` captured (formals plus `utils::removeSource()`d body) before and after
`pkgload::load_all(export_all = FALSE, quiet = TRUE)`, compared by name.

```
installed closures: 677
pin closures:       677
identical:          677
differing:            0
installed-only:       0
pin-only:             0
```

**The installed build is byte-identical in every closure to `PIN`'s `R/`.**

---

## 3. Part 1 — the `R/` changes

### 3.1 Anchors

| | SHA | date | subject | `DESCRIPTION` diff |
|---|---|---|---|---|
| **C0** | `f3975b99e5a625419b33a77f6504ff7ab793b751` | 2026-09-01 10:04:50 -0700 | *perf: unadjusted survival fits on subset vectors (0.3.5) — …* | `-Version: 0.3.4` / `+Version: 0.3.5` |
| **C1** | `66c4469b22d4c84be7ff8308650f29531d8d691c` | 2026-09-18 20:53:49 -0700 | *chore(DESCRIPTION): 0.3.5 -> 0.3.5.9000 development version* | `-Version: 0.3.5` / `+Version: 0.3.5.9000` |

**A structural fact that governs the "before/after C1" column and §A16:** `git log 66c4469b..$PIN -- R/` is
**empty**. Every one of the 55 `R/` commits in `C0..$PIN` landed *before* C1 — the last, `f9b794f6`, at
2026-09-18 20:02:11, 51 minutes before the version bump. So "the `R/` changes in 0.3.5.9000" cannot mean
"changes made after the version was bumped"; read as `C0..$PIN`, the range is the whole 0.3.5 → 0.3.5.9000
development line.

### 3.2 Inventory — 55 commits touching `R/` in `C0..$PIN`, oldest first

`before/after C1` is **before** for every row, so the column is omitted and stated once here.
`in §A` = named, directly or by its workstream, in the handoff inventory reproduced in §A of the task.

| # | SHA | date | `R/` files | functions touched | class | governing task document (basis) | §A |
|---|---|---|---|---|---|---|---|
| 1 | `9bd8f96c` | 09-02 | `run_subgroup_sims.R`, `summary_subgroup_sims.R` | `run_subgroup_sims`, `.run_subgroup_sims_one`, `subgroup_cox`, `validate_subgroups`, `benchmark_spec`, `print.subgroup_sims`, `summary.subgroup_sims`, `print.subgroup_sims_summary` | add-only | none; "Phase 1" (commit message only — no `dev/tasks/` document had been committed yet) | N |
| 2 | `2f301cdc` | 09-02 | `plot_subgroup_sims.R`, `run_subgroup_sims.R`, `summary_subgroup_sims.R` | `plot.subgroup_sims_summary`, `forest_height`, `.subgroup_sims_footnote`, `run_subgroup_sims`, `summary.subgroup_sims` | add-only | none; "Phase 2" | N |
| 3 | `f71f39f0` | 09-02 | `run_subgroup_sims.R` | `subgroup_cox`, `run_subgroup_sims`, `benchmark_spec` | add-only (new `lean` formal) | none; "Phase 2.1" | N |
| 4 | `b54d82d7` | 09-02 | `run_subgroup_sims.R` | `run_subgroup_sims` | messages/display only (scoped `options(future.rng.onMisuse="ignore")`, restored `on.exit`; no RNG state change) | none; "Phase 2.2" | N |
| 5 | `a8430e73` | 09-02 | `compare_subgroup_sims.R` | `compare_subgroup_sims` | add-only + moves existing code (replaces inline `.col_stats()`) | none; "Phase 3" | N |
| 6 | `14df006c` | 09-02 | `run_subgroup_sims.R` | `run_subgroup_sims`, `.run_subgroup_sims_one`, `subgroup_glm`, `subgroup_cox`, `validate_subgroups` | add-only | none; "Phase 4.1" | N |
| 7 | `310039f6` | 09-03 | `compare_subgroup_sims.R`, `summary_subgroup_sims.R` | `compare_subgroup_sims`, `summary.subgroup_sims`, `print.subgroup_sims_summary` | changes behaviour (effect-aware summaries, within the same-arc wrappers) | none; "Phase 4.2" | N |
| 8 | `f9b1e4cf` | 09-03 | `plot_subgroup_sims.R` | `plot.subgroup_sims_summary`, `.subgroup_sims_footnote*` | messages/display only | none; "Phase 4.3" | N |
| 9 | `e56afaf8` | 09-03 | `plot_subgroup_sims.R`, `run_subgroup_sims.R`, `summary_subgroup_sims.R` | `subgroup_glm`, `subgroup_cox`, `run_subgroup_sims`, `plot.subgroup_sims_summary` | add-only (binary gate-lift) | none; "Phase 4.4" | N |
| 10 | `9dd588b7` | 09-03 | `run_subgroup_sims.R` | `subgroup_glm`, `subgroup_cox`, `run_subgroup_sims` | add-only (count gate-lift) | none; "Phase 4.5" | N |
| 11 | `4066a953` | 09-03 | `generate_glm_dgm.R`, `run_subgroup_sims.R`, `simulate_from_glm_dgm.R` | `generate_glm_dgm`, `simulate_from_glm_dgm`, `run_subgroup_sims` | add-only (new `df_source` return element and `baseline="fixed"`; no RNG consumed, `df_super` byte-unchanged) | none; "Arc C' step C1" | N |
| 12 | `c1752b99` | 09-03 | `generate_glm_dgm.R` | — (roxygen) | docs/roxygen only | none | N |
| 13 | `736088e3` | 09-04 | `fs_mr_inference.R` | `fs_mr_inference`, `.fs_mr_se_from_ij` | add-only (`return_reselection`, D6) | `TASK_mr_vs_guohe_2026-09-04_addendum-A.md` (task doc committed same session, `276c3cee`) | **Y** |
| 14 | `87880a24` | 09-05 | `fs_mr_inference.R` | `fs_mr_inference`, `.fs_mr_se_from_ij` | add-only (`ci_method="field"` as a choice) | `TASK_mr_field_vs_guohe_2026-09-05.md` (commit message) | **Y** |
| 15 | `a6702fd8` | 09-05 | `fs_mr_inference.R` | `fs_mr_inference`, `.fs_mr_se_from_ij` | changes the method (complement SE Wald → IJ under `"field"`) | `TASK_mr_field_section7_2026-09-05.md` (commit message) | N |
| 16 | `b1316d22` | 09-05 | `forestsearch_main.R`, `fs_mr_field_uniform.R`, `fs_mr_inference.R` | `fs_mr_field_uniform`, `fs_mr_inference`, `forestsearch` | add-only | `TASK_mr_field_uniform_2026-09-05.md` (commit message) | N |
| 17 | `3ad6ee13` | 09-05 | `fs_mr_field_uniform.R` | `fs_mr_field_uniform` | changes behaviour (`M_cap` 12→40, `kappa_grid` [1,2]→[1,3]) | `TASK_mr_field_uniform_2026-09-05.md` (Gate 1 adjudication, commit message) | N |
| 18 | `c1b957b2` | 09-06 | `forestsearch_main.R`, `fs_mr_inference.R` | `forestsearch`, `.validate_outcome_threshold_config`, `fs_mr_inference` | docs/roxygen only (`df_source` in the return value) | `TASK_gbsg_frozen_intervals_2026-09-05.md` (commit message) | N |
| 19 | `72bd714a` | 09-06 | `fs_bias_coverage.R` | `fs_sim_bias_coverage`, `fs_plot_bias_coverage` | add-only (two new exports) | `TASK_bias_coverage_display_2026-09-06.md` (`64fb252f`, same session) | N |
| 20 | `ef3e609a` | 09-06 | `forestsearch_main.R`, `fs_bias_coverage.R`, `fs_mr_inference.R` | `fs_mr_inference`, `.fs_mr_field_complement`, `fs_sim_bias_coverage`, `forestsearch` | add-only (`field_complement`, `side`) | `TASK_mr_field_complement_2026-09-06.md` (`2c519337`, same session) | **Y** |
| 21 | `8fbb3bdc` | 09-06 | `forestsearch_main.R`, `fs_bias_coverage.R`, `fs_mr_inference.R` | `.fs_mr_ij_floor`, `.fs_mr_field_joint`, `.fs_mr_field_complement`, `fs_mr_inference`, `forestsearch` | add-only (`ij_residual`, `field$joint`) | `TASK_complement_refinements_2026-09-06.md` (`95e6c8c5`, same session) | N |
| 22 | `2f118042` | 09-07 | `fs_bias_coverage.R` | `fs_sim_bias_coverage` | add-only (`scale = c("log","identity")`, default byte-identical) | `TASK_continuous_field_mac_2026-09-07.md` / `HANDOFF_continuous_field_mac_2026-09-07.md` (`8ae83be9`) | N |
| 23 | `7245e898` | 09-08 | `forestsearch_main.R`, `fs_mr_inference.R` | `fs_mr_inference`, `.fs_mr_field_complement`, `forestsearch` | add-only (`field_decompose`, default `FALSE`) | `PROPOSAL_complement_field_scale_2026-09-08_v2.md` + `TASK_field_studentize_stage1_e0_2026-09-08.md` (`1bb8bdf1`, same session) | **Y** |
| 24 | `3880022e` | 09-08 | `forestsearch_main.R`, `fs_mr_inference.R` | `fs_mr_inference`, `.fs_mr_field_complement`, `forestsearch` | add-only (`field_scale_complement`, default `"none"` at the time) | `TASK_field_studentize_e1_2026-09-08.md` (`9a6e21cc`, same session) | **Y** |
| 25 | `c0f48a7c` | 09-08 | (merge) `fs_bias_coverage.R` | `fs_sim_bias_coverage` | merge commit — brings the Mac line's `scale="identity"`; no new `R/` content of its own | `TASK_continuous_field_mac_2026-09-07.md` / `TASK_actg175_continuous_intervals_2026-09-07.md` | N |
| 26 | `fb62705c` | 09-08 | `forestsearch_main.R`, `fs_mr_inference.R` | `fs_mr_inference`, `forestsearch` | **changes behaviour** — three defaults: `field_complement TRUE`, `field_scale_complement "selected"`, `return_reselection TRUE` | `TASK_cert20_2026-09-08.md` Part D (commit message + `c5044911`) | **Y** |
| 27 | `beda000b` | 09-09 | `forestsearch_main.R`, `fs_mr_inference.R` | `fs_mr_inference`, `forestsearch`, `.validate_outcome_threshold_config` | **changes behaviour** — `ci_method` default `"ij"` → `"field"` | `TASK_cimethod_note_2026-09-09.md` Part D2 (commit message + `ab9afb20`) | **Y** |
| 28 | `df251a96` | 09-09 | `forestsearch_methods.R` | `print.forestsearch`, `summary.forestsearch`, `.fs_print_mr_products`, `.fs_mr_products`, `.fs_mr_caveats`, `.fs_cat_caveat`, `.fs_sg_labels` | messages/display only (existing output byte-identical with no MR results) | `TASK_print_vignette_2026-09-09.md` (commit message) | N |
| 29 | `105186cd` | 09-09 | `forestsearch_methods.R` | `.fs_mr_caveats`, `.fs_mr_products`, `.fs_print_mr_products` | messages/display only (low-p̂ caveat) | `TASK_print_vignette_2026-09-09.md` (`8b9da196`) | N |
| 30 | `e819a58a` | 09-09 | `forestsearch_main.R`, `fs_mr_inference.R` | `.fs_mr_field_recovery`, `fs_mr_inference`, `forestsearch` | add-only (`field_recovery`, default `FALSE`) | `TASK_field_recovery_2026-09-09.md` (commit message + `b69fffab`) | N |
| 31 | `0f1eb18c` | 09-09 | `forestsearch_methods.R` | `print.forestsearch`, `.fs_mr_products`, `.fs_print_mr_products` | messages/display only (`sens_H` reported) | `TASK_field_recovery_2026-09-09.md` (commit message) | N |
| 32 | `1d9401cb` | 09-10 | `fs_mr_inference_methods.R` | `.fs_apply_mr`, `.fs_mr_reselection_from_focus` | add-only (ten dropped arguments forwarded; no default changes on any branch) | `TASK_o1_forwarding_2026-09-10.md` (`8fe49362`, same session) | **Y** (Q1) |
| 33 | `0cd33f7b` | 09-16 | `grf_helpers.R`, `grf_subg_harm_glm.R`, `grf_subgroup_labels.R` | `.grf_code_column` (new), `.build_grf_X`, `.grf_evaluate_subgroup`, `validate_grf_data` | **changes behaviour** (a selection path: factor cuts no longer evaluate to `NA` membership) | `TASK_grf_dina_fixes_2026-09-16.md` P1 (commit message) | **Y** (Q3) |
| 34 | `064fce91` | 09-16 | `dina_subgroup.R`, `forestsearch_helpers.R` | `dina_subgroup`, `.dina_tau_sign` (new), `.forestsearch_dina_select`, `.coerce_covariates_numeric` | **changes behaviour** (a selection path: `tau_sign` orients DINA's proposal floor) | `TASK_grf_dina_fixes_2026-09-16.md` P2 (commit message) | N |
| 35 | `cc1714e0` | 09-18 | `consistency_resample.R`, `forestsearch_main.R`, `fpr_approximation.R`, `fpr_calibration.R`, `fs_oc_predict.R`, `mrct_simulation.R`, `subgroup_consistency_helpers.R`, `subgroup_consistency_main.R`, `subgroup_search.R` | roxygen blocks of the exported entry points; `.validate_outcome_threshold_config`, `.consistency_glm_pieces`, `evaluate_subgroup_consistency`, `get_split_hr_fast`, `.make_eval_consistency_twostage`, `.make_candidate_rng_seeds` | docs/roxygen only | `TASK_threshold_naming_docs_2026-09-18.md` (`216f3405`, same session) | **Y** (A8) |
| 36 | `0671525f` | 09-18 | `forestsearch_main.R` | `.validate_outcome_threshold_config` (roxygen) | docs/roxygen only (`threshold_config` `@return`) | `TASK_threshold_naming_docs_2026-09-18.md` | **Y** (A8) |
| 37 | `458cef48` | 09-18 | `dina_subgroup.R`, `forestsearch_main.R` | `.dina_check_cap`, `.validate_outcome_threshold_config` (roxygen) | docs/roxygen only | `TASK_threshold_naming_docs_2026-09-18.md` | N |
| 38 | `5b1dac43` | 09-18 | `forestsearch_main.R` | `forestsearch` | **changes behaviour** (binary `effect_measure` default `"RD"` → `"OR"`) | `TASK_binary_default_or_2026-09-18.md` (`6d0ce4df`, same session) | **Y** (A3) |
| 39 | `1ba22ded` | 09-18 | `forestsearch_main.R` | `forestsearch` | dead-code removal (the unreachable second resolution site) | `TASK_binary_default_or_2026-09-18.md` | **Y** (A3) |
| 40 | `fda9dcb0` | 09-18 | `forestsearch_main.R` | `forestsearch` (roxygen) | docs/roxygen only | `TASK_binary_default_or_2026-09-18.md` | N |
| 41 | `7e4e4c87` | 09-18 | `subgroup_search.R` | `subgroup.search` (roxygen), `extract_idx_flagredundancy` | docs/roxygen only | `TASK_binary_default_or_2026-09-18.md` | N |
| 42 | `571b202c` | 09-18 | `forestsearch_main.R` | `forestsearch` (SECTION 2B-ii) | **changes behaviour** (a resolution path — on replay only; the parent fit is untouched) | `TASK_threshold_sync_2026-09-18.md` (`d0708011`, same session) | **Y** (A2) |
| 43 | `85847aa0` | 09-18 | `consistency_resample.R`, `glm_effect_estimators.R` | `make_effect_estimator`, `.consistency_glm_pieces`, `.consistency_adj_terms` | **changes behaviour** (binary default `"OR"` at the estimation entry points) | `TASK_threshold_sync_2026-09-18.md` (rider) | **Y** (A4) |
| 44 | `822cd69a` | 09-18 | `forestsearch_main.R` | `.fs_resolve_threshold_pair` (new), `forestsearch`, `.sync_args_call_all` | **changes behaviour** (a new error: `c2 > c1`, and disagreeing spellings) | `TASK_directive_A_2026-09-18.md` (`ed34ce13`, same session) | **Y** (A1) |
| 45 | `0037ef6d` | 09-18 | `forestsearch_main.R` | `.fs_resolve_threshold_pair`, `.sync_args_call_all` | **changes behaviour** (a new default derivation: `c2 = 0.80 * c1`) | `TASK_directive_A_2026-09-18.md` | **Y** (A1) |
| 46 | `70aab091` | 09-18 | `forestsearch_main.R` | `.validate_outcome_threshold_config` (roxygen) | docs/roxygen only | `TASK_directive_A_2026-09-18.md` | N |
| 47 | `787d138e` | 09-18 | `glm_effect_estimators.R` | `make_effect_estimator` | add-only / alignment (`"OR"` first in the binary `match.arg` choices; unreachable today) | `TASK_directive_A_2026-09-18.md` (rider) | **Y** (A4) |
| 48 | `c5bd3b5d` | 09-18 | `forestsearch_helpers.R`, `forestsearch_main.R` | `.dina_assert_ratio_estimand` (new), `.forestsearch_dina_select`, `.dina_tau_sign`, `forestsearch` | **changes behaviour** (a new error: DINA refuses `RD`/`IRD`) | `TASK_directive_C_2026-09-18.md` (`bef180ce`, same session) | **Y** (A5) |
| 49 | `b8243432` | 09-18 | `forestsearch_helpers.R`, `forestsearch_main.R` | `.forestsearch_dina_select`, `.resolve_dina_args`, `.map_dina_family`, `.validate_outcome_threshold_config` | **changes behaviour** (a new warning on inert frontier keys) | `TASK_directive_C_2026-09-18.md` | **Y** (A6) |
| 50 | `31255426` | 09-18 | `dina_subgroup.R`, `forestsearch_main.R` | `dina_frontier`, `.dina_check_cap`, `.dina_warn_cap_trim` (new), `forestsearch` | **changes behaviour** (two defaults `3L`/`10L` → `Inf`, plus a trim warning) | `TASK_directive_C_2026-09-18.md` | **Y** (A7) |
| 51 | `c7e95160` | 09-18 | `forestsearch_helpers.R` | `.forestsearch_dina_select` | messages/display only (frontier print retitled) | `TASK_directive_C_2026-09-18.md` | N |
| 52 | `c8b475b4` | 09-18 | `bootstrap_dofuture_main.R`, `forestsearch_cross_validation.R`, `forestsearch_helpers.R`, `forestsearch_main.R`, `forestsearch_methods.R`, `fs_family_report.R`, `interpret_search_config.R` | `.fs_c2_inert_note` (new), `forestsearch`, `forestsearch_bootstrap_dofuture`, `print_cv_params`, `summary.forestsearch`, `fs_family_report`, `interpret_search_config`, `.map_dina_family` | messages/display only (seven echo sites annotated) | `TASK_directive_C_2026-09-18.md` | **Y** (A8) |
| 53 | `b9a32705` | 09-18 | `subgroup_search.R` | `extract_idx_flagredundancy` (Rd link) | docs/roxygen only | `TASK_directive_C_2026-09-18.md` (remediation; see §3.5) | N |
| 54 | `1719056f` | 09-18 | `fs_family_report.R`, `glm_effect_estimators.R`, `subgroup_search.R` | `.fs_binary_cells`, `.fs_existence_reason`, `.fs_nonestimable` (all new), `make_effect_estimator`, `.make_cox_estimator`, `.make_glm_binary_estimator`, `.make_poisson_rate_estimator`, `fs_family_report`, `subgroup.search`, `search_combinations_parallel`, `evaluate_combination_with_status`, `fit_glm_for_subgroup` | **changes behaviour** (an estimate path: `NA` with a reason where the estimand does not exist) | `TASK_estimability_boundary_2026-09-18.md` (`cf97ec70`, same session) | **Y** (A9–A11) |
| 55 | `f9b794f6` | 09-18 | `fs_dgm_feasibility.R` (new file) | `fs_dgm_feasibility`, `print.fs_dgm_feasibility` | add-only (one new export) | `TASK_estimability_boundary_2026-09-18.md` | **Y** (A12) |

Scale: **31** distinct files under `R/`, **5,042** added-or-changed lines (`git diff C0 $PIN -- R/`).

### 3.3 `NEWS.md` — the development section at `PIN`, and its commit map

Header: **`NEWS.md:1`** — `# forestsearch (development version)`. The section runs to `:286`; `# forestsearch 0.3.5`
begins at `:288`. Nineteen entries, newest first as written:

| `NEWS.md` lines | entry (first phrase) | commit(s) that added the entry | `R/` commit(s) it describes |
|---|---|---|---|
| 3–27 | *New `fs_dgm_feasibility()`: a design-time check…* | `f9b794f6` | `f9b794f6` |
| 29–55 | *The estimator boundary now returns `NA` with a reason…* | `1719056f` | `1719056f` |
| 57–72 | *`dina_frontier()`'s display caps now default to `Inf`…* | `31255426` | `31255426`, `c7e95160` |
| 74–84 | *`dina_args` frontier keys now warn when they are ignored.* | `b8243432` | `b8243432` |
| 86–103 | *DINA now refuses the identity-scale estimands `"RD"` and `"IRD"`.* | `c5bd3b5d` | `c5bd3b5d` |
| 105–129 | *The two effect thresholds are now a pair…* | `7cabca90` (docs-only) | `822cd69a`, `0037ef6d`, `70aab091` |
| 131–152 | *Bootstrap replicates and cross-validation folds now resolve…* | `a389ebf8` (docs-only) | `571b202c` |
| 154–161 | *The estimation-layer entry points default binary to the odds ratio.* | `a389ebf8` (docs-only) | `85847aa0`, `787d138e` |
| 163–178 | *Binary outcomes now default to the odds ratio.* | `f391b816` (docs/tests-only) | `5b1dac43`, `1ba22ded`, `fda9dcb0` |
| 180–185 | *GRF membership on factor covariates.* | `323de083` (docs-only) | `0cd33f7b` |
| 187–193 | *DINA proposal floor orientation.* | `323de083` (docs-only) | `064fce91` |
| 195–211 | *The internal MR wrapper `.fs_apply_mr()` … now forwards the full argument set.* | `1d9401cb` | `1d9401cb` |
| 213–232 | *`fs_mr_inference()` gains `field_recovery`…* | `e819a58a` | `e819a58a`, `0f1eb18c` |
| 234–248 | *The `ci_method` default is now `"field"` (was `"ij"`)…* | `beda000b` | `beda000b` |
| 250–259 | *`fs_mr_inference()` defaults are now the recommended constructions…* | `fb62705c` | `fb62705c` |
| 261–265 | *`fs_mr_inference()` gains `field_decompose`…* | `7245e898` | `7245e898` |
| 266–272 | *`fs_mr_inference()` gains `field_scale_complement = c("none","selected")`…* | `3880022e` | `3880022e` |
| 274–286 | *Documentation: the threshold arguments.* | `cc1714e0` | `cc1714e0`, `0671525f`, `458cef48`, `fda9dcb0`, `7e4e4c87`, `70aab091` |

**`R/` commits in `C0..$PIN` with no `NEWS.md` entry — 24 of 55:**

- the twelve subgroup-simulation-wrapper commits `9bd8f96c`, `2f301cdc`, `f71f39f0`, `b54d82d7`, `a8430e73`,
  `14df006c`, `310039f6`, `f9b1e4cf`, `e56afaf8`, `9dd588b7`, `4066a953`, `c1752b99` — **ten new exported
  functions and five new S3 methods, none mentioned in `NEWS.md`** (`run_subgroup_sims`, `subgroup_cox`,
  `subgroup_glm`, `benchmark_spec`, `validate_subgroups`, `compare_subgroup_sims`, `forest_height`, plus
  `generate_glm_dgm`'s `df_source` and `simulate_from_glm_dgm(baseline="fixed")`);
- `736088e3` (`return_reselection` added — its later default flip is in `NEWS.md`, its introduction is not);
- `87880a24` (`ci_method="field"` added as a choice);
- `a6702fd8` (**changes the method**: the complement interval SE under `"field"` moves Wald → IJ);
- `b1316d22`, `3ad6ee13` (`fs_mr_field_uniform()` and its `M_cap`/`kappa_grid` changes);
- `c1b957b2` (roxygen);
- `72bd714a`, `2f118042`, `c0f48a7c` (`fs_sim_bias_coverage()` / `fs_plot_bias_coverage()` — **two new exports**,
  and the `scale` argument);
- `ef3e609a`, `8fbb3bdc` (`field_complement` and `ij_residual`/`field$joint` introduced);
- `df251a96`, `105186cd`, `0f1eb18c` (the print/summary reporting of certified products);
- `c8b475b4`, `c7e95160`, `b9a32705` (the display annotations and the Rd fix).

### 3.4 §A claims, checked against source at `PIN`

| claim | verdict | evidence |
|---|---|---|
| **A1** `.fs_resolve_threshold_pair()` at `forestsearch_main.R` ~line 90; `c2 > c1` a `stop()`; unset c2 derives 0.80·c1; 1.25→1.00, 1.00→0.80, 0.90→0.72 | **CONFIRMED**, line number exact | definition at `R/forestsearch_main.R:90`; `stop()` at `:140-143`; derivation `c2 <- 0.80 * c1` at `:135`; the three worked values in the roxygen at `:58`, and again in `NEWS.md:112-113` |
| **A2** SECTION 2B-ii: bootstrap and CV replicates resolve the parent fit's thresholds; parent-fit resolution unchanged | **CONFIRMED** | `R/forestsearch_main.R:2318` (section header) – `:2351`; the resolved naturals written into `effect.threshold`/`consistency.threshold` at `:2346-2348` and synced at `:2349-2351`; *"This runs AFTER resolution and writes only to `args_call_all`, so THIS fit's resolved thresholds are untouched; only a replay sees the difference"* (`:2337-2338`) |
| **A3** binary `effect_measure` default RD → OR at the live site; a second unreachable site deleted | **CONFIRMED** | live site `R/forestsearch_main.R:1614-1620` (`binary = "OR"` at `:1616`); the flip is `5b1dac43` (`-binary = "RD"` / `+binary = "OR"`); the second site deleted by `1ba22ded` (nine lines removed at the then-`:1831`) |
| **A4** `glm_effect_estimators.R` and `consistency_resample.R` binary default `"OR"`; `match.arg` order `c("OR","RD","RR","IRR","IRD")` | **CONFIRMED** | `R/glm_effect_estimators.R:112-120` (`binary = "OR"` at `:116`), `match.arg` choices at `:150-153`; `R/consistency_resample.R:245-249` (`binary = "OR"` at `:248`, inside `.consistency_glm_pieces()`, the resolution `consistency_resample()` hands it) |
| **A5** `.dina_assert_ratio_estimand(family, effect_measure)` in `forestsearch_helpers.R`, called at the `use_dina` derivation site; floor `if (family=="gaussian") hr.threshold else log(hr.threshold)` | **CONFIRMED**, with one addition | definition `R/forestsearch_helpers.R:1446-1462`; floor at `:1525-1526`. It is called at **two** sites, not one: `R/forestsearch_helpers.R:1524` (inside `.forestsearch_dina_select()`, the `subgroup_method="dina"` route) **and** `R/forestsearch_main.R:3050` (the `use_dina` + `selected_only` screening route). The roxygen says so at `:1434-1437` |
| **A6** `.DINA_FRONTIER_KEYS` and a frontier-key warning in `.forestsearch_dina_select()` | **CONFIRMED** | `R/forestsearch_helpers.R:1055-1056` (the seven keys); the single-per-fit `warning()` at `:1496-1507`, inside `.forestsearch_dina_select()` (`:1479`) |
| **A7** `dina_subgroup.R`: `max_per_covariate`/`max_subgroups` default `Inf`; `.dina_warn_cap_trim()` with class `dina_frontier_cap_trim` | **CONFIRMED** | `R/dina_subgroup.R:1195-1196` (both `Inf`, in `dina_frontier()`'s formals); `.dina_warn_cap_trim()` at `:1034`, condition class at `:1036`; raised at `:1306` and `:1327` |
| **A8** `.fs_c2_inert_note()` applied to both config banners; `threshold_config` `@return` documented | **CONFIRMED** | banners: GLM at `R/forestsearch_main.R:2227` used at `:2244-2246`; survival at `:2280` used at `:2295-2297`. `threshold_config` `@return` at `:1145-1166`. (Five further call sites exist: `interpret_search_config.R:181`, `bootstrap_dofuture_main.R:429`, `forestsearch_cross_validation.R:2062`, `forestsearch_methods.R:451`, plus `fs_family_report`) |
| **A9** existence condition — OR all four cells ≥ 1; RR/IRR/HR ≥ 1 event per arm; RD/IRD/MD untouched; non-estimable returns `estimate=NA`, `se=NA`, `converged=FALSE`, `reason`; counted in `filter_counts`; in `fs_family_report()`'s stage map; cannot rank | **CONFIRMED in full** | condition `R/glm_effect_estimators.R:261-275` (`.fs_existence_reason()`), cells `:251-258`, return shape `.fs_nonestimable()` `:237-248`; call sites HR `:295-299`, binary `:392-396`, rate `:742-751`; the written contract at `:217-234`. Counted: `R/subgroup_search.R:308-309` (`n_nonestimable`, `nonestimable_reasons`), incremented `:326-328`, printed under `details` `:228-235`. Stage row `R/fs_family_report.R:314-322`. Cannot rank: the candidate returns `status = 5L` (`R/subgroup_search.R:660-661`), the pre-existing fit-failure status |
| **A10** `converged` computed by every estimator and discarded at `subgroup_search.R:879` | **CONTRADICTED on the line reference; the substance holds** | at `PIN`, `R/subgroup_search.R:879` is `dimnames = list("Treat", names(hr.cox)))` — the one-row `conf.int` matrix coercion inside `fit_cox_for_subgroup()` (`:802`). The string `converged` appears **nowhere** in `R/subgroup_search.R`, at `PIN` or at the fix's parent `1719056f^`. What is true: `fit_glm_for_subgroup()` builds its return from `res$estimate` / `res$se` only (`:946-952`) and never reads `res$converged`; it reads `res$reason` instead (`:934-939`) |
| **A11** RD's tier-3 raw-proportions fallback returns `converged = FALSE` and is untouched | **CONFIRMED** | tier 3 returns `converged = FALSE` at `R/glm_effect_estimators.R:552`; the non-convergence catch is scoped to `OR`/`RR` at `:433-437` with the reason stated in the comment `:429-432`. `.estimate_rd()` is unchanged over `C0..$PIN` (`git log -S '.estimate_rd'` in range: no hits) |
| **A12** `R/fs_dgm_feasibility.R`, exported, commit `f9b794f6`; `fs_dgm_feasibility(dgm, n, n.min, d0.min, d1.min, n_rep, tolerance)`; draws through the DGM's own generator; does not change the RNG kind; GLM path only | **CONFIRMED**, formals list incomplete | file added by `f9b794f6`; `export(fs_dgm_feasibility)` at `NAMESPACE:162`, `S3method(print,fs_dgm_feasibility)` at `:25`. Signature at `R/fs_dgm_feasibility.R:104` has **ten** formals — the seven named plus `effect_measure`, `seed`, `rand_ratio`. Own generator: `simulate_from_glm_dgm()` at `:171-172`. RNG: `set.seed(s_i)` only, *"kind unchanged: seed only"* (`:170`), the caller's `.Random.seed` saved at `:151-152` and restored `on.exit` at `:153-155`. GLM only: `stop()` unless `inherits(dgm, "glm_dgm")` at `:116-121`, naming `setup_gbsg_dgm()` as unsupported |
| **A13** `n.min = 60` strict; `d0.min`/`d1.min` = 10/10 non-strict, skipped for continuous and count — unchanged over `C0..$PIN` | **CONFIRMED** | defaults `R/forestsearch_main.R:1526` (`n.min = 60`), `:1546-1547` (`d0.min = 10`, `d1.min = 10`). Strict: `if (nx <= n.min) return(status 4)` at `R/subgroup_search.R:646` (GLM) and `:696` (survival) — the size must *exceed* 60. Non-strict: `d0_sg < d0.min` at `:631`; `meets_event_criteria()` uses `>=` at `:788`. Skipped: `is_continuous` covers `c("continuous","count")` at `:621-622`, Status 3 bypassed at `:634` / `:642`. Unchanged: `git diff C0 $PIN` on those two files shows **no** changed line mentioning `n.min`, `d0.min` or `d1.min` except one added comment (`subgroup_search.R:701`) |
| **A14** the oracle helper keeps the pooled 5/5 criterion and adds the four-cell condition — where does it live? | **CONFIRMED; it lives in a template, not in `R/`** | `.logit_or_ci()` appears nowhere under `R/` at `PIN`. The in-scope copy is `quarto/simulations/actg175/binary_020/sim_fs_mr_field_or_template.qmd:552-580`: *"LEGACY guard, unchanged: at least 5 events and 5 non-events POOLED over arms"* at `:558` with the test at `:559`; *"FOUR-CELL condition (TASK_binary_study_redesign_2026-09-18 Step 3)"* at `:560-563` with the four tests at `:565-568`, citing `.fs_existence_reason(), R/glm_effect_estimators.R:261-276`. A second copy is in `quarto/simulations/actg175/binary_020/maxeffCons_mr_coverage_sweep_or075.qmd`; the template asserts the copies match at `:480-491` |
| **A15** the five test files exist | **CONFIRMED** (read, not run) | `tests/testthat/helper-threshold-sync.R` (479 lines), `test-threshold-sync.R` (184), `test-threshold-pair-directive-a.R` (432), `test-directive-c.R` (748), `test-binary-default-or-entry-points.R` (109) |
| **A16** *"Two workstreams (thresholds; admission floors) produced every `R/` change in 0.3.5.9000"* | **CONTRADICTED** | The §3.2 table names **at least eleven** distinct governing workstreams across `C0..$PIN`, and 32 of the 55 `R/` commits predate 2026-09-16 — the subgroup-simulation-wrapper arc (12), the MR-field / complement / studentization line (13 commits across 9 task documents, including the four Q1 default flips), the print/vignette and field-recovery reporting line (4), MR forwarding (1), and the GRF/DINA identifier fixes (2). The thresholds workstream (`cc1714e0`…`c8b475b4`, `b9a32705`) and the admission-floor / estimability workstream (`1719056f`, `f9b794f6`) account for **21** of the 55. Under the narrower reading — commits *after* C1 — the claim is vacuous: there are none |

### 3.5 Check coverage, from the committed records (nothing was run)

| | tree | record | `R/` commits after it, up to `$PIN` |
|---|---|---|---|
| last `--as-cran` run | **`b9a32705`** (2026-09-18 16:12:43) | `dev/reports/CHECK_ascran_2026-09-18.md:154-179` — the Larry-authorized verification run, `rcmdcheck(args="--as-cran")`, PDF manual and vignettes included, 10.2 min, **`0 errors \| 1 warning \| 2 notes`**, reproducing the reference set exactly (`:161-168`); `:208-209` *"Workstream record closed."* | **2**: `1719056f`, `f9b794f6` |
| last full-suite run | **`b9a32705`** (the `checking tests` phase of that same run) | `CHECK_ascran_2026-09-18.md:178` — `checking tests ... OK`, **`FAIL 0 \| WARN 21 \| SKIP 76 \| PASS 5063`** | **2**: `1719056f`, `f9b794f6` |
| last *standalone* `devtools::test()` | tree **not recorded**; the report was committed at `9de4d0c5` (2026-09-18 13:02:31) | `quarto/simulations/actg175/binary_020/REPORT_binary_default_or_2026-09-18.md:120-122` — 7.12 min, 361 files, **`FAIL 0 \| ERROR 0 \| PASS 5160 \| SKIP 3 \| WARN 32`**; quoted again at `dev/reports/STATUS_thresholds_workstream_2026-09-18.md:108-109` | **14**: `571b202c`, `85847aa0`, `822cd69a`, `0037ef6d`, `70aab091`, `787d138e`, `c5bd3b5d`, `b8243432`, `31255426`, `c7e95160`, `c8b475b4`, `b9a32705`, `1719056f`, `f9b794f6` |

An earlier `--as-cran` run at `4447a6c6` (2026-09-18 14:32:26) is recorded at `CHECK_ascran_2026-09-18.md:13-32`
with `1 error | 2 warnings | 2 notes`; `b9a32705` is the remediation commit that closed both new findings.

### 3.6 Static CRAN/style scan of the 5,042 added-or-changed `R/` lines in `C0..$PIN`

Report only; nothing was fixed.

**(a) RNG state.** Four `set.seed()` sites added; one restores the caller's stream, three do not.

| site | restores caller's RNG state? |
|---|---|
| `R/fs_dgm_feasibility.R:170` | **yes** — `.Random.seed` saved `:151-152`, restored `on.exit` `:153-155`; comment *"kind unchanged: seed only"* |
| `R/fs_mr_field_uniform.R:93` (`fs_mr_field_uniform()`, only when `seed` is supplied) | **no** — no `.Random.seed` save/restore and no `on.exit()` anywhere in the file |
| `R/fs_mr_inference.R:864` (`fs_mr_inference()`, the `ci_method=="field"` branch, only when `seed` is supplied) | **no** — same |
| `R/run_subgroup_sims.R:512` (benchmark-subgroup regeneration, only when `benchmarks` is supplied) | **no** — the file restores the `future` plan (`:818`, `:830`) and the `options()` (`:881`) on exit, but not the RNG stream |

No `RNGkind()` was added. One assignment to `.Random.seed` was added —
`R/fs_dgm_feasibility.R:154`, and it is the restoration itself.

**(b) `options()` / `par()` / `Sys.setenv()` / `setwd()` without `on.exit()`.** One `options()` pair added, and it
is restored: `R/run_subgroup_sims.R:880-881`. No `par()`, `Sys.setenv()` or `setwd()` was added.

**(c) `print()` / `cat()` outside print/format/summary methods and not behind a verbosity flag.** 48 sites added;
**every one is inside a print method or behind a verbosity flag.**

| file | enclosing function | sites | disposition |
|---|---|---|---|
| `forestsearch_methods.R` | `.fs_print_mr_products`, `.fs_cat_caveat` | 21 | helpers called only from `print.forestsearch()` (`:375`) and `summary.forestsearch()` (`:591`) |
| `fs_dgm_feasibility.R` | `print.fs_dgm_feasibility` | 13 | print method |
| `run_subgroup_sims.R` | `print.subgroup_sims` | 6 | print method |
| `summary_subgroup_sims.R` | `print.subgroup_sims_summary` | 2 | print method (one is the `print()` hit, `:332`) |
| `bootstrap_dofuture_main.R` | `forestsearch_bootstrap_dofuture` | 2 | inside `if (details)` (`:406`) |
| `forestsearch_cross_validation.R` | `print_cv_params` | 2 | called only under `if (details)` (`:451`, `:944`) |
| `subgroup_search.R` | `subgroup.search` | 2 | inside `if (details)` (`:194`) |

**(d) `<<-`:** none added. **Bare `T` / `F`:** none added. **`library()` / `require()` inside functions:** none
added. **`forestsearch:::` self-calls:** none added.

**(e) Non-ASCII** (`tools::showNonASCIIfile()` on each of the 31 changed files at `PIN`). 30 hit lines across
7 files; **3 of them were added in `C0..$PIN`**, all in one file and all added by `72bd714a`:

- `R/fs_bias_coverage.R:238`, `:239`, `:251` — plot-caption strings carrying `Φ`, `·`, `−`, `∈`.

This is the file the certification WARNING names: *"checking code files for non-ASCII characters —
`R/fs_bias_coverage.R`"* (`dev/reports/CHECK_ascran_2026-09-18.md:40`), carried over unchanged from the
reference set. The other 27 hits are pre-existing: `R/fpr_calibration.R` (16), `R/forestsearch_helpers.R` (3),
`R/subgroup_search.R` (3), `R/bootstrap_dofuture_main.R` (2), `R/generate_glm_dgm.R` (2),
`R/mrct_simulation.R` (1).

**(f) Functions newly exported in the range** — ten, from the `NAMESPACE` diff, plus five S3 methods
(`plot.subgroup_sims_summary`, `print.fs_dgm_feasibility`, `print.subgroup_sims`,
`print.subgroup_sims_summary`, `summary.subgroup_sims`). Every one has `@export` in its roxygen, consistent
with `NAMESPACE`, and **every formal of every one carries an `@param`**.

| function | file:line | formals | `@param` missing | `@return` | `@examples` |
|---|---|---|---|---|---|
| `benchmark_spec` | `run_subgroup_sims.R:49` | 3 | none | yes | **unwrapped** |
| `compare_subgroup_sims` | `compare_subgroup_sims.R:87` | 6 | none | yes | `\dontrun{}` |
| `forest_height` | `plot_subgroup_sims.R:50` | 4 | none | yes | `\dontrun{}` |
| `fs_dgm_feasibility` | `fs_dgm_feasibility.R:104` | 10 | none | yes | `\donttest{}` |
| `fs_plot_bias_coverage` | `fs_bias_coverage.R:208` | 5 | none | yes | **none** |
| `fs_sim_bias_coverage` | `fs_bias_coverage.R:73` | 7 | none | yes | **none** |
| `run_subgroup_sims` | `run_subgroup_sims.R:665` | 20 | none | yes | `\dontrun{}` |
| `subgroup_cox` | `run_subgroup_sims.R:109` | 2 | none | yes | **unwrapped** |
| `subgroup_glm` | `run_subgroup_sims.R:252` | 10 | none | yes | `\dontrun{}` |
| `validate_subgroups` | `run_subgroup_sims.R:432` | 4 | none | yes | **none** |

**(g) One further observation.** Five `# [delivery sentinel: …]` comment lines, added in this range, are still
present at `PIN` on line 2 of `R/generate_glm_dgm.R`, `R/plot_subgroup_sims.R`, `R/run_subgroup_sims.R`,
`R/simulate_from_glm_dgm.R` and `R/summary_subgroup_sims.R`. Recorded, not fixed.

---

## 4. Q1 evidence — the four MR defaults

### 4.1 Every site in `R/` at `PIN` that supplies a value when the caller does not

| # | default | function | path:line | value | identifier paths that reach it |
|---|---|---|---|---|---|
| 1 | `ci_method` | `fs_mr_inference()` formal (`match.arg` first element) | `R/fs_mr_inference.R:559` | `c("field","ij","wald")` → **`"field"`** | any direct call to `fs_mr_inference()` that omits it |
| 2 | `ci_method` | `forestsearch()` pass-through fallback (`.g_mr`) | `R/forestsearch_main.R:3775` | **`"field"`** | `subgroup_method = "consistency"` |
| 3 | `ci_method` | `.fs_apply_mr()` internal fallback (`.g`) | `R/fs_mr_inference_methods.R:171` | **`"ij"`** | `subgroup_method = "dina"` and `"grf"` |
| 4 | `ci_method` | `fs_fdr_report()` hard-coded `mr_inference_args` | `R/fs_fdr_report.R:241` | `"ij"` | `fs_fdr_report()` only |
| 5 | `field_complement` | `fs_mr_inference()` formal | `R/fs_mr_inference.R:566` | **`TRUE`** | direct calls |
| 6 | `field_complement` | `forestsearch()` fallback | `R/forestsearch_main.R:3790` | **`TRUE`** | consistency |
| 7 | `field_complement` | `.fs_apply_mr()` | `R/fs_mr_inference_methods.R:182-183` | `.d("field_complement")` → **`TRUE`** (read from `formals(fs_mr_inference)` at call time, `:154-155`) | dina, grf |
| 8 | `field_scale_complement` | `fs_mr_inference()` formal (`match.arg`) | `R/fs_mr_inference.R:568` | `c("selected","none")` → **`"selected"`** | direct calls |
| 9 | `field_scale_complement` | `forestsearch()` fallback | `R/forestsearch_main.R:3798` | **`"selected"`** | consistency |
| 10 | `field_scale_complement` | `.fs_apply_mr()` | `R/fs_mr_inference_methods.R:186-187` | `.d(...)` → **`"selected"`** | dina, grf |
| 11 | `return_reselection` | `fs_mr_inference()` formal | `R/fs_mr_inference.R:561` | **`TRUE`** | direct calls |
| 12 | `return_reselection` | `forestsearch()` fallback | `R/forestsearch_main.R:3784` | **`TRUE`** | consistency |
| 13 | `return_reselection` | `.fs_apply_mr()` | `R/fs_mr_inference_methods.R:176-177` | `.d(...)` → **`TRUE`** | dina, grf |

A fourteenth mention, `R/forestsearch_methods.R:122` (`ci_method = g$ci_method`), reads the resolved value for
display and supplies nothing. The chat record's lead is correct: **site 3 is the `"ij"` fallback for DINA/GRF.**

### 4.2 The commit that introduced each current value

| site(s) | commit | date | subject | value before |
|---|---|---|---|---|
| 1, 2 (`ci_method` → `"field"`) | `beda000b` | 2026-09-09 13:09:43 -0700 | *Part D2: the `ci_method` default becomes "field", plus the certified survival-products NOTE (TASK_cimethod_note_2026-09-09). Classification: changes behaviour. Gate D2 PASS on all three limbs.* | formal `c("ij","wald","field")` → `"ij"`; fallback `.g_mr(..., "ij")` |
| 3 (`.fs_apply_mr()` → `"ij"`) | `b22c743c` | 2026-06-09 23:23:59 -0700 | *feat(debias-gate): extend Tier-2 gate to dina/grf via per-method candidate families* (renamed `debias_gate_args` → `mr_inference_args` at `a909cece`, 2026-07-30) | — introduced there; **outside `C0..$PIN`, and never flipped** |
| 5, 6 (`field_complement` → `TRUE`) | `fb62705c` | 2026-09-08 20:50:48 -0700 | *Part D of TASK_cert20_2026-09-08 (classified changes-behaviour, decided by Larry F-2) …* | `FALSE` in both |
| 8, 9 (`field_scale_complement` → `"selected"`) | `fb62705c` | same | same | formal `c("none","selected")`; fallback `"none"` |
| 11, 12 (`return_reselection` → `TRUE`) | `fb62705c` | same | same | `FALSE` in both |
| 7, 10, 13 (`.fs_apply_mr()` inherits) | `1d9401cb` | 2026-09-10 23:09:07 -0700 | *O-1: `.fs_apply_mr()` forwards the full MR argument set (DINA + GRF operational)* | the three were **not forwarded at all** before; the wrapper dropped 10 of 25 arguments |

**Other MR defaults changed by the same commits.** `fb62705c` changed exactly the three above, and no others —
it states *"`field_decompose` stays `FALSE` and `ci_method` is untouched (open decision …)"*. `beda000b`
changed exactly `ci_method`. `1d9401cb` changed **no** default: *"Add-only: no default changes on any branch"*
(`NEWS.md:206-207`), each forwarded value read from `formals(fs_mr_inference)` rather than restated.

### 4.3 Authorization, quoted

**`field_complement` / `field_scale_complement` / `return_reselection` — `TASK_cert20_2026-09-08.md`**
(committed as received at `c5044911`, 2026-09-08 20:39:51):

> `:3` — *"Date: 2026-09-08. Author: chat (spec). Executor: Claude Code (Linux, unattended overnight; Larry
> offline until the morning). **Approver: Larry (F-1/F-2 decided 2026-09-08; kickoff paste = compute go).**
> Reviewer: the Linux MR-field chat."*

> `:13` — *"## Part D — Defaults become the recommendations (classified: changes behaviour; **decided by Larry
> F-2**)"*

> `:15` — *"D1. `R/fs_mr_inference.R`: `fs_mr_inference()` defaults `field_complement = TRUE`,
> `field_scale_complement = c("selected", "none")` (so `match.arg` yields `"selected"`),
> `return_reselection = TRUE`; `field_decompose` stays `FALSE`; **`ci_method` default untouched** (open
> decision, Larry: which two-sided interval is primary in user-facing output …)"*

Its report repeats the attribution: `quarto/simulations/gbsg_020/REPORT_defaults_flip_2026-09-08.md:3` —
*"Spec: `dev/tasks/TASK_cert20_2026-09-08.md`, Part D (classified *changes behaviour*; **decided by Larry,
F-2**)."*

The decision trail behind F-2, for `field_scale_complement`, is `dev/notes/REVIEW_E1_fields_2026-09-08.md:99`
(committed `7361ae68`, 2026-09-08 17:47:57):

> *"**D-1 Adoption of field-s.** Options: (a) promote field-s to the documented one-sided complement product …
> **Recommendation: (a), with the shortfall stated** … **Larry's call.**"*

and `dev/notes/HANDOFF_mr_field_linux_2026-09-08.md:20`, which records the adoption *and* that the package
default had not yet moved:

> *"**Complement Ĥᶜ, upper bound on β(Ĥᶜ): field-s** (`field_scale_complement = "selected"`, R1 studentized
> complement field; `dev/notes/NOTE_complement_product_2026-09-08.md`, **decided by Larry 2026-09-08**) … Campaign
> convention from here: `FS_S7_FIELD_SCALEC=selected`. **Package default stays `"none"`** so committed bundles
> stay byte-reproducible; **a default flip is a separate decision.**"*

That separate decision is F-2, three hours later the same evening.

**`ci_method` — `TASK_cimethod_note_2026-09-09.md`** (committed as received at `ab9afb20`,
2026-09-09 12:59:10):

> `:3` — *"Date: 2026-09-09. Author: chat (spec). Executor: Claude Code (Linux). **Approver: Larry
> (2026-09-09).** Reviewer: the Linux MR-field chat."*

> `:6` — *"**Why now.** Part D flipped `field_complement`, `field_scale_complement` and `return_reselection` to
> the recommendations but left `ci_method` at `"ij"`. … C-4 closed the open question. **Larry's constraint, and
> the point of Gate D2c below: the IJ two-sided interval must not be lost** — the IJ SE is computed
> unconditionally, before and independent of the field gate, so `"field"` adds the field block without removing
> anything."*

> `:32` — *"**D2c — IJ IS NOT LOST (the criterion Larry named).** …"*

Its report confirms the criterion was met: `quarto/simulations/gbsg_020/REPORT_cimethod_flip_2026-09-09.md:181`
— *"The side-by-side table above is the concrete form of **Larry's criterion** … **The default gains the
certified one-sided products and loses nothing.**"* The same report also records, at `:197` and `:231`, that
*"Larry put tests in scope for this fix"* on 2026-09-09 (`f4664aed`), including the stale assertion at
`test-mr-inference.R:85` that pinned the old default.

**The chat record's leads, each verified with `git log --all`:**

| lead | verdict |
|---|---|
| `dev/tasks/PROPOSAL_complement_field_scale_2026-09-08_v2.md` | **present**, added `1bb8bdf1` (2026-09-08 09:47:34). It is explicitly *not* a decision record: `:3` — *"**Status: proposal for Larry's decision — not a task document; nothing goes to CC and no code changes anywhere until approval.**"*; `:4` — *"`ci_method = "ij"` remains the reported two-sided default throughout"*; `:34` — the new argument's default is *`"none"` (byte-identical)*; `:49-62` list P-1…P-6 as *"Decisions for Larry (defaults in brackets)"* |
| `dev/notes/REVIEW_E1_fields_2026-09-08.md` | **present**, added `7361ae68` (2026-09-08 17:47:57). Contains D-1 (quoted above); `:106` — *"Constructions in the package: `field` (documented), `field-s` (add-beside, **defaults off**, validated) … Interim documented rule unchanged until D-1 is decided"* |
| `dev/tasks/claude_cc_task_guohe_supplement_2026-09-09.md` | **present**, added `179ae409` (2026-09-09 12:40:02). `:25` — *"defaults for `field_complement`, `field_scale_complement`, `return_reselection` … **If the defaults are not `TRUE` / `"selected"` / `TRUE`, that contradicts the recorded adoption — STOP and report.**"* — it treats the adoption as already on record, consistent with F-2 the previous evening. `:66` notes *"the open two-sided `ci_method` default"*, i.e. written before `beda000b` |
| `REVIEW_certification_2026-09-09.md` | **absent from every reachable commit** — `git log --all --diff-filter=A -- '*REVIEW_certification_2026-09-09.md'` is empty in forestsearch **and** in `~/Documents/GitHub/fs-glms-interpretable`. Confirms the 2026-09-13 report |
| `SUMMARY_survival_properties_2026-09-10.md` | **absent from every reachable commit**, in both repositories. Confirms the 2026-09-13 report |

A sweep of `dev/tasks/`, `dev/notes/` and `dev/reports/` at `PIN` for the four formal names returns 30 further
files; none records a decision on any of the four beyond the documents quoted above. `claude/` does not exist
as a directory (only `.claude/` and `CLAUDE.md`); neither names the four. In
`~/Documents/GitHub/fs-glms-interpretable` (read-only; nothing written there), the four names appear in
`dev/verification/REPORT_revalidation_2026-09-19.md`, `dev/tasks/TASK_gbsg_mr_2026-09-17.md`, the
`dev/writeup/identifiers/` reports and the five `quarto/` analysis documents — as settings passed explicitly,
not as decision records. `REPORT_revalidation_2026-09-19.md:65-68` states the consequence for the manuscript:
*"7–9 are default changes only; all five documents pass these arguments explicitly, so none is reachable."*

### 4.4 `NEWS.md`

| default | `NEWS.md` lines at `PIN` | commit that added the entry |
|---|---|---|
| `ci_method` | `:234-248` — *"The `ci_method` default is now `"field"` (was `"ij"`) …"* | `beda000b` |
| `field_complement`, `field_scale_complement`, `return_reselection` | `:250-259` — *"`fs_mr_inference()` defaults are now the recommended constructions … `field_complement = TRUE`, `field_scale_complement = "selected"`, `return_reselection = TRUE`, mirrored in the `forestsearch()` pass-through fallbacks."* | `fb62705c` |

The `.fs_apply_mr()` `"ij"` fallback is documented at `:210-211` — *"`ci_method` is unchanged and keeps the
wrapper's own `"ij"` default"* (added by `1d9401cb`).

### 4.5 Verdicts

| default | verdict |
|---|---|
| `ci_method` `"ij"` → `"field"` | **RECORDED DECISION** — `TASK_cimethod_note_2026-09-09.md:3` (*"Approver: Larry (2026-09-09)"*) and `:6` / `:32` (the constraint and criterion attributed to Larry) |
| `field_complement` → `TRUE` | **RECORDED DECISION** — `TASK_cert20_2026-09-08.md:13` (*"decided by Larry F-2"*), `:3`, `:15` |
| `field_scale_complement` → `"selected"` | **RECORDED DECISION** — same, plus `REVIEW_E1_fields_2026-09-08.md:99` (D-1) and `HANDOFF_mr_field_linux_2026-09-08.md:20` |
| `return_reselection` → `TRUE` | **RECORDED DECISION** — same as `field_complement` |

Facts only; Larry confirms each.

---

## 5. Q2 evidence — the scope of threshold rule (A)

### 5.1 The record

**`dev/tasks/TASK_directive_A_2026-09-18.md`** (committed as received at `ed34ce13`, 2026-09-18 13:43:43):

> `:3` — *"**Opened:** 2026-09-18 · **Repository:** forestsearch · **Authorized by: Larry, 2026-09-18.**"*

> `:14` — *"## Scope — dispositions already taken (these govern)"*

> `:16-17` — *"1. **HR and binary OR only.** RD, IRD, MD and IRR: no error, no derivation, resolution unchanged,
> byte-identical to baseline in every probe cell."*

> `:18-19` — *"2. **The error and the derivation act only where the consistency stage runs** —
> `subgroup_method = "consistency"`. Under `"dina"` and `"grf"` c2 is inert: no error, no derivation, exactly as
> today."*

> `:20-21` — *"3. **c2 > c1 is a `stop()`**, message naming values and the spellings the caller used … c2 = c1
> passes."*

> `:22-25` — *"4. **c1 supplied with c2 not supplied derives c2 = 0.80 × c1 on the ratio scale** — an additive
> shift on the log scale, never `0.8 × log(c1)`. Worked values: 1.25 → 1.00, 1.00 → 0.80, 0.90 → 0.72."*

> `:26` — *"5. **Explicit beats derived.** A supplied c2 (either spelling) is never overridden."*

> `:28` — *"7. **No method change.** Selection logic untouched."*

> `:112-116` — *"## Out of scope … RD / IRD / MD / IRR behaviour. Floors. fs-glms-interpretable.
> `fpr_calibration()`'s own c2 ≤ c1 check (already present and already documented)."*

**`dev/reports/STATUS_thresholds_workstream_2026-09-18.md`:**

> `:20` — *"`c2 > c1` is an error on the forest-search consistency path (`822cd69a`); a silent `c2` derives
> `0.80 * c1` on the ratio scale (`0037ef6d`). **Scope-gated to `subgroup_method = "consistency"` with survival
> `HR` or binary `"OR"`.** The default pair (1.25, 1.0) is the rule's fixed point. 86 acceptance checks; 36 of
> 175 probe cells move, all in scope."*

> `:35` — *"**`c2 > c1` stops, and a silent `c2` derives `0.80 * c1`** — **HR and binary OR, consistency path
> only**; disagreeing spellings of one threshold also stop; derivation on the ratio scale, announced once,
> carried into every replicate by the sync."*

**`dev/reports/REPORT_directive_A_2026-09-18.md`:**

> `:96` — *"**6a** … **PASS** — 36 cells move, all in scope and on the consistency path: 18 derive (6 of them at
> the fixed point, where only the announcement is new) and 18 error. **The other 139, including all 28 dina/grf
> cells and every RD / IRD / MD / IRR cell, are byte-identical.**"*

> `:141-142` — *"## Out of scope, untouched … RD / IRD / MD / IRR behaviour. Floors. `fs-glms-interpretable`.
> `fpr_calibration()`'s own c2 ≤ c1 check."*

**Date of the decision.** Every committed statement dates it **2026-09-18**: the task header's *"Authorized by:
Larry, 2026-09-18"* (`:3`), and its framing of the scope as *"dispositions already taken"* (`:14`) within the
same day's workstream, whose opening audit (`TASK_audit_criterion_and_defaults_2026-09-18.md:3`) is likewise
*"Opened: 2026-09-18"*. A repository-wide search for `c2 > c1` / `0.80 * c1` in `dev/` finds no 2026-09-17
record. **The chat record's date of 2026-09-17 is not corroborated in-repo.**

### 5.2 The implementation

`.fs_resolve_threshold_pair()`, `R/forestsearch_main.R:90-158`. The gate, verbatim:

```r
 96  # The pair rule governs the consistency path on a ratio estimand only.
 97  ratio_pair <- identical(as.character(subgroup_method)[1L], "consistency") &&
 98    (identical(outcome_type, "survival") ||
 99     (identical(outcome_type, "binary") && identical(effect_measure, "OR")))
100  if (!ratio_pair) return(c2)
```

The derivation and the stop:

```r
134  derived <- isTRUE(user_set_threshold) && !isTRUE(user_set_consistency)
135  if (derived) c2 <- 0.80 * c1
...
140  if (isTRUE(c2 > c1)) {
141    stop(sprintf("c2 > c1 not allowed for FS: %s = %.2f exceeds %s = %.2f",
142                 c2_name, c2, c1_name, c1), call. = FALSE)
143  }
```

The disagreeing-spellings errors are at `:113-124`; the announcement at `:149-155`; `derived` is `FALSE` inside a
replicate because SECTION 2B-ii supplies the resolved `consistency.threshold` explicitly (`:145-148`).

**Single call site:** `R/forestsearch_main.R:1652-1658`, inside SECTION 1A3 (`:1639-1641`), placed *"before the
Section 1B capture so the resolved value is captured with the rest, and re-synced at SECTION 2B-ii"* (`:1650-1651`).

### 5.3 Per-estimand × identifier, from the code

Each cell gives **validation** (does the `c2 > c1` stop fire?) / **derivation** (does an unset `c2` become
`0.80 · c1`?). The governing lines are the same for every cell: `:97-100` decides, `:135` derives, `:140-143` stops.

| estimand | FS consistency | dina | grf |
|---|---|---|---|
| **HR** (survival) | **yes / yes** — `:98` `identical(outcome_type,"survival")` makes `ratio_pair` TRUE; stop `:140-143`, derivation `:135` | no / no — `:97` requires `"consistency"`; `:100` returns `c2` unchanged | no / no — `:97`, `:100` |
| **OR** (binary) | **yes / yes** — `:99` `binary && effect_measure=="OR"`; `:135`, `:140-143` | no / no — `:97`, `:100` | no / no — `:97`, `:100` |
| **RR** (binary) | no / no — `:99` requires `effect_measure=="OR"`, so `ratio_pair` is FALSE; `:100` | no / no | no / no |
| **RD** (binary) | no / no — `:99`; `:100` | no / no | no / no |
| **IRR** (count) | no / no — `outcome_type=="count"` is in neither `:98` nor `:99`; `:100` | no / no | no / no |
| **IRD** (count) | no / no — `:98-99`; `:100` | no / no | no / no |
| **MD** (continuous) | no / no — `:98-99`; `:100` | no / no | no / no |

### 5.4 The tests that pin the out-of-scope behaviour (read, not run)

`tests/testthat/test-threshold-pair-directive-a.R`:

- `:18` — the scope comment: *"(HR) or binary "OR". dina, grf, RD, IRD, MD and IRR are inert."*
- `:60-83` — `test_that("nothing outside the consistency path on a ratio estimand moves", …)`:
  `:70` `expect_identical(sum(out), 112L)` (175 cells, 63 in scope), `:71`
  `expect_identical(b[out, ], p[out, ])`;
  `:73-77` the named dispositions — `for (est in c("binary-RD","continuous-MD","count-IRD","count-IRR"))`
  with `expect_identical(b[k, ], p[k, ], info = est)`;
  `:78-82` — `for (m in c("dina","grf"))` with `expect_identical(sum(k), 14L)` and
  `expect_identical(b[k, ], p[k, ], info = m)`.
- `:279-304` — `test_that("dina, grf and the non-OR estimands are inert", …)`: a degenerate pair
  `c1 = 0.90, c2 = 1.00` returns `1.00` untouched under `"dina"` and `"grf"` (`:283-290`), under
  `effect_measure` `"RD"`, `"RR"`, `"IRD"` (`:291-296`), and under `continuous`/`"MD"` and `count`/`"IRR"`
  (`:297-303`) — *"No error, no derivation, c2 returned untouched"* (`:282`).
- `:408-431` — the rider: `make_effect_estimator()`'s binary choice order.

### 5.5 One line

**Record and implementation agree.** The only divergence is in the chat record's *date* for the decision
(2026-09-17); every committed record dates it 2026-09-18, and no 2026-09-17 record of the rule exists in the
repository.

---

## 6. Q3 evidence — the GRF factor-covariate fix

### 6.1 The fix

| | |
|---|---|
| SHA | **`0cd33f7b99475995af26c23864dde1215e60ab0a`** |
| date | **2026-09-16 21:34:01 -0700** (the chat record says "landed 2026-09-17"; the commit is 2026-09-16) |
| subject | *GRF membership coding (TASK_grf_dina_fixes_2026-09-16 P1): `.grf_evaluate_subgroup()` codes each covariate with `.grf_code_column()`, the per-column coding moved out of `.build_grf_X()` so the forest matrix and the evaluator share one definition; a cut on a 0/1 factor covariate no longer evaluates to NA membership (on the MD design sim_id 1: 645 of 1,257 candidates were dropped, now 0; admitted 234 -> 474; selection unchanged); guard fits F1-F3, F5, F6 identical to the installed build; unit tests* |
| governing task | `dev/tasks/TASK_grf_dina_fixes_2026-09-16.md`, Part P1 (committed as received at `0a697abb`, 2026-09-16 21:19:06) |
| records | `quarto/simulations/actg175/continuous/REPORT_grf_dina_fixes_2026-09-16.md` (added `baaf0bc6`, 2026-09-16 21:50:21) and `quarto/simulations/actg175/continuous/REPORT_md_grf_2026-09-16.md` (added `94d9cf53`, 2026-09-17 00:15:45). **Both verified present.** The chat record's attribution to the md-field-rerun workstream is consistent: the defect was found in `REPORT_md_grf_stage1_2026-09-16.md` §1.6(c), and `TASK_grf_dina_fixes` was issued to fix it |

**The before/after hunk**, `R/grf_subgroup_labels.R`, inside `.grf_evaluate_subgroup()` (`:208`):

```diff
@@ -219,7 +219,10 @@
       v  <- cj$variable[r]; op <- cj$op[r]; val <- cj$value[r]
       if (!v %in% names(df))
         stop("Subgroup variable '", v, "' not found in data.", call. = FALSE)
-      x <- df[[v]]
+      # Code the column exactly as the forest's covariate matrix does
+      # (.grf_code_column(), shared with .build_grf_X()): a factor with levels
+      # "0"/"1" is compared as 0/1, not as a factor, which R cannot order.
+      x <- .grf_code_column(df[[v]])
       member <- switch(op,
                        "<=" = x <= val,
                        ">"  = x >  val,
```

The shared helper is new at `R/grf_helpers.R:719-731` (`.grf_code_column()`), lifted verbatim out of
`.build_grf_X()` (`R/grf_subg_harm_glm.R:882-895`, which now calls it at `:890`).

### 6.2 Reach

Callers of `.grf_evaluate_subgroup()` at the fix's parent, `0cd33f7b^` (identical set at `PIN`, with shifted
line numbers):

| caller | path:line at `0cd33f7b^` | reached by |
|---|---|---|
| `.forestsearch_grf_select()` — effect re-selection over the candidate family | `R/forestsearch_helpers.R:1612` | `subgroup_method = "grf"` |
| `.forestsearch_grf_select()` — the selected subgroup's membership on `df` / `df.predict` / `df.test` | `R/forestsearch_helpers.R:1850`, `:1853`, `:1856` | `subgroup_method = "grf"` |
| `.fs_mr_grf_members()` — MR's family membership | `R/fs_mr_inference_methods.R:29` | `subgroup_method = "grf"` with `mr_inference = TRUE` |
| `grf.subg.harm.glm()` — the frontier path's `treat.recommend` | `R/grf_subg_harm_glm.R:548` | `subgroup_method = "grf"`, and a standalone `grf.subg.harm.glm(grf_selection = "frontier")` call |
| `grf.subg.harm()` | `R/grf_main.R:308` | the survival GRF search |
| `.betaHhat_truth_*()` | `R/betaHhat_truth.R:88` | truth computation on GRF definitions |

**DINA does not reach the fixed code.** No DINA function appears in that list, at the fix's parent or at `PIN`.
The committed source read is `REPORT_grf_factor_exposure_2026-09-16.md:70` (**§5**): *"DINA never compares a
data-frame column: candidate membership is computed on a numeric matrix … **DINA's evaluator does not have the
defect:** factors are coerced (all-numeric levels) or rejected with an error (other levels) before any
comparison; there is no path on which a factor is compared with `<=` / `>=`."* Consistent with the record's
identification of `dinamr` as DINA and `grfmr` as GRF.

The `use_grf = TRUE` **screening** path also does not reach it: it consumes `grf_res$tree.cuts` as strings for
the consistency engine, whose `evaluate_comparison()` does its own factor coercion
(`REPORT_grf_factor_exposure_2026-09-16.md:56`, quoting `R/forestsearch_main.R:2548-2601` and
`R/forestsearch_helpers.R:200-205`).

### 6.3 The run set

Established from the catalogs (`quarto/simulations/gbsg_020/current_status.md`,
`quarto/simulations/actg175/continuous/current_status.md`), the ACTG175 binary campaign task, and the
companion repository's own scope statement.

| # | run | where | GRF? | named as manuscript input by |
|---|---|---|---|---|
| 1 | **`grfmr`** — 18 cells × 2,000 replicates, MR on | `quarto/simulations/gbsg_020/results/grf_effMaxSG_*_grfmr_*.rds` (54 bundles) | yes, `subgroup_method="grf"`, frontier | `gbsg_020/current_status.md:48-53` (*"### 2.3 `grfmr` — GRF, complete"*), `:157` |
| 1a | `grfmrsmk` — the 5-replicate Stage 1 smoke | `…/results/grf_effMaxSG_…_grfmrsmk_res_1_5.rds` | yes | same catalog entry |
| 2 | **`idsweep`** GRF rows — part of 288 cell-runs × 500 replicates, **MR off** | `…/results/grf_eff*_nomr_idsweep_res_1_500.rds` | yes | `gbsg_020/current_status.md:93-121`, `:159` (*"dina, consistency, grf"*) |
| 3 | **`mdgrf`** — 4 cells × 2,000 replicates, MR on | `quarto/simulations/actg175/continuous/mr_md_harm/grf_effMaxSG_mr_field_*_mdgrf_d5000/` (12 bundles) | yes | `actg175/continuous/current_status.md:24`, `:44`, `:67-74` |
| 4 | **`orgrf`** — 1 cell × 1,000 replicates, MR on | `quarto/simulations/actg175/binary_020/mr_or_harm/grf_effMaxSG_mr_field_or150_n500_nb20_orgrf_d5000/` (2 bundles) | yes | `dev/tasks/TASK_actg175_binary_campaign_2026-09-17.md:1`, `:7`, `:105-109` |
| 5 | `quarto/gbsg/analysis_gbsg_grf_mr.qmd` + payload (fs-glms-interpretable) | companion repo | yes, `subgroup_method="grf"`, frontier | `fs-glms-interpretable/dev/verification/REPORT_revalidation_2026-09-19.md:81` (document 3 of 6) |

**Not assessed, and why.** `pBoc` (5 GRF bundles, 30 replicates, `…_nomr_pBoc_res_1_30.rds`, added `96e743db`
2026-09-12) is a Part B OC smoke that no brief or catalog names as manuscript input. GRF results under
`quarto/applications/`: no brief or catalog names any `quarto/applications/` payload as manuscript input —
`REPORT_revalidation_2026-09-19.md:89-90` puts *"any multimethod vignette in the forestsearch repo"* explicitly
out of scope. Their exposure is nevertheless already settled by a committed record (§6.5, F3).

### 6.4 Vintage

`pkg_version` does not discriminate: the fix did not bump `DESCRIPTION`, so both pre- and post-fix bundles record
`forestsearch_version` **0.3.5**. Basis is therefore the run's own committed record where one exists, otherwise
ancestry of the commit that added its bundles relative to `0cd33f7b`.

| run | vintage | basis |
|---|---|---|
| `grfmr` (+ `grfmrsmk`) | **pre-fix** | ancestry — bundles added `6765c372` (2026-09-12 11:04:52) and `4e21ea2c` (2026-09-12 20:38:54); `0cd33f7b` is not an ancestor of either |
| `idsweep` GRF rows | **pre-fix** | ancestry — bundles added `d62e1391` (2026-09-13 12:07:59) |
| `mdgrf` | **post-fix** | the run's committed record — `actg175/continuous/current_status.md:24` names *"**the fix it waited for** `REPORT_grf_dina_fixes_2026-09-16.md`"*; corroborated by ancestry: `git merge-base --is-ancestor 0cd33f7b c4c49572` → **yes** (first bundle commit `c4c49572`, 2026-09-16 22:36:55, 62 minutes after the fix; last `432f69dc`, 2026-09-17 00:11:32) |
| `orgrf` | **post-fix** | ancestry — bundles added `5ef6b748` (2026-09-19 03:15:25) |
| companion `analysis_gbsg_grf_mr` payload | **post-fix** | ancestry in the companion repo — `8014cfa` (2026-09-17 17:40:48) |

### 6.5 The committed exposure record comes first

`quarto/simulations/gbsg_020/REPORT_grf_factor_exposure_2026-09-16.md`, commit **`35816257`**
(2026-09-16 17:53:36 -0700), from `dev/tasks/TASK_grf_factor_exposure_2026-09-16.md` (`e3e6f652`) — a read-only
check made **before** the fix, on the finding that became it.

> `:78` — *"**F1.** The survival `grfmr` results are **not exposed** to the §1.6(c) defect: no factor covariate
> reaches the evaluator, and the direct count on sim_id 1 is 0 NA memberships of 779 candidates, with the
> committed selection reproduced exactly."*

> `:80` — *"**F3.** The ACTG175 binary and GBSG applied documents coerce or carry numeric covariates and are not
> exposed; the three factor-carrying applied documents use GRF as a screening cut generator or the standalone
> tree, not the frontier evaluator; **no applied document runs `subgroup_method = "grf"` on factor covariates.**"*

> `:81` — *"**F4.** DINA coerces (all-numeric levels) or rejects (other levels) factor covariates before
> evaluation."*

**Which cells it covers, and on what basis.** Two bases, both quoted:

1. *Classes as passed*, `:27-38` (**§2.3**) — `str(df[confs])` on the regenerated `grfmr` replicate:
   *"Every covariate the identifier sees is `integer`; **no factor or character column among `confs`**. The
   simulated frame does carry factor columns (`v1`–`v7`, the DGM's internal dichotomizations), but they are not
   in `confounders_base` and are not passed. No chunk converts columns before the `forestsearch()` call."*
   This is a property of the **template and its `confounders_base`**, so it covers every cell of every campaign
   rendered from `sim_fs_maxeffCons_fb_mr_field_m1_template.qmd` — all 18 `grfmr` cells, the `grfmrsmk` smoke,
   and (see §6.6) `idsweep`'s GRF rows.
2. *A direct count*, `:58-64` (**§4**) — one `grfmr` cell, **A124_h150_n500 sim_id 1**, regenerated with the
   campaign's seed scheme and RNG kind: *"candidates enumerated 779; with NA membership on the data (direct
   `.grf_evaluate_subgroup()` over every candidate) **0 of 779**; `sel_effect` finite on 779 of 779; admitted
   115 … **selected `{age > 45} & {meno <= 0}`, n = 100, naive Cox HR 1.643501**"*, matching the committed
   bundle row exactly (`:64`: *"**They agree** on definition, n, naive effect, admitted count and family
   size"*). This covers that one cell/replicate directly.

Not re-derived here.

### 6.6 The two instruments, for the pre-fix runs the record does not name

`grfmr` is named by F1, so it is not re-derived. `idsweep`'s GRF rows are **not** named by F1, so both
instruments were applied to them.

**Instrument 1 — classes as passed.** `scripts_dinamr/idsweep.sh:51`, `:61-73` render the **same template** the
exposure record names as the one `render.sh` renders:
`quarto/simulations/gbsg_020/sim_fs_maxeffCons_fb_mr_field_m1_template.qmd`
(cf. `REPORT_grf_factor_exposure_2026-09-16.md:23`). At the producing commit `d62e1391`, at the `pBoc` commit
`96e743db`, and at the fix's parent `0cd33f7b^`, that template's two governing lines are **character-identical**:

```
:704   confounders_base <- c("er", "age", "meno", "pgr", "nodes", "size", "grade")
:1013  confs <- intersect(confounders_base, names(df))
```

No coercion step sits between the draw and the `forestsearch()` call. No committed DGM or super-population
object exists under `gbsg_020` to corroborate against (`git ls-tree` finds no `*dgm*.rds` / `*super*.rds`
there), so the corroboration step is not available — recorded as an open item. No generator was called.

**Instrument 2 — stored payloads.** One bundle,
`quarto/simulations/gbsg_020/results/grf_effMaxSG_fb_mr_field_m1_h150_knoise0_n500_nb20_nomr_idsweep_res_1_500.rds`,
read at `PIN` into a `mktemp` directory (`meta`: `subgroup_method grf`, `campaign idsweep`,
`forestsearch_version 0.3.5`, `hostname Mac-Studio-3.local`, **`mr_inference FALSE`**), 500 rows × 168 columns,
498 detected (`detected == 1L`; `detected` is an integer flag):

| statistic | value |
|---|---|
| `n_family`, `p_hat_H`, `nv_H_est` — NA count | **500 of 500 (all rows), 498 of 498 detected rows** |
| `admitted_n` on detected rows — NA count; min / median / max | **0 NA**; 3 / 110 / 600 |
| `n_sel` on detected rows — NA count; min / median / max | **0 NA**; 60 / 96 / 350 |
| covariates appearing in stored selected (`sg_def`) and candidate (`p_hat_top_labels`) labels | **`age`, `er`, `grade`, `meno`, `nodes`, `pgr`, `size` — all seven** |

Reading: the membership and re-selection columns are NA **because `idsweep` ran with MR off** — `n_family`,
`p_hat_H` and the naive/field columns are not populated on that campaign at all, so they carry no signal about
membership. The columns that *are* populated carry the signal instead: `admitted_n` is finite on every detected
row with a median of 110, and all **seven** covariates — including `grade`, the 1/2/3 integer column — appear in
the stored selected labels. A factor-covariate drop produces the opposite signature (a covariate absent from
every label, and a collapsed admitted count), which is what the pre-fix MD design showed.

**Positive control.** The pre-fix `mdgrf` Gate-1 smoke **artifacts are not committed**: the only `mdgrf` bundles
in the tree are the four post-fix cells (first added `c4c49572`, all with `0cd33f7b` as an ancestor), and
`scripts_mdgrf/smoke_identity.R` is a script, not a payload. The pre-fix numbers survive only in the committed
*report*, `quarto/simulations/actg175/continuous/REPORT_md_grf_stage1_2026-09-16.md`, which is the positive
control in narrative form:

> `:96` — *"On this design the analysis frame carries `hemo`, `homo`, `drugs`, `race`, `gender`, `symptom`,
> `str2` as **factors** (`sapply(df[confounders_analysis], class)` on sim_id 1: the six continuous covariates
> integer/numeric, the seven binary ones `factor`) … Every candidate whose `v1` or `v2` is one of them is
> evaluated to NA membership."*

> `:97` — *"1,257 candidates enumerated; **612 use continuous covariates only and all 612 are scorable
> (`sel_effect` finite); 645 use a binary covariate and none is scorable** … Direct check: `{gender <= 0}` on
> the factor column → **NA on 500 of 500 rows** (with the warning); on the same column coerced to numeric → 74
> members. In the 20 smoke replicates **no selected rule uses a binary covariate**; MR's `n_family` (604–680)
> coincides with the scorable continuous-only count (612 on sim_id 1), so MR's family excludes them too."*

That is the signature — a whole class of covariates absent from every selected rule, and `n_family` collapsed to
the continuous-only count. The `idsweep` GRF bundle shows neither.

### 6.7 Answer

**No run in the assessed set passed a factor-class covariate to the GRF membership evaluator before the fix.**
The classes, as evidence:

| run | vintage | covariate classes as passed | verdict |
|---|---|---|---|
| `grfmr` (18 cells) + `grfmrsmk` | pre-fix | seven `integer` columns (`er`, `age`, `meno`, `pgr`, `nodes`, `size`, `grade`), template `:704`/`:1013`, no coercion | **not exposed** — settled by `REPORT_grf_factor_exposure_2026-09-16.md:78` (F1) and `:27-38`, `:58-64` |
| `idsweep` GRF rows | pre-fix | the same seven `integer` columns, same template, character-identical at the producing commit | **not exposed** — instrument 1 (§6.6) plus a payload reading with no drop signature |
| `mdgrf` (4 cells) | **post-fix** | the MD design's seven binary covariates *are* factors — but every bundle was produced after `0cd33f7b`, which is precisely the campaign that waited for the fix | **not affected** — the evaluator codes them |
| `orgrf` (1 cell) | **post-fix** | `sim_fs_mr_field_or_template.qmd:719-725` coerces `bin_vars` to numeric before the call (*"DINA (and GRF) require NUMERIC covariates"*) | **not affected** |
| companion `analysis_gbsg_grf_mr.qmd` | **post-fix** | `survival::gbsg` with `grade3 <- ifelse(grade == "3", 1, 0)`; `confounders <- c("age","meno","size","grade3","nodes","pgr","er")` (`:185`, `:188`) — all numeric | **not affected** |
| companion `analysis_actg175_continuous_mr.qmd` | n/a | six binaries **are** `as.factor` (`:215`) — but the fit passes `use_grf = FALSE` (`:301`) and no `subgroup_method = "grf"`, so GRF never runs | **cannot reach the evaluator** |

No output is affected, and nothing was re-run.

---

## 7. OPEN ITEMS

1. **§A10's line reference does not hold at `PIN`.** `R/subgroup_search.R:879` is the Cox `conf.int` matrix
   coercion inside `fit_cox_for_subgroup()`; the string `converged` appears nowhere in that file, at `PIN` or at
   `1719056f^`. The substance — the search never reads the estimator's `converged` — holds at
   `R/subgroup_search.R:946-952`. Recorded, not corrected anywhere.
2. **§A12's formals list is incomplete.** `fs_dgm_feasibility()` has ten formals; the claim names seven. The
   three unnamed are `effect_measure`, `seed`, `rand_ratio`.
3. **§A5's call-site count is incomplete.** `.dina_assert_ratio_estimand()` is called at two sites, not one.
4. **§A16 does not hold** under either reading (see §3.4).
5. **Q2's decision date.** The chat record says 2026-09-17; every committed record says 2026-09-18. Searched:
   `dev/tasks/TASK_directive_A_2026-09-18.md`, `dev/tasks/TASK_audit_criterion_and_defaults_2026-09-18.md`,
   `dev/reports/STATUS_thresholds_workstream_2026-09-18.md`, `dev/reports/REPORT_directive_A_2026-09-18.md`, and
   `git grep -E '0\.80 \* c1|c2 > c1' $PIN -- dev/`. No 2026-09-17 record found.
6. **Q3's fix date.** The chat record says the fix landed 2026-09-17; the commit is 2026-09-16 21:34:01 -0700.
   Its second record, `REPORT_md_grf_2026-09-16.md`, was committed 2026-09-17 00:15:45 — which may be the source
   of the 2026-09-17 reading.
7. **`REVIEW_certification_2026-09-09.md` and `SUMMARY_survival_properties_2026-09-10.md` are absent** from every
   reachable commit in forestsearch **and** in `~/Documents/GitHub/fs-glms-interpretable`. Searched with
   `git log --all --diff-filter=A -- '*<name>'` in both. Confirms the 2026-09-13 report.
8. **No committed DGM or super-population object exists under `quarto/simulations/gbsg_020/`**, so instrument 1's
   corroboration step (`vapply(<frame>[<covariates>], class, "")` on a stored object) could not be run for
   `idsweep`. Searched `git ls-tree -r --name-only $PIN | grep -iE 'gbsg_020/.*(dgm|super).*\.rds'` — empty.
   The trace from the template is character-identical at the producing commit, which is the instrument the task
   specifies first.
9. **The pre-fix `mdgrf` Gate-1 smoke artifacts are not committed**, so the positive control could only be quoted
   from `REPORT_md_grf_stage1_2026-09-16.md` rather than recomputed. Searched
   `git ls-tree -r --name-only $PIN | grep -iE 'actg175/continuous.*(mdgrf|grf).*(smk|smoke|stage1)'` — only the
   two reports and `scripts_mdgrf/smoke_identity.R`.
10. **`idsweep` ran with MR off**, so its bundles carry no membership or re-selection columns to count `NA` in
    (`n_family`, `p_hat_H` and the naive/field columns are `NA` on all 500 rows by construction). Instrument 2
    was read on the columns that are populated (`admitted_n`, `n_sel`, the stored labels); the requested NA
    counts on the membership and re-selection columns do not exist for that campaign.
11. **Three `set.seed()` sites added in the range do not restore the caller's RNG stream**
    (`R/fs_mr_field_uniform.R:93`, `R/fs_mr_inference.R:864`, `R/run_subgroup_sims.R:512`). Reported in §3.6(a);
    nothing changed. Each fires only when the caller supplies `seed` (or `benchmarks`).
12. **Three non-ASCII lines were added in the range**, all in `R/fs_bias_coverage.R` (`:238`, `:239`, `:251`, by
    `72bd714a`) — the file named by the carried-over certification WARNING. Reported; not fixed.
13. **Five `# [delivery sentinel: …]` comments added in the range remain in `R/` at `PIN`** (§3.6(g)).
    Recorded; not fixed.
14. **Two `--as-cran`-uncovered `R/` commits, and fourteen full-suite-uncovered ones.** The last certification run
    was at `b9a32705`; `1719056f` and `f9b794f6` followed it (§3.5). The last standalone `devtools::test()`
    predates fourteen `R/` commits, and its tree SHA is not recorded in the report that carries its counts.
    Stated as coverage, not as a recommendation; nothing was run.
15. **Twenty-four `R/` commits have no `NEWS.md` entry** (§3.3), including ten newly exported functions, five new
    S3 methods, and one commit classified *changes the method* (`a6702fd8`). Recorded as a fact about the
    development section at `PIN`.

## 8. Commits landed by other sessions during this task

**None.** `git log $PIN..HEAD` at closeout contains only this task's own two commits. No `.git/index.lock` was
ever observed during the run, so the §1 wait loop never engaged.
