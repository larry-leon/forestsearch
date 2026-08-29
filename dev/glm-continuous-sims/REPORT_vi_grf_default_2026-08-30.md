# REPORT — `vi.grf.min = NULL` as the default: blast radius — **STOPPED AT §2 GATE, outcome (B); no edit made**

**Task:** `dev/tasks/cc_task_vi_grf_default_2026-08-30.md` (commit `5332f00f`)
**Repo / branch:** `~/Documents/GitHub/forestsearch` (pop-os), `feature/glm-extension`
**Ran after:** `cc_task_oc_wrapper_confs_2026-08-30.md` (commits `7e4a5dfd`, `58024dfe`, `49f6e1ce`), not concurrently.
**Category as executed:** read-only. **No file under `R/`, `man/`, `NAMESPACE`, `DESCRIPTION` or `NEWS.md` was edited.** Version stays `0.2.6`. Artefacts: this report; a scratch RNG check (not committed, reproduced verbatim in §2.5).

---

## ⚠ GATE FAILURE AT THE TOP — §2 verdict **(B)**

**Reference-producing callers rely on the default.** In the two application documents whose 0.2.0 reproduction was verified (`REPORT_repro_check_vs_0.2.0_2026-08-27.md`), three `forestsearch()` fits do **not** pass `vi.grf.min`:

| document | call | passes `vi.grf.min`? |
|---|---|---|
| `quarto/applications/gbsg/analysis_gbsg_survival_multimethod.qmd` | `fs <- forestsearch(` L693 | yes, `-0.2` (L716) |
| same | `fs_dina_screen <- forestsearch(` L1263 | yes, `-0.2` (L1286) |
| same | **`fs_dina <- forestsearch(` L1502–1534** | **no** |
| same | **`fs_grf <- forestsearch(` L1743–1777** | **no** |
| `quarto/applications/actg175/analysis_actg175_continuous_compare_all.qmd` | `compare_selection_rules(` L286 (pass-through `...`) | yes, `vi.grf.min = fs_vi_grf_min` (L346; `fs_vi_grf_min <- -0.2` L146) |
| same | **`fs_anchor <- forestsearch(` L404–444** | **no** |

And the **committed-payload simulation drivers** under `quarto/simulations/gbsg_redux/` (370 tracked `.rds` payloads) call `do.call(forestsearch, c(base_args, method_args))` (e.g. `sim_fs_maxeffCons_fb_mr_m1_h10_knoise0_n500_batch_1_50.qmd` L489) with **no `vi.grf.min` anywhere in the file** — every one of the 56 non-legacy `gbsg_redux` drivers relies on the default. (The ACTG175 MD40 drivers behind the OC-wrapper payloads pass `-0.2` explicitly, L215 / L568 — those are safe.)

Per §2's gate, the edit was **not made**. §3–§5 were not run. What the change would require before it could be made is in §4 below; that is a compute go/no-go for Larry.

---

## 1. Provenance (§1 — GATE passed)

```
pop-os · /home/larryleon/Documents/GitHub/forestsearch · feature/glm-extension
HEAD at start 49f6e1ce · git status --porcelain: empty
git log -6: 49f6e1ce 58024dfe 7e4a5dfd f2ca69e6 b8bd0f0a 69ccc84d   (confs task's three commits present)
packageVersion forestsearch 0.2.6
```

`ls ~/Downloads/cc_task_vi_grf_default_2026*`: exactly one match, `cc_task_vi_grf_default_20260830.md` (hyphens stripped) → `dev/tasks/cc_task_vi_grf_default_2026-08-30.md`, committed **`5332f00f`**.

`devtools::test()` parity baseline (unchanged tree, 0.2.6): `[ FAIL 0 | WARN 31 | SKIP 3 | PASS 4737 ]` (the same figures as the 0.2.6 close-out)

---

## 2. Blast radius (read-only)

### 2.1 Every caller of `forestsearch()` and whether it passes `vi.grf.min`

Method: every tracked file under `R/`, `tests/`, `vignettes/`, `man/`, `quarto/` (`.R`, `.qmd`, `.Rmd`, `.Rd`; `.html` and `quarto/.quarto` excluded), scanned for call sites (`forestsearch(`, `forestsearch::forestsearch(`, `do.call(forestsearch, …)` on non-comment lines) and for `vi.grf.min` in the call span or the argument object the call consumes. Prose and roxygen mentions were discarded by reading each hit.

**`R/` — programmatic callers (all forward a caller-built argument list; none injects `vi.grf.min` itself):**

| file:line | how | passes `vi.grf.min`? |
|---|---|---|
| `R/bootstrap_analysis_dofuture.R:614` | `do.call(forestsearch, args_FS_boot)` — inherited from the parent fit's `args_call_all` | as the parent did |
| `R/forestsearch_cross_validation.R:479`, `:1001` | `do.call(forestsearch, cv_args)` — from `args_call_all` | as the parent did |
| `R/compare_selection_rules.R:561` | `do.call(forestsearch, c(…, list(...)))` — pass-through | only if the document passes it |
| `R/calibrate_null_correction.R:355`, `R/mrct_simulation.R:466`, `R/fpr_calibration.R:297`, `R/fs_fdr_report.R:362` | `do.call(forestsearch, <fs_args>)` — caller-supplied list | as supplied |
| `R/run_simulation_analysis.R:824` | `do.call(forestsearch, fs_args)`, `fs_args` from `default_fs_params_*()` | **yes, `-0.2`** — `R/run_simulation_analysis.R:68` (§2.3), the one package-internal site that pins it |
| all other `R/` hits (≈ 40) | roxygen `@examples` / error-message text | n/a |

**`tests/testthat/`** — five files call `forestsearch()` directly; **none passes `vi.grf.min`** (all rely on the default): `test-glm-pipeline-integration.R` (L27, L80, L136), `test-input-validation.R` (L22, L239, L283, L321), `test-mr-inference.R` (L137), `test-return-shape-contract.R` (L159), `test-search-reproducibility.R` (L38). Others reach it through `forestsearch_bootstrap_dofuture()` / CV wrappers on those fits.

**`vignettes/`**: `vignettes/forestsearch.Rmd` L290 passes `-0.2` (L324, "default -0.2"); `vignettes/articles/methodology.Rmd` L941 and L1052 pass `-0.2` (L975, L1086); `vignettes/articles/extreme_subgroups.Rmd`'s only hit (L1286) is a prose table row, not a call.

**`man/` examples**: 58 `.Rd` files contain `forestsearch(` in examples; none passes `vi.grf.min`; all but a handful are `\dontrun{}` (the runnable ones — `dina_frontier`, `dina_subgroup`, `fs_oc_*`, `default_sim_params`, `interpret_search_config`, `mr_estimates_table`, `test_hte_from_forestsearch` — do not execute a search on the default path or take a fitted object).

**`quarto/` — active (non-`_archive`, non-`_broken`, non-`legacy`) files with a real call site:**

| area | files with a call | pass `vi.grf.min` at every call | rely on the default at ≥ 1 call |
|---|---:|---|---|
| `quarto/applications/gbsg/` | 5 | — | **all 5**: `analysis_gbsg_survival_multimethod` (2 of 4 fits), `_effMaxSG` (0 of 5), `_frozen_family` (0 of 4), `_sgfocus`, `_maxeff_mrconfirm` (`do.call(forestsearch, confirm_args)`, L285) |
| `quarto/applications/actg175/` | 7 | — | **all 7**: `analysis_actg175_continuous_compare_all` (anchor fit L404), `_binary_multimethod_{fixed_family,frontend,psi_v2_2,psi_v3a}` (fits 3–4 of 4: `fs_dina`, `fs_grf`), `_binary_sgfocus` (calls 1–2 of 3), `template_actg175_continuous` |
| `quarto/applications/count_data_demo.qmd`, `validation_glm_simulation_study.qmd`, `validation_hte_tests_crump.qmd` | 3 | — | all 3 (no `vi.grf.min` in the file) |
| `quarto/simulations/` (all 332 tracked non-legacy `.qmd` drivers, file-level scan for `vi.grf.min =`) | 332 | **13** — the ACTG175 MD40 family (`sim_fs_maxeffCons_mr_md40_*`, `sim_fs_maxeffCons_fb_mr_md40_*`, `mr_coverage_sweep_md_harm`, `actg175_continuous_simulations`), `binary_020/maxeffCons_mr_coverage_sweep_or075`, `binary_methods/…_template`, `survival/actg175_survival_benefit_simulations`, `survival/gbsg_poisson_simulations`, `gbsg/reproduce_table1_m1` | **319** — every `gbsg_redux/*` driver (`fs_t1_t2_*`, `grf_t1_t2_*`, `sim_fs_{maxeffCons,maxcons,effMaxSG,effMinSG}_fb_mr_*`, `mr_coverage_sweep_*`), `gbsg/_sim_mr_coverage*`, and the `actg175/binary/*mr_coverage_sweep*` drivers (`do.call(forestsearch, c(base_args, …))`, L424–451, `base_args` without it) |
| `quarto/GuoHe/` 6, `quarto/qc/` 5, `quarto/resampling/` 5, `quarto/dina/` 8, `quarto/calibration/` 6, `quarto/guides/` 15, `quarto/methodology/` 6, `quarto/extreme_subgroups/` 7, `quarto/grf/` 2, `quarto/smoke-tests/` 1 (files with a call) | 61 | 9 of them mention `vi.grf.min` somewhere in the file | the other 52 never mention it |

The full per-file table (518 rows, with call line numbers and every `vi.grf.min` line) was generated by the scan and is reproducible from the command in this report's history; it is not committed.

### 2.2 The applications specifically

Verified in `REPORT_repro_check_vs_0.2.0_2026-08-27.md` (both `.qmd` at `cf4d6432`, `.html` at `43b051b6`, "source-identical" reproduction): **`analysis_gbsg_survival_multimethod.qmd`** — `fs` and `fs_dina_screen` pass `-0.2`; **`fs_dina` (L1502–1534) and `fs_grf` (L1743–1777) do not** (their argument blocks carry `seedit`, `sg_focus`, engine arguments; the in-code comment lists "use_grf / use_lasso / use_twostage / fs.splits / maxk / …" as inherited defaults). **`analysis_actg175_continuous_compare_all.qmd`** — the `compare_selection_rules()` run passes `vi.grf.min = fs_vi_grf_min` (`-0.2`); **the `fs_anchor` fit (L404–444) does not.** Siblings: `analysis_gbsg_survival_effMaxSG.qmd` and `_frozen_family.qmd` (both with rendered `.html`) never mention `vi.grf.min`; the four `analysis_actg175_binary_multimethod_*.qmd` pass it in fits 1–2 (L555, L1018 in `psi_v3a`) and not in `fs_dina` / `fs_grf` (L1199–1232, L1421–1456).

Note on what "relies on the default" means for a fit with `use_grf = FALSE`, `use_dina = FALSE`: Section 5 still runs at `-0.2` (it is gated on `vi.grf.min`, not on `use_grf`), fits a causal forest per call, and **orders** the cut columns by importance. §2.5's fixture shows the ordering reaching the output even when the selected rule is the same set of clauses: `sg.harm` = `{biomarker_hi} {age <= 47}` at `-0.2` and `{age <= 47} {biomarker_hi}` at `NULL`.

### 2.3 Both `-0.2` sites, verbatim

```r
# R/forestsearch_main.R
712:  #' @param vi.grf.min Numeric. Minimum GRF variable importance. Default -0.2.
1216:                         vi.grf.min = -0.2,
# R/run_simulation_analysis.R  (inside default_fs_params_*(); explicit caller choice, NOT to be changed here)
68:    vi.grf.min               = -0.2,
```

### 2.4 The `max_n_confounders` coupling

`R/forestsearch_main.R` L2752–2832: the whole block is `if (!is.null(vi.grf.min)) { … } else { conf.screen <- FSconfounders.name }`; inside it,

```r
2814:      vi_max <- max(vi.cs2[, 3])
2815:      if (vi_max > 0) {
2817:        vi_ratio <- vi.cs2[, 3] / vi_max
2818:        selected.vars <- which(vi_ratio > vi.grf.min)
2819:        conf.screen <- conf.screen[selected.vars]
2820:        conf.screen <- conf.screen[seq_len(min(length(conf.screen), max_n_confounders))]
2821:      }
```

The truncation sits inside `if (vi_max > 0)`, inside `if (!is.null(vi.grf.min))`. **With `vi.grf.min = NULL` it never executes: the default change would make `max_n_confounders` inert by default.** Formal default `max_n_confounders = 1000` (L1167). Callers setting it to anything else: **none** in `R/`, `tests/`, `vignettes/`, `quarto/` (the only mentions outside `forestsearch_main.R` are the `.Rd` and `args_call_all` bookkeeping). Reported as a finding; its placement is not touched.

### 2.5 RNG — does Section 5 consume R's stream?

From source: Section 5 calls `apply(df[, …], 2, as.numeric)`, `fit_causal_forest_glm(X, Y_vi, Treat, is.RCT, seedit = seedit, tune_grf)` → `grf::causal_forest(…, seed = seedit, tune.parameters = "none")` (`R/grf_helpers.R` L115–140; survival path `grf::causal_survival_forest(…, seed = seedit)`), `grf::variable_importance()`, `order()`. grf seeds its own C++ RNG from `seed`; no `set.seed`, `sample`, `rnorm`, `runif` in the block.

Bounded check (scratch script; `.make_continuous_data(N = 300, seed = 42)` from `tests/testthat/helper-synthetic-dgm.R`; `forestsearch()` continuous / MD, `use_lasso = use_grf = FALSE`, `maxk = 2`, `n.min = 30`, `fs.splits = 20`, `seedit = 1`, sequential; `set.seed(20260830)` before each call):

```
RNG stream changed by the call:  vi=-0.2: TRUE  vi=NULL: TRUE
post-call .Random.seed identical between -0.2 and NULL: TRUE ; next runif identical: TRUE
sg.harm -0.2: {biomarker_hi} {age <= 47}  | NULL: {age <= 47} {biomarker_hi}
confounders.evaluated -0.2: age <= 55, age <= 47, age <= 62, biomarker <= -0.1, biomarker <= -1, biomarker <= 0.6, biomarker_hi, sex
confounders.evaluated NULL : (identical)
```

`forestsearch()` does advance R's stream (the consistency stage's splits), but **identically** under `-0.2` and `NULL`: Section 5 consumes none of it. **Outcome (C) does not apply.** The change is a reordering, not a stream shift — and the reordering is visible in the output (clause order of `sg.harm`) even on this tiny fixture.

### 2.6 Verdict

**(B).** Reference-producing callers rely on the default: two fits in the verified gbsg application, one in the verified actg175 application, both rendered sibling gbsg documents, and every `gbsg_redux` committed-payload driver. **Stopped before editing.**

---

## 3. §3 — not run (gate)

No change to `R/forestsearch_main.R`, roxygen, `man/`, `NAMESPACE`, `DESCRIPTION`, `NEWS.md`. `R/run_simulation_analysis.R:68` untouched, as instructed regardless.

## 4. What the change would require — for Larry's go/no-go

1. **Re-verification of the two 0.2.0-verified applications** after the default flips: `analysis_gbsg_survival_multimethod.qmd` (`fs_dina`, `fs_grf`) and `analysis_actg175_continuous_compare_all.qmd` (`fs_anchor`) would run Section 5 off instead of on. The reproduction check is ~19 min of render (per the repro report). The sibling `_effMaxSG` and `_frozen_family` gbsg documents, whose `.html` is tracked, would also change source-vs-render status.
2. **Or** pin `vi.grf.min = -0.2` explicitly at those three call sites first (an application-document edit, out of scope here) so the default change is inert for them; then the re-render is a no-op check rather than a re-verification.
3. **The `gbsg_redux` payloads** (370 `.rds`) were produced with the ordering on; re-running any of those drivers after the change would not reproduce them. Either their `base_args` gain `vi.grf.min = -0.2`, or the payloads are accepted as belonging to the pre-change default — a decision, not a compute item.
4. Tests: the five test files above rely on the default; after the change they exercise the `NULL` path. `test-search-reproducibility.R` in particular pins a fixed-seed result and may change (unknown until run).

## 5. Close-out

`git status --porcelain` before staging: this report only. Staged by explicit path; committed; hash in the follow-up line. **No push. No install** (nothing in `R/` changed; version stays 0.2.6).

Commits: `5332f00f` task doc; `06728783` this report; the hash line in the follow-up commit.

## 6. Verdict (ten lines)

1. §1 gate passed; task doc committed `5332f00f`; test baseline recorded above.
2. The formal is `vi.grf.min = -0.2` at `R/forestsearch_main.R:1216` (roxygen L712); the only package-internal explicit `-0.2` is `R/run_simulation_analysis.R:68` (left alone).
3. At `-0.2` the screen drops nothing (`vi_ratio ∈ [0,1]`, block guarded by `vi_max > 0`) but reorders cut columns by a per-replicate causal-forest importance; the reordering reaches `sg.harm`'s clause order on a 300-row fixture.
4. `max_n_confounders` is applied only inside that block (L2820, under `vi_max > 0`, under `!is.null(vi.grf.min)`): the change would make it inert by default; no caller sets it.
5. Section 5 consumes **no** R RNG: `.Random.seed` after the call and the next `runif()` are identical at `-0.2` and `NULL` — a reordering, not a stream shift; (C) does not apply.
6. **(B):** `fs_dina` and `fs_grf` in the verified gbsg document and `fs_anchor` in the verified actg175 document do not pass `vi.grf.min`; nor do the rendered `_effMaxSG` / `_frozen_family` gbsg documents; nor any of the 56 `gbsg_redux` drivers behind 370 committed payloads.
7. Callers that do pin `-0.2`: the ACTG175 MD40 drivers (the OC-wrapper payloads), `default_fs_params_*()`, the two vignettes, the actg175 `compare_selection_rules()` run, and fits 1–2 of the multimethod documents.
8. The five direct-calling test files rely on the default; the two vignettes pin `-0.2`.
9. **No edit made; §3–§5 not run; version 0.2.6; tree unchanged except this report and the task doc.**
10. Re-enabling this task needs Larry's decision on §4 items 1–3 (re-verify ~19 min, or pin `-0.2` at three application call sites first, and what to do about the `gbsg_redux` payloads).
