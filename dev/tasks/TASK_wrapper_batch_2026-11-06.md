# TASK: wrapper batch — P3, P4, P6, P7, P8 in one task (review proposals, bundled)

**Date:** 2026-11-06
**Origin:** REVIEW_wrapper_functions_2026-11-05.md, bundling approved by Larry 2026-11-06.
Standing constraint: functionality preserved; overall improvement.
**Compute:** no simulation campaign. Every fixture is seconds to ~2 minutes; the one full
suite + `check()` at the end is the dominant cost (~25 min). No compute go/no-go applies.
**Git:** commit the task doc on `master`, create work branch `fs-wrappers-batch` from it, do
everything else there. **One commit per item** (B1–B5 below), so review maps commits to items.
Commit as you go. **Do not push.**

## Testing cadence (binding for this task)

1. **No baseline run.** Verify `git diff <P1+P2 gated commit 84afa853>..HEAD -- R tests man
   vignettes inst DESCRIPTION NAMESPACE .Rbuildignore` is empty (the merge added only the
   report, the amended task doc and `_pkgdown.yml`); then reuse
   `REPORT_mrct_trial_analysis_2026-11-06.md`'s final figures as the baseline (check 0/0/0
   OK; suite 0 failures, 3 skips, 38 warnings, 7,962 expectations) and record the reuse.
   If the diff is not empty: STOP and report what changed.
2. **Inner loop:** only the new test files and each item's own fast fixture run between
   edits. Never the full suite or `check()` mid-task.
3. Each gate's fixture runs **once** at its gate.
4. **Full suite + `check()` + `pkgdown::check_pkgdown()`: once, at the final stage only.**

## R/ impact callout (per protocol)

- **New code only.** Five new `R/` files (or fewer, grouped sensibly), tests, roxygen, NEWS,
  `_pkgdown.yml` index entries.
- **Moves existing code: none. Changes behaviour: none. Changes method: none.**
  No existing `R/` file is edited — P8's build helper is a thin wrapper, not an argument on
  `inject_mrct_structure()`; `ahr_design()`, `mrct_region_sims()`, `mr_estimates_table()`
  and `fs_declaration_calibration()` are all untouched.

## Out of scope

Any edit to an existing `R/` file; any vignette or article edit (later pass); any
winner-variant quantity anywhere; any new statistic.

---

## Stage 0 — verify structure from source (one pass for all items)

Read and record, with file:line references (prior summaries and the KB are not evidence):

- **For B1 (P3):** `mrct_region_sims()`'s return — its class (if any), the per-trial columns
  (`sim`, `any_found`, `hr_itt`, `hr_test`, the `TP_/FP_/FN_/TN_{nonregion,region}` columns),
  and the `"estimators"` attribute's shape under `record = "full"` (estimator keys
  `itt_all`, `itt_region`, `region_Ghat`, `region_Gstar`; `hr`/`lo`/`hi` fields). Decide per
  D-B1.1 whether a `summary()` method or a plain function results.
- **For B2 (P4):** `ahr_population()` / `ahr_limit()` / `ahr_design()` signatures and return
  fields (`hr`, `marghr_inf`, `ahr`, `cde`; the `$prob` names `p_itt_ub_lt1`,
  `p_region_itt_gt1`, `p_trigger`, `p_region_Gstar_lt1`, `p_trigger_and_Gstar_lt1`), and
  the set-specification forms `ahr_limit()` accepts (named set strings and logical vectors).
- **For B3 (P6):** the `forestsearch_bootstrap_dofuture()` result fields the articles read
  (`H_estimates$H2`, `Hc_estimates$H2`) and whether the object stores its own wall time.
- **For B4 (P7):** `fs_declaration_calibration()`'s return (`c0$table` with `kappa_hat`,
  `fw_size`; `c0$source`; `kappa_hat`/`fw_size` at the claim threshold; `T_hat`;
  `admitted_calibrated`) and the fit's `mr_inference$declaration_field` fields
  (`family_id`, `meta$selected_index`), plus what the fit must carry
  (`keep_declaration_field`).
- **For B5 (P8):** `get_dfpred()`'s signature (the `version = 2` idiom), and the shipped
  extdata layout (`inst/extdata/mrct-consistency/`, `PROVENANCE.md`) with the dev-tree /
  installed resolution the articles use.

**Gate 0 (stop-on-failure):** every field named above verified or its absence reported; the
baseline-reuse diff empty.

---

## The items

### B1 (P3) — core summary of an `mrct_region_sims()` campaign

The committed core-table rows, computed from a `record = "full"` result; **transplant** the
`demo-summary` chunk of `vignettes/articles/mrct-consistency-simulations.qmd` as the
definitions — do not re-derive:

trial power P(ITT ub < 1.0) over all trials; trials reaching the trigger (count); subgroup
declared (% of triggered); Sens/Spec of Ĥ vs H*, non-region and region (averaged over the
declared triggered trials); declared & region(Ĝ) < 1.0 (% of triggered); region(G*) < 1.0
(% of triggered); declared & region(Ĝ) ≥ 1.0 (% of triggered); mean share of region patients
kept in Ĝ (over declared triggered). Counts and binomial SEs carried as the core table does;
`as.data.frame()`-ready.

- **D-B1.1:** if the return carries a dedicated S3 class, implement as a `summary()` method
  (+ a `print` for the summary); otherwise export `mrct_core_summary(res)`. Never modify
  `mrct_region_sims()` to add a class.
- **D-B1.2:** requires `record = "full"`; informative error otherwise. Denominators are fixed
  by the definitions above — all conditional rows are % of triggered.
- **Gate B1:** on a small synthetic campaign (small DGM, `n_sims` 20–40, small `n`,
  `record = "full"`; a couple of minutes at most), the function's every cell equals the
  transplanted `demo-summary` computation run on the same object, exactly.

### B2 (P4) — analytic-layer tables

```r
ahr_truths_table(x, sets, b, h, analysis_time, max_entry, ...)  # x: dgm/pop or named list
ahr_design_table(x, b, h, n, ...)                               # + the 7-row decomposition
```

Tidy data frames: cell × set × margHR(analysis time) / margHR(∞) / AHR / CDE; and cell ×
the five probabilities plus P(G* < 1 | trigger) plus the seven consistency rows of the
application article's `pre-table` chunk (identities among the rows asserted inside the
function). `h` scalar or named per cell; `sets` as `ahr_limit()` accepts them (named strings
or logical vectors). Pure loops over the existing `ahr_*` calls — no quadrature code of
their own.

- **Gate B2:** build S2 (and the other `r = −log 4` cells as needed) with the shipped
  `latent_R` via the simulations article's `build-fn` chunk (transplant; builds are ~1 s
  each), run both functions, and compare against the committed analytic records
  (`predictions_full.csv`, `predictions_S1.csv`, `predictions_n700.csv`) at 1e-9 — the
  article's own tolerance, its check code transplanted.

### B3 (P6) — naive / bootstrap / multiplier comparison, and cost

```r
mr_method_comparison(fits, boots, est.scale = "hr", digits = 2)
mr_cost_table(fits, boots, boot_seconds = NULL, project_B = 1000L)
```

`fits` a named list of `forestsearch()` fits with MR, `boots` the matching
`forestsearch_bootstrap_dofuture()` results. Output rows exactly as the committed
`three_columns()`/`side_by_side` build them (Method, SG, n (%), naive with CI, bootstrap
**point estimate only**, multiplier with CI); plain data frames, the caller styles them.
`mr_cost_table` takes `boot_seconds` as a named vector (auto-read only if Stage 0 found the
object stores wall time). Roxygen carries the committed captions' framing verbatim in
substance: the bootstrap point at small B is an illustration without an interval, and for
model-generated families FB and MR target different estimands — the table compares displays,
it does not adjudicate one quantity.

- **Gate B3:** on a fast fixture (one small synthetic survival `forestsearch()` fit with MR
  and a `nb_boots = 2` bootstrap; ~1–2 min), equality with the transplanted `three_columns()`
  and cost-chunk code from the committed `survival-post-selection` article, cell for cell.

### B4 (P7) — declaration-calibration table

```r
fs_declaration_table(fit, alpha = 0.10, c0 = NULL, c0_primary = NULL)
```

Wraps the two `fs_declaration_calibration()` calls (the `c0` grid and the claim threshold)
plus the `declaration_field` reads of the ACTG article's `s17_column()`: region N; T of the
region and the family maximum; per-`c0` κ̂, implied level 2Φ(κ̂)−1, family-wise size at the
executed screen; the claim-threshold κ̂ and size; the admits verdict. Tidy data frame (one
row per `c0` plus the claim-threshold row), from which the article's S17 entries are exactly
recoverable. Informative error when the fit lacks the declaration field.

- **Gate B4:** on a small fit carrying `keep_declaration_field` with a short `c0` grid,
  equality with the transplanted `s17_column()` entries on the same fit, exactly.

### B5 (P8) — minor helpers (one commit, three small functions)

1. `mrct_record(file)`: the shipped-record reader — dev-tree/installed resolution,
   `PROVENANCE.md` check, `read.csv` as the articles' `rd()` does.
   *Gate:* identical content to the articles' `rd()` on three shipped files; informative
   error when the directory or provenance is absent.
2. `fs_in_subgroup(fit, data)`: logical membership, retiring the
   `get_dfpred(..., version = 2)$treat.recommend == 0` idiom.
   *Gate:* equals the raw idiom on a small fit, and equals the fit's own
   `df.est`/`df.test` classification (the articles' stopifnot, transplanted).
3. `mrct_build_dgm(data, ..., m = 1)`: thin wrapper over `inject_mrct_structure()`
   (+ `scale_dgm_hazard()` when `m != 1`) with warnings captured to a `build_warnings`
   attribute — the `withCallingHandlers` pattern of the three committed build chunks,
   transplanted. **A wrapper, not an argument on the existing function.**
   *Gate:* at a fixed seed, output identical (attribute aside) to the direct
   `inject_mrct_structure()` + `scale_dgm_hazard()` calls; warnings captured on a fixture
   known to warn (the rep5-trial build of the application article's `pre-build` chunk,
   which records censoring-model warnings) and empty where the build is clean.

---

## Shared decisions

- **D1:** every output is a plain data frame (plus a light print method where an item says
  so); no gt/kable objects — articles style their own tables.
- **D2:** evaluated construction set only, everywhere MR quantities appear (naive, IJ
  two-term, field, field-s, Bonferroni); `mr_values()` used internally where convenient;
  bounds labeled by location, no significance framing in any printed output.
- **D3:** roxygen with runnable examples (`\donttest` over ~5 s); every new export added to
  `_pkgdown.yml` beside its family; one NEWS entry listing the batch.
- **D4:** tests committed per item (the gate's transplant-equality where the fixture is fast
  enough for the suite; `skip_on_cran()` over ~5 s; the B2 committed-record comparison and
  any slow gate run once as gate artifacts in the report, not as suite tests).
- **D5:** fixtures are built in-test (synthetic, seeded); no baseline RDS enters the tree.

## Final stage — confinement, preservation, report

1. **Diff confinement:** changed existing files are exactly `NAMESPACE`, `NEWS.md`,
   `_pkgdown.yml`; everything else in the diff is a new file. Any other existing file in the
   diff fails the gate.
2. Full suite once: green, skip list unchanged from the reused baseline, existing
   expectations intact.
3. `devtools::check()` once: 0 errors, 0 warnings, notes ≤ baseline.
   `pkgdown::check_pkgdown()` clean.
4. **Report** `REPORT_wrapper_batch_2026-11-06.md` beside the prior two: the baseline-reuse
   record; per-item sections B1–B5 with the Stage-0 field notes and **verbatim** gate
   outputs; the commit-per-item list; the final suite/check comparison; OPEN ITEMS for
   anything non-blocking. No push.
