# REPORT — the threshold naming collision: documentation only

**Task:** `dev/tasks/TASK_threshold_naming_docs_2026-09-18.md` (committed alone at `216f3405`).
**Date:** 2026-09-18. **Branch:** `feature/glm-extension`. **Machine:** `pop-os`.

**Commits produced by this task**

| Commit | Contents |
|---|---|
| `216f3405` | the task document, as received, committed alone before anything else |
| `cc1714e0` | **Part 1**, sections 3.1-3.3: the threshold vocabulary, the per-estimand resolved defaults, the `fpr_calibration()` bridge; plus `NEWS.md` |
| `0671525f` | **Section 3.4** (offered/separable), committed alone so it can be reverted without touching `cc1714e0` |
| `458cef48` | **Part 2** (separable), the DINA frontier caps |
| this report | filed beside the audit it follows from |

**Placement.** Filed at `quarto/simulations/actg175/binary_020/REPORT_threshold_docs_2026-09-18.md`,
beside `REPORT_criterion_and_defaults_audit_2026-09-18.md`, the audit this task follows from. The
audit recorded that this repository has no root-level `REPORT_*` location; its reports live in
campaign directories. Same reasoning, same directory.

**Category.** Documentation only. No behaviour change, no default change, no method change, no code
line altered, no install, no push. The only computation is `devtools::document()` and the two
`R CMD check --as-cran` runs required by post-conditions 2 and 4.

**Source-first rule (task section 1).** Every value below was read from the current source before it
was written into any roxygen block; nothing was carried over from the task document or from the audit
report. Two citations in the audit turned out to be off by one to two lines (section 8, finding 4),
which is the reason the rule exists.

---

## 1. Every roxygen block changed

Line numbers are at HEAD (`458cef48`), after the edits.

### Part 1, sections 3.1-3.3 (`cc1714e0`)

| # | File and line | Function (exported) | What the block now states |
|---|---|---|---|
| 1 | `R/forestsearch_main.R:492` | `forestsearch()` | `effect.threshold` is `c1`, the screening threshold on a candidate's own effect; names `c2` and `p*` and says how each differs |
| 2 | `R/forestsearch_main.R:504` | `forestsearch()` | `consistency.threshold` is `c2`, **an effect threshold**, applied within each split half; explicitly not a proportion, not the consistency rate; names `pconsistency.threshold` as its companion |
| 3 | `R/forestsearch_main.R:517` | `forestsearch()` | `hr.threshold` is the legacy name for `effect.threshold` (`c1`) |
| 4 | `R/forestsearch_main.R:521` | `forestsearch()` | `hr.consistency` is the legacy name for `consistency.threshold` (`c2`), the per-split **effect** threshold |
| 5 | `R/forestsearch_main.R:705` | `forestsearch()` | `pconsistency.threshold` is `p*`, a proportion in `[0, 1]`; a rate, never remapped per estimand |
| 6 | `R/forestsearch_main.R:1054` | `forestsearch()` | **new `@section Threshold vocabulary and resolved defaults`** — the three-way `describe` list, the one-line contrast, and the per-estimand table (section 2 below) |
| 7 | `R/subgroup_search.R:22` | `subgroup.search()` | `hr.threshold` is `c1`; states that neither `c2` nor `p*` is applied here, because this function performs the screening stage only |
| 8 | `R/subgroup_consistency_main.R:79` | `subgroup.consistency()` | `hr.threshold` is `c1`, re-applied to the supplied `hr.subgroups` table on the estimand's comparison scale |
| 9 | `R/subgroup_consistency_main.R:84` | `subgroup.consistency()` | `hr.consistency` is `c2`, an effect threshold required in **each** half; not a proportion |
| 10 | `R/subgroup_consistency_main.R:89` | `subgroup.consistency()` | `pconsistency.threshold` is `p*`, a proportion in `[0, 1]`; `c2` sets the bar, `p*` counts how often it is met |
| 11 | `R/consistency_resample.R:346` | `consistency_resample()` (and `consistency_resample_compare()` via `@inheritParams`) | `hr.consistency` is `c2`, an effect threshold on the half-sample effect; notes that this function returns the rate rather than thresholding it, so it takes no `p*` |
| 12 | `R/subgroup_consistency_helpers.R:228` | `run_single_consistency_split()` | `hr.consistency` is `c2`, an effect threshold; not a proportion |
| 13 | `R/subgroup_consistency_helpers.R:1295`, `:1300` | `evaluate_subgroup_consistency()` | `c2` and `p*`, both roles stated |
| 14 | `R/subgroup_consistency_helpers.R:1575`, `:1580` | `evaluate_consistency_twostage()` | `c2` and `p*`, both roles stated |
| 15 | `R/fpr_calibration.R:49`, `:55` | `fpr_calibration()` | the vocabulary bridge (section 3 below) |
| 16 | `R/fpr_approximation.R:27` | `fpr_approx()` | `c1` is `forestsearch()`'s `effect.threshold`; notes that `c2` and `p*` do not appear because this approximation covers the screening stage only |
| 17 | `R/fs_oc_predict.R:97`, `:101` | `fs_oc_predict()` | `c1` / `c2` named against the `forestsearch()` arguments; `c2` stated as an effect threshold, not a proportion |
| 18 | `R/fs_oc_predict.R:113` | `fs_oc_predict()` | `pconsistency` is `p*`, a consistency-*rate*, not an effect threshold |
| 19 | `R/mrct_simulation.R:40`, `:43`, `:47` | `mrct_region_sims()` | `c1`, `c2`, `p*` roles stated |

Plus `NEWS.md:97-109`, one entry under documentation (post-condition 7).

`fs_oc_grid()` inherits from `fs_oc_predict()` (`R/fs_oc_grid.R:49`,
`@inheritParams fs_oc_predict`) and `fs_oc_invert()` inherits from `fs_oc_grid()` (`:225`), so both
pick up row 18's `pconsistency` wording without `R/fs_oc_grid.R` being edited -- which is why
`man/fs_oc_grid.Rd` and `man/fs_oc_invert.Rd` appear in the `man/` diff. Neither is a separately
edited block.

### Section 3.4 (`0671525f`)

| # | File and line | What it states |
|---|---|---|
| 20 | `R/forestsearch_main.R:997` | new `\item{threshold_config}` in `forestsearch()`'s `@return`: the component was previously **undocumented**. Says to read `$scale` first, and that on the survival path `$screening` is `log(hr.threshold)` — the admission / multiplier-resampling scale — while the candidate search compares the fitted hazard ratio against `hr.threshold` on the **natural** scale (`$screening_natural`); on the GLM paths `$screening` is the value the search compares against directly. `$pconsistency` has no scale. |

### Part 2 (`458cef48`)

| # | File and line | What it states |
|---|---|---|
| 21 | `R/dina_subgroup.R:1040` | `dina_frontier()` `@details`, "Two caps": both caps act on this function's **return value** and nothing else; they trim a table of **single cuts**, never conjunctions; `dina_subgroup()` does not take them; no effect under `subgroup_method = "dina"`; they bite only where the output is used directly, including `use_dina = TRUE, dina_args = list(selected_only = FALSE)` |
| 22 | `R/dina_subgroup.R:1094` | `@param max_per_covariate`: a report-trimming cap on the returned table, not a search control |
| 23 | `R/dina_subgroup.R:1098` | `@param max_subgroups`: the "DINA analog of forestsearch's `max_subgroups_search`" phrase **removed** and replaced by the explicit contrast — `max_subgroups_search` truncates the pool forestsearch **evaluates** and defaults `Inf`; `max_subgroups` trims the rows this function **returns** and defaults the finite `10L` |
| 24 | `R/forestsearch_main.R:321` | `forestsearch()`'s `@param dina_args`: `m_diff`'s "ignored on that path" sentence extended to all **seven** frontier keys by name (section 4 below) |

### Blocks found and deliberately **not** changed

| File and line | Function | Why |
|---|---|---|
| `R/subgroup_search.R:537` | `evaluate_combination_with_status()` | `@keywords internal`, not exported. Task section 3.1 scopes the work to exported functions. |
| `R/subgroup_consistency_helpers.R:1173-1174` | `.make_eval_subgroup_consistency()` | `@noRd`; generates no documentation. |
| `R/forestsearch_helpers.R:2353-2359` | `.fs_resolve_admission()` | `@noRd`; generates no documentation. |
| `R/calibrate_null_correction.R:146-147`, `:217-219` | `predict_fpr_corrected()`, `run_null_calibration()` | Neither is in `NAMESPACE`; not exported. |
| `R/truefind_asymptotic_glm.R:176`, `:179` | — | Not exported. |
| `R/fpr_calibration.R:61` | `fpr_calibration()` | `Must satisfy \code{c2 <= c1}.` — left byte-identical and untouched. See section 3. |

---

## 2. Every value documented, and the source line it was verified against

Read from source at the time of writing, not taken from the task document or the audit.

| Value written into the documentation | Source line verified against | What the source says |
|---|---|---|
| `c1` default `1.25` | `R/forestsearch_main.R:1354` | `hr.threshold = 1.25` (formal) |
| `c2` default `1.0` | `R/forestsearch_main.R:1355` | `hr.consistency = 1.0` (formal) |
| `p*` default `0.90` | `R/forestsearch_main.R:1361` | `pconsistency.threshold = 0.90` (formal) |
| `RD` resolves `c1 = 0.05` | `R/forestsearch_main.R:1904` | `effect_threshold <- switch(effect_measure, RD = 0.05, IRD = 0.01)` |
| `IRD` resolves `c1 = 0.01` | `R/forestsearch_main.R:1904` | same line |
| `RD` / `IRD` resolve `c2 = 0.0` | `R/forestsearch_main.R:1907` | `consistency_threshold <- 0.0` |
| `MD` resolves `c1 = 0.0` | `R/forestsearch_main.R:1956` | `effect_threshold <- 0.0` |
| `MD` resolves `c2 = 0.0` | `R/forestsearch_main.R:1959` | `consistency_threshold <- 0.0` |
| `OR` / `RR` / `IRR` resolve `c1 = log(1.25)` | `R/forestsearch_main.R:1989` | `effect_threshold <- log(effect_threshold)` |
| `OR` / `RR` / `IRR` resolve `c2 = log(1.0)`, i.e. `0` | `R/forestsearch_main.R:1990` | `consistency_threshold <- log(consistency_threshold)` |
| survival `c1 = 1.25` compared on the **natural HR** scale by the search | `R/forestsearch_main.R:1829` (`effect_threshold <- NULL`, assigned only inside the `outcome_type != "survival"` branch), `:3110` (`hr.threshold = if (!is.null(effect_threshold)) effect_threshold else hr.threshold`), `R/subgroup_search.R:132-136` (`screen_threshold` falls through to `hr.threshold`), `R/subgroup_search.R:640` / `:681` (`hr <= hr.threshold`) | for survival `effect_threshold` is `NULL`, so the natural `1.25` reaches the comparison, and `fit_cox_for_subgroup()` returns `exp(beta)` |
| survival `c2 = 1.0` | `R/forestsearch_main.R:1355`, `:2109` | formal `1.0`; `threshold_config$consistency = log(max(hr.consistency, 0.001))` |
| `threshold_config$screening` is `log(hr.threshold)` on the survival path | `R/forestsearch_main.R:2108` | `screening = log(hr.threshold)` |
| `threshold_config$screening_natural` holds the search's value | `R/forestsearch_main.R:2110` | `screening_natural = hr.threshold` |
| `$pconsistency` is `pconsistency.threshold`, unremapped | `R/forestsearch_main.R:2065` (GLM branch), `:2113` (survival branch) | `pconsistency = pconsistency.threshold` in both, with no transform anywhere |
| a split is consistent when **both** halves clear `c2` | `R/subgroup_consistency_helpers.R:298-299` (GLM, `>=`), `:312` (survival, `>`) | `res1$estimate >= c && res2$estimate >= c`; `hr.split1 > c && hr.split2 > c` |
| `fpr_calibration()` passes `c1` / `c2` through under the `forestsearch()` legacy names, on the natural scale | `R/fpr_calibration.R:259-261`, `:307` | `modifyList(fs_params, list(hr.threshold = c1, hr.consistency = c2, ...))`, then `do.call(forestsearch, fs_call)` |
| `fpr_approx()` `c1` default `1.25` | `R/fpr_approximation.R:94` | `c1 = 1.25` (formal) |
| `mrct_region_sims()` `c1 = 0.90`, `c2 = 0.80`, `p* = 0.90` | `R/mrct_simulation.R:175-177` | `hr.threshold = 0.90`, `hr.consistency = 0.80`, `pconsistency.threshold = 0.90` |
| `max_subgroups_search` defaults `Inf` | `R/forestsearch_main.R:1380` | `max_subgroups_search = Inf` (formal) |
| `max_subgroups` defaults `10L`, `max_per_covariate` defaults `3L` | `R/dina_subgroup.R:1163-1164` | `max_per_covariate = 3L`, `max_subgroups = 10L` (formals) |
| the caps trim the returned table | `R/dina_subgroup.R:1268`, `:1281` | `ff <- ff[seq_len(max_per_covariate), ]`; `fr <- fr[seq_len(max_subgroups), ]` |
| `dina_subgroup()` does not take the caps | `R/dina_subgroup.R:318-329` | its formals are `fit, df, covariates, m_diff, n_min, n_min.frac, direction, max_depth, grid_probs, sg_focus, selection_rule, effect_neighborhood, alpha, tau_sign` — neither cap is among them |
| the seven frontier keys | `R/forestsearch_helpers.R:1044-1045` | `frontier_keys <- c("scope", "m_diff", "n_min", "direction", "max_per_covariate", "max_subgroups", "digits")` |
| all seven inert under `subgroup_method = "dina"` | `R/forestsearch_main.R:2342` (branch), `:2459` (`return(out)`), `:2761` (the `use_dina` block, never reached); `R/forestsearch_helpers.R:1461-1462` (the only `dina_frontier()` call on that path, hard-coded `scope = "wide", n_min = n.min`, inside `if (isTRUE(details))`) | `da$frontier` is never referenced on the `subgroup_method = "dina"` path |
| under `use_dina` + `selected_only = TRUE`, six of seven change nothing and only `digits` acts | `R/forestsearch_main.R:2790-2793` (frontier computed), `:2795` (`if (isTRUE(da$selected_only))` -- the selected cut is used and `fr` is discarded), `:2821` (`signif(sgsel$threshold, da$frontier$digits)`), `:2826` (`else fr$cut_expr`) | exactly as documented |

**Not documented because it could not be traced:** nothing. Every value written has a line above.

---

## 3. The `fpr_calibration()` bridge, and the one sentence left alone

`R/fpr_calibration.R:49` and `:55` now say, in both directions, that `c1` is this function's name for
`forestsearch()`'s `effect.threshold` (legacy `hr.threshold`) and `c2` its name for
`consistency.threshold` (legacy `hr.consistency`), both passed through under the legacy names on the
natural scale, and that the consistency *rate* is `pconsistency.threshold`, which travels separately
in `fs_params`.

`R/fpr_calibration.R:61` (was `:51` before the edit) -- `Must satisfy \code{c2 <= c1}.` -- **was already in the source and is left
byte-identical.** It appears in the diff as unchanged context, not as an addition. Task section 2
forbids *adding* a `c2 <= c1` statement, and post-condition 6 greps the diff for added text; nothing
was added. The existing sentence is also **correct for this function**: `fpr_calibration()` enforces
it itself, at `R/fpr_calibration.R:212` -- `"c2 must be <= c1" = c2 <= c1` inside the `stopifnot()` at
`:209`. It is a precondition of `fpr_calibration()`, not of `forestsearch()`, and nothing asserting it
was added to `forestsearch()` or to any other block.

---

## 4. Part 2: what was verified, and the one item not done

Verified from source, not from the cited report (which does not exist — section 8, finding 2):

1. **The caps trim the returned table.** `R/dina_subgroup.R:1268` and `:1281` subset `ff`
   and `fr`, the data frame this function returns. Nothing else reads them.
2. **Single cuts, never conjunctions.** The emitted expression is the canonical
   `"<covariate> <= <threshold>"` built per covariate; the frontier is computed per covariate and
   pooled. There is no conjunction path in `dina_frontier()`.
3. **Not consulted by `dina_subgroup()`.** Its formals (`R/dina_subgroup.R:318-329`) carry neither cap.
4. **No effect under `subgroup_method = "dina"`.** That branch opens at `R/forestsearch_main.R:2342`
   and returns at `:2459`, before the `use_dina` screening block at `:2761` where `da$frontier` is
   used. On the `dina` path the only `dina_frontier()` call is the `details`-time one at
   `R/forestsearch_helpers.R:1461-1462`, which passes hard-coded `scope = "wide", n_min = n.min` and
   not `da$frontier`.
5. **The `max_subgroups_search` contrast.** `Inf` at `R/forestsearch_main.R:1380` against `10L` at
   `R/dina_subgroup.R:1164`; one truncates the evaluated pool, the other trims a returned table.
6. **All seven keys.** `R/forestsearch_helpers.R:1044-1045` lists exactly seven. The extension at
   `R/forestsearch_main.R:321` names them and states a split finer than the task asked for, because
   the source supports the finer statement: all seven are inert under `subgroup_method = "dina"`;
   under `use_dina` screening with `selected_only = TRUE` (the default) six are inert while `digits`
   still rounds the emitted cut's threshold (`R/forestsearch_main.R:2821`); all seven govern the pool
   only under `selected_only = FALSE` (`:2826`).

**Item not done — Part 2's fourth bullet.** "Retitle the `details`-time frontier print as a display of
proposed single cuts shown beside the family counts — not as 'candidates'." That title is **not
roxygen**. It is a string literal in executable code:

- `R/forestsearch_helpers.R:1470` — `lines <- c(lines, "  DINA frontier candidates (per-covariate non-dominated):")`
- `R/forestsearch_main.R:2840` -- `else "frontier candidates"` (the `details` screening-mode label)

Changing either is a code change, which task section 2 forbids and post-condition 1 would fail on.
Per section 2 ("stop and report; do not make it") the item was stopped and is reported here. No
substitute wording was invented elsewhere. The remaining three Part 2 bullets are complete.

---

## 5. Post-condition 1 — the mechanical assertion

The assertion script is `assert_roxygen_only.sh` (written to the session scratchpad, not committed).
For every added or removed line under `R/` in the given diff range, the line body must match
`^[[:space:]]*#'`.

Over the whole task range, `216f3405..HEAD`:

```
changed lines under R/: 267
non-roxygen changed lines: 0
ASSERTION PASS: every changed line under R/ is a roxygen comment line
```

Per commit:

```
-- cc1714e0  (Part 1, 3.1-3.3)
changed lines under R/: 194
non-roxygen changed lines: 0
ASSERTION PASS: every changed line under R/ is a roxygen comment line
-- 0671525f  (section 3.4)
changed lines under R/: 23
non-roxygen changed lines: 0
ASSERTION PASS: every changed line under R/ is a roxygen comment line
-- 458cef48  (Part 2)
changed lines under R/: 50
non-roxygen changed lines: 0
ASSERTION PASS: every changed line under R/ is a roxygen comment line
```

**Post-conditions 2, 3, 6 and 7**

- **2.** `devtools::document()` was run after each part. Changed paths outside `R/` are `man/` only
  (15 `.Rd` files) plus `NEWS.md`, which post-condition 7 requires.
  Over `216f3405..458cef48`,
  `git diff --name-only | grep -v '^R/' | grep -v '^man/'` returns exactly one path: `NEWS.md`.
  This report adds one more, in its own commit.
- **3.** `NAMESPACE` is unchanged: `git diff 216f3405..HEAD -- NAMESPACE` is empty. Expected — no
  `@export` or `@importFrom` was touched.
- **6.** `git diff 216f3405..HEAD -- R/ | grep '^+' | grep -v '^+++' | grep -niE 'c2 *<=? *c1|deriv|binary default'`
  returns nothing. No added text asserts `c2 <= c1`, describes a derivation of `c2` from `c1`, or
  names `OR` as the binary default. The words "derivation" and "`OR` default" appear nowhere in the
  added lines. (`R/fpr_calibration.R:61` is unchanged context, section 3.)
  A separate check confirmed the added lines are **ASCII only**, per the repository's CRAN-hygiene
  convention.
- **7.** `NEWS.md:97-109`, one entry beginning `**Documentation: the threshold arguments.**`, stating
  that the roles and the per-estimand resolved defaults are now documented and that
  `consistency.threshold` is an effect threshold distinct from `pconsistency.threshold`.

---

## 6. Post-condition 4 -- the pre-change and post-change check finding sets

**Surface used.** The task names `devtools::check()`. This repository's `CLAUDE.md` records that
`devtools::check()` is the dev loop and **not** the certification surface, because its
`manual = FALSE` default passes `--no-manual` and so skips the PDF manual build -- which is exactly
the step that catches LaTeX-unsafe Rd content. This task writes Rd content and nothing else, so both
runs used the certification surface instead, `rcmdcheck::rcmdcheck(args = "--as-cran")`, which builds
the manual. That is a superset of what post-condition 4 asks for.

**Environment note.** A first pre-change attempt failed in `R CMD build` with "Pandoc is required to
build R Markdown vignettes but not available" -- `pandoc` is not on `PATH` in a non-interactive
`Rscript` session on this machine. Both recorded runs were made with
`RSTUDIO_PANDOC=/usr/lib/rstudio/resources/app/bin/quarto/bin/tools/x86_64` on `PATH`, so both build
and re-build all vignettes. This is an environment fact, not a package finding.

**Pre-change** (tree at `216f3405`, before any roxygen edit):

```
COUNTS: 0 errors | 1 warnings | 2 notes
```

| Kind | Finding |
|---|---|
| WARNING | `checking code files for non-ASCII characters` -- `R/fs_bias_coverage.R` |
| NOTE | `checking R code for possible problems` -- 9 "no visible binding for global variable" in `fs_plot_bias_coverage` (`b`, `cell`, `cov`, `cov1`, `cov2`, `estimator`, `obs`, `r`, `ref`) |
| NOTE | `checking HTML version of manual` -- HTML validation skipped, no `tidy` command on this machine (environmental) |

**Post-change** (tree at `458cef48`, all three commits applied):

```
COUNTS: 0 errors | 1 warnings | 2 notes
```

**Comparison.** Not merely equal in count -- **identical in text**. Both runs' ERRORS / WARNINGS /
NOTES blocks were extracted, timing annotations (`[NNs/NNs]`) stripped, and compared:

```
diff set_pre.txt set_post.txt   ->  no output
IDENTICAL: pre-change and post-change finding sets match exactly
```

Post-condition 4 requires only that the set does not grow. It did not grow and did not change. In
particular the manual built cleanly both times, so nothing added to the Rd is LaTeX-unsafe.

Neither pre-change finding is attributable to this task, and neither was touched by it. The non-ASCII
WARNING concerns `R/fs_bias_coverage.R`, which this task did not edit; see section 8, finding 3 for a
separate roxygen defect in that same file. A check confirmed independently that every line **added**
by this task is ASCII.

---

## 7. Section 3.4 and Part 2: kept or struck, and by whose decision

| Item | Outcome | Decided by |
|---|---|---|
| **Section 3.4** (the two meanings of "the screening threshold") | **Kept**, and committed alone as `0671525f` | The user's instruction opening this task said sections 3.4 and Part 2 "are separable and may be struck" and gave no strike. Read as: left in scope. The user's commit instruction named only Part 1 and Part 2 as separate commits; giving 3.4 its own commit as well is **my decision**, so that the item the user singled out as separable can be reverted alone -- which the task's stated purpose for separating them requires, and which a 3.4-inside-Part-1 commit would not allow. |
| **Part 2** (the DINA frontier caps) | **Kept**, committed alone as `458cef48`, three of its four bullets done | The user's instruction: "Commit Part 1 and Part 2 separately so either can be reverted alone" -- Part 2 in scope. The fourth bullet was stopped, not struck: it requires a code change, see section 4. That stop is required by task section 2, not a discretionary call. |

Reverting `0671525f` removes only the `threshold_config` `@return` item. Reverting `458cef48` removes
only the DINA caps work. Neither touches `cc1714e0`. Each would need a `devtools::document()` run
afterwards to bring `man/` back in step.

## 8. Findings -- wrong at HEAD, or otherwise worth recording

Recorded as findings. **No fix applied, no task attached**, per task section 6.6.

1. **`R/subgroup_search.R:25` says `hr.threshold` is on the log scale for `HR`. It is not.** The
   sentence reads "On the log scale for ratio measures (OR, HR), identity scale for difference
   measures (RD, MD)." For `OR` on the GLM path that is right -- `effect_threshold` arrives already
   logged. For `HR` it is **wrong**: on the survival path `effect_threshold` is `NULL`
   (`R/forestsearch_main.R:1829`, assigned only inside the non-survival branch), so `:3110` passes the
   **natural** `hr.threshold`, `screen_threshold` falls through to it
   (`R/subgroup_search.R:132-136`), and the comparison at `:681` is against
   `cox_result$hr`, which `fit_cox_for_subgroup()` returns as `exp(beta)`. A user reading this
   `@param` and passing `log(1.25)` for a survival run would screen at `HR >= 1.13`, not `1.25`.
   **The sentence is pre-existing and was not modified**: in the diff of `cc1714e0` it appears as
   unchanged context, and the text added around it makes no scale claim. Left as found, flagged here.
2. **Part 2's cited evidence does not exist in this repository.** The task names
   `quarto/gbsg/REPORT_dina_family_survey_2026-09-17.md` (`dba3073`). Neither the file nor any file
   matching `*dina_family_survey*` is present in the working tree or anywhere in
   `git log --all --diff-filter=A`, and `dba3073` resolves to no object
   (`git rev-parse --verify` fails). Every Part 2 claim in this report was therefore verified from
   source only -- which the task also required -- and nothing was taken from the missing report.
3. **`devtools::document()` errors on a pre-existing roxygen defect in `R/fs_bias_coverage.R:22`.**
   Every run emits:
   `@description failed to evaluate inline markdown code. Caused by error: Failed to parse the inline
   R code: `r = se_mean / sd_emp``. With `Roxygen: list(markdown = TRUE)` (DESCRIPTION:25), roxygen
   treats a code span whose content begins `r ` as **inline R to evaluate**, so
   `` `r = se_mean / sd_emp` `` is parsed as R and fails on the leading `=`. The companion
   `` `b = bias_log / sd_emp` `` on the same line is unaffected because it does not begin with `r `.
   Documentation still writes; the block's `@description` is the casualty. Pre-existing, untouched,
   unrelated to this task.
4. **Two citations in the predecessor audit are off by one to two lines.**
   `REPORT_criterion_and_defaults_audit_2026-09-18.md` cites the ratio log conversion at
   `forestsearch_main.R:1884-1885`; at the pre-edit tree the lines are `:1885-1886`. It also cites
   the `RD`/`IRD` consistency remap guard's `> 1.0` test at `:1829`; the actual line is `:1830`. The
   values in that report are correct; only those two line anchors are off. This is precisely the
   failure mode task section 1's read-from-source rule exists to prevent, and no value in this task's
   documentation was taken from that report.
5. **`dina_frontier()` is computed and discarded under the default screening mode.**
   `R/forestsearch_main.R:2790-2793` calls it unconditionally; `:2795` then takes the
   `selected_only = TRUE` branch (the default, `R/forestsearch_helpers.R:1075`) and uses
   `dina_subgroup()`'s cut instead, reading only `da$frontier$digits` from the resolved key set
   (`:2821`). The frontier table itself is dropped. This is now documented as the reason six of the
   seven keys change nothing on that path; the wasted computation is behaviour and is recorded here
   without a fix.
6. **`threshold_config` was an undocumented component of the `forestsearch()` return object.** It is
   produced at `R/forestsearch_main.R:2055` (GLM) and `:2105` (survival) and attached at `:2392`,
   `:2581` and `:3616`, but appeared nowhere in `@return`. Section 3.4 documents it; noted because
   its absence, not its wording, is what let "the screening threshold" name two numbers unremarked.

## 9. Scope

Documentation only. `R/` carries roxygen-comment changes and nothing else, asserted mechanically
(section 5). No behaviour, default, method or code change; no install; nothing pushed. Out of scope
and untouched, as task section 7 requires: Directive A's validation and the `c2 <= c1` sentence that
belongs with it; Directive B, the binary default estimand; Directive C parts 2 and 3; the
`args_call_all` sync defect; and the audit's other section 6 findings.
