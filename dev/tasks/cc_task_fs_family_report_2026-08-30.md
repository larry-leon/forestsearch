# CC task — `fs_family_report()`: what is and is not deterministic about the candidate family

**File:** `dev/tasks/cc_task_fs_family_report_2026-08-30.md` · **Issued:** 2026-08-29 by chat, commissioned by Larry
**Repo:** `~/Documents/GitHub/forestsearch` (pop-os) · **Branch:** `feature/glm-extension`
**Transport:** Larry places this file in `~/Downloads`; CC copies it into `dev/tasks/` and commits it first (§1).
**Run order:** after `cc_task_oc_wrapper_confs_2026-08-30.md` **and** after `cc_task_vi_grf_default_2026-08-30.md`. Never concurrently with either.

**Why.** A user intending a bona-fide fixed-family analysis can set `use_lasso = FALSE`, `use_grf = FALSE`, `vi.grf.min = NULL` and reasonably believe the candidate family is now fixed. It is not: cuts land at *sample* quantiles, `n.min` / `minp` / `rmin` are floors on *sample* counts, and near-duplicate removal is statistics-keyed with no off switch. Those are the method, not screening, and no argument disables them. This function says so, for a given set of arguments, in a form that cannot overpromise. It is a reporter, not a setter.

---

## ⚠ `R/` CALLOUT — **ADDITIVE ONLY**

**One new file, one new export, one new `print()` method. No existing `R/` file may be edited.** The function computes nothing that enters any result and is called by nothing else in the package. If you believe an existing file must change to make this work, **stop and report** — that is a proposal for Larry.

Files that may be created or regenerated: `R/fs_family_report.R`, `tests/testthat/test-fs-family-report.R`, `man/`, `NAMESPACE`, `DESCRIPTION` (patch bump from whatever HEAD carries), `NEWS.md`.

**Compute:** trivial. §2 is reading; §5's tests build one small fixture. No simulations, no replicates, no renders.

**Unattended.** Larry is away. Every gate is stop-on-failure: stop, write what you have into the report with the failure at the top, commit, end. **Never ask. Never work around.**

---

## 1. Provenance and first commit — GATE

```bash
hostname; pwd
git rev-parse --abbrev-ref HEAD; git rev-parse --short HEAD
git status --porcelain
git log --oneline -6
Rscript -e 'cat("forestsearch", as.character(packageVersion("forestsearch")), "\n")'
```

*GATE:* branch `feature/glm-extension`, clean tree, and the two prior tasks' commits present. Copy this document into `dev/tasks/` and commit it alone (the `~/Downloads` stem arrives with hyphens stripped — copy under the header name, report both). Run `devtools::test()` for the parity baseline.

---

## 2. Verify the stage inventory from source — GATE

The table below is **chat's claim, not established fact**. Verify every row from source and quote the governing line. Correct anything wrong; add any stage missing. The function's value is entirely in this table being right.

| stage | governed by | chat's claim |
|---|---|---|
| cut construction | `conf.cont_jcuts`, `cut_type`, `cont.cutoff`, `conf.cont_medians`, `conf.cont_medians_force`, `conf_force`, `defaultcut_names`, `exclude_cuts`, `collapse_cuts`, `collapse_cuts_args` | data-dependent — cuts at sample quantiles; fixable only by supplying cut *values* |
| LASSO screen | `use_lasso` | disableable |
| DINA / GRF subgroup paths | `use_dina`, `use_grf`, `subgroup_method` | disableable |
| GRF VI ordering | `vi.grf.min` | ordering only at values ≤ 0 (`variable_importance()` non-negative, block guarded by `vi_max > 0`); `NULL` skips Section 5 entirely |
| confounder cap | `max_n_confounders` | applied **only** inside the VI block, so inert when `vi.grf.min` is `NULL` |
| dummy expansion | — | deterministic |
| combination enumeration | `maxk` | deterministic given the columns |
| per-factor prevalence floor | `minp` | data-dependent — `colMeans(x) >= minp` on the sample |
| redundancy | `rmin` | data-dependent — subject counts, and sensitive to column order |
| subgroup size floor | `n.min`, `n.min.frac` | data-dependent — `sum(id.x) > n.min` on the sample |
| arm-count floors | `d0.min`, `d1.min` | skipped entirely for continuous outcomes |
| effect screen | `effect.threshold` / `hr.threshold` | data-dependent — per-candidate fit |
| near-duplicate removal | `sg_focus` | data-dependent, **no off switch**; `maxeff` gets exact-membership dedup, every other focus including `maxeffCons` gets the statistics-keyed one |
| consistency screen | `consistency_method`, `pconsistency.threshold`, `consistency.threshold`, `fs.splits`, `use_twostage`, `twostage_args` | data-dependent — the gate |
| early stopping | `stop_threshold`, `sg_focus` | meaningful for `maxeffCons` only, reset to `NULL` for other foci; truncates what is evaluated but the prefix winner is the global winner, so the answer is unchanged |
| winner selection | `sg_focus`, `selection_rule`, `effect_neighborhood` | deterministic given the qualifying set |
| time cap | `max.minutes` | **inert** — no path consults it |

*GATE:* if any row is wrong in a way that changes its status classification, correct it and **say so prominently** — chat's design rests on this table.

---

## 3. The function

`R/fs_family_report.R`, one export.

```r
fs_family_report(x, data = NULL, outcome_type = NULL)
```

- **`x`** — either a named list of `forestsearch()` arguments, or a fitted `forestsearch()` object, in which case read its `args_call_all`. Detect which; error clearly on anything else.
- **`data`** — optional. When supplied, ground the report in counts: the number of cut columns `get_FSdata()` produces and the number of `maxk` combinations that implies. Never fit anything; never run the search.
- **`outcome_type`** — needed because the `d0.min`/`d1.min` row differs; take it from `x` when present.

**It must mirror `forestsearch()`'s own argument resolution**, or it will misreport. At minimum: `stop_threshold` reset to `NULL` for foci other than `maxeffCons` (`forestsearch_main.R` ~L1575–1610); `max_n_confounders` inert when `vi.grf.min` is `NULL`; `d0.min`/`d1.min` inapplicable for continuous. §4 test 2 guards this against drift.

**Return** a classed `data.frame` (`fs_family_report`) with one row per stage:

| column | contents |
|---|---|
| `stage` | short name |
| `arguments` | the governing argument names |
| `values` | their resolved values, formatted |
| `status` | one of `deterministic`, `disabled`, `inert`, `data-dependent`, `data-dependent (not disableable)` |
| `note` | one line: why, and what would change it — `NA` where nothing would |

Attributes: the overall verdict, counts by status, and whether `data` was supplied.

**`print()`** — a verdict line first (e.g. *"Candidate family is data-dependent: 6 of 14 stages vary with the sample; 3 are disableable, 3 are intrinsic to the method."*), then the table, then a footer listing the intrinsic ones with their one-line notes. The footer is the point of the function: it is what a user setting `use_grf = FALSE` and expecting a fixed family needs to read.

Roxygen must state plainly that this reports and does not change anything, and that no combination of arguments makes the family deterministic while cuts are placed at sample quantiles.

---

## 4. Tests — `tests/testthat/test-fs-family-report.R`

1. **Shape and contract** — every row has a status from the allowed set; `arguments` names are all real `forestsearch()` formals; the object prints without error; `data = NULL` and `data` supplied both work.
2. **Resolution agreement — the drift guard.** Build a small fixture, call `forestsearch()` once with `sg_focus` set to something other than `maxeffCons` and `stop_threshold` left at its default, and assert the report's resolved `stop_threshold` matches the fitted object's `args_call_all$stop_threshold`. Same for `max_n_confounders` inertness under `vi.grf.min = NULL`. If the report and the engine ever disagree, this fails.
3. **Coverage guard.** `setdiff(names(formals(forestsearch)), <arguments the report classifies>)` must be a documented, explicitly-listed set of out-of-scope formals (data names, parallel args, plotting, MR arguments). A new formal that is neither classified nor on that list fails the test. Keep the out-of-scope list in the test, not in the function, so adding a formal forces a decision.
4. **Status sensitivity** — `vi.grf.min = NULL` reports the VI row `disabled` and the cap row `inert`; a numeric value reports the VI row active and the cap row applicable. `use_lasso = TRUE` flips its row. The dedup and cut-construction rows report `data-dependent` under **every** argument combination tried — assert that explicitly, since it is the claim the function exists to make.
5. **Fitted-object path** — passing a fitted `forestsearch()` object gives the same table as passing its `args_call_all`.

Per the v5 §9 principle, show at least one test failing against an injected defect (e.g. hard-code the dedup row to `disabled`) and say which you checked.

---

## 5. Close-out

- `devtools::document()`; `devtools::test()` — raw counts and **warning-count parity** against §1.
- `R CMD check` (`RSTUDIO_PANDOC` on Rscript's PATH). Status verbatim; 0/0/0 is the target. Examples must be runnable — a `@examples` block that a user can execute is part of the deliverable.
- Patch version bump from whatever HEAD carries; `NEWS.md` entry stating it reports and changes nothing.
- Commit, **explicit paths only**, never `git add -A`. **No push.** `devtools::install()`; confirm the version.
- Report `dev/glm-continuous-sims/REPORT_fs_family_report_2026-08-30.md`: §2's verified table with every correction to chat's claims called out prominently; the function's signature and return; the tests including the injected-defect check; test and check output raw; a worked `print()` output on the drivers' own argument set — which is the most useful single artifact in the report; commits; ten-line verdict.

---

## 6. Out of scope

- No edit to any existing `R/` file. No `screening = "none"` setter — that is a separate, later task.
- No change to any default, to any driver, to any document, or to the OC-wrapper files.
- No family-passing feature in `forestsearch()`.
- No simulations, no replicate runs, no renders, no push.
