# `forestsearch` Codebase Review

**Repository:** `larry-leon/forestsearch`  
**Version:** 0.1.0  
**Review scope:** All 46 R source files (~31,000 lines), DESCRIPTION, NAMESPACE  
**Focus:** CRAN compliance · Readability · Efficiency · Style

---

## Executive Summary

The package is impressively engineered — the methodology is sophisticated, the
roxygen2 documentation is thorough, and the overall architecture is well-structured.
The primary path to a clean CRAN submission involves addressing two groups of issues:
(1) **must-fix** items that will cause `R CMD check --as-cran` errors or notes, and
(2) **should-fix** items that are style, safety, or API-quality improvements.

---

## 1. CRAN — Must Fix

### 1.1 Hyphenated source filename

`R/forestsearch_cross-validation.R` contains a hyphen. CRAN strongly prefers
alphanumeric-plus-underscore filenames for R source files.

**Fix:** Rename to `forestsearch_cross_validation.R`.

```bash
git mv R/forestsearch_cross-validation.R R/forestsearch_cross_validation.R
```

---

### 1.2 `par()` changed without `on.exit` restore — `cox_spline_fit.R`

`cox_ahr_cde_wrapper.R` correctly saves and restores graphics parameters:
```r
oldpar <- par(no.readonly = TRUE)
on.exit(par(oldpar), add = TRUE)
```

But `cox_spline_fit.R` (lines 697–758) modifies `par(mfrow)` without saving first:
```r
par(mfrow = c(1, 2))   # line 697 – no save
...
par(mfrow = c(1, 1))   # line 758 – manual "reset" only
```

Resetting to `(1, 1)` is not the same as restoring the user's prior settings. This
is a **CRAN policy violation** (`par()` must be restored via `on.exit`).

**Fix:** At the start of the relevant function:
```r
oldpar <- par(no.readonly = TRUE)
on.exit(par(oldpar), add = TRUE)
```
Remove the manual `par(mfrow = c(1, 1))` reset at the end.

---

### 1.3 `options(warn = -1)` double-restore in `bootstrap_dofuture_main.R`

`forestsearch_bootstrap_dofuture()` (lines 185–270) has both an `on.exit` handler
*and* an explicit manual restore of `options(warn)`:

```r
old_warn <- getOption("warn")
options(warn = -1)         # line 187

on.exit({
  options(warn = old_warn) # line 191 — fires on function exit
  ...
}, add = TRUE)
...
options(warn = old_warn)   # line 270 — also called explicitly mid-function
```

The explicit call on line 270 is redundant with `on.exit` and leaves warnings globally
suppressed for any code *between* lines 189 and 270 (the entire bootstrap section).
CRAN dislikes blanket warning suppression.

**Fix:** Remove line 270. Let `on.exit` handle the restore exclusively. Prefer
`suppressWarnings()` on specific calls that are known to warn, rather than
`options(warn = -1)` across the whole function.

---

### 1.4 `seq_len()` / `seq_along()` instead of `1:n` patterns

Several files use `1:nrow()`, `1:length()`, `1:ncol()` which fail when `n = 0`.
R CMD check `--as-cran` now warns about these. Affected locations include:

| File | Lines |
|------|-------|
| `generate_aft_dgm_helpers.R` | 287, 688, 1355, 1454 |
| `get_fsdata.R` | 259, 328 |
| `plot_sg_results.R` | 820, 843 |
| `subgroup_search.R` | 371 |
| `summary_utility_functions.R` | 510 |

**Fix:** Replace `1:nrow(x)` → `seq_len(nrow(x))`, `1:length(x)` → `seq_along(x)`.

---

### 1.5 Remaining `eval(parse())` — `forestsearch_helpers.R` line 247

Although `safe_eval_expr()` is well-designed (restricted `baseenv()` parent), CRAN
reviewers flag all `eval(parse())` patterns. The function handles compound expressions
like `"BM > 1 & tmrsize > 19"` that `evaluate_comparison()` doesn't cover because
it only handles single comparisons.

**Options:**
- If expressions are always binary comparisons joined by `&`/`|`, split on those
  operators and evaluate each side with `evaluate_comparison()`, then combine.
- Otherwise, add a comment with `# nocov` and a prominent `@note` in the docs
  explaining why `eval(parse())` is unavoidable here, and be prepared to defend
  this to CRAN in `cran-comments.md`.

---

### 1.6 56 exported functions with no `\examples{}` section

Of the 113 exported functions, 56 lack any `\examples{}` entry in their `.Rd` files.
CRAN policy requires examples for all exported functions (or explicit `\dontrun{}`
wrapping when examples require unavailable resources).

**Recommendation:** For functions that are expensive to run, use `\donttest{}` (runs
in `R CMD check --run-donttest` but not on CRAN's servers by default) rather than
`\dontrun{}` (never runs). For truly internal helpers that are exported only for
technical reasons, consider moving them to `@keywords internal` with no example.

Functions currently missing examples include many helpers like `build_cox_formula()`,
`calculate_hazard_ratios()`, `compute_ahr()`, `format_CI()`, etc. Even a two-line
minimal example with synthetic data covers the requirement.

---

### 1.7 `@import survival` and `@import data.table` declared in multiple files

These `@import` directives appear in:
- `globals.R` (canonical location ✓)
- `find_k_inter_main.R` (duplicate)
- `oc_analyses_gbsg.R` (duplicate)

Duplicate `@import` declarations produce duplicate entries in `NAMESPACE` after
`roxygenise()`, which generates a NOTE.

**Fix:** Remove the redundant declarations from `find_k_inter_main.R` and
`oc_analyses_gbsg.R`. Keep only the ones in `globals.R`.

---

## 2. CRAN — Likely Notes

### 2.1 `<<-` in `get_FSdata_helpers.R` and `mrct_simulation.R`

`<<-` is used in two `tryCatch` error handlers (get_FSdata_helpers.R:384–385) and one
inner closure (mrct_simulation.R:1431). While technically correct, CRAN sometimes
queries `<<-` in non-obvious contexts.

**`get_FSdata_helpers.R` (lines 384–385):**  
The `<<-` modifies `has_error[i]` and `is_valid[i]` defined in the enclosing loop.
This is a legitimate closure pattern, but it can be replaced with an environment:

```r
# Preferred: use an environment or preallocate and index
results_env <- new.env(parent = emptyenv())
results_env$has_error <- logical(n_confs)
results_env$is_valid  <- logical(n_confs)
# Inside tryCatch error handler:
results_env$has_error[i] <- TRUE
results_env$is_valid[i]  <- FALSE
```

Or simply restructure so `tryCatch` returns a value and assigns normally:

```r
result_i <- tryCatch({
  r <- evaluate_comparison(thiscut, df)
  list(eval = as.logical(r), valid = length(unique(r)) > 1, error = FALSE)
}, error = function(e) {
  list(eval = NULL, valid = FALSE, error = TRUE, msg = e$message)
})
has_error[i] <- result_i$error
is_valid[i]  <- result_i$valid
```

---

### 2.2 `DESCRIPTION` missing `Date:` field

While not strictly required, CRAN policy recommends including a `Date:` field.

```
Date: 2025-01-01
```

---

## 3. API Design — Should Fix

### 3.1 Naming convention inconsistency: dot-style vs snake_case

Four exported functions still use dot-style names, creating an inconsistent public API:

| Dot-style name | Suggested snake_case alias |
|---|---|
| `subgroup.search()` | `subgroup_search()` |
| `subgroup.consistency()` | `subgroup_consistency()` |
| `grf.subg.harm.survival()` | `grf_subgroup_harm_survival()` |
| `grf.subg.eval()` | `grf_subgroup_eval()` |

**Recommended approach:** Add snake_case aliases that call the dot-style functions,
mark the dot-style names as deprecated with `.Deprecated()`, and document in
`NEWS.md`. This maintains backward compatibility while presenting a clean API.

```r
#' @rdname subgroup_search
#' @export
subgroup_search <- function(...) {
  .Deprecated("subgroup_search", old = "subgroup.search")
  subgroup.search(...)
}
```

---

### 3.2 Parameter naming inconsistency in `forestsearch()`

The main function mixes dot-style and snake_case parameter names within the
same function signature:

**Dot-style parameters:**
`is.RCT`, `est.scale`, `n.min`, `hr.threshold`, `hr.consistency`, `plot.sg`,
`plot.grf`, `by.risk`, `m1.threshold`, `vi.grf.min`, `dmin.grf`, `frac.tau`

**Snake_case parameters:**
`use_lasso`, `use_grf`, `grf_res`, `grf_cuts`, `parallel_args`, `use_twostage`,
`twostage_args`, `sg_focus`, `stop_threshold`, `max_subgroups_search`

This mixed convention is the single most noticeable style issue in the package.
Since this is a breaking change, it's best addressed now before CRAN submission.

**Suggested migration path:**
1. Rename all parameters to `snake_case` in `forestsearch()` and downstream functions.
2. During a transition period, accept both via `...` matching and warn on old names.

---

### 3.3 Hardcoded bootstrap seed in `forestsearch_bootstrap_dofuture()`

```r
boot_seed <- 8316951L  # hardcoded, not user-settable
```

The bootstrap seed is not exposed as a parameter, making exact reproducibility
impossible for users who want to vary the seed. This is unusual for a statistical
method function.

**Fix:** Add a `seed` parameter:
```r
forestsearch_bootstrap_dofuture <- function(fs.est,
                                            nb_boots,
                                            seed = 8316951L,
                                            details = FALSE,
                                            show_three = FALSE,
                                            parallel_args = list()) {
```

---

## 4. Readability / Efficiency

### 4.1 Very large files warrant splitting

| File | Lines | Functions |
|------|-------|-----------|
| `oc_analyses_gbsg.R` | 2,061 | 15 |
| `mrct_simulation.R` | 1,847 | 6 |
| `generate_aft_dgm_helpers.R` | 1,584 | 19 |
| `sim_aft_gbsg.R` | 1,280 | 8 |

`generate_aft_dgm_helpers.R` with 19 functions is the clearest candidate for
splitting: the DGM construction helpers, censoring calibration functions, and
comparison/output functions are logically distinct groups.

---

### 4.2 `acm.disjctif()` naming and visibility

`acm.disjctif()` is an internal helper (has `@keywords internal`, no `@export`)
named after an old SAS/FactoMineR convention. It has a man page (`acm.disjctif.Rd`)
generated from its roxygen block but is not in `NAMESPACE`. This is a minor anomaly
but may confuse users browsing the documentation.

**Recommendation:** Either:
- Add `@noRd` to suppress the `.Rd` file entirely (preferred for pure internals), or
- Rename to `dummy_encode_factors()` consistent with the public `dummy_encode()` wrapper.

---

### 4.3 `bootstrap_dofuture_main.R` global warning suppression scope

Related to issue 1.3: the `options(warn = -1)` on line 187 suppresses warnings for
the **entire** function body, including input validation and the bootstrap setup.
If a legitimate warning fires during setup (e.g., a bad `parallel_args` value), it
will be silently swallowed.

**Better pattern:** Suppress warnings only at the parallel dispatch point:
```r
results <- suppressWarnings(
  bootstrap_results(...)
)
```

---

### 4.4 Consolidate `simulate_from_dgm()` and `simulate_from_gbsg_dgm()`

Per known roadmap: these two simulation engines share substantial logic and the
consolidation is already planned. From a CRAN perspective, two very similar
top-level simulation functions with partially overlapping signatures increases
documentation burden and user confusion. Completing this consolidation before
submission is worthwhile.

---

## 5. Positive Observations

- **`globals.R`** is exemplary: a single, well-organised file with section headers
  covering all NSE variables. This is the right approach.
- **S3 methods** are properly registered in `NAMESPACE` via `S3method()`.
- **`stop_threshold` validation** is thorough — the guard against values > 1.0 with
  a clear error message is exactly right.
- **`safe_eval_expr()`** using `list2env(..., parent = baseenv())` is a sound,
  security-conscious design for expression evaluation.
- **`on.exit(par(oldpar), add = TRUE)`** in `cox_ahr_cde_wrapper.R` — the `add = TRUE`
  is critical and correctly used.
- **`do.call()` with merged param lists** for flexible parameter passing is a clean
  and CRAN-safe pattern.
- **`.Rbuildignore`** is comprehensive — `dev/`, `docs/`, `_pkgdown.yml` all correctly
  excluded.
- **`DESCRIPTION`** has a proper DOI reference and dual URL entries.

---

## Summary Checklist

| # | Issue | Severity | File(s) |
|---|-------|----------|---------|
| 1.1 | Hyphenated source filename | ⚠️ CRAN NOTE | `forestsearch_cross-validation.R` |
| 1.2 | `par()` not restored via `on.exit` | ⚠️ CRAN POLICY | `cox_spline_fit.R` |
| 1.3 | `options(warn=-1)` double-restore | ⚠️ CRAN POLICY | `bootstrap_dofuture_main.R` |
| 1.4 | `1:nrow()` / `1:length()` patterns | ⚠️ CRAN NOTE | multiple |
| 1.5 | `eval(parse())` remaining instance | ⚠️ CRAN QUERY | `forestsearch_helpers.R` |
| 1.6 | 56 exported fns missing `\examples` | ⚠️ CRAN NOTE | multiple `.Rd` |
| 1.7 | Duplicate `@import` declarations | ⚠️ CRAN NOTE | `find_k_inter_main.R`, `oc_analyses_gbsg.R` |
| 2.1 | `<<-` in tryCatch handlers | ℹ️ Style | `get_FSdata_helpers.R`, `mrct_simulation.R` |
| 2.2 | Missing `Date:` in DESCRIPTION | ℹ️ Minor | `DESCRIPTION` |
| 3.1 | Dot-style exported function names | ℹ️ API | `grf_main.R`, `subgroup_*.R` |
| 3.2 | Mixed param naming in `forestsearch()` | ℹ️ API | `forestsearch_main.R` |
| 3.3 | Hardcoded bootstrap seed | ℹ️ Reproducibility | `bootstrap_dofuture_main.R` |
| 4.1 | Very large files (1500–2000+ lines) | ℹ️ Readability | `oc_analyses_gbsg.R`, `generate_aft_dgm_helpers.R` |
| 4.2 | `acm.disjctif` naming / visibility | ℹ️ Minor | `get_FSdata_helpers.R` |
