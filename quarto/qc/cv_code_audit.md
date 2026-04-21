# ForestSearch Cross-Validation: Code Audit

**Scope.**  `R/forestsearch_cross_validation.R` (~1,360 lines), plus the
external dependencies it exercises: `get_dfpred()` and
`evaluate_comparison()` in `forestsearch_helpers.R`, and
`setup_parallel_SGcons()` / `resolve_cv_parallel_args()` for parallel
setup.

**Method.**  Top-to-bottom read; cross-referenced to León et al. (2024)
Section 4.1 definitions; checked invariants (fold coverage, seed
reproducibility, column carry-through, treatment-direction flips).  No
empirical testing — this is a code-only audit.  Items marked *(not
fully verifiable by inspection)* below would need a short R session to
confirm.

**Key fact.**  The two CV entry points (`forestsearch_Kfold`,
`forestsearch_tenfold`) share the same per-fold body, so most findings
apply symmetrically.  Where they diverge, I call that out explicitly.

---

## 1. What appears correct

These parts I traced through and am confident are implemented as the
paper and function docs say they should be.

### 1.1 Fold assignment (both functions)

- LOO case (`Kfolds = nrow(fs.est$df.est)`): scramble is skipped, so
  fold *k* holds out row *k* of `fs$df.est` by position.  Reproducible.
  Matches the paper's N-fold definition.
- K < N case: one `set.seed(seedit)` + `sample()` + `cut(..., breaks =
  Kfolds)`.  Deterministic given `seedit`, and every observation lands
  in exactly one fold.  Matches the standard `caret`-style fold
  construction.
- Post-loop validation (`forestsearch_Kfold` lines 284-292) explicitly
  checks that (a) every input row is represented exactly once in the
  returned `resCV`, and (b) fold sizes match the original `folds`
  assignment.  Good defensive checks.

### 1.2 Per-fold argument plumbing

- `cv_args` is built from `fs$args_call_all` and then mutated with four
  forced overrides: `parallel_args = sequential`, `details = FALSE`,
  `quiet = TRUE`, `plot.sg = FALSE`.  The sequential-inside-parallel
  pattern is correct (outer loop parallelizes across folds / sims).
- `cv_args$ps_hat <- NULL` is right — a pre-estimated propensity score
  has length `nrow(full data)`, not `nrow(training fold)`, so re-estimation
  on each fold is required.  `ps_method` and `ps_adjust_method` carry
  through via `args_call_all`, which is the intended behavior.

### 1.3 Training → test subgroup application

- `get_dfpred(df.predict = x.test, sg.harm = fs.train$sg.harm, version = 2)`
  applies the training-fold subgroup to the held-out test data.
  `get_dfpred` in `forestsearch_helpers.R` handles three label formats —
  brace-wrapped expressions (`"{er <= 0}"`), plain expressions
  (`"er <= 0"`), and bare column names (`"q3.1"`) — via the
  `evaluate_comparison()` operator-dispatch evaluator.  No `eval(parse())`
  anywhere, which is good.
- The `version = 2` argument is a no-op (maintained for backward
  compatibility, per the function's roxygen).

### 1.4 "No subgroup found" branch (the else in the foreach body)

When a fold's `fs.train$sg.harm` is `NULL`:

- `df.test$treat.recommend <- 1.0` assigns every held-out subject to
  the complement (Hc = recommend treatment).  This matches the paper's
  description: *"if FS does not identify a subgroup then Ĥ = ∅, Ĥc
  represents the ITT population, and the CV metrics are defined
  accordingly."*
- `sg1`, `sg2` set to `NA`, which makes `any_found` in
  `forestsearch_KfoldOut` correctly report the fold as "no subgroup
  found".

### 1.5 Sensitivity / PPV calculation (the dominant case)

In the `n_sgfound == 2` branch of `forestsearch_KfoldOut` (lines
781-786):

```r
tabit  <- with(df_CV, table(treat.recommend, treat.recommend.original))
sensH  <- tabit[1, 1] / sum(tabit[, 1])
sensHc <- tabit[2, 2] / sum(tabit[, 2])
ppvH   <- tabit[1, 1] / sum(tabit[1, ])
ppvHc  <- tabit[2, 2] / sum(tabit[2, ])
```

`table()` orders rows and columns by factor level, so `tabit[1,1]` is
`P(CV=0, Orig=0) × N`.  With the convention `treat.recommend = 0 ⇒ H`
and `treat.recommend = 1 ⇒ Hc`, these four lines compute exactly:

- `sensH` = P(CV classifies as H | originally H) — matches paper's
  sensCV(Ĥ).
- `sensHc` = P(CV classifies as Hc | originally Hc) — matches
  sensCV(Ĥc).
- `ppvH` = P(originally H | CV classifies as H) — matches ppvCV(Ĥ).
- `ppvHc` = P(originally Hc | CV classifies as Hc) — matches ppvCV(Ĥc).

Definitions align with León et al. (2024) Section 4.1.

### 1.6 Reproducibility across knits

- `forestsearch_Kfold`: single `set.seed(seedit)` at line 190 makes the
  scramble reproducible.  The `%dofuture%` loop uses
  `.options.future = list(seed = TRUE)`, which installs per-worker
  L'Ecuyer-CMRG streams.  Downstream random operations inside each
  `forestsearch()` call (bootstrap splits, etc.) are therefore
  reproducible.
- `forestsearch_tenfold`: `set.seed(seed + 1000 * ksim)` inside each
  sim-level iteration gives each simulation its own reproducible
  scramble.  The additive spacing of 1000 means sim 1 / sim 2 bootstrap
  streams are not independent in the strict statistical sense (they
  share a generator, just at different offsets), but at the spacing /
  generator-period relevant here this is negligible.  Defensible design
  choice, not a bug.

### 1.7 Parallel plumbing

- `setup_parallel_SGcons()` whitelists four plan types (`multisession`,
  `multicore`, `callr`, `sequential`), caps workers at `detectCores()`,
  and installs a future plan.  `old_plan <- future::plan()` +
  `on.exit(future::plan(old_plan), add = TRUE)` inside both CV entry
  points ensures the caller's plan is restored even on error.
- `resolve_cv_parallel_args()` correctly falls back to the fs call's
  `parallel_args` when the user doesn't pass its own, and caps workers
  to `detectCores() - 1` with a warning if overshot.

### 1.8 `CV_sgs` exact / one / any matching

For depth-1 `sg_analysis` (the GBSG case):

- `any_found`, `exact_match`, `one_match` all work correctly when the
  fold's subgroup labels have any recognizable string form.
- `cov1_match` is set equal to `exact_match` for depth 1, which is
  appropriate.

For depth-2 `sg_analysis`: logic is symmetric in `(sg1a, sg2a)`, and
`exact_match` correctly requires both factors to appear in the fold's
(sg1, sg2).

---

## 2. Potential issues, ranked by severity

### 2.1 MEDIUM — `grf_res` / `grf_cuts` not nullified in `cv_args`

The fs call signature exposes two caching parameters for re-use across
calls:

```
@param grf_res GRF results object (optional, for reuse).
@param grf_cuts List. Custom GRF cut points (optional).
```

These are stored in `args_call_all` and therefore carried into
`cv_args` by both CV entry points.  But unlike `ps_hat` (which is
explicitly nullified), `grf_res` is a **full-data** GRF fit whose
internal row indices point at the full dataset — reusing it on an
N-1-row training fold is incorrect.

Two sub-cases:

- **(a) User never passed `grf_res` / `grf_cuts`.**  Both are `NULL` in
  `args_call_all`, so the fold's `forestsearch()` call receives `NULL`
  and refits GRF on the training fold.  Correct.  This is the normal
  path and what's happening in the GBSG analysis.
- **(b) User passed `grf_res` to speed up a rerun of the main analysis
  by skipping GRF.**  CV would then reuse the full-data GRF fit on
  every training fold, silently biasing every fold toward the full-data
  subgroup.  Incorrect.

Recommendation for a future package patch: after the `cv_args$ps_hat <-
NULL` line in both functions, add `cv_args$grf_res <- NULL` and
`cv_args$grf_cuts <- NULL`.  Not blocking the paper reproduction
(case (a)), but a real hazard for any user who adopts the `grf_res`
optimization.

### 2.2 MEDIUM — `find_covariate_any_match` is dependent on brace format

`find_covariate_any_match()` in `forestsearch_cross_validation.R`
(lines 1139-1201) branches on whether `sg_target` starts with `"{"` or
`"!{"`:

```r
dda <- charmatch("{", sg_target, nomatch = 0)
ddb <- charmatch("!{", sg_target, nomatch = 0)

if (dda == 1)      { aa <- rep("{",  length(confs)) }
else if (ddb == 1) { aa <- rep("!{", length(confs)) }
else { return(rep(0, length(sg1))) }          # <<-- silent zero-out
```

If `sg.harm` is stored in any other format (bare expressions like
`"er <= 0"`, cut-name tokens like `"q3.1"` or `"z7.1"`), the function
returns a vector of zeros silently.  The consequence is that
`cov1_any` / `cov2_any` in the `CV_summary` output will be reported as
0% even when the fold did find a subgroup sharing the same underlying
covariate.

*(Not fully verifiable by inspection.)*  The actual format of
`fs$sg.harm` depends on which code path produced it in
`subgroup_consistency_main.R`.  The fallback comment in
`extract_fs_subgroup_definition()` — *"Fall back to sg.harm factor
names"* — suggests the typical production format is factor-token
names (`"q3.1"`), which would silently hit the zero-out branch.

**Impact on paper reproduction.**  Minimal — the paper reports
`exact_match` (70% for GBSG), not `cov1_any`.  The primary sens/ppv
metrics are unaffected.

Recommendation: once verified with an actual fs object, either
(a) normalize `sg.harm` to brace format at the point of construction,
or (b) rewrite `find_covariate_any_match` to be format-agnostic using
`evaluate_comparison`-style parsing.  Flag candidates for v0.3.0.

### 2.3 LOW — `n_sgfound == 1` branch computes the wrong quantity

`forestsearch_KfoldOut` lines 787-793, entered when `treat.recommend`
takes only one unique value across the entire `resCV`:

```r
tabit <- with(df_CV, table(treat.recommend, treat.recommend.original))
sensH <- if (nrow(tabit) > 0 && ncol(tabit) > 0)
           tabit[1, 1] / sum(tabit[, 1]) else NA
```

When does this branch fire?  Only if **every** CV observation gets the
same `treat.recommend` value.  The realistic trigger: no fold
identifies any subgroup, so every test-fold row defaults to
`treat.recommend = 1.0` via the else-branch of the foreach body.  In
that case `tabit` is 1×2 and its only row corresponds to CV = 1
(recommend), not CV = 0 (harm).  The formula `tabit[1, 1] / sum(tabit[
, 1])` computes `P(CV = 1 | Orig = 0)` — that is, one minus the true
sens_H — rather than the intended sens_H (which should be 0 in this
degenerate case).

**Impact on paper reproduction.**  None — Larry's GBSG analysis has
plenty of folds identifying the subgroup, so `n_sgfound == 2` and the
main branch runs.

Recommendation for a future patch: in the `n_sgfound == 1` branch,
detect which value is present and branch correctly, or simply return
sens_H = 0, ppv_H = NaN, sens_Hc = 1, ppv_Hc = (observed Hc
prevalence) since that's the meaningful summary when no fold finds any
subgroup.

### 2.4 LOW — `sg_analysis = NULL` input causes input-validation stop

`forestsearch_KfoldOut` line 697-701:

```r
required_elements <- c("cv_args", "sg_analysis", "sg0.name",
                       "sg1.name", "Kfolds", "resCV")
missing_elements <- required_elements[
  !sapply(required_elements, function(x) !is.null(res[[x]]))
]
if (length(missing_elements) > 0) {
  stop("The following elements are missing in 'res': ",
       paste(missing_elements, collapse = ", "))
}
```

This flags `sg_analysis = NULL` as "missing" and stops.  But the
downstream `CV_sgs` logic explicitly handles the `sg_depth == 0` case
(no subgroup in the full-data analysis), so the strict `!is.null` gate
is inconsistent with the rest of the function.  A user whose full-data
FS returned no subgroup cannot run CV diagnostics at all.

**Impact on paper reproduction.**  None — GBSG full-data FS identifies
`er <= 0`.

Recommendation: relax the validation to require presence-of-name but
allow NULL value for `sg_analysis`; the downstream code already handles
this correctly.

### 2.5 LOW — depth-1 "exact match" over-counts when a fold finds a
       deeper subgroup containing the full-data cut

`CV_sgs` depth-1 branch:

```r
exact_match <- ifelse((sg1 == sg1a | sg2 == sg1a), 1, 0)
```

If the full-data subgroup is depth 1 (`sg_analysis = "er <= 0"`) and a
fold finds a depth-2 subgroup `("er <= 0", "age <= 50")`, then sg1
matches sg1a → the fold is counted as an "exact match" even though
its subgroup is strictly smaller than the full-data subgroup.

This is a definitional choice, not a coding bug.  The paper's phrasing
("exact full analysis definition was reproduced") is ambiguous enough
that either interpretation is defensible.  Worth clarifying in the
package docs.

Suggested disambiguation: rename the current `exact_match` to
`contains_full_sg` or similar, and add a strict `same_depth_and_cuts`
metric that requires the fold's subgroup to have the same arity *and*
the same cuts as `sg_analysis`.

### 2.6 LOW — `resCV`-level `est.scale == "1/hr"` flip happens in
       `forestsearch_Kfold` but not in `forestsearch_tenfold`

`forestsearch_Kfold` lines 298-303 post-process `resCV` when
`est.scale = "1/hr"` (inverting `treat` for downstream display).
`forestsearch_tenfold` does *not* perform this flip before handing the
per-sim `resCV` to `forestsearch_KfoldOut`.

**Impact.**  The sens/ppv calculation uses `treat.recommend` /
`treat.recommend.original` (not `treat`), so `est.scale = "1/hr"` has
no effect on the numeric metrics.  The inconsistency only matters for
any downstream code that uses the `treat` column from a tenfold-era
`resCV` for analysis in the inverted scale.  Since `forestsearch_tenfold`
doesn't actually return `resCV` to the user (only summary statistics
via `sens_out` / `find_out`), practical impact is zero.

Recommendation: either (a) replicate the flip in `forestsearch_tenfold`
for symmetry, or (b) document that `est.scale = "1/hr"` post-processing
is only meaningful in `forestsearch_Kfold`.

### 2.7 NITPICK — redundant `suppressWarnings` nesting

`forestsearch_Kfold` lines 230 and 246 both call `suppressWarnings(...)`
around the same `try(do.call(forestsearch, cv_args), TRUE)` expression
(outer wraps the whole foreach, inner wraps the try).  Harmless but
not DRY.

---

## 3. Items not fully verifiable by inspection

These would need a short R session with an actual `fs` object to
confirm one way or the other:

1. **Actual format of `fs$sg.harm` for the GBSG run.**  The package
   has code paths producing `"q3.1"`-style factor tokens, brace-wrapped
   expressions, and readable bare expressions.  Whichever is the
   production output determines whether finding 2.2 is an active
   silent-zero bug or a dormant one.
   - Quick check: after running `fs <- forestsearch(...)`,
     `print(fs$sg.harm)` and `print(fs$grp.consistency$out_sg$sg.harm_label)`
     will reveal both formats side by side.

2. **Whether `get_dfpred` correctly handles factor tokens `"q3.1"` as
   produced by the current pipeline.**  The function's docs say yes
   (bare column names are treated as `df[[name]] == 1`), but this
   relies on the fold's `df.est` having a `"q3.1"` column populated
   consistently with training.  The CV code takes the raw
   `df_scrambled` (original confounders) not any post-processed
   encoding, so if `sg.harm` is in factor-token format but the test
   data doesn't have those columns, `get_dfpred` will fail with a
   missing-column error — and the CV will report every fold as
   `try-error` → no subgroup.  A quick single-fold manual call would
   confirm or refute this.

3. **Whether the GBSG CV results Larry observed (51% sens_H, 60% find
   rate) actually reflect algorithmic behavior vs. a format mismatch
   between `sg.harm` and `df_scrambled`.**  This ties back to item (2)
   above.

---

## 4. Summary verdict

The cross-validation machinery is, in its main execution path,
correctly implemented and faithful to León et al. (2024) Section 4.1:

- fold mechanics are right;
- sens/ppv metric formulas match the paper's definitions;
- the "no subgroup found" branch defaults to the paper-correct ITT
  convention;
- reproducibility is well-handled;
- the original `gbsg_analysis.qmd` discrepancies are driven by **user
  configuration**, not CV code defects, as we already established
  earlier in this review thread (LOO vs 7-fold; `use_twostage`; sim
  count).

The potential issues I identified above are all either (a) dormant in
Larry's current GBSG configuration, or (b) edge-case handling that
doesn't affect the primary CV metrics.  **No finding in this audit
invalidates the paper-reproduction plan.**

The two items most worth addressing in a future v0.3.0 package patch
are §2.1 (nullify `grf_res` / `grf_cuts` in `cv_args`) and §2.2
(make `find_covariate_any_match` format-agnostic).  Neither requires
any change to the CV-QC Quarto document we've built.
