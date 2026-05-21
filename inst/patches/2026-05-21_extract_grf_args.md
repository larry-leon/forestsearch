# PR: `extract_grf_args()` and internal-builder refactor

This patch introduces a public function `extract_grf_args()` that
recovers the exact arguments `forestsearch()` passed to its inner GRF
subgroup-identification call (`grf.subg.harm.glm()` for binary,
continuous, and count outcomes; `grf.subg.harm.survival()` for survival
outcomes).  It also refactors the `forestsearch()` call site to use the
same internal argument-list builders the extractor reads from, so the
contract cannot drift between the two locations.

## Motivation

Before this patch, `forestsearch_main.R` lines 1303–1346 inlined the
argument lists for the inner GRF calls.  Any reproducibility audit,
sensitivity analysis, or independent standalone GRF run required the
caller to hand-mirror those 14+ arguments — and silently broke whenever
the inline list shifted.  After this patch, both the call site and the
public extractor share a single source of truth.

## Files in this patch

| File | Purpose | Action |
|---|---|---|
| `R/grf_args.R` | **New file**: `.build_grf_glm_args()`, `.build_grf_survival_args()` (internal), `extract_grf_args()` (exported) | Add |
| `R/forestsearch_main.R` | **Revised file**: refactors lines 1300–1346 to call the internal builders.  Drop-in replacement for the existing file. | Overwrite |

## Apply the refactor

```sh
cd path/to/forestsearch
cp /mnt/user-data/outputs/R/grf_args.R          R/grf_args.R
cp /mnt/user-data/outputs/R/forestsearch_main.R R/forestsearch_main.R
```

Then in R:

```r
devtools::document()       # adds export(extract_grf_args) to NAMESPACE
devtools::load_all()
```

## One-shot round-trip sanity check

Run once after applying the patch.  Three `TRUE` results means the
refactor is behavior-preserving and the new function works:

```r
set.seed(20260520)
n  <- 300L
df <- data.frame(
  id    = seq_len(n),
  treat = rbinom(n, 1, 0.5),
  x1    = rnorm(n),
  x2    = rnorm(n),
  x3    = sample(c(0L, 1L), n, replace = TRUE)
)
df$y <- 0.2 * df$treat + 0.3 * df$x1 * df$treat + 0.5 * df$x2 +
        rnorm(n, sd = 0.5)

fs <- forestsearch(
  df.analysis      = df,
  confounders.name = c("x1", "x2", "x3"),
  outcome.name     = "y",
  treat.name       = "treat",
  id.name          = "id",
  outcome_type     = "continuous",
  effect_measure   = "MD",
  hr.threshold     = 0.10,
  hr.consistency   = 0.05,
  use_grf          = TRUE,
  use_lasso        = FALSE,
  is.RCT           = TRUE,
  seedit           = 1L,
  fs.splits        = 50L,
  n.min            = 50L,
  quiet            = TRUE
)

args <- extract_grf_args(fs)
fn   <- match.fun(attr(args, "grf_function"))
res  <- do.call(fn, args)

identical(res$tree.cuts,  fs$grf_res$tree.cuts)    # TRUE
identical(res$sg.harm.id, fs$grf_res$sg.harm.id)   # TRUE
identical(res$grf_varimp, fs$grf_res$grf_varimp)   # TRUE


```

Then optionally `devtools::check()` to confirm the package builds
cleanly.

## Public API change

One symbol added to the NAMESPACE: `extract_grf_args`.  No existing
function signatures change.  The internal call-site refactor is
behavior-preserving — the new builders construct the same lists the
old inline code constructed.

## Usage

```r
# Reproduce the inner GRF run that forestsearch() performed
args <- extract_grf_args(fs)
fn   <- match.fun(attr(args, "grf_function"))
grf_standalone <- do.call(fn, args)

# Sensitivity analysis: deeper trees, same data and seed
args_deep <- modifyList(args, list(maxdepth = 3L))
grf_deep  <- do.call(fn, args_deep)

# Pull standalone variable importance
grf_standalone$grf_varimp
```
