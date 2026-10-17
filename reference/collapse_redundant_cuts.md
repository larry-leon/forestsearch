# Collapse near-redundant continuous candidate cuts

Continuous covariates can generate many candidate threshold cuts that
are practically redundant – e.g. "age \<= 35.7" and "age \<= 35", or
"wtkg \<= 75.2" and "wtkg \<= 74.5" – because the candidate pool unions
quantile cuts, GRF / DINA splits, and calibrated-DGM cuts at different
precisions. This helper merges such near-duplicates to a single
representative threshold, leaving categorical and indicator cuts
untouched.

## Usage

``` r
collapse_redundant_cuts(
  cuts,
  df,
  confounders.name = NULL,
  c_band = 1,
  safety_tol = 0.05,
  digits = 0L,
  cont.cutoff = 4,
  details = FALSE,
  protect = character(0)
)
```

## Arguments

- cuts:

  Character vector of candidate cut expressions, e.g. `"age <= 35.7"`,
  `"wtkg >= 87.2"`, `"gender"`.

- df:

  Data frame supplying the covariate columns, used for sd(x), the sample
  size n, and the membership safety check.

- confounders.name:

  Character vector of confounder names (currently informational;
  variable identity is parsed from each cut expression).

- c_band:

  Numeric \>= 0. Band multiplier; `band = c_band * sd(x) / sqrt(n)`.
  Larger values merge more aggressively. Default 1.

- safety_tol:

  Numeric \> 0. Membership safety tolerance: a fraction of n when \< 1,
  an absolute subject count when \>= 1. Default 0.05.

- digits:

  Integer \>= 0. Decimal places for the representative threshold.
  Default 0 (nearest integer). Raise for variables measured on a
  sub-unit scale.

- cont.cutoff:

  Integer. A variable with at least this many unique values is treated
  as continuous (see
  [`is.continuous()`](https://larry-leon.github.io/forestsearch/reference/is.continuous.md)).
  Default 4.

- details:

  Logical. If TRUE, print a per-cluster collapse report.

- protect:

  Character vector of cut expressions to exempt from coarsening (matched
  verbatim against `cuts`). Intended for screening-engine candidates
  (GRF/DINA) that should survive the collapse intact. Default
  `character(0)` (collapse everything eligible).

## Value

Character vector of de-duplicated cut expressions (length \<= input).

## Details

For each (variable, operator) group of continuous-variable cuts,
distinct thresholds are single-linkage clustered using a per-variable
band `band = c_band * sd(x) / sqrt(n)` (the standard error of the
variable: thresholds finer than ~1 SE are not statistically resolvable).
Each cluster collapses to one representative, the cluster mean rounded
half-up to `digits` places (see
[`collapse_cuts_round()`](https://larry-leon.github.io/forestsearch/reference/collapse_cuts_round.md)).
Rounding alone can also merge singletons that round to the same value
(e.g. 75.2 and 74.5 both round to 75 at `digits = 0`); de-duplication
then drops the duplicate.

Safety check: a cluster is collapsed only if replacing each member
threshold by the representative changes subgroup membership for no more
than `safety_tol` subjects (a fraction of n when `safety_tol < 1`, an
absolute count when `safety_tol >= 1`). Clusters that would move more
than that are kept unchanged, so coarsening never silently redefines a
candidate subgroup by more than the tolerance.

Cut expressions that are categorical, indicator-valued (fewer than
`cont.cutoff` unique values), bare variable names, or equality tests
(`==`) pass through unchanged. Cuts listed in `protect` are also passed
through verbatim: they are model-identified screening candidates (e.g.
GRF/DINA cuts) rather than points on a default continuous grid, so
coarsening them would defeat the purpose of the screening engine.
