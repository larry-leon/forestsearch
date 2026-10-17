# Compare two extreme-subgroups simulation studies

Builds the per-subgroup comparison frame between two
[`run_subgroup_sims()`](https://larry-leon.github.io/forestsearch/reference/run_subgroup_sims.md)
results (or RDS payloads loaded from the vignettes – any list carrying
`sim_hrs` / `sim_ubs` / `sim_ns` matrices with identical column names).
For each input the per-design statistics replicate the vignettes'
Section 6.6.1 definitions exactly: tail probabilities are percentages
denominated over converged fits (`na.rm = TRUE`), medians are
unconditional, `N` is the across-trial mean subgroup size, and
structurally empty subgroups (all-`NA` columns) are mapped from `NaN` to
`NA`.

## Usage

``` r
compare_subgroup_sims(
  x,
  y,
  expect_designs = NULL,
  suffixes = c("_r", "_f"),
  est_thresholds = NULL,
  ub_thresholds = NULL
)
```

## Arguments

- x, y:

  The two studies to compare, in display order (the memo passes random-X
  then fixed-X).

- expect_designs:

  Optional length-2 character: required `design` labels for `x` and `y`
  respectively (e.g. `c("resample", "fixed")`); `NULL` skips the check.

- suffixes:

  Length-2 character appended to each statistic's column for `x` and
  `y`; the default `c("_r", "_f")` reproduces the memo's column names
  (`N_r`, `N_f`, `ub2_r`, ...).

- est_thresholds:

  Length-2 numeric `c(low, high)` for the `hr05`-style (`est < low`) and
  `hr1`-style (`est > high`) percentages; `NULL` resolves via `effect`
  metadata, then the HR legacy `c(0.5, 1.0)`.

- ub_thresholds:

  Length-2 numeric for the `ub2`-style and `ub3`-style (`UB >= t`)
  percentages; resolution as above (HR legacy `c(2, 3)`).

## Value

A plain `data.frame` (one row per subgroup, matching the memo's
static-mode frame type) with columns `subgroup`, then for each of `N`,
`ub2`, `ub3`, `mUB`, `hr05`, `hr1`, `mHR` the `x` and `y` columns
interleaved. The inputs' `n_sims` and `design` values are attached as
attributes `"n_sims"` and `"designs"`; when the inputs carry effect
metadata it is attached as attribute `"effect"` (informational;
attributes are dropped by most transformations).

## Details

Column names are structural and retained across outcome types (`hr05`,
`hr1`, `ub2`, `ub3`, `mHR`, `mUB`), exactly as `sim_hrs` serves generic
duty on GLM results: for
[`subgroup_glm()`](https://larry-leon.github.io/forestsearch/reference/subgroup_glm.md)
fits they hold estimate-scale statistics at the resolved thresholds.
Threshold resolution: explicit arguments win; otherwise the inputs'
`effect` metadata supplies them; otherwise the HR legacy `c(0.5, 1.0)` /
`c(2, 3)`. An `NA` threshold yields an all-`NA` column. The two inputs
must carry compatible metadata (both none, or the same effect measure) –
comparing an MD study against an HR study is a hard error, raised before
the panel-alignment guard so the categorical incompatibility is reported
first.

Two validations guard the row-wise alignment (formerly the memo's Guard
1 and Guard 2): the subgroup panels must be identical, and, when
`expect_designs` is supplied, each input's `design` label must match it
exactly.

The medians here use [`median()`](https://rdrr.io/r/stats/median.html),
matching the memo; they are equal to
[`summary.subgroup_sims()`](https://larry-leon.github.io/forestsearch/reference/summary.subgroup_sims.md)'s
type-7 `quantile(..., 0.50)` values up to floating-point formulation
(`(a + b) / 2` versus `a + 0.5 * (b - a)`), which is why the memo's
numbers and the vignettes' tables have always agreed to the printed
digit.

## See also

[`run_subgroup_sims()`](https://larry-leon.github.io/forestsearch/reference/run_subgroup_sims.md),
[`summary.subgroup_sims()`](https://larry-leon.github.io/forestsearch/reference/summary.subgroup_sims.md),
[`subgroup_glm()`](https://larry-leon.github.io/forestsearch/reference/subgroup_glm.md)

## Examples

``` r
if (FALSE) { # \dontrun{
pl_r <- readRDS("results/extreme_sims_resample_10000_payload.rds")
pl_f <- readRDS("results/extreme_sims_fixed_10000_payload.rds")
cmp <- compare_subgroup_sims(pl_r, pl_f,
                             expect_designs = c("resample", "fixed"))
} # }
```
