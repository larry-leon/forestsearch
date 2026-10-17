# ForestSearch Data Preparation and Feature Selection

Prepares a dataset for ForestSearch, including options for LASSO-based
dimension reduction, GRF cuts, forced cuts, and flexible cut strategies.
Returns a list with the processed data, subgroup factor names, cut
expressions, and LASSO selection results.

## Usage

``` r
get_FSdata(
  df.analysis,
  use_lasso = FALSE,
  use_grf = FALSE,
  grf_cuts = NULL,
  confounders.name,
  cont.cutoff = 4,
  conf_force = NULL,
  conf.cont_medians = NULL,
  conf.cont_medians_force = NULL,
  conf.cont_jcuts = NULL,
  dina_cuts = NULL,
  collapse_cuts = TRUE,
  collapse_cuts_args = list(),
  replace_med_grf = TRUE,
  defaultcut_names = NULL,
  cut_type = "default",
  exclude_cuts = NULL,
  outcome.name = "tte",
  event.name = "event",
  details = TRUE,
  outcome_type = "survival",
  offset.name = NULL
)
```

## Arguments

- df.analysis:

  Data frame containing the data.

- use_lasso:

  Logical. Whether to use LASSO for dimension reduction.

- use_grf:

  Logical. Whether to use GRF cuts.

- grf_cuts:

  Character vector of GRF cut expressions.

- confounders.name:

  Character vector of confounder variable names.

- cont.cutoff:

  Integer. Cutoff for continuous variable determination.

- conf_force:

  Character vector of forced cut expressions.

- conf.cont_medians:

  Character vector of continuous confounders to cut at median.

- conf.cont_medians_force:

  Character vector of additional continuous confounders to force median
  cut.

- conf.cont_jcuts:

  Named list of positive integers, one per continuous confounder for
  which J-quantile cuts are desired. For a variable `X` listed as
  `X = J`, the default cut set
  ([`cut_var()`](https://larry-leon.github.io/forestsearch/reference/cut_var.md):
  mean, median, Q1, Q3) is replaced by `J` binary cut points at the
  (k/(J+1))-th quantiles of X for `k = 1, ..., J`, defining J+1
  non-overlapping intervals \$\$\[\min(X), c_1),\\ \[c_1, c_2),\\
  \ldots,\\ \[c_J, \max(X)\].\$\$ Variables not listed here retain
  default behaviour. J-quantile cuts are unconditional (not filtered by
  LASSO), matching `defaultcut_names` semantics. Names must be in
  `confounders.name` and must not overlap with `defaultcut_names` or
  `conf.cont_medians_force`.

- dina_cuts:

  Character vector of DINA cut expressions (optional), merged into the
  candidate pool exactly like `grf_cuts`.

- collapse_cuts:

  Logical. If TRUE, collapse near-redundant continuous candidate cuts
  after the full pool is assembled and resolved to literal numerics.
  Cuts on the same variable and operator whose thresholds lie within a
  per-variable standard-error band are merged to a single
  rounded-centroid threshold, subject to a membership safety check; see
  [`collapse_redundant_cuts`](https://larry-leon.github.io/forestsearch/reference/collapse_redundant_cuts.md).
  Categorical, indicator, bare-name and equality cuts are untouched.
  Default TRUE; pass `collapse_cuts = FALSE` to recover the un-coarsened
  candidate pool used by earlier package versions.

- collapse_cuts_args:

  List of overrides for the coarsening, merged onto the defaults
  `list(c = 1.0, tol = 0.05, digits = 0L)`: `c` is the band multiplier
  (`band = c * sd(x)/sqrt(n)`), `tol` the membership safety tolerance
  (fraction of n when \< 1, absolute count when \>= 1), and `digits` the
  rounding for the representative threshold. Ignored when
  `collapse_cuts = FALSE`.

- replace_med_grf:

  Logical. If TRUE, removes median cuts that overlap with GRF cuts.

- defaultcut_names:

  Character vector of confounders to force default cuts.

- cut_type:

  Character. "default" or "median" for cut strategy.

- exclude_cuts:

  Character vector of cut expressions to exclude.

- outcome.name:

  Character. Name of outcome variable.

- event.name:

  Character. Name of event indicator variable.

- details:

  Logical. If TRUE, prints details during execution.

- outcome_type:

  Character. One of `"survival"` (default), `"binary"`, or
  `"continuous"`.

- offset.name:

  Character or `NULL`. Name of the follow-up time column for rate-based
  measures (IRR, IRD).

## Value

A list with components:

- `df`:

  Data frame with derived subgroup factor columns.

- `confs_names`:

  Character vector of factor column names added to `df`.

- `confs`:

  Named character vector of cut expressions (continuous cuts plus
  categorical indicators).

- `lassokeep`:

  Character vector of confounders retained by LASSO screening (empty if
  `use_lasso = FALSE`).

- `lassoomit`:

  Character vector of confounders dropped by LASSO screening (empty if
  `use_lasso = FALSE`).

## Examples

``` r
# \donttest{
library(survival)
df <- survival::gbsg
df$grade3 <- as.integer(df$grade == "3")
fs_data <- get_FSdata(df.analysis = df,
  confounders.name = c("age", "meno", "size", "grade3", "nodes", "pgr", "er"),
  outcome.name = "rfstime", event.name = "status",
  use_lasso = FALSE, use_grf = FALSE)
#> # of continuous/categorical characteristics 5 2 
#> Continuous characteristics: age size nodes pgr er 
#> Categorical characteristics: meno grade3 
#> Default cuts included (1st 20) age <= mean(age) age <= median(age) age <= qlow(age) age <= qhigh(age) size <= mean(size) size <= median(size) size <= qlow(size) size <= qhigh(size) nodes <= mean(nodes) nodes <= median(nodes) nodes <= qlow(nodes) nodes <= qhigh(nodes) pgr <= mean(pgr) pgr <= median(pgr) pgr <= qlow(pgr) pgr <= qhigh(pgr) er <= mean(er) er <= median(er) er <= qlow(er) er <= qhigh(er) 
#> Categorical: meno grade3 
#> 
#> Cut coarsening: 22 -> 21 candidates  (c=1, safety_tol=0.05, digits=0)
#>   age <= {53.0, 53.1} -> age <= 53  (0 subj moved)
#>   size <= {29.3} -> size <= 29  (0 subj moved)
#>   pgr <= {32.5} -> pgr <= 33  (4 subj moved)
#>   pgr <= {131.8} -> pgr <= 132  (2 subj moved)
#>   er <= {96.3} -> er <= 96  (0 subj moved) 
#> 
#> ===== CONSOLIDATED CUT EVALUATION (IMPROVED) =====
#> Evaluating 21 cut expressions once and caching...
#> Cut evaluation summary:
#>   Total cuts:  21 
#>   Valid cuts:  21 
#>   Errors:  0 
#> ✓ All 21 factors validated as 0/1
#> ===== END CONSOLIDATED CUT EVALUATION =====
#> 
#> # of candidate subgroup factors= 21 
#>  [1] "age <= 53"  "age <= 46"  "age <= 61"  "size <= 29" "size <= 25"
#>  [6] "size <= 20" "size <= 35" "nodes <= 5" "nodes <= 3" "nodes <= 1"
#> [11] "nodes <= 7" "pgr <= 110" "pgr <= 33"  "pgr <= 7"   "pgr <= 132"
#> [16] "er <= 96"   "er <= 36"   "er <= 8"    "er <= 114"  "meno"      
#> [21] "grade3"    
names(fs_data)
#> [1] "df"          "confs_names" "confs"       "lassokeep"   "lassoomit"  
# }
```
