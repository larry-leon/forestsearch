# Create Summary Tables from forestsearch_KfoldOut Results

Formats the detailed output from
[`forestsearch_KfoldOut`](https://larry-leon.github.io/forestsearch/reference/forestsearch_KfoldOut.md)`(outall=TRUE)`
into publication-ready gt tables. This includes ITT estimates, original
subgroup estimates, and K-fold subgroup estimates.

## Usage

``` r
cv_summary_tables(
  kfold_out,
  title = "Cross-Validation Summary",
  subtitle = NULL,
  show_metrics = TRUE,
  digits = 3,
  font_size = 12,
  use_gt = TRUE
)
```

## Arguments

- kfold_out:

  List. Result from `forestsearch_KfoldOut(res, outall = TRUE)`. Must
  contain `itt_tab`, `SG_tab_original`, `SG_tab_Kfold`, and optionally
  `tab_all`.

- title:

  Character. Main title for combined table. Default: "Cross-Validation
  Summary".

- subtitle:

  Character. Subtitle for table. Default: NULL (auto-generated).

- show_metrics:

  Logical. Include agreement and finding metrics in output. Default:
  TRUE.

- digits:

  Integer. Decimal places for numeric formatting. Default: 3.

- font_size:

  Integer. Font size in pixels. Default: 12.

- use_gt:

  Logical. Return gt table if TRUE, data.frame if FALSE. Default: TRUE.

## Value

If `use_gt = TRUE`, returns a list with gt table objects:

- `combined_table`: Combined ITT and subgroup estimates

- `itt_table`: ITT estimates only

- `original_table`: Original full-data subgroup estimates

- `kfold_table`: K-fold subgroup estimates

- `metrics_table`: Agreement and finding metrics (if
  `show_metrics = TRUE`)

If `use_gt = FALSE`, returns equivalent data.frames.

## Details

The estimate tables are those of `kfold_out`, formatted without
recomputation. For GLM outcomes
[`forestsearch_KfoldOut`](https://larry-leon.github.io/forestsearch/reference/forestsearch_KfoldOut.md)
builds their effect column with the fit's orientation (the fit's stored
`adverse_outcome` or, when none is stored,
[`forestsearch`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)'s
rule: `TRUE` for binary and count outcomes, `FALSE` otherwise), while,
for a fit with the default `est.scale = "hr"`, the arm columns
(`Rate(C)` or `Mean(C)`, `Rate(T)` or `Mean(T)`) and their raw
difference `Diff` (treated minus control) are the raw outcome. So for an
unadjusted fit where a lower outcome is adverse
(`adverse_outcome = FALSE`, the default for continuous outcomes) the
effect and `Diff` have opposite signs. Count outcomes and the rate
measures (IRR, IRD) are not reoriented. See
[`sg_tables`](https://larry-leon.github.io/forestsearch/reference/sg_tables.md)
for the details.

## See also

[`cv_metrics_tables`](https://larry-leon.github.io/forestsearch/reference/cv_metrics_tables.md)
for formatting
[`forestsearch_tenfold()`](https://larry-leon.github.io/forestsearch/reference/forestsearch_tenfold.md)
results

## Examples

``` r
if (FALSE) { # \dontrun{
# Run K-fold CV
cv_results <- forestsearch_Kfold(fs.est = fs_result, Kfolds = 10)

# Get detailed output
kfold_out <- forestsearch_KfoldOut(cv_results, outall = TRUE)

# Create summary tables
cv_tables <- cv_summary_tables(kfold_out)
cv_tables$combined_table
cv_tables$metrics_table
} # }
```
