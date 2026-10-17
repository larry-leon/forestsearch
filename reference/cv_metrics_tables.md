# Create Metrics Tables for Cross-Validation Results

Formats the `find_summary` and `sens_summary` outputs from
[`forestsearch_tenfold`](https://larry-leon.github.io/forestsearch/reference/forestsearch_tenfold.md)
or
[`forestsearch_Kfold`](https://larry-leon.github.io/forestsearch/reference/forestsearch_Kfold.md)
into publication-ready gt tables.

## Usage

``` r
cv_metrics_tables(
  cv_result,
  sg_definition = NULL,
  title = "Cross-Validation Metrics",
  show_percentages = TRUE,
  digits = 1,
  include_raw = FALSE,
  table_style = c("combined", "separate", "minimal"),
  use_gt = TRUE
)
```

## Arguments

- cv_result:

  List. Result from
  [`forestsearch_tenfold()`](https://larry-leon.github.io/forestsearch/reference/forestsearch_tenfold.md)
  or
  [`forestsearch_Kfold()`](https://larry-leon.github.io/forestsearch/reference/forestsearch_Kfold.md).
  Must contain `find_summary` and `sens_summary` elements.

- sg_definition:

  Character vector. Subgroup factor definitions for labeling (optional).
  If NULL, extracted from `cv_result$sg_analysis`.

- title:

  Character. Main title for combined table. Default: "Cross-Validation
  Metrics".

- show_percentages:

  Logical. Display metrics as percentages (0-100) instead of proportions
  (0-1). Default: TRUE.

- digits:

  Integer. Decimal places for formatting. Default: 1.

- include_raw:

  Logical. Include raw matrices (`sens_out`, `find_out`) in the output
  for detailed analysis. Default: FALSE.

- table_style:

  Character. One of "combined", "separate", or "minimal".

  - "combined": Single table with both agreement and finding metrics

  - "separate": Two separate gt tables

  - "minimal": Compact single-row summary

  Default: "combined".

- use_gt:

  Logical. Return gt table(s) if TRUE, data.frame(s) if FALSE. Default:
  TRUE.

## Value

Depending on `table_style`:

- "combined": A single gt table (or data.frame)

- "separate": A list with `agreement_table` and `finding_table`

- "minimal": A single-row gt table (or data.frame)

If `include_raw = TRUE`, also includes `sens_out` and `find_out`
matrices in the returned list.

## See also

[`cv_summary_tables`](https://larry-leon.github.io/forestsearch/reference/cv_summary_tables.md)
for formatting `forestsearch_KfoldOut(outall=TRUE)` results

## Examples

``` r
if (FALSE) { # \dontrun{
# After running forestsearch_tenfold
tenfold_results <- forestsearch_tenfold(
  fs.est = fs_result,
  sims = 100,
  Kfolds = 10
)

# Create combined metrics table
cv_tables <- cv_metrics_tables(tenfold_results)
cv_tables

# Create separate tables
cv_tables <- cv_metrics_tables(tenfold_results, table_style = "separate")
cv_tables$agreement_table
cv_tables$finding_table

# Minimal one-row summary
cv_metrics_tables(tenfold_results, table_style = "minimal")
} # }
```
