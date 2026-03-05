# Render Reference Simulation Table as gt

Converts a data frame of pre-computed reference simulation results
(e.g., digitized from a published LaTeX table) into a styled `gt` table.
This is useful for displaying published benchmark results alongside new
simulation output within vignettes or reports.

## Usage

``` r
render_reference_table(
  ref_df,
  title = "Reference Simulation Results",
  subtitle = NULL,
  bold_threshold = 0.05
)
```

## Arguments

- ref_df:

  `data.frame`. Must contain a `Metric` column and one column per
  analysis method, with a `Scenario` column for row grouping.

- title:

  Character. Table title.

- subtitle:

  Character. Table subtitle. Default: `NULL`.

- bold_threshold:

  Numeric. Values in `any(H)` rows exceeding this threshold are shown in
  bold. Set `NULL` to disable. Default: 0.05.

## Value

A `gt` table object.

## Examples

``` r
if (FALSE) { # \dontrun{
ref <- data.frame(
  Scenario = "M1 Null: N=700",
  Metric   = "any(H)",
  FS       = 0.02,
  FSlg     = 0.03,
  GRF      = 0.25
)
render_reference_table(ref, title = "Reference Results")
} # }
```
