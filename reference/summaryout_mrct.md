# Summary Tables for MRCT Simulation Results

Creates summary tables from MRCT simulation results using the gt
package. Summarizes hazard ratio estimates, subgroup identification
rates, and classification of identified subgroups. Optionally displays
two scenarios (e.g., alternative and null hypotheses) side by side.

## Usage

``` r
summaryout_mrct(
  pop_summary = NULL,
  mrct_sims,
  mrct_sims_null = NULL,
  scenario_labels = c("Alternative", "Null"),
  pop_summary_null = NULL,
  sg_type = 1,
  tab_caption = "Identified subgroups and estimation summaries",
  digits = 3,
  trim_threshold = 1000,
  trim_fraction = 0.01,
  table_width = 600,
  font_size = 11,
  showtable = TRUE
)
```

## Arguments

- pop_summary:

  List. Population summary from large sample approximation (optional).
  Default: NULL

- mrct_sims:

  data.table. Simulation results from
  [`mrct_region_sims`](https://larry-leon.github.io/forestsearch/reference/mrct_region_sims.md)
  (first / primary scenario).

- mrct_sims_null:

  data.table. Optional second set of simulation results (e.g., null
  hypothesis). When supplied, the table displays two value columns side
  by side. Default: NULL (single-scenario table).

- scenario_labels:

  Character vector of length 2. Column headers for the two scenarios.
  Only used when `mrct_sims_null` is supplied. Default:
  `c("Alternative", "Null")`.

- pop_summary_null:

  List. Population summary for the null scenario (optional). Default:
  NULL

- sg_type:

  Integer. Type of subgroup summary: 1 = basic summary (found,
  biomarker, age); 2 = extended summary (all subgroup types). Default: 1

- tab_caption:

  Character. Caption for the output table. Default: "Identified
  subgroups and estimation summaries"

- digits:

  Integer. Number of decimal places for numeric summaries. Default: 3

- trim_threshold:

  Numeric. When the raw mean of a metric exceeds this value in absolute
  terms, the summary switches to a symmetrically trimmed mean and SD
  (excluding the lower and upper `trim_fraction` of observations).
  Trimmed values are marked with `*` and a footnote is added to the
  table. Set to `NULL` to disable trimming entirely. Default: 1000.

- trim_fraction:

  Numeric between 0 and 0.5. Fraction of observations to trim from each
  tail when trimming is triggered. Default: 0.01 (1 percent from each
  tail, i.e., the central 98 percent of values).

- table_width:

  Numeric. Total table width in pixels. Column widths are allocated
  proportionally. Increase for HTML/wide displays (e.g., 750), decrease
  for beamer slides (e.g., 550). Default: 600.

- font_size:

  Numeric. Base font size in pixels. Title is `font_size + 2`, subtitle
  matches base. Reduce for beamer (e.g., 9 or 10). Default: 11.

- showtable:

  Logical. Print the table. Default: TRUE

## Value

List with components:

- res:

  List of summary statistics from population. When dual-scenario,
  contains `res_alt` and `res_null`.

- out_table:

  Formatted gt table object, or data.frame if gt is unavailable.

- data:

  Processed mrct_sims data.table with derived variables. When
  dual-scenario, also contains `data_null`.

- summary_df:

  Data frame of computed summary statistics.

## See also

[`mrct_region_sims`](https://larry-leon.github.io/forestsearch/reference/mrct_region_sims.md)
for generating simulation results

## Examples

``` r
if (FALSE) { # \dontrun{
# Single scenario (backward-compatible)
summaryout_mrct(
  mrct_sims = results_alt,
  tab_caption = "H1: Heterogeneous treatment effect"
)

# Dual scenario: alternative vs null side-by-side
summaryout_mrct(
  mrct_sims = results_alt,
  mrct_sims_null = results_null,
  scenario_labels = c("Alternative (HTE)", "Null (uniform)"),
  tab_caption = "Operating characteristics: Alternative vs Null"
)

# Custom trimming: 2% from each tail when mean > 500
summaryout_mrct(
  mrct_sims = results_alt,
  mrct_sims_null = results_null,
  trim_threshold = 500,
  trim_fraction = 0.02
)
} # }
```
