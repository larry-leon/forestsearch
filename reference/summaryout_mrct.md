# Summary Tables for MRCT Simulation Results

Creates summary tables from MRCT simulation results using the gt
package. Summarizes hazard ratio estimates, subgroup identification
rates, and classification of identified subgroups.

## Usage

``` r
summaryout_mrct(
  pop_summary = NULL,
  mrct_sims,
  sg_type = 1,
  tab_caption = "Identified subgroups and estimation summaries",
  digits = 3,
  showtable = TRUE
)
```

## Arguments

- pop_summary:

  List. Population summary from large sample approximation (optional).
  Default: NULL

- mrct_sims:

  data.table. Simulation results from
  [`mrct_region_sims`](https://github.com/larry-leon/forestsearch/reference/mrct_region_sims.md)

- sg_type:

  Integer. Type of subgroup summary:

  - 1: Basic summary (found, biomarker, age)

  - 2: Extended summary (all subgroup types)

  Default: 1

- tab_caption:

  Character. Caption for the output table. Default: "Identified
  subgroups and estimation summaries"

- digits:

  Integer. Number of decimal places for numeric summaries. Default: 3

- showtable:

  Logical. Print the table. Default: TRUE

## Value

List with components:

- res:

  List of summary statistics from population (if provided)

- out_table:

  Formatted gt table object (or data.frame if gt unavailable)

- data:

  Processed mrct_sims data.table with derived variables

- summary_df:

  Data frame of computed summary statistics

## See also

[`mrct_region_sims`](https://github.com/larry-leon/forestsearch/reference/mrct_region_sims.md)
for generating simulation results

## Examples

``` r
if (FALSE) { # \dontrun{
# After running simulations
summary_results <- summaryout_mrct(
  pop_summary = NULL,
  mrct_sims = results_alt,
  sg_type = 1,
  tab_caption = "Subgroups under strong biomarker effect"
)
} # }
```
