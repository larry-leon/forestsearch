# Enhanced Subgroup Summary Tables (gt output)

Returns formatted summary tables for subgroups using the gt package,
with search metadata and customizable decimal precision.

## Usage

``` r
sg_tables(
  fs,
  which_df = "est",
  est_caption = "Training data estimates",
  potentialOutcome.name = NULL,
  hr_1a = NA,
  hr_0a = NA,
  ndecimals = 3,
  include_search_info = TRUE,
  font_size = 12
)
```

## Arguments

- fs:

  ForestSearch results object.

- which_df:

  Character. Which data frame to use ("est" or "testing").

- est_caption:

  Character. Caption for estimates table.

- potentialOutcome.name:

  Character. Name of potential outcome variable (optional).

- hr_1a:

  Character. Adjusted HR for subgroup 1 (optional).

- hr_0a:

  Character. Adjusted HR for subgroup 0 (optional).

- ndecimals:

  Integer. Number of decimals for formatted numbers (default: 3).

- include_search_info:

  Logical. Include search metadata table (default: TRUE).

- font_size:

  Numeric. Font size in pixels for table text (default: 12).

## Value

List with gt tables for estimates, subgroups, and optionally search
info.
