# Enhanced Subgroup Summary Tables (gt output)

Returns formatted summary tables for subgroups using the gt package,
with search metadata and customizable decimal precision. Produces two
tables: a treatment effect estimates table and an identified subgroups
table, each with fully customizable titles and subtitles.

## Usage

``` r
sg_tables(
  fs,
  which_df = "est",
  est_title = "Treatment Effect Estimates",
  est_caption = "Training data estimates",
  sg_title = "Identified Subgroups",
  sg_subtitle = NULL,
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

- est_title:

  Character or NULL. Main title for the estimates table (default:
  "Treatment Effect Estimates"). Rendered as bold markdown. Set to NULL
  to suppress the title and display only `est_caption`.

- est_caption:

  Character. Subtitle for the estimates table (default: "Training data
  estimates").

- sg_title:

  Character or NULL. Main title for the identified subgroups table
  (default: "Identified Subgroups"). Rendered as bold markdown. Set to
  NULL to suppress the title and display only `sg_subtitle`.

- sg_subtitle:

  Character or NULL. Subtitle for the identified subgroups table. When
  NULL (default), an informative subtitle is auto-generated from `maxk`
  (e.g., "Two-factor subgroups (maxk=2)").

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

## Examples

``` r
if (FALSE) { # \dontrun{
library(survival)
fs <- forestsearch(
  gbsg,
  confounders.name = c("age", "meno", "size", "grade3", "nodes", "pgr", "er"),
  outcome.name     = "rfstime",
  treat.name       = "hormon",
  event.name       = "status",
  fs.splits        = 50,
  use_lasso        = FALSE,
  use_grf          = FALSE
)
tabs <- sg_tables(fs)
tabs$tab_estimates
} # }
```
