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
  subgroup_notation = NULL,
  font_size = 12
)
```

## Arguments

- fs:

  ForestSearch results object, or a `grf_glm_result` object from
  [`grf.subg.harm.glm`](https://larry-leon.github.io/forestsearch/reference/grf.subg.harm.glm.md).
  When a GRF object is supplied, only Table 1 (treatment effect
  estimates) is produced; Table 2 (identified subgroups with
  consistency) returns `NULL`.

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
  Controls precision in both `sg10_out` (formatted via
  [`gt::fmt_number`](https://gt.rstudio.com/reference/fmt_number.html))
  and `tab_estimates` (threaded as `digits` to
  [`SG_tab_estimates_glm`](https://larry-leon.github.io/forestsearch/reference/SG_tab_estimates_glm.md)
  for arm summaries, Diff, effect-estimate CI, and bias-corrected CI).
  GLM and GRF-GLM paths only; the survival path uses its own formatting
  in
  [`SG_tab_estimates`](https://larry-leon.github.io/forestsearch/reference/SG_tab_estimates.md).

- include_search_info:

  Logical. Include search metadata table (default: TRUE).

- subgroup_notation:

  Character or `NULL`. When set to `"harm"` or `"benefit"`, overrides
  the default "Recommend"/"Questionable" labels with the Unicode labels
  from
  [`get_sg_labels`](https://larry-leon.github.io/forestsearch/reference/get_sg_labels.md)
  (e.g., \\\hat{H}\\ / \\\hat{H}^c\\ for harm). Default `NULL` (use
  "Recommend"/"Questionable").

- font_size:

  Numeric. Font size in pixels for table text (default: 12).

## Value

List with gt tables for estimates, subgroups, and optionally search
info.

## Details

**Orientation of the GLM estimates table.** For a non-survival fit with
the default `est.scale = "hr"` (and for a `grf_glm_result`), the
estimates table shows the per-arm `Rate(C)`, `Rate(T)` (binary) or
`Mean(C)`, `Mean(T)` (continuous, count), their raw difference `Diff`
(treated minus control; in percentage points for binary outcomes) and
the effect column. The effect column follows the fit's orientation: its
estimator is built with the fit's stored `adverse_outcome` or, when none
is stored, with
[`forestsearch`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)'s
rule (`TRUE` for binary and count outcomes, `FALSE` otherwise), so it is
on the scale of the fit's own estimates. The arm columns and `Diff` do
not follow it: they are the raw outcome. So for an unadjusted fit where
a lower outcome is adverse (`adverse_outcome = FALSE`, the default for
continuous outcomes) the effect and `Diff` have opposite signs: an MD is
`-Diff`, an RD is `-Diff / 100`, and the log of an OR or RR has the sign
opposite to `Diff`. Count outcomes and the rate measures (IRR, IRD) are
not reoriented (see
[`make_effect_estimator`](https://larry-leon.github.io/forestsearch/reference/make_effect_estimator.md)),
so for them `adverse_outcome` does not change the effect column.

## Examples

``` r
if (FALSE) { # \dontrun{
gbsg <- survival::gbsg
gbsg$grade3 <- as.integer(gbsg$grade == 3)
fs <- forestsearch(
  gbsg,
  id.name          = "pid",
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
