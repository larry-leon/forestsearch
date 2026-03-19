# Plot Kaplan-Meier Survival Difference Bands for ForestSearch Subgroups

Creates Kaplan-Meier survival difference band plots comparing the
identified ForestSearch subgroup (sg.harm) and its complement against
the ITT population. This function wraps `plotKM.band_subgroups()` from
the weightedsurv package, automatically extracting subgroup definitions
from ForestSearch results.

## Usage

``` r
plot_km_band_forestsearch(
  df,
  fs.est = NULL,
  sg_cols = NULL,
  sg_labels = NULL,
  sg_colors = NULL,
  itt_color = "azure3",
  outcome.name = "tte",
  event.name = "event",
  treat.name = "treat",
  xlabel = "Time",
  ylabel = "Survival differences",
  yseq_length = 5,
  draws_band = 1000,
  tau_add = NULL,
  by_risk = 6,
  risk_cex = 0.75,
  risk_delta = 0.035,
  risk_pad = 0.015,
  ymax_pad = 0.11,
  show_legend = TRUE,
  legend_pos = "topleft",
  legend_cex = 0.75,
  ref_subgroups = NULL,
  verbose = FALSE
)
```

## Arguments

- df:

  Data frame. The analysis dataset containing all required variables
  including subgroup indicator columns.

- fs.est:

  A forestsearch object containing the identified subgroup, or `NULL` if
  using pre-defined subgroup indicators.

- sg_cols:

  Character vector. Names of columns in `df` containing subgroup
  indicators (0/1). These columns must already exist in `df`. If `NULL`
  and `fs.est` is provided, columns will be created automatically.
  Default: `NULL`

- sg_labels:

  Character vector. Subsetting expressions for each subgroup,
  corresponding to `sg_cols`. These are passed to
  `plotKM.band_subgroups()` which evaluates them as R expressions (e.g.,
  `"age < 65"`, `"er <= 0"`, `"Qrecommend == 1"`). Must be same length
  as `sg_cols`. Default: `NULL` (auto-generated as `"colname == 1"`).

- sg_colors:

  Character vector. Colors for each subgroup curve, corresponding to
  `sg_cols`. Must be same length as `sg_cols`. Default: `NULL` (uses
  default color palette).

- itt_color:

  Character. Color for ITT population band. Default: `"azure3"`.

- outcome.name:

  Character. Name of time-to-event column. Default: `"tte"`.

- event.name:

  Character. Name of event indicator column. Default: `"event"`.

- treat.name:

  Character. Name of treatment column. Default: `"treat"`.

- xlabel:

  Character. X-axis label. Default: `"Time"`.

- ylabel:

  Character. Y-axis label. Default: `"Survival differences"`.

- yseq_length:

  Integer. Number of y-axis tick marks. Default: `5`.

- draws_band:

  Integer. Number of bootstrap draws for confidence band. Default:
  `1000`.

- tau_add:

  Numeric. Time horizon for the plot. If `NULL`, auto-calculated from
  data. Default: `NULL`.

- by_risk:

  Numeric. Interval for risk table. Default: `6`.

- risk_cex:

  Numeric. Character expansion for risk table text. Default: `0.75`.

- risk_delta:

  Numeric. Vertical spacing for risk table. Default: `0.035`.

- risk_pad:

  Numeric. Padding for risk table. Default: `0.015`.

- ymax_pad:

  Numeric. Y-axis maximum padding. Default: `0.11`.

- show_legend:

  Logical. Whether to display the legend. Default: `TRUE`.

- legend_pos:

  Character. Legend position (e.g., "topleft", "bottomright"). Default:
  `"topleft"`.

- legend_cex:

  Numeric. Character expansion for legend text. Default: `0.75`.

- ref_subgroups:

  Named list. Optional additional reference subgroups to include. Each
  element should be a list with:

  subset_expr

  :   Character. R expression to define subgroup (e.g., "age \< 65")

  label

  :   Character. Display label (optional, defaults to subset_expr)

  color

  :   Character. Color for the curve

  The function automatically creates indicator columns from the
  expressions. Default: `NULL`.

- verbose:

  Logical. Print diagnostic messages. Default: `FALSE`.

## Value

Invisibly returns a list containing:

- df:

  The modified data frame with subgroup indicators

- sg_cols:

  Character vector of subgroup column names used

- sg_labels:

  Character vector of subgroup labels used

- sg_colors:

  Character vector of colors used

- sg_harm_definition:

  The subgroup definition extracted from fs.est

- ref_subgroups:

  The reference subgroups list (if provided)

## Details

This function simplifies the workflow of creating KM survival difference
band plots for ForestSearch-identified subgroups. It can work in two
modes:

**Mode 1: With ForestSearch result (`fs.est` provided)**

1.  Extracts the subgroup definition from the ForestSearch result

2.  Creates binary indicator columns (Qrecommend, Brecommend) in `df`

3.  Generates appropriate labels from the subgroup definition

4.  Calls `plotKM.band_subgroups()` with configured parameters

**Mode 2: With pre-defined columns (`sg_cols` provided)**

1.  Uses existing indicator columns in `df`

2.  Requires `sg_labels` and `sg_colors` to match `sg_cols`

The sg.harm subgroup (Qrecommend) represents patients with questionable
treatment benefit (where `treat.recommend == 0` in ForestSearch output).
The complement (Brecommend) represents patients recommended for
treatment.

## Note

This function requires the weightedsurv package, which can be installed
from GitHub: `devtools::install_github("larry-leon/weightedsurv")`

## Subgroup Extraction

When `fs.est` is provided, the subgroup definition is extracted from:

- `fs.est$grp.consistency$out_sg$sg.harm_label` - Human-readable labels

- `fs.est$sg.harm` - Technical factor names (fallback)

- `fs.est$df.est$treat.recommend` - Subgroup membership indicator

## See also

[`forestsearch`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
for running the subgroup analysis
[`plot_sg_weighted_km`](https://larry-leon.github.io/forestsearch/reference/plot_sg_weighted_km.md)
for weighted KM plots
[`plot_sg_results`](https://larry-leon.github.io/forestsearch/reference/plot_sg_results.md)
for comprehensive subgroup visualization
