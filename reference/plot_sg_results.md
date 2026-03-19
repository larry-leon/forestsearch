# Plot ForestSearch Subgroup Results

Creates comprehensive visualizations of subgroup results from
ForestSearch, including Kaplan-Meier survival curves, hazard ratio
comparisons, and summary statistics. This function is designed to work
with the output from
[`forestsearch`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md),
specifically the `df.est` component.

## Usage

``` r
plot_sg_results(
  fs.est,
  outcome.name = "Y",
  event.name = "Event",
  treat.name = "Treat",
  plot_type = c("combined", "km", "forest", "summary"),
  by.risk = NULL,
  conf.level = 0.95,
  est.scale = c("hr", "1/hr"),
  sg0_name = "Questionable (H)",
  sg1_name = "Recommend (H^c)",
  treat_labels = c(`0` = "Control", `1` = "Treatment"),
  colors = NULL,
  title = NULL,
  show_events = TRUE,
  show_ci = TRUE,
  show_logrank = TRUE,
  show_hr = TRUE,
  verbose = FALSE,
  ...
)
```

## Arguments

- fs.est:

  A forestsearch object or list containing at minimum:

  df.est

  :   Data frame with analysis data including `treat.recommend`

  sg.harm

  :   Character vector of subgroup-defining variable names

  grp.consistency

  :   Optional. Consistency results from `sg_consistency_out`

- outcome.name:

  Character. Name of time-to-event outcome column. Default: "Y"

- event.name:

  Character. Name of event indicator column (1=event, 0=censored).
  Default: "Event"

- treat.name:

  Character. Name of treatment column (1=treatment, 0=control). Default:
  "Treat"

- plot_type:

  Character. Type of plot to create. One of:

  "km"

  :   Kaplan-Meier survival curves

  "forest"

  :   Forest plot of hazard ratios

  "summary"

  :   Summary statistics panel

  "combined"

  :   All plots combined (default)

- by.risk:

  Numeric. Risk interval for KM survival curves. Default: NULL
  (auto-calculated)

- conf.level:

  Numeric. Confidence level for intervals. Default: 0.95

- est.scale:

  Character. Effect scale: "hr" (hazard ratio) or "1/hr" (inverse).
  Default: "hr"

- sg0_name:

  Character. Label for subgroup 0 (harm/questionable). Default:
  "Questionable (H)"

- sg1_name:

  Character. Label for subgroup 1 (recommend/complement). Default:
  "Recommend (H^c)"

- treat_labels:

  Named character vector. Labels for treatment arms. Default: c("0" =
  "Control", "1" = "Treatment")

- colors:

  Named character vector. Colors for plot elements. Default: uses
  package defaults

- title:

  Character. Main plot title. Default: auto-generated

- show_events:

  Logical. Show event counts on KM curves. Default: TRUE

- show_ci:

  Logical. Show confidence intervals. Default: TRUE

- show_logrank:

  Logical. Show log-rank p-value. Default: TRUE

- show_hr:

  Logical. Show hazard ratio annotation. Default: TRUE

- verbose:

  Logical. Print diagnostic messages. Default: FALSE

- ...:

  Additional arguments passed to plotting functions.

## Value

An object of class `fs_sg_plot` containing:

- plots:

  List of ggplot2 or base R plot objects

- summary:

  Data frame of subgroup summary statistics

- hr_estimates:

  Data frame of hazard ratio estimates

- call:

  The matched call

## Details

The function extracts subgroup membership from
`fs.est$df.est$treat.recommend`:

- `treat.recommend == 0`: Harm/questionable subgroup (H)

- `treat.recommend == 1`: Recommend/complement subgroup (H^c)

For `est.scale = "1/hr"`, treatment labels and subgroup interpretation
are reversed to maintain clinical interpretability.

## Kaplan-Meier Plots

When `plot_type = "km"`, creates side-by-side survival curves for:

1.  The identified subgroup (H) with treatment vs control

2.  The complement subgroup (H^c) with treatment vs control

## Forest Plot

When `plot_type = "forest"`, creates a forest plot showing hazard ratios
with confidence intervals for: ITT population, H subgroup, and H^c
complement.

## See also

[`forestsearch`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
for running the subgroup analysis
[`sg_consistency_out`](https://larry-leon.github.io/forestsearch/reference/sg_consistency_out.md)
for consistency evaluation
[`plot_subgroup_results_forestplot`](https://larry-leon.github.io/forestsearch/reference/plot_subgroup_results_forestplot.md)
for publication-ready forest plots
