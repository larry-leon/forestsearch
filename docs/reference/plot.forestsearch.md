# Plot ForestSearch Results

Dispatches to
[`plot_sg_results`](https://larry-leon.github.io/forestsearch/reference/plot_sg_results.md)
for Kaplan-Meier curves, hazard-ratio forest plots, or combined panels.

## Usage

``` r
# S3 method for class 'forestsearch'
plot(
  x,
  type = c("combined", "km", "forest", "summary"),
  outcome.name = "Y",
  event.name = "Event",
  treat.name = "Treat",
  ...
)
```

## Arguments

- x:

  A `forestsearch` object returned by
  [`forestsearch`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md).

- type:

  Character. Type of plot:

  `"combined"`

  :   KM curves + forest plot (default)

  `"km"`

  :   Kaplan-Meier survival curves only

  `"forest"`

  :   Hazard-ratio forest plot only

  `"summary"`

  :   Summary statistics panel

- outcome.name:

  Character. Name of time-to-event column. Default: `"Y"`.

- event.name:

  Character. Name of event indicator column. Default: `"Event"`.

- treat.name:

  Character. Name of treatment column. Default: `"Treat"`.

- ...:

  Additional arguments passed to
  [`plot_sg_results`](https://larry-leon.github.io/forestsearch/reference/plot_sg_results.md),
  such as `by.risk`, `conf.level`, `est.scale`, `sg0_name`, `sg1_name`,
  `treat_labels`, `colors`, `title`, `show_events`, `show_ci`,
  `show_logrank`, `show_hr`.

## Value

Invisibly returns the plot result from
[`plot_sg_results`](https://larry-leon.github.io/forestsearch/reference/plot_sg_results.md).

## See also

[`plot_sg_results`](https://larry-leon.github.io/forestsearch/reference/plot_sg_results.md)
for full control over appearance,
[`plot_sg_weighted_km`](https://larry-leon.github.io/forestsearch/reference/plot_sg_weighted_km.md)
for weighted KM curves,
[`plot_subgroup_results_forestplot`](https://larry-leon.github.io/forestsearch/reference/plot_subgroup_results_forestplot.md)
for publication-ready forest plots.

## Examples

``` r
if (FALSE) { # \dontrun{
fs <- forestsearch(df.analysis = mydata, ...)

# Combined KM + forest plot (default)
plot(fs)

# KM curves only
plot(fs, type = "km")

# Forest plot only
plot(fs, type = "forest")

# With non-standard column names
plot(fs, type = "km",
     outcome.name = "os_months",
     event.name = "os_event",
     treat.name = "treatment")

# With custom labels
plot(fs, sg0_name = "High Risk", sg1_name = "Standard Risk",
     treat_labels = c("0" = "Placebo", "1" = "Active Drug"))
} # }
```
