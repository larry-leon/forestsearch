# Plot GLM Subgroup Outcomes

Creates a publication-ready visualization of identified subgroup
outcomes for binary or continuous endpoints. This is the GLM counterpart
to
[`plot_sg_weighted_km`](https://larry-leon.github.io/forestsearch/reference/plot_sg_weighted_km.md)
for survival endpoints.

## Usage

``` r
plot_sg_glm_outcomes(
  fs.est,
  fs_bc = NULL,
  outcome.name,
  treat.name,
  offset.name = NULL,
  effect_measure = NULL,
  outcome_type = NULL,
  adjust_covariates = NULL,
  sg0_name = "Questionable",
  sg1_name = "Recommend",
  E.name = "Treatment",
  C.name = "Control",
  bar_colors = c("steelblue", "coral"),
  show_effect = TRUE,
  show_bc = TRUE,
  show_itt = TRUE,
  conf.level = 0.95,
  title = NULL,
  subtitle = NULL,
  verbose = FALSE
)

# S3 method for class 'fs_binary_plot'
print(x, ...)

# S3 method for class 'fs_binary_plot'
plot(x, ...)
```

## Arguments

- fs.est:

  A forestsearch object containing `df.est` with `treat.recommend`
  column.

- fs_bc:

  Optional. Bootstrap results from
  [`forestsearch_bootstrap_dofuture`](https://larry-leon.github.io/forestsearch/reference/forestsearch_bootstrap_dofuture.md)
  for bias-corrected annotations.

- outcome.name:

  Character. Name of outcome column.

- treat.name:

  Character. Name of treatment column.

- offset.name:

  Character or `NULL`. Name of follow-up time column for rate-based
  measures (IRR, IRD).

- effect_measure:

  Character or `NULL`. Effect measure. Auto-detected from
  `fs.est$effect_measure` if `NULL`.

- outcome_type:

  Character or `NULL`. One of `"binary"`, `"continuous"`. Auto-detected
  from `fs.est$outcome_type` if `NULL`.

- adjust_covariates:

  Character vector or `NULL`. Baseline covariates to adjust the plotted
  effect model for, matching the adjustment used during
  search/estimation. When `NULL`, resolved from
  `fs.est$args_call_all$adjust_covariates` so the plot reflects the same
  (adjusted or unadjusted) scale as the fit. Columns must be present in
  `fs.est$df.est`.

- sg0_name:

  Character. Label for H subgroup (treat.recommend == 0). Default:
  `"Questionable"`.

- sg1_name:

  Character. Label for Hc subgroup (treat.recommend == 1). Default:
  `"Recommend"`.

- E.name:

  Character. Label for experimental/treatment arm. Default:
  `"Treatment"`.

- C.name:

  Character. Label for control arm. Default: `"Control"`.

- bar_colors:

  Character vector of length 2. Colors for control and treatment bars.
  Default: `c("steelblue", "coral")`.

- show_effect:

  Logical. Annotate with effect estimate. Default `TRUE`.

- show_bc:

  Logical. Annotate with bias-corrected estimate when `fs_bc` is
  provided. Default `TRUE`.

- show_itt:

  Logical. Include an ITT (full population) panel. Default `TRUE`.

- conf.level:

  Numeric. Confidence level for intervals. Default `0.95`.

- title:

  Character or `NULL`. Plot title.

- subtitle:

  Character or `NULL`. Plot subtitle.

- verbose:

  Logical. Print diagnostic messages. Default `FALSE`.

- x:

  An fs_binary_plot object.

- ...:

  Additional arguments (unused).

## Value

A list of class `"fs_binary_plot"` containing:

- plot:

  A ggplot2 object.

- data:

  Data frame of per-arm summaries used for the plot.

- effects:

  Data frame of effect estimates per subgroup.

- sg_definition:

  Character. Subgroup definition label.

- figure_note:

  Character. Auto-generated figure note.

## Details

For binary outcomes, displays grouped bar charts of event rates by
treatment arm within each subgroup. For continuous outcomes, displays
point-and-errorbar charts of means (+/- SE) by treatment arm. Both are
annotated with the effect estimate and optional bias-corrected estimate
from the bootstrap.

## Examples

``` r
if (FALSE) { # \dontrun{
fs <- forestsearch(df, outcome_type = "binary", effect_measure = "OR",
                   ...)
result <- plot_sg_glm_outcomes(fs, outcome.name = "status",
                                treat.name = "hormon")
print(result)
} # }
```
