# Three-Panel Effect-Distribution Violin from Simulation Results

Builds a three-panel violin plot (ITT / identified \\\hat H\\ /
identified \\\hat H^c\\) of treatment-effect estimates across simulation
replicates produced by
[`run_simulation_analysis()`](https://larry-leon.github.io/forestsearch/reference/run_simulation_analysis.md).
Handles both GLM endpoints (OR, RD, IRR) and survival (HR), optionally
inverts effects to a benefit scale (`subgroup_notation = "benefit"`),
and overlays a dotted reference line at the DGM-implied true effect in
\\\hat H\\.

## Usage

``` r
plot_effect_distribution(
  results,
  analysis_method = "FS",
  effect_measure = c("OR", "HR", "RD", "IRR", "MD"),
  subgroup_notation = c("harm", "benefit"),
  reference_effect = NA_real_,
  title = NULL,
  subtitle = NULL,
  panel_labels = NULL,
  drop_undetected = TRUE,
  trim_threshold = 1000,
  trim_fraction = 0.01,
  ylim = NULL
)
```

## Arguments

- results:

  `data.frame` returned by
  [`run_simulation_analysis()`](https://larry-leon.github.io/forestsearch/reference/run_simulation_analysis.md)
  (optionally filtered). Must contain the columns `analysis`, `any.H`,
  `hr.itt`, `hr.H.hat`, and `hr.Hc.hat`. The column names end in `.hr`
  regardless of the underlying effect measure (OR / HR / RD / IRR); this
  follows the package's established naming convention.

- analysis_method:

  Character. Method label in the `analysis` column to plot (e.g. `"FS"`,
  `"GRF"`). Default: `"FS"`.

- effect_measure:

  Character. One of `"OR"`, `"HR"`, `"RD"`, `"IRR"`, `"MD"`. Used only
  for axis labelling and to decide the y-axis reference line (\\1\\ for
  ratio measures, \\0\\ for risk / mean differences). `"MD"` (mean
  difference) is treated additively like `"RD"`. Default: `"OR"`.

- subgroup_notation:

  Character. `"harm"` leaves effects on the original (switched) scale;
  `"benefit"` inverts multiplicative effects via \\1/\hat\theta\\ and
  flips the sign of additive effects so that larger values favour the
  benefit subgroup. Default: `"harm"`.

- reference_effect:

  Numeric or `NA`. Optional overlay showing the DGM-implied true effect
  in \\\hat H\\ (e.g. `1 / dgm$hazard_ratios$harm_subgroup` on a benefit
  scale). Drawn as a dotted darkgreen line. Default: `NA`.

- title:

  Character. Plot title. Default: auto-constructed.

- subtitle:

  Character or `NULL`. Plot subtitle. Default: `NULL`.

- panel_labels:

  Named character vector with names `"itt"`, `"H"`, `"Hc"` and the
  display labels for each violin. Default:
  `c(itt = "ITT", H = "Identified H-hat", Hc = "Identified Hc-hat")`.
  When `subgroup_notation = "benefit"`, the default `Hc` label is
  adjusted to `"Identified Hc-hat (benefit)"`.

- drop_undetected:

  Logical. Restrict the \\\hat H\\ and \\\hat H^c\\ panels to replicates
  where a subgroup was identified (`any.H == 1`). Default: `TRUE`. The
  ITT panel always uses all replicates.

- trim_threshold:

  Numeric or `NULL`. Trimming is applied only if any of the three groups
  has a raw mean exceeding this absolute value – i.e., extreme outliers
  that distort the violin or push annotated mean/SD into scientific
  notation. Default: `1000`, matching
  [`build_estimation_table`](https://larry-leon.github.io/forestsearch/reference/build_estimation_table.md)
  so the plot and the companion gt table auto-trim on the same
  condition. Set `NULL` to disable trimming entirely; set a smaller
  value (e.g., `100`) to be more aggressive on plots where reasonable OR
  / HR / IRR estimates sit well below 100 and even moderate outliers
  visibly distort the violin axis.

- trim_fraction:

  Numeric in (0, 0.5). Fraction of observations to trim from each tail
  of each group when trimming triggers. Default: `0.01` (1\\ of values
  per group is kept). Has no effect when `trim_threshold = NULL`.

- ylim:

  Numeric vector of length 2 or `NULL`. Explicit y-axis limits as
  `c(lower, upper)` applied via
  [`coord_cartesian`](https://ggplot2.tidyverse.org/reference/coord_cartesian.html)
  (data are *not* dropped, only the viewport is clipped). Useful when
  extreme outlier estimates stretch the axis beyond the interesting
  range — e.g., `ylim = c(0, 4)` for OR benefit-scale plots. When `NULL`
  (default), limits are determined automatically by ggplot2.

## Value

A `ggplot2` object. Has `attr(p, "panel_data")` set to the long-format
`data.table` used for plotting (for downstream summaries or diagnostic
inspection). When trimming is active, also has `attr(p, "trim_info")`
containing per-group diagnostics (`n_total`, `n_trimmed`, `n_flagged`,
`raw_mean`, `raw_sd`, `trimmed_mean`, `trimmed_sd`, `lower_bound`,
`upper_bound`).

## Details

The MRCT analogue is
[`SGplot_estimates`](https://larry-leon.github.io/forestsearch/reference/SGplot_estimates.md),
which works on the four-panel
[`mrct_region_sims()`](https://larry-leon.github.io/forestsearch/reference/mrct_region_sims.md)
schema (ITT / training / testing / subgroup). This function targets the
three-panel
[`run_simulation_analysis()`](https://larry-leon.github.io/forestsearch/reference/run_simulation_analysis.md)
schema that does not involve a train/test split.

## See also

[`run_simulation_analysis`](https://larry-leon.github.io/forestsearch/reference/run_simulation_analysis.md)
for the simulation pipeline,
[`SGplot_estimates`](https://larry-leon.github.io/forestsearch/reference/SGplot_estimates.md)
for the MRCT analogue,
[`build_estimation_table`](https://larry-leon.github.io/forestsearch/reference/build_estimation_table.md)
for the tabular counterpart.

## Examples

``` r
if (FALSE) { # \dontrun{
# Binary endpoint, benefit search, OR scale
true_or_benefit <- 1 / dgm_calibrated$hazard_ratios$harm_subgroup
p <- plot_effect_distribution(
  results_alt,
  analysis_method   = "FS",
  effect_measure    = "OR",
  subgroup_notation = "benefit",
  reference_effect  = true_or_benefit,
  title             = "ACTG175 Binary (HTE): OR Estimates Across Simulations",
  subtitle          = sprintf("n = 1000, 50 sims | truth OR(Q) = %.2f",
                              true_or_benefit)
)
print(p)

# Survival endpoint, harm search, HR scale (no inversion)
p_surv <- plot_effect_distribution(
  results_alt,
  effect_measure    = "HR",
  subgroup_notation = "harm",
  reference_effect  = dgm_alt$hazard_ratios$harm_subgroup
)

# Clip y-axis to focus on the interesting OR range
p_clipped <- plot_effect_distribution(
  results_alt,
  effect_measure    = "OR",
  subgroup_notation = "benefit",
  reference_effect  = true_or_benefit,
  ylim              = c(0, 4)
)
} # }
```
