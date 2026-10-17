# Plot Subgroup Analysis Results

Creates diagnostic plots for subgroup treatment effects from df_super
object

## Usage

``` r
plot_subgroup_effects(
  df_super,
  z,
  hrz_crit = 0,
  log.hrs = NULL,
  ahr_empirical = NULL,
  plot_type = c("both", "profile", "ahr"),
  add_rug = TRUE,
  zpoints_by = 1,
  ...
)
```

## Arguments

- df_super:

  A data frame containing subgroup analysis results with columns:
  loghr_po (log hazard ratios), and optionally theta_1 and theta_0
  (treatment-specific parameters)

- z:

  Character string specifying the column name to use as the subgroup
  score (e.g., "z_age", "z_size", "subgroup"). Required.

- hrz_crit:

  Critical hazard ratio threshold for defining optimal subgroup. Default
  is 1 (HR=1 on log scale is 0).

- log.hrs:

  Optional vector of reference log hazard ratios to display as
  horizontal lines. Default is NULL.

- ahr_empirical:

  Optional empirical average hazard ratio to display. If NULL,
  calculated from data. Default is NULL.

- plot_type:

  Character string specifying plot type: "both" (default), "profile", or
  "ahr".

- add_rug:

  Logical indicating whether to add rug plot of z values. Default is
  TRUE.

- zpoints_by:

  Step size for z-axis grid when calculating AHR curves. Default is 1.

- ...:

  Additional graphical parameters passed to plot()

## Value

A list containing:

- cut.zero:

  The minimum z value where loghr_po \< hrz_crit

- AHR_opt:

  Average hazard ratio for optimal subgroup (z \>= cut.zero)

- zpoints:

  Grid of z values used for AHR calculations

- HR.zpoints:

  AHR for population with z \>= zpoints

- HRminus.zpoints:

  AHR for population with z \<= zpoints

- HR2.zpoints:

  Alternative AHR calculation for z \>= zpoints

- HRminus2.zpoints:

  Alternative AHR calculation for z \<= zpoints

## Details

The function creates up to two plots:

1.  Treatment effect profile: Shows log hazard ratio as function of z

2.  Average hazard ratio curve: Shows AHR for subgroups z \>= threshold

The "optimal" subgroup is defined as patients with z \>= cut.zero, where
cut.zero is the minimum z value with favorable treatment effect (loghr
\< hrz_crit).

## Examples

``` r
if (FALSE) { # \dontrun{
# Using z_age as the subgroup score
results <- plot_subgroup_effects(dgm_spline$df_super, z = "z_age", hrz_crit = 0)

# Using subgroup identifier
results <- plot_subgroup_effects(dgm_spline$df_super, z = "subgroup", hrz_crit = 0)

# With reference lines
results <- plot_subgroup_effects(dgm_spline$df_super, z = "z_size",
                                  hrz_crit = 0,
                                  log.hrs = c(-0.5, 0, 0.5))

# Only AHR plot
results <- plot_subgroup_effects(dgm_spline$df_super, z = "z_pgr",
                                  plot_type = "ahr")
} # }
```
