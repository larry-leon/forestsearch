# Quick Plot KM Bands from ForestSearch

Convenience wrapper with sensible defaults for quick visualization.

## Usage

``` r
quick_km_band_plot(df, fs.est, outcome.name, event.name, treat.name, ...)
```

## Arguments

- df:

  Data frame with analysis data.

- fs.est:

  ForestSearch result object.

- outcome.name:

  Character. Time-to-event column name.

- event.name:

  Character. Event indicator column name.

- treat.name:

  Character. Treatment column name.

- ...:

  Additional arguments passed to
  [`plot_km_band_forestsearch()`](https://github.com/larry-leon/forestsearch/reference/plot_km_band_forestsearch.md),
  such as `ref_subgroups`, `tau_add`, `by_risk`, etc.

## Value

Invisibly returns the plot result.

## Examples

``` r
if (FALSE) { # \dontrun{
# Quick plot with minimal configuration
quick_km_band_plot(df_analysis, fs_result, "os_months", "os_event", "treat")

# With reference subgroups
ref_sgs <- list(
  age_young = list(subset_expr = "age < 65", color = "brown"),
  age_old = list(subset_expr = "age >= 65", color = "grey")
)
quick_km_band_plot(df_analysis, fs_result, "os_months", "os_event", "treat",
                   ref_subgroups = ref_sgs, tau_add = 60)
} # }
```
