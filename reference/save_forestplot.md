# Save ForestSearch Forest Plot to File

Saves a forest plot to a file (PDF, PNG, etc.) with explicit dimensions.

## Usage

``` r
save_forestplot(x, filename, width = 12, height = 10, dpi = 300, bg = "white")
```

## Arguments

- x:

  An fs_forestplot object.

- filename:

  Character. Output filename. Extension determines format.

- width:

  Numeric. Plot width in inches. Default: 12.

- height:

  Numeric. Plot height in inches. Default: 10.

- dpi:

  Numeric. Resolution for raster formats. Default: 300.

- bg:

  Character. Background color. Default: "white".

## Value

Invisibly returns the filename.

## Examples

``` r
if (FALSE) { # \dontrun{
# Create plot with custom theme
large_theme <- create_forest_theme(base_size = 14, row_padding = c(6, 4))

result <- plot_subgroup_results_forestplot(
  fs_results = list(fs.est = fs, fs_bc = fs_bc),
  df_analysis = df.analysis,
  outcome.name = "time",
  event.name = "status",
  treat.name = "treatment",
  theme = large_theme
)

# Save to file
save_forestplot(result, "forest_plot.pdf", width = 14, height = 12)
} # }
```
