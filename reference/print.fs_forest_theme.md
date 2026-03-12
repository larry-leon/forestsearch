# Print Method for ForestSearch Forest Theme

Print Method for ForestSearch Forest Theme

## Usage

``` r
# S3 method for class 'fs_forest_theme'
print(x, ...)
```

## Arguments

- x:

  An fs_forest_theme object

- ...:

  Additional arguments (ignored)

## Examples

``` r
theme <- create_forest_theme()
print(theme)
#> ForestSearch Forest Plot Theme
#> ==============================
#> Base size:        10 
#> Scale:            1 
#> Size factor:      1 x (relative to base_size=10)
#> 
#> Calculated values:
#>   Body font:      10 
#>   Header font:    11 
#>   CV font:        10 
#>   Footnote font:  9 
#>   Row padding:    4, 4  mm
#>   CI line width:  1.5 
#> 
#> Use with: plot_subgroup_results_forestplot(..., theme = x)
```
