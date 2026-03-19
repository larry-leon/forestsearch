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
