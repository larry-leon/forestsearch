# Default ForestSearch Parameters for GBSG Simulations

Returns a list of default parameters for ForestSearch analysis in
GBSG-based simulations.

## Usage

``` r
default_fs_params()
```

## Value

List of default ForestSearch parameters

## Details

Default parameters are optimized for GBSG simulation scenarios with
moderate sample sizes (300-1000) and typical event rates.

Variable selection defaults:

- use_lasso = TRUE: LASSO-based variable importance (default for FS)

- use_grf = FALSE: GRF-based variable importance (enable for FSlg)

The `use_twostage` parameter is set to FALSE by default for backward
compatibility. Set to TRUE for faster exploratory analyses.
