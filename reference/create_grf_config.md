# Helper Functions for GRF Subgroup Analysis

This file contains helper functions used by grf.subg.harm.survival() to
improve readability and modularity. Create GRF configuration object

## Usage

``` r
create_grf_config(
  frac.tau,
  n.min,
  dmin.grf,
  RCT,
  sg.criterion,
  maxdepth,
  seedit
)
```

## Arguments

- frac.tau:

  Numeric. Fraction of tau for GRF horizon

- n.min:

  Integer. Minimum subgroup size

- dmin.grf:

  Numeric. Minimum difference in subgroup mean

- RCT:

  Logical. Is the data from a randomized controlled trial?

- sg.criterion:

  Character. Subgroup selection criterion

- maxdepth:

  Integer. Maximum tree depth

- seedit:

  Integer. Random seed

## Value

List with configuration parameters

## Details

Creates a configuration object to organize GRF parameters

## Examples

``` r
cfg <- create_grf_config(frac.tau = 0.6, n.min = 60, dmin.grf = 6,
                         RCT = TRUE, sg.criterion = "mDiff",
                         maxdepth = 2, seedit = 42L)
str(cfg)
#> List of 10
#>  $ frac.tau                 : num 0.6
#>  $ n.min                    : num 60
#>  $ dmin.grf                 : num 6
#>  $ RCT                      : logi TRUE
#>  $ sg.criterion             : chr "mDiff"
#>  $ maxdepth                 : num 2
#>  $ seedit                   : int 42
#>  $ valid_criteria           : chr [1:2] "mDiff" "Nsg"
#>  $ max_tree_depth           : num 3
#>  $ return_selected_cuts_only: logi FALSE
```
