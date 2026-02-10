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
