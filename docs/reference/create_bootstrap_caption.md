# Calculate Bootstrap Table Caption

Generates an interpretive caption for bootstrap results table.

## Usage

``` r
create_bootstrap_caption(est.scale, nb_boots, boot_success_rate)
```

## Arguments

- est.scale:

  Character. "hr" or "1/hr"

- nb_boots:

  Integer. Number of bootstrap iterations

- boot_success_rate:

  Numeric. Proportion successful

## Value

Character string with caption
