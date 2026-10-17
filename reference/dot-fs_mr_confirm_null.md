# Default near-null harm-confirmation threshold for a measure

The value `t_confirm` takes when left `NULL`: the null effect on the
effect scale, i.e. no effect at all rather than the screening threshold.

## Usage

``` r
.fs_mr_confirm_null(effect_measure)
```

## Arguments

- effect_measure:

  Character measure label (`"OR"`, `"RR"`, `"IRR"`, `"HR"`, `"RD"`,
  `"MD"`, ...) or `NULL`.

## Value

`1` for ratio measures (effect scale), `0` for differences.
