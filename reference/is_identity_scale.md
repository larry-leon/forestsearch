# Is this Effect Measure on the Identity Scale?

Identity-scale measures (MD, RD, IRD) use additive differences;
ratio-scale measures (HR, OR, RR, IRR) use multiplicative ratios. The
distinction matters for benefit-search inversion: ratio-scale values are
reciprocated (`1/x`), identity-scale values are negated (`-x`).

## Usage

``` r
is_identity_scale(effect_measure)
```

## Arguments

- effect_measure:

  Character. One of `"HR"`, `"OR"`, `"RR"`, `"IRR"`, `"RD"`, `"IRD"`,
  `"MD"`, or `NULL` (treated as ratio-scale for backward compatibility).

## Value

Logical. `TRUE` for identity-scale measures.
