# Which rules produced NA `beta(Hhat)` targets, and why

Diagnostic companion to
[`fs_betaHhat_neff_parity()`](https://larry-leon.github.io/forestsearch/reference/fs_betaHhat_neff_parity.md).
Separates the two causes of a parity break: a rule the evaluation frame
cannot express, and a region small enough for a per-family guard to
fire.

## Usage

``` r
fs_betaHhat_neff_report(bundle, block = c("H", "Hc"))
```

## Arguments

- bundle:

  A results data frame, a bundle list with a `results` element, or a
  path to an `.rds` holding either.

- block:

  Character. `"H"` or `"Hc"`.

## Value

A data frame of distinct offending rules with `is_disjunction`,
`nH_eval` and `status` where available, or `NULL` when there are none.
