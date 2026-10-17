# Refuse to print coverage targets computed on different denominators

`.coverage_meta()`-style assembly computes
`ok <- is.finite(target) & is.finite(lo) & is.finite(hi)` and
`n_eff <- sum(ok)`, so a non-finite target is dropped with no error and
no warning. Coverage is then computed over the replicates whose targets
happened to score, which is not a random subset: whatever made a target
non-finite is a property of the realized rule.

## Usage

``` r
fs_betaHhat_neff_parity(
  cov_df,
  strict = TRUE,
  target_beta = "C_betaHhat",
  target_ref = "C_dagger"
)
```

## Arguments

- cov_df:

  Coverage table with `target`, `n_eff`, `block` and `estimator`
  columns.

- strict:

  Logical. `TRUE` stops; `FALSE` warns, for a run whose cause is known.

- target_beta, target_ref:

  Character. The target to check and the reference to check it against.

## Value

`cov_df`, invisibly.

## Details

Because a reference target that is a finite scalar is dropped only when
the interval itself is non-finite, `n_eff` parity between the two is the
exact check for a target-specific loss, and it is free – `n_eff` is
already in the table.

## Two causes, deliberately not distinguished

A parity break can mean an unresolvable rule, or a genuinely degenerate
region whose per-family guard returned `NA` by design. Both make the two
coverage figures incomparable as printed, so both stop. Use
[`fs_betaHhat_neff_report()`](https://larry-leon.github.io/forestsearch/reference/fs_betaHhat_neff_report.md)
on the bundle to tell them apart before deciding.
