# Assemble the influence matrix from candidate membership lists

For each candidate (a vector of row indices into `df`), fit the
treatment effect once and place its `dfbeta` into the candidate's rows
(zero elsewhere). Candidates whose fit fails or whose `dfbeta` length
does not match are dropped.

## Usage

``` r
.fs_mr_assemble(df, candidates, spec)
```

## Arguments

- df:

  Analysis data frame.

- candidates:

  Named list of integer row-index vectors.

- spec:

  See
  [`.fs_mr_pieces()`](https://larry-leon.github.io/forestsearch/reference/dot-fs_mr_pieces.md).

## Value

List with `B` (`N x S_kept`), `beta_hat`, `sigma_D`, `sizes`, `names`,
`keep`, `log_scale`.
