# Selected-covariate-pair frequencies

Bar-chart companion to
[`fs_rule_covariate_pairs()`](https://larry-leon.github.io/forestsearch/reference/fs_rule_covariate_pairs.md),
in the same idiom as the selected-covariate figure: horizontal bars,
ordered by share, coloured by how the signature relates to the true
region.

## Usage

``` r
plot_fs_rule_pairs(
  x,
  true_covariates = NULL,
  top_n = 12L,
  detected_only = TRUE
)
```

## Arguments

- x:

  An `fs_identification` object, or anything
  [`fs_rule_covariate_pairs()`](https://larry-leon.github.io/forestsearch/reference/fs_rule_covariate_pairs.md)
  accepts.

- true_covariates:

  Character vector; taken from `x` when it is an `fs_identification`
  object.

- top_n:

  Integer. Show only the most frequent `top_n` signatures; the remainder
  are pooled into a single "other" bar so the shares still sum to 1 and
  nothing is silently dropped.

- detected_only:

  Passed through.

## Value

A `ggplot` object.
