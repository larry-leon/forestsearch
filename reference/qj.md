# k-th J-Quantile

Returns the (k/J)-th quantile of a numeric vector. Used by
[`cut_var_jq()`](https://larry-leon.github.io/forestsearch/reference/cut_var_jq.md)
to emit deferred cut expressions of the form `"X <= qj(X, k, J)"`, which
are then resolved to literal numerics by
[`process_conf_force_expr()`](https://larry-leon.github.io/forestsearch/reference/process_conf_force_expr.md).

## Usage

``` r
qj(x, k, J)
```

## Arguments

- x:

  A numeric vector.

- k:

  Integer in `1, ..., J - 1L`. Index of the cut point (where `J` here is
  the third argument to `qj()`, not the user-facing `conf.cont_jcuts`
  value – see
  [`cut_var_jq()`](https://larry-leon.github.io/forestsearch/reference/cut_var_jq.md)
  for how these relate).

- J:

  Integer \>= 2. Total number of intervals implied by the probability
  `k/J`. In the cut expressions emitted by `cut_var_jq(x, J_user)`, this
  argument is `J_user + 1`.

## Value

Numeric value of the (k/J)-th quantile of `x`.
