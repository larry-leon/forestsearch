# Infinitesimal-jackknife variance of the de-biased (bagged) estimate

Implements Leon et al. (2024) Eq. (VInfJ)/(VInfJ_bc) in multiplier form.
The multiplier matrix `Xi` supplies the centered bootstrap
multiplicities \\K^\*\_{bi} - \bar K^\*\_i\\, and `r` is the per-draw
residual \\r_b = \hat\beta(\widehat H) - \eta^\*\_b(\widehat H^\*\_b) -
\eta^\*\_b(\widehat H) - \hat\beta^\*(\widehat H)\\, i.e.
`(selection_bias + fixed_bias) - D_{H*_b}(b) - D_{H}(b)`, evaluated only
on the draws in `ok` that produced a re-selected winner.

## Usage

``` r
.fs_mr_ij_var(Xi, r, ok)
```

## Arguments

- Xi:

  `N x draws` multiplier matrix (one column per draw).

- r:

  Length-`draws` residual vector (`NA` where no winner).

- ok:

  Integer indices of usable draws.

## Value

List with `tilde_V` (raw IJ), `hat_V` (Wager 2014 bias-corrected), and
`B_ok` (number of usable draws).
