# Multiply every hazard of a Weibull AFT DGM by a constant

Sets the AFT intercept of a DGM to `mu_fitted - tau * log(m)`, where
`mu_fitted` is the intercept as fitted when the DGM was built. Under the
Weibull AFT model the baseline cumulative hazard is \\\Lambda_0(s) = (s
e^{-\mu})^{1/\tau}\\, so every patient's hazard under either arm is
multiplied by `m` and every hazard ratio (individual, AHR, CDE) is
unchanged; only event probabilities in a finite window move. The scaling
is not cumulative: it is always applied to the fitted intercept, so
`scale_dgm_hazard(scale_dgm_hazard(dgm, 2), 3)` equals
`scale_dgm_hazard(dgm, 3)`.

## Usage

``` r
scale_dgm_hazard(dgm, m)
```

## Arguments

- dgm:

  An `aft_dgm_flex` object, typically from
  [`inject_mrct_structure()`](https://larry-leon.github.io/forestsearch/reference/inject_mrct_structure.md).

- m:

  Positive number: the hazard multiplier. `m = 1` restores the fitted
  intercept (for an unscaled DGM, the DGM is returned unchanged except
  for the record).

## Value

The DGM with `model_params$mu` set, the multiplier recorded in
`dgm$mrct$hazard_multiplier` and the fitted intercept in
`dgm$mrct$mu_fitted`, and the attribute `"m"` set to `m` (as the
mrct-consistency helper `ideal_set_m()` sets it, so scripts reading
`attr(dgm, "m")` keep working).

## Details

The fitted intercept is taken, in order, from `dgm$mrct$mu_fitted`
(stored by an earlier call), `attr(dgm, "mu_fitted")` (the attribute the
mrct-consistency builders set), or, for a DGM never scaled,
`dgm$model_params$mu`. A DGM that records a multiplier but no fitted
intercept is refused.

## See also

[`inject_mrct_structure()`](https://larry-leon.github.io/forestsearch/reference/inject_mrct_structure.md),
[`ahr_population()`](https://larry-leon.github.io/forestsearch/reference/ahr_population.md)

## Examples

``` r
if (FALSE) { # \dontrun{
dgm3 <- scale_dgm_hazard(dgm, 3)   # every hazard x 3, HRs unchanged
} # }
```
