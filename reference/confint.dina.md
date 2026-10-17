# Wald confidence intervals for DINA coefficients.

Computes normal-approximation confidence intervals \\\hat\beta_j \pm
z\_{1 - \alpha/2} \sqrt{\widehat{V}\_{jj}}\\ using the sandwich (or Cox
model-based) variance returned by
[`dina_fit()`](https://larry-leon.github.io/forestsearch/reference/dina_fit.md).

## Usage

``` r
# S3 method for class 'dina'
confint(object, parm, level = 0.95, ...)
```

## Arguments

- object:

  an object of class `"dina"`.

- parm:

  optional character or integer vector selecting which coefficients to
  return intervals for; defaults to all.

- level:

  confidence level (default `0.95`).

- ...:

  unused.

## Value

a matrix with two columns giving lower and upper limits, named in
standard [`confint()`](https://rdrr.io/r/stats/confint.html) style (e.g.
`"2.5 %"`, `"97.5 %"`).
