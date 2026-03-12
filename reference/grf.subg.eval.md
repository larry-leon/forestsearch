# GRF Subgroup Evaluation and Performance Metrics

Evaluates the performance of GRF-identified subgroups, including hazard
ratios, bias, and predictive values. This function is typically used in
simulation studies to assess the performance of the GRF subgroup
identification method.

## Usage

``` r
grf.subg.eval(
  df,
  grf.est,
  dgm,
  cox.formula.sim,
  cox.formula.adj.sim,
  analysis = "GRF",
  frac.tau = 1
)
```

## Arguments

- df:

  Data frame containing the analysis data.

- grf.est:

  List. Output from `grf.subg.harm.survival`.

- dgm:

  List. Data-generating mechanism (truth) for simulation.

- cox.formula.sim:

  Formula for unadjusted Cox model.

- cox.formula.adj.sim:

  Formula for adjusted Cox model.

- analysis:

  Character. Analysis label (default: "GRF").

- frac.tau:

  Numeric. Fraction of tau for GRF horizon (default: 1.0).

## Value

A data frame with evaluation metrics.

## Examples

``` r
if (FALSE) { # \dontrun{
# grf.subg.eval() is called internally to evaluate GRF subgroup quality.
# See grf.subg.harm.survival() for the standard entry point.
} # }
```
