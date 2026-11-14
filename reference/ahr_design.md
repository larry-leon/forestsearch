# Design quantities and trigger probabilities at a trial size

At trial size `n`, for the standard sets of an
[`ahr_population()`](https://larry-leon.github.io/forestsearch/reference/ahr_population.md)
with hazard ratio `b` outside H\* and `h` inside: each set's limit and
its model-based and robust SEs (using \\E\[1/n_S\]\\, \\n_S \sim Bin(n,
\pi_S)\\), expected events, the exact and SE-ratio correlations between
nested sets' estimates, and, from the trivariate normal approximation
(estimates \\\sim N(\beta^\*, SE\_{robust}^2)\\, exact correlations, the
Wald limit of the ITT with the model SE):

- `p_itt_ub_lt1`: P(ITT upper 95\\

- `p_region_itt_gt1`: P(region ITT estimate \> 1);

- `p_region_Gstar_lt1`: P(region(G\*) estimate \< 1);

- `p_trigger`: P(ITT ub \< 1 and region ITT \> 1);

- `p_trigger_and_Gstar_lt1`: P(trigger and region(G\*) \< 1), the oracle
  joint.

Probabilities use `mvtnorm::pmvnorm(algorithm = Miwa(steps = 1024))`.

## Usage

``` r
ahr_design(
  pop,
  b,
  h,
  n,
  sets = c("all", "region", "region_Gstar"),
  b_R = b,
  h_R = h
)
```

## Arguments

- pop:

  An
  [`ahr_population()`](https://larry-leon.github.io/forestsearch/reference/ahr_population.md).

- b, h:

  Hazard ratios outside and inside H\*. All sets use the population's
  cached grid (built for \\0.2 \le b, h \le 30\\; outside that range the
  grid's first panel is coarser than the library's).

- n:

  Trial size (both arms).

- sets:

  Standard set names to report (see
  [`ahr_limit()`](https://larry-leon.github.io/forestsearch/reference/ahr_limit.md));
  the trigger sets `"all"`, `"region"` and `"region_Gstar"` are always
  computed.

- b_R, h_R:

  Hazard ratios inside the region, outside and inside H\* (the four-cell
  model of
  [`ahr_limit()`](https://larry-leon.github.io/forestsearch/reference/ahr_limit.md));
  default `b_R = b`, `h_R = h`, the two-cell model.

## Value

A list: `sets` (a data.frame per set: `set`, `share`, `beta`, `marghr`
(at the analysis time), `se_model`, `se_robust`, `se_d4`
(\\\sqrt{4/events}\\), `info_ratio`, `events` (expected, of `n`), `ahr`,
`cde`), `corr` (a data.frame per nested pair: `big`, `small`, `exact`,
`se_ratio`, `pi`), `prob` (the named probabilities above), and `b`, `h`,
`n`, `b_R`, `h_R`.

## Limits

Weibull AFT DGMs from
[`inject_mrct_structure()`](https://larry-leon.github.io/forestsearch/reference/inject_mrct_structure.md)
only; entry uniform on \\(a_i, a_i + R)\\ for each entry group;
administrative censoring only (the random censoring model of the DGM,
e.g. `cens_adjust`, is not in the layer); correlations exact for nested
sets. Each distinct delay is one quadrature group with its own breaks,
so a delay vector with many distinct values is slow. Derivations: the
technical report (`docs/manuscript/techreport/analytic-layer.qmd`,
Chapter 6, of the mrct-consistency repository) and the theory notes
(`docs/theory/AHR_mrct_dgm.qmd`, `docs/theory/AHR_complete.qmd`).

## See also

[`ahr_limit()`](https://larry-leon.github.io/forestsearch/reference/ahr_limit.md),
[`ahr_solve()`](https://larry-leon.github.io/forestsearch/reference/ahr_solve.md)

## Examples

``` r
if (FALSE) { # \dontrun{
pop <- ahr_population(dgm, 40, 11)
ahr_design(pop, b = 0.55, h = 2, n = 1000)$prob
} # }
```
