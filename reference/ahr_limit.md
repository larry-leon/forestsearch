# The exact limit of the two-sample Cox estimate on a patient set

For a set of patients of an
[`ahr_population()`](https://larry-leon.github.io/forestsearch/reference/ahr_population.md),
the limit \\\beta^\*\\ of the log hazard ratio estimated by a Cox fit of
treatment alone on that set (1:1 randomisation, analysis at the
population's calendar time with each patient's uniform entry window,
administrative censoring only): the root of the limiting score
(Struthers and Kalbfleisch), with the information, the model-based and
robust (Lin-Wei) variances per patient, the event probabilities, the
type table and \\\omega\\ (the H\* share of events), and the set's AHR,
CDE and its margHR at infinity (no censoring). Computed by quadrature;
nothing is drawn.

## Usage

``` r
ahr_limit(
  pop,
  set,
  b = NULL,
  h = NULL,
  psi = NULL,
  horizon = NULL,
  variance = TRUE,
  b_R = NULL,
  h_R = NULL
)
```

## Arguments

- pop:

  An
  [`ahr_population()`](https://larry-leon.github.io/forestsearch/reference/ahr_population.md).

- set:

  A set name: `"all"`, `"nonregion"`, `"region"`, `"region_Gstar"`
  (region outside H\*), `"Hstar"`, `"Hstar_region"`,
  `"Hstar_nonregion"`; or a logical vector of length `pop$N`.

- b, h:

  Hazard ratios outside and inside H\*: the treated log hazard is
  \\\theta_0 + \log b + \log(h/b) 1\\H^\*\\\\. Both `NULL` (and `psi`
  `NULL`): the DGM's own per-patient log hazard ratios (`loghr_po`).

- psi:

  `NULL`, or a per-patient log hazard ratio vector of length `pop$N`
  (overrides `b`, `h`, `b_R`, `h_R`).

- horizon:

  Calendar time of the analysis; `NULL` (default) is `pop$x`; `Inf`
  gives the uncensored limit.

- variance:

  Logical: compute the robust variance (default `TRUE`).

- b_R, h_R:

  Hazard ratios inside the region, outside and inside H\*: with `b` and
  `h` they give the four-cell treated log hazard \\\theta_0 + \log b\\
  (outside H\*, outside the region), \\\log h\\ (H\*, outside), \\\log
  b_R\\ (outside H\*, region), \\\log h_R\\ (H\*, region). `NULL`
  (default) means `b_R = b` and `h_R = h`, the two-cell model. Used only
  with `b` and `h`. When they differ from `b`, `h`, a type group whose
  cells carry different hazard ratios has `psi = NA` in `types`.

## Value

A list: `beta` (the limit) and `hr` (margHR at the horizon), `i1`
(information per patient), `var_model`, `var_robust`, `var_robust_fixed`
(allocation fixed at n/2) (variances per patient, so divide by the set's
size), `sandwich`, `d1`, `d0`, `d` (event probabilities by arm and at
1:1, i.e. expected events per patient), `types` (per type group:
patients, log-HR, events, event share, own and pooled imbalance
moments), `omega`, `ahr`, `cde`, `band_share`, `omega0`, `marghr_inf`,
`n_pop`, `share` (\\n_S/N\\), `horizon`, and the node arrays and `grid`
used by
[`ahr_design()`](https://larry-leon.github.io/forestsearch/reference/ahr_design.md).

## Details

A standard set at given (`b`, `h`) with \\0.2 \le b, h \le 30\\ uses the
population's cached grid and Laplace sums (the route of
[`ahr_design()`](https://larry-leon.github.io/forestsearch/reference/ahr_design.md));
any other call (a logical set, `psi`, or the DGM's own effects) builds
the set's own grid.

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

## References

Struthers, C. A. and Kalbfleisch, J. D. (1986). Misspecified
proportional hazards models. *Biometrika* 73, 363-369.

## See also

[`ahr_population()`](https://larry-leon.github.io/forestsearch/reference/ahr_population.md),
[`ahr_design()`](https://larry-leon.github.io/forestsearch/reference/ahr_design.md)

## Examples

``` r
if (FALSE) { # \dontrun{
pop <- ahr_population(dgm, 40, 11)
e <- ahr_limit(pop, "region", b = 0.55, h = 2)
c(margHR40 = e$hr, AHR = e$ahr, CDE = e$cde, margHR_inf = e$marghr_inf)
} # }
```
