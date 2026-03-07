# Set Up a GBSG-Based AFT Data Generating Mechanism

Creates a GBSG-based data generating mechanism that is fully compatible
with
[`simulate_from_dgm`](https://larry-leon.github.io/forestsearch/reference/simulate_from_dgm.md).
This is the replacement for
[`create_gbsg_dgm()`](https://larry-leon.github.io/forestsearch/reference/create_gbsg_dgm.md):
it accepts exactly the same arguments and produces the same numeric
output, but returns an object of class `"aft_dgm_flex"` instead of
`"gbsg_dgm"`.

## Usage

``` r
setup_gbsg_dgm(
  model = c("alt", "null"),
  k_treat = 1,
  k_inter = 1,
  k_z3 = 1,
  z1_quantile = 0.25,
  n_super = 5000L,
  cens_type = c("weibull", "uniform"),
  use_rand_params = FALSE,
  seed = 8316951L,
  verbose = FALSE
)
```

## Arguments

- model:

  Character. Either "alt" for alternative hypothesis with heterogeneous
  treatment effects, or "null" for uniform treatment effect. Default:
  "alt"

- k_treat:

  Numeric. Treatment effect multiplier applied to the treatment
  coefficient from the fitted AFT model. Values \> 1 strengthen the
  treatment effect. Default: 1

- k_inter:

  Numeric. Interaction effect multiplier for the treatment-subgroup
  interaction (z1 \* z3). Only used when model = "alt". Higher values
  create more heterogeneity between HR(H) and HR(Hc). Default: 1

- k_z3:

  Numeric. Effect multiplier for the z3 (menopausal status) coefficient.
  Default: 1

- z1_quantile:

  Numeric. Quantile threshold for z1 (estrogen receptor). Observations
  with ER \<= quantile are coded as z1 = 1. Default: 0.25

- n_super:

  Integer. Size of super-population for empirical HR estimation.
  Default: 5000

- cens_type:

  Character. Censoring distribution type: "weibull" or "uniform".
  Default: "weibull"

- use_rand_params:

  Logical. If TRUE, modifies confounder coefficients using estimates
  from randomized subset (meno == 0). Default: FALSE

- seed:

  Integer. Random seed for super-population generation. Default: 8316951

- verbose:

  Logical. Print diagnostic information. Default: FALSE

## Value

An object of class `c("aft_dgm_flex", "gbsg_dgm", "list")` with all
fields from
[`create_gbsg_dgm()`](https://larry-leon.github.io/forestsearch/reference/create_gbsg_dgm.md)
plus:

- `df_super`:

  Super-population data frame with
  [`simulate_from_dgm()`](https://larry-leon.github.io/forestsearch/reference/simulate_from_dgm.md)-compatible
  column names.

- `model_params$tau`:

  Copy of `model_params$sigma`.

- `model_params$censoring`:

  Sub-list with `type`, `mu`, `tau` for the censoring model.

## Details

Internally the function calls
[`create_gbsg_dgm()`](https://larry-leon.github.io/forestsearch/reference/create_gbsg_dgm.md)
and then:

1.  Adds a `df_super` field with column names aligned to
    [`simulate_from_dgm()`](https://larry-leon.github.io/forestsearch/reference/simulate_from_dgm.md)
    conventions (`lin_pred_1`, `lin_pred_0`, `lin_pred_cens_1`,
    `lin_pred_cens_0`, `flag_harm`).

2.  Adds a `model_params$tau` field (= `model_params$sigma`) and a
    `model_params$censoring` sub-list.

3.  Sets class to `c("aft_dgm_flex", "gbsg_dgm", "list")`.

The original `df_super_rand` field is kept so that
[`compute_dgm_cde()`](https://larry-leon.github.io/forestsearch/reference/compute_dgm_cde.md)
and `print.gbsg_dgm` continue to work.

## See also

[`create_gbsg_dgm`](https://larry-leon.github.io/forestsearch/reference/create_gbsg_dgm.md),
[`simulate_from_dgm`](https://larry-leon.github.io/forestsearch/reference/simulate_from_dgm.md),
[`compute_dgm_cde`](https://larry-leon.github.io/forestsearch/reference/compute_dgm_cde.md)

## Examples

``` r
if (FALSE) { # \dontrun{
dgm <- setup_gbsg_dgm(model = "alt", k_inter = 2, verbose = FALSE)
dgm <- compute_dgm_cde(dgm)
print(dgm)
sim <- simulate_from_dgm(dgm, n = 400, seed = 1)
} # }
```
