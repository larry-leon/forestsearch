# Compute and Attach CDE Values to a DGM Object

Calculates Controlled Direct Effect (CDE) hazard ratios from the
super-population potential outcomes (`theta_0`, `theta_1`) and attaches
them to the DGM's `hazard_ratios` list. This enables automatic CDE
detection by
[`build_estimation_table`](https://larry-leon.github.io/forestsearch/reference/build_estimation_table.md).

## Usage

``` r
compute_dgm_cde(dgm, harm_col = NULL)
```

## Arguments

- dgm:

  A DGM object (e.g., from `create_gbsg_dgm` or
  `generate_aft_dgm_flex`). Must contain `df_super_rand` with columns
  `theta_0` and `theta_1`.

- harm_col:

  Character. Name of the subgroup indicator column in
  `dgm$df_super_rand`. If `NULL` (default), auto-detected from
  `flag.harm`, `flag_harm`, or `H`.

## Value

The DGM object with CDE values added to `dgm$hazard_ratios` (`CDE`,
`CDE_harm`, `CDE_no_harm`) and to top-level fields (`dgm$CDE`,
`dgm$cde_H`, `dgm$cde_Hc`).

## Details

The CDE for subgroup S is defined as:

`CDE(S) = mean(exp(theta_1[S])) / mean(exp(theta_0[S]))`

which is the ratio of average hazard contributions on the natural scale.
This differs from the AHR (`exp(mean(loghr_po))`) due to Jensen's
inequality. In the notation of Leon et al. (2024), CDE corresponds to
theta-ddagger.

The function detects the subgroup indicator column automatically,
checking for `flag.harm`, `flag_harm`, and `H` in the super-population
data frame.

## See also

[`build_estimation_table`](https://larry-leon.github.io/forestsearch/reference/build_estimation_table.md),
[`get_dgm_hr`](https://larry-leon.github.io/forestsearch/reference/get_dgm_hr.md)

## Examples

``` r
if (FALSE) { # \dontrun{
dgm <- create_gbsg_dgm(model = "alt", k_inter = 2.0)
dgm <- compute_dgm_cde(dgm)
dgm$hazard_ratios$CDE_harm   # theta-ddagger(H)
dgm$hazard_ratios$CDE         # theta-ddagger overall
} # }
```
