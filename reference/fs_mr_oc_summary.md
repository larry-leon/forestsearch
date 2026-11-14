# Operating-Characteristics Summary of a Multiplier-Resampling Payload

Summarises a simulation payload produced by the multiplier-resampling
(MR) harness into the estimation, coverage and classification quantities
used in simulation summaries, returning plain numbers rather than a
formatted table.

## Usage

``` r
fs_mr_oc_summary(
  payload,
  estimators = NULL,
  blocks = c("H", "Hc"),
  digits = NULL
)
```

## Arguments

- payload:

  A list with `results`, `truth` and `meta`, as written by the MR batch
  harness, or a path to such an `.rds` file.

- estimators:

  Character vector selecting estimators, from `"oracle"`, `"naive"`,
  `"FB"`, `"MR"`. Default `NULL` takes all present, dropping `"FB"` when
  the run carried no full bootstrap (`meta$nb_boots` absent or below 1).

- blocks:

  Character vector of blocks, from `"H"` and `"Hc"`. Default both.

- digits:

  Integer or `NULL`. `NULL` (default) returns unrounded values; a
  formatter should round, not the computation.

## Value

An object of class `c("fs_mr_oc", "list")`:

- `estimation`:

  Data frame, one row per block and estimator: `n`, `avg`, `sd_emp`,
  `se_hat`, the four `bias_*`, the four `cov_*`, and `width`.

- `identification`:

  Data frame with one row per convention: detection rate, mean subgroup
  size, and mean sens/spec/ppv/npv.

- `targets`:

  The oriented target values used.

- `meta`:

  Key run metadata, echoed.

## Orientation

Estimator columns (`or_`, `nv_`, `fb_`, `mr_`) are stored on the
ORIENTED scale, where positive means harm. `betaHhat_H` and
`betaHhat_Hc` are stored on the RAW outcome scale. The bridge is
`orient = if (adverse_outcome) 1 else -1`, applied here as the
simulation harness applies it. Getting this backwards silently negates
every bias.

## Four targets, deliberately

Every estimator is scored against every target in both blocks, so that
disagreement is visible rather than assumed away:

- `beta`:

  Per-replicate exact \\\beta(\hat H)\\ – the conditional estimand at
  the realised region.

- `oracle`:

  Per-replicate refit on the true region in that trial.

- `struct`:

  Scalar structural potential-outcome effect in the true region; zero
  Monte Carlo error.

- `marg`:

  Scalar fitted effect in the true region on one realised draw; carries
  sampling noise.

Bias is ABSOLUTE, in outcome units. On an identity scale a complement
effect can sit near zero, so a relative bias would explode.

## Two classification conventions

Both are returned, because the project contains both and they disagree
whenever detection is not near certain:

- `conditional`:

  Averaged over DETECTED replicates only. This is the detected-only
  summary.

- `unconditional`:

  Averaged over ALL replicates. A replicate with no identified subgroup
  has everyone in its estimated complement, so it is scored as sens = 0,
  ppv = 0, spec = 1 and npv = \|Hc\|/N, the share of that replicate's
  subjects truly in the complement (Leon et al. 2024, Section 3),
  computed as `1 - n_true / meta$n_sample`. Under the null this is 1.
  These are the spec and npv values
  [`run_simulation_analysis()`](https://larry-leon.github.io/forestsearch/reference/run_simulation_analysis.md)
  stores for such a replicate. The unconditional npv is `NA` when the
  payload has no numeric `n_true` column or no usable `meta$n_sample`,
  or when `n_true` is missing for any non-detected replicate (for
  example a replicate that failed before its data were simulated).

With detection at 999/1000 the two agree to three decimals; in a null
cell at 0.930 they do not. Quote which one you mean.

## See also

[`fs_dgm_scale`](https://larry-leon.github.io/forestsearch/reference/fs_dgm_scale.md)

## Examples

``` r
if (FALSE) { # \dontrun{
oc <- fs_mr_oc_summary("fs_maxeffCons_mr_md40_knoise0_n500_res_1_1000.rds")
oc$estimation[oc$estimation$block == "H", c("estimator", "avg", "bias_beta")]
oc$identification
} # }
```
