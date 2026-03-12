# Simulate Survival Data from AFT Data Generating Mechanism

Generates simulated survival data from a previously created AFT data
generating mechanism (DGM). Samples from the super population and
generates survival times with specified censoring.

## Usage

``` r
simulate_from_dgm(
  dgm,
  n = NULL,
  rand_ratio = 1,
  entry_var = NULL,
  max_entry = 24,
  analysis_time = 48,
  cens_adjust = 0,
  draw_treatment = TRUE,
  seed = NULL,
  strata_rand = NULL,
  hrz_crit = NULL,
  keep_rand = FALSE,
  time_eos = NULL
)
```

## Arguments

- dgm:

  An object of class `"aft_dgm_flex"` created by
  [`generate_aft_dgm_flex`](https://larry-leon.github.io/forestsearch/reference/generate_aft_dgm_flex.md).

- n:

  Integer specifying the sample size. If `NULL` (default), uses the
  entire super population without sampling.

- rand_ratio:

  Numeric randomisation ratio (treatment:control). Default `1` (1:1
  allocation).

- entry_var:

  Character string naming an entry-time variable in the super
  population. If `NULL`, entry times are drawn as
  `Uniform(0, max_entry)`. Default `NULL`.

- max_entry:

  Numeric maximum entry time for staggered entry simulation. Only used
  when `entry_var = NULL`. Default `24`.

- analysis_time:

  Numeric calendar time of analysis. Follow-up is
  `analysis_time - entry_time`. Must be on the same time scale as the
  DGM (i.e. the same units as `outcome_var` passed to
  `generate_aft_dgm_flex`). Default `48`.

- cens_adjust:

  Numeric log-scale adjustment to censoring distribution. Positive
  values increase censoring times; negative values decrease them.
  Default `0` (no adjustment).

- draw_treatment:

  Logical. If `TRUE` (default), reassigns treatment according to
  `rand_ratio`. If `FALSE`, retains original treatment assignments from
  the super population.

- seed:

  Integer random seed. Default `NULL`.

- strata_rand:

  Character string naming a column in the sampled data for
  within-stratum balanced treatment allocation. If `NULL`, marginal
  allocation is used. Default `NULL`.

- hrz_crit:

  Numeric log-HR threshold. If supplied, a column `hrz_flag` is added
  marking subjects with `lin_pred_1 - lin_pred_0 >= hrz_crit`. Default
  `NULL`.

- keep_rand:

  Logical. If `TRUE`, appends a `rand_order` column preserving the
  randomisation sequence. Default `FALSE`.

- time_eos:

  Numeric secondary administrative censoring cutoff (end-of-study time
  on the DGM scale). Applied after `follow_up` censoring. Default
  `NULL`.

## Value

A `data.frame` with columns:

- `id`:

  Subject identifier.

- `treat`:

  Original treatment from super population.

- `treat_sim`:

  Simulated treatment assignment.

- `flag_harm`:

  Subgroup indicator (1 = all subgroup conditions met).

- `z_*`:

  Covariate values.

- `lin_pred_1`, `lin_pred_0`:

  Counterfactual log-time linear predictors.

- `y_sim`:

  Observed survival time (`min(T, C)`).

- `event_sim`:

  Event indicator (1 = event, 0 = censored).

- `t_true`:

  Latent true survival time (pre-censoring).

- `c_time`:

  Effective censoring time (post admin-censoring).

- `hrz_flag`:

  (Optional) Individual harm-zone indicator.

- `rand_order`:

  (Optional) Randomisation sequence index.

## Details

### Time-scale consistency

All time parameters (`analysis_time`, `max_entry`, `time_eos`) must be
expressed in the same units as `outcome_var` supplied to
[`generate_aft_dgm_flex()`](https://larry-leon.github.io/forestsearch/reference/generate_aft_dgm_flex.md).
A common error is building the DGM on days (e.g. `rfstime`) and then
passing `analysis_time` in months, which causes follow-up windows far
shorter than the DGM event-time scale and produces universal
administrative censoring (`event_sim = 0` for all subjects).

Verify with: `exp(dgm$model_params$mu)` — the implied median event time
should be plausible given your `analysis_time`.

### n = NULL path

When `n = NULL` the entire super population is used as-is, with no
staggered entry and no administrative censoring (`follow_up = Inf`).
Treatment assignments and linear predictors already stored in
`dgm$df_super` are retained unchanged.

### Censoring adjustment

`cens_adjust` shifts the log-scale location parameter of the censoring
distribution:

- `cens_adjust = log(2)` doubles expected censoring times.

- `cens_adjust = log(0.5)` halves expected censoring times.

## See also

[`generate_aft_dgm_flex`](https://larry-leon.github.io/forestsearch/reference/generate_aft_dgm_flex.md),
[`check_censoring_dgm`](https://larry-leon.github.io/forestsearch/reference/check_censoring_dgm.md)

## Examples

``` r
if (FALSE) { # \dontrun{
dgm <- setup_gbsg_dgm(model = "null", verbose = FALSE)
sim_data <- simulate_from_dgm(dgm, n = 200, seed = 42)
dim(sim_data)
head(sim_data[, c("y_sim", "event_sim", "treat_sim")])
} # }
```
