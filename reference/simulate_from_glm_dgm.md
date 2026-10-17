# Simulate Trial Data from a GLM Data Generating Mechanism

Generates a simulated clinical trial dataset by sampling from the
super-population in a `"glm_dgm"` object. Assigns treatment and
generates outcomes from the individual-level potential outcomes stored
in the DGM.

## Usage

``` r
simulate_from_glm_dgm(
  dgm,
  n = NULL,
  rand_ratio = 1,
  seed = NULL,
  draw_treatment = TRUE,
  replace = TRUE,
  baseline = c("resample", "fixed")
)
```

## Arguments

- dgm:

  An object of class `"glm_dgm"` created by
  [`generate_glm_dgm`](https://larry-leon.github.io/forestsearch/reference/generate_glm_dgm.md).

- n:

  Integer. Sample size. If `NULL` (default), uses the entire
  super-population without resampling.

- rand_ratio:

  Numeric. Treatment:control randomisation ratio. Default `1` (1:1
  allocation).

- seed:

  Integer. Random seed. Default `NULL` (no explicit seed; use the
  calling environment's RNG state).

- draw_treatment:

  Logical. If `TRUE` (default), randomly assigns treatment. If `FALSE`,
  retains the treatment column from the super-population.

- replace:

  Logical sampling scheme for drawing the `n` subjects from
  `dgm$df_super`. `TRUE` (default) samples with replacement
  (bootstrap-style trial sampling). `FALSE` samples without replacement
  and so requires `n <= nrow(dgm$df_super)`; treatment is still
  re-randomised and the outcome re-drawn. Used to draw the full pool
  once (`n = nrow(dgm$df_super)`, `replace = FALSE`) for the beta(Hhat)
  evaluation frame. Ignored on the `n = NULL` path.

- baseline:

  \\resample\\ (default; draw the trial panel from the super-population)
  or \\fixed\\ (use dgm\$df_source – the observed panel, every patient
  exactly once, observed order; requires a DGM built by a version that
  stores df_source; n must be NULL or nrow(df_source)).

## Value

A data frame with columns:

- `id`:

  Integer subject identifier.

- `y_sim`:

  Simulated outcome.

- `treat_sim`:

  Treatment indicator (0/1).

- `flag_harm`:

  True subgroup membership indicator.

Plus all covariate columns from the super-population.

For binary outcomes, `y_sim` is a 0/1 Bernoulli draw from the potential
outcome probabilities (`p0` or `p1`). There is no `event_sim` column –
for ForestSearch, set `event.name = "y_sim"`.

## See also

[`generate_glm_dgm`](https://larry-leon.github.io/forestsearch/reference/generate_glm_dgm.md),
[`run_simulation_analysis`](https://larry-leon.github.io/forestsearch/reference/run_simulation_analysis.md)

## Examples

``` r
if (FALSE) { # \dontrun{
dgm <- generate_glm_dgm(...)
df  <- simulate_from_glm_dgm(dgm, n = 500, seed = 42)
table(df$y_sim, df$treat_sim)
} # }
```
