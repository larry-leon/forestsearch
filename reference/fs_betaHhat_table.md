# Deduplicated `beta(Hhat)` targets over distinct realized rules

Scores each **distinct** realized rule once. Deduplication is by rule,
not by replicate: two replicates landing on the same rule score the same
target, and the target does not depend on the trial.

## Usage

``` r
fs_betaHhat_table(
  sg_defs,
  frame,
  focus,
  outcome_type = c("survival", "binary", "continuous", "count"),
  effect_measure = NULL,
  outcome.name = "y_sim",
  event.name = "event_sim",
  treat.name = "treat_sim",
  sg_def_structs = NULL
)
```

## Arguments

- sg_defs:

  Character vector of realized rules, one entry per replicate. `NA` or
  `""` marks a replicate where no subgroup was identified.

- frame:

  Data frame. The fixed evaluation population.

- focus, outcome_type, effect_measure, outcome.name, event.name,
  treat.name:

  Passed to
  [`fs_betaHhat_one()`](https://larry-leon.github.io/forestsearch/reference/fs_betaHhat_one.md).

- sg_def_structs:

  Optional named list of structured GRF definitions, keyed by the rule
  string, used in preference to string parsing.

## Value

A data frame keyed by `sg_def`, one row per distinct rule, with the
[`fs_betaHhat_one()`](https://larry-leon.github.io/forestsearch/reference/fs_betaHhat_one.md)
schema. Resolution counters are attached as attributes and retrievable
with
[`fs_betaHhat_counts()`](https://larry-leon.github.io/forestsearch/reference/fs_betaHhat_counts.md).

## Resolution accounting

Six counters are attached to the result and are the fix for targets
being dropped silently. `n_eff` reported beside a target is what makes a
reduced denominator visible; a coverage figure computed on an unknown
fraction of replicates is not interpretable, and the fraction has to
travel with it.

A seventh, `n_reps_undetected`, counts replicates with no rule at all,
so that `n_reps_resolved + n_reps_unresolved + n_reps_undetected` closes
against `n_reps_total`.
