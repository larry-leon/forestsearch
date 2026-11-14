# Attach `beta(Hhat)` targets to a results frame

The one call each engine adds to its `run_cell()`, immediately before
the bundle is assembled. Runs in the main process, after the replicate
loop, so the evaluation frame stays in one place.

## Usage

``` r
fs_attach_betaHhat(
  results,
  frame,
  focus,
  outcome_type = c("survival", "binary", "continuous", "count"),
  effect_measure = NULL,
  outcome.name = "y_sim",
  event.name = "event_sim",
  treat.name = "treat_sim",
  sg_def_structs = NULL,
  method = c("frame", "exact"),
  dgm = NULL,
  analysis_time = 84,
  max_entry = 24,
  cens_adjust = log(1.5)
)
```

## Arguments

- results:

  Data frame carrying an `sg_def` column.

- frame:

  Data frame. The fixed evaluation population.

- focus, outcome_type, effect_measure, outcome.name, event.name,
  treat.name:

  Passed to
  [`fs_betaHhat_one()`](https://larry-leon.github.io/forestsearch/reference/fs_betaHhat_one.md).

- sg_def_structs:

  Optional named list of structured GRF definitions.

- method, dgm, analysis_time, max_entry, cens_adjust:

  Passed to
  [`fs_betaHhat_one()`](https://larry-leon.github.io/forestsearch/reference/fs_betaHhat_one.md).
  `method = "frame"` (default) is the route above, unchanged.
  `method = "exact"` (survival only, administrative censoring only)
  replaces the fit on `frame` by the analytic limit on `dgm`'s
  super-population, built once; the returned columns are the same. See
  *The exact route* in
  [`fs_betaHhat_one()`](https://larry-leon.github.io/forestsearch/reference/fs_betaHhat_one.md).
  Under `"exact"`, `frame` is not used and the counts are over
  `dgm$df_super`.

## Value

`results` with `betaHhat_H`, `betaHhat_Hc`, `betaHhat_status` and
`nH_eval` / `nHc_eval` columns added, and counters attached.

## Replicates with no realized rule

An undetected replicate is not a missing measurement – it is a run in
which the whole population is the complement. It therefore receives the
no-subgroup record
[`fs_betaHhat_one()`](https://larry-leon.github.io/forestsearch/reference/fs_betaHhat_one.md)
already returns: `nH_eval = 0`, `nHc_eval = nrow(frame)`, `betaHhat_Hc`
the ITT effect, `betaHhat_H` `NA` (an empty region has no target), and
`status = "ok"`.

This is the partition invariant reaching this layer:
`nH_eval + nHc_eval` now equals `nrow(frame)` on **every** row that is
not `"unresolved"`, not only on rows carrying a rule. Before this change
such rows were all-`NA`, which silently excluded them from any consumer
keying on a finite target. Their count is still reported separately as
`n_reps_undetected`.
