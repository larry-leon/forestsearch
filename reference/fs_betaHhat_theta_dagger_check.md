# `theta-dagger` at the true subgroup flag, on the scoring frame

The marginal target at the DGM's own harm flag, computed on the same
frame `beta(Hhat)` is scored on. Use it as a sanity gate: it should
reproduce the DGM's own subgroup effects.

## Usage

``` r
fs_betaHhat_theta_dagger_check(
  frame,
  outcome_type = c("survival", "binary", "continuous", "count"),
  harm.name = "flag_harm",
  outcome.name = "y_sim",
  event.name = "event_sim",
  treat.name = "treat_sim",
  effect_measure = NULL
)
```

## Arguments

- frame:

  Data frame. The evaluation population, normally from
  [`fs_build_eval_frame()`](https://larry-leon.github.io/forestsearch/reference/fs_build_eval_frame.md).

- outcome_type:

  Character. One of `"survival"`, `"binary"`, `"continuous"`, `"count"`.

- harm.name:

  Character. The true-subgroup flag column; `== 1L` is in.

- outcome.name, event.name, treat.name:

  Character. Column names in `frame`. `event.name` is survival only.

- effect_measure:

  Character. Required for `"continuous"` and `"count"`; ignored for
  `"survival"` and `"binary"`, which keep their fitted targets.

## Value

A named numeric vector, `thetaDagger_H` and `thetaDagger_Hc`.

## Details

For `"continuous"` and `"count"` this is an **exact identity**, not
agreement to Monte Carlo error – the frame is the super-population and
the arithmetic is the same
[`compute_aor()`](https://larry-leon.github.io/forestsearch/reference/compute_aor.md)
dispatch the DGM used. A tolerance there would hide a real defect.

This is a thin dispatch onto the same per-family effect used for every
region; it introduces no arithmetic of its own.
