# SPEC: `fs_build_eval_frame()` and `fs_betaHhat_theta_dagger_check()`

Two exported functions, each dispatching on `outcome_type` — mirroring
`.fs_region_effect()`, not six per-family exports. Both live in
`R/betaHhat_truth.R`. Export pattern `@keywords internal` + `@export`,
matching the existing six. Roxygen markdown is on: literal `% < > &`, plain
`@section` titles.

## `fs_build_eval_frame(dgm, outcome_type, eval_seed = 20260628L, analysis_time = 84, cens_adjust = log(1.5), n_eval = NULL)`

The frame `beta(Hhat)` is scored on, per family:

- `"survival"`: `simulate_from_dgm(dgm, n = nrow(dgm$df_super),
  replace = FALSE, analysis_time = analysis_time, cens_adjust = cens_adjust,
  seed = eval_seed)`. The `n_eval` trap carries over verbatim from the shim:
  non-NULL hard-errors with the full-pool message. Defaults must equal the
  shim's exactly — including `cens_adjust = log(1.5)`, the default restored
  at STOP 4 of the consolidation.
- `"binary"`: `simulate_from_glm_dgm(dgm, n = nrow(dgm$df_super),
  replace = FALSE, seed = eval_seed)`. `analysis_time`/`cens_adjust`/`n_eval`
  are rejected with an error if supplied non-default — silently ignoring a
  survival argument on the binary path is how conventions drift.
- `"continuous"`: returns `dgm$df_super` unchanged. **Decision, vetoable at
  STOP 1:** the MD target is an exact finite mean over `df_super`, so the
  scoring frame *is* the super-population; returning it gives a uniform
  producer contract (`fs_attach_betaHhat(results, fs_build_eval_frame(...))`)
  across families. The roxygen must state plainly that no simulation occurs
  and the target on this frame carries zero Monte Carlo error, so nobody
  reads uniformity as sameness. `eval_seed` is accepted and ignored (with a
  documented note), so generic harness code need not branch.

## `fs_betaHhat_theta_dagger_check(frame, outcome_type, harm.name = "flag_harm", outcome.name = "y_sim", event.name = "event_sim", treat.name = "treat_sim", effect_measure = NULL)`

Theta-dagger at the TRUE flag on the same frame:
`.fs_region_effect()` at `frame[[harm.name]] == 1L` and its complement,
returning `c(thetaDagger_H, thetaDagger_Hc)`. No new arithmetic — this is a
thin dispatch, and for continuous it is the exact identity gate (agreement
with the DGM's own `effect_Q`/`effect_Qc` to `0`, not to tolerance).

## Contract tests (append to the existing suite)

- T11 — frame identity: on the same seed, the survival frame is bitwise
  identical to the shim's `build_eval_frame()`, the binary frame to
  `build_eval_frame_glm()`, and the continuous return is `identical()` to
  `dgm$df_super`.
- T12 — the `n_eval` trap fires on survival; a survival-only argument
  supplied on the binary path errors.
- T13 — theta values: each family's check is `identical()` to the
  corresponding shim function's output on the same frame.
- T14 — negative control: the wrong `outcome_type` on the same frame must
  differ (or error), so T11/T13 cannot pass vacuously.

`R CMD check --as-cran` clean is the bar.
