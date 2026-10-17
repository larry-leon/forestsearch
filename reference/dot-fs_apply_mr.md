# Apply multiplier resampling with forestsearch's defaulting rules.

Thin wrapper around
[`fs_mr_inference()`](https://larry-leon.github.io/forestsearch/reference/fs_mr_inference.md)
shared by the consistency, DINA, and GRF branches so the call (and
`mr_inference_args` overrides) stays identical across methods. Wrapped
in tryCatch so a reconstruction failure yields `NULL` (the branch's
normal output is unaffected) rather than aborting.

## Usage

``` r
.fs_apply_mr(
  df,
  candidates,
  selected_members,
  spec,
  admission,
  effect_neighborhood,
  reselection_default,
  selection_rule_default = "neighborhood",
  mr_inference_args = list(),
  seedit = NULL
)
```

## Arguments

- reselection_default:

  Method-appropriate default re-selection rule (`"maxcons"` for
  consistency;
  [`.fs_mr_reselection_from_focus()`](https://larry-leon.github.io/forestsearch/reference/dot-fs_mr_reselection_from_focus.md)
  otherwise).

## Forwarded argument set

The wrapper forwards the *full*
[`fs_mr_inference()`](https://larry-leon.github.io/forestsearch/reference/fs_mr_inference.md)
argument set, so the DINA and GRF branches are controllable through
`mr_inference_args` exactly as the consistency branch is through its own
direct call. Before this, ten arguments were dropped here –
`return_reselection`, `field_R_out`, `field_R_in`, `field_uniform`,
`field_M_cap`, `field_complement`, `field_decompose`,
`field_scale_complement`, `ij_residual`, `field_recovery` – which left
`field_decompose` and `field_recovery` **unreachable** on those engines
and the other three **inert while their defaults happened to coincide**
with the intended values: a run's meta recorded them as set while they
controlled nothing.

**Defaults are unchanged on every branch.** Each newly forwarded
argument defaults to
[`fs_mr_inference()`](https://larry-leon.github.io/forestsearch/reference/fs_mr_inference.md)'s
own formal default, read from `formals(fs_mr_inference)` at call time
rather than restated here, so the two cannot drift and a call that asks
for nothing is byte-identical to the same call before the change.
[`match.arg()`](https://rdrr.io/r/base/match.arg.html)-style defaults
are forwarded as their whole candidate vector, which
[`fs_mr_inference()`](https://larry-leon.github.io/forestsearch/reference/fs_mr_inference.md)
resolves exactly as it would have resolved a missing argument.
`ci_method` keeps this wrapper's own `"ij"` default, which is
deliberately *not* the consistency branch's `"field"`; changing it is a
separate decision.
