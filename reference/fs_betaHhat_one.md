# Conditional estimand `beta(Hhat)` for one realized rule

Computes the population value of the standard within-subgroup analysis
on a realized selection region and its complement, on a fixed evaluation
frame. This is the conditional estimand of León and Anderson (2026,
Theorem 2), evaluated at the rule a detector actually returned – not the
marginal true-subgroup target.

## Usage

``` r
fs_betaHhat_one(
  rule,
  frame,
  focus,
  outcome_type = c("survival", "binary", "continuous", "count"),
  effect_measure = NULL,
  outcome.name = "y_sim",
  event.name = "event_sim",
  treat.name = "treat_sim",
  sg_def_struct = NULL,
  method = c("frame", "exact"),
  dgm = NULL,
  analysis_time = 84,
  max_entry = 24,
  cens_adjust = log(1.5)
)
```

## Arguments

- rule:

  Character. The realized rule: a named `sg.harm` vector, the
  `" & "`-joined `sg_def` string, or a GRF disjunction
  `"(a & b) | (c & d)"`. `NULL` or `NA` means no subgroup was
  identified.

- frame:

  Data frame. The fixed evaluation population.

- focus:

  Character, required, no default. `"harm"` reads `treat.recommend == 0`
  as in-region; `"benefit"` inverts.

- outcome_type:

  Character. One of `"survival"`, `"binary"`, `"continuous"`, `"count"`.

- effect_measure:

  Character. Passed to
  [`compute_aor()`](https://larry-leon.github.io/forestsearch/reference/compute_aor.md)
  for the continuous and count families; ignored for survival and
  binary, which keep their fitted targets.

- outcome.name, event.name, treat.name:

  Character. Column names in `frame`. `event.name` is used only by the
  survival family.

- sg_def_struct:

  Optional structured GRF subgroup definition. When supplied it takes
  precedence over `rule` and membership is computed by
  `.grf_evaluate_subgroup()`, which is the package's own answer for GRF
  membership.

- method:

  Character. `"frame"` (default) scores the rule on `frame` by the
  per-family fit, exactly as before this argument existed. `"exact"`
  (survival only) computes the same population quantity from the
  analytic layer; see the section *The exact route*.

- dgm:

  The DGM `frame` was drawn from (an
  [`inject_mrct_structure()`](https://larry-leon.github.io/forestsearch/reference/inject_mrct_structure.md)
  build). Required by `method = "exact"`; ignored by `"frame"`.

- analysis_time, max_entry, cens_adjust:

  The analysis time, the entry maximum and the censoring adjustment of
  the evaluation draw. Used by `method = "exact"` only; the defaults are
  those of the draw
  [`fs_build_eval_frame()`](https://larry-leon.github.io/forestsearch/reference/fs_build_eval_frame.md)
  makes (`analysis_time = 84`, `cens_adjust = log(1.5)`, and
  [`simulate_from_dgm()`](https://larry-leon.github.io/forestsearch/reference/simulate_from_dgm.md)'s
  default `max_entry = 24`, which
  [`fs_build_eval_frame()`](https://larry-leon.github.io/forestsearch/reference/fs_build_eval_frame.md)
  does not override). Pass the values the frame was drawn with.

## Value

A one-row data frame with columns `betaHhat_H`, `betaHhat_Hc`,
`nH_eval`, `nHc_eval`, `status` and `missing_cols`. The column names are
scale-agnostic and identical across families so downstream
`paste0("betaHhat_", suffix)` scoring is shared verbatim.

## Details

Membership is resolved once, by the internal `.fs_resolve_membership()`,
and the same resolution serves every outcome family: membership is a
property of the rule and the frame, not of the effect measure.

## The partition invariant

`nH_eval + nHc_eval` must equal `nrow(frame)` on every call where the
rule resolves. Hhat and its complement partition the frame by
construction, so a sum below `nrow(frame)` means some rows fell into
neither side – which happens when a membership vector carries `NA`. That
is a hard error, not a warning: if membership is incoherent, every
number downstream of it is meaningless.

## Unresolvable rules

A rule naming a column the frame lacks has no membership on that frame.
The record is all-`NA` with `status = "unresolved"` and `missing_cols`
naming the columns, rather than a partially-resolved region. Silently
dropping the offending clause returns a finite value for the wrong
region.

## The exact route

With `method = "exact"` the survival target is computed without a draw:
[`ahr_population()`](https://larry-leon.github.io/forestsearch/reference/ahr_population.md)`(dgm, analysis_time, max_entry)`
once, then `exp(beta*)` from
[`ahr_limit()`](https://larry-leon.github.io/forestsearch/reference/ahr_limit.md)`(pop, set = <membership>, psi = pop$psi)`
for the realized region and for its complement, where `beta*` is the
Struthers-Kalbfleisch limit of the two-sample Cox estimate (1:1
randomisation, entry U(0, `max_entry`), analysis at `analysis_time`).
Membership is resolved on `dgm$df_super` by the same internal
`.fs_resolve_membership()` (the
[`get_dfpred()`](https://larry-leon.github.io/forestsearch/reference/get_dfpred.md)
route) the frame route uses on `frame`, so a rule the frame route cannot
evaluate is `NA` on both routes, and the counts `nH_eval` / `nHc_eval`
are over `dgm$df_super`. An empty region is `NA`, as on the frame route.

The two routes estimate the same population quantity: the frame route
fits [`survival::coxph()`](https://rdrr.io/pkg/survival/man/coxph.html)
on one
[`simulate_from_dgm()`](https://larry-leon.github.io/forestsearch/reference/simulate_from_dgm.md)
draw of the whole super-population and so carries the Monte Carlo error
of that single evaluation draw; the exact route has none. The frame
route's small-region guard (fewer than five events or a single arm in
the draw) has no counterpart: the limit exists for any non-empty region.

The exact route applies when the censoring is administrative only: the
analytic layer has no random-censoring factor. It stops, rather than
approximating, unless the DGM's random censoring time falls before
`analysis_time` with probability below `.Machine$double.eps` for every
patient in both arms under `cens_adjust`. In
[`simulate_from_dgm()`](https://larry-leon.github.io/forestsearch/reference/simulate_from_dgm.md)
`cens_adjust = 0` means *no adjustment* of the DGM's censoring model,
not no random censoring, and the frame route's default `log(1.5)`
lengthens the random censoring times but keeps them; `cens_adjust = 50`
(the administrative-only convention of the MRCT campaigns) or `Inf`
removes it. Requires a DGM made by
[`inject_mrct_structure()`](https://larry-leon.github.io/forestsearch/reference/inject_mrct_structure.md),
as
[`ahr_population()`](https://larry-leon.github.io/forestsearch/reference/ahr_population.md)
does.
