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
  sg_def_struct = NULL
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
