# Sort Subgroups by Focus (post-consistency)

Sorts a data.table of subgroup results according to the specified focus.
For `"hrMaxSG"` and `"hrMinSG"`, candidates are first partitioned into a
*candidate inclusion set* (the "band"), then ranked within the band by
sample size. The `selection_rule` argument controls how the band is
defined:

- `"neighborhood"` – 1D effect-size band (default; legacy behaviour): a
  candidate is in-band iff its effect is within `effect_neighborhood` of
  the maximum.

- `"pareto"` – 2D Pareto-non-dominated set in (effect, N) space: a
  candidate is in-band iff no other candidate has both a larger effect
  and a larger N.

- `"both"` – intersection of the two: a candidate is in-band iff it
  satisfies the neighborhood test *and* is Pareto-non-dominated.

See Details.

## Usage

``` r
sort_subgroups(
  result_new,
  sg_focus,
  selection_rule = "neighborhood",
  effect_neighborhood = 0.1,
  effect_log_scale = FALSE
)
```

## Arguments

- result_new:

  A data.table of subgroup results with columns `Pcons`, `hr`, `N`, `K`.

- sg_focus:

  Character. Sorting focus. One of `"hr"`, `"maxSG"`, `"hrMaxSG"`,
  `"minSG"`, `"hrMinSG"`.

- selection_rule:

  Character. Rule defining the candidate inclusion set for `"hrMaxSG"` /
  `"hrMinSG"`. One of `"neighborhood"` (default; current behaviour),
  `"pareto"`, or `"both"`. Must be `"neighborhood"` for other `sg_focus`
  values.

- effect_neighborhood:

  Numeric in `[0, 1)`. Relative tolerance defining the effect-size
  neighborhood for `"hrMaxSG"` and `"hrMinSG"`. A candidate is in the
  neighborhood iff its (natural-scale) effect is at least
  `(1 - effect_neighborhood) * max(effect)`. Default `0.10` (i.e.,
  within 10\\ `selection_rule` `"neighborhood"` and `"both"`; ignored by
  `"pareto"` (and rejected if non-default with `"pareto"`).

- effect_log_scale:

  Logical. If `TRUE`, the `hr` column stores log-scale values (log-OR,
  log-IRR) and is exponentiated before applying the inclusion test.
  Default `FALSE` (Cox stores natural-scale HR).

## Value

A sorted data.table. The top row is the selected subgroup under
`sg_focus`; remaining rows are diagnostic.

## Details

**Sort keys by `sg_focus`:**

- `"hr"`:

  `(-Pcons, -hr, K)` – prefer high consistency and large effect.

- `"maxSG"`:

  `(-N, -Pcons, K)` – prefer large subgroups with high consistency.

- `"minSG"`:

  `(N, -Pcons, K)` – prefer small subgroups with high consistency.

- `"hrMaxSG"`:

  `(-in_band, -N, -Pcons, -hr, K)` – among candidates in the inclusion
  band defined by `selection_rule`, prefer the largest sample size.

- `"hrMinSG"`:

  `(-in_band, N, -Pcons, -hr, K)` – among candidates in the inclusion
  band defined by `selection_rule`, prefer the smallest sample size.

`"hrMaxSG"` and `"hrMinSG"` use a *lexicographic with candidate band*
rule: effect size is primary, but a bounded set of candidates near (or
trading off against) the maximum is accepted to optimize sample size.

**Selection rules:**

- `"neighborhood"`:

  1D \\\varepsilon\\-band on the effect scale:
  `in_band = (hr >= (1 - effect_neighborhood) * max(hr))`. Setting
  `effect_neighborhood = 0` reduces this to a strict max-effect filter.
  This is the legacy behaviour.

- `"pareto"`:

  2D Pareto-non-dominated set: a candidate is in-band iff no other
  candidate has both `hr_j >= hr_i` and `N_j >= N_i` with at least one
  strict inequality. Reuses the same dominance computation as
  [`compute_pareto_frontier()`](https://larry-leon.github.io/forestsearch/reference/compute_pareto_frontier.md).
  Removes the `effect_neighborhood` tuning parameter from the selection
  criterion.

- `"both"`:

  Intersection: in-band iff in the \\\varepsilon\\-band *and*
  Pareto-non-dominated. Strictest of the three.

## Examples

``` r
# \donttest{
library(data.table)
dt <- data.table(Pcons = c(0.92, 0.95, 0.88, 0.90),
                 hr    = c(2.5,  2.4,  2.3,  2.0),
                 N     = c(70,   100,  150,  200),
                 K     = c(1, 2, 1, 2))
# hrMaxSG with 10% neighborhood: HR 2.5, 2.4, 2.3 are in-band;
# among those, N = 150 is largest -> top row is N=150, hr=2.3.
sort_subgroups(dt, sg_focus = "hrMaxSG", effect_neighborhood = 0.10)
#>    Pcons    hr     N     K
#>    <num> <num> <num> <num>
#> 1:  0.88   2.3   150     1
#> 2:  0.95   2.4   100     2
#> 3:  0.92   2.5    70     1
#> 4:  0.90   2.0   200     2

# Same data with selection_rule = "pareto": the dominance set in
# (hr, N) is used instead of the 1D effect band.
sort_subgroups(dt, sg_focus = "hrMaxSG", selection_rule = "pareto")
#>    Pcons    hr     N     K
#>    <num> <num> <num> <num>
#> 1:  0.90   2.0   200     2
#> 2:  0.88   2.3   150     1
#> 3:  0.95   2.4   100     2
#> 4:  0.92   2.5    70     1
# }
```
