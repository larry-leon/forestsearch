# Invert the family declaration rate for a target

Finds the threshold – `c1` at fixed `c2`, or `c2` at fixed `c1` – at
which the family declaration rate equals `target`, on one fixed draw
set: "the `c1` giving 80\\ read as a type-I error under a null DGM.

## Usage

``` r
fs_oc_invert(
  dgm = NULL,
  forestsearch_args = list(),
  n,
  target,
  solve_for = c("c1", "c2"),
  c1 = NULL,
  c2 = NULL,
  family = NULL,
  consistency_method = c("resample", "split"),
  pconsistency = NULL,
  draws = 2e+05,
  seed = NULL,
  tol = 0.001,
  ...
)
```

## Arguments

- dgm:

  An object of class `"glm_dgm"`; used only when `family = NULL`.

- forestsearch_args:

  Named list of
  [`forestsearch`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
  arguments. Supplies the default `c1` (`effect.threshold`) and `c2`
  (`consistency.threshold`) and, when `family = NULL`, the cut
  specification and floors for
  [`fs_oc_family_enumerate`](https://larry-leon.github.io/forestsearch/reference/fs_oc_family_enumerate.md).

- n:

  Single trial size.

- target:

  Numeric in (0, 1); may be a vector when `solve_for = "c1"`. Target
  family declaration rate(s).

- solve_for:

  `"c1"` (at fixed `c2`) or `"c2"` (at fixed `c1`).

- c1, c2:

  The fixed threshold (the one not solved for; the solved-for argument
  is ignored). Defaults from `forestsearch_args`.

- family:

  `NULL` (enumerate per `n`) or an `fs_oc_family` object, accepted only
  when its `n` equals the single `n` requested – a family built at one
  `n` has the wrong prevalence floor for another.

- consistency_method:

  One gate.

- pconsistency:

  Numeric in (0, 1). Consistency-*rate* threshold for the `"resample"`
  gate – the `p*` of
  [`forestsearch`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)'s
  `pconsistency.threshold`: the fraction of splits that must clear `c2`,
  not an effect threshold. Default
  `forestsearch_args$pconsistency.threshold`, itself defaulting to
  [`forestsearch()`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)'s
  default (0.90); an explicit value overrides. Ignored by `"split"`.

- draws:

  Integer. Number of Monte-Carlo draws.

- seed:

  Integer or `NULL`. Passed to
  [`set.seed()`](https://rdrr.io/r/base/Random.html) before the draws
  when non-`NULL`.

- tol:

  Bisection tolerance on the threshold scale (`solve_for = "c2"` only).

- ...:

  Passed to
  [`fs_oc_family_enumerate`](https://larry-leon.github.io/forestsearch/reference/fs_oc_family_enumerate.md)
  when `family = NULL`.

## Value

For a single `target`, an object of class `c("fs_oc_invert", "list")`:
`value` (the threshold, or `NA`), `achieved` (rate at `value`),
`achieved_se`, `target`, `ceiling`, `ceiling_se`, `binding` (which
threshold sets the ceiling), `attainable`, `solve_for`, `fixed`, `k`
(the order-statistic rank, `solve_for = "c1"`), `next_step_rate`,
`iterations` (bisections, `solve_for = "c2"`), `settings`. For a vector
`target`, a list of such objects (class `fs_oc_invert_list`) with a
`table` attribute.

## Details

**Solving for `c1`: an order statistic.** At fixed `(n, gate, c2)`
neither the eligible set nor the winner depends on `c1`: with
`eligible_g = (Bhat_g - c2 >= z_p * se_g)` (`"resample"`) or
`(W1_g >= c2) & (W2_g >= c2)` (`"split"`), each draw's winner is
`w = argmax_{eligible} Bhat_g` and its statistic `T = Bhat_w` (`-Inf`
when nothing is eligible). Declaration at `c1` is exactly `T >= c1`, and
`w` is the maxeffCons winner whenever anything declares. So the
reduction is computed once (`.fs_oc_reduce()`), the declaration rate at
any `c1` is `mean(T >= c1)`, and the largest `c1` attaining a target is
the `k`-th largest `T` with `k = ceiling(target * draws)`: an order
statistic, no search. The **ceiling** is `mean(T > -Inf)`, the fraction
of draws with any eligible member, set by the consistency screen alone;
a `target` above it is unattainable and returns `NA` with the ceiling
and the binding threshold named – no extrapolation. `target` may be a
vector: one reduction serves every target.

**Solving for `c2`** at fixed `c1`: the eligible set moves with `c2`, so
the rate is bracketed from the draws' range and bisected to `tol` on
fixed draws (a monotone non-increasing step function); the returned
value is the largest threshold whose rate is still at least `target`.

Every result reports the achieved rate and its MC SE, so the resolution
the draw count supports is visible.

## See also

[`fs_oc_grid`](https://larry-leon.github.io/forestsearch/reference/fs_oc_grid.md)

## Examples

``` r
piQ <- 0.34
fam <- structure(list(
  lab = c("Q", "P", "D"), Pg = c(piQ, 0.45, 0.31),
  PQg = c(1, 0.28 / 0.45, 1), sens_g = c(1, 0.28 / piQ, 0.31 / piQ),
  spec_g = c(1, 1 - 0.17 / (1 - piQ), 1),
  ovl = matrix(c(piQ, 0.28, 0.31, 0.28, 0.45, 0.28 * 0.31 / piQ,
                 0.31, 0.28 * 0.31 / piQ, 0.31), 3, 3),
  M = 3L, PQ = piQ, n = 500), class = c("fs_oc_family", "list"))
fam$beta_g <- 26 + 14 * fam$PQg
fam$se_g   <- 13.7 * sqrt(2) * sqrt(piQ / fam$Pg)
fs_oc_invert(family = fam, n = 500, target = 0.5, solve_for = "c1", c2 = 10,
             consistency_method = "resample", draws = 2e4, seed = 1)
#> Declaration-rate inversion (fs_oc_invert)
#>   solve for c1 at fixed c2 = 10; gate = resample; n = 500; M = 3; draws = 20,000
#>   ceiling (threshold -> -Inf): 0.5890 (MC SE 0.0035), set by c2 + z_p * se_g (resample consistency screen)
#>   target 0.500: c1 = 44.6017, achieved rate 0.5000 (MC SE 0.0035), next step 0.5000; order statistic k = 10000 of 20,000
```
