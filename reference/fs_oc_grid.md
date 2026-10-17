# Sweep the declaration thresholds on one draw set per trial size

Evaluates every
[`fs_oc_predict()`](https://larry-leon.github.io/forestsearch/reference/fs_oc_predict.md)
quantity over the full crossing of `n`, `c1` and `c2`, enumerating the
family and drawing the candidates' fitted effects **once per `n`** (and
per gate) and sweeping the thresholds against those draws. Thresholds
enter only the gate, so the sweep costs arithmetic, not draws.

## Usage

``` r
fs_oc_grid(
  dgm = NULL,
  forestsearch_args = list(),
  n,
  c1,
  c2,
  family = NULL,
  consistency_method = c("resample", "split"),
  pconsistency = NULL,
  draws = 2e+05,
  block = 50000,
  seed = NULL,
  verbose = FALSE,
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

  Numeric vector of trial sizes.

- c1, c2:

  Numeric vectors of screening and consistency floors; the grid is their
  full crossing with `n`.

- family:

  `NULL` (enumerate per `n`) or an `fs_oc_family` object, accepted only
  when its `n` equals the single `n` requested – a family built at one
  `n` has the wrong prevalence floor for another.

- consistency_method:

  `"resample"`, `"split"`, or both.

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

- block:

  Integer or `Inf`. Draw-block size; `Inf` means one block (exact).
  Default `5e4`.

- seed:

  Integer or `NULL`. Passed to
  [`set.seed()`](https://rdrr.io/r/base/Random.html) before the draws
  when non-`NULL`.

- verbose:

  Logical. Report per-`n` timings.

- ...:

  Passed to
  [`fs_oc_family_enumerate`](https://larry-leon.github.io/forestsearch/reference/fs_oc_family_enumerate.md)
  when `family = NULL` (`max_M`, `verbose`).

## Value

An object of class `c("fs_oc_grid", "list")`:

- `table`:

  Data frame, one row per `(n, gate, c1, c2)`:
  `n, consistency_method, pconsistency, c1, c2, M, draws, block, seed, det_rate, det_rate_se, EnH, Esens, Espec, Eppv, Enpv, EbetaH, Enaive_bias, mass_below`.

- `results`:

  List, parallel to the rows of `table`, of the full per-point objects
  (class `fs_oc_predict`, with `P1`, `p_sel`, `sel_c` and their MC SEs).

- `families`:

  Per `n`: `M`, the floor, the stage counts and the family.

- `timing`:

  Per `n` and gate: seconds for enumeration, drawing and sweeping.

- `settings`:

  The call's settings.

## Details

For each `n` the family is enumerated with
[`fs_oc_family_enumerate`](https://larry-leon.github.io/forestsearch/reference/fs_oc_family_enumerate.md)
(or the supplied `family` is used when its `n` matches), the draw set is
generated under `seed` (the same seed for every `n`: common random
numbers across trial sizes), and every `(c1, c2)` is evaluated on it.
With `block = Inf` the whole draw set is held and the computation is
exactly
[`fs_oc_predict`](https://larry-leon.github.io/forestsearch/reference/fs_oc_predict.md)'s
– one grid point is
[`identical()`](https://rdrr.io/r/base/identical.html) to a call of
[`fs_oc_predict()`](https://larry-leon.github.io/forestsearch/reference/fs_oc_predict.md)
at the same settings and seed. With a finite `block` the draws are
generated in row blocks of that size and the quantities accumulated, so
memory is O(block x M) instead of O(draws x M); blocked results agree
with the one-block results to Monte-Carlo precision (the RNG stream is
laid out differently and proportions are formed as count / draws), not
bit-for-bit.

## See also

[`fs_oc_predict`](https://larry-leon.github.io/forestsearch/reference/fs_oc_predict.md),
[`fs_oc_invert`](https://larry-leon.github.io/forestsearch/reference/fs_oc_invert.md)

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
g <- fs_oc_grid(family = fam, n = 500, c1 = c(20, 30, 40), c2 = 10,
                consistency_method = "resample", draws = 2e4, seed = 1)
g
#> Threshold sweep on one draw set per n (fs_oc_grid)
#>   n: 500 | gates: resample | c1: 3 values | c2: 1 values | draws = 20,000, block = 50000, seed = 1
#>   3 grid rows; table columns: n, consistency_method, pconsistency, c1, c2, M, draws, block, seed, det_rate, det_rate_se, EnH, Esens, Espec, Eppv, Enpv, EbetaH, Enaive_bias, mass_below
```
