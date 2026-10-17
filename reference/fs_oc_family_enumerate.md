# Enumerate the candidate family on a DGM's population frame

Builds the family of candidate subgroups that
[`forestsearch`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
would search over, but on the DGM's super-population frame
`dgm$df_super` rather than on a sampled trial, so that every cut lands
at a population quantile and every prevalence, purity and overlap is a
population proportion. The result carries the nine quantities the
operating-characteristics prediction of
[`fs_oc_predict`](https://larry-leon.github.io/forestsearch/reference/fs_oc_predict.md)
consumes.

## Usage

``` r
fs_oc_family_enumerate(
  dgm,
  forestsearch_args,
  n,
  max_M = 2000L,
  verbose = FALSE
)
```

## Arguments

- dgm:

  An object of class `"glm_dgm"` from
  [`generate_glm_dgm`](https://larry-leon.github.io/forestsearch/reference/generate_glm_dgm.md),
  with the true-region indicator in `df_super$flag_harm`.

- forestsearch_args:

  Named list of
  [`forestsearch`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
  arguments. `confounders.name` is required; the cut, `maxk`, `minp`,
  `rmin`, `n.min` and `n.min.frac` entries are read, and any omitted
  entry takes
  [`forestsearch()`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)'s
  default.

- n:

  Integer. Trial size at which the family is evaluated: sets the size
  floor `n.min / n` and the `sqrt(1000 / n)` rescaling of the anchored
  standard errors.

- max_M:

  Integer. Size guard: if more than `max_M` candidates survive the
  floors, the function stops before allocating the M-by-M overlap matrix
  (the downstream symmetric root is O(M^3)).

- verbose:

  Logical. Print the enumeration counts at each stage.

## Value

An object of class `c("fs_oc_family", "list")` with elements

- `lab`:

  Character, length M: the rule of each candidate.

- `Pg`:

  Population prevalence `P(g)`.

- `PQg`:

  Purity `P(g & Q) / P(g)`.

- `beta_g`:

  Oriented mixture mean
  `s * (m_tau[Qc] + (m_tau[Q] - m_tau[Qc]) * PQg)`, with
  `s = sign(m_tau[Q])` the harm direction.

- `se_g`:

  Anchored standard error at `n`.

- `sens_g`:

  `P(g & Q) / P(Q)`.

- `spec_g`:

  `1 - P(g & Qc) / P(Qc)`.

- `ovl`:

  M-by-M matrix of `P(g_i & g_j)`.

- `M`:

  Number of candidates.

- `PQ`:

  `P(Q)`, the prevalence of the true region.

- `memb`:

  Logical matrix, `nrow(df_super)` by M, the population membership of
  each candidate.

- `null`:

  Logical: `TRUE` when the null branch was taken (Q empty).

- `orientation`:

  Alternative branch only: list with the harm-direction sign `s`, the
  signed region effects `m_tau_Q` and `m_tau_Qc`, and the oriented
  mixture coefficients `tauQc = s * m_tau_Qc` (may be negative for
  opposite-sign families) and `bint = s * (m_tau_Q - m_tau_Qc)`. Absent
  on the null branch.

- `scale`:

  The `fs_dgm_scale` object used.

- `n`:

  The trial size.

- `args_used`:

  The `forestsearch_args` entries consumed, with the resolved `n.min`.

- `cuts`:

  The cut expressions
  [`get_FSdata()`](https://larry-leon.github.io/forestsearch/reference/get_FSdata.md)
  produced.

- `counts`:

  Named integer vector: candidates at each stage.

## Details

**Cuts.** The cut-related entries of `forestsearch_args`
(`confounders.name`, `conf.cont_jcuts`, `cut_type`, `cont.cutoff`,
`conf.cont_medians`, `conf.cont_medians_force`, `conf_force`,
`defaultcut_names`, `exclude_cuts`, `collapse_cuts`,
`collapse_cuts_args`) are handed to the package's own
[`get_FSdata()`](https://larry-leon.github.io/forestsearch/reference/get_FSdata.md)
on `df_super`. Both directions of every cut are generated, as the search
does. LASSO, GRF and DINA cut sources are off: they are sample-fitted
screens, not part of the cut specification.

**Combinations.** All combinations of up to `maxk` indicator columns are
enumerated with
[`generate_combination_indices()`](https://larry-leon.github.io/forestsearch/reference/generate_combination_indices.md),
in the order
[`forestsearch()`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
uses when it composes the MR family.

**Structural floors** (all on the population frame, in the order of
[`evaluate_combination_with_status()`](https://larry-leon.github.io/forestsearch/reference/evaluate_combination_with_status.md)):

1.  combinations with a constant column or an empty pairwise
    intersection are skipped (the search's status 0);

2.  `minp`: every constituent factor must have population prevalence
    `>= minp`;

3.  `rmin`: each added factor must shrink the membership by more than
    `rmin` subjects *of a trial of size n*, i.e. by more than `rmin / n`
    in population proportion;

4.  size: `Pg >= n.min / n`, where `n.min` is resolved as
    [`forestsearch()`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
    resolves it – the supplied value (default 60), or
    `max(60, ceiling(n.min.frac * n))` when `n.min = NULL`.

The GRF variable-importance pre-screen, the effect screen, the
consistency screen and the near-duplicate removal are *not* applied.
Candidates whose population membership vectors are identical are
collapsed to the first one enumerated.

**Null DGMs.** A DGM whose `df_super$flag_harm` has no member (Q empty;
`generate_glm_dgm(model = "null")`) is detected structurally,
cross-checked against `dgm$model` when present. Under the null every
candidate has the same true effect, so the fields become: `beta_g` = the
common effect (the DGM's `effect_Qc` = `effect_ITT`, oriented as the
alternative path orients); `se_g` from the whole-population effective
variance (`fs_dgm_scale(dgm, regions = list(S = ...))`, the S row) at
`(n, Pg)` with the same prevalence scaling; `PQg = 0`; `sens_g = NA`
(0/0, undefined – not zero); `spec_g = 1 - Pg`; `PQ = 0` (from which NPV
= 1 follows downstream). Enumeration, floors, `ovl` and the covariance
are unchanged. The element `null` records the branch taken.

**Scale.** Every mean and standard error is derived from
[`fs_dgm_scale`](https://larry-leon.github.io/forestsearch/reference/fs_dgm_scale.md)`(dgm)`.
The orientation is the harm direction `s = sign(m_tau[Q])`, so the
planted effect is oriented positive: with `Q` the true harm region,
`tauQc = s * m_tau[Qc]`, `bint = s * (m_tau[Q] - m_tau[Qc])`,
`seQ1000 = sqrt(V_eff[Q] / (1000 * P(Q)))`, and for each candidate
`beta_g = tauQc + bint * PQg` – the signed mixture
`s * (m_tau[Qc] + (m_tau[Q] - m_tau[Qc]) * PQg)` – and
`se_g = seQ1000 * sqrt(1000 / n) * sqrt(P(Q) / Pg)` (sign-free, from
`V_eff[Q]`). Opposite-sign families (`sign(m_tau[Qc]) != s`) are
supported: benefit-direction candidates carry oriented-negative
`beta_g`, and `tauQc` may be negative for such families. When both
region effects share a sign the values coincide exactly with the former
oriented-absolute reading. A DGM with `m_tau[Q]` exactly zero is
rejected (no harm direction to orient by): plant a nonzero Q effect or
use the null path.

## See also

[`fs_oc_predict`](https://larry-leon.github.io/forestsearch/reference/fs_oc_predict.md),
[`fs_dgm_scale`](https://larry-leon.github.io/forestsearch/reference/fs_dgm_scale.md)

## Examples

``` r
# \donttest{
set.seed(1)
N <- 2000
age <- round(rnorm(N, 35, 9)); pre <- round(rexp(N, 1 / 500))
V   <- factor(rbinom(N, 1, 0.42), levels = 0:1)
inQ <- as.integer(age > 34 & pre <= 745)
mu0 <- 40 + 0.2 * age
dgm <- structure(list(
  df_super = data.frame(age = age, preanti = pre, V = V,
                        mu0 = mu0, mu1 = mu0 - 26 - 14 * inQ,
                        flag_harm = inQ),
  outcome_type = "continuous", effect_measure = "MD",
  model_params = list(sigma = 127.5)), class = c("glm_dgm", "list"))
fam <- fs_oc_family_enumerate(
  dgm, list(confounders.name = c("age", "preanti", "V"),
            conf.cont_jcuts = list(age = 4, preanti = 4), n.min = 60),
  n = 500)
fam
#> Population-enumerated candidate family (fs_oc_family)
#>   M            : 104 candidates (from 18 cut columns, maxk = 2)
#>   n            : 500
#>   floors       : n.min = 60 (Pg >= 0.1200), minp = 0.025, rmin = 5
#>   prevalence   : 0.1205 to 0.8215
#>   P(Q)         : 0.3980;  purity range 0.000 to 1.000
#>   dropped      : empty 21, minp 0, rmin 12, size 22, duplicate 12
# }
```
