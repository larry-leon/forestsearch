# Sampling Scale of the Difference-in-Means Estimator

Computes the finite-population sampling scale of the within-region
difference-in-means estimator directly from a data generating
mechanism's individual-level potential outcomes.

## Usage

``` r
fs_dgm_scale(
  dgm,
  regions = NULL,
  harm_col = "flag_harm",
  rand_ratio = 1,
  labels = c("Q", "Qc", "S")
)
```

## Arguments

- dgm:

  An object of class `"glm_dgm"` from
  [`generate_glm_dgm`](https://larry-leon.github.io/forestsearch/reference/generate_glm_dgm.md).

- regions:

  Optional named list of region specifications. Each element is either a
  logical vector of length `nrow(dgm$df_super)` or an integer vector of
  row indices. Defaults to `NULL`, which uses the true harm region, its
  complement, and the whole super-population.

- harm_col:

  Character. Name of the true-region indicator column in `df_super`.
  Default `"flag_harm"`.

- rand_ratio:

  Numeric. Treatment:control randomisation ratio, so the randomisation
  probability is `rand_ratio / (1 + rand_ratio)`. Default `1` (1:1),
  matching
  [`simulate_from_glm_dgm`](https://larry-leon.github.io/forestsearch/reference/simulate_from_glm_dgm.md).

- labels:

  Optional character vector of length 3 renaming the default regions.
  Ignored when `regions` is supplied. Default `c("Q", "Qc", "S")`.

## Value

An object of class `c("fs_dgm_scale", "list")`:

- `regions`:

  Data frame, one row per region, with columns `region`, `n_g`, `P_g`,
  `m_mu0`, `m_mu1`, `m_tau`, `V_mu0`, `V_mu1`, `V_tau`, `C_mu0_tau`,
  `v_cond0`, `v_cond1`, `V_arm0`, `V_arm1`, `bracket`, `V_eff`.

- `sigma`:

  Residual standard deviation for continuous outcomes, otherwise
  `NA_real_`.

- `outcome_type`, `effect_measure`:

  Copied from `dgm`.

- `rand_ratio`, `p_treat`:

  Allocation used.

- `n_super`:

  Super-population size.

## Details

For a region \\g\\ the estimator is the within-region difference in arm
means. Conditioning on the realized arm counts,
\$\$\mathrm{Var}\[\hat\beta(g)\] = V\_{\mathrm{eff}}(g) / (n P(g)),\$\$
with the region- and \\n\\-free constant \$\$V\_{\mathrm{eff}}(g) =
V_1(g)/p + V_0(g)/(1-p),\$\$ where \\p\\ is the randomisation
probability and \\V_w(g)\\ is the total within-arm outcome variance in
\\g\\ under arm \\w\\: \$\$V_w(g) = m_g\[v_w\] + V_g\[\mu_w\].\$\$ Here
\\\mu_w\\ is the arm mean surface, \\v_w(x) = \mathrm{Var}(Y \mid X = x,
W = w)\\ the conditional outcome variance, and \\m_g\[\cdot\]\\,
\\V_g\[\cdot\]\\ finite-population moments over the rows of the
super-population lying in \\g\\. The conditional variance is supplied by
the outcome family:

- continuous:

  \\v_w = \sigma^2\\, the residual variance (constant).

- binary:

  \\v_w = p_w (1 - p_w)\\.

- count:

  \\v_w = \mu_w\\ (Poisson).

## Balanced-arm reading

At 1:1 allocation \\V\_{\mathrm{eff}}(g) = 2\\V_0(g) + V_1(g)\\\\, and
`bracket` \\= V\_{\mathrm{eff}}(g)/4\\ admits the decomposition
\$\$\sigma^2 + V_g\[\mu_0\] + C_g\[\mu_0, \tau\] + \tfrac12
V_g\[\tau\],\$\$ with \\\tau = \mu_1 - \mu_0\\ the conditional average
treatment effect. The component columns are returned for every region
and every outcome type; the decomposition above is exact for continuous
outcomes at 1:1 allocation. Away from 1:1, `bracket` remains
\\V\_{\mathrm{eff}}/4\\ but no longer equals that sum, so prefer `V_eff`
for any downstream calculation.

## Idealisations

`V_eff` treats the region size and the arm split as fixed at their
expectations. A trial draws \\n_g \sim \mathrm{Bin}(n, P(g))\\ and \\n_1
\mid n_g \sim \mathrm{Bin}(n_g, p)\\, and by Jensen's inequality
\\E\[1/n_1 + 1/n_0\]\\ and \\E\[1/n_g\]\\ exceed their fixed-count
counterparts. Use
[`fs_scale_se`](https://larry-leon.github.io/forestsearch/reference/fs_scale_se.md)
with `jensen = TRUE` for the unconditional standard deviation.

## Scope

Exact for identity-scale effect measures (`"MD"`, `"RD"`, `"IRD"`),
where the estimator is a difference in means. Ratio measures (`"OR"`,
`"RR"`, `"IRR"`) require a delta-method layer and are rejected rather
than silently approximated.

## See also

[`fs_scale_se`](https://larry-leon.github.io/forestsearch/reference/fs_scale_se.md),
[`generate_glm_dgm`](https://larry-leon.github.io/forestsearch/reference/generate_glm_dgm.md)

## Examples

``` r
if (FALSE) { # \dontrun{
dgm <- generate_glm_dgm(..., outcome_type = "continuous",
                        effect_measure = "MD")
sc <- fs_dgm_scale(dgm)
sc$regions[, c("region", "P_g", "bracket", "V_eff")]

# Standard deviation of the estimator on the true region at n = 500
fs_scale_se(sc, n = 500, region = "Q")

# Arbitrary regions, e.g. a candidate family
fs_dgm_scale(dgm, regions = list(big = dgm$df_super$age > 40))
} # }
```
