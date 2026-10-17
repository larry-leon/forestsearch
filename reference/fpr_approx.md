# Approximate Procedure-Level False Positive Rate

Computes an analytical approximation to the ForestSearch false positive
rate under the null hypothesis of homogeneous treatment effect. The
calculation accounts for (i) the sampling variability of the effect
estimate in a subgroup of size `n_min`, (ii) the number of candidate
subgroups searched (driven by the covariate space and `maxk`), and (iii)
the screening threshold `c1`.

## Usage

``` r
fpr_approx(
  effect_true = 0.65,
  p_event = 0.4,
  n_min = 60L,
  c1 = 1.25,
  n_confounders_cont = 6L,
  n_confounders_bin = 6L,
  cont_cutoff = 4L,
  maxk = 2L,
  effect_type = c("ratio", "difference"),
  treat_alloc = 0.5,
  rho = 0.3,
  n_eff = NULL,
  thresholds = NULL,
  verbose = TRUE
)
```

## Arguments

- effect_true:

  Numeric. True effect under H0 on the ratio scale (e.g., OR or HR).
  Default 0.65 (protective treatment).

- p_event:

  Numeric in (0, 1). Overall event rate (control arm baseline rate for
  binary; used to compute SE). Default 0.40.

- n_min:

  Integer. Minimum subgroup size (`n.min` in
  [`forestsearch()`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)).
  Default 60.

- c1:

  Numeric. The **screening threshold** on a candidate subgroup's own
  effect –
  [`forestsearch`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)'s
  `effect.threshold` (legacy `hr.threshold`). Candidates with effect \>
  `c1` pass screening. Default 1.25. This approximation covers the
  screening stage only, so the per-split effect threshold `c2`
  (`consistency.threshold`) and the consistency rate `p*`
  (`pconsistency.threshold`) do not appear here.

- n_confounders_cont:

  Integer. Number of continuous confounders. Default 6.

- n_confounders_bin:

  Integer. Number of binary confounders. Default 6.

- cont_cutoff:

  Integer. Number of cut points per continuous covariate (`cont.cutoff`
  in
  [`forestsearch()`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)).
  Each continuous variable produces `cont_cutoff` candidate subgroups.
  Default 4.

- maxk:

  Integer. Maximum interaction depth. 1 = single-variable subgroups
  only; 2 = up to two-variable interactions. Default 2.

- effect_type:

  Character. `"ratio"` (OR, HR, IRR: log-scale SE) or `"difference"`
  (RD, MD: identity-scale SE). Default `"ratio"`.

- treat_alloc:

  Numeric in (0, 1). Treatment allocation fraction. Default 0.5 (1:1
  randomisation).

- rho:

  Numeric in \[0, 1). Assumed average pairwise correlation between
  candidate test statistics, used for the correlation-adjusted FPR
  estimate. Default 0.30. Set to 0 for the independence (upper-bound)
  estimate only. Ignored when `n_eff` is set.

- n_eff:

  Numeric or `NULL`. Pre-specified effective number of independent
  tests. When set, overrides the `rho`-based Sidak adjustment. Default
  `NULL` (compute from `rho`).

- thresholds:

  Numeric vector. Additional thresholds to evaluate besides `c1`.
  Default `c(1.0, 1.5, 2.0)`.

- verbose:

  Logical. Print a formatted summary. Default `TRUE`.

## Value

A list (class `"fpr_approximation"`) with components:

- `se_log_effect`:

  SE of log(effect) in a subgroup of size `n_min`.

- `n_candidates`:

  Total number of candidate subgroups searched.

- `n_single`:

  Number of single-variable candidates.

- `n_pairs`:

  Number of two-variable candidates (0 if `maxk == 1`).

- `per_candidate`:

  Named numeric vector: P(effect \> thr) for each threshold evaluated.

- `fpr_indep`:

  Named numeric vector: procedure-level FPR assuming independent
  candidates, for each threshold.

- `fpr_adjusted`:

  Named numeric vector: correlation-adjusted FPR, for each threshold.

- `params`:

  List of input parameters.

## Details

The approximation assumes:

- Normal approximation to the log-OR (or log-HR) distribution.

- Independent candidates (upper bound; true FPR is lower due to
  overlap). A correlation-adjusted estimate is also returned.

- Equal allocation to treatment arms within each subgroup.

- No consistency check adjustment (the reported FPR is the
  *screening-stage* FPR before `pconsistency.threshold` filtering).

## Examples

``` r
# ACTG175 binary harm search
fpr_approx(effect_true = 0.65, p_event = 0.40, n_min = 60, c1 = 1.25)
#> ============================================================
#>   ForestSearch FPR Approximation (screening stage)
#> ============================================================
#> 
#>   True effect (H0):   0.650  (log = -0.431)
#>   Event rate (ctrl):  40.0%
#>   Event rate (trt):   30.2%  (implied by OR = 0.65)
#>   Min subgroup (n):   60  (~30 per arm)
#>   SE(log effect):     0.545
#> 
#>   Confounders:        6 continuous (x4 cuts) + 6 binary
#>   Single candidates:  30
#>   Pair candidates:    435  (maxk = 2)
#>   Total candidates:   465
#>   Effective tests:    3  (rho = 0.30)
#> 
#>   Threshold   P(per candidate)   FPR (indep)   FPR (adj)
#>   ---------   ----------------   -----------   ---------
#>     1.25         11.5%             100.0%         33.3% <-- c1
#> 
#>   Note: Screening-stage only.  Consistency check (pconsistency)
#>   and bootstrap correction further reduce the operational FPR.
#>   Independence FPR is an upper bound; adjusted uses Sidak with
#>   effective tests (rho = 0.30).
#> ============================================================

# Tighter threshold
fpr_approx(effect_true = 0.65, p_event = 0.40, n_min = 60, c1 = 2.0)
#> ============================================================
#>   ForestSearch FPR Approximation (screening stage)
#> ============================================================
#> 
#>   True effect (H0):   0.650  (log = -0.431)
#>   Event rate (ctrl):  40.0%
#>   Event rate (trt):   30.2%  (implied by OR = 0.65)
#>   Min subgroup (n):   60  (~30 per arm)
#>   SE(log effect):     0.545
#> 
#>   Confounders:        6 continuous (x4 cuts) + 6 binary
#>   Single candidates:  30
#>   Pair candidates:    435  (maxk = 2)
#>   Total candidates:   465
#>   Effective tests:    3  (rho = 0.30)
#> 
#>   Threshold   P(per candidate)   FPR (indep)   FPR (adj)
#>   ---------   ----------------   -----------   ---------
#>     2.00          2.0%             100.0%          6.3% <-- c1
#> 
#>   Note: Screening-stage only.  Consistency check (pconsistency)
#>   and bootstrap correction further reduce the operational FPR.
#>   Independence FPR is an upper bound; adjusted uses Sidak with
#>   effective tests (rho = 0.30).
#> ============================================================

# Larger minimum subgroup
fpr_approx(effect_true = 0.65, p_event = 0.40, n_min = 100, c1 = 1.25)
#> ============================================================
#>   ForestSearch FPR Approximation (screening stage)
#> ============================================================
#> 
#>   True effect (H0):   0.650  (log = -0.431)
#>   Event rate (ctrl):  40.0%
#>   Event rate (trt):   30.2%  (implied by OR = 0.65)
#>   Min subgroup (n):   100  (~50 per arm)
#>   SE(log effect):     0.422
#> 
#>   Confounders:        6 continuous (x4 cuts) + 6 binary
#>   Single candidates:  30
#>   Pair candidates:    435  (maxk = 2)
#>   Total candidates:   465
#>   Effective tests:    3  (rho = 0.30)
#> 
#>   Threshold   P(per candidate)   FPR (indep)   FPR (adj)
#>   ---------   ----------------   -----------   ---------
#>     1.25          6.1%             100.0%         18.7% <-- c1
#> 
#>   Note: Screening-stage only.  Consistency check (pconsistency)
#>   and bootstrap correction further reduce the operational FPR.
#>   Independence FPR is an upper bound; adjusted uses Sidak with
#>   effective tests (rho = 0.30).
#> ============================================================
```
