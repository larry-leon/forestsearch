# Suggest Screening and Consistency Thresholds

Recommends `(c1, c2)` threshold pairs for ForestSearch based on the
corrected null approximation. Uses the per-subgroup detection
probability and an \\L\_{\text{eff}}\\ multiplicity correction to find
the smallest `c1` (screening threshold) such that the corrected FPR is
at most `fpr_target`.

## Usage

``` r
suggest_thresholds(
  d_eff,
  N,
  fpr_target = 0.1,
  n_min = 60,
  L_eff_C = 1,
  L_eff_alpha = 0,
  c1_range = c(1.05, 2.5),
  c1_step = 0.05,
  c2_candidates = seq(0.9, 2, by = 0.05),
  effect_scale = c("ratio", "difference")
)
```

## Arguments

- d_eff:

  Numeric. Effective information in the expected harm subgroup. Use
  [`d_eff_survival`](https://larry-leon.github.io/forestsearch/reference/d_eff_survival.md),
  [`d_eff_binary`](https://larry-leon.github.io/forestsearch/reference/d_eff_binary.md),
  [`d_eff_count`](https://larry-leon.github.io/forestsearch/reference/d_eff_count.md),
  or
  [`d_eff_continuous`](https://larry-leon.github.io/forestsearch/reference/d_eff_continuous.md).

- N:

  Integer. Total sample size.

- fpr_target:

  Numeric. Maximum acceptable procedure-level FPR. Default: 0.10.

- n_min:

  Integer. Minimum subgroup size. Default: 60.

- L_eff_C:

  Numeric. Calibration constant for \\L\_{\text{eff}} = C
  (N/n\_{\min})^\alpha\\. Default: 1.0 (conservative; see Details).

- L_eff_alpha:

  Numeric. Calibration exponent. Default: 0.0 (no N-dependence; see
  Details).

- c1_range:

  Numeric vector of length 2. Search range for c1. Default:
  `c(1.05, 2.50)`.

- c1_step:

  Numeric. Grid step for c1. Default: 0.05.

- c2_candidates:

  Numeric vector. Candidate c2 values to evaluate for each c1. Default:
  `seq(0.90, 2.0, by = 0.05)`.

- effect_scale:

  Character. `"ratio"` (default) or `"difference"`.

## Value

A data.frame with columns:

- c1:

  Screening threshold.

- min_c2:

  Minimum consistency threshold achieving the target.

- FPR_corr:

  Corrected FPR at this (c1, min_c2).

- P1:

  Per-subgroup detection probability.

- L_eff:

  Effective number of candidates at this N.

Sorted by c1 ascending. The first row is the smallest c1 that achieves
the target FPR.

## Details

The corrected FPR is: \$\$\text{FPR}\_{\text{corr}} = 1 - (1 -
P_1)^{L\_{\text{eff}}}\$\$

where \\P_1\\ = `compute_detection_probability_glm(1.0, d_eff, c1, c2)`
and \\L\_{\text{eff}} = C (N / n\_{\min})^\alpha\\.

**Default L_eff values (C=1, alpha=0)** assume \\L\_{\text{eff}} = 1\\,
which gives the per-subgroup approximation with no multiplicity
correction. This is appropriate when:

- No calibration data is available

- The analysis has few confounders and small N

For more accurate suggestions, calibrate \\L\_{\text{eff}}\\ from an H0
simulation using
[`calibrate_L_eff`](https://larry-leon.github.io/forestsearch/reference/calibrate_L_eff.md)
or
[`run_null_calibration`](https://larry-leon.github.io/forestsearch/reference/run_null_calibration.md).

**Reference calibrations** (from the forestsearch GLM extension analysis
documents):

- Binary, 4 confounders:

  C = 0.220, alpha = 1.298 (binary_threshold_calibration.qmd, 5000 sims)

- GBSG survival, 7 confounders:

  C = 0.029, alpha = 0.882 (cox_vs_glm_approximation.qmd, 500 sims per
  N)

## See also

[`compute_detection_probability_glm`](https://larry-leon.github.io/forestsearch/reference/compute_detection_probability_glm.md),
[`calibrate_L_eff`](https://larry-leon.github.io/forestsearch/reference/calibrate_L_eff.md),
[`predict_fpr_corrected`](https://larry-leon.github.io/forestsearch/reference/predict_fpr_corrected.md)

## Examples

``` r
# \donttest{
# Survival: paper defaults should be recovered
d <- d_eff_survival(n_sg = 60, prop_cens = 0.45)
suggest_thresholds(d_eff = d, N = 700)
#>      c1 min_c2 FPR_corr       P1 L_eff
#> 1  1.30   1.30 0.088652 0.088652     1
#> 2  1.35   1.30 0.087324 0.087324     1
#> 3  1.40   1.25 0.094620 0.094620     1
#> 4  1.45   1.20 0.096150 0.096150     1
#> 5  1.50   1.10 0.098494 0.098494     1
#> 6  1.55   0.90 0.098858 0.098858     1
#> 7  1.60   0.90 0.084792 0.084792     1
#> 8  1.65   0.90 0.072690 0.072690     1
#> 9  1.70   0.90 0.062144 0.062144     1
#> 10 1.75   0.90 0.052848 0.052848     1
#> 11 1.80   0.90 0.044966 0.044966     1
#> 12 1.85   0.90 0.038170 0.038170     1
#> 13 1.90   0.90 0.032550 0.032550     1
#> 14 1.95   0.90 0.027632 0.027632     1
#> 15 2.00   0.90 0.023324 0.023324     1
#> 16 2.05   0.90 0.019712 0.019712     1
#> 17 2.10   0.90 0.016704 0.016704     1
#> 18 2.15   0.90 0.014098 0.014098     1
#> 19 2.20   0.90 0.012010 0.012010     1
#> 20 2.25   0.90 0.010182 0.010182     1
#> 21 2.30   0.90 0.008516 0.008516     1
#> 22 2.35   0.90 0.007242 0.007242     1
#> 23 2.40   0.90 0.006178 0.006178     1
#> 24 2.45   0.90 0.005234 0.005234     1
#> 25 2.50   0.90 0.004434 0.004434     1

# Binary: needs higher thresholds
d <- d_eff_binary(n_sg = 150, p_event = 0.30)
suggest_thresholds(d_eff = d, N = 500,
                   L_eff_C = 0.220, L_eff_alpha = 1.298)
#>      c1 min_c2   FPR_corr       P1    L_eff
#> 1  1.65   1.65 0.08608888 0.025766 3.448621
#> 2  1.70   1.60 0.09827876 0.029552 3.448621
#> 3  1.75   1.60 0.09268443 0.027810 3.448621
#> 4  1.80   1.55 0.09444052 0.028356 3.448621
#> 5  1.85   1.45 0.09820185 0.029528 3.448621
#> 6  1.90   1.35 0.09671406 0.029064 3.448621
#> 7  1.95   1.00 0.09974552 0.030010 3.448621
#> 8  2.00   0.90 0.08671634 0.025960 3.448621
#> 9  2.05   0.90 0.07399799 0.022046 3.448621
#> 10 2.10   0.90 0.06339480 0.018812 3.448621
#> 11 2.15   0.90 0.05389984 0.015938 3.448621
#> 12 2.20   0.90 0.04632551 0.013660 3.448621
#> 13 2.25   0.90 0.03972054 0.011684 3.448621
#> 14 2.30   0.90 0.03388431 0.009946 3.448621
#> 15 2.35   0.90 0.02865819 0.008396 3.448621
#> 16 2.40   0.90 0.02459162 0.007194 3.448621
#> 17 2.45   0.90 0.02113811 0.006176 3.448621
#> 18 2.50   0.90 0.01803003 0.005262 3.448621

# Use calibration object from run_null_calibration()
# cal <- run_null_calibration(...)
# suggest_thresholds(d_eff = d, N = 700,
#                    L_eff_C = cal$C, L_eff_alpha = cal$alpha)
# }
```
