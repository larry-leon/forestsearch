# REPORT — ACTG175 continuous (MD) re-run under the current field constructions: campaign `mdsgnb20`

Date: 2026-09-16. Machine: `pop-os` (AMD Ryzen Threadripper PRO 5995WX, 64 physical / 128 logical cores, 251 GB; R 4.6.1, reference BLAS/LAPACK 3.12.0). Branch `feature/glm-extension`. Task: `dev/tasks/TASK_md_field_rerun_2026-09-15.md` (`bb84120e`), on Larry's approval of the Gate 0 recommendations (`REPORT_md_field_rerun_stage0_2026-09-15.md`, `b6e30ac7`). Records of this task: Stage 1 `REPORT_md_field_rerun_stage1_2026-09-15.md` (`e6abdb29`), Gate 2 `REPORT_md_field_rerun_gate2_2026-09-15.md` (per cell, committed by the runner), this record. Installed forestsearch 0.3.5, `Built: R 4.6.1; ; 2026-09-16 05:57:14 UTC; unix` (a rebuild of unchanged `R/`); no `R/` change in this task. Every render ran with `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1`.

## Gate 0 dispositions (Larry, 2026-09-15), as executed

- **D1, selection rule:** `effMaxSG`, ε = 0.20, `selection_rule = "neighborhood"` (the survival `effMaxSG` campaigns' rule, Stage 1 §1.3). The same-draws anchor to `mdf1` is data-level (`n_true` and the oracle columns), checked in both directions per cell (Gate 2).
- **D2, null cell:** presented and labelled by its truth (no subgroup; homogeneous +26 on the harm-oriented scale), with its declaration rate beside it.
- **D3, thresholds:** one ladder, τ = 0, 10, 20, 30, 40, 60, 80, 100 on the harm-oriented MD scale, both blocks, every cell, oracle beside: Ĥ lower bound ≥ τ ("harm of at least τ supported"); Ĥᶜ upper bound ≤ τ ("harm of at most τ supported").
- **D4, machine:** `pop-os`, 63 workers (Stage 1 §1.7).
- **D5, applied document:** not touched here (its own task).
- Also executed: the template's recorder / knob / `meta` edits (Stage 1 §1.4, template only); the display reads field-s through a renamed copy (no `R/` change); scripts, the progress log and the directory catalog committed; the extract is new code. **FB off** (the template default; Larry, 2026-09-15 overnight): no FB rows in any table or extract row.
- **Gate 1 go, given in advance** (Larry, 2026-09-15 overnight): on the condition that every Stage 1 gate was green and the §1.7 projection under 12 hours. Both held (Stage 1 record §1.8); Stage 2 launched at W = 63 with timeout 3,361 s per render and ceiling 19,071 s cumulative.

## Gate 1 and Gate 2 in brief

- **Gate 1 (Stage 1 record):** install rebuilt from HEAD (finding: `upgrade = "never"` is rejected by this devtools; `upgrade = FALSE` used); template edits E1–E5 transplanted from the survival m1 template (`2cffb95f`); scripts (`92d6ceb7`); smoke at defaults reproduced `mdf1` on sim_id 1–20 of all four cells with **0 selection flips** (17 / 19 / 14 / 20 of 20 rows identical to ≤ 4e-11 relative; 3 / 1 / 6 / 0 enumerated rows, all label ties in `mdf1`'s Gate 2 classes); field-s wiring green on every filled replicate; the campaign rule reproduced the draws (`n_true` identical, oracle columns ≤ 2.1e-12) with |Ĥ| grown on 17 of 20; calibration W = 63 (31.3 replicates/min; per-replicate `fit_mr_secs` 74.5 s mean, p90 95.5 s; peak summed RSS 72 GB), projection 212 min, ceiling 19,071 s, timeout 3,361 s.

- **Gate 2 (`REPORT_md_field_rerun_gate2_2026-09-15.md`, per cell):** every cell passed all 58 checks of `scripts_mdsgnb20/gate2.R`: 2,000 rows with `sim_id` exactly 1–2000 and the batch files matching the combined bundle on all 156 columns; `meta` carrying `effMaxSG` / ε 0.20 / `neighborhood`, `field_scale_complement = "selected"`, `pkg_version 0.3.5`, `hostname pop-os`, 63 workers; the same draws as `mdf1` in both directions (`n_true` identical on 2,000 rows; oracle columns within 2.8e-11 / 2.8e-11 / 2.5e-10 / 2.8e-11 relative for md40 n500 / md120 / null (complement only) / n700); field-s finite on every filled replicate with the interval invariants and the bound identities exact; MR failures on declared replicates 0. Declared: 1998 / 2000 / 1993 / 1999 of 2,000. Render walls (s): md40 n500 1045 + 1035 + 41; md120 1235 + 1245 + 41; null 1014 + 1010 + 41; md40 n700 1533 + 1503 + 40 — **9,783 s (163 min) in total** against the 212-min projection and the 19,071-s ceiling; peak summed RSS 72–73 GB at n = 500 and 82.5 GB at n = 700 (63 workers). Cell commits `72e6f129`, `b12f5983`, `dbeee517`, `f528d188`; campaign complete `2e4540c6`.
- **The one halt** (07:03:42Z, md40 n500): `gate2.R`'s own `sprintf` on the same-draws oracle label lacked its tolerance argument and threw after 52 of 53 checks had passed (no data check failed). The checker was fixed by explicit path (`34580fb0`, the check unchanged), re-run on the completed, un-rerendered cell output (58 / 58), the cell committed by hand with the runner's paths (`72e6f129`), the halt file removed (`76d6a8fd`) and the runner relaunched, skipping the committed cell. The relaunched runner's cumulative-wall counter restarted at 0 (its ceiling then applied to the remaining three cells, 7,662 s); the true cumulative, 9,783 s, is under the ceiling either way.

## Tables (pasted from the render of `summary_continuous_field_mdsgnb20.qmd`, `fece16e6`; every number is in `md_field_metrics.csv`, `4ab27707`)

Conventions: harm-oriented MD scale (positive = harm); targets β(Ĥ), β(Ĥᶜ) exact per replicate, the oracle against the structural true-region effect; one-sided coverage on the exposed side (Ĥ: LOWER bound; Ĥᶜ: UPPER bound); Wilson 95% intervals in parentheses; Monte Carlo SEs in parentheses in the ladder tables. The null cell is labelled by its truth (no subgroup; homogeneous +26); its oracle row on Ĥ is blank (Q is empty).

### Declaration

|cell                                     | replicates| declared|   rate|  mc_se| mdf1_rate|
|:----------------------------------------|----------:|--------:|------:|------:|---------:|
|md40 n500                                |       2000|     1998| 0.9990| 0.0007|    0.9990|
|md120 n500                               |       2000|     2000| 1.0000| 0.0000|    1.0000|
|null n500 (no subgroup; homogeneous +26) |       2000|     1993| 0.9965| 0.0013|    0.9965|
|md40 n700                                |       2000|     1999| 0.9995| 0.0005|    0.9995|

- Declaration is 0.997–1.000 in every cell, the null included, and identical to `mdf1`'s on the same seeds: the consistency screen at thresholds 30 / 10 declares on almost every replicate under either rule.
- The null cell declares on 1993 of 2000 replicates; every one of those Ĥ has β(Ĥ) = +26.26, so its coverage rows are read as coverage of a homogeneous harm, not as a false-claim rate.

### Ĥ — unadjusted, oracle, IJ two-term, field

|cell                                     |estimator  |    n| declaration rate| bias (MD)| bias (SD units)| SE/SD|one-sided LOWER coverage (Wilson) |two-sided coverage (Wilson) |
|:----------------------------------------|:----------|----:|----------------:|---------:|---------------:|-----:|:---------------------------------|:---------------------------|
|md40 n500                                |naive      | 1998|            0.999|    58.837|           3.536| 1.524|0.137 (0.122, 0.152)              |0.293 (0.273, 0.313)        |
|md40 n500                                |oracle     | 1998|            0.999|    -0.660|          -0.033| 0.991|0.952 (0.942, 0.961)              |0.943 (0.932, 0.952)        |
|md40 n500                                |MR (IJ)    | 1998|            0.999|    14.980|           0.765| 1.696|0.975 (0.968, 0.981)              |0.996 (0.992, 0.998)        |
|md40 n500                                |MR (field) | 1998|            0.999|     7.067|           0.333| 1.175|0.948 (0.937, 0.957)              |0.975 (0.967, 0.981)        |
|md120 n500                               |naive      | 2000|            1.000|    27.304|           1.491| 1.195|0.672 (0.652, 0.693)              |0.774 (0.755, 0.792)        |
|md120 n500                               |oracle     | 2000|            1.000|    -0.684|          -0.034| 0.991|0.953 (0.942, 0.961)              |0.943 (0.932, 0.952)        |
|md120 n500                               |MR (IJ)    | 2000|            1.000|    -4.147|          -0.181| 1.408|0.992 (0.986, 0.995)              |0.987 (0.980, 0.991)        |
|md120 n500                               |MR (field) | 2000|            1.000|    -6.403|          -0.255| 1.005|0.967 (0.958, 0.974)              |0.916 (0.903, 0.927)        |
|null n500 (no subgroup; homogeneous +26) |naive      | 1993|            0.997|    60.823|           3.705| 1.573|0.116 (0.103, 0.131)              |0.262 (0.244, 0.282)        |
|null n500 (no subgroup; homogeneous +26) |oracle     |    0|            0.997|       NaN|             NaN|    NA|NaN (NA, NA)                      |NaN (NA, NA)                |
|null n500 (no subgroup; homogeneous +26) |MR (IJ)    | 1993|            0.997|    16.196|           0.832| 1.724|0.974 (0.967, 0.980)              |0.994 (0.990, 0.997)        |
|null n500 (no subgroup; homogeneous +26) |MR (field) | 1993|            0.997|     7.973|           0.377| 1.187|0.947 (0.937, 0.956)              |0.975 (0.968, 0.981)        |
|md40 n700                                |naive      | 1999|            1.000|    60.861|           3.587| 1.469|0.085 (0.073, 0.098)              |0.208 (0.190, 0.226)        |
|md40 n700                                |oracle     | 1999|            1.000|     0.074|           0.004| 0.967|0.941 (0.930, 0.950)              |0.943 (0.932, 0.952)        |
|md40 n700                                |MR (IJ)    | 1999|            1.000|    15.806|           0.801| 1.582|0.958 (0.948, 0.966)              |0.986 (0.980, 0.990)        |
|md40 n700                                |MR (field) | 1999|            1.000|     7.725|           0.363| 1.154|0.939 (0.928, 0.949)              |0.959 (0.949, 0.967)        |

- The field's one-sided lower bound on β(Ĥ) covers at 0.948 / 0.967 / 0.947 / 0.939 (md40 n500 / md120 n500 / null / md40 n700); the Wilson interval contains 0.95 in the md40 n500 and null cells, sits above it at md120 (0.958–0.974) and just below at n = 700 (0.928–0.949). The oracle's is 0.941–0.953.
- Retained bias in SD units: unadjusted +3.5 / +1.5 / +3.7 / +3.6, IJ +0.77 / −0.18 / +0.83 / +0.80, field +0.33 / −0.26 / +0.38 / +0.36; at md120 both corrections now overshoot slightly (negative bias) where `mdf1` had +0.36 / +0.06. IJ's SE/SD is 1.41–1.72 (two-sided 0.986–0.996); the field's is 1.00–1.19.
- The field's two-sided interval covers at 0.975 / 0.916 / 0.975 / 0.959: at md120 it is below nominal while the one-sided lower bound is above it, the two sides of an asymmetric Λ* about β̃.

### Ĥᶜ — unadjusted, oracle, IJ two-term, field (paired before/after), field-s (evaluated)

|cell                                     |estimator    |    n| declaration rate| bias (MD)| bias (SD units)| SE/SD|one-sided UPPER coverage (Wilson) |two-sided coverage (Wilson) |
|:----------------------------------------|:------------|----:|----------------:|---------:|---------------:|-----:|:---------------------------------|:---------------------------|
|md40 n500                                |naive        | 1998|            0.999|   -16.798|          -1.374| 1.080|0.660 (0.639, 0.681)              |0.780 (0.762, 0.798)        |
|md40 n500                                |oracle       | 1998|            0.999|    -0.262|          -0.018| 1.005|0.957 (0.947, 0.965)              |0.953 (0.943, 0.962)        |
|md40 n500                                |MR (IJ)      | 1998|            0.999|    -5.995|          -0.473| 1.875|0.994 (0.990, 0.997)              |0.999 (0.996, 1.000)        |
|md40 n500                                |MR (field)   | 1998|            0.999|    -3.851|          -0.298| 0.998|0.914 (0.901, 0.926)              |0.933 (0.922, 0.944)        |
|md40 n500                                |MR (field-s) | 1998|            0.999|    -3.861|          -0.299| 1.005|0.914 (0.901, 0.925)              |0.936 (0.924, 0.946)        |
|md120 n500                               |naive        | 2000|            1.000|   -11.239|          -0.770| 0.970|0.804 (0.786, 0.821)              |0.875 (0.860, 0.889)        |
|md120 n500                               |oracle       | 2000|            1.000|    -0.301|          -0.021| 1.002|0.956 (0.946, 0.964)              |0.953 (0.942, 0.961)        |
|md120 n500                               |MR (IJ)      | 2000|            1.000|    -1.764|          -0.115| 1.639|0.997 (0.993, 0.999)              |0.999 (0.996, 1.000)        |
|md120 n500                               |MR (field)   | 2000|            1.000|    -0.472|          -0.030| 0.866|0.938 (0.927, 0.948)              |0.925 (0.913, 0.936)        |
|md120 n500                               |MR (field-s) | 2000|            1.000|    -0.508|          -0.032| 0.892|0.941 (0.929, 0.950)              |0.930 (0.919, 0.941)        |
|null n500 (no subgroup; homogeneous +26) |naive        | 1993|            0.997|   -16.505|          -1.367| 1.087|0.667 (0.646, 0.688)              |0.782 (0.764, 0.800)        |
|null n500 (no subgroup; homogeneous +26) |oracle       | 1993|            0.997|    -0.317|          -0.028| 1.026|0.953 (0.943, 0.962)              |0.956 (0.946, 0.964)        |
|null n500 (no subgroup; homogeneous +26) |MR (IJ)      | 1993|            0.997|    -5.806|          -0.463| 1.888|0.996 (0.992, 0.998)              |0.999 (0.997, 1.000)        |
|null n500 (no subgroup; homogeneous +26) |MR (field)   | 1993|            0.997|    -3.690|          -0.288| 1.004|0.915 (0.902, 0.927)              |0.942 (0.931, 0.952)        |
|null n500 (no subgroup; homogeneous +26) |MR (field-s) | 1993|            0.997|    -3.702|          -0.290| 1.008|0.917 (0.904, 0.929)              |0.942 (0.931, 0.951)        |
|md40 n700                                |naive        | 1999|            1.000|   -11.591|          -1.117| 1.040|0.720 (0.700, 0.740)              |0.822 (0.805, 0.839)        |
|md40 n700                                |oracle       | 1999|            1.000|     0.186|           0.015| 1.007|0.955 (0.945, 0.963)              |0.950 (0.940, 0.959)        |
|md40 n700                                |MR (IJ)      | 1999|            1.000|    -3.759|          -0.354| 1.876|0.994 (0.990, 0.997)              |0.999 (0.997, 1.000)        |
|md40 n700                                |MR (field)   | 1999|            1.000|    -2.161|          -0.200| 0.984|0.926 (0.914, 0.937)              |0.942 (0.931, 0.951)        |
|md40 n700                                |MR (field-s) | 1999|            1.000|    -2.168|          -0.201| 0.989|0.928 (0.916, 0.938)              |0.943 (0.932, 0.953)        |

- The field-s one-sided upper bound on β(Ĥᶜ) covers at 0.914 / 0.941 / 0.917 / 0.928; the unstudentized field beside it at 0.914 / 0.938 / 0.915 / 0.926. On this design the two constructions differ by at most 0.003 in coverage, 0.01 MD in bias and 0.03 in SE/SD: λ-SDᶜ / naive SEᶜ is 0.97–0.98 for the field and 0.98–1.00 for field-s, and SD(β̃ᶜ) / naive SEᶜ is 0.96–1.09, so the studentization has almost nothing to rescale here.
- Both are under 0.95 at n = 500 (Wilson upper limits 0.925–0.929 in md40 and null) and nearer it at md120 (0.941, Wilson 0.929–0.950) and n = 700 (0.928, 0.916–0.938); the oracle covers at 0.953–0.957. IJ's upper bound covers at 0.994–0.997 with SE/SD 1.64–1.89.
- Retained bias in SD units on Ĥᶜ: unadjusted −0.77 to −1.37, IJ −0.12 to −0.47, field-s −0.03 to −0.30.

### Bound location on Ĥ — one-sided 95% LOWER bound, the D3 ladder (field | oracle; IJ and unadjusted for orientation)

|cell                                     |estimator  |    n|  mean|   q05|   q25| median|   q75|   q95|P(L>=0)       |P(L>=10)      |P(L>=20)      |P(L>=30)      |P(L>=40)      |P(L>=60)      |P(L>=80)      |P(L>=100)     |
|:----------------------------------------|:----------|----:|-----:|-----:|-----:|------:|-----:|-----:|:-------------|:-------------|:-------------|:-------------|:-------------|:-------------|:-------------|:-------------|
|md40 n500                                |MR (field) | 1998|  -3.5| -35.2| -18.0|   -3.9|  10.0|  31.2|0.412 (0.011) |0.250 (0.010) |0.126 (0.007) |0.055 (0.005) |0.027 (0.004) |0.005 (0.002) |0.001 (0.001) |0.000 (0.000) |
|md40 n500                                |oracle     | 1998|   6.8| -27.0|  -5.9|    7.3|  20.1|  39.8|0.639 (0.011) |0.442 (0.011) |0.252 (0.010) |0.119 (0.007) |0.048 (0.005) |0.003 (0.001) |0.000 (0.000) |0.000 (0.000) |
|md40 n500                                |MR (IJ)    | 1998|  -7.8| -37.3| -21.5|   -8.3|   4.5|  23.7|0.334 (0.011) |0.176 (0.009) |0.067 (0.006) |0.027 (0.004) |0.009 (0.002) |0.001 (0.001) |0.000 (0.000) |0.000 (0.000) |
|md40 n500                                |naive      | 1998|  49.0|  24.9|  38.3|   48.1|  59.0|  74.2|1.000 (0.000) |0.999 (0.001) |0.978 (0.003) |0.899 (0.007) |0.706 (0.010) |0.230 (0.009) |0.026 (0.004) |0.003 (0.001) |
|md120 n500                               |MR (field) | 2000|  51.5|  11.7|  33.9|   50.3|  68.2|  93.2|0.991 (0.002) |0.961 (0.004) |0.903 (0.007) |0.806 (0.009) |0.667 (0.011) |0.357 (0.011) |0.129 (0.007) |0.029 (0.004) |
|md120 n500                               |oracle     | 2000|  86.7|  53.0|  74.1|   87.3| 100.1| 119.8|1.000 (0.000) |1.000 (0.000) |0.999 (0.001) |0.997 (0.001) |0.986 (0.003) |0.907 (0.006) |0.638 (0.011) |0.252 (0.010) |
|md120 n500                               |MR (IJ)    | 2000|  43.3|   8.4|  29.1|   42.5|  57.0|  79.0|0.982 (0.003) |0.939 (0.005) |0.872 (0.007) |0.737 (0.010) |0.555 (0.011) |0.207 (0.009) |0.045 (0.005) |0.003 (0.001) |
|md120 n500                               |naive      | 2000|  91.7|  63.9|  79.6|   91.0| 103.3| 121.7|1.000 (0.000) |1.000 (0.000) |1.000 (0.000) |1.000 (0.000) |0.999 (0.001) |0.972 (0.004) |0.743 (0.010) |0.305 (0.010) |
|null n500 (no subgroup; homogeneous +26) |MR (field) | 1993|  -8.6| -40.0| -22.9|   -9.8|   4.4|  27.5|0.325 (0.010) |0.179 (0.009) |0.087 (0.006) |0.040 (0.004) |0.015 (0.003) |0.003 (0.001) |0.000 (0.000) |0.000 (0.000) |
|null n500 (no subgroup; homogeneous +26) |oracle     |    0|    NA|    NA|    NA|     NA|    NA|    NA|NA            |NA            |NA            |NA            |NA            |NA            |NA            |NA            |
|null n500 (no subgroup; homogeneous +26) |MR (IJ)    | 1993| -12.7| -41.7| -25.8|  -13.7|  -0.5|  19.3|0.240 (0.010) |0.113 (0.007) |0.048 (0.005) |0.017 (0.003) |0.004 (0.001) |0.001 (0.001) |0.000 (0.000) |0.000 (0.000) |
|null n500 (no subgroup; homogeneous +26) |naive      | 1993|  44.6|  20.1|  34.1|   44.0|  54.4|  70.5|1.000 (0.000) |1.000 (0.000) |0.951 (0.005) |0.829 (0.008) |0.600 (0.011) |0.160 (0.008) |0.018 (0.003) |0.002 (0.001) |
|md40 n700                                |MR (field) | 1999|  -2.6| -33.4| -16.5|   -4.2|   9.4|  35.2|0.419 (0.011) |0.244 (0.010) |0.125 (0.007) |0.070 (0.006) |0.037 (0.004) |0.008 (0.002) |0.000 (0.000) |0.000 (0.000) |
|md40 n700                                |oracle     | 1999|  12.5| -15.5|   0.8|   12.0|  24.0|  41.5|0.763 (0.010) |0.552 (0.011) |0.322 (0.010) |0.165 (0.008) |0.059 (0.005) |0.004 (0.001) |0.000 (0.000) |0.000 (0.000) |
|md40 n700                                |MR (IJ)    | 1999|  -3.7| -31.7| -15.9|   -4.7|   7.0|  30.1|0.389 (0.011) |0.206 (0.009) |0.102 (0.007) |0.051 (0.005) |0.017 (0.003) |0.003 (0.001) |0.000 (0.000) |0.000 (0.000) |
|md40 n700                                |naive      | 1999|  51.7|  28.4|  41.5|   50.6|  60.4|  79.1|1.000 (0.000) |0.999 (0.001) |0.989 (0.002) |0.940 (0.005) |0.782 (0.009) |0.259 (0.010) |0.048 (0.005) |0.006 (0.002) |

- md40 cells (planted MD 40, prevalence 0.345): the field's lower bound lies at or above 0 on 41.2% / 41.9% of trials, at or above 30 on 5.5% / 7.0%, at or above 40 on 2.7% / 3.7%; the oracle's, which knows Q, at or above 0 on 63.9% / 76.3% and at or above 40 on 4.8% / 5.9%. A harm claim of any size on the selected Ĥ is supported on two trials in five; at the planted size on one in thirty.
- md120: the field's lower bound is at or above 40 on 66.7%, above 60 on 35.7%, above 80 on 12.9%, above 100 on 2.9%; the oracle's above 60 on 90.7% and above 80 on 63.8%. The unadjusted bound says "at least 60" on 97.2% of trials and "at least 100" on 30.5%.
- Null (no subgroup; homogeneous +26): the field's lower bound is at or above 0 on 32.5%, above 10 on 17.9%, above 20 on 8.7%, above 30 on 4.0% — read against a true +26 everywhere, this is the bound's location under a homogeneous harm, not a false-claim rate; the unadjusted bound says "at least 30" on 82.9%.
- Against `mdf1` (Stage 0 §7.3, `maxeffCons`): the field's lower bound at or above 0 moves from 30.7% to 41.2% in md40 n500 and from 24.0% to 32.5% in the null cell as Ĥ grows under `effMaxSG`; its 5–95% range is −35 to 31 (md40 n500) against −44 to 31 before.

### Bound location on Ĥᶜ — one-sided 95% UPPER bound, the D3 ladder (field-s | oracle; field, IJ and unadjusted for orientation)

|cell                                     |estimator    |    n| mean|  q05|  q25| median|  q75|   q95|P(U<=0)       |P(U<=10)      |P(U<=20)      |P(U<=30)      |P(U<=40)      |P(U<=60)      |P(U<=80)      |P(U<=100)     |
|:----------------------------------------|:------------|----:|----:|----:|----:|------:|----:|-----:|:-------------|:-------------|:-------------|:-------------|:-------------|:-------------|:-------------|:-------------|
|md40 n500                                |MR (field-s) | 1998| 48.2| 26.9| 39.8|   48.1| 56.5|  70.1|0.000 (0.000) |0.002 (0.001) |0.017 (0.003) |0.076 (0.006) |0.253 (0.010) |0.830 (0.008) |0.995 (0.002) |1.000 (0.000) |
|md40 n500                                |oracle       | 1998| 49.7| 27.1| 39.8|   50.0| 59.5|  72.5|0.001 (0.001) |0.006 (0.002) |0.020 (0.003) |0.080 (0.006) |0.252 (0.010) |0.758 (0.010) |0.984 (0.003) |1.000 (0.000) |
|md40 n500                                |MR (field)   | 1998| 48.1| 26.9| 39.7|   48.1| 56.7|  70.6|0.000 (0.000) |0.003 (0.001) |0.018 (0.003) |0.079 (0.006) |0.258 (0.010) |0.825 (0.008) |0.994 (0.002) |1.000 (0.000) |
|md40 n500                                |MR (IJ)      | 1998| 63.8| 43.6| 55.5|   63.9| 71.9|  85.6|0.000 (0.000) |0.000 (0.000) |0.001 (0.001) |0.006 (0.002) |0.033 (0.004) |0.371 (0.011) |0.897 (0.007) |0.999 (0.001) |
|md40 n500                                |naive        | 1998| 35.6| 16.0| 27.7|   35.7| 43.5|  56.6|0.002 (0.001) |0.020 (0.003) |0.094 (0.007) |0.315 (0.010) |0.646 (0.011) |0.978 (0.003) |1.000 (0.000) |1.000 (0.000) |
|md120 n500                               |MR (field-s) | 2000| 61.8| 37.1| 51.2|   62.0| 72.7|  86.7|0.000 (0.000) |0.003 (0.001) |0.007 (0.002) |0.019 (0.003) |0.071 (0.006) |0.450 (0.011) |0.881 (0.007) |0.997 (0.001) |
|md120 n500                               |oracle       | 2000| 49.7| 26.8| 39.8|   49.9| 59.5|  72.5|0.001 (0.001) |0.006 (0.002) |0.021 (0.003) |0.081 (0.006) |0.253 (0.010) |0.758 (0.010) |0.985 (0.003) |1.000 (0.000) |
|md120 n500                               |MR (field)   | 2000| 61.2| 35.9| 50.5|   61.5| 72.1|  86.9|0.001 (0.000) |0.003 (0.001) |0.007 (0.002) |0.024 (0.003) |0.080 (0.006) |0.466 (0.011) |0.887 (0.007) |0.997 (0.001) |
|md120 n500                               |MR (IJ)      | 2000| 78.8| 55.5| 68.6|   79.2| 89.0| 102.5|0.000 (0.000) |0.000 (0.000) |0.001 (0.000) |0.003 (0.001) |0.006 (0.002) |0.095 (0.007) |0.522 (0.011) |0.931 (0.006) |
|md120 n500                               |naive        | 2000| 51.2| 28.8| 41.4|   51.3| 61.0|  74.0|0.001 (0.000) |0.005 (0.002) |0.018 (0.003) |0.061 (0.005) |0.216 (0.009) |0.732 (0.010) |0.983 (0.003) |1.000 (0.000) |
|null n500 (no subgroup; homogeneous +26) |MR (field-s) | 1993| 43.7| 22.9| 35.5|   43.8| 52.2|  65.6|0.001 (0.001) |0.006 (0.002) |0.032 (0.004) |0.139 (0.008) |0.379 (0.011) |0.899 (0.007) |0.998 (0.001) |1.000 (0.000) |
|null n500 (no subgroup; homogeneous +26) |oracle       | 1993| 45.1| 26.6| 37.5|   45.1| 52.6|  64.0|0.000 (0.000) |0.001 (0.001) |0.015 (0.003) |0.090 (0.006) |0.332 (0.011) |0.897 (0.007) |0.999 (0.001) |1.000 (0.000) |
|null n500 (no subgroup; homogeneous +26) |MR (field)   | 1993| 43.7| 22.5| 35.3|   43.7| 52.3|  65.7|0.001 (0.001) |0.006 (0.002) |0.031 (0.004) |0.145 (0.008) |0.383 (0.011) |0.901 (0.007) |0.999 (0.001) |1.000 (0.000) |
|null n500 (no subgroup; homogeneous +26) |MR (IJ)      | 1993| 59.4| 39.1| 51.3|   59.5| 67.7|  81.2|0.000 (0.000) |0.000 (0.000) |0.001 (0.001) |0.012 (0.002) |0.057 (0.005) |0.517 (0.011) |0.937 (0.005) |1.000 (0.000) |
|null n500 (no subgroup; homogeneous +26) |naive        | 1993| 31.4| 12.2| 23.6|   31.1| 39.3|  52.0|0.005 (0.002) |0.034 (0.004) |0.175 (0.009) |0.456 (0.011) |0.769 (0.009) |0.992 (0.002) |1.000 (0.000) |1.000 (0.000) |
|md40 n700                                |MR (field-s) | 1999| 46.1| 28.9| 38.9|   45.9| 53.5|  63.5|0.000 (0.000) |0.001 (0.001) |0.010 (0.002) |0.063 (0.005) |0.290 (0.010) |0.902 (0.007) |0.999 (0.001) |1.000 (0.000) |
|md40 n700                                |oracle       | 1999| 46.5| 27.2| 38.4|   46.6| 54.6|  66.5|0.000 (0.000) |0.002 (0.001) |0.012 (0.002) |0.086 (0.006) |0.295 (0.010) |0.872 (0.007) |0.996 (0.001) |1.000 (0.000) |
|md40 n700                                |MR (field)   | 1999| 46.1| 28.7| 38.9|   46.0| 53.6|  63.5|0.000 (0.000) |0.001 (0.001) |0.010 (0.002) |0.066 (0.006) |0.294 (0.010) |0.902 (0.007) |0.999 (0.001) |1.000 (0.000) |
|md40 n700                                |MR (IJ)      | 1999| 59.8| 42.9| 52.6|   59.6| 67.2|  76.9|0.000 (0.000) |0.000 (0.000) |0.001 (0.001) |0.005 (0.001) |0.025 (0.003) |0.517 (0.011) |0.971 (0.004) |1.000 (0.000) |
|md40 n700                                |naive        | 1999| 36.9| 20.2| 30.2|   36.8| 44.0|  53.9|0.001 (0.001) |0.008 (0.002) |0.048 (0.005) |0.244 (0.010) |0.625 (0.011) |0.988 (0.002) |1.000 (0.000) |1.000 (0.000) |

- "Harm of at most 30" on Ĥᶜ is supported by the field-s upper bound on 7.6% / 1.9% / 13.9% / 6.3% of trials (md40 n500 / md120 / null / n700) and by the oracle on 8.0% / 8.1% / 9.0% / 8.6%; "at most 40" on 25.3% / 7.1% / 37.9% / 29.0% (oracle 25.2% / 25.3% / 33.2% / 29.5%); "at most 60" on 83.0% / 45.0% / 89.9% / 90.2%.
- Field-s and the unstudentized field agree to within 0.005 at every τ in every cell; IJ's upper bound sits about 15 MD higher (its "at most 60" share is 37.1% where field-s's is 83.0%, md40 n500).
- At md120 two thirds of Q sit in Ĥᶜ (β(Ĥᶜ) ≈ +52 under `mdf1`'s rule; larger recovery of Q here), so the complement bound is located higher there: median 62.0 against 48.1 / 43.8 / 45.9 in the other cells.

### Joint pair (Ĥ lower, Ĥᶜ upper)

|cell                                     |pair                                            | declared| both_bounds| share_both|joint                | joint_mc_se| cov_H| cov_Hc| margin_H| margin_Hc|
|:----------------------------------------|:-----------------------------------------------|--------:|-----------:|----------:|:--------------------|-----------:|-----:|------:|--------:|---------:|
|md40 n500                                |Bonferroni field-s (gamma = 0.025)              |     1998|        1998|          1|0.928 (0.916, 0.938) |       0.006| 0.976|  0.950|   59.768|    27.501|
|md40 n500                                |Bonferroni unstudentized (gamma = 0.025)        |     1998|        1998|          1|0.925 (0.913, 0.936) |       0.006| 0.976|  0.947|   59.768|    27.473|
|md40 n500                                |separate 95% bounds: field lower, field-s upper |     1998|        1998|          1|0.865 (0.850, 0.880) |       0.008| 0.948|  0.914|   50.387|    23.435|
|md120 n500                               |Bonferroni field-s (gamma = 0.025)              |     2000|        2000|          1|0.949 (0.938, 0.958) |       0.005| 0.985|  0.964|   53.855|    28.884|
|md120 n500                               |Bonferroni unstudentized (gamma = 0.025)        |     2000|        2000|          1|0.949 (0.938, 0.957) |       0.005| 0.985|  0.963|   53.855|    28.191|
|md120 n500                               |separate 95% bounds: field lower, field-s upper |     2000|        2000|          1|0.910 (0.897, 0.922) |       0.006| 0.967|  0.941|   44.747|    24.475|
|null n500 (no subgroup; homogeneous +26) |Bonferroni field-s (gamma = 0.025)              |     1993|        1993|          1|0.932 (0.920, 0.942) |       0.006| 0.976|  0.954|   60.475|    27.325|
|null n500 (no subgroup; homogeneous +26) |Bonferroni unstudentized (gamma = 0.025)        |     1993|        1993|          1|0.931 (0.919, 0.942) |       0.006| 0.976|  0.954|   60.475|    27.365|
|null n500 (no subgroup; homogeneous +26) |separate 95% bounds: field lower, field-s upper |     1993|        1993|          1|0.869 (0.853, 0.883) |       0.008| 0.947|  0.917|   51.053|    23.282|
|md40 n700                                |Bonferroni field-s (gamma = 0.025)              |     1999|        1999|          1|0.927 (0.915, 0.938) |       0.006| 0.963|  0.962|   59.880|    22.416|
|md40 n700                                |Bonferroni unstudentized (gamma = 0.025)        |     1999|        1999|          1|0.926 (0.914, 0.937) |       0.006| 0.963|  0.961|   59.880|    22.390|
|md40 n700                                |separate 95% bounds: field lower, field-s upper |     1999|        1999|          1|0.869 (0.854, 0.883) |       0.008| 0.939|  0.928|   50.285|    19.075|

|cell                                     |metric       |  value|  mc_se|    n|
|:----------------------------------------|:------------|------:|------:|----:|
|md40 n500                                |gamma_mean_s | 0.0251| 0.0000| 1998|
|md40 n500                                |corr_s       | 0.0883| 0.0014| 1998|
|md40 n500                                |gamma_mean   | 0.0251| 0.0000| 1998|
|md40 n500                                |corr         | 0.0920| 0.0014| 1998|
|md120 n500                               |gamma_mean_s | 0.0252| 0.0000| 2000|
|md120 n500                               |corr_s       | 0.0009| 0.0017| 2000|
|md120 n500                               |gamma_mean   | 0.0251| 0.0000| 2000|
|md120 n500                               |corr         | 0.0047| 0.0017| 2000|
|null n500 (no subgroup; homogeneous +26) |gamma_mean_s | 0.0251| 0.0000| 1993|
|null n500 (no subgroup; homogeneous +26) |corr_s       | 0.0884| 0.0014| 1993|
|null n500 (no subgroup; homogeneous +26) |gamma_mean   | 0.0251| 0.0000| 1993|
|null n500 (no subgroup; homogeneous +26) |corr         | 0.0919| 0.0014| 1993|
|md40 n700                                |gamma_mean_s | 0.0251| 0.0000| 1999|
|md40 n700                                |corr_s       | 0.0885| 0.0013| 1999|
|md40 n700                                |gamma_mean   | 0.0251| 0.0000| 1999|
|md40 n700                                |corr         | 0.0919| 0.0013| 1999|

- The field-s Bonferroni pair (γ = 0.025 each) covers (β(Ĥ), β(Ĥᶜ)) jointly on 0.928 / 0.949 / 0.932 / 0.927 of declared replicates (Wilson lower limits 0.915–0.938); both bounds exist on every declared replicate. The unstudentized pair is within 0.003 of it; the separate 95% pair covers at 0.865–0.910.
- The calibrated split returns the Bonferroni floor in every cell (mean γ 0.0251–0.0252) with corr(Λ*, Λ*ᶜ) 0.001–0.092: nothing is gained over Bonferroni here, as on `mdf1` and on survival.
- Compared with `mdf1`'s Bonferroni 0.940 / 0.943 / 0.940 / 0.936 (Stage 0 §7.2), the pair covers about 0.01 less in the md40 and null cells under `effMaxSG`, the complement side (cov_Hc 0.950–0.964) being the smaller component.

### Rule contrast with `mdf1`, paired by sim_id

|cell                                     |quantity                                                           |effMaxSG       |maxeffCons_mdf1 |
|:----------------------------------------|:------------------------------------------------------------------|:--------------|:---------------|
|md40 n500                                |declared replicates                                                |1998.0000      |1998.0000       |
|md40 n500                                |mean size of Hhat (n_harm)                                     |111.6241       |72.0781         |
|md40 n500                                |mean n_sel (gate)                                                  |111.6241       |72.0781         |
|md40 n500                                |sensitivity (mean over declared; NA where n_true = 0)              |0.2679         |0.1660          |
|md40 n500                                |PPV (mean over declared)                                           |0.4122         |0.3978          |
|md40 n500                                |size of Hhat grew / stayed / shrank vs mdf1 (paired by sim_id) |1881 / 117 / 0 |-               |
|md40 n500                                |field one-sided lower coverage on beta(Hhat)                       |0.9479         |0.9505          |
|md40 n500                                |unstudentized field one-sided upper coverage on beta(Hhat^c)       |0.9144         |0.9239          |
|md40 n500                                |field-s one-sided upper coverage on beta(Hhat^c)                   |0.9139         |NA              |
|md120 n500                               |declared replicates                                                |2000.0000      |2000.0000       |
|md120 n500                               |mean size of Hhat (n_harm)                                     |152.7385       |74.8425         |
|md120 n500                               |mean n_sel (gate)                                                  |152.7385       |74.8425         |
|md120 n500                               |sensitivity (mean over declared; NA where n_true = 0)              |0.7116         |0.3345          |
|md120 n500                               |PPV (mean over declared)                                           |0.7980         |0.7567          |
|md120 n500                               |size of Hhat grew / stayed / shrank vs mdf1 (paired by sim_id) |1988 / 12 / 0  |-               |
|md120 n500                               |field one-sided lower coverage on beta(Hhat)                       |0.9670         |0.9480          |
|md120 n500                               |unstudentized field one-sided upper coverage on beta(Hhat^c)       |0.9380         |0.9345          |
|md120 n500                               |field-s one-sided upper coverage on beta(Hhat^c)                   |0.9405         |NA              |
|null n500 (no subgroup; homogeneous +26) |declared replicates                                                |1993.0000      |1993.0000       |
|null n500 (no subgroup; homogeneous +26) |mean size of Hhat (n_harm)                                     |107.4922       |72.0768         |
|null n500 (no subgroup; homogeneous +26) |mean n_sel (gate)                                                  |107.4922       |72.0768         |
|null n500 (no subgroup; homogeneous +26) |sensitivity (mean over declared; NA where n_true = 0)              |NA             |NA              |
|null n500 (no subgroup; homogeneous +26) |PPV (mean over declared)                                           |0.0000         |0.0000          |
|null n500 (no subgroup; homogeneous +26) |size of Hhat grew / stayed / shrank vs mdf1 (paired by sim_id) |1844 / 149 / 0 |-               |
|null n500 (no subgroup; homogeneous +26) |field one-sided lower coverage on beta(Hhat)                       |0.9473         |0.9503          |
|null n500 (no subgroup; homogeneous +26) |unstudentized field one-sided upper coverage on beta(Hhat^c)       |0.9152         |0.9257          |
|null n500 (no subgroup; homogeneous +26) |field-s one-sided upper coverage on beta(Hhat^c)                   |0.9172         |NA              |
|md40 n700                                |declared replicates                                                |1999.0000      |1999.0000       |
|md40 n700                                |mean size of Hhat (n_harm)                                     |117.9905       |73.2571         |
|md40 n700                                |mean n_sel (gate)                                                  |117.9905       |73.2571         |
|md40 n700                                |sensitivity (mean over declared; NA where n_true = 0)              |0.2041         |0.1210          |
|md40 n700                                |PPV (mean over declared)                                           |0.4101         |0.3968          |
|md40 n700                                |size of Hhat grew / stayed / shrank vs mdf1 (paired by sim_id) |1848 / 151 / 0 |-               |
|md40 n700                                |field one-sided lower coverage on beta(Hhat)                       |0.9395         |0.9470          |
|md40 n700                                |unstudentized field one-sided upper coverage on beta(Hhat^c)       |0.9260         |0.9365          |
|md40 n700                                |field-s one-sided upper coverage on beta(Hhat^c)                   |0.9280         |NA              |

- Under `effMaxSG` at ε = 0.20 the selected subgroup is larger on 1881 / 1988 / 1844 / 1848 of the paired replicates and smaller on none; mean |Ĥ| is 112 / 153 / 107 / 118 patients against 72 / 75 / 72 / 73 under `maxeffCons` (`n_sel` and `n_harm` coincide in both campaigns).
- Sensitivity rises from 0.166 / 0.335 / – / 0.121 to 0.268 / 0.712 / – / 0.204 at a PPV of 0.412 / 0.798 / 0 / 0.410 against 0.398 / 0.757 / 0 / 0.397: the band rule recovers more of Q at about the same purity; in the null cell PPV is 0 by construction and sensitivity undefined (`n_true` = 0).
- The field's one-sided lower coverage on β(Ĥ) is 0.948 / 0.967 / 0.947 / 0.940 here against 0.951 / 0.948 / 0.950 / 0.947 under `maxeffCons`; the unstudentized complement's upper coverage 0.914 / 0.938 / 0.915 / 0.926 against 0.924 / 0.935 / 0.926 / 0.937, about 0.01 lower in three cells; field-s adds 0.000–0.003 to the unstudentized field.

### Regime diagnostics

|cell                                     | n_det| p_hat_mean| p_hat_lt05| sd_btc_naive| lamc_naive| lamc_s_naive| ij_sd_H| ij_sd_Hc|
|:----------------------------------------|-----:|----------:|----------:|------------:|----------:|------------:|-------:|--------:|
|md40 n500                                |  1998|      0.079|      0.999|        0.960|      0.977|        0.982|   1.696|    1.875|
|md120 n500                               |  2000|      0.089|      1.000|        1.087|      0.971|        1.000|   1.408|    1.639|
|null n500 (no subgroup; homogeneous +26) |  1993|      0.082|      0.999|        0.956|      0.980|        0.982|   1.724|    1.888|
|md40 n700                                |  1999|      0.080|      0.999|        0.985|      0.983|        0.987|   1.582|    1.876|

- p̂(Ĥ) is 0.079–0.089 with 99.9–100% of replicates below 0.5: a deeper tie regime than `mdf1`'s 0.155–0.189, on the same 1,842-candidate family with its duplicate-membership labels.
- SD(β̃ᶜ) / mean naive SEᶜ is 0.96–1.09 and λ-SDᶜ / naive SEᶜ 0.97–0.98 (field), 0.98–1.00 (field-s): the complement is in the survival n = 500 regime, not the moved regime where field-s was built to help.
- IJ SE / empirical SD is 1.41–1.72 on Ĥ and 1.64–1.89 on Ĥᶜ.

### The display (identity scale)

|cell                                     |block |estimator |    n| bias (MD)|     SD| mean SE|      b|     r| 1-sided cov| 1-sided ref| 2-sided cov| 2-sided ref|
|:----------------------------------------|:-----|:---------|----:|---------:|------:|-------:|------:|-----:|-----------:|-----------:|-----------:|-----------:|
|md40 n500                                |H     |naive     | 1998|    58.837| 16.640|  25.362|  3.536| 1.524|       0.137|       0.152|       0.293|       0.292|
|md40 n500                                |H     |mr        | 1998|    14.980| 19.590|  33.232|  0.765| 1.696|       0.975|       0.979|       0.996|       0.995|
|md40 n500                                |H     |fld       | 1998|     7.067| 21.235|  24.943|  0.333| 1.175|       0.948|       0.945|       0.975|       0.971|
|md40 n500                                |Hc    |naive     | 1998|   -16.798| 12.224|  13.208| -1.374| 1.080|       0.660|       0.657|       0.780|       0.771|
|md40 n500                                |Hc    |mr        | 1998|    -5.995| 12.674|  23.763| -0.473| 1.875|       0.994|       0.995|       0.999|       0.999|
|md40 n500                                |Hc    |fld       | 1998|    -3.851| 12.929|  12.902| -0.298| 0.998|       0.914|       0.910|       0.933|       0.939|
|md40 n500                                |Hc    |fld_s     | 1998|    -3.861| 12.903|  12.966| -0.299| 1.005|       0.914|       0.912|       0.936|       0.941|
|md120 n500                               |H     |naive     | 2000|    27.304| 18.310|  21.873|  1.491| 1.195|       0.672|       0.682|       0.774|       0.802|
|md120 n500                               |H     |mr        | 2000|    -4.147| 22.872|  32.214| -0.181| 1.408|       0.992|       0.994|       0.987|       0.993|
|md120 n500                               |H     |fld       | 2000|    -6.403| 25.088|  25.203| -0.255| 1.005|       0.967|       0.972|       0.916|       0.944|
|md120 n500                               |Hc    |naive     | 2000|   -11.239| 14.591|  14.157| -0.770| 0.970|       0.804|       0.795|       0.875|       0.867|
|md120 n500                               |Hc    |mr        | 2000|    -1.764| 15.395|  25.227| -0.115| 1.639|       0.997|       0.995|       0.999|       0.999|
|md120 n500                               |Hc    |fld       | 2000|    -0.472| 15.863|  13.745| -0.030| 0.866|       0.938|       0.919|       0.925|       0.910|
|md120 n500                               |Hc    |fld_s     | 2000|    -0.508| 15.862|  14.155| -0.032| 0.892|       0.941|       0.924|       0.930|       0.920|
|null n500 (no subgroup; homogeneous +26) |H     |naive     | 1993|    60.823| 16.415|  25.824|  3.705| 1.573|       0.116|       0.132|       0.262|       0.267|
|null n500 (no subgroup; homogeneous +26) |H     |mr        | 1993|    16.196| 19.461|  33.558|  0.832| 1.724|       0.974|       0.977|       0.994|       0.995|
|null n500 (no subgroup; homogeneous +26) |H     |fld       | 1993|     7.973| 21.177|  25.134|  0.377| 1.187|       0.947|       0.942|       0.975|       0.971|
|null n500 (no subgroup; homogeneous +26) |Hc    |naive     | 1993|   -16.505| 12.078|  13.135| -1.367| 1.087|       0.667|       0.664|       0.782|       0.778|
|null n500 (no subgroup; homogeneous +26) |Hc    |mr        | 1993|    -5.806| 12.551|  23.701| -0.463| 1.888|       0.996|       0.996|       0.999|       0.999|
|null n500 (no subgroup; homogeneous +26) |Hc    |fld       | 1993|    -3.690| 12.811|  12.866| -0.288| 1.004|       0.915|       0.914|       0.942|       0.942|
|null n500 (no subgroup; homogeneous +26) |Hc    |fld_s     | 1993|    -3.702| 12.786|  12.894| -0.290| 1.008|       0.917|       0.915|       0.942|       0.942|
|md40 n700                                |H     |naive     | 1999|    60.861| 16.969|  24.930|  3.587| 1.469|       0.085|       0.121|       0.208|       0.240|
|md40 n700                                |H     |mr        | 1999|    15.806| 19.727|  31.212|  0.801| 1.582|       0.958|       0.964|       0.986|       0.989|
|md40 n700                                |H     |fld       | 1999|     7.725| 21.252|  24.530|  0.363| 1.154|       0.939|       0.938|       0.959|       0.967|
|md40 n700                                |Hc    |naive     | 1999|   -11.591| 10.377|  10.789| -1.117| 1.040|       0.720|       0.723|       0.822|       0.821|
|md40 n700                                |Hc    |mr        | 1999|    -3.759| 10.629|  19.938| -0.354| 1.876|       0.994|       0.997|       0.999|       1.000|
|md40 n700                                |Hc    |fld       | 1999|    -2.161| 10.782|  10.611| -0.200| 0.984|       0.926|       0.922|       0.942|       0.941|
|md40 n700                                |Hc    |fld_s     | 1999|    -2.168| 10.769|  10.650| -0.201| 0.989|       0.928|       0.923|       0.943|       0.943|

- Field points on Ĥ sit within 0.006 of the Gaussian reference on the one-sided coverage (0.948 vs 0.945, 0.967 vs 0.972, 0.947 vs 0.942, 0.939 vs 0.938) with r 1.00–1.19 and b −0.26 to +0.38; the field's two-sided point at md120 (0.916 vs 0.944) is the largest departure on this block.
- Field-s points on Ĥᶜ: 0.914 vs 0.912, 0.941 vs 0.924, 0.917 vs 0.915, 0.928 vs 0.923 (b −0.03 to −0.30, r 0.89–1.01); the md120 point is 0.017 above its reference, the others within 0.005. IJ points are within 0.006 on both blocks.
- Figures: `fig_mdsgnb20_bias_coverage_display_H.png`, `fig_mdsgnb20_bias_coverage_display_Hc.png` (also embedded in the summary HTML).

## Scope

These are operating characteristics on one continuous design. They do not verify condition (A3) on the GLM paths, and no construction is promoted on this design's performance.

## Findings

1. **Checker halt, not a data halt** (Gate 2 brief above): `gate2.R:47`'s `sprintf` lacked an argument; fixed (`34580fb0`), cell 1 committed by hand with the runner's paths; the relaunched runner's ceiling counter restarted at zero (true cumulative 9,783 s, under the 19,071-s ceiling either way).
2. **Cost:** 9,783 s for 8,000 replicates at 63 workers (163 min against the 212-min projection); per replicate `fit_mr_secs` 56.2 / 68.0 / 53.9 / 84.9 s (field pass 31–42 s, complement 3.3–5.8 s) against `mdf1`'s 15.0 / 17.2 / 14.7 / 17.6 s at 13 workers on the M4 Max under `maxeffCons` — the band rule's field pass and the host both contribute (Stage 1 F6).
3. **Field-s equals the unstudentized field on this design** to within 0.003 in coverage at every cell and 0.005 at every ladder point: λ-SDᶜ / naive SEᶜ is already 0.97–0.98, and SD(β̃ᶜ) / naive SEᶜ 0.96–1.09, so the candidate-wise rescaling changes almost nothing. The moved-regime behaviour field-s was built for (survival p30sg, 1.08–1.12) does not occur here.
4. **Complement one-sided coverage is below nominal at n = 500** (0.914–0.917 for field-s in md40 and null, Wilson upper 0.925–0.929) and about 0.01 below `mdf1`'s on the same seeds; it is nearer nominal at md120 (0.941) and n = 700 (0.928). Reported, not interpreted.
5. **md120's field two-sided interval** covers at 0.916 while its one-sided lower bound covers at 0.967: the two-sided departure comes from the upper side of an asymmetric Λ*. The one-sided lower bound is the certified product; reported as a fact.
6. **The Stage 1 record listed `logs_mdsgnb20/` in full for deletion**, written before the campaign's raw render logs and Gate 2 outputs were placed there by the runner; §1.5 keeps raw render logs untracked under that directory, so the closeout deleted the smoke, calibration and dry-run items only and left the campaign's raw logs (untracked) in place.
7. **Summary extract, two rounds:** the first render's extract lacked the paired `mdf1` comparators, the grew/stayed/shrank counts and the joint margins that the record's tables print; the summary was amended and re-rendered (`fece16e6`, `4ab27707`) so the post-condition (every table number in the CSV) holds. An empty-bound branch (the null cell's harm oracle) was fixed in the same file before the first successful render.
8. The Gate 2 report's per-cell `meta` line shows `built_at` in local time (the template stamps `Sys.time()`), while the progress log is UTC.

## Commits of this task (`git log --oneline 0071c17e..HEAD`, oldest last; the catalog commits follow this record)

```
4ab27707 mdsgnb20 extract: 1,114 rows (1,086 mdsgnb20 + 28 paired mdf1 comparator rows); COLUMNS updated for campaign, margins and counts
fece16e6 mdsgnb20 summary: export the paired mdf1 comparators (campaign = mdf1), the grew/stayed/shrank counts and the joint margins so every number the record's tables print is in the extract; re-rendered
486249d2 mdsgnb20 extract: md_field_metrics.csv (1,050 rows: cell x block x estimator x metric, MC SEs, Wilson bounds, denominators, bundle commit) and COLUMNS_md_field.md
9bd7cb7a mdsgnb20 Stage 3 summary (transplant of summary_continuous_field_mdf1.qmd): Hhat and Hhat^c tables (field-s evaluated, unstudentized field beside), the D3 ladder with MC SEs and the oracle beside, the field-s Bonferroni pair, the paired rule contrast with mdf1, regime diagnostics, identity-scale display (field-s via a renamed copy); rendered HTML and the two display figures
2e4540c6 mdsgnb20: campaign complete (cumulative render wall 7662 s)
f528d188 mdsgnb20 md40_n700: md 40 n 700, 2000 replicates + combine; Gate 2 PASS (GATE_COUNTS run=58 passed=58 failed=0)
dbeee517 mdsgnb20 null_n500: md null n 500, 2000 replicates + combine; Gate 2 PASS (GATE_COUNTS run=58 passed=58 failed=0)
b12f5983 mdsgnb20 md120_n500: md 120 n 500, 2000 replicates + combine; Gate 2 PASS (GATE_COUNTS run=58 passed=58 failed=0)
76d6a8fd mdsgnb20: remove HALT (the halt was the checker's own error, fixed in 34580fb0; cell md40_n500 passed Gate 2 58/58 and is committed); resume
72e6f129 mdsgnb20 md40_n500: md 40 n 500, 2000 replicates + combine; Gate 2 PASS (GATE_COUNTS run=58 passed=58 failed=0) -- committed by hand after the checker fix, the runner's paths
34580fb0 gate2.R: the same-draws oracle label's sprintf lacked its tolerance argument ("too few arguments" after 52 of 53 checks passed on md40_n500); the check itself is unchanged
414467b3 mdsgnb20 md40_n500: HALT -- gate2.R failed: GATE_COUNTS run=53 passed=52 failed=1
e6abdb29 ACTG175 continuous field re-run Stage 1 record: install (upgrade = FALSE; Built 2026-09-16 05:57:14 UTC), the survival reference table, template edits E1-E5 and scripts (quoted sources), smoke -- zero flips on 4 cells (17/19/14/20 of 20 identical; 3/1/6/0 label ties), field-s wiring green, campaign rule reproduces the draws (n_true identical, oracle <= 2.1e-12; |Hhat| grew on 17 of 20); calibration W = 63 (31.3 reps/min, peak 72 GB), projection 212 min, ceiling 19071 s, timeout 3361 s; Gate 1 PASS, advance go noted
92d6ceb7 scripts_mdsgnb20 (TASK_md_field_rerun_2026-09-15 §1.5, transplants): mem_sampler.sh (scripts_mdf1, Linux ps); smoke_identity.R (gate2_check.R's pairing proof and classification, every column except *_secs, field-s wiring, rule mode); gate2.R (p12x20's gate2F.R pointed at the MD bundles with the §2.3 checks, same draws both directions vs mdf1); run_mdsgnb20.sh (run_p12x20.sh sequencing with run_cell.sh's render lines, GNU timeout, progress log, halt file, per-cell explicit-path commits)
2cffb95f MD template E1-E5 (TASK_md_field_rerun_2026-09-15 §1.4, transplanted from the survival m1 template): FS_MD_FOCUS / FS_MD_NBHD knobs with the six-focus guard, the band-foci eps guard and the _nb tag on the stem; selection_rule guard; FS_MD_FIELD_SCALEC knob (default selected) passed as field_scale_complement; nine fld_Hc_*_s and nine fld_joint_s_* recorder columns with their fill blocks; _s interval invariant pairs; effect_neighborhood / field_scale_complement in meta and effect_neighborhood / selection_rule / field_scale_complement in the poolability keys; knob echo. Defaults reproduce mdf1's stem and settings.
bb84120e Add TASK_md_field_rerun_2026-09-15 as received
<this record, then scripts_mdsgnb20 catalog files + status_curated.md, then current_status.md alone>
```
