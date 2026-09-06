# REPORT — Field interval, uniform calibration: Stage 3 (final)

**Task:** `dev/tasks/TASK_mr_field_uniform_2026-09-05.md` (1313f597). Records: Stage 0 2291ca02 · Stage 1 28022b54 · Gate 1 adjudication 5577b079 · H2/H3 adjustments 3ad6ee13 · 2b driver 8975af1c · Gate 2 record beside this file.
**Date:** 2026-09-06. Data: campaign `s7u` (h100/h175, sims 1–500, `field_uniform = TRUE`, one gate call per replicate yielding MR (IJ), MR (field) and the uniform block) and campaign `uniform` (t35 β₂ = 0.3/0.5, t6 k = 10, t7 β₂ = 0.0 vs the committed 2026-09-05 field bundles). Rendered documents: `fs_maxeffCons_fb_mr_field_m1_{h100,h175}_knoise0_n500_s7u_{batch,combine}_1_500.html`.

Conventions: log scale throughout (bias, SD, SE, half-widths); coverage of β(Ĥ) (s7u) or θ (Guo–He) with Wilson intervals; the four rows are the task's Stage 3 layout — MR (IJ), MR (field, one-sided), MR (field, plain two-sided), MR (field, uniform two-sided); the field rows share the est2 point estimate, so their retained bias/SE-SD lines repeat by construction. κ*, M, mass, minC₁, mcse(κ*) recorded per replicate (H3/H4 as adjudicated).

## Per-cell tables

### h100 (campaign s7u, sims 1–500): 329 detections; target β(Ĥ)

| Row | Coverage (Wilson) | Retained bias (log) | SE/SD | Half-width / z·SE (log) |
|---|---|---|---|---|
| MR (IJ) | 0.979 (0.957, 0.990) | +0.2382 | 1.66 | 0.771 |
| MR (field, one-sided) | 0.957 (0.930, 0.974) | +0.1362 | 1.32 | — |
| MR (field, plain 2s) | 0.988 (0.969, 0.995) | +0.1362 | 1.32 | 0.695 |
| MR (field, uniform 2s) | 1.000 (0.988, 1.000) | +0.1362 | 1.32 | 1.124 |

κ*: mean 1.62, SD 0.14, quartiles 1.52/1.61/1.72, range [1.15, 2.12], ceiling(3) hits 0 of 329; mcse(κ*) mean 0.096; M mean 39.8 (cap-40 hits 323), mass ≥ 0.99 on 8 of 329; κ-sweep secs mean 583.

### h175 (campaign s7u, sims 1–500): 477 detections; target β(Ĥ)

| Row | Coverage (Wilson) | Retained bias (log) | SE/SD | Half-width / z·SE (log) |
|---|---|---|---|---|
| MR (IJ) | 0.964 (0.944, 0.978) | -0.0081 | 1.17 | 0.765 |
| MR (field, one-sided) | 0.962 (0.941, 0.976) | -0.0725 | 0.94 | — |
| MR (field, plain 2s) | 0.881 (0.848, 0.907) | -0.0725 | 0.94 | 0.704 |
| MR (field, uniform 2s) | 0.985 (0.970, 0.993) | -0.0725 | 0.94 | 1.093 |

κ*: mean 1.56, SD 0.13, quartiles 1.47/1.55/1.63, range [1.07, 2.06], ceiling(3) hits 0 of 477; mcse(κ*) mean 0.093; M mean 39.0 (cap-40 hits 431), mass ≥ 0.99 on 51 of 477; κ-sweep secs mean 651.

### t35_beta2_03 (campaign uniform, 2000 replicates; K = 2–2; target θ)

| Row | Coverage (Wilson) | Retained bias (log) | SE/SD | Half-width / z·SE (log) |
|---|---|---|---|---|
| MR (IJ) | 0.998 (0.995, 0.999) | +0.0045 | 1.59 | 0.624 |
| MR (field, one-sided) | 0.936 (0.925, 0.946) | -0.0030 | 0.87 | — |
| MR (field, plain 2s) | 0.902 (0.888, 0.914) | -0.0030 | 0.87 | 0.358 |
| MR (field, uniform 2s) | 0.959 (0.949, 0.967) | -0.0030 | 0.87 | 0.452 |

κ*: mean 1.26, SD 0.08, quartiles 1.22/1.27/1.31, range [1.03, 1.54], ceiling(3) hits 0 of 2000; mcse(κ*) mean 0.075; M mean 1.9 (cap-40 hits 0), mass ≥ 0.99 on 2000 of 2000; κ-sweep secs mean 76.
Runtime: gate+field+κ per replicate mean 77.1 s; cell elapsed 26.1 min at 100 cores.

### t35_beta2_05 (campaign uniform, 2000 replicates; K = 2–2; target θ)

| Row | Coverage (Wilson) | Retained bias (log) | SE/SD | Half-width / z·SE (log) |
|---|---|---|---|---|
| MR (IJ) | 0.998 (0.994, 0.999) | -0.0129 | 1.64 | 0.662 |
| MR (field, one-sided) | 0.942 (0.930, 0.951) | -0.0113 | 0.89 | — |
| MR (field, plain 2s) | 0.900 (0.887, 0.913) | -0.0113 | 0.89 | 0.370 |
| MR (field, uniform 2s) | 0.944 (0.933, 0.953) | -0.0113 | 0.89 | 0.453 |

κ*: mean 1.23, SD 0.10, quartiles 1.13/1.24/1.30, range [1.03, 1.51], ceiling(3) hits 0 of 2000; mcse(κ*) mean 0.074; M mean 1.6 (cap-40 hits 0), mass ≥ 0.99 on 2000 of 2000; κ-sweep secs mean 74.
Runtime: gate+field+κ per replicate mean 75.9 s; cell elapsed 25.7 min at 100 cores.

### t6_k10 (campaign uniform, 2000 replicates; K = 10–10; target θ)

| Row | Coverage (Wilson) | Retained bias (log) | SE/SD | Half-width / z·SE (log) |
|---|---|---|---|---|
| MR (IJ) | 0.997 (0.993, 0.999) | +0.0884 | 1.57 | 0.437 |
| MR (field, one-sided) | 0.944 (0.933, 0.953) | +0.0382 | 1.09 | — |
| MR (field, plain 2s) | 0.970 (0.961, 0.976) | +0.0382 | 1.09 | 0.325 |
| MR (field, uniform 2s) | 0.998 (0.994, 0.999) | +0.0382 | 1.09 | 0.509 |

κ*: mean 1.56, SD 0.09, quartiles 1.50/1.56/1.62, range [1.21, 1.92], ceiling(3) hits 0 of 2000; mcse(κ*) mean 0.096; M mean 7.9 (cap-40 hits 0), mass ≥ 0.99 on 2000 of 2000; κ-sweep secs mean 237.
Runtime: gate+field+κ per replicate mean 240.6 s; cell elapsed 81.6 min at 100 cores.

### t7_beta2_00 (campaign uniform, 2000 replicates; K = 122–179; target θ)

| Row | Coverage (Wilson) | Retained bias (log) | SE/SD | Half-width / z·SE (log) |
|---|---|---|---|---|
| MR (IJ) | 1.000 (0.998, 1.000) | +0.0449 | 1.91 | 0.688 |
| MR (field, one-sided) | 0.940 (0.928, 0.949) | +0.0298 | 0.99 | — |
| MR (field, plain 2s) | 0.954 (0.944, 0.962) | +0.0298 | 0.99 | 0.359 |
| MR (field, uniform 2s) | 0.975 (0.967, 0.981) | +0.0298 | 0.99 | 0.414 |

κ*: mean 1.15, SD 0.06, quartiles 1.11/1.14/1.19, range [1.01, 1.39], ceiling(3) hits 0 of 2000; mcse(κ*) mean 0.068; M mean 39.5 (cap-40 hits 1872), mass ≥ 0.99 on 141 of 2000; κ-sweep secs mean 560.
Runtime: gate+field+κ per replicate mean 567.7 s; cell elapsed 192.7 min at 100 cores.


## Reading criteria (Larry's, from the task — factual mapping)

- *"Uniform two-sided coverage of β(Ĥ) ≥ 0.93 in every cell including the t35 dip cells and h175."* **Met in all six cells**: h100 1.000, h175 **0.985** (plain was 0.881), t35_03 **0.959** (plain 0.902), t35_05 **0.944** (plain 0.900), t6_k10 0.998, t7 0.975. Minimum 0.944.
- *"One-sided unchanged."* **Met**: 0.936–0.962 across all six cells — the same range as the plain field runs (the field block reproduced the committed bundles exactly, so the one-sided indicators are identical by construction on the Guo–He cells).
- *"κ* on forestsearch families reported as found."* Means 1.56–1.62 (s7u cells), 1.15–1.56 (Guo–He); no replicate reached the extended grid ceiling of 3; mcse(κ*) ≈ 0.07–0.10 at the H4 counts.

## Findings (in the record)

1. **The κ(Σ̂) widening repairs the two-sided dips it was built for**: h175 0.881 → 0.985, t35 dips 0.90 → 0.94–0.96, at the price of ~26–56% wider intervals (half-width 0.70 → 1.09 log on h175; 0.36 → 0.45 on t35). Where the plain interval already covered (h100, t6_k10), the widening lands conservative (1.000 / 0.998).
2. **Cost is the method's main liability on enumerated families**: the κ sweep at cap 40 runs 583–651 s per detected replicate under load on the forestsearch cells (vs ~15 s for the field block itself) and 560 s on t7's nested family; the t35/t6 families are cheap (74–237 s). The 6 h projection was missed (7 h 15 m actual) chiefly through this term.
3. **H3 remains cap-bound on forestsearch families**: even at M_cap 40, the ≥ 0.99 mass target is reached on only 8/329 (h100) and 51/477 (h175) replicates (mass means 0.93–0.95, min 0.78) and on 141/2000 of t7; the t35/t6 families cover ≥ 0.99 fully. minC₁ stayed ≥ 0.88 everywhere.
4. **Uniform-vs-IJ width**: on h100 the uniform interval (half-width 1.124) is wider than IJ (0.771) — in the tie regime IJ's conservatism already over-covers, and κ widening from the field's shorter base overshoots past it. On h175 the uniform half-width (1.093) is comparable to IJ (0.765 at coverage 0.964 vs 0.985). The uniform interval's value proposition is the *guarantee*, not width.
5. Default-path integrity held end to end: 2a non-uniform columns 64/64 identical to s7; 2b pairing proof exact on all 8,000 replicates.

## H7 (Larry's call on this record, per the task)

The one-sided bound as the standard reported quantity is supported everywhere tested (0.936–0.981 across this task and the OC map, with the OC map's c3 exception at 0.920 noted in its own record). For the two-sided choice — κ(Σ̂) at ~10 min/replicate compute with a uniform guarantee over the winner-profile family, vs IJ at zero extra cost with SE/SD 1.1–1.9 conservatism that under-covers in the moderate-harm large-n regime (OC map c5) — the record now carries both sides' numbers. No further task proposed.

## Done means (task §Done)

Stage 3 report + rendered documents committed; Gate 2 records beside them; `field_uniform` landed add-only with default-path identities recorded (Stage 1 + Gate 2); the four reference files in `dev/tasks/` (1313f597); template and harness updated under their own stems; branch unpushed.
