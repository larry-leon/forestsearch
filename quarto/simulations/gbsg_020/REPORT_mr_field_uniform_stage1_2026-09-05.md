# REPORT — Field interval, uniform calibration: Stage 1 (Implementation, identities, projection)

**Task:** `dev/tasks/TASK_mr_field_uniform_2026-09-05.md` (1313f597); Stage 0 record 2291ca02. Implementation commit **b1316d22**.
**Date:** 2026-09-05. Executor: Claude Code. H1–H5 at defaults. **Parked at Gate 1 per instruction — H6 (compute go) is Larry's decision on this record.**

---

## GATE 1 status: identities 46 of 47 within stated tolerances (the one exception argued below to be reference-noise-limited); κ-sweep cost and projection measured at 100 workers under load; mcse(κ*) at H4 counts reported. **Stopped here for H6, with H3/H4 adjustment data as requested.**

## 1a — Implementation (b1316d22)

- `R/fs_mr_field_uniform.R`: new internal `fs_mr_field_uniform()` implementing task §method steps 1–6 — top-M mass-carrying reduction (winner always kept), winner-profile family μ_δ, per-δ hypothetical trials with Gaussian multipliers of covariance exactly Σ̂ (via `crossprod(B_red, ξ)`, or a Cholesky root of a supplied `Sigma` for the reference checks), the gate's own selection map (vectorized thresholded argmax under `reselection = "maxeff"` — proven equal to `.fs_mr_select` on 4,000 draws — per-draw fallback otherwise), conditioning on detection, the field procedure per hypothetical trial, C₁/C₂ profiles, κ* by grid with mcse from the binding profile's binomial SE over the local slope.
- `fs_mr_inference()` gains add-only `field_uniform = FALSE`; when TRUE the κ block attaches as `field$uniform` (kappa, kappa_mcse, M, mass_covered, minC1, profiles, `lower_2u`/`upper_2u` = the κ*-widened interval for the real analysis per task step 6, timing) under seed `+ 910000L`, drawn strictly after the field's stream.
- `forestsearch()` forwarding line (Larry's classification: add-only pass-through, no behaviour change); template `FS_S7_UNIFORM` knob + the eight recorder columns (`fld_H_kappa/M/mass/minC1/lo2u/hi2u/kappa_mcse/uniform_secs`), meta records `field_uniform`; harness `mv_mr(field_uniform = FALSE)` pass-through.
- **One implementation lesson recorded:** per-trial fresh perturbation batches are load-bearing. A first draft shared the outer/inner batches across trials within a δ (an over-extension of the reference's inner-draw economy); the K = 1 identity check exposed it immediately (coverage wandering 0.88–0.97 where binomial error is ±0.006). The committed version draws all batches fresh per trial (inner batches still shared across that trial's outer draws, the reference economy), after which K = 1 is exactly nominal.

## 1b — Identities

**Default-path byte-identity (8/8 PASS):** the three fixed-seed Guo–He cases under `"ij"` are `identical()` to the pre-field-feature 736088e3 reference (timing excluded); the same three under `"field"` with `field_uniform = FALSE` are `identical()` to pre-uniform captures taken minutes before the change; under `field_uniform = TRUE` the field element is unchanged (uniform sub-block stripped, `identical()`); same-seed gate calls give identical uniform blocks. Suite at standing parity 0 fail / 4941 pass / 3 skip / 32 warn. Template level: default-path re-renders reproduce the committed smoke bundles on **73/73 shared columns, both cells**, and the eight new columns are all-NA there.

**Limit-experiment references (22/23 within stated tolerances):**

- K = 2, Σ = I, pure argmax (task's rule-A targets), R_rep = 4,000: C₂(δ;1) = 0.9475/0.9343/0.9050/0.9060/0.9277 at δ = 0/1/2/2.5/4 vs 0.945/0.93/0.905/0.90/0.935 (tol 0.015, all pass); C₁ ∈ [0.941, 0.952] at every δ (band [0.93, 0.96]); **κ\* = 1.25** vs 1.20 ± 0.05 (pass at the boundary; mcse 0.024); C₂(2.5;1.2) = 0.9487 vs 0.947 ± 0.015.
- K = 10, ρ = 0.7: all eight C₂ points (κ = 1 and 1.2) within 0.025 of the CSV; C₁ in band.
- K = 10 disjoint: δ = 0/1/2 within 0.025 (0.9625/0.9421/0.8562 vs 0.966/0.935/0.840); C₁ in band. **δ = 3: 0.7913 vs 0.766 ± 0.025 — outside by 0.0003.** Convergence analysis: three independent estimates (R_rep = 2,400/2,400/800, binomial SE 0.008/0.008/0.015) give 0.7829/0.7921/0.7975 — a converged value ≈ 0.785–0.79, which is ~1.5 combined-SE from the CSV's own R_rep = 800, shared-inner-draw estimate of 0.766. Every other disjoint point and the entire ρ = 0.7 sweep match; the flagged point sits on the steepest part of the coverage curve where the reference's own MC error spans the tolerance. Recorded as reference-noise-limited, not an implementation discrepancy — **Larry adjudicates at H6.**
- K = 1: C₁ 0.933–0.953 and C₂(·;1) 0.934–0.952 at every δ (nominal ± 0.02); κ* = 1.07 (within its coarse MC error of 1).
- Top-M: K = 12 with two candidates at 0.99 of p̂ — reduced M = 2, mass 0.990; κ* full 1.15 vs reduced 1.18, within MC SE (0.068).
- Seed reproducibility and the vectorized-S equivalence: exact.

**1c smoke (10/10 PASS):** 5 replicates per cell at the committed seeds, `field_uniform = TRUE` (campaign `usmoke`): every non-uniform column identical to the s7 bundles (64/64, FB-join and wall-clock columns excluded); uniform block finite on all detected reps (3 h100 / 5 h175); the κ*-widened interval contains the plain quantile interval on every replicate. Committed evidence bundles: `..._usmoke_res_1_5.rds` (both cells), `..._uburst_res_1_100.rds` (both cells).

## 1d — Cost and projection (100 workers, loaded, campaign `uburst`, 100 replicates/cell)

**κ-sweep cost, separately from the gate cost (Larry's request):**

| Cell | κ sweep per detected rep (mean / q50 / q90) | Gate ex-κ (search + MR-IJ + field) | Detections | Wall, 100 reps @ 100 workers |
|---|---|---|---|---|
| h100 | **159 / 159 / 165 s** | ~35 s (s7 anchor 28.3) | 66/100 | 220 s |
| h175 | **246 / 246 / 266 s** | ~40 s (s7 anchor 41.7) | 95/100 | 333 s |

Lightly loaded (5 workers, smoke): 78–83 s/rep — load inflation ≈ ×2–3, as with the field block. Peak host memory during the bursts: 119 GB (baseline 33; budget 176). Zero κ NAs on 161 detected replicates.

**mcse(κ\*) at H4 counts (R_rep = 300):** mean ≈ **0.10** on the real forestsearch families (both cells; synthetic test families 0.03–0.07). Against that noise, the observed per-replicate κ* spread is real but coarse: h100 mean 1.57 (SD 0.20, range 1.15–**2.00**), h175 mean 1.49 (SD 0.15, range 1.23–**2.00**) — some replicates hit the κ grid ceiling of 2.

**H3 flag:** the mass-carrying target (≥ 0.99) is unreachable under M ≤ 12 on these enumerated families — **M = 12 (the cap) on every detected replicate, mass covered 0.52–0.96 (means 0.75 / 0.79)**. The protection family is therefore evaluated on a 12-candidate core holding, on some replicates, barely half the re-selection mass. minC₁ stayed 0.893–0.939. H3/H4 adjustments to consider before H6: raise M_cap (cost is roughly linear in the big-matrix width; ~×2 for M = 24), raise R_rep (mcse ∝ 1/√R_rep; ×4 for mcse ≈ 0.05), and/or extend the κ grid past 2 for the ceiling-hitters.

**Projection for Stage 2 as scoped by H5 (at current H4 defaults):**
- 2a (h100 + h175, sim_id 1–2,000, campaign `s7u`): from measured throughput, h100 ≈ 2000 × 2.2 s ≈ 1.2 h, h175 ≈ 2000 × 3.3 s ≈ 1.85 h, + combines ≈ **~3.2 h wall**, ~310 core-h.
- 2b (four Guo–He cells, smaller families so cheaper per rep, est. **~1–1.5 h**) — total Stage 2 ≈ **4–5 h wall** at 100 workers. Within any plausible ceiling, but the H3/H4 choices above change it multiplicatively.

## Deviations

- The committed `limit_sweep_K2_rules_2026-09-05.csv` contains rules D/E/F only; the rule-A tolerance targets used are the task text's inline values (noted for the record).
- The K2 reference computes m(·) in closed form; this implementation (and the K10 reference) uses inner Monte Carlo — the stated tolerances absorbed it everywhere except the flagged K10-disjoint δ = 3 point discussed above.

**Parked at Gate 1.** Awaiting H6 (compute go, possibly with H3/H4 adjustments) and the K10-δ=3 adjudication.
