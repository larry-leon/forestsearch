# REPORT — MR gate vs Guo & He (2021): Stage 1 (adapter, identities, smoke test)

Date: 2026-09-04. Task: `dev/tasks/TASK_mr_vs_guohe_2026-09-04.md` + addendum A (`276c3cee`).
Repo state: branch `feature/glm-extension`, on top of Stage 0 (`46cb9602`), forestsearch 0.3.5.
Decisions applied: D2/D3 defaults; **D6 as directed** (add-only `return_reselection = FALSE` on `fs_mr_inference()`).

## GATE 1: GREEN on all checks — STOPPED FOR THE COMPUTE GO (stop-to-ask by design)

Every 1b and addendum-A1 identity passed with the measured values below; the 1c smoke test ran with the naive identity holding bit-exactly on all five replicates; the existing test suite passes with the D6 change (0 failures). **The Stage 2 projection is 20.0 core-hours → ~10.0 minutes wall at 120 workers.** Awaiting Larry's go on that projection before any Stage 2 compute.

---

## 1a — Adapter and the D6 change

**R/ change (D6, add-only).** `fs_mr_inference()` gains `return_reselection = FALSE` as its **last** parameter (`R/fs_mr_inference.R:376`; all production call sites pass named arguments, checked at `R/fs_mr_inference_methods.R:127` and `R/forestsearch_main.R:3367`, so trailing position also protects any positional caller). When `TRUE`, the already-computed per-draw `winner` vector and its tabulation `p_hat = tabulate(winner)/draws` are appended as a `reselection` element (`R/fs_mr_inference.R:618-623`); the list construction itself is unchanged (`out <- list(...)` at :597 is the former literal). Docs at :290-298 and :340-341; `man/fs_mr_inference.Rd` regenerated (only that Rd changed under `devtools::document()`). Nothing in the arithmetic moves; the default return object is the previous one exactly (verified V8a below by exact name-vector comparison).

**R/ classification:** adds code; moves nothing; changes no existing behaviour; changes no method.

**Test suite on the edited source:** `devtools::test()` → **0 fail / 4941 pass / 3 skip / 32 warnings** (log in session scratchpad). The memory-recorded parity point (31 warnings, 4864 passing) predates the subgroup-sims wrapper phases; the 32nd warning comes from `tests/testthat/test-subgroup-sims.R:46` ("no non-missing arguments to max"), a file last touched at `a8430e73` (2026-09-02) — before this change, unrelated to it. All MR fixture tests (`test-mr-closedform-fixture.R`) pass.

**Adapter.** `quarto/GuoHe/mr_vs_guohe_sim.R` (new, beside the replication harnesses, which are untouched). Per the task, the gate accepts an external family, so nothing else was added under `R/`. Following the §5.1/§5.2 repository convention the adapter lives in a sourced sim file rather than inline in the qmd — the Stage 3 comparison qmd will source it, targets never retyped (recorded as a convention choice, not a deviation of substance). Contents:

- Replicate regeneration from the bundles' seed schemes (`mv_gh51_regen`, `mv_gh52_regen`; bases re-derived and asserted against the stored `seed_base`).
- Candidate families as named row-index lists (`mv_cand_idx_51/52`) — the gate's input format.
- **Orientation (§5.1 only):** `treat_flip = 1 - treat` appended to the regenerated data; the treatment-only Cox partial likelihood is symmetric under the flip, so the raw coefficient on flipped data equals the replication's oriented effect and `maxeff` is its argmax. §5.2 needs no adapter (`orient = +1`, raw-coefficient argmax already). Fidelity measured in V6 below.
- **MR call (D2/D3 defaults):** `reselection = "maxeff"`, `admission = list(effect_floor = NULL, consistency = NULL)` (the engine's unrestricted path), `draws = 5000`, `multiplier = "poisson"` (centred), `ci_method = "ij"`, `return_reselection = TRUE`. MR seed = `seed_base + m + 700000L` (disjoint from the data seed `base+m` and the replication's bootstrap seed `base+m+500000`).
- Addendum-A quantities: `mv_sigma` (Σ̂ = B_effᵀB_eff via the engine's own `.fs_mr_assemble`, read-only internal use), `mv_a6`, `mv_c_of_M` (adaptive quadrature), `mv_m0_hat` (chunked MC, fixed seed, 2e5 draws, 1e-8 ridge on the near-singular nested correlation — recorded), `mv_M_eff` (uniroot on [1, K+1]), `mv_a_record` (the full A3 column set).

Note: `fs_mr_inference` is not exported; the harness calls `forestsearch:::fs_mr_inference` / `:::.fs_mr_assemble` (read-only use of internals, stated here per the engine-level-entry-point caveat in its docs).

## 1b — Identities, measured values (all PASS; log: `mr_vs_guohe_sim.R` executable block, exit 0)

| ID | Identity | Tolerance | Measured |
|---|---|---|---|
| V1sum (§5.2, 154 cands, t7_b2_00 m=1) | influence sums to zero within every candidate | \|Σdb\|/σ_D < 1e-6 | worst 1.25e-08 |
| V1var (§5.2) | σ_D² = sum(dfbeta²) equals robust (sandwich) coxph variance | rel < 1e-6 | worst 1.24e-15 |
| V1sum (§5.1, 12 cands, t6_k12 m=1) | same | < 1e-6 | worst 5.44e-09 |
| V1var (§5.1) | same | rel < 1e-6 | worst 1.05e-15 |
| V2a (K=1) | de-biased = naive up to MC error at B=5000 | 4 sd = 0.0203 (log) | −0.0018 |
| V2b (K=1) | V̂_IJ/σ_D² in the separated-regime band | [3.6, 4.4] | 4.130 (`ij_source = "ij"`) |
| V3 | exchangeable candidates → exchangeable selection | χ² p > 0.01 | counts 51/43/40/34/32 over 200 reps, p = 0.219 |
| V4 (ρ=0.3) | M_eff exchangeable target (A1) | \|Δ\| < 0.1 | m0 1.2866 (MC se 0.0016; analytic 1.2874), M_eff 6.21 vs 6.22 |
| V4 (ρ=0.6) | 〃 | < 0.1 | m0 0.9727 (se 0.0019; analytic 0.9732), M_eff 3.65 vs 3.65 |
| V4 (ρ=0.9) | 〃 | < 0.1 | m0 0.4872 (se 0.0022; analytic 0.4866), M_eff 1.80 vs 1.80 |
| V5 | R̂ = I, K=10 → M_eff = K | \|Δ\| < 0.2 | m0 1.5390 (c(10) = 1.5388), M_eff 10.00 |
| V6 | treat-flip negates β̂(g) and B_eff columns, Σ̂ invariant (A1) | rel < 1e-10 | 5.6e-15 / 2.1e-15 / 8.5e-16 |
| V7a (§5.1 disjoint) | off-diagonal Σ̂ exactly 0; A6_std = 1 + 2p̂_Ĥ | exact; < 1e-12 | max\|offdiag\| = 0; A6_std 1.551 |
| V7b (§5.2 nested) | A6_mass_std ≥ 1 (structural, measured) | ≥ 1 | 2.946; min winner cross-term +2.29e-02 (all non-negative on this replicate) |
| V8a (D6) | default return object unchanged | exact name vector | PASS (18 names, no `reselection`) |
| V8b (D6) | winner length = draws; Σp̂ = selection_rate | 1e-12 | 1.000000 = 1.000000 |

**Interpretation recorded for V3.** The task's "all K candidates identical" is implemented as identical **in distribution** (K = 5 disjoint equal-probability candidates, all β = 0), with exchangeability tested **across** 200 independent replicates on the modal MR winner. Literally identical membership vectors give identical perturbed effects on every draw, and `which.max` resolves ties to the first index deterministically — that variant tests tie-breaking, not exchangeability, and cannot pass; within one replicate the frequencies legitimately concentrate on the realized argmax. Stated in the harness comment at the V3 block.

**(A6) status, as the addendum directs:** structurally satisfied on both replication families rather than empirically tested — §5.1 off-diagonals are exactly zero (disjoint dfbeta supports), §5.2 winner cross-terms measured non-negative; V7 records both.

## 1c — Smoke test (t7_beta2_00, replicates 1–5, paired by seed)

- **S1a PASS:** recomputed naive `point`, `lower`, `c_hat`, `gamma_s` are `identical()` (bit-exact) to the stored bundle values on all 5 replicates.
- **S1b PASS:** recomputed full-sample argmax cutpoint equals stored `c_hat_gh` on all 5 — the selection MR receives is the replication's.
- **S1c PASS:** MR ran on every replicate with finite estimate and SE, `ij_source = "ij"` throughout.
- G&H per-replicate values were taken from the **stored** columns (debiased = `r4_bias + gamma_s`, lower = `gamma_s − r4_dist`); no G&H re-run, per Stage 0's D5 finding.

Paired values (log-HR scale, r = 1/30 column for G&H):

| m | γ_s | naive | G&H debiased | G&H lower | MR debiased | MR lower(1s) | MR se_ij | bias_sel | p̂_Ĥ | t_MR (s) |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 0 | 0.207 | 0.086 | −0.216 | 0.151 | −0.389 | 0.328 | 0.056 | 0.30 | 1.23 |
| 2 | 0 | −0.093 | −0.198 | −0.508 | −0.194 | −0.841 | 0.394 | 0.098 | 0.26 | 1.00 |
| 3 | 0 | 0.089 | −0.031 | −0.350 | 0.012 | −0.625 | 0.388 | 0.081 | 0.13 | 1.32 |
| 4 | 0 | 0.404 | 0.300 | −0.009 | 0.354 | −0.226 | 0.353 | 0.047 | 0.35 | 1.03 |
| 5 | 0 | 0.097 | −0.028 | −0.358 | −0.003 | −0.566 | 0.342 | 0.098 | 0.17 | 1.30 |

(Direction sanity only, no conclusions at n = 5: MR de-biases toward the null on every replicate, less aggressively than G&H at r = 1/30, with wider one-sided bounds.)

Addendum-A record columns computed on replicate 1 (schema viability): p̂_Ĥ 0.299; top-3 0.299/0.211/0.082 (cand_100/cand_098/cand_073); Σ̂_ĤĤ 0.02851; A6_mass 0.08403 (std 2.948); m̂(0) 0.6784 (MC se 0.00211 at 2e5 draws); **M_eff 2.34 of K = 154**; tie_resid_implied 0.0335. (Same machinery gives M_eff 11.94 of K = 12 on a t6_k12 replicate and 2.00 of K = 2 on t35 — the effective-competition reading behaves as the supplement intends on both family types.)

## Projection (the Gate 1 stop-to-ask object)

Per-replicate Stage 2 cost = regeneration + naive recompute + MR (B = 5000) + addendum-A record, measured single-core on this host:

| Cell type | per replicate | of which MR / A-record |
|---|---|---|
| t35 (k=2, n=400) — also t6_k02 | 0.11 s | 0.09 / 0.02 |
| t6_k12 (k=12, n=2400) | 0.97 s | 0.84 / 0.13 |
| t7 (~151 cands, n=400) | 5.50 s | 1.19 / 4.31 |

t6 k=6/k=10 interpolated linearly in k between the measured k=2 and k=12 (per-subgroup size is constant at 200, so cost is ~linear in k — the §5.1 RUN note's measured scaling).

**Full grid (D1 default, 16 cells × 2000 replicates): 20.0 core-hours → ~10.0 minutes wall at 120 workers.** The dominant term is the §5.2 A-record M_eff Monte Carlo (2e5 draws on a ~151×151 correlation; the addendum's floor is 1e5 draws, so this could halve to ~13 core-h if wanted — not assumed). No G&H compute is included (none is needed; D5 moot).

## Files and commits

- `R/fs_mr_inference.R` + `man/fs_mr_inference.Rd` — the D6 change (committed separately).
- `quarto/GuoHe/mr_vs_guohe_sim.R` — adapter + identities + smoke + timing, executable and fail-loud (exit 0 on this run).
- This report.

## Deviations

- V3 interpretation as recorded above (identical-in-distribution, across-replicate exchangeability).
- Adapter in a sourced sim file beside the qmd rather than inline in the qmd, per the §5.1/§5.2 convention; the Stage 3 qmd will source it.
- None otherwise: no replication file touched, no bundle re-run, no method change.

— End of Stage 1. **STOPPED at Gate 1 for the compute go on the 20.0 core-h / ~10 min projection.**
