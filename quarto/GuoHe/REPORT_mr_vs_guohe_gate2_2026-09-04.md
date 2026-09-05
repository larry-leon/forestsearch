# REPORT — MR gate vs Guo & He (2021): Gate 2 record

Date: 2026-09-04. Task: `dev/tasks/TASK_mr_vs_guohe_2026-09-04.md` + addendum A. Compute GO given on the Stage 1 projection.
Driver: `quarto/GuoHe/mr_vs_guohe_run.R` (sources the committed Stage 1 adapter `mr_vs_guohe_sim.R`). Settings as authorized: 16 cells × 2000 replicates, MR draws = 5000 centred Poisson, IJ interval, M_eff Monte Carlo at 2e5 draws (fixed seed 20260904).

## GATE 2: GREEN

## Precondition (ran before any cell; stop-on-failure)

5 paired replicates each of `t35_beta2_03` and `t6_k12`: regenerated naive column `identical()` to the stored bundle values (`naive_beta_s`, `naive_cover`, `naive_dist`, `naive_bias`) and recomputed argmax equal to the stored `naive_sel` on all 10 — **PASS** (selected indices 1,1,2,1,1 and 9,1,10,8,6, all matching). Transcript in the driver output; the check also re-runs at the head of every driver invocation.

## Per-cell tally (from the `mr_vs_guohe_<id>.rds` gate blocks)

All 16 cells identical on every criterion:

| Criterion | Result |
|---|---|
| Replicate completeness | 2000/2000 in every cell (32,000 total), 0 worker errors |
| Naive column `identical()` to stored replication values | 0 mismatches in 32,000 (asserted per replicate inside the worker; §5.2 additionally asserts the recomputed argmax cutpoint and `gh52_truth_at` value against stored `c_hat_gh` / `gamma_s`) |
| MR NA count | 0 in 32,000 (finite estimate, IJ SE, and lower bound everywhere; `ij_source = "ij"` share reported per scenario in the qmd) |
| MR-side argmax agreement with the replication's selection | 32,000/32,000 (the treat-flip/§5.2-raw argmax equals the replication's argmax on every replicate) |
| G&H re-run | none (per-replicate G&H values are exact arithmetic on the stored bundles; the "re-run matches committed" clause is therefore moot, per Stage 0/D5) |

The Stage 3 qmd refuses to render on any completeness/naive-identity violation (load-time guard in its setup chunk).

## Cost — actual vs authorized projection (deviation, reported plainly)

- **Authorized projection:** 20.0 core-h, ~10 min wall at 120 workers.
- **Actual:** 89.5 min wall, ≈179 core-h at 120 workers (per-cell wall in each bundle's `elapsed_sec`; t35 ≈ 0.8–1.1 min, t6 rising to ≈ 4 min at k = 12, t7 ≈ 12.8 min each).
- **Cause:** the single-core per-replicate measurement (5.5 s on t7) did not scale under 120-way concurrency — the addendum-A M_eff Monte Carlo is a memory-bandwidth-bound dense-matrix workload (chunks of 5e4 × ~151 products), and 120 concurrent instances contend for bandwidth; the worker-measured per-replicate cost inflated to a mean 41.4 s on t7 cells. Nothing about the run's *scope* changed: same cells, replicates, seeds, draws, and Monte Carlo size as authorized; the overrun is wall/core-hours only.
- The addendum's 1e5-draw floor for M_eff, or fewer workers for the t7 cells, would have roughly halved/improved this; neither was substituted mid-run since the GO named 2e5 draws.

## Outputs

`mr_vs_guohe_t35_beta2_00..05.rds`, `mr_vs_guohe_t6_k02/06/10/12.rds`, `mr_vs_guohe_t7_beta2_00..05.rds` (16 bundles, 7.1 MB total, tracked). Each carries: per-replicate results (seeds; selection and agreement flags; naive est/SE/lower/cover; G&H debiased/lower/cover for all four r; MR est/se_ij/se_wald/lower_1s/two-sided CI/cover/bias_sel/bias_fix/ij_source/selection_rate/mean_r; p̂_Ĥ and top-3 with labels; Σ̂_ĤĤ, A6_mass, A6_mass_std, m̂(0) with MC SE, M_eff, tie_resid_implied; timings), the gate tally, all seeds and settings, elapsed, and full `sessionInfo`.
