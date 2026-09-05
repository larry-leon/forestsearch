# REPORT — MR (field) vs Guo & He: Stage 1 (implementation, identities, smoke, projection)

Date: 2026-09-05. Task: `dev/tasks/TASK_mr_field_vs_guohe_2026-09-05.md` (Stage 0 report `a5320ef1`).
Provenance cited per instruction: `dev/tasks/POC_mr_interval_alternatives_2026-09-05.md` (the first proof of concept — plug-in field, floor, tie-field — whose "preferred" recommendation is this `ci_method = "field"`), and `dev/tasks/poc2_results_2026-09-05.csv` (PoC-2, the authoritative reference for the 1b identity values and the shrunk-field / second-order design implemented here).
Decisions: E1–E5 defaults. E6: **resolved by Larry's unattended pre-authorization** — compute pre-approved if all identities pass, byte-identity and suite hold, smoke holds, and the loaded projection is ≤ 12 h wall at ≤ 120 workers.

## GATE 1: GREEN — pre-authorization conditions met; Stage 2 started immediately

Projection (loaded-machine, the 1d requirement): **47 core-h → ~23 min wall at 120 workers** — well inside the 12 h ceiling.

## 1a — Implementation (add-only; acknowledged method change)

`R/fs_mr_inference.R`:
- `ci_method = c("ij", "wald", "field")` — `match.arg` default unchanged (`"ij"` first); new trailing arguments `field_R_out = 1000L, field_R_in = 500L` (E2 defaults) after `return_reselection`.
- The field block sits after the complement block, before the return assembly, and runs only under `ci_method == "field"`, drawing under the derived seed `seed + 900000L` **after** the main multiplier stream is fully consumed — the `"ij"`/`"wald"` paths are untouched, RNG stream included. Implementation follows the spec exactly: shrunk field `w = β̂, w[Ĥ] = β̃` (E1); `ζ = crossprod(B, ξ)`, `ξ ~ N(0, I_n)` (no Cholesky, no explicit Σ̂); shared inner draws; `Λ*_r = ζ*_{r,G_r} − m̂(v_r)` with no-winner draws skipped (the bias_sel convention); the `S(·)` map is the gate's own — a vectorized `which.max`/`max.col(ties="first")` fast path when `reselection = "maxeff"` with no admission floor (equivalence machine-checked, F4), `.fs_mr_select` per draw otherwise.
- Returns `out$field` (the D6 append pattern): `lambda_mean`, `lambda_sd`/`se_field` and quantiles `q05/q25/q50/q75/q95/q025/q975` on the working scale; `n_out_used`, `n_in_used_mean`, `R_out`, `R_in`, `seed_offset`, `timing_seconds`; and on the effect scale `est2 = β̃ − λ̄` (E3), primary `lower_1s = β̃ − q95` (E4), two-sided quantile `lower_2s`/`upper_2s`, supplementary SE-type `lower_se`/`upper_se`.
- One behavior-preserving line: the debiased-CI switch is now `se <- if (ci_method == "wald") se_wald else se_ij$se`, so under `"field"` the `debiased` element is computed exactly as under `"ij"` — one call yields both MR rows. `"ij"` and `"wald"` semantics unchanged.
- Rd regenerated (only `man/fs_mr_inference.Rd` changed under `devtools::document()`).

**Test suite on the edited source: 0 fail / 4941 pass / 3 skip / 32 warnings — totals identical to the pre-change run of 2026-09-04.**

Harness: `mv_mr()` in `mr_vs_guohe_sim.R` gains a `ci_method` pass-through (default `"ij"`; field args forwarded only under `"field"`, keeping the default call argument-identical to the pre-change capture calls). Checks harness: `mr_field_stage1_checks.R` (committed; fail-loud, exit 0 on this run).

## 1b — Identities, measured values (all PASS)

| ID | Identity | Tolerance | Measured |
|---|---|---|---|
| F0a/b/c | default-path byte-identity, three fixed-seed cases (§5.1 m=1, §5.2 m=1, K=1) vs objects captured with the **pre-change installed package** | exact `identical()` on every component except `timing_seconds` | PASS ×3 |
| F1a (K=1, R_out=20,000) | `lambda_mean` ≈ 0 | 4 sd = 0.0325 | +0.00635 |
| F1b | `lambda_sd`/σ̂_D | [0.95, 1.05] | 1.0018 |
| F1c | est2 = β̃ | 4 sd (log) = 0.0325 | −0.00635 |
| F1d | q95 / (1.645 σ̂_D) | ±5% | 1.0220 |
| F2 | cov(ζ) = B_effᵀB_eff, §5.2 family K=154, R=20,000 | max correlation-scale dev < 0.06 (MC sd 0.010/entry) | 0.0256 |
| F3 | seed reproducibility of `field` | exact (timing excluded) | PASS |
| F4 | fast path ≡ `.fs_mr_select(maxeff, unrestricted)`; `max.col("first")` ≡ column-wise `which.max` | exact, 1,000 vectors + 500-column matrix | PASS |
| F5a | tie sign, K=10 disjoint × 200 reps: mean `lambda_mean` > 0 | t > 3 | 0.0487 (t = 35.9) |
| F5b | est2 closer to θ than β̃ | strict | β̃ bias 0.1037 → est2 0.0551 (MCSE 0.0117) |

**Byte-identity exclusion, stated plainly:** the first F0 run compared the full objects and failed; the diagnosis (single component, `timing_seconds`, 0.086 vs 0.085 s) showed the only difference was the wall-clock measurement of the selection loop — a quantity that differs between *any* two runs, pre-change included, so the literal test is unsatisfiable. F0 now asserts the name vector matches and **every other component** is `identical()`; the exclusion is confined to that one field and is recorded here and in the check's comment.

**F1 gating configuration:** gated at R_out = 20,000 so the q95 Monte Carlo error (~1.5%) sits inside the 5% tolerance; at the production configuration (1000/500) the same quantities are reported ungated: `lambda_mean` −0.0017, ratio 1.026, q95 ratio 1.047.

**F5 reference (PoC-2, Cox here vs normal-means there):** K=10 retained 0.1037 → 0.0551 against the csv's 0.0634 → 0.0249; K=2 reference row 0.0292 → 0.0102 against 0.0177 → 0.0041. Direction and rough magnitude of the second-order removal agree; the Cox tie bias itself runs higher than the normal-means model's, as expected for a different DGM — Stage 3 measures this on the real paired grid.

## 1c — Smoke (5 paired replicates × t35_beta2_00 / t6_k12 / t7_beta2_00)

All PASS: naive `identical()` to the stored replication bundles 15/15; **MR (current) est / se_ij / bias_sel / bias_fix `identical()` to the 2026-09-04 `mr_vs_guohe_*` bundles 15/15** on the same seeds (the one-call design working as intended); field output finite 15/15 (mean field pass 0.09 / 0.24 / 0.36 s per replicate by cell type).

## 1d — Projection, measured loaded (the 2026-09-04 lesson applied)

Probes of the actual Stage 2 per-replicate work (regeneration + naive recompute + one `ci_method = "field"` gate call + pairing asserts; **no M_eff pass**, E5) run under real 120-way `mclapply` load — 360/240/360 replicates of t35/t6_k12/t7 respectively, every probe replicate also passing the smoke asserts:

- Loaded core-s per replicate: t35 3.8, t6_k12 6.4, t7 6.7.
- **Full grid 16 × 2000: 47 core-h → ~23 min wall at 120 workers.** (Pre-authorization ceiling: 12 h. Met with 30× margin.)

## Disposition

Per the pre-authorization: Gate 1 conditions all hold → Stage 2 launched immediately at these settings (16 cells × 2000, same seeds, E1–E5 defaults, 120 workers, 14 h hard timeout), driver `mr_field_vs_guohe_run.R`.
