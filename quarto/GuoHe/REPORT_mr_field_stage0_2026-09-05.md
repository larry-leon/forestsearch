# REPORT — MR (field) vs Guo & He: Stage 0 (discovery)

Date: 2026-09-05. Task: `dev/tasks/TASK_mr_field_vs_guohe_2026-09-05.md` (committed `1f771bf8`).
Repo state at discovery: branch `feature/glm-extension`, HEAD `1f771bf8`; `R/fs_mr_inference.R` last touched by `736088e3` (the D6 change), worktree clean for that file. Predecessor arc complete through `b0b5da5e`; the 16 `mr_vs_guohe_*.rds` bundles are tracked and untouched.
Decisions: E1–E5 taken at defaults (the kickoff left the placeholder); E6 deferred to Gate 1 as specified.

## GATE 0: GREEN

All of 0a–0c resolved and quoted below. **No change to existing arithmetic is needed to reach the field quantities**: every input the method consumes (β̂ over the kept family, B_eff, Ĥ, β̃, the re-selection routine, the seed) is in scope inside `fs_mr_inference()` at the point where an add-only `ci_method = "field"` block would sit, and the re-selection helper is a pure function callable on an arbitrary K-vector. The method-change acknowledgement is in the task and is not re-litigated here.

One provenance flag (not a stop condition): of the two proof-of-concept artifacts the task cites, only `poc2_results_2026-09-05.csv` exists on disk (delivered to `~/Downloads`, copied to `dev/tasks/` beside the task doc for durability); **`POC_mr_interval_alternatives_2026-09-05.md` is nowhere under `~/Downloads` or the repo** (case-insensitive find). The csv contains the 1b reference values the task quotes — S5 "K=10 disjoint, tie": `bias_bt` 0.0634 → `bias_bt2_shr` 0.0249 ("0.063 → 0.025"); S4 "K=2 disjoint, tie": 0.0177 → 0.0041 ("0.018 → 0.004") — so Stage 1's sign-check reference is verifiable without the md. If Larry wants the md in the record, deliver it and it will be added; nothing in Stages 0–3 requires it.

## 0a — The gate's internals (all in `R/fs_mr_inference.R` at `736088e3`)

**β̂ and B_eff after the fits.** `asm <- .fs_mr_assemble(df, candidates, spec)` (:394) builds the kept family; the working objects unpack at :409:

```r
B <- asm$B; bh <- asm$beta_hat; sdv <- asm$sigma_D; sz <- asm$sizes
```

`B` is the n × K influence matrix (candidate g's dfbeta in its member rows, zero elsewhere — `.fs_mr_assemble`, :88-105), `bh` the raw-coefficient vector β̂, `sdv` the per-candidate robust SEs, `sel` the winner index (`sel <- match(sel_lab, asm$names)`, :396). The two-term estimate is `beta_deb <- beta_naive - selection_bias - fb` (:498). Everything the field spec's step 1–2 needs (`w = β̂; w[Ĥ] = β̃`; `ζ = B_effᵀ ξ`) is therefore live in scope after :498, before the return assembly.

**Re-selection routine, callable on an arbitrary K-vector.** Per draw the gate calls (:467-469):

```r
s <- .fs_mr_select(bs, .zcons(bs), sz, pass, reselection,
                   effect_neighborhood, selection_rule, log_scale)
if (!is.na(s)) { sel_bias[b] <- P[s, b]; winner[b] <- s }
```

`.fs_mr_select(beta, zcons, sizes, passers, rule, nbhd, selection_rule, log_scale)` (:132) is a pure function of its arguments — no captured state — so `S(v)` on any perturbed K-vector `v` is exactly `.fs_mr_select(v, .zcons(v), sz, .admit(v), ...)`. Under this task's configuration the rule is `maxeff = passers[which.max(beta[passers])]` (:170) with the unrestricted admission path `.admit <- function(bs) seq_along(bs)` (:441), so the inner loop vectorizes as a row-wise argmax of an R_in × K matrix, exactly as the task's cost note anticipates (an implementation shortcut valid only for maxeff + unrestricted; the general path stays `.fs_mr_select` per draw — Stage 1 will take the vectorized route and assert its equality to `.fs_mr_select` on a sample of draws).

**Gaussian multiplier path.** `.fs_mr_multipliers` (:110) already has `gaussian = matrix(stats::rnorm(n * draws), n, draws)` (:113) — mean-zero unit-variance i.i.d. N(0,1), which is precisely the ξ the field spec wants; `ζ = B_effᵀξ = crossprod(B, xi)` matches the gate's own `P <- crossprod(B, Xi)` pattern (:459). No Cholesky and no explicit Σ̂ anywhere, as specified.

**Seeding.** `if (!is.null(seed)) set.seed(seed)` at :380, then the stream is consumed by `Xi <- .fs_mr_multipliers(nrow(B), draws, multiplier)` (:458) and the draw loop. A field block inserted after the IJ computation (:510) draws its ξ from the stream *after* Xi under a derived seed (Stage 1 plan: `set.seed(seed + 900000L)` when `seed` is non-NULL, recorded in the output; offset disjoint from the arc's data/boot/MR offsets 0 / 5e5 / 7e5). Because the field block is behind `ci_method == "field"`, the default path's RNG consumption is untouched.

**Return-list assembly (the D6 model).** `out <- list(...)` at :597, conditional add-only append at :618-623 (`out$reselection <- ...`). The field element follows the identical pattern: `out$field <- list(...)` behind the new `ci_method` value. `ci_method = c("ij", "wald")` (:374) is resolved by `match.arg`, so appending `"field"` leaves the default `"ij"` (first element) and both existing values byte-identical. Stage 1 design note, recorded now: under `ci_method = "field"` the `debiased` element will be computed exactly as under `"ij"` (IJ SE), so one call yields both the current MR row and the field row — which is what lets Gate 2 assert MR-current est/se_ij `identical()` to the 2026-09-04 bundles from the new run's own output.

## 0b — Harness and driver; minimal change

Regeneration fidelity is not merely "confirmed on 2026-09-04" by assertion — it is enforced per replicate: the Stage 2 driver `quarto/GuoHe/mr_vs_guohe_run.R` recomputes the naive column inside every worker and Gate 2 recorded `naive_mismatch 0` on all 32,000 replicates across the 16 cells (`REPORT_mr_vs_guohe_gate2_2026-09-04.md`), with the same seed bases this task will reuse (`mv_gh51_base`/`mv_gh52_base` in `mr_vs_guohe_sim.R`, asserted against each bundle's stored `seed_base`).

Minimal change for the new run (Stage 1 work, listed here per 0b):
- `mr_vs_guohe_sim.R`: `mv_mr()` gains a pass-through `ci_method` argument (default `"ij"`, so the committed Stage 1 harness behavior is unchanged); no other adapter change — family construction, orientation, seeds identical.
- A new driver `mr_field_vs_guohe_run.R` (the existing driver stays untouched, per "do not modify the existing bundles" — new bundles `mr_field_vs_guohe_<id>.rds`): per replicate it (i) recomputes naive and asserts against the replication bundle as before, (ii) calls the gate once with `ci_method = "field"` and the *same MR seed* `base + m + 700000L`, (iii) asserts MR-current `est`/`se_ij` `identical()` to the 2026-09-04 `mr_vs_guohe_<id>.rds` row (the new pairing proof — same seed → same Xi → same β̃/IJ), (iv) records the `field` fields and the two field coverage indicators, (v) **skips the M_eff pass** (E5) and joins `p_hat_H`/`p_top*`/`Sigma_HH`/`A6_*`/`m0_*`/`M_eff`/`tie_resid_implied` from the old bundle by `(id, m)` — the old results carry `m` and `seed_data` per row, so the join is exact and verifiable on `seed_data`.

Identity subtlety recorded in advance: the `identical()` in (iii) requires the field block to consume no RNG before Xi and the IJ computation — guaranteed by placing it after :510 under the derived seed, per 0a. If exact identity nonetheless fails at Stage 1's smoke, the fallback comparison is equality to 0 ulps via `==`, and a failure there is a Gate 1 stop.

## 0c — K, n, orientation per lineage (unchanged from the 2026-09-04 arc)

| Lineage | Cells | K | n | Orientation handling |
|---|---|---|---|---|
| §5.1 Tables 3–5 | t35_beta2_00..05 | 2 (disjoint `S1`,`S2`) | 400 | treat-flip adapter: `treat_flip = 1 - treat` in `mv_gh51_regen()`, spec `mv_spec51` (`treat.name = "treat_flip"`); flip exactness measured at ~1e-15 (Stage 1 V6, 2026-09-04) |
| §5.1 Table 6 | t6_k02/06/10/12 | 2/6/10/12 (disjoint) | 200·k | same treat-flip adapter |
| §5.2 Table 7 | t7_beta2_00..05 | data-determined ≈151 (nested order-statistic cutpoints, random per replicate) | 400 | raw scale (`orient = +1` equivalent): argmax of the raw coefficient, no adapter |

MR settings carried over unchanged for the current-MR row: `reselection = "maxeff"`, unrestricted admission, B = 5000 centred Poisson, IJ; field adds R_out = 1000 / R_in = 500 Gaussian (E2 default) under the derived seed.

## Decisions E1–E5 against source (defaults, as the kickoff indicated)

- **E1** shrunk field `w[Ĥ] = β̃`, no second unshrunk field — spec step 1 as written; both quantities in scope (0a).
- **E2** R_out = 1000, R_in = 500, Gaussian — the `gaussian` multiplier path exists (:113); cost note's vectorized argmax valid on this configuration (0a).
- **E3** point estimate `est2 = β̃ − lambda_mean`; the current MR row keeps β̃ + IJ — enabled by the one-call design in 0a/0b.
- **E4** primary one-sided `β̃ − q95`; two-sided quantile supplementary — both from the same Λ* sample.
- **E5** M_eff off, joined by `(id, m)` from the 2026-09-04 bundles — join keys verified present.

— End of Stage 0. Stopped at Gate 0 per instruction. Next (on go): Stage 1 — the add-only `ci_method = "field"` implementation, byte-identity check on three fixed-seed cases, the 1b identities, the 5×3-cell smoke, and the loaded-machine (≥30 concurrent workers) projection the 2026-09-04 overrun mandates.
