# REPORT — the null path: fast inversion, a null branch for the family, and the null-cell comparison

**Task:** `dev/tasks/cc_task_oc_wrapper_null_2026-08-29.md` (commit `6574beb5`)
**Repo / branch:** `~/Documents/GitHub/forestsearch` (pop-os), `feature/glm-extension`
**Follows:** `REPORT_oc_wrapper_grid_2026-08-29.md` and the two 08-28 reports (read first).
**`R/` category:** `R/fs_oc_grid.R` — the inversion and one-block sweep replaced by the exact order-statistic reduction (results unchanged; guards below); `R/fs_oc_family.R` — a null branch, additive. **`R/fs_dgm_scale.R` not edited**; §2 is read-only and its diff is a proposal.

---

## 0. Provenance (§1 — GATE passed)

```
pop-os · /home/larryleon/Documents/GitHub/forestsearch · feature/glm-extension
HEAD at start 399ccabe (descends from itself) · git status --porcelain: empty
packageVersion forestsearch 0.2.5
```

`~/Downloads/cc_task_oc_wrapper_null*`: one match, `cc_task_oc_wrapper_null_20260829.md` (hyphens stripped) → `dev/tasks/cc_task_oc_wrapper_null_2026-08-29.md`, committed **`6574beb5`**.

**Parity baseline** `devtools::test()`: `[ FAIL 0 | WARN 31 | SKIP 3 | PASS 4681 ]`.

---

## 1. §2 — read-only: why `fs_dgm_scale()` cannot serve the null as called, and what it would take (GATE: **opened** — the whole-population effective variance is computable without editing the file)

**(1) What fails.** Under a null DGM `df_super$flag_harm` is all zero (`R/generate_glm_dgm.R` L316–319: `if (model == "null") { data$flag_harm <- 0L }`). `fs_dgm_scale()` with the default `regions = NULL` builds `list(in_q, !in_q, rep(TRUE, n_super))` (`.fs_scale_regions()`, L317–329) and then, per region (L153–157):

```r
  rows <- lapply(names(regions), function(nm) {
    g <- regions[[nm]]
    if (!any(g)) {
      stop("region '", nm, "' is empty.", call. = FALSE)
    }
```

so it **errors**, cleanly, at the `Q` row: `Error: region 'Q' is empty.` — confirmed on the rebuilt null DGM. Nothing malformed is returned and nothing silently wrong is produced; the only step that depends on a non-empty Q is the default region construction. Everything before it (`.fs_scale_po()` L285, `.fs_scale_condvar()` L298, `sigma`) is Q-free.

**(2) What survives.** Every whole-population (`S`) quantity is computable, because the public `regions` argument accepts any named list of logical vectors (L317–350) and the per-region arithmetic (L158–178) never references Q. Run on the rebuilt null DGM with `regions = list(S = rep(TRUE, 5000))`:

| field | null value (S) | computable? |
|---|---:|---|
| `n_g`, `P_g` | 5000, 1 | yes (`sum(g)`, `mean(g)`) |
| `sigma` | 127.5001 (`sigma² = 16256.28`) | yes (`model_params$sigma`, region-free) |
| `v_cond0` | 16256.282 | yes (`sigma²`, constant for continuous) |
| `V_mu0` | 728.741 | yes (finite-population variance of `mu0` on S; equals the alternative's S row exactly — same frame) |
| `V_tau` | 2.2e-29 (= 0) | yes — under the null `tau` is constant (`mu1 − mu0 = −26.255236` on every row) |
| `bracket` (`V_eff/4`) | 16985.023 | yes |
| `V_eff` | **67940.092** | yes |
| `m_tau` | −26.255236 | yes; equals `effect_Qc = effect_ITT` of the payload's `truth` |

**(3) `V_eff` across regions.** Tracked n = 500 alternative payload, `$scale$regions`: `V_eff` = 67453.760 (Q), 68034.437 (Qc), 67891.079 (S) — spread **max/min − 1 = 0.861%**. It is **not** region-invariant by construction: `V_eff = V_arm1/p + V_arm0/(1−p)` with `V_arm = σ² + Var_g(mu)` (L158–160), and `Var_g(mu0)` is the finite-population variance of the control mean *within the region* (607.2 on Q, 752.3 on Qc, 728.7 on S against a common `σ² = 16256.3`), so regions differ by their spread of baseline means; the spread is small here only because `σ²` dominates `Var_g(mu0)` by 20×. The alternative wrapper reads `V_eff[Q]` (the document's anchor); under the null there is no Q, and the S row is the only region-agnostic choice — its `V_eff` (67940.1) sits inside the alternative's 0.86% band, so using it for `se_g` is sound at this fixture and changes `se_g` by ≤0.4% relative to the Q-row convention. For a DGM with large baseline-mean heterogeneity the S row would be the correct null scale in any case, since under the null every candidate is a mixture of S alone.

**(4) Proposal, not an edit — a null branch in `fs_dgm_scale()`.** Not needed for this task (the `regions` route is used), but if Larry wants `fs_dgm_scale(dgm)` to work as called on a null DGM:

```diff
--- a/R/fs_dgm_scale.R
+++ b/R/fs_dgm_scale.R
@@ .fs_scale_regions <- function(regions, df, harm_col, n_super, labels) {
   if (is.null(regions)) {
     if (!harm_col %in% names(df)) {
       stop("df_super has no column '", harm_col, "'.", call. = FALSE)
     }
     if (length(labels) != 3L) {
       stop("'labels' must have length 3.", call. = FALSE)
     }
     in_q <- df[[harm_col]] == 1
+    if (!any(in_q)) {
+      # Null DGM (Q empty): no Q/Qc partition exists.  Return the whole
+      # population only, under the third default label ("S").
+      out <- list(rep(TRUE, n_super))
+      names(out) <- labels[3L]
+      return(out)
+    }
     out <- list(in_q, !in_q, rep(TRUE, n_super))
     names(out) <- labels
     return(out)
   }
```

Effect: `fs_dgm_scale(null_dgm)` returns a one-row `regions` table (`S`) instead of erroring; `sigma`, `outcome_type`, `effect_measure`, `rand_ratio`, `p_treat`, `n_super` unchanged. **Callers that could be affected:** anything indexing `match("Q", scale$regions$region)` unguarded — the prediction document's `anchor` chunk (`stopifnot(!is.na(.iQ))`, so it would stop rather than misread), `fs_scale_se(scale, n, region = "Q")` (would return NA/error for a missing region name), and `fs_oc_family_enumerate()`'s alternative branch (guards `is.na(iQ)` and stops). The drivers guard the null cell before calling (`dgm_scale <- if (isTRUE(null_cell)) NULL else fs_dgm_scale(dgm)`), so they are unaffected either way. `print.fs_dgm_scale()` prints whatever rows exist. Nothing in `R/` calls `fs_dgm_scale()` on a null DGM today, so the change is additive; the two `stop()`s above it remain. **Not applied.**

---

## 2. §3 — `R/fs_oc_grid.R`: the exact order-statistic reduction (GATE passed)

At fixed `(n, gate, c2)`: `eligible_g = (Bhat_g − c2 >= z_p·se_g)` (resample) or `(W1_g >= c2) & (W2_g >= c2)` (split); per draw `w = argmax_{eligible} Bhat_g` (NA if none), `T = Bhat_w` (−Inf if none). Declaration at `c1` is `T >= c1` (conjoined with "any eligible" so that `c1 = −Inf` does not admit `T = −Inf`), and `w` is the maxeffCons winner whenever anything declares.

Implemented as `.fs_oc_reduce(dr, fam, c2, pc, gate)` (one `max.col` per `c2`; stores `T`, `w`, `any_el`, the winner's noise `Bhat_w − beta_g[w]`, and the eligible matrix for the per-member `P1`) and `.fs_oc_from_reduction(red, fam, n, c1, Rdraw)`. **Comparison operators unchanged** (`>=` for declaration and per-member pass, `<` for `mass_below`). The one-block sweep now loops `c2` → reduction → every `c1`; `fs_oc_invert(solve_for = "c1")` returns the `k`-th largest `T`, `k = ceiling(target·draws)`, with the ceiling `mean(T > −Inf)` and the binding threshold named above it; `target` may be a vector (one draw set serves all). `solve_for = "c2"` keeps the bisection (the eligible set moves with `c2`).

**On bit-identity.** The per-`c1` quantities are formed from the same elements in the same order as `.fs_oc_functionals()` — `det` is the same logical vector, `winner` the same integer vector, `P1 = colMeans(eligible & Bhat >= c1)`, and the winner-noise mean is `mean()` over the same values in row order — so a grid point remains `identical()` to `fs_oc_predict()` (test 2 of `test-fs-oc-grid.R`, both gates). The task's `cumsum`-over-sorted-prefix form for the conditional means was **not** used for those quantities precisely because a different summation order would break that identity; the sorted order is used where it is exact — the inversion's order statistic.

**Guard — the twelve stored inversions** (`invert_reduction_guard_2026-08-29.R`, same seed 20260825 and 2×10⁵ draws, one draw set per `(n, gate)`, vector targets):

| n | gate | target | c1 (stored, bisection) | c1 (reduction) | Δ | local step T(k−1)−T(k) / T(k)−T(k+1), max | old tol | within one step |
|---|---|---:|---:|---:|---:|---:|---:|---|
| 500 | resample | 0.80 | 91.9010 | 91.9016 | 6.4e−04 | 1.7e−03 | 1e−3 | TRUE |
| 500 | resample | 0.90 | 84.7186 | 84.7187 | 1.6e−04 | 1.1e−03 | 1e−3 | TRUE |
| 500 | resample | 0.95 | 78.9870 | 78.9876 | 5.9e−04 | 9.6e−04 | 1e−3 | TRUE |
| 500 | split | 0.80 | 91.8461 | 91.8470 | 9.4e−04 | 1.2e−03 | 1e−3 | TRUE |
| 500 | split | 0.90 | 84.5709 | 84.5714 | 4.9e−04 | 1.2e−03 | 1e−3 | TRUE |
| 500 | split | 0.95 | 78.8253 | 78.8260 | 6.9e−04 | 1.6e−03 | 1e−3 | TRUE |
| 700 | resample | 0.80 | 93.1123 | 93.1129 | 6.2e−04 | 4.8e−05 | 1e−3 | TRUE |
| 700 | resample | 0.90 | 86.5179 | 86.5179 | 7.6e−06 | 2.5e−04 | 1e−3 | TRUE |
| 700 | resample | 0.95 | 81.2837 | 81.2841 | 3.7e−04 | 1.9e−03 | 1e−3 | TRUE |
| 700 | split | 0.80 | 93.1398 | 93.1403 | 5.6e−04 | 1.9e−04 | 1e−3 | TRUE |
| 700 | split | 0.90 | 86.5218 | 86.5221 | 3.7e−04 | 4.4e−04 | 1e−3 | TRUE |
| 700 | split | 0.95 | 81.3058 | 81.3063 | 5.4e−04 | 7.6e−04 | 1e−3 | TRUE |

Comparison used: `|c1_new − c1_old| <= max(local order-statistic step around k, old bisection tol 1e−3)` on the same draw set (same seed); **largest discrepancy 9.41×10⁻⁴**, all twelve within. The stored bisection values sit just *below* the exact order statistic (the bisection returned the lower bracket end within `tol`), which is why the stored achieved rates exceed the target by up to 1×10⁻⁵ (one or two extra draws) while the reduction's achieved rate is exactly `k/draws`. Nothing exceeds one step. `REDUCTION GUARD: PASS`.

**Wall-clock.** Twelve inversions by the reduction: **2148 s** of `fs_oc_invert()` time (four draw sets of 322 / 641 / 409 / 776 s — the draw *is* the cost; the reduction and the seven-target read-off are seconds) against **8443 s (2.35 h)** stored for the slow path — 3.9×, and 12× per inversion once a draw set exists. (The guard script itself ran 4305 s because it draws each set a second time to measure the local step.) The one-block sweep benefits the same way: 21 `c1` values in 70 s instead of 376 s on M = 1601.

Existing GATE tests (grid/predict identity, monotonicity, round-trip, ceiling): `test-fs-oc-grid.R` `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 78 ]` after the reduction.

---

## 3. §4 — `R/fs_oc_family.R`: the null branch (GATE passed: alternative unchanged)

**Detection.** From source: `generate_glm_dgm(model = "null")` zeroes `flag_harm` (L316–319) and stores `model = model` in the object (L559); the null driver asserts `sum(inQ) == 0L`, `is.na(truth$effect_Q)`, `isTRUE(all.equal(truth$beta_inter, 0))`; the null payload's `truth` has `effect_Q = NA`, `prevalence_Q = 0`, `beta_inter = 0`. The wrapper tests **`!any(df_super$flag_harm == 1)`** — the structural fact the driver asserts, available on any `glm_dgm` (hand-built test fixtures included) — and, when `dgm[["model"]]` is present, requires agreement (either mismatch → `stop("… the object is inconsistent")`). An alternative DGM always has a non-empty calibrated Q, so the flag test cannot misfire on one; `effect_Q`/`beta_inter` were not chosen because they live in `hazard_ratios`/`model_params` sub-lists that a minimal object need not carry. (Found and fixed while testing: `dgm$model` partial-matches `model_params` when `model` is absent — indexing is `dgm[["model"]]`.)

**Fields under the null** (as the task's table): `beta_g = |m_tau[S]|` for every g (= |effect_Qc| = |effect_ITT|, oriented as the alternative orients); `se_g = sqrt(V_eff[S]/(1000·1))·sqrt(1000/n)·sqrt(1/Pg)` — the alternative's `seQ1000·sqrt(1000/n)·sqrt(piQ/Pg)` with `P(Q) → P(S) = 1` and the S row from `fs_dgm_scale(dgm, regions = list(S = …))`; `PQg = 0`; `sens_g = NA_real_`; `spec_g = 1 − Pg`; `PQ = 0` (from which `Enpv = 1` follows in `fs_oc_predict()`'s formula). Enumeration, floors, `ovl`, `Rho`, `Sg`, root untouched. New element `null` (logical); `print()` flags the branch.

**Propagation checked** (not assumed), synthetic null DGM, both gates: `Esens = NA`; `det_rate`, `EnH`, `Espec`, `Eppv` (= 0), `Enpv` (= 1), `EbetaH` (= the common effect), `Enaive_bias`, `mass_below` (= 1 at `c1 = 30` since 26 < 30) all finite; `fs_oc_grid()` rows carry `NA` only in `Esens`; `fs_oc_invert()` runs. Tests: `tests/testthat/test-fs-oc-null.R` `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 56 ]` — field definitions, `se_g` against the S row, same M/labels/`ovl` as the alternative on the same covariates, detection without a `model` field, both inconsistency errors, NA propagation, monotonicity in `c1`, vector-target inversion, and the alternative path against the frozen reference.

**Alternative unchanged** — after the branch: `prerefactor_reference_2026-08-29.R check` → `REFACTOR GUARD: PASS (identical to the 0.2.4 reference)`; `fidelity_fs_oc_predict_2026-08-28.R` → `FIDELITY GATE: PASS (bit-identical)`.

---

## 4. §5 — the null-cell comparison (`oc_wrapper_null_2026-08-29.R` → sibling `oc_wrapper_null_2026-08-29.rds`)

Null DGM rebuilt as the driver builds it under `null_cell = TRUE` (`generate_glm_dgm(model = "null")`, same frame/factors/subgroup spec, seed 8316951), gated against the null payload's committed `truth`: `Q_empty`, `effect_Q_NA`, `effect_Qc` (diff 0), `beta_inter` (0), `effect_ITT` (diff 0), `model_null` — all TRUE. The payload's `NA`-named third element is neutralised by renaming to `""` and indexing by name. Settings asserted against `$meta` (`resample`, 30 / 10, n = 500).

**Family under the null floor:** M = **1601**, the same 1601 candidates as the alternative at n = 500 (same super-population frame under the same seed; only `flag_harm` is zeroed) — stage counts 72 columns → 2628 → empty 119 / minp 0 / rmin 99 / size 631 → kept 1779 → duplicates 178. Common effect `26.255236`; S-row `V_eff = 67940.09`; `PQg = 0`, `sens_g = NA`, `spec_g = 1 − Pg` on every member.

**The null cell at n = 500, c1 = 30, c2 = 10** (2×10⁵ draws, seed 20260825; measured: null payload, 998 declaring of 1000):

| quantity | analytic `"resample"` | analytic `"split"` | simulation (null payload) |
|---|---:|---:|---:|
| false declaration rate | 0.9970 (0.0001) | 1.0000 (0.0000) | 0.998 (0.0014) |
| per-candidate rate range, min .. max `P1` | 0.1211 .. 0.3784 | 0.3460 .. 0.3778 | — |
| implied `L_eff = log(1−fam)/log(1−max P1)` | 12.2 | 25.7 (at fam = 1 − 5×10⁻⁶ resolution) | — |
| E\|Ĥ\| given declaration | 71.0 | 70.9 | 72.05 (0.44) |
| E[sens] | NA (0/0) | NA (0/0) | NA (payload: 0/0) |
| E[spec] | 0.858 | 0.858 | 0.856 (0.002) |
| E[PPV] / E[NPV] | 0 / 1 | 0 / 1 | 0 / 1 |
| E[β(Ĥ)] oriented | 26.255 | 26.255 | 26.255 (0) |
| naive bias vs β(Ĥ) | 76.22 | 76.05 | 74.85 (0.58) |
| selection mass on rules with mean < c1 | 1.000 | 1.000 | — |
| M | 1601 | 1601 | — |

Wall-clock: 337 s (resample) / 648 s (split) for the single-point evaluations.

**Type-I error curve and inversions** (`c1 = 20 … 120` by 5 at `c2 = 10`; one draw set per gate — 315 / 626 s — then 70 / 69 s for the 21-point sweep; inversions: seven targets from the same draw set in 334 / 628 s):

| c1 | resample | split | split − resample | | target | c1 resample | c1 split |
|---:|---:|---:|---:|---|---:|---:|---:|
| 30 | 0.9970 | 1.0000 | +0.0030 | | 0.05 | 133.11 | 133.24 |
| 50 | 0.9963 | 0.9993 | +0.0030 | | 0.10 | 125.76 | 125.71 |
| 60 | 0.9933 | 0.9946 | +0.0012 | | 0.20 | 117.11 | 117.01 |
| 70 | 0.9723 | 0.9715 | −0.0009 | | 0.50 | 101.57 | 101.54 |
| 80 | 0.8995 | 0.8980 | −0.0015 | | 0.80 | 87.15 | 87.07 |
| 100 | 0.5352 | 0.5349 | −0.0003 | | 0.90 | 79.96 | 79.81 |
| 120 | 0.1608 | 0.1601 | −0.0007 | | 0.95 | 74.21 | 74.04 |

All fourteen inversions attainable (ceilings 0.9970 resample, 1.0000 split); achieved rates equal their targets to the draw resolution; MC SEs 0.0005–0.0011.

**Do the two gates differ here?** Yes, but only at the ceiling, and only by 0.003: the production `"resample"` screen leaves 0.3% of null draws with no eligible candidate (`Bhat_g − c2 < z_p·se_g` for all 1601 members), so its false-declaration rate is 0.9970 for every `c1 ≤ 50`, while the single-split gate declares in every draw (1.0000). Once `c1` binds (`c1 ≥ 65`) the two curves agree within 0.0017 — one to two MC SEs — and the inverted `c1` at any type-I target differ by ≤ 0.17. So at the null the gates are distinguishable at the ceiling and not along the curve: the consistency screen's form matters only through *whether anything is eligible*, and under this null almost everything is, because the common effect (26.3) sits 0.7–2.4 candidate-SEs above `c2 = 10` at these prevalences. The per-candidate range makes the same point: `"resample"` has members at `P1 = 0.12` (small, high-SE candidates the `z_p·se_g` term screens out), `"split"` none below 0.35.

**Against the measured null cell:** false declaration 0.9970 / 1.0000 vs 0.998 (the measured 2 non-declarations of 1000 sit between the two gates, 0.7 MC SE from `"resample"`); E|Ĥ| 71.0 vs 72.05 ± 0.44 (2.4 SE); specificity 0.858 vs 0.856 ± 0.002; E[β(Ĥ)] exactly the common effect on both sides; naive bias 76.2 vs 74.9 ± 0.6 (2.3 SE, the same sign and size as at the alternative). The document's stylised M = 16 `worked-null` could not reach the measured 998/1000 because sixteen near-threshold candidates cannot produce a family-level false declaration of 0.998; the enumerated 1601 do, with `L_eff` 12–26 rather than the document's implied effective size — the discrepancy was the family, as at the alternative.

---

## 5. §6 — the document

`oc_wrapper_verification.qmd` gains "The null cell": the field table, the comparison table (analytic resample · analytic split · simulation, plus the `P1` range, `L_eff` and `M`), the type-I error curve in the base-graphics idiom, the inverted-`c1` table, and the two sentences on the M = 1601 result versus the document's M = 16. All numbers inline from the sibling `.rds`.

Rendered with RStudio's bundled Quarto 1.9.38: **exit 0 in 8 s**; `oc_wrapper_verification.html` (2.4 MB) committed beside it.

---

## 6. §7 — close-out

`devtools::document()`: `NAMESPACE` +`S3method(print, fs_oc_invert_list)` only; `man/fs_oc_family_enumerate.Rd`, `man/fs_oc_invert.Rd` rewritten.

`devtools::test()`:

```
[ FAIL 0 | WARN 31 | SKIP 3 | PASS 4737 ]
```

Parity against the baseline `[ FAIL 0 | WARN 31 | SKIP 3 | PASS 4681 ]`: WARN 31 = 31 (parity), SKIP 3 = 3, FAIL 0 = 0, PASS 4681 + 56 = 4737 — the new null tests and nothing else.

`R CMD check` (`devtools::check(document = FALSE, args = "--no-manual")`, `RSTUDIO_PANDOC` set):

```
── R CMD check results ───────────────────────────────── forestsearch 0.2.6 ────
Duration: 11m 24.8s

0 errors ✔ | 0 warnings ✔ | 0 notes ✔
```

Version `0.2.5 → 0.2.6`; `NEWS.md` names the null branch and the reduction and states alternative-cell results are unchanged. `devtools::install()`: `packageVersion("forestsearch")` → `0.2.6`.

Commits (explicit paths; no push): `6574beb5` task doc; __COMMITS__.

---

## 7. Verdict (ten lines)

1. §2 gate opened without touching `R/fs_dgm_scale.R`: the default regions stop at "region 'Q' is empty", but the public `regions` argument yields the whole-population row — `V_eff[S] = 67940.1`, inside the alternative's 0.86% band; the file-level null branch is written as a diff for Larry, not applied.
2. `V_eff` is not region-invariant by construction (it carries the within-region variance of baseline means); it is numerically close here because `σ²` dominates 20×.
3. The order-statistic reduction reproduces all twelve stored inversions to ≤9.4×10⁻⁴ in `c1` (within the old tolerance and the local step) and cuts twelve inversions from 2.35 h to 36 min of `fs_oc_invert()` time — the draw is now the only cost; grid points remain `identical()` to `fs_oc_predict()`.
4. The null branch detects Q-empty structurally (the driver's own assertion), cross-checked against `dgm[["model"]]` (exact indexing — `$` partial-matched `model_params`, found and fixed); fields as specified; `NA` sensitivity reaches only `E[sens]`.
5. Alternative path unchanged: 0.2.4 reference `identical()` on both gates and the document's chunk bit-identical, re-run after every `R/` edit.
6. Null cell, M = 1601: false declaration 0.9970 (production) / 1.0000 (split) against the measured 0.998; E|Ĥ|, specificity, E[β(Ĥ)] and the naive bias line up as at the alternative.
7. The two gates differ at the null — by 0.003, at the ceiling only, from the 0.3% of draws with no eligible member under the `z_p·se_g` screen; along the type-I curve (`c1 ≥ 65`) they agree within 1–2 MC SEs and their inverted `c1` within 0.17.
8. The document's stylised M = 16 null figure was a family artefact: the enumerated family reaches the measured 998/1000 with `L_eff` 12–26.
9. Tests +56 (parity 31 = 31), check 0/0/0 at 0.2.6, document extended and rendered in 8 s, all committed by explicit path, no push.
10. Open for Larry: whether to apply the `fs_dgm_scale()` null diff; and that a type-I error of 0.05 at this fixture needs `c1 ≈ 133` — the driver's 30 is a saturating screen under this null, which is the workstream's last-gap finding stated as a number.
