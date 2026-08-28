# REPORT — `sigma_D` verified, the closed-form gate added, n = 700 run, both documented

**Task:** `dev/tasks/cc_task_oc_wrapper_gate_and_n700_2026-08-29.md` (commit `aa83b3ea`)
**Repo / branch:** `~/Documents/GitHub/forestsearch` (pop-os), `feature/glm-extension`
**Follows:** `REPORT_oc_wrapper_build_2026-08-28.md`, `REPORT_oc_wrapper_fixture_run_2026-08-28.md` (both read first).
**`R/` category:** one conditional edit, `R/fs_oc_predict.R` — the `"resample"` branch — taken because §2's gate opened. No other `R/` file.

---

## 0. Provenance (§1 — GATE passed)

```
pop-os · /home/larryleon/Documents/GitHub/forestsearch · feature/glm-extension
HEAD at start 8ecec589 (descends from itself) · git status --porcelain: empty
git log -6: 8ecec589 1334a6dd 438a2e81 9e538ec2 47769eb3 c016df02
packageVersion forestsearch 0.2.3
```

`~/Downloads/cc_task_oc_wrapper_gate_and_n700*`: exactly one match, `cc_task_oc_wrapper_gate_and_n700_20260829.md` (hyphens stripped); copied to `dev/tasks/cc_task_oc_wrapper_gate_and_n700_2026-08-29.md`, committed **`aa83b3ea`**.

**Parity baseline** `devtools::test()`: `[ FAIL 0 | WARN 31 | SKIP 3 | PASS 4577 ]`.

Note: the task's "Follow `claude/quarto_applications_conventions.md`" — that file does not exist in the repo or under `~` (searched; the only mention is the task doc itself). The new document follows the prediction document's own YAML and style.

---

## 1. §2 — the `sigma_D` identification (GATE: **opened**, ratio ≈ 1, no factor)

### 2.1 From source

**(1) Half-effect convention.** `R/consistency_resample.R` L11–12 (header): "the random-split half effects are `beta_hat +/- D, D = sum_i G_i * dfbeta[i, treat]`". The split is 50/50 by construction of the Rademacher multiplier (L13: "Rademacher = the 50/50 split"). Neither half is named `W1`/`W2` in source; the document's `W1`, `W2` are the two half effects `beta_hat + D` and `beta_hat − D` in either order — the rate is symmetric in them.

**(2) `D` and the multipliers.** L12–14: `D = sum_i G_i * dfbeta[i, treat]` "with `G_i` any mean-zero, unit-variance multiplier (Rademacher = the 50/50 split, normal = Lin 1997, centred Poisson = Dobler-Beyersmann-Pauly)". The MC branch (L451–456) realizes exactly that: `rademacher = matrix(2 * rbinom(draws * n, 1, 0.5) - 1, …)`, `D <- as.numeric(G %*% pieces$dfbeta)`. Every package call site uses `method = "closed"`, so `D` is never realized in production — only its SD is used.

**(3) `sigma_D`.** L86 (Cox) and L318 (GLM): `sigma_D <- sqrt(sum(dfbeta^2))` / `sqrt(sum(dfb^2))` — confirmed. `dfb` is `.dfbeta_glm(fit)[, treat]` (L310–316), the exact leave-one-out influence. **No scaling factor** is applied anywhere between computation and use: `out$sigma_D <- pieces$sigma_D` (L439) and the rate uses `delta / pieces$sigma_D` (L446) directly.

**(4) Comparison scale `c`.** L430–436: `c_cmp <- comparison_threshold` if supplied, else `log(thr_nat)` for ratio measures, else `thr_nat` (identity). On the GLM path `forestsearch()` passes `comparison_threshold = consistency_threshold` (`forestsearch_main.R` L3071–3079), with `consistency_threshold <- hr.consistency` after `hr.consistency <- consistency.threshold` (L1712; L1290s) — for MD that is `consistency.threshold = 10` on the effect scale, no log. **Orientation:** `.consistency_glm_pieces()` L255–261 negates `Y` for continuous outcomes when `adverse_outcome = FALSE` (the driver's setting), so `beta_hat` is the oriented (harm-positive) coefficient and `c = 10` is compared on the oriented scale; no flip is applied to `c`.

**(5) Pass condition.** L446: `rate_closed <- max(0, 2 * pnorm(delta / sigma_D) - 1)`, `delta = beta_hat − c_cmp`. Call site `R/subgroup_consistency_helpers.R` L1457: `p.consistency <- round(rr$rate_closed, pconsistency.digits)` (digits = 2), then L1507: `if (isTRUE(p.consistency < pconsistency.threshold)) { … return(NULL) }` — a candidate **passes iff `p.consistency >= pconsistency.threshold`**, i.e. `>=`, after rounding to two decimals. (The two-stage path L1843 uses the same `<` → drop.)

### 2.2 The threshold in general form

From (4)–(5), with `p = pconsistency.threshold` symbolic:

```
2 * Phi((beta_hat - c) / sigma_D) - 1  >=  p
  <=>  Phi((beta_hat - c) / sigma_D)  >=  (1 + p) / 2
  <=>  beta_hat  >=  c + qnorm((1 + p) / 2) * sigma_D
```

`qnorm((1 + 0.90) / 2) = 1.644854` — confirmed numerically in the diagnostic. The source's left-hand side is exactly this (the `max(0, ·)` clamp cannot change a `>= p` decision for `p > 0`); the only source detail not in the inversion is the **two-decimal rounding** of the rate before comparison, which the wrapper does not replicate (it would move the threshold by at most `sigma_D · Δz` for a 0.005 change in rate — ~0.01 SE — and is not a population quantity).

### 2.3 DIAGNOSTIC — one simulated trial (identity check, not a result)

`dev/glm-continuous-sims/sigma_d_diagnostic_2026-08-29.R` (+ `.log`): MD40 fixture rebuilt and gated; **one** trial at n = 500, `simulate_from_glm_dgm(dgm, n = 500, seed = 20260829)`; ten candidates from the M = 1601 enumeration at evenly spaced prevalence ranks (population `Pg` 0.120–0.920, sample sizes 61–453). `sigma_D` from `consistency_resample()` called exactly as the GLM call site calls it; `se_model` from `summary(lm(-y ~ treat))`; `se_hc0 = sqrt(sum(dfbeta[, treat]^2))` from base R's `lm.influence()$coefficients` (not the package helper).

| candidate | Pg (pop) | n_g | sigma_D | se_model | se_hc0 | se_g (pop) | sigma_D/se_model | sigma_D/se_hc0 | sigma_D/se_g |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| age > 39 & cd40 > 349 | 0.120 | 61 | 34.466 | 35.907 | 34.466 | 33.530 | 0.960 | 1 | 1.028 |
| age <= 31 & karnof <= 95 | 0.150 | 83 | 27.739 | 26.983 | 27.739 | 29.970 | 1.028 | 1 | 0.926 |
| age > 31 & age <= 35 | 0.185 | 88 | 26.797 | 26.759 | 26.797 | 27.033 | 1.001 | 1 | 0.991 |
| wtkg > 83 & gender | 0.221 | 105 | 24.888 | 24.605 | 24.888 | 24.685 | 1.012 | 1 | 1.008 |
| age > 39 | 0.267 | 145 | 20.966 | 20.730 | 20.966 | 22.495 | 1.011 | 1 | 0.932 |
| age <= 37 & cd40 > 339 | 0.324 | 152 | 19.263 | 19.165 | 19.263 | 20.399 | 1.005 | 1 | 0.944 |
| age <= 42 & wtkg > 75 | 0.385 | 185 | 17.426 | 17.329 | 17.426 | 18.729 | 1.006 | 1 | 0.930 |
| preanti > 188 | 0.453 | 209 | 18.048 | 17.923 | 18.048 | 17.250 | 1.007 | 1 | 1.046 |
| preanti <= 875 & race != 1 | 0.563 | 292 | 14.340 | 14.300 | 14.340 | 15.482 | 1.003 | 1 | 0.926 |
| age <= 47 | 0.920 | 453 | 11.622 | 11.604 | 11.622 | 12.112 | 1.002 | 1 | 0.960 |

- `sigma_D / se_hc0`: **1.0000000000 on every row** (identical to 1e-10) — `sigma_D` *is* the dfbeta sum of squares, computed by an independent route.
- `sigma_D / se_model`: 0.960–1.028, **Spearman with `Pg` = −0.07** — no systematic dependence on prevalence or size. The largest deviation (0.96) is the n_g = 61 candidate, where the `1/(1−h)` inflation and heteroscedasticity are most visible; the plain HC0 sandwich (without `1/(1−h)`) sits below `se_model` there (33.35 vs 35.91).
- `sigma_D / se_g`: 0.926–1.046, **Spearman −0.12** — no trend; `se_g` is a population, fixed-count scale and the ratio scatters around 1 by the sampling variability of a single trial's residual variance.
- `rate_closed` reproduces `max(0, 2·pnorm((beta_hat − 10)/sigma_D) − 1)` to 0 difference.

**Gate outcome: first case — ratio ≈ 1.** No constant factor is carried. `se_model` and `se_hc0` agree within 4% (recorded, not a divergence).

---

## 2. §3 — `R/fs_oc_predict.R`: the `"resample"` branch

- `stop()` replaced. New argument `pconsistency = NULL` → `forestsearch_args$pconsistency.threshold` → `formals(forestsearch)$pconsistency.threshold` (0.90); explicit values override; validated in (0, 1); recorded in `settings$pconsistency` (`NA` for `"split"`).
- Branch body: `root_full <- fs_sym_root(Sg, scale = 1)`; `Bhat <- Z %*% t(root_full) + beta_g` — **one draw matrix** (`Bhat ~ N(beta_g, Sg)`), no `W1`/`W2`; `z_p <- qnorm((1 + pconsistency)/2)`; `pass <- (Bhat >= c1) & (Bhat - c2 >= z_p * se_g)`. Everything from `det` on is shared with `"split"`, which is untouched (its block is byte-for-byte the previous code inside an `else`).
- Roxygen documents both branches, the identification with its §2 evidence, and states that `"resample"` is production and `"split"` the document's historical gate. File header updated likewise.
- Tests added (`tests/testthat/test-fs-oc-predict.R`): one-candidate family — `"resample"` `det_rate` equals `pnorm(max(c1, c2 + z_p·se), beta, se, lower = FALSE)` within 4 MC SEs at `pconsistency ∈ {0.5, 0.9, 0.99}`; `pconsistency` default chain and validation; `"resample"` at `pconsistency = 1e-6` exceeds the `"split"` rate, lies closer to it than at 0.90, and does not equal it; finite outputs and print. File: `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 89 ]` (was 63).
- **§5 fidelity harness re-run after the edit: `FIDELITY GATE: PASS (bit-identical)`** — all 13 compared quantities `identical()`.

Memory / wall-clock difference between the branches: §3 below (from the grid run).

---

## 3. §4 — the grid: `dev/glm-continuous-sims/oc_wrapper_grid_2026-08-29.R` → `oc_wrapper_grid_2026-08-29.rds`

Seed `20260825` (the document's), `draws = 2e5`, `c1 = 30`, `c2 = 10`, `pconsistency = 0.90`, the driver's arguments (asserted against each payload's `$meta`). Measured columns from `payload$oc$identification[conditional]` (detection, mean_size_H, sens, spec, ppv, npv), `payload$oc$estimation[H, oracle]` (`EbetaH = avg − bias_beta`) and `[H, naive]` (`Enaive_bias = bias_beta`), with MC SEs from `payload$results` per-replicate columns — computed, never typed.

**Families (population enumeration, driver's arguments):**

| n | cut columns | enumerated | empty | minp | rmin | size (Pg < n.min/n) | kept | duplicate | **M** | floor n.min/n | secs |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 500 | 72 | 2628 | 119 | 0 | 99 | 631 | 1779 | 178 | **1601** | 0.1200 | 8.5 |
| 700 | 72 | 2628 | 119 | 0 | 94 | 435 | 1980 | 196 | **1784** | 0.0857 | 9.4 |

Only the size floor changes with n (and `rmin/n`: 0.0100 → 0.0071), admitting 183 more candidates at n = 700.

**Wall-clock and memory per evaluation** (2×10⁵ draws):

| n | gate | M | seconds | peak MB (gc max-used, cumulative in session) |
|---|---|---:|---:|---:|
| 500 | resample | 1601 | 328 | 11,384 |
| 500 | split | 1601 | 633 | 15,068 |
| 700 | resample | 1784 | 403 | 15,161 |
| 700 | split | 1784 | 783 | 16,233 |

The `"resample"` branch runs in **about half the wall-clock** of `"split"` at both n (one `Rdraw × M` product instead of two, and no `W1`/`W2`/`W1 + W2` temporaries). The memory column is R's cumulative max-used and was not reset between evaluations, so only the first row (a fresh session) is a clean per-branch peak: 11.4 GB for `"resample"` at M = 1601 against a `Bhat` matrix of 2.56 GB; the `"split"` branch holds two such matrices plus their sum and pushed the session maximum to 15.1 GB.

**n = 500** (M = 1601; `c1 = 30`, `c2 = 10`, `pconsistency = 0.90`; MC SE in parentheses where it exists on both sides):

| quantity | analytic `"resample"` | analytic `"split"` | simulation (payload, 1000 reps) |
|---|---:|---:|---:|
| det_rate | 0.9991 (0.0001) | 1.0000 (0.0000) | 1.0000 (0.0000) |
| EnH | 70.88 | 70.83 | 72.34 (0.41) |
| Esens | 0.1638 | 0.1637 | 0.1706 (0.0035) |
| Espec | 0.8698 | 0.8699 | 0.8686 (0.0019) |
| Eppv | 0.3994 | 0.3996 | 0.4055 (0.0080) |
| Enpv | 0.6643 | 0.6643 | 0.6650 (0.0016) |
| EbetaH | 31.75 | 31.75 | 31.76 (0.11) |
| Enaive_bias | 75.44 | 75.36 | 74.27 (0.58) |
| M | 1601 | 1601 | — |

**n = 700** (M = 1784):

| quantity | analytic `"resample"` | analytic `"split"` | simulation (payload, 999 detected of 1000) |
|---|---:|---:|---:|
| det_rate | 0.9998 (0.0000) | 1.0000 (0.0000) | 0.9990 (0.0010) |
| EnH | 73.03 | 73.07 | 73.63 (0.42) |
| Esens | 0.1217 | 0.1218 | 0.1229 (0.0028) |
| Espec | 0.9048 | 0.9048 | 0.9042 (0.0015) |
| Eppv | 0.4027 | 0.4029 | 0.4016 (0.0087) |
| Enpv | 0.6621 | 0.6621 | 0.6618 (0.0012) |
| EbetaH | 31.79 | 31.79 | 31.76 (0.12) |
| Enaive_bias | 75.53 | 75.49 | 77.09 (0.56) |
| M | 1784 | 1784 | — |

Measured `EbetaH` via `oc$estimation[H, oracle]$avg − bias_beta` equals the direct mean of `orient · results$betaHhat_H` over detected replicates (31.7645 / 31.7632) — the two routes agree, so the row is computed, not back-out-typed. Measured MC SEs are `sd/sqrt(n_detected)` of the per-replicate columns (`n_harm`, `sens`, `spec`, `ppv`, `npv`, `betaHhat_H`, `nv_H_est − orient·betaHhat_H`).

**Reading.** The two gates agree with each other to the third decimal on every functional at both n — as they must here, since every candidate near the floor is many SEs above `c2` and both gates saturate on declaration (`det_rate` 0.999–1.000 on all three sides: the gates **cannot be discriminated** on this fixture). Against the simulation, E|Ĥ|, sensitivity, PPV and NPV sit within about 1–2 measured MC SEs; specificity and `EbetaH` within one; the naive bias is the one row clearly outside (75.4 vs 74.3 ± 0.6 at n = 500; 75.5 vs 77.1 ± 0.6 at n = 700), in opposite directions at the two n, which is what a fixed population family scored against a sample-realized, pre-screened family would be expected to show first on the winner's-curse term.

---

## 4. §5 — the document

`quarto/simulations/actg175/continuous/oc_wrapper_verification.qmd` — "fs_oc_predict(): Verification and Methods Note". Reads the `.rds` and the two payloads; every number in prose is inline R; sections in the task's order (interface → family → means/SEs → covariance → the two gates with the derivation and identification → selection/functionals → the two tables → reading the comparison). Writes no payload.

Rendered once with RStudio's bundled Quarto (`/usr/lib/rstudio/resources/app/bin/quarto/bin/quarto`, v1.9.38; `quarto` is not on the shell PATH): **exit 0 in 7 s**; `oc_wrapper_verification.html` (2.1 MB, embedded resources) committed beside it, as the other rendered documents in that directory are tracked.

---

## 5. §6 — close-out

`devtools::document()`: `man/fs_oc_predict.Rd` rewritten; `NAMESPACE` unchanged (no new export).

`devtools::test()` after the edit:

```
[ FAIL 0 | WARN 31 | SKIP 3 | PASS 4603 ]
```

Warning-count parity against the baseline `[ FAIL 0 | WARN 31 | SKIP 3 | PASS 4577 ]`: WARN 31 = 31 (parity), SKIP 3 = 3, FAIL 0 = 0, PASS 4577 + 26 = 4603 — the new resample tests and nothing else.

`R CMD check` (`devtools::check(document = FALSE, args = "--no-manual")`, `RSTUDIO_PANDOC` set):

```
── R CMD check results ───────────────────────────────── forestsearch 0.2.4 ────
Duration: 10m 51.4s

0 errors ✔ | 0 warnings ✔ | 0 notes ✔
```

Version `0.2.3 → 0.2.4`; `NEWS.md` entry names the `"resample"` branch and the identification. `devtools::install()`: `packageVersion("forestsearch")` → `0.2.4`.

Commits (explicit paths; no push): `aa83b3ea` task doc; `4a17641d` build, grid, document and this report; the hash line itself in the follow-up commit.

---

## 6. Verdict (ten lines)

1. §2's gate opened: `sigma_D` is `sqrt(sum(dfbeta²))` with no scaling, equal to an independent dfbeta computation to 1e-10 and within 4% of the model SE and 7% of the population `se_g` on one trial, with no prevalence trend — ratio ≈ 1, no factor carried.
2. The source's pass condition is `p.consistency >= pconsistency.threshold` after two-decimal rounding; the inversion `beta_hat >= c2 + qnorm((1+p)/2)·sigma_D` is what the wrapper implements, `z_0.90 = 1.6449`.
3. `R/fs_oc_predict.R`'s `"resample"` branch is implemented on one draw matrix (`Bhat ~ N(beta_g, Sg)`), with a `pconsistency` argument; `"split"` is untouched and the fidelity harness is still bit-identical.
4. Tests +26 (89 in the file), warning parity 31 = 31, `R CMD check` 0/0/0, version 0.2.4 installed.
5. Families: M = 1601 at n = 500 and M = 1784 at n = 700 from the same 2628 combinations; only the size floor moves.
6. `"resample"` runs in about half the wall-clock of `"split"` (328 vs 633 s; 403 vs 783 s).
7. The two gates agree to three decimals on every functional at both n; `det_rate` saturates on all sides, so this fixture cannot discriminate them.
8. Against the payloads (computed from `$oc`/`$results`, never typed), E|Ĥ|, sens, spec, PPV, NPV and E[β(Ĥ)] sit within 1–2 measured MC SEs at both n; the naive bias is off by ~1–1.5 in opposite directions.
9. The verification document reads the `.rds` and payloads, types no number, and renders in 7 s.
10. What is next is Larry's: a null/near-null cell where the gates separate, the grid/inversion on this seam, and whether the naive-bias residual is family (sample cuts, pre-screen) or gate.
