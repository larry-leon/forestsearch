# REPORT — `fs_oc_grid()` and `fs_oc_invert()`: one draw set per n, thresholds swept, declaration rate inverted

**Task:** `dev/tasks/cc_task_oc_wrapper_grid_2026-08-29.md` (commit `07aaab2a`)
**Repo / branch:** `~/Documents/GitHub/forestsearch` (pop-os), `feature/glm-extension`
**Follows:** `REPORT_oc_wrapper_gate_and_n700_2026-08-29.md` (read first).
**`R/` category:** new `R/fs_oc_grid.R` (additive); `R/fs_oc_predict.R` **move-only refactor** — draw / gate / functional / assembly blocks factored into `.fs_oc_draw()`, `.fs_oc_gate()`, `.fs_oc_functionals()`, `.fs_oc_result()`, called in the original order. No other `R/` file.

---

## 0. Provenance (§1 — GATE passed)

```
pop-os · /home/larryleon/Documents/GitHub/forestsearch · feature/glm-extension
HEAD at start 69c2e045 (descends from itself) · git status --porcelain: empty
packageVersion forestsearch 0.2.4
```

`~/Downloads/cc_task_oc_wrapper_grid*`: one match, `cc_task_oc_wrapper_grid_20260829.md` (hyphens stripped) → `dev/tasks/cc_task_oc_wrapper_grid_2026-08-29.md`, committed **`07aaab2a`**.

**Parity baseline** `devtools::test()`: `[ FAIL 0 | WARN 31 | SKIP 3 | PASS 4603 ]`.

**Pre-refactor reference** (frozen at 0.2.4, before any edit): `dev/glm-continuous-sims/prerefactor_reference_2026-08-29.R write` → `prerefactor_reference_2026-08-29.rds` — `fs_oc_predict()` at `seed = 20260829`, `draws = 2e4`, `c1 = 30`, `c2 = 10`, both gates, on the hand-built M = 3 family and on the synthetic-DGM enumeration (M = 104); the returned objects minus `$family`.

---

## 1. §2 — the refactor and its guards (GATE passed)

`R/fs_oc_predict.R`: the statements from `set.seed(seed)` through `mass_below` and the object assembly were moved, unchanged and in order, into four internal helpers:

| helper | body |
|---|---|
| `.fs_oc_draw(Sg, beta_g, M, Rdraw, gate)` | the `"resample"` (`fs_sym_root(Sg, 1)`, one `rnorm(R*M)` matrix) and `"split"` (`fs_sym_root(Sg, 2)`, `W1`, `W2`, `Bhat = (W1+W2)/2`) draw blocks — RNG consumption order unchanged |
| `.fs_oc_gate(dr, c1, c2, se_g, pcons, gate)` | the two gate expressions, verbatim |
| `.fs_oc_functionals(pass, Bhat, family, n, c1, Rdraw)` | `det`, `Bmask`/`max.col`, `P1`, `p_sel`, `det_rate`, `sel_c`, `EnH … mass_below`, verbatim |
| `.fs_oc_result(fx, family, …)` | the returned list and class, verbatim |

`fs_oc_predict()` now reads: `Rho/Sg → set.seed → .fs_oc_draw → .fs_oc_gate → .fs_oc_functionals → .fs_oc_result`. Thresholds enter only `.fs_oc_gate()`, which is what lets the grid share the gate.

**Guard 1 — 0.2.4 reference:** `prerefactor_reference_2026-08-29.R check` after the refactor:

```
hand_resample  identical: TRUE
hand_split     identical: TRUE
fam_resample   identical: TRUE
fam_split      identical: TRUE
REFACTOR GUARD: PASS (identical to the 0.2.4 reference)
```

**Guard 2 — the build task's fidelity harness** (`fidelity_fs_oc_predict_2026-08-28.R`, the document's M = 16 chunk at seed 20260825, 2×10⁵ draws) after the refactor: `FIDELITY GATE: PASS (bit-identical)` — all 13 compared quantities `identical()`.

Both guards are also tests (`test-fs-oc-grid.R` test 1 re-runs guard 1 from the frozen `.rds`; it skips under `R CMD check`, where `dev/` is not shipped).

---

## 2. §3–§4 — `R/fs_oc_grid.R`

**`fs_oc_grid(dgm, forestsearch_args, n, c1, c2, family = NULL, consistency_method = c("resample","split"), pconsistency = NULL, draws = 2e5, block = 5e4, seed = NULL, verbose = FALSE, ...)`** — signature as suggested (plus `...` to the enumerator, as `fs_oc_predict()` has). `n`, `c1`, `c2` vectors; `consistency_method` may be both; one row per `(n, gate, c1, c2)`. Per `n`: enumerate once (or use `family` if `family$n` matches the single `n` — otherwise `stop()` naming both n's; a family with no recorded `n` is refused as "unrecorded"), then per gate: `set.seed(seed)` (common random numbers across `n`), draw once, sweep every `(c1, c2)`. Returns class `fs_oc_grid`: `table` (scalars + MC SE of `det_rate`), `results` (the full per-point `fs_oc_predict` objects), `families` (M, floor, stage counts), `timing` (enumerate / draw / sweep seconds per `(n, gate)`), `settings`. `print()` and `summary()`.

- **`block = Inf`** — one block: the sweep calls `.fs_oc_gate` / `.fs_oc_functionals` / `.fs_oc_result` on the whole draw set, so a grid point **is** `fs_oc_predict()`'s computation.
- **finite `block`** — draws generated in row blocks, accumulating declaration counts, per-member pass and win counts and the winner-noise sum; every quantity is recovered from the accumulators. Memory O(block × M). **Where it diverges from one-block, precisely:** (i) the same RNG stream is laid out differently across the `Rdraw × M` matrix (a single `rnorm(R*M)` fills column 1 first; blocks fill `block × M` at a time), so blocked and one-block results are two Monte-Carlo estimates of the same quantity — they agree to MC precision, not bit-for-bit; (ii) proportions are `count/draws` rather than `mean()`/`colMeans()` (long-double accumulation), a last-bit matter. Stated in the file header, the roxygen and test 3.

**`fs_oc_invert(dgm, forestsearch_args, n, target, solve_for = c("c1","c2"), c1, c2, family, consistency_method, pconsistency, draws, seed, tol = 1e-3, ...)`** — one draw set; `rate_at(x) = mean(rowSums(gate) > 0)` is a monotone non-increasing step function of the solved-for threshold. **Ceiling first:** the rate at threshold `−Inf` (other threshold binding: `c2 + z_p·se_g` for resample, both halves ≥ `c2` for split, or the `c1` screen when solving for `c2`); `target > ceiling` → `value = NA`, `attainable = FALSE`, `ceiling` and `binding` reported, no extrapolation. Otherwise bracket from the draws' range and bisect to `tol`; the returned threshold is the largest with rate ≥ target (the step the target falls on), with `achieved`, `achieved_se = sqrt(p(1−p)/draws)`, `next_step_rate` and the iteration count. No `uniroot()`.

`devtools::document()`: `NAMESPACE` +`export(fs_oc_grid)`, +`export(fs_oc_invert)`, +`S3method(print, fs_oc_grid)`, +`S3method(summary, fs_oc_grid)`, +`S3method(print, summary.fs_oc_grid)`, +`S3method(print, fs_oc_invert)`; `man/fs_oc_grid.Rd`, `man/fs_oc_invert.Rd`.

---

## 3. §5 — tests: `tests/testthat/test-fs-oc-grid.R`

1. **Refactor guard (GATE)** — identical to the frozen 0.2.4 reference, both gates, both families. Pass.
2. **Grid/predict identity (GATE)** — `fs_oc_grid(block = Inf)` at one point `identical()` to `fs_oc_predict()`, both gates. Pass.
3. **Blocking invariance** — `block = 4000` vs `Inf` at 2×10⁴ draws: `det_rate` within 4 MC SE, `EnH` within 2%, `Eppv` within 0.02. **Agreement to MC precision holds; identity does not and is not claimed** (RNG layout differs, see §2).
4. **Monotonicity** — `det_rate` non-increasing along `c1` at every `c2` and along `c2` at every `c1`, both gates, 48-row grid.
5. **Inversion round-trip** — `fs_oc_invert(target = 0.4)`'s `c1` re-evaluated by `fs_oc_grid()` on the same draws gives `det_rate` `identical()` to `achieved`; `achieved ≥ 0.4`, within 0.01 above, `next_step_rate < 0.4`; and the `c2` solve at fixed `c1`.
6. **Ceiling** — `target = 0.999` → `NA`, `attainable = FALSE`, `binding` names `c2`, and `ceiling` `identical()` to the grid rate at `c1 = −Inf`.
7. **Family reuse guard** — a family built at n = 500 refused at n = 700, at `c(500, 700)`, in `fs_oc_invert()`, and when its `n` is unrecorded; per-`n` enumeration inside the grid gives different M and floor `60/700`.

File: `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 78 ]`.

**v5 §9 — a GATE test shown to fail against a defect.** Injecting "seed ignored" into `.fs_oc_sweep()` (the `set.seed(seed)` line replaced by a comment): `[ FAIL 8 | PASS 70 ]` — test 2's two `expect_identical` (L71, L73) fail on both gates, and tests 5–6's same-draw identities (L135, L148, L165) with them. Restored (`grep` confirms) → `[ FAIL 0 | PASS 78 ]`. (A first injection, `c1 + 1e-9`, flipped no comparison at 8000 draws and is recorded as a non-demonstrating defect.) Checked: test 2.

---

## 4. §3/§4 precompute — `oc_wrapper_grid_sweep_2026-08-29.R` (extends `oc_wrapper_grid_2026-08-29.rds`)

`fs_oc_grid(dgm, fs_args, n = c(500, 700), c1 = seq(20, 120, by = 5), c2 = 10, both gates, pconsistency = 0.90, draws = 2e5, block = Inf, seed = 20260825)` on the MD40 fixture (rebuilt and gated as before), then `fs_oc_invert()` at targets 0.80 / 0.90 / 0.95 per `(n, gate)`. The sweep's point at `c1 = 30` is asserted `identical()` to the 08-29 `fs_oc_predict()` objects (all four: `TRUE` for n = 500 resample / split and n = 700 resample / split, asserted by `stopifnot()` in the script).

**Wall-clock split** (2×10⁵ draws, `block = Inf`, 21 values of `c1`):

| n | gate | M | enumerate (s) | draw (s) | sweep (s) | per c1 (s) |
|---|---|---:|---:|---:|---:|---:|
| 500 | resample | 1601 | 8.6 | 312 | 376 | 17.9 |
| 500 | split | 1601 | 8.6 | 619 | 370 | 17.6 |
| 700 | resample | 1784 | 10.3 | 385 | 417 | 19.8 |
| 700 | split | 1784 | 10.3 | 762 | 414 | 19.7 |

The 84-row sweep took 3675 s in total: drawing once per `(n, gate)` costs 5–13 min, after which each additional `c1` costs ~18 s of arithmetic (gate, argmax, functionals on a 2×10⁵ × M matrix) and no draws. The design point holds as stated: the 21-point curve costs 6–7 min of sweeping against 5–13 min of drawing — a second draw set per threshold would have cost 21× the draw column.

**The curves.** Under the alternative every curve sits at 0.998–0.999 at `c1 = 60`; the rate first drops below 0.95 at `c1 = 80` (n = 500) and `c1 = 85` (n = 700); at `c1 = 100` it is 0.64 / 0.66 and at `c1 = 120` 0.23 / 0.22. The two gates track each other to the third decimal along the whole curve (e.g. n = 500, `c1 = 100`: 0.6414 vs 0.6410).

**Inversions** (`solve_for = "c1"` at `c2 = 10`; 2×10⁵ draws, seed 20260825; each inversion draws its own set — 490–960 s apiece, ~2.4 h for the twelve; 18–20 bisections each):

| n | gate | target | c1 (inverted) | achieved rate (MC SE) | ceiling (threshold → −∞) |
|---|---|---:|---:|---:|---:|
| 500 | resample | 0.80 | 91.90 | 0.8000 (0.0009) | 0.9991 |
| 500 | resample | 0.90 | 84.72 | 0.9000 (0.0007) | 0.9991 |
| 500 | resample | 0.95 | 78.99 | 0.9500 (0.0005) | 0.9991 |
| 500 | split | 0.80 | 91.85 | 0.8000 (0.0009) | 1.0000 |
| 500 | split | 0.90 | 84.57 | 0.9000 (0.0007) | 1.0000 |
| 500 | split | 0.95 | 78.83 | 0.9500 (0.0005) | 1.0000 |
| 700 | resample | 0.80 | 93.11 | 0.8000 (0.0009) | 0.9998 |
| 700 | resample | 0.90 | 86.52 | 0.9000 (0.0007) | 0.9998 |
| 700 | resample | 0.95 | 81.28 | 0.9500 (0.0005) | 0.9998 |
| 700 | split | 0.80 | 93.14 | 0.8000 (0.0009) | 1.0000 |
| 700 | split | 0.90 | 86.52 | 0.9000 (0.0007) | 1.0000 |
| 700 | split | 0.95 | 81.31 | 0.9500 (0.0005) | 1.0000 |

All twelve targets are attainable (ceilings 0.9991 / 0.9998 for `"resample"`, where `c2 + z_p·se_g` binds, and 1.0000 for `"split"`); the achieved rate equals the target to the draw resolution (one draw = 5×10⁻⁶), so the "next step" rate differs from the target below the fourth decimal. The inverted `c1` for 80% declaration is 92–93, for 95% it is 79–81 — the interesting range of `c1` at this fixture is 2.5–3× the driver's 30, and the two gates give the same `c1` to within 0.2.

---

## 5. §6 — the document

`oc_wrapper_verification.qmd` gains the section "Threshold sweeps and the inversion" after the tables: the design point, the wall-clock split table, the `c1` → declaration-rate curves at `c2 = 10` (base graphics, 7×5, the existing document's idiom: labelled axes, `points()` at the driver's (30, 10), dotted `abline`s at the targets), the inverted-`c1` table, and the two interpretation sentences. Every number inline from the `.rds`.

Rendered with RStudio's bundled Quarto 1.9.38 (`RSTUDIO_PANDOC` set): **exit 0 in 7 s** (both before and after the extension); `oc_wrapper_verification.html` (2.3 MB) committed beside it, as before.

---

## 6. §7 — close-out

`devtools::test()`:

```
[ FAIL 0 | WARN 31 | SKIP 3 | PASS 4681 ]
```

Parity against the baseline `[ FAIL 0 | WARN 31 | SKIP 3 | PASS 4603 ]`: WARN 31 = 31 (parity), SKIP 3 = 3, FAIL 0 = 0, PASS 4603 + 78 = 4681 — the new grid tests and nothing else.

`R CMD check` (`devtools::check(document = FALSE, args = "--no-manual")`, `RSTUDIO_PANDOC` set):

```
── R CMD check results ───────────────────────────────── forestsearch 0.2.5 ────
Duration: 10m 54.1s

0 errors ✔ | 0 warnings ✔ | 0 notes ✔
```

Version `0.2.4 → 0.2.5`; `NEWS.md` names both exports and states `fs_oc_predict()`'s results are unchanged. `devtools::install()`: `packageVersion("forestsearch")` → `0.2.5`, `fs_oc_grid` and `fs_oc_invert` exported.

Commits (explicit paths; no push): `07aaab2a` task doc; `d6b7d1ae` build, precompute, document and this report; the hash line in the follow-up commit.

---

## 7. Verdict (ten lines)

1. `fs_oc_predict()` was refactored move-only into `.fs_oc_draw` / `.fs_oc_gate` / `.fs_oc_functionals` / `.fs_oc_result`; both guards held — `identical()` to the frozen 0.2.4 outputs on both gates, and the document's chunk still bit-identical.
2. `fs_oc_grid()` sweeps `(n, c1, c2)` on one draw set per `(n, gate)`; at one point with `block = Inf` it is `identical()` to `fs_oc_predict()` (GATE passed, both gates), and the fixture sweep's `c1 = 30` points are `identical()` to the 08-29 objects.
3. Blocked draws agree with the one-block path to MC precision, not bit-for-bit — the RNG layout differs across blocks; stated in code, roxygen and test 3, not hidden.
4. `fs_oc_invert()` bisects the monotone step function on fixed draws, reports the ceiling and the binding threshold, and returns `NA` above the ceiling — no `uniroot()`, no extrapolation.
5. On the MD40 fixture the draw costs 5–13 min per `(n, gate)` and each further `c1` ~18 s: the 21-point curves cost 6–7 min of arithmetic each.
6. Under the alternative the rate saturates to `c1 ≈ 60` and falls below 0.95 only at `c1 = 80–85`; 80% declaration needs `c1 ≈ 92–93`, 95% needs 79–81; the two gates agree to 0.2 in `c1`.
7. All twelve inversions are attainable (ceilings 0.9991–1.0000), achieved rates equal their targets to the draw resolution.
8. Tests +78 with warning parity 31 = 31; `R CMD check` 0/0/0 at 0.2.5; a GATE test shown to fail against an injected seed defect.
9. The document gained the sweep section in its own figure idiom, types no number, and renders in 7 s.
10. Cost note for Larry: `fs_oc_invert()` redraws per call (~2.4 h for twelve inversions here); accepting a draw set or a grid object would make inversion free after one sweep — a design choice for a later task, not made here.
