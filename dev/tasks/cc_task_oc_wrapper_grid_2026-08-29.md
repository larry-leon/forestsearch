# CC task — `fs_oc_grid()` and `fs_oc_invert()`: sweep (c1, c2) on one draw set, and invert for a target rate

**File:** `dev/tasks/cc_task_oc_wrapper_grid_2026-08-29.md` · **Issued:** 2026-08-29 by chat, commissioned by Larry
**Repo:** `~/Documents/GitHub/forestsearch` (pop-os) · **Branch:** `feature/glm-extension`
**Transport:** Larry places this file in `~/Downloads`; CC copies it into `dev/tasks/` and commits it first (§1).
**Follows:** `REPORT_oc_wrapper_build_2026-08-28.md`, `REPORT_oc_wrapper_fixture_run_2026-08-28.md`, `REPORT_oc_wrapper_gate_and_n700_2026-08-29.md`. Read the third before starting.

---

## ⚠ `R/` CALLOUT

Two things, both confined to files this workstream created:

1. **New file `R/fs_oc_grid.R`** — adds `fs_oc_grid()` and `fs_oc_invert()`. Additive.
2. **`R/fs_oc_predict.R`: MOVES EXISTING CODE.** The gate and the selection/functional block are factored into internal helpers so the grid uses the same code rather than a second copy. **This must change no result.** The build task's fidelity harness is the guard: re-run it and it must still be bit-identical, and `fs_oc_predict()`'s own output at fixed seed must be `identical()` to 0.2.4's on both gates. Category: *moves existing code; changes no behaviour and no method.*

**No other `R/` file may be edited.** If factoring cleanly seems to require changing what `fs_oc_predict()` computes, stop and report — do not adjust results to make the refactor tidy.

**Compute:** the same four evaluations as the 08-29 grid (~25 min total), plus threshold sweeps that cost arithmetic, not draws. No simulation study, no replicate runs.

---

## 1. Provenance and first commit — GATE

```bash
hostname; pwd
git rev-parse --abbrev-ref HEAD; git rev-parse --short HEAD
git status --porcelain
git log --oneline -6
Rscript -e 'cat("forestsearch", as.character(packageVersion("forestsearch")), "\n")'
```

*GATE:* branch `feature/glm-extension`; HEAD descends from `69c2e045`; installed `0.2.4`. Dirty tree is not a failure.

```bash
ls ~/Downloads/cc_task_oc_wrapper_grid*
cp ~/Downloads/cc_task_oc_wrapper_grid*.md dev/tasks/cc_task_oc_wrapper_grid_2026-08-29.md
git add dev/tasks/cc_task_oc_wrapper_grid_2026-08-29.md
git commit -m "task doc — fs_oc_grid() and fs_oc_invert()"
```

*GATE:* exactly one match (stem may have hyphens stripped; copy under the name above, report both). Then `devtools::test()` for the parity baseline, and record `fs_oc_predict()`'s output at a fixed seed on both gates as the **pre-refactor reference** for §2.

---

## 2. The design — why this is cheap

At fixed `n` the family, `Pg`, `beta_g`, `se_g`, `Rho` and `Sg` are all fixed, and **both gates are thresholds on draws that do not depend on `c1` or `c2`**:

```
resample:  (Bhat >= c1) & (Bhat - c2 >= z_p * se_g)
split:     (W1 + W2 >= 2*c1) & (W1 >= c2) & (W2 >= c2)
```

So one draw set serves the entire `(c1, c2)` grid and any inversion at that `n`. Only `n` forces re-enumeration and re-drawing, because the prevalence floor is `n.min/n` (M was 1601 at n = 500 and 1784 at n = 700).

**Accumulate over draw blocks rather than materialising the full matrix.** For each block, compute the gate at every grid point and accumulate per grid point: declaration count; per-member selection counts; running sums for `E|Ĥ|`'s prevalence weighting, `E[β(Ĥ)]`, and the winner's-noise mean. Every reported quantity is a function of those accumulators, so memory becomes O(block × M) instead of O(draws × M) — the 2.6 GB matrices become a knob. Block size is an argument with a stated default; `Inf` means one block.

*Consistency gate (§5):* `fs_oc_grid()` at a single grid point with one block must be `identical()` to `fs_oc_predict()` at the same seed and settings. If chunking cannot preserve that, keep the single-block path exact and document precisely where multi-block draws diverge — do not silently accept a difference.

---

## 3. `R/fs_oc_grid.R` — `fs_oc_grid()`

Suggested signature; adjust to package convention and say so if you do:

```r
fs_oc_grid(dgm, forestsearch_args, n, c1, c2,
           family = NULL, consistency_method = c("resample", "split"),
           pconsistency = NULL, draws = 2e5, block = 5e4, seed = NULL,
           verbose = FALSE)
```

- `n`, `c1`, `c2` each accept a vector. The result is the full crossing, **one row per `(n, c1, c2)`**, as a data frame or a classed object carrying one.
- Enumerate and draw **once per `n`** (and once per `consistency_method` if a vector is given); sweep all `(c1, c2)` against those draws. Report, per `n`: `M`, the floor used, the enumeration stage counts, and the wall-clock split between enumeration, drawing and sweeping.
- Columns: every quantity `fs_oc_predict()` returns, plus its MC SE where it is a proportion, plus `n`, `c1`, `c2`, `M`, `consistency_method`, `draws`, `seed`.
- `family` supplied explicitly is used for every `n` **only if** its `n` matches; otherwise stop with a clear message — a family built at one `n` has the wrong prevalence floor for another.
- A `print()` and a `summary()` method are enough; no plotting in `R/`.

---

## 4. `fs_oc_invert()`

```r
fs_oc_invert(dgm, forestsearch_args, n, target, solve_for = c("c1", "c2"),
             c1 = NULL, c2 = NULL, ...)
```

Find the threshold at which the family declaration rate equals `target` — "the `c1` giving 80% power" in the document's language, the same quantity read as type-I error under a null DGM.

- On a fixed draw set the rate is a **monotone decreasing step function** of `c1` (and of `c2`). Solve on those fixed draws: bracket first, then bisect to a stated tolerance, or evaluate a fine grid and refine. Do not call `uniroot()` on a step function without checking it converged to a genuine crossing; report the achieved rate at the returned threshold alongside the target.
- **Attainability.** As `c1 → -Inf` the rate tends to a ceiling below 1 — the other threshold still binds (`c2 + z_p * se_g` for resample; both halves clearing `c2` for split). If `target` exceeds that ceiling, return `NA` **with the ceiling value and which threshold caused it**, exactly as the document's `c1_for()` labels "infeasible (c2 ceiling …)". Do not extrapolate.
- Same treatment solving for `c2` at fixed `c1`.
- Report the MC SE of the achieved rate so the user can see the resolution the draw count supports.

---

## 5. Tests — extend `tests/testthat/test-fs-oc-predict.R` or add `test-fs-oc-grid.R`

1. **Refactor guard — GATE.** `fs_oc_predict()` output at a fixed seed, both gates, `identical()` to the §1 pre-refactor reference. And the build task's fidelity harness re-run: still bit-identical to the document's chunk.
2. **Grid/predict identity — GATE.** `fs_oc_grid()` at one `(c1, c2)`, `block = Inf`, same seed → `identical()` to `fs_oc_predict()`.
3. **Blocking invariance.** Same grid point at `block = Inf` and at a small block: agreement to within the MC SE, and identical if the RNG order is preserved. State which holds.
4. **Monotonicity.** `det_rate` non-increasing along `c1` and along `c2`, on a small family.
5. **Inversion round-trip.** `fs_oc_invert()` for a target inside the attainable range returns a `c1` whose `fs_oc_grid()` rate matches the target within tolerance.
6. **Ceiling.** A target above the ceiling returns `NA` and reports the ceiling; no silent extrapolation.
7. **Family reuse guard.** A family built at one `n` passed with a different `n` errors.

Keep them fast — small synthetic DGM, a few thousand draws. Per the v5 §9 principle, show at least one of the two GATE tests failing against an injected defect, and say which you checked.

---

## 6. Document — extend `oc_wrapper_verification.qmd`

Add one section, after the comparison tables:

- The design point of §2 stated plainly: one draw set per `n`, thresholds swept for free.
- A `c1` → declaration-rate curve at each `n`, both gates, at the driver's `c2`. **Match the existing document's figure style** (it already carries a contour plot); do not introduce a new visual idiom.
- A small table of inverted `c1` at declaration targets 0.80 / 0.90 / 0.95, per `n` and gate, with the achieved rate and `NA` where the ceiling binds.
- Two sentences on interpretation: that this is the design calculation the wrapper now supports at family level, and that under the alternative at this fixture the rate saturates so the interesting range of `c1` lies well above the driver's 30.

The document must keep typing no number — everything inline R from the stored `.rds`, which the §3/§4 precompute extends. Re-render and report exit code and time.

---

## 7. Close-out

- `devtools::document()`; `devtools::test()` — raw counts and warning-count parity against §1.
- `R CMD check` (`RSTUDIO_PANDOC` on Rscript's PATH). Status verbatim; 0/0/0 is the target.
- Version → `0.2.5`; `NEWS.md` naming both new exports and stating that `fs_oc_predict()`'s results are unchanged.
- Commit on `feature/glm-extension`, **explicit paths only**, never `git add -A`. **No push.**
- `devtools::install()`; confirm the version.
- Report `dev/glm-continuous-sims/REPORT_oc_wrapper_grid_2026-08-29.md`: the refactor guard result, the grid/predict identity result, wall-clock split (enumerate vs draw vs sweep) per `n`, the inversion table, the document's render time, commits, ten-line verdict.

---

## 8. Out of scope

- No `R/` file but the new `R/fs_oc_grid.R` and the §2 refactor of `R/fs_oc_predict.R`.
- No null-DGM path, no binary/OR path, no second family constructor, no changes to the prediction document or any driver.
- No new inversion targets beyond the declaration rate — not `E|Ĥ|`, not the classification metrics.
- No simulation study, no replicate runs.
- If the refactor cannot preserve bit-identical results, stop and report rather than accepting a difference.
