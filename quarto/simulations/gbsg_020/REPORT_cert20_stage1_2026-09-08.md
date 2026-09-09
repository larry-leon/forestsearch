# REPORT — Campaign `cert20`, Stage 1 and Gate 1

Task: `dev/tasks/TASK_cert20_2026-09-08.md`, Part T (committed as received at c5044911). Part D committed at **fb62705c** (Gate D PASS; `REPORT_defaults_flip_2026-09-08.md`). Executor: Claude Code, Linux (pop-os, 128 logical cores, 100 workers), unattended. forestsearch 0.3.5 installed from fb62705c.

## Stage 1 — the two new-n smokes

Committed template driven by `FS_S7_*` env only, 5 replicates each, `FS_S7_WORKERS=5`, campaign tag `cert20smoke` (a separate namespace, so `combine_glob` can never pool a smoke into a production stem):

```
FS_S7_FOCUS=effMaxSG FS_S7_NBHD=0.20 FS_S7_Z1Q=0.60 FS_S7_FIELD_COMPLEMENT=TRUE
FS_S7_FIELD_SCALEC=selected FS_S7_FIELD_DECOMP=TRUE FS_S7_IJ_RESIDUAL=two_term
FS_S7_FB=none FS_S7_CAMPAIGN=cert20smoke FS_S7_HR=1.50 FS_S7_NSIMS=5 FS_S7_START=1
```

| n | wall | detections | truth `marg_H` / `marg_Hc` | design prevalence | realized \|H\| | realized prevalence |
|---|---|---|---|---|---|---|
| 1000 | 79 s | 5 / 5 | 1.499 / 0.7206 | 0.30655 | 304, 316, 313, 328, 300 | 0.300–0.328 |
| 1500 | 90 s | 5 / 5 | 1.499 / 0.7206 | 0.30655 | 461, 468, 475, 492, 467 | 0.3073–0.328 |

**The 31% prevalence scales at n = 1500.** `harm_z1_quantile = 0.60` produces the same super-population prevalence 0.30655 at every n (it is a quantile rule on the covariate template, not a count), and the realized region tracks it: 0.3073–0.328 at n = 1500 against 0.300–0.328 at n = 1000 and the design 0.3065. The truth targets are identical across n — `marg_H` 1.499, `marg_Hc` 0.7206, `cde_H` 1.7087, `cde_Hc` 0.6564 — as they must be, since the DGM is the same and only the draw count changes. Realized `|Ĥ|` ranges 215–328 (n = 1000) and 350–668 (n = 1500), so the planted region scales and the identifier is not saturating.

Every construction finite on every replicate of both smokes: the harm field block, the complement field block, the field-s companions, the decompose scalars, the joint and `joint_s` pairs, and the p̂ columns. Bundles carry the full 158 columns at both n. Invariants (interval ordering, γ ∈ [0.025, 0.05], the bound↔quantile identities, p̂ validity) are gated per cell at Gate 2 rather than on 5 replicates; the checker `gate2_checks_cert20.R` was exercised end to end on the committed `e1stud` HR 1.50 n 500 bundle against `p30sgnb20j20` and returned PASS on every substantive check (the only non-passes were the three campaign-identity fields, which correctly reject a non-`cert20` bundle).

## Gate 1 — projection and the amended ceiling

**Amendment (Larry, 2026-09-08, mid-run message):** the Gate 1 ceiling is raised from **10 h to 13 h** wall and the hard timeout from **12 h to 15 h**, applied to the projection and to every per-cell re-projection, so a cell is deferred only if it would cross 13 h. Nothing else changes: same cells in the same order, same knobs, same defer order (HR 1.00 n 1500 first, then HR 1.00 n 1000), same gates. No cell had been deferred under the 10 h ceiling at the time of the amendment — the projection was already inside it — so there was nothing to re-instate.

Projection basis, anchored on a realized 100-worker campaign wall rather than on a worker-count model. The light-load per-replicate cost is `fit_mr_secs + fld_H_secs + fld_Hc_secs`:

| basis | per-replicate cost | source |
|---|---|---|
| n 500, HR 1.50, 5 workers | 55.1 s | `dflt_on` gate bundle (Part D) |
| n 1000, HR 1.50, 5 workers | 70.5 s | this Stage 1 smoke |
| n 1500, HR 1.50, 5 workers | 84.3 s | this Stage 1 smoke |
| n 500, HR 1.50, 100 workers | 111.5 s | committed `e1stud` bundle, realized wall **30 m 58 s** for 2,000 replicates |

The 100-worker contention factor at n = 500 is 111.5 / 55.1 = **2.02**, matching the handoff's independent figure (73.1 s/rep under campaign load against 36.1 s light, ×2.02). Taking it as roughly constant in n, a cell's wall projects as 1858 s × T_light(n) / 55.1. Cross-checked against arm B (J = 20, n = 1000): the same relation predicts 90 min against a realized 77–85 min, so it over-predicts by about 10% — the projection is conservative. The null/harm cost ratio comes from arm B at matched n: 134.1/165.2 = 0.81 (n 500), 281.6/325.3 = 0.87 (n 1000); 0.85 is used for n 1500.

| # | cell | projected wall |
|---|---|---|
| 1 | HR 1.50, n 1000 | 40 min |
| 2 | HR 1.75, n 1000 | 41 min |
| 3 | HR 1.50, n 1500 | 48 min |
| 4 | HR 1.75, n 1500 | 49 min |
| 5 | HR 1.00, n 500 | 26 min |
| 6 | HR 1.00, n 1000 | 35 min |
| 7 | HR 1.00, n 1500 | 41 min |
| | 2,000 replicates each | **280 min** |
| | plus one combine render per cell (~4 min) | +28 min |
| | **total** | **≈ 5 h 8 m** |

**Gate 1: GO.** The projection is 5.1 h against the amended 13 h ceiling (and was already inside the original 10 h one), so **no cell is deferred**. The driver still re-projects before each cell from the realized wall of a completed cell at the same n — nulls at 0.85 of it — and defers, listing, any cell whose projected finish would cross 13 h; the hard timeout is 15 h.

## Run configuration

Every cell: `FS_S7_FOCUS=effMaxSG FS_S7_NBHD=0.20 FS_S7_Z1Q=0.60 FS_S7_FIELD_COMPLEMENT=TRUE FS_S7_FIELD_SCALEC=selected FS_S7_FIELD_DECOMP=TRUE FS_S7_IJ_RESIDUAL=two_term FS_S7_FB=none FS_S7_CAMPAIGN=cert20 FS_S7_RETURN_RESEL=TRUE FS_S7_WORKERS=100`, J = 10 (default), seeds 8316951 + sim_id, sim_id 1–2000 as two seed-disjoint batches of 1,000 (`FS_S7_START` 1 and 1001) then `FS_S7_MODE=combine`. Stems `fs_effMaxSG_fb_mr_field_m1_h{100,150,175}_knoise0_n{500,1000,1500}_z1q60_nb20_cert20`. `.refuse_if_tracked()` live on every save; `devtools::install()` done before the run; no `load_all()`. Driver `cert20_driver.sh`, checker `gate2_checks_cert20.R`, render wrapper `render_cell.sh` (session scratchpad); per-render logs and per-cell Gate 2 logs beside them.

Comparators for the same-draws assertions, where one exists: nb20 arm B (J = 20, same DGM draws) at HR 1.50 n 1000, HR 1.00 n 500 and HR 1.00 n 1000. Cells 2, 3, 4 and 7 have none on record, as the task states.
