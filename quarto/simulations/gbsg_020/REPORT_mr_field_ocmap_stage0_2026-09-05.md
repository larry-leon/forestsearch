# REPORT — OC map: Stage 0 (Discovery)

**Task:** `dev/tasks/TASK_mr_field_ocmap_2026-09-05.md` (0ef90878). Unattended per the task header; H-Q1–H-Q6 at defaults; started after the uniform task parked at Gate 1 (28022b54), whose state is untouched.
**Date:** 2026-09-05. Template at b1316d22 (`FS_S7_UNIFORM` defaults FALSE — every OC-map render runs the plain `ci_method = "field"` path, byte-identical to the s7-era template per the Stage 1 uniform record). No `R/`, template, or committed-document changes anywhere in this task.

## GATE 0: PASS — no cell dropped (all five calibrate; all five stems new)

## 0a — Identity anchors per cell

The consistency-detector coverage-sweep bundles (`mr_sweep/m1_h*_*_maxeffCons_mr_seedtab_s500/fs_mr_n*_res.rds`) are **not** seed-compatible: their meta records `seed_scheme = "pre-generated table indexed by global sim_id"` with a per-row `seed` column (e.g. 1530735852, 66010672, …), not the template's `seed_base + sim_id`. Per H-Q5 default, cells whose only committed relative is a sweep bundle get **no identity gate beyond completeness and invariants**.

| Cell | Anchor | Basis |
|---|---|---|
| 1 — h1.5, n=500, knoise0 | **`results/fs_maxeffCons_fb_mr_m1_h15_knoise0_n500_combined_1_500.rds`** (500 rows, sim_id 1–500, nb_boots 300, mr_draws 5000, seed_base 8316951, maxeffCons/resample, fs 0.2.2) | fb_mr family; seed lines quoted from `sim_fs_maxeffCons_fb_mr_m1_h15_knoise0_n500_batch_1_200.qmd` — `seed_base <- 8316951L` :294, DGM `seed = seed_base + sim_id` :469, noise `set.seed(seed_base + sim_id + 1000000L)` :485, search `seedit = seed_base + sim_id` :505 — line-identical to the template's scheme |
| 2 — h0.75, n=500 | none (sweep `m1_h08_knoise0_..._seedtab` only) | seedtab scheme |
| 3 — h1.0, n=1000 | **`results/fs_maxeffCons_fb_mr_m1_h10_knoise0_n1000_res_{1_20,21_100,101_110,111_500}.rds`** (pool = sim_id 1–500; res_1_20 meta: n_sample 1000, nb_boots 300, seed_base 8316951, fs 0.2.2) | fb_mr family; same seed lines quoted from `..._n1000_batch_1_20.qmd` :294/:469/:485/:505 |
| 4 — h1.75, knoise3, n=500 | none (nearest committed: sweep h10/h15/h20 knoise3, seedtab) | seedtab scheme |
| 5 — h1.5, n=1500 | none (sweep `fs_mr_n1500_res.rds` only, seedtab) | seedtab scheme |

Anchored comparisons (cells 1, 3) use tolerance ≤ 1e-12 on the MR (IJ) and naive columns, memberships not strings, with any cross-vintage flips enumerated and reported individually (the 2026-09-05 convention; anchors are 0.2.2-vintage, so knife-edge flips of the sims-393/780 class are possible and will be handled exactly as at Gate 2 of the s7 arc).

## 0b — Calibration and stems

`calibrate_k_inter(model = "alt", use_ahr = FALSE)` handles every cell without edits: k_inter = 0.1923 (hr 0.75, protective), 0.5679 (1.0), 1.1115 (1.5), 1.3203 (1.75). `n_super = 100000L` is a fixed knob independent of `n_sample`; `FS_S7_N` feeds `n_sample` directly; `FS_S7_KNOISE=3` engages the template's own noise-covariate block (seeded `+ 1000000L`). Stems (template rule `..._fb_mr_field_m1_h%03d_knoise%d_n%d_<campaign>`, campaign `map1`):

`h150_knoise0_n500_map1` · `h075_knoise0_n500_map1` · `h100_knoise0_n1000_map1` · `h175_knoise3_n500_map1` · `h150_knoise0_n1500_map1` — all five distinct and absent from `results/` (no collision with `s7`/`s7u`/`smoke`/`uburst`/`diag` stems; the three-digit h-tag keeps h075 ≠ h75-style ambiguity out).

## 0c — Cost anchors

- s7 campaign at 100 workers (fit + MR-IJ + field): h100 28.3 s, h175 41.7 s per replicate.
- Sweep bundles (115 workers, MR-IJ only, no field — add the s7 field increment ~15–20 s at n=500, more at larger n): h15 **n=1000: 79.2 s** (q90 94), **n=1500: 133.7 s** (q90 155); h10 **knoise3 n=500: 57.4 s**; h15 n=500 sweep elapsed 188 s / 500 reps.
- Indicative Stage 2 projection at 100 workers, 2,000 reps/cell (refined at Gate 1 from smoke): cell 1 ≈ 0.4–0.5 h, cell 2 ≈ 0.4 h, cell 3 ≈ 1.0 h, cell 4 ≈ 0.8 h, cell 5 ≈ 1.5–1.7 h — cumulative ≈ **4.1–4.6 h** against the 5 h ceiling: cell 5 is the marginal one and will be the first deferred if Gate 1 timings run high.

Gate 0 complete; proceeding to Stage 1 smoke.
