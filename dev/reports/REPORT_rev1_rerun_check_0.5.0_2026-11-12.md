# REPORT — rev1 simulations and analyses re-run under forestsearch 0.5.0 (check only), 2026-11-12

Task: `dev/tasks/TASK_rev1_rerun_check_0.5.0_2026-11-12.md`, committed as received (`c63398bb`) on `cran-prep`. Host `pop-os`; the run started 2026-11-12 21:35 and the driver finished 2026-11-13 04:15:07 (6 h 40 min 06 s; the 12 h safety cap was not reached).

The five `re-done` payloads of `manuscript_rev1` (A1, A2, S4, S1, S2) were re-generated with forestsearch 0.5.0 from the same generators, seeds and settings as the rev1 runs, then compared with the payloads at fs-post-selection HEAD. Nothing was written into fs-post-selection, and no rev1 payload was replaced. This is a check of the package; it changes no manuscript result.

## Top table

| Item | Classification | Sims run | Wall clock now | Wall clock recorded by rev1 | Displayed cells changed | Largest displayed change |
|---|---|---|--:|--:|---|---|
| A1 | SELECTION SAME (against the Linux rev1 payload: IDENTICAL) | — (one analysis) | 320 s (5 min 20 s) | 493 s (Mac, the reference); 300 s (Linux rev1-reruns) | 6 of 210; results-table cells changed 0; changed = 5 timing phrases (§6.1) + manifest `Built` | no estimate; timing: roughly 355 to 1,929 times faster → roughly 14 to 117 times faster (`speedup_phrase`) |
| A2 | SELECTION SAME (against the Linux rev1 payload: IDENTICAL) | — (one analysis) | 400 s (6 min 40 s) | 680 s (Mac, the reference); 401 s (Linux rev1-reruns) | 12 of 275; results-table cells changed 0; changed = 5 timing phrases (S7) + 6 in standalone §6.2 (timing, workers/host) + manifest `Built` | no estimate; timing: roughly 331 to 1,222 times faster → roughly 7 to 54 times faster (`speedup_phrase`) |
| S4 | IDENTICAL | 500 (1–500) | 7,402 s (2 h 03 min 22 s) | 7,495 s (1,581 + 5,907 + 7) | 0 of 209; figure data 0 of 240 | none |
| S1 | IDENTICAL | 500 (1–500) | 6,742 s (1 h 52 min 22 s) | 6,781 s (1,380 + 5,381 + 20) | 0 of 113 | none |
| S2 | IDENTICAL | 400 (1–400) | 8,302 s (2 h 18 min 22 s) | 8,242 s (4,181 + 4,041 + 20) | 0 of 135 | none |

- **Reference** = the payload in `manuscript_rev1/_payloads/` at fs-post-selection HEAD `5d165acc`. A1 and A2 there come from the Mac Studio (`REPORT_rev1_mac_timing_2026-11-11.md`), so both were also compared with the Linux payloads of the rev1-reruns commit `3a61f05` (read with `git show` into scratch).
- **Displayed cells** = every table cell, caption and inline-R value that the reader fragments listed in `REVISION_LOG.md` display, computed by each fragment's own code on a scratch copy with only the payload directory changed. Figure-data values are counted separately.
- **Wall clock now** = `quarto render` start to end, all pieces of the item in sequence (the driver's `timing.log`). Rev1 figures are as the rev1 reports record them.

## Gates

| Gate | Result | Evidence |
|---|---|---|
| 1. forestsearch on `cran-prep`, `DESCRIPTION` at HEAD `Version: 0.5.0` | PASS | branch `cran-prep`; `git show HEAD:DESCRIPTION` → `Version: 0.5.0` (HEAD `c63398bb`, the task commit) |
| 2. 0.5.0 from `/tmp/chk050/lib` in the main process and one multisession worker (`workers = 2`) | PASS | main pid 270044 and worker pid 270097: `packageVersion("forestsearch")` 0.5.0, `find.package()` and the namespace path `/tmp/chk050/lib/forestsearch`, `.libPaths()[1]` `/tmp/chk050/lib` in both |
| 3. post-conditions | FAIL as written (condition b only: the concurrent part-7 session's four commits); (a), (c) and this session's own changes (b′) PASS | see "Gate 3" below |

## Provenance

| Repository | Branch | HEAD at Step 0 | `git status --porcelain` at Step 0 |
|---|---|---|---|
| forestsearch | `cran-prep` | `05f9ec68b54bcb42bafab0ed4e89833b9baf3872` before the task commit; `c63398bb029f6b8b90cdd2b9bd4062c8cc63b414` after it | empty (before and after) |
| fs-post-selection | `main` | `5d165acce104a336e11093f1881cfdb572348462` | empty |
| fs-glms-interpretable | `main` | `dc94a17e0e93653d5d4a9fb83bacca93bb8f87ea` | empty |

- **Package.** `git archive c63398bb DESCRIPTION NAMESPACE LICENSE NEWS.md R man inst vignettes` → `/tmp/chk050/src`; `R CMD INSTALL --no-test-load --library=/tmp/chk050/lib` (16 s, "version 0.5.0"). Every archived path is as at `e3479b82` ("chore: release 0.5.0"); the two later commits touch `dev/` only. Nothing was installed into the user or site library.
- **Other packages** came from the existing libraries, read only: forestploter 1.2.0 (the A1/A2 forest-plot chunk needs ≥ 1.2.0) in the user library; future 1.80.0.
- **R / Quarto / host.** R 4.6.2, Quarto 1.10.18, `pop-os` (128 physical cores, 251 GiB).
- **Reference payloads** (md5 at fs-post-selection HEAD; each equals the md5 its rev1 report records): A1 `7b5f125a…`, A2 `97e16b63…`, S4 `3eac501b…`, S1 `28a9ea14…`, S2 `cf41cceb…`. Linux rev1 (`3a61f05`): A1 `2b1e94c7…`, A2 `64698536…`.
- **Item set.** `REVISION_LOG.md` marks exactly these five payloads `re-done` (rows at lines 881, 882, 886, 887, 889); no other payload was added.

### Generators

Each generator was taken from forestsearch at the blob its rev1 report names (`git cat-file -p <blob>`), copied to `/tmp/chk050/<item>/` and rendered there. The copies' code equals the rev1 scratch copies that produced the current payloads (`/tmp/rev1/…`, `/tmp/a002_tableS4/run/…`) line for line, once the library path and the edit markers are set aside; for S4's batch 51–500 the two also differ in the comment text kept after the marker.

| Item | Generator (forestsearch) | Blob | Report that produced the reference | Pieces run |
|---|---|---|---|---|
| A1 | `quarto/applications/gbsg/analysis_gbsg_survival_multimethod.qmd` | `fded7718` (HEAD blob is `114aa11d`; both the rev1-reruns and the Mac timing reports name `fded7718`) | `REPORT_rev1_mac_timing_2026-11-11.md` | one render |
| A2 | `quarto/applications/actg175/analysis_actg175_binary_multimethod_fixed_family.qmd` | `fc59934b` (= HEAD; both reports name it) | `REPORT_rev1_mac_timing_2026-11-11.md` | one render |
| S4 | `quarto/simulations/gbsg_020/sim_fs_maxeffCons_fb_mr_m1_h10_knoise0_n500_batch_1_50.qmd` + `…_n500_combine_1_500.qmd` | `d9446e92` / `52de6afa` (= HEAD and `808e284d`; the a002 report names the paths, not the blobs) | `REPORT_a002_table_s4_2026-11-10.md` | batch 1–50, batch 51–500, combine |
| S1 | `…/gbsg_020/sim_fs_maxeffCons_fb_mr_m1_h10_knoise0_n1000_batch_1_20.qmd` + `…_n1000_combine_1_500.qmd` | `f0121563` / `e9ba13a8` (= HEAD) | `REPORT_rev1_reruns_2026-11-10.md` | batch 1–50, batch 51–500, combine |
| S2 | `…/gbsg_020/sim_fs_maxeffCons_fb_mr_m1_h15_knoise0_n1000_batch_1_200.qmd` + `…_h15_knoise0_n1000_combine_1_400.qmd` | `74b0e8eb` / `8dc86a4b` (= HEAD) | `REPORT_rev1_reruns_2026-11-10.md` | batch 1–200, batch 201–400, combine |

### `# CHK050:` edits (all of them)

Every line that differs from the blob carries `# CHK050:` (checked with `diff` per file).

- **`.libPaths()`**, all 14 copies: `.libPaths(c("/tmp/chk050/lib", .libPaths()))` as the first line of the first R chunk (`setup`): A1 line 94, A2 line 41, every sim document line 67.
- **Worker count ≤ 100.**
  - A1, A2: `n_cores <- min(100L, n_cores)` added after the formula (A1 line 128, A2 line 58); the formula `floor(0.80 × 128)` gives 102. Both renders print "Using 100 of 128 total cores".
  - S4 (all four copies, as the a002 run did it): line 296 `n_workers <- min(100L, ceiling(0.90*max(1L, min(parallel::detectCores(logical = FALSE) - 1L))))`.
  - S1, S2 (all copies, as the rev1-reruns run did it): `n_workers <- min(100L, n_workers)` added at line 297. The template formula gives 115.
- **Batch settings** (line 93 `n_sims`, line 115 `sim_id_start`), as the rev1 runs used them:
  - S4: batch 1–50 keeps the template's 50 / 1 (no edit, as in the a002 run); batch 51–500 450 / 51.
  - S1: 50 / 1 and 450 / 51 (template 20 / 1).
  - S2: 200 / 1 (the template's values, marked) and 200 / 201.
  - Smoke copies (`/tmp/chk050/<item>/smoke/`): 2 / 1.
- **Output paths:** no edit. Every generator writes relative to its own directory (`results/…` for the sims, `_payloads/<stem>/` for A1/A2), so all output landed under `/tmp/chk050/<item>/`. The smoke copies sit in their own `smoke/` folder so that the combiners' `results/<stem>_res_*.rds` glob cannot pick up their sims 1–2.
- Nothing else changed: seeds, DGM, sample sizes, `nb_boots` 300, `mr_draws` 5,000, thresholds, `run_cv`/`run_loo` `FALSE`, B = 1,000 for A1/A2.

## Run record

| Item | Piece | Start → end | Exit | Wall clock now (s) | Rev1 recorded (s) |
|---|---|---|---|--:|--:|
| A1 | analysis | 21:35:01 → 21:40:21 | 0 | 320 | 300 (Linux rev1-reruns); 493 (Mac, reference) |
| A2 | analysis | 21:40:21 → 21:47:01 | 0 | 400 | 401 (Linux rev1-reruns); 680 (Mac, reference) |
| S1 | smoke, sims 1–2 | 21:47:01 → 21:58:21 | 0 | 680 | — |
| S4 | smoke, sims 1–2 | 21:47:01 → 22:00:01 | 0 | 780 | — |
| S2 | smoke, sims 1–2 | 21:47:01 → 22:01:01 | 0 | 840 | — |
| S4 | sims 1–50 | 22:01:01 → 22:27:02 | 0 | 1,561 | 1,581 |
| S4 | sims 51–500 | 22:27:02 → 00:04:03 | 0 | 5,821 | 5,907 |
| S4 | combine | 00:04:03 → 00:04:23 | 0 | 20 | 7 |
| S1 | sims 1–50 | 00:04:23 → 00:26:43 | 0 | 1,340 | 1,380 |
| S1 | sims 51–500 | 00:26:43 → 01:56:25 | 0 | 5,382 | 5,381 |
| S1 | combine | 01:56:25 → 01:56:45 | 0 | 20 | 20 |
| S2 | sims 1–200 | 01:56:45 → 03:07:46 | 0 | 4,261 | 4,181 |
| S2 | sims 201–400 | 03:07:46 → 04:14:47 | 0 | 4,021 | 4,041 |
| S2 | combine | 04:14:47 → 04:15:07 | 0 | 20 | 20 |

- One driver (`/tmp/chk050/driver.sh`) ran every render in task order, logged `timing.log` and held the 12 h safety cap (deadline 2026-11-13 09:35:01). The cap was not reached: driver done 2026-11-13T04:15:07-08:00 capped=0.
- **Smoke.** The three smoke renders ran side by side (each with the 100-worker plan; two replicates each). Each passed: exit 0, both replicates `DETECTED`, `err_msg` empty, and sims 1–2 equal to the reference payload in every non-timing column.
- **Payload copies** in `~/fs050_check/payloads/`, under the original names:

| Item | File | Bytes | md5 | `forestsearch_version` |
|---|---|--:|---|---|
| A1 | `analysis_gbsg_survival_multimethod_payload.rds` | 1,910 | `2a27fc2661d375431e6c61992b3f771e` | 0.5.0 |
| A2 | `analysis_actg175_binary_multimethod_fixed_family_payload.rds` | 1,917 | `4c23540129d8e432fb459d8786ae7085` | 0.5.0 |
| S4 | `fs_maxeffCons_fb_mr_m1_h10_knoise0_n500_combined_1_500.rds` | 114,010 | `810c2b2ec3b9bfb541cd3d8e7d547c50` | 0.5.0 |
| S1 | `fs_maxeffCons_fb_mr_m1_h10_knoise0_n1000_combined_1_500.rds` | 113,434 | `d52d8d7a133ea3a4dd9d4d8164ffe36e` | 0.5.0 |
| S2 | `fs_maxeffCons_fb_mr_m1_h15_knoise0_n1000_combined_1_400.rds` | 116,942 | `82f995115d65495fb2637a06dce01f0b` | 0.5.0 |

## Comparisons

### A1 — GBSG survival application: SELECTION SAME (against the Linux rev1 payload: IDENTICAL)

- Against the reference (Mac Studio payload): the selected subgroups (FS, DINA, GRF, DINA standalone) and the subgroup sizes are identical; estimates differ in the last bits only (largest 4.0e-15). Against the Linux rev1 payload: every estimate field and every subgroup definition is identical at full precision.
- The differences from the Mac payload are, field for field, the Mac-versus-Linux differences that `REPORT_rev1_mac_timing_2026-11-11.md` (estimate gate) records between its Mac payload and the Linux one. They are a platform difference, not a code difference.
- Fields moved against the reference: naive limits, FB limits, MR estimate and limits (last bits). No displayed estimate changed. The changed displayed entries are the timing phrases computed from `meta$timings` (host and worker count) and the manifest's build stamp.

**Displayed table (`_supp_gbsg_table.qmd`, supplement Table S2), new beside the reference.** A cell shows one value when new = reference and **new** (reference) when it changed.

| Method | SG | n (%) | Naive | FB | MR | changed |
|---|---|---|---|---|---|---|
| FS | Ĥ | 61 (9) | 2.54 (1.25, 5.17) | 1.45 (0.59, 3.58) | 1.44 (0.62, 3.36) | 0 |
|  | Ĥᶜ | 625 (91) | 0.61 (0.47, 0.79) | 0.65 (0.38, 1.09) | 0.64 (0.40, 1.03) | 0 |
| DINA | Ĥ | 89 (13) | 1.60 (0.88, 2.91) | 0.89 (0.38, 2.07) | 1.11 (0.45, 2.75) | 0 |
|  | Ĥᶜ | 597 (87) | 0.63 (0.48, 0.82) | 0.66 (0.40, 1.11) | 0.66 (0.40, 1.09) | 0 |
| GRF | Ĥ | 82 (12) | 1.95 (1.04, 3.67) | 1.04 (0.49, 2.24) | 1.05 (0.47, 2.34) | 0 |
|  | Ĥᶜ | 604 (88) | 0.61 (0.47, 0.80) | 0.66 (0.39, 1.13) | 0.66 (0.41, 1.06) | 0 |

- `_supp_gbsg_table.qmd`: 37 displayed entries (36 table cells), 0 changed. Against the Linux rev1 payload: 0 of 37 changed.
- `_gbsg_section6_1.qmd`: 17 displayed entries (0 table cells), 5 changed. Against the Linux rev1 payload: 3 of 17 changed.

  | Part | Object / expression | New | Reference | |
  |---|---|---|---|---|
  | inline | `boot_min` | about 42 seconds | about 74 seconds | changed |
  | inline | `boot_max` | about 1.8 minutes | about 4.1 minutes | changed |
  | inline | `mr_min` | 0.36 seconds | 0.11 seconds | changed |
  | inline | `mr_max` | about 5 seconds | 0.39 seconds | changed |
  | inline | `speedup_phrase` | roughly 14 to 117 times faster | roughly 355 to 1,929 times faster | changed |
- `payload_manifest.qmd`: 156 displayed entries (156 table cells), 1 changed.
  Changed: `analysis_gbsg_survival_multimethod_payload.rds` Built: 2026-11-12 21:39 (reference 2026-11-11 09:55)

Per-field differences:

| Field | vs reference (Mac): rows differing of 6 | max abs diff | vs Linux rev1: rows differing | max abs diff |
|---|--:|--:|--:|--:|
| `n` | 0 | 0 | 0 | 0 |
| `pct` | 0 | 0 | 0 | 0 |
| `naive_est` | 0 | 0 | 0 | 0 |
| `naive_lo` | 4 | 4.4e-16 | 0 | 0 |
| `naive_hi` | 4 | 1.8e-15 | 0 | 0 |
| `fb_est` | 0 | 0 | 0 | 0 |
| `fb_lo` | 6 | 1.1e-16 | 0 | 0 |
| `fb_hi` | 5 | 8.9e-16 | 0 | 0 |
| `mr_est` | 1 | 1.1e-16 | 0 | 0 |
| `mr_lo` | 4 | 6.7e-16 | 0 | 0 |
| `mr_hi` | 5 | 4e-15 | 0 | 0 |

- Selected-subgroup definitions (`labels$fs`, `$dina`, `$grf`, `$dina_standalone`): 4 of 4 identical to the reference, 4 of 4 to the Linux rev1 payload. New: fs `{er <= 0} & {size <= 35}`; dina `{grade3 >= 1} & {pgr <= 10}`; grf `{er <= 0}`; dina_standalone `TRUE; grade3 >=  1.0; pgr >= 32.5; 127`.
- Other content (`est_scale`, `extras`, `meta$B`, `meta$mr_draws`, `meta$n_total`, covariates, outcome, event): 8 of 8 identical to the reference.

- `_gbsg_section6_1_legacy.qmd` (standalone) does not run, on either payload: `Error in vapply(P$meta$timings, function(m) m$gate_sec, numeric(1)): values must be length 1, but FUN(X[[1]]) result is length 0`. It reads `gate_sec`, which neither the reference nor the new payload records. It is not tabulated.

### A2 — ACTG175 binary application: SELECTION SAME (against the Linux rev1 payload: IDENTICAL)

- Against the reference (Mac Studio payload): the selected subgroups and the subgroup sizes are identical; FB differs by up to 1.07e-4 and the rest by up to 1.4e-14. Against the Linux rev1 payload: every estimate field and every subgroup definition is identical at full precision.
- As for A1, these are the Mac-versus-Linux differences recorded in the Mac timing report (FB `fb_est` 1.07e-04, `fb_lo` 7.1e-05, `fb_hi` 1.04e-04): a platform difference, not a code difference.
- Fields moved against the reference: naive, FB and MR estimates and limits (FB at the fourth decimal, the rest in the last bits). No displayed estimate changed. The changed displayed entries are the timing phrases, the standalone §6.2 worker/host phrase, and the manifest's build stamp.

**Displayed table (`_actg175_supp.qmd`, supplement Table S3), new beside the reference.** A cell shows one value when new = reference and **new** (reference) when it changed.

| Method | SG | n (%) | Naive | FB | MR | changed |
|---|---|---|---|---|---|---|
| FS | Ĥ | 72 (7) | 3.58 (1.31, 9.77) | 1.85 (0.53, 6.50) | 1.79 (0.63, 5.09) | 0 |
|  | Ĥᶜ | 1011 (93) | 0.59 (0.45, 0.76) | 0.62 (0.38, 1.02) | 0.62 (0.37, 1.03) | 0 |
| DINA | Ĥ | 98 (9) | 1.94 (0.84, 4.47) | 0.87 (0.25, 2.97) | 0.95 (0.29, 3.10) | 0 |
|  | Ĥᶜ | 985 (91) | 0.60 (0.46, 0.78) | 0.63 (0.40, 1.02) | 0.63 (0.38, 1.06) | 0 |
| GRF | Ĥ | 95 (9) | 2.25 (0.97, 5.20) | 1.09 (0.40, 3.01) | 1.10 (0.43, 2.80) | 0 |
|  | Ĥᶜ | 988 (91) | 0.58 (0.45, 0.76) | 0.63 (0.38, 1.04) | 0.62 (0.37, 1.04) | 0 |

- `_actg175_supp.qmd`: 59 displayed entries (36 table cells), 5 changed. Against the Linux rev1 payload: 4 of 59 changed.

  | Part | Object / expression | New | Reference | |
  |---|---|---|---|---|
  | inline | `boot_min` | about 39 seconds | about 64 seconds | changed |
  | inline | `boot_max` | about 2.4 minutes | about 6.3 minutes | changed |
  | inline | `mr_min` | 0.72 seconds | 0.15 seconds | changed |
  | inline | `mr_max` | about 12 seconds | 0.56 seconds | changed |
  | inline | `speedup_phrase` | roughly 7 to 54 times faster | roughly 331 to 1,222 times faster | changed |
- `_actg175_section6_2.qmd`: 60 displayed entries (36 table cells), 6 changed. Against the Linux rev1 payload: 4 of 60 changed.

  | Part | Object / expression | New | Reference | |
  |---|---|---|---|---|
  | inline | `cores_phrase` | 100 workers (AMD Ryzen Threadripper PRO 5995WX 64-Cores) | 11 workers (Apple M4 Max, Accelerate BLAS) | changed |
  | inline | `boot_min` | 38.8 seconds | 63.7 seconds | changed |
  | inline | `boot_max` | 2.4 minutes | 6.3 minutes | changed |
  | inline | `mr_min` | 0.72 seconds | 0.15 seconds | changed |
  | inline | `mr_max` | 12.0 seconds | 0.56 seconds | changed |
  | inline | `speedup_phrase` | roughly 7 to 54 times faster | roughly 331 to 1,222 times faster | changed |
- `payload_manifest.qmd`: 156 displayed entries (156 table cells), 1 changed.
  Changed: `analysis_actg175_binary_multimethod_fixed_family_payload.rds` Built: 2026-11-12 21:46 (reference 2026-11-11 10:06)

Per-field differences:

| Field | vs reference (Mac): rows differing of 6 | max abs diff | vs Linux rev1: rows differing | max abs diff |
|---|--:|--:|--:|--:|
| `n` | 0 | 0 | 0 | 0 |
| `pct` | 0 | 0 | 0 | 0 |
| `naive_est` | 6 | 7e-15 | 0 | 0 |
| `naive_lo` | 5 | 5.6e-15 | 0 | 0 |
| `naive_hi` | 6 | 1.2e-14 | 0 | 0 |
| `fb_est` | 6 | 0.00011 | 0 | 0 |
| `fb_lo` | 6 | 7.1e-05 | 0 | 0 |
| `fb_hi` | 6 | 0.0001 | 0 | 0 |
| `mr_est` | 5 | 7.3e-15 | 0 | 0 |
| `mr_lo` | 5 | 4.9e-15 | 0 | 0 |
| `mr_hi` | 5 | 1.4e-14 | 0 | 0 |

- Selected-subgroup definitions (`labels$fs`, `$dina`, `$grf`, `$dina_standalone`): 4 of 4 identical to the reference, 4 of 4 to the Linux rev1 payload. New: fs `!{cd40 <= 380} & !{wtkg <= 86}`; dina `{preanti >= 849.40000000000009} & {cd40 >= 338}`; grf `{wtkg > 84.369600000000005} & {cd40 > 368.20000000000005}`; dina_standalone `TRUE; cd40 >= 500; drugs >=   0; 100`.
- Other content (`est_scale`, `extras`, `meta$B`, `meta$mr_draws`, `meta$n_total`, covariates, outcome, event): 8 of 8 identical to the reference.

### S4 — Table S4 simulation, n = 500, h = 1.0, sims 1–500: IDENTICAL

- Simulations paired by `sim_id`: 500 of 500 (reference 500). Detection differs in 0; the selected subgroup definition (`sg_def`) differs in 0; `n_sel`/`n_harm` differ in 0.
- Every compared value is identical, in every simulation; `truth` is identical; the only `meta` difference is `forestsearch_version` (0.4.1.9003 → 0.5.0).

**Displayed table (`_sim_estimation_fs_h10.qmd`), new beside the reference.** A cell shows one value when new = reference and **new** (reference) when it changed.

| Block | Estimator | Avg | SD_β | SD̂_β | Length | b^oracle | b^θ | b^‡ | b^† | C^oracle | C^θ | C^‡ | C^† | changed |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Ĥ | oracle | 1.13 | 0.31 | 0.33 | 1.55 | +0.0 | — | +13.3 | +13.2 | 1.00 | — | 0.96 | 0.96 | 0 |
| Ĥ | naive | 1.63 | 0.17 | 0.32 | 2.22 | +55.7 | +129.3 | +62.9 | +62.8 | 0.78 | 0.13 | 0.86 | 0.86 | 0 |
| Ĥ | FB | 0.95 | 0.21 | 0.39 | 1.65 | -9.6 | +33.5 | -4.9 | -5.0 | 0.96 | 0.99 | 1.00 | 1.00 | 0 |
| Ĥ | MR | 0.94 | 0.24 | 0.39 | 1.64 | -11.5 | +31.2 | -6.4 | -6.5 | 0.97 | 0.98 | 0.99 | 0.99 | 0 |
| Ĥᶜ | oracle | 0.66 | 0.13 | 0.14 | 0.36 | +0.0 | — | +12.4 | +0.0 | 1.00 | — | 0.89 | 0.96 | 0 |
| Ĥᶜ | naive | 0.60 | 0.13 | 0.14 | 0.33 | -7.9 | -8.5 | +3.3 | -8.0 | 1.00 | 0.91 | 0.95 | 0.92 | 0 |
| Ĥᶜ | FB | 0.66 | 0.14 | 0.25 | 0.68 | +0.1 | -0.5 | +12.4 | +0.1 | 1.00 | 1.00 | 0.99 | 1.00 | 0 |
| Ĥᶜ | MR | 0.65 | 0.14 | 0.25 | 0.67 | -1.2 | -1.9 | +10.9 | -1.3 | 1.00 | 1.00 | 0.99 | 1.00 | 0 |

- `_sim_estimation_fs_h10.qmd`: 113 displayed entries (104 table cells), 0 changed.
- `_sim_synthesis_h10.qmd`: 16 displayed entries (0 table cells), 0 changed. Figure data (`fig_df`, 120 plotted values): 0 changed.
- `postselection_subgroup_inference_main.qmd`: 5 displayed entries (0 table cells), 0 changed.
- `_sim_additional_results_intro.qmd`: 21 displayed entries (0 table cells), 0 changed.
- `_supp_sim_synthesis.qmd`: 33 displayed entries (0 table cells), 0 changed. Figure data (`fig_df`, 120 plotted values): 0 changed.
- `_survey_2026-08/_sim_additional_results_intro.qmd`: 21 displayed entries (0 table cells), 0 changed.

Per-simulation differences:

Columns compared: 54 (every shared column except `sim_id`, `fb_secs`, `fit_mr_secs`); columns with any difference: 0.

### S1 — n = 1,000, h = 1.0, sims 1–500: IDENTICAL

- Simulations paired by `sim_id`: 500 of 500 (reference 500). Detection differs in 0; the selected subgroup definition (`sg_def`) differs in 0; `n_sel`/`n_harm` differ in 0.
- Every compared value is identical, in every simulation; `truth` is identical; the only `meta` difference is `forestsearch_version` (0.4.1.9003 → 0.5.0).

**Displayed table (`_sim_estimation_fs_h10_n1000.qmd`), new beside the reference.** A cell shows one value when new = reference and **new** (reference) when it changed.

| Block | Estimator | Avg | SD_β | SD̂_β | Length | b^oracle | b^θ | b^‡ | b^† | C^oracle | C^θ | C^‡ | C^† | changed |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Ĥ | oracle | 1.10 | 0.22 | 0.23 | 1.01 | +0.0 | — | +10.2 | +10.2 | 1.00 | — | 0.94 | 0.94 | 0 |
| Ĥ | naive | 1.39 | 0.13 | 0.25 | 1.45 | +30.8 | +84.9 | +38.9 | +38.8 | 0.88 | 0.30 | 0.91 | 0.91 | 0 |
| Ĥ | FB | 0.90 | 0.17 | 0.30 | 1.16 | -15.7 | +19.5 | -10.0 | -10.1 | 0.94 | 0.98 | 0.98 | 0.98 | 0 |
| Ĥ | MR | 0.92 | 0.18 | 0.30 | 1.16 | -14.1 | +21.8 | -8.3 | -8.4 | 0.95 | 0.99 | 0.99 | 0.99 | 0 |
| Ĥᶜ | oracle | 0.65 | 0.10 | 0.10 | 0.25 | +0.0 | — | +10.4 | -1.7 | 1.00 | — | 0.83 | 0.95 | 0 |
| Ĥᶜ | naive | 0.63 | 0.09 | 0.10 | 0.24 | -2.5 | -4.1 | +7.5 | -4.3 | 1.00 | 0.95 | 0.91 | 0.93 | 0 |
| Ĥᶜ | FB | 0.66 | 0.10 | 0.18 | 0.48 | +2.3 | +0.6 | +12.8 | +0.4 | 1.00 | 1.00 | 0.98 | 1.00 | 0 |
| Ĥᶜ | MR | 0.66 | 0.10 | 0.18 | 0.47 | +2.0 | +0.4 | +12.5 | +0.1 | 1.00 | 1.00 | 1.00 | 1.00 | 0 |

- `_sim_estimation_fs_h10_n1000.qmd`: 113 displayed entries (104 table cells), 0 changed.

Per-simulation differences:

Columns compared: 54 (every shared column except `sim_id`, `fb_secs`, `fit_mr_secs`); columns with any difference: 0.

### S2 — n = 1,000, h = 1.5, sims 1–400: IDENTICAL

- Simulations paired by `sim_id`: 400 of 400 (reference 400). Detection differs in 0; the selected subgroup definition (`sg_def`) differs in 0; `n_sel`/`n_harm` differ in 0.
- Every compared value is identical, in every simulation; `truth` is identical; the only `meta` difference is `forestsearch_version` (0.4.1.9003 → 0.5.0).

**Displayed table (`_sim_estimation_fs_h15_n1000.qmd`), new beside the reference.** A cell shows one value when new = reference and **new** (reference) when it changed.

| Block | Estimator | Avg | SD_β | SD̂_β | Length | b^oracle | b^θ | b^‡ | b^† | C^oracle | C^θ | C^‡ | C^† | changed |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Ĥ | oracle | 1.61 | 0.21 | 0.22 | 1.42 | +0.0 | — | -3.6 | +6.8 | 1.00 | — | 0.94 | 0.94 | 0 |
| Ĥ | naive | 1.65 | 0.17 | 0.24 | 1.59 | +4.3 | +42.4 | -1.1 | +9.6 | 1.00 | 0.74 | 0.98 | 0.97 | 0 |
| Ĥ | FB | 1.16 | 0.23 | 0.31 | 1.51 | -27.5 | -0.9 | -30.4 | -22.9 | 0.91 | 0.93 | 0.77 | 0.86 | 0 |
| Ĥ | MR | 1.20 | 0.25 | 0.30 | 1.52 | -25.5 | +1.8 | -28.4 | -20.7 | 0.95 | 0.94 | 0.76 | 0.85 | 0 |
| Ĥᶜ | oracle | 0.63 | 0.10 | 0.10 | 0.25 | +0.0 | — | +8.3 | -3.6 | 1.00 | — | 0.88 | 0.94 | 0 |
| Ĥᶜ | naive | 0.64 | 0.10 | 0.10 | 0.24 | +0.6 | -3.0 | +8.9 | -3.1 | 1.00 | 0.95 | 0.87 | 0.94 | 0 |
| Ĥᶜ | FB | 0.66 | 0.10 | 0.18 | 0.49 | +4.9 | +1.1 | +13.6 | +1.1 | 1.00 | 1.00 | 0.98 | 1.00 | 0 |
| Ĥᶜ | MR | 0.66 | 0.10 | 0.18 | 0.48 | +4.5 | +0.8 | +13.2 | +0.8 | 1.00 | 1.00 | 0.99 | 1.00 | 0 |

- `_sim_estimation_fs_h15_n1000.qmd`: 114 displayed entries (104 table cells), 0 changed.
- `_sim_additional_results_intro.qmd`: 21 displayed entries (0 table cells), 0 changed.

Per-simulation differences:

Columns compared: 54 (every shared column except `sim_id`, `fb_secs`, `fit_mr_secs`); columns with any difference: 0.

## Errors

- **Renders:** none failed. All 14 renders exited 0 (A1, A2, three smokes, and three pieces each for S4, S1, S2); no item was skipped or capped.
- **Replicates:** `err_msg` and `fb_err` are empty in every simulation of the three new payloads (1,400 rows), as in the references; FB returned an estimate on every detected replicate (S4 329, S1 329, S2 384).
- **Reader fragment `_gbsg_section6_1_legacy.qmd`** (standalone, reads the A1 payload) stops on the reference payload and on the new one alike:
  `Error in vapply(P$meta$timings, function(m) m$gate_sec, numeric(1)): values must be length 1, but FUN(X[[1]]) result is length 0` (call: `vapply(P$meta$timings, function(m) m$gate_sec, numeric(1))`, fragment line 124). Neither payload records `gate_sec`.
- **Messages that are not errors**, identical on both payloads: `_supp_sim_synthesis.qmd` warns that `geom_point()` removed 6 rows (the oracle's conditional-estimand entries, which the fragment sets to `NA` on purpose); the two h = 1.0 estimation fragments print "beta(H) source: fallback constant (fs sim doc, DGM eval frame) (1.002890)".

## Gate 3

Run by `/tmp/chk050/gate3.sh staged` on the staged report immediately before the report commit. After the commit it was re-run as `gate3.sh committed <report commit>`; the closing message prints that run. Staged output:

```text
a. fs-post-selection: branch main, HEAD 5d165acce104a336e11093f1881cfdb572348462, porcelain lines 0 (Step 0: main, 5d165acce104a336e11093f1881cfdb572348462, 0) -> PASS
a. fs-glms-interpretable: branch main, HEAD dc94a17e0e93653d5d4a9fb83bacca93bb8f87ea, porcelain lines 0 (Step 0: main, dc94a17e0e93653d5d4a9fb83bacca93bb8f87ea, 0) -> PASS
b. forestsearch (staged): branch cran-prep; files changed since 05f9ec68 (parent of the Step 0 commit): dev/reports/REPORT_cran_prep_policy_fixes_2026-11-12.md dev/reports/REPORT_rev1_rerun_check_0.5.0_2026-11-12.md dev/tasks/TASK_cran_prep_policy_fixes_2026-11-12.md dev/tasks/TASK_rev1_rerun_check_0.5.0_2026-11-12.md NEWS.md R/mrct_simulation.R tests/testthat/test-mrct-case-study.R tests/testthat/test-mrct-handlers-restored.R ; other porcelain lines 0 -> FAIL
   commits since the Step 0 commit: da43a29a dev: CRAN preparation report, part 7;f6bb74a7 fix: mrct_region_sims() restores the caller's progressr.handlers (part 7, fix 2);4482a0d9 test: cap test-mrct-case-study.R at 2 workers (part 7, fix 1);9c94727b dev: task part 7, CRAN policy fixes (as received);
b'. files changed by this session (Step 0 commit c63398bb + the staged report): dev/reports/REPORT_rev1_rerun_check_0.5.0_2026-11-12.md dev/tasks/TASK_rev1_rerun_check_0.5.0_2026-11-12.md  -> PASS
c. forestsearch: refs/remotes + refs/tags equal to the Step 0 snapshot (9 refs) -> PASS
c. fs-post-selection: refs/remotes + refs/tags equal to the Step 0 snapshot (4 refs) -> PASS
c. fs-glms-interpretable: refs/remotes + refs/tags equal to the Step 0 snapshot (3 refs) -> PASS
c. forestsearch: merge commits in 05f9ec68..HEAD: 0 -> PASS
GATE 3: FAIL
```

- **(a) fs-post-selection and fs-glms-interpretable: PASS.** Branch, HEAD and `git status --porcelain` (empty) are as at Step 0.
- **(b) forestsearch: FAIL as written.** Since the Step 0 commit, files other than this task file and the report have changed. They come from the four part-7 commits of the concurrent session (`9c94727b`, `4482a0d9`, `f6bb74a7`, `da43a29a`; see Notes): `R/mrct_simulation.R`, `NEWS.md`, `tests/testthat/test-mrct-case-study.R`, `tests/testthat/test-mrct-handlers-restored.R`, and that task's own task and report files.
- **(b′) Restricted to this session's own changes: PASS.** The Step 0 commit `c63398bb` holds only the task file, and the report commit holds only the report. Nothing else in forestsearch was written, and the working tree was otherwise clean.
- **(c) Nothing pushed, merged or tagged: PASS.** `refs/remotes` and `refs/tags` of all three repositories equal the Step 0 snapshots, and `05f9ec68..HEAD` holds no merge commit.
- **Overall: FAIL as written, on condition (b) only.** Gates stop on failure, so Step 5.4 (deleting `/tmp/chk050/lib` and `/tmp/chk050/src`) was not run.

## Notes

- **Classification rule.** Each item is classified against its reference payload at fs-post-selection HEAD. Timing columns and fields (`fb_secs`, `fit_mr_secs`, `meta$timings`), host fields (`meta$machine`, `meta$cores`), `built_at` and the package version are not compared values.
- **Comparison tools** (scratch, `/tmp/chk050/compare/`): `extract.R` is rev1's `extract_tables.R` extended to capture `knitr::kable()` inputs (so qualified calls in `payload_manifest.qmd` count), inline values, captions and the synthesis fragments' `fig_df`, and to skip `eval: false` chunks. On the Linux rev1 payloads it reproduces rev1's own extractions of Tables S2, S7 and S8 byte for byte. The fragments are read from a `git archive` of fs-post-selection HEAD, so nothing ran inside the repository.
- **Standalone `_survey_2026-08/_sim_additional_results_intro.qmd`** reads `.pdir <- "_payloads/"` relative to its own folder, where no payloads exist. Its scratch copies point at the payload directory instead, for both runs.
- **`payload_manifest.qmd`** (metadata only) is compared per application item: its one payload-driven cell per application is the `Built` timestamp.
- **Concurrent session.** While this task ran, the session for `TASK_cran_prep_policy_fixes_2026-11-12.md` (CRAN preparation, part 7) committed four commits on `cran-prep`: `9c94727b` 21:41:16, `4482a0d9` 21:46:16, `f6bb74a7` 21:51:58, `da43a29a` 21:53:46. They change `R/mrct_simulation.R`, `NEWS.md`, two test files and that task's own task and report files. That task file says another session is running simulations here (`/tmp/chk050`) and tells its session to leave this session's commits and files alone. This task's package build was archived from `c63398bb` at 21:35, before those commits, so they are not in the tested build. Those commits fall within this task's A2 render and smoke runs (21:40–22:01); that task file limits its tests and reproducers to 2 workers and 2 threads.
- **Rendered generator HTML, logs and scratch** stay in `/tmp/chk050/` (not committed). `/tmp/chk050/lib` and `/tmp/chk050/src` are **kept**: Gate 3 failed as written (condition b), and a failed gate stops the task before Step 5.4. `~/fs050_check/payloads/` is kept.
