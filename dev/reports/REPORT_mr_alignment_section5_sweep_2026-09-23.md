# REPORT — Section 5 sweep: the MR alignment's shift across all 18 cells

- **Task:** `dev/tasks/TASK_mr_alignment_section5_sweep_2026-09-23.md` (`60bf95b1`); extends the probe
  (`REPORT_mr_alignment_section5_probe_2026-09-23.md`).
- **Result: complete, 18 of 18 cells, not partial. Refresh, not revision.** Declaration is identical to the
  committed rows on every replicate of every cell. The largest single shift anywhere in the design is 0.0152, on
  a Bonferroni lower bound sitting 0.28 below the nearer threshold. The only threshold readings that change are 7
  replicates whose before value was within 0.0025 of 0.75 or 1.25.

## 0. Record

| item | value |
|---|---|
| HEAD at start | `60bf95b1` (task doc); `R/` last changed at `7713942e` |
| R, platform | R 4.6.1, pop-os, Linux 7.1.5, 128 physical cores |
| Workers | 64 per render |
| Library the workers load | `/tmp/claude-1000/-home-larryleon-Documents-GitHub-forestsearch/ab0a4121-64a0-4c61-972b-998577a64b16/scratchpad/Rlib_head/forestsearch` (HEAD build, first on `R_LIBS`); installed build (stale, 0.3.5.9000 of 02:32 UTC) not used, not reinstalled |
| Untracked before the task (never staged) | the two `actg175/.../_d5000/` directories, `actg175/binary_020/smoke_redes.html`, `smoke_relaunch.html`, `gbsg_020/scripts_dinamr/logs/nullmr_findings.err` |
| Sweep window | 17:05:15 → 18:43:05 PDT, **98 min** of the 150 min budget |

## 1. The 18 cells, from the committed record

The committed `p12x20` record holds **nine** cells (12.4%). The other nine (31%) are the `cert20` campaign's
seven plus the two n = 500 harm cells that `cert20` itself took from `e1stud` rather than re-running
(`TASK_cert20_2026-09-08.md` line 39). All 18 share `effMaxSG`, ε 0.20, J 10, field complement on, field-s
scale `selected`, `two_term` IJ residual, FB none, seed `8316951 + sim_id`, 2,000 replicates; the metas differ
only in n, planted HR, prevalence, and two meta fields (`mr_inference`, `field_recovery`) that the older campaigns
predate.

| cell | source (committed bundle) | prevalence | HR | n | committed declaration | committed wall (workers) |
|---|---|---|---|---|---|---|
| A7 | `p12x20` (`p12x20_2026-09-12/`) | 12.4% | 1.00 | 500 | 0.6805 | 1309 s (64) |
| A8 | `p12x20` | 12.4% | 1.00 | 1000 | 0.6595 | 1763 s (64) |
| A9 | `p12x20` | 12.4% | 1.00 | 1500 | 0.6240 | 2396 s (64) |
| A1 / A2 / A3 | `p12x20` | 12.4% | 1.50 | 500 / 1000 / 1500 | 0.9110 / 0.9740 / 0.9880 | 1575 / 2191 / 3352 s (64) |
| A4 / A5 / A6 | `p12x20` | 12.4% | 1.75 | 500 / 1000 / 1500 | 0.9500 / 0.9950 / 0.9990 | 1618 / 2213 / 3371 s (64) |
| P31_h100_n500 / n1000 / n1500 | `cert20` (`results/`) | 31% | 1.00 | 500 / 1000 / 1500 | 0.9205 / 0.9545 / 0.9590 | 1616 / 2427 / 3743 s (100) |
| P31_h150_n500 | `e1stud` (`results/`) | 31% | 1.50 | 500 | 0.9995 | 1858 s (100) |
| P31_h175_n500 | `e1stud` | 31% | 1.75 | 500 | 0.9995 | 1901 s (100) |
| P31_h150_n1000 / n1500 | `cert20` | 31% | 1.50 | 1000 / 1500 | 1.0000 / 1.0000 | 2668 / 4118 s (100) |
| P31_h175_n1000 / n1500 | `cert20` | 31% | 1.75 | 1000 / 1500 | 1.0000 / 1.0000 | 2694 / 4147 s (100) |

Full stems: `mrs5sweep/cells.txt`. Declaration rates computed from the committed bundles (`detected == 1L`);
walls from `STATUS_p12x20.md`, `REPORT_cert20_2026-09-08.md`, `REPORT_field_studentize_e1_2026-09-08.md`.

**Before-arm validity for the 31% family.** The route's premise was established on `p12x20` only (A7, 163/163).
`cert20` and `e1stud` ran on older builds (0.3.5, 2026-09-08), so the same parent-build check was repeated
before the sweep: 10 replicates on a `git archive` export of `ba595f4b` installed to a separate library
(`.../Rlib_pre`; worker assertion: fix absent in all 64 workers).

| committed bundle | columns identical | max abs diff |
|---|---|---|
| `cert20` HR 1.00 n 500 (`..._h100_knoise0_n500_z1q60_nb20_cert20_combined_1_2000.rds`) | 153 / 153 | 0 |
| `e1stud` HR 1.50 n 500 (`..._h150_knoise0_n500_z1q60_nb20_e1stud_combined_1_2000.rds`) | 153 / 153 | 0 |

So all three source campaigns are valid before arms, and every shift below is the alignment alone.

## 2. Run

Order: n 500, then 1000, then 1500; within each, HR 1.00 before harm; 12.4% before 31%. Projection =
committed wall × 0.1 × (committed workers / 64) × 1.44 (A7's measured/projected ratio). Driver
`mrs5sweep/sweep.sh`, per-cell `mrs5sweep/runcell.sh` (the unchanged `render.sh` + template with the campaign
knob set, `FS_S7_Z1Q=0.60` for the 31% cells, `FS_S7_CAMPAIGN=mrs5sweep`, `FS_S7_START=1 FS_S7_NSIMS=200`).

| cell | projection (s) | measured (s) |
|---|---|---|
| A7 | — | 188 (probe, reused) |
| P31_h100_n500 | 364 | 228 |
| A1 | 227 | 219 |
| A4 | 233 | 221 |
| P31_h150_n500 | 418 | 258 |
| P31_h175_n500 | 428 | 262 |
| A8 | 254 | 239 |
| P31_h100_n1000 | 546 | 313 |
| A2 | 316 | 297 |
| A5 | 319 | 301 |
| P31_h150_n1000 | 600 | 340 |
| P31_h175_n1000 | 606 | 353 |
| A9 | 345 | 324 |
| P31_h100_n1500 | 842 | 434 |
| A3 | 483 | 415 |
| A6 | 485 | 425 |
| P31_h150_n1500 | 927 | 482 |
| P31_h175_n1500 | 933 | 488 |

No budget stop, no overrun. The 31% projections ran high because the 100-worker committed walls carry more
contention than 64 workers do.

## 3–4. Read-out (200 replicates per cell, sim_id 1–200, paired on seed; shift = after − before, HR scale, over declaring replicates)

| cell | prev | HR | n | declaring | decl. rate before / after | est shift mean / med / p05 / p95 | field lower shift mean / med / p05 / p95 | field-s upper shift mean / med / p05 / p95 | largest single shift (sim, qty) | region chg | MR top-label chg | bound flips at 0.75/1.25 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| A7 | 12.4% | 1.00 | 500 | 126/200 | 0.6300 / 0.6300 | +0.00095 / +0.00069 / -0.00024 / +0.00323 | +0.00052 / +0.00041 / -0.00023 / +0.00184 | +0.00047 / +0.00017 / -0.00080 / +0.00273 | +0.00630 (sim 149, bonf_up) | 0 | 1 | 1 |
| P31_h100_n500 | 31% | 1.00 | 500 | 182/200 | 0.9100 / 0.9100 | +0.00032 / +0.00014 / -0.00017 / +0.00146 | +0.00019 / +0.00010 / -0.00010 / +0.00099 | +0.00011 / +0.00002 / -0.00144 / +0.00164 | -0.01522 (sim 18, bonf_lo) | 0 | 0 | 0 |
| A1 | 12.4% | 1.50 | 500 | 178/200 | 0.8900 / 0.8900 | +0.00039 / +0.00017 / -0.00038 / +0.00175 | +0.00020 / +0.00009 / -0.00028 / +0.00122 | +0.00028 / +0.00007 / -0.00055 / +0.00187 | -0.00977 (sim 87, bonf_lo) | 0 | 0 | 0 |
| A4 | 12.4% | 1.75 | 500 | 188/200 | 0.9400 / 0.9400 | +0.00017 / +0.00000 / -0.00049 / +0.00123 | +0.00014 / +0.00000 / -0.00040 / +0.00098 | +0.00020 / +0.00002 / -0.00071 / +0.00162 | +0.00934 (sim 102, bonf_up) | 0 | 0 | 0 |
| P31_h150_n500 | 31% | 1.50 | 500 | 200/200 | 1.0000 / 1.0000 | +0.00003 / +0.00000 / -0.00014 / +0.00032 | +0.00002 / +0.00000 / -0.00007 / +0.00018 | -0.00002 / +0.00000 / -0.00015 / +0.00032 | -0.00574 (sim 197, up) | 0 | 0 | 0 |
| P31_h175_n500 | 31% | 1.75 | 500 | 200/200 | 1.0000 / 1.0000 | +0.00000 / +0.00000 / -0.00007 / +0.00007 | +0.00000 / +0.00000 / -0.00004 / +0.00003 | +0.00001 / +0.00000 / -0.00002 / +0.00002 | +0.00296 (sim 61, up) | 0 | 0 | 0 |
| A8 | 12.4% | 1.00 | 1000 | 128/200 | 0.6400 / 0.6400 | +0.00076 / +0.00060 / -0.00020 / +0.00235 | +0.00073 / +0.00055 / -0.00016 / +0.00241 | +0.00036 / +0.00015 / -0.00033 / +0.00193 | +0.00832 (sim 197, lo) | 0 | 0 | 1 |
| P31_h100_n1000 | 31% | 1.00 | 1000 | 190/200 | 0.9500 / 0.9500 | +0.00022 / +0.00006 / -0.00011 / +0.00113 | +0.00018 / +0.00004 / -0.00007 / +0.00092 | +0.00005 / +0.00002 / -0.00069 / +0.00083 | +0.00591 (sim 157, bonf_up) | 0 | 0 | 0 |
| A2 | 12.4% | 1.50 | 1000 | 194/200 | 0.9700 / 0.9700 | +0.00010 / +0.00000 / -0.00037 / +0.00095 | +0.00007 / +0.00000 / -0.00031 / +0.00091 | +0.00006 / +0.00001 / -0.00075 / +0.00105 | -0.00891 (sim 199, bonf_lo) | 0 | 0 | 0 |
| A5 | 12.4% | 1.75 | 1000 | 199/200 | 0.9950 / 0.9950 | +0.00002 / +0.00000 / -0.00037 / +0.00054 | -0.00003 / +0.00000 / -0.00031 / +0.00040 | +0.00008 / +0.00000 / -0.00015 / +0.00069 | -0.01114 (sim 199, bonf_lo) | 0 | 0 | 0 |
| P31_h150_n1000 | 31% | 1.50 | 1000 | 200/200 | 1.0000 / 1.0000 | +0.00000 / +0.00000 / +0.00000 / +0.00000 | +0.00000 / +0.00000 / +0.00000 / +0.00000 | +0.00001 / +0.00000 / +0.00000 / +0.00000 | -0.00127 (sim 199, bonf_up) | 0 | 0 | 0 |
| P31_h175_n1000 | 31% | 1.75 | 1000 | 200/200 | 1.0000 / 1.0000 | -0.00000 / +0.00000 / +0.00000 / +0.00000 | -0.00000 / +0.00000 / +0.00000 / +0.00000 | -0.00000 / +0.00000 / +0.00000 / +0.00000 | -0.00004 (sim 176, est) | 0 | 0 | 0 |
| A9 | 12.4% | 1.00 | 1500 | 134/200 | 0.6700 / 0.6700 | +0.00080 / +0.00070 / -0.00012 / +0.00206 | +0.00080 / +0.00060 / -0.00019 / +0.00255 | +0.00030 / +0.00016 / -0.00044 / +0.00149 | +0.00590 (sim 95, up) | 0 | 2 | 1 |
| P31_h100_n1500 | 31% | 1.00 | 1500 | 194/200 | 0.9700 / 0.9700 | +0.00017 / +0.00005 / -0.00007 / +0.00077 | +0.00014 / +0.00004 / -0.00005 / +0.00064 | +0.00010 / +0.00000 / -0.00026 / +0.00097 | -0.00444 (sim 26, bonf_up) | 0 | 1 | 0 |
| A3 | 12.4% | 1.50 | 1500 | 198/200 | 0.9900 / 0.9900 | +0.00002 / +0.00000 / -0.00042 / +0.00056 | +0.00002 / +0.00000 / -0.00039 / +0.00053 | +0.00008 / +0.00000 / -0.00027 / +0.00076 | +0.00432 (sim 129, bonf_lo) | 0 | 3 | 0 |
| A6 | 12.4% | 1.75 | 1500 | 200/200 | 1.0000 / 1.0000 | -0.00003 / +0.00000 / -0.00026 / +0.00013 | -0.00002 / +0.00000 / -0.00023 / +0.00013 | +0.00003 / +0.00000 / -0.00013 / +0.00027 | +0.00456 (sim 187, bonf_lo) | 0 | 0 | 0 |
| P31_h150_n1500 | 31% | 1.50 | 1500 | 200/200 | 1.0000 / 1.0000 | -0.00000 / +0.00000 / +0.00000 / +0.00000 | -0.00000 / +0.00000 / +0.00000 / +0.00000 | -0.00000 / +0.00000 / +0.00000 / +0.00000 | -0.00018 (sim 25, up) | 0 | 0 | 0 |
| P31_h175_n1500 | 31% | 1.75 | 1500 | 200/200 | 1.0000 / 1.0000 | +0.00000 / +0.00000 / +0.00000 / +0.00000 | +0.00000 / +0.00000 / +0.00000 / +0.00000 | -0.00000 / +0.00000 / +0.00000 / +0.00000 | +0.00001 (sim 128, bonf_up) | 0 | 0 | 0 |

"Largest single shift" is over all five quantities (the three above plus the Bonferroni pair
`fld_joint_s_bonf_loH` / `fld_joint_s_bonf_upHc`). The Bonferroni pair's own mean/percentiles are in
`mrs5sweep/readout.rds`. "Region chg" is the search's selection (`label` or `n_sel`); "MR top-label chg" is the
first entry of `p_hat_top_labels` (MR's re-selection winner). The declaration rate in the committed 2,000 differs
from the 200 here only because these are the first 200 draws.

**Every threshold crossing in the design** (est = corrected estimate; the reading changes only at the knife edge):

| cell | quantity | threshold | sim | before | after |
|---|---|---|---|---|---|
| A7 | est | 1.25 | 119 | 1.24961 | 1.25049 |
| A7 | Bonferroni upper | 0.75 | 66 | 0.74996 | 0.75507 |
| P31_h100_n500 | est | 0.75 | 105 | 0.74797 | 0.75002 |
| A1 | est | 1.25 | 71 | 1.24958 | 1.25096 |
| A8 | field-s upper | 0.75 | 156 | 0.74800 | 0.75098 |
| A9 | est | 0.75 | 139 | 0.74943 | 0.75126 |
| A9 | field-s upper | 0.75 | 168 | 0.74984 | 0.75131 |

The design-wide largest shift: P31_h100_n500 sim 18 (region q18.0 & q27.0, N 69), Bonferroni lower
0.47132 → 0.45611; its other quantities moved by at most 0.0008.

**Reading.**

- **Systematic pattern: the shift scales with how marginal the admitted set is.** It is largest in the
  attenuated-benefit cells (HR 1.00: mean corrected-estimate shift +0.0008 to +0.0010 at 12.4%, +0.0002 to
  +0.0003 at 31%). It shrinks under harm and with prevalence: exactly zero at the 95th percentile in the 31%
  harm cells at n ≥ 1000, where true harm candidates clear the threshold by a wide margin. n has little effect
  at HR 1.00 (A7/A8/A9: +0.00095 / +0.00076 / +0.00080), and a shrinking effect under harm. Direction: upward
  (towards larger HR) on all three quantities in all six HR 1.00 cells and in every cell whose mean shift exceeds
  0.0001; the harm cells' means below 0.0001 in magnitude take either sign. A plausible mechanism, not tested
  here: the aligned threshold is the rounded one, slightly lower than the exact cutoff, so it admits marginal
  candidates the exact cutoff excluded.
- **Small relative to the thresholds.** The 95th-percentile shift never exceeds 0.0033. The median distance from
  a bound to the nearer of 0.75 / 1.25 is at least 0.039 in every cell. The largest single shift anywhere (0.015)
  is on a bound 0.28 from the nearer threshold. The 7 readings that change are all within 0.0025 of the line
  before the change: 4 bound readings and 3 corrected-estimate readings, out of roughly 3,300 declaring
  replicates × 5 quantities.
- **Re-running Section 5 at full replicates is a refresh, not a revision.** The published figures would move
  in the third or fourth decimal. No declaration, region or N changes anywhere.

## 5. Gates

- **Gate 1 — search untouched:** `detected` identical to the committed rows on 200 / 200 replicates in all 18
  cells; `label` and `n_sel` identical on every declaring replicate (0 region changes). PASS.
- **Gate 2 — right build:** before each of the 17 new cells, 64 multisession futures under the render's own
  environment reported one library path (`.../Rlib_head/forestsearch`), `pconsistency.digits` in
  `formals(fs_mr_inference)` and `.fs_pcons_eff` in its body: 17 / 17 `OK=TRUE`
  (`mrs5sweep/logs/*.buildcheck`). A7's check is in `mrs5probe/logs/`. PASS.
- **Gate 3 — no `R/` file modified:** `git status --short -- R/` empty. PASS.
- **Gate 4 — errors:** 0 `err_msg`, 0 MR failures among declaring replicates, 0 non-finite shifts, across
  3,600 replicates. PASS.

## 6. Files (under `quarto/simulations/gbsg_020/`)

| path | count / size |
|---|---|
| `results/fs_effMaxSG_fb_mr_field_m1_h{100,150,175}_knoise0_n{500,1000,1500}[_z1q60]_nb20_mrs5sweep_res_1_200.rds` | 17 after-arm bundles, 118–166 KB each (A7's is the probe's `..._mrs5probe_res_1_200.rds`) |
| `results/fs_effMaxSG_fb_mr_field_m1_h{100,150}_knoise0_n500_z1q60_nb20_mrs5pre_res_1_10.rds` | 2 parent-build check bundles (`cert20`, `e1stud`) |
| `mrs5sweep_<cell>_batch_1.html` | 17 renders, 73 MB total, each under 5 MB |
| `mrs5pre_C31_h100_n500_gateP_batch_1.html`, `mrs5pre_E31_h150_n500_gateP_batch_1.html` | parent-build check renders |
| `mrs5sweep/logs/` | driver log, 17 render logs (`WALL_SECONDS`), 17 worker build checks |
| `mrs5pre/logs/` | the two parent-build check logs and build checks |
| `mrs5sweep/{sweep.sh,runcell.sh,pair.R,readout.R,cells.txt}` | driver, runner, paired comparison, read-out, cell list |
| `mrs5sweep/pair/<cell>.rds`, `mrs5sweep/readout.{rds,md}` | per-cell paired statistics and the table above |

Nothing over 50 MB. Nothing written outside the simulation directory, `dev/tasks/` and `dev/reports/`.
