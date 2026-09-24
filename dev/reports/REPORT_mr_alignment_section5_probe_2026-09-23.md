# REPORT — Section 5 probe: does the MR alignment move the simulation results?

- **Task:** `dev/tasks/TASK_mr_alignment_section5_probe_2026-09-23.md` (`ebfc4548`).
- **Result: refresh, not revision.** In cell A7, every declaring replicate shifts, by a mean of +0.0010 on the
  corrected estimate and +0.0005 on each bound; the largest single shift is 0.0063. Declaration, region and N
  are identical in both arms on all 200 replicates.

## 0. Record

| item | value |
|---|---|
| HEAD at Gate T start / at probe run | `57807ec2` / `ebfc4548` (task doc only; `R/` last changed at `7713942e`) |
| R, platform | R 4.6.1, pop-os, Linux 7.1.5, 128 physical cores |
| Workers | 64 (`FS_S7_WORKERS=64`, the `p12x20` production count) |
| Installed forestsearch | 0.3.5.9000, built 2026-09-23 02:32 UTC, **stale**: no `pconsistency.digits` in `fs_mr_inference`. Not used, not reinstalled. |
| Library the workers load (after arm) | `/tmp/claude-1000/-home-larryleon-Documents-GitHub-forestsearch/ab0a4121-64a0-4c61-972b-998577a64b16/scratchpad/Rlib_head/forestsearch` (HEAD working tree, `R CMD INSTALL`, set first with `R_LIBS`) |
| Worker assertion (after arm, 200-rep run) | Before the render, 64 multisession futures under the same env: 64 distinct PIDs, one library path (above), `pconsistency.digits` in `formals(fs_mr_inference)` TRUE in all, `.fs_pcons_eff` in its body TRUE in all. `mrs5probe/logs/mrs5probe_A7_batch_1.buildcheck`. |
| Untracked before the task (never staged) | `actg175/binary_020/mr_or_harm/..._redes_d5000/`, `..._relaunch_d5000/`, `actg175/binary_020/smoke_redes.html`, `smoke_relaunch.html`, `gbsg_020/scripts_dinamr/logs/nullmr_findings.err` |

## 1. Driver, payload, route

- **Driver:** `quarto/simulations/gbsg_020/scripts_p12x20/campaign_p12x20.sh` (`118a66ce`) → `scripts_dinamr/render.sh` →
  `sim_fs_maxeffCons_fb_mr_field_m1_template.qmd`. The template and `render.sh` are reused unchanged. `cell()`
  hard-codes 1,000-replicate batches, so `render.sh` is called directly with `cell()`'s knob set
  (`FS_S7_FOCUS=effMaxSG FS_S7_NBHD=0.20 FS_S7_FIELD_COMPLEMENT=TRUE FS_S7_FIELD_SCALEC=selected
  FS_S7_FIELD_DECOMP=TRUE FS_S7_IJ_RESIDUAL=two_term FS_S7_FB=none FS_S7_RETURN_RESEL=TRUE FS_S7_WORKERS=64`,
  every other `FS_S7_*` unset), plus `FS_S7_CAMPAIGN=mrs5probe`, `FS_S7_MODE=batch FS_S7_START=1
  FS_S7_NSIMS=200`. Wrapper: `gbsg_020/mrs5probe/runcell.sh`.
- **Payload (before arm):** `gbsg_020/p12x20_2026-09-12/fs_effMaxSG_fb_mr_field_m1_h100_knoise0_n500_nb20_p12x20_combined_1_2000.rds`,
  committed in `bc5e6a69`. Cell **A7** (HR 1.00, n 500, prevalence 12.4%), 2,000 replicates, seed `8316951 + sim_id`.
- **Cell note.** The task named n = 500 and cited 0.624 as the design's lowest declaration rate. In the
  committed record 0.624 is **A9** (n 1500); A7 declares at 0.6805. The probe ran A7 as named (Larry's go).
- **Route: after arm only.** The payload retains every per-replicate quantity the task lists: `detected`,
  `label`/`n_sel`, `mr_H_est`, `fld_H_lo1s`, `fld_Hc_up1s_s`, `fld_joint_s_bonf_loH`/`fld_joint_s_bonf_upHc`.

## 2. The committed payload is a valid before arm (parent-build check)

The Gate T after-arm rows showed ~1e-3 MR shifts. To rule out drift from other commits since the `p12x20` base
(`0ab5d1c5`), the same 10 replicates were run on a `git archive` export of the fix's parent **`ba595f4b`**,
installed to a separate scratch library (`.../Rlib_pre`). That tree has the exact cutoff at
`R/fs_mr_inference.R:660` and `pconsistency.digits` on `forestsearch()`; `git diff ba595f4b 7713942e -- R/`
touches only the MR admission (`R/fs_mr_inference.R`, and one line in `R/forestsearch_main.R` forwarding
`pconsistency.digits` to MR).

**Result: 163 / 163 non-timing columns identical to the committed rows, maximum absolute difference 0**
(timing columns excluded: `fb_secs fit_mr_secs fld_H_secs fld_H_uniform_secs fld_Hc_secs`). Re-verified from
the files on disk for this report. So the committed payload *is* the before arm, and any shift is the
alignment alone. Payload: `results/fs_effMaxSG_..._n500_nb20_mrs5pre_res_1_10.rds`; render
`mrs5pre_A7_gateT_batch_1.html`; log `mrs5probe/logs/mrs5pre_A7_gateT_batch_1.log` (61 s).

## 3. Gate T

10 replicates, 64 workers: **61 s** wall, mean 18.5 s per replicate (declaring 27–31 s, non-declaring ~2.7 s),
projected 2–3 min for 200. **Measured for 200: 188 s** (`WALL_SECONDS=188`). The 200-replicate run reproduces
its own Gate T rows 1–10 exactly (173 / 173 columns).

## 4–5. Probe read-out (A7, 200 replicates, sim_id 1–200, paired on seed)

Declaring replicates: 126 / 200 in both arms. Shift = after − before, HR scale.

| quantity (column) | before mean | after mean | shift mean | median | p05 | p95 | largest (sim_id; region, N) | flips at 0.75 / 1.25 |
|---|---|---|---|---|---|---|---|---|
| corrected estimate (`mr_H_est`) | 0.8504 | 0.8514 | +0.00095 | +0.00069 | −0.00024 | +0.00323 | +0.00548 (188; q8.1 & q22.1, 87) | 0 / 1 |
| field lower bound (`fld_H_lo1s`) | 0.4312 | 0.4317 | +0.00052 | +0.00041 | −0.00023 | +0.00184 | −0.00498 (87; q4.1 & q25.1, 67) | 0 / 0 |
| field-s upper bound (`fld_Hc_up1s_s`) | 0.8302 | 0.8307 | +0.00047 | +0.00017 | −0.00080 | +0.00274 | +0.00484 (66; q15.1 & q26.0, 61) | 0 / 0 |
| Bonferroni lower (`fld_joint_s_bonf_loH`) | 0.3733 | 0.3738 | +0.00050 | +0.00037 | −0.00015 | +0.00185 | −0.00605 (87; q4.1 & q25.1, 67) | 0 / 0 |
| Bonferroni upper (`fld_joint_s_bonf_upHc`) | 0.8685 | 0.8689 | +0.00042 | +0.00007 | −0.00077 | +0.00310 | +0.00630 (149; q4.1 & q11.1, 116) | 1 / 0 |

| replicate-level count | value |
|---|---|
| declaration rate, before / after | 0.630 / 0.630 (identical on every replicate) |
| region changed (rule or N), of 126 declaring | 0 |
| MR re-selection top label changed, of 126 declaring | 1 (sim 183) |
| declaring replicates with any nonzero shift | 126 / 126 |
| errors / MR failures | 0 / 0 of 200 |
| median distance of the before value to the nearer of 0.75 / 1.25 | est 0.071, lo 0.349, up 0.084, Bonf lo 0.397, Bonf up 0.113 |

The two threshold flips are knife-edge cases: sim 119, corrected estimate 1.24961 → 1.25049; sim 66,
Bonferroni upper 0.74996 → 0.75507. Both sat within 0.0005 of the threshold before the change.

**Reading.** The shift is two orders of magnitude smaller than the typical distance between a bound and the
threshold it is read against: the 95th percentile shift is at most 0.0032, against median distances of
0.07–0.40. It moves a threshold reading only where the before value was already within 0.0005 of the line
(2 of 126 × 5 readings). The only columns that change are the MR/field outputs (`mr_*`, `fld_*`, `p_hat_*`);
every search-side column is byte-identical. **Re-running Section 5 is a refresh, not a revision.**

## 6. Gates

- **Gate T:** reported, and the task stopped for the go (session `ab0a4121`). PASS.
- **Gate 1:** the route reuses the committed payload, so the task needs no check. The parent-build check (§2)
  confirms it anyway: 163 / 163, max diff 0. PASS.
- **Gate 2:** declaration, `label` and `n_sel` identical on 200 / 200 replicates; 0 region changes. PASS.
- **Gate 3:** `git status --short -- R/` empty at the end. PASS.

## 7. Files

| path (under `quarto/simulations/gbsg_020/`) | size |
|---|---|
| `results/fs_effMaxSG_fb_mr_field_m1_h100_knoise0_n500_nb20_mrs5probe_res_1_200.rds` (probe, after arm) | 117 KB |
| `results/fs_effMaxSG_fb_mr_field_m1_h100_knoise0_n500_nb20_mrs5probe_res_1_10.rds` (Gate T) | 8 KB |
| `results/fs_effMaxSG_fb_mr_field_m1_h100_knoise0_n500_nb20_mrs5pre_res_1_10.rds` (parent-build check) | 8 KB |
| `mrs5probe_A7_batch_1.html`, `mrs5probe_A7_gateT_batch_1.html`, `mrs5pre_A7_gateT_batch_1.html` | 4.3–4.5 MB each |
| `mrs5probe/logs/*` (render logs with `WALL_SECONDS`, worker build check) | < 4 KB each |
| `mrs5probe/runcell.sh`, `mrs5probe/pair.R` (runner and paired read-out) | < 3 KB each |
| `mrs5probe/pair/A7_probe.rds` (the read-out table) | 2 KB |

Nothing over 50 MB. No `R/` file modified; nothing written outside the simulation directory, `dev/tasks/`
and `dev/reports/` (the scratch libraries are outside the repo and are not deliverables).
