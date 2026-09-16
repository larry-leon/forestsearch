# REPORT — ACTG175 continuous (MD) re-run under the current field constructions: Stage 1 (edits, smoke, calibration → Gate 1)

Date: 2026-09-15/16 (UTC). Machine: `pop-os` (AMD Ryzen Threadripper PRO 5995WX, 64 physical / 128 logical cores, 251 GB; reference BLAS/LAPACK 3.12.0; R 4.6.1). Branch `feature/glm-extension`. Task: `dev/tasks/TASK_md_field_rerun_2026-09-15.md` (committed as received, `bb84120e`). Stage 0 record: `REPORT_md_field_rerun_stage0_2026-09-15.md` (`b6e30ac7`). Every render in this stage ran with `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1`. Stages 2–3 were **not** started.

## 1.1 Provenance and first commit — GATE PASS

```
pop-os
feature/glm-extension
0071c17e
0071c17e Merge branch 'feature/glm-extension' of github.com:larry-leon/forestsearch into feature/glm-extension
b6e30ac7 ACTG175 continuous field re-run Stage 0 (read-only): record -- ...
70975bc4 Add TASK_md_field_rerun_stage0_2026-09-15 as received
[tracked modifications: none]   [pending under R/ DESCRIPTION NAMESPACE: none]
[R / Rscript / quarto processes: none]
Stage 0 record in HEAD
```
- HEAD at kickoff was `0071c17e`, one merge ahead of the Stage 0 record: `git diff --stat b6e30ac7..0071c17e` is a single binary file (`FEASIBILITY_extreme_subgroups_HTA.docx`, `d9c0b4ba`); no change under `R/`, `DESCRIPTION`, `NAMESPACE` or to the MD template.
- First commit: `bb84120e Add TASK_md_field_rerun_2026-09-15 as received` (the task document alone; `~/Downloads/TASK_md_field_rerun_2026-09-15.md`, the exact name).

## 1.2 Install — GATE PASS (with one finding)

- `Rscript -e 'devtools::install(quick = TRUE, upgrade = "never")'` **failed** before installing: `cli::cli_abort("{.arg upgrade} must be a single TRUE, FALSE, or NA")` — this devtools/remotes does not accept `upgrade = "never"`. **Finding:** the task's literal command does not run here; the equivalent `devtools::install(quick = TRUE, upgrade = FALSE)` was run instead (rc 0, `* DONE (forestsearch)`), a rebuild of unchanged `R/` at HEAD `bb84120e`.
- `packageDescription("forestsearch")$Built` = **`R 4.6.1; ; 2026-09-16 05:57:14 UTC; unix`** (was `2026-09-13 06:52:22 UTC` before the rebuild).
- Two `doFuture` multisession workers each returned the same `Built` (pids 1224337 / 1224338): **agree = TRUE**.
- Installed `fs_mr_inference()` formals: `ci_method = c("field", "ij", "wald")`, `field_scale_complement = c("selected", "none")` — as Stage 0 §3 recorded.

## 1.3 The survival reference (read only)

Sources: **MD** = `sim_fs_maxeffCons_mr_field_md_template.qmd` at `0071c17e` (line numbers before the §1.4 edits); **m1** = `quarto/simulations/gbsg_020/sim_fs_maxeffCons_fb_mr_field_m1_template.qmd` at HEAD; **runner** = `quarto/simulations/gbsg_020/scripts_p12x20/campaign_p12x20.sh:24–27` (`KN=(FS_S7_FOCUS=effMaxSG FS_S7_NBHD=0.20 FS_S7_FIELD_COMPLEMENT=TRUE FS_S7_FIELD_SCALEC=selected FS_S7_FIELD_DECOMP=TRUE FS_S7_IJ_RESIDUAL=two_term FS_S7_FB=none FS_S7_CAMPAIGN=p12x20 FS_S7_RETURN_RESEL=TRUE FS_S7_WORKERS=$P12X20_WORKERS)`) and `run_p12x20.sh:14` (`WORKERS=64`); **meta** = `p12x20_2026-09-12/fs_effMaxSG_fb_mr_field_m1_h100_knoise0_n500_nb20_p12x20_{combined_1_2000,res_1_1000}.rds` and `results/fs_effMaxSG_fb_mr_field_m1_h100_knoise0_n1000_z1q60_nb20_cert20_combined_1_2000.rds` (computed here).

| Setting | MD template | FS survival `effMaxSG` campaigns (m1 / runner / meta) | Campaign takes |
|---|---|---|---|
| `sg_focus` | `"maxeffCons"` literal (`:135`) | `.env_chr("FS_S7_FOCUS", "maxeffCons")` (m1 `:318`), guard `stopifnot(sg_focus %in% c("effMaxSG","maxeffCons","effMinSG","maxSG","minSG","maxeff"))` (`:333`); runner `FS_S7_FOCUS=effMaxSG`; meta `sg_focus effMaxSG` (p12x20, cert20) | **survival's: `effMaxSG`** |
| `effect_neighborhood` | `0.10` literal (`:214`) | `.env_num("FS_S7_NBHD", 0.10)` (`:342`) with the band-foci guard (`:350–353`) and the `_nb%02d` stem tag (`:443–444`); runner `FS_S7_NBHD=0.20`; meta `0.2` | **survival's: 0.20** |
| `selection_rule` | `"neighborhood"` literal (`:213`) | `"neighborhood"` literal (`:547`) with the package-mirroring guard (`:555–558`); meta (batch) `neighborhood` | **survival's: `"neighborhood"`** (same value; no `FS_MD_RULE` knob needed) |
| `stop_threshold` | `NULL` (`:215`) | `NULL` (`:559`); meta `"NULL"` | MD's (same) |
| `consistency_method` | `"resample"` (`:210`) | `"resample"` (`:542`); meta `resample` | MD's (same) |
| `pconsistency.threshold` | `0.90` (`:220`) | `0.90` (`:570`) | MD's (same) |
| `fs.splits` | `400L` (`:221`) | `400L` (`:571`) | MD's (same) |
| `maxk` | `2L` (`:221`) | `2L` (`:572`) | MD's (same) |
| `n.min` | `60L` (`:221`) | `NULL` (`:573`) | MD's — **difference (finding F1)** |
| arm minima `d0.min` / `d1.min` | `12L` / `12L` (`:221`) | `10L` / `10L` (`:574–575`) | MD's — **difference (F1)** |
| effect / consistency thresholds | `effect.threshold = 30`, `consistency.threshold = 10` (MD scale; `:216–218`) | `hr_threshold = 0.90`, `hr_consistency = 0.80` (HR scale; `:568–569`) | MD's — design-specific (F1) |
| seed scheme | `seed_base 8316951L` (`:131`); `sd_i <- seed_base + sim_id` with `RNGkind("L'Ecuyer-CMRG")` (`:565–568`); `seedit = sd_i` (`:606`) | `seed_base 8316951L` (`:707`); `seedit = seed_base + sim_id` (`:1044`); meta `seed_base 8316951` | MD's (same base and offset) |
| MR `draws` | `mr_draws 5000L` (`:122`) | `5000L` (`:488`); meta `5000` | MD's (same) |
| `multiplier` | not passed → `fs_mr_inference()` default `"poisson"` | not passed → `"poisson"` | same |
| `ci_method` | `.env_chr("FS_MD_CI", "field")` (`:263`) | `"field"` literal (`:615`); meta `field` | same |
| `ij_residual` | `.env_chr("FS_MD_IJ_RESIDUAL", "two_term")` (`:274`) | `.env_chr("FS_S7_IJ_RESIDUAL", "two_term")` (`:639`); runner `two_term` | same |
| `confirm_rule` | `"point"` (`:266`) | `"point"` (`:617`) | same |
| `t_confirm` | `NULL` (`:265`) → near-null default (0 for MD) | `NULL` (`:616`) → 1 for HR | same rule, measure-specific value |
| `field_complement` | `FS_MD_FIELD_COMPLEMENT` default TRUE (`:271`) | `FS_S7_FIELD_COMPLEMENT` default TRUE (`:632`); runner TRUE; meta TRUE | same |
| `field_scale_complement` | **not passed** (inherited `"selected"`) | `.env_chr("FS_S7_FIELD_SCALEC", "selected")` passed (`:658`, `:688`); runner `selected`; meta `selected` | E2 adds `FS_MD_FIELD_SCALEC` (default `"selected"`) |
| `field_decompose` | not passed (FALSE) | `FS_S7_FIELD_DECOMP` (`:648`); runner TRUE; meta TRUE | MD's (FALSE) — **difference (F2)**: the survival campaigns recorded the scale diagnostics (`rho^c`); this campaign does not |
| `return_reselection` | `FS_MD_RESELECTION` default TRUE (`:278`) | `FS_S7_RETURN_RESEL` default TRUE (`:682`); runner TRUE | same |
| field `R_out` / `R_in` | not forwarded; package defaults 1000 / 500 (meta `field_R`) | same (m1 `:612–614`; batch meta `field_R`) | same |
| workers | `FS_MD_WORKERS`, cap physical − 1 (`:116`) | `FS_S7_WORKERS` default 60, cap physical − 1 (`:711`); runner 64; batch meta `n_workers 64` | from Gate 1 |

Rule (task §1.3): the campaign takes `sg_focus = "effMaxSG"`, `effect_neighborhood = 0.20`, `selection_rule = "neighborhood"`; everything else stays the MD template's. Differences from the survival campaigns, as findings for Gate 1: **F1** the MD design's own identification floors (`n.min 60`, `d0.min/d1.min 12`, thresholds 30/10 on the MD scale, J = 10 cut grids on `age`/`preanti`) vs survival's (`n.min NULL`, 10/10, 0.90/0.80, J = 10 on `er`); **F2** `field_decompose` FALSE here vs TRUE in p12x20/cert20 (diagnostic-only; excluded from every table by convention 7 anyway). Both are by the task's rule, not changes.

## 1.4 Template edits — commit `2cffb95f` (template only; `R/` untouched)

Applied bottom-up by line anchor at `0071c17e`'s numbering (Stage 0's), so every Stage 0 citation held while editing; the file went from 1,803 to 1,888 lines (`git diff --stat`: 91 insertions, 7 deletions). All 32 R chunks parse. Sources are the survival template `sim_fs_maxeffCons_fb_mr_field_m1_template.qmd` (m1) at HEAD.

- **E1 — focus and band knobs.** `:135` `sg_focus <- "maxeffCons"` replaced by m1's focus knob, guard and band block (m1 `:318` `sg_focus <- .env_chr("FS_S7_FOCUS", "maxeffCons")`; `:333` the six-focus `stopifnot`; `:334` `.band_foci <- c("effMaxSG", "effMinSG")`; `:342–343` `effect_neighborhood <- .env_num("FS_S7_NBHD", 0.10)` + `stopifnot(is.finite(...), >= 0, < 1)`; `:350–353` the "FS_S7_NBHD is set but sg_focus does not consult effect_neighborhood" stop), renamed `FS_MD_FOCUS` / `FS_MD_NBHD` with defaults `"maxeffCons"` / `0.10`; the band is read there, before the stem, as in m1. `:214` `effect_neighborhood <- 0.10` replaced by m1's pointer comment (`:565–567`). `:213` `selection_rule <- "neighborhood"` kept, followed by m1's guard (`:555–558`: `selection_rule` other than `"neighborhood"` only with a band focus). Stem tagging: m1's `nbhd_tag` (`:443–444`, `_nb%02d` only when non-default) inserted before `rds_stem` and the stem format gains `%s` before the campaign tag (m1 `:450–453`), so defaults reproduce `mdf1`'s stems exactly and the campaign stem is `fs_effMaxSG_mr_field_md40_knoise0_n500_nb20_mdsgnb20` (the survival form `fs_effMaxSG_..._nb20_p12x20`). No `FS_MD_RULE` knob: §1.3 shows the survival `selection_rule` is the literal `"neighborhood"`.
- **E2 — complement-scale knob.** After `:271`, m1 `:655–659` (`mr_field_scalec <- .env_chr("FS_S7_FIELD_SCALEC", "selected")`; `stopifnot(mr_field_scalec %in% c("none", "selected"))`) as `FS_MD_FIELD_SCALEC`, default `"selected"`; `field_scale_complement = mr_field_scalec,` added to `mr_inference_args` (`:285`; m1 `:688`).
- **E3 — recorder.** Stage 0 §5.2 (a)–(d), verbatim from m1: (a) after `:518` the nine `fld_Hc_*_s` NA columns (m1 `:928–936`); (b) after `:530` the nine `fld_joint_s_*` NA columns (m1 `:973–978`); (c) after `:733` the nine `rec$fld_Hc_*_s <- fc$*_s %||% NA_real_` fills (m1 `:1259–1268`); (d) after `:751` the `joint_s` fill block (m1 `:1305–1316`).
- **E4 — invariants.** After `:827`, the pairs `c("fld_Hc_lo2s_s", "fld_Hc_hi2s_s")`, `c("fld_Hc_lo_se_s", "fld_Hc_hi_se_s")`, `c("fld_Hc_lo1s_s", "fld_Hc_up1s_s")` in the existing pair form. **Finding F3:** m1's own `.ci_check` (`:1482–1487`) carries no `_s` pair (p12x20's `_s` invariants live in `gate2F.R:168, :172`), so these three lines follow the MD template's pattern rather than a survival source line.
- **E5 — `meta` and poolability.** `effect_neighborhood = effect_neighborhood,` after `selection_rule` (`:895`; m1 `:1631–1632`) and `field_scale_complement = mr_field_scalec,  # record only (E1)` after `return_reselection` (`:912`; m1 `:1616`); the combine-mode keys (`:945–948`) gain `"effect_neighborhood", "selection_rule", "field_scale_complement"` (m1 `:1665–1667` gates `sg_focus` and `effect_neighborhood`; `selection_rule` and `field_scale_complement` added here so batches of different constructions never pool). The batch banner (`:307–311`) echoes `scalec`, `focus`, `nbhd`, `rule` (m1 `:728` echoes them).

## 1.5 Scripts — commit `92d6ceb7`, under `scripts_mdsgnb20/`

- `mem_sampler.sh` ← `scripts_mdf1/mem_sampler.sh`; the one change is Linux `ps -eo rss,comm` for macOS `ps -Ao rss,comm`.
- `smoke_identity.R` ← `scripts_mdf1/gate2_check.R:22–53` (the script that wrote `mdf1`'s `gate2_flips.txt`; `identity_postchange.R` compares `fs_sim_bias_coverage()` tables, not bundles). Pairs every `mdf1` recorder column except `*_secs` (numeric ≤ 1e-8 relative, NA ≡ NA; `sg_def`/`covs` as term sets; other character columns identical); message (`mr_msg`, `err_msg`, `fb_err`) and joined FB columns are reported apart, per `mdf1`'s Gate 2 convention (`REPORT_continuous_field_gate2_2026-09-07.md:5`: "minus timings, messages and the FB columns"); classification exactly as `gate2_check.R:41–51`; the §1.6(b) field-s checks; a `rule` mode for §1.6(c).
- `gate2.R` ← `quarto/simulations/gbsg_020/scripts_p12x20/gate2F.R`, pointed at `mr_md_harm/fs_effMaxSG_mr_field_<cell>_nb20_mdsgnb20_d5000/` with the §2.3 checks (batch files match the combined bundle on every column; `meta` rule / `field_scale_complement` / `pkg_version` / `hostname`; same draws vs `mdf1` in both directions with `n_true` identical and the oracle columns ≤ 1e-8, complement only in the null cell; §1.6(b) on filled replicates; MR failures ≤ max(20, 2 × `mdf1`'s)); identities on the identity scale (`est2_s + lam_mean_s == est2 + lam_mean`).
- `run_mdsgnb20.sh` ← `quarto/simulations/gbsg_020/scripts_p12x20/run_p12x20.sh` (preflight, sequencing, heartbeats, halt, explicit-path commits) with `scripts_mdf1/run_cell.sh:18–31`'s three renders per cell, GNU `timeout` (`MDSG_TIMEOUT`), threads exported, `mem_sampler.sh` per render, `LOG_mdsgnb20_progress.txt`, `HALT_mdsgnb20.md`, the go's ceiling (`MDSG_CEILING`), cells in `mdf1`'s order, per-cell commits of the result directory, the combine HTML, the Gate 2 report section and the progress log. Raw logs under `logs_mdsgnb20/` (untracked).

## 1.6 Smoke — GATE PASS

Renders (`logs_mdsgnb20/smoke.sh`, sequential, `FS_MD_CAMPAIGN=mdsmoke FS_MD_START=1 FS_MD_NSIMS=20 FS_MD_WORKERS=20 FS_MD_CI=field`, template at `2cffb95f`, installed Built `2026-09-16 05:57:14 UTC`): md40 n500 77 s, md120 n500 76 s, null n500 73 s, md40 n700 89 s; the rule render (md40 n500, `FS_MD_FOCUS=effMaxSG FS_MD_NBHD=0.20`) 85 s. Every render exit 0; every replicate `DETECTED`; no `CONFIG-ERROR`.

**(a) Defaults, all four cells** (`scripts_mdsgnb20/smoke_identity.R ... identity`; 119 paired `mdf1` recorder columns = 108 numeric ≤ 1e-8 relative + 11 character; `*_secs` excluded; messages and FB columns reported apart and identical / joined-only):

| cell | rows compared | within tolerance | enumerated: pure label tie | label tie, target moves | MR-numerics | selection / detection flips | max rel diff among identical rows |
|---|---|---|---|---|---|---|---|
| md40 n500 | 20 | 17 | 2 (sims 5, 8) | 1 (sim 9) | 0 | **0** | 2.85e-12 |
| md120 n500 | 20 | 19 | 1 (sim 19) | 0 | 0 | **0** | 2.11e-12 |
| null n500 | 20 | 14 | 6 (sims 4, 6, 8, 9, 11, 13) | 0 | 0 | **0** | 4.00e-11 |
| md40 n700 | 20 | 20 | 0 | 0 | 0 | **0** | 2.15e-11 |

Every enumerated row is a `mdf1` Gate 2 class. Sims 8 and 9 (md40 n500, null) are the rows `mdf1`'s own Gate 2 enumerated as label ties (`gate2_flips.txt`); sims 4 and 6 (null) differ only by term order inside the `label` string; sims 5, 11, 13, 19 differ only in a re-selection top-3 label (`p_lab1..3`), the same duplicate-membership ties (`{str2}` ≡ `!{preanti <= 0}`). `mr_msg`, `err_msg`, `fb_err` identical on every row; the 13 FB columns are finite on `mdf1`'s md40 n500 rows (joined) and NA here (`FS_MD_FB` unset), as designed.

**(b) Field-s wiring**, on every smoke replicate whose complement field block was filled (20 of 20 in each of the five renders): the nine `fld_Hc_*_s` and nine `fld_joint_s_*` columns finite — PASS; `fld_Hc_lo1s_s ≤ fld_Hc_up1s_s` and `fld_Hc_lo2s_s ≤ fld_Hc_hi2s_s` — PASS; joint and joint_s draw counts agree on 20/20 rows and their Bonferroni harm bounds are identical (max |diff| 0) — PASS; field-s inverted around the same β̃ᶜ (`est2_s + lam_mean_s` vs `est2 + lam_mean`, max |diff| ≤ 7.1e-15) — PASS; `meta$field_scale_complement == "selected"` — PASS (and `meta` now carries `effect_neighborhood 0.1`, `selection_rule neighborhood`, `pkg_version 0.3.5`, `hostname pop-os`, `r_version 4.6.1`).

**(c) The campaign rule, md40 n500, sim_id 1–20** (`... effMaxSG 0.20 mdsmoke 20 rule`): stem `fs_effMaxSG_mr_field_md40_knoise0_n500_nb20_mdsmoke`, `meta` `sg_focus effMaxSG`, `effect_neighborhood 0.2`, `selection_rule neighborhood` — PASS; `n_true` identical on every row — PASS; the eight oracle columns within 1e-8 relative (max 2.11e-12) — PASS; (b)'s checks — PASS. **Fact:** against `mdf1`, |Ĥ| (`n_harm`) grew on 17 replicates, stayed on 3, shrank on 0 (the same rule string on those 3); mean |Ĥ| 120.3 here vs 74.2 under `maxeffCons` (`n_sel` moves identically).

Smoke timing at 20 workers (one round of 20; loaded; not a calibration figure): `fit_mr_secs` mean 28.3 / 30.4 / 27.8 / 39.4 s (md40 n500 / md120 / null / n700; `fld_H_secs` 12.1–14.8, `fld_Hc_secs` 1.7–2.6); under `effMaxSG` at md40 n500 38.7 s (`fld_H_secs` 21.8, `fld_Hc_secs` 2.1) — the band rule's field pass is longer than `maxeffCons`'s on the same draws.

Smoke outputs, untracked (listed for §3.5 deletion): `mr_md_harm/fs_maxeffCons_mr_field_{md40_knoise0_n500,md120_knoise0_n500,mdnull_knoise0_n500,md40_knoise0_n700}_mdsmoke_d5000/`, `mr_md_harm/fs_effMaxSG_mr_field_md40_knoise0_n500_nb20_mdsmoke_d5000/`, `logs_mdsgnb20/mdsmoke_*.html`, `logs_mdsgnb20/mdsmoke_*.log`, `logs_mdsgnb20/smoke_identity_*.txt`, `logs_mdsgnb20/smoke_rule_md40_n500.txt`, `logs_mdsgnb20/smoke.sh`, `logs_mdsgnb20/smoke_driver.log`.

## 1.7 Calibration — md40 n700 under the campaign knobs (`FS_MD_FOCUS=effMaxSG FS_MD_NBHD=0.20 FS_MD_FIELD_SCALEC=selected FS_MD_CI=field FS_MD_FB=none`), 3 × W replicates, memory sampled every 5 s (`logs_mdsgnb20/calib.sh`; summary computed by `calib_summary.R`, pasted below)

| W | replicates | render wall (s) | fit_mr_secs mean / median / p90 / max (s) | ratio to 16-worker mean | fld_H_secs / fld_Hc_secs mean (s) | peak summed RSS (MB) | replicates per minute |
|---|---|---|---|---|---|---|---|
| 16 | 48 | 207 | 50.76 / 50.10 / 57.21 / 60.37 | 1.000 | 26.22 / 3.15 | 19305 | 13.91 |
| 32 | 96 | 233 | 54.35 / 54.90 / 60.80 / 63.43 | 1.071 | 28.06 / 3.43 | 36929 | 24.72 |
| 63 | 189 | 362 | 74.47 / 73.68 / 95.48 / 109.38 | 1.467 | 36.78 / 4.92 | 71835 | 31.33 |

Fixed render overhead (16-worker render, wall minus 3 rounds x mean): 54.7 s. Wall-based loop cost per replicate (render wall minus overhead, over replicates): W 16: 3.172 s; W 32: 1.857 s; W 63: 1.626 s.

Projection at W = 63 (loop cost 1.626 s per replicate at md40 n700, scaled per cell by mdf1 fit_mr_secs ratios md40 n500 0.8523, md120 n500 0.9773, null n500 0.8352, md40 n700 1.0000; overhead 55 s per batch render; combine render 90 s):

| cell | ratio | batch of 1,000 (s) | batch (min) | cell: 2 batches + combine (min) |
|---|---|---|---|---|
| md40 n500 | 0.8523 | 1440 | 24.0 | 49.5 |
| md120 n500 | 0.9773 | 1644 | 27.4 | 56.3 |
| null n500 | 0.8352 | 1413 | 23.5 | 48.6 |
| md40 n700 | 1.0000 | 1681 | 28.0 | 57.5 |

**Projection: 12714 s = 211.9 min = 3.53 h.  Ceiling (1.5x): 19071 s = 317.9 min = 5.30 h.  Per-render timeout (2 x the longest projected batch, at least 20 min): 3361 s = 56.0 min.**
Busiest-worker model for comparison (ceil(1000/W) rounds x fit_mr_secs mean + overhead, the template's own projection form): at W = 63, md40 n700 batch 20.8 min; total 158.9 min. It ignores per-round straggler loss (p90/mean 1.28 at W = 63) and is the lower figure.
MDSG_WORKERS=63 MDSG_TIMEOUT=3361 MDSG_CEILING=19071

- **W = 63.** It minimizes projected wall: 31.3 replicates per minute against 24.7 at 32 and 13.9 at 16. **What limits it:** the template's worker cap (physical cores − 1 = 63, `sim_fs_maxeffCons_mr_field_md_template.qmd:115–116`), i.e. cores; scaling falloff is visible (per-replicate `fit_mr_secs` 1.47× the 16-worker mean, p90/mean 1.28) but throughput still rises to the cap; memory is not binding (peak summed RSS 72 GB of 251 GB, ≈ 1.1 GB per worker).
- **Projection for Stage 2 at W = 63: 212 min (3.53 h)** for four cells × (two 1,000-replicate batches + combine), per-cell cost scaled by `mdf1`'s `fit_mr_secs` ratios (Stage 0 §8.1: 15.0 / 17.2 / 14.7 / 17.6 s) and the fixed render overhead from these renders (54.7 s per batch render; 90 s allowed per combine render, above `mdf1`'s 20–24 s on the Mac). The wall-based per-replicate loop cost (1.626 s at W = 63) carries the observed straggler loss; the template's busiest-worker model gives 159 min and is the lower figure.
- **Limits:** projection **12,714 s (211.9 min)**; ceiling 1.5 × projection = **19,071 s (317.9 min, 5.30 h)**; per-render timeout 2 × the longest projected batch (28.0 min) = **3,361 s (56.0 min)**, above the 20-min floor. Runner environment: `MDSG_WORKERS=63 MDSG_TIMEOUT=3361 MDSG_CEILING=19071`.
- Calibration outputs, untracked (listed for §3.5 deletion): `mr_md_harm/fs_effMaxSG_mr_field_md40_knoise0_n700_nb20_mdcal{16,32,63}_d5000/`, `logs_mdsgnb20/mdcal{16,32,63}.{html,log,peak_mb}`, `logs_mdsgnb20/calib.sh`, `logs_mdsgnb20/calib_driver.log`, `logs_mdsgnb20/calib_summary.txt`.
- Stage 3 preparation done in the waiting time, no compute beyond a read-only render: `summary_continuous_field_mdsgnb20.qmd` (the §3.1 transplant of `summary_continuous_field_mdf1.qmd`, untracked until Stage 3) dry-rendered against the md40 n500 smoke bundle (`MDSG_SUMMARY_TAG=mdsmoke MDSG_SUMMARY_GLOB=res_1_20`, exit 0, 266 extract rows) under `logs_mdsgnb20/dryrun/` (untracked; listed for §3.5 deletion).

## 1.8 Gate 1 — PASS; the go was given in advance

Every Stage 1 gate is green (§1.1, §1.2, §1.6) and the §1.7 projection (3.53 h) is under 12 hours. **Larry's instruction of 2026-09-15 (overnight), given in advance on exactly this condition:** "if every Stage 1 gate is green and the §1.7 projection for Stage 2 is under 12 hours, do not stop at Gate 1 … run Stage 2 at the worker count chosen in §1.7, with the ceiling and per-render timeout stated in the Gate 1 record, and run Stage 3 once Stage 2 is green"; FB off (the template default; no FB rows in Stage 3); no reduction of cells, replicates, gates or knobs. Stage 2 therefore launches after this record's commit with `MDSG_WORKERS=63 MDSG_TIMEOUT=3361 MDSG_CEILING=19071` and the §1.3 rule knobs, cells in `mdf1`'s order.

## Findings

- **F0 (install):** `devtools::install(quick = TRUE, upgrade = "never")` is rejected by this devtools (`upgrade` must be TRUE/FALSE/NA); `upgrade = FALSE` was used. Same intent, different literal.
- **F1 (§1.3):** the MD design's identification floors differ from the survival campaigns' (`n.min 60` vs `NULL`; `d0.min/d1.min 12` vs `10`; thresholds 30/10 on the MD scale vs 0.90/0.80 on the HR scale; cut grids on `age`/`preanti` vs `er`). By the task's rule they stay the MD template's.
- **F2 (§1.3):** `field_decompose` is FALSE here and was TRUE in `p12x20`/`cert20` (diagnostic only; excluded from tables by convention 7).
- **F3 (§1.4 E4):** the survival template's `.ci_check` carries no `_s` invariant pairs (they live in `gate2F.R`); the three `_s` pairs were added in the MD template's own pair form.
- **F4 (§1.6):** the `mdf1` bundle `meta` and poolability keys never carried `field_scale_complement`, `effect_neighborhood` (Stage 0 findings 7, 11); the smoke `meta` now does (`effect_neighborhood 0.1`, `selection_rule neighborhood`, `field_scale_complement selected`).
- **F5 (§1.6(c), fact):** under `effMaxSG` at ε = 0.20 the selected subgroup is larger on 17 of 20 md40 n500 replicates and unchanged on 3 (mean |Ĥ| 120 vs 74); under the calibration's md40 n700 the mean |Ĥ| is 118–121 of 700.
- **F6 (§1.7):** the per-replicate cost on this machine is well above the Mac's (`fit_mr_secs` 50.8 s at 16 workers vs 17.6 s at 13 on the M4 Max for the same cell under `maxeffCons`; the field pass 26 s vs 12.8 s). Part is the rule (§1.6: 38.7 vs 28.3 s at 20 workers on md40 n500) and part is the host; the projection uses this machine's measurements only.
- **F7 (§1.7):** at W = 63 the p90/mean of `fit_mr_secs` is 1.28 and the maximum 109 s, so the busiest-worker projection (159 min) undershoots; the ceiling and timeout are set from the wall-based projection (212 min).

## Untracked outputs of this stage (for §3.5 deletion)

`logs_mdsgnb20/` in full (smoke, calibration, dry-run logs and HTML, `smoke.sh`, `calib.sh`); `mr_md_harm/*_mdsmoke_d5000/` (5 directories) and `mr_md_harm/*_mdcal{16,32,63}_d5000/` (3 directories). The Stage 3 summary draft `summary_continuous_field_mdsgnb20.qmd` is untracked until Stage 3 commits it.

## git log --oneline for this stage

```
92d6ceb7 scripts_mdsgnb20 (TASK_md_field_rerun_2026-09-15 §1.5, transplants): mem_sampler.sh (scripts_mdf1, Linux ps); smoke_identity.R (gate2_check.R's pairing proof and classification, every column except *_secs, field-s wiring, rule mode); gate2.R (p12x20's gate2F.R pointed at the MD bundles with the §2.3 checks, same draws both directions vs mdf1); run_mdsgnb20.sh (run_p12x20.sh sequencing with run_cell.sh's render lines, GNU timeout, progress log, halt file, per-cell explicit-path commits)
2cffb95f MD template E1-E5 (TASK_md_field_rerun_2026-09-15 §1.4, transplanted from the survival m1 template): FS_MD_FOCUS / FS_MD_NBHD knobs with the six-focus guard, the band-foci eps guard and the _nb tag on the stem; selection_rule guard; FS_MD_FIELD_SCALEC knob (default selected) passed as field_scale_complement; nine fld_Hc_*_s and nine fld_joint_s_* recorder columns with their fill blocks; _s interval invariant pairs; effect_neighborhood / field_scale_complement in meta and effect_neighborhood / selection_rule / field_scale_complement in the poolability keys; knob echo. Defaults reproduce mdf1's stem and settings.
bb84120e Add TASK_md_field_rerun_2026-09-15 as received
<this record: the next commit>
```
