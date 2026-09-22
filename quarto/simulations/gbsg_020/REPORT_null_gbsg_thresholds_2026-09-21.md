# REPORT — campaign (A): the structural null at c1 = 1.25, c2 = 1.00 (identification only)

**Task:** `dev/tasks/TASK_null_gbsg_thresholds_2026-09-21_v2.md` (committed `8a318f26`) ·
**Opened:** 2026-09-21 · **Branch:** `feature/glm-extension` · **Study directory:**
`quarto/simulations/gbsg_020` · **Campaign tag:** `nullc125` · **Machine:** pop-os.

**What this is.** The six `nullid` cells — no planted region, uniform treatment benefit at
super-population marginal Cox HR 0.657 and 0.721, n 500 / 1000 / 1500, 2,000 replicates, FS / DINA /
GRF on identical draws, `effMaxSG` at ε = 0.20 — rerun at a stricter screen: **(c1, c2, p⋆) =
(1.25, 1.00, 0.90)**, c1 and c2 passed explicitly. DINA's proposal floor follows c1
(`m_diff = log(1.25) = 0.2231`), as does GRF's effect-scale admission floor; `dmin.grf` stays 0.0.
Everything else is `nullid`'s. The 0.90 / 0.80 rows below are read, read-only, from `nullid`'s
committed bundles.

**Scope held.** No MR anywhere (no multiplier resampling, field, field-s, IJ or Bonferroni product),
no bootstrap, no cross-validation. No edit to `R/`, no install. `nullid` and every committed cell
untouched; nothing committed was re-run. No R CMD check, vignette build or test suite. Nothing
fetched, pulled or pushed. Every `git add` named its paths.

**Outcome.** 18 of 18 runs; Gate A PASS on every run, Gate C PASS on every cell; no halt record.
Grid wall 4,002 s of render (driver 4,005 s end to end, 1.11 h) against a 1.12 h projection.

---

## 0. Preconditions and build

| check | result |
|---|---|
| host / branch | `pop-os` / `feature/glm-extension` |
| HEAD contains `cabbd3f3` | yes |
| tracked modifications / R, Rscript, quarto processes | none / none |
| last commit touching `R/` | 2026-09-18T20:02:11-07:00 (= 2026-09-19 03:02 UTC) |
| installed forestsearch | **0.3.5.9000**, `Built: R 4.6.1; ; 2026-09-19 03:54:02 UTC; unix` — newer than the last `R/` commit; no install needed |
| R | 4.6.1 (2026-06-24) |
| machine | AMD Ryzen Threadripper PRO 5995WX, 64 physical cores, 251 GB; 64 workers |

## 1. The threshold knobs, and that they are inert at the default

Template `sim_fs_maxeffCons_fb_mr_field_m1_template.qmd`, add-only (commit `6292e8c2`):

- **`FS_S7_C1` / `FS_S7_C2`**, default 0.90 / 0.80 — the literals formerly at `:601-602` — fed to the
  same `hr_threshold` / `hr_consistency` the `forestsearch()` call passes. They are read beside the
  stem (they tag it), so `:601-602` now read `hr_threshold <- thr_c1`, `hr_consistency <- thr_c2`.
  Setting one alone stops the render, and so does C2 > C1 (both verified: rc 1 with the stated
  message, no output written). p⋆ stays the literal 0.90.
- **Stem token** `_c%03dc%03d` (`_c125c100`) after `_nomr`, **empty at the default**, so no committed
  stem changes and no stricter-screen bundle can share a stem or `combine_glob` with `nullid`.
- **Meta:** `c1`, `c2`, `pstar`, `dmin_grf`, `dina_m_diff` (= `log(c1)`). No earlier bundle records
  its thresholds; there they were 0.90 / 0.80 / 0.90.
- **Recorder:** `c1_resolved` / `c2_resolved` from `fs.est$args_call_all$hr.threshold` /
  `$hr.consistency`, read before the no-detection return; **recorded on 2,000 / 2,000 replicates on
  all three engines in all 18 runs — no NA anywhere.** `itt_est` / `itt_se`: `.cox_hr_ci()` on the
  whole trial, computed right after the draw, before the search.
- The `Thresholds:` audit line is printed in every render.

**Inertness, same machine** (`scripts_dinamr/logs/nullc125_inertness.txt`,
`scripts_dinamr/nullthr_inertness.R`): `null0657_n500`, `sim_id` 1–20, all three identifiers,
rendered from HEAD before the edit (`nullc125inertpre`), after it with the knobs unset
(`nullc125inertunset`) and with them set to 0.90 / 0.80 (`nullc125inertexpl`), each quickrun-tagged.
**Identical on 169 of 169 shared non-timing columns and on the truth object, for every identifier and
both comparisons;** the four added columns agree between the two post-edit renders. The pre-edit
declarations, 6 / 9 / 20 of 20, reproduce `nullid`'s Mac smoke exactly.

## 2. Driver and gates

- **`scripts_dinamr/nullthr.sh`** — bash port of `nullid.sh` (pop-os has no zsh; `nullid.sh` also
  calls `sysctl`, `memory_pressure` and BSD `date -r`). It reads `thresholds <c1> <c2>` and the cells
  in run order from `nullthr.cells`; invokes `render.sh` as `campaign_p12x20.sh` does
  (`bash "$RENDER"`), at 64 workers (`run_p12x20.sh`); unsets every inherited `FS_S7_*` first
  (which includes `nullid`'s `FS_S7_Z1Q` and `FS_S7_ER_JCUTS`); keeps `nullid.sh`'s pinned knob set,
  `FS_S7_DGM=null`, `FS_S7_MR=FALSE`, `FS_S7_FB=none` and halt-and-continue. Timeouts: 90 min per
  render (`timeout(1)`), 6 h for the campaign (watchdog). It commits each completed cell by named
  paths.
- **Gate A** (`nullthr_gateA.R`, 42 checks FS / 36 DINA / 37 GRF): `nullid_gateA.R`'s invariants,
  plus meta `c1`/`c2` equal the run's knobs, `pstar` 0.90, `dmin_grf` 0, `dina_m_diff = log(c1)`;
  resolved thresholds equal the knobs on every replicate where recorded; `itt_est` finite and
  `itt_se` finite and positive on every replicate; `status == "DETECTED"` exactly where
  `detected == 1`; `p_sel >= p⋆` on every FS declaration; `er_jcuts == 10` in the meta.
- **Two of `nullid_gateA.R`'s checks were counts or shares, not invariants, and were replaced** —
  the same mistake §3.1 of the `nullid` report records:
  1. "unadjusted within-region estimate present on ≥ 95% of declaring replicates" → **the estimate,
     its SE and its one-sided bound are finite together on declaring rows, and absent on
     non-declaring rows**; coverage is printed, not gated. (It was 100% in every run of both
     campaigns.)
  2. GRF "`admitted_n` recorded on some replicate (> 0)" → **`admitted_n ≥ 1` on every declaring
     GRF replicate** (the pick was admitted). A first draft asserting `admitted_n` finite on every
     replicate was checked against the `nullid` bundles before use and is false there (NA on 30–819
     non-declaring replicates per cell), so it was never run.
  Every new invariant was checked against the 18 `nullid` bundles before the campaign.
- **Gate C** (`nullthr_gateC.R`, 16 checks): `itt_est` and `itt_se` **identical** (`identical()`)
  across the three identifiers on all 2,000 rows, same `sim_id`, same `k_treat`, same truth object.
  Unlike `nullid`'s `or_Hc_*` fingerprint it needs no declaration, so it has no row floor.
- Driver, gates and cells: commit `072a052a`.

## 3. Smoke, projection, grid

**Smoke** (`scripts_dinamr/logs/nullc125_smoke_projection.txt`, commit `156c6c6b`), through the real
driver, `null0657_n500` at (1.25, 1.00), 20 replicates, 64 workers: per-replicate search FS 2.427 s,
DINA 0.146 s, GRF 3.964 s (6.537 s across the three; light load); declarations 2 / 2 / 10 of 20;
all gates PASS. **Projection 1.12 h** (n-scaling as `nullid`; load factor 1.17–2.57 from the
measured 63-worker contention; 18 × 60 s render overhead) — under the 4 h bound, so the grid ran
straight away.

**Grid:** 6 cells × 3 identifiers × 2,000 replicates, `nullid`'s cell order, driver log
`scripts_dinamr/logs/nullc125.driver.log`. Cell commits `f7bbfa40`, `f16f7257`, `cbfd0f34`,
`a589b879`, `083033b5`, `eaff6ab0`. No failure, no halt, no retry.

---

## 4. Results

All tables are produced by `scripts_dinamr/nullthr_findings.R` (captured output
`scripts_dinamr/logs/nullc125_findings.md`). **Every declaration is a false declaration.**

### Table 1 — screen x cell x identifier

Every declaration is false here. `k` is the number of declaring replicates a conditional summary rests on. Lower-bound shares: *cond.* over declaring replicates with a finite bound, *uncond.* over all 2,000 (a replicate that declares nothing scores 0).

| cell | identifier | screen c1/c2 | declared / 2000 | rate [Wilson 95%] | mean \|H\|/n (MC SE) | spec uncond. (MC SE) | spec cond. (MC SE) | med HR(H) unadj. | med true β(Ĥ) (HR scale) [populated] | LB ≥ 1.00 cond. | LB ≥ 1.25 cond. | LB ≥ 1.00 uncond. | LB ≥ 1.25 uncond. |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| null0657_n500 | FS | 0.90/0.80 | 894 | 0.4470 [0.4253, 0.4689] | 0.2107 (0.0024) [k=894] | 0.9058 (0.0026) | 0.7893 (0.0024) [k=894] | 1.390 [k=894] | 0.616 [894/894] | 93/894 = 0.1040 | 7/894 = 0.0078 | 93/2000 = 0.0465 [0.0381, 0.0566] | 7/2000 = 0.0035 [0.0017, 0.0072] |
| null0657_n500 | FS | 1.25/1.00 | 190 | 0.0950 [0.0829, 0.1086] | 0.1883 (0.0039) [k=190] | 0.9821 (0.0013) | 0.8117 (0.0039) [k=190] | 1.772 [k=190] | 0.616 [190/190] | 153/190 = 0.8053 | 7/190 = 0.0368 | 153/2000 = 0.0765 [0.0656, 0.0890] | 7/2000 = 0.0035 [0.0017, 0.0072] |
| null0657_n500 | DINA | 0.90/0.80 | 1116 | 0.5580 [0.5361, 0.5796] | 0.1913 (0.0019) [k=1116] | 0.8933 (0.0024) | 0.8087 (0.0019) [k=1116] | 1.197 [k=1116] | 0.618 [1116/1116] | 93/1116 = 0.0833 | 16/1116 = 0.0143 | 93/2000 = 0.0465 [0.0381, 0.0566] | 16/2000 = 0.0080 [0.0049, 0.0130] |
| null0657_n500 | DINA | 1.25/1.00 | 234 | 0.1170 [0.1036, 0.1318] | 0.1586 (0.0025) [k=234] | 0.9814 (0.0012) | 0.8414 (0.0025) [k=234] | 1.455 [k=234] | 0.611 [234/234] | 31/234 = 0.1325 | 6/234 = 0.0256 | 31/2000 = 0.0155 [0.0109, 0.0219] | 6/2000 = 0.0030 [0.0014, 0.0065] |
| null0657_n500 | GRF | 0.90/0.80 | 1884 | 0.9420 [0.9309, 0.9514] | 0.2074 (0.0016) [k=1884] | 0.8046 (0.0019) | 0.7926 (0.0016) [k=1884] | 1.206 [k=1884] | 0.612 [1884/1884] | 113/1884 = 0.0600 | 20/1884 = 0.0106 | 113/2000 = 0.0565 [0.0472, 0.0675] | 20/2000 = 0.0100 [0.0065, 0.0154] |
| null0657_n500 | GRF | 1.25/1.00 | 1369 | 0.6845 [0.6638, 0.7045] | 0.1697 (0.0012) [k=1369] | 0.8839 (0.0019) | 0.8303 (0.0012) [k=1369] | 1.380 [k=1369] | 0.608 [1369/1369] | 113/1369 = 0.0825 | 20/1369 = 0.0146 | 113/2000 = 0.0565 [0.0472, 0.0675] | 20/2000 = 0.0100 [0.0065, 0.0154] |
| null0721_n500 | FS | 0.90/0.80 | 1356 | 0.6780 [0.6572, 0.6981] | 0.2212 (0.0020) [k=1356] | 0.8500 (0.0027) | 0.7788 (0.0020) [k=1356] | 1.419 [k=1356] | 0.687 [1356/1356] | 225/1356 = 0.1659 | 21/1356 = 0.0155 | 225/2000 = 0.1125 [0.0994, 0.1271] | 21/2000 = 0.0105 [0.0069, 0.0160] |
| null0721_n500 | FS | 1.25/1.00 | 465 | 0.2325 [0.2145, 0.2515] | 0.1941 (0.0030) [k=465] | 0.9549 (0.0020) | 0.8059 (0.0030) [k=465] | 1.745 [k=465] | 0.686 [465/465] | 382/465 = 0.8215 | 21/465 = 0.0452 | 382/2000 = 0.1910 [0.1744, 0.2088] | 21/2000 = 0.0105 [0.0069, 0.0160] |
| null0721_n500 | DINA | 0.90/0.80 | 1480 | 0.7400 [0.7203, 0.7588] | 0.2022 (0.0019) [k=1480] | 0.8504 (0.0024) | 0.7978 (0.0019) [k=1480] | 1.334 [k=1480] | 0.688 [1480/1480] | 269/1480 = 0.1818 | 46/1480 = 0.0311 | 269/2000 = 0.1345 [0.1202, 0.1502] | 46/2000 = 0.0230 [0.0173, 0.0305] |
| null0721_n500 | DINA | 1.25/1.00 | 426 | 0.2130 [0.1956, 0.2315] | 0.1643 (0.0019) [k=426] | 0.9650 (0.0016) | 0.8357 (0.0019) [k=426] | 1.514 [k=426] | 0.683 [426/426] | 106/426 = 0.2488 | 21/426 = 0.0493 | 106/2000 = 0.0530 [0.0440, 0.0637] | 21/2000 = 0.0105 [0.0069, 0.0160] |
| null0721_n500 | GRF | 0.90/0.80 | 1969 | 0.9845 [0.9781, 0.9891] | 0.2156 (0.0018) [k=1969] | 0.7878 (0.0019) | 0.7844 (0.0018) [k=1969] | 1.323 [k=1969] | 0.684 [1969/1969] | 268/1969 = 0.1361 | 49/1969 = 0.0249 | 268/2000 = 0.1340 [0.1198, 0.1496] | 49/2000 = 0.0245 [0.0186, 0.0322] |
| null0721_n500 | GRF | 1.25/1.00 | 1670 | 0.8350 [0.8181, 0.8506] | 0.1798 (0.0012) [k=1670] | 0.8499 (0.0018) | 0.8202 (0.0012) [k=1670] | 1.427 [k=1670] | 0.681 [1670/1670] | 268/1670 = 0.1605 | 49/1670 = 0.0293 | 268/2000 = 0.1340 [0.1198, 0.1496] | 49/2000 = 0.0245 [0.0186, 0.0322] |
| null0657_n1000 | FS | 0.90/0.80 | 637 | 0.3185 [0.2984, 0.3392] | 0.1668 (0.0023) [k=637] | 0.9469 (0.0019) | 0.8332 (0.0023) [k=637] | 1.240 [k=637] | 0.616 [637/637] | 26/637 = 0.0408 | 4/637 = 0.0063 | 26/2000 = 0.0130 [0.0089, 0.0190] | 4/2000 = 0.0020 [0.0008, 0.0051] |
| null0657_n1000 | FS | 1.25/1.00 | 65 | 0.0325 [0.0256, 0.0412] | 0.1362 (0.0045) [k=65] | 0.9956 (0.0006) | 0.8638 (0.0045) [k=65] | 1.721 [k=65] | 0.610 [65/65] | 53/65 = 0.8154 | 4/65 = 0.0615 | 53/2000 = 0.0265 [0.0203, 0.0345] | 4/2000 = 0.0020 [0.0008, 0.0051] |
| null0657_n1000 | DINA | 0.90/0.80 | 593 | 0.2965 [0.2769, 0.3169] | 0.1578 (0.0022) [k=593] | 0.9532 (0.0017) | 0.8422 (0.0022) [k=593] | 1.029 [k=593] | 0.615 [593/593] | 12/593 = 0.0202 | 4/593 = 0.0067 | 12/2000 = 0.0060 [0.0034, 0.0105] | 4/2000 = 0.0020 [0.0008, 0.0051] |
| null0657_n1000 | DINA | 1.25/1.00 | 36 | 0.0180 [0.0130, 0.0248] | 0.1335 (0.0040) [k=36] | 0.9976 (0.0004) | 0.8665 (0.0040) [k=36] | 1.312 [k=36] | 0.603 [36/36] | 3/36 = 0.0833 | 0/36 = 0.0000 | 3/2000 = 0.0015 [0.0005, 0.0044] | 0/2000 = 0.0000 [0.0000, 0.0019] |
| null0657_n1000 | GRF | 0.90/0.80 | 1662 | 0.8310 [0.8139, 0.8468] | 0.1895 (0.0017) [k=1662] | 0.8425 (0.0021) | 0.8105 (0.0017) [k=1662] | 1.002 [k=1662] | 0.613 [1662/1662] | 38/1662 = 0.0229 | 6/1662 = 0.0036 | 38/2000 = 0.0190 [0.0139, 0.0260] | 6/2000 = 0.0030 [0.0014, 0.0065] |
| null0657_n1000 | GRF | 1.25/1.00 | 576 | 0.2880 [0.2686, 0.3082] | 0.1325 (0.0013) [k=576] | 0.9618 (0.0014) | 0.8675 (0.0013) [k=576] | 1.333 [k=576] | 0.606 [576/576] | 42/576 = 0.0729 | 6/576 = 0.0104 | 42/2000 = 0.0210 [0.0156, 0.0283] | 6/2000 = 0.0030 [0.0014, 0.0065] |
| null0721_n1000 | FS | 0.90/0.80 | 1234 | 0.6170 [0.5955, 0.6381] | 0.1930 (0.0021) [k=1234] | 0.8809 (0.0025) | 0.8070 (0.0021) [k=1234] | 1.213 [k=1234] | 0.688 [1234/1234] | 85/1234 = 0.0689 | 9/1234 = 0.0073 | 85/2000 = 0.0425 [0.0345, 0.0523] | 9/2000 = 0.0045 [0.0024, 0.0085] |
| null0721_n1000 | FS | 1.25/1.00 | 222 | 0.1110 [0.0980, 0.1255] | 0.1470 (0.0029) [k=222] | 0.9837 (0.0011) | 0.8530 (0.0029) [k=222] | 1.598 [k=222] | 0.684 [222/222] | 176/222 = 0.7928 | 9/222 = 0.0405 | 176/2000 = 0.0880 [0.0764, 0.1012] | 9/2000 = 0.0045 [0.0024, 0.0085] |
| null0721_n1000 | DINA | 0.90/0.80 | 1152 | 0.5760 [0.5542, 0.5975] | 0.1841 (0.0022) [k=1152] | 0.8940 (0.0024) | 0.8159 (0.0022) [k=1152] | 1.062 [k=1152] | 0.688 [1152/1152] | 75/1152 = 0.0651 | 7/1152 = 0.0061 | 75/2000 = 0.0375 [0.0300, 0.0468] | 7/2000 = 0.0035 [0.0017, 0.0072] |
| null0721_n1000 | DINA | 1.25/1.00 | 112 | 0.0560 [0.0467, 0.0670] | 0.1325 (0.0029) [k=112] | 0.9926 (0.0007) | 0.8675 (0.0029) [k=112] | 1.308 [k=112] | 0.676 [112/112] | 13/112 = 0.1161 | 1/112 = 0.0089 | 13/2000 = 0.0065 [0.0038, 0.0111] | 1/2000 = 0.0005 [0.0001, 0.0028] |
| null0721_n1000 | GRF | 0.90/0.80 | 1886 | 0.9430 [0.9320, 0.9523] | 0.2154 (0.0020) [k=1886] | 0.7969 (0.0022) | 0.7846 (0.0020) [k=1886] | 1.063 [k=1886] | 0.685 [1886/1886] | 97/1886 = 0.0514 | 11/1886 = 0.0058 | 97/2000 = 0.0485 [0.0399, 0.0588] | 11/2000 = 0.0055 [0.0031, 0.0098] |
| null0721_n1000 | GRF | 1.25/1.00 | 1028 | 0.5140 [0.4921, 0.5359] | 0.1418 (0.0012) [k=1028] | 0.9271 (0.0017) | 0.8582 (0.0012) [k=1028] | 1.329 [k=1028] | 0.681 [1028/1028] | 103/1028 = 0.1002 | 11/1028 = 0.0107 | 103/2000 = 0.0515 [0.0426, 0.0621] | 11/2000 = 0.0055 [0.0031, 0.0098] |
| null0657_n1500 | FS | 0.90/0.80 | 337 | 0.1685 [0.1527, 0.1855] | 0.1600 (0.0031) [k=337] | 0.9730 (0.0014) | 0.8400 (0.0031) [k=337] | 1.144 [k=337] | 0.616 [337/337] | 6/337 = 0.0178 | 0/337 = 0.0000 | 6/2000 = 0.0030 [0.0014, 0.0065] | 0/2000 = 0.0000 [0.0000, 0.0019] |
| null0657_n1500 | FS | 1.25/1.00 | 14 | 0.0070 [0.0042, 0.0117] | 0.1287 (0.0072) [k=14] | 0.9991 (0.0002) | 0.8713 (0.0072) [k=14] | 1.579 [k=14] | 0.603 [14/14] | 12/14 = 0.8571 | 0/14 = 0.0000 | 12/2000 = 0.0060 [0.0034, 0.0105] | 0/2000 = 0.0000 [0.0000, 0.0019] |
| null0657_n1500 | DINA | 0.90/0.80 | 278 | 0.1390 [0.1245, 0.1549] | 0.1498 (0.0026) [k=278] | 0.9792 (0.0012) | 0.8502 (0.0026) [k=278] | 0.966 [k=278] | 0.612 [278/278] | 1/278 = 0.0036 | 0/278 = 0.0000 | 1/2000 = 0.0005 [0.0001, 0.0028] | 0/2000 = 0.0000 [0.0000, 0.0019] |
| null0657_n1500 | DINA | 1.25/1.00 | 7 | 0.0035 [0.0017, 0.0072] | 0.1322 (0.0089) [k=7] | 0.9995 (0.0002) | 0.8678 (0.0089) [k=7] | 1.270 [k=7] | 0.596 [7/7] | 0/7 = 0.0000 | 0/7 = 0.0000 | 0/2000 = 0.0000 [0.0000, 0.0019] | 0/2000 = 0.0000 [0.0000, 0.0019] |
| null0657_n1500 | GRF | 0.90/0.80 | 1129 | 0.5645 [0.5427, 0.5861] | 0.1821 (0.0019) [k=1129] | 0.8972 (0.0023) | 0.8179 (0.0019) [k=1129] | 0.946 [k=1129] | 0.613 [1129/1129] | 6/1129 = 0.0053 | 1/1129 = 0.0009 | 6/2000 = 0.0030 [0.0014, 0.0065] | 1/2000 = 0.0005 [0.0001, 0.0028] |
| null0657_n1500 | GRF | 1.25/1.00 | 138 | 0.0690 [0.0587, 0.0810] | 0.1230 (0.0019) [k=138] | 0.9915 (0.0007) | 0.8770 (0.0019) [k=138] | 1.317 [k=138] | 0.601 [138/138] | 10/138 = 0.0725 | 1/138 = 0.0072 | 10/2000 = 0.0050 [0.0027, 0.0092] | 1/2000 = 0.0005 [0.0001, 0.0028] |
| null0721_n1500 | FS | 0.90/0.80 | 969 | 0.4845 [0.4626, 0.5064] | 0.1982 (0.0028) [k=969] | 0.9040 (0.0026) | 0.8018 (0.0028) [k=969] | 1.107 [k=969] | 0.688 [969/969] | 17/969 = 0.0175 | 1/969 = 0.0010 | 17/2000 = 0.0085 [0.0053, 0.0136] | 1/2000 = 0.0005 [0.0001, 0.0028] |
| null0721_n1500 | FS | 1.25/1.00 | 65 | 0.0325 [0.0256, 0.0412] | 0.1301 (0.0038) [k=65] | 0.9958 (0.0005) | 0.8699 (0.0038) [k=65] | 1.500 [k=65] | 0.684 [65/65] | 58/65 = 0.8923 | 1/65 = 0.0154 | 58/2000 = 0.0290 [0.0225, 0.0373] | 1/2000 = 0.0005 [0.0001, 0.0028] |
| null0721_n1500 | DINA | 0.90/0.80 | 748 | 0.3740 [0.3531, 0.3954] | 0.1731 (0.0023) [k=748] | 0.9353 (0.0021) | 0.8269 (0.0023) [k=748] | 0.982 [k=748] | 0.687 [748/748] | 14/748 = 0.0187 | 1/748 = 0.0013 | 14/2000 = 0.0070 [0.0042, 0.0117] | 1/2000 = 0.0005 [0.0001, 0.0028] |
| null0721_n1500 | DINA | 1.25/1.00 | 27 | 0.0135 [0.0093, 0.0196] | 0.1238 (0.0046) [k=27] | 0.9983 (0.0003) | 0.8762 (0.0046) [k=27] | 1.280 [k=27] | 0.673 [27/27] | 2/27 = 0.0741 | 0/27 = 0.0000 | 2/2000 = 0.0010 [0.0003, 0.0036] | 0/2000 = 0.0000 [0.0000, 0.0019] |
| null0721_n1500 | GRF | 0.90/0.80 | 1689 | 0.8445 [0.8280, 0.8597] | 0.2281 (0.0022) [k=1689] | 0.8074 (0.0026) | 0.7719 (0.0022) [k=1689] | 0.957 [k=1689] | 0.687 [1689/1689] | 29/1689 = 0.0172 | 2/1689 = 0.0012 | 29/2000 = 0.0145 [0.0101, 0.0207] | 2/2000 = 0.0010 [0.0003, 0.0036] |
| null0721_n1500 | GRF | 1.25/1.00 | 380 | 0.1900 [0.1734, 0.2078] | 0.1288 (0.0013) [k=380] | 0.9755 (0.0012) | 0.8712 (0.0013) [k=380] | 1.301 [k=380] | 0.676 [380/380] | 39/380 = 0.1026 | 2/380 = 0.0053 | 39/2000 = 0.0195 [0.0143, 0.0265] | 2/2000 = 0.0010 [0.0003, 0.0036] |

### Table 2 — FS: the candidate family, the consistency screen, and max_g T_g against 1.645

`clearing the floor` is at the screen's own c1, so the family `max_g T_g` is taken over differs by screen.

| cell | screen | enumerated Q1/med/Q3 | clearing the floor Q1/med/Q3 | empty floor family | consistency-qualifying Q1/med/Q3 [k] | floor > 0 but no declaration [Wilson] | max_g T_g Q1/med/Q3 | 90% / 95% / 99% | share > 1.645 [Wilson] |
|---|---|---|---|---|---|---|---|---|---|
| null0657_n500 | 0.90/0.80 | 1711 / 1711 / 1830 | 16 / 43 / 89 | 31/2000 | 1 / 3 / 8 [k=894] | 1075/1969 = 0.5460 [0.5239, 0.5678] | 0.483 / 0.833 / 1.225 | 1.577 / 1.807 / 2.176 | 166/1969 = 0.0843 [0.0728, 0.0974] |
| null0657_n500 | 1.25/1.00 | 1711 / 1711 / 1830 | 0 / 1 / 5 | 745/2000 | 1 / 2 / 3 [k=190] | 1065/1255 = 0.8486 [0.8277, 0.8674] | 0.868 / 1.114 / 1.414 | 1.723 / 1.941 / 2.254 | 166/1255 = 0.1323 [0.1146, 0.1521] |
| null0721_n500 | 0.90/0.80 | 1711 / 1711 / 1830 | 44 / 100 / 185 | 7/2000 | 2 / 7 / 19 [k=1356] | 637/1993 = 0.3196 [0.2995, 0.3404] | 0.767 / 1.159 / 1.564 | 1.887 / 2.093 / 2.575 | 414/1993 = 0.2077 [0.1905, 0.2261] |
| null0721_n500 | 1.25/1.00 | 1711 / 1711 / 1830 | 1 / 5 / 15 | 384/2000 | 1 / 2 / 4 [k=465] | 1151/1616 = 0.7123 [0.6897, 0.7338] | 0.995 / 1.299 / 1.660 | 1.944 / 2.164 / 2.592 | 414/1616 = 0.2562 [0.2355, 0.2780] |
| null0657_n1000 | 0.90/0.80 | 1711 / 1711 / 1830 | 6 / 16 / 37 | 93/2000 | 1 / 2 / 5 [k=637] | 1270/1907 = 0.6660 [0.6445, 0.6868] | 0.139 / 0.522 / 0.877 | 1.238 / 1.465 / 1.975 | 56/1907 = 0.0294 [0.0227, 0.0379] |
| null0657_n1000 | 1.25/1.00 | 1711 / 1711 / 1830 | 0 / 0 / 1 | 1424/2000 | 1 / 1 / 2 [k=65] | 511/576 = 0.8872 [0.8587, 0.9105] | 0.894 / 1.093 / 1.332 | 1.637 / 1.839 / 2.321 | 56/576 = 0.0972 [0.0756, 0.1242] |
| null0721_n1000 | 0.90/0.80 | 1711 / 1711 / 1830 | 23 / 53 / 104 | 13/2000 | 2 / 5 / 12 [k=1234] | 753/1987 = 0.3790 [0.3579, 0.4005] | 0.488 / 0.874 / 1.253 | 1.637 / 1.840 / 2.270 | 188/1987 = 0.0946 [0.0825, 0.1083] |
| null0721_n1000 | 1.25/1.00 | 1711 / 1711 / 1830 | 0 / 1 / 3 | 940/2000 | 1 / 1 / 3 [k=222] | 838/1060 = 0.7906 [0.7650, 0.8140] | 1.007 / 1.224 / 1.526 | 1.832 / 2.020 / 2.495 | 188/1060 = 0.1774 [0.1555, 0.2015] |
| null0657_n1500 | 0.90/0.80 | 1711 / 1711 / 1830 | 1 / 5 / 14 | 330/2000 | 1 / 2 / 3 [k=337] | 1333/1670 = 0.7982 [0.7783, 0.8168] | -0.137 / 0.146 / 0.489 | 0.831 / 1.037 / 1.579 | 13/1670 = 0.0078 [0.0046, 0.0133] |
| null0657_n1500 | 1.25/1.00 | 1711 / 1711 / 1830 | 0 / 0 / 0 | 1871/2000 | 1 / 1 / 2 [k=14] | 115/129 = 0.8915 [0.8261, 0.9342] | 0.940 / 1.111 / 1.373 | 1.628 / 1.865 / 2.391 | 13/129 = 0.1008 [0.0598, 0.1648] |
| null0721_n1500 | 0.90/0.80 | 1711 / 1711 / 1830 | 11 / 26 / 53 | 51/2000 | 1 / 3 / 8 [k=969] | 980/1949 = 0.5028 [0.4806, 0.5250] | 0.172 / 0.511 / 0.892 | 1.237 / 1.478 / 1.915 | 58/1949 = 0.0298 [0.0231, 0.0383] |
| null0721_n1500 | 1.25/1.00 | 1711 / 1711 / 1830 | 0 / 0 / 0 | 1612/2000 | 1 / 1 / 2 [k=65] | 323/388 = 0.8325 [0.7921, 0.8663] | 1.059 / 1.234 / 1.479 | 1.721 / 1.917 / 2.312 | 58/388 = 0.1495 [0.1174, 0.1884] |

### Table 3 — per cell: the realized whole-trial ITT HR (`itt_est`), and the cross-machine draw check

`itt_est` exists on the nullc125 bundles only (Gate C: identical across the three identifiers on all 2,000 rows). nullid has no `itt_est`, but its `or_Hc_est` is the same Cox fit of treatment alone on the whole trial, written on its declaring rows; the last column compares the two on the rows where nullid's FS declared.

| cell | target marginal HR | itt_est Q1 / median / Q3 | mean (MC SE) | nullid or_Hc_est vs nullc125 itt_est (FS rows) |
|---|---|---|---|---|
| null0657_n500 | 0.657 | 0.5713 / 0.6287 / 0.6877 | 0.6321 (0.0019) | 894 rows; max |diff| 7.994e-15; identical FALSE |
| null0721_n500 | 0.721 | 0.6365 / 0.6987 / 0.7602 | 0.7017 (0.0020) | 1356 rows; max |diff| 3.331e-16; identical FALSE |
| null0657_n1000 | 0.657 | 0.5922 / 0.6282 / 0.6687 | 0.6305 (0.0013) | 637 rows; max |diff| 5.551e-16; identical FALSE |
| null0721_n1000 | 0.721 | 0.6566 / 0.6966 / 0.7387 | 0.6994 (0.0014) | 1234 rows; max |diff| 2.220e-16; identical FALSE |
| null0657_n1500 | 0.657 | 0.5967 / 0.6266 / 0.6593 | 0.6286 (0.0010) | 337 rows; max |diff| 8.882e-16; identical FALSE |
| null0721_n1500 | 0.721 | 0.6611 / 0.6960 / 0.7305 | 0.6974 (0.0011) | 969 rows; max |diff| 1.110e-16; identical FALSE |

### Table 4 — the composition of Ĥ, per identifier and screen


**FS** — share of declaring replicates whose rule contains the term, pooled over the six cells (declaring replicates pooled: 0.90/0.80 5427, 1.25/1.00 1021; one-factor rules: 157/5427 = 0.0289 [0.0248, 0.0337] / 17/1021 = 0.0167 [0.0104, 0.0265]):

| term | 0.90/0.80 replicates (share) | 1.25/1.00 replicates (share) |
|---|---|---|
| `NOT er <=` | 1520 (0.2801) | 313 (0.3066) |
| `NOT pgr <=` | 1163 (0.2143) | 261 (0.2556) |
| `er <=` | 1377 (0.2537) | 244 (0.2390) |
| `size <=` | 1014 (0.1868) | 195 (0.1910) |
| `NOT size <=` | 986 (0.1817) | 195 (0.1910) |
| `nodes <=` | 1039 (0.1915) | 191 (0.1871) |
| `age <=` | 850 (0.1566) | 156 (0.1528) |
| `NOT age <=` | 816 (0.1504) | 151 (0.1479) |
| `NOT nodes <=` | 628 (0.1157) | 104 (0.1019) |
| `pgr <=` | 614 (0.1131) | 95 (0.0930) |
| `NOT meno (indicator)` | 305 (0.0562) | 53 (0.0519) |
| `meno (indicator)` | 196 (0.0361) | 31 (0.0304) |
| `grade (indicator)` | 73 (0.0135) | 21 (0.0206) |
| `NOT grade (indicator)` | 116 (0.0214) | 15 (0.0147) |

**DINA** — share of declaring replicates whose rule contains the term, pooled over the six cells (declaring replicates pooled: 0.90/0.80 5367, 1.25/1.00 842; one-factor rules: 983/5367 = 0.1832 [0.1730, 0.1937] / 259/842 = 0.3076 [0.2774, 0.3396]):

| term | 0.90/0.80 replicates (share) | 1.25/1.00 replicates (share) |
|---|---|---|
| `pgr >=` | 1564 (0.2914) | 433 (0.5143) |
| `er >=` | 997 (0.1858) | 235 (0.2791) |
| `age >=` | 872 (0.1625) | 121 (0.1437) |
| `nodes <=` | 857 (0.1597) | 107 (0.1271) |
| `size <=` | 712 (0.1327) | 98 (0.1164) |
| `size >=` | 830 (0.1546) | 89 (0.1057) |
| `age <=` | 827 (0.1541) | 84 (0.0998) |
| `nodes >=` | 538 (0.1002) | 57 (0.0677) |
| `grade >=` | 459 (0.0855) | 49 (0.0582) |
| `pgr <=` | 617 (0.1150) | 37 (0.0439) |
| `grade <=` | 308 (0.0574) | 37 (0.0439) |
| `er <=` | 574 (0.1069) | 29 (0.0344) |
| `meno <=` | 356 (0.0663) | 27 (0.0321) |
| `meno >=` | 240 (0.0447) | 22 (0.0261) |

**GRF** — share of declaring replicates whose rule contains the term, pooled over the six cells (declaring replicates pooled: 0.90/0.80 10219, 1.25/1.00 5161; one-factor rules: 657/10219 = 0.0643 [0.0597, 0.0692] / 202/5161 = 0.0391 [0.0342, 0.0448]):

| term | 0.90/0.80 replicates (share) | 1.25/1.00 replicates (share) |
|---|---|---|
| `nodes <=` | 2927 (0.2864) | 1831 (0.3548) |
| `pgr >` | 2418 (0.2366) | 1479 (0.2866) |
| `er >` | 1892 (0.1851) | 1082 (0.2096) |
| `size <=` | 2006 (0.1963) | 1047 (0.2029) |
| `age <=` | 1889 (0.1849) | 933 (0.1808) |
| `size >` | 1729 (0.1692) | 769 (0.1490) |
| `age >` | 1639 (0.1604) | 713 (0.1382) |
| `er <=` | 1439 (0.1408) | 652 (0.1263) |
| `pgr <=` | 1036 (0.1014) | 409 (0.0792) |
| `nodes >` | 880 (0.0861) | 332 (0.0643) |
| `meno <=` | 671 (0.0657) | 313 (0.0606) |
| `grade <=` | 503 (0.0492) | 217 (0.0420) |
| `grade >` | 406 (0.0397) | 204 (0.0395) |
| `meno >` | 346 (0.0339) | 139 (0.0269) |

### Table 5 — wall per cell, and the machine, R version and build of each campaign

| cell | 0.90/0.80 wall s (FS / DINA / GRF) | 0.90/0.80 cell s | 1.25/1.00 wall s (FS / DINA / GRF) | 1.25/1.00 cell s |
|---|---|---|---|---|
| null0657_n500 | 316 / 167 / 508 | 991 | 186 / 105 / 418 | 709 |
| null0721_n500 | 382 / 238 / 519 | 1139 | 216 / 131 / 428 | 775 |
| null0657_n1000 | 320 / 99 / 637 | 1056 | 190 / 72 / 353 | 615 |
| null0721_n1000 | 399 / 165 / 676 | 1240 | 200 / 75 / 400 | 675 |
| null0657_n1500 | 329 / 73 / 750 | 1152 | 196 / 68 / 323 | 587 |
| null0721_n1500 | 394 / 116 / 821 | 1331 | 196 / 65 / 380 | 641 |
| **total (18 renders)** | | **6909 s = 1.919 h** | | **4002 s = 1.112 h** |

| campaign | screen | hostnames | R | forestsearch | workers | bundles |
|---|---|---|---|---|---|---|
| nullid | 0.90/0.80 | Mac-Studio-3.local | 4.5.2 | 0.3.5.9000 | 12 | 18 |
| nullc125 | 1.25/1.00 | pop-os | 4.6.1 | 0.3.5.9000 | 64 | 18 |

Wall per campaign is **not a screen comparison**: `nullid` ran on the Mac at 12 workers, `nullc125`
on pop-os at 64. The strict screen is cheaper per render on every engine, but the two totals
confound machine, worker count and screen.

**The draws are the same across the two campaigns, not only within one.** Table 3's last column
compares `nullid`'s `or_Hc_est` — under the structural null, the same whole-trial Cox fit — with
`nullc125`'s `itt_est` on every row where `nullid`'s FS declared (337–1,356 rows per cell): they agree
to **≤ 8.0e-15**, not bit for bit (Mac R 4.5.2 against Linux R 4.6.1). So the side-by-side rows of
Table 1 are paired on the same 2,000 trials per cell, and differences between the screens are
differences of the screen, not of the draw.

---

## 5. Findings

- **The stricter screen cuts the false-declaration rate on every identifier in every cell, but
  does not bring it to any nominal level.** FS 0.1685–0.6780 → **0.0070–0.2325**; DINA
  0.1390–0.7400 → **0.0035–0.2130**; GRF 0.5645–0.9845 → **0.0690–0.8350**. The worst cell is
  still `null0721_n500`: FS 0.2325 [0.2145, 0.2515], DINA 0.2130 [0.1956, 0.2315], GRF 0.8350
  [0.8181, 0.8506]. As at 0.90 / 0.80, the rate falls with n and is higher at the weaker uniform
  benefit.
- **GRF remains the outlier.** It declares on 68.5% and 83.5% of replicates at n 500 under
  c1 = 1.25; its admission floor follows c1, but the frontier still finds a candidate that clears it
  far more often than FS's consistency screen or DINA's proposal floor lets through.
- **Where the floor bites (FS, Table 2).** At c1 = 1.25 the median number of candidates clearing the
  floor is 0–5 (against 5–100 at 0.90), and the family is **empty** on 384–1,871 of 2,000 replicates.
  Where it is non-empty, the consistency screen at c2 = 1.00 declines it on **71.2%–89.2%** of
  replicates (against 32.0%–79.8%).
- **Fewer declarations, but not fewer unadjusted "harm" claims — for FS the unconditional rate of
  a declaration whose unadjusted one-sided 95% lower bound reaches HR 1.00 *rises* in every cell**:
  0.0465 → 0.0765, 0.1125 → 0.1910, 0.0130 → 0.0265, 0.0425 → 0.0880, 0.0030 → 0.0060,
  0.0085 → 0.0290 (cell order of Table 1). The mechanism, read row by row at `null0657_n500`: every
  one of the 190 strict-screen declarations is a replicate that also declared at 0.90; all 93 of the
  0.90-screen claims are kept, and 60 more appear. Removing the low-effect candidates from the family
  moves the `effMaxSG` band — "the largest subgroup within ε of the maximal effect" — onto smaller,
  higher-effect regions (80 of the 190 shared declarations change their label; median |Ĥ| 99 → 89.5).
  Conditionally, **79%–89% of FS's strict-screen declarations carry an unadjusted lower bound ≥ 1.00**.
  The unadjusted bound is not selection-adjusted, and this is what selection does to it.
- **The same quantity does not rise on DINA, and on GRF it rises only at n ≥ 1000.** DINA falls in
  every cell (e.g. 0.1345 → 0.0530 at `null0721_n500`). GRF is **unchanged exactly** at n 500
  (113 and 268 of 2,000 on both screens) and rises at n 1000 and 1500 (38 → 42, 97 → 103, 6 → 10,
  29 → 39).
- **At HR 1.25 the unconditional count is identical across screens for FS and GRF in every cell**
  (FS 7, 21, 4, 9, 0, 1; GRF 20, 49, 6, 11, 1, 2 of 2,000). A declaration whose unadjusted lower
  bound already reaches 1.25 survives the stricter screen, as it should — its estimate clears 1.25. On
  DINA the count falls (16 → 6, 46 → 21, …).
- **`max_g T_g` over the whole replicate set does not depend on the screen.** The number of
  replicates with `max_g T_g > 1.645` is identical on both screens in every cell (166, 414, 56, 188,
  13, 58 of 2,000): a candidate with T > 1.645 has an HR well above 1.25 and is in both floor
  families. The *conditional* share in Table 2 differs only because the denominator — replicates with
  a non-empty floor family — shrinks. Neither screen changes the finding of `nullid` §5 that a fixed
  1.645 is uncalibrated across this grid (0.0065–0.2070 of all replicates).
- **Size and specificity.** Mean |Ĥ|/n falls to 0.123–0.194 (from 0.150–0.228), so conditional
  specificity — exactly 1 − |Ĥ|/n here — rises to 0.806–0.877; unconditional specificity is
  0.850–0.9995, driven by the declaration rate.
- **The truth in Ĥ is the uniform benefit, whatever the estimate says.** Median true β(Ĥ) (HR scale,
  super-population marginal Cox HR in the selected region) is 0.596–0.618 at the 0.657 cells and
  0.673–0.688 at the 0.721 cells, on both screens, against median unadjusted within-region HRs of
  1.27–1.77 at c1 = 1.25.
- **Composition of Ĥ.** FS and GRF shift little between screens; FS's ER terms still lead because
  of how candidates are built (`nullid` §5), not because of signal. **DINA concentrates:** `pgr >=`
  appears in 51.4% of strict-screen declarations (29.1% at 0.90), and one-factor rules rise from
  18.3% to 30.8%.
- **The realized trial ITT HR sits below the super-population target** — median 0.627–0.629 at the
  0.657 cells, 0.696–0.699 at the 0.721 cells. The target is the marginal Cox HR on the uncensored,
  stacked potential outcomes; the trial fits one arm-mixed sample with 84-month administrative and
  random censoring. Recorded, not investigated.

## 6. What is not reported, and why (as `nullid` §6)

- **Sensitivity and PPV are undefined with an empty planted region** and are not reported.
  Sensitivity is NA on all 36,000 replicates (Gate A asserts it per run); PPV is 0 by construction on
  every declaring replicate and carries no information.
- **NPV is 1 by construction** on every declaring replicate; Gate A asserts it.
- **The identified-to-planted size ratio has no denominator** and is not reported.
- **No MR product of any kind exists in these bundles**; Gate A asserts every `mr_*` / `fld_*` /
  `fb_*` product column is NA and `mr_ok == 0`, per run.
- `n_family` stays NA and the maximum consistency rate over the full screened family is not
  recoverable (`nullid` §0.6); the family counts, `maxT` and `p_sel` are FS-only.
- Rates carry Wilson 95% intervals and means carry Monte Carlo standard errors; every conditional
  summary carries its `k`.

## 7. Files

| what | path |
|---|---|
| bundles (18) | `results/{fs,dina,grf}_effMaxSG_fb_mr_field_m1_h0{66,72}_knoise0_n{500,1000,1500}_null{657,721}_nb20_nomr_c125c100_nullc125_res_1_2000.rds` |
| smoke / inertness bundles (12) | `results/*_nullc125smoke_quickrun_res_1_20.rds`, `results/*_nullc125inert{pre,unset,expl}_quickrun_res_1_20.rds` |
| driver, cells, gates, smoke, inertness, findings | `scripts_dinamr/nullthr.sh`, `nullthr.cells`, `nullthr_smoke.cells`, `nullthr_gateA.R`, `nullthr_gateC.R`, `nullthr_smoke.sh`, `nullthr_inertness.R`, `nullthr_findings.R` |
| logs | `scripts_dinamr/logs/nullc125*` (driver, per-render, Gate A, Gate C, inertness, smoke/projection, findings) |
| renders (not committed; not required by the task) | `nullc125_*.html`, `nullc125smoke_*.html`, `nullc125inert*_*.html` in the study directory |
