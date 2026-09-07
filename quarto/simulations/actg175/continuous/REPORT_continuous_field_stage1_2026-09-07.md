# REPORT — Continuous/MD field port, Stage 1 (port, identities, smoke, projection)

Date: 2026-09-07. Machine: Mac Studio (Apple M4 Max, 14 physical cores = 10 P + 4 E, 36 GB). Branch: `feature/glm-extension-mac`. Package: forestsearch 0.3.5 installed from this tree (with the `scale` change below). Task: `dev/tasks/TASK_continuous_field_mac_2026-09-07.md`. Stage 0 record: `REPORT_continuous_field_stage0_2026-09-07.md` (Gate 0 PASS).

## 1a. The port and the R/ change

**R/ change (the one permitted; add-only).** `fs_sim_bias_coverage()` in `R/fs_bias_coverage.R` gains `scale = c("log", "identity")` (commit `2f118042`): `match.arg`; a working-scale transform `.tr` (log or identity) applied to the target and the estimates; the normal-based one-sided bound is the original `exp(log e ± z·se)` with its `e > 0` guard under `"log"` and `e ± z·se` with no positivity guard under `"identity"`. Output columns unchanged (`bias_log` holds the identity-scale bias under `"identity"`; documented in the Rd). `fs_plot_bias_coverage()` untouched (reads only `b`, `r`, coverages). Classification: adds code; no behaviour change under the default.
- **Default byte-identical:** `fs_sim_bias_coverage()` on the seven committed `s7`/`map1` combined bundles, both blocks × both sides × all five estimators (28 tables) — `identical()` to the pre-change output, 28 of 28.
- **14-point fixture:** unchanged from the adjudicated state (`REPORT_bias_coverage_display_2026-09-06.md`): all 56 observed quantities within 0.005 (max 0.0045), 27 of 28 references within 0.001, the known `harm 1.75, n=500, field` c2_pred row at 0.00104 (fixture rounding, not computation).
- **Suite:** new `tests/testthat/test-bias-coverage-scale.R` (26 assertions: default ≡ `"log"`; identity-scale hand computations for bias, SD, SE, b, r, both coverages and both references on both sides; non-positive estimates kept under `"identity"` and dropped under `"log"`); `test-mr-field-complement.R` still 38/38. `devtools::document()` changed only `man/fs_sim_bias_coverage.Rd`. Following the 09-06 precedents (display functions, complement field), no version bump and no NEWS entry.

**The template.** New document `quarto/simulations/actg175/continuous/sim_fs_maxeffCons_mr_field_md_template.qmd` (1,803 lines, 32 chunks), beside the twin; the committed batch documents and bundles are untouched. Every Stage 0 port-map row landed:

| addition | in the template |
|---|---|
| `FS_MD_*` knobs (`NSIMS START MODE N MD BUILD CAMPAIGN CI FIELD_COMPLEMENT IJ_RESIDUAL RESELECTION UNIFORM WINNER_ROWS FB FB_PATH JOIN_SKIP WORKERS KNOISE QUICKRUN SAVE_COMBINED`) with the twin's values as defaults; `FS_MD_MD=null` is the null cell; `FS_MD_MD=120` builds directly at the locked `k_inter` (the committed md120 route) | `setup-knobs` |
| stem `fs_maxeffCons_mr_field_md%02d_knoise%d_n%d_<campaign>` (token `mr_field`; `_mdnull_`; `_quickrun`); results dir `mr_md_harm/<stem>_d<draws>` (the twin's `_s<n_sims>` piece dropped so batches of any size pool) | `setup-knobs` |
| `.refuse_if_tracked()` before both saves | `setup-knobs`, `run-batch` |
| `mr_inference_args` = `ci_method` (default `"field"`), `draws`, `include_complement = TRUE`, `confirm_rule`, `field_uniform` (off), `field_complement` (default TRUE), `ij_residual` (`"two_term"`), `return_reselection` (TRUE) | `setup-knobs` → `record_replicate()` |
| recorder: the twin's 55 columns first, in its order and types; then `label`, `mr_harm_flag`, `fld_H_*` (23), uniform (8), `fld_Hc_*` (24), winner variants (12), `fld_joint_*` (9), `p_hat_H`, `p_top1..3`, `p_lab1..3` — 138 columns before `fs_attach_betaHhat()` adds three | `machinery` |
| `.safe_record()` firewall; `.ci_check()` interval invariant incl. field pairs; FB `"join"` with the 1e-12 naive-identity licence and enumerated skips; provenance meta; combine-mode poolability keys (`ci_method campaign_tag target_md_harm null_cell dgm_build field_complement ij_residual k_random_noise effect_threshold consistency_threshold` added) | `run-batch` |
| `mr-settings-readout`; counts with field/complement/FB tallies; four-target estimation tables with the MR (field) row; **Table-2 layout** (bias in MD and SD units, SD, SE, SE/SD, two-sided and exposed-side one-sided coverage with Wilson intervals, half-widths, margins) in the full and the compact format; **bound-location tables** against MD thresholds 0/10/20/30/40 (harm lower bound) and 0/10/20/30 (complement upper bound); field coverage with Wilson intervals against all four targets; the display on the identity scale for both blocks; retained-bias and draw-usage diagnostics; p̂(Ĥ) diagnostics; complement regime diagnostics (SD(β̃ᶜ)/naive SE, λ-SDᶜ/naive SE, corr(Λ*, Λ*ᶜ)); the joint pair; identification and timing (twin's) | summary layer |

Smoke render (`FS_MD_NSIMS=3 FS_MD_QUICKRUN=TRUE`, field + complement): exit 0, 25 s wall, every section present, no warnings; smoke outputs deleted.

## 1b. Identities

Eight 5-replicate renders of the template (`sim_id` 1–5; campaigns `idij` / `idfield`; 14–29 s each), compared with the committed bundles' rows 1–5 on every pre-existing column (the twin's 59 minus timings and `mr_msg`; 49 numeric columns; NA ≡ NA; relative tolerance 1e-8; rule strings as term **sets**):

| cell (committed pkg) | `ci_method = "ij"`: max rel diff, pre-existing numerics | `"field"`: same | rule sets / status / n_harm / ij_source / betaHhat_status | truth |
|---|---|---|---|---|
| md40 n = 500 (0.2.2) | 2.1e-12 | 2.1e-12 | identical | equal |
| md40 n = 700 (0.2.2) | 2.2e-13 | 2.2e-13 | identical | equal |
| md120 n = 500 (0.3.1) | 2.1e-12 | 2.1e-12 | identical | equal |
| null n = 500 (0.2.2) | 6.1e-13 | 6.1e-13 | identical as sets (sim 4: `{str2} & {cd80 <= 623}` committed vs `{cd80 <= 623} & {str2}` here — term order only; same members, same n_harm, same numerics) | equal |

Under `"field"`, on all four cells: the MR (IJ) columns equal the `"ij"` render's to 0.0 relative; field and complement-field columns finite on 20 of 20 detected replicates; **bound identities** (`lo1s = β̃ − q95`, `lo2s = β̃ − q975`, `hi2s = β̃ − q025`, `est2 = β̃ − mean Λ*`, `lo_se = est2 − 1.96·λ-SD`; complement `up1s = β̃ᶜ − q05`, `lo1s = β̃ᶜ − q95`, `lo2s/hi2s`, `est2`) hold to **0.0** absolute (≤ 1e-12 required); γ ∈ {0.025, 0.026} (Bonferroni or one grid step above; joint probability ≥ 0.95); p̂(Ĥ) finite in [0, 1]; outer draws used ≥ 999 of 1000; complement fits ≈ 520–590 per replicate.

**Re-selection map (replicate 1, md40 n = 500, `"field"`).** The family rebuilt exactly as `forestsearch_main.R:3351–3365` (Z = `dummy()` over the 36 screened cut factors, ≤ 2 conjunctions, `n.min = 60`) has 1,842 candidates = the gate's `n_family`; the gate's `selected_label` (`q25.0 & q26.0` = `!{cd40 <= 415} & !{cd80 <= 1022}`) is in the family and its members equal the observed Ĥ as a set (|Ĥ| = 69); the assembled unperturbed β̂ at that candidate equals the gate's naive estimate (103.1366726933) to 1e-10; the gate's admission set on β̂ (effect floor 30, consistency floor 10 at p* = 0.90 → 294 passers) and its `maxeff` rule return the same label: **reproduced**. p̂(Ĥ) = 0.147 (top-3: 0.147 / 0.065 / 0.062).

**Save guard.** `.refuse_if_tracked()` stops on the committed md40 bundle path and allows an untracked path under the new stem (verified directly).

**Cross-vintage note.** Bundles built with 0.2.2 and 0.3.1 reproduce under 0.3.5 to ≤ 2.1e-12 on 20 replicates across the four cells; the only difference seen is the term order inside one rule string. Rule strings are therefore compared as term sets in Gate 2.

## 1c. Worker calibration and projection (measured under load)

| run (`"field"`, complement on) | workers | reps | run-loop wall | throughput | fit+MR+field per replicate: mean / max | field / complement per replicate |
|---|---|---|---|---|---|---|
| md40 n = 500 | 13 | 26 | 40.5 s | 1.56 s/rep | 14.8 / 19.5 s | 10.8 / 0.36 s |
| md40 n = 500 | 10 | 20 | 32.9 s | 1.65 s/rep | 12.7 / 16.2 s | 9.3 / 0.30 s |
| md40 n = 700 | 13 | 26 | 46.0 s | 1.77 s/rep | 17.2 / 22.7 s | 12.3 / 0.37 s |

Unloaded single replicate: 2.4 s (`"ij"`), 9.4 s (`"field"` + complement). Under 13-worker load a replicate costs ≈ 1.3× the unloaded time (the E-cores and memory bandwidth); 13 workers still out-throughput 10 (1.56 vs 1.65 s/rep), so **M-4 resolves to 13 workers**. The two-round calibration runs carry ≈ 10 s of worker spawn; for a 1,000-replicate batch (77 rounds) the busiest-worker model is ceil(N/13) × per-replicate mean:

| cell | per-replicate (13 workers) | 1,000 replicates | 2,000 replicates |
|---|---|---|---|
| md40 n = 500 | 14.8 s | ≈ 19 min (model) → ≈ 21 min with overhead | ≈ 42 min |
| md40 n = 700 | 17.2 s | ≈ 22 min → ≈ 24 min | ≈ 48 min |
| md120 n = 500 | ≈ 15 s (identity renders: 13.2 s unloaded-ish) | ≈ 21 min | ≈ 42 min |
| null n = 500 | ≈ 14 s (10.7 s) | ≈ 20 min | ≈ 40 min |
| **all four cells** | | **≈ 1.5 h** | **≈ 3 h** |

Plus ≈ 10–20 s of summary-layer rendering per document and the DGM build (< 1 s). Memory was not the binding constraint at 13 workers (all calibration runs completed; no swapping observed).

## Gate 1

- Identities: **PASS** on all four cells (≤ 2.1e-12 relative; memberships equal as sets; field invariants exact; IJ columns bit-identical between `"ij"` and `"field"`; re-selection map reproduced; guard verified).
- Projection reported: 13 workers; ≈ 1.5 h for the committed count (1,000 per cell), ≈ 3 h for 2,000 per cell.

**Gate 1: PASS. Stopping for M-5 (compute go / pre-authorization with a ceiling).**

Proposal for M-5 / M-2: 2,000 replicates per cell (the task's M-2 rule: the committed count is 1,000 and cost allows), run as two seed-disjoint batches per cell (`sim_id` 1–1000, then 1001–2000; the first batch carries the pairing proof against the committed bundle on all 1,000 anchored replicates) and combined; campaign `contfieldmac` (the stem's campaign tag is alphanumeric-only, so the task's `cont_field_mac` is written without underscores); `FS_MD_FB=join` on md40 n = 500 batch 1 (the committed FB bundle covers `sim_id` 1–100); ceiling ≈ 4 h wall, hard timeout 1.5 h per render.

## Files

- Committed with this record: the template; the R/ change (`2f118042`). Identity and calibration bundles stay on disk, untracked, under `mr_md_harm/*_id{ij,field}_d5000/` (8 × 5 replicates) and `*_cal{13n500,10n500,13n700}_d5000/` (26 / 20 / 26 replicates); scratch scripts in the session scratchpad.
