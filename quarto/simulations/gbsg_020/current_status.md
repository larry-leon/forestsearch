# current_status — `quarto/simulations/gbsg_020`

- **Pin:** `ca9b0a62` on `feature/glm-extension` — HEAD at the time this file was committed; the closeout commit that adds this file is its child.
- **Updated:** 2026-09-12
- **Purpose:** a catalog of what has been run in this directory and where the payloads are, so a chat or workstream on another machine can be brought up to speed by attaching this one file. It points at authoritative files; it does not restate their numbers.
- **Maintenance:** regenerated as the closeout step of every task that touches this directory. The pin above must equal HEAD at commit time.

---

## 1. The DGM, in one block

- GBSG-based survival simulation. One template: `sim_fs_maxeffCons_fb_mr_field_m1_template.qmd`.
- Cells are (prevalence × target HR × n): prevalence 12.4% (`FS_S7_Z1Q` unset) or 31% (`FS_S7_Z1Q=0.60`); HR 1.50, 1.75 or 1.00; n = 500, 1000, 1500. 2,000 replicates per cell unless stated.
- **HR 1.00 is not a global null.** `FS_S7_HR` calibrates `k_inter` to a target Cox HR inside the planted region; the region rule depends on `FS_S7_Z1Q` alone; `dgm_model <- "alt"` is a literal, so the harness's global-null path is unreachable. At HR 1.00 the planted region carries HR 1.00 against a benefiting complement (HR 0.657 at 12.4%, 0.721 at 31%). Describe it as **differentially null against a benefiting complement**, never with harm vocabulary.
- Seeds are `8316951 + sim_id` throughout, so cells at matched coordinates share DGM draws across campaigns. Verified per cell: `n_true` identical on all 2,000 rows, truth agreeing to ≤ 8.9e-15 (cross-machine BLAS; the largest observed is 8.882e-15, on the 31% cells).

## 2. Campaigns

### 2.1 FS comparators (the reference grid)

| Campaign | Prevalence | Criterion | Cells |
|---|---|---|---|
| `p12ext` | 12.4% | `maxeffCons` ε 0.10 | HR 1.50 (all n); HR 1.00 at n 1000, 1500 |
| `tier2` | 12.4% | `maxeffCons` ε 0.10 | HR 1.75 (all n); HR 1.00 at n 500 |
| `e1stud` | 31% | `effMaxSG` ε 0.20 | HR 1.50 and 1.75 at n 500 |
| `cert20` | 31% | `effMaxSG` ε 0.20 | all others |

- All 18 grid cells covered; all 18 designated comparator bundles verified present. The table names the **designated** comparator — the one `gate2G.R` and `fs_extraction.R` resolve to. Earlier FS campaigns (`map1`, `s7`, and others) also hold bundles at some of these cells; those are not the comparator and must not be substituted.
- **Criterion is matched to DINA and GRF at 31% and not at 12.4%.** Any 12.4% cross-identifier gap carries a criterion confound on top of identifier, family construction and detection set.
- FS one-sided products, recomputed from the committed bundles across **all 12 harm cells**, both prevalences:
  - **field lower on β(Ĥ): 0.941–0.975.** The certification record's 0.944–0.974 is the same quantity on a **smaller cell set** — it predates `p12ext`, which supplies the three 12.4% HR 1.50 cells, and 12.4% HR 1.50 n 1000 (0.9410) is the only harm cell below 0.944. Both are right for their set; quote the cell set with the range.
  - **field-s upper on β(Ĥᶜ): 0.9125–0.9605.** Identical to the record's "0.912–0.960" — the endpoints are the same numbers under a different rounding convention, not a disagreement.
  - **Bonferroni joint: unscaled `joint` 0.932–0.964; studentized `joint_s` 0.940–0.964.** The record's "0.939–0.963" is **`joint_s`**, whose 9-cell range is 0.9395–0.9640. `joint` and `joint_s` are different constructions and must be named when quoted.
- Classification and bound-location extraction: `REPORT_fs_extraction_2026-09-11.md` (all 18 cells; reading of committed bundles, no re-run).

### 2.2 `dinamr` — DINA, complete

- Engine `dina`, `effMaxSG`, ε 0.20, effect floor log(0.90), no consistency term. 18 of 18 cells, 2,000 replicates.
- Reports: `REPORT_dinamr_*`, `REPORT_dinamr_blockC_2026-09-11.md`. Summary: `summary_dinamr.qmd` / `.html` (18 cells). Review: `REVIEW_dinamr_blockC_2026-09-11.md`.
- DINA's `n_family` is the **surface-proposed candidate set** — small, volatile and shrinking with n. CV 0.296–1.458 across all 18 cells; the 1.008–1.458 band quoted in `REPORT_dinamr_blockC` is the **six HR 1.00 cells only**, where it is most extreme.
- **Detection falls with n**, sharply at 12.4%: 0.875 → 0.839 → 0.799 on HR 1.50 harm, 0.908 → 0.907 → 0.898 on HR 1.75 harm, and 0.714 → 0.525 → 0.344 at the null. Overall range across the 18 cells is 0.344–1.000. Every DINA n-trend conditions on a shrinking detection set; say so when reading one.
- On DINA, `admitted_n` is not recorded (the column postdates this campaign) and `n_family` with p̂ are strongly confounded (ρ = −0.47), so the joint count table is needed to read either stratification.

### 2.3 `grfmr` — GRF, partial

- Engine `grf`, `effMaxSG`, ε 0.20, `grf_selection = "frontier"`, `grf_select_statistic = "effect"`, `dmin.grf = 0.0`. **10 of 12 harm cells**, 2,000 replicates.
- **Deferred, not dropped:** 31% HR 1.75 at n 1000 and n 1500. **Never run:** all six HR 1.00 null cells. Together ≈ 6.7 h at the realized rate (1.9 h for the two deferred, ≈ 4.9 h for the six null).
- Reports: `REPORT_grfmr_2026-09-11.md`, `TABLES_grfmr_percell_2026-09-12.md`. Summary: `summary_grfmr.qmd` / `.html` (10 cells, 6 guarded skips). Review: `REVIEW_grfmr_2026-09-12.md`.
- **GRF's `n_family` is the outcome-independent enumerated pool**, not the qualified set. It is identical quantile-for-quantile across both prevalences and both hazard ratios at each n. The outcome-dependent quantity is **`admitted_n`**, which spans 115–449 across those same cells with ρ(`admitted_n`, `n_family`) ≤ +0.042.
  - **Stratify GRF on `admitted_n`, never on `n_family`.**
  - `admitted_n` and p̂ carry independent information on GRF (joint counts near-uniform) — the opposite of DINA.
- Detection is 1.0000 at all four 31% cells and 0.9975–1.0000 at 12.4%, flat in n, so no GRF n-trend carries a detection-conditioning caveat.
- Source trace (Gate 0a, 13 sites): GRF's candidates are enumerated from quantiles of X subject to `n_min` with DR scores entering only the effect column; MR re-evaluates the full enumerated pool, not the forest-qualified subset, with qualification re-applied per draw. Recorded as a fact about the construction. It sits against handoff §3's description of GRF; that discrepancy is noted, unresolved, and blocks nothing.
- Two floors on two scales: `dmin.grf` is a DR-score pre-filter in RMST units; the binding effect-scale floor on the re-selection path is `hr.threshold = 0.90`, the same floor DINA carries. GRF is not "unfloored" relative to DINA.

### 2.4 `grfprobe` — GRF cost and mechanism

- Five 36-replicate probes. Cost 13.5–19.8 s per replicate, rising with n and prevalence, **flat in family size** (|ρ| ≤ 0.171). Frontier band empty on 0 of 180; `admitted_n` never 0 across 20,000 campaign replicates either.

## 3. Payload inventory

| path pattern (first match wins) | tracked/disk | total | largest single file | what it is |
|---|---|---|---|---|
| `results/*dinamr*.rds` | 66/66 | 50.54 MB | `dina_effMaxSG_fb_mr_field_m1_h150_knoise0_n500_z1q60_nb20_dinamr_combined_1_2000.rds` 1.54 MB | `dinamr` per-replicate bundles + metas (batch and combined) |
| `results/*grfmr*.rds` | 30/30 | 30.78 MB | `grf_effMaxSG_fb_mr_field_m1_h150_knoise0_n500_nb20_grfmr_combined_1_2000.rds` 1.54 MB | `grfmr` per-replicate bundles + metas (batch and combined) |
| `results/*grfprobe*.rds` | 5/5 | 155 KB | `grf_effMaxSG_fb_mr_field_m1_h150_knoise0_n500_nb20_grfprobe_res_1_36.rds` 31 KB | `grfprobe` cost probes, 36 replicates each |
| `results/fs_*.rds` | 315/315 | 149.12 MB | `fs_maxeffCons_fb_mr_field_m1_h150_knoise0_n1000_p12ext_combined_1_2000.rds` 1.51 MB | FS bundles — the comparator grid plus every earlier FS campaign |
| `results/*.rds` | 14/14 | 960 KB | `grf_eff_fb_mr_m1_h10_knoise0_n500_combined_1_500.rds` 150 KB | other bundles in `results/` |
| `mr_sweep/**` | 210/210 | 15.19 MB | `grf_mr_n500_res.rds` 130 KB | `mr_sweep/` — an earlier seed-table sweep, superseded, kept for provenance |
| `scripts_dinamr/logs/*` | 40/40 | 122 KB | `grfmr_A124_h150_n1000_batch_1001.log` 3 KB | per-render and driver logs — `WALL_SECONDS` / `CELL DONE wall=` |
| `scripts_dinamr/*.R` | 20/20 | 148 KB | `gate2G.R` 16 KB | drivers, checkers, projections, extractions (R) |
| `scripts_dinamr/*.sh` | 8/8 | 12 KB | `grfprobe.sh` 3 KB | render/campaign drivers and the closeout checker (shell) |
| `scripts_dinamr/*.py` | 2/2 | 13 KB | `transplant_grfmr.py` 12 KB | transplant / chunk-diff helpers (Python) |
| `scripts_dinamr/*.cells` | 6/6 | 1 KB | `grfmr.cells` 0 KB | cell lists, one line per cell |
| `scripts_dinamr/*.rds` | 5/5 | 44 KB | `grfmr_tables.rds` 35 KB | saved derived objects (projections, extracted tables) |
| `scripts_dinamr/*.md` | 1/1 | 7 KB | `README.md` 7 KB | `README.md` — the standing rules, and what each script is |
| `scripts_dinamr/*.txt` | 1/1 | 5 KB | `grf_mechanism_output.txt` 5 KB | captured script output (GRF mechanism probe) |
| `dinamr_*.html` | 54/54 | 235.69 MB | `dinamr_C31_h100_n500_combine_1.html` 4.39 MB | `dinamr` batch and combine renders |
| `grfmr_*.html` | 30/30 | 130.99 MB | `grfmr_A124_h150_n500_combine_1.html` 4.38 MB | `grfmr` batch and combine renders |
| `grfprobe_*.html` | 5/5 | 21.56 MB | `grfprobe_g_p124_h150_n1500.html` 4.33 MB | `grfprobe` renders |
| `probe_*.html` | 10/10 | 43.14 MB | `probe_p31_h100_n500.html` 4.33 MB | Gate 1 cost-probe renders (`dinamr` era) |
| `summary_*.html` | 14/14 | 49.89 MB | `summary_dinamr.html` 7.62 MB | summary rendered outputs, all campaigns |
| `summary_*.qmd` | 14/14 | 371 KB | `summary_grfmr.qmd` 87 KB | summary sources, all campaigns |
| `sim_fs_maxeffCons_fb_mr_field_m1_template.qmd` | 1/1 | 159 KB | `sim_fs_maxeffCons_fb_mr_field_m1_template.qmd` 159 KB | **the** template `dinamr` / `grfmr` / the FS grid all render |
| `sim_*.qmd` | 53/53 | 4.19 MB | `sim_fs_maxeffCons_fb_mr_m1_h10_knoise0_n500_batch_1_100.qmd` 98 KB | other simulation templates (earlier campaigns and variants) |
| `sim_*.html` | 52/52 | 150.72 MB | `sim_fs_maxeffCons_fb_mr_m1_h10_knoise0_n500_batch_1_1000.html` 3.33 MB | renders of those other templates |
| `fs_*.html` | 245/245 | 973.62 MB | `fs_maxeffCons_fb_mr_field_m1_h175_knoise0_n1000_tier2_combine_1_2000.html` 4.34 MB | FS campaign batch/combine renders (`p12ext`, `tier2`, `e1stud`, `cert20`, earlier) |
| `*.html` | 29/29 | 87.87 MB | `smoke_p124_h150_n500_batch_1_5.html` 4.27 MB | remaining renders (smoke, gate, dflt, compare) |
| `REPORT_*.md` | 64/64 | 1.07 MB | `REPORT_fixedphat_ij2s_2026-09-09.md` 79 KB | REPORT documents |
| `TABLES_*.md` | 1/1 | 30 KB | `TABLES_grfmr_percell_2026-09-12.md` 30 KB | TABLES documents |
| `REVIEW_*.md` | 2/2 | 23 KB | `REVIEW_dinamr_blockC_2026-09-11.md` 13 KB | REVIEW documents |
| `current_status.md` | 1/1 | 18 KB | `current_status.md` 18 KB | this file — the directory's catalog at a pin |
| `*.md` | 1/1 | 2 KB | `payload_runbook_mr_only_20260819.md` 2 KB | other notes in the directory |
| `*.qmd` | 21/21 | 1.44 MB | `gate_d2_cim_unset.qmd` 156 KB | remaining `.qmd` |
| `*.R` | 4/4 | 63 KB | `p12ext_findings.R` 24 KB | top-level ad-hoc R scripts |
| `*` | 1/1 | 5 KB | `compare_1_20_vs_1_500.csv` 5 KB | everything else |
| **total** | **1325/1325** | **1948 MB** | `summary_dinamr.html` 7.62 MB | every file, each counted once |

Sizes are **apparent size** (`st_size`), not disk usage; `du` reports block-allocated size and reads larger for many small files. **Every file is counted exactly once** — the rules are applied first-match-wins and the rows sum to the total. **Files over 50 MB: 0; over 100 MB: 0.** The gitignored `_gateT_pre_template_files/` (1.7 MB) and `.DS_Store` are excluded throughout.

**Where to start, by question**

| If you need | Read |
|---|---|
| Per-cell coverage, bias, SDs, misses — DINA | `summary_dinamr.html`, or the Block reports |
| Per-cell numbers — GRF | `TABLES_grfmr_percell_2026-09-12.md` |
| FS classification and bound location | `REPORT_fs_extraction_2026-09-11.md` |
| Why a quoted FS range differs from the certification record | `REPORT_fs_products_reconciliation_2026-09-12.md` |
| How a number was computed | `scripts_dinamr/` — the summary's own chunks are authoritative over any re-implementation |
| Cost and wall provenance | `scripts_dinamr/logs/`, `projectC.R`, `walls.R` |
| Raw per-replicate rows | `results/*<campaign>*.rds` |

## 4. Reading conventions that must travel with these numbers

- Read every bound **by its location** against clinically meaningful effect sizes. Never frame a result as significance at HR = 1.00.
- **Selection rate is not an error rate.** At the differentially-null cells both FS and DINA select frequently and admissibly; what speaks to an unsupported claim is the share of lower bounds reaching HR 1.00 and 1.25.
- Coverage and bias as one comparative table (cells × estimators), then a short plain-language reading — not narrative prose.
- Marginal SD and error SD side by side. A Gaussian reference on the marginal SD understates the prediction by up to 7.5 points where the target moves with the estimate; both forms are in the summaries with their formulas stated.
- Wilson intervals on every rate. A Wilson interval belongs on pooled subject-level counts; a replicate-mean rate goes beside it without one.
- Cross-identifier comparisons are descriptive, never a ranking, and carry the confound: identifier, family construction, detection set, and at 12.4% the selection criterion.
- Evaluated estimator set: naive, oracle, IJ two-term, field, and Bonferroni for two-subgroup claims. Winner-only and winner-floor IJ variants are closed.
- `n_cons_qual`, `band_n` and `p_star` are structurally NA on DINA and GRF — not failures.

## 5. Superseded or known-wrong — do not quote

The first three entries name files that live **outside this repository** (the `fs_glms_interpretable` workstream and the manuscript build); they are listed because their numbers circulate alongside this directory's, not because they are here. Only the last entry is a file in `gbsg_020`.

- `claude/forestsearch_audit_spec.md` conditional-coverage figures: built from an earlier submission (different title, 35 + 85 pp.). The current build is 30 + 75 pp., dated 2026-08-20.
- Supplement §8.3 prose describing DINA/GRF conditional coverage as lower than FS's or recovering with n: does not agree with the adjacent tables and Figures S5/S7. Quote the tables and figures, not the prose.
- `BRIEF_dinamr_for_fs_glms_interpretable_2026-09-11.md` (v1): §8 listed FS classification metrics and bound-location shares as unavailable; both exist. Superseded by `BRIEF_fs_identifier_for_fs_glms_interpretable_2026-09-12.md`.
- `REVIEW_grfmr_2026-09-12.md` v1 framed its §2 as a blocking decision; withdrawn in the committed version.
- **A reconciliation of the FS one-sided product ranges (2026-09-12; per-cell values in `REPORT_fs_products_reconciliation_2026-09-12.md`).** An earlier pass of this file "corrected" the certification figures; **two of those corrections were wrong and are withdrawn**. The certification records are correct as written and were not edited. What actually differs:
  - **field lower on β(Ĥ)** — a genuine **cell-set** difference. `NOTE_survival_products_2026-09-09.md` reports 0.944–0.974 over its harm cells; its evidence list is `cert20` / `tier2` / `fixedphat_ij2s` / `field_studentize_e1` / `cimethod_flip` and does **not** include `REPORT_p12ext_2026-09-09`, the campaign that supplies the three 12.4% HR 1.50 cells. Excluding those, the committed bundles give exactly 0.944–0.974; including them gives 0.941–0.975, the 0.9410 coming from 12.4% HR 1.50 n 1000. *(The note says "ten harm cells"; nine is what reproduces the range, and no tenth bundle on disk carries field columns. Unreconciled, and it does not move the range.)*
  - **field-s upper on β(Ĥᶜ)** — **no difference at all.** The value is 0.9125–0.9605; "0.912–0.960" and "0.913–0.961" are the same endpoints rounded differently. The note's per-cell quotes reproduce exactly: 0.912 / 0.942 / 0.947 at 31% HR 1.50, 0.919 / 0.942 / 0.946 at HR 1.75, 0.941 / 0.956 / 0.961 at 12.4% HR 1.75.
  - **Bonferroni joint** — **different construction, not a different number.** The record's 0.939–0.963 is the **studentized** pair `fld_joint_s_bonf_*` (9-cell range 0.9395–0.9640). The 0.932–0.964 quoted against it was the **unscaled** pair `fld_joint_bonf_*`. Always name which.
  - Neither `SUMMARY_survival_properties_2026-09-10.md` nor `REVIEW_certification_2026-09-09.md` is in this repository, so the "0.941–0.980" attributed to them could not be checked; **0.980 does not reproduce from any field-lower computation on the committed bundles** (the nearest 0.98 in the record is the IJ two-sided at 31%, 0.971–0.981).

## 6. Not derivable from the committed columns

- **DINA non-detection causes** — the recorder returns before `n_family` is written, so an empty proposal and a proposal with nothing admitted are indistinguishable.
- **Per-subject membership** — so no classification metric finer than the four rates, and no agreement between two identifiers' Ĥ on the same draw.
- **2×2 counts** — not stored, but exactly recoverable from the recorded rates with `n_sel` and `n_true` (checked to 5.68e-14).
- **`dinamr` render logs** — written to a session scratchpad that no longer exists. The compute wall is rebuildable from bundle timing columns (`walls.R`); the per-render overhead for that campaign is not.

## 7. Open work in this directory

- GRF at the six HR 1.00 null cells, and 31% HR 1.75 at n 1000/1500 — ≈ 6 h.
- A criterion-matched FS comparator at 12.4%. None is committed; this needs compute and is the largest gap for a like-for-like low-prevalence comparison.
- Whether DINA and GRF should default to the field constructions — now informed by both grids.
- Whether to pin a commit in `forestsearch_version`.
