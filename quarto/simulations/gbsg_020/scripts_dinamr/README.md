# `dinamr` campaign drivers and checkers

> ## Standing rules for `quarto/simulations/gbsg_020`
>
> These govern every task that touches this directory, not only the campaigns below.
>
> **1. Artifact tracking.** Every artifact carrying payload information or analysis/summary
> content is **committed**, so other repositories can read this directory **at a pin**, build
> derived payloads and run their own summaries. In scope: per-replicate bundles and their metas;
> every driver, checker, projection and extraction script; the summary `.qmd` **and** its
> rendered `.html`; batch and combine renders; `REPORT_*`, `TABLES_*`, `REVIEW_*`; every saved
> derived `.rds`; and `logs/`, whose `WALL_SECONDS` and `CELL DONE ... wall=` lines are the
> provenance for the Gate 1 cost accounting. **Nothing another repo would need stays in a
> session scratchpad.** Report every file's path, tracked status and size before finishing;
> **flag anything over 50 MB before committing it, and treat 100 MB as a hard stop**
> (`REPORT_push_size_fix_2026-08-31.md`: a 105 MB blob forced a history rewrite).
>
> **2. Closeout.** Every task regenerates **`current_status.md`** as its **last action** — the
> directory's catalog at a pin, so another machine can be brought up to speed by attaching that
> one file. Refresh the pin and date, regenerate the §3 payload inventory **from the directory,
> never from chat records**, re-verify the curated sections against the files, and report every
> correction. **Post-condition, machine-checked: the stated pin equals HEAD at commit time.**
> Commit `current_status.md` **alone**, as the child of the commit it describes, so its parent
> *is* that HEAD. Enforce with:
>
> ```sh
> ./scripts_dinamr/check_current_status.sh            # before: stated pin == HEAD
> ./scripts_dinamr/check_current_status.sh --commit   # after:  stated pin == HEAD~1, and the
>                                                     #         commit touched only that file
> ```


Everything needed to reproduce the DINA MR campaign (`TASK_dinamr_campaign_2026-09-10`,
`REPORT_dinamr_2026-09-10.md`) from the repository. These ran from a session scratchpad during the
campaign and are committed here afterwards so a fresh session does not have to reconstruct them.

**Paths were session-absolute and are now portable.** Each script resolves two locations, both
overridable:

- `DINAMR_QMD_DIR` — the directory holding `sim_fs_maxeffCons_fb_mr_field_m1_template.qmd`
  and `results/`. Defaults to this script's parent, i.e. `quarto/simulations/gbsg_020`.
- `DINAMR_SCRATCH` — where logs and intermediate `.rds` go. Defaults to this directory.

| file | role |
|---|---|
| `render.sh` | One render. **Exports `VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1`** (Apple Accelerate is not fork-safe; tier2's choice on this host) and times the render. |
| `campaign.sh` | Block driver: for each cell, batch 1 → batch 1001 → combine. **Pins `FS_S7_WORKERS=12`** and the full campaign knob set (`FS_S7_METHOD=dina FS_S7_FOCUS=effMaxSG FS_S7_NBHD=0.20 FS_S7_FIELD_COMPLEMENT=TRUE FS_S7_FIELD_SCALEC=selected FS_S7_FIELD_DECOMP=TRUE FS_S7_FIELD_RECOV=TRUE FS_S7_IJ_RESIDUAL=two_term FS_S7_FB=none FS_S7_CAMPAIGN=dinamr`). **`FS_S7_ER_JCUTS` is deliberately unset** — inert on DINA. |
| `blockA.cells` / `blockB.cells` | The cell lists as run. `blockB.cells` holds **five** cells: HR 1.75 n 1500 is the checkpoint's deferred cell. |
| `probe.sh` | The ten 36-replicate Gate 1 cost probes (tag `dinamrprobe`). |
| `project.R` | Gate 1 projection from the probe corners. |
| `checkpoint.R` | The Amendment 2 post-Block-A re-projection, calibrated by n, against the 13 h ceiling. |
| `stage1_checks.R` | The Stage 1 battery: finiteness, the interval invariants, γ, realized prevalence, family size, and the presence of p̂ / ρᶜ / the nine recovery columns. Its `bonf vs raw` line is the **wrong-pair** comparison described below, kept and relabelled as what it is; the corrected identity is computed beside it. |
| `check_current_status.sh` | **The closeout post-condition**: stated pin in `current_status.md` == HEAD at commit time. |
| `gate2G.R` | Gate 2 for the **GRF** cells: `gate2.R`'s content plus the `admitted_n` distribution, with `n_family` relabelled as the outcome-independent enumerated pool. |
| `gate3.R` | GRF alignment per batch: `grf_select_statistic`, `grf_selection`, `dmin.grf`, resolved three ways. |
| `grfmr.sh`, `grfmr.cells`, `grfmr_deferred.cells` | The `grfmr` Part A driver and the ten-plus-two split. |
| `projectG.R`, `grfmr_numbers.R`, `grfmr_tables.R`, `fs_extraction.R` | `grfmr` Gate 1 projection, the per-cell accumulator, the table extractor that re-executes the summary's own chunks, and the Part C FS extraction. |
| `grfmrC_smoke.sh`, `stage1G.R` | `grfmr` completion (`TASK_grfmr_completion_2026-09-12`): the Stage 1 GRF null smoke and its battery. |
| `projectGC.R`, `grfmrC.sh`, `grfmrC.cells`, `wallsGC.R` | The completion's Gate 1 projection (calibrated ×0.8962), the eight-cell driver (Gate 3 per batch, stop-on-failure, 10 h watchdog), its cell list, and realized walls against the projection. `gate2G.R C` gates the six HR 1.00 cells. |
| `partB_stage0_readout.R` | `TASK_partB_measurement_2026-09-12`, Stage 0 readout from committed bundles only: MR-on per-replicate cost at the three comparators, the field-excluded bound, FS cost across `sg_focus`, the pop-os/Mac host factor, 18-cell bounds, and Wilson resolution. **Not a measurement:** Stage 0c stopped the task before Stage M. |
| `t3gate.sh`, `t3gate.R` | `TASK_partB_enabling_2026-09-12` Gate T3. Standing identity cell, 5 replicates: pre- vs post-change with `FS_S7_MR` unset (every non-timing column and truth `identical()`), then MR on vs MR off (identification and classification `identical()`). Also runs a supplementary DINA/GRF on/off pair. |
| `partBoc.sh`, `partBoc_table.R` | Part B OC smoke: 16 runs (consistency × 6 criteria, DINA and GRF × 5) at 12.4% HR 1.50 n 500, 30 replicates, MR off, tag `pBoc`. Produces the OC table, within-engine same-subgroup agreement, and the 288-cell-run projection. |
| `partBoc_checks.R` | Beside the OC table. (1) Compares the smoke's MR-off selections with the committed MR-on bundles at the same cell and seeds (`dinamr`, `grfmr`, `p12ext`, sim_id 1–30). (2) Gives an approximate cell-profile-adjusted projection from `partB_stage0_readout.rds`. |
| `idsweep.sh`, `idsweep.cells` | `TASK_idsweep_2026-09-12`, the Part B identification sweep: 18 cells × 16 runs (consistency × 6 criteria, DINA and GRF × 5), 500 replicates, one batch per run, MR off, tag `idsweep`. Gate A per run and Gate I per cell (both stop-on-failure), a Gate 1 re-projection before every cell (defer from the tail past 13 h), and a 16 h watchdog. `IDSWEEP_T0` / `IDSWEEP_FROM` resume a stopped run without re-running completed cells, with the clock kept from the original start. |
| `idsweep_gateA.R`, `idsweep_gateI.R` | Gate A (alignment: MR off / `_nomr`, `sg_focus`, ε, GRF `grf_selection` and `frontier_rule` resolved from the installed package and from source) and Gate I (integrity, structural NA, realized prevalence, same draws against committed bundles). Gate I carries **Amendment 1**: a classification rate may be NA on a detected replicate only where its denominator is zero, e.g. NPV when the selection is the whole trial. |
| `idsweep_project.R`, `idsweep_findings.R` | `gate1` (the reconciled projection), `next` (the per-cell re-projection from realized cost), `walls` (realized against projected); and the across-cell readout `REPORT_idsweep` quotes. |
| `status_inventory.R` | Regenerates `current_status.md` §3, the payload inventory, **from the directory** (first-match-wins rules, tracked/disk counts, apparent sizes, 50 MB / 100 MB flags). |
| `gate2.R` | **The Gate 2 checker.** Usage: `Rscript gate2.R A` / `B` / `C`. Carries the **corrected** bound↔quantile identity (see below) and Amendment 3 at `TOL_TRUTH <- 1e-8`. |
| `blockA_numbers.R`, `blockA_rest.R`, `chunkdiff.py` | Ad-hoc readouts that produced report tables (the Block A standard tables, the stratified tables, and the per-chunk transplant accounting). Not part of the run. |

## The corrected bound↔quantile identity, and what it replaced

The first version of `gate2.R` compared `fld_joint_bonf_loH` against `fld_joint_loH`. **That is the
wrong pair.** They are not the same quantity: the Bonferroni bound is built at the replicate's own
γ and the raw one at 0.025, so they coincide **exactly when γ sits at its 0.025 floor** and differ
otherwise. The check therefore passed on the 5-replicate smokes (γ at the floor throughout) and read
`0.0213` on the first 2,000-replicate cell.

`gate2.R` now gates the identity `REPORT_cert20_2026-09-08` actually names — **field-s is inverted
around the same `bdc` as the unscaled complement field**:

```
log(fld_Hc_est2_s) + fld_Hc_lam_mean_s  ==  log(fld_Hc_est2) + fld_Hc_lam_mean
```

which lands at **2.22e-16** on every cell of the campaign. The `bonf == raw` identity is kept, gated
on the γ-at-floor rows only, with the share at the floor reported beside the γ range as `cert20` did.

## Verified to run from this directory

All three `.sh` pass `zsh -n`; all seven `.R` parse; `chunkdiff.py` parses. Run from here against
the committed bundles, every script reproduces its recorded result:

| script | reproduces |
|---|---|
| `gate2.R A` | **204 passes, 0 failures** |
| `gate2.R B` | **170 passes, 0 failures** |
| `project.R` | the Gate 1 projection verbatim — A 2.12 h, B 8.53 h, C 2.17 h, A+B 10.65 h, all-18 12.83 h |
| `stage1_checks.R` | the Stage 1 battery on both smokes — 9/9 recovery columns, all finite, γ in range |

No path under this directory contains an absolute session path.

## Reproducing

```sh
cd quarto/simulations/gbsg_020/scripts_dinamr
./probe.sh                                   # Gate 1 probes
Rscript project.R                            # Gate 1 projection
./campaign.sh A blockA.cells                 # Block A, six cells
Rscript gate2.R A                            # Gate 2
Rscript checkpoint.R                         # Amendment 2 checkpoint
./campaign.sh B blockB.cells                 # Block B, five cells
Rscript gate2.R B
cd .. && quarto render summary_dinamr.qmd --output summary_dinamr.html
```

`checkpoint.R` reads the campaign driver's own stdout for the realized walls; point it at that log
(it defaults to the path the campaign used) or pass the walls in directly.
