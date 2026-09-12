# `dinamr` campaign drivers and checkers

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
| `stage1_checks.R` | The Stage 1 battery: finiteness, the interval invariants, γ, realized prevalence, family size, and the presence of p̂ / ρᶜ / the nine recovery columns. |
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
