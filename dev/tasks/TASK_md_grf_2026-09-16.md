# TASK — GRF on the ACTG175 continuous (MD) design: campaign `mdgrf`, Stages 1–3

**File:** `dev/tasks/TASK_md_grf_2026-09-16.md` · **Issued:** 2026-09-16 by chat, on Larry's "proceed as you recommend" after the DINA/GRF Stage 0
**Machine:** `pop-os` · **Repo:** `~/Documents/GitHub/forestsearch` · **Branch:** `feature/glm-extension`, HEAD containing `700bdcdd`
**Transport:** `~/Downloads` → `dev/tasks/`, committed first and alone (§1.1)
**Directory:** `quarto/simulations/actg175/continuous/`, written `<dir>`
**Campaign tag:** `mdgrf`
**References, all committed:**
- `<dir>/REPORT_md_dina_grf_stage0_2026-09-16.md`, written `S0` — its sections S0.1 (identifiers), S0.2 (orientation and floors), S0.3 (one replicate), S0.4 (rule and MR), S0.5 (what a campaign needs)
- `dev/tasks/TASK_md_field_rerun_2026-09-15.md` — the FS re-run task, the transplant source for every stage of this task
- `<dir>/REPORT_md_field_rerun_stage1_2026-09-15.md`, `REPORT_md_field_rerun_gate2_2026-09-15.md`, `REPORT_md_field_rerun_2026-09-15.md`, `scripts_mdsgnb20/`, `summary_continuous_field_mdsgnb20.qmd`, `md_field_metrics.csv`, `COLUMNS_md_field.md` — the FS campaign's records, scripts, summary and extract

**What this is.** Campaign `mdsgnb20` evaluated the FS identifier on the ACTG175 continuous design under `effMaxSG`, ε = 0.20. This task runs the GRF identifier on the same four cells, replicates and seeds, recording the same constructions — unadjusted, oracle, IJ two-term, the field on Ĥ, field-s on Ĥᶜ with the unstudentized complement field beside it, the Bonferroni pair, p̂(Ĥ) — paired with `mdsgnb20` at the data level (`n_true` and the oracle columns). DINA is not run; it is blocked by a floor defect that needs an `R/` change, which comes as a separate proposal.

**Governing constraint.** GRF's candidate family is generated from a fitted surface, so the fixed-family condition does not hold. Every coverage figure in this task is coverage of the estimand conditional on the proposed family, and every table, figure caption and extract row says so. Comparisons with FS are descriptive, not a contest, and state the confound.

## Dispositions (Larry, 2026-09-16)

- **Identifier:** GRF only.
- **Rule:** `effMaxSG`, `effect_neighborhood = 0.20`, `selection_rule` as `mdsgnb20`'s `meta` records; the MD template's `FS_MD_FOCUS` and `FS_MD_NBHD` knobs.
- **GRF admission floor:** `dmin.grf = 30` on the harm-oriented MD scale, aligned with FS's effect threshold. S0.2's alignment statement names the argument value that does this; if S0.2 states that no argument places GRF's floor on that scale, that is a Gate 0 stop (§1.3).
- **Other GRF arguments:** as the survival `grfmr` campaign passed them, transplanted and quoted; a survival-only argument is omitted with the reason stated.
- **MR:** `ci_method = "field"` passed explicitly, with every MR argument the FS branch passes (`draws = 5000L`, `include_complement`, `field_complement`, `field_scale_complement = "selected"`, `return_reselection = TRUE`, and the `ij_residual`, `confirm_rule` and `t_confirm` values the MD template sets). Any of these that S0.4 lists as hard-coded on the GRF path to a value different from the FS branch's is a Gate 0 stop; no `R/` change.
- **Re-selection alignment:** S0.4 must state that MR's re-selection uses the same focus and band as the GRF identifier on this outcome. "Not aligned" is a Gate 0 stop: a misaligned correction is not run.
- **Labelling:** the conditional-family sentence on every table and figure, the confound sentence once at the top of the summary and the record (§3.1).
- **Gate 1 advance go (Larry offline):** if every Stage 1 gate is green and the §1.7 projection for Stage 2 is under 8 hours, do not stop at Gate 1. Note the advance go and its condition in the Gate 1 record, run Stage 2 at the worker count chosen in §1.7 with the record's ceiling and per-render timeout, and run Stage 3 once Stage 2 is green. Otherwise stop at Gate 1 and report. Never reduce cells, replicates, gates or knobs to fit under the threshold.

## ⚠ CATEGORY

- **No `R/` change.**
  - Edited: the MD template (identifier knob, GRF argument block, recorder fields, `meta` and poolability keys).
  - New: `scripts_mdgrf/`, a summary document, records, the extract, and the catalog update.
- **No install.** §1.2 asserts that `R/` is unchanged since the install of 2026-09-16.
- **Compute authorized by this kickoff:** Stage 1's smoke (20 FS + 20 GRF replicates) and calibration (333 GRF replicates), ceiling 2 h. Stage 2 and Stage 3 under the advance go above.
- **Unattended.** Gates stop on failure, never to ask. On a stop: commit what is green, record the failure, and stop. A statement in this document that does not hold is a finding, not a stop.

## Conventions

1. Verify from source; quote `path:line` at a stated commit.
2. Transplant, don't author: copy the committed precedent, change named lines, quote the source lines.
3. `git add` by explicit path only; untracked files never staged; no `fetch`, `pull` or `push`.
4. Every render runs with `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1`.
5. Numbers in records are computed and pasted, never typed.
6. `fit_mr_secs` contains `fld_H_secs` and `fld_Hc_secs`; nested timing columns are never summed.
7. Excluded from every table, figure and extract row: IJ winner-only and winner-floor, κ variants, `field_uniform`, covariate adjustment, tuned inflation factors.
8. Bounds are read by location against the D3 ladder (τ = 0, 10, 20, 30, 40, 60, 80, 100, harm-oriented MD), no significance language.
9. **Identification columns.** Sensitivity and PPV are computed from the recorder's columns as the recorder defines them (quote the lines). The true-positive count is the number of truly harmed patients in Ĥ; it is never equated with |Ĥ|. The `mdsgnb20` extract's `mean_n_harm` equals its `mean_n_sel` — a labelling defect noted as a finding and not transplanted; the FS extract is not re-run.
10. The submitted parent paper is not a topic.

---

## Stage 1 — edits, smoke, calibration → Gate 1

### 1.1 Provenance and first commit — GATE

```bash
cd ~/Documents/GitHub/forestsearch
hostname; git branch --show-current; git rev-parse --short HEAD; git log --oneline -3
git status --porcelain --untracked-files=no
git merge-base --is-ancestor 700bdcdd HEAD && echo "DINA/GRF Stage 0 in HEAD"
ps -eo pid,etime,args | grep -E '[e]xec/R|[R]script|[q]uarto' | head
```

*GATE:* the host is `pop-os`, the branch is `feature/glm-extension`, `700bdcdd` is in HEAD, there are no tracked modifications, and no R, Rscript or quarto process is running. Otherwise stop.

Copy this document from `~/Downloads` to `dev/tasks/TASK_md_grf_2026-09-16.md` (exact name, else the single match for `*TASK_md_grf_2026-09-16*.md`) and commit it alone.

### 1.2 Package check — GATE

- `git diff --quiet 0071c17e..HEAD -- R/ DESCRIPTION NAMESPACE` succeeds.
- `packageDescription("forestsearch")$Built` is `2026-09-16 05:57:14 UTC`, and two doFuture workers report the same.

*GATE:* both hold. No install is authorized.

### 1.3 Gate 0 from the Stage 0 record — GATE

Quote from `S0` and decide:
- **Alignment (S0.2):** the argument value that places GRF's admission floor on the harm-oriented MD scale. If S0.2 states there is none, stop.
- **Forwarding (S0.4):** for each MR argument in the dispositions, forwarded or hard-coded on the GRF path. A hard-coded value that differs from the FS branch's for any of them: stop, naming the argument and line.
- **Re-selection (S0.4):** aligned or not aligned for GRF on this outcome. Not aligned: stop.
- **Recorder fields (S0.5):** the list of DINA/GRF fields the MD template lacks.
- **Warnings (S0.3):** the two factor-comparison warnings, verbatim, for §1.6(c).

### 1.4 The survival reference — read only

One table, three columns — setting / MD template at HEAD / the survival GRF campaign — from `quarto/simulations/gbsg_020/sim_fs_maxeffCons_fb_mr_field_m1_template.qmd` at HEAD, the committed `grfmr` runner and campaign scripts (locate with `git ls-files`), and a `grfmr` bundle `meta`. Rows: the identifier knob and how each identifier is selected; every GRF argument (`dmin.grf`, forest and honesty settings, tree count, GRF seeding, the frontier filter and its band use, `vi.grf.min`); focus, band, rule; the MR arguments as passed on the GRF path. Every difference is listed as a finding. **Rule:** the campaign takes the survival campaign's GRF arguments except `dmin.grf = 30`; everything else stays the MD template's.

### 1.5 Template edits — transplant, no `R/`

Edit `<dir>/sim_fs_maxeffCons_mr_field_md_template.qmd` bottom-up, copying each block from the survival template and quoting its source lines:
- **E1 — identifier knob.** The survival template's identifier knob, under the `FS_MD_` prefix, default `consistency`; it drives `subgroup_method`, `use_dina` and `use_grf` as the survival template's does. The stem tag carries the identifier as the survival template's does.
- **E2 — GRF argument block**, applied only when the knob selects GRF: the §1.4 arguments, `dmin.grf = 30` as a knob `FS_MD_DMIN_GRF` with that default, GRF seeding as the survival template does it (quote).
- **E3 — recorder.** The §1.3 recorder fields and their fills, transplanted.
- **E4 — `meta` and poolability keys** carry the identifier, `dmin.grf` and the GRF arguments.

Then a default-path check: with the knob at its default, the template's rendered code path is unchanged (`git diff` shows additions guarded by the knob and the new recorder fields only). Commit the template.

### 1.6 Smoke — GATE

Transplant `scripts_mdsgnb20/` to `scripts_mdgrf/` (sampler, identity checker, `gate2.R`, runner), adding to the checker and `gate2.R` the GRF fields and the `meta` items below. Commit the scripts.

**(a) FS regression, md40 n500, sim_id 1–20,** knob at its default, campaign knobs as `mdsgnb20`'s, tag `mdgrfsmokefs`: every column except `*_secs` against `mdsgnb20`'s combined bundle, identical within 1e-8 relative, label ties enumerated, zero selection flips. This protects the template edits.

**(b) GRF, md40 n500, sim_id 1–20,** knob `grf`, tag `mdgrfsmoke`, 20 workers:
- runs without error; warnings captured verbatim;
- `n_true` identical and the oracle columns within 1e-8 of `mdsgnb20`'s on all 20 sim_ids;
- `meta` carries identifier `grf`, `dmin.grf` 30, focus `effMaxSG`, band 0.20, the rule, `ci_method` field, `field_scale_complement` selected, `pkg_version` 0.3.5, host `pop-os`;
- the floor as applied, read from the returned object or the recorded family, is 30 on the harm-oriented scale;
- the E3 fields are filled on every declared replicate; the nine `fld_Hc_*_s` and nine `fld_joint_s_*` columns are finite where the complement block is filled; `fld_Hc_lo1s_s ≤ fld_Hc_up1s_s`; the Bonferroni harm bound is identical between `joint` and `joint_s` where draw counts agree; p̂(Ĥ) is recorded;
- report, as facts: the declared count; sim_id 1's selection beside S0.3's; `fit_mr_secs` and the GRF fit time (mean, median, max).

**(c) The warnings.** Locate the source line of each factor-comparison warning; state what is compared and why the candidate evaluation is correct. If the comparison could mis-evaluate or drop a candidate on this design, stop.

*GATE:* (a) and (b) hold and (c) is resolved.

### 1.7 Calibration

md40 n700, knob `grf`, campaign knobs, `FS_MD_WORKERS` = 16, 32 and 63 with 3 × workers replicates each (tags `mdgrfcal16/32/63`), memory sampled every 5 s. Report per worker count: render wall; `fit_mr_secs` and GRF fit time (mean, median, 90th percentile); peak summed RSS; replicates per minute. Choose W (minimum projected wall; name what limits it). Project Stage 2 at W: four cells × (two 1,000-replicate batches + combine), scaling per-cell cost by `mdsgnb20`'s cell ratios and taking fixed overhead from these renders. State the projection, a ceiling of 1.5 × the projection, and a per-render timeout of 2 × the longest projected batch (at least 20 min).

### 1.8 Gate 1 record

`<dir>/REPORT_md_grf_stage1_2026-09-16.md`: provenance and the package check; §1.3 with quotations; the §1.4 table; the edits with quoted sources; the smoke results and the warning resolution; the calibration table, W, projection, ceiling, timeout; the advance-go condition and whether it is met; findings; `git log --oneline` for the stage. Commit it; smoke and calibration outputs stay untracked and are listed. Then apply the advance go, or stop with the closing message.

---

## Stage 2 — campaign `mdgrf`

### 2.1 Preconditions — GATE

HEAD descends from the Stage 1 record commit with no change since under `R/` or to the template; `Built` equals §1.2's; no R, Rscript or quarto process is running.

### 2.2 Launch

The transplanted runner, detached, with knob `grf`, `FS_MD_DMIN_GRF=30`, the rule knobs, `FS_MD_FIELD_SCALEC=selected FS_MD_CI=field FS_MD_CAMPAIGN=mdgrf`, `FS_MD_WORKERS=W` and the timeout; progress log `<dir>/LOG_mdgrf_progress.txt`, halt file `<dir>/HALT_mdgrf.md`, raw logs untracked under `<dir>/logs_mdgrf/`. Cells in `mdsgnb20`'s order: md40 n500, md120 n500, null n500, md40 n700.

### 2.3 Per-cell Gate 2

`gate2.R` checks: 2,000 rows with `sim_id` exactly 1–2000, batch files matching the combined bundle on every column; `meta` as in §1.6(b); **same draws as `mdsgnb20`, both directions:** `n_true` identical and the oracle columns within 1e-8 on all 2,000 sim_ids (complement oracle only in the null cell); §1.6(b)'s field-s and recorder checks on every declared replicate; MR failures on declared replicates at most 40. On pass, commit by explicit path the cell's result directory, its combine HTML, its section of `<dir>/REPORT_md_grf_gate2_2026-09-16.md`, and the progress log. On failure, write the halt file, commit what is green, stop. Ceiling: finish the current cell, then halt and stop.

---

## Stage 3 — summary, extract, record, catalog, closeout

### 3.1 Summary — transplant

Copy `summary_continuous_field_mdsgnb20.qmd` to `summary_continuous_field_mdgrf.qmd` (and `stage3_aggregate.R` into `scripts_mdgrf/` if used). Inputs: the four `mdgrf` bundles; `mdsgnb20`'s four for pairing; the committed `md_field_metrics.csv` for FS comparator rows, not recomputed. At the top, once, verbatim:

"GRF's candidate family is generated from a fitted surface, so the fixed-family condition does not hold: every coverage figure below is coverage of the estimand conditional on the proposed family. Comparisons with FS are descriptive, not a contest: the identifier and the family construction differ, and each summary conditions on a different set of detected replicates."

Every table and figure caption ends with "(conditional on the proposed family)". Tables:
- **Ĥ and Ĥᶜ** — cells × unadjusted, oracle, IJ two-term, field (Ĥ) / field-s (Ĥᶜ, the unstudentized field once beside): declaration rate, bias in MD and SD units, SE/SD, one-sided coverage, two-sided coverage.
- **Bound location** — the ladder with Monte Carlo SEs: field | oracle on Ĥ, field-s | oracle on Ĥᶜ.
- **Joint** — the field-s pair's joint coverage and share of declared replicates with both bounds; the unstudentized pair once beside; γ and correlation as diagnostics.
- **Identification** — GRF: declaration rate, mean |Ĥ|, sensitivity, PPV, mean family and admitted sizes; FS beside, from the CSV; paired by sim_id: replicates on which GRF's |Ĥ| is larger, equal, smaller than FS's; per convention 9.
- **Regime diagnostics** — p̂(Ĥ) mean and share below 0.5; SD(β̃ᶜ)/naive SEᶜ; λ-SDᶜ/naive SEᶜ for field and field-s.
- **Display** — `fs_sim_bias_coverage(scale = "identity")`, both blocks, field-s via the renamed copy.
- **Null cell** labelled by its truth in every table.

Render; commit the `.qmd`, HTML and figures.

### 3.2 Extract

`<dir>/md_grf_metrics.csv` and `<dir>/COLUMNS_md_grf.md`, written from the same objects the tables print, with `md_field_metrics.csv`'s schema plus an `identifier` column (`grf`; FS comparator rows carry `fs` and are copied from the committed CSV with their commit). Monte Carlo SEs and Wilson intervals as in `COLUMNS_md_field.md`. The COLUMNS file states the conditional-family reading in its scale-convention section.

### 3.3 Record

`<dir>/REPORT_md_grf_2026-09-16.md`: provenance and dispositions; Gate 1 and Gate 2 in brief; each §3.1 table pasted from the render with a two-to-four-bullet reading (numbers from the CSV, bounds by location); the confound paragraph from §3.1; this scope sentence verbatim: "These are operating characteristics of the GRF identifier on one continuous design, conditional on the proposed family. They do not verify condition (A3), no construction is promoted on this design's performance, and DINA was not run."; findings and commits.

### 3.4 Catalog

Add to `<dir>/status_curated.md`: campaign `mdgrf` with payload locations and the conditional reading convention; open work: "DINA on this design: blocked pending the proposal-floor orientation fix (S0.2)". Commit it, regenerate `current_status.md` as the last commit, `check_current_status.sh --commit` passes.

### 3.5 Closeout

Delete the Stage 1 smoke and calibration outputs and nothing else. Copy `md_grf_metrics.csv`, `COLUMNS_md_grf.md` and the record to `~/Downloads`. Post-conditions, printed: no tracked modification; `git diff --quiet <§1.1 HEAD>..HEAD -- R/` succeeds; `Built` unchanged; the four `mdgrf` combined bundles tracked with 2,000 rows; every number in the record's tables found in the CSV at printed precision; the catalog check passes; the `~/Downloads` copies `cmp`-identical. Closing message: the commit range to push, the post-conditions, the findings, and the paths of the record, the summary HTML and the CSV. Then stop.

## Out of scope

- `R/`; DINA; the applied documents; any threshold, search or MR setting beyond the dispositions; the excluded constructions; pushing.
