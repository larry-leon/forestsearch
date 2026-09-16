# TASK — ACTG175 continuous (MD) re-run under the current field constructions: Stages 1–3

**File:** `dev/tasks/TASK_md_field_rerun_2026-09-15.md` · **Issued:** 2026-09-15 by chat, on Larry's approval of the Gate 0 recommendations
**Machine:** `pop-os` · **Repo:** `~/Documents/GitHub/forestsearch` · **Branch:** `feature/glm-extension`, HEAD containing `b6e30ac7` (the Stage 0 record)
**Transport:** `~/Downloads` → `dev/tasks/`, committed first and alone (§1.1)
**Directory:** `quarto/simulations/actg175/continuous/`, written `<dir>` below; every record, script and output of this task lives there
**Campaign tag:** `mdsgnb20` · **Stage 0 record:** `<dir>/REPORT_md_field_rerun_stage0_2026-09-15.md`

**What this is.** Campaign `mdf1` (2026-09-07) evaluated the MR intervals on the ACTG175 continuous design under `maxeffCons`, before field-s existed. This task re-runs its four cells × 2,000 replicates on the same seeds under the survival grid's selection rule, `effMaxSG` at ε = 0.20. It records the current constructions: the field one-sided lower bound on β(Ĥ), the field-s one-sided upper bound on β(Ĥᶜ), and their Bonferroni pair, with unadjusted, oracle and IJ two-term as references. It ends with a summary, the extract `md_field_metrics.csv`, and a record.

## Gate 0 dispositions (Larry, 2026-09-15)

- **D1, selection rule:** `effMaxSG`, ε = 0.20, with the `selection_rule` of the survival `effMaxSG` campaigns (§1.3). The same-draws anchor to `mdf1` is data-level: `n_true` and the oracle columns.
- **D2, null cell:** presented, labelled by its truth (no subgroup; homogeneous +26 on the harm-oriented scale), with its declaration rate beside it.
- **D3, thresholds:** one bound-location ladder, τ = 0, 10, 20, 30, 40, 60, 80, 100 on the harm-oriented MD scale, for both blocks in every cell, with the oracle beside.
  - Ĥ: lower bound ≥ τ, read as "harm of at least τ supported".
  - Ĥᶜ: upper bound ≤ τ, read as "harm of at most τ supported".
- **D4, machine:** `pop-os`.
- **D5, applied document:** gains field-s rows in its own task after this one; not here.
- **Also:**
  - the recorder, knob and `meta` edits of §1.4 (template only);
  - the display reads field-s through a renamed copy (§3.1), with no `R/` change;
  - scripts, the progress log and the directory catalog are committed (§1.5, §2.3, §3.4);
  - the extract is new code (§3.2).

## ⚠ CATEGORY

- **No `R/` change.**
  - Edited: the MD template (knobs, recorder, `meta`).
  - New: scripts under `<dir>/scripts_mdsgnb20/`, a summary document, records, the extract, and the directory catalog.
- **Install:** `devtools::install()` from HEAD at §1.2, a rebuild of unchanged `R/`.
- **Compute authorized by the Stage 1 kickoff:** the install, the smoke (§1.6, 100 replicates) and the calibration (§1.7, 333 replicates). That is about 30–45 min of compute, with a ceiling of 1.5 h.
- **Stages 2–3 run only on Larry's Gate 1 go**, which states the worker count and the wall-clock ceiling for 4 cells × 2,000 replicates.
- **Unattended.**
  - Gates stop on failure, never to ask.
  - On a stop: commit what is green, record the failure, and stop.
  - A statement in this document that does not hold is a finding, not a stop.

## Conventions

1. **Verify from source.** Quote `path:line` at a stated commit.
2. **Transplant, don't author.** Where a precedent exists, copy it, change named lines, and quote the source lines.
3. **Git.** `git add` by explicit path only; untracked files are never staged; no `fetch`, `pull` or `push`.
4. **Threads.** Every render runs with `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1`.
5. **Numbers.** Numbers in records are computed and pasted, never typed.
6. **Timing.** `fit_mr_secs` contains `fld_H_secs` and `fld_Hc_secs`; nested timing columns are never summed.
7. **Excluded** from every table, figure and extract row: IJ winner-only and winner-floor, κ variants, `field_uniform`.
8. **Reading bounds.** Bounds are read by location against the D3 ladder, with no significance language anywhere.
9. **Names and scope.** Record names carry this task's date. The submitted parent paper is not a topic.

---

## Stage 1 — edits, smoke, calibration → Gate 1 (then STOP)

### 1.1 Provenance and first commit — GATE

```bash
cd ~/Documents/GitHub/forestsearch
hostname; git branch --show-current; git rev-parse --short HEAD; git log --oneline -3
git status --porcelain --untracked-files=no
git status --porcelain -- R DESCRIPTION NAMESPACE
ps -eo pid,etime,args | grep -E '[e]xec/R|[R]script|[q]uarto' | head
git merge-base --is-ancestor b6e30ac7 HEAD && echo "Stage 0 record in HEAD"
```

*GATE:* all of the following hold, or stop.
- The host is `pop-os` and the branch is `feature/glm-extension`.
- `b6e30ac7` is in HEAD.
- There are no tracked modifications, and nothing tracked or untracked is pending under `R/`, `DESCRIPTION` or `NAMESPACE`.
- No R, Rscript or quarto process is running.

Copy this document from `~/Downloads` to `dev/tasks/TASK_md_field_rerun_2026-09-15.md` and commit it alone. Use the exact name; failing that, the single match for `*TASK_md_field_rerun_2026-09-15*.md`.

### 1.2 Install — GATE

- Run `Rscript -e 'devtools::install(quick = TRUE, upgrade = "never")'`.
- Record `packageDescription("forestsearch")$Built`.
- Assert that two doFuture workers each see that same `Built`.
- Assert the installed `fs_mr_inference()` formals as Stage 0 §3 recorded them: `ci_method` defaults to `"field"` and `field_scale_complement` to `"selected"`.

*GATE:* all hold.

### 1.3 The survival reference — read only

Build a side-by-side table with three columns.
- **Setting.**
- **MD template:** `<dir>/sim_fs_maxeffCons_mr_field_md_template.qmd`.
- **The FS survival `effMaxSG` campaigns**, read from three sources:
  - `quarto/simulations/gbsg_020/sim_fs_maxeffCons_fb_mr_field_m1_template.qmd` at HEAD;
  - the knobs of the committed `p12x20` runner (`run_p12x20.sh`, `campaign_p12x20.sh`; locate with `git ls-files`);
  - the `p12x20` and `cert20` bundle `meta`.

Rows:
- **Identification:** `sg_focus`, `effect_neighborhood`, `selection_rule`, `stop_threshold`, `consistency_method`, `pconsistency.threshold`, `fs.splits`, `maxk`, `n.min`, arm minima, seed scheme.
- **MR:** `draws`, `multiplier`, `ci_method`, `ij_residual`, `confirm_rule`, `t_confirm`, `field_complement`, `field_scale_complement`, `return_reselection`, field `R_out`/`R_in`.

**Rule:** the campaign takes the survival campaigns' `sg_focus`, `effect_neighborhood` and `selection_rule`. Every other setting stays the MD template's, and each difference is listed as a finding for Gate 1.

### 1.4 Template edits — transplant, no `R/`

Edit the MD template bottom-up so the Stage 0 line numbers hold. Copy each block from the survival template and quote its source lines.
- **E1 — focus and band knobs.**
  - Replace `sg_focus <- "maxeffCons"` (`:135`) and `effect_neighborhood <- 0.10` (`:214`) with the survival template's focus and band knob blocks, including any guard and stem tagging they carry.
  - Rename them `FS_MD_FOCUS` / `FS_MD_NBHD`, with defaults `"maxeffCons"` / `0.10`.
  - If §1.3 shows a survival `selection_rule` other than `"neighborhood"`, also transplant its knob as `FS_MD_RULE`, default `"neighborhood"`.
- **E2 — complement-scale knob.** Add the survival template's complement-scale knob as `FS_MD_FIELD_SCALEC`, default `"selected"`, and pass it in `mr_inference_args` (`:282–287`).
- **E3 — recorder.** Add the four recorder and fill blocks listed in Stage 0 §5.2 (a)–(d).
- **E4 — invariants.** Add the `_s` interval invariant pairs at `:824–827`.
- **E5 — `meta`.** Add `effect_neighborhood`, `selection_rule` and `field_scale_complement` to `meta` (`:892–921`) and to the poolability keys (`:945–948`).

Commit the template.

### 1.5 Scripts — transplant

Create these under `<dir>/scripts_mdsgnb20/`, quoting their sources.
- **`mem_sampler.sh`**, from `scripts_mdf1/mem_sampler.sh`, using Linux `ps`.
- **`smoke_identity.R`**, from `scripts_mdf1/identity_postchange.R` (or from whatever wrote `mdf1`'s `gate2_flips.txt`, if that script does not compare bundles column by column).
  - It applies `mdf1`'s Gate 2 classification as `REPORT_continuous_field_gate2_2026-09-07.md` defines it: tolerance 1e-8 relative, label-tie and MR-numerics rows enumerated, selection flips counted.
- **`gate2.R`**, from the committed `gate2F.R` used by `p12x20` (combine checks, same draws in both directions), pointed at MD bundles, with the checks of §2.3.
- **`run_mdsgnb20.sh`**, from the committed `p12x20` runner.
  - It runs detached, writes the progress log `<dir>/LOG_mdsgnb20_progress.txt` (each render's wall and peak memory) and the halt file `<dir>/HALT_mdsgnb20.md`, and commits per cell by explicit path.
  - Its render lines come from `scripts_mdf1/run_cell.sh`: batch 1–1000, batch 1001–2000, combine.
  - It uses GNU `timeout` in place of `tmo.sh`. Raw render logs stay untracked under `<dir>/logs_mdsgnb20/`.

Commit the scripts.

### 1.6 Smoke — GATE

**(a) Defaults, all four cells.**
- Run `FS_MD_CAMPAIGN=mdsmoke FS_MD_START=1 FS_MD_NSIMS=20` at 20 workers.
- `smoke_identity.R` compares every column of `mdf1`'s recorder except `*_secs` against sim_id 1–20 of `mdf1`'s combined bundles.

**(b) Field-s wiring.** On every smoke replicate where the complement field block was filled:
- the nine `fld_Hc_*_s` and nine `fld_joint_s_*` columns are finite;
- `fld_Hc_lo1s_s ≤ fld_Hc_up1s_s` and `fld_Hc_lo2s_s ≤ fld_Hc_hi2s_s`;
- where the `joint` and `joint_s` draw counts agree, their Bonferroni harm bounds are identical (both come from the same harm draws);
- `meta` records `field_scale_complement = "selected"`.

**(c) The campaign rule, md40 n500, sim_id 1–20.**
- Run with the §1.3 rule: `FS_MD_FOCUS=effMaxSG FS_MD_NBHD=0.20`, plus `FS_MD_RULE` if set.
- The stem and `meta` carry the rule.
- Against `mdf1`, `n_true` is identical and the oracle columns agree within 1e-8 relative.
- (b)'s checks pass.
- Report on how many of the 20 replicates |Ĥ| grew, stayed or shrank. This is a fact, not a gate.

*GATE:*
- In (a), zero selection flips, and every other difference falls in a `mdf1` Gate 2 class.
- (b) and (c) hold.

### 1.7 Calibration

Run md40 n700 under the campaign knobs at `FS_MD_WORKERS` = 16, 32 and 63, with 3 × workers replicates each (tags `mdcal16`, `mdcal32`, `mdcal63`), sampling memory every 5 s. For each worker count, report:
- the render wall;
- `fit_mr_secs` mean, median and 90th percentile, and its ratio to the 16-worker mean;
- peak summed RSS;
- replicates per minute.

Then:
- **Choose W.** Pick the worker count that minimizes projected wall, and name what limits it: cores, memory, or scaling falloff.
- **Project Stage 2** at W: four cells × (two 1,000-replicate batches + combine), scaling per-cell cost by `mdf1`'s recorded `fit_mr_secs` ratios and taking fixed render overhead from these renders.
- **State the limits:** the projection, a ceiling of 1.5 × the projection, and a per-render timeout of 2 × the longest projected batch (at least 20 min).

### 1.8 Gate 1 record, then STOP

Write `<dir>/REPORT_md_field_rerun_stage1_2026-09-15.md`. It contains:
- provenance and install;
- the §1.3 table;
- the edits, with quoted sources, and the scripts;
- the smoke results — per cell: rows compared, rows within tolerance, rows enumerated by class, flips — and the (b)/(c) checks;
- the calibration table, W, the projection, the ceiling and the timeout;
- findings, and `git log --oneline` for this stage.

Commit it. The smoke and calibration outputs stay untracked and are listed in the record.

**Closing message:** the commit range to push, W, the projection, the ceiling, and the findings. Then stop.

---

## Stage 2 — campaign `mdsgnb20` (only on Larry's Gate 1 go)

### 2.1 Preconditions — GATE

- HEAD descends from the Stage 1 record commit, with no change since under `R/` or to the template.
- The installed `Built` equals §1.2's.
- No R, Rscript or quarto process is running.

### 2.2 Launch

Launch `run_mdsgnb20.sh` detached, with:
- the §1.3 rule knobs;
- `FS_MD_FIELD_SCALEC=selected FS_MD_CI=field FS_MD_CAMPAIGN=mdsgnb20`;
- `FS_MD_WORKERS` and the timeout from the go.

Run the cells in `mdf1`'s order: md40 n500, md120 n500, null n500, md40 n700.

### 2.3 Per-cell Gate 2

`gate2.R` checks:
- the combined bundle has 2,000 rows with `sim_id` exactly 1–2000, and the batch files match it on every column;
- `meta` carries the campaign rule, `field_scale_complement = "selected"`, `pkg_version` 0.3.5 and host `pop-os`;
- **same draws as `mdf1`, in both directions:** `n_true` is identical and the oracle columns agree within 1e-8 relative on all 2,000 sim_ids. In the null cell, only the complement oracle is checked, since the harm oracle is empty.
- §1.6(b)'s field-s checks hold on every replicate where the complement field block was filled;
- MR failures on declared replicates number no more than the larger of 20 and twice `mdf1`'s count in that cell.

**On pass**, commit by explicit path:
- the cell's result directory;
- its combine HTML, as `mdf1` tracked it;
- its section of `<dir>/REPORT_md_field_rerun_gate2_2026-09-15.md`;
- the progress log.

**On failure**, write the halt file naming the failing check, commit what is green, and stop.

**Ceiling:** if cumulative wall exceeds the go's ceiling, finish the current cell, then halt and stop.

---

## Stage 3 — summary, extract, record, closeout (after Stage 2 is green)

### 3.1 Summary document — transplant

Copy `<dir>/summary_continuous_field_mdf1.qmd` to `<dir>/summary_continuous_field_mdsgnb20.qmd`, together with `scripts_mdf1/stage3_aggregate.R` into `scripts_mdsgnb20/` if the summary relies on it. Named changes:
- **Inputs:** the four `mdsgnb20` combined bundles, plus `mdf1`'s four for the paired sections.
- **Ĥ table** — cells × unadjusted, oracle, IJ two-term, field. Columns: declaration rate, bias on the MD scale and in SD units, SE/SD, one-sided lower coverage, two-sided coverage.
- **Ĥᶜ table** — the same, with field-s as the evaluated row and the unstudentized field shown once beside it as the paired before/after.
- **Bound location** — the D3 ladder, with Monte Carlo SEs: field | oracle on Ĥ, field-s | oracle on Ĥᶜ.
- **Joint** — the field-s Bonferroni pair's joint coverage and the share of declared replicates on which both bounds exist, with the unstudentized pair once beside it. The calibrated split's mean γ and correlation appear as diagnostics.
- **Rule contrast with `mdf1`, paired by sim_id:**
  - identification: mean |Ĥ|, sensitivity and PPV from `n_sel`, `n_harm` and `n_true` as the recorder defines them (quote the lines; sensitivity is undefined where `n_true` is 0), and replicates where |Ĥ| grew, stayed or shrank;
  - one-sided coverage: field on Ĥ and unstudentized field on Ĥᶜ under both rules, plus field-s under `effMaxSG`.
- **Regime diagnostics** — p̂(Ĥ) mean and its share below 0.5; SD(β̃ᶜ) / mean naive SEᶜ; mean λ-SDᶜ / naive SEᶜ for field and for field-s.
- **Display** — `fs_sim_bias_coverage(scale = "identity")` for both blocks. The Ĥᶜ field-s points come from a copy of the results with `fld_Hc_*_s` renamed to `fld_Hc_*`, labelled field-s.
- **Null cell** — labelled per D2 in every table.

Every metric definition is taken from `mdf1`'s summary and quoted. Render, then commit the `.qmd`, the HTML and the figures.

### 3.2 The extract

The summary document writes `<dir>/md_field_metrics.csv` and `<dir>/COLUMNS_md_field.md` from the same objects its tables print.
- **Rows:** one per cell × block (H, Hc, joint) × estimator × metric.
- **Columns:** `campaign`, `cell`, `block`, `estimator`, `metric`, `tau` (bound-location rows), `value`, `mc_se`, `wilson_lo` and `wilson_hi` (proportions), `n` (the denominator), `commit` (the commit that added the bundles).
- **Monte Carlo SEs:**
  - proportions: √(p(1−p)/n);
  - bias: SD/√n;
  - empirical SD: SD/√(2(n−1));
  - mean SE: SD(SE)/√n;
  - SE/SD: the delta method.
- **`COLUMNS_md_field.md`** defines every column and metric, maps each estimator code to its plain-language label, and states the scale convention.

### 3.3 The record

`<dir>/REPORT_md_field_rerun_2026-09-15.md` contains:
- provenance and the Gate 0 dispositions;
- Gate 1 and Gate 2 results in brief;
- each §3.1 table, pasted from the render and followed by a plain-language reading of two to four bullets, with numbers from the CSV and bounds read by location;
- this scope sentence, verbatim: "These are operating characteristics on one continuous design. They do not verify condition (A3) on the GLM paths, and no construction is promoted on this design's performance.";
- findings and commits.

### 3.4 Directory catalog

Create `<dir>/current_status.md` by transplanting `quarto/simulations/gbsg_020/`'s generator, curated file and check script (`current_status_regen.R`, `status_curated.md`, `check_current_status.sh`; locate them with `git ls-files`). The curated content covers:
- **Campaigns:** `mdf1` and `mdsgnb20`, with payload locations.
- **Reading conventions:** the harm-oriented scale, bound location, and the null cell's truth.
- **Superseded:** `mdf1`'s unstudentized complement bound as the reported complement product.
- **Open work:** the applied document's field-s rows.

Commit the generator, curated file and check script, then regenerate `current_status.md` as the last commit.

### 3.5 Closeout

- Delete the smoke and calibration outputs listed in the Stage 1 record, and nothing else.
- Copy `md_field_metrics.csv` and `COLUMNS_md_field.md` to `~/Downloads`.
- **Post-conditions**, printed in the closing message:
  - no tracked modification;
  - `git diff --quiet <§1.1 HEAD>..HEAD -- R/` succeeds;
  - the installed `Built` equals §1.2's;
  - the four `mdsgnb20` combined bundles are tracked, with 2,000 rows each;
  - every number in the record's tables is found in the CSV at the printed precision;
  - the transplanted catalog check passes;
  - the `~/Downloads` copies are byte-identical (`cmp`) to the committed files.
- **Closing message:** the commit range to push, the post-conditions, the findings, and the paths of the record, the summary HTML and the CSV. Then stop.

## Out of scope

- `R/`.
- The applied document (it gets its own task).
- DINA and GRF.
- Any threshold, search or MR setting other than §1.3's rule.
- The excluded constructions.
- Pushing.
