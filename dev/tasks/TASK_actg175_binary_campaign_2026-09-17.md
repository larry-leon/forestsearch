# TASK — ACTG175 binary (OR) design under the current constructions: campaigns `orfs`, `orgrf`, `ordina`, Stages 1–3

**File:** `dev/tasks/TASK_actg175_binary_campaign_2026-09-17.md` · **Issued:** 2026-09-17 by chat, on Larry's decisions below
**Machine:** `pop-os` · **Repo:** `~/Documents/GitHub/forestsearch` · **Branch:** `feature/glm-extension`, HEAD containing `98475a05`
**Transport:** `~/Downloads` → `dev/tasks/`, committed first and alone (§1.1)
**Directory:** `quarto/simulations/actg175/binary_020/`, written `<dir>`, where the study's driver and payloads live
**Campaign tags:** `orfs`, `orgrf`, `ordina`

**References (all committed):**
- `quarto/simulations/actg175/binary/REPORT_actg175_binary_stage0_2026-09-17.md`, written `S0` — the study's design (§3), the package on the binary path (§4), the gaps (§5), the costs (§6), the six fits (§7).
- `<dir>/maxeffCons_mr_coverage_sweep_or075.qmd` (`1d42f6da`), written `the study driver` — the producer of the supplement's Figures S9 and S10, and this task's source for the data recipe, seeds, truths and detector settings.
- `<dir>/mr_sweep/maxeffCons_actg175_or075_seedtab_s1000/` — its 21 cell bundles and grid.
- `quarto/simulations/actg175/binary/mr_coverage_sweep_or10.qmd` and `mr_coverage_sweep_or15.qmd` — the committed OR 1.0 and OR 1.5 drivers, for their `target_or_h` and calibration settings only.
- `quarto/simulations/actg175/continuous/sim_fs_maxeffCons_mr_field_md_template.qmd`, `scripts_mddina/`, `summary_continuous_field_mddina.qmd`, `md_dina_metrics.csv`, `COLUMNS_md_dina.md`, and `dev/tasks/TASK_md_dina_campaign_2026-09-17.md` — the MD machinery this task transplants.

**What this is.** The supplement's binary study ran `maxeffCons` at ε 0.10 with IJ intervals only, so it carries no field columns. This task re-runs that design under the current constructions and adds two design points, then reports the three identifiers together.

## Dispositions (Larry, 2026-09-17)

- **Rule:** `effMaxSG`, `effect_neighborhood = 0.20`, `selection_rule = "neighborhood"`, replacing the study's `maxeffCons` at ε 0.10.
- **Thresholds:** the study's, unchanged and in every cell: `effect.threshold = 0.90`, `consistency.threshold = 0.80`, `pconsistency.threshold = 0.90`, with `adverse_outcome = TRUE` so OR > 1 is harm.
- **Design points,** each built by the study driver's calibration with only `target_effect` changed:
  - **protective**, `target_or_h = 0.75` — the supplement's;
  - **harm**, `target_or_h = 1.5`;
  - **borderline null**, `target_or_h = 1.0` — the planted region at the null against a protective complement, as the committed OR 1.0 driver has it. The driver's homogeneous `dgm_model = "null"` branch is not used.
- **Sizes:** n = 500 and n = 2000, the ends of the study's sweep. Six cells.
- **Replicates:** 2,000 per cell, against the study's 1,000. The seed table holds 5,000.
- **Identifiers:** FS, GRF and DINA, one campaign each over all six cells.
- **Constructions:** the field one-sided lower bound on Ĥ, field-s one-sided upper on Ĥᶜ, their Bonferroni pair, with unadjusted, oracle and IJ two-term as references, on the OR scale.
- **Bound-location ladder:** τ = 0.7, 0.8, 0.9, 1.0, 1.25, 1.5, 2.0. On Ĥ, a lower bound ≥ τ reads "harm of at least τ"; on Ĥᶜ, an upper bound ≤ τ reads "harm of at most τ".
- **Naming:** the recorder follows the MD naming (`nv_*`, `mr_*`, `or_*`, `fld_*`), not the study's (`t2_*`, `ora_*`), so `fs_sim_bias_coverage()` reads it without renames beyond the field-s copy.
- **Gate 1 advance go:** if every Stage 1 gate is green and the §1.7 projection for Stage 2 is under 12 hours, do not stop at Gate 1: note the advance go in the Gate 1 record, then run Stages 2 and 3. Otherwise stop at Gate 1 and report. Never reduce cells, replicates, gates or knobs to fit under the threshold.

**Governing constraint.** GRF's and DINA's candidate families are generated from fitted surfaces, so the fixed-family condition does not hold for them: their coverage is coverage of the estimand conditional on the proposed family. FS's family is the prespecified cut grid. Every table, caption and extract row says which, and comparisons are descriptive.

## ⚠ CATEGORY

- **No `R/` change. No install.** §1.2 asserts `R/` is unchanged since `064fce91` and the build is `2026-09-17 04:47:31 UTC`.
- **New:** a binary campaign template and `<dir>/scripts_or/`, a summary, records, the extract, and the directory's catalog.
- **Not edited:** the study driver, its bundles, and everything under `quarto/simulations/actg175/binary/`.
- **Compute authorized by this kickoff:** Stage 1's smoke (four renders of 20 replicates) and calibration (three renders), ceiling 3 h. Stages 2 and 3 under the advance go.
- **Unattended.** Gates stop on failure, never to ask. On a stop: commit what is green, record the failure, and stop. A statement here that does not hold is a finding.

## Conventions

The MD DINA task's conventions 1–10 apply, with these changes:
- the scale is the odds ratio, computed on log-OR and reported as OR;
- "harm-oriented MD" reads as "OR above 1";
- the ladder is the dispositions'.

---

## Stage 1 — template, smoke, calibration → Gate 1

### 1.1 Provenance and first commit — GATE

Assert: host `pop-os`; branch `feature/glm-extension`; `98475a05` in HEAD; no tracked modifications; no R, Rscript or quarto process running. Then copy this document from `~/Downloads` to `dev/tasks/TASK_actg175_binary_campaign_2026-09-17.md` (exact name, else the single match for `*TASK_actg175_binary_campaign_2026-09-17*.md`) and commit it alone.

### 1.2 Package check — GATE

`git diff --quiet 064fce91..HEAD -- R/ DESCRIPTION NAMESPACE`; `packageDescription("forestsearch")$Built` is `2026-09-17 04:47:31 UTC`, and two doFuture workers report the same.

### 1.3 The template — transplant

Create `<dir>/sim_fs_mr_field_or_template.qmd` by copying the MD template and changing only the named blocks. Quote both sources for every block.

- **Knobs.** The MD template's, renamed `FS_OR_*`, including `FS_OR_METHOD` (`consistency` | `dina` | `grf`), `FS_OR_FOCUS`, `FS_OR_NBHD`, `FS_OR_CI`, `FS_OR_FIELD_SCALEC`, `FS_OR_CAMPAIGN`, `FS_OR_WORKERS`, `FS_OR_START`, `FS_OR_NSIMS`, `FS_OR_MODE`, plus a design knob `FS_OR_TARGET` (`0.75` | `1.0` | `1.5`) and `FS_OR_N`.
- **Data and DGM** from the study driver, verbatim except `target_effect = FS_OR_TARGET`: the trial and arm coding, the outcome `y_neg`, the covariate pools, the factor-to-numeric coercion before every fit and in the evaluation frame, H's quantile cuts, `calibrate_glm_interaction()` with its `k_inter_range`, `grid_step`, `n_super` and `seed`, and the evaluation frame with `eval_seed`.
- **Seeds:** the study's pre-generated table indexed by `sim_id`, and its per-replicate `RNGkind("L'Ecuyer-CMRG")` and `seedit`, verbatim.
- **Identification:** the study's thresholds and search settings, with `sg_focus`, `effect_neighborhood` and `selection_rule` from the knobs; the study's GRF and DINA arguments.
- **MR:** the MD template's argument list — `ci_method = "field"`, `draws = 5000L`, `include_complement = TRUE`, `field_complement = TRUE`, `field_scale_complement = "selected"`, `return_reselection = TRUE`, `ij_residual`, `confirm_rule` — and MR requires `consistency_method = "resample"`, which the study already sets.
- **Recorder:** the MD template's, unchanged in naming, plus the study's targets per replicate: `C_dagger` and `C_ddagger` for both blocks beside `betaHhat_*`, via the study's `fs_attach_betaHhat()` call and its truth table.
- **`meta` and poolability keys:** the MD template's, plus `target_or_h`, the truths, and the thresholds.
- **Stem and results directory:** the MD template's convention, with the design and n in the stem. `.refuse_if_tracked()` stays live, so no committed study path can be written.

Commit the template.

### 1.4 Scripts — transplant

Copy `scripts_mddina/` to `<dir>/scripts_or/`: the sampler, the identity checker, `gate2.R` and the runner, with the OR fields and `meta` items. The runner writes no hard-coded trailer. Commit them.

### 1.5 Smoke — GATE

Four renders of 20 replicates each (`sim_id` 1–20), at `target_or_h = 0.75`, n = 500, tag `orsmoke`:

- **(a) The data recipe.** With the knobs at the study's rule, compare against the committed study bundle `…_mr_n500_res.rds` for `sim_id` 1–20: the truths, the oracle estimates and every other rule-independent data-level column agree within 1e-8 relative. A difference means the port is unfaithful: stop.
- **(b) FS under the campaign rule.** It runs, declares, and its `meta` carries the rule, the thresholds, `target_or_h`, the truths, `field_scale_complement = "selected"`, `pkg_version` and host.
- **(c) GRF and DINA under the campaign rule.** Both run with zero factor-comparison warnings and zero NA-membership candidates; DINA's proposal and admission floors, as applied, are the OR-scale threshold on the harm side.
- **(d) The constructions**, on every declared replicate of (b) and (c): the nine `fld_Hc_*_s` and nine `fld_joint_s_*` columns are finite; `fld_Hc_lo1s_s ≤ fld_Hc_up1s_s`; the Bonferroni harm bound agrees between `joint` and `joint_s` where the draw counts agree; every bound is an OR, so positive; p̂(Ĥ) is recorded.

Report, as facts: the declared counts; each identifier's selection on `sim_id` 1 beside S0 §7's fit; `fit_mr_secs` and field seconds.

### 1.6 Calibration

The most expensive configuration, FS at `target_or_h = 1.5`, n = 2000, at `FS_OR_WORKERS` = 16, 32 and 63, with 3 × workers replicates (tags `orcal16`, `orcal32`, `orcal63`), memory sampled every 5 s. Report render wall, `fit_mr_secs` mean, median and 90th percentile, peak summed RSS, and replicates per minute per worker count. Choose W and name what limits it.

### 1.7 Projection and Gate 1 record

Project Stage 2 at W: three identifiers × six cells × 2,000 replicates, scaling each configuration's per-replicate cost from S0 §7's six fits and this calibration. State the projection, a ceiling of 1.5 × it, and a per-render timeout of 2 × the longest projected batch, at least 20 min.

Write `<dir>/REPORT_actg175_or_stage1_2026-09-17.md` with provenance, the package check, the template and script transplants with quoted sources, the smoke results, the calibration table, W, the projection, the ceiling and the advance-go decision. Commit it, then apply the advance go or stop.

---

## Stage 2 — campaigns `orfs`, `orgrf`, `ordina`

Run the MD DINA task's §2.1–§2.3 mechanics, with these changes.

- **Order:** `orfs` first, then `orgrf`, then `ordina`; within each, the six cells in the order 0.75/500, 0.75/2000, 1.5/500, 1.5/2000, 1.0/500, 1.0/2000.
- **Run files:** `<dir>/LOG_or_progress.txt`, `<dir>/HALT_or.md`, raw logs untracked under `<dir>/logs_or/`, Gate 2 record `<dir>/REPORT_actg175_or_gate2_2026-09-17.md`.
- **Per-cell Gate 2:**
  - 2,000 rows with `sim_id` exactly 1–2000, and the batch files matching the combined bundle on every column;
  - `meta` carries the rule, thresholds, `target_or_h`, the truths, `field_scale_complement`, `pkg_version` and host;
  - **same draws across identifiers:** for the two later campaigns, the data-level columns — the truths, the oracle columns and the seed — are identical to `orfs`'s in the same cell, in both directions, within 1e-8 relative;
  - the §1.5(d) checks on every declared replicate;
  - MR failures on declared replicates at most 40.
- **Pass, failure and ceiling rules:** as in the MD DINA task.

---

## Stage 3 — summary, extract, record, catalog, closeout

### 3.1 Summary

`<dir>/summary_actg175_or.qmd`, transplanted from `summary_continuous_field_mddina.qmd`, reading the eighteen combined bundles. Changes: `fs_sim_bias_coverage(scale = "log")`; the ladder of the dispositions; bounds and biases reported as ORs, with bias in SD units on the log-OR scale; the conditional-family sentence at the top and on every GRF and DINA caption.

Tables, each cell labelled by its design and n:
- **Ĥ and Ĥᶜ**, per identifier: declaration rate, bias, SE/SD, one-sided and two-sided coverage, for unadjusted, oracle, IJ two-term, and field on Ĥ / field-s on Ĥᶜ with the unstudentized field once beside it.
- **Bound location** on the ladder, with Monte Carlo SEs: field | oracle on Ĥ, field-s | oracle on Ĥᶜ.
- **Joint:** the field-s pair's joint coverage and the share of declared replicates carrying both bounds, the unstudentized pair once beside, and γ and correlation as diagnostics.
- **Three identifiers:** one table over the six cells — declaration rate, field lower coverage on Ĥ, field-s upper coverage on Ĥᶜ, Bonferroni joint, mean |Ĥ|, sensitivity and PPV.
- **Against the supplement's study**, quoted from the committed grid, not recomputed: detection rate and the naive and IJ coverage of θ† at n = 500 and 2000 under the study's rule, beside this campaign's at the same cells. State that the rule, replicates and constructions differ.
- **Regime diagnostics:** p̂(Ĥ) and its share below 0.5; SD(β̃ᶜ)/naive SEᶜ; λ-SDᶜ/naive SEᶜ for field and field-s.
- **Targets:** every coverage row names its target — β(Ĥ) for the conditional rows, with θ† and θ‡ reported beside for the harm block, as the study did.

Render; commit the `.qmd`, HTML and figures.

### 3.2 Extract

`<dir>/or_metrics.csv` and `<dir>/COLUMNS_or.md`, in `md_dina_metrics.csv`'s schema with its `identifier` column and a `design` column (`or075`, `or150`, `or100`). Monte Carlo SEs and Wilson intervals as in `COLUMNS_md_dina.md`. The COLUMNS file states the OR scale, the ladder's reading and the conditional-family reading.

### 3.3 Record

`<dir>/REPORT_actg175_or_2026-09-17.md`, in the MD records' structure: provenance and dispositions, Gate 1 and Gate 2 in brief, each §3.1 table with a two-to-four-bullet reading by bound location, findings and commits.

Scope sentence, verbatim: "These are operating characteristics on one binary design at three effect sizes, with GRF's and DINA's figures conditional on the proposed family. They do not verify condition (A3), no construction is promoted on this design's performance, and they supersede nothing in the committed study, which ran a different rule and different constructions."

### 3.4 Catalog

Create `<dir>/current_status.md` by transplanting `quarto/simulations/actg175/continuous/`'s generator, curated file and check script (in `scripts_mdsgnb20/`). The curated content covers the committed study and its payloads, the three new campaigns, the OR reading conventions, and open work. Regenerate as the last commit; the check must pass.

### 3.5 Closeout

Delete the Stage 1 smoke and calibration outputs and nothing else. Write `~/Downloads/bundle_or_2026-09-17.zip` from HEAD with `or_metrics.csv`, `COLUMNS_or.md`, the Stage 1 record, the Gate 2 record and the record. Post-conditions: no tracked modifications; `git diff --quiet <§1.1 HEAD>..HEAD -- R/`; nothing under `quarto/simulations/actg175/binary/` changed; the eighteen combined bundles tracked with 2,000 rows each; every number in the record's tables found in the CSV at printed precision; the catalog check passes; the bundle listed. Closing message: the commit range to push, the post-conditions, the three-identifier table, the findings, and the bundle path.

## Out of scope

- `R/`; the applied documents; the committed study's files; the MD campaigns.
- Any threshold, search or MR setting beyond the dispositions.
- The excluded constructions; pushing.
