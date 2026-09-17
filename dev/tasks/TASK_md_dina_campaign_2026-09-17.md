# TASK — DINA on the ACTG175 continuous (MD) design: campaign `mddina`, Stages 1–3

**File:** `dev/tasks/TASK_md_dina_campaign_2026-09-17.md` · **Issued:** 2026-09-17 by chat, on Larry's yes to the sequence "P1 and P2, then the GRF campaign, then the DINA campaign"
**Machine:** `pop-os` · **Repo:** `~/Documents/GitHub/forestsearch` · **Branch:** `feature/glm-extension`, HEAD containing `853c43c6`
**Transport:** `~/Downloads` → `dev/tasks/`, committed first and alone (§1.1)
**Directory:** `quarto/simulations/actg175/continuous/`, written `<dir>` below
**Campaign tag:** `mddina`

**References (all committed):**
- `dev/tasks/TASK_md_grf_2026-09-16.md` ("the GRF task") and `dev/tasks/TASK_md_grf_resume_2026-09-16.md` — the mechanics this task transplants.
- `<dir>/REPORT_grf_dina_fixes_2026-09-16.md` ("the fix record") — P2 landed; F5 and F6 are the DINA fits on md40 n500, sim_id 1.
- `<dir>/REPORT_md_dina_grf_stage0_2026-09-16.md` ("S0").
- `<dir>/REPORT_md_grf_2026-09-16.md`, `scripts_mdgrf/`, `summary_continuous_field_mdgrf.qmd`, `md_grf_metrics.csv`, `COLUMNS_md_grf.md` — the GRF campaign.
- `md_field_metrics.csv` — the FS campaign `mdsgnb20`.

**What this is.** Campaigns `mdsgnb20` (FS) and `mdgrf` (GRF) evaluated the MR intervals on the ACTG175 continuous design under `effMaxSG`, ε = 0.20. This task runs the DINA identifier on the same four cells, replicates and seeds, now that its proposal floor is oriented (P2).

It records the same constructions:
- unadjusted, oracle and IJ two-term;
- the field on Ĥ;
- field-s on Ĥᶜ, with the unstudentized complement field beside it;
- the Bonferroni pair;
- p̂(Ĥ).

It pairs with `mdsgnb20` at the data level (`n_true` and the oracle columns). It ends with a three-identifier comparison.

**Governing constraint.** DINA's candidate family is generated from a fitted surface, so the fixed-family condition does not hold. Every DINA coverage figure is coverage of the estimand conditional on the proposed family, and every table, caption and extract row says so. Comparisons with FS and GRF are descriptive, not a contest, and state the confound.

## Dispositions (Larry, 2026-09-17)

- **Identifier:** DINA only, through the MD template's committed identifier knob (commit `894da993`; quote its name and its DINA value). The call must run DINA's identifier path, the one P2 fixed, and not the `use_dina` screening path under FS, whose floor is still unoriented (§1.3).
- **Rule:** `effMaxSG`, `effect_neighborhood = 0.20`, and `selection_rule` as `mdsgnb20`'s `meta` records.
- **DINA floors:** the proposal floor and the admission floor as the package applies them after P2. Both sit on the harm-oriented MD scale at the effect threshold, 30. The smoke checks them as applied (§1.6(c)).
- **Other DINA arguments:** as the survival campaign `dinamr` passed them, transplanted and quoted. A survival-only argument is omitted, with the reason stated.
- **MR:** exactly as the GRF task's dispositions: `ci_method = "field"` passed explicitly, with the same argument set. Gate 0 stops on a hard-coded difference or an unaligned re-selection for DINA (§1.3).
- **Labelling:** the conditional-family sentence on every table and figure, and the confound sentence once at the top of the summary and the record (§3.1).
- **Two carried fixes to committed files, not in `R/`:**
  - the runner's commit messages carry no hard-coded model or co-author trailer;
  - the directory's catalog generator gains inventory rows for `mdgrf` and `mddina`.
- **Gate 1 advance go (Larry offline):**
  - If every Stage 1 gate is green and the §1.7 projection for Stage 2 is under 8 hours, do not stop at Gate 1. Note the advance go and its condition in the Gate 1 record, run Stage 2 at the worker count chosen in §1.7 with the record's ceiling and per-render timeout, and run Stage 3 once Stage 2 is green.
  - Otherwise, stop at Gate 1 and report.
  - Never reduce cells, replicates, gates or knobs to fit under the threshold.

## ⚠ CATEGORY

- **No `R/` change. No install.**
- **Edited:**
  - the MD template: a DINA argument block, any DINA recorder fields not yet present, and `meta` and poolability keys;
  - the catalog generator.
- **New:** `scripts_mddina/`, a summary document, records, the extract.
- **Compute authorized by this kickoff:** Stage 1's smoke (20 FS, 20 GRF and 20 DINA replicates) and calibration (333 DINA replicates), with a 2 h ceiling. Stages 2 and 3 run under the advance go.
- **Unattended.**
  - Gates stop on failure, never to ask.
  - On a stop, commit what is green, record the failure, and stop.
  - A statement here that does not hold is a finding, not a stop.

## Conventions

The GRF task's conventions 1–10 apply unchanged, with "GRF" read as "DINA" where they name the identifier.

---

## Stage 1 — edits, smoke, calibration → Gate 1

### 1.1 Provenance and first commit — GATE

```bash
cd ~/Documents/GitHub/forestsearch
hostname; git branch --show-current; git rev-parse --short HEAD; git log --oneline -3
git status --porcelain --untracked-files=no
git merge-base --is-ancestor 853c43c6 HEAD && echo "mdgrf closeout in HEAD"
ps -eo pid,etime,args | grep -E '[e]xec/R|[R]script|[q]uarto' | head
```

*GATE:*
- the host is `pop-os` and the branch is `feature/glm-extension`;
- `853c43c6` is in HEAD;
- there are no tracked modifications;
- no R, Rscript or quarto process is running.

Copy this document from `~/Downloads` to `dev/tasks/TASK_md_dina_campaign_2026-09-17.md` (exact name, or else the single match for `*TASK_md_dina_campaign_2026-09-17*.md`) and commit it alone.

### 1.2 Package check — GATE

- `git diff --quiet 019be60f..HEAD -- R/ DESCRIPTION NAMESPACE` succeeds.
- `packageDescription("forestsearch")$Built` is `2026-09-17 04:47:31 UTC`, and two doFuture workers report the same.

*GATE:* both hold.

### 1.3 Gate 0 from the records — GATE

Quote each item from its source, then apply the stop rules.

1. **Identifier path.**
   - Quote the committed template's knob lines that map the DINA value onto the `forestsearch()` arguments, and the call F5 used in the fix record.
   - Show from source that the DINA value runs DINA's identifier path — the lines P2 changed — and not the `use_dina` screening path under FS.
   - If it runs the screening path, stop.
2. **Floors.** Quote the fix record's F5 after-P2 floors, proposal and admission, as applied.
3. **Forwarding (S0.4).**
   - For each MR argument in the dispositions, state whether the DINA path forwards it or hard-codes it.
   - Stop on a hard-coded value that differs from the FS branch's, naming the argument and line.
4. **Re-selection (S0.4).**
   - State whether MR's re-selection is aligned for DINA on this outcome, quoting both sides.
   - If it is not aligned, stop.
5. **Recorder fields.** List the DINA recorder fields from S0.5 that the template still lacks after the GRF task's edit E3.

### 1.4 The survival reference — read only

Build the GRF task's §1.4 table for DINA, from:
- the survival template at HEAD;
- the committed `dinamr` runner and campaign scripts (locate them with `git ls-files`);
- a `dinamr` bundle's `meta`.

The rows are:
- DINA's selection through the knob;
- every DINA argument;
- how DINA uses the band;
- focus, band and rule;
- the MR arguments as passed on the DINA path.

List every difference as a finding.

**Rule:** the campaign takes `dinamr`'s DINA arguments; everything else stays the MD template's.

### 1.5 Template and generator edits — transplant, no `R/`

Edit `<dir>/sim_fs_maxeffCons_mr_field_md_template.qmd` bottom-up. Copy each block from the survival template and quote its source lines.

- **E1 — DINA argument block.** Apply the §1.4 arguments only when the knob selects DINA. Add no knob for the floors: they follow the effect threshold as the package applies them.
- **E2 — recorder.** Add the §1.3.5 fields and their fills, transplanted.
- **E3 — `meta` and poolability keys.** They must carry the identifier and the DINA arguments.

Then check the default path: `git diff` shows only additions guarded by the DINA value, plus the new recorder fields. Commit the template.

In `<dir>/current_status_regen.R`, add inventory rows for `mdgrf` and `mddina`, copied from the `mdsgnb20` rows with the tag and paths changed. Commit it.

### 1.6 Smoke — GATE

Transplant `scripts_mdgrf/` to `scripts_mddina/`: the sampler, the identity checker, `gate2.R` and the runner.
- Add the DINA fields and `meta` items to the checker and to `gate2.R`.
- Remove the runner's hard-coded commit-message trailer.

Commit the scripts.

**(a) FS regression.** md40 n500, sim_id 1–20, knob at its default, tag `mddinasmokefs`.
- Compare every column except `*_secs` against `mdsgnb20`'s combined bundle.
- Values must agree within 1e-8 relative, with label ties enumerated and zero selection flips.

**(b) GRF regression.** md40 n500, sim_id 1–20, knob at its GRF value with `mdgrf`'s knobs, tag `mddinasmokegrf`.
- Compare every column except `*_secs` against `mdgrf`'s combined bundle, with the same tolerance, enumeration and zero flips.

**(c) DINA.** md40 n500, sim_id 1–20, knob at its DINA value, tag `mddinasmoke`, 20 workers. Check that:
- it runs without error, with warnings captured verbatim and zero factor-comparison warnings;
- `n_true` is identical and the oracle columns agree within 1e-8 of `mdsgnb20`'s on all 20 sim_ids;
- `meta` carries:
  - identifier DINA and its arguments;
  - focus `effMaxSG`, band 0.20, and the rule;
  - `ci_method` field and `field_scale_complement` selected;
  - `pkg_version` 0.3.5 and host `pop-os`;
- on sim_id 1, the proposal floor and admission floor as applied both sit at 30 on the harm-oriented scale. Read them from the returned object, or from a direct identification fit on that replicate in a temporary script;
- every proposed candidate satisfies oriented `tau_hat ≥ 30`;
- the E2 fields are filled on every declared replicate;
- the nine `fld_Hc_*_s` and nine `fld_joint_s_*` columns are finite wherever the complement block is filled;
- `fld_Hc_lo1s_s ≤ fld_Hc_up1s_s`;
- the Bonferroni harm bound is identical between `joint` and `joint_s` wherever the draw counts agree;
- p̂(Ĥ) is recorded.

Report as facts:
- the declared count;
- sim_id 1's selection, beside the fix record's F5 after-P2 selection. If the campaign's DINA arguments equal F5's, the two selections must agree, and that is a gate;
- `fit_mr_secs` and DINA fit time (mean, median, maximum);
- the proposed and admitted family sizes.

*GATE:* (a), (b) and (c) all hold.

### 1.7 Calibration

Run the GRF task's §1.7 with these changes:
- knob at its DINA value;
- tags `mddinacal16`, `mddinacal32` and `mddinacal63`;
- DINA fit time reported where GRF fit time was.

Choose W by projected wall under the memory limit, and name what limits it. Then state:
- the projection;
- a ceiling of 1.5 × the projection;
- a per-render timeout of 2 × the longest projected batch, and at least 20 min.

### 1.8 Gate 1 record

Write `<dir>/REPORT_md_dina_stage1_2026-09-17.md`, with the contents of the GRF task's §1.8 for this run.

Commit it. The smoke and calibration outputs stay untracked and are listed in the record. Then apply the advance go, or stop with the closing message.

---

## Stage 2 — campaign `mddina`

Run the GRF task's §2.1–§2.3 with these changes:
- **Launch:** knob at its DINA value; `FS_MD_FIELD_SCALEC=selected FS_MD_CI=field FS_MD_CAMPAIGN=mddina`; the rule knobs; `FS_MD_WORKERS=W` and the timeout.
- **Run files:** progress log `<dir>/LOG_mddina_progress.txt`, halt file `<dir>/HALT_mddina.md`, raw logs untracked under `<dir>/logs_mddina/`, and Gate 2 record `<dir>/REPORT_md_dina_gate2_2026-09-17.md`.
- **Per-cell Gate 2:**
  - `meta` as in §1.6(c);
  - same draws as `mdsgnb20`, checked in both directions;
  - the §1.6(c) field-s and recorder checks on every declared replicate;
  - MR failures on declared replicates at most 40.
- **Order, pass, failure and ceiling rules:** as in the GRF task.

---

## Stage 3 — summary, extract, record, catalog, closeout

### 3.1 Summary — transplant

Copy `summary_continuous_field_mdgrf.qmd` to `summary_continuous_field_mddina.qmd`.

**Inputs:**
- the four `mddina` bundles;
- `mdsgnb20`'s four, for pairing;
- FS and GRF comparator rows copied from the committed `md_field_metrics.csv` and `md_grf_metrics.csv`, not recomputed.

**Opening statement.** At the top, once, verbatim:

"DINA's and GRF's candidate families are generated from fitted surfaces, so the fixed-family condition does not hold for them: every DINA and GRF coverage figure below is coverage of the estimand conditional on the proposed family. FS's family is the prespecified cut grid. The three identifiers are compared descriptively, not as a contest: the identifier, the family construction and the set of detected replicates each summary conditions on all differ."

**Captions.** Every DINA table and figure caption ends with "(conditional on the proposed family)".

**Tables:** those of the GRF summary, for DINA, plus one **three-identifier table**. Its rows are the cells; for each identifier (FS, GRF, DINA) it gives:
- declaration rate;
- field lower-bound coverage on Ĥ;
- field-s upper-bound coverage on Ĥᶜ;
- Bonferroni joint coverage;
- mean |Ĥ|, taken from `mean_n_sel` for FS;
- sensitivity and PPV.

**Null cell:** labelled by its truth in every table.

Render, then commit the `.qmd`, the HTML and the figures.

### 3.2 Extract

Write `<dir>/md_dina_metrics.csv` and `<dir>/COLUMNS_md_dina.md` with the GRF extract's schema, including its `identifier` column (`dina`).
- FS and GRF comparator rows are copied from the committed CSVs, with their commits.
- The COLUMNS file states the conditional-family reading, and that FS's `mean_n_harm` equals its `mean_n_sel`.

### 3.3 Record

Write `<dir>/REPORT_md_dina_2026-09-17.md` as the GRF record's structure, including the three-identifier table and its reading.

The scope sentence, verbatim: "These are operating characteristics of the DINA identifier on one continuous design, conditional on the proposed family. They do not verify condition (A3), no construction is promoted on this design's performance, and DINA ran with the proposal-floor orientation fix landed in the commits recorded in REPORT_grf_dina_fixes_2026-09-16.md."

### 3.4 Catalog

In `<dir>/status_curated.md`:
- add campaign `mddina`, with its payload locations and the conditional reading convention;
- replace DINA's open-work line with "DINA on this design: campaign `mddina` complete";
- add under open package work, as not blocking: "the `use_dina` screening path under FS still applies the unoriented floor".

Commit it. Then regenerate `current_status.md` as the last commit; `check_current_status.sh --commit` must pass, and the regenerated inventory must list `mdsgnb20`, `mdgrf` and `mddina`.

### 3.5 Closeout

1. Delete the Stage 1 smoke and calibration outputs, and nothing else.
2. Copy `md_dina_metrics.csv`, `COLUMNS_md_dina.md` and the record to `~/Downloads`.
3. Check the post-conditions:
   - no tracked modifications;
   - `git diff --quiet <§1.1 HEAD>..HEAD -- R/` succeeds;
   - `Built` is unchanged;
   - the four `mddina` combined bundles are tracked, with 2,000 rows each;
   - every number in the record's tables is found in the CSV at printed precision;
   - the catalog check passes;
   - the `~/Downloads` copies are `cmp`-identical to the committed files.
4. **Closing message:**
   - the commit range to push;
   - the post-conditions;
   - the three-identifier table;
   - the findings;
   - the paths of the record, the summary HTML and the CSV.

   Then stop.

## Out of scope

- `R/`, including the `use_dina` screening path and the survival GRF forest matrix.
- The applied documents.
- Re-running `mdsgnb20` or `mdgrf`.
- Any threshold, search or MR setting beyond the dispositions.
- The excluded constructions.
- Pushing.
