# TASK — ACTG175 binary (OR) simulation study under the current constructions: Stage 0, read-only (v2)

**File:** `dev/tasks/TASK_actg175_binary_stage0_2026-09-17_v2.md` · **Issued:** 2026-09-17 by chat, on Larry's direction to start with the binary ACTG175 simulations described in the post-selection supplement · **Supersedes** `TASK_actg175_binary_stage0_2026-09-17.md`, which is not to be run
**Machine:** `pop-os` · **Repo:** `~/Documents/GitHub/forestsearch` · **Branch:** `feature/glm-extension`
**Runs after:** `TASK_md_dina_campaign_2026-09-17.md` has closed out or stopped
**Transport:** `~/Downloads` → `dev/tasks/`, committed first and alone (§1)
**Manuscript copy (read-only):** `dev/reference/post_selection/manuscript_jrssb_initialsubmit/`, written `<ms>`
**Record:** `REPORT_actg175_binary_stage0_2026-09-17.md`, at the path §2 chooses

**What this is.** The supplement's section S8 ("Binary outcomes on the odds-ratio scale", Figures S9 and S10) reports a binary study. Its design, as the supplement states it:
- **Data:** a data-generating mechanism on ACTG175 covariates, arms ZDV+ddI vs ddI, outcome "no CD4 improvement".
- **Subgroup:** H = {wtkg > 70th percentile} ∩ {cd40 > 70th percentile}.
- **Truths:** marginal OR θ†(H) = 0.75, a protective region, against θ†(Hᶜ) ≈ 0.66; CDE θ‡(H) = 0.73, θ‡(Hᶜ) = 0.63.
- **Detectors:** consistency (FS), DINA and GRF, at n = 500–2000.
- **Estimators:** MR and naive, with the oracle in the text; 5,000 multiplier draws; targets the conditional subgroup OR and the marginal OR; summaries conditional on detection.
- **Two internal conflicts:** the text says 500 simulations per cell while the Figure S9/S10 captions say 1,000, and the captions print the MR draws as "NULL".

Larry's direction is to start with this study and update it analogously to the continuous re-runs:
- `effMaxSG` at ε = 0.20, with the survival campaigns' `selection_rule`;
- the field on Ĥ, field-s on Ĥᶜ and their Bonferroni pair, beside unadjusted, oracle and IJ two-term;
- FS, GRF and DINA, with GRF and DINA labelled conditional on the proposed family;
- summaries and extracts.

This stage establishes what the study ran and what a re-run needs. It decides nothing.

## ⚠ CATEGORY

- **Read-only.**
  - No `R/`, template, script, document or payload edit; no install; no render; no campaign.
  - `<ms>` is never modified.
- **Allowed computation:** §7's data check and fits, under their condition, with a 1.5 h ceiling.
- **Writes:** this task document and the record.
  - If the record lands in a directory that has a `current_status.md` generator, regenerate the catalog as §8 allows.
  - If the directory is new, catalog creation waits for the first campaign's closeout.
- **Unattended.** A gate stops on failure, never to ask. A statement here that does not hold is a finding, not a stop.

## Conventions

1. Verify from source: quote `path:line` at a stated commit. Read other branches with `git show <ref>:<path>`; never check them out.
2. `git add` by explicit path only. Untracked files are never staged. No `fetch`, `pull` or `push`.
3. Numbers in the record are computed and pasted, never typed. Nested timing columns are never summed.
4. `<ms>` and any `fs-post-selection` clone are read for design facts only. Nothing here assesses that manuscript or proposes changes to it.
5. No significance language. Excluded constructions: IJ winner-only and winner-floor, κ variants, `field_uniform`.
6. Every R process runs with `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1`.

---

## 1. Provenance and first commit — GATE

```bash
cd ~/Documents/GitHub/forestsearch
hostname; git branch --show-current; git rev-parse --short HEAD; git log --oneline -5
git status --porcelain --untracked-files=no
git cat-file -e HEAD:dev/tasks/TASK_md_dina_campaign_2026-09-17.md && echo "DINA task started"
ps -eo pid,etime,args | grep -E '[e]xec/R|[R]script|[q]uarto' | head
ls -la dev/reference/post_selection/manuscript_jrssb_initialsubmit | head -40
test -d ~/Documents/GitHub/fs-post-selection/.git && git -C ~/Documents/GitHub/fs-post-selection log --oneline -1 || echo "fs-post-selection not cloned here"
Rscript -e 'cat(packageDescription("forestsearch")$Built, "\n")'
```

*GATE:*
- the host is `pop-os` and the branch is `feature/glm-extension`;
- the DINA task document is at HEAD;
- there are no tracked modifications;
- no R, Rscript or quarto process is running, meaning the DINA campaign is finished or stopped;
- `<ms>` exists.

Record in the record how the DINA task ended, from its record or halt file. Copy this document from `~/Downloads` to `dev/tasks/TASK_actg175_binary_stage0_2026-09-17_v2.md` (exact name, else the single match for `*TASK_actg175_binary_stage0_2026-09-17_v2*.md`) and commit it alone.

## 2. S0.1 — Where the binary study lives

1. **In `<ms>`:** locate the supplement source that renders S8's binary paragraph and Figures S9 and S10. Quote the includes and chunks, the fragment files (earlier sessions named `_sim_mr_coverage_or075_actg175_*.qmd`), the payload paths they read, and the README entry for those figure labels.
2. **Across forestsearch** (all refs: `git log --all --name-only`, `git grep` across refs for the payload names, `wtkg` with `cd40` quantile cuts, `or075`, `or15`), and in `fs-post-selection` if cloned: locate the sweep driver (earlier sessions named `mr_coverage_sweep_or15.qmd` and an OR-0.75 variant), the DGM build (`build_actg175_glm_dgm.R` or inline), the payload files, and any records.
3. **Report one table**, one row per file:
   - repo, branch or ref, path;
   - last commit and date;
   - tracked or not, and present on this machine;
   - role: driver, DGM, fragment, payload, README or record.
4. **State which files produced Figures S9 and S10.**
5. **Record location:** the committed forestsearch simulation directory holding this study's drivers or records, if §2.2 finds one; otherwise `quarto/simulations/actg175/binary/`.

## 3. S0.2 — The design as committed

From the driver, DGM code, fragments and payload, quote:
- **Data and subgroup:** covariates and their coding (factor or numeric binary indicators, and any coercion loop); arms and treatment coding; outcome and `adverse_outcome`; H's definition; calibration (`k_inter`); prevalence; super-population size.
- **Truths:** θ†(H), θ†(Hᶜ), θ‡(H), θ‡(Hᶜ), and the conditional target column (earlier sessions named `C_betaHhat`).
- **Designs in the payloads:** the OR-0.75 protective grid; any OR-1.5 harm or null branch, stated only as present or absent.
- **The n grid,** as run.
- **Replicates per cell actually run,** from payload row counts. This settles 500 against 1,000.
- **MR draws and multiplier,** from payload `meta` or the driver. This settles "NULL".
- **Seeds and RNG kind.**
- **Settings per detector:** `sg_focus`, `effect_neighborhood`, `selection_rule`, thresholds on the OR scale, cut grid, `maxk`, `n.min`, `consistency_method`, and the DINA and GRF arguments.
- **Estimators recorded:** naive, MR, oracle (prefix `ora_`), full bootstrap if any.
- **Detection and targets:** the detection definition, and the target of each bias and coverage figure.
- **Per-replicate cost** by detector and n, from recorded timing columns, with the machine and workers used.

## 4. S0.3 — The current package on the binary path

From source at HEAD, quote:
- **The field for a binomial GLM:**
  - the per-candidate pieces for logistic outcomes;
  - the scale of the bounds (log-OR, reported as OR);
  - any outcome-type guard;
  - MR's `consistency_method` requirement.
- **Field-s and `joint_s`:** whether they are returned for binary with the same fields as continuous.
- **`adverse_outcome` for binary:** the default, and how FS, DINA and GRF orient their surfaces under each value, after the P2 fix.
- **Factor covariates:** GRF membership after the P1 fix, and DINA's coercion.
- **The `effMaxSG` band for binary:**
  - the expression and its scale, with `effect_log_scale` where it applies;
  - the identifier side and MR's re-selection side, and whether they are aligned.
- **`fs_sim_bias_coverage()` on the OR scale:** the scale argument, and whether it reads field and field-s columns directly or needs the renamed copy used on the continuous path.

## 5. S0.4 — What a re-run needs

- **Mapping.** Compare the committed MD campaign machinery against the binary driver: `quarto/simulations/actg175/continuous/sim_fs_maxeffCons_mr_field_md_template.qmd` with its identifier knob and field-s recorder, `scripts_mdgrf/`, `scripts_mddina/` if committed, the summaries and the extracts.
- **The gaps to list:**
  - the outcome family and DGM builder;
  - θ† and θ‡ targets beside β(Ĥ);
  - the oracle prefix (`ora_` vs `or_`);
  - recorder columns;
  - the bound-location threshold scale.
- **Not made.** List these gaps only; make no changes.
- **Harm and null designs.** Whether the DGM builder can produce them, and with which arguments. State this as a fact, not a recommendation.

## 6. S0.5 — Cost references

Quote, without running anything:
- this study's recorded costs;
- per-replicate costs and walls from `REPORT_md_field_rerun_2026-09-15.md`, `REPORT_md_grf_2026-09-16.md` and, if committed, the DINA campaign's record.

## 7. S0.6 — Data check and one replicate per detector, conditional

**Condition:**
- the installed `Built` postdates the last commit touching `R/`, and a namespace comparison against HEAD's source finds no differences (closures compared after `utils::removeSource()`, source extracted to a temporary directory);
- and §2 located the DGM recipe and the committed payload.

If either fails, skip this section and say which.

1. **Data check — GATE for the fits.**
   - In `$(mktemp -d)`, rebuild the OR-0.75 design at the grid's smallest and largest n, sim_id 1, with the committed seeds and RNG kind.
   - Compare against the payload's sim_id 1 rows: the data-level columns present (the oracle estimate, n in H, the targets) must match within 1e-8 relative.
   - If they don't, the recipe isn't faithful: skip the fits and record the difference.
2. **Fits.**
   - **What runs:** on those two replicates, `forestsearch()` for FS, GRF and DINA — six fits — on one worker each.
   - **Settings:**
     - `sg_focus = "effMaxSG"` and `effect_neighborhood = 0.20`;
     - `selection_rule` as a committed `p12x20` or `mdsgnb20` bundle's `meta` records;
     - the committed thresholds, cut grid and detector arguments;
     - MR on with `ci_method = "field"`, `draws = 5000L`, `include_complement = TRUE`, `field_complement = TRUE`, `field_scale_complement = "selected"` and `return_reselection = TRUE`.
   - **Limits:** each fit has a 30-min timeout; a timeout is recorded as a cost fact.
   - **Report for each fit:**
     - warnings, verbatim;
     - family and admitted sizes;
     - the selection, its n and its oriented log-OR;
     - whether `field$lower_1s`, `field$complement$upper_1s_s` and the `joint_s` Bonferroni bounds are present and finite;
     - wall time;
     - the committed payload's selection for the same replicate, beside it.
3. Remove the temporary directories.

## 8. Record, catalog, closeout

1. **Record,** at the §2.5 path. It contains:
   - provenance, including how the DINA task ended;
   - S0.1–S0.6;
   - a **Facts for Larry's decisions** block, facts only and no recommendation:
     - the supplement's design as committed: cells, replicates, draws, detector settings, targets;
     - what the package supports on binary today;
     - the rule and settings previously used, beside the current campaign standard;
     - template gaps;
     - measured and recorded cost per replicate by detector and n;
     - harm and null designs as available options.

   Commit the record by explicit path.
2. **Catalog.** If the record's directory already has `current_status.md` and its generator, regenerate it as the last commit and run its check. If not, say so.
3. **Post-conditions,** printed in the closing message:
   - `git diff --name-only <§1 HEAD>..HEAD` lists only this task document, the record and any catalog file;
   - `<ms>` is unchanged;
   - there are no tracked modifications;
   - the installed `Built` is unchanged;
   - the temporary directories are gone.
4. **Closing message:**
   - the commit range to push;
   - where the study lives;
   - replicates and draws as actually run;
   - whether the package runs field, field-s and the pair on binary;
   - the six fits' walls;
   - the findings.

   Then stop.

## Out of scope

- Any edit beyond §8; installs; renders; campaigns.
- Changes to `<ms>` or `fs-post-selection`.
- `R/` proposals.
- Pushing.
