# TASK — ACTG175 continuous (MD) re-run under the current field constructions: Stage 0, read-only

**File:** `dev/tasks/TASK_md_field_rerun_stage0_2026-09-15.md` · **Issued:** 2026-09-15 by chat, commissioned by Larry
**Repo:** `~/Documents/GitHub/forestsearch` · **Machine:** whichever this is run on (the record states it) · **Branch:** as checked out; never switched
**Transport:** `~/Downloads` → `dev/tasks/`, committed first and alone (§1)
**Record:** `quarto/simulations/actg175/continuous/REPORT_md_field_rerun_stage0_2026-09-15.md`, beside `mdf1`'s REPORT

**What this is.** Campaign `mdf1` (2026-09-07) evaluated the MR intervals on the ACTG175 continuous design before the studentized complement field (field-s) existed and before the MR defaults moved. A re-run is planned under the current constructions: the field one-sided lower bound on β(Ĥ), the field-s one-sided upper bound on β(Ĥᶜ), and their Bonferroni pair, with unadjusted, oracle and IJ two-term as references, paired replicate by replicate with `mdf1`. This stage establishes, from source and committed records, what that re-run would execute, save, reproduce and cost, so Larry can take decisions D1–D5 at Gate 0 (§10). It decides nothing and runs nothing.

---

## ⚠ CATEGORY

- **No `R/` change.** No edit to any template, script, document, payload or bundle. **No install, render, simulation or timing run**; no package function is called on data.
- **Writes:** this task document, the record, and `current_status.md` only as §10.2 allows.
- R is used only to read namespaces, formals and function bodies (§2.6, §3), and to summarise committed `mdf1` bundles with base R where §7.3 allows.
- **Unattended.** A gate stops on failure, never to ask. A claim below that does not hold is a **finding** for the record, never a stop. Estimate: under an hour.

## Conventions

1. **Verify from source.** Quote `path:line` at a stated commit. Read refs other than HEAD with `git show <ref>:<path>`. If a file to be read has tracked modifications, read `HEAD:<path>` and note the modification.
2. **Leave the machine as found.** No branch switch, `fetch`, `pull`, `push` or install, and nothing sent to a running process. A campaign running here is reported and left alone.
3. `git add` by explicit path only; untracked files are never staged.
4. **Quote, don't recompute.** `mdf1`'s figures are quoted from its committed record with line numbers. A figure is computed from its bundles only where §7.3 asks for one the record lacks, and is labelled as computed here.
5. **Timing.** State which per-replicate timing columns nest inside others; never sum nested columns.
6. **Excluded everywhere:** IJ winner-only and winner-floor, κ variants, `field_uniform`.
7. Bounds are described by location against thresholds; no significance language.
8. The submitted parent paper is not a topic.

---

## 1. Provenance and first commit — GATE

```bash
cd ~/Documents/GitHub/forestsearch
hostname; pwd; git branch --show-current; git rev-parse --short HEAD
git status --porcelain --untracked-files=no; git log --oneline -5
echo "behind ahead (vs upstream, as of last fetch):"; git rev-list --left-right --count @{u}...HEAD 2>/dev/null || echo "no upstream"
ls -l .git/FETCH_HEAD
uptime; ps -eo pid,etime,args | grep -E '[e]xec/R|[R]script|[q]uarto' | head -20
getconf _NPROCESSORS_ONLN; free -g 2>/dev/null || sysctl -n hw.memsize
Rscript -e 'd <- packageDescription("forestsearch"); si <- sessionInfo(); cat(d$Version, d$Built, find.package("forestsearch"), R.version.string, si$BLAS, si$LAPACK, sep = "\n")'
```

*GATE:* stop if `dev/tasks/TASK_md_field_rerun_stage0_2026-09-15.md` or the record already exists at HEAD, since this stage has then already run. Everything else is recorded as found.

Copy this document from `~/Downloads` (exact name; failing that, the single file matching `*TASK_md_field_rerun_stage0_2026-09-15*.md`; none or several: stop) to `dev/tasks/TASK_md_field_rerun_stage0_2026-09-15.md` and commit it alone.

## 2. S0.1 — Branch, version, and what has changed since `mdf1`

1. **Carriers.** Every local and remote-tracking ref containing `quarto/simulations/actg175/continuous/sim_fs_maxeffCons_mr_field_md_template.qmd`.
2. **The Mac branch.** Which refs named `feature/glm-extension-mac` exist here, their tips, and whether each tip is an ancestor of HEAD and of `feature/glm-extension` (local and remote-tracking), by `git merge-base --is-ancestor`. If none exists, report the newest commit on any ref touching `quarto/simulations/actg175/continuous/` or `R/fs_bias_coverage.R`, and whether HEAD contains it.
3. **`mdf1`'s run commit and platform.** From its record's provenance or its bundle meta; failing both, the parent of the commit that added its bundles. Say which. Include the machine, R version and BLAS it ran on, where recorded.
4. **`R/` since `mdf1`.** `git log --oneline <mdf1 commit>..HEAD -- R/`, one line per commit: whether it can change anything the template computes on the continuous FS path (the simulated data, identification and selection, the unadjusted and oracle quantities, IJ, the harm field, the unstudentized complement field) or cannot (for example DINA/GRF-only, survival-only, or add-only with defaults preserved). Quote hunks only for commits classed "can change". This fixes which `mdf1` columns a re-run identity gate can require exactly.
5. **Loading.** Quote how the template, and the scripts under `quarto/simulations/actg175/continuous/scripts_mdf1/`, load forestsearch (`library()`, `load_all()`, `devtools::install()`).
6. **Installed package against source.** `git archive HEAD DESCRIPTION NAMESPACE R` (plus `src` if the package has one) into `$(mktemp -d)`. In one R session: capture every closure in `asNamespace("forestsearch")` (the installed build), then `pkgload::load_all(<extract>, export_all = FALSE, quiet = TRUE)` and capture again; compare each common name on `formals()` and `body()` after `utils::removeSource()`. Report the counts identical, differing, and present on one side only, with the differing names. Repeat against the tip of `feature/glm-extension` if it differs from HEAD. Remove the extract. If `pkgload` is not installed, say so and skip; do not install it.

## 3. S0.2 — The loaded engine

Confirm or contradict each claim, from the installed namespace (what parallel workers load) and from source at HEAD with `path:line`.

- `fs_mr_inference()` formals: `ci_method = c("field", "ij", "wald")`; `field_complement = TRUE`; `field_scale_complement = c("selected", "none")`; `return_reselection = TRUE`; `field_R_out = 1000L`; `field_R_in = 500L`; `include_complement = FALSE`. Quote also the multiplier-draw and seed formals.
- `forestsearch()` passes, on the FS branch, `include_complement = TRUE`, `ci_method = "field"`, `field_complement = TRUE`, `field_scale_complement = "selected"`. MR on a GLM outcome requires `consistency_method = "resample"`.
- The field code has no outcome-type guard. Each candidate's coefficient and influence contributions come from `.consistency_glm_pieces()`. MD bounds are on the identity scale. Under `adverse_outcome = FALSE` the refit negates Y, so every bound is on the harm-oriented MD scale.
- Returned: `field$lower_1s`; `field$complement$upper_1s_s`, `se_field_s`, `est2_s`, `lower_2s_s`, `upper_2s_s`; `field$complement$upper_1s`, computed beside `upper_1s_s` from the same draws; `field$joint_s$bonf_lower_H`, `bonf_upper_Hc`, `bonf_joint_prob`, `gamma`, `joint_prob`, `lower_H`, `upper_Hc`, `corr`. Say whether the unstudentized `field$joint` is still returned.

## 4. S0.3 — The template's MR call

1. Quote the identification call and the MR call verbatim with `path:line`. State whether MR runs inside `forestsearch()` or through a direct `fs_mr_inference()` call.
2. For each of `ci_method`, `include_complement`, `field_complement`, `field_scale_complement` and `return_reselection`: passed (quote the value) or inherited (from which function's default). State plainly which complement construction the template, as committed, would run today.
3. List every `FS_MD_*` knob the template reads, with its default and what it sets. Quote the campaign-tag guard.

## 5. S0.4 — What it saves per replicate

1. Quote the per-replicate recorder and the batch save. Say whether it saves a fixed set of named columns or a whole result object, and list what it saves.
2. If `upper_1s_s`, `se_field_s`, `est2_s`, `lower_2s_s`, `upper_2s_s` and the `joint_s` fields are not saved, list the template edit that would add them, as named-line additions with insertion points, in transplant form. As the source, quote the committed lines that record the same fields in the survival template under `quarto/simulations/gbsg_020/` (locate it). This is a template change, not an `R/` change, and it is listed, not made.
3. Quote `.refuse_if_tracked()` on the batch save and confirm it is live.
4. From source: can `fs_sim_bias_coverage()` read the field-s complement columns under the recorder's naming with no `R/` change? Its expected signature is `(results, block = c("H", "Hc"), estimators, level, target, side = c("lower", "upper"), scale = c("log", "identity"))`. Quote how `estimators` maps to columns. If it cannot, state exactly what is missing; that is a finding, and no proposal is drafted here.
5. Paths only: the committed document(s) that summarised `mdf1` (its tables and bias–coverage display), and any committed document that writes a long-format metrics file (one row per cell × block × estimator × metric, with Monte Carlo SEs). These are the transplant sources for the re-run's summary and extract.

## 6. S0.5 — Settings, estimators, seeds, and the selection-rule question

1. **Identification**, as the template sets it (quote): `sg_focus`; `effect_neighborhood`, and whether that focus consults it; the effect and consistency thresholds (argument names as the template uses them) and `pconsistency.threshold`; the cut grid and `cut_type`; `maxk`, `n.min` and the arm minima; `fs.splits`; `consistency_method`; `adverse_outcome`; the effect screen; `stop_threshold`.
2. **The four cells** (MD 40 n 500, MD 120 n 500, null n 500, MD 40 n 700): the DGM arguments per cell, and the recorded truths on the harm-oriented scale for the planted region, its complement, and overall.
3. **Estimators.** Quote the lines computing the unadjusted, oracle, IJ two-term, field and complement-field quantities per replicate. For each, say whether it depends on Ĥ.
4. **Seeds.** Data seed, MR seed and any offset, the `sim_id` range and the batch structure, from the template and from `mdf1`'s bundle meta.
5. **`effMaxSG` at ε = 0.20 on this path** (facts for D1; nothing is acted on):
   - (i) whether `sg_focus` and `effect_neighborhood` are set by `FS_MD_*` knobs, or the template edit that would make them so;
   - (ii) the expression computing the effMaxSG band for a continuous outcome, and the scale it acts on (quote);
   - (iii) whether the MR gate's re-selection applies the same focus and band for a GLM outcome; quote both sides and state aligned or not aligned.

## 7. S0.6 — The `mdf1` record

1. **Locations.** The record, the summary document(s), the per-replicate bundles (tracked or untracked, size, rows per cell, column count, present on this machine or not), and the logs.
2. **Claims.** Confirm or contradict each by quotation with line numbers:
   - harm field one-sided lower coverage 0.947–0.950 in every cell;
   - p̂(Ĥ) ≈ 0.16 on a family of about 1,842 candidates with many duplicate-membership labels;
   - unstudentized complement field upper coverage 0.924–0.936;
   - IJ SE/SD 1.5–1.75;
   - Bonferroni joint coverage 0.936–0.943;
   - retained bias in SD units: unadjusted +4, IJ +0.9, field +0.4;
   - field two-sided coverage 0.942–0.978;
   - every display point within 0.008 of the Gaussian reference;
   - the pairing proof enumerated 55 label-tie rows and no selection flips;
   - 13 workers, 172.8 min, peak 17.5 GB, about 15 s per replicate with field and complement at n = 500.
3. **Per cell, both blocks.** Quote these where the record has them; otherwise compute them here with base R from the committed bundles and label them so; if a needed column is absent, say so.
   - declaration rate;
   - mean true β(Ĥ) and β(Ĥᶜ) over declared replicates;
   - the 5, 25, 50, 75 and 95% quantiles of the harm field lower bound and the oracle lower bound on Ĥ;
   - the same quantiles of the unstudentized complement field upper bound and the oracle upper bound on Ĥᶜ;
   - the proportion of field and oracle lower bounds on Ĥ at or above 0, 10, 20, 30 and 40, with Monte Carlo SEs;
   - the complement thresholds the record used, if any, quoted.
4. **The null cell.** Its truth as recorded; its declaration rate; how coverage of β(Ĥ) is conditioned; and every bound-location figure the record gives for it, with the figure's exact definition.

## 8. S0.7 — Estimated cost

1. From `mdf1`'s logs and record: wall, workers, replicates and peak memory per cell; the per-replicate time distribution where logged (median, 90th percentile, maximum), with the nesting of timing columns stated.
2. From source (§2.4, §3): whether today's defaults do per-replicate work `mdf1` did not, such as the field-s scale step or any change to draw counts or re-selection bookkeeping. Quote; no timing run.
3. An estimate, labelled as such, on `mdf1`'s machine at its worker count: 4 cells × 2,000 replicates, and double that for a run under both rules (assuming the second rule costs the same per replicate, which is unmeasured).
4. Quote the render environment and timeout wrapper the `mdf1` scripts used.

## 9. S0.8 — The applied document

For `quarto/applications/actg175/analysis_actg175_continuous_oc.qmd`:

1. Quote the gate call's MR arguments, and the payload element that carries the intervals with its fields. Say whether any field-s field is present.
2. Quote the lines that build the intervals section's complement and joint rows, naming the fields read. Say whether re-rendering under today's defaults would change a number it displays now.
3. Say whether the intervals section can be re-rendered without re-running the document's OC evaluation (quote any flag or cache). Give the last recorded render time and memory per worker, and the cross-platform tolerance its record states for payload identity.

## 10. Record, closeout, post-conditions

1. **Record**, at the header path: §1's output verbatim; S0.1–S0.8 in order, every claim marked confirmed, differs, or not determinable here, with its evidence; then **Facts for Gate 0**, one bullet list per decision, with facts copied from the sections and no recommendation:
   - **D1**, the re-run's selection rule (`mdf1`'s, `effMaxSG` at ε = 0.20, or both): from §6, §2.4 and §8.3.
   - **D2**, whether the null cell is presented: from §7.4.
   - **D3**, bound-location thresholds on the harm-oriented MD scale, including whether the MD 120 cell needs higher ones: from §7.3.
   - **D4**, the machine: from §1, §2.2, §2.3, §2.6, §7.1 and §8.
   - **D5**, the optional applied-document update: from §9.

   Commit the record by explicit path.
2. **`current_status.md`.** If `quarto/simulations/actg175/continuous/` has one and a generator for it, regenerate it with that generator as the last action and commit it alone. If either is missing, leave the directory as it is and say so in the closing message.
3. **Post-conditions**, asserted after the last commit and printed in the closing message:
   - `git diff --name-only <§1 HEAD>..HEAD` lists only the task document, the record and, if regenerated, `current_status.md`;
   - `git status --porcelain --untracked-files=no` shows the same tracked modifications as §1, and no others;
   - `packageDescription("forestsearch")$Built` equals §1's value;
   - the record has non-empty sections S0.1–S0.8 and a Facts for Gate 0 block covering D1–D5;
   - if regenerated, the pin stated in `current_status.md` resolves to the same commit as `HEAD~1`.
4. **Closing message:** the commit range for Larry to push, the post-conditions, and the findings (each claim that differs). Then stop.

## 11. Out of scope

Any edit beyond the files in §10; installs, renders, simulations and timing runs; drafting the template edit beyond listing it; any `R/` proposal; DINA and GRF; pushing.
