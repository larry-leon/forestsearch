# TASK — GRF factor-membership exposure on the survival path: read-only check

**File:** `dev/tasks/TASK_grf_factor_exposure_2026-09-16.md` · **Issued:** 2026-09-16 by chat, on Larry's go
**Machine:** `pop-os` · **Repo:** `~/Documents/GitHub/forestsearch` · **Branch:** `feature/glm-extension`, HEAD containing `a4c063bf`
**Transport:** `~/Downloads` → `dev/tasks/`, committed first and alone (§1)
**Record:** `quarto/simulations/gbsg_020/REPORT_grf_factor_exposure_2026-09-16.md`
**Finding under test:** `quarto/simulations/actg175/continuous/REPORT_md_grf_stage1_2026-09-16.md`, §1.6(c): `.grf_evaluate_subgroup()` (`R/grf_subgroup_labels.R`) compares raw data columns with numeric cuts without coercing factors as the forest's covariate matrix does, so candidates on factor covariates get NA membership and are dropped from re-selection and from MR's family.

**What this is.** The same evaluator served campaign `grfmr` (GRF on the GBSG survival design) and any applied document that runs GRF. Whether those results dropped candidates depends on one fact: whether their data carry factor covariates. This task answers that, read-only, before Larry decides on the evaluator fix and before the identifiers brief relies on `grfmr`. It decides nothing.

## ⚠ CATEGORY

- **Read-only.** No `R/`, template, script, document or payload edit; no install; no campaign.
- **Allowed computation:** data regeneration and at most two GRF identification fits with `mr_inference = FALSE` (§2, §4), from a temporary script outside the repo. Hard timeout 20 min.
- **Writes:** this task document, the record, and the `gbsg_020` catalog as §6.2 allows.
- **Unattended.** Gates stop on failure, never to ask; a statement here that does not hold is a finding.

## Conventions

1. Verify from source; quote `path:line` at HEAD.
2. `git add` by explicit path only; untracked files never staged; no `fetch`, `pull` or `push`.
3. Numbers in the record are computed and pasted, never typed.
4. `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1` for every R process.

---

## 1. Provenance and first commit — GATE

```bash
cd ~/Documents/GitHub/forestsearch
hostname; git branch --show-current; git rev-parse --short HEAD; git log --oneline -3
git status --porcelain --untracked-files=no
git merge-base --is-ancestor a4c063bf HEAD && echo "GRF Stage 1 record in HEAD"
ps -eo pid,etime,args | grep -E '[e]xec/R|[R]script|[q]uarto' | head
Rscript -e 'cat(packageDescription("forestsearch")$Built, "\n")'
```

*GATE:* `pop-os`; branch `feature/glm-extension`; `a4c063bf` in HEAD; no tracked modifications; no R, Rscript or quarto process running. Otherwise stop.

Copy this document from `~/Downloads` to `dev/tasks/TASK_grf_factor_exposure_2026-09-16.md` (exact name, else the single match for `*TASK_grf_factor_exposure_2026-09-16*.md`) and commit it alone.

## 2. The survival simulation's covariates

1. Locate the survival GRF campaign: its template (`quarto/simulations/gbsg_020/sim_fs_maxeffCons_fb_mr_field_m1_template.qmd` or the template `grfmr` ran; quote the runner), its runner and campaign scripts, and its bundles (`git ls-files`). State whether the `grfmr` bundles are present on this checkout.
2. In a temporary directory outside the repo, extract the template with `knitr::purl()` and run only its setup, DGM and data chunks for one harm cell `grfmr` ran (name it) at sim_id 1, with the template's seed scheme. Name the chunks run.
3. Report `str()` of the simulated covariates the identifier sees, and list every factor or character column. If a chunk converts columns before the `forestsearch()` call, quote it and report the classes after conversion.

## 3. The applied data

1. From `quarto/applications/actg175/REPORT_actg175_applied_effmaxsg_stage0_2026-09-16.md` (its S0.1 inventory), list every applied document that runs GRF or DINA (`use_grf`, `use_dina`), with path.
2. For each, quote the lines that construct its analysis data and report the class of every covariate passed to the identifier.
3. Say, per document, whether any covariate is a factor.

## 4. Direct count on the survival replicate — the deciding fact

On the §2 replicate, in the temporary script, run GRF identification once with `mr_inference = FALSE` and the arguments `grfmr` passed (from its runner and a bundle `meta`, quoted; `dmin.grf` and focus/band as recorded), with warnings captured (`withCallingHandlers`). Report:
- every warning verbatim, with its count;
- the number of candidates enumerated, the number with NA membership, and the number admitted;
- the selected subgroup's definition, n and effect;
- if the `grfmr` bundle for that cell is present: the committed sim_id 1 selection beside it (definition, `n_sel`, effect), and whether they agree; if absent, say so.

If §2.3 finds no factor or character column, run this fit anyway and confirm zero NA memberships.

Remove the temporary directory.

## 5. DINA's evaluator — source read

Quote the lines where DINA evaluates candidate membership on the data (`R/dina_subgroup.R` or the helper it calls). State whether they coerce factor covariates before comparison. No fit is run.

## 6. Record, catalog, closeout

1. **Record**, at the header path: provenance; §2–§5 with quotations and the fits' output; then **The answer**, one paragraph: whether the survival simulation carries factor covariates; whether the `grfmr` replicate dropped candidates (count of NA-membership candidates); which applied documents are exposed; whether DINA's evaluator has the same defect. Facts only, no recommendation. Commit by explicit path.
2. **Catalog.** If `quarto/simulations/gbsg_020/` has `status_curated.md` and its generator and check script: add one line under open work stating the answer in a sentence, commit it, regenerate `current_status.md` as the last commit, and run the check script. If any of the three is missing, leave the catalog as it is and say so.
3. **Post-conditions**, printed in the closing message: `git diff --name-only <§1 HEAD>..HEAD` lists only the task document, the record and, if updated, the two catalog files; no tracked modifications; the installed `Built` equals §1's; the temporary directory is gone.
4. **Closing message:** the commit range to push; the answer paragraph; the findings. Then stop.

## Out of scope

- Any edit beyond §6, `R/`, installs, campaigns, the fix itself, pushing.
