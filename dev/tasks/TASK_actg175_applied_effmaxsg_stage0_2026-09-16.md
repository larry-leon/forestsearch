# TASK — Applied analyses under effMaxSG, ε = 0.20: Stage 0, read-only

**File:** `dev/tasks/TASK_actg175_applied_effmaxsg_stage0_2026-09-16.md` · **Issued:** 2026-09-16 by chat, on Larry's direction that the analysis run `effMaxSG` with ε = 0.20
**Machine:** `pop-os` · **Repo:** `~/Documents/GitHub/forestsearch` · **Branch:** `feature/glm-extension`
**Runs after:** `TASK_actg175_continuous_field_s_2026-09-16.md` has closed out
**Transport:** `~/Downloads` → `dev/tasks/`, committed first and alone (§1)
**Record:** `quarto/applications/actg175/REPORT_actg175_applied_effmaxsg_stage0_2026-09-16.md`

**What this is.** Larry wants the applied analysis under `effMaxSG` at ε = 0.20. That is the rule of the survival grid and of the MD simulation re-run (`mdsgnb20`).

The ACTG175 continuous applied document currently runs `maxeffCons`, as the field-s task confirmed in its §2.1. Switching the rule changes the anchor subgroup and everything built on it. Before any edit or render, this stage establishes:
- which applied documents exist and which rule each runs;
- what the anchor subgroup becomes under the new rule;
- what depends on the anchor or on the rule, and whether the OC functions carry the rule;
- what each re-render would cost.

It decides nothing.

## ⚠ CATEGORY

- **Read-only.** No `R/` change, no edit to any document, script or payload, no install, no render, no simulation.
- **One allowed computation:** the two identification fits of §3 on the applied data, with `mr_inference = FALSE`. They run from a temporary script outside the repo.
- **Writes:** this task document and the record.
- **Unattended.** A gate stops on failure, never to ask. A statement here that does not hold is a finding, not a stop. Estimate: under an hour.

## Conventions

1. Verify from source; quote `path:line` at HEAD.
2. `git add` by explicit path only; untracked files are never staged; no `fetch`, `pull` or `push`.
3. Numbers in the record are computed and pasted, never typed.
4. No significance language. The submitted parent paper is not a topic.

---

## 1. Provenance and first commit — GATE

```bash
cd ~/Documents/GitHub/forestsearch
hostname; git branch --show-current; git rev-parse --short HEAD; git log --oneline -3
git status --porcelain --untracked-files=no
git cat-file -e HEAD:quarto/applications/actg175/REPORT_actg175_continuous_field_s_2026-09-16.md && echo "field-s task closed"
ps -eo pid,etime,args | grep -E '[e]xec/R|[R]script|[q]uarto' | head
Rscript -e 'cat(packageDescription("forestsearch")$Built, "\n")'
```

*GATE:* all of the following hold, or stop.
- The host is `pop-os` and the branch is `feature/glm-extension`.
- The field-s task's record exists at HEAD.
- There are no tracked modifications.
- No R, Rscript or quarto process is running.

Copy this document from `~/Downloads` to `dev/tasks/` (exact name, else the single match for `*TASK_actg175_applied_effmaxsg_stage0_2026-09-16*.md`) and commit it alone.

## 2. S0.1 — Inventory of applied documents

For every `.qmd` under `quarto/applications/`, list:
- its path and last commit;
- each `forestsearch()` call, with its `sg_focus`, `effect_neighborhood`, `selection_rule`, identifiers used (`use_dina`, `use_grf`), and whether it runs MR;
- its last recorded render: wall, machine and workers, quoted from its record where one exists.

Present this as one table.

## 3. S0.2 — The anchor under both rules

Extract the headline document `quarto/applications/actg175/analysis_actg175_continuous_oc.qmd` with `knitr::purl()` into a temporary directory outside the repo. Run only its setup and data chunks and its anchor `forestsearch()` call; name the chunks run. Change nothing but the arguments below, and set `mr_inference = FALSE` in both fits.

- **(a) As committed (`maxeffCons`).**
  - *GATE:* it reproduces the committed anchor subgroup, `{age ≤ 37 & cd40 > 507}` with n = 66.
  - A label with identical membership counts as a reproduction.
  - Anything else means the extraction is not faithful: stop.
- **(b) `effMaxSG`, `effect_neighborhood = 0.20`.** Use the `selection_rule` recorded in the `mdsgnb20` bundle `meta`; quote it.

For each fit, report:
- the selected subgroup's definition, n, oriented MD estimate and consistency proportion;
- the complement's n;
- the family size.

For (b), also report:
- the band floor and the number of candidates in the band;
- the five largest in-band candidates (definition, n, MD);
- |Ĥ(a) ∩ Ĥ(b)|, and whether Ĥ(b) contains Ĥ(a).

Remove the temporary directory at the end.

## 4. S0.3 — What depends on the anchor or the rule

For every ACTG175 continuous applied document, quote the lines where the anchor feeds a later section. Classify each by what it depends on:
- the selected subgroup itself;
- a quantity derived from it (T̂, thresholds, the anchored truth Q);
- the selection rule.

Name each section affected.

For the OC functions the documents call (`fs_oc_*` and `fs_family_report`; list them from the namespace), quote the formals and the selection lines. State plainly, from source:
- whether they take `sg_focus`, `effect_neighborhood` or `selection_rule`, directly or through `forestsearch_args`;
- whether their model of the selected subgroup implements `effMaxSG` (largest size in the band) or only the `maxeffCons` pick.

If only the latter, that is a finding. Draft no proposal.

## 5. S0.4 — Cost

For each affected document, quote from its record the last render's wall time, CPU time where recorded, machine, workers and memory per worker.

## 6. Record, closeout

1. **Record,** at the header path:
   - provenance;
   - S0.1–S0.4 in order;
   - a **Facts for Larry's decision** block, facts only and no recommendation, covering:
     - the anchor under each rule;
     - the documents and sections a switch changes;
     - whether the OC sections can run under `effMaxSG` without an `R/` change;
     - the cost of each re-render.

   Commit the record by explicit path.
2. **Post-conditions,** printed in the closing message:
   - `git diff --name-only <§1 HEAD>..HEAD` lists only the task document and the record;
   - there are no tracked modifications;
   - the installed `Built` equals §1's value;
   - the temporary directory is gone.
3. **Closing message:**
   - the commit range to push;
   - the anchor under both rules, side by side;
   - the affected documents and sections;
   - whether the OC functions carry the rule;
   - the costs and findings.

   Then stop.

## Out of scope

- Any edit, render or install.
- Choosing which documents switch; that is Larry's decision.
- `R/` proposals.
- Pushing.
