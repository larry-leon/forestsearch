# TASK — DINA and GRF on the ACTG175 continuous (MD) design: Stage 0, read-only

**File:** `dev/tasks/TASK_md_dina_grf_stage0_2026-09-16.md` · **Issued:** 2026-09-16 by chat, on Larry's go to scope the three-identifier set on the continuous path
**Machine:** `pop-os` · **Repo:** `~/Documents/GitHub/forestsearch` · **Branch:** `feature/glm-extension`, HEAD containing `b77ca649`
**Runs after:** `TASK_actg175_applied_effmaxsg_stage0_2026-09-16.md` has closed out
**Transport:** `~/Downloads` → `dev/tasks/`, committed first and alone (§1)
**Directory:** `quarto/simulations/actg175/continuous/`, written `<dir>` below
**Record:** `<dir>/REPORT_md_dina_grf_stage0_2026-09-16.md`

**What this is.** The MD re-run (`mdsgnb20`) ran the FS identifier only. The survival grid runs FS, DINA and GRF, all at `effMaxSG` with ε = 0.20. Matching it on this design would mean DINA and GRF campaigns on the same four cells.

This stage establishes, from source and from one identification fit per identifier on a single replicate:
- whether DINA and GRF can run on this design as the package stands;
- on which floors and orientation each would run;
- which MR arguments a campaign must pass;
- what the template needs;
- cost references.

It decides nothing and drafts no `R/` proposal.

**Governing constraint** (from the survival DINA/GRF workstream):
- DINA and GRF fail the fixed-family condition, so any interval they produce is conditional on the proposed family.
- Comparisons with FS are descriptive, with the confound named.
- The record uses this language wherever it touches their intervals.

## ⚠ CATEGORY

- **Read-only.** No `R/`, template, script, document or payload edit. No install, render or simulation.
- **One allowed computation:** the two identification fits of §4, run from a temporary script outside the repo.
- **Writes:** this task document, the record, and the directory catalog (§8).
- **Unattended.** A gate stops on failure, never to ask. A claim below that does not hold is a finding, not a stop. Estimate: under an hour.

## Conventions

1. Verify from source; quote `path:line` at HEAD.
2. `git add` by explicit path only. Untracked files are never staged. No `fetch`, `pull` or `push`.
3. Numbers in the record are computed and pasted, never typed.
4. **Excluded everywhere:** IJ winner-only and winner-floor, κ variants, `field_uniform`, covariate adjustment, tuned inflation factors.
5. No significance language. The submitted parent paper is not a topic.

---

## 1. Provenance and first commit — GATE

```bash
cd ~/Documents/GitHub/forestsearch
hostname; git branch --show-current; git rev-parse --short HEAD; git log --oneline -3
git status --porcelain --untracked-files=no
git merge-base --is-ancestor b77ca649 HEAD && echo "mdsgnb20 closeout in HEAD"
ps -eo pid,etime,args | grep -E '[e]xec/R|[R]script|[q]uarto' | head
Rscript -e 'cat(packageDescription("forestsearch")$Built, "\n")'
```

*GATE:* all of the following hold, or stop.
- The host is `pop-os` and the branch is `feature/glm-extension`.
- `b77ca649` is in HEAD.
- There are no tracked modifications.
- No R, Rscript or quarto process is running.

Copy this document from `~/Downloads` to `dev/tasks/` (exact name, else the single match for `*TASK_md_dina_grf_stage0_2026-09-16*.md`) and commit it alone.

## 2. S0.1 — The identifiers on a continuous outcome

1. **Selecting DINA and GRF.** Quote the argument that makes `forestsearch()` run each of them.
2. **Continuous outcomes.** Say whether each accepts `outcome_type = "continuous"` with `effect_measure = "MD"` and `adverse_outcome = FALSE`, and quote any guard, error or warning.
3. **Transplant source.** Quote the survival template's identifier knob and the arguments it passes per identifier: `quarto/simulations/gbsg_020/sim_fs_maxeffCons_fb_mr_field_m1_template.qmd` at HEAD.

## 3. S0.2 — Orientation and floors

**Claims to confirm or contradict.** These are the MD re-run handoff's code readings, not runs.
- DINA's proposal floor is applied to the raw surface with no adverse-outcome orientation, so under `adverse_outcome = FALSE` its proposal floor and its re-ranking floor point in opposite directions.
- DINA's continuous floor is `hr.threshold` itself (default 1.25), while FS screens at MD > 0.

**Side-by-side table** for FS, DINA and GRF on this path, one row per identifier:
- the surface it ranks;
- that surface's orientation under `adverse_outcome = FALSE`;
- every floor applied, with name, default, scale and line;
- the effect screen.

**For GRF,** also cover `dmin.grf` (its scale on a continuous outcome and its default), the frontier filter, and `vi.grf.min`.

**Alignment.** For DINA and for GRF, state which argument values put their floors on FS's harm-oriented MD scale. If no argument can, say so; that makes it an `R/` behaviour change. Record it as a finding and draft no proposal.

## 4. S0.3 — One replicate, observed

In a temporary script outside the repo, rebuild md40 n500 at sim_id 1 with the MD template's data generation. Use its seed scheme (8316951 + sim_id) and its identification settings.

Run `forestsearch()` twice with `mr_inference = FALSE`, once for DINA and once for GRF:
- `sg_focus = "effMaxSG"` and `effect_neighborhood = 0.20`;
- `selection_rule` as recorded in the `mdsgnb20` bundle `meta`;
- every DINA- or GRF-specific argument at its package default, with the defaults recorded.

For each fit, record:
- whether it ran, with any error or warning verbatim (an error is a finding);
- the size of the proposed family;
- the floors as applied, from the returned object where it carries them;
- the selected subgroup's definition, n and oriented MD;
- beside these, FS's result for sim_id 1 from the committed `mdsgnb20` md40 n500 bundle.

Remove the temporary directory at the end.

## 5. S0.4 — Rule and MR for DINA and GRF on this outcome

1. **Focus and band.**
   - Quote how `effMaxSG` resolves for DINA and GRF (`fs_focus_tag()`).
   - Quote where each uses the band: DINA as a sort key, GRF as a frontier filter.
   - Name the scale the band acts on for MD.
2. **MR eligibility.** Quote the conditions `forestsearch()` checks before running MR for DINA and GRF (`.mr_dina_ok`, `.mr_grf_ok`), and say whether they admit a GLM outcome.
3. **Arguments.**
   - For each MR argument the MD template passes, say whether `.fs_apply_mr()` forwards it or hard-codes a value. The arguments are `ci_method`, `draws`, `include_complement`, `field_complement`, `field_scale_complement`, `return_reselection`, `ij_residual`, `confirm_rule` and `t_confirm`.
   - *Claim to confirm or contradict:* `.fs_apply_mr()` hard-codes `ci_method = "ij"`, so DINA and GRF get the field only when `ci_method = "field"` is passed.
4. **Re-selection.** For DINA and GRF on this outcome, say whether MR's re-selection uses the same focus and band as the identifier: aligned or not, quoting both sides.

## 6. S0.5 — What a campaign would need

1. **Template.**
   - List the edits to `<dir>/sim_fs_maxeffCons_mr_field_md_template.qmd` that would run DINA and GRF: an identifier knob and the per-identifier arguments.
   - Write them in transplant form from the survival template, quoting source lines and insertion points.
   - List any recorder fields that the survival DINA and GRF campaigns record and this template lacks.
   - List only; make no edits.
2. **Same draws.**
   - Quote the template's order of seeding, data simulation and identification.
   - State whether the choice of identifier can change the simulated data. Pairing with `mdsgnb20` rests on `n_true` and the oracle columns.
3. **Drivers and gates.** Give the paths of the survival DINA and GRF campaign scripts (`scripts_dinamr/`, including `gate2.R` and `gate2G.R`; locate with `git ls-files`) as transplant sources.
4. **Cost references,** quoted from records; no timing run.
   - Per-replicate `fit_mr_secs` for DINA (`dinamr`) and GRF (`grfmr`) against FS on the same survival cells.
   - FS on this design, from the `mdsgnb20` Gate 2 record.

## 7. Facts for Larry's decision

A block in the record, facts only and no recommendation:
- whether DINA can run on this design without an `R/` change, and whether GRF can;
- the floors and orientation each would run on, beside FS's;
- the MR arguments a campaign must pass explicitly;
- the template edits;
- the cost references.

## 8. Record, catalog, closeout

1. **Record** at the header path, with provenance, §2–§7 in order, findings and commits. Commit it by explicit path.
2. **Catalog.**
   - Add one open-work line to `<dir>/status_curated.md`: "DINA and GRF on this design: Stage 0 recorded, decisions pending". Commit it.
   - Regenerate `<dir>/current_status.md` with the directory's generator, as the last commit.
   - `check_current_status.sh --commit` passes.
3. **Post-conditions,** printed in the closing message:
   - `git diff --name-only <§1 HEAD>..HEAD` lists only the task document, the record, `status_curated.md` and `current_status.md`;
   - there are no tracked modifications;
   - the installed `Built` equals §1's value;
   - the temporary directory is gone.
4. **Closing message:**
   - the commit range to push;
   - one line each for DINA and GRF on whether it can run as the package stands, naming the blocking lines if it cannot;
   - the one-replicate results beside FS;
   - the findings.

   Then stop.

## Out of scope

- Any edit beyond §8.
- `R/` proposals.
- Campaigns, renders and installs.
- The applied documents.
- Pushing.
