# TASK — ACTG175 continuous applied analysis: the field-s complement bound and pair

**File:** `dev/tasks/TASK_actg175_continuous_field_s_2026-09-16.md` · **Issued:** 2026-09-16 by chat, on Larry's D5 disposition (Gate 0 of the MD field re-run, 2026-09-15)
**Machine:** `pop-os` · **Repo:** `~/Documents/GitHub/forestsearch` · **Branch:** `feature/glm-extension`, HEAD containing `b77ca649` (the MD re-run's closeout)
**Transport:** `~/Downloads` → `dev/tasks/`, committed first and alone (§1)
**Document:** `quarto/applications/actg175/analysis_actg175_continuous_oc.qmd`, written `<doc>` below
**Line numbers:** from §9 of `quarto/simulations/actg175/continuous/REPORT_md_field_rerun_stage0_2026-09-15.md` (at `b15c8b7a`); re-read before editing
**Record:** `quarto/applications/actg175/REPORT_actg175_continuous_field_s_2026-09-16.md`

**What this is.** The applied analysis still reports the unstudentized complement field bound and joint pair, while the package now defaults to the studentized complement field (field-s) and the MD simulation re-run recorded field-s as the reported complement product. This task:
- makes the gate call pass `field_scale_complement = "selected"` explicitly;
- adds the field-s complement bound and the field-s pair beside the existing rows;
- points the interpretation at them;
- re-renders the document once.

Every number the document displays today must reproduce.

## ⚠ CATEGORY

- **No `R/` change.** One file is edited, `<doc>`: its gate call, its intervals payload element, two tables and the intervals interpretation. Nothing else in the document changes.
- **One render** of `<doc>`. It re-runs the document's OC loop, since the intervals section cannot render alone (Stage 0 §9.3). Estimate: 15–30 min on pop-os, hard timeout 90 min. Authorized by this kickoff.
- **Unattended.**
  - Gates stop on failure, never to ask.
  - On a stop, write the record, commit it with the task document, leave everything else uncommitted as it stands, and stop.
- **Self-contained.** The document stays self-contained: its wrapper calls stay visible, and nothing is read from `dev/` or from a simulation directory.

## Conventions

1. Verify from source; quote `path:line` at a stated commit.
2. `git add` by explicit path only. Untracked files are never staged. No `fetch`, `pull` or `push`.
3. The render runs with `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1` and the document's committed `params`.
4. Numbers in the record are computed and pasted, never typed.
5. Excluded everywhere: IJ winner-only and winner-floor, κ variants, `field_uniform`.
6. Bounds are read by location against the thresholds the document already uses, with no significance language.
7. The submitted parent paper is not a topic.

---

## 1. Provenance and first commit — GATE

```bash
cd ~/Documents/GitHub/forestsearch
hostname; git branch --show-current; git rev-parse --short HEAD
git status --porcelain --untracked-files=no
git merge-base --is-ancestor b77ca649 HEAD && echo "MD re-run closeout in HEAD"
ps -eo pid,etime,args | grep -E '[e]xec/R|[R]script|[q]uarto' | head
Rscript -e 'cat(packageDescription("forestsearch")$Built, "\n")'
```

*GATE:* all of the following hold.
- The host is `pop-os` and the branch is `feature/glm-extension`.
- `b77ca649` is in HEAD.
- There are no tracked modifications.
- No R, Rscript or quarto process is running.

Copy this document from `~/Downloads` to `dev/tasks/TASK_actg175_continuous_field_s_2026-09-16.md`, using the exact name or, failing that, the single match for `*TASK_actg175_continuous_field_s_2026-09-16*.md`. Commit it alone.

## 2. Read — no edit

1. **The anchor's focus.**
   - Quote the anchor call (`:92–132`), including its `sg_focus`, `selection_rule` and MR arguments.
   - Quote the line in `R/fs_mr_inference_methods.R` where `.fs_mr_reselection_from_focus()` maps that focus.
   - State what the payload's `settings$reselection = maxeff` names: the identifier's focus, or MR's label for it.
2. **The identity reference.**
   - List every payload file `<doc>` writes, quoting the save lines, and say whether each is tracked.
   - Copy each existing one to a temporary directory outside the repo. These copies are the identity reference for §4.
3. **The interpretation.** Locate the intervals interpretation. Say whether it uses inline R or literal numbers, and quote each sentence that reads the complement bound or the joint pair.

## 3. Edits — `<doc>` only

- **E1 — gate call (`:128–131`).** Add `field_scale_complement = "selected"` to `mr_inference_args`.
- **E2 — intervals payload element `iv` (`:207–239`).**
  - Add `Hc$field_s`, using the names `Hc$field` uses, filled from the complement's `_s` fields: `est2_s`, `upper_1s_s`, `lower_1s_s`, `lower_2s_s`, `upper_2s_s`, `se_field_s` (as `lambda_sd`), `lambda_mean_s`.
  - Add `joint_s`, with the fields `joint` carries, from `f$joint_s`.
  - Add `field_scale_complement` to `settings`.
  - Existing fields stay exactly as they are.
- **E3 — complement table (`tab2`, `:277–290`, with `up1` at `:203–204`).**
  - Add a row "field-s": point `est2_s`, SE `se_field_s`, two-sided `lower_2s_s`/`upper_2s_s`, one-sided upper `upper_1s_s`.
  - Relabel the existing field row "field (unstudentized)".
- **E4 — joint table (`tab3`, `:302–311`).**
  - From `joint_s`, add the field-s Bonferroni pair (`bonf_lower_H`, `bonf_upper_Hc`, `bonf_joint_prob`) and the field-s calibrated split (`lower_H`, `upper_Hc`, `joint_prob`).
  - Relabel the existing rows "(unstudentized)".
- **E5 — interpretation.**
  - The sentences quoted in §2.3 now read the field-s complement bound and the field-s Bonferroni pair, each by location against the document's thresholds.
  - The unstudentized value is named once, as the paired before-and-after.
  - Use inline R from `iv` wherever the document already does.
  - No other prose changes.

## 4. Render — GATE

Render `<doc>` once under GNU `timeout` 90 min, with the 5-second summed-RSS sampler `quarto/simulations/actg175/continuous/scripts_mdsgnb20/mem_sampler.sh`. Record wall and peak memory.

*GATE — identity against the §2.2 copies,* field by field over everything the copies contain:
- **Tolerance.** Character, integer and logical fields are identical. Numeric fields agree within 1e-8 relative, or 1e-10 absolute near zero. The copies were built on the Mac under R 4.5.2, so bit-identity is not expected.
- **One allowed class: a label tie.** Only label-valued fields differ, while the selected subgroup's size (66) and every numeric field agree within tolerance. Enumerate it.
- **Anything else fails the gate:** a different subgroup, or a numeric field outside tolerance.
- **New fields.** The new `Hc$field_s` and `joint_s` fields are finite, and `joint_s$bonf_lower_H` equals `joint$bonf_lower_H` wherever the two draw counts agree.

## 5. Record and closeout

1. **Record,** at the header path, containing:
   - provenance and the §2 facts;
   - the edits, as the quoted diff;
   - render wall and peak memory;
   - the identity result: fields compared, largest relative difference, any enumerated label tie;
   - a small table of the field-s complement bound and field-s pair beside the unstudentized values, followed by a two-to-four-bullet reading by bound location;
   - findings and commits.
2. **Commit** by explicit path: `<doc>`, its rendered HTML if tracked, each payload file only if tracked, and the record.
3. **Copy** this record and `quarto/simulations/actg175/continuous/REPORT_md_field_rerun_2026-09-15.md` to `~/Downloads`.
4. **Post-conditions,** printed in the closing message:
   - there are no tracked modifications;
   - `git diff --quiet <§1 HEAD>..HEAD -- R/` succeeds;
   - `git diff --name-only <§1 HEAD>..HEAD` lists only the task document, `<doc>`, its HTML and tracked payloads, and the record;
   - the `~/Downloads` copies are `cmp`-identical to the committed files.
5. **Closing message:** the commit range to push, the identity result, the field-s bound and pair beside the unstudentized values, and the findings. Then stop.

## Out of scope

- `R/`.
- Any other document, or any other section of this one.
- The simulation directories.
- Thresholds.
- The excluded constructions.
- Pushing.
