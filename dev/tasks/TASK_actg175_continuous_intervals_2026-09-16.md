# TASK — ACTG175 continuous applied analysis under effMaxSG, ε = 0.20: the intervals document

**File:** `dev/tasks/TASK_actg175_continuous_intervals_2026-09-16.md` · **Issued:** 2026-09-16 by chat, on Larry's decision that the applied analysis runs `effMaxSG` at ε = 0.20, the rule of the survival grid and of campaign `mdsgnb20`
**Machine:** `pop-os` · **Repo:** `~/Documents/GitHub/forestsearch` · **Branch:** `feature/glm-extension`, HEAD containing `700bdcdd`
**Transport:** `~/Downloads` → `dev/tasks/`, committed first and alone (§1)
**Source document:** `quarto/applications/actg175/analysis_actg175_continuous_oc.qmd` at HEAD, written `<oc>`; read only, never edited
**New document:** `quarto/applications/actg175/analysis_actg175_continuous_intervals.qmd`, written `<doc>`
**Facts:** `quarto/applications/actg175/REPORT_actg175_applied_effmaxsg_stage0_2026-09-16.md` (the anchor under both rules) and `quarto/applications/actg175/REPORT_actg175_continuous_field_s_2026-09-16.md` (the maxeffCons intervals)
**Record:** `quarto/applications/actg175/REPORT_actg175_continuous_intervals_2026-09-16.md`

**What this is.** The applied analysis reports its intervals from `<oc>`, whose anchor is selected under `maxeffCons` and hard-wired into its operating-characteristics (OC) sections, and whose OC functions model the `maxeffCons` pick only. Larry's decision is that the analysis runs `effMaxSG` at ε = 0.20. This task creates a self-contained intervals document under that rule — the data, one `forestsearch()` call with MR, the interval table, the interpretation — and leaves `<oc>` untouched. The OC evaluation is not switched; that needs an `R/` method change and is deferred.

## Dispositions (Larry, 2026-09-16)

- A new document; `<oc>` is not edited.
- One identification call; FS only.
- Rule: `sg_focus = "effMaxSG"`, `effect_neighborhood = 0.20`, `selection_rule` as a `mdsgnb20` bundle's `meta` records (§2.2).
- MR arguments exactly as `<oc>`'s gate call passes them.
- Products: the field one-sided lower bound on Ĥ; the field-s one-sided upper bound on Ĥᶜ; their Bonferroni pair; unadjusted, IJ two-term and the unstudentized complement field beside them; p̂(Ĥ).
- Reading: by bound location against `<oc>`'s thresholds, on the harm-oriented scale and on raw `cd4_change`, with no significance language.
- The `maxeffCons` result is cited once, for the contrast, from the field-s record, stated as quoted.

## ⚠ CATEGORY

- **No `R/` change.** One new `.qmd`, its render, its payload and the record. `<oc>` is unchanged.
- **One render** of `<doc>`. There is no OC loop, so minutes; hard timeout 30 min. Authorized by this kickoff.
- **Unattended.** Gates stop on failure, never to ask. On a stop, commit the task document and the record, leave the rest uncommitted, and stop.
- **Self-contained.** Nothing is read from `dev/` or from a simulation directory; the wrapper calls are visible and few.

## Conventions

1. Verify from source; quote `path:line` at HEAD.
2. `git add` by explicit path only. Untracked files are never staged. No `fetch`, `pull` or `push`.
3. The render runs with `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1`.
4. Numbers in the record are computed and pasted, never typed.
5. Excluded everywhere: IJ winner-only and winner-floor, κ variants, `field_uniform`.
6. No significance language. The submitted parent paper is not a topic.

---

## 1. Provenance and first commit — GATE

```bash
cd ~/Documents/GitHub/forestsearch
hostname; git branch --show-current; git rev-parse --short HEAD; git log --oneline -3
git status --porcelain --untracked-files=no
git merge-base --is-ancestor 700bdcdd HEAD && echo "DINA/GRF Stage 0 in HEAD"
git diff --quiet 0071c17e..HEAD -- R/ DESCRIPTION NAMESPACE && echo "no R/ change since the install"
test -e quarto/applications/actg175/analysis_actg175_continuous_intervals.qmd && echo "EXISTS" || echo "new path free"
ps -eo pid,etime,args | grep -E '[e]xec/R|[R]script|[q]uarto' | head
Rscript -e 'cat(packageDescription("forestsearch")$Built, "\n")'
```

*GATE:* all of the following hold, or stop.
- The host is `pop-os` and the branch is `feature/glm-extension`.
- `700bdcdd` is in HEAD, and there is no `R/`, `DESCRIPTION` or `NAMESPACE` change since `0071c17e`.
- The installed `Built` is `2026-09-16 05:57:14 UTC` (no install is authorized).
- The new path is free.
- There are no tracked modifications, and no R, Rscript or quarto process is running.

Copy this document from `~/Downloads` to `dev/tasks/TASK_actg175_continuous_intervals_2026-09-16.md` (exact name, else the single match for `*TASK_actg175_continuous_intervals_2026-09-16*.md`) and commit it alone.

## 2. Read — no edit

1. **Transplant sources.** Quote from `<oc>` at HEAD with line numbers: the YAML header and `params`; the setup and data chunks; the anchor `forestsearch()` call with every argument; the intervals chunks, from the guards and `lo1`/`up1` through `tab1`–`tab3` and the `iv` payload element; the payload save lines, with the payload's path and whether it is tracked; the intervals interpretation.
2. **Rule.** Quote `selection_rule` from a `mdsgnb20` combined bundle's `meta` under `quarto/simulations/actg175/continuous/`.
3. **Reproduction target.** Quote the effMaxSG anchor from the effMaxSG Stage 0 record: definition, n = 79, MD = 71.10, consistency 0.91, band floor, and the candidates in band.
4. **Contrast citation.** Quote from the field-s record the maxeffCons results: Ĥ's definition and n; the unadjusted, IJ and field one-sided lower bounds on Ĥ; the field-s and unstudentized upper bounds on Ĥᶜ; the Bonferroni pair and its joint probability.

## 3. Build `<doc>` — transplant

Create `<doc>` from the §2.1 chunks, in this order. Every chunk is a copy with only the named changes.

1. **Header.** Title: "ACTG175 continuous (CD4 change): post-selection intervals under effMaxSG, ε = 0.20". `params` keeps `seed` and any parameter the anchor or intervals chunks read; OC-only parameters (`draws`, `n_workers`) are dropped. The header note that everything is computed in the document and nothing is read from disk stays, reworded for this document.
2. **Setup and data chunks** unchanged.
3. **The anchor call**, changed only in: `sg_focus = "effMaxSG"`, `effect_neighborhood = 0.20`, `selection_rule` per §2.2. The MR arguments stay exactly as quoted (`ci_method = "field"`, `draws = 5000L`, `include_complement = TRUE`, `field_complement = TRUE`, `field_scale_complement = "selected"`, `return_reselection = TRUE`). `<oc>`'s literal anchor assertion is replaced by an assertion of the §2.3 anchor (its definition and n = 79), so the document stops if the selection changes.
4. **The intervals chunks** unchanged in code: `lo1`, `up1` with `field_s`, `tab1`–`tab3`, the `iv` element. Captions state the rule. The payload file name gains the suffix `_effmaxsg`, in the same directory and with the same tracking status as `<oc>`'s. If `<oc>` does not print p̂(Ĥ), add one line printing it from the re-selection object.
5. **Interpretation section**, inside the document, written for this anchor:
   - the subgroup, its size, its MD on both scales, its consistency proportion, and the in-band candidates it was chosen among;
   - the bounds on Ĥ by location against the ladder, on both scales;
   - Ĥᶜ: what benefit the field-s upper bound supports, with the unstudentized bound named once beside it;
   - the pair;
   - p̂(Ĥ) and what regime it indicates;
   - one paragraph contrasting with the maxeffCons analysis, citing §2.4 at commit `570bf8a5` as quoted values, not recomputed;
   - the scope sentence, verbatim: "These are selection-adjusted bounds on one trial under one selection rule; the operating-characteristics evaluation in analysis_actg175_continuous_oc.qmd remains under maxeffCons."
6. Nothing from `<oc>`'s OC sections is transplanted.

Commit `<doc>`.

## 4. Render — GATE

Render `<doc>` under `timeout 30m` with the thread variables. Record wall and peak memory.

*GATES:*
- The render succeeds.
- The anchor assertion passed: Ĥ is §2.3's definition with n = 79; MD within 1e-6 of 71.10; consistency within 1e-8 of the Stage 0 value. Quote the printed values.
- Every numeric field of `iv` is finite; `joint_s$bonf_lower_H` equals `joint$bonf_lower_H`; the payload file exists at the stated path.

Report, as facts: `Hc$field_s$upper_1s` beside `Hc$field$upper_1s`; p̂(Ĥ).

## 5. Record and closeout

1. **Record**, at the header path: provenance; the §2 quotations; the transplant, as the diff against `<oc>`'s chunks; render wall and peak memory; the anchor check; the intervals table pasted from the render on both scales; a reading of three to five bullets by bound location; the maxeffCons contrast, quoted; findings; commits.
2. **Commit** by explicit path: `<doc>`; its HTML if `<oc>`'s is tracked; its payload if `<oc>`'s is tracked; the record.
3. **Copy** the record and the HTML to `~/Downloads`.
4. **Post-conditions**, printed in the closing message:
   - no tracked modifications;
   - `git diff --quiet <§1 HEAD>..HEAD -- R/` succeeds;
   - `git diff --name-only <§1 HEAD>..HEAD` lists only the task document, `<doc>`, its HTML and payload where tracked, and the record;
   - the `~/Downloads` copies are `cmp`-identical to the committed files.
5. **Closing message:** the commit range to push; the anchor; the field lower bound on Ĥ, the field-s upper bound on Ĥᶜ and the pair, on both scales; p̂(Ĥ); the findings. Then stop.

## Out of scope

- `<oc>` and every other document.
- `R/`.
- The OC evaluation.
- DINA and GRF.
- Pushing.
