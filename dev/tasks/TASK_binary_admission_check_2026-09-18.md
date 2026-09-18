# TASK — The minimum-events admission criterion on the binary path: read-only check

**File:** `dev/tasks/TASK_binary_admission_check_2026-09-18.md` · **Issued:** 2026-09-18 by chat, after Larry identified the intended criterion as at least 10 events in **each** arm, mirrored from the survival path
**Machine:** `pop-os` · **Repo:** `~/Documents/GitHub/forestsearch` · **Branch:** `feature/glm-extension`
**Transport:** `~/Downloads` → `dev/tasks/`, committed first and alone (§1)
**Directory:** `quarto/simulations/actg175/binary_020/`, written `<dir>`
**Record:** `<dir>/REPORT_binary_admission_check_2026-09-18.md`

**Why.** Stage 2 of the OR campaign halted on a diverged logistic fit: the true region on one replicate held 43 subjects with the treated arm at 17 events and 0 non-events. Two facts about the design stand behind it:
- the planted region averages about 48 subjects at n = 500 (prevalence 9.632%), below the search's own `n.min = 60`;
- the study's `.logit_or_ci()` guard counts at least 5 events and 5 non-events **pooled over arms**, not per arm.

Larry's criterion is at least 10 events in **each** arm, as the survival path applies it, governing whether a subgroup may be evaluated as a candidate. This task establishes, read-only, whether the binary path implements it, whether declared subgroups satisfy it, and how often the true region would fail it. **It decides nothing and relaunches nothing.**

## ⚠ CATEGORY

- **No `R/` change. No campaign. No relaunch. No install.**
- **Writes:** this task document, the record, and §2's settlement of the working tree.
- **Allowed computation:** §4 regenerates simulated data and recomputes counts. No `forestsearch()` call, no MR, no render other than §2.1's.
- **Unattended.** Gates stop on failure, never to ask. Ceiling 1.5 h.

## Conventions

1. Verify from source; quote `path:line` at a stated commit.
2. `git add` by explicit path only; untracked files never staged; no `fetch`, `pull` or `push`.
3. Numbers are computed and pasted, never typed.
4. The study's data recipe is frozen; nothing here proposes changing it.

---

## 1. Provenance — GATE

Assert: host `pop-os`; branch `feature/glm-extension`; HEAD contains the halt commit `98d1e3bb`; no R, Rscript or quarto process running. Record `git status --porcelain` in full, including untracked, and `git log --oneline -6`. A dirty tree is expected here and is not a stop.

Copy this document from `~/Downloads` to `dev/tasks/TASK_binary_admission_check_2026-09-18.md` (exact name, else the single match for `*TASK_binary_admission_check_2026-09-18*.md`) and commit it alone.

## 2. Settle the working tree

### 2.1 The template and checker patches — assertion-only, then commit

The patches to `<dir>/sim_fs_mr_field_or_template.qmd`, `<dir>/scripts_or/gate2.R` and `<dir>/scripts_or/smoke_identity.R` split the OR invariant: estimates must be strictly positive, a bound fails only when negative, and zero or infinite bounds are counted and reported. That split stands.

- **Check:** re-render `sim_id` 1–20 of a committed green cell with the patched template, tag `admsmoke`, and compare every column except `*_secs` against the committed rows. They must agree within 1e-12 relative, with identical NA patterns.
- *GATE:* if anything differs, the patch is not assertion-only. Stop and report; commit nothing from this section.
- **On pass:** commit the three files together, with a message stating the split and that it changes no recorded value. Leave the smoke output untracked and list it.

### 2.2 The summary's oracle masking — revert

`<dir>/summary_actg175_or.qmd` was patched to mask non-convergent oracle fits out of that estimator's statistics. That convention is superseded by Larry's decision pending this check, and masking one estimator only is not the convention that will be adopted. Restore the file to HEAD (`git checkout -- <path>`) and record that it was reverted and why.

## 3. The criterion, from source

For the **survival** path and the **binary / GLM** path separately, quote the lines and answer:

1. **Is there a minimum-events rule?** Where it is applied, its value, and whether it counts events **per arm** or pooled.
2. **What does it govern?** Candidate admission — whether a subgroup may enter the evaluated family or be declared — or only whether an interval is computed.
3. **The arm minima.** What `d0.min` and `d1.min` count on each path: subjects or events. Likewise `n.min`.
4. **The study's helper.** Quote `.logit_or_ci()`'s guard and state which of the above it is.
5. **The mirror.** State plainly whether the binary path implements the survival path's rule, and if not, what it implements instead.

Name every function and line; do not infer from argument names.

## 4. What the current design produces

No `forestsearch()` calls. Regenerate data with the committed recipe and count.

### 4.1 The true region against the floors

For each design point (`target_or_h` = 0.75, 1.0, 1.5) and each of n = 500, 750, 1000, 2000, regenerate 200 replicates (`sim_id` 1–200, the study's seed scheme) and report for the true region H:
- the mean, 5th and 95th percentiles of |H|;
- the mean and minimum of events and non-events **in each arm**;
- the share of replicates failing each of: `n.min = 60`; 10 subjects in each arm; **10 events in each arm**; 5 events and 5 non-events pooled;
- the share with complete separation in an arm, meaning zero events or zero non-events.

Report the same for Hᶜ where a floor could bind.

### 4.2 Declared subgroups against the criterion

For the two committed green cells and the committed rows of the halted cell, take 200 declared replicates each. Regenerate the data, evaluate the recorded subgroup definition, and report for Ĥ and Ĥᶜ:
- events and non-events in each arm: mean and minimum;
- the share failing 10 events in each arm;
- for any replicate that fails, its `sim_id`, definition, size and arm counts.

This says whether the search declared subgroups the intended criterion would have barred.

### 4.3 Reach

State whether the same question arises for the committed continuous campaigns (`mdsgnb20`, `mdgrf`, `mddina`), which have no events, and for the committed survival campaigns, from §3's answer alone. One sentence each, no computation.

## 5. Record and closeout

1. **Record** at the header path: provenance; §2's outcome; §3 with quotations; §4's tables; then **Facts for Larry's decisions**, facts only:
   - whether the binary path implements at least 10 events in each arm, and what it implements instead;
   - whether any declared subgroup in the committed cells violates that criterion;
   - how often the true region fails each floor, by design point and n;
   - whether n = 500 can ever declare a subgroup as small as the planted region.

   Commit the record by explicit path.
2. **Bundle:** write `~/Downloads/bundle_binary_admission_2026-09-18.zip` from HEAD with the record.
3. **Post-conditions**, printed in the closing message: `git diff --name-only <§1 HEAD>..HEAD` lists only this document, the record, and §2.1's three files; `R/` unchanged; `<dir>/summary_actg175_or.qmd` matches HEAD as of §1; no campaign was launched; the installed build is unchanged.
4. **Closing message:** the commit range to push; the answer to §3.5 in one line; whether any declared subgroup fails the criterion; the §4.1 failure shares at n = 500; the findings. Then stop.

## Out of scope

- `R/` changes and proposals; relaunching Stage 2; re-running any committed cell; the applied documents; pushing.
