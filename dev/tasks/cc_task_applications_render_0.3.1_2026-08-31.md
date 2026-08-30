# CC task — render the four applications under 0.3.1, compare against baseline HTML, stage dated payloads

**Date:** 2026-08-31 · **Repo:** `~/Documents/GitHub/forestsearch` · **Branch:** `feature/glm-extension` · **Package under test:** installed forestsearch, must be **0.3.1** and match `DESCRIPTION` at HEAD
**Type:** renders + one knob-level `.qmd` edit + commits. **Nothing under `R/` is touched. No LOO, no K-fold, no cross-validation of any kind executes. `quarto/GuoHe/` must be byte-identical before and after.**
**Report:** `dev/tasks/REPORT_applications_render_0.3.1_2026-08-31.md`
**Compute (approved by Larry via kickoff):** ≈ 35–55 min total, ceiling ~60 (multimethod ≈ 13 min and compare_all ≈ 6 min per the 08-27 measurements; template minutes; frozen_family without LOO unmeasured, B = 1000 GH bootstrap, expect 10–30 min).

## Scope and standing decisions (Larry, 2026-08-31)

Four documents, cheap-first render order:

1. `quarto/applications/actg175/template_actg175_continuous.qmd`
2. `quarto/applications/actg175/analysis_actg175_continuous_compare_all.qmd`
3. `quarto/applications/gbsg/analysis_gbsg_survival_frozen_family.qmd`  *(after the flag edit below)*
4. `quarto/applications/gbsg/analysis_gbsg_survival_multimethod.qmd`

- **Core analyses only.** frozen_family's `run_loo` is flipped `TRUE → FALSE` as a **standing, committed** configuration change. All CV-type flags across the four must verify FALSE before any render. **MR stays ON** (`run_mr <- TRUE`, `mr_draws <- 5000L` in multimethod) — MR is core inference, not cross-validation.
- **Comparison primitive:** text-layer HTML diff against the pinned baseline per document; every hunk classified; any SUBSTANTIVE hunk ⇒ renders and payloads stay uncommitted.
- **Dated payloads:** on a passing triage, this render's payload outputs are staged under `quarto/applications/gbsg/_payloads_2026-08-31/` and `quarto/applications/actg175/_payloads_2026-08-31/` and **committed** (force-add if gitignored, with the rule quoted). Template's tracked 0.2.0 payload at the default path is restored so both vintages survive.
- **Commits are separate and droppable pre-push** (Larry pushes; CC never pushes).

The read-only inventory behind this task is committed at `dev/tasks/REPORT_applications_inventory_2026-08-31.md` (its §6 quotes the 08-27 method and baseline). Consult it if a baseline question arises; do not consult anything outside the repository.

## Gates — stop semantics

On a STOP: write the report with everything gathered, banner `STOPPED AT <gate>: <reason>` as the first line, commit **only the report** (and the flag edit if already committed), leave all other new/changed files uncommitted in the worktree, print their paths, end. Do not ask questions.

- **G1 — repo state.** Branch is `feature/glm-extension`; `git status --porcelain` is empty. Dirty tree ⇒ STOP. Record `git log -1 --oneline`.
- **G2 — version.** `Rscript -e 'cat(as.character(packageVersion("forestsearch")))'` equals the `Version:` field of `DESCRIPTION` at HEAD **and** equals `0.3.1`. If installed ≠ DESCRIPTION: install from the tree (`R CMD INSTALL .` or `Rscript -e 'devtools::install(upgrade = "never")'`), re-verify, then proceed. If DESCRIPTION ≠ 0.3.1 ⇒ STOP (wrong tree state for this task).
- **G3 — flags.** After step 1's edit, `grep -n -E '^run_(loo|cv|cv_kfold)\s*<-' ` across the four documents shows **every** such flag `FALSE`, and no other `run_*` flag that gates a CV/LOO/K-fold section is TRUE (read the flag block comments to judge; `run_mr` is exempt and stays TRUE). Any CV-type flag TRUE ⇒ STOP.
- **G4 — per-render success.** Each render exits 0, regenerates its `.html` (mtime after the sentinel), and writes its expected fresh output(s) (payload `.rds` newer than the sentinel). Any failure ⇒ STOP.
- **G5 — invariants.** After all renders: `git status --porcelain quarto/GuoHe` is empty AND `find . -name 'cv_out_*.rds' -newer /tmp/fs_render_sentinel_20260831 -not -path './.git/*'` returns nothing. Any hit means LOO ran somewhere ⇒ STOP.
- **G6 — triage.** Any hunk classified SUBSTANTIVE in step 5 ⇒ triage FAIL: skip steps 7–8, write the report with banner `TRIAGE FAILED`, commit the report only, leave renders and payload outputs uncommitted, print `git checkout --` restore hints without executing them, end.
- **G7 — payload identity.** Every staged payload reads back `forestsearch_version == 0.3.1`. Mismatch ⇒ STOP.

## Steps

### 0. Environment
Record: hostname; branch; `git log -1 --oneline`; `git status --porcelain | wc -l`; installed forestsearch version and `DESCRIPTION` Version (G2 applied here, including the install-if-stale path and, if installed, the re-verification output); R version; the quarto binary used (below) and its `--version`.

### 1. Flag audit and the single edit
1. `grep -n -E 'run_(loo|cv|cv_kfold|mr)\s*<-' <each of the four .qmd>` — record all hits verbatim (repo copies are the source of truth; expected: multimethod `run_cv <- FALSE`, `run_loo <- FALSE`, `run_mr <- TRUE`; frozen_family `run_loo <- TRUE`, `run_cv_kfold <- FALSE`; none in compare_all or template).
2. Edit `quarto/applications/gbsg/analysis_gbsg_survival_frozen_family.qmd`: the line reading `run_loo    <- TRUE      # leave-one-out (N-fold) honest out-of-sample subgroup estimate` becomes `run_loo    <- FALSE     # leave-one-out (N-fold) honest out-of-sample subgroup estimate` (change only TRUE→FALSE; preserve alignment and comment). Re-grep to confirm exactly one changed line and no other diff (`git diff --stat` shows only this file, 1 insertion 1 deletion).
3. Commit: `gbsg frozen_family: run_loo -> FALSE (core analyses only; LOO/CV off per decision 2026-08-31)`.
4. Apply G3 across all four.

### 2. Baselines pinned per document, then the sentinel
Create scratch dirs `dev/tasks/_baseline_html_2026-08-31/` and `dev/tasks/_diffs_2026-08-31/` (both stay uncommitted). For each of the four documents' `.html` sibling:
- If tracked (`git ls-files -- <path>` non-empty): record `git log -1 --format='%h %ad %s' --date=short -- <path>`, then snapshot the HEAD copy: `git show HEAD:<path> > dev/tasks/_baseline_html_2026-08-31/<name>.html`.
- If untracked but present in the worktree: copy it into the snapshot dir and record that the baseline is a local, uncommitted render.
- Either way, record the in-file stamp: `grep -o -E 'forestsearch[^<]{0,40}[0-9]+\.[0-9]+\.[0-9]+' <snapshot> | sort -u` (may be empty for the actg175 pair — record "no legible stamp").
- If no `.html` exists at all: record "no baseline; this render becomes the baseline" and exclude the document from triage counts (not a failure).

Then `touch /tmp/fs_render_sentinel_20260831` — every fresh output must be `-newer` than this file.

### 3. Renders — cheap first, each timed
Use the RStudio-bundled binary: `QUARTO=/usr/lib/rstudio/resources/app/bin/quarto/bin/quarto` (if absent, `command -v quarto` and record which binary was used). From the repo root, for each document in the scope order:

```
/usr/bin/time -f "elapsed %E  maxRSS %MkB" $QUARTO render <path-to-qmd> --to html
```

After each render apply G4: exit code; `.html` newer than sentinel; expected outputs newer than sentinel —
- template → `quarto/applications/actg175/_payloads/template_actg175_continuous/template_actg175_continuous_payload.rds` (this **overwrites the tracked 0.2.0 payload in the worktree** — expected; it is restored in step 7);
- compare_all → `.../_payloads/analysis_actg175_continuous_compare_all/{selected_subgroups_continuous.rds,comparison_continuous.rds}`;
- frozen_family → `quarto/applications/gbsg/_payloads/analysis_gbsg_survival_frozen_family/analysis_gbsg_survival_frozen_family_payload.rds`;
- multimethod → `quarto/applications/gbsg/_payloads/analysis_gbsg_survival_multimethod/analysis_gbsg_survival_multimethod_payload.rds`.

Record elapsed per document. Render all four before triage (a later triage failure does not retroactively invalidate the timing/flag evidence).

### 4. Invariants
Apply G5 (GuoHe untouched; no new `cv_out_*` anywhere). Also record `git status --porcelain | head -40` — the changed/untracked set should consist only of: the four `.html`, the payload outputs above, and the scratch dirs.

### 5. Diff and triage
For each document with a baseline, text layer first:

```
sed -E 's/data:image[^"]*//g' dev/tasks/_baseline_html_2026-08-31/<name>.html > /tmp/old_<name>.txt
sed -E 's/data:image[^"]*//g' <worktree .html>                                > /tmp/new_<name>.txt
diff -u /tmp/old_<name>.txt /tmp/new_<name>.txt > dev/tasks/_diffs_2026-08-31/<name>.diff ; echo "exit $?"
```

Classify **every** hunk into exactly one class:
- **VOLATILE** — dates, wall-clock/timing tables, `built_at`, version stamps (0.2.0 → 0.3.1 is expected everywhere versions print), session/provenance echoes, worker counts, paths.
- **CLAUSE-ORDER** — a subgroup definition string (`sg.harm` / `sg_harm` and kin) whose clause **set** is unchanged but whose order differs. Quote old and new verbatim.
- **SUBSTANTIVE** — anything else: any estimate, CI, p-value, count, table value, membership size, or a definition whose clause set changed. Quote verbatim.

Note image-byte differences separately (count of image blocks; not compared — cross-OS baselines make byte drift expected) — images are outside the triage classes.

Report rules: VOLATILE summarized (count per document + up to 3 examples); every CLAUSE-ORDER and SUBSTANTIVE hunk quoted in full. Dedicated subsection **5.x — vi.grf.min NULL path, measured**: multimethod's `fs_dina` and `fs_grf` selected-subgroup lines, old vs new, whatever the classification — these two fits are the only calls in scope that exercise the new default. Then apply G6.

### 6. (Reserved for triage failure — handled by G6's stop semantics.)

### 7. Payload staging (triage passed)
1. `mkdir -p quarto/applications/gbsg/_payloads_2026-08-31 quarto/applications/actg175/_payloads_2026-08-31`
2. Move the untracked fresh outputs (they are this render's products): `git status --porcelain` must show them as `??` before moving —
   - `mv quarto/applications/gbsg/_payloads/analysis_gbsg_survival_multimethod quarto/applications/gbsg/_payloads_2026-08-31/`
   - `mv quarto/applications/gbsg/_payloads/analysis_gbsg_survival_frozen_family quarto/applications/gbsg/_payloads_2026-08-31/`
   - `mv quarto/applications/actg175/_payloads/analysis_actg175_continuous_compare_all quarto/applications/actg175/_payloads_2026-08-31/`
3. Template (tracked path — copy then restore): `cp -r quarto/applications/actg175/_payloads/template_actg175_continuous quarto/applications/actg175/_payloads_2026-08-31/` then `git checkout -- quarto/applications/actg175/_payloads/template_actg175_continuous/` and confirm `git status` is clean for that path (0.2.0 vintage restored).
4. Read back every staged `*_payload.rds`: `Rscript -e 'p <- readRDS("<path>"); cat(as.character(p$forestsearch_version), "|", format(p$built_at), "|", p$est_scale, "\n")'` — table in the report; apply G7. For compare_all's two `.rds` (no version field expected), record `file.info()$mtime` and top-level `names()` instead.
5. `git check-ignore -v quarto/applications/gbsg/_payloads_2026-08-31 quarto/applications/actg175/_payloads_2026-08-31` — record output; then `git add` the two dated dirs (`-f` if ignored, quoting the matching rule in the report).

### 8. Commits (each separate; Larry may drop any before pushing) and close-out
- b) `applications: 0.3.1 payloads under _payloads_2026-08-31 (dated; committed for cross-machine provenance)` — the two dated dirs.
- c) `applications: 0.3.1 reference renders (actg175 template+compare_all; gbsg frozen_family+multimethod)` — the four `.html` (this **adds** any previously-untracked HTML, making it a tracked reference going forward, and replaces tracked 0.2.0 baselines at HEAD; prior vintages remain at their commits).
- d) `applications: 0.3.1 render-and-compare report` — the report file.
Scratch dirs (`_baseline_html_2026-08-31/`, `_diffs_2026-08-31/`) stay uncommitted; state their paths in the report. Print the final `git log --oneline -6` and stop. **No push.**

## Report skeleton — use exactly these sections

```
# REPORT — applications render-and-compare under 0.3.1 — 2026-08-31
[banner if stopped: STOPPED AT <gate> / TRIAGE FAILED]
## 0 Environment and version gate
## 1 Flag audit and the run_loo edit   (verbatim grep hits before/after; edit commit sha)
## 2 Baselines   (per doc: tracked? · baseline commit or "local uncommitted render" · in-file stamp)
## 3 Render timings   (per doc elapsed; total)
## 4 Invariants   (GuoHe clean; no new cv_out_*; worktree change set)
## 5 Diff triage   (per doc: hunk counts by class; all CLAUSE-ORDER and SUBSTANTIVE hunks verbatim; image-block note)
### 5.x vi.grf.min NULL path, measured   (multimethod fs_dina and fs_grf lines, old vs new)
## 6 Triage verdict
## 7 Payload staging   (moves/copies; template restore confirmed; readback table; check-ignore result)
## 8 Commits   (shas and messages; scratch paths left uncommitted)
## 9 Deviations and notes
```

## Hard rules
- No file under `R/` is read-modified in any way; the only source edit is frozen_family's `run_loo` line.
- No LOO, K-fold, or other cross-validation executes; no file under `quarto/GuoHe/` is created, modified, or removed.
- New artifacts limited to: the four rendered `.html`, their payload outputs (staged as specified), the dated `_payloads_2026-08-31/` dirs, the report, and the two uncommitted scratch dirs under `dev/tasks/`.
- Verify from source: every report claim traces to a quoted line, a command output, or a readback captured in this task. No push.
