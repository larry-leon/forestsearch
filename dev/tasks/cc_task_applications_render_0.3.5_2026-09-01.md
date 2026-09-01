# CC task — applications re-render at 0.3.5: `_payloads_2026-09-01` versus the committed 0.3.1 reference

**File:** `dev/tasks/cc_task_applications_render_0.3.5_2026-09-01.md` · **Issued:** 2026-09-01 by chat, commissioned by Larry
**Repo:** `~/Documents/GitHub/forestsearch` (pop-os) · **Branch:** `feature/glm-extension`
**Transport:** `~/Downloads` → `dev/tasks/`, committed first (§1). Hyphen-stripped stem `cc_task_applications_render_0.3.5_20260901.md` is the same file.
**Reads (in-repo):** the 08-31 workstream pair — `dev/tasks/cc_task_applications_render_0.3.1_2026-08-31.md` and its committed report (locate it; state the path) — which are **the procedure this task transplants**: document list, render commands and toolchain, restore list, triage definitions, LOO/CV flag decisions. Also `dev/glm-continuous-sims/HANDOFF`-era reports as context only: the bit-identity chain 0.3.2 → 0.3.5 covered the gbsg and ACTG175 anchor configurations, so identity is the expectation here, and 0.3.1 → 0.3.2 (the orientation commit) is the one link this measurement newly pronounces on.

**Why.** Larry wants direct confirmation that the applications do not change under the 0.3.1 → 0.3.5 revisions. The 08-31 pass rendered at 0.3.1, compared HTML against the tracked (Mac-era) baseline, and committed `_payloads_2026-09-01`'s reference set: `_payloads_2026-08-31/`. This task renders the same documents at **0.3.5**, stages **`_payloads_2026-09-01/`** (directory name per Larry's instruction), and compares **payload-to-payload at full precision** against the 08-31 set — the field-by-field comparison the workstream queued, now with a committed reference — plus the established HTML triage as the display-layer check.

---

## ⚠ CATEGORY — renders and comparisons only; nothing existing is modified or deleted

- **Zero `R/` edits.** No package, driver, or document source change of any kind. No version bump.
- **Existing payloads are read-only.** `_payloads_2026-08-31/` (both applications), any `_payloads_mac_*` directories, the tracked template payloads, and every other `.rds` in the tree are read, never written, never moved, **never deleted**. No `rm` of any `.rds` occurs anywhere in this task.
- **In-place render outputs are restored.** Rendering writes HTML and payload files at the documents' default paths, overwriting working-tree copies of tracked files. Before rendering, enumerate exactly which tracked paths the renders will touch (§2.3); after staging copies into the dated directories, restore every one of them with `git checkout -- <path>` so the working tree's tracked files finish **byte-identical** to HEAD.
- **Writes, exhaustively:** this task doc; the comparison/render-driver scripts; `quarto/applications/gbsg/_payloads_2026-09-01/` and `quarto/applications/actg175/_payloads_2026-09-01/` (payloads plus an `html/` subdirectory holding the new renders, so the verification bundle is self-contained); the report. *Final gate (§6): `git status --porcelain` shows only these paths.*

**Compute:** four document renders — the 08-31 pass ran under a ~60-minute ceiling at 0.3.1; the search portions are 2–4× faster at 0.3.5 while DINA/GRF/MR costs are unchanged, so **estimate 30–50 minutes, ceiling 60**. Comparisons add seconds. Toolchain per the 08-31 task verbatim: the RStudio-bundled quarto binary, plain Linux terminal.

**Unattended.** Gates stop the task; never ask, never work around. Provenance gate dirt-tolerant; stop only if dirt touches this task's own paths.

---

## 1. Provenance and first commit — GATE

Standard block; `packageVersion("forestsearch")` expect **0.3.5**. *GATE:* branch `feature/glm-extension`; the efficiency-workstream close (`e31e612d`, `af254e18`) and the 08-31 render-task commits in the log (identify them from the log and name them in the report). Copy this document into `dev/tasks/`, commit **alone**, report filenames and hash. Record the `vi.grf.min` default in force (expect `NULL`).

## 2. Transplant the 08-31 procedure — verify from source, author nothing

1. **Document list and render commands:** quote from the committed 08-31 task/report the exact four documents, the render invocations, and the toolchain lines. This task runs the same four — no additions (the OC applied document has its own record and is out of scope).
2. **Flag state:** confirm at HEAD that every LOO / K-fold / CV flag the 08-31 task set FALSE is still FALSE (including frozen_family's `run_loo`, flipped and committed on 08-31). Quote each. If any has drifted TRUE, STOP and report — do not flip anything in this task.
3. **The restore list:** from the 08-31 report and the documents' output declarations, enumerate every tracked path a render overwrites (HTML files, in-place payload files such as the template's tracked payload). This list drives the §3 restore step and the §6 gate.
4. **Reference-set census:** list the contents of both `_payloads_2026-08-31/` directories — filenames, sizes, and each payload's self-described `forestsearch_version` and `built_at`. *GATE:* the reference set exists and self-describes as the 0.3.1 vintage; if either directory is absent or the vintage is not 0.3.1, STOP and report what is actually there.
5. **Triage definitions:** quote the 08-31 VOLATILE / CLAUSE-ORDER / SUBSTANTIVE definitions verbatim for reuse in §5.

Commit the render-driver and comparison scripts as the second commit before any render runs.

## 3. Render at 0.3.5 and stage

Render the four documents with the 08-31 invocations. Then:

1. Copy every produced payload `.rds` into the matching `_payloads_2026-09-01/` directory; copy the four rendered HTML into `_payloads_2026-09-01/html/`.
2. Restore every §2.3 path with `git checkout -- <path>`; verify each is byte-identical to HEAD (`git diff --stat` empty for those paths).
3. Record wall-clock per document beside the 08-31 report's figures — the realized render-time effect of 0.3.2→0.3.5 on full documents is itself a finding for the record (the first full-document measurement of the efficiency workstream's changes; DINA/GRF/MR portions dilute the search speedup — say by how much).

## 4. Comparison A — payloads, field by field, full precision

Pair each `_payloads_2026-09-01/` file with its `_payloads_2026-08-31/` counterpart by name. For each pair, recursively compare the loaded objects and classify every leaf:

- **EXCLUDED-VOLATILE:** version strings (`forestsearch_version` 0.3.1 vs 0.3.5 — expected), `built_at`/timestamps, wall-clock/timing fields, session info. Enumerate every excluded field by path and print both values — the exclusion list is a deliverable, and an over-broad exclusion is a gate failure in spirit.
- **STRUCTURAL-VERSION:** fields present in one vintage only (result-object fields added or removed between 0.3.1 and 0.3.5 — e.g., any early-stop bookkeeping that postdates 0.3.1). List each by path with which side carries it. These are package-evolution facts, not result changes — but they must be listed, not silently skipped.
- **COMPARED:** everything else — `identical()` required. Tables, estimates, CIs, subgroup definitions, membership, `meta`, `extras`, MR outputs, consistency values.

*GATE:* any COMPARED leaf that is not `identical()` → report its path and both values verbatim at full precision, for every such leaf, then proceed to §5 (the HTML layer localizes the same difference in display terms) but trigger the §6 no-commit rule. The gbsg/actg175 selections themselves (harm-subgroup definitions) get an explicit named line each in the report regardless of outcome.

## 5. Comparison B — HTML against the tracked baseline

Diff each new render against the tracked HTML with the 08-31 triage taxonomy verbatim. Expected: VOLATILE hunks only (version strings, timestamps, render dates). Report the hunk census per document — VOLATILE / CLAUSE-ORDER / SUBSTANTIVE counts — and quote every SUBSTANTIVE hunk in full if any exists.

## 6. Commit rule and final gate

- **Clean** (§4: no COMPARED differences; §5: no SUBSTANTIVE hunks): commit both `_payloads_2026-09-01/` directories (payloads + `html/`) and the report as the final commit. The 0.3.5 payload set becomes the current-codebase reference, sitting beside — never replacing — the 08-31 set.
- **Not clean:** leave both `_payloads_2026-09-01/` directories **staged but uncommitted**; commit only the report (and scripts, already in). The differences go to Larry for adjudication; nothing is deleted either way — the staged directories remain on disk untouched.
- *Final gate, both branches:* `git status --porcelain` shows only this task's declared paths; every §2.3 restored path is byte-identical to HEAD; **no existing payload directory has been modified** (`git diff --stat HEAD -- '*_payloads*'` empty apart from the new dated directories).

## 7. Report

Beside the 08-31 render report (same directory; state the path): `REPORT_applications_render_0.3.5_2026-09-01.md`:

1. Provenance; the §2 transplant quotes, flag confirmations, restore list, and reference-set census with vintages.
2. Render walls beside the 08-31 figures, with the dilution statement.
3. Comparison A: the exclusion list with values; the STRUCTURAL-VERSION list; the COMPARED verdict per file — and the named lines for the gbsg and actg175 selected subgroups.
4. Comparison B: the triage census per document; SUBSTANTIVE hunks quoted in full if any.
5. Which §6 branch was taken and why; the final-gate outputs proving the tree state.
6. Ten-line verdict, opening with the one sentence Larry asked this task to produce: whether the analyses changed under the 0.3.1 → 0.3.5 revisions.

Commits, in order: (1) this task doc alone; (2) scripts; (3) report + the dated directories iff clean. Explicit paths; tree clean at close apart from any deliberately-staged-uncommitted directories under the not-clean branch. **No push. No `R/` edit. No deletion of anything, anywhere.**
