# CC task — the complete applications at 0.3.5: LOO and K-fold exercised, main analysis guarded

**File:** `dev/tasks/cc_task_applications_complete_render_2026-09-01.md` · **Issued:** 2026-09-01 by chat, commissioned by Larry
**Repo:** `~/Documents/GitHub/forestsearch` (pop-os) · **Branch:** `feature/glm-extension`
**Transport:** `~/Downloads` → `dev/tasks/`, committed first (§1). Hyphen-stripped stem `cc_task_applications_complete_render_20260901.md` is the same file.
**Reads (in-repo):** `dev/tasks/cc_task_applications_render_0.3.5_2026-09-01.md` and `REPORT_applications_render_0.3.5_2026-09-01.md` (the just-closed comparison pass — toolchain, restore mechanics, size precedents, and the committed comparator at its revised form `7fef766b`); the 08-31 task/report pair (the flag decisions being reversed here, and the inventory's LOO/CV cost estimates); `dev/tasks/render_driver_applications_0.3.5_2026-09-01.sh` and `compare_payloads_applications_0.3.5_2026-09-01.R` (transplant sources).

**Why.** The 08-31 pass set every LOO and K-fold flag FALSE so comparisons would be like-for-like and cheap, and the 09-01 pass inherited that. Larry now wants the **complete** documents — cross-validation content included — rendered at the current codebase as a refresh, explicitly **not** for comparison with prior work. Two things must still hold: the compute is projected before it is spent (LOO cost is unknown until the census — an iteration may be a full search or a cheap membership recompute), and the headline analysis must not change silently — the just-verified `_payloads_2026-09-01` set is the reference for a shared-content invariance guard, so any RNG-ordering effect of enabling the flags surfaces as a reported finding, never as a quiet drift in the canonical documents.

---

## ⚠ CATEGORY

- **Zero `R/` edits.** No package, driver, or comparator-logic change. No version bump.
- **Edits `.qmd` flag lines only** — knob-level, per the direct-edit rule: the exact LOO/K-fold/CV flags the §2 census identifies, flipped by `sed` with before/after quoted, committed as a configuration change (the symmetric reversal of the 08-31 flip `d1f0e88e`, rationale in the commit message: complete-document source must reproduce its rendered output).
- **Existing payloads read-only; nothing deleted.** `_payloads_2026-08-31/`, `_payloads_2026-09-01/`, parked copies, manifests, tracked template payloads: read, never written, never removed. No `rm` of any `.rds`. New outputs land only in **`_payloads_2026-09-01_complete/`** directories (plus `html/` subdirs), the parked-oversize location, and the report.
- **Tracked HTML disposition is gate-controlled (§6):** on the clean branch the complete renders become the tracked documents (and unflagged documents get their verified 0.3.5 HTML promoted by copy from the committed `_payloads_2026-09-01/html/` bundle, so all four finish at a coherent vintage); on the dirty branch every tracked path is restored via `git checkout` and nothing canonical changes.

**Compute:** unknown until §4's probe — **ceiling 3 hours of renders, enforced by a projection gate before anything long launches.** If the projection exceeds the ceiling, STOP and report the table; Larry raises or declines. Probes themselves are minutes.

**Unattended.** Gates stop the task; never ask, never work around. Provenance gate dirt-tolerant; stop only if dirt touches this task's own paths.

---

## 1. Provenance and first commit — GATE

Standard block; `packageVersion("forestsearch")` expect **0.3.5**. *GATE:* branch `feature/glm-extension`; the 09-01 comparison commits (`84d34a08`, `1ada23eb`, `7fef766b`, `86cbd5bb`) and the park/manifest follow-on commit in the log (name its hash). Copy this document into `dev/tasks/`, commit **alone**, report filenames and hash. Record the `vi.grf.min` default in force.

## 2. Census from source — what "complete" actually means, per document

1. **Flag census.** For each of the four documents, quote every LOO / K-fold / CV flag line with its current value and file:line — `run_loo`, any `run_kfold`/`run_cv`/`K`-fold knobs, and anything the 08-31 task's flag list named. Classify each document: **FLAGGED** (gains content when enabled) or **UNFLAGGED** (no such machinery — renders identically regardless).
2. **What one unit costs.** For each flagged mechanism, quote the implementation entry point and establish from source what a single LOO iteration and a single fold actually execute — full `forestsearch()` search, consistency-only recompute, membership refit, or other. Quote the loop bounds (n per document; K).
3. **Seeding.** Quote how the LOO/CV sections seed — document-level `set.seed`, derivation from `seedit`, or unseeded. If unseeded, record it: the complete render is then a snapshot rather than a reproducible artifact, and the report must say so plainly.
4. **Output paths.** Enumerate from source every file the enabled machinery writes — `cv_out_*` payloads, anything under `quarto/GuoHe/`, additions to the document payload bundles — with expected locations. These drive §5's staging and the size policy.
5. **The 08-31 estimates.** Quote the inventory's LOO/CV render-cost estimates verbatim; they anchor §4's projection alongside the fresh probes (divided by the measured 0.3.5 speedups where the mechanism is search-bound).
6. **Ordering relative to the main analysis.** Establish from source whether the LOO/CV sections execute after the main analysis completes and whether they share the RNG stream — this is the mechanism question behind the §6 guard, and the report states the answer either way.

Transplant the render driver and the comparator (from their committed 09-01 forms, changing only the named lines: document set, output names, `_complete` paths, and a NEW-CONTENT classification for payload fields present only on the complete side). *GATE:* show the diffs; anything beyond the named lines is a stop. Commit scripts as the second commit.

## 3. Flag flips — committed configuration change

Flip every §2.1 flagged knob to its enabled value by `sed`, quote before/after for each, run nothing else, and commit as the third commit (message citing the 08-31 flip being reversed and Larry's completeness instruction). *GATE:* `git diff` for this commit touches only the quoted flag lines.

## 4. Probe and projection — the compute gate

For each flagged document: time **one** LOO iteration and **one** fold in isolation (a scratch driver invoking the §2.2 unit, or a truncated-bounds dry run if the document structure allows nothing cleaner — say which), then project: units × per-unit + the document's known flags-FALSE render wall. Present the projection table (per document and total). *GATE:* total ≤ 3 hours → proceed. Total > 3 hours → STOP, commit the projection table in a short interim report, leave the flag flips committed (they are the intended end-state configuration regardless), and end the task — Larry decides whether to raise the ceiling.

## 5. Render and stage

Render the flagged documents with the transplanted driver (unflagged documents are **not** re-rendered). Then:

1. Copy every produced payload — including the new `cv_out_*` / LOO artifacts — into the matching `_payloads_2026-09-01_complete/` directory; copy rendered HTML into its `html/` subdir.
2. **Size policy, mirroring the standing precedents:** any single file > 50 MB is parked (moved, never deleted) to `~/fs_parked_2026-09-01_complete/` with a manifest entry (filename, bytes, sha256 computed before the move, `forestsearch_version`, `built_at`) committed in the `_complete` directory. Verify the new `.gitignore` rule covers `comparison_continuous.rds` under the `_complete` path too; if its pattern misses, extend it by one line in the same style and say so.
3. Leave the fresh tracked-path overwrites **in place for now** — §6 decides their fate.

## 6. The main-analysis invariance guard — GATE, and the branch

Run the transplanted comparator: `_payloads_2026-09-01_complete/` against the committed `_payloads_2026-09-01/` reference, with the taxonomy extended by **NEW-CONTENT** (fields present only on the complete side — LOO tables, CV results, their meta; expected and listed, not alarming). The standing classes keep their meanings; EXCLUDED-VOLATILE per the revised comparator.

- **Clean** (every shared COMPARED leaf `identical()`): the complete render leaves the headline analysis untouched. Commit, as the final commit: both `_complete` directories (+ manifests), the complete HTML **in place** as the tracked documents, the promoted 0.3.5 HTML for unflagged documents (copied from the committed `_payloads_2026-09-01/html/`), and the report. All four tracked documents now render-reproduce from committed source at 0.3.5.
- **Dirty** (any shared COMPARED leaf differs): restore **every** tracked path via `git checkout` (documents keep their current canonical state), commit the `_complete` directories and the report only, quote every differing leaf with both values at full precision, name the §2.6 mechanism if it explains the drift, and STOP — whether the complete configuration becomes canonical with changed headline numbers is Larry's adjudication, made with the values in hand.

*Final gate, both branches:* `git status --porcelain` shows only this task's declared paths (plus the deliberately-untracked oversize file(s) if the ignore rule needed no extension); every pre-existing payload directory unmodified (`git diff --stat HEAD~1 -- '*_payloads_2026-08-31*' '*_payloads_2026-09-01/*'` empty); nothing deleted anywhere.

## 7. Report

`dev/tasks/REPORT_applications_complete_render_2026-09-01.md`, beside its siblings:

1. Provenance; the §2 census — flags quoted, unit semantics, loop bounds, seeding verdict, output-path enumeration, the 08-31 estimates, and the §2.6 ordering/RNG answer.
2. The flag flips (before/after), and the transplant diffs.
3. §4's probe timings and projection table, and which side of the ceiling it landed.
4. Render walls per document; the staging census with sizes; park/manifest entries; the ignore-rule coverage finding.
5. §6's guard result: the NEW-CONTENT inventory (what completeness added, per document), the shared-leaf verdict, and which branch was taken — with full quoted values on the dirty branch.
6. The tracked-HTML end state, per document, with vintage.
7. Ten-line verdict, opening with whether the complete documents are now canonical and whether the headline analysis moved.

Commits, in order: (1) this task doc alone; (2) transplanted scripts; (3) the flag flips; (4) the §6-branch final commit (or the §4 interim report if the ceiling stopped it). Explicit paths; tree clean at close apart from deliberately-ignored oversize files. **No push. No `R/` edit. No deletion of anything, anywhere.**
