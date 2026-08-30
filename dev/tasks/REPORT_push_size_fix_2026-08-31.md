# REPORT — push size fix — 2026-08-31

All gates passed. `713cd93e` rewritten in place to `8ad2c3f6`; the two later commits replayed clean; history now contains no blob over 50 MB in the dated payload dirs (and none over 100 MB anywhere in the range). The excluded file is intact on disk, untracked, byte-identical to the blob that was removed from history. Backup branch kept. **No push** — GitHub Desktop's queue should now clear the size check.

## 0 Environment, topology (G1), remote safety (G2)

- Host `pop-os`, branch `feature/glm-extension`.
- **G1 pass** — `git log --oneline -8`, newest first, exactly as specified:
```
4ca38eec applications: 0.3.1 render-and-compare report
e872763d applications: 0.3.1 reference renders (actg175 template+compare_all; gbsg frozen_family+multimethod)
713cd93e applications: 0.3.1 payloads under _payloads_2026-08-31 (dated; committed for cross-machine provenance)
d1f0e88e gbsg frozen_family: run_loo -> FALSE (core analyses only; LOO/CV off per decision 2026-08-31)
728ccf5a tasks: add applications render task (2026-08-31)
2767197c applications: read-only inventory ahead of 0.3.1 re-run (report)
dfc8f3c9 tasks: add applications inventory task (2026-08-31)
d25d7354 report — record commit hash
```
  `BASE = d25d7354` (`d25d7354fbe3f1a190fc7aa8af9e0c18f25ccefa`). `git status --porcelain` showed only the expected untracked scratch (`dev/tasks/_baseline_html_2026-08-31/`, `dev/tasks/_diffs_2026-08-31/`); no tracked modifications.
- **G2 pass** — `git fetch origin` completed; `origin/feature/glm-extension = d25d7354` — **equal to `BASE`**, so origin was not ahead and step 5 was skipped. Per-sha ancestor checks: all seven of `4ca38eec e872763d 713cd93e d1f0e88e 728ccf5a 2767197c dfc8f3c9` are **not** reachable from `origin/feature/glm-extension` (`git merge-base --is-ancestor` failed for each). Rewriting is safe; nothing local is on the remote.

## 1 Size census (G3)

`git rev-list --objects d25d7354..HEAD | git cat-file --batch-check … | awk '$3 > 26214400'` — every blob ≥ 25 MB in the range:

| size | blob | path | introduced by |
|---|---|---|---|
| **105.3 MB** (110,436,997 bytes) | `15156d3e152ed4c83f454f5e3275de7c71be67ae` | `quarto/applications/actg175/_payloads_2026-08-31/analysis_actg175_continuous_compare_all/comparison_continuous.rds` | `713cd93e` (`git log --find-object=15156d3e` → that commit only) |

That is the **only** blob ≥ 25 MB — exactly one blob > 100 MB, introduced by `713cd93e`. **G3 pass**; the expected diagnosis (compare_all's full-comparison object) is confirmed. The `<stem>_payload.rds` class is indeed small (largest 0.3 MB).

**`EXCLUDE`** = { `quarto/applications/actg175/_payloads_2026-08-31/analysis_actg175_continuous_compare_all/comparison_continuous.rds` } (one path).

`du -sh`: `gbsg/_payloads_2026-08-31` = **20K**; `actg175/_payloads_2026-08-31` = **106M** (of which `analysis_actg175_continuous_compare_all/` 106M, `template_actg175_continuous/` 280K; the gbsg subdirs 8K each). The dated dirs minus the excluded file total well under 1 MB.

## 2 Rewrite log

- Safety branch: **`backup/pre-size-fix-2026-08-31`** created at `4ca38eec` (did not previously exist). Larry deletes it once the push succeeds; this task never deletes it.
- Rebase, scripted non-interactively:
```
GIT_SEQUENCE_EDITOR="sed -i 's/^pick 713cd93e/edit 713cd93e/'" git rebase -i d25d7354
git rm --cached quarto/applications/actg175/_payloads_2026-08-31/analysis_actg175_continuous_compare_all/comparison_continuous.rds
git commit --amend -m 'applications: 0.3.1 payloads under _payloads_2026-08-31 (dated; committed for cross-machine provenance; files >50 MB excluded from tracking: comparison_continuous.rds)'
git rebase --continue
```
- Transcript summary: rebase stopped at `713cd93e` for edit; `git rm --cached` removed the one path from the index only (`rm '…/comparison_continuous.rds'`); amend succeeded; `git rebase --continue` replayed the remaining two commits without conflict — `Successfully rebased and updated refs/heads/feature/glm-extension.` (The later commits touch only the four `.html` and the report — no path overlap with `EXCLUDE`, as expected.)

## 3 Post-rewrite verification (G4)

Census re-run over `d25d7354..HEAD` after the rewrite:
- Blobs > 100 MB: **0**.
- Blobs > 50 MB under `quarto/applications/*/_payloads_2026-08-31/`: **0**.
- Largest blobs now in the range (all fine): `analysis_gbsg_survival_multimethod.html` 17.9 MB; `analysis_actg175_continuous_compare_all.html` 3.1 MB; `analysis_gbsg_survival_frozen_family.html` 2.3 MB; `template_actg175_continuous.html` 2.0 MB; `template_actg175_continuous_payload.rds` 0.3 MB.
- Excluded file on disk and untracked: `ls -la` → `110436997` bytes at the original path; `git status --porcelain` shows it as `??`. **Content byte-identical** to the removed blob: `git hash-object` on the working-tree file returns `15156d3e152ed4c83f454f5e3275de7c71be67ae`, the same blob id (the file's mtime moved 12:12 → 13:04 during the rebase's tree checkout; contents unchanged).
- Branch back on `feature/glm-extension`; 7 commits with original messages (the payload commit's message amended as specified). Sha map:

| old | new | message |
|---|---|---|
| `dfc8f3c9`, `2767197c`, `728ccf5a`, `d1f0e88e` | unchanged | (before the edited commit) |
| `713cd93e` | **`8ad2c3f6`** | `applications: 0.3.1 payloads under _payloads_2026-08-31 (dated; committed for cross-machine provenance; files >50 MB excluded from tracking: comparison_continuous.rds)` |
| `e872763d` | **`2c2406a9`** | `applications: 0.3.1 reference renders (actg175 template+compare_all; gbsg frozen_family+multimethod)` |
| `4ca38eec` | **`0e9d48fa`** | `applications: 0.3.1 render-and-compare report` |

**G4 pass.**

## 4 Notes (no action taken)

- **Excluded-file manifest** (for cross-machine awareness — this file exists only on `pop-os`, untracked):
  - `quarto/applications/actg175/_payloads_2026-08-31/analysis_actg175_continuous_compare_all/comparison_continuous.rds` — 110,436,997 bytes (105.3 MB), built 2026-08-30 12:12 by the 0.3.1 compare_all render; blob/content id `15156d3e152ed4c83f454f5e3275de7c71be67ae`; top-level names `combos,fs,ci_tab,plots,plot_grid,plot_combined,plot_combined_subsets,combined_skip_reason,console,diagnostics,errors` (readback in `REPORT_applications_render_0.3.1_2026-08-31.md` §7). Its small tidy companion `selected_subgroups_continuous.rds` (648 bytes) **is** tracked in `8ad2c3f6`.
- **Open decision, not acted on:** slimming compare_all's saved `comparison` object so future runs are trackable — the object carries eight full `forestsearch` fits (`fs`), per-combo plots/plot grids, and captured console text; dropping or thinning the heavy elements (e.g. `fs`, `plots`, `console`) at the `saveRDS` site would put it under the 50 MB policy. Left for Larry; no file was edited.
- The working-tree copies of everything else are exactly as the render task left them; no file contents were modified by this task (verified for the excluded file by blob-id equality; nothing else was touched).
