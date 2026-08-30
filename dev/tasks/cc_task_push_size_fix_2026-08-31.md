# CC task — unblock the push: rewrite `713cd93e` to exclude files over the size limit

**Date:** 2026-08-31 · **Repo:** `~/Documents/GitHub/forestsearch` · **Branch:** `feature/glm-extension`
**Type:** local git history rewrite (nothing is on the remote — verified by G2 before any rewrite). **No file contents change; no renders; nothing under `R/` or `quarto/` is edited.** Working-tree copies of excluded files are preserved on disk, untracked.
**Policy (Larry, 2026-08-31):** track only files **≤ 50 MB** in the dated payload directories; exclude larger files from git, keep them on disk, manifest them in the report.
**Report:** `dev/tasks/REPORT_push_size_fix_2026-08-31.md`
**Compute:** negligible.

## Context — complete

GitHub rejected the push of 7 local commits because at least one blob exceeds the 100 MB per-file server limit. The offending commit is expected to be `713cd93e` (`applications: 0.3.1 payloads under _payloads_2026-08-31 ...`). Deleting files in a new commit cannot fix this — the blobs remain in history — so `713cd93e` must be rewritten in place and the three later commits replayed. The excluded-file class is expected to be the compare_all full-comparison object(s); the `<stem>_payload.rds` files are designed small and the template's tracked 0.2.0 payload proves that class pushes fine.

## Gates — stop semantics

On any STOP: write the report with what was gathered, banner `STOPPED AT <gate>: <reason>` first, commit only the report (on top of the *unmodified* history if the rewrite has not begun; otherwise state the exact rewrite state and how to return to `backup/pre-size-fix-2026-08-31`), end. Do not ask questions.

- **G1 — topology.** Branch `feature/glm-extension`. `git log --oneline -8` shows, newest first: `4ca38eec`, `e872763d`, `713cd93e`, `d1f0e88e`, `728ccf5a`, `2767197c`, `dfc8f3c9`, then the pre-workstream base (record its sha as `BASE`). Any extra, missing, or reordered commit ⇒ STOP. Untracked scratch (`dev/tasks/_baseline_html_2026-08-31/`, `_diffs_2026-08-31/`, render logs) is expected; tracked modifications ⇒ STOP.
- **G2 — remote safety.** `git fetch origin`. None of the 7 shas may be reachable from `origin/feature/glm-extension` (`git merge-base --is-ancestor <sha> origin/feature/glm-extension` must fail for each) ⇒ otherwise STOP (something was pushed; rewriting is unsafe). Record whether `origin/feature/glm-extension == BASE`; if origin is ahead (e.g. the Mac capture landed), note it — handled in step 5.
- **G3 — diagnosis confirmed.** The size census (step 2) must find **at least one blob > 100 MB introduced by `713cd93e`**. If none ⇒ STOP: the push failure has a different cause; report every blob ≥ 25 MB across `BASE..HEAD` with paths, and the exact remote error text if available in Desktop, so the diagnosis can be redone.
- **G4 — post-rewrite verification.** Across `BASE..HEAD` after the rewrite: zero blobs > 100 MB; zero blobs > 50 MB under `quarto/applications/*/_payloads_2026-08-31/`; the excluded files still exist on disk and show as untracked; `git log` shows 7 commits with the original messages (the rewritten one's message amended as specified); working tree otherwise unchanged. Any failure ⇒ STOP (leave state; name the backup branch).

## Steps

### 1. Environment + G1 + G2
Record hostname, branch, the 8-line log with `BASE`, fetch result, and the per-sha ancestor checks.

### 2. Size census (drives everything; G3 here)
```
git rev-list --objects BASE..HEAD \
  | git cat-file --batch-check='%(objecttype) %(objectname) %(objectsize) %(rest)' \
  | awk '$1=="blob" && $3 > 26214400 {printf "%.1f MB  %s\n", $3/1048576, substr($0, index($0,$4))}' \
  | sort -rn
```
Report every blob ≥ 25 MB with its path and which commit introduced it (`git log --oneline --find-object=<sha>`). Define `EXCLUDE` = the set of paths in `713cd93e` whose blob size > 50 MB. Also record `du -h` of both `_payloads_2026-08-31` directories.

### 3. Safety branch
`git branch backup/pre-size-fix-2026-08-31` (must not already exist; if it does ⇒ STOP). Note in the report that Larry deletes it once the push succeeds.

### 4. Rewrite `713cd93e`
Interactive rebase scripted non-interactively:
```
GIT_SEQUENCE_EDITOR="sed -i 's/^pick 713cd93e/edit 713cd93e/'" git rebase -i BASE
git rm --cached <each path in EXCLUDE>
git commit --amend -m 'applications: 0.3.1 payloads under _payloads_2026-08-31 (dated; committed for cross-machine provenance; files >50 MB excluded from tracking: <basenames of EXCLUDE>)'
git rebase --continue
```
The later commits (`e872763d`, `4ca38eec` equivalents) replay clean — no path overlap with `EXCLUDE`. Any conflict ⇒ STOP (state, backup branch named).

### 5. If origin was ahead at G2
`git rebase origin/feature/glm-extension` (paths disjoint from the Mac capture; conflict ⇒ STOP). Otherwise skip.

### 6. G4 verification
Re-run the step 2 census on the new `BASE..HEAD` (or `origin/..HEAD` if step 5 ran); confirm the zero-blob conditions; `ls -la` each `EXCLUDE` path (on disk, untracked per `git status --porcelain`); `git log --oneline -8` for the new shas.

### 7. Task doc + report + close-out
Copy this task document into `dev/tasks/` and commit it, then write and commit the report. Print the final `git log --oneline -9` — the new branch tip sha is what the roxygen task's G1 note will be updated to reference. **No push** (Larry pushes via GitHub Desktop; the queue should now show these commits and clear the size check).

## Report skeleton
```
# REPORT — push size fix — 2026-08-31
[banner if stopped]
## 0 Environment, topology (G1), remote safety (G2)
## 1 Size census (every blob ≥ 25 MB, path, size, introducing commit; EXCLUDE set; du -h)
## 2 Rewrite log (backup branch; rebase transcript summary; amended message)
## 3 Post-rewrite verification (G4: censuses, on-disk/untracked confirmation, new shas old→new table)
## 4 Notes (no action taken): excluded-file manifest for cross-machine awareness; the open decision on slimming compare_all's saved object so future runs are trackable
```

## Hard rules
- No file contents are modified anywhere; the only changes are git history (which paths `713cd93e` tracks) and the two committed documents.
- Excluded files are never deleted from disk.
- The backup branch is created before the rewrite and never deleted by this task.
- Verify from source; no push; no follow-up work on §4 notes.
