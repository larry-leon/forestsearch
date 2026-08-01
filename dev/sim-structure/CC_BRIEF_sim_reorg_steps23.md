# CC BRIEF — proposal Steps 2 and 3: relocate artifacts, untrack batch renders

```
claude "Read dev/sim-structure/CC_BRIEF_sim_reorg_steps23.md and execute it."
```

Continues `dev/sim-structure/SIM_REORG_PROPOSAL.md`. Steps 0 and 1 are complete
(manifest written; the κ drift resolved by putting κ in the `sim_fs_maxcons_*`
template and archiving the six κ-less combine documents). This brief executes
**Step 2 and Step 3** as written in the proposal.

**No file is deleted.** Every operation is `git mv`, `git rm --cached`
(untrack, file stays on disk), or a path-string edit. If a step appears to
require a real deletion, stop and report.

---

## 1. Scope, stated plainly

| proposal step | this brief |
|---|---|
| Step 0 — manifest | done |
| Step 1 — reconcile κ drift | done (κ in the maxcons template; six archived) |
| **Step 2 — move `.rds` into `results/`** | **yes** |
| **Step 3 — untrack batch renders** | **yes** |
| Step 4 — introduce the folded template | no |
| Step 5 — retire the 60 `t1_t2` copies | no |
| Step 6 — fold the `mr_coverage_sweep_*` family | no |
| Step 7 — resolve the `_sim_*` fragments | no |

Steps 4–6 are excluded for one reason: the maintainer is keeping the `t1_t2`
family in place until the `sim_fs_maxcons_*` setup is settled, and those steps
all end in retiring it. Step 7 needs a maintainer answer. Steps 2 and 3 touch
neither.

---

## 2. Step 2 — move the artifacts into `results/`

Create `quarto/simulations/gbsg_redux/results/` and `git mv` **every top-level
`.rds`** into it. **Basenames do not change** — `results/` is a prefix and
nothing more. Derive the file list from the directory rather than from a list in
this brief; the directory has changed twice since the proposal was measured.

Do **not** touch `mr_sweep/`. It has its own `<run_tag>/` structure and is out
of scope.

### The path edits

Three classes of reader, all updated in the same commit:

1. **The `sim_fs_maxcons_*` template.** `rds_path`, `combine_glob` and
   `combined_path` are built from `rds_stem`. Add one `file.path()` component so
   they resolve into `results/`, and make the directory on write with
   `dir.create(..., showWarnings = FALSE, recursive = TRUE)` so a fresh clone
   works. **This is required, not optional** — the maintainer is about to run a
   production baseline from this template, and it must write into the new
   location rather than the old one.
2. **The `t1_t2` documents.** Same edit, same reason. They are staying in place
   and must keep working.
3. **The `dev/` consumers** — the three lines the proposal identified in
   `phase0_cells.R`, `phase0_summarise.R` and `test_combine_guard.R`, plus
   whatever `dev/sim-check/fixcheck/` uses to read
   `fs_maxcons_fb_mr_m1_h10_knoise0_n500_res_1_100.rds` (proposal §2.5). Grep
   for the full set rather than trusting that count; the proposal measured it
   before two rounds of moves.

Report the exact edits made, per file.

### Verification

* Re-run `dev/sim-structure/rds_manifest.csv` and require an **exact match** on
  content — same file count, same `nrow(results)`, same `meta` for every
  artifact. Only the path column changes. Any other difference is a fault.
* Confirm `git status` shows renames and zero deletions.
* Confirm every `.rds` referenced anywhere in `quarto/` or `dev/` resolves to a
  file that exists. A stale path here fails at `readRDS`, not at parse.

### Fixed points

These are preserved by the move, never regenerated (proposal §2.5):

* `fs_t1_t2_m1_h20_knoise3_n1500_combined_1_250.rds` — no source reproduces 250 rows
* `fs_maxcons_fb_mr_m1_h10_knoise0_n500_res_1_100.rds` — its only source now sits
  at smoke-test values
* the eight `*_combined_1_500.rds` — reproducible only while the glob set is
  unchanged, which this step preserves

Confirm all ten are present in `results/` after the move.

---

## 3. Step 3 — untrack the batch renders

200 MB of `.html` is the largest item in the directory, and `embed-resources:
true` means size is nearly independent of content.

1. Add `quarto/simulations/gbsg_redux/renders/` to `.gitignore`.
2. `git mv` the **combine-mode renders** into `renders/combined/` and
   `git add -f` them — these stay tracked, because a combine render is the
   record of a pooled result.
3. `git rm --cached` the remaining batch renders. **They stay on disk.** Confirm
   that: the working tree must still contain every `.html` afterwards.

Determine combine-vs-batch by the `run_mode` **inside** each source document,
never by its filename. Six documents named `..._combine_*` in the
`h20_knoise3_n1500` family are `run_mode <- "batch"`; that mismatch is why they
were archived and then moved back. Report the split you derive and the rule you
used.

Renders in `legacy/` stay where they are — that directory is a closed archive.

This step is fully reversible and touches no artifact. Note in the report that
untracking stops growth but does not reclaim history.

---

## 4. Deliverable

Report directly, no separate document:

1. Step 2: files moved, path edits per file, manifest match, fixed points
   confirmed, `git status` rename/deletion counts.
2. Step 3: combine/batch split and how it was derived, files untracked, working
   tree still complete, `.gitignore` line added.
3. Top-level entry count before and after, and tracked size before and after.
4. Any `.rds` path reference you found that this brief did not anticipate.
5. Anything that did not go as described.

Commit Steps 2 and 3 as **two separate commits** — the artifact move and the
render untracking are independently revertible, and should stay that way.
Report and stop.
