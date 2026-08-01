# CC BRIEF — scoped reorganization of `quarto/simulations/gbsg_redux`

```
claude "Read dev/sim-structure/CC_BRIEF_sim_reorg_execute.md and execute it."
```

This executes **part** of `dev/sim-structure/SIM_REORG_PROPOSAL.md`. Most of that
proposal is deliberately deferred — see §5.

**No file is deleted in this brief.** Every operation is `git mv` or a drop-in
replacement. If any step seems to require a deletion, stop and report instead.

---

## 1. Answer to the proposal's §6.1 (the blocker)

The proposal correctly identified the κ / β(H) split as blocking the fold. It is
now settled, and the answer is not "backfill the six":

* **`sim_fs_maxcons_fb_mr_*` becomes the model template**, and it now carries the
  κ treatment. The maintainer has already applied it — see §2.
* **The combine documents lacking κ are archived**, not backfilled. They belong
  to the `t1_t2` generation, which is being superseded rather than repaired.

The count in the proposal was off: the split is **4 with κ / 12 without**, not
4/6. The proposal treated the `h20_knoise3_n1500` family as one cell; it has
eight combine documents of its own. Verified by
`grep -ci kappa` across all 16 combine documents — the split is exactly binary,
1291 lines (κ) vs 1231 lines (no κ).

---

## 2. Task A — drop in the revised template. **Do not re-implement it.**

The maintainer has revised
`quarto/simulations/gbsg_redux/sim_fs_maxcons_fb_mr_m1_h10_knoise0_n500_batch_1_100.qmd`
to carry the κ treatment: **1347 → 1415 lines, 11 → 12 chunks**, parses clean.

**Replace the repo copy with the maintainer's copy. Write no κ code yourself.**
Earlier in this project a parallel edit forked this exact file and the
reconciliation was expensive.

Confirm after the drop-in — report the count for each:

| marker | expected |
|---|---|
| chunk `beta-true-region` | 1 |
| `betaH_H` | ≥ 3 |
| `is_or <- blk$Estimator == "oracle"` | 2 (mean and median tables) |
| `kappa_txt` | ≥ 3 |
| `sub_missing(columns = c(b_betaHhat, Cov_betaHhat)` | 2 |
| `mr_avail` | 3 (from the previous round) |

Then verify statically only — **do not render, do not run a simulation**:

1. Extract the chunks and `parse()` them.
2. Source the setup chunk in a clean session (stub `library()`/`source()`); it
   must complete with no `object ... not found`.
3. `sprintf` arity across the file, including the two new κ footnotes.
4. Confirm the four κ pieces are symmetric across both tables — anything present
   for the mean table and absent for the median table is a finding.

One design point to check rather than assume: the new `beta-true-region` chunk
guards on `!identical(run_mode, "batch")` and rebuilds the eval frame only in
combine mode, reusing `.td_chk` from `build-dgm` in batch mode. Confirm
`.td_chk`, `harm_col`, `outcome_name`, `event_name`, `treat_name`, `n_super`,
`analysis_time` and `cens_adjust` are all in scope at that point in both modes.

---

## 3. Task B — archive the 12 combine documents without κ

Create `quarto/simulations/gbsg_redux/legacy/` and `git mv` **24 files** into it
— each `.qmd` *and its `.html` render*, so source and render stay paired:

```
fs_t1_t2_m1_h10_knoise0_n1000_combine_1_400
fs_t1_t2_m1_h10_knoise0_n1000_combine_1_500
fs_t1_t2_m1_h15_knoise0_n500_combine_1_500
fs_t1_t2_m1_h15_knoise3_n1500_combine_1_500
fs_t1_t2_m1_h20_knoise3_n1500_combine_1_100
fs_t1_t2_m1_h20_knoise3_n1500_combine_101_150
fs_t1_t2_m1_h20_knoise3_n1500_combine_151_200
fs_t1_t2_m1_h20_knoise3_n1500_combine_1_250
fs_t1_t2_m1_h20_knoise3_n1500_combine_201_250
fs_t1_t2_m1_h20_knoise3_n1500_combine_251_350
fs_t1_t2_m1_h20_knoise3_n1500_combine_351_400
fs_t1_t2_m1_h20_knoise3_n1500_combine_1_500
```

Verify each has both extensions before moving. Confirm `git status` shows 24
renames and **zero deletions**.

### Do not move

* **The 4 combine documents that carry κ** —
  `fs_t1_t2_m1_h10_knoise0_n500_combine_1_500`,
  `fs_t1_t2_m1_h15_knoise0_n1000_combine_1_500`,
  `dina_t1_t2_m1_h10_knoise0_n500_combine_1_500`,
  `grf_t1_t2_m1_h10_knoise0_n500_combine_1_500`.
* **Any `.rds`.** In particular the nine `*_combined_*.rds` are the *outputs* of
  these documents and stay at the top level. Moving them would break the
  cross-directory consumers the proposal identified.
* **Any `*_batch_*` document.** The maintainer is keeping the `t1_t2` batch
  family in place until the maxcons setup is settled.

### Write a `legacy/README.md`

Short, stating: these are `t1_t2`-generation combine documents superseded by the
`sim_fs_maxcons_*` template; they lack the κ / β(H) treatment; they are archived
for reference, **not runnable in place** — they `source("betaHhat_truth.R")` and
read `.rds` by relative path, so re-rendering from `legacy/` would need path
adjustment; and their combined `.rds` outputs remain at the parent level.

Note in it that this leaves four cells (`h10_knoise0_n1000`, `h15_knoise0_n500`,
`h15_knoise3_n1500`, `h20_knoise3_n1500`) with no combine document at top level.
That is intended.

---

## 4. Protection: manifest first

Before any move, write `dev/sim-structure/rds_manifest.csv` — proposal Step 0.
For every `.rds` under `quarto/simulations/gbsg_redux/` (including `mr_sweep/`):
path, `nrow(results)`, and every `meta` field.

Re-run it after the moves and require an exact match. No `.rds` is moved by this
brief, so any difference is a fault.

---

## 5. Explicitly deferred — do not do these

From `SIM_REORG_PROPOSAL.md`:

* **Step 2, moving `.rds` into `results/`.** It changes `rds_path` in the
  template the maintainer is still settling. Revisit after the maxcons baseline
  re-run, so the template is not moving underneath them.
* **Step 3, untracking the batch renders.**
* **The 60 → 1 fold of the `t1_t2` family.** The maintainer is keeping those
  files until the maxcons setup is settled.
* `cells.yml`, `run_batch.sh`, `params:` blocks, any `renders/` directory.
* Retiring `fsNEW_`, the `_sim_*` fragments, or the `mr_sweep/*.zip`.

---

## 6. Deliverable

Report directly, no separate document:

1. Task A: the marker counts, and the four static checks.
2. Task B: files moved, `git status` rename/deletion counts, README written.
3. Manifest: before/after match.
4. Top-level entry count before and after.
5. Anything that did not go as described.

Commit the drop-in, the moves and the manifest together. Report and stop.
