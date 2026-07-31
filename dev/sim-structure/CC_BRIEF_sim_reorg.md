# CC BRIEF — review and propose a structure for `quarto/simulations/gbsg_redux`

```
claude "Read dev/sim-structure/CC_BRIEF_sim_reorg.md and execute it."
```

**Propose, do not restructure.** Produce a design document and a migration plan.
Move no file, delete no file, edit no `.qmd`. Write only under
`dev/sim-structure/`.

The point is to get the design right before touching 204 files, several of which
carry results that cannot be regenerated cheaply.

---

## 1. The problem, measured

`quarto/simulations/gbsg_redux/` holds **204 entries**: 73 `.qmd`, 69 `.html`,
59 `.rds`, plus `betaHhat_truth.R` and two loose notes. Disk: 4.6M of `.qmd`,
**182M of `.html`**, 2.5M of `.rds`, all tracked in git.

The 73 `.qmd` files are not 73 analyses. Strip the `_batch_a_b` /
`_combine_a_b` suffix and they collapse to **22 distinct stems**, of which
roughly 11 are the León Table 2 design cells. One cell,
`fs_t1_t2_m1_h15_knoise3_n1500`, has **11 sibling files**.

Those siblings are near-copies. Measured pairwise:

* `..._n1000_batch_1_20.qmd` vs `..._batch_21_100.qmd` — **4 differing lines of
  1231**, all in `n_sims` and `sim_id_start`
* `h15_knoise3_n1500` siblings — **2 differing lines of 1231**

So the same ~1230-line document is stored dozens of times to vary two integers
that the document already derives its output filename from
(`sim_id_end <- sim_id_start + n_sims - 1L`).

**This is the cost.** Every fix to the simulation logic must be applied N times
or the copies diverge. That has already happened in this codebase: three copies
of one simulation drifted during recent work, one carrying an edit set the
others lacked, one carrying a render-blocking bug. The batch-per-file layout
makes that outcome structural rather than accidental.

---

## 2. What to review

Read the directory and the machinery, not just the file names.

1. **The cell space.** Which stems are live design cells, which are supersedes
   (`fsNEW_`, `_sim_*`, `mr_coverage_sweep_*`), and which are one-off probes.
   `fsNEW_` in particular looks like a name that meant something once.
2. **The batching machinery.** How `sim_id_start` / `n_sims` / `rds_stem` /
   `rds_path` / `combine_glob` / `combine_files` interact; what `run_mode`
   values exist; how a combine run finds and validates its inputs.
3. **The agreement check.** Which `meta` fields are compared, which are recorded
   but not compared, and what happens to batches written before a field existed.
   Note the deliberate design decision already made: provenance fields
   (`sg_focus`, `forestsearch_version`, `n_workers`, …) are **recorded but not
   gated**, because gating would `stop()` on every pre-existing batch.
4. **Result storage.** `_res_a_b.rds` per batch, `_combined_a_b.rds` per cell.
   Which are inputs to something downstream, which are terminal.
5. **Cross-directory consumers.** Results here are read from elsewhere —
   `quarto/GuoHe/guohe_sec52.qmd` reads `n_sel`,
   `quarto/simulations/gbsg/fs_subgroup_id_sweep_n500.qmd:266` reads `covs`.
   Find the full set. **Any proposal that moves or renames an `.rds` must
   account for these**, and they are the main constraint on the design.
6. **The `.html` renders.** 182M, tracked. Establish whether any is referenced
   by anything, or whether they are byproducts kept by default.

---

## 3. What the proposal must address

Give a concrete target layout with real paths, not a sketch. It must answer:

**One document per cell, or one per batch?** The measured duplication argues for
one parameterised document per cell, batches driven by parameters rather than by
copies. If you propose that, say exactly how a batch is launched (Quarto `-P`,
a driver script, a params block) and how the rendered `.html` is kept per batch
without overwriting — Quarto writes `X.qmd` → `X.html`, which is the constraint
that made per-batch copies attractive in the first place.

If you conclude per-batch files are right after all, say why the drift cost is
acceptable and how it is mitigated.

**Where results live.** Batch `.rds` and combined `.rds` currently sit beside the
sources. Consider separating source from output, and say what that costs the
cross-directory consumers in §2.5.

**What git tracks.** 182M of `.html` is the single largest item. Options: keep,
`.gitignore` and regenerate, or keep only combined-run renders. State the
tradeoff — a tracked render is a record of what was actually produced, and this
project has already needed exactly that twice.

**What retires.** Name the stems you believe are superseded and the evidence.
Do not delete them; list them for the maintainer.

**Naming.** The recently agreed convention for this file family is
`sim_<engine>_<focus>_<procedures>_<design>_<cell>`, with batch range **out** of
the document name because the `.rds` already carries it
(`..._res_1_100.rds`). Apply it, and flag anything it does not cover.

**Migration.** An ordered sequence of steps that leaves the repo working at each
one. Identify anything unrecoverable — an `.rds` that no longer has a source
that can regenerate it — and treat those as fixed points to be preserved, not
migrated.

---

## 4. Constraints

* Existing `.rds` files are the expensive artifacts. Some represent many hours of
  compute. Any step that would orphan or invalidate one must be called out.
* Pre-existing batches carry `meta$gate_draws` rather than `mr_draws`, so they
  already cannot pool with new ones. Note where that bites.
* `betaHhat_truth.R` is `source()`d from the working directory by the simulation
  documents. Any layout that moves the documents must keep that resolvable.
* The recorder schema is shared across roughly 180 simulation documents
  repo-wide. Do not propose schema changes here.

---

## 5. Deliverable

`dev/sim-structure/SIM_REORG_PROPOSAL.md`:

1. **Inventory** — what is actually in the directory, grouped by role, with the
   duplication quantified.
2. **Findings** — what the current structure costs, with evidence.
3. **Proposed layout** — real paths, and how a batch is run under it.
4. **Migration plan** — ordered, each step leaving a working repo.
5. **Retirement candidates** — listed, not acted on.
6. **Open questions for the maintainer** — anything where the right answer
   depends on intent rather than evidence.

Where a recommendation is a judgement call rather than a measurement, say so and
give the alternative.

Report and stop.
