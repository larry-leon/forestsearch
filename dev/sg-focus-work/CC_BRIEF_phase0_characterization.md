# CC BRIEF — Phase 0: characterise the existing simulation runs (READ-ONLY)

**Run first, before any implementation work.**

```
git checkout -b feature/sg-focus-maxeffcons
claude "Read dev/sg-focus-work/CC_BRIEF_phase0_characterization.md and execute it."
```

This is a **measurement** task. Its output decides whether the 500-replicate
simulation runs can be re-described in the new vocabulary or must be rerun.
Nothing is implemented here.

---

## 1. Why this must run first

After the implementation brief lands, `sg_focus = "hr"` no longer early-stops
in **any** configuration -- `forestsearch_main.R:1200` resets
`stop_threshold <- NULL` unconditionally for its listed foci and only *warns*
when the user set it explicitly. So the behaviour being characterised here
ceases to be reachable. **The opportunity does not come back.**

---

## 2. Hard constraints

**Do not modify:** `R/`, `NAMESPACE`, `DESCRIPTION`, `man/`, `tests/`,
`vignettes/`, `quarto/`.

Write **only** under `dev/sg-focus-work/`.

Hash every file under the installed package's source directory and the git
worktree's `R/` before and after the run; `stop()` on any difference and report
it. (`dev/efficiency-eval/R/00_guard.R` from the earlier evaluation has a
working `fs_hash_sources()` / `fs_guard_verify()` pair -- reuse it rather than
rewriting.) Record the guard verdict in the findings.

`git status` at the end must show changes only under `dev/sg-focus-work/`.

Install with `devtools::install(dependencies = FALSE)`, not `load_all()` --
`future` workers are separate processes that see only installed packages.

---

## 3. What the simulations actually ran

`sg_focus = "eff"` (normalised to `"hr"` at `forestsearch_main.R:1081`) with
`stop_threshold = 0.95` **inherited, never passed**, against
`pconsistency.threshold = 0.90`. Candidates are enumerated HR-descending
(`sort_subgroups_preview()`, `setorder(-HR, K)`) and, with
`batch_size_parallel = 1L`, evaluated one at a time.

That yields a **per-replicate hybrid**, not a single rule:

* **branch 1** -- some candidate reaches `Pcons >= 0.95` -> the loop halts and
  returns the first candidate in HR-descending order clearing 0.95;
* **branch 2** -- none does -> all candidates are evaluated, then
  `sort_subgroups()` applies `(-Pcons, -hr, K)` and returns the
  **max-consistency** candidate.

The intended rule is *argmax effect among candidates with `Pcons >= 0.90`*.
Branch 1 approximates it but skips qualifying candidates in `[0.90, 0.95)`.
Branch 2 is a different rule entirely.

**Verify this account before relying on it.** It comes from source reading in a
chat session, not execution:

| # | Claim | Location |
|---|-------|----------|
| V1 | `sg_focus` normalised before the line-1200 reset and before `subgroup.consistency()` | `forestsearch_main.R:1081` |
| V2 | line-1200 reset does **not** list `"hr"`, so `stop_threshold` stays 0.95 | `forestsearch_main.R:1200-1219` |
| V3 | batch guard keys on `sg_focus %in% c("hr","minSG")` -> `batch_size_parallel <- 1L` | `subgroup_consistency_main.R:836-841` |
| V4 | `evaluate_subgroup_consistency()` returns `NULL` below `pconsistency.threshold` | `subgroup_consistency_helpers.R:1476` — condition read, **body not** |

If any is false, stop and report; the framing above is then wrong.

---

## 4. Location and coverage

Simulation `.qmd` files and their `.rds` outputs both live in
`quarto/simulations/gbsg_redux/`. The DGM builder is a `forestsearch` package
function; the `.qmd` files show its usage.

**Run every cell used in the manuscript**, or at minimum both extremes. The
branch mix is expected to differ across the grid: `Pcons = P(HR > hr.consistency)`
rises with both signal strength and sample size, so high-`h` / large-`n` cells
(e.g. `h20_knoise3_n1500`) should clear 0.95 at candidate 1 almost always,
while low-`h` / small-`n` cells (e.g. `h10_knoise0_n500`) are where candidates
land in `[0.90, 0.95)` and branch 2 can fire.

**Characterising one cell does not settle the others.** Report per cell. The
low-`h` / small-`n` cells matter most -- likeliest to diverge, and they carry
the type-I-error and FDR claims.

---

## 5. Method

Re-run the **search stage only** on the same seeds (`seedit = seed_base + sim_id`,
`seed_base = 8316951`), skipping the Tier-1 bootstrap -- the wall-clock
bottleneck. Identify the toggle from the `.qmd`; `nb_boots` and the
`debias_gate` block are the relevant controls. Set `debias_gate = FALSE`.

**Seed-fidelity check first.** On ~10 replicates, confirm the search-only pass
reproduces the `sg_def` values already stored in the saved bundles. (Confirmed
present: each combined `.rds` holds a `results` data.frame with `sg_def`,
`n_sel`, `detected` and effect columns, 500 rows per cell.) **If it does not
reproduce, stop and report** -- everything downstream is meaningless without
seed fidelity. Note that the DGM builder link is unverified; if `sg_def` values
mismatch, suspect the DGM call before suspecting the seeds.

Record per replicate:

* `early_stop_triggered`, `n_candidates_evaluated`
* `Pcons` of the selected candidate
* max `Pcons` among candidates clearing `pconsistency.threshold`
* HR-rank of the selected candidate among qualifiers

The last two together identify replicates where a **higher-effect qualifying
candidate was passed over** -- exactly the cases where the intended rule and
the as-run behaviour diverge.

---

## 6. The decision this drives

Per cell:

* **branch 2 never fired, and the selected candidate was always HR-rank 1
  among qualifiers** -> the runs already implement the intended rule.
  `sg_focus = "maxeffCons", stop_threshold = 0.95` will reproduce them exactly.
  **Footnote, no rerun.**
* **otherwise** -> report the divergence rate and let the maintainer decide.

---

## 7. Deliverable

`dev/sg-focus-work/PHASE0_FINDINGS.md`:

* guard verdict (files hashed, unchanged);
* V1-V4 verification results;
* seed-fidelity result;
* per-cell table: branch-1 rate, branch-2 rate, divergence rate
  (selected candidate not HR-rank 1 among qualifiers), and the distribution of
  selected `Pcons`;
* a per-cell **footnote / rerun** recommendation;
* what could not be determined, and why.

Do not implement anything. Do not open a PR. Stop after writing the findings.
