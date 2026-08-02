# Handoff — simulation structure work, `mr_coverage_sweep` family

**Purpose of the next session:** evaluate whether the `mr_coverage_sweep_*`
files need revision, and if so, revise and/or fold them. This is Step 6 of the
reorganization proposal, brought forward.

---

## Who and where

* Larry F. León, biostatistician at Merck; `forestsearch` R package, v0.2.0.
* Repo `larry-leon/forestsearch`, working branch **`feature/mr-in-replicates`**,
  base branch `master` (not `main`).
* Files under discussion:
  `quarto/simulations/gbsg_redux/mr_coverage_sweep_*.{qmd,html}` and the
  `quarto/simulations/gbsg_redux/mr_sweep/<run_tag>/` result tree.

## How the work is split (important)

Two agents, and they must never edit the same file:

* **This chat assistant** — edits one file at a time in a sandbox with R 4.3.3,
  verifies (parse + execute the setup chunk + behavioural tests), and returns
  the file in the artifact panel. Every response that creates or edits a file
  ends with `present_files`. Larry downloads from the panel.
* **CC (Claude Code)** — runs, renders, commits, and does multi-file / repo
  operations on Larry's machine (Pop!_OS, ~127 cores).

**The recurring failure mode this session was both agents editing the same
`.qmd`, which forks it.** Reconciling forks was expensive and happened more than
once. Keep a single source of truth per file: if CC is going to run/commit it,
the assistant hands over a finished copy and then stops touching it.

## Working conventions Larry holds firmly

* **Recommend, then implement.** Present options A/B/C; "proceed" triggers the
  work. Do not fold new decisions into a task silently.
* **Transparency.** Do not defer sub-decisions into a "deferred" list without
  saying so. Do not surface new issues *after* answering a question as though
  they were part of the answer. If something new turns up, say "new finding"
  explicitly. Larry called this out repeatedly; treat it as a hard rule.
* **Lose no file.** Prefer `git mv` to deletion. Larry keeps a manual backup, but
  the default is still move-not-delete.
* **`eff*` spelling** for `sg_focus`, never the internal `hr*`, because `eff*`
  generalises to GLM outcomes.
* **Verify against the code, not comments or filenames.** Both lie in this
  directory (see naming trap below). Where a comment and the code disagree,
  trust the code; a stale comment is a documentation fix — unless the code it
  describes has an unwanted effect.

---

## Fixed vocabulary (do not substitute descriptive forms in results/prose)

* **FS / DINA / GRF** — identifiers (engines).
* **FB** = full bootstrap, **MR** = multiplier resampling — inference procedures,
  *not* estimators. The identifier is the estimator; FB/MR correct its selection
  bias.
* **MR does not "re-select" a subgroup.** It de-biases against a fitted candidate
  family. Language implying MR identifies or re-selects a subgroup is wrong and
  Larry has flagged it specifically.
* **β(Ĥ)** conditional estimand; **θ†(H)** marginal causal HR in the true
  subgroup. Never conflate, and never conflate the log-HR (β) and HR scales.
* **κ = β(H)/β(Ĥ)** identification factor.

---

## State of the world (as of this handoff)

Completed and committed on `feature/mr-in-replicates` (local; **confirm these are
pushed — at last check the remote was several commits behind**):

* `analysis_gbsg_survival_multimethod.qmd` — vocabulary purge (Tier→FB, gate→MR),
  provenance tables reading `fit$args_call_all`, 12-finding bug review fully
  applied, numeric invariance proven via a **same-machine control render**.
* `sim_fs_maxcons_fb_mr_m1_h10_knoise0_n500_batch_1_100.qmd` — the **model
  template** for the simulation family. Carries: FB/MR vocabulary, `sg_focus`
  reordering fix, `mr_in_replicates = FALSE` pinned, record-only provenance in
  the payload `meta`, and now the **κ / β(H) treatment** in both table captions.
* Reorg Steps 0–2 done: manifest at `dev/sim-structure/rds_manifest.csv`; six
  κ-less combine documents archived to `legacy/`; **all 59 top-level `.rds`
  moved into `results/`** (basenames unchanged), with the template, the `t1_t2`
  documents, and the `dev/` consumers repathed.
* Reorg Step 3 (untrack renders) was **cancelled by Larry** — anything tracked
  stays tracked. Do not revisit untracking.

Deferred, still open (proposal Steps 4–7):

* Fold the 60 `t1_t2` documents → one parameterised template. **Held until Larry
  says the maxcons setup is settled.** Do not start it.
* **Step 6 — the `mr_coverage_sweep_*` fold. This is the next session's subject.**
* Resolve the `_sim_*` include-fragments (Step 7).

Larry's own next action, independent of the assistant: run the maxcons
production baseline (template currently at smoke values n_sims=50/nb_boots=20/
mr_draws=500 — set deliberately before running; move the old
`fs_maxcons_..._res_1_100.rds` aside first so it doesn't collide in the glob).

---

## The `mr_coverage_sweep` family — what is actually there

**9 documents, `mr_coverage_sweep_<cell>.qmd`**, cells:
`h075, h10, h10_knoise0, h10_knoise3, h10_knoise6, h15, h20, h25` (and one more).
Each has a paired `.html`. They are near-identical: an `h10_knoise0` vs `h15`
diff is **4 lines**. So this is the same 60→1 duplication as the `t1_t2` family,
one order smaller.

**They are already the better design.** Each is driven by a `run_tag` at the top
of the file:

```r
run_tag     <- "m1_h10_knoise0new_s1000"          # bump for a fresh isolated run
results_dir <- file.path("mr_sweep", run_tag)      # bundles + grid live here
grid_path   <- file.path(results_dir, sprintf("mr_coverage_grid_%s.rds", run_tag))
```

So results already live in a per-run subtree, `mr_sweep/<run_tag>/`, **not** at
top level — the sweep family is unaffected by the Step 2 `.rds` move. The
proposal calls Step 6 "parameter extraction with no logic change," and that
looks right: the `run_tag` machinery is the fold mechanism already; the cells
differ only in the DGM parameters the tag encodes.

**But they carry the same stale vocabulary the other files had, unfixed:**

* `gate` — 5 occurrences in `mr_coverage_sweep_h10_knoise0.qmd`, all real
  (prose and comments: "MR de-biased gate", "gate re-selection", "MR gate"), none
  in a protected word.
* `t1`/`t2` — 4 occurrences, including a **code** one:
  `.EST_KEYS <- c(oracle = "or", naive = "nv", FB = "t1", MR = "t2")` at ~L447,
  and `t1 <- bundles[[1]]$truth` at ~L629.
* `sg_focus <- "eff"` with a comment mentioning gate re-selection (L128) — the
  same `"eff"`→`"maxcons"` and re-selection-language issues already fixed in the
  maxcons template.

So the honest answer to "do they need revision" is **yes on two axes**, and the
next session should treat them separately:

1. **Vocabulary / correctness parity** with the maxcons template — the gate→MR,
   t1/t2→fb/mr, `eff`→`maxcons`, and MR-re-selection-language fixes. These are
   the same edits already validated on two other files; the risk is the `t2m`/
   `focus_tag` class of near-miss (see below), not the edits themselves.
2. **Structural fold** — 9 → 1 via a `params:` block, the same move the proposal
   defines for `t1_t2`. This is optional and larger; it can wait or be done in
   the same pass, Larry's call.

---

## Traps proven real this session — check for each

* **Word-boundary regex misses.** `\bt2\b` does not match `t2m`; `\b` does not
  match between `_` and `t` (so `_t1_t2_` survives). A residual check of
  `Tier|tier1` missed lowercase `tier`. **Search case-insensitively and include
  substrings**; protected words are `frontier`, `delegate*`, `gatekeeper`,
  `ungated`, `Negate`.
* **Definition order.** A symbol can be referenced before it is defined and still
  `parse()` cleanly — it only fails at execution. **Source the setup chunk in a
  clean session** after any reorder; parsing is not enough.
* **Filenames and `run_mode` lie.** In this directory, files named `..._combine_*`
  were `run_mode <- "batch"`; ranges in names (`combine_1_400`) did not match the
  rows produced (500). **The `run_mode` / `sim_id_start` / `n_sims` triple inside
  a file is the only reliable descriptor.** Never key logic on a filename.
* **Cross-machine numeric drift.** The pipeline is bit-reproducible *within* a
  machine but not across — worker count feeds the search batch size, changing
  summation order. `tolerance = 0` is not a cross-machine test. Any invariance
  check must be a **same-machine control render** of the unedited source. See
  `dev/replication-check/REPLICATION_FINDINGS.md` §2.2.
* **`.rds` provenance.** A batch `.rds` records what produced it only if the doc
  writes provenance into `meta`. The maxcons template now does (`sg_focus`,
  `forestsearch_version`, `n_workers`, host, …); **check whether the sweep does**,
  and consider adding the same record-only block if not. It is record-only —
  never add provenance fields to the combine agreement-check key vector, which
  would `stop()` on every pre-existing bundle.

---

## Suggested first moves for the next session

1. Read one sweep document end to end (`mr_coverage_sweep_h10_knoise0.qmd` is a
   fair representative) and confirm the vocabulary/config inventory above against
   the current file — it may have moved since this handoff.
2. Diff two cells to confirm the 4-line delta and identify exactly which lines
   carry the cell parameters — that set is the eventual `params:` block.
3. Decide with Larry: **parity-only pass first, or parity + fold in one.** Parity
   is the same validated edit set; the fold is new structural work with the
   `t2m`/order traps in play. Recommend, then implement.
4. Whatever the scope, verify by: parse, source the setup chunk clean, behavioural
   test any changed expression, and — if anything that could touch a number
   changes — a same-machine control render diffed at `tolerance = 0` on
   non-timing columns.

## Key paths

* Sweep docs: `quarto/simulations/gbsg_redux/mr_coverage_sweep_*.qmd`
* Sweep results: `quarto/simulations/gbsg_redux/mr_sweep/<run_tag>/`
* Model template: `quarto/simulations/gbsg_redux/sim_fs_maxcons_fb_mr_m1_h10_knoise0_n500_batch_1_100.qmd`
* Proposal + manifest: `dev/sim-structure/SIM_REORG_PROPOSAL.md`, `dev/sim-structure/rds_manifest.csv`
* `betaHhat_truth.R` is `source()`d from the working directory by ~70 docs — any
  layout that moves a `.qmd` into a subdirectory breaks that. Documents stay at
  top level; only `.rds` moved.
* Reference framework: `references.bib`.
