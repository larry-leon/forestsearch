# Phase 0 — characterisation of the as-run `sg_focus = "eff"` selection behaviour

Executed per `dev/sg-focus-work/CC_BRIEF_phase0_characterization.md`.
Measurement only; nothing was implemented and no package file was modified.

**Bottom line.** The brief's account of the mechanism is correct in every
particular. The as-run behaviour is the described hybrid, and it **does**
diverge from the intended rule (*argmax effect among candidates with
`Pcons >= 0.90`*) — in **166 of 2642** selections across the six cells.
Exactly one cell, `h20_knoise3_n1500`, meets the brief's "footnote, no rerun"
condition. The other five do not.

---

## 1. Guard verdict

`fs_hash_sources()` / `fs_guard_verify()` reused verbatim from
`dev/efficiency-eval/R/00_guard.R` (sourced, not copied).

| | |
|---|---|
| Files hashed | **132** (installed package tree + git worktree `R/`) |
| Before snapshot | `dev/sg-focus-work/out/guard_before.rds`, taken after the install |
| Verdict | **`ok = TRUE` — 132 source files verified unchanged** |
| `git status` | changes confined to `dev/sg-focus-work/` (see §8) |

The package was installed from source with
`devtools::install(dependencies = FALSE, upgrade = FALSE)` before the guard
snapshot, per the brief. A 10-replicate pass run *before* and *after* that
install returned bit-identical instrumentation, so the install changed nothing
— the previously installed copy was already built from the current `R/`.

Environment: R 4.6.1, x86_64-pc-linux-gnu, forestsearch 0.2.0, survival 3.8.9,
data.table 1.18.4, future 1.75.0, reference BLAS.

---

## 2. V1–V4 verification

All four claims hold. V1–V3 were checked **twice**: by reading the source, and
by observing the values at runtime.

| # | Claim | Verdict | Evidence |
|---|-------|---------|----------|
| **V1** | `sg_focus` normalised before the line-1200 reset and before `subgroup.consistency()` | **TRUE** | `forestsearch_main.R:1081` calls `.normalize_sg_focus()`, defined at `forestsearch_helpers.R:902-912`, which maps `eff -> "hr"`. Runtime: a live `forestsearch(sg_focus = "eff", ...)` call returns `args_call_all$sg_focus == "hr"` **and** `grp.consistency$sg_focus == "hr"`. |
| **V2** | line-1200 reset does **not** list `"hr"`, so `stop_threshold` stays 0.95 | **TRUE** | `forestsearch_main.R:1200` gates on `sg_focus %in% c("hrMaxSG", "hrMinSG")` only; the default is `stop_threshold = 0.95` at `forestsearch_main.R:887`. Runtime: `grp.consistency$stop_threshold == 0.95`, and no reset warning is emitted. |
| **V3** | batch guard keys on `sg_focus %in% c("hr","minSG")` -> `batch_size_parallel <- 1L` | **TRUE** | `subgroup_consistency_main.R:837-841`. Runtime: across 44 early-stopped pilot replicates, `n_candidates_evaluated == early_stop_candidate` in **44/44** — including one that stopped at candidate 2, the case that discriminates batch size 1 from anything larger (a bigger batch would have evaluated the whole batch and reported a larger `n_evaluated`). |
| **V4** | `evaluate_subgroup_consistency()` returns `NULL` below `pconsistency.threshold` | **TRUE** | `subgroup_consistency_helpers.R:1476-1482`. The brief flagged that only the condition had been read; **the body was read for this report** and it is a bare `return(NULL)` (preceded by a `details`-gated `cat()`). Two earlier `is.na(p.consistency)` paths also return `NULL`. |

### 2a. One correction to the brief's branch-1 description

The brief says branch 1 "returns the first candidate in HR-descending order
clearing 0.95". That is not literally what happens, though the distinction does
not change the conclusion.

What happens: the loop halts at the first candidate reaching 0.95
(`subgroup_consistency_main.R:944-968`), and then `sort_subgroups()` — with key
`(-Pcons, -hr, K)` for `sg_focus = "hr"` — is applied to **all** qualifiers
found in the evaluated prefix, and row 1 is selected
(`subgroup_consistency_helpers.R:1002-1006`). So the winner is the
**max-`Pcons`** candidate of the prefix, which is usually but not always the
stopping candidate. Branch 1 is therefore branch 2 restricted to a prefix,
not a different rule. The measured divergence below already accounts for this.

---

## 3. Method, and why one pass per replicate suffices

Search stage only: `forestsearch()` with `debias_gate = FALSE` and the Tier-1
bootstrap simply never called. Everything else is verbatim from the cell's
`.qmd` (§4). Runner: `dev/sg-focus-work/R/phase0_run_cell.R`.

**No package instrumentation was needed.** `subgroup.consistency()` already
returns `early_stop_triggered`, `early_stop_candidate`,
`n_candidates_evaluated`, `n_candidates_total`, `n_passed`, and
`out_sg$result` — the post-`sort_subgroups()` table of every qualifier found,
carrying each candidate's global index `m`, `Pcons`, `hr`, `N`, `K`. All of it
surfaces on `fs$grp.consistency`.

**No second run was needed either.** `sort_subgroups_preview()` orders
candidates `setorder(-HR, K)` for `sg_focus = "hr"`
(`subgroup_consistency_helpers.R:632-635`), so ascending `m` **is** descending
effect. The intended rule's winner is therefore the smallest-`m` qualifier —
and it always lies inside the evaluated prefix, because every candidate before
the stop *was* evaluated and every candidate after it has lower effect. So:

* **branch 1** = `early_stop_triggered == TRUE`
* **branch 2** = `early_stop_triggered == FALSE` (whole family evaluated)
* **divergence** = selected candidate's `m` != smallest qualifying `m`

This is an argument, so it was **checked against execution** (§6).

---

## 4. Cells and configuration

The six `consistency` (`fs_*`) cells in `quarto/simulations/gbsg_redux/` that
have a pooled 500-replicate bundle. Their six batch `.qmd` files are identical
apart from three fields:

| cell | `target_hr_harm` | `n_sample` | `k_random_noise` |
|---|---|---|---|
| `h10_knoise0_n500`  | 1.0 | 500  | 0 |
| `h10_knoise0_n1000` | 1.0 | 1000 | 0 |
| `h15_knoise0_n500`  | 1.5 | 500  | 0 |
| `h15_knoise0_n1000` | 1.5 | 1000 | 0 |
| `h15_knoise3_n1500` | 1.5 | 1500 | 3 |
| `h20_knoise3_n1500` | 2.0 | 1500 | 3 |

Shared: `seed_base = 8316951`, `sg_focus = "eff"`,
`consistency_method = "resample"`, `pconsistency.threshold = 0.90`,
`hr.threshold = 0.90`, `hr.consistency = 0.80`, `maxk = 2`,
`d0.min = d1.min = 10`, `use_twostage = TRUE`, `use_lasso = use_grf = FALSE`,
`conf_force = c("meno == 0", "er <= 0", "pgr <= 0")`,
`conf.cont_jcuts = list(er = 10)`, DGM `setup_gbsg_dgm(model = "alt", ...)` at
`n_super = 1e5`, `analysis_time = 84`, `cens_adjust = log(1.5)`.
Per replicate: `simulate_from_dgm(seed = seed_base + sim_id)` and
`seedit = seed_base + sim_id`; noise covariates seeded at
`seed_base + sim_id + 1e6`. Full spec in `dev/sg-focus-work/R/phase0_cells.R`.

3000 replicates were run (6 cells × 500), sharded 100 at a time.

---

## 5. Seed fidelity — **PASSES**

Compared against the `results` frame of each `*_combined_1_500.rds`.

| cell | n | `sg_def` exact | `sg_def` normalised | `detected` | `n_sel` (membership) |
|---|---|---|---|---|---|
| `h10_knoise0_n500`  | 500 | 494 | 499 | 500 | 345 / 346 |
| `h10_knoise0_n1000` | 500 | 463 | 500 | 500 | 349 / 349 |
| `h15_knoise0_n500`  | 500 | 466 | 500 | 500 | 462 / 462 |
| `h15_knoise0_n1000` | 500 | 378 | 500 | 500 | 487 / 487 |
| `h15_knoise3_n1500` | 500 | 415 | 500 | 500 | 498 / 498 |
| `h20_knoise3_n1500` | 500 | 398 | 500 | 500 | 500 / 500 |

**Detection reproduces 3000/3000. Subgroup membership reproduces 2999/3000.**

Two things need explaining.

**(a) The "exact `sg_def`" column is not a fidelity measure.** The gap is
entirely a *label rendering* change: a cut on a forced binary is now printed
`{meno == 0}` where the June-2026 package printed `!{meno}` (and
`!{meno == 0}` vs `{meno}`). Same subjects — `n_sel` matches on every one of
these rows. The normalised column collapses the two spellings; it reaches 500
for five cells and 499 for the sixth. Downstream code that parses `sg_def`
should be aware of this, but it is not a reproduction failure.

**(b) One genuine, isolated discrepancy** — `h10_knoise0_n500`, `sim_id = 349`:

```
saved (Jun 2026) : !{er <= 117} & !{age <= 46}   n = 121
re-run  (HEAD)   :  {nodes <= 3} & !{age <= 46}  n = 202
```

A different subgroup, not a relabelling. It is a **branch-2** replicate with
**10 qualifiers** — a max-`Pcons` contest among ten near-tied candidates,
exactly where a small change in candidate enumeration or duplicate removal
flips the winner. The `.rds` batch files date from 2026-06-29/06-30, while `R/`
last changed 2026-07-21, so roughly three weeks of package drift separates the
two runs. 1 in 3000 (0.03%). I did not bisect which commit caused it.

The brief's stop-condition ("if it does not reproduce, stop and report") is not
triggered: reproduction is 99.97% at the membership level and 100% on detection.

---

## 6. Cross-check: does the evaluated prefix contain the global winner?

The divergence measurement rests on the §3 argument. To check it against
execution rather than reasoning, the first 100 replicates of every cell were
re-run with `stop_threshold = NULL`, which forces the **entire** candidate
family to be evaluated, and the smallest-`m` qualifier was compared.

| cell | compared | full eval confirmed | `first_qual_m` agrees | `first_qual_hr` agrees | qualifiers seen (as-run, median) | qualifiers in full family (median) | max |
|---|---|---|---|---|---|---|---|
| `h10_knoise0_n500`  | 63  | yes | **63**  | 63  | 1 | 6   | 72  |
| `h10_knoise0_n1000` | 70  | yes | **70**  | 70  | 1 | 4   | 82  |
| `h15_knoise0_n500`  | 90  | yes | **90**  | 90  | 1 | 13  | 140 |
| `h15_knoise0_n1000` | 95  | yes | **95**  | 95  | 1 | 25  | 168 |
| `h15_knoise3_n1500` | 99  | yes | **99**  | 99  | 1 | 34  | 195 |
| `h20_knoise3_n1500` | 100 | yes | **100** | 100 | 1 | 111 | 200 |

**517/517 agreement, every cell.** The argument holds: early stopping hides a
great many qualifiers (median 1 seen against 4–111 present) but never hides a
*higher-effect* one, so a single as-run pass measures the divergence exactly.

---

## 7. Per-cell results

### 7a. Branch mix and divergence

Rates are over replicates where a subgroup was actually selected
(`n_selected`); a replicate with no qualifying candidate exercises no
selection rule. `diverge` = the selected candidate was **not** HR-rank 1 among
qualifiers, i.e. a higher-effect qualifying candidate was passed over.

| cell | detected /500 | selections | branch-1 (early stop) | branch-2 (full eval) | **divergence** | div. within branch 1 | div. within branch 2 |
|---|---|---|---|---|---|---|---|
| `h10_knoise0_n500`  | 346 | 346 | 263 (**76.0%**) | 83 (**24.0%**) | 47 (**13.6%**) | 32 (12.2%) | 15 (18.1%) |
| `h10_knoise0_n1000` | 349 | 349 | 250 (**71.6%**) | 99 (**28.4%**) | 57 (**16.3%**) | 37 (14.8%) | 20 (20.2%) |
| `h15_knoise0_n500`  | 462 | 462 | 417 (**90.3%**) | 45 (**9.7%**)  | 34 (**7.4%**)  | 21 (5.0%)  | 13 (28.9%) |
| `h15_knoise0_n1000` | 487 | 487 | 465 (**95.5%**) | 22 (**4.5%**)  | 16 (**3.3%**)  | 9 (1.9%)   | 7 (31.8%)  |
| `h15_knoise3_n1500` | 498 | 498 | 489 (**98.2%**) | 9 (**1.8%**)   | 12 (**2.4%**)  | 9 (1.8%)   | 3 (33.3%)  |
| `h20_knoise3_n1500` | 500 | 500 | 500 (**100%**)  | 0 (**0%**)     | **0 (0%)**     | 0 (0%)     | —          |

The brief's expectation about the grid is confirmed: branch 2 fires in ~24–28%
of low-`h`/small-`n` selections and never at `h20`/`n1500`. But note that
**branch 1 diverges too** — 108 of the 166 divergences are branch-1 events,
where candidate 1 qualified in `[0.90, 0.95)`, failed to trigger the stop, and
a later candidate with higher `Pcons` but *lower* effect won the
`(-Pcons, -hr, K)` sort. Divergence is not a branch-2 phenomenon.

### 7b. Distribution of the selected candidate's `Pcons`

| cell | n | min | p10 | p25 | median | p75 | p90 | max | mean | `>= 0.95` | in `[0.90, 0.95)` |
|---|---|---|---|---|---|---|---|---|---|---|---|
| `h10_knoise0_n500`  | 346 | 0.90 | 0.92 | 0.95 | 0.97 | 0.99 | 1.00 | 1.00 | 0.9636 | 76.0% | **24.0%** |
| `h10_knoise0_n1000` | 349 | 0.90 | 0.92 | 0.94 | 0.97 | 0.99 | 1.00 | 1.00 | 0.9605 | 71.6% | **28.4%** |
| `h15_knoise0_n500`  | 462 | 0.90 | 0.95 | 0.97 | 0.99 | 1.00 | 1.00 | 1.00 | 0.9787 | 90.3% | **9.7%**  |
| `h15_knoise0_n1000` | 487 | 0.90 | 0.96 | 0.99 | 1.00 | 1.00 | 1.00 | 1.00 | 0.9897 | 95.5% | **4.5%**  |
| `h15_knoise3_n1500` | 498 | 0.91 | 0.99 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 0.9953 | 98.2% | **1.8%**  |
| `h20_knoise3_n1500` | 500 | 0.98 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 0.9998 | 100%  | **0%**    |

(`>= 0.95` and the branch-1 rate coincide by construction: a selection with
`Pcons >= 0.95` is exactly one where the stop fired.)

### 7c. How much effect the as-run rule gave up

Among divergent replicates only, comparing the selected candidate's `hr` to the
passed-over HR-rank-1 qualifier's:

| cell | n div. | median HR gap | median gap (%) | max gap (%) | median HR-rank of selected | max rank | median `Pcons` of passed-over |
|---|---|---|---|---|---|---|---|
| `h10_knoise0_n500`  | 47 | 0.082 | 5.8% | 32.1% | 2 | 9 | 0.930 |
| `h10_knoise0_n1000` | 57 | 0.067 | 5.8% | 22.2% | 2 | 6 | 0.920 |
| `h15_knoise0_n500`  | 34 | 0.083 | 5.9% | 31.4% | 2 | 5 | 0.920 |
| `h15_knoise0_n1000` | 16 | 0.050 | 3.9% | 25.5% | 2 | 4 | 0.915 |
| `h15_knoise3_n1500` | 12 | 0.016 | 1.4% | 10.3% | 2 | 4 | 0.925 |
| `h20_knoise3_n1500` | 0  | —     | —    | —     | — | — | —     |

Typically the runner-up is taken (median rank 2) at a ~5% effect cost, but the
tail is not negligible: the selected candidate ranked as low as **9th** among
qualifiers, and single replicates gave up **>30%** of the hazard ratio.

---

## 8. Recommendation per cell

Applying the brief's rule — *branch 2 never fired **and** the selected candidate
was always HR-rank 1 among qualifiers* -> footnote; otherwise report and let the
maintainer decide.

| cell | branch 2 fired? | always HR-rank 1? | **verdict** |
|---|---|---|---|
| `h20_knoise3_n1500` | **no** (0/500) | **yes** (0 divergences) | **FOOTNOTE — no rerun.** `sg_focus = "maxeffCons", stop_threshold = 0.95` reproduces these runs exactly. |
| `h15_knoise3_n1500` | yes (1.8%) | no — 2.4% diverge | **maintainer's call.** Divergence is small and the effect cost is small (median 1.4%). A footnote may well be defensible. |
| `h15_knoise0_n1000` | yes (4.5%) | no — 3.3% diverge | **maintainer's call.** |
| `h15_knoise0_n500`  | yes (9.7%) | no — 7.4% diverge | **maintainer's call.** |
| `h10_knoise0_n1000` | yes (28.4%) | no — **16.3%** diverge | **maintainer's call — this is the weakest case for a footnote.** |
| `h10_knoise0_n500`  | yes (24.0%) | no — **13.6%** diverge | **maintainer's call — likewise.** Also carries the single genuine reproduction discrepancy (§5b). |

The gradient is monotone in signal strength and sample size, and the two
low-`h` cells — which the brief identifies as carrying the type-I-error and FDR
claims — are precisely the ones where the as-run behaviour departs most from
the rule it is described as implementing. In `h10_knoise0_n1000`, roughly one
detection in six selected a subgroup that was not the effect maximiser among
qualifying candidates.

---

## 9. What could not be determined

1. **Downstream impact.** This is a search-stage characterisation. It
   establishes *which subgroup* was selected, not what the divergence does to
   the Tier-1/Tier-2 estimates, coverage, type-I error, or FDR reported in the
   manuscript. Answering that requires re-running the affected replicates
   through the full pipeline (bootstrap included) under both rules, which was
   out of scope. Note the divergence is not obviously conservative in either
   direction: the as-run rule systematically prefers *higher consistency, lower
   effect*, which will tend to select more stable but less extreme subgroups.

2. **Global maximum `Pcons` under the as-run configuration.** Under branch 1
   the candidates after the stop are never evaluated, so the maximum `Pcons`
   over the whole family is unobservable from an as-run pass. §6 measures it
   on a 100-replicate subsample per cell (median 4–111 qualifiers exist where
   the as-run pass sees a median of 1). This does not affect the divergence
   measure, which depends only on the prefix.

3. **Cause of the `sim_id = 349` discrepancy** (§5b). Not bisected. It is a
   ten-way near-tie in a branch-2 replicate under ~3 weeks of package drift;
   identifying the responsible commit would need a bisect over
   `2026-06-30..2026-07-21`.

4. **The cross-check covers replicates 1–100 per cell**, not all 500. It is a
   verification of a structural argument, and it agreed 517/517, so extending
   it was not judged necessary.

5. **DINA and GRF cells** (`dina_*`, `grf_*` bundles in the same directory)
   were not characterised. They use different candidate-generation engines and
   the brief scopes the question to the `consistency` runs.

---

## 10. Artefacts

Everything written by this task lives under `dev/sg-focus-work/`:

```
PHASE0_FINDINGS.md              this report
R/phase0_cells.R                cell + shared-knob definitions (read from the .qmd files)
R/phase0_guard.R                thin wrapper over dev/efficiency-eval/R/00_guard.R
R/phase0_run_cell.R             search-stage runner + selection-rule instrumentation
R/phase0_summarise.R            seed fidelity, branch mix, divergence, Pcons distribution
R/phase0_crosscheck.R           prefix-vs-full-family verification
out/guard_preinstall.rds        hash snapshot taken before devtools::install()
out/guard_before.rds            hash snapshot used for the guard verdict
out/run_<cell>_<a>_<b>.rds      30 shards, 3000 replicates, as-run configuration
out/xchk_<cell>_<a>_<b>.rds     12 shards, 600 replicates, stop_threshold = NULL
out/xchk_probe_<cell>.rds       2 cost-probe runs (10 reps each), same stop_threshold = NULL
                                config; phase0_crosscheck.R's glob reads all 14 files, but
                                these repeat sim_id 1-10 of the *_1_50 shards and are dropped
                                by its dedup -- verified a no-op (600 replicates and
                                bit-identical per-cell frames with or without them)
out/pilot_*.rds                 10-replicate seed-fidelity pilots
out/summary_full.rds            pooled tables + all 166 divergent replicates
out/summary_crosscheck.rds      cross-check table
logs/                           per-shard stdout
```

To reproduce a cell:

```sh
Rscript dev/sg-focus-work/R/phase0_run_cell.R h10_knoise0_n500 1 500 4 out.rds
Rscript dev/sg-focus-work/R/phase0_summarise.R "dev/sg-focus-work/out/run_*.rds" full
```

Worker-count invariance was verified (8 vs 4 workers -> identical
instrumentation), consistent with the per-candidate RNG streams that
`subgroup_consistency_main.R:894-904` generates from the global candidate index.
