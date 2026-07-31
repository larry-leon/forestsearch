# CC BRIEF — pre-merge replication check (VERIFICATION ONLY)

```
claude "Read dev/replication-check/CC_BRIEF_replication.md and execute it."
```

Run from a worktree on `feature/mr-in-replicates`. **Change no package code, no
`.qmd`, and no test.** Write only under `dev/replication-check/`. Hash the
package sources before and after and report the verdict, reusing
`fs_hash_sources()` / `fs_guard_verify()` from `dev/efficiency-eval/R/00_guard.R`.

This is the gate before merging into `feature/glm-extension`.

---

## 1. Expectations differ by artifact --- read this first

Do not treat "replication" as one criterion. The `sg_focus` work deliberately
changed behaviour for one class and deliberately did not for the other.

| artifact | `sg_focus` | expected |
|---|---|---|
| GBSG multimethod / maxeff analyses | `effMaxSG` | **EXACT replication.** Early stopping was already disabled for `hrMaxSG` by the pre-existing reset at `forestsearch_main.R:1200`, so none of the `sg_focus` work touches this path. |
| `gbsg_redux` simulation cells | `eff` -> `hr` | **DIVERGENCE EXPECTED.** `"hr"` no longer early-stops, so the full candidate pool is now always evaluated where previously a truncated prefix often was. |

For the simulations, **a difference is the correct outcome and a match would be
suspicious.** The task there is to confirm the difference has the predicted
mechanism, not to make it go away. Under no circumstances adjust anything to
force agreement.

---

## 2. Part 0 --- the RNG question (do this first; everything depends on it)

`mr_in_replicates = FALSE` means MR no longer runs inside bootstrap replicates or
CV folds. **If MR consumed from the ambient RNG stream, removing it shifts every
subsequent draw** and the GBSG analyses cannot replicate exactly --- for a reason
unrelated to selection.

Establish, by reading `fs_mr_inference()` and then by execution:

1. Does MR seed its own multiplier draws (e.g. from `mr_inference_args$seed`),
   or does it draw from the ambient stream?
2. If it seeds itself: does it restore the prior RNG state on exit, or leave the
   stream advanced?
3. **Decisive test.** With a fixed outer seed, run a short bootstrap twice ---
   once with `mr_in_replicates = TRUE`, once `FALSE` --- and compare the
   *bootstrap* outputs (selected subgroups per replicate, `Ystar_mat`, the
   bias-corrected estimates). If they are identical, MR is RNG-neutral and Part 1
   can expect exact replication. If not, report the discrepancy: exact
   replication of any analysis using bootstrap or CV becomes impossible, and
   that is a finding the maintainer needs before merging.

**Stop and report if MR is not RNG-neutral.** Do not proceed to Part 1 on the
assumption that it is.

---

## 3. Part 1 --- the GBSG analyses (exact replication expected)

Target: `quarto/applications/gbsg/analysis_gbsg_cox_maxeff.qmd` and the
multimethod notebook the maintainer identified
(`analysis_gbsg_cox_multimethod_psi_v2_2new.qmd`, or the `_v3a` variant if that
is the current one --- confirm which before running, since several near-identical
copies exist).

Method: render or source under the current branch, and compare against the
**existing rendered `.html`**, which predates all of this work.

Compare, and report each separately:

* the **selected subgroup** for FS, DINA and GRF --- the primary criterion;
* subgroup size `n`;
* the naive and bias-corrected effect estimates and their intervals;
* the MR de-biased estimate and `mr_harm_confirmed`;
* the CV/LOO summary metrics.

Note the notebook sets `run_debias_gate <- TRUE` and `gate_draws <- 5000L` and
calls **both** `forestsearch_tenfold()` (line ~1595, `Kfolds = 10`) and
`forestsearch_Kfold()` (line ~1711, LOO). Those are the arms `mr_in_replicates`
now affects, so they are where an RNG-related difference would surface first.

**Also report wall-clock** against whatever the notebook's own timing table
recorded previously. Part A measured +221% for MR inside replicates at
`draws = 5000`, so a substantial reduction is expected and is the practical
benefit of the change.

If any reported number differs, **stop and report it with the mechanism** ---
do not adjust the notebook.

---

## 4. Part 2 --- the simulations (characterise, do not replicate)

The six `gbsg_redux` cells characterised in
`dev/sg-focus-work/PHASE0_FINDINGS.md` ran `sg_focus = "eff"` with
`stop_threshold = 0.95` inherited. Under the current branch `"hr"` no longer
early-stops.

**Do not re-run all 500 replicates per cell.** Take the first 100 of two cells
--- `h10_knoise0_n500` (highest branch-2 rate, 24.0%) and `h20_knoise3_n1500`
(branch 1 fired in 100% of selections) --- and compare against the stored
`*_combined_1_500.rds` bundles.

For each replicate report: whether `sg_def` matches the stored value (membership,
using the label normaliser from `dev/sg-focus-work/R/phase0_summarise.R`), and
where it differs, the stored subgroup, the new one, and both `Pcons` values.

**The predicted mechanism**, to be confirmed or refuted:

> Previously the loop halted at the first candidate reaching 0.95 and
> `sort_subgroups()` applied `(-Pcons, -hr, K)` to that **prefix**. Now the full
> pool is evaluated and the same key is applied to **all** qualifiers. So the new
> selection should equal the *old branch-2 answer* in every replicate --- and
> Phase 0 §6 found a median of 1 qualifier visible in the truncated prefix
> against 4--111 present in the full family, so divergence should be **common**,
> not rare.

Report the divergence rate per cell. **A rate near zero would refute the
mechanism and needs investigating.** Also report whether `h20_knoise3_n1500`
diverges less than `h10_knoise0_n500`, as the branch mix predicts.

**No configuration reproduces the old behaviour** --- the `stop_threshold` reset
for `"hr"` is unconditional. Do not attempt to construct one.

---

## 4b. Part 3 --- pin the defaults explicitly (ONLY after Parts 0-2 pass)

**Gated.** Do not start until Parts 0 and 1 have completed and the GBSG analyses
replicate exactly. Pinning values changes the input, so doing it first would
destroy the comparison.

The maintainer's aim: a reader of a `.qmd` should see what the analysis actually
used, without having to consult `?forestsearch`. Most of the trouble traced in
`dev/sg-focus-work/` came from a value that was never written down --- the
simulations inherited `stop_threshold = 0.95` against a `pconsistency.threshold`
of 0.90, and nothing in the notebook said so.

### Pin these, at their current effective values

`sg_focus`, `pconsistency.threshold`, `hr.threshold`, `hr.consistency`,
`consistency_method`, `maxk`, `n.min`, `d0.min`, `d1.min`,
`max_subgroups_search`, `mr_inference` --- and `mr_in_replicates` on every
`forestsearch_bootstrap_dofuture()`, `forestsearch_Kfold()` and
`forestsearch_tenfold()` call.

Each must be pinned to the value that run **currently resolves to**, so pinning
is a no-op. Verify that: re-run and confirm the results are unchanged from
Part 1. If any pinned value changes a result, the value was wrong --- stop and
report.

### Do NOT pin these

* **`stop_threshold`.** It is reset to `NULL` for `hr`, `maxSG`, `minSG`,
  `hrMaxSG` and `hrMinSG`, and the reset **warns when the value was supplied
  explicitly**. Pinning it would emit a warning on every run. Add a comment
  instead, e.g. `# stop_threshold: not applicable -- reset to NULL for this
  sg_focus; see ?forestsearch`.
* **`selection_rule` and `effect_neighborhood` on non-band foci.** They are
  consulted only by `hrMaxSG`/`hrMinSG` (equivalently `effMaxSG`/`effMinSG`).
  Pinning `selection_rule = "neighborhood"` under `sg_focus = "eff"` implies the
  band applies when it does not. Pin them **only** in the `effMaxSG` analyses,
  where they are operative, and comment their inertness elsewhere.

### Scope

The GBSG application notebooks and the `gbsg_redux` simulation cells. If several
near-identical variants exist, ask which are current rather than editing all of
them --- the maintainer has flagged this duplication before.

**Do not change any value.** This is annotation only: every pinned argument gets
the value already in force.

**Acceptance:** re-running one GBSG analysis and one simulation cell after
pinning reproduces the Part 1 / Part 2 results exactly, and emits no new
warnings.

## 5. Deliverable

`dev/replication-check/REPLICATION_FINDINGS.md`:

* guard verdict; `git status` confined to `dev/replication-check/`;
* **Part 0**: whether MR is RNG-neutral, with the evidence;
* **Part 1**: side-by-side table, old vs new, per reported quantity, with a
  clear exact / not-exact verdict and wall-clock comparison;
* **Part 2**: per-cell divergence rate, whether the mechanism is confirmed, and
  a handful of worked examples;
* **Part 3**: which arguments were pinned in which files, confirmation that
  pinning was a no-op, and the list of arguments deliberately left unpinned with
  the reason;
* a **merge recommendation**: is anything blocking, and what does the maintainer
  need to decide first.

---

## 6. Out of scope

* Fixing anything. Parts 0-2 are verification; Part 3 is annotation that must
  provably change no result.
* Re-running full 500-replicate cells.
* The `quarto/` prose sweep, F1, F3, Phase 6.
* Any judgement about which `sg_focus` the manuscript should use --- that is the
  maintainer's, and it is informed by Part 2 rather than settled by it.
