---
bibliography: []
---

# CC prompt

Paste everything below the rule into Claude Code, from the repo root, on
`mr-in-replicates`.

---

Read `dev/identifier-alignment/CC_HANDOFF_identifier_alignment.md` in full
before doing anything. It is the specification for this task; everything below
is scheduling. Its "Directory layout" section is verified against the actual
tree — note that the subdirectory is `sim_analyses` with an underscore.

Before Phase C you will need `betaHhat_truth.R`, which the three simulation
documents `source()` with a bare relative path. It is not in
`dev/identifier-alignment/`; find it in the repo and either copy it beside the
documents or repoint the call.

Work on branch `mr-in-replicates`. Do not merge to `feature/glm-extension`
first — two of the enforcement sites are `mr_in_replicates`, which only exists
on this branch. Commit after each phase separately, so Phase A can be reviewed
and reverted independently of the analyses.

Two standing rules for the whole task.

**Nothing writes into the reference files.** Three documents write filenames
that collide with the bundles they are compared against — see B0 and C0. Copy
the references somewhere safe, or render elsewhere, before running anything.
Treat `sim_analyses/` as read-only.

**Comparisons are not expected to match exactly.** Report what matched, what
did not, by how much, and the likely cause. Do not tune settings to force
agreement. Do not assert pass/fail thresholds. Where a difference follows from
a decision recorded in the handoff, report it as the expected consequence of
that decision.

## Phase A — package changes

Implement A1 through A7 as written. The six default changes carry exact line
numbers; verify each line still says what the handoff claims before editing,
since the file may have moved. Three of the six are `match.arg()` vector
re-orderings that change the default while leaving the permitted set intact —
do not delete any option.

Write `.validate_mr_configuration()` once and call it from all three code
enforcement sites in A3, including `forestsearch_Kfold()`. Errors only: no
warnings, no silent coercion.

Read A3b first. The per-engine re-selection mapping already exists in
`.fs_mr_reselection_from_focus()` — do not duplicate it.

Run `devtools::document()` and `R CMD check`. Do not create testthat files.

**Stop here and report** the diff summary and the check result before Phase B.

## Phase B — applications

Read B0 and B1 before touching either document.

Revise, render, and compare `analysis_gbsg_survival_multimethod.qmd` against
`gbsg_table_payload.rds`, and `analysis_actg175_binary_multimethod.qmd` against
`actg175_table_payload.rds`.

ACTG175's main fit sets `use_lasso = TRUE` and `use_grf = TRUE` with MR on,
which the Phase A validator now rejects. Set both `FALSE`, matching GBSG.
Expect the FS subgroup to change — if the label differs from the payload,
report old and new side by side as a finding. Do not try to recover the old
label. GBSG needs no configuration change.

## Phase C — simulations

Read C0 before rendering.

Install the package with `devtools::install()`, not `load_all()` — the doFuture
workers are separate R processes and only see the installed package.

Render three documents, one per engine: `sim_fs_maxcons_*`,
`sim_dina_maxcons_*`, `sim_grf_maxcons_*`. Each has `subgroup_method` pinned,
so do not edit between runs. Their configuration is already correct —
`sg_focus = "maxcons"`, `sim_id_start = 1L`, `n_sims = 20L`,
`run_mode = "batch"`. Do not regenerate them from the older
`sim_fs_maxeffCons_*` file, which sets a different focus.

They contain no validation code by design. Add the pre-loop check per C3 to
each, reusing the Phase A validator.

**Twenty replicates each. Do not run 1–500.** Compare each against its
reference per the C1 table, extracting the matching rows with
`subset(x$results, sim_id <= 20)`.

## Reporting

One summary at the end covering all five runs: what matched, what did not, the
likely cause of each difference, and anything that could not be run.

**Ask rather than guessing** if a file's contents contradict the handoff, if a
comparison is off by more than resampling noise would explain, if the validator
rejects a configuration you expected to pass, or if the reference bundles are not
where the handoff's verified layout says they are.
