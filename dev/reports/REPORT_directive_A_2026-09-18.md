# REPORT — Directive A: c2 > c1 fails loudly; a silent c2 derives 0.80 · c1

**Task:** `dev/tasks/TASK_directive_A_2026-09-18.md`, committed alone at `ed34ce13`.
**Date:** 2026-09-18. **Branch:** `feature/glm-extension`. **Machine:** pop-os. **Evidence base
re-verified from source:** the audit (§3.5, §3.6, finding 11) and the sync report — pointers hold, line
numbers had drifted by the sync's own insertion. **Compute:** two micro-fits inside the acceptance tests.
**Tests:** this task's one acceptance file, nothing else. No `R CMD check`, no full suite. Not pushed.

| SHA | What |
|---|---|
| `ed34ce13` | the task document, as received |
| `5d2db58f` | the extended probe + the 175-cell baseline at the pre-edit SHA |
| `822cd69a` | the validation (`stop()` on c2 > c1, and on disagreeing spellings) |
| `0037ef6d` | the derivation (silent c2 → 0.80 · c1), announced and synced |
| `70aab091` | roxygen: the pair rule, `document()` |
| `4e1b6dd4` | acceptance tests |
| `787d138e` | the rider: `"OR"` first in `make_effect_estimator()`'s binary choices |
| this one | NEWS + this report |

## 1. Gate 1 — the prerequisite (Step 1)

Every site located and matching the sync report's description exactly; only the line numbers had moved, and
only by the sync's own 36 inserted lines. **Line numbers below are at the pre-edit SHA `a389ebf8`**, which
is where the gate was taken; §3 onward cites current HEAD. `forestsearch_main.R`: `effect_measure` default
`:1446-1452`; `user_set_*` detection + alias merge `:1462-1465`; `args_call_all <- mget(...)` `:1512-1513`,
after the merge; GLM resolution `:1896-1990`; **SECTION 2B-ii present** `:2121`;
`.sync_args_call_all()` `:21-30`. `subgroup_method <- match.arg()` at `:1424` and `outcome_type` at
`:1421` — both resolved **before** the merge, so the method and the estimand are known where the pair rule
must run. "c2 not supplied" is detectable there for both spellings, as the sync's analysis said:
`is.null()` on the new names, `missing()` on the legacy ones, evaluated in `forestsearch()`'s own frame.
**Gate 1: PASS.**

## 2. Baseline (Step 2)

`tests/testthat/helper-threshold-sync.R` extended — one copy, still source-driven: it lifts the new
statements out of `body(forestsearch)` under a `validate` switch on the same footing as the existing `sync`
switch, so the probe cannot drift from the source and the pre-edit tree remains reproducible in-session.
**175 cells** = the sync's 63, unchanged in identity and reproducing
`postsync_threshold_sync_2026-09-18.csv` column for column, plus four Directive A families: c1 alone at a
value whose 0.80 multiple is not the old default; c2 alone against the default c1; c2 > c1 under each of
`consistency` / `dina` / `grf`; and both spellings of one quantity disagreeing.
`dev/reports/baseline_directive_A_2026-09-18.csv`, committed at `5d2db58f` before any source edit: **no
cell errors, no cell announces anything, and c2 = 1.00 against c1 = 0.90 resolves silently.**

## 3. The change (Step 3)

`.fs_resolve_threshold_pair()` (`forestsearch_main.R:90`, beside `.sync_args_call_all()`), called once from
the new **SECTION 1A3** (`:1637-1665`), immediately after the alias merge and ahead of every fit. In order:
disagreeing spellings → `stop()`; c1 supplied and c2 not → `c2 <- 0.80 * c1`; `c2 > c1` → `stop()`;
announce the derivation. Deriving *before* the check is what makes the check the single guarantor of the
invariant on both routes, and holding the announcement until after it means a rejected call never reports a
threshold it will not use.

- **Scope gate**, first statement in the helper: `subgroup_method == "consistency"` **and**
  (survival **or** binary `"OR"`). Everything else returns `c2` untouched — no error, no derivation.
- **The derivation is on the ratio scale.** `0.80 * c1` before `log()` at `:2179-2180`, i.e. an additive
  `log(0.80)` shift, never `0.80 * log(c1)`. 1.25 → 1.00 (exactly, in doubles), 1.00 → 0.80, 0.90 → 0.72.
  The default pair (1.25, 1.0) is the rule's **fixed point**.
- **The sync carries it.** SECTION 2B-ii already writes the resolved natural-scale c2 into
  `consistency.threshold`; on a ratio estimand that is `hr.consistency`, now the derived value. A replay
  supplies that spelling, `is.null()` is FALSE, `user_set_consistency` is TRUE, and the derivation branch is
  not re-entered — so the message cannot fire per replicate. Asserted, not argued (§5).
- **`isTRUE()` on the comparison**, not a bare `if`: an `NA` or non-scalar threshold reaches the branches
  that already handle it rather than failing here with "missing value where TRUE/FALSE needed".

**Deliberate choice, recorded:** the announcement is gated on `!quiet`, matching every other resolution
announcement in `forestsearch()` (the adaptive `n.min` message at `:1686-1689`). `quiet = TRUE` suppresses the
message; it does not suppress the derivation.

## 4. Step 4 — the deletion gate: **DO NOT DELETE**

The audit's finding 11 says the entry condition *binds* only when c2 > c1, and that is right — with
c2 ≤ c1 enforced, screening's strict `hr > c1` (`subgroup_search.R:651`, `:692`) makes `hr > c2` automatic.
But *binding* and *reachable* are different, and the skip branch remains reachable by legal calls on two
independent routes:

1. **An empty candidate family.** `format_search_results()` returns `out.found = NULL` when nothing passed
   (`subgroup_search.R:997-1009`), so the guard at `forestsearch_main.R:3414-3419` fails and
   `has_subgroups` stays FALSE whatever the thresholds are. This is the no-subgroup contract every
   downstream consumer depends on.
2. **`sg_focus = "maxeff"` on the consistency engine.** `.fs_admission_applies("maxeff", "consistency")`
   returns `effect = FALSE` (`forestsearch_helpers.R:2308-2310`), so `disable_effect_floor = TRUE` at
   `forestsearch_main.R:3345` and candidate effects are **not** bounded below by c1. `any(hr_values > c2)`
   can then be FALSE with c2 ≤ c1, and the condition still binds.

`dina` and `grf` do not reach it at all — they `return()` from their own sections (`:2685`, `:2870`), well
before the stage boundary — so the task's phrasing of the gate is satisfied a fortiori. Nothing deleted;
both source facts are pinned by a test.

## 5. Gates (Step 6) and tests (Step 7)

`dev/reports/postA_directive_A_2026-09-18.csv`, 175 cells against the committed baseline.

| Gate | Result |
|---|---|
| **6a** every cell matches the Step 2 expectation table | **PASS** — 36 cells move, all in scope and on the consistency path: 18 derive (6 of them at the fixed point, where only the announcement is new) and 18 error. The other 139, including all 28 dina/grf cells and every RD / IRD / MD / IRR cell, are byte-identical. |
| **6b** parent-resolution invariance off the derivation/error cells | **PASS** — 139 cells unchanged; the 6 fixed-point cells are unchanged in every resolution column too |
| **6c** replicate == parent everywhere, derived cells included | **PASS** — 0 violations in 175; `n_rep_messages` is 0 in every cell |

`tests/testthat/test-threshold-pair-directive-a.R`: **86 pass, 0 fail, 0 skip, 4.2 s** including
`load_all()`, against a 3-minute abort. The two micro-fits are inside that: a binary-OR fit (N = 150,
20 splits) with `hr.threshold = 0.90` and no consistency spelling resolves
`threshold_config$consistency_natural = 0.72` and `$consistency = log(0.72)`, announces exactly once, and
carries 0.72 into `args_call_all`; and a c2 > c1 call on a 3-row frame raises the pair error, while the
same call with a legal pair gets past it and fails on the data instead. `test-threshold-sync.R` (16) and
`test-binary-default-or-entry-points.R` (25) re-run green, since the helper they share was extended.

## 6. Findings (no tasks attached)

**A1 — the package's own MRCT default already pairs the two.** `mrct_region_sims()`
(`mrct_simulation.R:175-176`) defaults `hr.threshold = 0.90`, `hr.consistency = 0.80` and passes both
explicitly at `:424-425`. It is unaffected (0.80 ≤ 0.90, both supplied), and it is the one place in the
package where the pair was already chosen rather than inherited.

**A2 — no committed caller is broken, and the exposure is smaller than the rule's reach.** Every in-package
`do.call(forestsearch, …)` site was checked: the bootstrap and CV replays supply both spellings (post-sync)
and inherit a validated pair; `fpr_calibration()` supplies both and has its own c2 ≤ c1 check;
`run_simulation_analysis.R:55-56` is (1.25, 1.0). `tests/testthat/helper-synthetic-dgm.R`'s `.fs_args_for()`
supplies (1.25, 1.00) for every outcome type, so the synthetic suite neither errors nor derives. A
file-level scan of tracked `.R` / `.qmd` / `.Rmd` for a c1 spelling without a c2 spelling returned only
`dina` / `grf` drivers, MD campaigns, and non-`forestsearch()` hits. The scan is file-level, not
call-level; it is weaker than the sync task's parsed inventory and is reported as such.

**A3 — a committed GRF campaign runs exactly the degenerate pair, deliberately untouched.**
`quarto/simulations/gbsg_020/summary_grfmr.qmd:26` documents `hr.threshold = 0.90` on the GRF re-selection
path against the default c2 = 1.0. Under `subgroup_method = "grf"` c2 is inert, so Directive A leaves it
alone by construction — which is the scope disposition working, not an oversight. It is also the clearest
illustration of why Directive C's dina/grf display annotation matters: nothing in that run's output says
c2 was never consulted.

**A4 — the derivation makes an explicit default-valued c1 newly chatty.** A driver that passes
`hr.threshold = 1.25` and no c2 resolves exactly what it resolved before (the fixed point) but now emits
one message per fit. Worth knowing when diffing run logs across this change.

**A5 — deviation from the task's Step 5.** The roxygen landed in its own commit (`70aab091`) rather than
inside the two behaviour commits, because the behaviour was already committed when the documentation was
written. The content is as specified; only the commit boundary differs.

## Out of scope, untouched

Directive C (DINA guard, frontier warnings, display defaults, dina/grf display annotation). RD / IRD / MD /
IRR behaviour. Floors. `fs-glms-interpretable`. `fpr_calibration()`'s own c2 ≤ c1 check.
