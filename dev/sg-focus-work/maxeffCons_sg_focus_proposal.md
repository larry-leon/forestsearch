# Proposal: `sg_focus = "maxeffCons"`

**Status:** design specification for maintainer review. Companion to
`dev/sg-focus-work/CC_BRIEF_sg_focus_implementation.md`, which carries the
implementation phases and acceptance criteria. This document carries the
rationale.

---

## 1. The rule

> the candidate with the **highest effect** among those meeting at least
> `pconsistency.threshold` consistency.

This is the rule the León Table-2 simulations were built to exercise. No
current `sg_focus` implements it as its selection contract.

| focus | final sort key | why it is not the rule |
|---|---|---|
| `maxeff` | `(-hr, K)` | **no consistency gate at all** -- see §4 |
| `hr` / `eff` | `(-Pcons, -hr, K)` | consistency-primary; effect is only a tiebreak |
| `maxSG` / `minSG` | `(-N, -Pcons, K)` / `(N, -Pcons, K)` | size-primary |
| `hrMaxSG` / `hrMinSG` | `(-in_band, -N, -Pcons, -hr, K)` | size-primary within an effect band |

---

## 2. How the simulations obtained it anyway

By accident, through candidate-level early stopping. Enumeration is
HR-descending (`sort_subgroups_preview()`, `setorder(-HR, K)`), and the loop
halts at the first candidate reaching `stop_threshold` -- short-circuiting
`sort_subgroups()` before its consistency-primary key can apply.

Two problems with relying on that:

1. `stop_threshold` defaults to **0.95** while `pconsistency.threshold`
   defaults to **0.90** (`forestsearch_main.R:886-887`, adjacent lines). A
   candidate in `[0.90, 0.95)` qualifies but does not halt the scan, so a
   *lower*-effect candidate can be returned.
2. If nothing reaches `stop_threshold`, the run falls through to
   `(-Pcons, -hr, K)` and returns the **max-consistency** candidate instead.

The effective rule is therefore a per-replicate hybrid, not a rule at all.
Phase 0 of the CC brief measures which branch fired, per simulation cell.

---

## 3. Naming: the current convention is inverted

The package supplies its own evidence. `forestsearch()` derives the debias
gate's re-selection from `sg_focus`, mapping **`"eff"` -> `"maxcons"`**. The
same rule carries two names, and only the gate's is honest.

| sort key | what it is | honest name | current name |
|---|---|---|---|
| `(-hr, K)`, ungated | effect argmax | `maxeff` | `maxeff` (correct) |
| `(-hr, K)`, gated | effect argmax, consistency-gated | `maxeffCons` | *absent* |
| `(-Pcons, -hr, K)` | consistency argmax | `maxcons` | `eff` / `hr` (misleading) |

The bare stem `eff` denotes the *consistency*-primary rule, while a compound
name would denote the *effect*-primary one -- the reverse of the compositional
reading any user applies.

**`maxeffCons`, not `effCons`.** It derives from `maxeff`, which genuinely is
effect-primary, so the relation reads correctly: ungated effect argmax ->
gated effect argmax. `effCons` would inherit from a stem that does not mean
what it says.

**Do not reuse `"eff"`.** Existing code passing `sg_focus = "eff"` -- including
the 500-replicate simulation runs -- would silently switch from
consistency-primary to effect-primary. No error, no warning, a different
subgroup. Silent semantic drift on a live name is worse than the present
inconsistency.

### Proposed vocabulary

| name | rule | status |
|---|---|---|
| `maxeff` | effect argmax, **no** consistency gate | unchanged |
| `maxeffCons` | effect argmax among consistency-qualifying candidates | **new** |
| `maxcons` | consistency argmax, effect as tiebreak | **new alias** for `hr` / `eff` |
| `maxSG` / `minSG` | size argmax / argmin, `Pcons` tiebreak | unchanged |
| `hrMaxSG` / `effMaxSG` | band-filtered size argmax | unchanged |
| `hrMinSG` / `effMinSG` | band-filtered size argmin | unchanged |

`maxeffCons` takes no `hr*` canonical form, matching `maxeff`; the `hr*`
canonicalisation exists only for the legacy Cox vocabulary.

**Recommendation: adopt `maxcons` in the same release.** It is purely
additive -- `eff` and `hr` keep working and keep their current meaning -- and
it makes the `sg_focus` vocabulary agree with the `reselection` vocabulary that
already uses the name. Deferring it means shipping `maxeffCons` into a
vocabulary that still reads backwards.

---

## 4. Specification

**Contract.** Among candidates with `Pcons >= pconsistency.threshold`, select
the largest effect; ties by `K`.

Sub-threshold candidates already return `NULL` from
`evaluate_subgroup_consistency()` and never reach the sort, so the gate is
implicit and the key is simply:

```r
if (sg_focus == "maxeffCons") {
  data.table::setorder(result_new, -hr, K)
  return(result_new)
}
```

* **Preview sort** -- `setorder(-HR, K)`, as for `"hr"` and `"maxeff"`.
* **Do not** route through the `maxeff` winner-only fast path
  (`subgroup_consistency_main.R` SECTION 8). Every candidate must be evaluated
  or the gate is meaningless.
* `selection_rule` stays restricted to `"neighborhood"` and inert;
  `.validate_selection_rule()` already enforces this outside
  `{hrMaxSG, hrMinSG}`.

### Why `maxeff` and `maxeffCons` genuinely differ

Stronger than "no filter on non-winners". In the fast path
(`subgroup_consistency_main.R:717-721`):

```r
res_m <- .maxeff_eval(m, skip = (m != 1L))
if (m == 1L && is.null(res_m)) res_m <- .maxeff_eval(m, skip = TRUE)
```

Non-winners are never evaluated. The winner *is* -- but if it returns `NULL`,
the code re-runs it with `skip = TRUE` and keeps it with `Pcons = NA`. `NULL`
is returned both when a candidate is inestimable **and** when it simply fails
the threshold, and the `is.null()` branch does not distinguish them. So a
`maxeff` winner that fails `pconsistency.threshold` is retained regardless.

The roxygen's "No consistency filter" for `maxeff` is therefore literally
accurate. The two foci share a sort key and differ **only** in whether the gate
binds -- a load-bearing distinction, and precisely what the name encodes.

### The property worth having

Enumeration (`-HR`) and selection key (`-hr`) **agree**. Therefore:

* the first qualifying candidate in evaluation order *is* the winner;
* early stopping is **provably valid** -- unlike `"hr"`, whose key is
  `Pcons`-primary and cannot be pre-sorted, since `Pcons` is the quantity being
  computed;
* additional evaluated candidates can only be lower-effect and cannot displace
  the incumbent, so **`parallel_args$batch_size` is provably irrelevant to the
  result**.

That last point matters: the observed batch-size sensitivity is not patched for
this focus, it is structurally absent.

---

## 5. Threshold coupling

Replace the numeric default with a lazily-coupled formal:

```r
pconsistency.threshold = 0.90,
stop_threshold         = pconsistency.threshold,
```

R's lazy evaluation makes this track automatically if the floor moves. A
hardcoded `0.90` re-breaks the moment a user sets
`pconsistency.threshold = 0.85`.

---

## 6. Guard membership

Lowering the `stop_threshold` default **amplifies** the existing `"hr"` defect:
the stop fires sooner, at a higher-effect / lower-`Pcons` candidate, deviating
further from `"hr"`'s documented rule. The default change must therefore land
together with the guard correction.

First-passer validity follows from whether enumeration order agrees with the
selection key's primary term:

| `sg_focus` | enumeration | key primary | valid? | action |
|---|---|---|---|---|
| `hr` | `-HR` | `-Pcons` | no | **add to the `forestsearch_main.R:1200` disable list** |
| `minSG` | `n` asc | `N` asc | yes | keep in batch guard |
| `maxSG` | `-n` desc | `-N` | yes | **add to batch guard** |
| `maxeffCons` | `-HR` | `-hr` | yes | **add to batch guard** |
| `hrMaxSG` / `hrMinSG` | band-filtered | `-in_band` | no | already disabled |

`subgroup_consistency_main.R:836-841` becomes
`c("minSG", "maxSG", "maxeffCons")`. The line-1200 comment already states the
principle -- *"early stopping is incompatible with compound sg_focus criteria
that require comparing HR and size across all candidates"* -- and `"hr"`
requires comparing `Pcons` across all candidates: the same structure, applied
to two of the three foci that need it.

**`"hr"` is retained**, not deprecated. Consistency-primary selection is a
coherent estimand; the defect was only that early stopping implemented a
different one. Its roxygen should state plainly that it is consistency-primary
with effect as tiebreak, since users read the name as the opposite.

---

## 7. Relationship to the existing simulations

The 500-replicate runs used `sg_focus = "eff"` with `stop_threshold = 0.95`
against `pconsistency.threshold = 0.90`. If Phase 0 shows branch 2 never fired
and the selected candidate was always HR-rank 1 among qualifiers, then
`sg_focus = "maxeffCons", stop_threshold = 0.95` reproduces them exactly and a
footnote suffices. Otherwise the affected cells need rerunning.

Indirect evidence from the `h10_knoise0_n500` bundle points toward
effect-maximising behaviour -- naive HR median 1.574 against a truth of 1.0,
and `corr(log naive HR, n_sel) = -0.429`, i.e. small extreme subgroups winning,
which max-`Pcons` selection would not favour. But that is inference from a
single cell, and the branch mix should differ across the grid: low-`h` /
small-`n` cells are where candidates land in `[0.90, 0.95)`.

**Manuscript code should pass `stop_threshold` explicitly** regardless of the
default, so published results never depend on a package default.

---

## 8. Release notes

A **behaviour-changing default**, not a tweak. `forestsearch` 0.1.0 is on CRAN;
an existing analysis re-run under 0.2.0 can return a different subgroup. It
warrants a prominent `NEWS.md` entry covering: `"hr"` no longer early-stops;
`stop_threshold` now follows `pconsistency.threshold`; new `maxeffCons`; new
`maxcons` alias.

---

## 9. Gate alignment (resolved -- and mandatory)

The Tier-2 gate is about the MR de-biased estimate, not about choosing the
reported subgroup. But the correction works by **replaying the selection step
inside each MR draw** (`.fs_dg_select()`, `fs_debias_gate.R:130-153`), because
resampling a *fixed* subgroup would ignore that the subgroup was itself
selected -- the source of the winner's curse. The replayed rule must therefore
match the analysis rule. This is a validity requirement, not a preference: a
mismatch de-biases against a selection event that never occurred.

An earlier draft of this document offered "the tiers deliberately diverge" as
an option. That is incoherent and is withdrawn.

**Resolution.** `.fs_dg_reselection_from_focus()`
(`fs_debias_gate_methods.R:134-145`) needs one added case:

```r
maxeffCons = "maxeff",
```

The gate's `"maxeff"` is `passers[which.max(beta[passers])]` -- argmax effect
**among passers**, already restricted to the qualifying set. That is precisely
the `maxeffCons` rule. It differs from `sg_focus = "maxeff"`, which is ungated.
No new gate mode is needed.

**Why this cannot be deferred.** The `switch()` ends in a bare `hr_rule`
fallback, `"maxcons"` under the consistency engine. Omitting the case makes the
gate silently replay consistency-argmax selection under an effect-argmax
search -- no error, no warning, invalid Tier-2 estimates throughout. A
regression test should assert that every accepted `sg_focus` maps to an
explicitly named rule and never to the fallback.

Remaining to verify: that `passers` under the consistency engine is the
`Pcons >= pconsistency.threshold` set. The mapping is correct only if so.
