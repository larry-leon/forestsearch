# CC BRIEF — selection subsystem evaluation (`sg_focus` x `selection_rule`)

**Run from the package root:**

```
claude "Read dev/selection-eval/CC_BRIEF.md and execute it."
```

---

## 1. Objective

Establish **empirically** how `sg_focus` and `selection_rule` interact, and
whether the Pareto machinery behaves as documented. This is a
**measurement and verification** task. Do not modify `R/`.

The motivating question: `selection_rule = "pareto"` / `"both"` filter to the
Pareto-non-dominated set with **both** effect and $N$ maximized, but
`hrMinSG` **minimizes** $N$. Source reading suggests these are in tension.
That needs confirming or refuting against real data.

---

## 2. Hard constraints

- Do not modify `R/`, `NAMESPACE`, `DESCRIPTION`, `man/`, `tests/`,
  `vignettes/`, `quarto/`.
- Work only under `dev/selection-eval/`.
- Install with `devtools::install(dependencies = FALSE)`, not `load_all()`.
- Hash `R/` before and after; report any change as a contract violation.

---

## 3. Phase 1 — validate the premise (do this first, it is cheap)

Render `pareto_frontier_explainer.qmd` with `forestsearch` installed. Its
validation chunk tests the exact assumption everything else rests on:

```r
identical(forestsearch:::pareto_dominated_flags(dt_hn), cands$dominated)
```

**If that FAILS, stop and report.** The document's local `dominated_xy()`
mirror would then not match `.pareto_dominated_xy()`, and the whole analysis
below is built on a false premise.

Also record whether the other three `.ok()` checks pass.

---

## 4. Phase 2 — the legal combination matrix

`.validate_selection_rule()` (`subgroup_consistency_helpers.R:712-737`) permits
`selection_rule != "neighborhood"` **only** for `sg_focus` in
`{"hrMaxSG", "hrMinSG"}`, and forbids a non-default `effect_neighborhood`
under `"pareto"`.

Enumerate every `sg_focus` value (`"hr"`, `"eff"`, `"maxeff"`, `"maxSG"`,
`"minSG"`, `"hrMaxSG"`, `"effMaxSG"`, `"hrMinSG"`, `"effMinSG"`) crossed with
every `selection_rule` (`"neighborhood"`, `"pareto"`, `"both"`), and record for
each cell: **accepted / errored**, and if accepted, the selected subgroup on
GBSG.

Report the matrix. Confirm the alias pairs (`hr`/`eff`, `hrMaxSG`/`effMaxSG`,
`hrMinSG`/`effMinSG`) give **identical** selections.

---

## 5. Phase 3 — the `hrMinSG` question (the headline)

### 3a. Reproduce the analytic claim on synthetic data

Using `forestsearch::sort_subgroups()` directly (it is exported), evaluate:

```r
dt <- data.table(id = c("P", "Q"), Pcons = c(0.90, 0.90),
                 hr = c(2.30, 2.45), N = c(50L, 100L), K = c(1L, 1L))
```

Analytic prediction: `hrMinSG` + `"neighborhood"` selects **P**;
`hrMinSG` + `"pareto"` selects **Q**. Confirm or refute.

Then test the general claim: **the strict-minimum-$N$ candidate survives
`"pareto"` only if it also has the strictly largest effect.** Generate random
candidate sets, compute `pareto_dominated_flags()`, and check whether the
min-$N$ row is ever non-dominated while not being the max-`hr` row. Report a
counterexample if one exists.

### 3b. Quantify the degeneracy

Over many random candidate sets, record how often
`hrMinSG` + `"pareto"` selects the same row as `maxeff` (pure effect argmax).
If that rate is near 1, `"pareto"` effectively nullifies the `MinSG` half of
the focus. Report the rate, not an interpretation.

### 3c. Check against the package's own documentation

`compute_pareto_frontier()`'s roxygen
(`subgroup_consistency_helpers.R:787-791`) states the frontier is a
*"post-hoc reporting artifact, not a selection criterion"* and that
`"hrMinSG"` *"may select an N-dominated point by design (preferring small
subgroups)."*

But under `"pareto"` / `"both"`, `in_band` excludes dominated candidates, so
`hrMinSG` **cannot** select an N-dominated point. Determine whether these two
statements are actually in conflict, and report which describes observed
behaviour.

---

## 6. Phase 4 — the preview sort's extra key

`sort_subgroups_preview()` (`subgroup_consistency_helpers.R:677-683`) applies a
**five-key** sort for the effect-anchored foci:

```r
on_frontier <- as.integer(!.pareto_dominated_xy(hr_vec, n_vec))
order(-in_band, -on_frontier, -n_vec, -hr_vec, K_vec)   # hrMaxSG
```

`on_frontier` is computed **unconditionally** — it does not depend on
`selection_rule`. The inline comment (lines 665-676) says this was added so the
full frontier is guaranteed inclusion in the top-`stop_Kgroups` candidates
passed to consistency, citing an issue documented in NEWS.

Establish:

1. Find and quote the relevant `NEWS.md` entry.
2. Does preview truncation actually bind on GBSG? Compare the number of
   candidates surviving the search against `max_subgroups_search`. If
   truncation never binds, this key is inert in practice and the whole concern
   is theoretical.
3. Where truncation does bind, does removing frontier-preservation change which
   candidates reach consistency? **Do not patch the package to test this** —
   call `forestsearch:::sort_subgroups_preview()` directly on a candidate
   table and compare against a locally-defined four-key variant.
4. The comment claims *"selection semantics are preserved"* because the final
   sort is unchanged. Assess precisely: the **rule** is unchanged, but the
   **candidate set** reaching it is not. State whether the selected subgroup
   can differ.

---

## 7. Phase 5 — documentation alignment

Compare `pareto_frontier_explainer.qmd` against observed behaviour and report
every discrepancy. Known candidates to check:

- The doc's §"Where the rules genuinely diverge" claims the top pick differs
  *only* on a max-band-$N$ tie. Verify whether that holds for `hrMaxSG` (source
  reading suggests yes) and for `hrMinSG` (suggests no).
- The doc does not mention the preview sort's `-on_frontier` key at all.
- The doc says *"Only `in_band` depends on `selection_rule` — the sort itself
  is fixed."* True for `sort_subgroups()`; check the preview.
- §"Practical guidance" directs users to `compute_pareto_frontier()`, which is
  `@keywords internal` and not exported. Confirm it needs `:::`.
- The doc never states that for non-effect-anchored foci, `selection_rule`
  is (a) restricted to `"neighborhood"` by the validator and (b) inert, since
  `sort_subgroups()` returns before `.compute_inclusion_band()` is reached.
  Confirm (b) by checking that `effect_neighborhood` has no effect on the
  selection under `sg_focus = "hr"`.

---

## 8. Deliverable

`dev/selection-eval/FINDINGS.md` containing:

- the Phase 1 verdict (premise valid or not);
- the full `sg_focus` x `selection_rule` matrix;
- the `hrMinSG` degeneracy rate, with a counterexample if the analytic claim
  is refuted;
- the preview-truncation verdict (does it bind; does it change outcomes);
- a documentation discrepancy list with proposed corrections;
- a section on **what could not be determined and why**.

Do not propose or apply package changes. Do not open a PR.

---

## 9. Priors to test, not confirm

All of the following come from **source reading only** — no execution. Report
refutations prominently.

| # | Prior | Confidence |
|---|-------|-----------|
| 1 | `.pareto_dominated_xy()` maximizes both axes as the doc states | high (read directly) |
| 2 | `.compute_inclusion_band()` is unreachable for `hr`/`eff`/`maxeff`/`maxSG`/`minSG` | high (read directly) |
| 3 | min-$N$ candidate survives `"pareto"` only if it has the strict max effect | medium — analytic, untested |
| 4 | `hrMinSG` + `"pareto"` degenerates toward max-effect selection | medium — follows from 3, untested |
| 5 | Preview truncation binds on GBSG | **low — pure speculation** |
| 6 | Alias pairs select identically | medium |
