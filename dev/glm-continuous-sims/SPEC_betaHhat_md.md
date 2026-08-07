# SPEC: `betaHhat_truth_md.R`

Mean-difference sibling of
`quarto/simulations/actg175/binary/betaHhat_truth_glm.R`.

**Target path** `quarto/simulations/actg175/continuous/betaHhat_truth_md.R`

Basis: `NOTE_target_is_collapsibility.md`. Read it first — this spec depends on
its §2 and §3 and contradicts the originating handoff on both.

## What differs from the two existing siblings

The survival and OR modules both evaluate a fitted model on an evaluation frame,
because the hazard ratio and the odds ratio are non-collapsible. The mean
difference is collapsible, so the target is a finite mean over `dgm$df_super` and
is exact.

Three structural consequences:

- **No evaluation frame.** `build_eval_frame_glm()` has no counterpart. Delete
  the concept rather than port it; `eval_seed` disappears with it.
- **No model fit.** No `lm()`, no `glm()`. Separation and minimum-cell guards
  have no analogue; the only degenerate case is an empty region.
- **No Monte Carlo error in the target.** Repeated evaluation is bit-identical.

## Suffix convention

Every exported name carries `_md`, for the reason stated in
`betaHhat_truth_glm.R`'s header: a mistaken co-source of two modules must fail
loudly with object-not-found rather than silently score an MD target with Cox or
GLM code. Column names stay scale-agnostic and identical across all three modules
(`betaHhat_H`, `betaHhat_Hc`, `nH_eval`, `nHc_eval`) so the downstream
`paste0("betaHhat_", suffix)` scoring code is shared verbatim.

## The region computation

Delegate to the package. `compute_aor()` is exported, accepts an arbitrary
`subset_indicator`, and dispatches `"MD"` to `mean(mu1) - mean(mu0)`:

```r
.beta_region_md <- function(df_super, idx) {
  if (!length(idx) || !any(idx)) {
    return(NA_real_)
  }
  compute_aor(df_super,
              subset_indicator = as.integer(idx),
              effect_measure   = "MD")
}
```

Verified against a hand-computed target on an arbitrary region overlapping $Q$:
agreement to `0.00e+00` for both logical and integer indicators, and on the
no-subset path. Passing `effect_measure = "OR"` against the same columns returns
0.9939 rather than 21.6737, confirming the check discriminates.

Do **not** restate `mean(mu1) - mean(mu0)` inline. Under the code–theory
alignment standard a quantity reconstructed at a second site is a defect, and the
package already carries two sites for this arithmetic (`compute_aor()` and the
local `.effect_one()` closure inside `.compute_glm_effects()`). Do not add a
third, and do not refactor the existing two under this spec.

## Call surface

Unchanged from the engines' point of view, so `attach_betaHhat_md()` remains "the
one call each engine adds to `run_cell()`". Deduplication is by distinct
`sg_def`, not by replicate.

    .beta_region_md(df_super, idx)
    betaHhat_one_md(rule, df_super, focus)
    betaHhat_table_md(sg_defs, df_super, focus)
    attach_betaHhat_md(results, df_super, focus)
    betaHhat_theta_dagger_check_md(df_super)

`betaHhat_one_md()` accepts either the named `sg.harm` character vector or the
`" & "`-joined `sg_def` string, and resolves membership through `get_dfpred()`.

**Corrected.** This section previously read "GRF disjunctions live inside a
single element and survive the split", and instructed mirroring
`betaHhat_one_or()`. Both were wrong. A GRF disjunction carries `" & "` inside
each conjunct, so a `" & "` split shreds it:

    "(age > 34 & preanti <= 744.5) | (wtkg <= 60)"
      -> c("(age > 34", "preanti <= 744.5) | (wtkg <= 60)")

`get_dfpred()` then fails to resolve it and the target comes back `NA` rather
than erroring. The fix is dispatch order — test the disjunction form **first**,
on the unsplit string, exactly as `get_dfpred()` itself does
(`R/forestsearch_helpers.R:101`). Both `betaHhat_truth_md.R` and
`betaHhat_truth_glm.R` now do this.

The defect was **latent, not live**: `.grf_build_subgroup_definition()` emits
the `" | "` form only for a multi-leaf policy tree, and 0 of 250,025 non-empty
`sg_def` rows in committed artifacts under `quarto/` contain `"|"`. Multi-leaf
selection is reachable at `maxdepth >= 2`, so no data repair was needed but the
correction had to land before any sweep that produces one.

The deeper fix remains open: the package's own answer is
`.grf_evaluate_subgroup()` on the structured `sg_def`
(`R/forestsearch_helpers.R:1778, 1781, 1784`), chosen precisely because
`get_dfpred()` cannot evaluate the disjunction string (`:1772`). All three
`betaHhat_truth_*.R` modules re-parse the label string instead, which is why
this defect class exists at all.

## The membership convention is the one open line

`betaHhat_one_or()` reads `pred$treat.recommend == 0L` as in-region, which is the
**harm** convention. The continuous template searches for *benefit* via treatment
switching, so its realized regions are $\widehat G/\widehat G^{c}$ and the sense
inverts.

Implement with an explicit argument rather than a hard-coded comparison, so the
choice is visible at the call site and in the bundle:

```r
in_region <- if (focus == "harm") {
  pred$treat.recommend == 0L
} else {
  pred$treat.recommend == 1L
}
```

Do not default this silently. Until the harm/benefit question is settled, a
missing `focus` should be an error, not a guess.

## Acceptance criteria

Integration-style, as a Quarto living document. No `testthat` scaffolding. Every
assertion must be an exact identity, and each must be paired with a
plausible-but-wrong form asserted to fail — otherwise the check cannot
distinguish a correct implementation from one that happens to be close.

1. **Exactness.** Repeated evaluation of `.beta_region_md()` on one region is
   bit-identical. Assert `max - min == 0`, not a tolerance.
2. **Agreement with a fitted MD.** On a large simulated trial the fitted
   $\bar Y_1 - \bar Y_0$ on the region agrees with the target to Monte Carlo
   error, and the difference shrinks at $n^{-1/2}$ across at least three trial
   sizes. A single comparison cannot distinguish Monte Carlo error from a small
   systematic offset.
3. **The wrong form fails.** Assert that $\beta_{\mathrm{inter}} \cdot
   \mathbb P(Q \mid \widehat H)$ — the handoff's formula — differs from the
   fitted MD by $\delta$ on every region, including a region disjoint from $Q$
   where it gives $0$ against a truth of $\delta$.
4. **$\theta^{\dagger}$ identity.** `betaHhat_theta_dagger_check_md()` at
   `flag_harm` reproduces the DGM's own `effect_Q` and `effect_Qc` to
   `0.00e+00`. This is an exact identity, not agreement to Monte Carlo error; a
   tolerance here would hide a real defect.
5. **Rule shapes.** At least one disjunction, one negation, and one region
   disjoint from $Q$, since the closed form is shape-free and the check should
   demonstrate that rather than assume it. Every region must have a target
   distinguishable from $0$, from $\delta$, and from
   $\delta + \beta_{\mathrm{inter}}$ — a rule returning one of those by
   coincidence makes the assertion vacuous.
6. **Reachability.** At least one assertion must run on regions produced by an
   actual `forestsearch()` call rather than hand-constructed masks. A control
   that passes because the quantity is easy is not evidence the code path is
   reachable.

`R CMD check --as-cran` clean is the bar for the commit.
