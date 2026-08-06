---
title: "Task: expose mean_r, and close the fixture's residual-centering gap"
subtitle: "For CC. Branch mr-in-replicates, after b473257 and c12c75f."
date: 2026-08-06
bibliography: []
---

## Why this is worth doing

`dev/closedform-fixture/closedform_fixture.qmd` records the gap itself, in
`@sec-rbar-gap`: the $\bar r = 0$ assertion cannot be made because
`fs_mr_inference()` returns neither $r_b$, nor $\Xi$, nor the winner vector,
and rebuilding $D_g(b)$ inside the fixture is the first thing the fixture spec
forbids. That reasoning is correct and the substitute you chose --- the
$\widetilde V$ shift-invariance check against the shipped `.fs_mr_ij_var()` ---
was the strongest available. This task removes the need for a substitute.

The stake is specific. $\bar r = 0$ **is** what F13 was: the defect was that the
two bias terms carried different denominators, and its signature was a non-zero
residual mean. The fixture's inverted `fixed_bias` pair catches a reversion to
the unconditional convention, but that is one particular way the denominators
can diverge. Any other way produces $\bar r \ne 0$ and nothing currently
notices. Given that this single quantity carried three separate defects ---
`84203ed`, `8d60b1e`, `dad0415` --- a direct assertion on it is worth more than
most.

## The change

`R/fs_mr_inference.R`. **Read the current lines before editing**; the return
block below is from a pre-`dad0415` snapshot and the file has moved since.

The residual already exists at the point of the IJ call, on both paths: `r_H`
for the selected subgroup and `r_c` for the complement. Both are computed over
the full draw vector and then indexed to the shared set --- `ok_H` and `use_c`
respectively --- inside `.fs_mr_ij_var()`. Add to the returned object:

- `mean_r` --- the mean of `r_H` over `ok_H`, i.e. over exactly the draws the IJ
  uses. **Not** over all $B$: the property being asserted is that the residual
  is centered on the draws entering the variance, and on the excluded draws
  `r_H` is not a meaningful quantity.
- `mean_r_c` --- the same for `r_c` over `use_c`. `NA_real_` when no complement
  was fit, matching how the other complement fields degrade.

Place them beside `selection_rate` in the return list rather than inside
`debiased`; they are diagnostics of the correction's internal consistency, not
properties of the estimate. Roxygen for both should state the invariant --- zero
by construction when both bias terms share a denominator --- and name that as
the thing they exist to detect.

**Do not change any computation.** This is exposure only. If adding it requires
moving where `r_H` or `r_c` is formed, stop and report rather than
restructuring; the F13 change is recent and its shape should not shift for a
diagnostic.

Confirm while you are there that `ok_H` is genuinely the index set passed to
`.fs_mr_ij_var()` on the selected path and `use_c` on the complement path,
rather than inferring it from the variable names --- a formal's name has
described behaviour it no longer had at least once in this project.

## The fixture assertion

`dev/closedform-fixture/closedform_fixture.qmd`, Level 3. Replace the
substitute check --- keep the shift-invariance one as well, it tests a different
property --- with a direct assertion:

$$|\,\overline r\,| \;<\; 10^{-10}
\qquad\text{at every scenario, and } |\,\overline{r}_c\,| < 10^{-10}
\text{ wherever the complement was fit.}$$

Two things about where it must bind.

**S2 is the load-bearing case.** At `selection_rate` $= 1$ the shared-denominator
property holds trivially, since the two candidate draw sets coincide. The
assertion is only informative where the floor binds, so S2 must carry it and its
selection rate must be reported next to it. A pass at rate 1 alone is vacuous,
which is the failure mode the fixture's own closing section warns about.

**Assert that the old convention would fail it.** Compute what $\bar r$ would
have been with `fixed_bias` averaged over all $B$ while `selection_bias` stays
over the winner draws --- the pre-`dad0415` arrangement --- and assert that value
is **not** within tolerance of zero at S2. Without this, the assertion cannot
distinguish a correct implementation from one where the quantity happens to be
small. `verify_conditional_convention.R` in
`dev/identifier-alignment/verification/` establishes the arithmetic; do not
re-derive it, and do not edit that script.

Update `@sec-rbar-gap` to record that the gap is closed and how, rather than
deleting it. The section is the reason this task exists and the record is worth
keeping.

## Contract test

Add the direct assertion to `tests/testthat/test-mr-closedform-fixture.R` as
well. It is a static property of the arithmetic --- exactly the standing
exception to the no-testthat rule --- and it is the assertion that would have
caught F13 in CI rather than by a bisect three commits later.

The test must construct a configuration with `selection_rate` strictly between
roughly 0.3 and 0.9. Assert the rate falls in that band **before** asserting the
residual, so a future change that inadvertently drives the rate to 1 fails
loudly rather than passing vacuously.

## Measurement

Report, for one binding configuration and one rate-1 control:

`selection_rate`, `mean_r`, `mean_r_c`, and the counterfactual $\bar r$ under
the old denominator arrangement.

Expected: `mean_r` at machine epsilon in both; the counterfactual near zero at
rate 1 and materially non-zero at the binding configuration. If `mean_r` is not
at machine epsilon, stop and report --- that would mean the shared denominator
is not actually shared somewhere, which is a finding rather than a fixture
problem.

## Scope

In scope: the two returned fields and their roxygen, the fixture assertion, the
contract test, the `@sec-rbar-gap` update, a `NEWS.md` entry.

Out of scope: F4, F9, F10, F11, F14; `selection_criteria.qmd`; any analysis
document; any change to how the residual or either bias term is computed.

`R CMD check --as-cran` clean is the bar. Regenerating the fixture document is
required since its assertions change; regenerating anything else is not.

The `NEWS.md` entry is short --- two new diagnostic fields, no change to any
estimate. Say what they are for: the invariant is zero by construction when the
two bias terms share a denominator, and a non-zero value indicates they do not.
