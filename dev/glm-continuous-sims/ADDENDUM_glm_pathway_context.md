---
title: "Addendum: what the audit session established about the GLM pathway"
subtitle: "Companion to HANDOFF_glm_continuous_simulations.md. No new claims."
date: 2026-08-06
bibliography: []
---

## Why this exists

The handoff describes what the closed forms can verify. This records what the
audit session that preceded it *found*, specifically about the GLM and
continuous pathways, because three things there are counterintuitive enough
that assuming otherwise would cost time.

Nothing here is new work. It is context the parallel workstream is operating
downstream of without necessarily knowing it.

## 1. The GLM pathway is the most-audited part of the release, not the least

The natural assumption runs the other way --- survival is the manuscript
pathway and looks like the mature one. That is backwards for defect coverage.

Survival got the regeneration: five documents re-run, attribution verified
against `selection_rate`, numbers confirmed. **The GLM pathway got the
defect-finding.**

| finding | where it lives | how it was measured |
|---|---|---|
| F1 | `.dfbeta_glm()` on `lm()` fits | MD/continuous only --- the sole measure fitted with `lm()` |
| F2/F3 | identifier vs MR admission floors | measured on ACTG175 binary/OR |
| F12 | CV reusing cached fits | measured on ACTG175 binary/OR |
| F4 | propensity adjustment | GLM-only by construction: the estimator-closure rebuild is `is_glm`-gated |
| `ps_hat` misalignment | bootstrap replicates | **only observable on GLM** --- see §3 |
| the closed-form fixture | `fs_mr_inference()` | built entirely on MD |

So when the parallel workstream scrutinises the GLM pathway, it is largely
*confirming* rather than discovering. That should raise the prior on any
anomaly it finds being real rather than a misunderstanding --- the obvious
defects have been swept.

## 2. `effect_measure = "MD"` has an unusually short history

Two facts that belong together.

**Until `411a448`, MD produced no MR correction at all.** `.dfbeta_glm()` read
`fit$weights` unconditionally; an unweighted `lm()` has `NULL` there, so the
computation was non-conformable and errored on every call. MD is the only
measure fitted with `lm()`. The failure was silent at every level:
`.consistency_glm_pieces()` returned `NULL`, `consistency_resample()` returned
`NA` rates, the consistency engine fell back to literal splitting while
`consistency_method = "resample"` was in force, and `fs_mr_inference()` could fit
no candidate --- returning the short variant with `n_family = 0` and no
`debiased` element.

**Any MD result predating that commit is empty**, and the `NEWS.md` entry says
so.

**And MD is now the only outcome family with closed-form assertions on its
influence path.** The fixture asserts $\widehat\beta$, $\mathrm{db}_{g,i}$,
$\sigma_{D,g}$ and $\Sigma$ against exact algebraic forms at $10^{-10}$.

The reason it is MD-only is worth internalising rather than treating as a
scheduling accident. For logistic the one-step influence is
$(X'WX)^{-1}x_ie_i/(1-h_i)$ with IRLS weights --- which is what the code
computes. Writing that down and comparing is re-implementation. MD is the only
family where the reference comes from a genuinely different computation, arm
sums of squares.

So the simulation workstream is building on a path that went from broken to
best-verified in a single day. Both halves matter: the verification is real, and
it is recent enough that nothing has stress-tested it at scale. That is exactly
what a simulation would do.

## 3. Two measurement traps, both of which would recur in a harness

Each cost time in the audit session and each has a general form.

### A defect in an input is only observable where the input is consumed

The bootstrap carried a user-supplied `ps_hat` into replicates **unaligned** ---
a resample has the same row count but different subjects, so the length check
passed while every score belonged to a different person. Measured misalignment
of ~0.52 on scores constructed over a 0.60-wide range: not slightly off,
effectively unrelated.

The natural place to measure it was the survival frame. **That would have shown
zero change**, because no Cox path reads `sw`, `ps_hat` or `ips_covar` --- the
score is estimated, attached, and consumed by nothing. A survival measurement
would have read as a clean refutation of a real defect.

Moved to ACTG175 binary/OR, the full bootstrap estimate moved 2.0919 to 2.2434
(+7.2%) with a materially wider interval.

**General form:** before concluding a defect is absent, confirm the quantity it
corrupts is actually read on the path being measured.

### A degenerate configuration makes assertions vacuous

Three instances, all found by someone noticing the check could not fail:

- **`selection_rate = 1`.** The two denominator conventions coincide
  identically there, so every assertion about the convention passes without
  testing it. The fixture requires `0.3 < selection_rate < 0.9` and asserts the
  rate *before* asserting the residual, so a change driving the rate to 1 fails
  loudly rather than passing silently.
- **An empty complement.** A single candidate spanning the whole frame has no
  complement, so `mean_r_c` degrades to `NA` --- correctly. But the obvious way
  to measure it then asserts nothing. The fixture uses a strict-subset
  candidate, with a separate test pinning the `NA` case so the degradation stays
  deliberate.
- **A family of one.** At $M = 1$ the winner is the same candidate on every
  draw, so selection-bias assertions are trivially satisfied. It is a useful
  *boundary* check --- $\widetilde V \to 4\sigma^2_D$ is the only assertion that
  catches a dropped $\mathrm{bias_{fix}}$ --- but useless as a general one.

**General form:** every scenario needs a value distinguishable from zero, from
$\sigma_D$, and from what the naive estimate would give. A control that passes
because the quantity is easy is not evidence the code path is reachable.

## 4. Current MR behaviour, which any simulation must target

Three changes landed during the audit and all three move numbers.

**Both bias terms share one denominator, conditional on identification.** Each
averages over the draws that produced a winner, and the infinitesimal jackknife
runs on that same set. `mean_r` and `mean_r_c` are returned and sit at machine
epsilon when this holds. Note the consequence: `bias_fix` is **not** mean-zero
under conditioning --- conditioning on admission selects draws with upward
perturbation, so it is a genuine selection effect, monotone in the selection
rate and approaching zero as that rate approaches 1.

**MR's family is no longer restricted** to a neighbourhood of the best native
statistic. On one application this moved a GRF family from 2 candidates to
1289; on another, from 1 to 858. Estimated selection bias roughly doubled.

**The identifier and MR admit on the same floor**, read from one resolved
object, and an empty admitted set returns no subgroup rather than the native
winner.

Any simulation comparing against pre-audit results will see all three at once.
Attribution needs `selection_rate` recorded alongside every MR result --- that
was the single most useful diagnostic in disentangling the regeneration.

## 5. What is still open, and why none of it should block

Four audit findings remain open: F5 (MR's family is a superset of the
identifier's), F6 (truncation warned but recorded nowhere), F7 (the inclusion
band empties on a negative maximum effect, with three consumers diverging), F8
(display sites re-implement the band). Plus F4's Q2 and Q3, deferred by
decision.

**F7 is the only one that touches the continuous pathway**, and only in a
degenerate case: a negative maximum effect means no candidate in the family
shows harm, which is a run that should report nothing regardless. The three
consumers disagree about *how* to say nothing. It is reachable for RD and MD
because ratio measures are bounded below by zero.

Judged not worth fixing. Recorded here so that if a simulation produces an empty
band it is recognised rather than investigated from scratch.
