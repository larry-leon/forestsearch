---
title: "Handoff: GLM/continuous simulations with closed-form checks"
subtitle: "Opens a new chat. Design first; nothing has been built."
date: 2026-08-06
bibliography: []
---

## What this is

A new workstream, not a continuation. The preceding chat completed an
audit-and-fix cycle on `forestsearch`'s multiplier-resampling correction and
built a closed-form fixture for it. This one asks a different question: **can
simulations over the GLM/continuous pathway be run against closed-form
references rather than only against each other?**

The short answer, established below with numbers rather than argument: yes, and
the mean-difference case is uniquely well-suited to it in a way no other outcome
family is.

**Design comes before implementation.** The preceding cycle produced three spec
errors that were caught by measurement rather than review --- a stated
asymptotic identity asserted as exact, an envelope that contradicted a decision
made later in the same document, and a claim that an argument could pin cut
locations when it takes variable names. All three were cheap to fix because they
were caught against small code. A simulation harness built on unverified closed
forms would not be.

## State of the codebase

Branch `mr-in-replicates` on `github.com/larry-leon/forestsearch`. Everything
below is pushed and `R CMD check --as-cran` clean.

The MR correction was corrected in three places during the preceding cycle and
its current behaviour is what any simulation must target:

- The identifier and MR now admit on the same floor, $\widehat\beta(g) \ge t_g$,
  read from one resolved admission object. An empty admitted set returns **no
  subgroup** rather than the native winner.
- Both bias terms share one denominator, conditional on identification: each
  averages over the draws that produced a winner, and the infinitesimal
  jackknife runs on that same set. `mean_r` and `mean_r_c` are returned and are
  zero at machine epsilon when this holds.
- MR's candidate family is no longer restricted to a neighbourhood of the best
  native statistic. On one application this moved a GRF family from 2
  candidates to 1289 and roughly doubled the estimated selection bias.

`dev/closedform-fixture/closedform_fixture.qmd` is the existing fixture, and
`dev/identifier-alignment/verification/` holds four reference scripts. Read
both before designing anything --- they define what is already checkable, and
the new workstream should extend rather than duplicate.

## The infrastructure already exists

**This is not a greenfield workstream.** Three sets of simulation documents are
already in the repository, and one of them is already continuous:

| location | pathway |
|---|---|
| `quarto/simulations/gbsg_redux/` | survival --- manuscript simulations |
| `quarto/simulations/actg175/binary/` | binary --- manuscript simulations |
| `quarto/simulations/actg175/continuous/` | **a single file, believed an early template** |
| `dev/identifier-alignment/rerun/` | the three `maxcons` re-run documents |

**Read the continuous template first.** Whatever it is, it is the closest thing
to a starting point, and designing around it is cheaper than designing past it.

More important is `dev/identifier-alignment/rerun/betaHhat_truth.R`, which
already implements the conditional estimand and defines the architecture any new
work should follow:

- **Producer/consumer split.** Each engine computes $\beta(\widehat H)$ once per
  *distinct realized rule* on a fixed evaluation population; the manuscript
  fragments are pure consumers. `attach_betaHhat()` is "the one call each engine
  adds to `run_cell()`".
- **Deduplication by rule**, not by replicate --- `betaHhat_table()` computes
  each distinct `sg_def` once.
- **A sanity gate.** `betaHhat_theta_dagger_check()` evaluates
  $\theta^\dagger$ on the same frame at the true harm flag and should reproduce
  `dgm$hr_H_true` to Monte Carlo error, confirming the evaluation frame reads as
  "population".
- **One file sourced verbatim by both engines**, so the MR grid and the FB
  batches score an identical target.

That architecture is sound and should be reused. What changes for the continuous
pathway is what sits inside `.beta_region()` --- and, more consequentially, how
the target is computed at all.

## What closed forms can and cannot confirm

Three layers, in decreasing strength. This is the load-bearing section.

### Layer 1 --- exact, no Monte Carlo

For the mean difference the per-candidate objects are algebraic identities,
assertable at $10^{-10}$:

$$
\widehat\beta(g) = \bar Y_1 - \bar Y_0, \qquad
\mathrm{db}_{g,i} = \frac{e_i}{n_1 - 1}\ \text{or}\ -\frac{e_i}{n_0 - 1},
$$
$$
\sigma^2_{D,g} = \frac{\mathrm{SS}_1}{(n_1-1)^2} + \frac{\mathrm{SS}_0}{(n_0-1)^2},
\qquad
\Sigma = \mathtt{crossprod(B)}.
$$

Note the $(n_a - 1)$: the package's one-step influence carries a $1/(1-h_i)$
leverage factor and $h_i = 1/n_a$ in a two-group model. The Neyman form
$\widehat s_1^2/n_1 + \widehat s_0^2/n_0$ is the correct **asymptotic**
statement and the wrong **exact** one --- it differs by 1.0% at $n_a = 100$,
which is $10^{8}$ times the tolerance. This was a spec error in the preceding
cycle, caught by the fixture on first contact. Do not reintroduce it.

**Why this layer is MD-only.** For logistic the one-step influence is
$(X'WX)^{-1}x_ie_i/(1-h_i)$ with IRLS weights --- which is what the code
computes. Writing that down and comparing is re-implementation, not
verification. MD is the only family where the reference comes from a genuinely
different computation, arm sums of squares, rather than the same matrix algebra
restated.

### Layer 2 --- closed form as a target, Monte Carlo in the comparison

$\tau\varphi(d/\tau)$ for two candidates; $\sigma\,\mathbb{E}[\max]$ for $M$
orthogonal; $\sigma\varphi(c/\sigma)/\overline\Phi(c/\sigma)$ under an active
floor. Exact as targets. The estimate converges at $B^{-1/2}$, and that rate is
itself a check --- a systematic offset that does not shrink is a defect no
single comparison reveals.

### Layer 3 --- boundaries only

$\widetilde V \to 4\sigma^2_{D}$ at $M = 1$; $\bar r = 0$ always;
$\widetilde V$ invariant to a constant shift in the residual, always.

**The honest limit: for $M \ge 2$ with general overlap, $\widetilde V$ has no
closed form.** The interval's *width* is checkable at boundaries and not in
general. Any design that implies otherwise is wrong.

## What simulation adds that the fixture cannot

This is the reason the workstream is worth running, and it rests on one property
peculiar to MD under Gaussian multipliers.

$D \sim \mathcal{N}_M(0, \Sigma)$ **exactly**, at every $n$ --- no asymptotics,
because $D$ is a fixed linear combination of independent normals. And for MD,
$\Sigma$ is known exactly. So the $B \to \infty$ target can be evaluated
independently, by sampling the *law* rather than running the pipeline harder.

Measured, at pipeline scale $B = 20{,}000$:

| $M$ | $\rho$ | oracle SE | pipeline SD | ratio |
|---|---|---|---|---|
| 2 | 0.0 | 2.1e-04 | 5.4e-03 | 26x |
| 2 | 0.7 | 1.8e-04 | 4.8e-03 | 27x |
| 5 | 0.0 | 1.1e-04 | 3.5e-03 | 32x |
| 5 | 0.7 | 2.3e-04 | 5.7e-03 | 24x |
| 20 | 0.0 | 1.5e-04 | 4.0e-03 | 26x |
| 20 | 0.7 | 2.2e-04 | 7.0e-03 | 31x |

The reference is 24--32x sharper than the quantity it judges, across family
sizes and correlation structures. **The comparison is limited by the pipeline's
$B$, not by the reference.**

That buys an error decomposition nothing else provides: per replicate, split
the discrepancy into finite-$B$ Monte Carlo versus estimator behaviour, because
the first is now measurable rather than confounded with the second. In survival
this is impossible --- the influence is a one-step approximation, so $\Sigma$ is
itself approximate and the "oracle" would inherit the error it is meant to
isolate.

Two supporting facts, both checked:

**$\Sigma = \mathtt{crossprod(B)}$ stays exact under overlapping memberships.**
With three candidates sharing ~100 subjects pairwise, the analytic $\Sigma$ and
a 400,000-draw empirical estimate agree to 5.2e-05 against a maximum entry of
2.1e-02 --- Monte Carlo error, not bias. Overlap does not degrade the property.

**The conditional estimand becomes analytic, replacing a numerical oracle.**
This is the sharpest gain and it is specific to the continuous pathway.

`betaHhat_truth.R` computes $\beta(\widehat H)$ *numerically*: it builds a fixed
super-population, evaluates `coxph(Surv ~ treat)` on the realized region, and
treats the result as the population value. It has to --- for survival there is
no closed form. The file's own comment concedes the cost: "residual Monte-Carlo
noise is negligible at this size."

For the mean difference with an additive DGM there is no such concession to
make. The target is

$$\beta(\widehat H) \;=\; \tau \cdot \mathbb{P}\bigl(H_{\mathrm{true}} \mid \widehat H\bigr),$$

a ratio of rectangle probabilities --- and under independent covariates an
axis-aligned rectangle's probability is a product of marginals, so this is
closed form. Verified against 4,000,000 draws across the rule shapes the engines
actually realize:

| realized rule | closed form | empirical |
|---|---|---|
| exact rule | 0.80000 | 0.80039 |
| too narrow | 0.80000 | 0.79811 |
| too wide | 0.42857 | 0.42806 |
| one cut only | 0.48000 | 0.47915 |
| spurious third covariate | 0.80000 | 0.79996 |
| negation form | 0.80000 | 0.80039 |
| disjoint from truth | 0.00000 | 0.00065 |

Every diff is Monte Carlo error. **Mis-specification is handled exactly** --- the
target is correct whether the selected rectangle is too narrow, too wide,
missing a covariate, carrying a spurious one, written with negations, or
disjoint from the truth. That last row matters most: a rule capturing none of
the true harm region has $\beta(\widehat H) = 0$ exactly, and a numerical oracle
would return noise around zero.

Three consequences for the architecture. No super-population is needed for the
target, so the fixed evaluation frame and its `eval_seed` disappear. The residual
Monte Carlo noise in the target is not merely negligible but **zero**. And
`betaHhat_theta_dagger_check()` upgrades from an approximate agreement test to
an exact identity, which makes it a much stronger gate.

**The limits, stated honestly.** Product-of-marginals requires *independent*
covariates and *conjunction* rules. Correlated covariates need multivariate
normal rectangle probabilities --- `mvtnorm::pmvnorm`, numerical but to far
higher precision than a super-population. GRF's disjunctions need
inclusion-exclusion, which is closed form but grows with the number of terms. A
design assuming closed form unconditionally would be wrong.

## Design questions to settle before any code

Each of these changes what gets built. None is settled.

**What is the simulation for?** Three plausible aims and they need different
designs: verifying the correction is unbiased for $\beta(\widehat H)$;
characterising coverage of the IJ interval; or measuring how the correction
degrades as the family grows. The oracle helps all three but the DGM and
replicate counts differ.

**Which DGM?** The rectangle-subgroup Gaussian model above has an analytic
conditional estimand, which is its main virtue. Whether the candidate family
should be *constructed* (memberships supplied directly, so $\mathcal{F}$ is
exactly known) or *searched* (`forestsearch()` generates it, so the simulation
exercises the real pipeline) is the central trade: the first isolates the
correction, the second tests the thing that ships.

**Extend `betaHhat_truth.R` or write a sibling?** The architecture is worth
reusing verbatim --- producer/consumer, dedupe by distinct rule, one file
sourced by both engines. But the continuous version does something structurally
different: it computes the target from the DGM's parameters rather than
evaluating a model on an evaluation frame, so `build_eval_frame()` has no
counterpart. A sibling file with the same call surface (`attach_betaHhat()`
unchanged from the engines' point of view) is probably cleaner than branching
inside the existing one. Worth deciding before either is written.

**Does the closed form cover the rules the engines actually produce?** The
verification above used conjunctions and negations on independent uniform
covariates. Whether the DGM's covariates are independent, and whether GRF's
disjunctions appear in the realized rules, decides whether the target is closed
form, `pmvnorm`-numerical, or inclusion-exclusion. Check against the existing
DGM before assuming.

**Where does the error decomposition go?** Reporting it per replicate is more
informative and much larger; reporting it in aggregate is cheaper. The
finite-$B$ component should shrink at $B^{-1/2}$ and the estimator component
should not, which is itself a check on the decomposition.

**What is the relationship to `fs-glms-interpretable`'s GLM-native DGM?** That
was already on the horizon before this workstream. If the two want the same
generator, deciding so now avoids a second one later.

## Working conventions

Carried from the preceding cycle, all earned rather than stipulated.

- **The chat writes specifications; CC writes code.** Code in a handoff
  duplicates what the executor is about to build.
- **Commit spec documents when they arrive, before editing**, so corrections
  show as real diffs. Adopted late in the preceding cycle after four spec
  corrections landed as an untracked file's first commit.
- **Never infer file contents.** Ask for the file. Trace call sites, not
  formals or documentation --- a function's roxygen described behaviour it did
  not have at least twice in the preceding cycle.
- **A control that passes because the quantity is easy is not evidence the code
  path is reachable.** A gaussian-family control passed for weeks while the
  `lm()` path errored on every call.
- **Assert against the exact form and assert the plausible-but-wrong form
  fails.** Otherwise a check cannot distinguish a correct implementation from
  one where the quantity happens to be close.
- **A rate-1 or degenerate configuration makes most assertions vacuous.** Every
  scenario needs a value distinguishable from zero, from $\sigma_D$, and from
  what the naive estimate would give.
- Deliverables are `.qmd` with a PDF target, `scrartcl`, `\mathrm`, LaTeX in
  `$...$`. Testing means integration-style Quarto living documents, not
  `testthat` scaffolding --- with a standing exception for contract tests on
  static properties.
- `R CMD check --as-cran` clean is the bar for a commit; `devtools::test()`
  alone is not.
- Larry's machine is Pop!_OS, ~127 cores, R 4.6.1. The chat sandbox has R 4.3.3
  and no Quarto.

## What has not been done

Nothing. No DGM specified, no harness written, no replicate counts chosen.

The numbers in this document come from standalone scripts run in a sandbox
against the *law*, not against `forestsearch()`. **No claim here has been
checked against the package**, and the preceding cycle's lesson is that the gap
between a correct closed form and what the code returns is exactly where the
errors live. The first implementation step should be small and should measure,
not build.
