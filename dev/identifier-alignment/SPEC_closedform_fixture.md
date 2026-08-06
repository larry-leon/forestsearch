---
title: "Spec: the mean-difference closed-form fixture"
subtitle: "For CC. Branch mr-in-replicates, after 411a448 and 0302e8c."
date: 2026-08-05
bibliography: []
---

## Scope, and what is already done

Build a fixture that binds `fs_mr_inference()`'s output to closed-form
expressions. **The reference side already exists and is verified** --- do not
re-derive it. Four scripts are supplied:

| Script | What it establishes |
|---|---|
| `verify_closedform_fixture.R` | V1--V6: the closed forms themselves, checked by Monte Carlo against the multiplier construction |
| `verify_residual_centering.R` | V7: the residual-centering identity, run against the package's own `.fs_mr_ij_var()` |
| `verify_eq9_alignment.R` | V8: as-coded vs Eq. 9 literal vs half-aligned; establishes the as-coded arm is not centered |
| `verify_conditional_convention.R` | V9: the conditional option V8 did not consider --- the one `dad0415` adopted |

Three scripts, not four. `ij_var_verbatim.R` was **deliberately not filed**: all
three package-dependent scripts bind `forestsearch:::.fs_mr_ij_var` directly
rather than a snapshot, and the verification README records why. Keep that
property in anything you write.

`check_dina_actg.R` --- the closed form checked against the pre-fix
DINA/ACTG175 numbers --- sits one level up, in `dev/identifier-alignment/`, not
in `verification/`.

Those scripts verify that **the mathematics is right**. They do not touch the
package pipeline. This task is the missing half: does `fs_mr_inference()`, run
end to end, reproduce those values?

The reference scripts are already filed at
`dev/identifier-alignment/verification/`, with their own README. They bind
`forestsearch:::.fs_mr_ij_var` directly rather than copying it; keep that
property in anything you write. Do not edit them as part of this task.

**Do not change any package code.** The MR correction is settled as of
`dad0415` and this task asserts against it, not about it. If an assertion
fails, that is a finding to report --- not a licence to adjust the code until
it passes.

## Why the mean difference

Three reasons, all load-bearing.

The MD score is linear in $(\alpha, \beta)$, so Lemma 1's mean-value curvature
term vanishes *identically*. The dfbeta linearization carries no curvature error
at all; the only remainder is the perturbed information, entering at
$O_p(n_g^{-1})$ through the arm counts. Every other outcome type has a curvature
term that must be bounded rather than computed.

Every per-candidate object is closed form. **Corrected after measurement** ---
the forms below previously omitted the leverage factor and would have failed on
contact at the stated `1e-10`. `.dfbeta_glm()`
(`R/consistency_resample.R:168-176`) returns the one-step influence
$(X'X)^{-1}x_i e_i/(1-h_i)$, and for a two-group model $h_i = 1/n_a$, so

$$\mathrm{db}_{g,i} = \frac{e_i}{n_1 - 1}\ \text{(treated)},\quad
  -\frac{e_i}{n_0 - 1}\ \text{(control)},\qquad
  \sigma^2_{D,g} = \frac{\mathrm{SS}_1}{(n_1-1)^2} + \frac{\mathrm{SS}_0}{(n_0-1)^2}.$$

The Neyman form $\widehat s_1^2/n_1 + \widehat s_0^2/n_0$ is the correct
*asymptotic* statement and the wrong *exact* one: it differs by
$O(n_a^{-1})$ relative --- 0.2669415 against the package's 0.2696379 at
$n_a = 100$, a 1.0% gap. The fixture asserts the exact form, and additionally
asserts the Neyman form does **not** match, so a future change toward it is
caught rather than absorbed.

And it is the path F1 broke. `effect_measure = "MD"` is the only measure fitted
with `lm()`; an MD fixture asserting $\sigma^2_{D,g}$ fails **on contact** when
that path is down, rather than degrading silently to `n_family = 0`. The
gaussian control that passed while the path was broken could not have done this,
because it exercised the family and not the path.

## Configuration the fixture must set explicitly

**`multiplier = "gaussian"`.** The package default is `"poisson"`. Under
Gaussian multipliers $D_g(b)$ is exactly multivariate normal at every $n$ --- no
asymptotics, because the linearization is not being invoked, only the law of a
fixed linear combination of independent normals. Every closed form below assumes
this. A fixture left on the Poisson default will miss by a few percent and the
miss will look like a defect.

**Synthetic data, constructed not sampled.** Build the analysis frame so the arm
means, arm variances and membership overlaps are known before any fit runs.
Candidate memberships should be supplied directly rather than produced by a
search, so $\mathcal{F}$ is exactly what the fixture intends.

**$B$ large enough to separate signal from Monte Carlo.** The comparisons below
are $O(B^{-1/2})$; $B \ge 20{,}000$ keeps the envelope well inside the
quantities being checked at these sample sizes.

## The four scenarios

Named, because each isolates a different failure.

**S1 --- two candidates, floor inactive.** Both admitted on every draw. This is
the main selection-bias check.

**S2 --- $M = 1$, floor active.** One candidate, `hr.threshold` set so the floor
binds on roughly a third of draws. This is the shape of the pre-fix
DINA/ACTG175 row, and it is the only scenario that can distinguish the adopted
conditional convention from the unconditional one --- see @sec-convention. It is
therefore the load-bearing scenario, not an edge case.

**S3 --- $M = 1$, floor inactive.** Every draw has a winner. The degenerate IJ
case.

**S4 --- $M$ orthogonal candidates, disjoint membership, no separation.** For
$M = 2, 3, 4$. Assert $\min_{g \ne h}\Sigma_{gh} \ge 0$ *before* comparing to the
$\mathbb{E}[\max]$ reference, since the bound is Slepian's and holds under
non-negative pairwise covariance rather than universally.

## What to assert, in three levels

Level 1 must pass before Level 2 means anything; same for 2 before 3. Report
them in that order and stop on the first failure.

### Level 1 --- the per-candidate objects

Read from `.fs_mr_assemble()`, which returns `B`, `beta_hat`, `sigma_D`.

- $\widehat\beta(g)$ against the arm-mean difference. Exact.
- `sigma_D[g]` against
  $\sqrt{\mathrm{SS}_1/(n_1-1)^2 + \mathrm{SS}_0/(n_0-1)^2}$, the one-step
  influence norm. Exact. (Corrected --- see "Why the mean difference".)
- Column $g$ of `B` against $e_i/(n_1-1)$ / $-e_i/(n_0-1)$, and its within-arm
  sums against zero. Exact.
- $\Sigma = $ `crossprod(B)` against $\sum_i \mathrm{db}_{g,i}\mathrm{db}_{h,i}$
  computed from the memberships. Exact, and it is the same object.

These are algebraic identities, so assert equality at `1e-10`, not an envelope.
A failure here is a broken influence path, which is F1's signature.

### Level 2 --- the correction

- **S1**: `selection_bias` against $\tau\varphi(d/\tau)$, with
  $d = \widehat\beta(g_1) - \widehat\beta(g_2)$ and
  $\tau^2 = \sigma_1^2 + \sigma_2^2 - 2\sigma_{12}$, all from Level 1. Envelope
  $4\max_j\sigma_j B^{-1/2}$.
- **S1 rate check**: repeat at $B$, $4B$, $16B$ and confirm the error falls at
  $B^{-1/2}$. A systematic offset that does not shrink is the failure a single
  comparison cannot see.
- **S1 tie sanity**: at $d = 0$ with orthogonal unit-scale influences the value
  is $\sigma/\sqrt\pi$, the expected maximum of two normals.
- **`fixed_bias`, the inverted pair.** The mean-zero envelope
  $|\widehat{\mathrm{bias}}_{\mathrm{fix}}| \le z_{0.975}\,\sigma_{D,\widehat H}B^{-1/2}$
  holds at `selection_rate` $= 1$ **and only there**. Assert the envelope at
  S1, S3 and S4; at S2 assert it is **breached**, and that the conditional
  value is monotone in $p$, approaching the unconditional value as $p \to 1$.
  Vary the floor on a **fixed** draw matrix --- redrawing per floor confounds
  the trend with sampling noise.

  This corrects an earlier version that asked for the envelope in every
  scenario, written before the F13 decision and not reconciled with
  @sec-convention, which states that `fixed_bias` is *not* mean-zero under
  conditioning. The inverted pair is the strongest assertion in the fixture: an
  implementation that reverted to the unconditional convention passes at rate 1
  and fails S2, and nothing else detects that.
- **S4**: `selection_bias` against $\sigma\mathbb{E}[\max(Z_1,\dots,Z_M)]$,
  computed by quadrature ($\pi^{-1/2}$, $\tfrac32\pi^{-1/2}$, $1.02938$ for
  $M = 2, 3, 4$), after the $\Sigma \ge 0$ gate.

### Level 3 --- the interval

- **S3**: `var_ij` / $\sigma^2_{D,\widehat H} \to 4$. This is the assertion that
  certifies **both** terms of the residual reach the variance --- an
  implementation dropping $\widehat{\mathrm{bias}}_{\mathrm{fix}}$ returns
  $1$ here and passes every other check in this document. Report the ratio at
  three values of $B$ and confirm it approaches 4 rather than asserting a point
  value at one $B$.
- **Every scenario**: $\bar r$ over the draws the IJ uses, asserted **against
  zero** at $10^{-10}$.

  This bullet previously said to assert against the two-bracket identity
  $(\widehat{\mathrm{bias}}_{\mathrm{sel}} - \overline{D_{\widehat H^*}}|_\mathcal{B}) +
   (\widehat{\mathrm{bias}}_{\mathrm{fix}} - \overline{D_{\widehat H}}|_\mathcal{B})$
  rather than against zero, because $\bar r$ was zero only at
  `selection_rate` $= 1$. That was true of the pre-`dad0415` code and is false
  now: under the conditional convention each bracket vanishes identically, so
  the identity reduces to $0 = 0$ and asserting it **passes vacuously** while
  detecting nothing. Assert $|\bar r| < 10^{-10}$ directly, at every selection
  rate. Record `selection_rate` next to it in every row so a vacuous rate-1 pass
  is visible as such.
- **S2 and S3**: confirm $\widetilde V$ is unchanged when a constant is added to
  every $r_b$. Row-centering makes this exact, so it is a wiring check on the
  centering rather than a statistical one.

## The denominator convention: assert it {#sec-convention}

This section previously said the convention was undecided and told you to report
rather than assert. **It is decided.** The package implements the **conditional**
convention as of `dad0415`: both bias terms average over the draws that produced
a winner, and the IJ runs over that same set. Assert against it.

The defect that change fixed was never a denominator on its own --- it was that
the two terms carried *different* denominators, so the residual mixed
differently normalised quantities and $\bar r \ne 0$. Two repairs were coherent
and both centered; they differ in estimand, and conditioning was chosen because
the reported analysis exists only when a subgroup was identified. The reasoning
is in `dev/identifier-alignment/NOTE_f13_open_questions.md` and is not to be
relitigated here.

**S2 asserts the conditional closed form.** With $M = 1$, $c = t_g -
\widehat\beta(g)$ and $\sigma = \sigma_{D,g}$:

$$
\widehat{\mathrm{bias}}_{\mathrm{sel}}
\;\longrightarrow\;
\frac{\sigma\,\varphi(c/\sigma)}{\overline\Phi(c/\sigma)},
\qquad
p = \overline\Phi(c/\sigma).
$$

Envelope $4\sigma/\sqrt{Bp}$ --- the average is over $Bp$ draws, not $B$.

**Assert the unconditional form fails.** $\sigma\varphi(c/\sigma)$ is the
divide-by-$B$ value; the two differ by exactly $1/p$, which is a factor of 3.24
at $c = \sigma/2$. A fixture that only checks the value the code produces cannot
tell a convention change from a regression, so S2 must also assert the observed
value is **not** within tolerance of the unconditional form. Report $p$ alongside
both, so the separation is visible rather than implicit.

**Three further assertions this convention makes checkable**, all exact and all
cheap:

- $\bar r = 0$ over the winner draws, to $10^{-10}$. This is the property whose
  absence *was* the defect, and no other check in the fixture detects a
  denominator mismatch.
- $\widehat V$ equals Wager–Hastie–Efron's centered form exactly, which follows
  from $\bar r = 0$ and holds only under a shared denominator.
- $\widehat{\mathrm{bias}}_{\mathrm{fix}}$ is **not** mean-zero under
  conditioning. Unconditionally $\mathbb{E}[D_{\widehat H}] = 0$; conditioning on
  admission selects draws with upward perturbation, so the conditional value is
  a genuine selection effect. It should be monotone in $p$ and approach the
  unconditional value as $p \to 1$. Assert the monotonicity across three or four
  floors on a **fixed** draw matrix --- redrawing per floor confounds the trend
  with sampling noise, which is a mistake already made once.

**A rate-1 control is mandatory.** At $p = 1$ the two conventions coincide
identically, so every assertion above is vacuous there. S2 must bind on a real
minority of draws, and a separate rate-1 scenario should confirm the conditional
and unconditional forms agree to every digit. A fixture that only ever runs at
$p = 1$ asserts nothing about the convention at all.

## Deliverable

Recommendation, for confirmation before you build:

A Quarto living document under `dev/closedform-fixture/`, rendering one table
per scenario with the closed form, the observed value, the envelope and a
pass/fail column --- **plus** a small testthat file carrying only the Level 1
identities and the S3 ratio. The split is deliberate: Level 1 and the $M = 1$
ratio are static properties that should fail in CI, which is the standing
exception to the no-testthat rule; the Level 2 and 3 comparisons are Monte Carlo
and belong in a rendered document Larry reads rather than in a suite that would
go flaky.

If you would rather put everything in one place, say so before building.

`R CMD check --as-cran` clean is the bar, as usual, and `devtools::test()` alone
is not.

## What would make this fixture worthless

Recorded because each has already happened once in this project.

**Re-implementing the pipeline inside the fixture.** If the fixture computes its
own $D_g(b)$ and compares that to a closed form, it tests arithmetic that is
already verified in the supplied scripts and says nothing about
`fs_mr_inference()`. The comparison must be against values the package returned.

**Asserting at the leaf.** Reading `.fs_mr_assemble()` directly and checking
`sigma_D` proves the assembler works. It does not prove `fs_mr_inference()`
calls it, or that the value reaches `debiased$se`. At least one assertion per
level must run the top-level function and read its returned object.

**A control that passes because the quantity is easy.** S1 at large separation
has $\mathrm{bias}_{\mathrm{sel}} \approx 0$, and so does a pipeline that
computes nothing. Every scenario needs a value that is distinguishable from
zero, from $\sigma_D$, and from what the naive estimate would give.
