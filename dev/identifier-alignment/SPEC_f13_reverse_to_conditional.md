---
title: "Task: reverse F13 -- conditional repair, not Eq. 9 literal"
subtitle: "For CC. Supersedes SPEC_f13_eq9_alignment.md."
date: 2026-08-05
bibliography: []
---

## What changed, and why

`SPEC_f13_eq9_alignment.md` gave you a two-way choice. There were three, and the
one it omitted is the one Larry wants. **Your implementation is correct against
the spec you were given; the spec was wrong.** Nothing you did needs
re-diagnosing --- the defect you fixed is real and your centering analysis
carries straight over.

The defect was never the denominator on its own. It was that the two bias terms
carried **different** denominators: `selection_bias` over winner draws,
`fixed_bias` over all $B$. The residual then mixed differently normalised
quantities and `mean(r) != 0`. There are two coherent repairs, not one:

| | `bias_sel` | `bias_fix` | IJ draws | `mean(r)` |
|---|---|---|---|---|
| shipped (broken) | winners | all $B$ | winners | $-0.169$ |
| **(iii) Eq. 9 literal** --- what you built | all $B$, $D := 0$ | all $B$ | all $B$ | $+2.3\times10^{-17}$ |
| **(iv) conditional** --- what to build | winners | **winners** | winners | $+9.7\times10^{-18}$ |

Both are centered to machine precision, so centering does not choose between
them. They differ in **estimand**.

Larry's decision is (iv), conditional on identification. The reported analysis
exists only because a subgroup was identified, and $\beta(\widehat H)$ is
already a conditional estimand; estimating the correction on draws where no
selection occurred does not mirror what is reported. Three further reasons, in
his judgement and recorded so they are not relitigated: Theorem 1 has MR
reproducing the full bootstrap on a fixed family, and a bootstrap replicate that
identifies nothing contributes nothing to FB, so (iii) makes the two diverge
exactly on the no-winner draws; the variance is already conservative, and (iii)
widened it a further 74% on your measurement; and (iii) forces an arbitrary
convention on the complement, since an empty $\widehat H^*_b$ has *everyone* as
its complement, so $D_{\widehat H^{*c}}(b)$ is the ITT perturbation rather than
zero or undefined.

That last point is worth seeing clearly, because it is the internal
inconsistency in (iii) rather than a preference: (iii) sets the primary term to
zero on a no-winner draw while excluding the same draw from the complement. Under
(iv) the question does not arise --- the draw is excluded on both sides and the
complement is always a genuine complement.

## The change

Revert `6f97d7c`'s successor --- the F13 commit --- if it has landed; otherwise
revert the working-tree edits. Then apply (iv), which is **one line** from the
original shipped code:

- `selection_bias` --- unchanged from the original: `mean(sel_bias, na.rm = TRUE)`.
- `fixed_bias` --- **the only change**: average over the winner draws, not all
  $B$. `mean(P[sel, ok])` where `ok <- which(is.finite(sel_bias))`.
- The residual and the `.fs_mr_ij_var()` call --- unchanged from the original,
  both on `ok`.
- `sel_bias0` --- delete. Not needed under (iv).
- `selection_rate` --- unchanged.

Complement block: the same principle, and your two-NA-source distinction was
right and should survive. Both complement bias terms already need a shared
denominator; under (iv) that denominator is the draws where the complement's
perturbation exists --- your `use_c`. Keep it. What changes is only that
`fixed_c` and `selbias_c` no longer need the `sel_bias0` substitution, since
no-winner draws are excluded rather than zeroed. Read the block and confirm this
reduces to a smaller edit than the one you made; say so if it does not.

Keep `.fs_mr_ij_var()` untouched.

## The contract test

`tests/testthat/test-mr-eq9-residual-centering.R` should be **kept and
adjusted**, not deleted. Its core assertion --- $|\overline{r}| < 10^{-10}$ over
the draws the IJ uses --- holds under (iv) and is exactly the assertion that
would have caught the original defect. Adjust the draw set to `ok`, keep the
check that the *old* convention fails it, and rename the file to drop the `eq9`
token since (iv) is not Eq. 9 literal. Retain the equality of the coded
$\widehat V$ with Wager's centered form; that also holds under (iv), since
$\overline r = 0$ is what it needs.

## Measurement

Re-run the same two configurations. Expected, and worth checking as exact
identities:

- `selection_bias` and the de-biased point estimate: **unchanged** from the
  original shipped values, to every digit. (iv) does not move them.
- `fixed_bias`: moves, since its denominator changed. It has expectation zero,
  so the move should be small relative to `selection_bias`.
- `tilde_V`, `se_ij`: unchanged or nearly so. On the synthetic reference (iv)
  gave $+0.0\%$ against the shipped value.
- `mean(r)`: from $-8.55\times10^{-2}$ to zero at machine precision.
- The rate-1 control: identical on every field, before and after, as before.

If the point estimate moves at all under (iv), stop and report --- that would
mean something other than `fixed_bias` changed.

## Also record

Add to the findings document, or a short note in
`dev/identifier-alignment/`, the open question this surfaced: **what does the
full bootstrap record for a replicate that identifies no subgroup?** If it drops
the replicate, FB conditions on identification and (iv) aligns MR with it, which
strengthens the choice. If it substitutes something, that is worth knowing. Do
not chase it now --- record it.

Also record, for the manuscript file rather than for action: Eq. 9 divides by
$B$, which appears to assume every draw admits a candidate. If so the manuscript
is **silent** on this case rather than contrary to (iv), and no amendment is
implied. Flag it as unverified --- it needs a read of §3.1, which is out of scope
here.

## Scope

In scope: the revert, the one-line change plus its complement counterpart, the
test adjustment, the measurement, a corrected `NEWS.md` entry.

Out of scope: F4, F9, F10, F11, F14; the fixture; any analysis document.

The `NEWS.md` entry must now say something different from the one you wrote.
Under (iv) point estimates and intervals are **materially unchanged**; what is
fixed is an incoherence between the two bias terms' denominators, which affected
the residual and hence the finite-$B$ correction wherever `selection_rate < 1`.
That is a much smaller user-facing claim than the reverted version made --- do
not carry the 74% figure over.

`R CMD check --as-cran` clean is the bar. Stage your own hunks only ---
`NEWS.md` and `R/forestsearch_main.R` still carry uncommitted
`max_subgroups_search` work.
