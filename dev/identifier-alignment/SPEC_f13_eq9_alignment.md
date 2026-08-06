---
title: "Task: F13 -- implement Eq. 9 as written"
subtitle: "For CC. Branch mr-in-replicates."
date: 2026-08-05
bibliography: []
---

> **SUPERSEDED — DO NOT FOLLOW.** Replaced by
> `SPEC_f13_reverse_to_conditional.md`. This spec framed the choice as two
> conventions when there were three; the omitted one (conditional on
> identification) is what was implemented, in `dad0415`. The defect this spec
> diagnoses is real and its centering analysis is correct — only the prescribed
> repair is wrong. Retained for the record.

## The decision

**Align with manuscript Eq. 9.** Both bias terms divide by $B$. This is settled;
do not re-open it.

Eq. 9 is:

$$
\mathrm{bias}_{\mathrm{sel}} = \frac{1}{B}\sum_{b=1}^{B} D_{\widehat H^*_b}(b),
\qquad
\mathrm{bias}_{\mathrm{fix}} = \frac{1}{B}\sum_{b=1}^{B} D_{\widehat H}(b).
$$

The code implements the second correctly and the first over the draws that
produced a winner (`mean(sel_bias, na.rm = TRUE)`).

## What "as written" requires -- read this before editing

Dividing by $B$ is **two** changes, not one, and making only the first is worse
than making neither.

Dividing by $B$ means a draw admitting no candidate contributes **zero** to the
sum rather than being removed from the average. That fixes what the residual is
on such a draw: $D_{\widehat H^*_b}(b) := 0$ there. $D_{\widehat H}(b)$ is always
defined, since $\widehat H$ is the observed selection. So the residual of Eq. 12
is defined on **every** draw, and the IJ of Eq. 13 must therefore be evaluated
over all $B$ draws rather than the winner subset.

Do both, and

$$
\bar r = (\mathrm{bias}_{\mathrm{sel}} + \mathrm{bias}_{\mathrm{fix}})
        - \mathrm{bias}_{\mathrm{sel}} - \mathrm{bias}_{\mathrm{fix}} = 0
$$

identically. Verified numerically at $-1.2\times10^{-16}$, with the coded
$\widehat V$ then matching Wager's centered $\widehat v$ to **exactly** zero
difference. Doing only the denominator, leaving the IJ on the winner subset,
made $\bar r$ worse: $-0.367$ against $-0.169$ as coded.

**Consequence worth stating in the write-up:** no amendment to Eq. 13 is needed.
Its uncentered $\overline{r^2}$ *is* Wager's centered quantity once Eq. 9 is
implemented literally. The apparent gap was an implementation artefact.

## The changes

All in `R/fs_mr_inference.R`. **Read the current lines before editing** --- the
numbers below are from a snapshot and the file has changed since.

Selected subgroup, around `:447-457`:

1. `selection_bias` --- denominator $B$, not the count of winner draws.
2. `sel_bias` --- substitute `0` for `NA` when forming the residual, so `r_H` is
   defined on every draw. Keep the `NA`-carrying vector for `selection_rate`.
3. The `.fs_mr_ij_var()` call --- pass all $B$ draw indices, not
   `which(is.finite(sel_bias))`.

Complement, around `:492-505`: the same three changes. Read that block rather
than assuming it mirrors the first --- it handles an unfit complement separately
(`vals[is.na(bh_c[...])] <- NA_real_`), and a complement that could not be *fit*
is not the same event as a draw that admitted *nobody*. Those must stay
distinct. Say so if the distinction makes any of the three changes inapplicable
there.

Leave `.fs_mr_ij_var()` itself alone. Its `ok` argument becomes the full index
set at the call site; the function needs no change, and its `length(ok) < 2L`
guard still stands.

`selection_rate` keeps its current meaning --- the proportion of draws admitting
a candidate --- and stays in the returned object. It is now the diagnostic that
tells a user whether this change moved anything at all.

## Measurement

The two conventions coincide exactly at `selection_rate == 1`, so a run at
rate 1 demonstrates nothing. Construct a configuration where the floor binds on
a substantial minority of draws and report, before and after:

`selection_rate`, `selection_bias`, `fixed_bias`, `debiased$est`, `var_ij`,
`tilde_V`, `se_ij`, `ij_source`.

Expected: `selection_bias` falls by exactly the factor `selection_rate` ---
check that identity, it is exact and a clean wiring test. The point estimate
moves by `selection_bias_old * (1 - selection_rate)`. `tilde_V` will also move,
since the IJ now runs over more draws with a different residual vector; in one
synthetic configuration at $p = 0.82$ it rose 29%. Report the direction and size
rather than predicting them.

Also run a rate-1 configuration and confirm **nothing** changes, to any digit.
That bounds the blast radius, and both the shipped DINA config and GRF on
ACTG175 ran at rate 1.

Add a contract test asserting $|\bar r| < 10^{-10}$ over the full draw set. This
is a static property and exactly the standing exception to the no-testthat rule
--- it is the assertion that would have caught this.

## Reference

`dev/identifier-alignment/verification/` holds the closed-form scripts. The new
`verify_eq9_alignment.R` accompanying this spec contains the three-arm
comparison the numbers above come from. File it alongside the other two and note
it in that README; it is reference material, not a package test.

## Scope

In scope: the three changes at both blocks, the contract test, the measurement,
a `NEWS.md` entry.

Out of scope: F4, F9, F10, F11, F14; the fixture; any analysis document;
`.fs_mr_ij_var()`'s body; the manuscript.

**This changes reported MR point estimates and intervals wherever
`selection_rate < 1`.** The `NEWS.md` entry must say so plainly and name
`selection_rate` as the field that tells a user whether they are affected.
Regeneration is still held, but this is now a reason it will be needed.

`R CMD check --as-cran` clean is the bar. Stage your own hunks only --- `NEWS.md`
and `R/forestsearch_main.R` still carry uncommitted `max_subgroups_search` work.
