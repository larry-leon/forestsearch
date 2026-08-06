---
title: "Task: F4 -- the propensity-score effect gap"
subtitle: "For CC. Branch mr-in-replicates, after 411a448, 0302e8c, 1856cba, 917e15d."
date: 2026-08-05
bibliography: []
---

## This is not a fix task

**Do not change any behaviour.** F2/F3 had a manuscript answer that told us
which side to move. F4 does not yet, and the plausible resolutions differ in
scope rather than in implementation detail --- one of them makes a whole
configuration unsupported. Producing evidence and stopping is the deliverable.

The one exception is at the end, and it is documentation only.

## What has been established from source

Verified by reading, not inferred from the audit summary. Line numbers are
current as of `917e15d`.

`R/forestsearch_main.R:1892-1893` resolves the default:
`ps_method <- if (is.RCT) "none" else "grf"`. So a non-randomised analysis gets
propensity adjustment **by default**; this is not an exotic path.

`:1903` sets `ps_adjust_resolved`, and `:1966-1973` **rebuilds the estimator
closure** with `ps_adjust_method = ps_adjust_resolved` for GLM outcomes. From
that point the identifier ranks candidates on the IPTW or g-computation
adjusted effect.

`.consistency_glm_pieces()` (`R/consistency_resample.R:221`) has formals `df`,
`outcome_type`, `effect_measure`, `treat.name`, `outcome.name`, `offset.name`,
`adjust_covariates`, `adverse_outcome`. **There is no `ps_method` or
`ps_adjust_method` among them**, and no `ps_hat`, `sw` or `ips_covar` is read.
It fits the plain conditional model.

`:1913` makes `adjust_covariates` mutually exclusive with
`ps_method != "none"`, so in this branch that formal is necessarily `NULL`.
`.consistency_glm_pieces()` therefore fits the **entirely unadjusted** model,
not merely a differently-adjusted one. Confirm this; it is the sharpest form of
the gap and it follows from the guard rather than from anything I ran.

`.fs_mr_pieces()` (`R/fs_mr_inference.R:61-73`) is a thin dispatcher to that
function and passes no propensity information either, because it has none.
So under `ps_method != "none"` the entire influence apparatus --- every
$\hat\beta(g)$, every $\mathrm{db}_{g,i}$, every $\sigma_{D,g}$, and hence $P$,
both bias terms and the infinitesimal-jackknife variance --- is computed for a
functional the identifier did not rank on.

Finally, the docstring at `R/consistency_resample.R:200-202` states that
propensity-adjusted (`iptw` / `dr_gcomp`) effects "return `NULL` (the caller
falls back to splitting)". The body's only such guard is
`if (effect_measure == "IRD") return(NULL)`. **The claim is not merely false but
unfulfillable**: the function has no argument from which it could learn it is in
the propensity-adjusted case. That is a defect independent of how the rest is
resolved.

## What to establish

Four questions. Answer each with a run, not an argument.

**1. Does MR execute at all under `ps_method != "none"`?** Run a GLM outcome
with `is.RCT = FALSE` (so the default `"grf"` applies) and `mr_inference = TRUE`.
Report whether `fs_mr_inference()` returns the full variant or the short one,
`n_family`, and whether `.validate_mr_configuration()` fires. It does not appear
to inspect `ps_method`; confirm rather than assume.

**2. How far apart are the two effects?** For the selected subgroup, report side
by side: the identifier's effect (from the rebuilt closure, PS-adjusted) and
`.consistency_glm_pieces()`'s `beta_hat` (unadjusted), plus `sigma_D`. Do this
for both `ps_adjust_method` values. A large gap makes this urgent; a small one
on this data does not make it safe, so report the magnitude without editorialising.

**3. Does the ranking itself change?** Compute the full candidate family's
effects under both definitions and report Spearman correlation, whether the
argmax coincides, and how many candidates change admission status against
`c_screen`. This is the quantity that decides whether F4 is a variance-estimation
problem or a selection problem --- they need different fixes.

**4. Which paths reach here?** Enumerate every route by which
`ps_method != "none"` can coexist with the influence path: MR, the resample
consistency engine, the bootstrap, cross-validation. Enumerate exhaustively ---
a `grep | head` reported four MR enforcement sites when there were five.

## Report, do not resolve

Write findings into `dev/identifier-alignment/`, in the same shape as
`code_theory_audit.qmd`. Set out the three resolutions and what the evidence
says about each, without recommending one:

**(a) Extend the influence path.** Give `.consistency_glm_pieces()` the
propensity information and return the adjusted influence. Note in the write-up
that IPTW has a standard weighted-GLM influence function, whereas `dr_gcomp`
additionally requires the propensity model's own estimation uncertainty --- the
two are not equally easy and may not both be feasible.

**(b) Make the docstring true.** Thread `ps_adjust_method` down solely so the
function can return `NULL`, letting the caller fall back to splitting and MR
decline. Report what MR then returns and whether that failure is visible to a
user or silent.

**(c) Hard-error.** Refuse `mr_inference = TRUE` together with
`ps_method != "none"` in `.validate_mr_configuration()`, as the three
misaligned `sg_focus` configurations already are.

State plainly which of these the evidence rules out, if any. The choice among
the survivors is Larry's and turns on whether the manuscript's scope covers
non-randomised data at all --- both applications are randomised trials, so it may
not.

## The one permitted change

The docstring at `R/consistency_resample.R:200-202` may be corrected to describe
what the function does: it returns `NULL` for `"IRD"`, and it has no propensity
argument, so a propensity-adjusted effect is fitted **as if unadjusted** rather
than declined. Do not add the guard --- only stop the documentation asserting a
behaviour the function cannot have.

If you would rather leave even that until the resolution is picked, say so and
leave it. Either is defensible; silently keeping the false claim is not.

## Scope

Out of scope: F9, F10, F11, F13, F14; the fixture; any analysis document;
`.fs_mr_ij_var()` and the bias denominators.

If the documentation change is made, `R CMD check --as-cran` clean is the bar and
stage your own hunks only --- `NEWS.md` and `R/forestsearch_main.R` still carry
uncommitted `max_subgroups_search` work. A findings-only session needs no commit
beyond the report.
