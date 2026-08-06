---
title: "Task: F12 -- the cross-validation cache leak"
subtitle: "For CC. Branch mr-in-replicates, after 411a448 and 0302e8c."
date: 2026-08-05
bibliography: []
---

## Before anything

One shell was left running at the end of the last session. Identify it and kill
it before starting. A stray process plus a buffered background check is how a
clean `R CMD check` gets attributed to the wrong tree state, which nearly
happened once already.

## The decision

**Cached fits are not reusable across folds. Null them.** This is settled --- do
not re-open it in the write-up.

The reasoning: cross-validation estimates the out-of-sample behaviour of the
procedure. A candidate family estimated on all the data, including the held-out
fold, is not the procedure that ships, so `sensCV`, `ppvCV` and exact-match
would be optimistic by construction. The bootstrap already made this call. CV
nulls `ps_hat` on identical reasoning and leaves the other four. **The asymmetry
is the defect, not the nulling.**

Cost is the only argument the other way, and it is a performance concern rather
than a correctness one. If refitting a forest per fold proves painful, the
honest form is an explicit opt-in argument defaulting to off and warning when
set --- not silent reuse. Do not add that argument as part of this task.

## The change

`R/forestsearch_cross_validation.R`, both entry points. The `ps_hat` nulling
sits at `:396` (`forestsearch_Kfold()`) and `:867` (`forestsearch_tenfold()`);
the four additions go alongside it, in both.

- `cv_args$grf_res <- NULL`
- `cv_args$grf_cuts <- NULL`
- `cv_args$dina_res <- NULL`
- `cv_args$dina_cuts <- NULL`

Mirror the bootstrap's comment rather than writing a new one. Its block at
`R/bootstrap_analysis_dofuture.R:573-582` is labelled `CATEGORY 2: VARIABLE
RE-SELECTION` and already states the reason for the DINA half in full. Reuse
that language so the two sites read as one decision, and say explicitly in the
CV comment that this matches the bootstrap --- the absence of such a note is
part of what made the gap invisible.

Confirm while you are there that `args_call_all`'s bulk `mget()` is in fact what
carries these four into `cv_args`, rather than assuming it from the bootstrap's
behaviour. If they arrive by some other route the nulling may need to go
elsewhere.

## Measurement

The leak only bites when a cached fit is supplied, so a run at
`grf_res = NULL` proves nothing. Construct the failing case deliberately.

Fit once with `subgroup_method = "grf"`, keep the returned `grf_res`, then run
CV twice on the same data and seed --- once passing the cached fit, once not ---
before and after the change. Report `sensCV`, `ppvCV` and exact-match for all
four cells.

The expected pattern: before the change, the cached run is optimistic relative
to the uncached one; after, the two agree. If they already agree before the
change, the cached fit is not reaching the folds and the finding needs
re-diagnosing rather than the fix declaring victory. Say so if that happens.

Repeat for `subgroup_method = "dina"` with a cached `dina_res`. Keep `Kfolds`
small; this is a wiring measurement, not an operating-characteristics run.

## Bundled: the shipped-config DINA re-run

Cheap, and it settles a gap in the audit record.

CC's F2/F3 before/after table ran `subgroup_method = "dina"` at
`hr.threshold = 1.25`. The analysis document
`analysis_actg175_binary_multimethod.qmd` runs that mode at
`fs_dina_or_threshold <- 1.0`, and its `fs_dina` call passes no
`max_subgroups_search`. So the table describes a configuration the document does
not use, on two axes at once.

Re-run the `fs_dina` call exactly as the .qmd has it --- `or_threshold = 1.0`,
`mr_inference = TRUE`, `draws = 5000` --- before and after `0302e8c`, and report
`found`, the selected label, `n_family`, `selection_rate`, `selection_bias` and
the naive estimate. No bootstrap needed.

Prediction on record, so the result is informative either way: **nothing
changes.** At a floor of `log(1.0) = 0` the beta-hat admission is close to
vacuous, so the F2/F3 fix should be inert at the shipped threshold. If it is
inert, say so plainly --- an alignment fix that binds at 1.25 and not at 1.0 is
a fact about the threshold worth recording, not a weak result. If something does
change, that is more interesting and needs a look before regeneration.

Do not regenerate any analysis document. Regeneration remains held.

## Scope

In scope: the four nullings at two sites, their comment, the measurement, the
DINA re-run.

Out of scope: F9, F10, F11, F13, F14; the opt-in reuse argument; DINA's
tau-bar/beta-hat ranking behaviour; anything touching `.fs_mr_ij_var()` or the
bias denominators, since the F13 convention is undecided.

`R CMD check --as-cran` clean is the bar. `NEWS.md` and
`R/forestsearch_main.R` still carry uncommitted `max_subgroups_search` work ---
stage your own hunks only, as you did for `0302e8c`. A `NEWS.md` entry is
warranted here: this changes CV metrics for anyone passing a cached fit.
