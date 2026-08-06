---
title: "Task: the two `ps_hat` findings from F4 Q4"
subtitle: "For CC. Branch mr-in-replicates, after 8bb3ed8. Establish before fixing."
date: 2026-08-06
bibliography: []
---

## What these are

Two findings you surfaced while enumerating the paths that reach the influence
apparatus for F4, and recorded as "not chased". Neither has been measured. They
are grouped here because both concern `ps_hat` handling and both are cheap ---
**not** because they are the same defect.

Neither is propensity-theory work. F4's Q2 and Q3 remain deprioritised and are
out of scope.

**Establish each before fixing it.** F12 became trustworthy because the leak was
constructed and measured rather than argued from the code's shape; the same
standard applies here. If a measurement contradicts the finding as recorded,
that is the result --- report it and stop rather than fixing something that
isn't there.

## Finding A --- the bootstrap carries `ps_hat` into replicates unaligned

**Severity: this produces wrong values, not merely optimistic ones.** F12
inflated a metric; this attaches a propensity score to the wrong subject.

As recorded: `forestsearch()` validates a user-supplied `ps_hat` at
`R/forestsearch_main.R:1920-1948` --- length against `nrow(df.analysis)`, no
`NA`, strictly inside $(0,1)$ --- then assigns `df.analysis$ps_hat <- ps_hat`
positionally and derives the IPTW weights from it. A bootstrap resample has the
**same row count** but different subjects, so the length check at `:1922` passes
while every value belongs to a different person.

Cross-validation nulls `ps_hat` for exactly this reason, with a comment saying
so. The bootstrap's `CATEGORY 2: VARIABLE RE-SELECTION` block
(`R/bootstrap_analysis_dofuture.R:573-582`) nulls `grf_res`, `grf_cuts`,
`dina_res` and `dina_cuts`, and does not null `ps_hat`.

### Establish it first

Construct the failing case deliberately; a run at `ps_hat = NULL` proves
nothing. Fit once with `is.RCT = FALSE` and a supplied `ps_hat`, then run the
bootstrap and confirm, inside a replicate, that the `ps_hat` column attached to
the resampled frame does **not** correspond to the resampled subjects. The
cleanest evidence is a `ps_hat` vector constructed to be a deterministic
function of a covariate --- then a replicate's `ps_hat` should equal that
function of *its own* rows and will not.

Report what you find. Three outcomes are possible and they need different
responses: the values are misaligned as described; something upstream
re-derives `ps_hat` per replicate so the finding does not hold; or the path is
unreachable because a supplied `ps_hat` never survives into `args_call_all`.

### If established

Null `ps_hat` in the bootstrap's CATEGORY 2 block, mirroring CV's comment and
saying explicitly that it matches CV. Then confirm PS is re-estimated per
replicate rather than silently dropped --- nulling a *supplied* score is only
correct if the replicate estimates its own, which is what CV relies on.

Measure before and after on the same seed: the selected subgroup, the naive and
FB estimates. State the size of the change rather than only its direction.

**Do not add an opt-in reuse argument.** If someone wants to pass a precomputed
score they can, per replicate, through a route that does not exist yet; that is
a separate design question.

## Finding B --- survival plus `ps_method != "none"` may be a no-op

As recorded: the estimator-closure rebuild at
`R/forestsearch_main.R:1966` is gated on `if (is_glm)`, and no Cox path reads
`sw`, `ps_hat` or `ips_covar`. If that holds, propensity scores are estimated,
attached to the frame, and never consumed --- a user asking for adjustment on a
survival outcome gets none, with no message.

This is a **different** defect from A. A is wrong values; B is a silently
ignored argument. Do not fix them in one commit.

### Establish it first

Run a survival outcome with `is.RCT = FALSE`, so `ps_method` defaults to
`"grf"`, at both `ps_adjust_method` values and at `ps_method = "none"`. If B
holds, all three give **bit-identical** estimates and the same selected
subgroup. Report the comparison either way.

Then confirm the negative directly rather than by absence of evidence:
enumerate exhaustively --- unbounded `grep -rn`, not `| head` --- every read of
`sw`, `ps_hat` and `ips_covar` across `R/`, and state which are inside GLM-only
closures. A `grep | head` reported four MR enforcement sites when there were
five.

### If established

**Do not add PS support to the Cox path.** That is a feature, and it is
propensity-theory work you have been asked to leave alone.

The fix is to stop the argument being silently ignored. Two options; present
both with what the evidence supports and let Larry choose:

1. **Warn** at resolution time that `ps_method` has no effect for survival
   outcomes, naming the argument and the outcome type.
2. **Error**, on the grounds that a silently ignored adjustment argument is the
   same class of defect as a silently ignored `dmin.grf` --- which was F2 and was
   treated as a defect rather than a documentation matter.

Whichever is chosen, the roxygen for `ps_method` and `ps_adjust_method` must say
which outcome types they affect. Today's `conf_force` work is the precedent:
four parameter descriptions pointed away from the only argument that did the
job, and that is what made a wrong spec reasonable to write.

## Scope

In scope: establishing both findings, fixing A if established, and whichever of
B's two options Larry picks --- in **separate commits**.

Out of scope: F4 Q2/Q3; F9, F10, F11, F14; adding PS support to the Cox path;
any opt-in reuse argument; any analysis document; the GLM/continuous simulation
work, which is running in a separate workstream.

`R CMD check --as-cran` clean is the bar. A `NEWS.md` entry is warranted for A
--- it changes bootstrap results for anyone supplying `ps_hat` --- and for B
whichever option is taken.

Stop after establishing both if either measurement contradicts the finding.
