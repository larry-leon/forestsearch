---
title: "Task: F9, F10, F11 --- closing the audit"
subtitle: "For CC. Branch mr-in-replicates, after c88e1c2."
date: 2026-08-06
bibliography: []
---

## Scope

The last three open findings from `code_theory_audit.qmd`. All three are low
severity and none changes any estimate. **F14 is already settled as
not-a-defect** and needs nothing; note it as closed in the audit record and move
on.

This is a sweep, not an investigation. If any of the three turns out to be
larger than described, stop and report rather than widening the work --- the
findings were rated low on a reading, and a reading can be wrong.

**Three commits, one per finding.** They are unrelated to each other and a
combined diff would be harder to read than three small ones.

Line references below are from a pre-`dad0415` snapshot and several files have
moved since. **Read the current lines before editing.**

## F9 --- the dispatch assert reads three of five sites

`.assert_sg_focus_dispatch_complete()` (`R/forestsearch_helpers.R:2063`) exists
to catch a `sg_focus` value gaining a dispatch site that some other site does not
handle. The audit found it reads three of five such sites: the two it misses are
`.dina_reselect_on_effect()`, where `maxeff`/`maxeffCons` route through a
default branch, and `subgroup.consistency()`, which re-inlines the canonical
vector as a literal.

**Behaviour agrees today.** This is a guard that would miss a future divergence,
not a current defect --- which is exactly why it is worth fixing: a check that
appears to cover five sites and covers three is worse than one that covers three
and says so.

Note that `selection_criteria.qmd` now documents the DINA fall-through as
deliberate: `maxeff` and `maxeffCons` are synonyms outside the consistency
engine precisely *because* neither matches a branch there. So the fix is to make
the assert read that site and confirm the fall-through, **not** to add branches
that would break the documented synonymy. Check `@sec-vocab` in that document
before editing, and if your change contradicts it, one of the two is wrong ---
say which rather than changing both to agree.

For `subgroup.consistency()`'s re-inlined literal: the canonical vector should
have one definition. `forestsearch_main.R:1329` already carries a comment saying
not to re-inline it, so the second site is a violation of a stated rule rather
than an oversight.

Enumerate the dispatch sites exhaustively --- unbounded `grep -rn`, not
`| head`. A bounded grep reported four MR enforcement sites when there were
five.

## F10 --- dead assignments

`.mr_cscr` (`R/forestsearch_main.R:2105` and `:2282`), and `c_screen_mr` /
`c_consistency_mr` (`:3148-3149` and `:3156-3157`). All assigned, none
consulted. They are leftovers of the reconstruction that
`.fs_resolve_admission()` removed.

Delete them. Confirm by unbounded grep that each is genuinely unread across all
of `R/` **and** `tests/` before deleting --- a test asserting on one of these
would make it live, and deleting it would be the same "reconstructed at a second
site" defect in reverse.

If any is read anywhere, stop and report. That would mean the resolver's
consolidation is incomplete, which is a finding rather than a cleanup.

## F11 --- an inert formal

`max.minutes` is a formal of `subgroup.search()` (`R/subgroup_search.R:77`) and
of `search_combinations_parallel()`, threaded through to `:150` and `:234`, and
no path consults it. `forestsearch()`'s roxygen already documents it as inert
(`R/forestsearch_main.R:610`, "Currently inert; scheduled for..."), but
`subgroup.search()`'s does not (`R/subgroup_search.R:25`, "Maximum minutes for
search").

**Do not remove the formal.** It is part of a public signature and the
`forestsearch()` roxygen says it is scheduled rather than abandoned. The fix is
that `subgroup.search()`'s roxygen should say what `forestsearch()`'s says ---
that it is currently inert --- so the two agree.

Check whether `search_combinations_parallel()` also documents it, and whether
`run_simulation_analysis.R:63` passes a value that a reader would take as
meaningful. If a caller supplies `max.minutes = 5` and nothing consumes it, that
call is misleading even though it is harmless; a comment at the call site is
enough.

This is the same defect class as today's `conf_force` work: a parameter
description that promises a behaviour the code does not have. That one was worse
because it pointed away from the argument that worked, but the shape is
identical.

## Then close the audit record

`dev/identifier-alignment/code_theory_audit.qmd` carries all fourteen findings
with their status. Update it so each of the fourteen shows its final state ---
fixed, settled, or deferred --- with the commit that closed it where one exists.

F4's Q2 and Q3 remain **deferred by decision**, not open by oversight. Record
them that way and note that Q3 is the load-bearing one: it decides whether F4 is
a variance-estimation problem or a selection problem, which need different
fixes.

Do not re-litigate any finding while updating the record. If the current state
contradicts what a finding says, that is a report, not an edit.

## Scope

In scope: F9, F10, F11, and the audit-record update. Four commits.

Out of scope: F4 Q2/Q3; the `forestsearch_Kfold` silent-failure hardening; any
analysis document; the GLM/continuous simulation work, which is running in a
separate workstream; anything in the manuscript.

`R CMD check --as-cran` clean is the bar. No `NEWS.md` entry is warranted for
F10 --- deleting unread variables is invisible to users. F9 and F11 are
judgement calls: F11 corrects user-facing documentation, so probably yes; F9 is
internal, so probably not. Say which you took and why.
