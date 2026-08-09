---
title: "forestsearch 0.2.0 --- release status"
subtitle: "What is settled, what is pending, what a submission still needs"
date: 2026-08-06
bibliography: []
---

## How to read this

Two kinds of entry. **Settled** items are recorded facts with the evidence
attached. **Pending** items are headings with nothing under them yet --- they
are listed so the gap is visible, not because anyone has looked.

Nothing here recommends a submission date. That decision depends on the GLM
pathway scrutiny in §4, which has not reported.

## 1. Branch position {#sec-branch}

```
feature/mr-in-replicates   <- all audit work; pushed, R CMD check clean
        |
        v  merge, then check the MERGED state
feature/glm-extension
        |
        v  separate decision
master                     <- CRAN target
```

**Settled.** The audit cycle landed entirely on `mr-in-replicates`, pushed and
`R CMD check --as-cran` clean at 0/0/0 on the final state of each change.

**Do next, and it is independent of everything else here:** merge to
`glm-extension` and run `R CMD check --as-cran` **on the merged result**. Every
check in the audit cycle validated `mr-in-replicates` in isolation. If
`glm-extension` has moved since the branch point, the merge produces a state
neither branch has been checked in. Conflicts, if any, will most likely be in
`NEWS.md` and `R/forestsearch_main.R`; both took heavy edits.

The merge to `master` is a separate decision and should not ride along with
this one.

## 2. What the audit cycle closed {#sec-audit}

**Settled.** Fourteen findings in `dev/identifier-alignment/code_theory_audit.qmd`.
Its `@sec-final-status` table is authoritative; the summary here is for release
purposes.

| | count | which |
|---|---|---|
| fixed | 8 | F1, F2, F3, F9, F10, F11, F12, F13 |
| settled, not a defect | 1 | F14 |
| partly fixed | 1 | F4 --- Q1/Q4 done, Q2/Q3 deferred by decision |
| open, judged not worth fixing | 4 | F5, F6, F7, F8 |

The eight fixed are the ones that produced incorrect results. In outline: an
entire outcome type silently producing no correction; identifiers selecting on
statistics the analysis does not report; held-out data leaking into training
folds; propensity scores attached to the wrong subjects; and a correction
computed against a family of one or two where the real family held hundreds.

**The four open findings were judged, not skipped.** F6 and F8 are bookkeeping
--- a truncation not recorded, display sites re-implementing a band --- and
neither changes an estimate. F5 is a scope question rather than a validity one:
multiplier resampling is not required to mimic the full bootstrap, so its family
being a superset is the design rather than a defect. F7 is a degenerate case
where no candidate in the family shows harm, and the three consumers that
diverge all arrive at "no subgroup" by different routes.

None is user-facing. No release note is warranted for any of them.

## 3. Breaking changes in this cycle {#sec-breaking}

**Settled.** Three, all in `NEWS.md`, all deliberate. A first CRAN release can
carry them; they should read as decisions rather than accidents.

| change | effect |
|---|---|
| `max_subgroups_search` defaults to `Inf` | was 200. Truncation applied a preview ordering rather than the selection criterion, so the candidate optimal under the actual criterion could be discarded before evaluation |
| `ps_method != "none"` errors for survival | was a structural no-op --- the score was estimated, attached, and read by no Cox path. Note the default resolves to `"grf"` when `is.RCT = FALSE`, so observational survival runs that never mentioned the argument now error |
| MR estimates move wherever `selection_rate < 1` | the two bias terms now share one denominator. Runs at `selection_rate = 1` are unaffected to every digit |

A fourth is user-visible without being breaking: **any `effect_measure = "MD"`
analysis run with `mr_inference = TRUE` before this cycle has no MR result at
all** and should be re-run. The influence path errored on every `lm()` fit and
failed silently at every level.

## 4. GLM pathway scrutiny --- PENDING {#sec-glm}

A parallel workstream in **`dev/glm-continuous-sims/`** is designing simulations
for the GLM/continuous pathway against closed-form references. Two documents
there record what is already established:

* `HANDOFF_glm_continuous_simulations.md` --- what closed forms can and cannot
  confirm, in three layers, and the open design questions.
* `ADDENDUM_glm_pathway_context.md` --- what this cycle's audit found about the
  GLM pathway specifically, including two measurement traps that would recur in
  a simulation harness.

**Nothing to record here until it reports.** What it finds bears directly on the
release decision, because the GLM pathway is where most of this cycle's defects
were found and where the only closed-form assertions on an influence path now
live.

Fill in: what the simulations exercised, what they found, whether anything
found blocks a release.

## 5. CRAN policy checklist --- NOT ASSESSED {#sec-cran}

`R CMD check --as-cran` clean is necessary and not sufficient. **None of the
following has been looked at during this cycle**, and each is a plausible
submission blocker.

- [ ] `DESCRIPTION` --- Title and Description fields per CRAN's wording rules; Authors@R; version bumped to 0.2.0
- [ ] LICENSE --- present, and matching what `DESCRIPTION` declares
- [ ] Examples --- runtime within CRAN limits; anything slow wrapped in `\donttest{}` and confirmed to still run
- [ ] Vignettes --- build from a clean checkout on a machine without the development environment
- [ ] URLs --- all resolving; CRAN checks these and fails on redirects
- [ ] `NEWS.md` --- parses as CRAN expects, and the breaking changes in §3 read as intentional
- [ ] Package size and check timing within limits
- [ ] Reverse dependencies --- none expected for a first release, worth confirming
- [ ] A check on a second platform, given the Mac/Linux floating-point differences already documented in this project

## 6. What to do next {#sec-next}

In order, with the first two independent of everything pending.

1. **Merge to `glm-extension` and check the merged state** (§1).
2. **Work the CRAN policy checklist** (§5). It needs no decisions and is
   entirely mechanical --- there is no reason for it to wait on §4.
3. **Wait for the GLM pathway scrutiny** (§4) before deciding on the merge to
   `master`.
4. **Then decide.** Not before.

## 7. A note on scope

The wait for a "code will evolve" stopping point has a failure mode worth
naming: this cycle found eight defects that mattered, and they were found by
building assertions, not by waiting. A tagged 0.2.0 gives something to bisect
against. Several investigations in this cycle began by reconstructing which
package state produced a number --- one of them cost an hour and turned on a
payload whose provenance was ambiguous.

A perfect codebase is not the bar. A codebase whose state is *identifiable* is.
