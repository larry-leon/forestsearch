---
title: "Spec: selection_criteria.qmd --- engine scope and the missing keys"
subtitle: "For CC. Branch mr-in-replicates, after 8694e3d."
date: 2026-08-06
bibliography: []
---

## What this document is, and what is wrong with it

`quarto/methodology/selection_criteria.qmd` is a careful, package-validated
account of `sg_focus` and `selection_rule`. Nothing in it needs rewriting. The
problem is narrower and specific: **it describes the consistency engine and
reads as though it describes all three.**

Measured on the current file: DINA appears once, GRF twice, MR once, in 677
lines. Yet "The identifier is the estimator" opens by stating the identifier has
two parts, engine and focus, and tables three `subgroup_method` values. A reader
arrives believing @sec-foci and @sec-inert cover all three engines. They cover
one.

That is the same confusion that cost real time in the audit --- a Pareto
frontier of 8 read against an MR family of 1, a "family of one" read as a
property of DINA rather than of a threshold. The fix is not more content but a
declared axis.

**Do not restructure the document.** Every change below is additive or a
scoped qualification. The existing prose, worked examples and validation block
stay as they are.

## Two corrections --- these are wrong as read, not merely incomplete

### C1. The scope statement

Add to `# Purpose`, before "The identifier is the estimator", a statement that
the document's sort keys, band behaviour and inert-argument table describe
`subgroup_method = "consistency"`. The other two engines rank on their own
statistics and are covered in the new section required by A1 below.

Keep the engine/focus table where it is --- it is correct --- but make the
sentence after it say that what follows develops the focus axis for one engine,
not all three.

### C2. `@sec-inert` needs the engine axis

Its `hr.threshold` row reads `all but maxeff` / `disabled`. That is true of the
consistency engine and **false of DINA and GRF as of `0302e8c`**: both now admit
on $\widehat\beta(g) \ge t_g$ with $t_g = c_{\mathrm{screen}} =
\log(\texttt{hr.threshold})$, read from the resolved admission set. The
threshold is active for them in a way it was not when the table was written.

Two options; pick one and say which in your report. Either add an `engine`
column so each row states where it holds, or keep the table consistency-only
under an explicit heading and add a second table for DINA and GRF. I lean to the
second --- the rows differ enough that one table would carry mostly caveats ---
but it is a presentation call and you can see the material.

Either way `@sec-inert`'s closing line, that arguments "degrade silently, which
is why a run's resolved settings are worth recording", should now point at
`admission_resolved` as the object that records them.

## Five additions

### A1. The DINA and GRF sort keys

A new section covering how the other two engines rank. The essential facts, all
verifiable against the source rather than to be taken from here:

Neither engine has a consistency screen, so **no key contains a `Pcons`
term**. Each ranks natively first --- DINA on mean $\bar\tau$, GRF on the DR
score --- and then, when the effect-path re-selection applies,
re-ranks on $\widehat\beta$ among the admitted set. `.dina_reselect_on_effect()`
and `.grf_reselect_on_effect()` in `R/forestsearch_helpers.R` are the two sites;
read them rather than describing them from this spec.

State plainly that an empty admitted set returns **no subgroup** rather than
falling back to the native winner, and why: returning the native winner would
report a subgroup selected on a statistic the analysis does not report.

### A2. `maxeff` and `maxeffCons` are synonyms outside the consistency engine

`@sec-foci` distinguishes them by whether the ranking is gated on consistency
qualification. DINA and GRF have no consistency term, so the distinction
collapses --- the two values name the same rule there. Say so, and say it in
`@sec-vocab` as well, which is where a reader checks whether two names mean one
thing.

### A3. MR reuses `maxeff` and `maxcons` for different rules

This is the addition most worth getting right, because it is a name collision
across two layers of the same package. `fs_mr_inference()`'s `reselection`
argument takes values spelled the same as two `sg_focus` values and does not
mean the same thing: it governs how a *perturbed* family is re-ranked on each
draw, not how the observed candidates are ranked once.

Give it its own short section, state the two meanings side by side, and say
which argument carries which. A reader who assumes one meaning transfers to the
other will misread both `sg_focus = "maxeff"` runs and MR's `reselection`
settings.

### A4. The band, per engine

`.compute_inclusion_band()` is now the only band on any selection path. State
where it applies per engine, and that the GRF effect path takes its floor from
the resolved admission set rather than the hard-coded `dmin = 1` it previously
used. One paragraph; the detail lives in the audit record, not here.

### A5. The `conf_force` freeze recipe

The document does not mention `conf_force`, `cut_type` or
`conf.cont_medians` at all. Under `cut_type = "default"` the cut locations are
derived from observed quantiles and therefore move with a resample, which is why
`.fs_family_status()` reports `"no-front-end"` rather than `"fixed"`. Give the
recipe for genuinely fixed cuts --- explicit `conf.cont_medians` /
`conf.cont_jcuts`, or `conf_force` with literal cut expressions --- and state
which configuration each produces.

Do **not** editorialise about whether fixed cuts are correct. That is settled and
out of scope here: the full bootstrap refits everything per resample by design.

## Validation is not optional

`@sec-validate` currently checks every substantive claim in the document against
the shipped package --- the dominance flags, both bands, the alias normaliser,
and each focus's pick against `sort_subgroups()`. That is why the document can be
trusted.

**Every new claim above that is checkable must gain a check there**, in the same
`.ok(label, cond)` form. At minimum:

- the DINA and GRF keys contain no `Pcons` term
- `maxeff` and `maxeffCons` resolve to the same ranking outside the consistency
  engine
- the GRF effect path's floor equals the resolved `effect_floor`, not `1`
- an empty admitted set returns `found = FALSE` / `sg_def = NULL`

If a claim cannot be checked without running a full `forestsearch()` fit, say so
in the text rather than asserting it unchecked. A document that mixes validated
and unvalidated assertions in one voice is worse than one that marks the
difference --- an unchecked claim in a commit message is exactly how `84203ed`'s
false statement survived.

## Scope

In scope: `quarto/methodology/selection_criteria.qmd` and its rendered output.

Out of scope: `pareto_frontier_explainer.qmd`; any package code; any analysis
document; the deferred audit findings.

Render must be clean and every `.ok` line must PASS. Report any that do not
rather than adjusting the claim to match --- a failing check is a finding.

`R CMD check` is not the bar here since no package code changes, but confirm
the document renders against the installed build at HEAD rather than a stale
one; a stale install has already caused trouble once in this project.
