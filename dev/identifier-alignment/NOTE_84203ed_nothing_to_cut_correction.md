---
title: "Correction: 84203ed's \"nothing to cut\" claim is false for GBSG"
subtitle: "Measured across six package states; the restriction was binding on both engines"
date: 2026-08-06
bibliography: []
---

## The claim being corrected

`84203ed` ("Remove the MR family restriction; MR uses the identifier's
qualifying family") removed `.fs_mr_restrict_native()`, which had narrowed MR's
re-selection family to candidates within a multiplicative band of the best
NATIVE statistic. Its commit message reports the size of what was excluded:

```
    GRF    2 -> 1289     (0.16% retained)
    DINA   1 -> 1        (not binding here)
```

and then states:

> For DINA on ACTG175 -- and for both engines on GBSG -- the qualifying
> frontier holds a single candidate, so there was nothing to cut.

**The GBSG half of that sentence is false.** The restriction was binding on
both GBSG engines, and more severely than on the ACTG175 GRF case the commit
did report.

## What was measured

`n_family`, `selection_bias` and `selection_rate` for the two GBSG fits that
`analysis_gbsg_survival_multimethod.qmd` performs -- `fs_dina` (`:1516`) and
`fs_grf` (`:1757`) -- reproduced call-for-call in a standalone script, run
against six installed package states via `git worktree` + `R CMD INSTALL`.

| package state | commit | when | DINA `n_family` | DINA bias | DINA rate | GRF `n_family` | GRF bias | GRF rate |
|---|---|---|---|---|---|---|---|---|
| baseline render | `daa1e0a` | 08-05 10:30 | 1 | 0.2632 | 0.4684 | 1 | 0.1771 | 0.6796 |
| | `e38acf1` | 08-05 11:11 | 1 | 0.2632 | 0.4684 | 1 | 0.1771 | 0.6796 |
| | `8d60b1e` | 08-05 12:03 | 1 | 0.0272 | 0.9688 | 1 | 0.0604 | 0.9186 |
| `84203ed^` | `c68fe70` | 08-05 14:08 | 1 | 0.0272 | 0.9688 | 1 | 0.0604 | 0.9186 |
| `0302e8c^` | `411a448` | 08-05 16:53 | **84** | **0.3573** | 1.0000 | **858** | **0.6088** | 0.9990 |
| HEAD | `b398661` | 08-05 22:02 | 84 | 0.3573 | 1.0000 | 858 | 0.6088 | 0.9990 |

Selected labels are `{grade3 >= 1} & {pgr <= 10}` (DINA) and `{er <= 0}` (GRF)
at **every** state. Selection never changed; only the family and the
correction computed over it did.

The HEAD row reproduces the regenerated payload
(`rerun/gbsg_survival_multimethod_payload.rds`) exactly, and the `daa1e0a` row
reproduces the previous payload exactly, so the series is anchored at both ends
to real rendered output rather than to the script alone.

## What the table shows

**1. The family jump is entirely `84203ed`.** `n_family` goes 1 -> 84 and
1 -> 858 across that single commit boundary (`c68fe70` -> `411a448`). So the
GBSG frontier did NOT hold a single candidate: it held 84 and 858, and the
restriction cut both to one. That is 100% of the competition removed, against
the 99.84% the commit reported for ACTG175 GRF.

**2. `0302e8c` is inert on GBSG.** `411a448` and HEAD agree to every digit on
every field for both engines. This also follows from the mechanism: `0302e8c`
adds a beta-hat admission floor, which is strictly more restrictive and cannot
grow a family. Its own commit message shows GRF `n_family` unchanged at 1289
before and after, with only `admitted` shrinking (62 -> 16).

**3. `8d60b1e` moved the bias term without touching the family.** Bias falls
0.2632 -> 0.0272 (DINA) and 0.1771 -> 0.0604 (GRF) while `n_family` stays 1 and
the labels stay fixed. The driver is `selection_rate`: 0.4684 -> 0.9688 and
0.6796 -> 0.9186. Before `8d60b1e`, MR reconstructed the admission set itself
and returned a selection rate below 0.5 -- correcting a selection event it
judged had usually not occurred.

## Why the reasoning failed

The claim inferred a property of the frontier from a number that was the
restriction's own output. `.fs_mr_restrict_native()` received
`effect_neighborhood = 0.10` on every DINA and GRF run ever made -- the same
commit establishes this ("THE DEFAULT WAS NEVER OFF"). So the observed
`n_family = 1` on GBSG was what the restriction had already reduced the family
to, not evidence about how many candidates the frontier produced. The
before-value and the after-value were read from the same restricted path, and
`1 -> 1` was recorded as "not binding" when the correct reading required
measuring the family with the restriction absent.

The ACTG175 GRF cell (`2 -> 1289`) was measured that way and is correct. The
GBSG cells were not measured; they were inferred from the pre-removal value
alone. The commit's own closing line -- "Checked rather than assumed" -- applies
to the `effect_neighborhood` audit in that commit, not to this claim.

## Consequence

**Every pre-`84203ed` GBSG DINA and GRF MR result was corrected against a
family of one.** The selection-bias term is the winner's curse over the family
MR re-selects within; with one candidate there is no competition and
essentially no correction. The bias values make this concrete: 0.0272 and
0.0604 against a family of one, versus 0.3573 and 0.6088 once the real family
of 84 and 858 is present -- an order of magnitude in both engines.

This affects the manuscript-era GBSG payload
(`sim_analyses/gbsg_table_payload.rds`) and every `rerun/` GBSG payload built
before `84203ed`, in the DINA and GRF rows only. The FS row is unaffected:
`.fs_mr_restrict_native()` was called only from the DINA and GRF branches, and
FS `n_family` holds at 1744 across the same boundary.

## The GBSG baseline was compromised three ways, not one

The three defects are independent -- different commits, different mechanisms,
and each would have been present had the other two never existed. All three
were live simultaneously in the `daa1e0a` baseline.

| | commit | mechanism | evidence at baseline |
|---|---|---|---|
| 1 | `84203ed` | MR's family restricted to a 0.10 band of the best native statistic, cutting 84 and 858 candidates to 1 | `n_family` = 1 for both engines |
| 2 | `8d60b1e` | MR reconstructed the admission set instead of reading the resolved one | `selection_rate` 0.4684 (DINA), 0.6796 (GRF) |
| 3 | `dad0415` | the two MR bias terms did not share one denominator conditional on identification | bites wherever `selection_rate` < 1, which both engines were, heavily |

Defect 3 is the reason the low rates in defect 2 are not merely a diagnostic
curiosity. The F13 denominator mismatch moves the MR estimate only where
`selection_rate < 1`; at the baseline both engines sat far below 1, so the
mismatch was near its maximum effect exactly where the misreconstructed
admission set had put them. At HEAD the same two engines run at 1.0000 and
0.9990, where F13 is inert or nearly so -- so the F13 exposure is a property of
the baseline, not of the current configuration, and it disappears once defect 2
is fixed rather than because F13 itself was applied.

That interaction is worth stating plainly: fixing `8d60b1e` alone would have
removed most of defect 3's effect without touching `dad0415`, and fixing
`dad0415` alone would have left the family of one and the misreconstructed
admission set untouched. The three had to be read together, which is why no
single before/after comparison against the baseline identifies any of them on
its own.

Not affected: labels, `n`, `pct`, and every `naive_*` and `fb_*` column, all of
which are identical across all six states.

## Reproducing

`measure_gbsg_nfamily.R` reproduces both document calls standalone and prints
the four fields. Driver: `git worktree add <dir> <commit>`,
`R CMD INSTALL <dir>`, run, `git worktree remove --force`. Reinstall HEAD with
`devtools::install()` afterwards and verify behaviourally -- checking
`formals()` alone is insufficient, since a stale install can carry correct
defaults while missing later changes.

## Status

Nothing regenerated for this note; the five documents were regenerated
separately at HEAD and already reflect the post-`84203ed` families. The
correction is to the record in `84203ed`'s commit message, which cannot be
amended in place.

## The ACTG175 DINA cell is correct, but unrepresentative

`84203ed` also reported `DINA 1 -> 1 (not binding here)` for ACTG175. Measured
the same way, reproducing the shipped `fs_dina` call
(`analysis_actg175_binary_multimethod.qmd:1230`) at both the threshold that
cell used and the threshold the document actually runs:

| threshold | package state | `n_family` | bias | rate | label |
|---|---|---|---|---|---|
| `hr.threshold = 1.25` | `84203ed^` (`c68fe70`) | 1 | 0.5320 | 0.3214 | `{cd40 >= 500} & {homo >= 1}` |
| `hr.threshold = 1.25` | `0302e8c^` (`411a448`) | 1 | 0.5320 | 0.3214 | `{cd40 >= 500} & {homo >= 1}` |
| `hr.threshold = 1.25` | HEAD | -- | -- | -- | no subgroup (`0302e8c`) |
| `or_threshold = 1.0` | `84203ed^` (`c68fe70`) | **5** | 0.3894 | 0.9948 | `{preanti >= 849.4} & {cd40 >= 338}` |
| `or_threshold = 1.0` | `0302e8c^` (`411a448`) | **147** | 0.7081 | 1.0000 | `{preanti >= 849.4} & {cd40 >= 338}` |
| `or_threshold = 1.0` | HEAD | 147 | 0.7081 | 1.0000 | `{preanti >= 849.4} & {cd40 >= 338}` |

**At 1.25 the `1 -> 1` is genuine.** The two states are identical on every
field, so the restriction really was inert there -- the qualifying family at
that threshold holds one candidate, as `check_dina_actg.R` records (147
qualifying at floor 0, 89 at 1.05, 1 at 1.25). This cell does NOT have the
GBSG defect, and unlike the GBSG cells it was a true measurement.

**But it is not the configuration the document runs.** At
`or_threshold = 1.0`, which `analysis_actg175_binary_multimethod.qmd` uses, the
restriction cut 147 to 5 -- 96.6% of the family removed. That is the 5 -> 147
jump seen in the regenerated ACTG175 payload, and it lands on `84203ed`
exactly, since `411a448` and HEAD agree.

So the cell is accurate and misleading at once: "not binding here" is true of
`hr.threshold = 1.25` and false of the shipped threshold, and the commit does
not say which one it measured. The same conflation is recorded in
`HANDOFF_2026-08-05_session2.md` -- DINA's "family of one" on ACTG175 is a
property of the threshold, not of DINA -- and it recurred here in a different
form.

Note also the `selection_rate` of 0.3214 at 1.25: a subgroup selected on two
draws in three that it was not, the signature `0302e8c` describes for a winner
sitting below `t_g`. That commit reports 0.338 for the comparable cell; the
small gap is a configuration difference between the two measurements
(`sg_focus`/`selection_rule`) and was not chased.
