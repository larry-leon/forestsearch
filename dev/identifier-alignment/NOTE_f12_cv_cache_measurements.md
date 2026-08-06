---
title: "F12 — the CV cache-leak measurements"
subtitle: "Supporting numbers for the fix in 917e15d"
date: 2026-08-05
bibliography: []
---

# Why this note exists

`NEWS.md` describes the F12 fix — cross-validation reusing a cached `grf_res` /
`dina_res` inside every training fold — and states the **direction** of the
effect without quoting magnitudes. The measured numbers live here instead.

They are a **wiring demonstration**, not an operating-characteristics result,
and NEWS is the wrong place to put figures a reader would take as stable. This
note records them with the caveats attached.

# What was measured

The leak only bites when a cached fit is supplied, so the failing case was built
deliberately: fit once, harvest the returned `grf_res` / `dina_res`, then make a
**second** `forestsearch()` call passing it back in. That puts the cached object
into `args_call_all`, which is what `forestsearch_Kfold()` reconstructs
`cv_args` from — so a fold that does not null it ranks over a family estimated
on all the data, including the fold being held out.

`before` is the pre-fix `R/forestsearch_cross_validation.R`; `after` is the
same file with the four nullings added at both entry points.

ACTG175, binary outcome, OR, `Kfolds = 3`, `seedit = 8316951`.

| method | cached | mode | cache in `args_call_all` | `sens_H` | `ppv_H` | `sens_Hc` | `ppv_Hc` | exact |
|---|---|---|---|---:|---:|---:|---:|---:|
| GRF | yes | **before** | TRUE | 0.55814 | **0.64865** | 0.97392 | 0.96234 | **0.33333** |
| GRF | yes | after | TRUE | 0.54651 | 0.48958 | 0.95085 | 0.96049 | 0.00000 |
| GRF | no | before | FALSE | 0.54651 | 0.48958 | 0.95085 | 0.96049 | 0.00000 |
| GRF | no | after | FALSE | 0.54651 | 0.48958 | 0.95085 | 0.96049 | 0.00000 |
| DINA | yes | **before** | TRUE | **0.46875** | **0.30000** | 0.93131 | 0.96541 | 0.00000 |
| DINA | yes | after | TRUE | 0.21875 | 0.15385 | 0.92444 | 0.94960 | 0.00000 |
| DINA | no | before | FALSE | 0.21875 | 0.15385 | 0.92444 | 0.94960 | 0.00000 |
| DINA | no | after | FALSE | 0.21875 | 0.15385 | 0.92444 | 0.94960 | 0.00000 |

The pattern is the one the leak predicts, in both engines: the cached run is
optimistic before the fix, and after it the cached and uncached runs agree
**exactly** on every metric. Rows with `grf_res = NULL` / `dina_res = NULL`
(the default) are unchanged by the fix, because there was nothing to reuse.

`cache in args_call_all` is `TRUE` in exactly the cached cells, which confirms
empirically — rather than by inference from the bootstrap's behaviour — that
`args_call_all`'s bulk `mget()` is the route by which the cached object reaches
the folds.

# Caveats, which is the point of keeping these out of NEWS

**The two engines were run at different thresholds.** GRF at
`hr.threshold = 1.25`, DINA at **1.05**. The DINA arm could not be run at 1.25:
since `0302e8c`, no candidate clears the beta-hat floor there, DINA selects
nothing, CV has no original subgroup to cross-validate, and every metric is
`NA` — that measures the F2/F3 fix rather than this one. So the table mixes two
configurations and the two engines' rows are not comparable with each other.

**`Kfolds = 3` makes exact-match granular to 1/3.** The GRF exact-match figure
of `0.33333 → 0.00000` is **one fold out of three**. It is a real difference in
that run and it is the right direction, but it is not an estimate of a rate.

**One dataset, one seed.** Nothing here supports a claim about the size of the
effect in general — only that the leak exists, that it is optimistic in
direction, and that the fix closes it exactly.

# What would be needed to state magnitudes

A sweep over datasets and seeds at a realistic `Kfolds` (10), with both engines
at a common threshold, reporting distributions rather than point values. That
was out of scope for the fix and has not been run.
