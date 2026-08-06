---
title: "Handoff — session 2, 2026-08-05"
subtitle: "State document. Branch mr-in-replicates, pushed."
date: 2026-08-05
bibliography: []
---

Reasoning behind settled decisions is in the specs and notes, not repeated here.

## What landed

Branch `mr-in-replicates`, **pushed** (`e956d0d..972f8fc`). Ten commits this
session, oldest first:

| commit | |
|---|---|
| `411a448` | `.dfbeta_glm()` errored on every `lm()` fit — MR was dead for `effect_measure = "MD"` |
| `0302e8c` | DINA and GRF admit on β̂ ≥ t_g, from the resolved admission set, not a literal or nothing |
| `1856cba` | filed the closed-form verification scripts and the DINA re-run script |
| `917e15d` | nulled cached `grf_res`/`dina_res` in CV folds, matching the bootstrap |
| `7d94636` | corrected `verify_residual_centering.R`'s header after the `forestsearch:::` change |
| `6f97d7c` | F4 findings (Q1, Q4) and the `.consistency_glm_pieces()` docstring correction |
| `dad0415` | one shared denominator across both MR bias terms, conditional on identification |
| `1ed6d38` | `max_subgroups_search` now defaults to `Inf` |
| `972f8fc` | uncapped the candidate pool in both analysis documents |

Plus the handoff commit carrying this file.

`R CMD check --as-cran` was clean (0/0/0) before each commit touching package
code.

## Still uncommitted or untracked

**`actg175_realignment_report.qmd`** — present at the start of this session, so
it predates all of the above. Never opened beyond confirming it is a
PDF-format ACTG175 re-analysis report. Left uncommitted because it is not a
spec, a note, or the audit, and was not part of any instruction this session.
**It is the only deliberate omission.** Decide whether it belongs in the repo.

Rendered `.html` outputs are gitignored and are not commit candidates:
`code_theory_audit.html`, `f4_ps_effect_gap_findings.html`,
`dfbeta_glm_validation.html`, and the `rerun/` bundles.

Working tree is otherwise clean.

## The two "also record" notes from F13

Both written, both in **`NOTE_f13_open_questions.md`** (committed in `dad0415`):

1. **What does the full bootstrap record for a replicate that identifies
   nothing?** If it drops the replicate, FB is itself conditional on
   identification and independently supports the convention chosen in
   `dad0415`. **Unverified** — `bootstrap_analysis_dofuture.R` was not read for
   this. Cheap to settle.
2. **Is the manuscript silent on the no-winner case, or contrary to it?** Eq. 9
   divides by `B`, which appears to assume every draw admits a candidate; if so
   the manuscript is silent rather than contrary, and no amendment is implied.
   **Unverified** — needs a read of §3.1.

That note also records the `fixed_bias` convergence sweep, and that Eq. 13
needs no amendment under either convention.

## Open findings not addressed

From `code_theory_audit.qmd` (committed this session):

| | |
|---|---|
| **F4 Q2** | how far apart the PS-adjusted and unadjusted effects are for the selected subgroup |
| **F4 Q3** | whether the *ranking* changes — Spearman over the family, argmax agreement, admission-status changes. **Load-bearing**: decides whether F4 is a variance-estimation problem or a selection problem, which need different fixes |
| **F9** | two `sg_focus` dispatch sites `.assert_sg_focus_dispatch_complete()` does not read |
| **F10** | dead `.mr_cscr` / `c_screen_mr` / `c_consistency_mr` variables |
| **F11** | `max.minutes` is a formal no code path consults |
| **F14** | both replay paths re-derive quantile cut locations per replicate/fold, so ℱ is not fixed across them |

F4's three candidate resolutions were never assessed; the choice turns on
whether the manuscript's scope covers non-randomised data, which is unresolved.

**Two further items surfaced by F4 Q4**, recorded in
`f4_ps_effect_gap_findings.qmd` and not chased:

- **The bootstrap carries a user-supplied `ps_hat` into replicates unaligned.**
  CV nulls it as "wrong length for training fold"; a bootstrap resample has the
  *same* row count, so the length check at `forestsearch_main.R:1922` passes and
  the original subjects' scores are applied to resampled rows. Bites only when
  `ps_hat` is passed explicitly.
- **Survival + `ps_method != "none"` appears to be a complete no-op.** The
  estimator-closure rebuild is guarded by `if (is_glm)` and no Cox path reads
  weights, so PS is estimated, the columns are attached, and nothing consumes
  them — while the object still reports `ps_method = "grf"`.

## Regeneration — still held, and now further behind

Regeneration was held throughout and **was not performed**. Five documents need
re-running, now for three separate reasons:

1. **`0302e8c`** — DINA and GRF admission changed. DINA on ACTG175 selects
   *nothing* at `hr.threshold = 1.25` where it previously reported a subgroup.
2. **`dad0415`** — MR point estimates move wherever `selection_rate < 1`
   (−8.2% on the binding configuration measured). Unaffected at rate 1.
3. **`1ed6d38` + `972f8fc`** — the caps are gone, and **all four previously
   bound**, so FS and DINA-screen labels may differ.

**`411a448` and `917e15d` may be inert for these five documents** — the first
affects only `effect_measure = "MD"`, the second only runs passing a cached
`grf_res`/`dina_res`. Neither condition was checked against the documents.
**This is unverified**; do not assume it.

Two measurements already on record bound part of the blast radius: the shipped
DINA configuration (`or_threshold = 1.0`, `sg_focus = "effMaxSG"`) is
**unchanged** by `0302e8c` on every field, and both it and GRF on ACTG175 ran at
`selection_rate = 1`, so `dad0415` does not move them either. The uncapping is
the change most likely to alter what those documents report.

## Things not in a commit message or a note

**Corrections made to my own earlier statements, so they are not repeated:**

- I reported the first `R CMD check` of the session as having "died when
  backgrounded and produced no output". It had not — the output was buffered
  and it completed normally. A second, redundant check was launched on that
  mistaken basis.
- I stated I had reverted my Eq. 9 `NEWS.md` entry and that the file carried
  only the `max_subgroups_search` work. **Wrong** — the entry was still present
  and was replaced later. Check `NEWS.md` state directly rather than trusting a
  claim about it.
- I described DINA's "family of one" on ACTG175 as an oddity in DINA. It is a
  property of `hr.threshold = 1.25`: 147 qualifying at floor 0, 89 at 1.05, 1 at
  1.25. `check_dina_actg.R` warns specifically against this conflation.
- I reported the `SPEC_f12` filing instruction as missing from that spec. It
  was — the spec has no such section, confirmed against every copy on disk. The
  filing scheme came from a direct answer, not the spec.

**Operational:**

- A stray `until`-loop shell survived session 1, spinning 1h37m on a log whose
  check had been killed. Kill background waiters when killing what they wait on.
- `devtools::check(args = c("--no-build-vignettes"))` does **not** skip vignette
  *building* — that flag applies to the check stage. Each run is ~15 minutes.
- Three `NEWS.md` commits this session needed hunk-level staging (revert to
  HEAD, re-apply own edit, stage, restore) because the file simultaneously
  carried uncommitted work belonging to someone else. That is now resolved —
  `NEWS.md` is fully committed.
- `verify_residual_centering.R` and `verify_eq9_alignment.R` cite line ranges in
  `R/fs_mr_inference.R`. `dad0415` shifted them once already. They were updated;
  they will go stale again on the next edit to that file.

**Deliberately left alone:**

- The F12 `NEWS.md` entry presents a four-row cached/uncached table under "the
  effect is not small". Those numbers mix GRF at `hr.threshold = 1.25` with DINA
  at 1.05, and exact-match `0.3333 → 0.0000` is one fold of three at
  `Kfolds = 3`. A correction was drafted and **deprioritised**; the entry shipped
  as-is. The caveats are in `NOTE_f12_cv_cache_measurements.md`.
- `verify_residual_centering.R`'s two arms describe *pre-repair* behaviour and
  no longer match shipped code. Kept deliberately — it is the record of the
  diagnosis — and its header says so.
