# CC BRIEF — implement all twelve bug-review findings, then verify

```
claude "Read dev/review/CC_BRIEF_qmd_fixes.md and execute it."
```

**Supersedes any earlier fix brief.** This covers every finding in
`dev/review/QMD_BUG_REVIEW.md` — F1–F7 and N1–N5, both files — with an explicit
disposition for each. Nothing from that report is left unaddressed; where the
disposition is "no action", the reason is stated.

Targets:

* `quarto/applications/gbsg/analysis_gbsg_survival_multimethod.qmd` — 2507 lines, 69 chunks
* `quarto/simulations/gbsg_redux/sim_fs_maxcons_fb_mr_m1_h10_knoise0_n500_batch_1_100.qmd` — 1329 lines, 11 chunks

**Confirm you are on the reviewed content before editing.** The multimethod file
should contain `mr_in_replicates = FALSE` exactly 3 times and a `fs_param_table()`
filtering on `.fs_selection_knobs`; the sim file should contain
`mr_in_replicates` exactly twice (one comment, one argument) and a
`RECORD ONLY` provenance block. If either differs, stop and report.

---

## 1. Two principles

**(a) The code is the source of truth; a stale comment is a documentation fix.**
The code in both files has been reviewed and executed. The comments have not.
Where they disagree, correct the comment — *unless* the code it describes is
itself producing an unwanted effect, in which case fix both. F3 is the only such
case.

**(b) A rendered claim must be generated from the object the code produced, or
gated by the same condition as that code.** F2, F5, F6 and F7 are one defect in
four places: a hand-typed parameter table beside pinned code; a timing table
listing rows by name; narrative describing analyses a toggle skipped; a
stability table for a fit that did not run. Prefer generating over correcting a
hardcoding — a corrected hardcoding drifts again.

**No fix in this brief may change a computed value.** §4 verifies that.

---

## 2. Disposition table

| # | file | disposition |
|---|---|---|
| F1 | sim | **No action** — filename/batch-range mismatch, deferred by the maintainer |
| F2 | multimethod | **Delete** the hardcoded parameter table |
| F3 | multimethod | **Code + comment** — the only code change here |
| F4 | multimethod | **Prose** — seven lowercase `tier`/`tiers` |
| F5 | multimethod | **Generate** the timing table |
| F6 | multimethod | **Gate** the GRF-cut stability block |
| F7 | multimethod | **Gate** the CV narrative |
| N1 | multimethod | **Pin** `mr_in_replicates` at two sites |
| N2 | sim | **No action** — see §3.11 |
| N3 | sim | **Comment** |
| N4 | sim | **Ask** — see §3.10 |
| N5 | multimethod | **Comment**, two sites |

---

## 3. The fixes

### 3.1 F3 — `stop_threshold = 0.95` warns on every render (multimethod ~753, ~1286)

The only finding where code changes.

`sg_focus = "effMaxSG"` normalises to `"hrMaxSG"`, which is in the reset list
`c("hrMaxSG", "hrMinSG", "hr", "maxSG", "minSG")`, and 0.95 is supplied
explicitly, so `user_explicit` is `TRUE` and the reset warns on every render.

Change both sites to `stop_threshold = NULL`. Then `!is.null(stop_threshold)` is
`FALSE`, the reset branch is never entered, no warning fires, and the setting
stays visible. Replace the adjacent comment — it states what the default is
without saying the value is reset — with one saying `NULL` is what this
configuration resolves to, and why.

**Verify the no-op rather than assume it:** both routes reach `NULL`.

### 3.2 F5 — timing table omits the GRF rows (multimethod ~2220–2237)

`timing_df` lists ten rows by hand. `timings$fs_grf_select` and
`timings$fs_grf_bootstrap` are absent, but `Total` is
`proc.time() - t_vignette_start` and includes them, so the column does not
reconcile.

Drive the rows from `names(timings)` with a name→label map for readable row
names. **Any name in `timings` with no label must `stop()`**, not be silently
omitted — that is what makes this a fix to the class rather than the instance.
`Total` stays last.

Report row count before/after and the reconciliation gap both ways.

### 3.3 F7 — CV narrative renders while CV is skipped (multimethod ~1940–2022)

Every CV chunk carries `eval = run_cv`; ~60 lines of surrounding prose do not.
With `run_cv <- FALSE` the document describes a ten-fold analysis it never ran,
including a hardcoded `50 × 10 = 500`.

Move the narrative under the same gate — a `results = "asis"` chunk with
`eval = run_cv` — and compute that arithmetic inline from the knobs. Keep the
existing `cv-skipped-note` (`eval = !run_cv`) and make the two complementary:
exactly one renders, never both, never neither.

Apply the same treatment to any LOO narrative. Report which blocks you moved.

### 3.4 F6 — GRF-cut stability renders a degenerate row (multimethod ~1124–1156)

The main fit sets `use_grf = FALSE`, so this renders a single
`(no GRF cut) — 100%` row presented as a result. Gate prose and chunk on
`use_grf`, the same pattern as `run_cv`, with a skipped-note counterpart. Read
the guard from the fit if it exposes it, otherwise from the knob the call passes.

### 3.5 F2 — hardcoded parameter table (multimethod ~612, 617, 1177, 1351)

Four typed values disagree with what the code pins.

**Delete the table.** It is redundant with the per-fit provenance tables, which
read `fit$args_call_all` at render time, and it is the half that drifted. If a
summary is wanted there, render it with `fs_param_table()`.

Check the other three sites individually — each is either a value that should
come from the fit, or prose that should not state a number. Report each.

### 3.6 F4 — seven lowercase `tier`/`tiers` (multimethod 1058, 1572, 1583, 1646, 1801, 1811, 1878)

The Tier→FB/MR sweep replaced `Tier` and `tier1` but not lowercase standalone
forms. These are prose — read each in context; replace with FB/MR or a neutral
word as the sentence requires.

**Protected words that legitimately contain the substring:** `frontier`,
`delegate`/`delegated`/`delegates`, `gatekeeper`, `ungated`, `Negate`. Confirm
their counts are unchanged afterwards.

### 3.7 N1 — `mr_in_replicates` unpinned at two sites (multimethod ~1524, ~1755)

The DINA and GRF bootstraps call `forestsearch_bootstrap_dofuture()` without the
argument while the main call at ~914 pins it with a five-line comment. Pin both,
with a one-line comment pointing at the main site rather than repeating it.

Verified no-op: `formals()` returns `FALSE` and
`bootstrap_analysis_dofuture.R:597` strips the flag regardless.

### 3.8 N5 — two comments describe behaviour the code does not have (multimethod 119–121, 1520–1522)

Both confirmed against the code:

* **119–121** — the comment says *"limited to 2 cores for CRAN checks"*; the code
  is `max(1L, floor(0.80 * parallel::detectCores(logical = FALSE)))`, i.e. 80% of
  physical cores. Correct the comment to describe that.
* **1520–1522** — the comment says *"Diagnostic run: sequential plan"*; the call
  immediately below passes `plan = "multisession", workers = n_cores`. It does
  set `details = TRUE`, so keep that part and drop the sequential claim.

Comment-only. Do not change the core count or the plan.

### 3.9 N3 — the FB-guard comment mis-states the failure mode (sim 566–567)

The comment says a missing `H2` column would leave a short record *"silently
dropped by the run-loop's rbind/.errorhandling = 'remove'"*. The review found it
**errors** rather than being dropped.

The guard itself is correct and stays. Correct the comment to state the actual
failure mode. Confirm the mechanism yourself before rewording — the run loop is
at sim ~610 (`.combine = rbind, .errorhandling = "remove"`).

### 3.10 N4 — recorded but never read (sim 399–402) — **ASK, do not act**

`mr_ok`, `n_sel`, `label` and `covs` are written into every record and persisted
to the `.rds`, but nothing in the document reads them.

Two readings: they are deliberate payload for downstream analysis of the `.rds`,
or they are dead weight. **Do not remove them.** Report which columns are
affected, confirm nothing in either document reads them, and ask the maintainer.
Removing persisted columns would change the `.rds` schema, which is not a
display-only change.

### 3.11 N2 — twelve unguarded assignments (sim 515–526) — **no action**

The MR block has 11 unguarded `rec$… <- g$…` assignments against 4 guarded; the
FB block below it guards 6 of 6. The asymmetry is real, but the review proved the
unguarded ones safe against the current package.

Adding guards that cannot fire would be noise, and the comment above the FB
block already explains why *that* one needs guarding. Leave it. Note it in the
deliverable as a known asymmetry with a rationale, so it is not rediscovered.

### 3.12 F1 — batch range vs filename (sim) — **no action**

`n_sims <- 50L` with `sim_id_start <- 1L` writes `..._res_1_50.rds` while the
document is named `..._batch_1_100.qmd`. Deferred by the maintainer: the file is
a template and will be renamed once the structural work is settled.

Do **not** rename anything or change `n_sims`. Record in the deliverable that
this remains open, including the review's finding that combining the on-disk
`_res_1_100.rds` with a batch written now would `stop()` on
`meta$nb_boots` (300 vs 20) and on duplicated `sim_id`.

---

## 4. Verification

### 4.1 Static, both files

1. **Symbol resolution.** Extract and source each setup chunk in a clean
   session; each must complete with no `object ... not found`. Then list any
   symbol read before it is assigned earlier in the document.
2. **`sprintf` arity.** Every `sprintf()`/`cat(sprintf())`: format specifiers vs
   arguments supplied. F5 and F7 both touch format strings.
3. **Residual vocabulary, case-insensitive.** Search `[Tt]ier` and `[Gg]ate`
   across both files, excluding the protected words. The earlier check used
   `Tier|tier1`, which is exactly why F4 survived — do not repeat that pattern.
4. **Chunk gating.** For every `eval =`, confirm the guard is defined before the
   chunk and that nothing after a gated block reads an object it creates.
5. **Fences and labels.** Count before/after; confirm no duplicate chunk labels.

### 4.2 Rendered — multimethod

Render once, then:

* **Numeric invariance is the hard requirement.** Compare the emitted
  `gbsg_survival_multimethod_payload.rds` against the previous one at
  `tolerance = 0`, matched on method and region. **Every numeric column must be
  identical.** Any difference means a fix was not display-only — report and stop.
* **No reset warning** in the render log (F3).
* **Timing table reconciles** and carries the GRF rows (F5).
* **CV and GRF blocks:** with `run_cv <- FALSE` and `use_grf = FALSE`, neither
  narrative renders and each skipped-note appears exactly once (F6, F7).
* **Diff the rendered HTML** against the previous render and account for every
  difference — each should map to one of F2–F7. Anything unaccounted for is a
  finding.

### 4.3 Sim

The sim changes are comment-only (N3), so no render is required. Confirm the
extracted chunks still parse and that the file is otherwise byte-identical apart
from the comment.

---

## 5. Deliverable

`dev/review/QMD_FIXES_APPLIED.md`:

1. **Verdict** — numeric invariance held or not; anything blocking.
2. **A row for all twelve findings** F1–F7, N1–N5, each marked fixed / no action
   / awaiting maintainer, with evidence or the reason. Nothing omitted.
3. §4.1–4.3 results per check.
4. The rendered-HTML diff accounting.
5. The N4 question for the maintainer.
6. Anything found while editing that the review did not anticipate.

Commit the two `.qmd` files and the multimethod render together. Report and stop.
