# QMD fixes applied — all twelve findings

Targets (both confirmed byte-identical to the reviewed copies in `dev/review/`
before any edit):

* `quarto/applications/gbsg/analysis_gbsg_survival_multimethod.qmd` — 2507 → 2632 lines, 69 → 79 chunks
* `quarto/simulations/gbsg_redux/sim_fs_maxcons_fb_mr_m1_h10_knoise0_n500_batch_1_100.qmd` — 1329 → 1335 lines, 11 chunks (comment-only)

---

## 1. Verdict

**Numeric invariance HELD. Every numeric value in the payload is bit-identical
at tolerance 0.** No fix changed a computed value.

Establishing that required a control, because the stored baseline was not
comparable:

| comparison | source | machine | result |
|---|---|---|---|
| **AFTER vs CONTROL** | edited vs **unedited** | both this box (x86_64, 128 cores) | **11 / 11 numeric columns bit-identical; 0 differences** |
| CONTROL vs BEFORE | unedited vs **unedited** | x86_64 vs arm64 | 8 / 11 columns differ, ≤ 2.2 ULP |
| AFTER vs BEFORE | edited vs unedited | x86_64 vs arm64 | 8 / 11 columns differ, ≤ 2.2 ULP — **the same 8 columns at the same magnitudes** |

The second and third rows are identical, and the first is empty. The drift is
therefore 100% attributable to the machine change and 0% to the edits.

**The brief's §4.2 test as written could not have passed on any render from this
machine, including one with no source changes at all.** The stored
`gbsg_survival_multimethod_payload.rds` was produced on an **Apple M4 Max
(14 physical cores, arm64, 11 workers)**; this box is an **AMD Threadripper PRO
5995WX (128 cores, x86_64, 102 workers)**. Cross-architecture floating point is
not bit-comparable, and the worker count feeds the candidate-search batch size
(`batch_size = 22` there, `30`/`50` here), which changes summation order. The
payload's own `.machine_info()` block — recorded so "timings pooled across
different machines stay attributable" — is what surfaced this.

Nothing is blocking.

Supporting evidence beyond the payload:

* All analysis result lines in the rendered HTML are **byte-identical** between
  the control and edited renders: `naive HR (95% CI)` (12 lines), `de-biased HR
  (95% CI)` (12), the three `*-selected subgroup:` lines (13), Cohen's kappa (8),
  `Candidate family:` (2), `selection bias` (11).
* Bit-identical across *all three* renders regardless of machine: subgroup
  labels (`{er <= 0} & {size <= 35}`, `{grade3 >= 1} & {pgr <= 10}`,
  `{er <= 0}`), `n_family` (1744 / 1 / 1), `extras$concordance`,
  `extras$harm_rates`, `meta$B`, `mr_draws`, `n_total`, `est_scale`,
  `forestsearch_version`.
* Even cross-architecture, every numeric column agrees to 6 decimal places, and
  `naive_est` (plain `coxph`, no matrix/parallel path) is bit-identical.

---

## 2. All twelve findings

| # | file | disposition | status | evidence |
|---|---|---|---|---|
| **F1** | sim | no action (deferred) | **open, unchanged** | §3.12 below — recorded, nothing renamed |
| **F2** | multi | delete hardcoded table | **fixed** | table deleted; 3 further sites handled individually (§3.5) |
| **F3** | multi | code + comment | **fixed** | `stop_threshold = NULL` ×2; reset warning 2 → **0** in the render; A/B proves no-op |
| **F4** | multi | prose | **fixed** | 7 replaced; rendered prose 5 → **0**; protected words unchanged |
| **F5** | multi | generate table | **fixed** | 10 → **12** rows; gap 202.4 s (30%) → **5.8 s (1%)** |
| **F6** | multi | gate GRF block | **fixed** | degenerate `(no GRF cut) 100.00` row gone; skipped-note renders once |
| **F7** | multi | gate CV narrative | **fixed** | 8 blocks moved; `50 × 10 = 500` now computed |
| **N1** | multi | pin ×2 | **fixed** | DINA + GRF bootstraps pinned |
| **N2** | sim | no action | **not acted on, by design** | §3.11 below |
| **N3** | sim | comment | **fixed** | mechanism re-confirmed empirically first |
| **N4** | sim | ask | **awaiting maintainer** | §5 below — nothing removed |
| **N5** | multi | comment ×2 | **fixed** | both corrected |
| — | multi | *not in the brief* | **fixed** | "This document never sets it" (§6.1) |

---

## 3. Per-fix detail

### 3.1 F3 — `stop_threshold` (the only code change)

Both sites (`forestsearch-main`, `fs-dina-screen`) changed `0.95` → `NULL`, and
the misleading comment replaced with one stating that `NULL` is what the
configuration *resolves* to and why.

**No-op verified, not assumed.** Two fits differing only in this argument, same
session, same seed:

```
reset warning fired?            old: TRUE    new: FALSE
args_call_all$stop_threshold    old: NULL    new: NULL
selected subgroup identical     TRUE   ({er <= 0} & {size <= 35})
naive est identical             TRUE   (2.5369175910 / 2.5369175910)
debiased est identical          TRUE   (1.4308462042 / 1.4308462042)
debiased CI identical           TRUE
n_family identical              TRUE
df.est treat.recommend          TRUE
grp.consistency out_sg          TRUE
```

Both routes reach `NULL`: `forestsearch_main.R:1402` guards on
`!is.null(stop_threshold)`, so passing `NULL` skips the reset branch entirely,
while any non-`NULL` value enters it, warns (`user_explicit <- !missing(...)`),
and is then set to `NULL` unconditionally at `:1435`.

Rendered result: the reset warning appears **twice** in the control render and
**zero** times in the edited one.

### 3.2 F5 — timing table generated

Rows are now driven by `names(timings)` through a name→label map. An unlabelled
timing **stops the render** — exercised directly:

```
STOP -> Computational Timing: no label for timings$brand_new_thing
        -- add it to .timing_labels in the timing-summary chunk.
```

Rendered, before → after:

| | rows | components summed | Total | unattributed |
|---|---|---|---|---|
| control (hand-written) | 10 | 479.5 s | 681.9 s | **202.4 s (29.7%)** |
| edited (generated) | 12 | 677.4 s | 683.3 s | **5.8 s (0.9%)** |

The 202.4 s gap was almost exactly the two missing rows, now present:
`FS (subgroup_method=grf)` 3.4 s and `FB bootstrap (GRF)` 193.5 s. Labels were
also made readable (`Bootstrap` → `FB bootstrap (main)`, `kfold` → `K-fold CV`).
Skipped analyses keep a row showing `NA` rather than vanishing, and the
reconciliation is stated in a footnote so a future gap is visible in the
document itself.

### 3.3 F7 — CV narrative gated

Eight prose blocks moved into `results = "asis"` chunks under `eval = run_cv`:
`cv-diag-heading` (which also carries the `### Diagnostic Tabulations` heading,
so no empty section remains), `cv-identification-note`, `cv-grf-cut-note`,
`cv-no-subgroup-note`, `cv-pcons-note`, `cv-fold-numeric-note`,
`cv-original-agreement-note`, plus `loo-note` under `eval = run_loo`.

The hardcoded `50 × 10 = 500` is now `sprintf("%d x %d = %d", Ksims, Kfolds_cv,
Ksims * Kfolds_cv)`. `Kfolds_cv <- 10` was added to the setup chunk and feeds
**both** the `forestsearch_tenfold(Kfolds = Kfolds_cv)` call and the prose, so
the two cannot diverge. The `0.90` threshold in the pconsistency prose now reads
`fs$args_call_all$pconsistency.threshold`.

Complementary gates confirmed: `cv-skipped-note` (`!run_cv`) /
`loo-skipped-note` (`!run_loo`) against everything else — exactly one side
renders. In the edited render the CV narrative is absent and each skipped note
appears once.

### 3.4 F6 — GRF-cut block gated

Guard read **from the fit**, not the call: `fs_used_grf <-
isTRUE(fs$args_call_all$use_grf)`, defined in `forestsearch-main` immediately
after the fit-integrity `stopifnot`. Prose and table both gated on it, with a
`grf-cuts-skipped-note` counterpart on `!fs_used_grf`.

Rendered: the degenerate `(no GRF cut) — 100.00%` row is gone; the skipped note
renders once. (Two `no GRF cut` strings remain in the *echoed comment* that
explains the guard — text, not a result.)

### 3.5 F2 — hardcoded parameters, all four sites

The table is **deleted**, replaced by a pointer to the provenance tables that
read `args_call_all` at render time. The other three, each checked individually:

| site | was | actual | action |
|---|---|---|---|
| table row (~612) | `hr.threshold` 1.25 | 1.0 | deleted with the table |
| table row (~617) | `d0.min`,`d1.min` 12 | 10, 10 | deleted with the table |
| DINA intro (~1177) | "mirror the main analysis … `hr.threshold = 1.25`" | screening 0.9, selection 0.90 — and it does *not* mirror the main fit's 1.0 | rewritten: names the knobs (`dina_screen_hr_threshold`, `fs_dina_hr_threshold`), states that the harm threshold is the one thing not inherited, and why |
| standalone DINA (~1351) | "harm floor `m_diff = log(1.25)`" | `log(fs_dina_hr_threshold)` = `log(0.90)` | rewritten to name the symbol, not a value |

The rest of the deleted table was in fact correct (`hr.consistency` 1.0,
`pconsistency.threshold` 0.90, `maxk` 2, `n.min` 60); it was deleted as a class
per the brief rather than corrected, since a corrected hardcoding drifts again.

### 3.6 F4 — seven `tier`/`tiers`

| line | was | now |
|---|---|---|
| 1058, 1646, 1878 | "descriptive columns and naive effect are **tier-independent**" | "…are **the same under either correction**" |
| 1572, 1801 | "the two **tiers** are not expected to agree" | "**FB and MR** are not expected to agree" |
| 1583, 1811 | "so it agrees across **tiers** regardless" | "so **FB and MR agree on it** regardless" |

Protected-word counts unchanged: `frontier` 33, `delegate` 3, `delegated` 1,
`delegates` 2, `gatekeeper` 1, `ungated` 1, `Negate` 0.

### 3.7 N1 — `mr_in_replicates` pinned at two sites

Added to `fs-dina-bootstrap` and `fs-grf-bootstrap` with a one-line comment
pointing at the main site. The document now pins it at all five call sites.
No-op reconfirmed: `formals()` returns `FALSE` for all three entry points, and
`bootstrap_analysis_dofuture.R:597` / `forestsearch_cross_validation.R:402,:866`
strip the flag regardless. Payload numerically identical, as above.

### 3.8 N5 — two comments

* `~119` — "limited to 2 cores for CRAN checks" → describes what the code does
  (80% of physical cores) and notes that the CRAN cap lives in the package's
  `.default_parallel_workers()`, not here. Core count unchanged (102).
* `~1520` — dropped the false "sequential plan" claim, kept the `details = TRUE`
  part, and noted the plan is `multisession`. Plan unchanged.

### 3.9 N3 — the FB-guard comment (sim)

Mechanism confirmed before rewording, in all three relevant paths:

```
do.call(rbind, lapply(...))                   -> ERROR: numbers of columns of arguments do not match
foreach .combine=rbind .errorhandling=remove  -> ERROR: Failed to combine foreach() %dofuture% results
control: body-level stop() with same settings -> combined ok, 4 rows  <- body error WAS removed
```

So `.errorhandling = "remove"` genuinely removes loop-body errors but does not
cover a `.combine` failure. The comment now states that a short record **aborts
the whole batch** rather than being silently dropped, and names both loop paths.
The guard itself is unchanged.

### 3.10 N4 — recorded but never read — **no action, see §5**

### 3.11 N2 — twelve unguarded assignments (sim 515–526) — no action

Left exactly as-is, per the brief. Recording the rationale so it is not
rediscovered: the MR block has 12 unguarded `rec$… <- g$…` assignments against
the FB block's 6-of-6 guarded. The asymmetry is real but cannot fire —
`fs_mr_inference.R:432` builds `naive` and `debiased` unconditionally, and the
complement is either fully populated (`:415`) or `list(note = …)` with no
`debiased`, which the `!is.null(g$complement$debiased)` guard at sim L518
already excludes. Guards that cannot fire would be noise. Should
`fs_mr_inference()`'s return shape ever become conditional, these 12 lines are
where it would bite — and per N3, the failure would abort the batch, not lose a
cell.

### 3.12 F1 — batch range vs filename (sim) — no action, still open

Nothing renamed, `n_sims` untouched. Recorded for the maintainer:

`n_sims <- 50L` with `sim_id_start <- 1L` resolves `rds_path` to
`fs_maxcons_fb_mr_m1_h10_knoise0_n500_res_1_50.rds`, while the document is named
`..._batch_1_100.qmd` and `..._res_1_100.rds` already sits in the same
directory, inside the same `combine_glob`. Executed against the real on-disk
bundle plus the batch a render would now write:

```
agreement check   -> STOP: combine mode: batches disagree on meta$nb_boots (300, 20).
duplicated sim_id -> 50   (would STOP the pool at the next check)
```

Both guards fire loudly rather than pooling silently, but a render gives no
warning that it has just created that state.

---

## 4. Verification

### 4.1 Static, both files

1. **Symbol resolution.** Both setup chunks extracted and evaluated in clean
   `--vanilla` sessions with `library()`/`source()` stubbed — both completed, no
   `object … not found`. Multimethod resolves 39 knobs (`run_cv` FALSE,
   `run_loo` FALSE, `Ksims` 50, **`Kfolds_cv` 10**, `n_cores` 102, `NB` 1000,
   `run_mr` TRUE, `mr_draws` 5000, `fs_sg_focus` "effMaxSG"); sim resolves 63.
   Forward-reference scan over every chunk, chunk option and inline `` `r ` ``
   in document order: **18 flags before, 18 after** — the same set, all
   previously classified false positives (`$`-field names harvested by
   `all.vars` from `eval =` guards, `within()` columns, gt/ggplot NSE). No new
   forward reference. New guards are defined far ahead of first use:
   `fs_used_grf` L795 → first used L1143; `Kfolds_cv` L147 → L1977.
2. **`sprintf` arity.** 50 calls in the multimethod (was 47), 57 in the sim,
   **zero mismatches**. Three new format strings are `paste()`-built and so
   invisible to the static checker; each was **executed directly** and its
   output inspected — `50 x 10 = 500`, the `0.90` threshold, and the timing
   footnote all render correctly. The checker itself was validated against a
   synthetic file with known too-few/too-many/`%%`/non-literal cases.
3. **Residual vocabulary, case-insensitive.** `[Tt]ier` and `[Gg]ate` across
   both files, excluding `frontier`/`delegate*`/`gatekeeper`/`ungated`/`Negate`:
   one hit, the pre-existing ordinary-English "not to gate anything" in a sim
   comment (L685). During editing my own `fs_used_grf` comment initially read
   "Gate for …" and was reworded to "Guard for …" so this scan stays clean.
4. **Chunk gating.** Every `eval =` guard is defined before its chunk. Nothing
   after a gated block reads an object it creates — the new note chunks create
   nothing, and the only consumer of CV/LOO objects outside their gates is the
   forest plot's `if (exists("fs_OOB"))` / `if (exists("fs_kfold"))`.
   `timings$kfold`/`$oob` remain pre-seeded `NA_real_`.
5. **Fences and labels.** Multimethod 138 → 158 fences (both even), 69 → 79
   chunks (+10: 2 GRF, 7 CV, 1 LOO); sim 22 → 22. **No duplicate chunk labels**
   in either file. Every chunk parses in both files (0 failures).

### 4.2 Rendered — multimethod

Render completed clean (`Output created: analysis_gbsg_survival_multimethod.html`),
no errors, no warnings in the log.

* **Numeric invariance** — see §1. AFTER vs CONTROL: 11/11 columns bit-identical
  at tolerance 0, `selection_bias` bit-identical for all three engines.
* **No reset warning** — 2 in control, **0** in edited.
* **Timing table reconciles and carries the GRF rows** — §3.2.
* **CV and GRF blocks** — with `run_cv <- FALSE` and `use_grf = FALSE`, neither
  narrative renders and each skipped note appears exactly once. Rendered-only
  counts (embedded source excluded), control → edited: `Diagnostic Tabulations`
  2 → 0, "Which subgroups emerged" 1 → 0, GRF skipped-note 0 → 1, cv/loo
  skipped-notes 2 → 2.

### 4.3 Sim

No render required (comment-only). `diff` against the pre-fix copy shows
**exactly one hunk**: two comment lines replaced by eight. All 11 chunks parse.
Nothing else in the file changed.

---

## 5. The rendered-HTML diff, accounted for

Compared **control vs edited on the same machine**, so timing is the only
expected numeric mover. 55 hunks / 355 added / 214 removed lines. Note the HTML
embeds the full `.qmd` via `code-tools`, so every source edit appears twice —
once as folded code, once in the appended source.

| difference | maps to |
|---|---|
| `stop_threshold = 0.95` → `NULL` + comment (`0.95` ×11 lost) | F3 |
| 2 × `Warning: stop_threshold = 0.95 reset to NULL` gone | F3 |
| parameter table rows gone (`1.25` ×6, `0.90`, `1.0`, `60`, `12`, `2`, `10`) | F2 |
| DINA intro + standalone-DINA prose reworded | F2 |
| 7 `tier`/`tiers` phrases reworded | F4 |
| timing table: 2 new rows, relabelled rows, new footnote | F5 |
| `(no GRF cut)` table row + its narrative gone; skipped note added (`100.00` lost) | F6 |
| CV narrative + `### Diagnostic Tabulations` heading gone (`500`, `1000`, `50`, `0.85`, `0.9` lost) | F7 |
| **section numbers shift** (e.g. `8.1.1` gone) | F7 — the gated heading changes auto-numbering downstream |
| `mr_in_replicates = FALSE` ×2 added to folded code | N1 |
| core-count and DINA-bootstrap comments reworded | N5 |
| "never sets it" → "pins it to `FALSE` at every … call site" | §6.1 |
| `1403` ×4 gained | my new F3 comment cites `forestsearch_main.R:1403` (2 folded + 2 source) |
| `100` ×4 gained | `100 * .resid_s / .total_s` in the new F5 code (2 folded + 2 source) |
| wall-clock values (`289.6`, `155.7`, `681.9`, `26.7`, `8.1`, `6.3`, MR seconds, speedups) | timing, expected — differs between any two renders |

**Nothing unaccounted for.** No estimate, CI, subgroup definition, kappa,
candidate-family size or selection bias appears anywhere in the diff.

---

## 6. Found while doing this, not anticipated by the review or the brief

### 6.1 A callout claim that N1 falsified (fixed)

The "Where MR runs" callout said the document *"never sets"* `mr_in_replicates`.
It already pinned it at three sites before this work and now pins five, so the
sentence was false before N1 and more visibly so after. Reworded to "pins it to
`FALSE` at every one of those call sites — the package default, so the pin is a
no-op that only makes the setting visible — and never sets it to `TRUE`."
Display-only.

### 6.2 The baseline payload was made on different hardware

Covered in §1. Beyond this task, it means **`tolerance = 0` is not a usable
acceptance criterion for this document across machines**. A same-machine control
render is the only version of that test that can discriminate. If cross-machine
reproducibility is wanted, the levers are pinning `n_cores` to a fixed number
rather than deriving it from `detectCores()`, and pinning `grf`'s `num.threads`
(never set anywhere in `R/`, so `grf` defaults to all cores).

### 6.3 The review under-called N4 (see §5 below)

My review noted `covs` had no consumer *in this document* and implied it was
vestigial. It has one — in a different document reading the same `.rds`. Detail
below.

### 6.4 The brief's precondition count was off by one

The brief expects `mr_in_replicates` **twice** in the sim ("one comment, one
argument"). It occurs **once** (the argument, L558); the adjacent seven-line
comment explains it but uses `mr_inference` and `mr_replicates`, never the
literal token. Every other precondition matched exactly (multimethod: 3
occurrences, `fs_param_table()` filtering on `.fs_selection_knobs`; sim:
`RECORD ONLY` provenance block), and both files were verified byte-identical to
the reviewed copies in `dev/review/`. I treated the check's purpose — "confirm
you are on the reviewed content" — as satisfied and proceeded rather than
stopping.

---

## 7. The N4 question for the maintainer

`mr_ok`, `n_sel`, `label` and `covs` are written into every record and persisted
to the `.rds`. **Nothing was removed.** Confirmed: none of the four is read
anywhere in either target document (`sg_def`, by contrast, *is* read — by
`attach_betaHhat()` in `betaHhat_truth.R`).

Looking beyond the two documents changes the picture, and this is the part the
review got wrong:

| column | read anywhere in `quarto/`? | where |
|---|---|---|
| `n_sel` | **yes** | `quarto/GuoHe/guohe_sec52.qmd` (`mean(t7[[i]]$results$n_sel)`); `quarto/smoke-tests/sg_focus_smoke_test.qmd` |
| `covs` | **yes** | `quarto/simulations/gbsg/fs_subgroup_id_sweep_n500.qmd:266` — `strsplit(det$covs, "+", fixed = TRUE)` → subgroup-structure classification and true-recovery rate |
| `mr_ok` | **no** — zero reads anywhere | written by ~180 sim documents, read by none |
| `label` | no read of the *record* column found | (the name is generic; other `$label` hits are unrelated) |

So the recorder is a **shared schema across ~180 sim documents**, and `covs`
does have the "id figures" consumer its comment refers to — just in a different
document reading the same bundle. Removing any of these would fork this file
from the family and change the `.rds` schema.

**The question:** given that, is `mr_ok` worth *surfacing* rather than removing?
Its stated purpose (sim L498) is that "an MR failure never masks a true
detection", but no table reports it, so an MR failure is currently invisible in
the output — `ij_source` is reported, `mr_ok` is not. Three options:

1. **Leave all four as-is** (schema parity with the family; `mr_ok` stays latent).
2. **Leave the schema, add one number** — report `mean(mr_ok)` alongside the
   existing "IJ variance held in X%" footnote in the estimation table. Display-only.
3. **Drop `mr_ok` and `label`** from the record — changes the `.rds` schema and
   diverges from the other ~180 documents. Not recommended.

I have not acted on any of these.

---

## Appendix — what was executed

| check | method |
|---|---|
| F3 no-op | two fits differing only in `stop_threshold`, same session/seed; 9 quantities compared with `identical()` |
| Numeric invariance | full control render of the **unedited** source on this machine; payload compared at tolerance 0, plus ULP quantification against the arm64 baseline |
| Rendered-value invariance | 8 classes of analysis result line compared by checksum between control and edited renders |
| F5 generator | executed standalone: normal path (12 rows, Total last, NA rows preserved, reconciliation) and the unlabelled-timing `stop()` path |
| F7/F6 note chunks | every new `asis` chunk executed and its markdown output inspected |
| N3 mechanism | ragged `rbind` via both `do.call` and `foreach %dofuture%`, plus a body-error control proving `.errorhandling = "remove"` covers only the body |
| Static suite | setup sourcing, forward references, `sprintf` arity, vocabulary scan, gate pairing, fence/label counts |
| N4 | read-vs-write discrimination across all `.qmd` in `quarto/` |

R 4.6.1, forestsearch 0.2.0, quarto render on pop-os (AMD Threadripper PRO
5995WX, 128 cores, 102 workers).
