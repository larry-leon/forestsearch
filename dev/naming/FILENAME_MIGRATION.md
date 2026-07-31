# FILENAME_MIGRATION — GBSG application filename migration

Executed from `dev/naming/CC_BRIEF_gbsg_filename_migration.md` on branch
`feature/mr-in-replicates`, base commit `4b642b9`. Scope: `quarto/applications/gbsg/`,
the two referrer edits in `quarto/GuoHe/` named in the brief's §3, and — by
maintainer decision after the brief's premises were checked — three path fixes
in `dev/replication-check/R/`.

Convention now applied: `analysis_<dataset>_<outcome_type>_<scope>.qmd`, with
`<outcome_type>` drawn from the package's own `outcome_type` vocabulary.

Nothing was re-rendered.

---

## 1. Renames executed

All via `git mv`, so rename detection holds (`git status` shows `R`, not `D`+`A`).

### From the brief

| from | to |
|---|---|
| `analysis_gbsg_cox_multimethod_psi_v2_2new.qmd` | `analysis_gbsg_survival_multimethod.qmd` |
| `analysis_gbsg_cox_maxeff.qmd` | `analysis_gbsg_survival_maxeff.qmd` |
| `analysis_gbsg_frozen_family.qmd` | `analysis_gbsg_survival_frozen_family.qmd` |
| `analysis_gbsg_cox_maxeff_rerun1_mrconfirm.qmd` | `analysis_gbsg_survival_maxeff_mrconfirm.qmd` |
| `analysis_gbsg_cox_maxeff_rerun1_mrconfirm.html` | `analysis_gbsg_survival_maxeff_mrconfirm.html` |
| `gbsg_table2new_payload.rds` | `gbsg_survival_multimethod_payload.rds` |

### Added by maintainer decision — orphaned renders repaired

The brief renamed the `.html` only for the mrconfirm pair, leaving three renders
under old names while their sources had moved. All three were repointed:

| from | to |
|---|---|
| `analysis_gbsg_cox_maxeff.html` | `analysis_gbsg_survival_maxeff.html` |
| `analysis_gbsg_frozen_family.html` | `analysis_gbsg_survival_frozen_family.html` |
| `analysis_gbsg_cox_maxeff_rerun1.html` | `analysis_gbsg_survival_maxeff_rerun1.html` |

The two baseline renders (`_v2_2A_linux.html`, `_v2_2new.html`) keep their old
names — see §3.

### Content of the renamed notebooks

The three revised notebooks supplied with the brief were copied over the renamed
paths after the `git mv`, so history follows *and* the revised content landed.
No convention from §5b/§5c was re-applied, and the `.payload_file` literal was
not edited a second time — the supplied multimethod notebook already carries
`gbsg_survival_multimethod_payload.rds` (line 2343). Writer and reader moved in
the same commit.

**Stale self-names fixed.** Line 12 of three notebooks is a header comment
naming the file it lives in, and each was stale after the rename. Corrected:

| file | was | now |
|---|---|---|
| `analysis_gbsg_survival_maxeff.qmd` | `analysis_gbsg_cox_maxeff.qmd` | `analysis_gbsg_survival_maxeff.qmd` |
| `analysis_gbsg_survival_frozen_family.qmd` | `analysis_gbsg_frozen_family.qmd` | `analysis_gbsg_survival_frozen_family.qmd` |
| `analysis_gbsg_survival_maxeff_mrconfirm.qmd` | `analysis_gbsg_cox_maxeff.qmd` | `analysis_gbsg_survival_maxeff_mrconfirm.qmd` |

The mrconfirm one was already wrong before the migration — it named the maxeff
sibling, not itself. Otherwise `analysis_gbsg_survival_maxeff_mrconfirm.qmd` was
renamed only; its body was not touched.

---

## 2. Referrer edits made

### In `quarto/GuoHe/` — the two from brief §3

Both are comment lines naming the notebook whose configuration they mirror:

| file | line | change |
|---|---|---|
| `quarto/GuoHe/eval_consistency_split_vs_resample.qmd` | 62 | `analysis_gbsg_cox_multimethod_psi_v2_2new.qmd` → `analysis_gbsg_survival_multimethod.qmd` |
| `quarto/GuoHe/run_guohe_gbsg_maxeff.R` | 50 | same |

### In `dev/replication-check/R/` — three live paths, by maintainer decision

See §5 for why these were held for a decision rather than done up front.

| file | line | change |
|---|---|---|
| `01_precondition_rng.R` | 16 | `QMD` → `quarto/applications/gbsg/analysis_gbsg_survival_multimethod.qmd` |
| `04_compare.R` | 19 | `new$payload` → `gbsg_survival_multimethod_payload.rds` |
| `04_compare.R` | 20 | `new$html` → `analysis_gbsg_survival_multimethod.html` |
| `04_compare.R` | 21 | `B$payload` → `gbsg_survival_multimethod_payload.rds` |

Lines 19 and 20 are the *fresh-render* outputs in `NEW_DIR`; a re-render of the
renamed notebook now emits both under the new names. Lines 22 and 24 (the
baseline-B and baseline-A `.html`s) were deliberately left alone: those renders
keep their old names.

A note was added to the top of `dev/replication-check/REPLICATION_FINDINGS.md`
recording that these paths were repointed post-hoc, and that
`analysis_gbsg_survival_multimethod.qmd` is **not** the file that produced those
findings — it additionally carries provenance tables, fit guards, and a
corrected `stop_threshold` comment, all intended to be inert.

### The three GuoHe `.md` briefs needed no edit

Every reference to `analysis_gbsg_cox_maxeff.qmd` and
`analysis_gbsg_frozen_family.qmd` in `guohe_promotion_CC_BRIEF.md`,
`gbsg_maxeff_guohe_CC_BRIEF.md` and `frozen_family_CC_BRIEF.md` is explicitly
path-qualified as `quarto/GuoHe/...`, or is a continuation line inside a
`quarto/GuoHe/` context (`guohe_promotion_CC_BRIEF.md:77-78`). All six hits mean
the GuoHe copies, which the brief says to leave in place and named as they are.
None was ambiguous, so nothing was changed and nothing needed asking.

`quarto/GuoHe/analysis_gbsg_cox_maxeff.qmd` and
`quarto/GuoHe/analysis_gbsg_frozen_family.qmd` are untouched.

---

## 3. Files retired

`git rm`. All had zero referrers in tracked source.

### From brief §5

| notebook | payload | render |
|---|---|---|
| `analysis_gbsg_cox_multimethod_psi_v2_2.qmd` | `gbsg_table2_payload.rds` | `_v2_2.html` removed |
| `analysis_gbsg_cox_multimethod_psi_v2_2A.qmd` | `gbsg_table2A_payload.rds` | `_v2_2A.html` removed |
| `analysis_gbsg_cox_multimethod_psi_v2_2A_linux.qmd` | **kept — see §5** | `_v2_2A_linux.html` **kept** |
| `analysis_gbsg_cox_multimethod_psi_v3a.qmd` | `gbsg_table_payload.rds` | `_v3a.html` removed |

### Added by maintainer decision

`analysis_gbsg_cox_maxeff_rerun1.qmd` — verified **byte-identical** to the
pre-rename `analysis_gbsg_cox_maxeff.qmd` against `HEAD`, so a pure source
duplicate. Its render survives as `analysis_gbsg_survival_maxeff_rerun1.html`.

`quarto/GuoHe/temp/` — deleted (plain `rm -rf`; the directory was untracked, so
nothing entered git history). It held a 2026-07-17/21 snapshot of twelve files,
not the two stale migration copies alone. Each was checked against its live
counterpart before deleting, and all twelve are superseded:

| files | status |
|---|---|
| `guohe_comparator.R`, `guohe_scaling_diagnostic.R`, `guohe_supplement_proofs.qmd`, `patch_spec_sg_focus_maxeff.md`, `test_guohe_algorithm3.qmd`, `test_guohe_comparator.qmd` | byte-identical to the tracked `quarto/GuoHe/` originals |
| `eval_consistency_split_vs_resample.qmd`, `run_guohe_gbsg_maxeff.R`, `guohe_from_forestsearch.R`, `maxeff_tests.rds` | older than their tracked `quarto/GuoHe/` counterparts (the first two are the stale copies of the files edited in §2) |
| `guohe_algorithm3.R`, `guohe_adaptive_r.R` | pre-promotion drafts; both now live in the package as `R/guohe_algorithm3.R` (533 vs 485 lines) and `R/guohe_adaptive_r.R` (320 vs 275) |

The two promoted files were diffed line-by-line against their package versions:
`guohe_algorithm3.R` has **zero** non-comment lines present only in the temp
copy, and `guohe_adaptive_r.R`'s seven are all pre-promotion scaffolding — a
`library(survival)` call superseded by `@importFrom`, a
`exists("guohe_algorithm3")` source-check made moot by packaging, and line
re-wrapping. The `fast` argument survives in both package versions
(`R/guohe_adaptive_r.R:194`, `R/guohe_algorithm3.R:244`). No unique work lost.

### Kept deliberately

* `analysis_gbsg_cox_multimethod_psi_v2_2A_linux.html` — baseline A in
  `dev/replication-check/REPLICATION_FINDINGS.md`.
* `analysis_gbsg_cox_multimethod_psi_v2_2new.html` — baseline B, same.
* `gbsg_table2Alinux_payload.rds` — baseline A's payload; see §5.

The `_v2_2A_linux.qmd` that was removed is **not** the source that produced the
`.html` kept beside it: the terminology purge rewrote it (81 lines different). A
faithful reconstruction is preserved at
`dev/replication-check/legacy_v2_2A_reconstructed.qmd`. Removing the misleading
`.qmd` while keeping the `.html` is the intended asymmetry.

Brief §5 says "the other **two** `.html` renders" and then lists three
(`_v2_2.html`, `_v2_2A.html`, `_v3a.html`). The list was followed; all three were
removed.

---

## 4. Out-of-repo coupling (report only — not acted on)

The manuscript fragment `_gbsg_section6_1_v5_8.qmd` lives outside this repo and
reads the payload by literal name at line 21:

```r
P <- readRDS(file.path(.pdir, "gbsg_table2new_payload.rds"))
```

**This still needs a manual edit** to `gbsg_survival_multimethod_payload.rds`, or
its next render fails with a `readRDS` file-not-found. No grep in this repo would
have caught it.

A reference copy of the fragment sits at `dev/naming/_gbsg_section6_1_v5_8.qmd`
(supplied with the brief, untracked). It was **not** edited — editing the copy
would create a silent divergence from the authoritative out-of-repo original
without fixing anything. Fix the original.

---

## 5. Two brief premises that did not hold

Both concern `dev/replication-check/R/`, which the brief's §6 verification
excludes from the grep on the grounds that it "legitimately records the old
names as historical baselines." That is true of the `.md` findings, which are
prose records. It is **not** true of the two `.R` files there, which are live,
executable, and read these paths at run time. Both were surfaced before acting
and resolved by the maintainer.

**(a) `gbsg_table2Alinux_payload.rds` — kept, not deleted.**

Brief §5 states each of the four retired notebooks "writes a payload referenced
only by itself." For this one that is false:

```r
# dev/replication-check/R/04_compare.R:23
A = list(payload = file.path(GB, "gbsg_table2Alinux_payload.rds"),
         html    = file.path(GB, "analysis_gbsg_cox_multimethod_psi_v2_2A_linux.html")))
```

It is baseline A's *payload*, the exact counterpart of the `_v2_2A_linux.html`
that §5 explicitly protects for the same reason; deleting it would break
`04_compare.R` at line 26 (`readRDS`). **Decision: keep.** The file is tracked
and unchanged. Its `.qmd` was still retired, per the brief.

**(b) Three live paths repointed.**

The notebook and payload renames broke `01_precondition_rng.R:16` and
`04_compare.R:19`/`:21`. **Decision: update all three** (done — see §2), so the
replication harness still runs, with the `REPLICATION_FINDINGS.md` note above
recording that the notebook it now points at is not the one that produced the
findings.

---

## 6. Verification

`git status` shows only intended paths: 9 renames, 9 deletions, 5 modifications
(2 in `quarto/GuoHe/`, 3 in `dev/replication-check/`), and untracked
`dev/naming/`.

The §6 grep

```bash
grep -rn "psi_v2_2\|psi_v3a\|table2new_payload\|table2A_payload" \
     --include=*.qmd --include=*.R --include=*.md .
```

returns nothing in `quarto/applications/gbsg/`, and a grep for
`analysis_gbsg_cox_maxeff|analysis_gbsg_frozen_family` across that directory's
`.qmd` files returns nothing. Remaining hits elsewhere, all outside scope:

* `dev/replication-check/` prose — `REPLICATION_FINDINGS.md`,
  `CC_BRIEF_replication.md`, `CC_BRIEF_replication_v2.md`,
  `v2_2new_rendered_source_prerename.qmd`,
  `legacy_v2_2A_reconstructed.qmd`. The historical record the brief exempts;
  correct as written.
* `dev/replication-check/R/04_compare.R:22` and `:24` — the baseline-B and
  baseline-A `.html` paths, kept under old names by design. Every live path in
  that harness now resolves.
* `dev/naming/CC_BRIEF_gbsg_filename_migration.md`,
  `dev/naming/_gbsg_section6_1_v5_8.qmd`, and this file — the brief, its supplied
  fragment, and the record of the migration.
* `quarto/applications/actg175/` and `quarto/simulations/actg175/` — the ACTG175
  application carries the same `psi_v2_2new` / `psi_v2_2A` / `psi_v3a` and
  `actg175_table2*_payload.rds` vocabulary. A parallel migration applies there
  but was not in scope.

The three staging copies of the revised notebooks were removed from `dev/naming/`
after being dropped into place, leaving the brief, this record, and the
manuscript fragment.

Nothing was re-rendered. Confirming the renamed notebooks still render is a
separate task.

---

## 7. Report-only observations

### 7.1 The already-conforming files carry stale focus vocabulary

`analysis_gbsg_survival_hrMaxSG-both.qmd` and
`analysis_gbsg_survival_hrMaxSG-pareto.qmd` use `hrMaxSG`, whose GLM-neutral
alias is `effMaxSG`, and use hyphens where the rest of the directory uses
underscores. Both cosmetic. Flagged, not changed.

### 7.2 The larger legacy tail

Still tracked in `quarto/applications/gbsg/`, out of scope here. The maintainer
should decide what is live:

**`gbsg_analysis_*` prefix (21 files):**

```
gbsg_analysis.qmd
gbsg_analysis_2sim-approaches.qmd
gbsg_analysis_cox_depth2_DINA.qmd
gbsg_analysis_cox_effMaxSG.qmd
gbsg_analysis_cox_effMinSG.qmd
gbsg_analysis_cox_hrMaxSG_b5k.qmd
gbsg_analysis_cox_maxSG.qmd
gbsg_analysis_cox_maxSG_b5k.qmd
gbsg_analysis_cox_multimethod.qmd      (+ .html)
gbsg_analysis_cox_poisson.qmd
gbsg_analysis_cox_poisson_1.qmd
gbsg_analysis_cox_poisson_check.qmd
gbsg_analysis_cox_poisson_checking-bootstrap.qmd
gbsg_analysis_cox_poisson_checking-bootstrap_new.qmd
gbsg_analysis_cox_poisson_checking-bootstrap_new-500.qmd
gbsg_analysis_cox_poisson_maxSG.qmd
gbsg_analysis_hrMaxSG.qmd
gbsg_analysis_maxSG.qmd
gbsg_analysis_working.qmd
gbsg_analysis_working2.qmd
gbsg_analysis_working_checking.qmd
```

**Other legacy:**

```
gbsg_poisson_simulations.qmd
gbsg_survival_cox_poisson.qmd
gbsg_survival_cox_poisson.rmarkdown     (stray render intermediate)
gbsg_survival_dina_cox.qmd
gbsg_survival_dina_cox_poisson.qmd
mainvignette_gbsg-analysis.Rmd
example_gbsg_pareto.R
```

### 7.3 Provenance-helper duplication

`.fs_fmt_val`, `.fs_fmt_default` and `fs_param_provenance` now appear in three
notebooks (multimethod, maxeff, frozen_family). Deliberate — each document must
render standalone, and a `source()`d shared file would add a path dependency
that breaks silently. The right long-term home is the package itself, as an
exported utility taking any fitted object. Flagged, not refactored.

### 7.4 Header cross-references still point at the old family glob

Line 14 of `analysis_gbsg_survival_maxeff.qmd` and
`analysis_gbsg_survival_maxeff_mrconfirm.qmd` reads "a simplified single-track
version of `analysis_gbsg_cox_multimethod_*.qmd`". The glob no longer matches
anything — that family is now the single `analysis_gbsg_survival_multimethod.qmd`.
Only the line-12 self-names were in scope. Flagged, not changed.
