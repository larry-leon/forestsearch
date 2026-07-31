# CC BRIEF — GBSG application filename migration

```
claude "Read dev/naming/CC_BRIEF_gbsg_filename_migration.md and execute it."
```

Scope is `quarto/applications/gbsg/` only. Use `git mv` throughout so history
follows. Do not re-render anything. Do not touch `quarto/GuoHe/` except for the
two referrer edits listed in §3.

---

## 1. The convention

```
analysis_<dataset>_<outcome_type>_<scope>.qmd
```

`<outcome_type>` uses the package's own vocabulary —
`outcome_type = c("survival", "binary", "continuous", "count")` — so the slot
maps 1:1 to a `forestsearch()` formal. `<scope>` carries `multimethod`,
`maxeff`, `frozen_family`, and so on.

**This is not a new scheme.** Three files in the directory already follow it:
`analysis_gbsg_survival_hr.qmd`, `analysis_gbsg_survival_hrMaxSG-both.qmd`,
`analysis_gbsg_survival_hrMaxSG-pareto.qmd`. This finishes the migration.

Version tokens (`v2_2`, `A`, `new`, `v3a`) and host tokens (`_linux`) come out
entirely. Git holds history; the provenance table in the notebook holds commit
SHA, package version, host and render date.

---

## 2. Renames

| from | to |
|---|---|
| `analysis_gbsg_cox_multimethod_psi_v2_2new.qmd` | `analysis_gbsg_survival_multimethod.qmd` |
| `analysis_gbsg_cox_maxeff.qmd` | `analysis_gbsg_survival_maxeff.qmd` |
| `analysis_gbsg_frozen_family.qmd` | `analysis_gbsg_survival_frozen_family.qmd` |
| `analysis_gbsg_cox_maxeff_rerun1_mrconfirm.qmd` | `analysis_gbsg_survival_maxeff_mrconfirm.qmd` |
| `analysis_gbsg_cox_maxeff_rerun1_mrconfirm.html` | `analysis_gbsg_survival_maxeff_mrconfirm.html` |

The last pair arrived on the branch at `784aa7f4` ("maxeff_rerun1 GBSG") after
this brief was first drafted. `cox` and `rerun1` both come out: the outcome slot
takes `survival`, and `rerun1` is an attempt token that git already records.
`mrconfirm` stays — it is a real scope distinction (the document carries an
independent MR confirmation fit that the plain maxeff analysis does not).

Rename the `.html` alongside the `.qmd` here, unlike the baseline renders in
§5: this render has no role in `REPLICATION_FINDINGS.md`, so nothing references
it by name.

Payload written by the first, renamed to match its notebook:

| from | to |
|---|---|
| `gbsg_table2new_payload.rds` | `gbsg_survival_multimethod_payload.rds` |

Update the `.payload_file` literal inside
`analysis_gbsg_survival_multimethod.qmd` in the same commit as the `.rds`
rename. Writer and reader must move together.

**Already done in the supplied file.** The revised
`analysis_gbsg_survival_multimethod.qmd` provided with this brief already
carries the new payload literal, so `git mv` the `.rds` and drop the revised
`.qmd` into place — do not edit the literal a second time.

---

## 3. Referrers to update

Verified by exact-token grep (`basename` + extension), so these are real, not
substring artifacts:

**`analysis_gbsg_cox_multimethod_psi_v2_2new`** — 2 source referrers:

* `quarto/GuoHe/eval_consistency_split_vs_resample.qmd`
* `quarto/GuoHe/run_guohe_gbsg_maxeff.R`

(A third hit is `eval_consistency_split_vs_resample.html`, a render — it
regenerates and needs no edit.)

**`analysis_gbsg_cox_maxeff`** and **`analysis_gbsg_frozen_family`** — the
referrers are `quarto/GuoHe/` briefs (`guohe_promotion_CC_BRIEF.md`,
`gbsg_maxeff_guohe_CC_BRIEF.md`, `frozen_family_CC_BRIEF.md`) plus the GuoHe
directory's **own separate copies** of those two analyses.

> The GuoHe copies are separate investigations, not duplicates to reconcile.
> Leave `quarto/GuoHe/analysis_gbsg_cox_maxeff.qmd` and
> `quarto/GuoHe/analysis_gbsg_frozen_family.qmd` where they are and named as
> they are. Update only the prose references in the three `.md` briefs, so they
> point at the renamed `applications/gbsg/` files where that is what they meant.
> Where a brief means the GuoHe copy, leave it. **If a reference is ambiguous,
> stop and ask rather than guessing.**

---

## 4. Out-of-repo coupling — the one that fails silently

The manuscript fragment `_gbsg_section6_1_v5_8.qmd` reads the payload by literal
name at line 21:

```r
P <- readRDS(file.path(.pdir, "gbsg_table2new_payload.rds"))
```

It lives **outside this repo**, so no grep here will find it. Renaming the
payload without updating that line gives a `readRDS` file-not-found at render.
Flag it in the findings; do not attempt to edit files outside the repo.

---

## 5. Retire

These four have **zero referrers** anywhere in the repo, and each writes a
payload referenced only by itself:

| notebook | its payload |
|---|---|
| `analysis_gbsg_cox_multimethod_psi_v2_2.qmd` | `gbsg_table2_payload.rds` |
| `analysis_gbsg_cox_multimethod_psi_v2_2A.qmd` | `gbsg_table2A_payload.rds` |
| `analysis_gbsg_cox_multimethod_psi_v2_2A_linux.qmd` | `gbsg_table2Alinux_payload.rds` |
| `analysis_gbsg_cox_multimethod_psi_v3a.qmd` | `gbsg_table_payload.rds` |

`git rm` the `.qmd` and `.rds` for all four.

### Keep the two baseline renders

**Do not delete `analysis_gbsg_cox_multimethod_psi_v2_2A_linux.html` or
`analysis_gbsg_cox_multimethod_psi_v2_2new.html`.** They are baselines A and B
in `dev/replication-check/REPLICATION_FINDINGS.md`, and renaming them would
break those references.

Note the asymmetry: the `_v2_2A_linux.qmd` being removed is *not* the source
that produced the `.html` beside it — the terminology purge rewrote it, 81 lines
different. A faithful reconstruction is preserved at
`dev/replication-check/legacy_v2_2A_reconstructed.qmd`. Removing the misleading
`.qmd` while keeping the `.html` is the correct outcome; say so in the commit
message.

The other two `.html` renders (`_v2_2.html`, `_v2_2A.html`, `_v3a.html`) have no
such role — remove them with their sources.

---

## 5b. Conventions applied to the revised multimethod notebook

The revised `analysis_gbsg_survival_multimethod.qmd` supplied with this brief
already carries the conventions lifted from
`analysis_gbsg_cox_maxeff_rerun1_mrconfirm.qmd`. Do not re-apply them:

* header comment block (what the document is, what it adds and drops relative
  to the maxeff sibling, and a VENUE note on installation and machine);
* per-method parallelism notes in `parallel-setup` (what inherits a plan, that
  MR is a single `crossprod` and must not be wrapped in a future, that CV/LOO
  pin their inner refits sequential);
* toggle documentation for `run_cv`, `run_loo` and `run_mr`, each stating what
  TRUE and FALSE do;
* fit-integrity guards after all four fits;
* parameter-provenance tables per fit plus a consolidated appendix, read from
  `fit$args_call_all`;
* the corrected `stop_threshold` annotation (`default = pconsistency.threshold`,
  0.90 here, not 0.95).

Two conventions were checked and found **already satisfied**: `plan("sequential")`
follows every parallel block (7 occurrences), and there are no positional reads
of result tables.

One was deliberately **adapted rather than copied**. The maxeff document asserts
`stopifnot(!is.null(fs$sg.harm))`. That would break the multimethod document,
which handles a NULL subgroup as a legitimate outcome via its `*-no-subgroup`
chunks. The guard used instead asserts fit integrity:
`stopifnot(!is.null(fs), !is.null(fs$args_call_all))`.

## 5c. Conventions applied to the maxeff and frozen-family notebooks

`analysis_gbsg_survival_maxeff.qmd` and
`analysis_gbsg_survival_frozen_family.qmd` are supplied with this brief, already
renamed and revised. Drop them into place rather than `git mv`-ing the originals.

Both were audited first, and most conventions were **already present**: header
comment block, `bibliography: []`, skip-note chunks for `run_cv_kfold` and
`run_loo`, `plan("sequential")` after every parallel block (4 each), and
`stopifnot()` guards. Nothing there was touched.

Added to each:

* a parameter-provenance table after the `forestsearch()` fit, reading
  `fs$args_call_all`;
* TRUE/FALSE documentation for the `run_loo` and `run_cv_kfold` toggles.

Note the contrast with the multimethod document: here
`stopifnot(!is.null(fs$sg.harm))` is **correct and was left as-is**, because
`sg_focus = "maxeff"` structurally guarantees a selected subgroup. Only the
multimethod notebook, which handles a NULL subgroup as a real outcome, needed
the weaker fit-integrity form.

### Known duplication

The provenance helper (`.fs_fmt_val`, `.fs_fmt_default`,
`fs_param_provenance`) now appears in **three** notebooks. That is deliberate --
each document must render standalone, and a shared `source()`d file would add a
path dependency that breaks silently. The right long-term home is the package
itself, as an exported utility taking any fitted object. Flag it; do not
refactor it here.

## 6. Verify

```bash
git status                      # only the intended paths
grep -rn "psi_v2_2\|psi_v3a\|table2new_payload\|table2A_payload" \
     --include=*.qmd --include=*.R --include=*.md . | grep -v dev/replication-check
```

The second command should return nothing outside `dev/replication-check/`
(which legitimately records the old names as historical baselines).

Do not re-render. Confirming the renamed notebook still renders is a separate
task and a long one.

---

## 7. Report, do not act

Two observations to record in the findings without acting on them:

1. **The already-conforming files carry stale focus vocabulary.**
   `analysis_gbsg_survival_hrMaxSG-both.qmd` and `-pareto.qmd` use `hrMaxSG`,
   whose GLM-neutral alias is `effMaxSG`. They also use hyphens where the rest
   of the directory uses underscores. Both are cosmetic; flag, do not change.

2. **There is a much larger legacy tail.** Roughly twenty files under an older
   `gbsg_analysis_*` prefix, several marked `_working`, `_working2`,
   `_checking`, `_check`, `_new`, `_new-500`, plus
   `mainvignette_gbsg-analysis.Rmd` and a stray `.rmarkdown`. Out of scope
   here. List them in the findings so the maintainer can decide what is live.

---

## 8. Deliverable

`dev/naming/FILENAME_MIGRATION.md`: the rename table as executed, referrer edits
made, files retired, the out-of-repo coupling from §4, and the two §7
observations.
