# DISPOSITIONS — declcal consumers v3, resuming after the Gate 2 stop (2026-09-24)

**Resumes** `dev/tasks/TASK_declcal_consumers_2026-09-24_v3.md` (`dd3b8331`) from the working tree, where the fixes
made before the stop are unstaged. This document gives Larry's answers to the four decisions in §5 of
`dev/reports/REPORT_declcal_consumers_2026-09-24_v3.md`. Everything in v3 that is not changed here stands,
including the rule of no `R/` change. **Hard abort: 1 h** for this resumption.

---

## 0. First action

1. Copy this file to `dev/tasks/DISPOSITIONS_declcal_consumers_2026-09-24.md` and `git add` that path only.
   Commit it with the message `docs(tasks): declcal consumers v3 dispositions (2026-09-24)`.
2. One background shell was still running at the stop. Identify it, and if it belongs to this task, stop it.
   Nothing it wrote may be staged.

---

## 1. The six scripts: archived, not fixed

These are the six scripts the report lists as "not fixed: §5.1":

- `fdr_mr_inference.R`
- `fdr_family_multiplier.R`
- `reselection_check.R`
- `mr_mechanism_A1.R`
- `mr_mechanism_probe_superset.R`
- `mr_mechanism_probe_restrict.R`

- **Reclassify all six as archived findings, class (a).** They verify pre-alignment behaviour of closed work,
  and rewriting them would change what their recorded outputs mean. The two scripts that `source()`
  `mr_mechanism_A1.R` (`…_earlystop.R` and `…_b1_price.R`) are treated the same way.
- **Add a header line only; make no code change.** On these eight files, use this wording: "Verifies
  pre-alignment (exact-cutoff) behaviour, before the 2026-09-23 alignment (96f84ad8, 7713942e); not valid
  against the aligned package. See dev/tasks/TASK_declcal_consumers_2026-09-24_v3.md."
- **Add no `R/` helper.** In the report's out-of-scope entry for the OC gate, record two things:
  - the forward map from (p\*, digits) to the rounded z cutoff that these scripts would need (the suggested
    `.fs_z_eff(p, digits)`) is the same map the OC gate lacks;
  - the third OC site found, `R/fs_oc_grid.R:383`.

  This is a record entry, not a task.

---

## 2. The three `quarto/resampling` Pcons-definition lines: new Gate 2 class (d)

Class (d) covers a line that computes the consistency proportion from T (the Pcons definition) and does not
compare it with p\* or any threshold. List these lines in the report. Make no change and add no header.

---

## 3. The dropped exact-cutoff columns: confirmed

The drivers no longer record `alpha_FW_hat_1645`, `declared_conv_exact`, `n_band` or `fw_1645_<c0>`. This
follows Larry's decision that the exact-threshold version is replaced outright. The committed payloads keep those
columns as records.

**Addition.** The findings scripts still read those columns. On a payload written by the new drivers, a missing
column reaches `sprintf()` as a zero-length or NA value, so rows drop out or print NA without any error. At the
top of each findings script, assert that every column it reads is present, and stop with a message naming any
missing column. Demonstrate this once on a scratch copy of one payload with one such column removed.

---

## 4. The two `qnorm(0.95)` comparisons in the findings scripts: use the executed cutoff

`qnorm(0.95)` is the exact cutoff for p\* = 0.90 written another way.

- Replace both comparisons with the executed cutoff as the payload recorded it (`z_exec()` or the stored
  `z_pstar`). Relabel any column whose name refers to 1.6449.
- Then grep the nine fixed consumers for `qnorm(0.95)`, `1.6449` and `1.644854`. Treat any **threshold** use
  the same way. A threshold use is one compared with T, M\* or κ̂.
- Leave interval critical values alone. For example, the Wilson intervals at alpha 0.10 use the same number
  legitimately.
- In the report, list each hit and what was done with it.

---

## 5. Finish

1. **Headers.** Add the header line to every class (a) file:
   - the inventory's archived files;
   - the committed outputs of the fixed findings scripts (per v3 §2b);
   - the eight files in §1.

   In `.qmd` and `.Rmd` files, place the line immediately after the YAML front matter, never above it.
2. **Gate 1.** Re-run it only on the files this document changes. Do not re-run Gate 3 unless the test file
   changes again.
3. **Gate 2.** Run the final Gate 2 under classes (a)–(d).
4. **Test 1, recorded only.** Record whether `test-declaration-calibration.R` runs under `skip_on_cran()`, and
   the wall clock of the band fit.
5. **Report.** Update it as follows:
   - keep the stop section as history;
   - add a Dispositions section that points to this file;
   - change the title and the outcome line to the final state;
   - make the Gate 2 sentence of record read: "Outside `R/`, no run script derives the threshold; remaining
     hits: n(a) archived, n(b) test oracles, n(c) OC pin, n(d) Pcons definitions."
6. **Commits.** Use explicit paths, in this order: run-script fixes; archived-file header lines; test fix;
   report.

---

## POST-CONDITIONS (machine-checkable)

1. v3's post-conditions 1–8 hold, with Gate 2 read under classes (a)–(d).
2. For each of the eight files in §1, `git diff dd3b8331 -- <file>` shows only the added header line.
3. Every findings script contains the column check, and the check is shown to stop on the scratch payload.
4. No threshold use of `qnorm(0.95)`, `1.6449` or `1.644854` remains in the nine fixed consumers.
5. Nothing is staged beyond the named paths, and the pre-existing untracked files remain unstaged.
