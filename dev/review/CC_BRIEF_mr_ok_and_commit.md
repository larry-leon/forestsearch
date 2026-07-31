# CC BRIEF — surface `mr_ok`, commit the review record, close open items

```
claude "Read dev/review/CC_BRIEF_mr_ok_and_commit.md and execute it."
```

Follow-on to `dev/review/QMD_FIXES_APPLIED.md` (committed as `1ec0bc1`). Three
small tasks and one question answered.

---

## 1. Answer to the N4 question: surface `mr_ok`, remove nothing

Your revised finding is the right one, and it changes the disposition. Of the
four columns, three are settled and need no action:

| column | status |
|---|---|
| `n_sel` | read by `quarto/GuoHe/guohe_sec52.qmd` and `quarto/smoke-tests/sg_focus_smoke_test.qmd` — **keep, no change** |
| `covs` | read by `quarto/simulations/gbsg/fs_subgroup_id_sweep_n500.qmd:266` — **keep, no change** |
| `label` | no reader found, but it is part of a schema shared across ~180 documents — **keep, no change** |

`mr_ok` is different, and not because it is unread. Its stated purpose (sim
L499) is that *"an MR failure never masks a true detection"* — it exists to
separate the two events. It succeeds at that internally: `detected` and `mr_ok`
are recorded independently at L508. But nothing reports it, so **an MR failure
is currently invisible in every rendered summary.** A run in which MR failed on
a third of the detected replicates would look identical to one where it
succeeded on all of them.

### The change

The summary footnote at sim **L1060** already reports the neighbouring
quantities:

> `Detection rate **%.0f%%** of %d sims; MR flags harm in **%.0f%%** of detected
> reps; IJ variance held (not fallback) in **%.0f%%**. ...`

Add MR availability to that same sentence — the share of **detected** replicates
for which MR returned an object, i.e. `mean(mr_ok[detected == 1])`. Place it
before the harm-flag clause, since harm-flagging is conditional on MR having run
at all, and the reader needs to know the denominator before the rate.

Match the existing style: same `%.0f%%` formatting, same bolding, same
detected-reps denominator. Check whether the second footnote at **L1166** should
carry it too, and treat both consistently.

Two things to get right:

* **Denominator.** `mr_ok` is `0L` in the all-NA failure record (L399), so a
  replicate that never detected also carries `mr_ok = 0`. Conditioning on
  `detected == 1` is what makes the number mean "MR failed" rather than "nothing
  was found." Verify that against a real batch `.rds` rather than assuming it.
* **Degenerate case.** With zero detections the mean is `NaN`. Confirm the
  surrounding footnote already guards that case; if it does not, guard the new
  clause the same way the existing clauses are guarded.

This is display-only: it adds a computed figure to a footnote and touches no
recorded value. Confirm no `.rds` schema change — the column set written must be
identical.

### Verification

The sim has no committed baseline render to compare against, so verify by
computation instead. Load the existing on-disk batch, compute the new figure
directly, and confirm the footnote as rendered would print that same number.
Report the value it yields on the real batch: if MR succeeded on every detected
replicate it reads 100%, which is the expected and uninteresting case — say so
rather than treating it as a null result.

---

## 2. Commit the review record

`dev/review/QMD_FIXES_APPLIED.md` and `dev/review/CC_BRIEF_qmd_fixes.md` are
untracked. Commit both.

`QMD_FIXES_APPLIED.md` carries the hardware finding, which is worth having in
the repo rather than only in a chat log: **`tolerance = 0` is not a usable
acceptance criterion for this document across machines**, and a same-machine
control render is the only form of that test that works. Any future replication
brief needs that constraint, and the report is where it is written down.

Also commit `dev/review/QMD_BUG_REVIEW.md` and
`dev/review/CC_BRIEF_qmd_bug_review.md` if they are still untracked, so the
review and its response sit together.

**Do not commit `quarto/methodology/selection_criteria.qmd`** — it is the
maintainer's and is still in draft.

---

## 3. Add the hardware constraint where it will be found

The finding in §2 is only useful if someone hits it again. Add a short note —
three or four sentences — to `dev/replication-check/REPLICATION_FINDINGS.md`,
which is where a future replication attempt will start:

* the payload records `.machine_info()`; check it before comparing payloads;
* cross-architecture renders are not bit-comparable, and the worker count feeds
  the candidate-search batch size, so summation order changes;
* the valid test is a same-machine control render of the unedited source;
* point at `dev/review/QMD_FIXES_APPLIED.md` §1 for the worked case.

Do not restate the whole analysis. A pointer plus the constraint is enough.

---

## 4. Deliverable

No separate document. Report:

1. The `mr_ok` value on the real batch, and where the footnote change landed
   (one site or two).
2. Confirmation the `.rds` schema is unchanged.
3. What was committed, and the SHA.
4. Anything that did not go as described here.

Commit the sim `.qmd`, the review documents, and the `REPLICATION_FINDINGS.md`
note together. Report and stop.
