# TASK — Restore the curated sections of `current_status.md`

- **Date:** 2026-09-13
- **Machine:** Pop!_OS
- **Branch:** `feature/glm-extension`, currently at `ef8bb985`
- **Cause:** the merge task's §5 forbade sourcing from the prior status file, to
  stop a stale pin propagating. It also stripped the curated prose. This restores
  it without reintroducing that risk.

The inventory stays generated. The curated prose becomes a hand-maintained
tracked file that the generator includes verbatim.

---

## 0. Standing constraints

- **No `R/` change.** If anything appears to need one, STOP and report.
- CC never fetches, pulls or pushes. Explicit named paths on every `git add`.
- Re-run no campaign, re-verify no committed cell.
- **Carry the curated text verbatim.** Do not paraphrase it, do not update it, do
  not re-derive any of it from memory, from a report, or from this document.
  A superseded value stays written as superseded; it is not re-quoted, corrected
  or refreshed.

## 1. Assert the starting point

- On `feature/glm-extension`, clean tracked tree, HEAD `ef8bb985`. STOP otherwise.

## 2. Extract

- Read `quarto/simulations/gbsg_020/current_status.md` at `0ab5d1c5`
  (`git show 0ab5d1c5:...`).
- Identify the sections that are **curated** — written by hand, not derivable
  from the directory or from git. Expected: DGM notes, reading conventions,
  superseded and do-not-quote claims, open work.
- Report the section headings you classified as curated, and those you
  classified as generated, before writing anything.
- Anything ambiguous: classify it as curated. A generated section wrongly carried
  is harmless duplication; a curated section wrongly dropped is the failure this
  task exists to fix.

## 3. Write the curated source

- Create `quarto/simulations/gbsg_020/status_curated.md` holding those sections,
  **byte-identical** to their text at `0ab5d1c5`.
- Head it with one line stating it is hand-maintained, that the generator
  includes it verbatim, and that it is edited directly rather than regenerated.
- Commit it with an explicit path.

## 4. Wire it into the generator

- Edit `scripts_dinamr/current_status_regen.R` to read `status_curated.md` and
  emit it verbatim into `current_status.md`, alongside the generated inventory.
- The generator must **fail loudly** if `status_curated.md` is missing — not
  silently emit a file without it.
- Place the curated sections where they sat at `0ab5d1c5` relative to the
  inventory, if that is recoverable; otherwise ahead of the inventory, and say so.

## 5. One addition to the do-not-quote list

Add exactly this line, verbatim, and nothing else of your own composition:

> The `cert20` Stage 1 report's per-cell timing sums field-block seconds into
> `fit_mr_secs`, which already contains them. Those totals are double-counted.
> The `p12x20` Gate 1 report uses `fit_mr_secs` alone.

## 6. Regenerate and assert

- Regenerate `current_status.md`. Commit it alone.
- Run `check_current_status.sh`. The stated pin must equal HEAD at commit time.
- **Idempotence:** regenerate a second time without committing and assert the
  output is unchanged apart from the pin.
- **Fidelity:** assert the curated text in the regenerated `current_status.md` is
  byte-identical to the corresponding text at `0ab5d1c5`. Report the comparison.
- `git status --porcelain --untracked-files=no` empty at the end.
- Any failure: STOP and report. Do not push.

## 7. Report

- Bullet form, one item per bullet.
- Curated vs generated classification, with headings.
- Fidelity and idempotence results.
- `check_current_status.sh` result.
- Commit range for Larry to push.
