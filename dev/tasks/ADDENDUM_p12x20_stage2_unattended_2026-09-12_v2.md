# ADDENDUM v2 — Part A Stage 2, unattended execution

- **Date:** 2026-09-12
- **Supersedes** the v1 addendum, which re-asserted what Stage 0 and Stage 1
  had already established. v2 keeps two preconditions and drops the rest.
- **Amends:** `TASK_p12x20_partA_2026-09-12_v2.md` §5–§8. Everything not restated
  there is unchanged. This document does not alter the cells, the knobs, the
  seeds or the gates.
- **Applies only after** Stage 1 has reported and Larry has given the compute go.

Stage 2 runs with nobody watching. The terminal, the CC session and the SSH
connection are all assumed to die at some point during the run; none of them may
be load-bearing.

---

## 1. Preconditions

Two, both because their outcome changes what happens next. Everything Stage 0
and Stage 1 already asserted stands; do not re-assert it.

- **Worker reap.** The Gate 1 probe spawned workers. `pkill -f workRSOCK`
  (ignore rc=1), then assert `pgrep -fc workRSOCK` returns zero.
- **Disk headroom.** Take the per-cell payload size from the Gate 1 probe and
  assert available space on the repo filesystem exceeds 9 x that plus 20 GB.
  Report both numbers. Running out mid-campaign is silent and unrecoverable.

## 2. The runner

- Write `quarto/simulations/gbsg_020/scripts_p12x20/run_p12x20.sh` and **commit it
  before launching**. The run must be reproducible from the repo, not from a
  scratchpad.
- **Transplant, do not author.** The script calls the committed `campaign.sh`
  driver with the knobs fixed in Stage 1 §1b. It does not reimplement the run,
  the combine, or the gates.
- Worker count and thread pinning are the values proposed in the Gate 1 report
  and approved with the go. State them at the top of the script as named
  variables.
- Cells run **sequentially**, one at a time. Workers are already saturated inside
  a cell; running cells concurrently only adds contention and makes the
  per-replicate distribution uninterpretable.
- Per cell, in order: batch 1 (`sim_id` 1–1000) → batch 2 (`sim_id` 1001–2000) →
  combine → combine assertions (§5 of the task) → Gate 2 → same-draws both
  directions → write the cell report → commit → heartbeat.
- Launch detached: `setsid nohup ./run_p12x20.sh > logs/p12x20_runner.log 2>&1 &`.
  Report the PID and the log path.
- Per-cell logs go to `logs/p12x20_<cell>.log`. Logs are scratch and are not
  committed; the reports are the record.

## 3. Resume semantics

- A cell is **done** when its combined payload and its gate record are both
  committed. The script checks that condition and skips done cells.
- Restarting the script after any interruption therefore resumes rather than
  re-runs. Never re-run a committed cell.
- Before any restart, repeat the worker reap from §1.
- An interrupted cell restarts from its first batch. Partial batch output is
  discarded, not salvaged.

## 4. Commit-as-you-go

- Each cell commits on completion, with explicit named paths: payload, cell
  report, gate record, heartbeat.
- A crash loses at most the in-flight cell.
- Commit messages carry the cell identifier, e.g. `p12x20 A4: HR 1.75 n 500`.
- Never `git add -A`. Never stage the untracked directories present on this
  machine. CC never pushes.

## 5. Heartbeat

- Append one timestamped line to
  `quarto/simulations/gbsg_020/LOG_p12x20_progress.txt` at every cell boundary:
  UTC timestamp, cell, event (`start` / `done` / `halt`), elapsed minutes, and
  the commit SHA for `done`.
- Committed at each cell boundary, so the file is readable from GitHub Desktop
  without attaching to the session.

## 6. Failure policy

Applying the standing rule, not asking.

- **Halt the campaign** on anything that makes results untrustworthy: a Gate 2
  identity failure, a same-draws mismatch in either direction, a combine
  assertion failure, or a payload over the 100 MB hard stop.
  - Write `HALT_p12x20.md` naming the cell, the assertion and the observed
    values; commit it; append a `halt` heartbeat; exit non-zero.
  - Completed cells stay committed. Do not roll back.
- **Continue** on anything documentary: a missing line number, an absent
  citation, a render glitch, a metadata field not found. Append to
  `OPEN_ITEMS_p12x20.md` and carry on.
- Flag any single artifact over 50 MB in the cell report and continue.

## 7. Closeout

- After the ninth cell, regenerate `STATUS_p12x20.md` per §8 of the task and
  commit it.
- **Do not regenerate `current_status.md`** and do not run
  `check_current_status.sh`. Both belong to the merge task on
  `feature/glm-extension`.
- Report the full commit range for Larry to push. Do not push.

## 8. Reading the state without CC

For Larry, at the repo root:

```
tail -40 quarto/simulations/gbsg_020/LOG_p12x20_progress.txt
git log --oneline -20
pgrep -fc workRSOCK
ls quarto/simulations/gbsg_020/HALT_p12x20.md 2>/dev/null && echo HALTED
```

## 9. Report on return

- Bullet form, one item per bullet.
- Per cell: wall-clock, replicates combined, gate counts (run / passed / failed),
  payload path and size, commit SHA.
- Realized total wall against the Gate 1 projected range.
- Contents of `OPEN_ITEMS_p12x20.md` if non-empty.
- Numbers verbatim. No interpretation and no recommendation.
