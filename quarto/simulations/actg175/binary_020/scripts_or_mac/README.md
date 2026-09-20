# scripts_or_mac — how the ten GRF/DINA cells ran on this Mac (2026-09-19)

Committed as a record of what ran, not as a toolkit, in the manner of
`../../continuous/scripts_mdf1/`. Cited by `REPORT_binary_stage2_2026-09-19.md`
(`dev/tasks/TASK_binary_stage2_v3_2026-09-19.md`, Step 2).

These files lived in `~/.orsg_shim_2026-09-19` during the pass — outside the repository, prepended to
`PATH` — and are copied here verbatim afterwards. **`run_or.sh` and `gate2.R` were not edited**: every
pop-os-specific element of the runner was resolved through this directory instead (Step 0d).

- `timeout` — Perl `alarm` wrapper standing in for GNU `timeout`, which macOS does not ship and for which
  no Homebrew `coreutils`/`gtimeout` was installed. The child runs in its own process group; on expiry the
  group gets TERM, a 10 s grace, then KILL, and the wrapper exits **124** as GNU `timeout` does, otherwise
  propagating the child's status (128+N for a signal death). A stricter relative of
  `../../continuous/scripts_mdf1/tmo.sh`, this Mac's earlier wrapper, which killed a bare pid with no grace.
  It also carries the mid-pass supplement: a `quarto render` belonging to DINA is re-capped from 52,339 s to
  **157,017 s**, identified by the render's own `FS_OR_METHOD` or its `_ordina_` output stem. The
  `quarto render` requirement keeps any other use of the wrapper on the seconds it was given.
- `pgrep` — macOS `pgrep(1)` has no `-c`, which `run_or.sh`'s preflight uses (`pgrep -fc '[w]orkRSOCK'`).
  This counts the matching lines, mirrors GNU's exit status, and delegates every other invocation to
  `/usr/bin/pgrep`. Without it the preflight fails before any render.
- `drive_stage2_mac.sh` — the sequencing driver. It adds nothing to the campaign's logic: it invokes the
  **unchanged** `../scripts_or/run_or.sh` once per cell, in the task's table order, under the Step 0d
  environment (shim `PATH`, W = 13, the scaled caps, `ORSG_HOST`, `VECLIB_MAXIMUM_THREADS=1`). Per-cell
  invocation is what produces the task's order — `run_or.sh` iterates its own `CELLS_ALL`, which puts
  `or150_n2000` before `or100_n2000`, and `ORSG_CELLS` subsets without reordering — and it is how "a cell
  that fails writes its halt record and the run continues to the next cell" is honoured, since `halt()`
  exits the invocation and the preflight refuses to start while `HALT_or.md` exists. No cell halted in this
  pass, so that path did not run. The driver keeps the true cumulative wall, `run_or.sh`'s `CUM` having
  reset at each invocation.
- `pairing_proof.R` — the per-replicate cross-machine check behind §4 of the report: a Mac `orgrf`/`ordina`
  bundle against the committed pop-os `orfs` bundle at the same design point and n, over the data-only
  columns (`seed`, `n_true`, the 8 oracle columns), tolerance-based at 1e-8 rather than `identical()`
  because the machines differ in BLAS and R version. Prints the matched fraction per column; if the row
  sets fail to line up it says the pairing is by construction and not verified, rather than reporting a
  fraction. Usage: `Rscript pairing_proof.R <qmd-dir> <orgrf|ordina> <or075|or100|or150> <n>`.
- `STEP0_record.txt` — the Step 0e launch record printed before the pass: machine, OS, R, BLAS/LAPACK, the
  installed build and C, the design of record, W, the caps, each 0d adaptation, and the template lines
  showing the seeds do not depend on the worker count.
- `timeout_applied.log` — one line per render, with the cap passed and the cap applied (10 GRF at 52,339 s,
  10 DINA at 157,017 s, plus two pre-launch test lines). Cell 1's batch render predates the log.
- `pairing_proof_all.txt` — `pairing_proof.R`'s output for all ten cells.
