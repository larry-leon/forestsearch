# REPORT — MR (field) vs Guo & He: Gate 2 record

Date: 2026-09-05. Task: `dev/tasks/TASK_mr_field_vs_guohe_2026-09-05.md`. Run under Larry's unattended pre-authorization (Gate 1 green: identities, byte-identity, suite, smoke, projection 47 core-h ≪ 12 h ceiling).
Driver: `mr_field_vs_guohe_run.R` at the task settings — 16 cells × 2000 replicates, the replication's own seeds, E1–E5 defaults (MR B = 5000 centred Poisson; field R_out/R_in = 1000/500 Gaussian under seed offset 900000; no M_eff pass, columns joined from the 2026-09-04 bundles by (id, m) with row-wise seed equality asserted).

## GATE 2: GREEN

All 16 cells: **2000/2000 replicates, 0 worker errors, 0 naive mismatches, 0 MR (current) mismatches, 0 MR NA, 0 field NA** (totals over 32,000 replicates: 0/0/0/0/0). The two pairing proofs held everywhere:

- naive column `identical()` to the stored replication bundles (as on 2026-09-04);
- **MR (current) est / se_ij / bias_sel / bias_fix / p̂_Ĥ `identical()` to the 2026-09-04 `mr_vs_guohe_*` bundles on every replicate** — the same gate call produced both MR rows, and the field addition changed nothing upstream of it.

The Stage 3 qmd refuses to render on any completeness / naive / cur violation (load-time guard); it rendered.

## Interruption and resume (reported plainly; no settings changed)

At ~34 min into the run, with **15 of 16 cells complete and every tally zero**, the driver process was killed externally by the session's low-memory watchdog while `t7_beta2_05` was in progress (the 14 h timeout did not trip; no in-band error occurred; host inspection immediately after showed 180 GB of 251 GB free, with desktop applications the top consumers — a transient spike meeting a conservative watchdog, not an exhaustion). Disposition under the pre-authorization: this was neither a cell failure nor a settings problem, and the driver is restartable by design (skip-if-exists; per-replicate seeding makes results worker-count-invariant, the property proven in the replication drivers). The single missing cell was re-invoked with **identical statistical settings** — same seeds, draws, field parameters — at 60 workers (within the authorized "≤ 120 workers" envelope, chosen to halve fork memory pressure). It completed in 2.3 min with all tallies zero. Cumulative wall ≈ 40 min, far inside the 14 h budget; nothing about scope or statistical settings changed at any point.

## Per-cell tally

Every row reads `reps 2000/2000, errored 0, naive_mm 0, cur_mm 0, mr_na 0, fld_na 0` — cells t35_beta2_00..05 (0.8–1.2 min each at 120 workers), t6_k02/06/10/12 (1.1–2.3 min), t7_beta2_00..05 (2.1–2.4 min; the final one 2.3 min at 60 workers). Full tally table in `mr_field_vs_guohe.qmd` §Gate 2 record.

## Outputs

`mr_field_vs_guohe_<id>.rds` × 16 (tracked): per-replicate records (identifiers/seeds; naive; G&H four-r columns from stored arithmetic; MR (IJ) columns identical to 2026-09-04; the full `field` block — λ̄, λ-sd, seven quantiles, est2, one-sided and two-sided bounds, SE-type bounds, draw-usage counts — with the two field coverage indicators; joined addendum-A diagnostics; timings), gate tallies, settings, elapsed, sessionInfo. The 2026-09-04 bundles were read-only inputs and are unmodified.
