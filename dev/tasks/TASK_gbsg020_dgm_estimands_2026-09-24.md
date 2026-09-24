# TASK — The `gbsg_020` DGM's estimands: the ITT effect and the prevalence pairing (read-only)

**Repo:** `forestsearch`, branch `feature/glm-extension`.
**Kind:** read-only. No simulation is re-run and no rate or estimate is computed. R may be
used only to read a committed object's structure and stored constants — `names()`, classes,
and the value of a stored truth field — in a session that writes nothing.
**Writes exactly two files:** this task document into `dev/tasks/`, and one report into
`dev/reports/`.
**Why:** the `fs-glms-interpretable` manuscript describes the `gbsg_020` design in two places —
Section 5 with its supplement, and the new Supplementary Section S1.8 on the declaration
calibration. Both currently name only the complement hazard ratios, 0.657 and 0.721. The
trial-wide (ITT) estimand those complements sit inside is not recorded anywhere in the
manuscript, and without it the design cannot be described coherently. It should be in the DGM.

**First action:** copy this file into `dev/tasks/` and commit it.

---

## Step 0 — Gates (stop on failure; do not ask)

- **G0.1** No **tracked** file is modified before any step. Untracked files in the tree are not
  a failure; list them in the report and leave them alone.
- **G0.2** The repository is `forestsearch` and the branch is `feature/glm-extension`. If not,
  stop and report.
- **G0.3** The `gbsg_020` campaign directory and at least one committed results payload exist.
  If not, stop and report.

## Step 1 — Read

Read-only: whichever script or record defines the `gbsg_020` data-generating processes and
their truth quantities, and the structure and stored constants of one committed bundle per
cell type. Name every path you read in the report.

## Step 2 — Answer four questions

Each answer carries **a file path and a line number**, or the exact words
**`not established from source`**. Do not infer, do not compute an estimate, and do not
reconcile a disagreement — record it.

**Q1 — the ITT estimand.** For each cell type of the campaign — the **uniform-benefit** cells,
the **attenuated-benefit** cells (planted region at a hazard ratio of one) and the **harm**
cells (planted region at 1.50 and 1.75) — give the **trial-wide marginal Cox hazard ratio** the
DGM targets or records, at each prevalence. Name the field or the script line that carries it.
State whether it is a **design constant** or a **per-replicate realized value**; if the latter,
give what the record itself quotes as its summary and say so.

**Q2 — the prevalence pairing.** Which prevalence, **12.4%** or **31%**, gives
`truth$marg_Hc` = **0.657**, and which gives **0.721**? State also for which cell types
`marg_Hc` is defined — the harm cells, the attenuated-benefit cells, or both — and whether the
two values are the same object in each.

**Q3 — what the DGM holds fixed.** Across the two prevalences, which quantity is held and
which is derived: the trial-wide effect, the complement's effect, or the planted region's? The
handoff records that the harm cells' complement "runs at `k_treat = 1`". Establish from source
whether the complement's effect is therefore the base treatment effect left unmodified — so
that `marg_Hc` differs between prevalences only because the complement is a different
subpopulation — or whether the trial-wide effect is pinned and the complement solved for.

**Q4 — the truth object.** List every field of the stored `truth` object with a one-line gloss
of what it holds and on which scale (marginal Cox, patient-level, or other). This is so the
manuscript can describe the design on one consistent scale.

## Step 3 — Write the report

`dev/reports/REPORT_gbsg020_dgm_estimands_2026-09-24.md` — **a new file; it does not exist.**
One section per question, each with its answer, its source path and line, and any disagreement
found between the source and the manuscript's or the handoff's wording. No recommendations and
no follow-up tasks: a question the record does not settle is recorded as
`not established from source` and left there.

## Step 4 — Bundle

Copy the report to `~/Downloads/gbsg020_dgm_estimands_2026-09-24/`.

## Step 5 — Post-conditions (stop on failure)

- **P1** `git status` shows exactly two **added** files — this task document in `dev/tasks/`
  and the report in `dev/reports/` — and **no modified tracked file**. If any tracked file
  differs, stop and report without committing.
- **P2** No simulation was re-run, no estimate or rate was computed, and no payload was written
  or modified. Any payload read leaves its checksum unchanged.
- **P3** Every one of the four answers carries either a path-and-line citation or the exact
  words `not established from source`.
- **P4** The report exists in `~/Downloads/gbsg020_dgm_estimands_2026-09-24/`.
- **P5** Commit the two files. **Do not push.**
