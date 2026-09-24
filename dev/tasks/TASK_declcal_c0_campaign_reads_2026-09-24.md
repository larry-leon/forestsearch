# TASK — Read the c₀ campaign record and answer five questions (read-only)

**Repo:** `forestsearch`, branch `feature/glm-extension`.
**Kind:** read-only. This task runs no R, re-runs no simulation, and re-verifies nothing that
is already committed. It writes exactly two files: this task document into `dev/tasks/`, and
one report into `dev/reports/`.
**Why:** the `fs-glms-interpretable` manuscript's Supplementary Section S1.8 is being drafted
from `HANDOFF_declaration_calibration_c0_2026-09-22.md`. Five statements in that handoff are
ambiguous or given only as a range, and the drafting chat will not guess at them. The answers
go in the paper.

**First action:** copy this file into `dev/tasks/` and commit it.

---

## Step 0 — Gates (stop on failure; do not ask)

- **G0.1** The working tree is clean before any step. If not, stop and report.
- **G0.2** The repository is `forestsearch` and the branch is `feature/glm-extension`. If not,
  stop and report.
- **G0.3** These files exist at HEAD:
  - `dev/reports/REPORT_declcal_c0_campaign_2026-09-22.md`
  - `quarto/simulations/gbsg_020/scripts_dinamr/logs/declcalc0_findings.txt`
  - `quarto/simulations/gbsg_020/current_status.md`

  If any is absent at HEAD, read it at `10b56b80` instead and say so in the report. If it is
  absent there too, stop and report.

## Step 1 — Read

Read-only. The three files of G0.3, and whichever script or record defines the `gbsg_020`
uniform-benefit data-generating processes — locate it yourself and name its path in the
report.

## Step 2 — Answer five questions

Each answer carries **a file path and a line number**, or the exact words
**`not established from source`**. Do not infer, do not compute anything the record does not
state, and do not reconcile a disagreement — record it.

**Q1 — the scale of the uniform-benefit designs.** The handoff lists the uniform-benefit cells
as "uniform HR 0.657 / 0.721", and separately says "their patient-level HR is 0.657 at the
0.721 marginal point". Those readings are not the same. Establish from source: are 0.657 and
0.721 **two designs**, or **one design stated on two scales**? For each uniform-benefit design
the campaign ran, give its target hazard ratio **on the marginal Cox scale** and, if the record
states it, on the patient-level scale. Name the file that defines it.

**Q2 — the conventional screen's uniform-benefit rates, cell by cell.** The handoff gives only
a range, 0.7% to 23%. Give the declaration rate of the `p* = 0.90` screen **as executed**, for
**every** uniform-benefit cell (design × sample size), with the replicate count per cell. A
table is fine; copy it as the record has it.

**Q3 — which cell is worst.** For each row of the handoff's §4.1 (α = 0.10) and §4.2
(α = 0.05) tables, say which uniform-benefit cell attains the worst false-declaration rate,
**described by its design and its sample size**, not by an internal cell label. State whether
it is the same cell in every calibrated row.

**Q4 — the cutoff the campaign's Eq. 8 used.** Were the campaign's `FŴ_α(c0)` values computed
against the exact cutoff `z_{(1+p*)/2} = 1.645`, or against the executed cutoff `1.621`
(`Pcons ≥ 0.895`)? Quote the line of code or of the report that shows it.

**Q5 — the free check's margins.** The handoff's §4.4 says Eq. 8 at c₀ overstates the realized
rate "smallest +0.04 to +0.09 at c₀ 0.70 in the HR 0.721 cells". Confirm that these are
differences between `FŴ_α(c0)` and the realized declaration rate **in the same cell**, and give
the c₀ values and cells they belong to.

## Step 3 — Write the report

`dev/reports/REPORT_declcal_c0_campaign_reads_2026-09-24.md`. Structure: one section per
question, each with its answer, its source path and line, and any disagreement found between
the record and the handoff's wording. No recommendations, no follow-up tasks. A question that
the record does not settle is recorded as `not established from source` and left there.

## Step 4 — Bundle

Copy to `~/Downloads/declcal_c0_reads_2026-09-24/`:

- the report written in Step 3;
- `quarto/simulations/gbsg_020/scripts_dinamr/logs/declcalc0_findings.txt`.

## Step 5 — Post-conditions (stop on failure)

- **P1** `git status` shows exactly two added files: this task document in `dev/tasks/` and the
  report in `dev/reports/`. **No other tracked file is modified.** If any other file differs,
  stop and report without committing.
- **P2** No `.R` script was executed and no `.rds` payload was written or rewritten.
- **P3** Every one of the five answers carries either a path-and-line citation or the exact
  words `not established from source`.
- **P4** Both files in Step 4 exist in `~/Downloads/declcal_c0_reads_2026-09-24/`.
- **P5** Commit the two files. **Do not push.**
