# TASK — Read the c₀ campaign record and answer six questions (read-only)

**v2, 2026-09-24. Supersedes the v1 of the same date — do not run v1.** The only change is
Step 2's new Q6 and the wording of post-conditions P2 and P3 that follow from it.

**Repo:** `forestsearch`, branch `feature/glm-extension`.
**Kind:** read-only. This task runs no R, re-runs no simulation, and re-verifies nothing that
is already committed. It writes exactly two files: this task document into `dev/tasks/`, and
one report into `dev/reports/`.
**Why:** the `fs-glms-interpretable` manuscript's Supplementary Section S1.8 is being drafted
from `HANDOFF_declaration_calibration_c0_2026-09-22.md`. Five statements in that handoff are
ambiguous or given only as a range, and a sixth question decides whether a further number can
be had without a simulation run, and the drafting chat will not guess at them. The answers
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

## Step 2 — Answer six questions

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

**Q6 — what a screening level other than 0.90 is recoverable from, without re-running.** The
campaign's committed payloads
(`quarto/simulations/gbsg_020/results/declcalc0_{inull,power}_*_res_1_2000.rds`) are recorded
as carrying `max_T_pre`, `Mstar_c0`, `Mstar_q90/95/99`, per-replicate `κ̂` and implied `p*`.
Establish, **from the campaign report and from the payload's own structure**, whether the
**executed** declaration rate — the quantity the handoff's `p* = 0.90 as executed` row reports
as 0.2325 — can be re-evaluated at another screening level from what is already committed, and
if so which column gives it. `max_T_pre` is understood to be the maximum over the
**pre-reduction** family; if that is what it is, say so plainly, and say whether any committed
column carries the maximum over the family the executed search actually admits from. **Do not
compute any rate.** Report what is recoverable and from which column; that is all.

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
- **P2** No simulation was re-run, no rate was computed, and no `.rds` payload was written or
  modified. R may be used only to read a payload's structure for Q6 — `names()` and the like —
  in a session that writes nothing.
- **P3** Every one of the six answers carries either a path-and-line citation or the exact
  words `not established from source`.
- **P4** Both files in Step 4 exist in `~/Downloads/declcal_c0_reads_2026-09-24/`.
- **P5** Commit the two files. **Do not push.**
