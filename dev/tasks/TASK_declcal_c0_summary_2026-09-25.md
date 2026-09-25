# TASK — Build the S1.8 payload: every number in the manuscript's declaration-calibration tables, from the committed campaign payloads

**Repo:** `forestsearch`, branch `feature/glm-extension`.
**Kind:** read-only on everything that exists. It reads committed payloads, computes summaries,
and writes one script, one CSV, one provenance file and one report. **No simulation is run, no
replicate is generated, no payload is modified, and nothing under `R/` is edited.** R is used to
read `.rds` payloads and summarise them; that is the whole of the compute.

**Why.** Supplementary Section S1.8 of the `fs-glms-interpretable` manuscript carries four tables
of campaign results — about 130 literal values — typed from the handoff and the campaign
report, under no check. This task produces the same numbers from the payloads themselves, as a
CSV the manuscript will draw its tables from at render, with a render-time check on the S1.7
pattern. **A disagreement between the payload and a typed value is a finding, not a failure:
the payload is right, and the typed table is what gets corrected afterwards.**

**First action:** copy this file into `dev/tasks/` and commit it.

---

## Step 0 — Gates (stop on failure; do not ask)

- **G0.1** No **tracked** file is modified before any step. Untracked files already in the tree
  are not a failure; list them in the report and leave them alone.
- **G0.2** The repository is `forestsearch` and the branch is `feature/glm-extension`. **The
  checkout is current:** `git fetch origin` and confirm the local branch is not behind
  `origin/feature/glm-extension`; if it is behind, **stop and report** — the pull is Larry's,
  through GitHub Desktop, not this task's. (This machine is Pop!_OS; the `s5rerun` closeout was
  committed from the Mac Studio, and the closeout regen in Step 5 must pin to a HEAD that
  includes it.) No campaign is running: `s5rerun` (§2.14 of
  `quarto/simulations/gbsg_020/status_curated.md`) is complete.
- **G0.3** These payloads exist under `quarto/simulations/gbsg_020/results/`:
  `declcal_{bnull,inull,power}_*_res_1_2000.rds` (13 cells) and
  `declcalc0_{inull,power}_*_res_1_2000.rds` (10 cells). Checksum each before and after; every
  checksum must be unchanged at the end (P2).

## Step 1 — The script

`dev/analysis/declcal_c0_summary/declcal_c0_summary.R`, runnable from the repository root with
no arguments. It reads the payloads' `results` data frames and computes the quantities below,
each with the column and the cell set it came from. **It types no number.** Every value in the
output is computed from a payload column.

The column semantics are those of `scripts_dinamr/declcalc0_run.R` at `a46bf9b7` and of
`dev/reports/REPORT_declcal_c0_campaign_reads_2026-09-24.md`, Q6:

- `declared_conv` — the conventional screen **as executed** (rounded rate, post-reduction family);
  matches the search's own indicator 2,000 / 2,000 in every cell.
- `max_T_post` — the maximum over the post-reduction family the executed screen evaluated;
  `NA` counts as not declared.
- `max_T_pre` — the maximum over the pre-reduction family the calibrated rule reads.
- `declared_cal05_<c0>`, `declared_cal10_<c0>` — the calibrated rule at α = 0.05 / 0.10 and
  protected level `<c0>` ∈ {c070, c075, c080, c085}; `kappa_hat_05_<c0>`, `kappa_hat_10_<c0>`;
  `n_admitted_cal05_<c0>`.
- In the `declcal` (claim-threshold) payloads: the corresponding unshifted `declared_cal05`,
  `declared_cal10`, `kappa_hat_05`, `kappa_hat_10`, and `n_admitted_cal05`. Read their exact names
  from the payload; do not assume.

**Cells.** Uniform-benefit B1–B6 (design 0.657 at $n$ 500 / 1000 / 1500 = B1–B3; design 0.721 =
B4–B6); planted-harm C1–C4 (HR 1.5 and 2.0 × $n$ 1000 / 1500). Read each cell's design and $n$
from its `meta`, not from its label, and record both.

## Step 2 — What to compute

Every rate is a mean of a 0/1 column over the cell's 2,000 replicates, with the count of ones
recorded beside it so Wilson intervals can be formed downstream.

**(a) Conventional screen, cell by cell (manuscript Table S3).** `declared_conv` in each of B1–B6.

**(b) Calibrated screen at α = 0.10 and at α = 0.05, by protected level (Tables S4 and S5).**
For each α and each `<c0>`: the rate of `declared_cal<α>_<c0>` in every B cell, and its maximum
over B1–B6 with the cell(s) attaining it — **report ties**; the rate in each C cell (power).
Add three comparator rows: the conventional screen (`declared_conv`, worst B cell and each C
cell); a fixed cutoff at $k$ = 2.0 on the executed family (`max_T_post ≥ 2.0`, `NA` as not
declared — the record's "tuned fixed p\* = 0.9545"); and the claim threshold, from the `declcal`
payloads' unshifted `declared_cal<α>` on the same B and C cells.

**(c) The cutoff by sample size (Table S6).** For each α and `<c0>`, and for the claim
threshold from `declcal`: a summary of `kappa_hat` at each $n$. The manuscript prints "medians
across cells at each sample size". Compute **both** candidate definitions — the median of the
pooled per-replicate `kappa_hat` over every cell at that $n$, and the median of the per-cell
medians — and report which reproduces the typed values; record the winner as the definition.
Also the implied level $2\Phi(\hat\kappa) - 1$ of each reported cutoff.

**(d) The re-selection footprint.** Over the C cells at α = 0.05: among replicates with
`declared_cal05_<c0>` = 1, the share with `n_admitted_cal05_<c0>` > 1, per `<c0>` and pooled —
the "~97%" of `status_curated.md` §2.11, computed at the α its column exists for.

**(e) The plug-in and free-check figures are not in scope.** They are prose, not a table, and
they stay as typed.

## Step 3 — The comparison against the typed tables (report every mismatch)

The typed values are in `~/Downloads/S18_typed_tables_2026-09-25.csv`, one row per printed
value with its table, row and column labels. For every value, compare the computed one at the
printed precision (half a unit in the last printed digit). **Write every mismatch to the report,
with both values and the column and cell set the computed one came from. Do not stop on a
mismatch; do not adjust the computed value.** A mismatch of definition (item (c)) is reported
as such, with both candidates.

## Step 4 — Outputs

- `quarto/simulations/gbsg_020/results/declcalc0_summary_tables.csv` — long format: one row per
  computed value, with columns `table`, `row`, `col`, `alpha`, `c0`, `cell`, `design_hr`, `n`,
  `value`, `count`, `denominator`, `source_column`, `source_payload`. Full precision.
- `quarto/simulations/gbsg_020/results/declcalc0_summary_tables.provenance.md` — the definition
  of every quantity in words, the cell sets, and the payload checksums read.
- `dev/reports/REPORT_declcal_c0_summary_2026-09-25.md` — **a new file**: what was computed;
  the comparison of Step 3 in full, mismatches first; the item (c) definition finding; the item
  (d) share; and anything the payloads could not supply, recorded as
  `not established from source`. No recommendations and no follow-up tasks.

## Step 5 — Closeout of the directory's record

`quarto/simulations/gbsg_020/status_curated.md` is hand-maintained: add **one line** under
§2.11 naming this task, the CSV and the report. Then run `scripts_dinamr/current_status_regen.R`
as the directory's closeout step requires, so `current_status.md` regenerates with the pin equal
to HEAD at commit. Edit nothing else in either file.

## Step 6 — Bundle

Copy the CSV, the provenance file and the report to `~/Downloads/declcal_c0_summary_2026-09-25/`.

## Post-conditions (stop on failure)

- **P1** `git status` shows exactly: this task document, the script, the CSV, the provenance file
  and the report **added**; `status_curated.md` and `current_status.md` **modified**; nothing
  else. Nothing under `R/` differs.
- **P2** Every payload checksum from G0.3 is unchanged.
- **P3** The CSV has a row for every value in the typed tables, and the report lists every
  mismatch found (possibly none).
- **P4** The three files exist in `~/Downloads/declcal_c0_summary_2026-09-25/`.
- **P5** Commit. **Do not push.**
