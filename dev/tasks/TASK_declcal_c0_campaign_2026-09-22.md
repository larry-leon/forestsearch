# TASK — declcal re-run with the clinically specified null level c0

**Repo:** `larry-leon/forestsearch`
**Branch:** `feature/glm-extension`
**Machine:** Pop!_OS, 64 workers (as the `declcal` campaign)
**Prerequisite:** `TASK_declcal_c0_rchange_2026-09-22.md` committed and green, and the package installed with
`devtools::install()` — workers see only the installed package.
**Kind:** simulation campaign. **This is compute.** Nothing runs until Larry's go; the go carries the ceiling.
**R/ call-out:** no change under `R/`. Any need for one is a STOP.

---

## 0. First action

1. Copy this file verbatim to `dev/tasks/TASK_declcal_c0_campaign_2026-09-22.md`; `git add` that path; commit
   `docs(tasks): add c0 campaign task document (2026-09-22)`.
2. Assert against the **installed** package: `declaration_c0` is a formal of `fs_mr_inference` with default
   `NULL`, and `fs_declaration_calibration` accepts `c0`. If either fails, run `devtools::install()` from HEAD
   and re-assert; STOP only if it still fails.
3. Record HEAD and the installed build in the report.

No `git fetch` / `pull` / `push`. Explicit paths on every `git add`.

---

## 1. What this campaign answers

For each clinically specified benefit level `c0` in **{0.70, 0.75, 0.80, 0.85}**: the false-declaration rate of
the calibrated screen under uniform benefit, and its power against a planted harm region, when the claim criterion
stays at `c1 = c2 = 1.0`. All against the same replicates, seeds and draws as `declcal`, so every difference from
the committed results is attributable to the null level alone.

"Declaration" is the family-wise event only, as before.

---

## 2. Cells — the `declcal` B and C blocks, unchanged

Block A (complete null) is **not** re-run. Larry has set it aside.

| block | cells | DGM | c1 = c2 | n |
|---|---|---|---|---|
| B | B1–B3 | uniform benefit HR 0.657 | 1.0 | 500 / 1000 / 1500 |
| B | B4–B6 | uniform benefit HR 0.721 | 1.0 | 500 / 1000 / 1500 |
| C | C1–C2 | planted harm HR 1.5, prevalence 0.1242 | 1.0 | 1000 / 1500 |
| C | C3–C4 | planted harm HR 2.0, prevalence 0.1242 | 1.0 | 1000 / 1500 |

2,000 replicates per cell. **Same seeds** (`8316951 + rep` for data, `seedit`, multipliers), same `B_cal = 500`,
same centred-Poisson law, same floors, same `pconsistency.digits` handling, same one-search-per-replicate design.
The Block C pre-flight firewall result stands (none excluded); do not re-run it.

---

## 3. The only change to the per-replicate engine

The field capture passes `declaration_c0 = c(0.70, 0.75, 0.80, 0.85)` (natural HR scale, as `hr.consistency`).
The capture now returns `Mstar_c0` (`B x 4`) alongside the unshifted `Mstar`. From these, per replicate and per
`c0`:

- `kappa_hat_05_c0`, `kappa_hat_10_c0` — `type = 1` quantiles of `Mstar_c0[, c0]`;
- `pstar_implied_05_c0`;
- `declared_cal05_c0 = (max_T_pre >= kappa_hat_05_c0)`, `declared_cal10_c0` likewise;
- `n_admitted_cal05_c0` — count over the pre-reduction family under the calibrated rule at that `c0`;
- `fw_1645_c0`, `fw_1621_c0` — `mean(Mstar_c0 > cutoff)` at the nominal and effective p\* = 0.90 cutoffs;
- `Mstar_c0_q90`, `_q95`, `_q99`.

The §5 schema of the `declcal` task is kept in full and these columns are appended, suffixed by the `c0` value
(e.g. `kappa_hat_05_c070`).

---

## 4. Identity gate — the campaign's central check

Because seeds and search are unchanged, the re-run must reproduce the committed `declcal` payloads exactly on the
quantities that do not depend on `c0`. Assert, replicate by replicate, against
`results/declcal_{inull,power}_<cell>_res_1_2000.rds` at `81752681`:

- `max_T_pre`, `max_T_post`, `G_pre`, `G_post`, `declared_conv`, `declared_conv_exact` — **identical**;
- the unshifted `Mstar` quantiles `Mstar_q90/q95/q99` and `kappa_hat_05`, `kappa_hat_10` — **identical**;
- `declared_cal05`, `declared_cal10` — identical.

Any disagreement is a **STOP** with the disagreeing replicate indices recorded. This gate is what makes the new
columns interpretable: it proves the search and the unshifted field are the same computation, so `Mstar_c0` is
the only new information.

Also assert per replicate: `Mstar_c0[, c0]` is non-increasing as `c0` decreases (larger shift, smaller maximum),
and `kappa_hat_05_c0 <= kappa_hat_05` for every `c0` in the grid.

---

## 5. Stages and the go

**Stage 1 — smoke, then continue without stopping.** Cell B5, 5 replicates. Run the identity gate of §4 on those
5 against the committed payload. Green → straight into Stage 2. Red → STOP and report.

**Stage 2 — the campaign.** 10 cells × 2,000. Blocks in order B then C, committing each block's payloads as it
completes, checkpointing per chunk as before.

**Measured pricing (from `declcal`, pop-os, 64 workers, B_cal 500):** Block B cells summed 2,213 s, Block C
1,743 s, driver overhead ~217 s. The added shift-and-max work is negligible. **Projection ≈ 4,200 s ≈ 1.2 h.**
The kickoff sets the ceiling; **3 h is recommended** (2.5× headroom).

Abort discipline as `declcal`: per-replicate cap 10× the `declcal` median (10 × 11.6 s = 116 s); per-cell gate at
1% aborts/errors; campaign hard cap 1.5 × projection = 6,300 s; do not compete with another live campaign.

---

## 6. Report — `dev/reports/REPORT_declcal_c0_campaign_2026-09-22.md`

### 6.1 Primary table — cells × screen, proportions with Wilson 95% intervals

| cell | DGM | n | p\* = 0.90 as executed | cal α 0.05, c0 0.70 | c0 0.75 | c0 0.80 | c0 0.85 | c0 = c2 (committed) |
|---|---|---|---|---|---|---|---|---|

Same table at α 0.10. B rows then C rows; no Block A rows. Short plain-language reading after each table — a
few bullets, one claim each.

### 6.2 The calibration's own quantities per `c0`

| cell | c0 | kappa_hat_05 (median, IQR) | implied p\* (median, IQR) | mean fw_1645 | mean fw_1621 | n_admitted_cal05 (median; share = 1) |
|---|---|---|---|---|---|---|

### 6.3 Readings the report must state

- For each `c0`: the worst uniform-benefit false-declaration rate across B1–B6 and which cell it is; the HR 1.5
  and HR 2.0 power at n 1000 / 1500. One table with `c0` as rows so the trade reads top to bottom.
- Beside it, the same three quantities for the committed fixed-p\* choice (k = 2.0, p\* 0.9545) and for the
  committed `c0 = c2` calibration, from `declcal_fixedk_practical.txt` and the committed payloads.
- Whether `kappa_hat_c0` is configuration-invariant as `kappa_hat` was: its median across B and C cells at each n.
- The share of declaring replicates with `n_admitted_cal05_c0 == 1`, per `c0` — this decides whether the
  re-selection question (which subgroup the calibrated screen would select) is live at the lower cutoffs.
- `fw_1645_c0` beside the executed p\* = 0.90 rate in the B cells whose uniform benefit is nearest the `c0`
  (B4–B6 for `c0 = 0.70`): whether the diagnostic now tracks a realized rate.
- Every gate with its measured value. OPEN ITEMS for documentation gaps.

### 6.4 Closeout

Regenerate `quarto/simulations/gbsg_020/current_status.md` (§2.11 gains the `c0` re-run, campaign tag
`declcalc0`), with the pin-equals-HEAD post-condition. Track every payload and summary file.

---

## 7. Build posture

Transplant: copy the committed `declcal_run.R`, `declcal.sh`, `declcal_findings.R` and change named lines only —
the `declaration_c0` argument, the appended columns, the identity gate, the campaign tag `declcalc0`, and the cells
files without Block A. Record the lines changed. Payloads to
`results/declcalc0_{inull,power}_<cell>_res_1_2000.rds`. Installed package only. No R CMD check, no full suite.

---

## 8. Commit plan

1. task doc; 2. scripts (before any run); 3. smoke record; 4. payloads per block; 5. report; 6. closeout.
Explicit paths; no push.

---

## 9. Out of scope

- No `R/` change. No Block A. No GRF / DINA. No MR post-selection correction, bounds or complement.
- No `sg_focus` re-run on the calibrated set — `n_admitted_cal05_c0` is recorded so the question can be priced.
- No change to c1, c2, p\*, `B_cal`, seeds, floors or the DGMs.
