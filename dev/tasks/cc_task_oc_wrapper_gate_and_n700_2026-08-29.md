# CC task — verify `sigma_D`, add the closed-form gate, run n = 700, and document both

**File:** `dev/tasks/cc_task_oc_wrapper_gate_and_n700_2026-08-29.md` · **Issued:** 2026-08-29 by chat, commissioned by Larry
**Repo:** `~/Documents/GitHub/forestsearch` (pop-os) · **Branch:** `feature/glm-extension`
**Transport:** Larry places this file in `~/Downloads`; CC copies it into `dev/tasks/` and commits it first (§1).
**Follows:** `dev/glm-continuous-sims/REPORT_oc_wrapper_build_2026-08-28.md` and `REPORT_oc_wrapper_fixture_run_2026-08-28.md`. Read both before starting.

---

## ⚠ `R/` CALLOUT

**One conditional edit to one file we wrote two days ago: `R/fs_oc_predict.R`, to implement the `consistency_method = "resample"` branch that currently `stop()`s.** Category: *changes behaviour of a new function by adding an already-declared argument value*. It changes no search behaviour and no method in the package's estimation path. **It happens only if §2's gate opens.** No other `R/` file may be edited.

**Compute:** minutes, not cells. §2 fits a handful of subgroup models on **one** simulated dataset (diagnostic, §2.3). §4 evaluates `fs_oc_predict()` four times at ~2×10⁵ draws — estimate 20–40 min total, single machine, no replicates. No simulation study, no renders except the new document.

---

## 1. Provenance and first commit — GATE

```bash
hostname; pwd
git rev-parse --abbrev-ref HEAD; git rev-parse --short HEAD
git status --porcelain
git log --oneline -6
Rscript -e 'cat("forestsearch", as.character(packageVersion("forestsearch")), "\n")'
```

*GATE:* branch `feature/glm-extension`; HEAD descends from `8ecec589`; installed version `0.2.3`. Dirty tree is not a failure — list it and work around it.

```bash
ls ~/Downloads/cc_task_oc_wrapper_gate_and_n700*
cp ~/Downloads/cc_task_oc_wrapper_gate_and_n700*.md dev/tasks/cc_task_oc_wrapper_gate_and_n700_2026-08-29.md
git add dev/tasks/cc_task_oc_wrapper_gate_and_n700_2026-08-29.md
git commit -m "task doc — verify sigma_D, add the closed-form gate, run n = 700, document both"
```

*GATE:* exactly one match in `~/Downloads` (stem may have hyphens stripped; copy under the name above and report both). Run `devtools::test()` now and record the counts as the parity baseline.

---

## 2. The `sigma_D` identification — GATE

Chat's claim, to be confirmed or refuted: **in the document's own draw model, `sigma_D` equals `se_g`**, so the closed-form screen needs no new quantity. The reasoning is that the document draws `W1`, `W2` independently with covariance `2*Sg` each, so `Bhat = (W1+W2)/2` has covariance `Sg` and `D = (W1-W2)/2` has covariance `Sg` also. Do not take this on trust; it is chat's derivation, not source.

### 2.1 From source — `R/consistency_resample.R` and its call site

Quote verbatim and state each answer:

1. The half-effect convention. Is it `beta_hat ± D`? Which of `W1`, `W2` is which, and is the split 50/50?
2. `D`'s definition: the multipliers `G_i` — their distribution and variance — and whether `D = sum_i G_i * dfbeta[i, treat]` as the header states.
3. `sigma_D`'s computed form (report says `sqrt(sum(dfbeta[i, treat]^2))` at `:86, :318` — confirm) and any scaling factor applied anywhere between its computation and its use.
4. The comparison scale `c` in the rate: is it `consistency.threshold` on the effect scale, and is any orientation flip applied?
5. The exact pass condition, and whether `pconsistency.threshold` is compared to the rate with `>=` or `>`.

### 2.2 The threshold in general form

From (4)–(5), derive and state the inversion, keeping `pconsistency.threshold` symbolic rather than substituting 0.90:

```
2 * Phi((beta_hat - c) / sigma_D) - 1  >=  p
    <=>   beta_hat  >=  c + qnorm((1 + p) / 2) * sigma_D
```

Confirm `qnorm((1+0.90)/2) = 1.6449`. If the source's condition differs in any way from the left-hand side, the wrapper must implement what the source does, not this.

### 2.3 DIAGNOSTIC — the numerical identification

**Scope: one simulated dataset. This is an identity check, not a result, and nothing from it is a reported number.**

Build the MD40 DGM deterministically as `fixture_run_fs_oc_predict_2026-08-28.R` does. Draw **one** trial at `n = 500` with a stated seed. Choose ~10 subgroups spanning the prevalence range (some near the `n.min` floor, some large). For each, using the package's own code path:

| quantity | how |
|---|---|
| `sigma_D` | the value `consistency_resample`'s machinery computes |
| `se_model` | classical SE of the treatment coefficient from the subgroup `lm` |
| `se_hc0` | `sqrt(sum(dfbeta[, treat]^2))` computed directly |

Report the table and the two ratios `sigma_D/se_model`, `sigma_D/se_hc0` per subgroup.

*GATE, three outcomes:*

- **Ratio ≈ 1** (to numerical tolerance) for `se_hc0`, and close for `se_model`: the identification holds. Proceed to §3.
- **Ratio ≈ a constant** other than 1 (e.g. `sqrt(2)`, `1/2`): the identification holds up to that factor. **Report the factor, carry it explicitly and name it in the code and the document**, then proceed to §3.
- **Ratio varies systematically with prevalence or subgroup size:** the identification fails. **Stop.** Do not implement §3. Write the report with the table, commit, and end — this is then a method question for Larry.

Also state whether `se_model` and `se_hc0` themselves agree; a large divergence is worth recording even where the gate opens, because `se_g` in the wrapper is a model-based population quantity.

---

## 3. `R/fs_oc_predict.R` — the `"resample"` branch (only if §2 opened)

Replace the `stop()` with the closed-form screen. Under the identification of §2 the gate is a single threshold on the full-sample draw:

```r
pass <- (Bhat >= c1) & (Bhat - c2 >= z_p * sigma_D_g)
```

with `z_p <- stats::qnorm((1 + pconsistency) / 2)` and `sigma_D_g` the family's `se_g` (times §2's constant factor if one was found).

Requirements:

- Add `pconsistency` as an argument, defaulting to `forestsearch_args$pconsistency.threshold`, itself defaulting to the package default. Explicit values override, as `c1`/`c2` do.
- **This branch needs only `Bhat`.** Do not materialise `W1`/`W2` when `consistency_method = "resample"` — one draw matrix, not two. Report the memory and wall-clock difference.
- The `"split"` branch is untouched. The §5 fidelity harness of the build task must still pass bit-identically; re-run it and say so.
- Roxygen: document both branches, state the `sigma_D = se_g` identification with its §2 evidence, and say plainly that `"resample"` is the package's production screen while `"split"` is the document's historical one.
- Extend `tests/testthat/test-fs-oc-predict.R`: on a one-candidate family, the `"resample"` declaration probability equals a direct `pnorm()` evaluation of the threshold; and `"resample"` at a very small `pconsistency` approaches, without equalling, the `"split"` rate. Keep them fast.

---

## 4. Precompute — four evaluations to one `.rds`

Extend `dev/glm-continuous-sims/fixture_run_fs_oc_predict_2026-08-28.R` into a new script `dev/glm-continuous-sims/oc_wrapper_grid_2026-08-29.R` (leave the old one untouched). It writes **one** `.rds` under `dev/glm-continuous-sims/` holding, for each of `n ∈ {500, 700}` × `gate ∈ {split, resample}`:

- the returned `fs_oc_predict()` object with its MC SEs;
- the family's `M`, the enumeration counts at each stage, and the floor used (`n.min/n` differs between the two n);
- the settings and seed.

Plus, read from the **tracked payloads**, the measured column for each n. **Take measured values from `payload$oc` (and `payload$results` where `oc` does not carry a quantity) — computed, never typed.** The prediction document's `m500` / `meas` / `measd` blocks are typed transcriptions; do not use them as the measured source. If a quantity genuinely has no machine-readable counterpart, record `NA` with a note rather than transcribing.

Seeds stated in the script. Report wall-clock per evaluation.

---

## 5. The document

New Quarto document beside the prediction document:
`quarto/simulations/actg175/continuous/oc_wrapper_verification.qmd`, titled for what it is (a verification and methods note for `fs_oc_predict()`).

It **reads the `.rds` from §4 and the two payloads; it computes nothing heavy and types no number.** Every figure in prose is inline R. Render must be seconds.

Contents, in this order:

1. **What the wrapper computes** — the family's nine-field interface, and the statement that everything above it is family-agnostic.
2. **The family** — population enumeration under the `forestsearch()` cut specification: cuts at population quantiles from `conf.cont_jcuts` and the `cut_type` defaults, both directions, combinations to `maxk`, then the structural floors (`minp`, `rmin`, `Pg >= n.min/n`). State what is deliberately excluded (GRF pre-screen, statistics-keyed dedup) and why. Give the realized `M` and the stage counts for both n.
3. **Means and standard errors** — `beta_g = tauQc + bint * PQg` as the mixture identity, `se_g` as the prevalence-scaled anchor SE, with the MD-only caveat stated.
4. **The covariance** — `Rho = ovl / sqrt(outer(Pg, Pg))`, `Sg = Rho * outer(se_g, se_g)`, the symmetric root and why the symmetric one (continuity at rank deficiency).
5. **The two gates** — the derivation from §2, written out: half effects `beta_hat ± D`; `P(|D| <= beta_hat - c) = 2*Phi((beta_hat-c)/sigma_D) - 1`; production thresholds that probability at `p`, the document realizes one `D`. Give the inversion, the `sigma_D = se_g` identification with its §2 evidence, and state plainly that `"resample"` is production and `"split"` is the document's historical gate.
6. **Selection and the functionals** — the maxeffCons argmax, `sel_c`, and each reported quantity's formula.
7. **The tables** — one per n, columns: analytic `"resample"` · analytic `"split"` · simulation. Rows as in the fixture-run report, plus `M`. MC SEs shown where they exist on both sides.
8. **Reading the comparison** — brief and honest: what agreement does and does not establish; that `det_rate` saturates at the alternative so the two gates cannot be discriminated there; and that the measured side is a sample-realized, screened family while the analytic side is a fixed population enumeration.

No payload is written. Follow `claude/quarto_applications_conventions.md` for anything it does write.

Render it once; confirm exit 0 and report the time.

---

## 6. Close-out

- `devtools::document()`; `devtools::test()` — raw counts and **warning-count parity** against §1's baseline.
- `R CMD check` if §3 ran (`RSTUDIO_PANDOC` must be on Rscript's PATH for the vignette build). Report Status verbatim.
- If §3 ran: version → `0.2.4`, `NEWS.md` entry naming the `"resample"` branch and the identification.
- Commit on `feature/glm-extension`, **explicit paths only**, never `git add -A`. **No push.**
- `devtools::install()` if `R/` changed; confirm the version.
- Report: `dev/glm-continuous-sims/REPORT_oc_wrapper_gate_and_n700_2026-08-29.md` — §2's source answers and diagnostic table with which gate outcome you took; §3's fidelity re-run result; §4's wall-clocks and `M` for both n; the two comparison tables; the document's render time; commits. Ten-line verdict at the end.

---

## 7. Out of scope

- No `R/` file but `R/fs_oc_predict.R`, and that only if §2 opened.
- No edit to the prediction document, to any driver, or to the earlier scripts and reports.
- No simulation study; no replicate runs; §2.3's single dataset is the only data generated.
- No grid over `c1`/`c2`, no inversion, no null-DGM path, no binary/OR path, no second family constructor.
- If §2's gate closes, no approximate `sigma_D` — report and stop.
