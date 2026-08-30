# CC task — re-point `oc_wrapper_verification.qmd` at the corrected run, and add the limitation and breadth sections

**File:** `dev/tasks/cc_task_oc_verification_repoint_2026-08-31.md` · **Issued:** 2026-08-31 by chat, commissioned by Larry
**Repo:** `~/Documents/GitHub/forestsearch` (pop-os) · **Branch:** `feature/glm-extension`
**Transport:** `~/Downloads` → `dev/tasks/`, committed first (§1).
**Reads:** `REPORT_oc_wrapper_confs_2026-08-30.md` §5 (the list of superseded numbers) + `oc_wrapper_grid_corrected_2026-08-30.rds`; `REPORT_residual_quantiles_2026-08-30.md` (the closure); `REPORT_oc_breadth_ladder_2026-08-30.md` + its `.rds` (the locked forecast, commit `c9cb0ca2`); `REPORT_oc_breadth_stage2_2026-08-31.md` + its `.rds` (the measured cell).

**What this is.** `oc_wrapper_verification.qmd` still reads the 08-29 `.rds`, built on the 12-confounder family (`M = 1601 / 1784 / 1601`). The corrected family is `M = 1696 / 1890 / 1696`; every functional moved by ≤ 0.1 and every inverted `c1` by < 0.25, so the document's conclusions stand but its numbers do not. This task re-points it, converts typed numbers to computed ones where they are typed, and adds two sections the record now requires: the residual's closure as a documented limitation, and the breadth result — the locked forecast beside the measured cell.

**Wording is Larry's.** The two limitation paragraphs in §4 are to be inserted **verbatim**; do not reword, tighten or soften them. Numbers inside them are the record's and are not to be recomputed here.

---

## ⚠ CATEGORY

**Document change only. No `R/` change of any kind**, no driver, no payload, no `.rds` produced by another task is modified. Writes: the `.qmd`, its rendered output (tracked under the repository's convention for `quarto/simulations/` HTML), and this task document.

**Compute:** the render only. The document reads stored `.rds` files and draws nothing new; if any chunk re-runs a draw set or an enumeration, that is a defect to fix by reading the stored object, not a cost to pay. Estimate **30–60 minutes** including the render and the checks. **Renders need** `/usr/lib/rstudio/resources/app/bin/quarto/bin` on PATH.

**Unattended.** Gates stop the task; never ask, never work around.

---

## 1. Provenance and first commit — GATE

Standard block plus the installed version (**0.3.1**). *GATE:* branch `feature/glm-extension`, clean tree, and the Stage 2 report's commit in the log (read its hash from `REPORT_oc_breadth_stage2_2026-08-31.md`'s ten-line summary). Copy this document into `dev/tasks/` and commit it alone; the `~/Downloads` stem arrives with hyphens stripped, so report both names.

## 2. Inventory — before any edit

Locate the document (expected under `quarto/simulations/actg175/continuous/`; say where it is). Then list, with line numbers:

1. every chunk that reads an `.rds`, with the path it reads;
2. every number in the prose or in a typed table that `REPORT_oc_wrapper_confs_2026-08-30.md` §5 marks as superseded — `M`, functionals, inverted `c1`, the null rates — and whether each is typed or computed;
3. every chunk that computes anything rather than reading it, and what it computes (expected: the fidelity checks against the prediction document's chunks, which stay).

*GATE:* the inventory must find the 08-29 `.rds` read(s). If the document already reads the corrected file, stop and report — the premise has moved.

## 3. Re-point and de-type

- Change every 08-29 `.rds` read to `oc_wrapper_grid_corrected_2026-08-30.rds`. Leave the 08-29 file where it is; do not delete or rename it.
- Where a superseded number is **typed** in prose, replace it with inline code that reads the value from the loaded object (`` `r ...` ``), formatted to the document's existing precision. Where a table is typed, rebuild it from the object. **Computed, never typed** is the rule the reports follow and the document should too.
- The fidelity chunks that compare the wrapper against the prediction document's own `worked-*` chunks are unchanged: they are the standing guard and do not depend on which family the grid used.
- Add one sentence near the top stating that the family is the corrected 13-confounder enumeration and giving `M` per cell as inline code, with a one-line note that the 08-29 document used 12 and that every functional moved by ≤ 0.1.

## 4. The limitation section — verbatim

Add a section titled **Limitations** after the results and before any closing summary, containing exactly these two paragraphs (numbers are the record's; the `E|Ĥ|` figures are read from the corrected `.rds` as inline code where the object carries them, otherwise typed as given):

> **The size of the selected rule.** The prediction conditions on a fixed candidate family enumerated at population quantiles on the DGM's super-population. Within every covariate signature the realized subgroup size matches its population value, but the search selects rules that are on average larger than the analytic selection weights predict — by +2.11 subjects at n = 500 under the alternative, +0.61 at n = 700, and +1.65 at the n = 500 null (corrected 13-confounder family, 2 × 10⁵ draws, seed 20260825). `E|Ĥ|`, PPV, sensitivity and the naive bias inherit this gap; the declaration rate, specificity, `E[β(Ĥ)]` and the null false-declaration rate do not. Five mechanisms were tested and eliminated: sample realization within signature; the confounder list; the prevalence scaling of `se_g`, measured within [0.992, 1.015] across 5,282 family members; statistics-keyed near-duplicate removal, at most 5% of the gap; and cut placement at replicate rather than population quantiles, under which a replicate-averaged family moves the prediction away from measurement, widening the gap to +3.03 / +1.99 / +2.55. The realized family also differs from the enumerated one in structure: at this DGM `wtkg` retains a fourth default cut in most n = 500 samples that the population enumeration collapses, `karnof` falls below four distinct values and is handled as categorical in about 40% of samples, 71% of the selection mass sits on affected candidates, and only about a quarter of replicates realize the population label set exactly — variation the fixed-family framework cannot represent. The gap is recorded as a limitation of the fixed-family prediction and is not explained.

> **The true rule need not be a family member.** The enumeration cuts continuous covariates at a specified population quantile grid, so a true subgroup defined by thresholds off that grid is represented only by its nearest members. At this DGM, Q — defined by thresholds at age 34 and preanti 744.5 (its definition as recorded in the DGM object, inserted as inline code) — is not a member of the enumerated family; the closest candidates have purity `PQg ≥ 0.95`, and at the breadth forecast rung 0.345 of the selection mass sits on them. Sensitivity and PPV of the selected rule are therefore bounded by the grid before any sampling enters, and "power" throughout this document is the probability of declaring some rule, whose composition the selection distribution reports.

## 5. The breadth section — read from the stored objects

Add a section titled **Breadth: the effect regime**, entirely computed from `oc_breadth_ladder_2026-08-30*.rds` and the Stage 2 `.rds`, containing:

1. One paragraph stating the design in the record's terms: same population and Q, `k_inter` moved so that the oriented Q effect runs 60–160 against MD40's 40; only `beta_g` moves across rungs (state that the family, `se_g` and the draws are identical, so the ladder is monotone by construction); the forecast rung, `c1*`, and the two predicted rates, each as inline code; the lock commit `c9cb0ca2`.
2. **The ladder table:** `q`, `k_inter`, power at `c1_05` and at `c1_10`, and the inverted `c1` for 80 / 90 / 95%, from the ladder `.rds`.
3. **The forecast-versus-measured table** at `c1*`, from the Stage 2 `.rds`: predicted (both gates) beside measured with SEs, for power, null type-I, `E|Ĥ|`, PPV, sensitivity, specificity, `E[β(Ĥ)]`, naive bias, and the between-rule size gap at this harm beside +2.11 — with the Stage 2 report's verdict word (within noise / population-versus-sample / gap) in a final column, transcribed from that report.
4. **The figure:** declaration rate against `c1` — predicted curves for MD40 and the forecast rung from the grid tables, the measured curve from Stage 2 run 1, and the null's predicted and measured curves; vertical guides at `c1_05` and `c1*`; a horizontal guide at 0.80. Both gates need not be drawn; the production gate is enough, labelled.
5. **The `se_g` band** per rung, one small table from the ladder `.rds`, followed by one sentence stating that the direct-`V_eff` sensitivity at the forecast rung moved power to 0.810 and PPV to 0.763 and that the constructor is unchanged — the decision is recorded as open in the handoff.

If Stage 2 named run 2 (the direct run) as the measured cell rather than run 1's post-hoc scoring, use that and say so in the paragraph.

## 6. Render and check — GATE

Render. Then, programmatically, compare every number the rendered document shows for the alternative and null cells against the corrected `.rds` (read the rendered markdown or HTML, extract the values at the document's precision, compare) — *GATE:* every superseded number now agrees with the corrected object; none of the 08-29 values remains anywhere in the rendered output (grep the render for `1601`, `1784` and each superseded functional at its old precision, expecting zero hits outside the one-line historical note of §3). Confirm the fidelity chunks still pass. Report the render's wall-clock.

Commit the `.qmd` and the rendered output by explicit path. **No push. No `devtools::install()`.**

Report `REPORT_oc_verification_repoint_2026-08-31.md`: the §2 inventory · the diff of the `.qmd` · §6's comparison table and the grep result · commits · ten-line summary.

## 7. Out of scope

- No `R/` change, no driver, no payload, no other document; the prediction document and its `worked-*` chunks are not touched.
- No new draws, enumerations, inversions or replicates — everything is read.
- No rewording of §4's paragraphs. No mechanism for the residual. No change to `se_g`.
- No push, no install.
