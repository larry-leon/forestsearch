# REPORT — `oc_wrapper_verification.qmd` re-pointed at the corrected run; Limitations and Breadth sections added

**Task:** `dev/tasks/cc_task_oc_verification_repoint_2026-08-31.md` (arrived in `~/Downloads` as `cc_task_oc_verification_repoint_20260831.md`, copied as `cc_task_oc_verification_repoint_2026-08-31.md`, committed alone as `4baa0be4`).
**Document:** `quarto/simulations/actg175/continuous/oc_wrapper_verification.qmd` (+ its tracked render `oc_wrapper_verification.html`, embedded resources, 2.7 MB).
**Document change only. No `R/` change, no driver, no payload, no `.rds` of another task modified; nothing drawn or enumerated.** Writes: the `.qmd`, its render, `oc_verification_repoint_check_2026-08-31.R` (+ `.log`) and `oc_verification_repoint_2026-08-31_qmd.diff` under `dev/glm-continuous-sims/`, this report.

## 1. Provenance — GATE

`pop-os`; `/home/larryleon/Documents/GitHub/forestsearch`; `feature/glm-extension`; HEAD `38922fee`; `git status --porcelain` empty; log `38922fee 6504e0ea b2d2a2de 8ff85f9c c9cb0ca2 f5388486` — the Stage 2 report's commit **`6504e0ea`** (read from `REPORT_oc_breadth_stage2_2026-08-31.md` line 10 of the summary) is in the log; `packageVersion("forestsearch") = 0.3.1`; `/usr/lib/rstudio/resources/app/bin/quarto/bin/quarto` v1.9.38 present (not on PATH by default; put on PATH for the render). PASS. Task doc committed alone as `4baa0be4`.

## 2. Inventory (before any edit; line numbers of the committed 502-line file)

**Location:** `quarto/simulations/actg175/continuous/oc_wrapper_verification.qmd`.

**(1) Chunks reading an `.rds`:**
- `setup`, L23–26: `.grid_path <- file.path(…, "dev", "glm-continuous-sims", "oc_wrapper_grid_2026-08-29.rds")`; `G <- readRDS(.grid_path)` — **the superseded alternative object.** L28: `.pay <- lapply(G$payloads, …readRDS…)` — the two tracked MD40 payloads (n = 500 / n = 700), paths carried by the object (unchanged in the corrected one).
- `null-setup`, L343–346: `.null_path <- …"oc_wrapper_null_2026-08-29.rds"`; `N0 <- readRDS(.null_path)` — **the superseded null object.**
*GATE:* both 08-29 reads found → the premise stands. PASS.

**(2) Superseded numbers (confs report §4, "What is superseded"):** the document's own header comment (L21–22) is accurate — *every* M, functional, inverted `c1` and null rate in prose or table is **computed** by inline R or a `kable` from `G`/`N0`: `M` at L71–76, L105–106, L373–374, L444, L450 (`G$families[[…]]$M`, `N0$family$M`); the comparison tables (`tables`, `tables-700`, L216–247); the timing tables (L249–255, L277–286); the sweep figure and inversion table (L288–322, `.iv$c1`); the null table, figure, inversion table and prose (L378–456: `.nr`, `.ns`, `.nm`, `L_eff`, `p1_range`, `.gap30`, `.gap_max`); the "Reading" paragraphs (L465–483). **No superseded number is typed.** Typed numbers that remain are not superseded: `M = 16` (the analytic document's stylised family, L56, L107, L444, L455), the `sigma_D` diagnostic figures (4 %, 7 %, Spearman −0.07/−0.12, L180–184; the 08-29 diagnostic on one trial), "0.86 %" (the alternative's three-region `V_eff` spread, L366; a property of the DGM scale, unchanged), "0.13–0.28 SEs" (the analytic document's `worked-null`, L446).

**(3) Chunks that compute rather than read:** none draws, enumerates or inverts. `sweep-reading` (L324–328) computes the first `c1` below 0.95 from the stored sweep; `null-reading` (L430–434) picks the largest gate gap from the stored `gate_compare`; everything else formats. **There are no fidelity chunks in this document**: the comparison against the prediction document's `worked-*` chunks lives in the precompute scripts (the corrected script's `stopifnot(identical(…))` guards, re-passed in the confs task) and is only *described* here (L57–58, "bit-identical to the `worked-predictions` chunk at its own seed"). Nothing to keep or break on that count; stated, not worked around.

## 3. Re-point and de-type (the diff: `oc_verification_repoint_2026-08-31_qmd.diff`, 190 lines; 140 insertions, 9 deletions)

- `setup`: reads `oc_wrapper_grid_corrected_2026-08-30.rds` as `C`; asserts `C$superseded[["alt"]]` names the 08-29 grid; `G <- C$alt` with `fs_args`, `seed`, `draws` copied from the top level — the corrected object's `$alt` mirrors the old object's fields (verified field by field: `runs`, `families`, `measured`, `sweep`, `sweep_timing`, `invert`, `scale`, `payloads`; the one renamed column is `invert$table$c1 → $value`). `.dev` (the `dev/glm-continuous-sims` path) is defined once for the new sections.
- `null-setup`: `N0 <- C$null` (with `fs_args`, `draws`, `seed`), asserting `C$superseded[["null"]]` names the 08-29 null. The 08-29 files are left in place.
- `sweep-invert-table`: `.ivt$c1` → `.ivt$value`.
- One sentence added at the end of §1 ("What the wrapper computes"): the family is the corrected 13-confounder enumeration, `M` per cell as inline code (`r .fm(G$families[["500"]]$M)` / `…"700"…` / `r .fm(C$null$family$M)` → 1,696 / 1,890 / 1,696), with the one-line note that the 08-29 build used 12 variables (M = 1601 / 1784 / 1601) and that every functional moved by ≤ 0.1 and every inverted `c1` by < 0.25.
- De-typing: nothing to de-type (§2(2)); no table was typed.
- Everything else in the existing document is unchanged.

## 4. Limitations — verbatim

Section **Limitations** (§11 of the render) inserted after the new Breadth section and before "Reading the comparison" (the closing reading), containing the two paragraphs exactly as issued. Inline-code substitutions: **none were possible** — the corrected `.rds` does not carry the between-rule gaps (+2.11 / +0.61 / +1.65 live in `oc_wrapper_confs_compare_2026-08-30.rds`, not in the corrected grid object), so they are typed as given; Q's thresholds (age 34, preanti 744.5) are recorded in no stored object the document reads (the DGM records Q as `z1 = 1 & z2 = 1`; the payload `truth` carries effects and prevalence only; the family's `cuts` are the population grid, which does not contain them), so they are typed as given and the parenthetical "(its definition as recorded in the DGM object, inserted as inline code)" is not reproduced because there is nothing to insert; the 0.345 is typed as issued. No other word changed.

## 5. Breadth: the effect regime (§10 of the render) — computed from `oc_breadth_ladder_2026-08-30_{gate,rung60…rung160,forecast120}.rds` and `oc_breadth_stage2_2026-08-31_score.rds`

1. Design paragraph: `q` 60–160 against MD40's 40; `k_inter(q) = k40 + s·(q − 40)` with `k40` and `s` inline; family/`se_g`/draws identical across rungs (monotone by construction); forecast rung, `c1*` = 135.741, predicted power 0.8000 (0.0009), predicted null 0.0387 (0.0004), lock commit `c9cb0ca2` — all inline. It states that the measured cell is Stage 2's run 1 scored post hoc, reproduced by run 2 (1000/1000, 786/786), via `S2$measured_cell` (the alternative wording is coded for the run-2 case).
2. Ladder table: `q`, `k_inter`, power at `c1_05` / `c1_10` (with MC SE), inverted `c1` at 80/90/95 % — from the six rung `.rds`.
3. Forecast-versus-measured table at `c1*`: predicted (both gates) beside measured (SE) for power, null type-I, E|Ĥ| (sample and population), PPV, sensitivity, specificity (sample and population each), E[β(Ĥ)], naive bias, the between-rule gap beside +2.11 — from the Stage 2 `.rds`; the verdict column is transcribed from the Stage 2 report §5.
4. Figure: declaration rate against `c1` on 0–300 — predicted q = 120 (rung grid), measured q = 120 (run 1 post hoc), predicted MD40 (corrected sweep, 20–120), null predicted (corrected null sweep, 20–120) and null measured (payload, post hoc); guides at `c1_05`, `c1*`, 0.80; production gate only, labelled.
5. `se_g` band per rung (min / median / max of `se_g / se_direct`, count outside ±2 %) from the gate and rung `.rds`, followed by the one sentence on the direct-`V_eff` sensitivity (0.800 → 0.810, 0.775 → 0.763, both read inline) and the constructor being unchanged, decision open in the handoff.

## 6. Render and check — GATE

Render: `PATH=/usr/lib/rstudio/resources/app/bin/quarto/bin:$PATH quarto render oc_wrapper_verification.qmd` → `Output created: oc_wrapper_verification.html`, **8.1 s wall** (8.5 s on the first clean render); 21 chunks, no warning, no error. Two syntax slips on the way (an escaped quote inside a `sprintf` caption; `c1*` inside an inline-R string rendering as emphasis) were fixed before this render; both are in the committed text's final form only.

**The comparison** (`oc_verification_repoint_check_2026-08-31.R`, log beside it): the rendered HTML is stripped of scripts, styles and the `code-tools` embedded source, and every value the document shows for the alternative and null cells is recomputed from the corrected object at the document's own precision (`.f`/`.fm`) and searched as a whole number token (not inside a longer number, not as a longer decimal):

| group | values checked | found in render |
|---|---:|---:|
| alternative predictions, both gates × n = 500/700 × 8 functionals + M | 36 | 36 |
| family M and stage counts (enumerated, kept, duplicate, size) per n | 10 | 10 |
| alternative inversions (12 rows, `.f(value, 2)`) | 12 | 12 |
| null predictions, both gates (det_rate, EnH, Espec, Enpv, EbetaH, Enaive_bias), M, L_eff | 15 | 15 |
| null inversions (14 rows) | 14 | 14 |
| measured columns (8 × 2) | 16 | 16 |
| **total** | **103** | **103 (0 missing)** |

The same keys formatted from the 08-29 objects give 46 superseded strings that differ from their corrected counterparts *and* equal no corrected value elsewhere in the document; **0 of the 46 appear in the render.** Raw greps: `1601` 2 hits and `1784` 1 hit — exactly the "M = 1601 / 1784 / 1601" of the §3 historical note and nothing else; `1,601`, `1,784`, `0.9970`, `133.11`, `grid_2026-08-29`, `null_2026-08-29` 0 hits; `1,696` 13, `1,890` 5, `corrected_2026-08-30` 0 in the body (the path appears only in folded code). (`91.90` has one hit: it is the corrected n = 500 `"split"` 80 % inversion, 91.90, not the old n = 500 `"resample"` value, which is now 92.04.) *GATE:* PASS. Fidelity chunks: none exist in this document (§2(3)); the precompute's guards are the standing check and were not touched.

## 7. Commits

`4baa0be4` task doc · the `.qmd`, the render, the check script and log, the diff and this report by explicit path as `606a9352`. **No push. No install.**

## 8. Ten-line summary

1. Task doc committed alone as `4baa0be4`; §1 gate passed (Stage 2's `6504e0ea` in the log, clean tree, 0.3.1, quarto 1.9.38).
2. Inventory: two superseded reads (`setup` L23–26, `null-setup` L343–346; payloads via `G$payloads` L28); no superseded number typed anywhere; no fidelity chunk in the document.
3. Re-pointed both reads at `oc_wrapper_grid_corrected_2026-08-30.rds` (`$alt` / `$null` aliased to the chunks' `G` / `N0`, top-level `fs_args`/`seed`/`draws` copied; `invert$c1 → $value`), 08-29 files left in place.
4. Added the family sentence: M = 1,696 / 1,890 / 1,696 as inline code with the 12-variable (1601 / 1784 / 1601) note and the "≤ 0.1 / < 0.25" line.
5. Limitations section inserted verbatim; nothing inside it could be read from the corrected `.rds` (gaps and Q's thresholds are not carried there), so those figures are typed as issued.
6. Breadth section: design paragraph, ladder table, forecast-vs-measured table with transcribed verdicts, the `c1` figure (five curves, three guides), the `se_g` band table and the sensitivity sentence — all read from the ladder and Stage 2 `.rds`; measured cell = run 1 post hoc (run 2 reproduced it).
7. Render 8.1 s; 21 chunks; no warnings.
8. §6 check: 103 / 103 corrected values present; 0 / 46 distinguishable superseded values present; `1601`/`1784` only in the historical note.
9. Nothing drawn, enumerated or inverted; no `R/`, driver, payload or foreign `.rds` touched.
10. Committed as `606a9352`; this line was added in the immediately following commit.
