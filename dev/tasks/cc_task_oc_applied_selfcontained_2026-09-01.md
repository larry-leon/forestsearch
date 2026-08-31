# CC task — the self-contained applied OC document: wrapper-forward, `c2 = 0.8·c1`, everything in `quarto/applications/actg175/`

**File:** `dev/tasks/cc_task_oc_applied_selfcontained_2026-09-01.md` · **Issued:** 2026-09-01 by chat, commissioned by Larry
**Repo:** `~/Documents/GitHub/forestsearch` (pop-os) · **Branch:** `feature/glm-extension`
**Transport:** `~/Downloads` → `dev/tasks/`, committed first (§1).
**Follows:** `REPORT_oc_applied_evaluation_2026-08-31.md` and its document. **Larry's three revisions, 2026-08-31/09-01:** (1) the evaluation must be **self-contained at the analysis stage** — it lives in `quarto/applications/actg175/`, computes everything itself, and reads nothing from `dev/`; (2) it must **illustrate the fundamental use of the wrappers** — the key calls visible and few, not buried in scratch scripts; (3) **`c2` grows with `c1`**: the practice recommendation is `c2 = 0.8·c1`, and a fixed `c2 = 5` against `c1 = 50` makes no sense — it is kept only as the liberal comparator column reconciling with the archived deep run.

---

## ⚠ CATEGORY

**No `R/` change.** Writes: **one new application document** `quarto/applications/actg175/analysis_actg175_continuous_oc.qmd`, its rendered HTML (per the directory's tracking convention) and its payload under `_payloads/`; **one bounded edit** to the existing `analysis_actg175_continuous_oc_evaluation.qmd` — a two-sentence note at the top pointing to the new self-contained document and stating that the archived run used fixed `c2 = 5` at `2e5` draws; nothing else in it changes and it is **not re-rendered**. Plus this task document.

**Compute — the go/no-go, stated.** The render is the compute: at the default parameters below (draws = 2e4, twelve (variant, rung) jobs concurrent via `parallel::mclapply`), per job ≈ enumerate 43 s + one draw set ≈ 400 s + ≈ 160 grid points × ≈ 6 s ≈ 22–25 min, run twelve-wide → **projected render ≈ 35–55 min, ≈ 10 GB peak**. *GATE:* time a single job first (a probe run of the same chunk code); if the twelve-job projection exceeds **2.5 h**, stop and report. State the measured render wall-clock in the report.

**Unattended.** Gates stop the task; never ask, never work around. **Provenance gates are dirt-tolerant:** list unrelated untracked/modified paths, leave them alone; stop only if dirt touches this task's own paths.

---

## 1. Provenance and first commit — GATE

Standard block plus the installed version (expect **0.3.2**). *GATE:* branch `feature/glm-extension`; commits `209a8e85` (stage 0), `eb136a35` (0.3.2) and `d597f289` (the evaluation report) in the log. Copy this document into `dev/tasks/` and commit it alone; the `~/Downloads` stem arrives with hyphens stripped, so report both names.

## 2. The document — design spec

**YAML `params`** (so Larry can re-render at other precisions without editing code): `draws` (default `2e4`), `n_workers` (default 12), `c2_ratio` (default 0.8), `seed` (default 8316951), and the rung and ladder vectors and the `Q_variants` list below as visible setup-chunk defaults.

**Visible, minimal wrapper chain — this is revision (2) and the document's point.** Each of these is its own short, echoed chunk with a sentence of prose, in this order:

1. **Data prep** exactly as the compare-all document (arms 1/3, `y_decline`, twelve confounders; N = 1083) — the template document's form.
2. **The analysis**: the fixed-family `maxeffCons` anchor `forestsearch()` call, run in-document (≈ 8 s, `parallel_args = list(plan = "sequential")`), with one sentence on why the OC evaluation attaches to this configuration (the MR-anchor rationale). Print Ĥ, n(Ĥ), `T̂_obs`, `p.consistency`. *GATE:* Ĥ = `{age <= 37} & !{cd40 <= 507}`, n = 66, `T̂_obs = 87.916667` — matching stage 0.
3. **The anchored truth, and the Q knob — Larry's revision (4), 2026-09-01.** A visible chunk defining the true subgroup as a **named list of cutpoint sets** — this is the knob, and the document says so:

   ```r
   Q_variants <- list(
     primary = list(age = 37, cd40 = list(type = "greater", value = 507)),  # Hhat itself
     wider   = <same factors, more permissive cuts, prevalence ~2x primary>,
     broad   = <same factors, more permissive still, prevalence ~3x primary>)
   ```

   The **primary** is Ĥ's exact clauses (the primary evaluation). The two supersets keep the same factors with more permissive thresholds, **chosen from the family's own population cut grid** to land nearest 2× and 3× Ĥ's prevalence — pick them, report the thresholds and realized prevalences, and *GATE:* strict nesting `Q_primary ⊂ Q_wider ⊂ Q_broad` on `df_super` (with both thresholds relaxed this holds by construction; verify anyway). One `generate_glm_dgm()` call shown for the primary at the top rung (`k_inter = q + 26.978725`, which is Q-independent since the within-Q effect is `beta_treat + k_inter` for any Q) and one `fs_dgm_scale()` print — three lines of prose on the two-point mixture stylization, plus one stating plainly: *the truth is a modeling choice; the primary plants the found rule, the variants ask what the same analysis means if the truly harmed population is broader than the sharp rule the search returned.*
4. **The family**: one `fs_oc_family_enumerate()` call with its own `print()` (M = 4508, the stage counts) — *GATE:* M = 4508 — and two sentences on Q's representability (purity 0.9969, not forced).
5. **The fundamental call**: one plain `fs_oc_predict(family = fam, n = 1083, c1 = 10, c2 = params$c2_ratio * 10, ...)` at the analyst's screening floor with the scaled consistency floor, its `print()` shown — "this one call is the wrapper; everything below is this call swept over truths and thresholds."
6. **The evaluation loop**: ~15 visible lines — `mclapply` over the twelve (variant, rung) jobs, each one `fs_oc_family_enumerate()` (43 s; the Q-referencing fields are what change) plus one `fs_oc_grid()` call whose `c1` and `c2` vectors are below; no helper files, no scratch, results held in the session. `n_workers` default **12**.

**Rungs.** Primary Q: `q ∈ {0.01, 10, 20, 40, 60, 87.916667}` — the zero-plus structural null (with the same one-paragraph orientation explanation as the archived document, condensed), the analyst's `c1` floor, two intermediates, and the observed statistic. Each superset variant: `q ∈ {20, 40, 87.916667}` — twelve rung-jobs in all. **The structural null is shared across variants**: at `q = 0.01` every candidate's mean is within `0.01·PQg` of the background for any Q — verify numerically (max |`beta_g` − background| ≤ 0.01 for each variant's family) and state it, rather than re-running the null three times.

**Thresholds:** `c1_ladder = {5, 10, 15, 20, 30, 40, 50, 60, 80, 100, 120, 87.916667}` (sorted; `T̂_obs` included for the calibration read-off); `c2` vector = `{5, params$c2_ratio × each c1_ladder value}` (deduplicated) — one `fs_oc_grid()` per rung evaluates the full crossing on one draw set (~160–170 points). From the crossing, three views:

- **The diagonal `c2 = 0.8·c1`** — the headline everywhere: the type-I curve (zero-plus rung) along the diagonal with `c1_05^diag` and `c1_10^diag` at ladder resolution; the ladder table (power at the diagonal fixed-type-I thresholds; rate at the analyst's `(10, 8)`); the declared-rule table (`E|Ĥ|`, PPV, sens, spec, `E[β(Ĥ)]`, naive bias, `mass_below`) at `(10, 8)`; the calibration column `P(T ≥ T̂_obs | q)` at the diagonal point nearest `T̂_obs`.
- **The `c2 = 5` column** — one comparator table row-set at the analyst's original gate, with one sentence reconciling to the archived deep run (numbers will differ at MC precision and ladder resolution; say so).
- **The type-I surface** — at the zero-plus rung only, a small heatmap or contour of the false-declaration rate over the `(c1, c2)` crossing, with the `c2 = 5` row and the `0.8·c1` diagonal marked: the picture of revision (3). Expected and worth one sentence: along the diagonal the consistency screen binds (for noisy candidates `0.8·c1 + z_p·se_g` overtakes `c1`), so type-I control arrives at much lower `c1` and selection shifts toward larger rules — report what the numbers actually show, not this expectation.

**Turning the true-subgroup knob — the new section.** One table across the three Q-variants (prevalence and thresholds; power at the diagonal `c1_05` for each shared rung; PPV and sensitivity of the declared rule *scored against that variant's own Q* at the analyst's `(10, 8)`; and the calibration crossings — the `q` where `P(T ≥ T̂_obs | q, Q)` crosses 0.05 and 0.5, interpolated at ladder resolution) plus one figure overlaying the three calibration curves. The reading this section exists for, reported from the numbers rather than asserted: **the observed statistic constrains a trade between breadth and severity** — a narrow Q needs a large planted harm to make 87.9 unsurprising, a broad Q a smaller one — so "what harm is consistent with what we saw" has no single answer without fixing who is harmed, and the table says what each choice implies.

**Reading and limitations**, condensed from the archived document (the declaration event vs the magnitude vs the calibration; the mixture stylization; Q off-grid; fixed-family conditioning; the ~2-subject `E|Ĥ|` machinery offset) plus one new sentence on precision: **all rates carry MC SE ≈ `sqrt(p(1−p)/draws)` ≈ 0.004 at the default draws; the archived companion holds the fixed-`c2` analysis at 2e5 draws.**

**Payload** per the applications conventions (`table` = the diagonal ladder; `extras` = type-I diagonal + surface summary, declared, calibration, anchor, purity, the `c2` policy; `est_scale = "md"`).

**Self-containment — GATE:** the qmd contains no path into `dev/` (grep the file for `dev/` — zero hits) and no `readRDS` of anything it did not itself write; every number in prose is inline R.

## 3. Render and checks — GATE

Probe one rung for the compute gate (§CATEGORY), then render with the repo's quarto convention, recording wall-clock and peak memory. *GATES:* the §2 anchor and family gates; the zero-plus diagonal type-I at the analyst's `c1 = 10` is finite and its declared-rule row shows the oriented-negative `E[β(Ĥ)]` (the honest false-positive reading); ladder power columns monotone in `q` (within 2 MC SEs); the payload written; the rendered HTML free of `NA ±` and error text in the body.

## 4. Close-out

- Commit by explicit path: the new qmd, its HTML per the directory convention, the payload, the two-sentence pointer edit to the archived qmd (its HTML **not** re-rendered — say so), and the report. **No push. No install. No `R/` change.**
- Report `REPORT_oc_applied_selfcontained_2026-09-01.md`: provenance · the §2 gates · the diagonal headline numbers beside the archived run's fixed-`c2` values (type-I at the analyst's gate, `c1_05^diag` vs the archived 86, the calibration crossings) with one paragraph on what the `c2` policy changed · render wall-clock and memory · commits · ten-line summary.

## 5. Out of scope

- No `R/` change; no re-render or re-run of the archived evaluation or its precompute; no change to any other document; no 2e5-draw run (Larry re-renders via `params` when wanted); no gbsg; no push.
