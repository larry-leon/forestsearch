# CC task — the self-contained applied OC document: wrapper-forward, `c2 = 0.8·c1`, everything in `quarto/applications/actg175/`

**File:** `dev/tasks/cc_task_oc_applied_selfcontained_2026-09-01.md` · **Issued:** 2026-09-01 by chat, commissioned by Larry
**Repo:** `~/Documents/GitHub/forestsearch` (pop-os) · **Branch:** `feature/glm-extension`
**Transport:** `~/Downloads` → `dev/tasks/`, committed first (§1).
**Follows:** `REPORT_oc_applied_evaluation_2026-08-31.md` and its document. **Larry's three revisions, 2026-08-31/09-01:** (1) the evaluation must be **self-contained at the analysis stage** — it lives in `quarto/applications/actg175/`, computes everything itself, and reads nothing from `dev/`; (2) it must **illustrate the fundamental use of the wrappers** — the key calls visible and few, not buried in scratch scripts; (3) **`c2` grows with `c1`**: the practice recommendation is `c2 = 0.8·c1`, and a fixed `c2 = 5` against `c1 = 50` makes no sense — it is kept only as the liberal comparator column reconciling with the archived deep run.

---

## ⚠ CATEGORY

**No `R/` change.** Writes: **one new application document** `quarto/applications/actg175/analysis_actg175_continuous_oc.qmd`, its rendered HTML (per the directory's tracking convention) and its payload under `_payloads/`; **one bounded edit** to the existing `analysis_actg175_continuous_oc_evaluation.qmd` — a two-sentence note at the top pointing to the new self-contained document and stating that the archived run used fixed `c2 = 5` at `2e5` draws; nothing else in it changes and it is **not re-rendered**. Plus this task document.

**Compute — the go/no-go, stated.** The render is the compute: at the default parameters below (draws = 2e4, six rungs concurrent via `parallel::mclapply`), per rung ≈ enumerate 43 s + one draw set ≈ 400 s + ≈ 160 grid points × ≈ 6 s ≈ 22–25 min, run six-wide → **projected render ≈ 30–45 min, ≈ 5 GB peak**. *GATE:* time a single rung first (a one-rung probe render or a scratch call of the same chunk code); if the six-rung projection exceeds **2 h**, stop and report. State the measured render wall-clock in the report.

**Unattended.** Gates stop the task; never ask, never work around. **Provenance gates are dirt-tolerant:** list unrelated untracked/modified paths, leave them alone; stop only if dirt touches this task's own paths.

---

## 1. Provenance and first commit — GATE

Standard block plus the installed version (expect **0.3.2**). *GATE:* branch `feature/glm-extension`; commits `209a8e85` (stage 0), `eb136a35` (0.3.2) and `d597f289` (the evaluation report) in the log. Copy this document into `dev/tasks/` and commit it alone; the `~/Downloads` stem arrives with hyphens stripped, so report both names.

## 2. The document — design spec

**YAML `params`** (so Larry can re-render at other precisions without editing code): `draws` (default `2e4`), `n_workers` (default 6), `c2_ratio` (default 0.8), `seed` (default 8316951), and the rung and ladder vectors below as defaults in the setup chunk.

**Visible, minimal wrapper chain — this is revision (2) and the document's point.** Each of these is its own short, echoed chunk with a sentence of prose, in this order:

1. **Data prep** exactly as the compare-all document (arms 1/3, `y_decline`, twelve confounders; N = 1083) — the template document's form.
2. **The analysis**: the fixed-family `maxeffCons` anchor `forestsearch()` call, run in-document (≈ 8 s, `parallel_args = list(plan = "sequential")`), with one sentence on why the OC evaluation attaches to this configuration (the MR-anchor rationale). Print Ĥ, n(Ĥ), `T̂_obs`, `p.consistency`. *GATE:* Ĥ = `{age <= 37} & !{cd40 <= 507}`, n = 66, `T̂_obs = 87.916667` — matching stage 0.
3. **The anchored truth**: one `generate_glm_dgm()` call (the clause mapping visible: `age = 37`; `cd40 = list(type = "greater", value = 507)`; `k_inter = q + 26.978725`) and one `fs_dgm_scale()` print for the top rung — three lines of prose on the two-point mixture stylization.
4. **The family**: one `fs_oc_family_enumerate()` call with its own `print()` (M = 4508, the stage counts) — *GATE:* M = 4508 — and two sentences on Q's representability (purity 0.9969, not forced).
5. **The fundamental call**: one plain `fs_oc_predict(family = fam, n = 1083, c1 = 10, c2 = params$c2_ratio * 10, ...)` at the analyst's screening floor with the scaled consistency floor, its `print()` shown — "this one call is the wrapper; everything below is this call swept over truths and thresholds."
6. **The evaluation loop**: ~15 visible lines — `mclapply` over the rungs, each rung one `fs_oc_grid()` call whose `c1` and `c2` vectors are below; no helper files, no scratch, results held in the session.

**Rungs (six):** `q ∈ {0.01, 10, 20, 40, 60, 87.916667}` — the zero-plus structural null (with the same one-paragraph orientation explanation as the archived document, condensed), the analyst's `c1` floor, two intermediates, and the observed statistic.

**Thresholds:** `c1_ladder = {5, 10, 15, 20, 30, 40, 50, 60, 80, 100, 120, 87.916667}` (sorted; `T̂_obs` included for the calibration read-off); `c2` vector = `{5, params$c2_ratio × each c1_ladder value}` (deduplicated) — one `fs_oc_grid()` per rung evaluates the full crossing on one draw set (~160–170 points). From the crossing, three views:

- **The diagonal `c2 = 0.8·c1`** — the headline everywhere: the type-I curve (zero-plus rung) along the diagonal with `c1_05^diag` and `c1_10^diag` at ladder resolution; the ladder table (power at the diagonal fixed-type-I thresholds; rate at the analyst's `(10, 8)`); the declared-rule table (`E|Ĥ|`, PPV, sens, spec, `E[β(Ĥ)]`, naive bias, `mass_below`) at `(10, 8)`; the calibration column `P(T ≥ T̂_obs | q)` at the diagonal point nearest `T̂_obs`.
- **The `c2 = 5` column** — one comparator table row-set at the analyst's original gate, with one sentence reconciling to the archived deep run (numbers will differ at MC precision and ladder resolution; say so).
- **The type-I surface** — at the zero-plus rung only, a small heatmap or contour of the false-declaration rate over the `(c1, c2)` crossing, with the `c2 = 5` row and the `0.8·c1` diagonal marked: the picture of revision (3). Expected and worth one sentence: along the diagonal the consistency screen binds (for noisy candidates `0.8·c1 + z_p·se_g` overtakes `c1`), so type-I control arrives at much lower `c1` and selection shifts toward larger rules — report what the numbers actually show, not this expectation.

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
