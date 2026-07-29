# FINDINGS — forestsearch efficiency evaluation

**Date:** 2026-07-28  **Machine:** 128 physical cores, R 4.6.1,
**reference `libblas`** (unoptimised), forestsearch 0.2.0 installed from
`feature/glm-extension` @ `75987a9`.
**Read-only contract:** held — 138 package source files verified unchanged
before/after every stage; `git status` shows changes only under
`dev/efficiency-eval/`.

All six candidates are **numerically equivalent** to the package functions
they would replace (every gate, including G4 and G5 which had never run, and
G6, passed exactly). So the decision is **not** about correctness of the
candidates — it is about whether the runtime they touch is worth touching.
The profile says: **mostly not.**

---

## 0. Headline (a correctness finding that outranks every timing here)

**`batch_size` — a parallelisation tuning knob — silently changes the selected
subgroup.** This is exactly the P2 invariance check the brief flags as
taking precedence over all timing, and it **failed**:

| `batch_size` | candidates evaluated | selected subgroup |
|---|---|---|
| **1** (the package default for `sg_focus="hr"` + `stop_threshold`) | 1 | **`{er ≤ 0} & {size ≤ 35}`** |
| 12, 51, 102, 204 | 8 | **`{er ≤ 0} & {pgr ≤ 33}`** |

Deterministic — identical within each batch size across all 3 repeats, differing
across. Not noise. Mechanism (verified by inspecting `grp.consistency`):

- Candidates are **HR-sorted**, then evaluated. Candidate #1 is
  `{er≤0}&{size≤35}` (HR = 2.537, the **maximum** HR; Pcons = 0.96).
- `stop_threshold = 0.95` fires on candidate #1 (0.96 ≥ 0.95). At
  `batch_size = 1` the loop stops **inside the first batch**, so
  `select_best_subgroup()` only ever sees candidate #1 and returns it.
- At any `batch_size` large enough to hold the whole pool (here ≥ 8), all 8
  candidates are evaluated and stored **before** the early-stop `break`. Six
  pass. `select_best_subgroup()` then runs over the full set and — under the
  package-default `selection_rule = "neighborhood"` — returns the
  **highest-consistency** subgroup, `{er≤0}&{pgr≤33}` (Pcons = 1.00, HR = 2.222),
  which is *not* candidate #1.

So the package has two selection semantics — "first HR-sorted passer"
(early-stop path) and "best over the full passing set" (`select_best_subgroup`)
— and `batch_size` decides which one the user gets. The default
(`batch_size = 1`) short-circuits the very `selection_rule` the function
otherwise advertises as its default.

This is a **latent bug in the package**, independent of any efficiency
candidate. It is not "fixed" here (out of scope, and the brief forbids it).
Two internal references frame the intent — `subgroup_consistency_main.R:838`
("For hr or minSG, first passing candidate is optimal") and `:943` ("process
in order to respect HR sorting") — but `select_best_subgroup()` with
`selection_rule = "neighborhood"` (`:506`, `in_band = hr ≥ (1-ε)·max(hr)`)
does **not** implement "first passing candidate," so the two paths genuinely
disagree. **Recommend the maintainer reconcile them** (either have the
early-stop path also apply `select_best_subgroup` over what it has seen, or
document that `stop_threshold` overrides `selection_rule`).

---

## 1. Per-candidate decision table

| Candidate | Gate | Unit-cost speedup | End-to-end effect (measured) | **Call** |
|---|---|---|---|---|
| **C1 Fast Cox engine** (`coxph.fit`) | G1,G2 exact | **12.1–14.5×** / fit | `coxph` is **8.9%** of runtime → **≈1.09×** overall (E7) | **Defer** — adopt only bundled into C3 |
| **C2 Exact futility stopping** | G3a,G3b exact | 2.6× *fewer fits* (prior, synthetic) | consistency is **6.75%** of runtime and the default already stops at candidate 1 → **no measurable win** | **Reject** (as a perf change) |
| **C3 Scorer w/o `summary.survfit`** | G4a,G4b exact | **3.46×** (scorer, incl. C1) | scorer ≈ **27%** of runtime → **≈1.2×** overall bundled with C1 | **Adopt** — best safe non-GRF lever |
| **C4 Vectorised pre-screen** | G5 exact | 2.3× / 3.2× (filters) | combinatorial filters are a **negligible** profile share (cost is the Cox fits, not the filters) | **Reject** |
| **C5 GLM via `glm.fit`** | G6 exact | **2.38×** / fit | GLM end-to-end path **not run** (no ACTG175); share unmeasured | **Defer** |
| **batch_size default** | — | see §0 | default is *faster* (8.24 vs 9.40 s) but **changes the answer** | **Escalate** (correctness, §0) |
| **`plan()` hygiene** in `subgroup_search.R` | — | — | not triggered from `forestsearch()` (see §3.4) | **Report only** |

**Bottom line for the maintainer:** none of C1–C5 moves the end-to-end number
by much, because **candidate construction (GRF ≈ 37%, LASSO ≈ 7%) and the
per-candidate search scorer (≈ 27%) dominate — the consistency stage the
candidates were built around is only ~7% of runtime on GBSG.** The one worth
doing is **C3 (bundling C1)**: a single rewrite of `fit_cox_for_subgroup` to
`coxph.fit` + a direct `survfit` median read, worth ~1.2× overall and exact.
Everything else is either a rounding error end-to-end (C2, C4) or unmeasured
(C5). The **highest-value change is not in this candidate set at all** — it is
GRF construction, and separately the correctness issue in §0.

---

## 2. Stage-by-stage evidence

### Gates (all pass, exact)
`coxph.fit == coxph`: max |Δβ| = 0 over 40 fits incl. tied times. Pcons
bit-identical under engine swap. Futility: decision identical and the reported
value is provably an upper bound (never wrongly rejects). Scorer HR/medians/CI
match to 0 (the package's `summary.coxph` uses `qnorm(0.975)`, so the
"1.96 vs qnorm" caveat produced **no** difference). Pre-screen selects the
identical combination set. `glm.fit == glm + vcov` exact.

### Micro (M) — reference BLAS
- **M1 fast Cox:** 12.1–14.5× (N = 60–500). *Confirms* prior #1; if anything a
  touch higher — this box runs unoptimised `libblas`, so the dispatch overhead
  `coxph.fit` removes is proportionally larger. On an optimised BLAS expect the
  low end of the range.
- **M2 scorer (C1+C3 combined):** 3.46×. **M3 pre-screen:** 2.3× (k=2) / 3.2×
  (k=3). **M4 GLM:** 2.38×.
- No `NaN`/`Inf`/`0 ms` — the auto-calibrating timer held.

### Parallel (P) — the decision-relevant stage
- **P2 (headline):** wall-clock median of 3 — `bs=1`: **8.24 s** (evaluates 1
  candidate); `bs≥12`: **~9.40 s** (evaluates 8). The default is *faster* here,
  because early stopping fires immediately. **But the selected subgroup is not
  invariant (§0).**
- **P1 pool spawn** grows with W: 0.02 s (W=1) → **1.57 s (W=102)**.
- **P3 scaling:** wall-clock **increases** with workers — 2.42 s (W=1) → 9.40 s
  (W=102). With only **8 candidates**, a multisession pool is pure overhead;
  sequential is fastest. The consistency stage does **not** scale on GBSG — it
  anti-scales.
- **P4 serialization:** 0.04 MB/dispatch, ≈ 0 s. The per-candidate closure is
  **not** the bottleneck — **pool spawn is** (P1). This refutes the framing that
  `batch_size=1` is expensive because of serialization.

### End-to-end (E) — sequential `Rprof`, honest attribution
Original patterns attributed **0%** to consistency and left ~50% unattributed;
corrected on-machine (real functions are `evaluate_consistency_twostage`,
`search_combinations_parallel`, and GRF/LASSO were missing rows entirely):

| stage | share |
|---|---|
| GRF candidate construction (`causal_survival_forest`) | **37.3%** |
| Subgroup search / candidate screening | **27.0%** |
| coxph (within screening) | 8.9% |
| consistency evaluation | **6.75%** |
| LASSO (`cv.glmnet`) | 6.9% |
| survfit / summary.survfit | 6.75% |
| model.frame/matrix/terms | 1.6% |

E2+E3+E3b+E3c = **78%** — a plausible majority. Projected fast-Cox-only
end-to-end speedup (Amdahl on the 8.9% coxph share at ~14×): **1.09×**.

---

## 3. What I could not measure, and why

1. **The GLM end-to-end path (P/E for C5).** `run_glm = TRUE` only exercises the
   self-contained synthetic gates/micro (G6, M4). The ACTG175 dataset needed for
   a real binary/OR `forestsearch()` run was not available, so C5's *runtime
   share* is unmeasured. Its unit cost (2.38×) is confirmed; its end-to-end
   value is a guess (likely the same story: dominated by GRF/LASSO).
2. **W = 128 (= `detectCores()`).** A 128-worker `multisession` (PSOCK) pool
   **cannot be created** on this machine: R's default connection limit is 128
   and ~4 are already in use, so node 125 fails ("only 124 connections left").
   The auto-derived sweep hit this and aborted stage P. I pinned the sweep to
   feasible values topping out at **102** = `.default_parallel_workers()`
   (`floor(0.80 × 128)`), the package's own default, which stays under the
   ceiling. **I did not extrapolate to 128.**
3. **The adjusted Cox path, DINA, GRF identifiers, continuous/count GLM,
   `use_twostage = FALSE`, bootstrap, CV** — all out of scope per the brief.
   Note the E/P runs use the **default `use_twostage = TRUE`**, so C2's futility
   numbers (measured against the single-stage split loop) do not transfer
   directly to the default pipeline, which already has its own screening.
4. **Unknown #4 (report-only):** `subgroup.search()` defaults
   `parallel_workers = parallel::detectCores()` (`subgroup_search.R:80`), which
   **bypasses** `.default_parallel_workers()`. But `forestsearch()` hard-codes
   `parallel_workers = 1` when it calls the search (`forestsearch_main.R:2410`),
   so the search stage always runs **sequentially** and the `detectCores()`
   default is never reached *from `forestsearch()`*. A **direct**
   `subgroup.search()` caller on this box would hit the same W=128 PSOCK failure
   as in (2). Worth changing the default to `availableCores()`; not changed here.

---

## 4. Every prior in §6 of the brief, adjudicated

| # | Prior | Verdict |
|---|---|---|
| 1 | `coxph.fit` ≈ 11–14× faster, identical coeffs | **Confirmed** (12.1–14.5×, exact; reference BLAS inflates it slightly) |
| 2 | Consistency is 95–99% of runtime | **Refuted, decisively.** Measured **6.75%**. GRF (37%) + search (27%) dominate. This inverts the whole premise. |
| 3 | Futility ≈ 2.6× fewer fits, decisions identical | Decisions identical **confirmed** (G3a/b). The 2.6× fit reduction is **untested on the real pool** and **immaterial**: the default already stops at candidate 1, and consistency is 6.75% of runtime. |
| 4 | `summary.survfit` medians ~33% of scorer, unused at `m1.threshold=Inf` | **Confirmed** (unused for selection; survfit ≈ 6.75/27 ≈ 25% of the scorer — close). |
| 5 | Pre-screen prunes 43–58% of combos | Pruning *fraction* **confirmed on the gate data** (~45–59% of combos, 25% of columns dropped) — but the end-to-end value is **refuted**: the filters are a negligible profile share; the screening cost is the Cox fits. |
| 6 | `batch_size = 1` wastes the pool at high W | **Refuted — and it is worse than "wasteful."** At `bs=1` early stopping fires at candidate 1, so it is the **fastest** option (pool simply isn't needed), **and** it **changes the selected subgroup** (§0). The default is optimal for speed and wrong for reproducibility. |

---

## 5. Recommendation

1. **Fix the `batch_size`/selection inconsistency (§0) first** — it is a
   correctness issue that the efficiency question is irrelevant beside.
2. If any efficiency work is done, do **C3 bundling C1** (one exact
   `fit_cox_for_subgroup` rewrite, ~1.2× overall). Skip **C2** and **C4**
   (no measurable end-to-end benefit). Hold **C5** until a real GLM run can
   attribute its share.
3. The real runtime is in **GRF candidate construction (37%)** and the
   **per-candidate search scorer (27%)** — neither is in this candidate set.
   Any serious speedup effort should start there (and with **not** spinning up a
   102-worker pool for 8 candidates: on GBSG, `plan = "sequential"` is 4× faster
   than W=102).
