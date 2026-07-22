---
title: "Guo & He (2021) §5.2 reproduction (Table 7) — Claude Code brief"
bibliography: []
---

# Claude Code brief: reproduce Guo & He (2021) §5.2, Table 7

## Purpose

Reproduce Table 7 of Guo & He (2021): empirical coverage of the 95% lower
bound of $\gamma_s$ (the best *selected* subgroup effect) in the post-hoc
identified, nested-family setting, across $\beta_2 \in \{0, 1/10, \ldots,
5/10\}$ and $r \in \{1/3, 1/12, 1/21, 1/30\}$ plus Naive. The §5.1
reproduction (Tables 3–6, disjoint families) is complete and committed; §5.2
is the nested regime closest to the `forestsearch` maxeff family on GBSG and
is the remaining validation target.

Reference: Guo, X. and He, X. (2021). *Inference on Selected Subgroups in
Clinical Trials.* JASA 116(535), 1498–1506. doi:10.1080/01621459.2020.1740096

This brief supersedes the session hand-off of 2026-07-21 where they differ.
Three items were settled in that session and are carried here as ground truth:

1. **Paper text verified.** §5.2 of `GuoHe_2021.pdf` was extracted and read
   verbatim. The DGM, the censoring inheritance from §5.1, the family
   $S(c) = \{W \le c\}$ on $c \in [30, 60]$, the step $b(\cdot)$, the 2000
   Monte Carlo samples, the absence of an Adaptive column, and every Table 7
   cell were confirmed against the PDF. Do not re-derive these; the targets
   in §"Published Table 7" below are transcription-verified.
2. **Orientation resolved** (§"Orientation" below). `orient = +1`. This is
   the single highest-risk design point; read that section before writing any
   code.
3. **Truth-curve module delivered and gate-cleared** at smoke scale
   (`guohe_sec52_truth.R`; evidence in §"Truth curve"). The main
   implementation risk named in the hand-off is therefore already retired.

## Files

```
quarto/GuoHe/guohe_sec52_truth.R      # DELIVERED, gate-cleared. Do not modify.
quarto/GuoHe/guohe_sec52_sim.R        # NEW: harness (this brief, §Harness)
quarto/GuoHe/guohe_sec52_run.R        # NEW: 127-core driver (§Compute)
quarto/GuoHe/guohe_sec52.qmd          # NEW: published-vs-reproduced report
```

Pattern sources, unchanged: `guohe_reproduction_sim.R`, `_run.R`, `.qmd`
(the §5.1 trio — mirror their structure, seeding, restartability, and
reporting conventions). Package API under test, unchanged:
`guohe_algorithm3()` via `library(forestsearch)` (exported at `d516a92` on
`feature/glm-extension`). `guohe_adaptive_r()` is **not** used — Table 7 has
no Adaptive column. No package code changes of any kind; if
`guohe_sec52_truth.R` appears to need a change, recommend it to the owner
first (recommend-then-implement), do not edit it unprompted.

Outputs: `guohe_sec52_truth_beta2_XX.rds` (truth caches, written by the
truth script), `guohe_repro_t7_beta2_XX.rds` (scenario results,
`XX = 00..05`, continuing the §5.1 `t35_`/`t6_` naming family).

## Orientation — read before coding

§2.1 of the paper sets larger $\beta$ = better effect *without loss of
generality*. §5.2 pins the direction for the nested family:

> "As $\beta_2$ increases, the subgroups are farther away from homogeneity,
> and **the best subgroup, S(30)**, is more distinctive from the others."

Under the DGM as printed (hazard $\lambda_0(t)e^{b(W)D}$), the within-$S(c)$
Cox coefficient attains its maximum $\beta(30) = \beta_2$ at $c = 30$ and
dilutes monotonically toward zero as $c$ grows. For the paper's sentence to
hold, selection must be the argmax of the **raw** coefficient. Operationally:

- `guohe_algorithm3(..., orient = +1)` in every call.
- The truth curve, the selection, the bounds, and the coverage indicator all
  live on the raw coefficient scale. **No negation appears anywhere in
  §5.2.**
- The score/report map at `orient = +1` is `report = exp(score)`; its exact
  inverse is `gh52_to_score <- function(report) log(report)`. Verify the
  round trip to 1e-12 in the harness verification block (the §5.1 V4
  pattern).
- **Do not copy `orient = -1` from `gh_one_rep()`**, and do not "reuse"
  `gh_to_score()` at its default. §5.1 used the negated convention; for
  *disjoint* families the sign is unidentified by the published tables
  (relabeling symmetry — either sign reproduces Tables 3–6), which is
  exactly why §5.1's success could never adjudicate it. The nesting breaks
  the symmetry: under negation, $S(30)$ becomes the *worst* candidate and
  selection drifts to $S(60)$, inverting the design. Any earlier session
  note framing §5.2 on the negated scale is superseded by this section.

A correct implementation shows selection concentrating near $c = 30$ as
$\beta_2$ grows (verification gate V4 below). If it concentrates near 60,
the sign is wrong.

## Design (stated in the paper; verified 2026-07-21)

- $n = 400$; hazard $\lambda_0(t)e^{b(W)D}$, $\lambda_0$ = Weibull(1,1)
  (unit exponential); $D \sim \text{Bernoulli}(0.5)$, $D \perp W$,
  $W \sim \text{Unif}[0, 80]$.
- Censoring inherited verbatim from §5.1: $\log C \sim \text{Unif}(-1.25,
  1.00)$; "about 40% across different choices of $b(\cdot)$".
- $b(w) = \beta_1 = 0$ for $w > 30$; $b(w) = \beta_2$ for $w \le 30$;
  $\beta_2 \in \{0, 1, 2, 3, 4, 5\}/10$.
- Family $S(c) = \{W \le c\}$, $c \in [30, 60]$; $\beta(c)$ the Lin–Wei
  estimand of the treatment-only Cox fit on $S(c)$ (well defined under
  misspecification; "a weighted average of $b(\cdot)$ in the range
  $[0, c]$").
- 2000 Monte Carlo samples; $r$ grid $\{1/3, 1/12, 1/21, 1/30\}$; Naive
  comparator; **no Adaptive column**.

The DGM and the analytic censoring mixture are implemented in
`guohe_sec52_truth.R` (`gh52_sim_data()`, `gh52_censoring_analytic()`);
the harness must use those, not re-implementations.

## Published Table 7 (transcription-verified targets)

Empirical coverage of the 95% lower bound of $\gamma_s$:

| $\beta_2$ | $r=1/3$ | $1/12$ | $1/21$ | $1/30$ | Naive |
|---|---|---|---|---|---|
| 0    | 0.947 | 0.961 | 0.962 | 0.962 | 0.872 |
| 1/10 | 0.960 | 0.972 | 0.972 | 0.972 | 0.879 |
| 2/10 | 0.958 | 0.966 | 0.967 | 0.967 | 0.890 |
| 3/10 | 0.959 | 0.969 | 0.970 | 0.970 | 0.895 |
| 4/10 | 0.962 | 0.968 | 0.968 | 0.968 | 0.906 |
| 5/10 | 0.964 | 0.972 | 0.973 | 0.973 | 0.901 |

Embed these in `guohe_sec52_sim.R` as `GH52_TABLE_7` (the qmd sources them
from there — never retyped, the §5.1 principle). Table 7 publishes coverage
only; distance and bias are computed and reported as ours-only diagnostics.
Note the proposed columns run **conservative** (0.947–0.973 vs nominal 0.95)
— expected from the nesting, since adjacent candidates differ by a single
subject and the effective number of independent comparisons is far below the
family size. A reproduction landing near-nominal is a red flag for the
family construction (most likely a hard-coded grid).

## Candidate family — data-determined, never a fixed grid

Per replicate, the distinct subgroups are indexed by cutpoints
$\{30\} \cup \{W_i : W_i \in (30, 60]\}$ — the indicator
$\mathbb{1}\{W_i \le c\}$ changes only as $c$ crosses an observed $W_i$, so
the supremum over $c \in [30,60]$ is an exact maximum over these
order-statistic cutpoints. Materialize each as a 0/1 membership column
(`cand_001`, …) appended to the replicate's data frame; pass the column
names as `candidates`; keep a cutpoint vector aligned to the names for
scoring. Expected family size $1 + 400 \times 3/8 = 151$, **random across
replicates**. Do NOT substitute a fixed integer grid — a ~31-point grid is a
~5× smaller family and understates selection bias, most visibly in the Naive
column (Table 6: Naive 0.900 at $k=2$ → 0.543 at $k=12$). All candidates are
large ($|S(30)| \approx 150$ to $|S(60)| \approx 300$); estimability
failures should be rare and are handled by the engine's `min_events` guard.

## Estimand and scoring

$\gamma_s = \beta(\hat c)$: the truth curve evaluated at the replicate's
selected cutpoint, via `gh52_truth_at(truth, c_hat)` (default `"smooth"`,
the isotonic projection). Coverage indicator per column: `lower <= gamma_s`.

- The engine's selection is $r$-invariant (it happens on full-sample scores
  before $r$ enters); assert this once per replicate and record the selected
  cutpoint once.
- The **naive comparator is an independent calculation** (the §5.1
  discipline): its own `coxph` per-candidate estimate/SE loop on the full
  sample, its own argmax, one-sided lower bound
  $\hat\beta(\hat c) - z_{0.95}\,\mathrm{SE}(\hat c)$. Record both the
  engine's and the independent selection plus an agreement flag
  (`sel_agree`); score each against $\gamma_s$ at its own $\hat c$. Verify
  engine-vs-independent agreement of the naive point to 1e-10 in the
  verification block (§5.1 V3 pattern).

## Harness (`guohe_sec52_sim.R`)

Self-contained apart from `library(survival)`, `library(forestsearch)`, and
`source("guohe_sec52_truth.R")` (path-robust, the `--file=` trick from the
§5.1 run script). Mirror the §5.1 harness structure:

- `GH52_R_GRID <- c(1/3, 1/12, 1/21, 1/30)`; `GH52_TABLE_7` as above.
- `gh52_candidates(df)` — cutpoints and membership columns as specified.
- `gh52_subgroup_fits(df, cand)` — per-candidate raw coefficient and SE.
- `gh52_naive(...)` — independent comparator.
- `gh52_one_rep(beta2, n, truth, r_grid, B, level, min_events, seed)`:
  simulate; build candidates; independent fits and naive; then one
  `guohe_algorithm3()` call per $r$ with **common random numbers**:
  `boot_seed <- seed + 500000L` passed as `seed =` to every $r$ (the engine
  draws all resample indices under `seed` before $r$ enters — property
  verified in §5.1; it removes between-column Monte Carlo noise, which
  matters because adjacent $r$ columns differ by as little as 0.000–0.001).
  Record per replicate: `beta2`, `n_cand`, `cens_rate`, `c_hat_gh`,
  `c_hat_naive`, `sel_agree`, `n_sel`, `gamma_s`, `naive_{point, lower,
  cover, dist, bias}`, and `r{i}_{cover, dist, bias}` with
  `bias = debiased - gamma_s`, `dist = gamma_s - lower`.
- `gh52_run_scenario()`, `gh52_summarise()` — §5.1 patterns (per-column
  means plus `mcse_cover`).
- Executable verification block (`sys.nframe() == 0L`), fail-loud PASS/FAIL
  lines, non-zero exit on failure:
  - V1 DGM: censoring vs `gh52_censoring_analytic()` at $\beta_2 = 0$ and
    $0.5$, large $n$ (tolerance 0.005).
  - V2 family size: mean `n_cand` $\approx 151$ over many draws, with
    sample-to-sample variation (fail if constant or $\approx 31$).
  - V3 engine agreement, one replicate: `selected` equals the independent
    argmax; engine naive point equals the independent `coxph` maximum to
    1e-10; report/score round trip exact at `orient = +1`.
  - V4 orientation direction: mean selected $\hat c$ at $\beta_2 = 0.5$
    strictly below mean selected $\hat c$ at $\beta_2 = 0$ (modest
    replicate counts suffice; also print the quantiles). Concentration near
    60 at $\beta_2 = 0.5$ means the sign is wrong.
  - V5 timing: one replicate at the pilot `B` and at `B = 2000`, with the
    single-core extrapolation printout (the §5.1 V6 pattern).

## Truth curve (`guohe_sec52_truth.R`) — delivered; contract and evidence

Functions: `gh52_sim_data(beta2, n)`; `gh52_censoring_analytic(beta2)`;
`gh52_truth_curve(beta2, n_big = 2e6, c_step = 0.25, seed, ...)` returning a
`gh52_truth` object (grids, `beta_raw`, `beta_smooth`, per-point `se` and
event counts, censoring diagnostics, gate table, bookkeeping) and
hard-stopping on any gate failure; `gh52_truth_at(truth, c_hat,
use = c("smooth", "raw"))`. Gates: G1 censoring vs analytic; G2 anchor
$|\hat\beta(30) - \beta_2| \le 4\,\mathrm{SE}(30)$; G3 max single positive
increment $\le 0.02$; G4 end-to-end rise $\le \max(0.01, 2\,\mathrm{SE}(30))$
(catches the gently-increasing curve of an accidental negation, which
per-step checks miss); G5 ($\beta_2 = 0$) null band; G6 ($\beta_2 > 0$)
dilution drop $\ge 0.15\,\beta_2$.

**Production step 0** (before any coverage run):
`Rscript quarto/GuoHe/guohe_sec52_truth.R` — defaults $n_{big} = 2\times10^6$,
`c_step = 0.25`, seeds $20260721 + 100\beta_2$; writes per-scenario caches,
skip-if-exists, `--force` to recompute; exits non-zero on any gate failure.
Roughly 6 min per scenario serial at ~0.3 s/fit-per-2e5-rows scaling;
embarrassingly parallel if wanted. The run driver loads the caches and must
refuse to proceed without them.

Smoke evidence (container, `--nbig=200000 --cstep=1`, all gates PASS, all
six scenarios):

| $\beta_2$ | $\hat\beta(30)$ | SE(30) | anchor $z$ | drop 30→60 | max rise | cens emp | cens analytic |
|---|---|---|---|---|---|---|---|
| 0.0 | −0.0140 | 0.0095 | −1.48 | −0.005 | 0.0035 | 0.4081 | 0.4095 |
| 0.1 | 0.0882 | 0.0094 | −1.26 | 0.039 | 0.0013 | 0.4058 | 0.4038 |
| 0.2 | 0.1863 | 0.0093 | −1.47 | 0.098 | 0.0000 | 0.4000 | 0.3982 |
| 0.3 | 0.2967 | 0.0092 | −0.36 | 0.148 | −0.0011 | 0.3945 | 0.3926 |
| 0.4 | 0.4053 | 0.0092 | 0.58 | 0.201 | −0.0023 | 0.3899 | 0.3873 |
| 0.5 | 0.4921 | 0.0091 | −0.87 | 0.248 | −0.0030 | 0.3821 | 0.3821 |

Dispositions: the dilution drop is $\approx 0.49\,\beta_2$, linear —
consistent with the weighted-average reading; single-step rises are all
$\le 0.0035$, an order of magnitude under the structural gate; censoring
tracks the analytic mixture within 0.003 everywhere. The six anchor
$z$-values are individually unremarkable but combine to $z \approx -2.1$
across the six independent draws — consistent with chance, and
re-adjudicated at production $n_{big}$ where SE(30) $\approx 0.003$. If the
production anchors also sit systematically low, stop and report; do not
proceed on a failed or suspicious truth curve.

## Compute protocol — pilot before production

Per replicate the engine evaluates all ~151 candidates on the original
sample and on each of $B$ resamples, once per $r$ value (CRN duplicates the
bootstrap work across the four $r$'s; §5.1 accepted this for fidelity — the
engine is the object under test). Roughly $4 \times (B+1) \times 151$ lean
Cox fits per replicate on subgroups of 150–300 subjects.

1. Pilot: one scenario ($\beta_2 = 0$, the maximal-bias case), `--reps=20`,
   the intended `B`, timed. Print the projected wall-clock for
   $6 \times 2000$ at the chosen core count. **Report the projection to the
   owner before launching production.**
2. Levers, in order: reduce `B` (paper silent; 2000 is the §5.1 inference;
   500 mainly costs bootstrap resolution); reduce `--reps` (500 gives
   coverage MCSE ≈ 0.010 — report the MCSE used). **Never thin the cutpoint
   family** — family size is a design parameter, not a tuning knob.
   (Fallback only if the projection is impractical and the owner approves:
   the §5.1-verified single-score-matrix reconstruction of all four $r$
   bounds from one engine call — agreement 0.000e+00 was demonstrated there
   — but prefer engine-per-$r$ fidelity.)
3. Driver mechanics, mirroring `guohe_reproduction_run.R`: parallelism over
   replicates only (`mclapply`, `mc.preschedule = FALSE`), engine
   `parallel = FALSE`; `RNGkind("Mersenne-Twister", "Inversion",
   "Rejection")`; scenario seed base from the id hash, replicate seed
   `base + m` so serial and parallel runs agree bit-for-bit; per-scenario
   `.rds` on completion, skip-if-exists, `--force`; record `B`, `r_grid`,
   seed base, truth-cache identity ($n_{big}$, `c_step`, seed), elapsed,
   `sessionInfo`.
4. `devtools::install()` (not `load_all()`) before any parallel run —
   workers see only the installed package.

## Report (`guohe_sec52.qmd`)

Mirror `guohe_reproduction.qmd`: source the sim file (targets never
retyped); per-cell comparison with the **combined** Monte Carlo SE,
$\sqrt{\mathrm{SE}^2_{\text{repro}} + \hat p_{\text{pub}}(1-\hat
p_{\text{pub}})/2000}$; `*`/`***` flags at 2/3 combined SE; the collected
miss table; the "B is an inference, not a quotation" callout. Add,
§5.2-specific: the orientation resolution with the S(30) quote; truth-curve
figures (raw vs smoothed with SE band, anchors marked) read from the
caches; the selected-$\hat c$ distribution by $\beta_2$ (diffuse at
$\beta_2 = 0$, concentrating toward 30 — the paper's design narrative made
visible); the conservatism check (proposed columns above nominal; flag
near-nominal); ours-only distance/bias tables; the stated-vs-inferred
ledger.

## Stated vs inferred — record in the report

**Stated (PDF-verified):** DGM, censoring, $D = [30,60]$, $S(c)$, step
$b(\cdot)$, $\beta_1 = 0$, $\beta_2$ grid, 2000 samples, $r$ grid, no
Adaptive column, Table 7 values.

**Inferred / chosen (record as such):** `B`; the order-statistic candidate
construction (our reading of "$c \in D$ compact" plus finiteness of the
empirical family); the truth-curve numerical method, $n_{big}$, `c_step`,
isotonic smoothing, and gate tolerances (calibration provenance in the
truth file's documentation); the `orient = +1` resolution (anchored to the
paper's S(30) sentence; the sign-flipped protective DGM with `orient = -1`
is the symmetric alternative, differing only through arm-level
event-rate asymmetry — we take the printed form); any reduction of
replicates from 2000.

## Conventions

Formal academic register; LaTeX math in `$...$`; "consistent with" /
"corroborated by" for simulation findings — never stronger. Quarto-first;
no `.qmd` ↔ `.Rmd` conversion. Recommend before implementing; show the plan
and the pilot projection before long runs. Every response that creates or
modifies a file ends by presenting the files. Conventional Commits on
`feature/glm-extension`; the owner commits and pushes via GitHub Desktop —
no PRs. Scope: the three new analysis files only; no package changes; do
not modify `guohe_sec52_truth.R` or any §5.1 file. Report failures loudly;
never fill a missing number with a plausible one.

## Open items for the owner (surface, do not decide silently)

1. Production `B` (recommend after the pilot; 2000 for §5.1-consistency vs
   500 for tractability).
2. Replicates: 2000 (published) vs a reported-MCSE reduction if the
   projection demands it.
3. Whether the qmd's ours-only distance/bias diagnostics should also be
   summarised in the companion-paper hand-off (recommend yes).
