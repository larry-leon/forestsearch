---
title: "Guo & He Section 5 reproduction — run instructions"
bibliography: []
---

# Run instructions (Claude Code / 127-core host)

## Verdict

**This is a compute-host job.** Authoring and single-replication verification
are complete and were done in the chat sandbox; the full study is 184–875
core-hours depending on one open decision (below), which is 1.7–8 hours wall
clock on the 127-core box and is not feasible in a sandbox limited to one core.

## Files

Place all of these in `quarto/GuoHe/` alongside the existing
`guohe_algorithm3.R` and `guohe_adaptive_r.R`, which are **not modified**:

| File | Role |
|---|---|
| `guohe_reproduction_sim.R` | DGM, naive comparator, one replication, scoring, published targets, and the fail-loud verification/timing block |
| `guohe_reproduction_run.R` | Parallel driver; one `.rds` per scenario; restartable |
| `guohe_reproduction_RUN.md` | This file |

## Step 1 — verify before spending anything

```bash
cd <repo>/quarto/GuoHe
Rscript guohe_reproduction_sim.R
```

Runs V1–V6 and exits non-zero on any failure. Expected output:

```
V1a PASS -- censoring 0.4080 vs analytic 0.4095 (paper: 'about 40%')
V1b PASS -- P(S1) = 0.5010, P(treat) = 0.5002, both target 0.5
V2  PASS -- fitted log-HR = (0.0005, 0.4989) vs specified (0.000, 0.500)
V3a PASS -- gamma_max: engine 0.3233635997 vs independent coxph 0.3233635997
V3b PASS -- selected S1; independent argmax S1
V4  PASS -- bound: report 1.01677606 <-> score -0.01663689
V5  PASS -- naive 0.3234 >= debiased 0.2075 > lower -0.0166
```

V1 and V2 confirm the DGM against an analytic censoring probability and against
the specified log hazard ratios. V3 confirms orientation and reproduces the
engine's $\hat\gamma_{\max}$ from an independent `coxph()` call. These are the
checks that guard against the failure mode the session directive identifies as
most likely: a mis-specified DGM rather than a broken estimator.

If any block fails on the host, **stop and diagnose** — do not tune.

## Step 2 — the run

```bash
# Tables 3, 4, 5 (k = 2, six beta_2 values)
Rscript guohe_reproduction_run.R --tables=35 --cores=120 --B=2000 --adaptive --adaptive-B=200

# Table 6 (k = 2, 6, 10, 12; n = 200k)
Rscript guohe_reproduction_run.R --tables=6  --cores=120 --B=2000 --adaptive --adaptive-B=200
```

Run them as two jobs, not one, so Tables 3–5 land while Table 6 is still going.
Each scenario writes `guohe_repro_<id>.rds` on completion and is skipped on
re-invocation unless `--force` is passed, so an interrupted job resumes.

Each `.rds` stores the per-replication results, the summary, the seed base,
`sessionInfo()`, and the elapsed time. Replications are individually
reproducible from `(scenario id, m)`.

## Measured cost

Timings are single-core in the authoring sandbox (R 4.3.3, survival 3.5.8).
Per-core speed on the host will differ; treat these as ±50%.

| Configuration | s / replication |
|---|---|
| $k=2$, $n=400$, $B=2000$, four fixed $r$ | 8.28 |
| $k=2$, $n=400$, $B=200$, adaptive marginal cost | 5.90 |
| $k=12$, $n=2400$, $B=2000$, four fixed $r$ | 60.57 |

Cost scales essentially linearly in $k$: each subgroup holds $n/k = 200$
subjects regardless of $k$, so the per-replication work is $k$ Cox fits per
resample. The $k=12$ measurement is 22% above the linear prediction from $k=2$,
which is dispatch overhead.

| Component | core-hours |
|---|---|
| Tables 3–5, fixed $r$ | 27.6 |
| Tables 3–5, adaptive at $B = 200$ | 19.7 |
| Tables 3–5, adaptive at $B = 2000$ | ~197 |
| Table 6, fixed $r$ | ~80 |
| Table 6, adaptive at $B = 200$ | ~57 |
| Table 6, adaptive at $B = 2000$ | ~570 |

| Whole study | core-hours | wall clock, 127 cores @ 85% |
|---|---|---|
| Adaptive at $B = 200$ | 184 | ~1.7 h |
| Adaptive at $B = 2000$ | 875 | ~8.1 h |

## The one open decision: $B$ inside Algorithm 2

The paper states $B$ nowhere in Section 5. The authors' own MONET1 script uses
$B = 2000$, which is the default for the fixed-$r$ columns here and is recorded
in each `.rds` as an inference rather than a quotation.

For Algorithm 2, `guohe_adaptive_r()` uses one $B$ for both the inner
cross-validation fits and the final refit, and the inner fits need only the
bias-reduced point estimate — a bootstrap *mean*, which stabilises far faster
than the 95th percentile the bounds need. Setting `--adaptive-B=200` is
therefore defensible and costs 184 core-hours; `--adaptive-B=2000` is more
conservative and costs 875.

**Recommendation: start at `--adaptive-B=200`.** If the Adaptive column is the
only one that fails to reproduce, re-run that column at 2000 before concluding
anything, because at that point $B$ becomes a live explanation.

## Known implementation choices, recorded in advance

Two decisions are recorded here so that neither is discovered after the fact and
rationalised.

**Common random numbers across the $r$ grid.** `guohe_algorithm3()` draws its
resample indices under `seed` before `r` is used (its lines 352–353), so passing
one seed to every $r$ within a replication gives all columns identical bootstrap
draws. This was verified by execution: reconstructing each bound from a single
score matrix agrees with the engine to `0.000e+00` across the whole grid. CRN is
applied because Table 3's $r$ columns differ by as little as 0.000–0.002, which
independent draws would bury in noise.

**Algorithm 2 uses independent draws across $r$.** `guohe_adaptive_r()` line 219
passes no `seed` to its inner `guohe_algorithm3()` calls, so the CV objective is
compared across $r$ on different Monte Carlo noise. Guo & He do not specify
either way. This is left as-is — it is the existing validated behaviour — but it
inflates the variance of $\hat r$, and the Adaptive column is precisely a
measurement of $\hat r$. **If Adaptive under-reproduces, examine this before
questioning the DGM.**

## Acceptance criteria

Judge against Monte Carlo error, not exact agreement. At 2000 replications the
paper reports coverage standard errors around 0.005; the same applies here.

- **Coverage**: reproduced within roughly $\pm 0.015$ (3 MCSE) of published.
- **Bias**: reproduced within roughly 3 MCSE, with the *ordering* across $r$ and
  against Naive preserved. The ordering is the stronger signal, since it is
  insensitive to any unstated $B$.
- **Table 6 trend**: the qualitative pattern must hold — naive coverage
  collapsing with $k$ (0.900 → 0.543), naive bias growing (0.105 → 0.302),
  $r = 1/3$ degrading while small $r$ holds near nominal, and Adaptive drifting
  down (0.939 → 0.925).

If the fixed-$r$ columns reproduce and Adaptive does not, that is still a
successful decoupling: the estimator and the shrinkage machinery are confirmed,
and the failure localises to Algorithm 2.

## Preliminary evidence already in hand

A 40-replication pilot at $\beta_2 = 0$, $B = 500$ (the maximal-bias
configuration, $|H| = 2$) gave, against published values in brackets:

| | $r=1/3$ | $1/12$ | $1/21$ | $1/30$ | Naive |
|---|---|---|---|---|---|
| bias | 0.045 [0.028] | 0.021 [0.008] | 0.019 [0.007] | 0.019 [0.006] | 0.120 [0.107] |
| dist | 0.223 [0.248] | 0.241 [0.265] | 0.242 [0.266] | 0.243 [0.266] | 0.188 [0.213] |
| cover | 0.875 [0.933] | 0.925 [0.950] | 0.925 [0.952] | 0.925 [0.952] | 0.850 [0.896] |

MCSE is approximately 0.025 on bias and 0.04 on coverage at 40 replications.
Every entry lies within roughly one MCSE of its published counterpart, and the
ordering across $r$ and against Naive is reproduced exactly. This is consistent
with the extracted design. It is not a reproduction; 40 replications and
$B = 500$ cannot support one.

## Deferred, deliberately

- **Table 7 (§5.2, post-hoc identified).** The truth target $\beta(c)$ is a
  weighted average of $b(\cdot)$ over $[0,c]$ rather than $b(c)$, so it requires
  separate derivation. Not needed for the tuning-sensitivity argument.
- **Table 8 (§5.3, synthetic model).** Requires supplement Tables B.1–B.2 and
  uses a different $r$ grid ($1/3, 1/9, 1/30$).
- **Zhao, Ivanova & Fine (2023).** Independent-group corroboration; a separate
  exercise once the primary reproduction lands.
- **A single-pass $r$-sweep wrapper.** Verified exact (`0.000e+00`) and worth
  roughly 4× on the fixed-$r$ columns, but not built: the compute is affordable
  without it, and adding a second code path to a validation exercise trades a
  known-good engine for an optimisation we do not need.

## After the run

The remaining deliverable is `guohe_reproduction.qmd` — the living document that
reads the `.rds` files and lays reproduced against published side by side, with
MCSE bands and a PASS/FAIL per cell. It is deliberately not written yet, since
it would render nothing until results exist.
