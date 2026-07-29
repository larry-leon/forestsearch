# forestsearch efficiency evaluation (measurement only)

> **Start here: [`CC_BRIEF.md`](CC_BRIEF.md).** That is the task specification
> Claude Code executes. Everything else in this directory is *scaffolding it
> validates and runs* — parse-checked and partly smoke-tested, but not trusted
> code. Run it with:
>
> ```
> claude "Read dev/efficiency-eval/CC_BRIEF.md and execute it."
> ```
>
> The file you are reading is the human-facing summary of what gets measured.

A **read-only** evaluation harness. It quantifies candidate efficiency changes
to `forestsearch()` **without modifying a single line of the package**, so the
incorporate / don't-incorporate decision can be made from measured numbers
rather than from argument.

## The read-only contract

`R/00_guard.R` SHA-256 hashes every file under the installed package's `R/`
source directory (and the git worktree, if the eval is run from inside a clone)
**before and after** the run, and hard-`stop()`s if anything changed. The guard
result is recorded in the output so a reviewer can confirm the contract held.

Candidate implementations live entirely inside this harness
(`R/01_candidates.R`). They are compared against the installed package's own
functions, reached via `forestsearch:::`. Nothing is written back.

## What it measures

| ID | Question | Why it needs your hardware |
|----|----------|---------------------------|
| **G** | Are the candidates numerically equivalent? | Gate. Speedups are meaningless until this passes. |
| **M1** | `coxph()` vs `coxph.fit()` per fit | Sandbox-measurable, but confirm on your BLAS. |
| **M2** | Search scorer with/without `summary.survfit` medians | As above |
| **M3** | Combination filtering: per-combination vs precomputed | As above |
| **M4** | GLM: `glm()` vs `glm.fit()` | As above |
| **P1** | Worker-pool spawn cost as a function of W | **Cannot be measured below ~127 cores** |
| **P2** | `batch_size` sweep under the default `stop_threshold = 0.95` | **The key unknown.** See below. |
| **P3** | Consistency-stage scaling curve vs worker count | Reveals whether the stage actually scales |
| **P4** | Closure serialization size / cost per batch | Scales with `nrow(df)`; matters at `batch_size = 1` |
| **E1** | Real `forestsearch()` wall-clock on GBSG, `Rprof` attributed by stage | Replaces the analytic cost model with measurement |

### Why P2 is the headline

`forestsearch()` defaults `stop_threshold = 0.95` (`forestsearch_main.R:887`).
With the default `sg_focus = "hr"`, `subgroup_consistency_main.R:837-840` then
sets `batch_size_parallel <- 1L`. Each `future_lapply()` call therefore carries
a **single** candidate: on a 127-core machine 126 workers idle, and the
evaluation closure — which captures `df`, `found.hrs`, `index.Z`,
`confs_labels` — is serialized once per candidate.

This is a deliberate tradeoff, not a bug: `batch_size = 1` gives the finest
possible early stopping, and if the threshold fires at candidate 3 it is
clearly correct. The open question is empirical — **where does early stopping
actually fire on your data, and does the serialization + idle-worker cost
exceed the savings?** P2 answers that by sweeping `batch_size` through the
package's own public `parallel_args`, changing no code.

## Running it

```r
# from this directory
source("run_eval.R")            # full run
source("run_eval.R"); run_eval(stages = c("G", "M"))   # gates + micro only
```

Stages are independent and each writes `results/<stage>.rds`, so long stages
can be run separately. `report.qmd` assembles whatever is present.

Expected runtime: gates + micro a few minutes; P and E stages depend on your
`forestsearch()` settings — budget an hour for a full sweep.

## Configuration

Edit `config.R`. Defaults target GBSG (survival/Cox). Set `run_glm = TRUE` to
add the ACTG175 binary/OR path if `speff2trial` is installed.

## What this harness deliberately does NOT do

- It does not modify, patch, or shim the package.
- It does not claim the candidates are correct on paths it did not test.
  The gates cover the **unadjusted survival path** and, optionally, the
  unadjusted GLM path. `adjust_covariates` (especially `strata()`) is **out of
  scope** — the proposed fast Cox path does not apply there and the harness
  does not pretend to have checked it.
- It does not measure DINA or GRF. Those identifiers have their own machinery;
  the findings here concern the FS search and consistency stages only.
- It reports what it could not measure rather than extrapolating.
