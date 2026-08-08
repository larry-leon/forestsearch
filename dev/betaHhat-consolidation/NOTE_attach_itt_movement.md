# ITT-complement semantics for undetected replicates: what moved


`fs_attach_betaHhat()` previously gave an undetected replicate — `sg_def` `NA`
or empty — an all-`NA` record. It now gives it the no-subgroup record that
`fs_betaHhat_one()` has always returned: `nH_eval = 0`,
`nHc_eval = nrow(frame)`, `betaHhat_Hc` the ITT effect, `betaHhat_H` `NA`, and
`status = "ok"`.

An undetected replicate is not a missing measurement. It is a run in which the
whole population is the complement, so the complement target exists and is the
ITT effect. This is the partition invariant reaching layer 4: `nH_eval +
nHc_eval` now equals `nrow(frame)` on **every** row that is not
`"unresolved"`, not only on rows carrying a rule.

### Movement, measured on three bundles

Read through the new attach in scratch; no bundle was overwritten.

| bundle | rows | moved | share |
|---|---|---|---|
| knoise3 pilot (30 reps, n = 500) | 30 | 9 | 30.0% |
| knoise6 pilot (30 reps, n = 500) | 30 | 5 | 16.7% |
| committed knoise0 `fs_mr_n500` (s1000) | 1000 | **520** | **52.0%** |

Every moved row changes identically:

```
betaHhat_Hc   NA  ->  0.6719670492      (the ITT effect on the eval frame)
nHc_eval      NA  ->  100000
nH_eval       NA  ->  0
status        NA  ->  "ok"
betaHhat_H    NA  ->  NA                (unchanged: an empty region has no target)
```

The ITT value is the same in all three because the DGM is the same; the noise
columns do not enter the outcome.

**Detected rows are untouched.** `betaHhat_H` and `betaHhat_Hc` are
`identical()` to their previous values in all three bundles. The change is
confined to rows that previously had no target at all.

**Partition after the change**: 30/30, 30/30 and 1000/1000 rows satisfy
`nH_eval + nHc_eval == N`, where before it held only on rows with a rule.

### The change is not contained to re-pointed documents

The shims strip *columns* (`betaHhat_status`, `nH_eval`, `nHc_eval`) but pass
*values* through unchanged. A consumer still on `attach_betaHhat()` therefore
sees the new ITT complement on undetected rows as soon as it is re-run:

```
shim attach_betaHhat() columns: sim_id, sg_def, betaHhat_H, betaHhat_Hc
undetected rows' betaHhat_Hc via the SHIM: 0.671967, 0.671967
```

**From this commit forward, ANY consumer re-run gets the new semantics** --
the 74 batch/combine documents in `gbsg_redux/`, the frozen
`dev/identifier-alignment/` copies, and anything else still on a shim. The
change is *not* gated by the migration: it lands on re-run, not on re-point.

At committed scale that is **52% of rows in a bundle moving from `NA` to a
finite value**, and any coverage computed over `betaHhat_Hc` now counts those
rows in `n_eff` where it previously dropped them.

This is the intended semantics, deliberately chosen and documented here rather
than discovered later. An undetected replicate contributes a real complement
observation; excluding it was the defect.

### Verification

- Contract tests: **145 assertions pass** (was 132). T5 and T6 updated; the
  discriminating assertion requires undetected rows' `betaHhat_Hc` to be
  `identical()` to `fs_betaHhat_one(NULL, ...)` on the same frame, and a
  negative control requires a detected row **not** to carry the no-subgroup
  record.
- `R CMD check --as-cran`: **0 errors, 0 warnings, 0 notes** (14 m 15.6 s).


## Why this is a separate NOTE, not an appendix to the T10 gate

`T10_GATE_RESULT.md` records a specific claim: that consolidating `betaHhat`
into the package moved **no** value except the one sanctioned D1 fix. That
record must stay exactly what it is.

This is a deliberate **post-gate semantics change**, made after T10 passed and
with its movement measured on purpose. Folding it into the gate document would
blur the distinction between "the consolidation changed nothing" and "we then
chose to change something", which is the distinction the gate exists to make.
