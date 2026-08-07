# T10 — the migration gate: result

**Verdict: PASS, with one sanctioned movement.**

15 of 16 rows bitwise identical (`diff == 0`, not a tolerance). The 16th is
D1's fix on the one family where D1 was still live.

Gate script: `T10_migration_gate.R`. Raw result: `T10_result.rds`.

## What the gate is for

The spec's scope line is binding: **no estimand changes.** Every family must
reproduce the target its module computes today, so any movement is a defect in
the consolidation until shown otherwise. T10 is the proof, not a review.

## Fixtures

| family | fixture |
|---|---|
| continuous | calibrated ACTG175 DGM from `acceptance_betaHhat_md.qmd`; `df_super` n = 5,000, MD(Q) = 40.000000 |
| binary | the DGM behind the `actg175_or075_*` `mr_sweep` bundles; `build_eval_frame_glm(eval_seed = 20260628)`, n = 100,000, OR(H) = 0.750000 |
| survival | `setup_gbsg_dgm(model = "alt", n_super = 100000, seed = 8316951)` + `build_eval_frame(analysis_time = 84, cens_adjust = log(1.5), eval_seed = 20260628)`, n = 100,000 |

Each family is scored on a rule set including a conjunction, a single cut, a
**disjunction**, a **negation**, and a region disjoint from the planted
subgroup.

## The table

### Continuous — old `betaHhat_one_md()` vs `fs_betaHhat_one()`

| rule | old_H | new_H | diff_H | old_Hc | new_Hc | diff_Hc |
|---|---|---|---|---|---|---|
| exact_Q | 40.00000000000 | 40.00000000000 | 0 | −26.25523587604 | −26.25523587604 | 0 |
| narrow | 40.00000000000 | 40.00000000000 | 0 | −11.06265474265 | −11.06265474265 | 0 |
| one_cut | 21.42971458947 | 21.42971458947 | 0 | −26.25523587604 | −26.25523587604 | 0 |
| disjunction | 26.15880058604 | 26.15880058604 | 0 | −26.25523587604 | −26.25523587604 | 0 |
| negation | −4.24382329013 | −4.24382329013 | 0 | −3.07252695054 | −3.07252695054 | 0 |
| disjoint_Q | −26.25523587604 | −26.25523587604 | 0 | 21.42971458947 | 21.42971458947 | 0 |

### Binary — old `betaHhat_one_or()` vs `fs_betaHhat_one()`

| rule | old_H | new_H | diff_H | old_Hc | new_Hc | diff_Hc |
|---|---|---|---|---|---|---|
| exact_Q | 0.740848939902 | 0.740848939902 | 0 | 0.649980364249 | 0.649980364249 | 0 |
| one_cut | 0.699410077131 | 0.699410077131 | 0 | 0.645716475038 | 0.645716475038 | 0 |
| disjunction | 0.666183731691 | 0.666183731691 | 0 | 0.655787840126 | 0.655787840126 | 0 |
| negation | 0.662529524560 | 0.662529524560 | 0 | 0.649862251202 | 0.649862251202 | 0 |
| disjoint_Q | 0.645716475038 | 0.645716475038 | 0 | 0.699410077131 | 0.699410077131 | 0 |

### Survival — old `betaHhat_one()` vs `fs_betaHhat_one()`

| rule | old_H | new_H | diff_H | old_Hc | new_Hc | diff_Hc |
|---|---|---|---|---|---|---|
| conjunction | 0.632695709116 | 0.632695709116 | 0 | 0.678488876601 | 0.678488876601 | 0 |
| one_cut | 0.669103318852 | 0.669103318852 | 0 | 0.631140528292 | 0.631140528292 | 0 |
| **disjunction** | **NA** | **0.688255552881** | **—** | **NA** | **0.641156539653** | **—** |
| negation | 0.639671477945 | 0.639671477945 | 0 | 0.720999350454 | 0.720999350454 | 0 |
| small | 0.621498832480 | 0.621498832480 | 0 | 0.671914864483 | 0.671914864483 | 0 |

Counts (`nH_eval`, `nHc_eval`) are identical on every comparable row.

## The sanctioned movement

`survival / disjunction`, rule `(er > 125 & size > 20) | (nodes > 5)`:

```
old:  betaHhat_H = NA              betaHhat_Hc = NA
      nH_eval    = NA              nHc_eval    = NA
new:  betaHhat_H = 0.688255552881  betaHhat_Hc = 0.641156539653
      nH_eval    = 42793           nHc_eval    = 57207
```

The old module emits, verbatim:

```
Warning: Column '(er' not found in data frame
Warning: Could not parse expression and column 'size > 20) | (nodes > 5)' not found in data frame
```

Those warnings **are the evidence of D1**. `betaHhat_truth.R:73-74` splits the
rule on `" & "` before testing for `"|"`, shredding
`(er > 125 & size > 20) | (nodes > 5)` into three fragments — `'(er > 125'`,
`'size > 20) | (nodes > 5'`, `'nodes > 5)'` — none of which names a real
column. The old module cannot express a disjunction at all.

### Why this is sanctioned rather than a failure

T10's invariant is **"no estimand changes, no unintended movement."**

- **15 of 16 rows bitwise identical** shows nothing moved unintentionally.
- The 16th is D1's fix on the one family where D1 was **still live**. The spec
  lists D1 as a defect being fixed and records it as "corrected in the binary
  and continuous modules at `e6f6024`; **still live in both survival copies**".
- **The binary and continuous disjunction rows at exactly 0 confirm the
  mechanism.** Both families dispatch correctly on both sides because
  `e6f6024` already fixed them, so their disjunctions agree to the bit. Only
  the family that never received that fix moves. That is the signature of a
  defect fix, not of a consolidation error.
- **Reproducing `NA` would mean shipping the defect** the spec exists to fix.
- **The new value is corroborated independently.** The same resolver is
  verified by T3 against `.grf_evaluate_subgroup()` on a realized GRF
  disjunction — identical membership, 921 of 1000 — and the partition closes
  here: `42793 + 57207 = 100000`.

## The `na.rm` masking error, and its correction

The first version of the gate script computed

```r
max(abs(c(T10$diff_H, T10$diff_Hc)), na.rm = TRUE)
```

`na.rm = TRUE` **silently dropped the one row that moved**, because that row's
diff is `NA` (old `NA` minus new finite). The script then printed:

```
max |diff| across all families: 0.000000e+00
GATE: PASS -- every value bitwise identical
```

which was false. The `counts identical: TRUE` lines were masked the same way,
by `all(..., na.rm = TRUE)`. The movement was found only by inspecting the
per-row table by hand.

**This is the same vacuous-check class as a `pgrep` that matches its own
command line: a check that cannot fail proves nothing.** A gate whose summary
drops the rows it exists to catch is worse than no gate, because it reports
confidence.

`T10_migration_gate.R` now classifies every row into exactly one of

| class | meaning | verdict |
|---|---|---|
| `identical` | both finite, `diff == 0` | ok |
| `moved` | both finite, `diff != 0` | **FAIL** unless sanctioned |
| `na_mismatch` | exactly one side `NA` | **FAIL** unless sanctioned |

with **no `na.rm` anywhere in the verdict path**, and sanctioned rows listed
explicitly in `.SANCTIONED` with their reason. A row that moves and is not on
that list fails the gate. The script also prints the old module's warnings
inline, so the D1 evidence appears in the run rather than only in this
document.

## Standing

The gate passes. Step 4 (shims) may proceed. The sanctioned row is the
expected behavioural change for the survival family and should be repeated in
the release notes when the survival modules are re-pointed, since any consumer
that previously received `NA` for a disjunction will now receive a value.

## Negative control on the gate itself

A gate that cannot fail proves nothing, so the corrected gate was checked
against its own failure mode. Re-classifying the same 16 rows with
`.SANCTIONED` emptied:

```
FAIL rows: 1  -> gate says: FAIL -- STOP
   family        rule       class verdict
 survival disjunction na_mismatch    FAIL
```

The gate does fail when the movement is not explicitly sanctioned. For
contrast, on the identical data the old summary line still reports

```
max(abs(diff), na.rm = TRUE) = 0.000000e+00  -> would have printed PASS
```

which is precisely the failure being corrected.
