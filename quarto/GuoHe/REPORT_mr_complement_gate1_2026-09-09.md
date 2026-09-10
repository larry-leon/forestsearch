# REPORT — Gate 1a (T1 pilot), Mac Studio, 2026-09-09

Governing documents: v1 (`179ae409`) §4, as amended by v2 (`669c9ef5`) A4 and
**v4 (`fb65b4a3`) N1/N3/N4**. Driver: `quarto/GuoHe/mr_field_complement_vs_guohe_run.R`.

**VERDICT: WITHIN the pre-authorized envelope. Production authorized and proceeding.**

## 1. Provenance

```
git log -1 --oneline  ->  78887c82 v4 addendum appended to the Stage 0 v2 record: probe scored PASS under N1, plus the two N6 corrections
git status -sb        ->  ## feature/glm-extension...origin/feature/glm-extension [ahead 11]
```

No `git fetch`, `git pull` or `git push` at any point. Machine: Mac Studio,
14 physical cores, 36 GiB (`hw.physicalcpu` 14, `hw.memsize` 38654705664).
Run with `VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1` (Accelerate is not
fork-safe; `mclapply` forks).

## 2. Pilot command and result

```
Rscript quarto/GuoHe/mr_field_complement_vs_guohe_run.R --pilot --cores=10
```

Pilot is `t7_beta2_00` (the maximal-bias cell) at reps = 20, per the A4
scaffolding transplanted from `guohe_sec52_run.R:70, 81-84, 144, 208-225`.

```
[run ] t7_beta2_00  20 reps ...
[done] t7_beta2_00  20/20 reps in 0.1 min  naive_mm 1  cur_mm 20  disc_mm 0  sel_mm 0
       worst_abs 6.66e-16  frac_id 0.424  mr_na 0  fld_na 0  c_na 0
```

## 3. Projection against the Gate 1a envelope (N4)

| quantity | measured | envelope | verdict |
|---|---|---|---|
| mean per-replicate (serial, gate + field) | **1.67 s** | — | — |
| full study, 6 × 2000 replicates, single core | **5.6 core-h** | ≤ 40 core-h | **WITHIN** (7× margin) |
| projected wall-clock at 10 cores | **33.4 min** | ≤ 90 min Mac wall | **WITHIN** (2.7× margin) |

Both conditions hold, so production is authorized without a further decision.
Ten workers were chosen against 14 physical cores, leaving headroom; the job is
compute-bound and small in memory (per replicate the largest objects are
`Xi` 400 × 5000 and `Xo` 400 × 1000).

## 4. N1 / N3 gate status on the pilot

| gate | requirement | measured | verdict |
|---|---|---|---|
| **N1 discrete** | integer / selection / flag / seed columns `identical()`, every replicate | `disc_mm` = **0** | **PASS** |
| **N1 float** | `all.equal()` at 1e-8 | no column outside tolerance; worst absolute deviation **6.66e-16**; 42.4% bit-identical | **PASS** |
| **N3 selection** | selected cutpoint matches the stored row, every replicate | `sel_mm` = **0/20** | **PASS** |

Selection margin on the pilot: the recorded `sel_gap` (top-1 minus top-2
oriented score) is well clear of the deviation scale throughout — at the one
replicate discussed below it is 1.93e-02, some 3.5 × 10¹³ times the largest
float deviation measured there.

### The two legacy counters, and why neither is a STOP

`naive_mm` and `cur_mm` are the **pre-v4** flags: both are built on plain
`identical()` over columns that include floats, which is exactly the standard
N1 replaced. They are retained in the bundle for continuity with the committed
16-cell campaign and are **not** gates under v4 (N9 replaced "any `identical()`
failure" with the N1 standard).

- **`cur_mm` = 20/20** — `cur_ok` compares the recomputed MR row against the
  Linux-built `mr_vs_guohe_t7_beta2_00.rds` with plain `identical()`. Failing on
  every replicate is the expected cross-platform result and is precisely the
  "provenance measurement, not a gate" N1 describes.
- **`naive_mm` = 1/20**, at m = 19 — diagnosed rather than assumed:

| component | `identical()` | \|diff\| | `all.equal` 1e-8 |
|---|---|---|---|
| `naive_point` | FALSE | 5.55e-17 | TRUE |
| `naive_lower` | FALSE | 5.55e-17 | TRUE |
| `c_hat_naive` | TRUE | 0 | TRUE |
| `gamma_s_naive` | TRUE | 0 | TRUE |
| `naive_cover` | TRUE | 0 | TRUE |
| selected `c_hat` | TRUE | 0 | TRUE |
| `gamma_s` at the selection | TRUE | 0 | TRUE |

  The failure is confined to the two float components, at 5.55e-17. Every
  discrete component — including the selection itself and both truth lookups —
  is bit-identical. Under N1 this replicate passes.

## 5. Authorization

- Gate 1a envelope: **WITHIN** on both conditions.
- N1 discrete: **PASS**. N1 float: **PASS**. N3 selection: **PASS**.
- **Production launched:** all six cells, `--cores=10`. The driver enforces the
  N1 and N3 STOPs at each cell boundary and aborts the run on any breach, so a
  violation cannot be written to a bundle and passed over.

The pilot bundle `mr_field_complement_vs_guohe_t7_beta2_00_pilot.rds` is a
scratch artifact, not part of the commit plan; it is removed after this record
is written and is reproducible from the command in §2.
