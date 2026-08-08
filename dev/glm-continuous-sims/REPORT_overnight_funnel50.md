# REPORT: overnight — 50-replicate funnel characterization (O1–O3)

Deliverable of `~/Downloads/CC_TASK_md_overnight_funnel50.md`, **as amended**
by Larry's authorizing message, recorded verbatim in §0 below. Where that
message conflicts with the task document, the message wins.

---

## 0. The authorizing message, recorded verbatim

> O1–O3 definitions and go — from Larry. The overnight document is
> ~/Downloads/CC_TASK_md_overnight_funnel50.md; read it now. This message
> supersedes any earlier O1–O3 message; where it conflicts with that file,
> this message wins. Record it verbatim in REPORT_overnight_funnel50.md's
> header.
>
> 1. If the addendum with §C–E and the provenance record is already committed
> and pushed, proceed; otherwise finish that first. Provenance stands: the
> harness qmd as-is pre-fix and the 1707-byte pilot bundle committed, and both
> follow-up verification scripts committed with them.
>
> 2. O1, amended: all 50 pilot replicates under the MR-COMPATIBLE
> configuration — use_lasso = FALSE, use_grf = FALSE, use_dina = FALSE — all
> else pilot-verbatim including the seed table, each replicate executed under
> RNGkind("L'Ecuyer-CMRG") with set.seed(seed_for(i)) set explicitly per
> replicate (settled by D.2: CMRG reproduces the bundle's n_true = 341). MR
> stays OFF in O1. Per replicate, using the corrected q-column extraction: the
> funnel; the age/preanti grids and bracketing of age > 34 and
> preanti <= 744.5; the TRUE region Q (flag_harm == 1) row — oriented MD, SE,
> pass/fail vs the 30 floor; the nearest grid-analogue conjunction's row
> likewise; detected yes/no; the selected rule with its n and oriented MD.
> Tabulate across 50: n/50 where Q clears the floor; n/50 where the grid
> analogue clears; n/50 bracketing each boundary; first-emptying gate if any;
> selected-rule size and oriented-MD distributions; family size
> min/median/max. This turns D.3's one-replicate caveat into measured rates.
> Parallel is fine given the explicit per-replicate RNG kind and seed;
> PID-watch per the document.
>
> 3. Then O2 and O3 as the document defines them. Gates unchanged: no harness
> edits, no fix application — the flags, consistency_method, tryCatch/quiet,
> and RNG-kind items all stay stated-only — no re-pilot, no STOP A, no DGM
> changes, the two Unicode man pages untouched. Finish by committing and
> pushing REPORT_overnight_funnel50.md — transport, not approval — then stop
> for the night.

---

## Header — preconditions and commit hashes

| precondition (task document §Queueing) | status |
|---|---|
| `REPORT_blockers_1-3.md` committed and pushed | **yes** — `1697b5a5` |
| Tasks 1–3 all completed; Task 1 green and committed | **yes, with a recorded history**: Task 1 initially halted on a check disagreement (`0 errors \| 1 warning \| 2 notes`). Larry amended the acceptance criterion after review, and the closeout committed green at `ef504707`. §C of the addendum then established the disagreement was a check-surface difference, not a regression. |
| Task 2 and Task 3 sections present with concrete numbers | **yes** — `REPORT_blockers_1-3.md` §2, §3, plus `REPORT_blockers_1-3_ADDENDUM.md` §D.1–D.4 |
| Addendum §C–E and provenance committed and pushed (message item 1) | **yes** — addendum `7d2ef0ed`, provenance `169847d5`; both follow-up verification scripts committed in `7d2ef0ed` |

Commit chain on `feature/mr-in-replicates`:

```
1697b5a5  REPORT_blockers_1-3.md
b86664b5  diagnose_md_harm_pilot_zero_detection.R (as written, unexecuted)
ef504707  closeout ONE-COMMIT (five files)
169847d5  provenance: harness qmd as-is pre-fix + 1707-byte pilot bundle
7d2ef0ed  REPORT_blockers_1-3_ADDENDUM.md + both verification scripts
```

Run provenance for O1: `pkg_commit 7d2ef0ed`, `forestsearch 0.2.0`, R 4.6.1,
25 workers, PID-watched by recorded PID (`kill -0 $PID`), never by a name
string the watcher's own command line contains.

---

## O1 — the funnel on all 50 pilot replicates

### O1.0 Configuration actually executed

```
use_lasso = FALSE   use_grf = FALSE   use_dina = FALSE   mr_inference = FALSE
everything else pilot-verbatim, including the pre-generated seed table
RNGkind("L'Ecuyer-CMRG") + explicit set.seed(seed_for(i)) per replicate
n = 1000, 50 replicates, seed_base = 8316951
```

Script: `dev/glm-continuous-sims/verification/overnight_funnel50_O1.R`.

**Amendment rationale, on the record.** The task document's O1 specifies "the
SAME configuration the pilot ran, verbatim". Under that configuration the MR
guard fires in ~2 ms on every replicate (addendum §D.1) and there is no funnel
to characterize. The amendment replaces it with the MR-compatible
configuration. This is a **measurement instrument, not an applied fix**: the
harness on disk is untouched and still carries `use_lasso <- TRUE` at `:125`
and no `consistency_method` argument.

**Fidelity confirmed:** replicate 1 gives `n_true = 341`, matching the pilot
bundle exactly.

**Cost:** 436.5 s wall for 50 replicates on 25 workers; per-replicate
47.8 / 154.8 / 247.7 s (min/median/max). **0 errored replicates.**

### O1.1 Per-replicate table (all 50)

`enum` = combinations enumerated; `floor` = surviving the effect floor;
`ageB`/`preB` = brackets that boundary; `ageNr`/`preNr` = nearest cut;
`Qmd`/`Qse` = TRUE region Q oriented MD and SE; `GAmd` = nearest grid-analogue
conjunction; `selN`/`selMD` = selected rule size and oriented MD. All MD values
are on the `adverse_outcome = FALSE` oriented scale where positive means harm,
compared against the floor of **30**.

```
id nTrue enum floor ageB preB ageNr preNr Qmd Qse Qok GAmd GAok det selN selMD
1 341 1485 609 Y Y 34 777  25.20 13.68 REJ   23.74 REJ  Y 85  105.76
2 322 1596 864 Y n 34 716  68.85 14.72 PASS  66.55 PASS Y 72  113.67
3 335 1485 229 Y n 34 690  26.97 14.42 REJ   26.92 REJ  Y 105   77.36
4 352 1378 955 Y Y 35 768  57.09 13.90 PASS  51.38 PASS Y 96  103.09
5 331 1485 777 Y n 34 659  49.82 14.98 PASS  48.40 PASS Y 138  108.96
6 355 1485 880 Y Y 34 802  56.41 13.96 PASS  56.01 PASS Y 103  109.02
7 347 1596 710 Y n 34 700  52.79 13.09 PASS  52.59 PASS Y 84  103.93
8 324 1485 640 Y Y 34 792  43.76 14.52 PASS  42.11 PASS Y 69  103.61
9 339 1596 919 Y Y 34 751  67.82 15.74 PASS  64.81 PASS Y 65  120.50
10 342 1596 983 Y n 34 715  26.22 13.45 REJ   24.41 REJ  Y 81   98.82
11 379 1596 975 Y n 35 735  44.77 13.62 PASS  43.44 PASS Y 66  120.37
12 349 1596 992 Y Y 34 762  59.87 14.52 PASS  59.28 PASS Y 63  101.38
13 382 1275 331 Y n 35 676  33.79 13.51 PASS  37.37 PASS Y 62   87.99
14 306 1485 224 Y Y 33 753  27.96 15.78 REJ   24.40 REJ  Y 111   84.72
15 350 1596 373 Y n 34 733  23.25 14.25 REJ   21.72 REJ  Y 85   86.90
16 360 1485 629 Y Y 35 755  49.28 13.67 PASS  49.95 PASS Y 84  102.42
17 358 1485 580 Y Y 35 754  15.75 14.54 REJ   17.25 REJ  Y 86  105.55
18 359 1596 582 Y n 34 699  16.30 13.51 REJ   18.79 REJ  Y 67  111.55
19 304 1596 1033 Y Y 33 756  34.36 15.18 PASS  36.36 PASS Y 63  105.74
20 336 1596 357 Y Y 34 764  40.38 13.73 PASS  38.20 PASS Y 68   85.05
21 341 1485 1046 Y Y 34 761  61.61 13.19 PASS  62.83 PASS Y 93  114.16
22 361 1596 1070 Y n 34 689  41.65 13.67 PASS  40.19 PASS Y 147   88.36
23 356 1596 309 Y n 34 740  18.56 14.04 REJ   18.56 REJ  Y 93   74.51
24 346 1485 1020 Y n 34 740  65.14 13.58 PASS  65.14 PASS Y 76  133.70
25 315 1596 691 Y Y 34 760  47.13 15.00 PASS  46.67 PASS Y 111   66.44
26 339 1485 851 Y n 34 722  58.51 14.25 PASS  57.18 PASS Y 72  104.83
27 327 1378 164 Y n 34 740  15.18 14.10 REJ   15.18 REJ  Y 88   88.93
28 332 1275 403 Y n 34 722  29.10 14.05 REJ   30.32 PASS Y 62   93.77
29 340 1485 479 Y Y 34 745  20.65 13.77 REJ   18.89 REJ  Y 63   98.71
30 341 1596 984 Y n 34 736  62.29 13.56 PASS  61.97 PASS Y 103  100.39
31 355 1378 857 Y n 34 740  47.73 13.98 PASS  47.73 PASS Y 75   98.30
32 322 1485 642 Y n 34 722  56.44 14.91 PASS  56.21 PASS Y 61  103.99
33 346 1596 964 Y n 34 719  37.51 14.10 PASS  36.16 PASS Y 106   89.15
34 327 1378 923 Y n 34 717  48.48 14.53 PASS  45.69 PASS Y 72  110.09
35 353 1378 324 Y n 34 736  31.51 14.06 PASS  31.28 PASS Y 62   97.73
36 346 1596 249 Y n 34 722  23.52 14.21 REJ   23.88 REJ  Y 70   87.69
37 350 1596 623 Y n 34 719  31.89 12.98 PASS  33.24 PASS Y 106  103.80
38 368 1378 1002 Y n 34 678  30.97 13.09 PASS  30.94 PASS Y 94  111.72
39 339 1596 578 Y Y 35 787  41.55 13.76 PASS  37.42 PASS Y 67   84.54
40 342 1485 465 Y Y 34 751  25.26 14.13 REJ   26.93 REJ  Y 71  101.62
41 342 1596 1003 Y n 34 718  40.03 14.22 PASS  38.86 PASS Y 106  103.67
42 351 1596 1115 Y n 34 731  61.66 13.25 PASS  61.45 PASS Y 70  120.05
43 359 1596 1018 Y Y 34 755  26.82 13.35 REJ   28.90 REJ  Y 68  102.38
44 337 1596 605 Y Y 34 772  53.15 14.56 PASS  53.66 PASS Y 61  105.61
45 346 1596 1032 Y n 34 692  35.39 14.01 PASS  36.39 PASS Y 62  111.37
46 357 1485 421 Y n 34 709  30.17 14.60 PASS  25.99 REJ  Y 63  120.77
47 359 1485 471 Y n 35 744  50.99 13.74 PASS  46.26 PASS Y 71  106.63
48 349 1485 760 Y n 34 661  71.98 13.82 PASS  75.33 PASS Y 222   83.08
49 317 1596 753 Y Y 34 754  42.98 14.00 PASS  41.40 PASS Y 115  101.05
50 358 1485 931 Y Y 34 760  41.59 13.41 PASS  41.11 PASS Y 63   97.86
```

### O1.2 (c) Resolved threshold and orientation — constant across all 50

Verified, not assumed. Every field has exactly one distinct value over 50
replicates:

```
  thr_screening    : 30             (n distinct = 1)
  thr_consistency  : 10             (n distinct = 1)
  thr_scale        : identity       (n distinct = 1)
  adverse_exec     : FALSE          (n distinct = 1)
  admission_floor  : 30             (n distinct = 1)
  family_status    : no-front-end   (n distinct = 1)
  sg_focus         : hr             (n distinct = 1)
```

**No variation to report.** `scale: identity` confirms the MD path takes no log
on any replicate, and the value reaching the engines is `+30` throughout.

### O1.3 (a) Candidate family size

```
  cut labels    min / median / max :   25 /   27.0 /   28
  enumerated    min / median / max : 1275 / 1485.0 / 1596
```

Cuts are quantile-derived, so both vary per replicate — the direct
demonstration of the `family_status` roxygen's point (`R/forestsearch_helpers.R:1939-1945`)
that `"no-front-end"` is deliberately weaker than a manuscript §2.1 fixed
family: no front end fits a model, yet the cut locations still move with the
sample.

### O1.4 (b) The funnel, and the first-emptying gate

Across 50 replicates, at each gate (min / median / max, and how often the gate
leaves the qualifying set empty):

```
  n_evaluated      min  1275  median  1485.0  max  1596   zero-count 0/50
  n_variance       min  1223  median  1429.0  max  1536   zero-count 0/50
  n_prevalence     min  1223  median  1429.0  max  1536   zero-count 0/50
  n_redundancy     min  1191  median  1396.0  max  1501   zero-count 0/50
  n_events         min  1191  median  1396.0  max  1501   zero-count 0/50
  n_sample_size    min  1090  median  1298.0  max  1404   zero-count 0/50
  n_model          min  1090  median  1298.0  max  1404   zero-count 0/50
  n_effect_floor   min   164  median   731.5  max  1115   zero-count 0/50
```

**First-emptying gate: NONE. No gate empties the qualifying set on any of the
50 replicates.** The distribution the task document asks for is therefore
degenerate — there is nothing to distribute.

`detected` is **50/50**. Against the pilot's 0/50 under the identical seeds,
this is the decisive confirmation that the pilot's zero detection was the MR
guard and nothing else. The effect floor is the largest gate (median 1298 →
731.5, a ~44% cut) but never an emptying one.

### O1.5 (d) Bracketing of the true boundaries

```
  brackets age > 34         : 50/50
  brackets preanti <= 744.5 : 20/50
  brackets BOTH             : 20/50

  nearest-cut gap, age     (|nearest - 34|)    min/median/max : 0.00 /  0.00 /  1.00
  nearest-cut gap, preanti (|nearest - 744.5|) min/median/max : 0.50 / 22.50 / 85.50

  age     nearest-below range : [29, 34]      nearest-above range : [35, 39]
  preanti nearest-below range : [371, 744]    nearest-above       : NA in 30/50
```

**`age > 34` is bracketed on every replicate**, and the nearest age cut is
*exactly* 34 on the median replicate (gap 0.00 at both min and median).

**`preanti <= 744.5` is bracketed on only 20/50.** In the other 30 the grid
has no cut above 744.5 at all, so the boundary sits above the whole preanti
grid. The nearest cut is still close in absolute terms — median gap 22.5 on a
variable ranging to 2851 — but it is a cut *below* the true boundary, so the
grid analogue is a strict subset of the true region rather than a bracketing
pair.

This is a genuine, previously unmeasured asymmetry between the two boundaries.
It is recorded, not acted on: no `conf.cont` change is proposed here, and O1.6
shows it costs little in floor clearance.

### O1.6 (e) Truth-adjacent clearance of the effect floor

```
  TRUE region Q clears the 30 floor        : 36/50
  grid-analogue conjunction clears         : 36/50

  Q oriented MD  min/q1/med/q3/max : 15.18 / 28.25 / 41.57 / 53.06 / 71.98
  Q SE           min/q1/med/q3/max : 12.98 / 13.63 / 14.00 / 14.49 / 15.78
  Q n            min/median/max    : 304 / 346 / 382
  mean Q oriented MD = 41.4015     (population value 40)

  GA oriented MD min/q1/med/q3/max : 15.18 / 27.43 / 39.53 / 52.28 / 75.33
  GA n           min/median/max    : 315 / 341 / 367

  per-boundary single-cut analogues:
    !{age <= nearest}      clears : 37/50   median MD 38.76
    {preanti <= nearest}   clears : 34/50   median MD 34.14
```

**The mean of Q's oriented MD over 50 replicates is 41.40 against a population
value of 40** — the DGM calibration is sound and the oriented-MD measurement is
unbiased for it.

**Q clears the floor on 36/50 = 72% of replicates**, not the ~50% that addendum
§D.3's single replicate might have suggested. That caveat was explicitly
hedged as "one replicate cannot establish a rate"; the measured rate is 72%,
and replicate 1 (Q = 25.20, rejected) was an unlucky draw in the lower tail.
**This supersedes the pessimistic reading of §D.3, and the correction runs in
the reassuring direction.**

The grid-analogue conjunction clears at exactly the same rate, 36/50, and its
MD distribution tracks Q's closely (median 39.53 vs 41.57). **The cut grid is
not costing detection of the truth** — including on the 30 replicates where
`preanti <= 744.5` is unbracketed. §O1.5's asymmetry is real but not
consequential for floor clearance.

### O1.7 The selected rule — selection optimism, measured

```
  detected                                : 50/50
  selected n           min/q1/med/q3/max  : 61.00 / 66.25 /  73.50 /  95.50 / 222.00
  selected oriented MD min/q1/med/q3/max  : 66.44 / 90.30 / 102.75 / 108.38 / 133.70
  mean selected oriented MD               : 100.83
  selected MD exceeds Q's MD on           : 50/50 replicates
  mean(selected MD - Q MD)                : 59.42

  selected rule names age     : 12/50
  selected rule names preanti : 16/50
  selected rule names both    :  2/50
  distinct selected rules     : 50/50  (never the same rule twice)
```

This is the headline number of O1. **The selected region's in-trial oriented MD
exceeds the true region's on every single replicate, by a mean of 59.42 on a
scale whose true value is 40** — the selected effect is roughly 2.5× the truth.
Selected regions are small (median n = 73.5 against Q's 346) and never repeat:
50 replicates produce 50 distinct rules, and only 2/50 name both true
boundary variables.

That is winner's-curse optimism over ~1485 candidates at n = 1000, and it is
exactly the quantity multiplier resampling exists to correct — which makes
addendum §D.4 (MR never runs, because `consistency_method` resolves to
`"split"`) the more consequential of the two configuration defects. **Stated,
not applied**, as are all four items: the flags, `consistency_method`, the
`tryCatch`/`quiet` suppression, and the RNG-kind dependence.

Note these are *in-trial* effects, not bias against the exact β(Ĥ). The bias
and coverage readouts the harness exists to produce still require MR, and are
out of scope here.

---

## O2 — the harness's own checks, on record

Quoted from `quarto/simulations/actg175/continuous/mr_coverage_sweep_md_harm.qmd`
(28,563 B, sha256 `1083c0b8…df1098`, committed as-is pre-fix at `169847d5`).
Reading and quoting only; no edits.

The document carries **five** internal checks. Two are guards inside the
recorder/readouts, three are the named closed-form assertions.

### O2.1 Seed-table bounds guard — `:94–98`

```r
seed_for <- function(sim_id) {
  if (sim_id < 1L || sim_id > MAX_SIMS)
    stop("sim_id ", sim_id, " outside the pre-generated seed table.")
  SEED_TABLE[sim_id]
}
```

Asserts every requested global replicate id lies inside the pre-generated
table, so a seed is never silently derived from loop position. **Did not fire**
in the pilot (`sim_id` 1–50, `MAX_SIMS` 5000).

### O2.2 CHECK B — θ† against the DGM's own effects — chunk `check-B`, `:218–231`

```r
  ok_B <- identical(chk$harness, chk$dgm)
  cat(sprintf("\nCHECK B: identical() = %s   max |difference| = %.3e\n",
              ok_B, max(abs(chk$difference))))
  if (!ok_B) stop("CHECK B FAILED: theta-dagger does not reproduce the DGM's own effects.")
```

Asserts `fs_betaHhat_theta_dagger_check()` reproduces `dgm$hazard_ratios$harm_subgroup`
/ `$no_harm_subgroup` under `identical()` — an exact identity, not a tolerance.
**PASSED**, and the pilot bundle proves it: `truth$marg_H == truth$effect_Q ==
-40.000000000000` and `truth$marg_Hc == truth$effect_Qc == -26.255235876036`.
This ran before the replicates and is why the pilot proceeded at all.

### O2.3 PARTITION guard — `:429–437`

```r
  s <- r$nH_eval + r$nHc_eval
  ...
  if (!all(s == N)) stop("PARTITION FAILED — machinery defect, phase stops.")
```

Asserts `nH_eval + nHc_eval == N` on every row including undetected ones.
**Would have PASSED**: all 50 bundle rows carry `0 + 5000 = 5000 = N`, the ITT
complement of `ab53239` behaving as designed. It is not reached, because
check-A halts the render first.

### O2.4 CHECK A — the two-valued CATE identity — chunk `check-A`, `:534–565`

Asserts β(R) = δ + β_inter · P(Q | R) on the harness's own realized rules, with
all three ingredients derived rather than typed:

```r
  mx <- max(abs(A$diff), na.rm = TRUE)
  cat(sprintf("\nCHECK A: distinct rules = %d;  max |exact - closed form| = %.3e\n",
              nrow(A), mx))
  cat(sprintf("  exactly 0 on %d of %d rules\n", sum(A$diff == 0), nrow(A)))
  # The identity is exact in real arithmetic; in floating point the two orders
  # of summation differ in the last bits.  Anything above ~1e-10 is a defect.
  if (!(mx == 0 || mx < 1e-10))
    stop(sprintf("CHECK A FAILED: max |diff| = %.3e — machinery defect.", mx))
```

**This is the check that halted the pilot render — and it did NOT fail its
assertion.** It errored on an empty input, upstream of any comparison.

The chain, with line numbers:

- `:537` — `rules <- unique(r$sg_def[!is.na(r$sg_def) & nzchar(r$sg_def)])`.
  All 50 rows have `sg_def = NA`, so `rules` is `character(0)`.
- `:544` — `A <- do.call(rbind, lapply(rules, function(g) {...}))` over an
  empty vector returns **`NULL`**.
- `:553` — `print(A[, c("nH","P_Q","exact","closed_form","diff")], ...)`.
  `NULL[, cols]` is `NULL`, so this prints `NULL` **without error**.
- `:555` — `mx <- max(abs(A$diff), na.rm = TRUE)`. `A$diff` is `NULL` and
  `abs(NULL)` raises.

Reproduced verbatim:

```
Error in abs(A$diff) : non-numeric argument to mathematical function
```

**"The exact failing values" the task document asks for do not exist**: check-A
never computed a `diff`, an `exact`, or a `closed_form` for any rule, because
there were no rules. The failure is an empty realized-rule set — a downstream
symptom of the MR guard — not a numeric disagreement between the exact target
and its closed form. Neither the pilot bundle nor any render log carries a
failing value; no HTML and no render log survive (the results directory
contains only `fs_md_harm_n1000_res.rds`, confirming the handoff's "no grid,
no HTML").

### O2.5 CHECK C — naive bias vs mean selection effect — chunk `check-C`, `:577–590`

```r
  if (!(abs(nb - sm) < 1e-8))
    stop(sprintf("CHECK C FAILED: %.3e — the two computations disagree.", nb - sm))
```

Asserts the gate's `naive` estimate and the recorder's independently computed
in-trial MD give the same mean deviation from β(Ĥ), to `1e-8`. **Not reached**
(check-A halts first). It would also have had nothing to compare: it subsets to
`detected == 1L`, which is empty, so `nb` and `sm` would both be `NaN` and
`abs(NaN - NaN) < 1e-8` is `NA` — itself an error in `if`. **A second latent
empty-input fragility, of the same class as check-A's.** Recorded as an
observation; no edit proposed and none made.

### O2.6 Summary

| check | location | status in the pilot |
|---|---|---|
| seed-table bounds | `:94–98` | passed (not triggered) |
| CHECK B (θ† identity) | `:218–231` | **PASSED**, exactly |
| PARTITION | `:429–437` | would have passed (0 + 5000 = N on all 50) |
| CHECK A (CATE identity) | `:534–565` | **ERRORED at `:555` on an empty rule set** — not an assertion failure |
| CHECK C (bias vs selection effect) | `:577–590` | not reached; latent same-class fragility |

The two checks that ran both passed. The render died on emptiness, not on
disagreement — so no harness check ever contradicted the machinery.

---

## O3 — shim re-point inventory (read-only)

Inventory of engines still `source()`-ing shims for the eval-frame/theta
pieces. No changes, no preparation of changes.

### O3.1 The shims and what they define

| shim | eval-frame / theta functions it defines |
|---|---|
| `quarto/simulations/gbsg_redux/betaHhat_truth.R` | `build_eval_frame()` `:51`, `betaHhat_theta_dagger_check()` `:117` |
| `quarto/simulations/actg175/binary/betaHhat_truth_glm.R` | `build_eval_frame_glm()` `:30`, `betaHhat_theta_dagger_check_or()` `:81` |
| `quarto/simulations/actg175/continuous/betaHhat_truth_md.R` | `betaHhat_theta_dagger_check_md()` `:74` |

These are the shim counterparts of the two functions committed at `ef504707`:
`fs_build_eval_frame()` and `fs_betaHhat_theta_dagger_check()`.

### O3.2 Counts

```
  betaHhat_truth.R      (survival) : 105 engine files
  betaHhat_truth_glm.R  (binary)   :   9 engine files
  betaHhat_truth_md.R   (continuous):  0 engine files

  TOTAL distinct engine files      : 114
```

By directory:

```
  85  quarto/simulations/gbsg_redux
   9  quarto/simulations/actg175/binary
   6  quarto/simulations/gbsg_redux/legacy
   3  dev/sim-check/run_post
   3  dev/sim-check
   3  dev/identifier-alignment/rerun
   1  dev/sim-check/setupcheck
   1  dev/sim-check/run_pre
   1  dev/sim-check/fixcheck
   1  dev/review
   1  dev/identifier-alignment/sim_analyses
```

**The continuous shim has zero consumers.** `betaHhat_truth_md.R` exists on
disk but no engine sources it: the md/harm harness already calls the package
entry points `fs_build_eval_frame()` and `fs_betaHhat_theta_dagger_check()`
directly (qmd `:178`, `:186`). The continuous pathway is therefore **already
re-pointed**; the deferred work is the 114 survival and binary engines.

### O3.3 Full `file:line` inventory

```
dev/identifier-alignment/rerun/sim_dina_maxcons_fb_mr_m1_h10_knoise0_n500.qmd:97 : source("betaHhat_truth.R")
dev/identifier-alignment/rerun/sim_fs_maxcons_fb_mr_m1_h10_knoise0_n500.qmd:97 : source("betaHhat_truth.R")
dev/identifier-alignment/rerun/sim_grf_maxcons_fb_mr_m1_h10_knoise0_n500.qmd:97 : source("betaHhat_truth.R")
dev/identifier-alignment/sim_analyses/sim_fs_maxeffCons_fb_mr_m1_h10_knoise0_n500_combine_1_500.qmd:85 : source("betaHhat_truth.R")
dev/review/sim_fs_maxcons_fb_mr_m1_h10_knoise0_n500_batch_1_100.qmd:85 : source("betaHhat_truth.R")
dev/sim-check/diag_probe.R:10 : source("betaHhat_truth.R")
dev/sim-check/fixcheck/setup_chunk.R:19 : source("betaHhat_truth.R")
dev/sim-check/fs_t1_t2_m1_h10_knoise0_n500_batch_1_20.qmd:82 : source("betaHhat_truth.R")
dev/sim-check/run_post/attr_check.R:3 : source("betaHhat_truth.R")
dev/sim-check/run_post/diag_probe.R:10 : source("betaHhat_truth.R")
dev/sim-check/run_post/post.qmd:85 : source("betaHhat_truth.R")
dev/sim-check/run_pre/pre.qmd:82 : source("betaHhat_truth.R")
dev/sim-check/setupcheck/setup_chunk.R:19 : source("betaHhat_truth.R")
dev/sim-check/sim_fs_maxcons_fb_mr_m1_h10_knoise0_n500_batch_1_100.qmd:85 : source("betaHhat_truth.R")
quarto/simulations/actg175/binary/effMaxSG_mr_coverage_sweep_or075.qmd:73 : source("betaHhat_truth_glm.R")
quarto/simulations/actg175/binary/maxcons_mr_coverage_sweep_or075.qmd:73 : source("betaHhat_truth_glm.R")
quarto/simulations/actg175/binary/maxeffCons_mr_coverage_sweep_or075.qmd:73 : source("betaHhat_truth_glm.R")
quarto/simulations/actg175/binary/mr_coverage_sweep_or075.qmd:73 : source("betaHhat_truth_glm.R")
quarto/simulations/actg175/binary/mr_coverage_sweep_or075_s100.qmd:73 : source("betaHhat_truth_glm.R")
quarto/simulations/actg175/binary/mr_coverage_sweep_or075_s1k.qmd:73 : source("betaHhat_truth_glm.R")
quarto/simulations/actg175/binary/mr_coverage_sweep_or075_s500.qmd:73 : source("betaHhat_truth_glm.R")
quarto/simulations/actg175/binary/mr_coverage_sweep_or20.qmd:73 : source("betaHhat_truth_glm.R")
quarto/simulations/actg175/binary/mr_coverage_sweep_or35_s1k.qmd:73 : source("betaHhat_truth_glm.R")
quarto/simulations/gbsg_redux/dina_t1_t2_m1_h10_knoise0_n1000_batch_1_500.qmd:82 : source("betaHhat_truth.R")
quarto/simulations/gbsg_redux/dina_t1_t2_m1_h10_knoise0_n500_batch_1_20.qmd:82 : source("betaHhat_truth.R")
quarto/simulations/gbsg_redux/dina_t1_t2_m1_h10_knoise0_n500_batch_21_500.qmd:82 : source("betaHhat_truth.R")
quarto/simulations/gbsg_redux/dina_t1_t2_m1_h10_knoise0_n500_combine_1_500.qmd:82 : source("betaHhat_truth.R")
quarto/simulations/gbsg_redux/effMaxSG_mr_coverage_sweep_h10_knoise0.qmd:73 : source("betaHhat_truth.R")
quarto/simulations/gbsg_redux/effMinSG_mr_coverage_sweep_h10_knoise0.qmd:73 : source("betaHhat_truth.R")
quarto/simulations/gbsg_redux/fsNEW_t1_t2_m1_h10_knoise0_n500_batch_1_100.qmd:82 : source("betaHhat_truth.R")
quarto/simulations/gbsg_redux/fs_t1_t2_m1_h10_knoise0_n1000_batch_101_150.qmd:82 : source("betaHhat_truth.R")
quarto/simulations/gbsg_redux/fs_t1_t2_m1_h10_knoise0_n1000_batch_1_20.qmd:82 : source("betaHhat_truth.R")
quarto/simulations/gbsg_redux/fs_t1_t2_m1_h10_knoise0_n1000_batch_151_200.qmd:82 : source("betaHhat_truth.R")
quarto/simulations/gbsg_redux/fs_t1_t2_m1_h10_knoise0_n1000_batch_201_250.qmd:82 : source("betaHhat_truth.R")
quarto/simulations/gbsg_redux/fs_t1_t2_m1_h10_knoise0_n1000_batch_21_100.qmd:82 : source("betaHhat_truth.R")
quarto/simulations/gbsg_redux/fs_t1_t2_m1_h10_knoise0_n1000_batch_251_350.qmd:82 : source("betaHhat_truth.R")
quarto/simulations/gbsg_redux/fs_t1_t2_m1_h10_knoise0_n1000_batch_351_400.qmd:82 : source("betaHhat_truth.R")
quarto/simulations/gbsg_redux/fs_t1_t2_m1_h10_knoise0_n1000_batch_401_500.qmd:82 : source("betaHhat_truth.R")
quarto/simulations/gbsg_redux/fs_t1_t2_m1_h10_knoise0_n500_batch_101_200.qmd:82 : source("betaHhat_truth.R")
quarto/simulations/gbsg_redux/fs_t1_t2_m1_h10_knoise0_n500_batch_1_20.qmd:82 : source("betaHhat_truth.R")
quarto/simulations/gbsg_redux/fs_t1_t2_m1_h10_knoise0_n500_batch_201_500.qmd:82 : source("betaHhat_truth.R")
quarto/simulations/gbsg_redux/fs_t1_t2_m1_h10_knoise0_n500_batch_21_100.qmd:82 : source("betaHhat_truth.R")
quarto/simulations/gbsg_redux/fs_t1_t2_m1_h10_knoise0_n500_combine_1_500.qmd:82 : source("betaHhat_truth.R")
quarto/simulations/gbsg_redux/fs_t1_t2_m1_h15_knoise0_n1000_batch_101_150.qmd:82 : source("betaHhat_truth.R")
quarto/simulations/gbsg_redux/fs_t1_t2_m1_h15_knoise0_n1000_batch_1_100.qmd:82 : source("betaHhat_truth.R")
quarto/simulations/gbsg_redux/fs_t1_t2_m1_h15_knoise0_n1000_batch_1_20.qmd:82 : source("betaHhat_truth.R")
quarto/simulations/gbsg_redux/fs_t1_t2_m1_h15_knoise0_n1000_batch_151_200.qmd:82 : source("betaHhat_truth.R")
quarto/simulations/gbsg_redux/fs_t1_t2_m1_h15_knoise0_n1000_batch_201_250.qmd:82 : source("betaHhat_truth.R")
quarto/simulations/gbsg_redux/fs_t1_t2_m1_h15_knoise0_n1000_batch_251_300.qmd:82 : source("betaHhat_truth.R")
quarto/simulations/gbsg_redux/fs_t1_t2_m1_h15_knoise0_n1000_batch_301_500.qmd:82 : source("betaHhat_truth.R")
quarto/simulations/gbsg_redux/fs_t1_t2_m1_h15_knoise0_n1000_combine_1_500.qmd:82 : source("betaHhat_truth.R")
quarto/simulations/gbsg_redux/fs_t1_t2_m1_h15_knoise0_n500_batch_101_200.qmd:82 : source("betaHhat_truth.R")
quarto/simulations/gbsg_redux/fs_t1_t2_m1_h15_knoise0_n500_batch_1_100.qmd:82 : source("betaHhat_truth.R")
quarto/simulations/gbsg_redux/fs_t1_t2_m1_h15_knoise0_n500_batch_201_300.qmd:82 : source("betaHhat_truth.R")
quarto/simulations/gbsg_redux/fs_t1_t2_m1_h15_knoise0_n500_batch_301_400.qmd:82 : source("betaHhat_truth.R")
quarto/simulations/gbsg_redux/fs_t1_t2_m1_h15_knoise0_n500_batch_401_500.qmd:82 : source("betaHhat_truth.R")
quarto/simulations/gbsg_redux/fs_t1_t2_m1_h15_knoise3_n1500_batch_101_150.qmd:82 : source("betaHhat_truth.R")
quarto/simulations/gbsg_redux/fs_t1_t2_m1_h15_knoise3_n1500_batch_1_50.qmd:82 : source("betaHhat_truth.R")
quarto/simulations/gbsg_redux/fs_t1_t2_m1_h15_knoise3_n1500_batch_151_200.qmd:82 : source("betaHhat_truth.R")
quarto/simulations/gbsg_redux/fs_t1_t2_m1_h15_knoise3_n1500_batch_201_250.qmd:82 : source("betaHhat_truth.R")
quarto/simulations/gbsg_redux/fs_t1_t2_m1_h15_knoise3_n1500_batch_251_300.qmd:82 : source("betaHhat_truth.R")
quarto/simulations/gbsg_redux/fs_t1_t2_m1_h15_knoise3_n1500_batch_301_350.qmd:82 : source("betaHhat_truth.R")
quarto/simulations/gbsg_redux/fs_t1_t2_m1_h15_knoise3_n1500_batch_351_400.qmd:82 : source("betaHhat_truth.R")
quarto/simulations/gbsg_redux/fs_t1_t2_m1_h15_knoise3_n1500_batch_401_450.qmd:82 : source("betaHhat_truth.R")
quarto/simulations/gbsg_redux/fs_t1_t2_m1_h15_knoise3_n1500_batch_451_500.qmd:82 : source("betaHhat_truth.R")
quarto/simulations/gbsg_redux/fs_t1_t2_m1_h15_knoise3_n1500_batch_51_100.qmd:82 : source("betaHhat_truth.R")
quarto/simulations/gbsg_redux/fs_t1_t2_m1_h20_knoise3_n1500_batch_401_500.qmd:82 : source("betaHhat_truth.R")
quarto/simulations/gbsg_redux/fs_t1_t2_m1_h20_knoise3_n1500_combine_101_150.qmd:82 : source("betaHhat_truth.R")
quarto/simulations/gbsg_redux/fs_t1_t2_m1_h20_knoise3_n1500_combine_1_100.qmd:82 : source("betaHhat_truth.R")
quarto/simulations/gbsg_redux/fs_t1_t2_m1_h20_knoise3_n1500_combine_151_200.qmd:82 : source("betaHhat_truth.R")
quarto/simulations/gbsg_redux/fs_t1_t2_m1_h20_knoise3_n1500_combine_201_250.qmd:82 : source("betaHhat_truth.R")
quarto/simulations/gbsg_redux/fs_t1_t2_m1_h20_knoise3_n1500_combine_251_350.qmd:82 : source("betaHhat_truth.R")
quarto/simulations/gbsg_redux/fs_t1_t2_m1_h20_knoise3_n1500_combine_351_400.qmd:82 : source("betaHhat_truth.R")
quarto/simulations/gbsg_redux/grf_t1_t2_m1_h10_knoise0_n500_batch_101_200.qmd:82 : source("betaHhat_truth.R")
quarto/simulations/gbsg_redux/grf_t1_t2_m1_h10_knoise0_n500_batch_1_20.qmd:82 : source("betaHhat_truth.R")
quarto/simulations/gbsg_redux/grf_t1_t2_m1_h10_knoise0_n500_batch_201_400.qmd:82 : source("betaHhat_truth.R")
quarto/simulations/gbsg_redux/grf_t1_t2_m1_h10_knoise0_n500_batch_21_100.qmd:82 : source("betaHhat_truth.R")
quarto/simulations/gbsg_redux/grf_t1_t2_m1_h10_knoise0_n500_batch_401_500.qmd:82 : source("betaHhat_truth.R")
quarto/simulations/gbsg_redux/grf_t1_t2_m1_h10_knoise0_n500_combine_1_500.qmd:82 : source("betaHhat_truth.R")
quarto/simulations/gbsg_redux/legacy/fs_t1_t2_m1_h10_knoise0_n1000_combine_1_400.qmd:82 : source("betaHhat_truth.R")
quarto/simulations/gbsg_redux/legacy/fs_t1_t2_m1_h10_knoise0_n1000_combine_1_500.qmd:82 : source("betaHhat_truth.R")
quarto/simulations/gbsg_redux/legacy/fs_t1_t2_m1_h15_knoise0_n500_combine_1_500.qmd:82 : source("betaHhat_truth.R")
quarto/simulations/gbsg_redux/legacy/fs_t1_t2_m1_h15_knoise3_n1500_combine_1_500.qmd:82 : source("betaHhat_truth.R")
quarto/simulations/gbsg_redux/legacy/fs_t1_t2_m1_h20_knoise3_n1500_combine_1_250.qmd:82 : source("betaHhat_truth.R")
quarto/simulations/gbsg_redux/legacy/fs_t1_t2_m1_h20_knoise3_n1500_combine_1_500.qmd:82 : source("betaHhat_truth.R")
quarto/simulations/gbsg_redux/maxcons_mr_coverage_sweep_h10_knoise0.qmd:73 : source("betaHhat_truth.R")
quarto/simulations/gbsg_redux/mr_coverage_sweep_h075.qmd:73 : source("betaHhat_truth.R")
quarto/simulations/gbsg_redux/mr_coverage_sweep_h10_knoise0.qmd:73 : source("betaHhat_truth.R")
quarto/simulations/gbsg_redux/mr_coverage_sweep_h10_knoise3.qmd:73 : source("betaHhat_truth.R")
quarto/simulations/gbsg_redux/mr_coverage_sweep_h10_knoise6.qmd:73 : source("betaHhat_truth.R")
quarto/simulations/gbsg_redux/mr_coverage_sweep_h10.qmd:73 : source("betaHhat_truth.R")
quarto/simulations/gbsg_redux/mr_coverage_sweep_h15.qmd:73 : source("betaHhat_truth.R")
quarto/simulations/gbsg_redux/mr_coverage_sweep_h20.qmd:73 : source("betaHhat_truth.R")
quarto/simulations/gbsg_redux/mr_coverage_sweep_h25.qmd:73 : source("betaHhat_truth.R")
quarto/simulations/gbsg_redux/sim_fs_effMaxSG_fb_mr_m1_h10_knoise0_n500_batch_1_20.qmd:85 : source("betaHhat_truth.R")
quarto/simulations/gbsg_redux/sim_fs_effMaxSG_fb_mr_m1_h10_knoise0_n500_batch_201_400.qmd:85 : source("betaHhat_truth.R")
quarto/simulations/gbsg_redux/sim_fs_effMaxSG_fb_mr_m1_h10_knoise0_n500_batch_21_200.qmd:85 : source("betaHhat_truth.R")
quarto/simulations/gbsg_redux/sim_fs_effMaxSG_fb_mr_m1_h10_knoise0_n500_batch_401_500.qmd:85 : source("betaHhat_truth.R")
quarto/simulations/gbsg_redux/sim_fs_effMaxSG_fb_mr_m1_h10_knoise0_n500_combine_1_500.qmd:85 : source("betaHhat_truth.R")
quarto/simulations/gbsg_redux/sim_fs_effMinSG_fb_mr_m1_h10_knoise0_n500_batch_101_200.qmd:85 : source("betaHhat_truth.R")
quarto/simulations/gbsg_redux/sim_fs_effMinSG_fb_mr_m1_h10_knoise0_n500_batch_1_100.qmd:85 : source("betaHhat_truth.R")
quarto/simulations/gbsg_redux/sim_fs_effMinSG_fb_mr_m1_h10_knoise0_n500_batch_1_200.qmd:85 : source("betaHhat_truth.R")
quarto/simulations/gbsg_redux/sim_fs_effMinSG_fb_mr_m1_h10_knoise0_n500_batch_201_400.qmd:85 : source("betaHhat_truth.R")
quarto/simulations/gbsg_redux/sim_fs_effMinSG_fb_mr_m1_h10_knoise0_n500_batch_401_500.qmd:85 : source("betaHhat_truth.R")
quarto/simulations/gbsg_redux/sim_fs_effMinSG_fb_mr_m1_h10_knoise0_n500_combine_1_500.qmd:85 : source("betaHhat_truth.R")
quarto/simulations/gbsg_redux/sim_fs_maxcons_fb_mr_m1_h10_knoise0_n500_batch_1_100.qmd:85 : source("betaHhat_truth.R")
quarto/simulations/gbsg_redux/sim_fs_maxcons_fb_mr_m1_h10_knoise0_n500_batch_1_20.qmd:85 : source("betaHhat_truth.R")
quarto/simulations/gbsg_redux/sim_fs_maxcons_fb_mr_m1_h10_knoise0_n500_batch_21_300.qmd:85 : source("betaHhat_truth.R")
quarto/simulations/gbsg_redux/sim_fs_maxcons_fb_mr_m1_h10_knoise0_n500_batch_301_500.qmd:85 : source("betaHhat_truth.R")
quarto/simulations/gbsg_redux/sim_fs_maxcons_fb_mr_m1_h10_knoise0_n500_combined_1_500.qmd:85 : source("betaHhat_truth.R")
quarto/simulations/gbsg_redux/sim_fs_maxeffCons_fb_mr_m1_h10_knoise0_n500_batch_1_50.qmd:85 : source("betaHhat_truth.R")
quarto/simulations/gbsg_redux/sim_fs_maxeffCons_fb_mr_m1_h10_knoise0_n500_batch_301_500.qmd:85 : source("betaHhat_truth.R")
quarto/simulations/gbsg_redux/sim_fs_maxeffCons_fb_mr_m1_h10_knoise0_n500_batch_51_300.qmd:85 : source("betaHhat_truth.R")
quarto/simulations/gbsg_redux/sim_fs_maxeffCons_fb_mr_m1_h10_knoise0_n500_combine_1_500.qmd:85 : source("betaHhat_truth.R")
```

---

## What changed in the picture, and what did not

**Superseded by measurement.** Addendum §D.3 reported, from one replicate, that
the true region Q scored 25.20 against the floor of 30 and was rejected, and
flagged — explicitly as a one-replicate observation that "cannot establish a
rate" — that rejecting the truth looked like an ordinary outcome. O1 measures
the rate: **Q clears the floor on 36/50 = 72%**, and the mean oriented MD over
50 replicates is **41.40 against a population value of 40**. Replicate 1 was a
lower-tail draw. The correction runs in the reassuring direction.

**Confirmed and sharpened.** The pilot's 0/50 was the MR guard: under identical
seeds with the front ends off, detection is **50/50** and **no funnel gate
empties on any replicate**. Selection optimism is now quantified: the selected
region's in-trial oriented MD exceeds Q's on **50/50** replicates by a mean of
**59.42**, with 50 distinct rules in 50 replicates and a median selected size
of 73.5 against Q's 346.

**New, previously unmeasured.** `preanti <= 744.5` is bracketed on only
**20/50** replicates — in the other 30 the quantile grid has no cut above the
boundary. It costs nothing measurable in floor clearance (the grid analogue
clears at 36/50, identical to Q), but it is a real asymmetry between the two
true boundaries and it is on the record now. Recorded, not acted on.

**Unchanged and still stated-only.** All four items remain proposals, none
applied: (1) `use_lasso`/`use_grf` at qmd `:125`; (2) `consistency_method =
"resample"`, without which MR never runs at all; (3) the `tryCatch`/`quiet`
suppression that turned two well-worded diagnostics into a silent null result;
(4) the RNG-kind dependence of replicate data. The harness on disk is byte-identical
to the committed pre-fix provenance copy.

## Gates honoured

No harness edits. No fix application. No re-pilot. No STOP A. No DGM changes.
No `R/` modifications. No new NOTE files. No design decisions. The two Unicode
man pages (`man/fpr_calibration.Rd`, `man/sg_tables.Rd`) untouched. Compute was
pilot-scale and PID-watched by recorded PID, never by `pgrep -f` on a string
the watcher's own command line contains.

Committing this report is transport, not approval.
