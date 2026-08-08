# REPORT: MR mechanism diagnostic (A)

Deliverable A3 of `CC_TASK_mr_mechanism_and_fb_pilot.md`. Baseline `c489c574`.
Investigation only: no `R/` edits, no harness edits, no DGM changes. The
superset gap is **not** touched; its pre-registration is recorded in §5.

Authority: `REPORT_stopB_md_harm_grid.md` (the finding),
`REPORT_blockers_1-3.md` §2.b (the MR consumption trace).

## Conclusion class: **RE-SELECTION FIDELITY GAP**

Stated up front, each clause tied to its evidence line below.

MR's simulated selection event is not the selection event that executed, and
the divergence widens with n exactly as the residual does. The gap is **not**
in the gates — every executing gate is replayed faithfully, with zero failures
across 250 replicates and 1.25 million draws (§2.3). It is in the **ranking
statistic**, and its cause is identified in source and confirmed against the
engine's own recorded table (§2.4):

- the identifier rounds the consistency rate to **2 decimals**
  (`p.consistency <- round(rr$rate_closed, pconsistency.digits)`,
  `R/subgroup_consistency_helpers.R:1457`, default
  `pconsistency.digits = 2`, `R/subgroup_consistency_main.R:357`), which puts
  `Pcons` on an 11-point grid above the 0.90 threshold and ties **31%** of
  candidates at exactly 1.00;
- `sort_subgroups()` for `sg_focus = "hr"` orders `(-Pcons, -hr, K)`
  (`R/subgroup_consistency_helpers.R:575`), so within that large tie set the
  winner is the **effect argmax** — a small, extreme subgroup;
- MR's `maxcons` ranks on **unrounded, untied** `zcons`
  (`.fs_mr_select`, `R/fs_mr_inference.R:169`), never reaches the `-hr`
  tiebreak, and therefore selects on consistency significance — which favours
  **large** subgroups (smaller `sigma_D`).

Consequence: MR simulates winners 2.5x the realized size at n = 1000 and 10x at
n = 4000 (§2.1). Bigger subgroups have smaller perturbations, so the simulated
selection bias understates the real optimism, and understates it worse as n
grows.

This is *not* a finding of linearization inadequacy. A2 (§3) shows the
correction decaying at n^-0.60 — close to the n^-1/2 a first-order term should
have — i.e. the correction is behaving like a well-formed leading-order term.
It is aimed at the wrong selection event.

**Not claimed:** that closing this gap would remove the residual. That is a
design question and is not settled here.

---

## 1. Method

Two cells, harness-verbatim config (resample consistency, MR args as
committed), the same seed table indexed by global replicate id, `RNGkind`
pinned to `L'Ecuyer-CMRG`:

| cell | replicates | draws | wall | ok |
|---|---|---|---|---|
| n = 1000 | 200 | 5000 | 216.1 s | 200/200 |
| n = 4000 | 50 | 5000 | 294.5 s | 50/50 |

Script: `dev/glm-continuous-sims/verification/mr_mechanism_A1.R`. For each
replicate it runs `forestsearch()` harness-verbatim, rebuilds the cut matrix
from the fit's own `$confounders.evaluated` labels, re-enumerates MR's family
with the package's own combination helpers, and replays the assembly + draw
pipeline through `.fs_mr_assemble`, `.fs_mr_multipliers` and `.fs_mr_select`.

**The replay is distributionally faithful with its own recorded draw seed. It
does not reproduce the grid run's individual draws bit-for-bit, and no such
claim is made.**

---

## 2. A1 results

### 2.1 Simulated-winner size vs realized — the discriminating quantity

```
############ n = 1000, s = 200 ############
  realized winner size      min/q1/med/q3/max :   61.0   69.0   85.5  117.2  386.0
  simulated winner MEDIAN   min/q1/med/q3/max :   76.0  177.8  223.5  302.0  562.0
  median simulated / realized size ratio       : 2.46
  mean frac of draws with size <= 100 : 0.1386    <= 200 : 0.4394

############ n = 4000, s = 50 ############
  realized winner size      min/q1/med/q3/max :   63.0  123.2  144.0  187.5  309.0
  simulated winner MEDIAN   min/q1/med/q3/max :  732.5 1414.0 1690.0 2054.0 2779.0
  median simulated / realized size ratio       : 10.29
  mean frac of draws with size <= 100 : 0.0017    <= 200 : 0.0154
```

**The ratio grows 2.46 → 10.29 across the two cells**, tracking the residual's
growth. At n = 4000 only 1.5% of draws select a winner even as large as 200,
while the realized winner has median size 144.

### 2.2 Re-selection frequency

```
  n = 1000  reselect_freq  min/q1/med/q3/max : 0.0  0.0  0.0  0.1  0.6
  n = 4000  reselect_freq  min/q1/med/q3/max : 0.0  0.0  0.0  0.0  0.1
```

The realized winner is re-selected on a **median of 0%** of draws in both
cells. The selection event MR averages over is almost disjoint from the one
that occurred.

### 2.3 Gate replay — traced from source, then measured

From source, the gates that appear inside MR's selection path:

- **Effect floor and consistency screen: PRESENT.** Resolved once by
  `.fs_resolve_admission` and carried as `admission`, then applied as
  `t_g <- pmax(admission$effect_floor, c_cons + z * sdv)` with
  `.admit <- function(bs) which(bs >= t_g)` (`R/fs_mr_inference.R:417-420`).
  The comment there is explicit that the domain "is NOT reconstructed from raw
  thresholds at this site, because reconstruction is what let MR's admission
  set drift from the identifier's."
- **`n.min`: PRESENT**, applied at family construction
  (`if (length(mem) >= n.min)`, `R/forestsearch_main.R:3195-3196`).
- **`d0.min` / `d1.min` and `max_subgroups_search`: ABSENT.** Recorded in
  source as a known gap (`R/forestsearch_main.R:3176-3181`).
- **No size gate inside `.fs_mr_select`** — the rules are pure argmax over
  `passers` (`R/fs_mr_inference.R:168-176`).

Measured fraction of draws whose simulated winner would FAIL each real gate:

```
                                       n = 1000        n = 4000
  per-arm minima d0.min/d1.min = 12    mean 0.00000    mean 0.00000   (max 0)
  effect floor 30                      mean 0.00000    mean 0.00000   (max 0)
  consistency screen                   mean 0.00000    mean 0.00000   (max 0)
```

**Zero gate failures anywhere**, across 250 replicates x 5000 draws. The two
gates MR does not replay are not being violated in practice: simulated winners
satisfy the per-arm minima regardless. **The fidelity gap is not a gate-replay
gap.**

### 2.4 The ranking statistic — where the gap actually is

Zero-perturbation probe: does MR's ranking statistic, evaluated on
**unperturbed** data, pick the candidate the identifier actually selected? It
must, if MR is linearizing the right selection event.

```
                                                       n = 1000      n = 4000
  obs_match (MR argmax == realized winner)             30/200 = 0.150   0/50 = 0.000
  rank of realized winner under MR's statistic
    min/q1/med/q3/max                        1.0 3.0 11.5 42.8 637.0   23 121 270.5 428 822
  admitted candidates  median                          477.5            794.0
```

At n = 4000 MR's statistic **never once** picks the realized winner, and ranks
it 270th of 794 on the median replicate.

Two candidate explanations were tested and separated:

**(i) Superset — partially responsible, not the driver.** MR's family is about
twice the identifier's qualifying set (median 1334 vs 666 at n = 1000). MR's
argmax lies inside the identifier's qualifying set on 70% (n = 1000) and 55%
(n = 4000) of replicates; the realized winner lies in MR's family on **100%**
of replicates in both cells. Restricting MR's argmax to the identifier's own
qualifying set barely helps:

```
  argmax-z restricted to the qualifying set == realized winner
      n = 1000 : 5/40  = 0.125        n = 4000 : 0/20 = 0.000
  restricted argmax size (median)  334  vs realized 86   |  2420 vs 138
```

So even on a common candidate set the two rules disagree almost always. The
superset is a contributor, not the mechanism.

**(ii) Two-stage early stopping — ruled out by measurement.**

```
  early_stop_triggered : 0/40 (n = 1000)   0/20 (n = 4000)
  candidates evaluated : median 608 of 608 total (fraction 1.000)
```

Every candidate is evaluated; the identifier is a genuine argmax over its
qualifying set, not a first-past-the-post rule.

**(iii) The rounding-induced tie set — this is the mechanism.** Read from the
engine's own recorded result table on a replicate where `obs_match` is FALSE
(n = 1000, `sim_id = 2`):

```
 K   N       hr Pcons
 2  72 113.6683     1
 2 124 110.2476     1
 2 113 110.0663     1
 2 105 109.4399     1
 2  75 108.9971     1
 2  88 106.1751     1

 Pcons: range 0.90000000 .. 1.00000000   ties at max: 198 of 634
 fraction of rows with Pcons == 1 exactly: 0.3123
 winner (row 1) : N = 72   hr = 113.668   Pcons = 1.00000000
 argmax-Pcons   : N = 72   hr = 113.668   Pcons = 1.00000000
 realized winner size : 72
```

198 of 634 candidates sit at `Pcons == 1.00` exactly. The sort key
`(-Pcons, -hr, K)` therefore falls through to `-hr` for a third of the field,
and the winner is the **effect argmax among the tied set** — `N = 72`, the
smallest and most extreme.

The rounding is explicit in source:

```r
# R/subgroup_consistency_helpers.R:1457 (and :1485, :1840)
p.consistency <- round(rr$rate_closed, pconsistency.digits)
```

with `pconsistency.digits = 2` (`R/subgroup_consistency_main.R:357`). Above the
0.90 threshold that leaves only 11 attainable values, so ties are pervasive by
construction, not by coincidence.

**A correction to an intermediate step of this investigation, recorded because
it changes what one would conclude.** An earlier check recomputed `Pcons` as
`2*pnorm(zcons) - 1` from MR's own `sigma_D` and found *no* saturation
(max 0.999780735718094, zero ties), which appeared to refute the tie-set
hypothesis. That check was measuring the wrong quantity: the engine stores the
**rounded** rate, and `round(0.999780735718094, 2) == 1`. Reading the engine's
recorded table rather than recomputing it is what settled the question.

MR, for its part, applies no rounding: `.zcons` is
`(bs - c_cons)/sdv` (`R/fs_mr_inference.R:442`) and `maxcons` is
`passers[which.max(zcons[passers])]` (`R/fs_mr_inference.R:169`) — continuous,
untied, and therefore never reaching an effect tiebreak.

### 2.5 Correction vs naive bias

```
                        n = 1000        n = 4000
  mean naive bias        66.3028         53.9172
  mean correction        32.1441         11.8236
    selection component  32.1350         11.8383
    fixed component       0.0091         -0.0147
  mean MR residual       34.1587         42.0936
  family size (median)     1342            1527
```

The correction is almost entirely the selection-bias term; the fixed-bias term
is ~0 at both n. The correction formula, quoted rather than assumed
(`R/fs_mr_inference.R`):

```r
:456    if (!is.na(s)) { sel_bias[b] <- P[s, b]; winner[b] <- s }
:481    selection_bias <- mean(sel_bias, na.rm = TRUE)   # over `ok_H`
:482    fixed_bias     <- mean(P[sel, ok_H])             # over `ok_H`
:485    beta_deb   <- beta_naive - selection_bias - fb
```

---

## 3. A2 — scaling decomposition (read-only, from the four committed bundles)

```
    n naive_bias correction mr_bias win_size se_wald   se_ij
  500    69.6892    39.8808 29.8085  89.9218 28.3173 34.9970
 1000    65.4906    32.2885 33.2021 108.1890 26.4876 31.0001
 2000    58.6018    21.4375 37.1643 137.6450 23.9693 27.0573
 4000    48.7735    11.3851 37.3884 195.0660 20.2589 21.9883
```

`win_size` is the mean realized-winner size (`n_harm`); `se_wald` is the winner
SE recovered from the recorded naive interval as
`(nv_H_hi - nv_H_lo) / (2 * qnorm(0.975))`; `se_ij` is the recorded
`t2_H_se_ij`.

Log-log slopes vs n:

```
  naive_bias  slope -0.1705   R^2 0.9495
  correction  slope -0.6017   R^2 0.9513
  mr_bias     slope +0.1143   R^2 0.9071
  win_size    slope +0.3699   R^2 0.9796
  se_wald     slope -0.1594   R^2 0.9586
  se_ij       slope -0.2208   R^2 0.9831
```

**This formalizes the STOP B reading with data.** The correction decays at
**n^-0.60**, close to (slightly steeper than) the n^-1/2 rate a first-order
term should have — the correction is well-formed. The quantity it must cancel,
the naive bias, decays at only **n^-0.17**. A term shrinking at n^-0.60 cannot
cancel one shrinking at n^-0.17, so the residual necessarily grows: measured
slope **+0.1143**.

The winner SE decays at n^-0.16 to n^-0.22, not n^-1/2, because the winner's
size grows only as n^0.37 — the selected region does not scale with the sample,
so its SE does not enjoy the root-n rate. That is the same reason the naive
bias decays slowly.

---

## 4. Evidence-to-conclusion map

| clause | evidence |
|---|---|
| Not a gate-replay gap | §2.3 — zero gate failures over 250 replicates x 5000 draws, all three gates |
| Not two-stage early stopping | §2.4(ii) — 0/60 replicates stopped early; all candidates evaluated |
| Not primarily the superset | §2.4(i) — restricted to the identifier's own qualifying set, agreement is still 0.125 / 0.000 |
| Is a ranking-statistic gap | §2.4 — obs_match 0.150 / 0.000; realized winner median rank 11.5 / 270.5 |
| Caused by rounding-induced ties + the `-hr` tiebreak | §2.4(iii) — 198/634 tied at `Pcons == 1.00`; `round(..., 2)` at `R/subgroup_consistency_helpers.R:1457`; sort key at `:575` |
| Correction is well-formed, aimed wrong | §3 — correction slope -0.60 vs naive-bias slope -0.17 |
| Gap widens with n, as the residual does | §2.1 — size ratio 2.46 → 10.29; §2.4 — obs_match 0.150 → 0.000 |

---

## 5. Superset gap — pre-registration, NOT acted on

`R/forestsearch_main.R:3176-3181` records that MR's reconstructed family does
not replay `d0.min`/`d1.min` or `max_subgroups_search`, so it is a superset of
the family the identifier chose among. Measured here: median 1334 vs 666 at
n = 1000 (§2.4(i)).

**Directional prediction, registered before any attempt to close it:** closing
the gap should make the residual **slightly LARGER**, not smaller. A superset
gives the simulated re-selection more candidates to win with, which inflates
the simulated optimism `selection_bias`; removing them lowers the correction,
and since the correction is subtracted, the residual rises.

The magnitude should be small, because §2.3 measured zero per-arm-minima
failures among simulated winners — the extra candidates are largely not the
ones winning. This is recorded so the prediction cannot be retrofitted.

---

## 6. Scope

Established only for: this DGM, `sg_focus = "hr"`, `consistency_method =
"resample"`, `focus = "harm"`, n in [500, 4000]. The tie-set mechanism depends
on `pconsistency.digits`, whose default is 2; nothing here measures behaviour
at other values, and no change to it is proposed.

Nothing here proposes a remedy. Whether MR's reselection rule should replay the
rounding and the `-hr` tiebreak, whether the rounding should exist at all, and
whether either change would reduce the residual, are design questions this
diagnostic does not settle.

Committing this report is transport, not approval.

---

## 7. B1 — full-bootstrap implementation check and pricing

Run because B1 was permitted in parallel with A. **B2 is NOT started**: the
price is above the stated gate.

### 7.1 Implementation check — the machinery DOES support this path

Verified from source and then confirmed by execution:

- **GLM outcomes are supported.** `forestsearch_bootstrap_dofuture()` branches
  on `is_glm <- outcome_type != "survival"` (`R/bootstrap_dofuture_main.R:295`)
  and builds `estimator_fn_boot` via `make_effect_estimator()` (`:302-308`).
- **MD is correctly treated as identity-scale.**
  `est_loghr <- effect_measure %in% c("OR", "RR", "IRR")`
  (`R/bootstrap_dofuture_main.R:561`) excludes MD, so no log transform is
  applied.
- **`consistency_method = "resample"` threads through unfiltered.**
  `args_FS_template <- fs.est$args_call_all`
  (`R/bootstrap_analysis_dofuture.R:406`) with no argument filtering, and each
  replicate re-runs `do.call(forestsearch, args_FS_boot)` (`:614`). Confirmed
  by execution: the run printed
  `consistency_method in args_call_all: resample`.
- **MR inside replicates is stripped by default** (`mr_in_replicates = FALSE`,
  `R/bootstrap_dofuture_main.R:218`; `args_FS_boot$mr_inference <- FALSE` at
  `R/bootstrap_analysis_dofuture.R:607`), which is what a FB pilot wants.

One FB run completed successfully and returned the expected structure:

```
FB OK: nb_boots=500, wall 361.5 s
returned names: results, SG_CIs, FSsg_tab, Ystar_mat, H_estimates,
                Hc_estimates, nb_boots, original_sg, outcome_type,
                effect_measure, est.scale, mr_replicates
```

Nothing is missing. **The path is implemented.**

### 7.2 Pricing — B2 is above the gate

`nb_boots` has no formal default; the documented recommendation is 500-1000
with a warning below 100 (`R/bootstrap_dofuture_main.R:29, :258-259`). Priced
at the lower end, **500**.

Measured, ONE replicate at n = 1000:

```
  nb_boots        500
  wall clock      361.5 s   (100 workers, multisession)
  => ~36,150 core-seconds per replicate, ~72.3 core-seconds per bootstrap fit
```

B2 as specified is n = 1000, s = 100:

| accounting | estimate |
|---|---|
| replicates run sequentially, each on 100 workers | 100 x 361.5 s = **10.0 hours** |
| total work repacked across all 128 cores at perfect efficiency | 100 x 500 = 50,000 fits x 72.3 core-s = 3,615,000 core-s / 128 = **7.8 hours** |

Both accountings land far above the **~3 hour** gate the task sets, and the
second is a floor, not a forecast — it assumes perfect packing and ignores
worker startup and the 100 outer fits.

**Per the task's instruction, B2 is not started and the price is on record.**
Halving to `nb_boots = 250` would still price at ~4-5 hours and would sit below
the documented recommendation; reducing `s` would weaken the comparison against
the MR n = 1000 row it exists to make. Neither trade is taken here — the
decision is the reviewer's.

For reference, the row B2 would have been scored against
(`REPORT_stopB_md_harm_grid.md`, n = 1000, s = 1000): MR bias **33.20**,
coverage **0.949**, width **121.5**.

### 7.3 What B2 would still be worth doing for

Recorded so the price can be judged against the value, not in isolation. §2-§4
establish that MR's residual comes from re-selecting on a different statistic
than the identifier uses. The full bootstrap re-runs the *actual* identifier
inside each replicate, so it has no re-selection fidelity gap by construction.
It is therefore the direct test of whether the residual is attributable to that
gap: if FB centers on beta(Hhat) where MR does not, the diagnosis in §2.4 is
confirmed end-to-end.

That is a strong reason to run it eventually. It is not a reason to overrun a
stated gate without approval.
