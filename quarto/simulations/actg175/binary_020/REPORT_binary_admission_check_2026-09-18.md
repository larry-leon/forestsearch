# REPORT — The minimum-events admission criterion on the binary path: read-only check

Date: 2026-09-18 (UTC). Machine: `pop-os` (64 physical cores, 251 GB; R 4.6.1, forestsearch 0.3.5).
Branch `feature/glm-extension`. Task: `dev/tasks/TASK_binary_admission_check_2026-09-18.md` (`633e15b4`).

**Read-only apart from §2's settlement.** No `R/` change, no install, no campaign, no relaunch, no
`forestsearch()` call, no MR. The only render is §2.1's 20-replicate assertion check.

---

## 1. Provenance — GATE PASS

```
pop-os
feature/glm-extension
98d1e3bb IS an ancestor of HEAD          (the halt commit)
[R / Rscript / quarto processes: none]
```

`<§1 HEAD>` = `ab89f11a`. `git status --porcelain` carried **no tracked modifications**; the untracked
entries were the Stage 1 smoke and calibration artefacts and the three continuous `logs_*`
directories:

```
?? binary_020/assert_check.html                 ?? binary_020/cal_orcal16.html
?? binary_020/cal_orcal32.html                  ?? binary_020/cal_orcal32b.html
?? binary_020/cal_orcal63.html                  ?? binary_020/cal_orcal63n500.html
?? binary_020/logs_or/                          ?? binary_020/smoke_{a_recipe,consistency,dina,grf}.html
?? binary_020/mr_or_harm/{10 smoke, assert and calibration bundle directories}
?? continuous/logs_{mddina,mdgrf,mdsgnb20}/
```

```
ab89f11a actg175 or: heartbeat records the second aborted launch; Stage 2 held at the user's instruction
adafc72c actg175 or: one non-convergence convention for every estimator; Gate 2 record gains its header and the halt section
6f80f292 actg175 or: clear HALT_or.md -- the cell-3 halt is diagnosed and fixed (1f99dfe4), Stage 2 resumes
1f99dfe4 actg175 or: separation in a true-region arm made a legitimate exp() underflow fatal -- guard corrected, recipe untouched
98d1e3bb actg175 or orfs_or150_n500: HALT -- render batch_1001_2000 failed (rc=1; 124 = timeout 9915s)
50c25059 actg175 or orfs_or075_n2000: consistency, target_or_h 0.75, n 2000, 2000 replicates + combine; Gate 2 PASS
```

## 2. The working tree

**The section assumes uncommitted changes; they were already committed.** All four files named in §2
match HEAD as of §1. `sim_fs_mr_field_or_template.qmd`, `scripts_or/gate2.R` and
`scripts_or/smoke_identity.R` were committed in `1f99dfe4`; `summary_actg175_or.qmd` in `adafc72c`.
This is recorded rather than worked around.

### 2.1 The invariant split — GATE PASS, and the commit already exists

The check was run as specified, tag `admsmoke`, on the committed green cell `orfs_or075_n500`,
`sim_id` 1–20, with the patched template:

```
columns: 180 total, 6 excluded as *_secs
         (fb_secs, fit_mr_secs, fld_H_secs, fld_H_uniform_secs, fld_Hc_secs, id_secs)
columns compared: 174
column sets identical both ways: TRUE
NA patterns identical on every compared column: TRUE
max relative difference over all compared numeric columns: 0

GATE 2.1: PASS -- assertion-only, no recorded value changes
```

The split is: estimates strictly positive and fatal; a bound fatal only when **negative**; zero or
infinite bounds counted and reported. It changes no recorded value. The three files were already
committed together in `1f99dfe4`, whose message states the split and that it changes nothing
recorded, so §2.1's commit step was already satisfied; nothing new was staged. The smoke output is
untracked: `mr_or_harm/fs_effMaxSG_mr_field_or075_n500_nb20_admsmoke_d5000/` (one bundle),
`admsmoke.html`, `logs_or/admsmoke.log`.

### 2.2 The summary — `git checkout` executed, and it was a no-op

`git checkout -- quarto/simulations/actg175/binary_020/summary_actg175_or.qmd` was run and exited 0.
The file already matched HEAD, so **the revert was a no-op**.

**The section's premise no longer holds, and this is flagged rather than silently accepted.** §2.2
describes the file as masking non-convergent **oracle** fits only, and reverts it because "masking
one estimator only is not the convention that will be adopted". The file has not carried that
version since `adafc72c`: on 2026-09-18, at the user's explicit instruction, the oracle-only masking
was **replaced** by the general convention — a replicate is non-convergent *for an estimator* when
that estimator's own point estimate or either two-sided bound is non-finite or its estimate is ≤ 0;
such rows are excluded from **that** estimator's statistics only; the count is reported per cell ×
block × estimator as `n_nonconvergent_fits` and as a `non-convergent` column in every affected
table; no row is dropped silently. The data recipe is the study's verbatim; the convention is
consumer-side.

So the intent behind §2.2 — "not the convention that will be adopted" — is already met, and
reverting to HEAD cannot remove the convention because the convention *is* HEAD. **If the intent is
that the summary should carry no non-convergence handling at all pending Larry's decision, that is a
one-commit revert of the summary hunk of `adafc72c`, and it was deliberately not done here**: it
would destroy work explicitly requested one turn earlier, on the strength of a premise this record
shows to be out of date.

---

## 3. The criterion, from source

All quotations at HEAD (`R/` unchanged since `064fce91`).

### 3.1 Is there a minimum-events rule? Where, what value, per arm or pooled?

**Yes, on both paths, and it is per arm.** Both are Status 3 of
`evaluate_subgroup_combination()` in `R/subgroup_search.R`.

*Binary / GLM* — `R/subgroup_search.R:591-600`:

```r
    # Status 3: Per-arm filter
    #   Binary:     minimum EVENTS (Y=1) per arm — analogous to survival
    #   Continuous:  d0.min/d1.min not meaningful — skip, rely on n.min
    if (is_binary) {
      # Events per arm: sum of binary outcome within each arm of the subgroup
      d0_sg <- sum(yy[id.x == 1 & tt == 0])
      d1_sg <- sum(yy[id.x == 1 & tt == 1])
      if (d0_sg < d0.min || d1_sg < d1.min) {
        return(list(status = 3L, result = NULL))
      }
```

*Survival* — `R/subgroup_search.R:652-656`, through two helpers:

```r
  event_counts <- calculate_event_counts(dd, tt, id.x)
  if (!meets_event_criteria(event_counts, d0.min, d1.min)) {
    return(list(status = 3L, result = NULL))
  }
```

`calculate_event_counts()` (`R/subgroup_search.R:727-733`) and `meets_event_criteria()`
(`R/subgroup_search.R:738-740`):

```r
calculate_event_counts <- function(dd, tt, id.x) {
  list(d0 = sum(dd[id.x == 1 & tt == 0]),
       d1 = sum(dd[id.x == 1 & tt == 1]),
       total = sum(dd[id.x == 1]))
}
meets_event_criteria <- function(event_counts, d0.min, d1.min) {
  return(event_counts$d0 >= d0.min && event_counts$d1 >= d1.min)
}
```

Both count **events in each arm separately**, never pooled. *Continuous / count*: Status 3 is skipped
entirely (`R/subgroup_search.R:609`, "Continuous: skip Status 3 entirely — d0.min/d1.min do not apply").

**Value in this campaign:** `d0.min = d1.min = 10`, the package default
(`R/forestsearch_main.R:1266-1267`) and the study's literal, carried by the template at
`sim_fs_mr_field_or_template.qmd:280` (`n_min <- 60L; d0_min <- 10L; d1_min <- 10L`) and passed at
`:664`. **This is exactly Larry's criterion: at least 10 events in each arm.**

### 3.2 What does it govern?

**Candidate admission.** Status 3 returns `list(status = 3L, result = NULL)`, so the candidate never
reaches the model fit (Status 5), the effect screen (Status 6) or the returned family (Status 7). It
is not a guard on whether an interval is computed.

### 3.3 What `d0.min` / `d1.min` and `n.min` count

- `d0.min` / `d1.min`: **events per arm** — control and treatment respectively — on survival (deaths,
  `dd`) and binary (`Y = 1`, `yy`). Ignored for continuous and count. Documented at
  `R/subgroup_search.R:17-21` and `R/forestsearch_main.R:687-694`, and the code above matches the
  documentation.
- `n.min`: **total subjects** in the candidate subgroup, both arms, on every path —
  `nx <- sum(id.x); if (nx <= n.min) return(status 4L)` (`R/subgroup_search.R:611-615` on the GLM
  branch, `:658-662` on the survival branch). Note the comparison is `<=`, so a subgroup must
  **exceed** `n.min`, not merely reach it.

### 3.4 The study's helper

`.logit_or_ci()`, defined at `maxeffCons_mr_coverage_sweep_or075.qmd:363` with its guard at
`:364-367`, mirrored in the template at `sim_fs_mr_field_or_template.qmd:471` / `:473-475`:

```r
  if (!nrow(df) || length(unique(df[[treat_name]])) < 2L) return(na4)
  y <- df[[outcome_name]]
  if (sum(y == 1L) < 5L || sum(y == 0L) < 5L) return(na4)
```

This counts **5 events and 5 non-events pooled over both arms**. It is **none of the above**: it is
not candidate admission — it is the recorder's own guard on whether the **oracle** interval is
computed on the **true** region, and it returns an all-NA quadruple when it fails. The true region is
not a candidate and never passes through Status 3, so `d0.min` / `d1.min` never touch it.

### 3.5 The mirror

**The binary path does implement the survival path's rule — at least 10 events in each arm, at
candidate admission — but only on the FS (consistency) path, and it counts events only.** Three
qualifications, each from source:

1. **It does not reach GRF or DINA.** `d0.min` / `d1.min` do not appear in
   `R/forestsearch_helpers.R` at all, so they are never forwarded to those identifiers. DINA's
   collector filters on subjects and effect only — `if (n_S < n_min) next` and
   `if (mean_tau < m_diff) next` (`R/dina_subgroup.R:746-749`). GRF's enumeration filters on
   subjects only — `if (nS < n_min || nS > n - 1L) next` (`R/grf_subgroup_labels.R:269`, `:302`) —
   with `dmin` an effect floor, not a count. **So for the `orgrf` and `ordina` campaigns the
   criterion is not applied at any point.**
2. **It counts events, not non-events, so it does not prevent separation.** An arm with 17 events and
   0 non-events passes a floor of 10 events and is completely separated. That is precisely the halt
   case. Preventing separation requires a floor on `min(events, non-events)` per arm.
3. **It never applies to the oracle.** The oracle refits on the *true* region in the recorder, guarded
   only by §3.4's pooled 5/5.

---

## 4. What the current design produces

No `forestsearch()` calls. Data regenerated with the committed recipe: the study's pre-generated seed
table indexed by global `sim_id`, `RNGkind("L'Ecuyer-CMRG")` per replicate, and the DGM calibrated
under the session-default generator before any kind switch, as the template builds it. 200
replicates (`sim_id` 1–200) per cell.

*Method note, recorded because it invalidated a first pass.* `gen()` sets `L'Ecuyer-CMRG`; a second
or third `calibrate_glm_interaction()` call in the same session therefore runs under that kind and
returns a **different** super-population. The first pass calibrated the OR 1.00 and OR 1.50 designs
that way, and its §4.2 n = 2000 arm counts likewise. Both were discarded and recomputed with the kind
reset before every DGM build. The corrected tables are internally consistent in a way the first pass
was not: |H| and every **control**-arm count are identical across the three design points, which they
must be, because the design points differ only in the treated-arm interaction.

#### 4.1a — the true region H: size and per-arm counts

| design | n | size_mean | size_q05 | size_q95 | e0_mean | e0_min | e1_mean | e1_min | ne0_mean | ne0_min | ne1_mean | ne1_min |
| ---|---|---|---|---|---|---|---|---|---|---|---|--- |
| OR 0.75 |   500 | 48.51 |    39 | 59.05 | 13.03 |     2 | 11.76 |     3 | 10.84 |     3 | 12.88 |     5 |
| OR 0.75 |   750 | 72.78 | 60.95 | 86.05 | 19.55 |     9 | 17.25 |     7 | 16.68 |     6 | 19.31 |    11 |
| OR 0.75 |  1000 | 96.64 |    83 |   113 | 26.36 |    15 | 23.46 |    12 | 21.75 |    10 | 25.08 |    11 |
| OR 0.75 |  2000 |   193 |   173 |   214 | 52.42 |    37 | 46.02 |    31 | 43.67 |    28 |  50.9 |    34 |
| OR 1.00 |   500 | 48.51 |    39 | 59.05 | 13.03 |     2 | 13.43 |     4 | 10.84 |     3 | 11.21 |     4 |
| OR 1.00 |   750 | 72.78 | 60.95 | 86.05 | 19.55 |     9 | 20.19 |    11 | 16.68 |     6 | 16.36 |     7 |
| OR 1.00 |  1000 | 96.64 |    83 |   113 | 26.36 |    15 | 26.88 |    16 | 21.75 |    10 | 21.67 |     8 |
| OR 1.00 |  2000 |   193 |   173 |   214 | 52.42 |    37 | 52.73 |    31 | 43.67 |    28 | 44.19 |    27 |
| OR 1.50 |   500 | 48.51 |    39 | 59.05 | 13.03 |     2 |  15.7 |     7 | 10.84 |     3 | 8.945 |     2 |
| OR 1.50 |   750 | 72.78 | 60.95 | 86.05 | 19.55 |     9 | 23.48 |    12 | 16.68 |     6 | 13.07 |     5 |
| OR 1.50 |  1000 | 96.64 |    83 |   113 | 26.36 |    15 |  31.3 |    18 | 21.75 |    10 | 17.25 |     7 |
| OR 1.50 |  2000 |   193 |   173 |   214 | 52.42 |    37 | 62.27 |    40 | 43.67 |    28 | 34.65 |    16 |

#### 4.1b — share of replicates failing each floor (true region H)

| design | n | fail_nmin | fail_10subj | fail_10events | fail_pooled55 | separation |
| ---|---|---|---|---|---|--- |
| OR 0.75 |  500 | 0.96 |    0 | 0.445 |    0 |    0 |
| OR 0.75 |  750 | 0.05 |    0 | 0.03 |    0 |    0 |
| OR 0.75 | 1e+03 |    0 |    0 |    0 |    0 |    0 |
| OR 0.75 | 2e+03 |    0 |    0 |    0 |    0 |    0 |
| OR 1.00 |  500 | 0.96 |    0 | 0.33 |    0 |    0 |
| OR 1.00 |  750 | 0.05 |    0 | 0.005 |    0 |    0 |
| OR 1.00 | 1e+03 |    0 |    0 |    0 |    0 |    0 |
| OR 1.00 | 2e+03 |    0 |    0 |    0 |    0 |    0 |
| OR 1.50 |  500 | 0.96 |    0 | 0.245 |    0 |    0 |
| OR 1.50 |  750 | 0.05 |    0 | 0.005 |    0 |    0 |
| OR 1.50 | 1e+03 |    0 |    0 |    0 |    0 |    0 |
| OR 1.50 | 2e+03 |    0 |    0 |    0 |    0 |    0 |

#### 4.1c — the complement Hc: no floor binds anywhere

| design | n | size_mean | e0_min | e1_min | ne0_min | ne1_min | fail_nmin | fail_10events | separation |
| ---|---|---|---|---|---|---|---|---|--- |
| OR 0.75 |   500 | 451.5 |    74 |    53 |   103 |   120 |     0 |     0 |     0 |
| OR 0.75 |   750 | 677.2 |   123 |    90 |   159 |   193 |     0 |     0 |     0 |
| OR 0.75 |  1000 | 903.4 |   158 |   122 |   219 |   258 |     0 |     0 |     0 |
| OR 0.75 |  2000 |  1807 |   353 |   260 |   460 |   545 |     0 |     0 |     0 |
| OR 1.00 |   500 | 451.5 |    74 |    53 |   103 |   120 |     0 |     0 |     0 |
| OR 1.00 |   750 | 677.2 |   123 |    90 |   159 |   193 |     0 |     0 |     0 |
| OR 1.00 |  1000 | 903.4 |   158 |   122 |   219 |   258 |     0 |     0 |     0 |
| OR 1.00 |  2000 |  1807 |   353 |   260 |   460 |   545 |     0 |     0 |     0 |
| OR 1.50 |   500 | 451.5 |    74 |    53 |   103 |   120 |     0 |     0 |     0 |
| OR 1.50 |   750 | 677.2 |   123 |    90 |   159 |   193 |     0 |     0 |     0 |
| OR 1.50 |  1000 | 903.4 |   158 |   122 |   219 |   258 |     0 |     0 |     0 |
| OR 1.50 |  2000 |  1807 |   353 |   260 |   460 |   545 |     0 |     0 |     0 |

#### 4.2 — declared subgroups

| cell | block | n_eval | size_mean | size_min | e0_mean | e0_min | e1_mean | e1_min | ne0_mean | ne0_min | ne1_mean | ne1_min | share_fail_10ev | share_separated |
| ---|---|---|---|---|---|---|---|---|---|---|---|---|---|--- |
| orfs_or075_n500 | Hhat |   200 |   103 |    61 | 16.64 |    10 | 24.51 |    12 | 35.31 |    15 | 26.56 |     6 |     0 |     0 |
| orfs_or075_n500 | Hhat^c |   200 |   397 |   230 | 92.55 |    61 | 64.94 |    36 | 105.3 |    45 | 134.2 |    88 |     0 |     0 |
| orfs_or075_n2000 | Hhat |   200 | 141.5 |    61 | 25.52 |    10 |  33.6 |    13 | 45.83 |    12 | 36.54 |     5 |     0 |     0 |
| orfs_or075_n2000 | Hhat^c |   200 |  1859 |  1418 | 418.1 |   327 |   317 |   248 | 510.1 |   374 | 613.4 |   467 |     0 |     0 |

**Reading 4.1.** Column `e0`/`e1` are events in the control and treated arms, `ne0`/`ne1` non-events.
`fail_nmin` is the share with |H| ≤ 60, the search's own rejection test. `fail_10events` is the share
with fewer than 10 events in **either** arm — Larry's criterion. `separation` is the share with a
zero cell in either arm.

**Reading 4.2.** `share_fail_10ev` is the share of **declared** Ĥ with fewer than 10 events in either
arm. Membership was recomputed from the recorded `sg_def` with
`forestsearch:::.fs_resolve_membership()` and reproduced the recorded `n_sel` (= `n_harm`) **exactly
on all 400 rows**, so the arm counts are the counts the run actually had.

**The halted cell has no committed rows.** `orfs_or150_n500` never completed a Gate 2; the halt
commit `98d1e3bb` added only `HALT_or.md` and the heartbeat, and the partial batch-1 bundle was
discarded by the runner's own restart cleanup. `git log --all` over that bundle path is empty. So
§4.2 covers the two committed green cells; the halted cell's declared subgroups cannot be checked
without re-running it, which is out of scope here.

### 4.3 Reach

- **Continuous campaigns** (`mdsgnb20`, `mdgrf`, `mddina`): the question does not arise — Status 3 is
  skipped entirely for continuous and count outcomes (`R/subgroup_search.R:609`), so there is no
  per-arm events criterion on that path and only `n.min` binds.
- **Survival campaigns**: the criterion is the survival path's own and is already applied at
  candidate admission through `meets_event_criteria()`, so committed survival campaigns enforce
  exactly this rule at whatever `d0.min` / `d1.min` they were run with.

---

## 5. Facts for Larry's decisions

**1. Does the binary path implement at least 10 events in each arm, and if not what instead?**

**Yes — on the FS (consistency) path, at candidate admission, at exactly 10 and 10 in this campaign.**
`R/subgroup_search.R:591-600` counts events per arm and rejects the candidate outright. Three
qualifications: it is **not applied on the GRF or DINA paths** (`d0.min`/`d1.min` are never forwarded
to them; they filter on subjects and effect only), it counts **events only and so does not prevent
separation** from an all-event arm, and it **never applies to the oracle**, which refits on the true
region under the study's pooled 5-events/5-non-events guard.

**2. Does any declared subgroup in the committed cells violate the criterion?**

**No.** Over 200 declared replicates in each of the two committed green cells — 400 in total — the
minimum events in either arm of Ĥ was **10** at n = 500 and **10** at n = 2000, and the share failing
was **0.000** in both. The floor binds exactly where it should: the smallest observed count is the
floor itself. Ĥᶜ is far from any floor (minimum 36 events in an arm). No declared subgroup showed
separation.

**3. How often does the true region fail each floor?**

| design | n | fails `n.min = 60` | fails 10 events per arm | fails 10 subjects per arm | fails pooled 5/5 | separated |
|---|---|---|---|---|---|---|
| OR 0.75 | 500 | **96.0%** | **44.5%** | 0 | 0 | 0 |
| OR 1.00 | 500 | **96.0%** | **33.0%** | 0 | 0 | 0 |
| OR 1.50 | 500 | **96.0%** | **24.5%** | 0 | 0 | 0 |
| any | 750 | 5.0% | 0.5–3.0% | 0 | 0 | 0 |
| any | 1000 | 0 | 0 | 0 | 0 | 0 |
| any | 2000 | 0 | 0 | 0 | 0 | 0 |

The pooled 5/5 guard — the only one the oracle is subject to — **never fails anywhere**, which is why
a separated fit reaches the Wald interval. Separation was not observed in these 200-replicate samples
at any cell; the halt's replicate is `sim_id` 1179 at OR 1.50 / n = 500, i.e. of order 1 in 1000.

**4. Can n = 500 ever declare a subgroup as small as the planted region?**

**No.** The search requires `nx > n.min`, so at `n.min = 60` the smallest declarable subgroup is 61
subjects. The planted region at n = 500 averages **48.5** subjects with a 95th percentile of **59**,
and **96.0%** of replicates put it at or below 60. Observed: the smallest **declared** Ĥ across 400
declared replicates was **61** at n = 500 — the floor exactly. So at n = 500 the search is
structurally unable to recover a region the size of the planted one; what it declares is a larger
region that overlaps it, which is consistent with the sensitivity of 0.245 and PPV of 0.126 recorded
for that cell.

## 6. Findings

1. **The intended criterion is implemented on the FS path and is being met**: no declared subgroup in
   either committed cell violates it, and the observed minimum sits exactly at the floor.
2. **It is absent from the GRF and DINA paths.** `d0.min` / `d1.min` are never forwarded there. Two of
   the three campaigns in this study would run without it.
3. **An events-only floor cannot prevent separation.** The halt case passed a 10-event floor in the
   arm that caused it (17 events) and failed only on non-events (0). A floor on `min(events,
   non-events)` per arm is what would bar it.
4. **The floor that binds at n = 500 is `n.min`, not the events floor**: 96% of planted regions are
   too small to be declarable at all, so that cell measures recovery of a *larger* overlapping
   region, not of the planted one.
5. **The oracle is governed by neither floor** — only by the study's pooled 5/5 — which is why the
   diverged fit reached the recorder at all.
6. **A latent defect in my own method, recorded**: calibrating a DGM after a replicate has switched
   the RNG kind yields a different super-population. It invalidated a first pass of §4 and was caught
   only because the recomputed membership stopped matching the recorded `n_sel`. The `n_sel`
   cross-check is worth keeping in any future consumer-side script.
7. The Gate 2 record had **no header**: the runner's `[ -f ... ]` test runs inside the `>> record`
   redirection that has already created the file. The committed `mdgrf` and `mddina` Gate 2 records
   have the same gap.

## 7. Scope

This is a read-only check of an admission criterion against source and against regenerated data. It
decides nothing, changes no `R/` code, proposes no change to the study's frozen data recipe, and
relaunches nothing.
