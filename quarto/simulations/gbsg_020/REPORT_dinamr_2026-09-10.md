# REPORT — DINA MR campaign at the FS-analogous criterion, both prevalences (campaign `dinamr`)

**Date:** 2026-09-11. **Executor:** Claude Code. **Task:** `dev/tasks/TASK_dinamr_campaign_2026-09-10.md`
(committed as received, with the kickoff amendments appended).
**Host:** Mac-Studio-3 (Apple M4 Max, 14 physical cores, 36 GB), R 4.5.2, arm64, `forestsearch` 0.3.5
installed from HEAD. **Not** the Linux host the task document was written for — see @Gate 1.
**Predecessors:** `REPORT_dinamr_stage0_2026-09-10.md` (the STOP), `REPORT_o1_forwarding_2026-09-10.md`
(which cleared it), `REPORT_cert20_2026-09-08.md`, `REPORT_tier2_2026-09-08.md`,
`REPORT_p12ext_2026-09-09.md`.
**Report and wait.** No acceptance criteria, no recommendation, no certification language.

---

## Framing — read this before any number

**This does not certify DINA.** DINA's candidates are read off a cross-fit surface that a bootstrap
would regenerate, so the manuscript's fixed-family condition (§2.1) **does not hold**. After the
alignment repair the estimate and interval are first-order exact **conditional on the proposed
family**. **Every coverage number in this report and in `summary_dinamr.qmd` is coverage of that
conditional-on-proposed-family estimand**, and each table says so in its own caption.

**The criterion is FS-analogous by construction.** `effMaxSG` at ε = 0.20 reuses
`.compute_inclusion_band()` — the single band helper shared with `forestsearch()`, DINA, GRF and MR —
with DINA's effect floor **log(0.90) = −0.105361** and **no consistency term**. That missing
consistency term is the one structural difference from FS, because neither DINA nor GRF has a
consistency floor (`REPORT_o1_forwarding_2026-09-10` §V2, measured). DINA applies the band as a
**sort key**, so it can never empty.

**`dina_select_statistic = "effect"` — confirmed to be the template default**, and confirmed to be
what the *package* resolves to: `formals(forestsearch)$dina_select_statistic` is
`c("effect", "dina")`, so `match.arg()` gives `"effect"` with or without the template's pin.

**`FS_S7_ER_JCUTS` was deliberately NOT set.** It is **inert on DINA**: it feeds
`fs_conf.cont_jcuts`, which reaches the consistency-only `method_args`; DINA's candidate grid does
not read it. The `er_jcuts = 10` that appears in every `dinamr` bundle's meta is the template default
and governs nothing on this path. (`REPORT_o1_forwarding` §V5 records the same for GRF.)

**GRF is out of scope**, as are any `R/` change, any recommendation change, and the certification
language.

---

## Protocol actions

**Archive (first action).** Eight stale campaign documents moved from `~/Downloads` to
`~/Downloads/cc_archive/`: `claude_cc_task_guohe_supplement_2026-09-09.md` and its `_v2` / `_v4`
variants, `claude_proposal_guohe_supplement_2026-09-09_v2.md`,
`claude_OC_summary_guohe_t7_2026-09-09.md`, `REPORT_mr_field_vs_guohe_2026-09-05.md`,
`REPORT_mr_vs_guohe_2026-09-04.md`, and `TASK_tier2_mac_2026-09-08.md`. Every one is already
committed under `dev/tasks/`, which is why each was safe to move. Only Markdown campaign documents
were touched; the `.R`, `.sh`, archive and installer files in `~/Downloads` were left alone.

> **`HANDOFF_guohe_comparison_2026-09-09.md`, which the protocol protects, is not present in
> `~/Downloads` at all** (nor anywhere under it; the only `HANDOFF_*` files there are in
> `cc_archive/` and `manuscript_jasa/`). There was nothing to protect, and nothing resembling it was
> moved.

**Task document** copied to `dev/tasks/TASK_dinamr_campaign_2026-09-10.md` verbatim, with an appendix
recording the four kickoff amendments. Committed. **Commit only; not pushed.**

**No `R/` change.** Nothing in the campaign needed one, and nothing was made.

**Seven pre-existing untracked files** left alone; `git status` is otherwise clean at every commit.

**Install.** The installed `forestsearch` did **not** match HEAD: `.fs_apply_mr()` was the
pre-O-1 version (it did not forward `field_recovery`, so the nine recovery columns could not have
been produced on the DINA path). `devtools::install(dependencies = FALSE)` was run, as the protocol
authorizes.

> **The package was reinstalled from commit `50bdeb8a`** (`dinamr: task document as received, with
> the kickoff amendments appended`) — HEAD at the moment the install ran. That commit touches only
> `dev/tasks/`, so **the `R/` tree it installed is byte-identical to the O-1 commit `1d9401cb`**:
> `git rev-parse 1d9401cb^{tree}:R` and `git rev-parse 50bdeb8a^{tree}:R` both give
> `aa9c3393eebcfb908895727f1eb8abee37ac487b`, and `git diff 1d9401cb 50bdeb8a -- R/ DESCRIPTION
> NAMESPACE` is empty. Every campaign cell was produced by that build; no commit since has touched
> `R/`.

Verified afterwards on the **installed** package:
`identical(deparse(installed .fs_apply_mr), deparse(source))` is **TRUE**, and all **25** of
`fs_mr_inference()`'s formals are forwarded (**none** missing).

---

## Part T — the `FS_S7_METHOD` engine knob

`subgroup_method` was a hard-coded literal at `sim_fs_maxeffCons_fb_mr_field_m1_template.qmd:300`.
Three touch points changed, all document-level and add-only:

1. `subgroup_method <- .env_chr("FS_S7_METHOD", "consistency")` with
   `stopifnot(subgroup_method %in% c("consistency","dina","grf"))` — verified to fire on a bad value.
2. The knob audit line now leads with `method=%s` (`Template knobs: method=dina hr=1.50 …`).
3. The `@sec-counts` settings readout gained a `subgroup_method (engine)` row.

**Nothing else changed.** In particular, `meta$subgroup_method` needed **no** edit: it already
travelled in **both** the batch and the combined bundle and was already in combine mode's poolability
key vector. That was verified against the committed `cert20` and `e1stud` bundles rather than
assumed. `method_tag`, `focus_tag` and the output stem already key on `subgroup_method`, so a DINA
run cannot overwrite a consistency bundle.

### Gate T — **PASS as a gate substitution, ratified by Larry**

**Recorded as a substitution, not as a literal pass.** The literal criterion — every non-timing
column and `truth` `identical()` to the committed `e1stud` rows 1–5 — is **unevaluable across
architectures**: 80 of 153 shared non-timing columns differ, all numerically, at **≤ 8.34e-15**,
while **every discrete and identification column is `identical()`**. The **same-Mac pre/post form**
was used instead — the HEAD template against the edited template, same host, same environment —
giving **162 of 162 non-timing columns, `truth`, and all 34 meta keys `identical()`**. Larry
ratified this substitution. The evidence for both halves follows.

5 replicates at the standing identity cell (`effMaxSG` ε 0.20, HR 1.50, n 500, `FS_S7_Z1Q=0.60`,
seeds 8316951 + sim_id, sim_id 1–5, tag `methknob`), `FS_S7_METHOD` **unset**.

**The decisive check — the HEAD template against the edited template, rendered back to back on this
host with an identical environment:**

| | |
|---|---|
| non-timing columns compared | **162** |
| `identical()` | **162 of 162** |
| `truth` `identical()` | **TRUE** |
| meta keys compared (`campaign_tag`, `built_at` excluded) | **34**, all `identical()` |

**Part T is byte-identical when the knob is unset.**

**The literal check against the committed `e1stud` rows 1–5 is NOT evaluable here, and the reason is
host arithmetic, not Part T.** `e1stud` was built on pop-os (x86-64, R 4.6.1, 100 workers); this is
Mac-Studio-3 (arm64, R 4.5.2, Accelerate). Measured:

- 80 of 153 shared non-timing columns disagree — **every one of them numerically**, at a **maximum
  relative difference of 8.34e-15**. That is last-bit BLAS/FMA, three orders of magnitude *inside*
  the campaign's own 1e-12 bound-identity tolerance, and `all.equal()` at default tolerance is
  **TRUE** over the whole frame.
- **Every discrete and identification column is exactly `identical()`**: `sim_id`, `n_true`,
  `detected`, `status`, `sg_def`, `n_sel`, `n_family`, `n_cons_qual`, `band_n`, `mr_ok`,
  `mr_harm_flag`, `covs`.
- **`truth`, with the magnitudes beside the verdicts:** `all.equal` at tolerance 1e-8 **TRUE**,
  `all.equal` at default tolerance **TRUE**, `identical()` **FALSE**. **Maximum absolute
  discrepancy 5.55e-15; maximum relative discrepancy 3.25e-15.** Only **2 of the 5** truth
  components move at all — `cde_H` (abs 5.55e-15, rel 3.25e-15) and `cde_Hc` (abs 1.44e-15, rel
  2.20e-15), the two that come from numerical integration. The three marginal Cox targets
  (`hr_causal`, `marg_H`, `marg_Hc`) are **exactly `identical()`** to all 17 digits.
- The bundle has **167** columns against `e1stud`'s 158: the nine `fld_recov_*` columns postdate
  `e1stud`, exactly as `REPORT_o1_forwarding_2026-09-10` §V4 records.

Substituting the decisive check for the unevaluable literal one follows the precedent
`REPORT_o1_forwarding` set for its own **Fb** and **Fe**. **The gate was not stopped on a
measurement artifact that was falsified in the same breath.**

> **Incidental finding, against `REPORT_o1_forwarding_2026-09-10` side issue 2.** That report could
> not reproduce a committed bundle's DGM draws on Linux (`n_true` differed at h175). **Here the draws
> reproduce exactly**: `n_true` is `154 172 162 171 155` on both the fresh Mac render and the
> committed `e1stud` bundle. Whatever blocked the h175 reproduction there does not bind at this cell
> on this host. Recorded, not chased.

---

## Stage 1 — the two mandated smokes

5 replicates each, full campaign knob set
(`FS_S7_METHOD=dina FS_S7_FOCUS=effMaxSG FS_S7_NBHD=0.20 FS_S7_FIELD_COMPLEMENT=TRUE
FS_S7_FIELD_SCALEC=selected FS_S7_FIELD_DECOMP=TRUE FS_S7_FIELD_RECOV=TRUE
FS_S7_IJ_RESIDUAL=two_term FS_S7_FB=none`), tag `dinamrsmk`.

| | **12.4%, HR 1.50, n 500** | **31%, HR 1.50, n 1500** |
|---|---|---|
| realized super-population prevalence | **0.12418** | **0.30655** |
| realized trial prevalence (mean `n_true`/n) | 0.1340 | 0.3151 |
| detected | 5 / 5 | 5 / 5 |
| **proposed-family size per replicate** | **112, 16, 17, 386, 54** | **1701, 3587, 3408, 2181, 1691** |
| `n_sel` | 84, 78, 115, 159, 87 | 522, 639, 261, 481, 569 |
| seconds per replicate | 5.53, 3.56, 3.13, 8.53, 4.38 | 28.10, 49.55, 47.83, 32.58, 27.87 |
| finiteness, 41 constructions | **ALL FINITE** | **ALL FINITE** |
| 13 interval / range invariants | **all OK** | **all OK** |
| γ joint / joint-s | 0.025 / 0.025, in range | 0.025–0.026 / 0.025, in range |
| bound ↔ quantile max \|diff\| | 0 | 0 |
| **nine recovery columns** | **9/9 PRESENT, none all-NA** | **9/9 PRESENT, none all-NA** |
| **p̂ block** (`p_hat_H/sum/top1`) | **PRESENT, populated** | **PRESENT, populated** |
| **ρᶜ** (`fld_Hc_scale_ratio`) | **PRESENT, populated** | **PRESENT, populated** |
| ρᶜ values | 0.9623 1.0002 1.0348 1.1057 1.0060 | 1.1292 1.1870 0.9821 1.1083 1.0655 |
| p̂(Ĥ) values | 0.3106 0.9488 0.4480 0.0320 0.2822 | 0.0790 0.0128 0.0180 0.1512 0.0256 |

**The O-1 fix took on this path.** Nothing was absent; no STOP was triggered.

`n_cons_qual` and `band_n` are present and all-NA — **structural** on DINA (both read
`grp.consistency$out_sg$result`, which the branch never produces). `p_star` is **not a recorder
column** at all; it is the admission set's consistency term, `NULL` on DINA by design. All three are
reported as structural, never as failures.

---

## Gate 1 — the projection, and the host

### Amendments carried

- **AMENDMENT 1** (kickoff, superseding the task document's §Stage 1): the ceiling is **9 h wall,
  not 10 h**; hard timeout 12 h. **Source, quoted:** FS projections have run **+16% (`cert20`) to
  −29% (`p12ext`)** against realized, and DINA's ~100× family-size skew is less predictable than
  FS's, so the wider margin is deliberate.
- **AMENDMENT 4** (Larry, mid-run, superseding Amendment 1): the ceiling is raised to **13 h wall**
  and the hard timeout to **16 h**, applied to the projection and to every per-cell re-projection
  including the post-Block-A checkpoint. **Stated intent:** Blocks A and B complete unattended and
  **Block C defers to a follow-up session**; the stated defer order stands, so Block C goes first if
  anything must be dropped. **The replicate count is unchanged at 2,000 per cell** — comparability
  with the 2,000-replicate FS grid is the point of this campaign and is not to be traded for cell
  count. Any cell deferred under the 9 h ceiling was to be re-instated if it now fits; **none had
  been deferred at the time the amendment arrived**, so there was nothing to re-instate.
- **AMENDMENT 3** (kickoff, as corrected by Larry mid-run): at Gate 2, where a committed FS
  comparator shares the DGM draws, assert **`n_true` `identical()` on all 2,000 rows** but compare
  **`truth` with `all.equal` at tolerance `1e-8`, NOT `identical()`** — `cert20`, `tier2` and
  `p12ext` are Linux-produced and cross-machine BLAS moves `truth` at ~1e-16, as tier2's own Stage 1
  recorded. **The tolerance used is 1e-8**, and `identical()` is reported beside it for information.
  A mismatch is a **finding about the DGM path, not a cell failure**.

### The host, and why the ceiling needed re-projecting

**The task document names a Linux executor and 100 workers. This ran on Mac-Studio-3.** No Linux
host is reachable from this session (no `~/.ssh/config`, no remote). The template caps workers at
`min(FS_S7_WORKERS, physical cores − 1)` = **13** here, so `FS_S7_WORKERS=100` would silently have
become 13. **Workers were set to 12**, matching tier2's headroom choice on this machine, with
`VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1` as tier2 did. Memory was sampled
during the run: **13 R processes at 5.5–16.5 GB total**, against tier2's measured 19.3 GB over 13
processes and the standing 24 GB rule — comfortable, and never close.

**The 9 h ceiling of Amendment 1 was calibrated for Linux at 100 workers.** At 12 workers this
machine is ~8× narrower, so the ceiling had to be read against *this machine's* measured walls rather
than the document's. Amendment 4 raised it to 13 h with the explicit expectation that **a partial
campaign is the intended outcome, not a failure.**

### The cost surface — measured, not assumed

The task forbids projecting from a mean of a 5-replicate smoke, and forbids using the FS walls. Ten
**36-replicate probes at 12 workers** were run at the corners instead, tagged `dinamrprobe` so they
can never pool with the campaign:

| block | HR | n | detection | K med | K q90 | K max | s/rep med | s/rep mean | s/rep q90 |
|---|---|---|---|---|---|---|---|---|---|
| A (12.4%) | 1.00 | 500 | 0.778 | 76.0 | 947.9 | 1195 | 4.421 | 5.269 | 11.886 |
| A (12.4%) | 1.00 | 1500 | 0.361 | 35.0 | 208.4 | 480 | 0.113 | 1.938 | 6.417 |
| A (12.4%) | 1.50 | 500 | 0.917 | 193.0 | 1275.0 | 1695 | 6.446 | 8.366 | 20.576 |
| A (12.4%) | 1.50 | 1000 | 0.833 | 246.0 | 757.4 | 1523 | 8.021 | 7.535 | 13.648 |
| A (12.4%) | 1.50 | 1500 | 0.778 | 128.5 | 435.1 | 1001 | 5.423 | 5.410 | 9.851 |
| B (31%) | 1.00 | 500 | 0.944 | 359.5 | 1634.1 | 3053 | 9.479 | 11.767 | 24.064 |
| B (31%) | 1.00 | 1500 | 0.944 | 258.5 | 1203.8 | 1699 | 8.450 | 10.165 | 18.955 |
| B (31%) | 1.50 | 500 | 1.000 | 1335.5 | 2468.0 | 3332 | 22.013 | 24.689 | 38.461 |
| B (31%) | 1.50 | 1000 | 1.000 | 2008.0 | 3086.5 | 3370 | 30.970 | 31.412 | 45.785 |
| B (31%) | 1.50 | 1500 | 1.000 | 1780.5 | 3123.5 | 3638 | 33.906 | 34.395 | 51.505 |

Two things in that table matter more than the timing.

1. **DINA's detection rate falls with n at 12.4%** — 0.917 → 0.833 → 0.778 at HR 1.50, and 0.778 →
   0.361 at HR 1.00. At 31% it is 1.000 everywhere at HR 1.50. This is the opposite of FS's
   behaviour and it is why the detection rate is stated before every coverage number.
2. **The family is an order of magnitude larger at 31% than at 12.4%** (median ~1300–2000 against
   ~130–250) and it is where the cost is. The pilot's "median 91, max 1745" at 12.4%/n500/HR 1.75 is
   a *lower* corner of this surface, not a typical one: **max 3638** is reached at 31%/n1500.

Projection per cell = `2000 × mean(fit_mr_secs) / 12` worker-seconds plus a measured ~30 s render
overhead × 3 renders (two batches plus a combine). HR 1.75 is costed at its HR 1.50 corner (both
detect near-always; the family is the admission-set candidate table, which the effect target barely
moves). The HR 1.00 n = 1000 cells are interpolated between the measured n 500 and n 1500 corners.

| campaign block | projected wall |
|---|---|
| **A** (12.4%, HR 1.50 & 1.75, n 500/1000/1500) | **2.12 h** |
| **B** (31%, HR 1.50 & 1.75, n 500/1000/1500) | **8.53 h** |
| **C** (HR 1.00, both prevalences, n 500/1000/1500) | **2.17 h** |
| **A + B** | **10.65 h** |
| **all 18 cells** | **12.83 h** |

### Gate 1 decision

**PROCEED with Blocks A and B (12 cells, 10.65 h projected, 18% headroom under the 13 h ceiling).
DEFER Block C (all six cells), per the stated defer order — Block C first, 31% before 12.4%.**

The full 18-cell campaign projects at **12.83 h against a 13 h ceiling — a 1.3% margin**, which is
not a margin at all against a projection error band measured at **+16% to −29%** on FS and a DINA
family-size skew that Amendment 1 itself calls less predictable. Deferring Block C is also the
**explicitly stated intent of Amendment 4**. Replicate count stays at **2,000 per cell**; no cell
was thinned.

Hard timeout: **16 h**.


---

## Where each instruction is recorded

Every instruction that governed this campaign, and the artifact that carries it. Nothing was
followed from conversation alone.

| Instruction | Source | Recorded in |
|---|---|---|
| The campaign as specified (Parts T and C, blocks, defer order, Stage 1/Gate 1/Gate 2/Stage 3, "Done means") | Task document `TASK_dinamr_campaign_2026-09-10.md` | `dev/tasks/TASK_dinamr_campaign_2026-09-10.md`, committed verbatim at `50bdeb8a` |
| **Amendment 1** — Gate 1 ceiling 9 h (not 10 h), hard timeout 12 h, with its source (FS projections ran +16% `cert20` to −29% `p12ext`; DINA's ~100× family skew less predictable) | Kickoff paste | Task-file **Appendix**, committed at `50bdeb8a`; and @Gate 1 of this report |
| **Amendment 2** — one post-Block-A checkpoint re-projecting B and C from A's realized walls; once, not per cell | Kickoff paste | Task-file **Appendix** (`50bdeb8a`); @Block A checkpoint below |
| **Amendment 3** — same-draws assertions vs the committed FS comparator; a mismatch is a DGM-path finding, not a cell failure | Kickoff paste | Task-file **Appendix** (`50bdeb8a`); @Gate 1 and the Gate 2 records |
| **Mac correction 1** — Amendment 3 refined: `n_true` stays `identical()`, but `truth` by `all.equal` at ~1e-8, **not** `identical()`, because the FS comparators are Linux-produced and cross-machine BLAS moves `truth` at 1e-16 (as tier2's own Stage 1 recorded); report the tolerance used | Larry, mid-run message | @Gate 1 "Amendments carried"; implemented in `gate2.R` as `TOL_TRUTH <- 1e-8`; **the tolerance is printed in every Gate 2 cell record** |
| **Mac correction 2** — 12 workers, not 100, matching tier2's headroom choice here (19.3 GB RSS over 13 processes); `VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1` as tier2 did | Larry, mid-run message | @Gate 1 "The host"; the three vars are exported by the render driver; `n_workers = 12` is in **every batch meta** and gated at Gate 2 |
| **Mac correction 3** — note that the 9 h ceiling was calibrated for Linux at 100 workers; re-project from this machine's measured walls; expect B and C to defer per the stated order; **a partial campaign is the expected outcome, not a failure** | Larry, mid-run message | @Gate 1 "The host, and why the ceiling needed re-projecting" |
| **Amendment 4** — ceiling raised to 13 h, hard timeout to 16 h, applied to the projection and every per-cell re-projection including the checkpoint; intent that A and B complete and C defers to a follow-up; stated defer order stands; **2,000 replicates per cell is not tradeable for cell count**; re-instate any cell deferred under 9 h if it now fits | Larry, mid-run message | Task-file **Appendix**; @Gate 1 "Amendments carried" and the Gate 1 decision |
| **Checkpoint pin** — at the post-Block-A checkpoint, **Block C stays deferred whatever the re-projection shows**; the checkpoint decides Block B's cells only | Larry, mid-run message | @Block A checkpoint below |
| **Gate T substitution ratified** — the literal `e1stud` identity is unevaluable across architectures; the same-Mac pre/post form stands as the gate | Larry, mid-run message | @Gate T, headed "PASS as a gate substitution, ratified by Larry" |
| Additional record items (install SHA; truth discrepancy magnitudes; non-detection split; FS family beside DINA's; the `summary_dinamr` diffstat; this table) | Larry, mid-run message | This section and the four that follow |

---

## The `summary_dinamr.qmd` transplant, diffstat against `summary_cert20.qmd`

```
git diff --no-index --stat summary_cert20.qmd summary_dinamr.qmd
 summary_cert20.qmd => summary_dinamr.qmd | 794 ++++++++++++++++++-----------
 1 file changed, 542 insertions(+), 252 deletions(-)
```

391 lines → **681 lines**. Chunk inventory, 17 → **25**:

| | chunks |
|---|---|
| **Retained** (13, same label and same machinery) | `setup` `cov` `cov-across` `wilson-fn` `vsn-h` `vsn-h-plot` `vsn-c` `ident` `ident-plot` `tert` `tert-plot` `null` `record` |
| **Dropped** (3) | `accept` — **the acceptance-criteria section, removed because this campaign has none**; `cov-plot-c` and `cov-plot-h`, merged into one guarded `cov-plots` |
| **Added** (12) | `inventory` `detection` `detection-plot` `structna` `cov-fns` `cov-plots` `twosided` `twosided-plot` `family` `family-plot` `fsgrid` `recov` |

The transplant is globs, labels and comparator names, plus: the grid widened from 9 cells at one
prevalence to **18 across two**; every coverage caption relabelled as the **conditional-on-proposed-
family** estimand; the acceptance-criteria section removed; and four sections the task asks for that
`cert20` has no counterpart to — detection stated before any coverage number (@sec-detect), the
structurally-NA columns named as such (@sec-structna), the proposed-family size distribution and its
stabilization (@sec-family), the FS grid beside with the confound stated (@sec-fsgrid), and the
two-sided decay question with the miss split by side (@sec-twosided).

**The absent-cell guard.** Every chunk that could render empty is wrapped in `have()` /
`skip_note()`. Verified by rendering `summary_dinamr.qmd` with **zero** campaign cells on disk: all
20 guarded chunks skipped with a named note, **none rendered empty** — the A3 behaviour seen on the
`e1stud` set does not recur — while the committed FS comparator grid still rendered, because those
bundles do exist.

---

## Non-detections: the split the task asks for is **not available from the bundle**

The instruction is to split non-detections into *empty proposed family* versus *candidates proposed
but none admitted*, **where the bundle allows**. It does not allow it, and the mechanism is exact
rather than inferred:

```r
  if (!found) { rec$status <- "NO-DETECTION"; return(rec) }   # template line 1049
  ...
  rec$n_family <- g$n_family %||% NA_integer_                  # template line 1055
```

The recorder **returns its all-NA record before `n_family` is ever written**, and `n_family` is read
off `g`, the MR gate object, which does not exist when nothing is selected. So on a non-detection
`n_family` is **NA — never 0, never positive** — and the two causes are indistinguishable.

Measured, not assumed. Across the ten 36-replicate probes: **52 non-detections, `n_family` NA on all
52, `== 0` on none, `> 0` on none, `err_msg` set on none** (so these are genuine no-selections, not
errors), `status` `NO-DETECTION` throughout. The same holds on every campaign cell — e.g. Block A
HR 1.50 n 500: 250 non-detections, `n_family` NA on all 250. Each Gate 2 cell record prints this
split with the reason.

**Separating the two would need an `R/` or recorder change, which this task forbids.** Flagged, not
fixed.

---

## FS's enumerated family beside DINA's — the sharpest contrast in the campaign

Both engines record `n_family`, so the comparison is read off committed bundles with no new
machinery. FS's numbers are from `p12ext` / `tier2` / `cert20` / `e1stud` **as they stand**.

| block | HR | n | FS campaign / focus / ε | FS K: min / med / q90 / max, **CV** | FS detection | DINA K: min / med / q90 / max, **CV** | DINA detection |
|---|---|---|---|---|---|---|---|
| A (12.4%) | 1.50 | 500 | p12ext, maxeffCons, 0.10 | 1037 / 1223 / 1318 / 1380, **0.053** | 0.9110 | 1 / 253.5 / 1111 / 3229, **1.201** | 0.8750 |
| A (12.4%) | 1.50 | 1000 | p12ext, maxeffCons, 0.10 | 1119 / 1299 / 1398 / 1418, **0.041** | 0.9740 | 1 / 156 / 628 / 2596, **1.185** | 0.8385 |
| A (12.4%) | 1.50 | 1500 | p12ext, maxeffCons, 0.10 | 1195 / 1297 / 1396 / 1412, **0.036** | 0.9880 | 1 / 103 / 394 / 1954, **1.142** | 0.7985 |
| A (12.4%) | 1.75 | 500 | tier2, maxeffCons, 0.10 | 1037 / 1223 / 1318 / 1380, **0.054** | 0.9500 | *(pending)* | |
| A (12.4%) | 1.75 | 1000 | tier2, maxeffCons, 0.10 | 1119 / 1299 / 1398 / 1418, **0.041** | 0.9950 | *(pending)* | |
| A (12.4%) | 1.75 | 1500 | tier2, maxeffCons, 0.10 | 1195 / 1297 / 1396 / 1412, **0.036** | 0.9990 | *(pending)* | |
| B (31%) | 1.50 | 500 | e1stud, effMaxSG, 0.20 | 1036 / 1223 / 1318 / 1380, **0.054** | 0.9995 | *(pending)* | |
| B (31%) | 1.50 | 1000 | cert20, effMaxSG, 0.20 | 1119 / 1299 / 1398 / 1418, **0.041** | 1.0000 | *(pending)* | |
| B (31%) | 1.50 | 1500 | cert20, effMaxSG, 0.20 | 1195 / 1297 / 1396 / 1412, **0.036** | 1.0000 | *(pending)* | |
| B (31%) | 1.75 | 500 | e1stud, effMaxSG, 0.20 | 1036 / 1223 / 1318 / 1380, **0.054** | 0.9995 | *(pending)* | |
| B (31%) | 1.75 | 1000 | cert20, effMaxSG, 0.20 | 1119 / 1299 / 1398 / 1418, **0.041** | 1.0000 | *(pending)* | |
| B (31%) | 1.75 | 1500 | cert20, effMaxSG, 0.20 | 1195 / 1297 / 1396 / 1412, **0.036** | 1.0000 | *(pending)* | |

Three facts, stated without interpretation beyond the record:

1. **FS's family is nearly deterministic and prevalence-invariant.** Its CV is 0.036–0.054, and its
   quantiles are **identical across the two prevalence blocks at matched n** (1037/1223/1318/1380 at
   n 500 in both A and B; the one-unit min difference at B is the only movement). It is a
   combinatorial enumeration off the covariate grid, so only **n** moves it — and it **grows**
   slightly with n.
2. **DINA's family is data-adaptive, ~25× more variable, and moves the other way.** CV 1.14–1.20, a
   floor of **1**, a ceiling of 3229, and a median that **shrinks with n** (253.5 → 156 → 103 at
   HR 1.50 in Block A) while FS's grows.
3. **Detection moves in opposite directions with n.** FS at 12.4%/HR 1.50: 0.911 → 0.974 → 0.988,
   **rising**. DINA at the same cells: 0.875 → 0.839 → **0.799, falling**.

FS additionally records `n_cons_qual` — the consistency-qualifying count the identifier sorts —
which is **structurally NA on DINA** (@sec-structna). For reference it is 17 / 25 / 31 at Block A
HR 1.50 and 174 / 292 / 357 at Block B HR 1.50. There is no DINA counterpart to place beside it.

---

## `summary_dinamr.qmd`, chunk by chunk against `summary_cert20.qmd`

### The 13 retained chunks — globs, labels and guards, **plus the following**

They do **not** differ only in globs, labels and absent-cell guards. Every additional change is
listed; there are no others.

| chunk | globs / labels / guards | anything else |
|---|---|---|
| `setup` | ✓ | Adds `BLK`, `A_labs`, `B_labs` — the two-prevalence-block grid the 9-cell original had no need of. **The inventory table was moved out** of this chunk into the new `inventory` chunk. |
| `cov` | ✓ | **The three helpers `subst_s` / `err_sd` / `cov_block` were moved out** into the new `cov-fns` chunk, unchanged in body. Adds a `block_grp` column; the single table becomes a per-block `for` loop over A and B. |
| `cov-across` | ✓ | `split()` now also keys on `block_grp`; a `prevalence_block` column is added to the output. |
| `wilson-fn` | ✓ | **Adds `miss_below` and `miss_above`** — the two-sided miss split by side, which Stage 3 item 5 asks for and `cert20` does not compute. Adds a `block` column; adds `if (!nrow(d)) return(NULL)`. |
| `vsn-h` | ✓ | Adds `block`, `miss_below`, `miss_above` to the displayed columns. |
| `vsn-h-plot` | ✓ | `facet_wrap(~ construction)` → `facet_grid(blk ~ construction)`, so both prevalence blocks show. |
| `vsn-c` | ✓ | **Nothing else.** |
| `tert` | ✓ | `need` is intersected with `names(r)` and five columns dropped from the requirement (`fld_H_lam_mean`, `mr_H_se_ij`, `mr_Hc_se_ij`, `nv_H_se`, `nv_Hc_se`) — they are *required* by `cert20`'s filter but never *used* in the chunk, and requiring them would drop replicates for no reason. Adds an "all (p-hat not splittable)" fallback for a degenerate p̂ distribution. Drops the unused intermediate `d$ef`. `sd()` → `stats::sd()`. |
| `tert-plot` | ✓ | Local `block` renamed `blk` (it collided with the prevalence block); filters out the "not splittable" row; x-axis labels angled, since 12 cells do not fit where 9 did. |
| `ident` | ✓ | Adds `block` and `K_med`. **Drops `qual_mean`** (mean `n_cons_qual`) — *deliberate and stated in the caption*: it is structurally NA on DINA, so the column would be a row of NAs. |
| `ident-plot` | ✓ | **Adds PPV and NPV panels** (4 panels → 6), per the standing convention that NPV sits beside sens/spec/PPV; adds a prevalence-block linetype. |
| `null` | ✓ | Adds `block` and `K_med`; adds `if (!nrow(d)) return(NULL)`. |
| `record` | ✓ | Adds a **detection-and-family-size headline** block first; adds the **two-sided miss split by side**; adds a **family-size / stabilization** block; adds the closing "no acceptance criteria" paragraph. The `cert20` acceptance-criteria loop is removed with the `accept` chunk. |

> **Two unforced losses found by this accounting and restored.** The first draft of `tert` had
> dropped `cov_field_C` (the *unscaled* field's upper coverage, which `cert20` reports beside
> field-s) and `ident` had dropped `absHtrue_mean` (mean `n_true`). Neither had a reason — they were
> transplant slips, not decisions. Both are **restored**, so the only column `dinamr` drops relative
> to `cert20` is `qual_mean`, and that one has a reason.

### The 3 dropped chunks, with reasons

| chunk | reason |
|---|---|
| `accept` | **The acceptance-criteria section.** The task states "No acceptance criteria. Exploratory. FS certified ranges as reference lines only; no recommendation." None was pre-registered for this campaign, so there is nothing to evaluate and the chunk is removed rather than left to score against `cert20`'s bars. |
| `cov-plot-c` | Merged into the single guarded `cov-plots`, which loops over the complement and harm sides. Same two figures, one guard. |
| `cov-plot-h` | As above. |

### The 12 added chunks, each mapped to the Stage 3 item it serves

| chunk | Stage 3 item served |
|---|---|
| `inventory` | Gate 2 record carried into the document — bundle provenance, knobs from the pooled meta, realized prevalence. (Split out of `setup`.) |
| `detection` | **"the detection rate … recorded prominently"** — Gate 2, and the framing that every summary conditions on it. |
| `detection-plot` | Same, against n and prevalence. |
| `structna` | **"the structurally-NA columns (`n_cons_qual`, `band_n`, `p_star`) reported as such rather than as failures."** |
| `cov-fns` | Mechanical: the helpers `cov` used, split out so the per-block loop and the guard read cleanly. Serves item 1. |
| `cov-plots` | Item 1, the standard tables' figures (merge of `cov-plot-c` / `cov-plot-h`). |
| `twosided` | **Item 5** — "does the harm-block two-sided decay with n at 12.4% appear on DINA as it does on FS … with the miss split by side", with the committed FS trajectory (0.977/0.932/0.901; 0.961/0.917/0.913) beside it. |
| `twosided-plot` | Item 5, the same comparison as a figure. |
| `family` | **Item 3** — "the proposed family: size distribution per cell, and whether it stabilizes with n — the mechanism Supplementary S8.3 offers." |
| `family-plot` | Item 3, the distribution and the relative-spread-versus-n read. |
| `fsgrid` | **Item 4** — "beside the FS grid (`cert20` / `tier2` / `p12ext` read as they stand), with the confound stated." |
| `recov` | **Item 6** — "p̂ and the recovery diagnostics on a model-generated family: distributions per cell … `sens_H` and containment quantiles." |

Item 2 (classification vs n **and prevalence**) is served by the retained `ident` / `ident-plot`
widened to two blocks; item 7 (no acceptance criteria) is served by the removal of `accept`.

---

## Gate 2 — Block A, all six cells: **PASS**, 204 checks, **zero failures**

Checker `gate2.R`; full per-cell record in the campaign log. Gated on the pooled 2,000-row bundle:
completeness (2,000 rows, `sim_id` 1–2000 with no duplicates, no `CONFIG-ERROR`, every campaign knob
as set in the pooled meta including `subgroup_method = "dina"`, two seed-disjoint batches,
`n_workers = 12` and `forestsearch_version = 0.3.5` in **each batch** meta); detection rate; family
distribution; realized prevalence; finiteness over 38 constructions on every detected replicate; ten
interval invariants; γ ∈ [0.025, 0.05]; the bound↔quantile identities; p̂ validity; the
structurally-NA columns reported as such; and the Amendment 3 same-draws assertions.

| cell | rows | detection | family K: min / med / q90 / max (CV) | realized prevalence | γ range (share at floor) | field-s bdc identity | same draws |
|---|---|---|---|---|---|---|---|
| HR 1.50 n 500 | 2000 | **0.8750** | 1 / 253.5 / 1111 / 3229 (1.201) | 0.12363 | 0.0250–0.0280 (0.925) | 2.22e-16 | **yes** |
| HR 1.50 n 1000 | 2000 | **0.8385** | 1 / 156 / 628 / 2596 (1.185) | 0.12388 | at floor 0.93 | ≤ 2.3e-16 | **yes** |
| HR 1.50 n 1500 | 2000 | **0.7985** | 1 / 103 / 394 / 1954 (1.142) | 0.12401 | at floor 0.93 | ≤ 2.3e-16 | **yes** |
| HR 1.75 n 500 | 2000 | **0.9080** | 1 / 332 / 1273 / 3246 (1.105) | 0.12363 | at floor 0.92 | ≤ 2.3e-16 | **yes** |
| HR 1.75 n 1000 | 2000 | **0.9065** | 1 / 232 / 807 / 2721 (1.059) | 0.12388 | at floor 0.92 | ≤ 2.3e-16 | **yes** |
| HR 1.75 n 1500 | 2000 | **0.8980** | 1 / 177 / 526 / 2315 (1.001) | 0.12401 | at floor 0.92 | ≤ 2.3e-16 | **yes** |

**Amendment 3 holds on all six.** `n_true` `identical()` on all 2,000 rows of every cell; `truth`
`all.equal` at **tolerance 1e-8** TRUE on every cell. `identical()` is FALSE on the three HR 1.50
cells at **max abs 6.66e-16 / max rel 3.99e-16** and **exactly TRUE (0 / 0) on all three HR 1.75
cells** — the HR 1.75 comparators are `tier2`, whose `truth` reproduces bit-for-bit here, while the
HR 1.50 comparators are `p12ext`, whose two CDE components do not. **No mismatch anywhere**, so
nothing had to be recorded as a DGM-path finding.

**Classification, and the finding that matters most.**

| cell | sens | spec | PPV | NPV | mean \|Ĥ\| | non-detections |
|---|---|---|---|---|---|---|
| HR 1.50 n 500 | 0.530 | 0.852 | 0.349 | 0.929 | 97.9 | 250 |
| HR 1.50 n 1000 | 0.650 | 0.895 | 0.492 | 0.948 | 173.1 | 323 |
| HR 1.50 n 1500 | 0.727 | 0.903 | 0.546 | 0.959 | 263.1 | 403 |
| HR 1.75 n 500 | 0.593 | 0.867 | 0.408 | 0.939 | 95.7 | 184 |
| HR 1.75 n 1000 | 0.726 | 0.916 | 0.582 | 0.960 | 164.5 | 187 |
| HR 1.75 n 1500 | 0.804 | 0.928 | 0.643 | 0.971 | 245.0 | 204 |

**Classification improves with n on every axis while detection falls.** At HR 1.50 sensitivity goes
0.530 → 0.650 → 0.727 and PPV 0.349 → 0.492 → 0.546, yet detection goes 0.875 → 0.839 → **0.799**.
More data makes DINA *better at describing the region it finds* and *less likely to return one at
all*. Every non-detection is a genuine no-selection — `err_msg` is set on none of the 1,551 across
the six cells.

---

## The Block A checkpoint (Amendment 2) — run once

> **Pinned:** Block C stays deferred to a follow-up session **whatever this re-projection shows**
> (Larry). This checkpoint decides **Block B's cells only**.

| Block A cell | projected (Gate 1) | **realized** | ratio |
|---|---|---|---|
| HR 1.50 n 500 | 24.7 min | **30.7 min** | 1.240 |
| HR 1.50 n 1000 | 22.4 min | **23.3 min** | 1.040 |
| HR 1.50 n 1500 | 16.5 min | **19.8 min** | 1.196 |
| HR 1.75 n 500 | 24.7 min | **34.4 min** | 1.389 |
| HR 1.75 n 1000 | 22.4 min | **28.2 min** | 1.257 |
| HR 1.75 n 1500 | 16.5 min | **25.7 min** | 1.556 |
| **total** | **2.123 h** | **2.701 h** | **1.272** |

Gate 1 under-projected Block A by **27%** overall, with a per-cell spread of 1.040–1.556 — outside
the +16%/−29% FS band that Amendment 1 quoted, and in the costly direction. The re-projection is
calibrated **by n**, since Block B's cost is n-driven: K₅₀₀ = 1.314, K₁₀₀₀ = 1.149, K₁₅₀₀ = 1.376.

| Block B cell | Gate 1 | **re-projected** |
|---|---|---|
| HR 1.50 n 500 | 70.1 min | **92.1 min** (1.535 h) |
| HR 1.50 n 1000 | 88.8 min | **102.0 min** (1.699 h) |
| HR 1.50 n 1500 | 97.0 min | **133.5 min** (2.226 h) |
| HR 1.75 n 500 | 70.1 min | **92.1 min** (1.535 h) |
| HR 1.75 n 1000 | 88.8 min | **102.0 min** (1.699 h) |
| HR 1.75 n 1500 | 97.0 min | **133.5 min** (2.226 h) |
| **total** | 8.529 h | **10.919 h** |

```
Ceiling 13 h  −  Block A realized 2.701 h  =  remaining budget 10.299 h
Block B re-projected in full                                  10.919 h
A + B in full                                                 13.620 h   -> OVER by 0.620 h (6.0%)
```

### Checkpoint decision

**Block B does not fit in full. One cell is deferred, per the stated order (B's n = 1500 first).**

- **RUN (5 cells, 8.694 h re-projected):** HR 1.50 at n 500 / 1000 / **1500**, HR 1.75 at n 500 /
  1000. With Block A's 2.701 h that is **11.395 h against the 13 h ceiling — 12% headroom**.
- **DEFER (1 cell):** **HR 1.75, n 1500** at 31%.
- **DEFER (6 cells, pinned):** all of **Block C**.

**Which n = 1500 cell to defer is a judgment the stated order does not settle**, so it is recorded
here. Deferring **HR 1.75** rather than HR 1.50 keeps the **complete HR 1.50 n-trajectory
(500 / 1000 / 1500) at both prevalences**, which is exactly what the Stage 3 two-sided decay question
needs — it is asked first at HR 1.50, and Block A already supplies the full HR 1.75 trajectory at
12.4%. Deferring HR 1.50 instead would have left the primary trajectory broken at 31% while
preserving a secondary one.

Dropping the **minimum** that makes it fit, rather than both n = 1500 cells, follows Amendment 4's
stated intent that Block B complete as far as it can; both would have cost a further 2.226 h of
science for no budget reason.

> **Ratified (Larry, after the checkpoint):** *"Keep HR 1.50 n1500 as you have it, continue."* The
> deferred Block B cell is therefore **HR 1.75, n 1500**, and the complete HR 1.50 n-trajectory at
> 31% is retained.

Block B launched immediately on this result, without waiting for the ratification. Hard timeout
stands at 16 h.

---

## The p̂ overlay: DINA against its FS comparator, described per cell

`@sec-fsphat` of `summary_dinamr.qmd` computes the p̂ tertile table **identically** on DINA and on
each cell's committed FS comparator — the bundle that shares the DGM draws — and overlays the two
retained-bias curves on **median p̂**. Each bundle is binned on **its own** within-cell empirical
tertiles (the FS `tert()` definition), so neither engine is forced onto the other's p̂ cut-points;
FS's p̂ distribution sits materially higher than DINA's, which is why the overlay must be read
against median p̂ rather than tertile index.

**The relationship changes with n, so it is described per cell rather than as one shape offset by a
constant.** At HR 1.50, 12.4%:

| | | p̂ T1 | p̂ T2 | p̂ T3 |
|---|---|---|---|---|
| **n = 500** | DINA median p̂ | 0.049 | 0.134 | 0.356 |
| | DINA retained bias | **+0.034** | **+0.145** | **+0.275** |
| | FS median p̂ | 0.072 | 0.166 | 0.325 |
| | FS retained bias | **−0.073** | **−0.118** | **+0.127** |
| | DINA − FS, by tertile index | +0.107 | +0.263 | +0.148 |
| **n = 1500** | DINA median p̂ | 0.079 | 0.213 | 0.504 |
| | DINA retained bias | **+0.027** | **+0.047** | **+0.056** |
| | FS median p̂ | 0.195 | 0.362 | 0.595 |
| | FS retained bias | **−0.239** | **−0.089** | **+0.069** |
| | DINA − FS, by tertile index | +0.266 | +0.136 | −0.013 |

- **At n = 500**, DINA's retained bias sits **above** FS's in **every** tertile, by **+0.11 to
  +0.26** read by tertile index.
- **At n = 1500**, DINA's curve is **nearly flat** — **+0.03 to +0.06** across median p̂ **0.08 to
  0.50** — while FS's **rises from −0.24 to +0.07** across median p̂ **0.20 to 0.60**. Read against
  median p̂, the gap is about **+0.29 near p̂ ≈ 0.2** and **closes near p̂ ≈ 0.5–0.6**.

> **At 12.4% this overlay is criterion-confounded**, and cannot be read as engine behaviour. The FS
> side is `maxeffCons` at ε 0.10 against DINA's `effMaxSG` at ε 0.20 — a different selection
> functional at half the band width, on top of the identifier, family-construction and
> detection-set differences that always apply (@sec-fsgrid).
>
> **Block B's overlay is the comparison to read.** There the FS comparators (`cert20` / `e1stud`)
> run **DINA's own criterion** — `effMaxSG` at ε 0.20 — so the criterion drops out and only the
> identifier, the family construction and the detection set remain.

---

## Gate 2 — Block B, all five cells run: **PASS**, 170 checks, **zero failures**

The sixth cell (HR 1.75, n 1500) is the checkpoint's deferred cell and correctly reports
**NOT ON DISK**, not a failure.

| cell | rows | detection | family K: min / med / q90 / max (**CV**) | realized prevalence | γ range (share at floor) | field-s bdc identity | same draws |
|---|---|---|---|---|---|---|---|
| HR 1.50 n 500 | 2000 | **0.9985** | 3 / 1495 / 2881 / 3556 (**0.564**) | 0.30591 | 0.0250–0.0270 (0.915) | 2.22e-16 | **yes** |
| HR 1.50 n 1000 | 2000 | **1.0000** | 48 / 1709.5 / 2984 / 3814 (**0.473**) | 0.30625 | 0.0250–0.0270 (0.907) | 2.22e-16 | **yes** |
| HR 1.50 n 1500 | 2000 | **1.0000** | 55 / 1743.5 / 2850 / 3794 (**0.433**) | 0.30646 | 0.0250–0.0270 (0.910) | 2.22e-16 | **yes** |
| HR 1.75 n 500 | 2000 | **1.0000** | 37 / 1954.5 / 2982 / 3557 (**0.449**) | 0.30591 | 0.0250–0.0270 (0.902) | 2.22e-16 | **yes** |
| HR 1.75 n 1000 | 2000 | **1.0000** | 241 / 2257.5 / 3160 / 3821 (**0.337**) | 0.30625 | 0.0250–0.0270 (0.904) | 2.22e-16 | **yes** |
| HR 1.75 n 1500 | — | — | — | — | — | — | **deferred** |

**Amendment 3 holds on all five.** `n_true` `identical()` on all 2,000 rows; `truth` `all.equal` at
**tolerance 1e-8** TRUE. `identical()` FALSE throughout, at **max abs 5.55e-15 / max rel 3.25e-15**
(HR 1.50, comparator `cert20` / `e1stud`) and **8.88e-15 / 4.25e-15** (HR 1.75) — the cross-machine
signature, confined as before to the two CDE components. **No mismatch anywhere.**

**Classification, Block B:** sens 0.497 / 0.642 / 0.788 and PPV 0.674 / 0.799 / 0.850 at HR 1.50
across n = 500 / 1000 / 1500, with NPV 0.808 / 0.863 / 0.918 and mean |Ĥ| 110 / 243 / 426.
**Detection is ~100% at every Block B cell**, against 0.80–0.91 in Block A — so the
conditioning-on-detection caveat, which bites hard at 12.4%, barely bites at 31%. Only 3
non-detections in the entire block (all at HR 1.50 n 500, all genuine no-selections).

### The proposed family stabilizes at 31% and does **not** at 12.4%

This is the mechanism Supplementary S8.3 offers for conditional coverage recovering with n, and the
two blocks answer it oppositely.

| block | HR | n 500 | n 1000 | n 1500 |
|---|---|---|---|---|
| **A (12.4%)** | 1.50 | med 254, CV **1.201**, IQR/med 2.026 | med 156, CV **1.185**, 1.917 | med 103, CV **1.142**, 2.000 |
| A (12.4%) | 1.75 | med 332, CV 1.105, 1.795 | med 232, CV 1.059, 1.595 | med 177, CV 1.001, 1.516 |
| **B (31%)** | 1.50 | med 1495, CV **0.564**, IQR/med 0.992 | med 1710, CV **0.473**, 0.790 | med 1744, CV **0.433**, 0.679 |
| B (31%) | 1.75 | med 1954, CV 0.449, 0.697 | med 2258, CV 0.337, 0.512 | *(deferred)* |

At **31%** the family **grows** with n (1495 → 1710 → 1744) while its relative spread **shrinks**
(CV 0.564 → 0.433; IQR/median 0.992 → 0.679) — it stabilizes. At **12.4%** the family **shrinks**
with n (254 → 156 → 103) and its relative spread **does not** (CV 1.20 → 1.14; IQR/median stays
≈ 1.9–2.0). The S8.3 mechanism is visible at 31% and absent at 12.4%. Reported, not interpreted
further.

### The matched-criterion comparison — the one to read

At 31% the FS comparators run **DINA's own criterion** (`effMaxSG`, ε 0.20), so the criterion drops
out and only identifier, family construction and detection set remain.

| cell | engine | criterion | field lower [Wilson] | field-s upper | IJ 2-sided | miss lo / hi | joint_s | sens | PPV |
|---|---|---|---|---|---|---|---|---|---|
| HR 1.50 n 500 | **DINA** | dinamr/effMaxSG/ε0.2 | **0.9159** [0.9029, 0.9273] | 0.9099 | 0.9905 | 0.0095 / 0.0000 | 0.9119 | 0.497 | 0.674 |
| | FS | e1stud/effMaxSG/ε0.2 | **0.9745** [0.9666, 0.9805] | 0.9125 | 0.9810 | 0.0035 / 0.0155 | 0.9420 | 0.590 | 0.706 |
| HR 1.50 n 1000 | **DINA** | dinamr/effMaxSG/ε0.2 | **0.9515** [0.9412, 0.9601] | 0.9365 | 0.9865 | 0.0065 / 0.0070 | 0.9455 | 0.642 | 0.799 |
| | FS | cert20/effMaxSG/ε0.2 | **0.9585** [0.9488, 0.9664] | 0.9420 | 0.9710 | 0.0055 / 0.0235 | 0.9500 | 0.708 | 0.804 |
| HR 1.50 n 1500 | **DINA** | dinamr/effMaxSG/ε0.2 | **0.9610** [0.9516, 0.9686] | 0.9420 | 0.9840 | 0.0045 / 0.0115 | 0.9535 | 0.788 | 0.850 |
| | FS | cert20/effMaxSG/ε0.2 | **0.9615** [0.9521, 0.9691] | 0.9465 | 0.9755 | 0.0030 / 0.0215 | 0.9505 | 0.839 | 0.837 |
| HR 1.75 n 500 | **DINA** | dinamr/effMaxSG/ε0.2 | 0.9155 [0.9025, 0.9269] | 0.9195 | 0.9870 | 0.0095 / 0.0035 | 0.9080 | 0.538 | 0.733 |
| | FS | e1stud/effMaxSG/ε0.2 | 0.9700 [0.9616, 0.9766] | 0.9195 | 0.9720 | 0.0040 / 0.0240 | 0.9395 | 0.660 | 0.779 |
| HR 1.75 n 1000 | **DINA** | dinamr/effMaxSG/ε0.2 | 0.9505 [0.9401, 0.9592] | 0.9435 | 0.9770 | 0.0070 / 0.0160 | 0.9440 | 0.689 | 0.848 |
| | FS | cert20/effMaxSG/ε0.2 | 0.9525 [0.9423, 0.9610] | 0.9420 | 0.9740 | 0.0040 / 0.0220 | 0.9485 | 0.769 | 0.860 |

**With the criterion matched, the two engines converge as n grows.** The harm field lower bound
differs by **−0.059** at n 500 (0.916 vs 0.975, Wilson intervals disjoint) but by only **−0.007** at
n 1000 and **−0.001** at n 1500 (Wilson intervals overlapping). The complement field-s upper bound
is **within 0.005 at every cell**. FS identifies the planted region better throughout (sens 0.590 vs
0.497 at n 500, 0.839 vs 0.788 at n 1500).

### The matched p̂ overlay — and it closes with n

| cell | | p̂ T1 | p̂ T2 | p̂ T3 | DINA − FS by tertile index |
|---|---|---|---|---|---|
| HR 1.50 n 500 | DINA bias (med p̂) | −0.140 (0.037) | +0.027 (0.092) | +0.341 (0.230) | **+0.089 / +0.124 / +0.182** |
| | FS bias (med p̂) | −0.229 (0.023) | −0.097 (0.069) | +0.160 (0.187) | |
| HR 1.50 n 1500 | DINA bias (med p̂) | −0.119 (0.035) | −0.058 (0.085) | +0.044 (0.182) | **+0.009 / +0.003 / −0.026** |
| | FS bias (med p̂) | −0.129 (0.038) | −0.061 (0.096) | +0.070 (0.199) | |
| HR 1.75 n 1000 | DINA bias (med p̂) | −0.166 (0.038) | −0.057 (0.095) | +0.118 (0.210) | **+0.018 / +0.023 / +0.003** |
| | FS bias (med p̂) | −0.184 (0.045) | −0.081 (0.109) | +0.116 (0.243) | |

At **n 500** DINA's retained bias sits above FS's in every tertile by **+0.09 to +0.18**. By
**n 1500 the two curves essentially coincide** (+0.009 / +0.003 / −0.026), and their p̂ distributions
coincide too (medians 0.035/0.085/0.182 against 0.038/0.096/0.199) — unlike at 12.4%, where FS's p̂
ran far higher. **At a matched criterion the engines differ at n 500 and are indistinguishable on
this curve by n 1500.**

### The two-sided question, answered at both prevalences

| | n 500 | n 1000 | n 1500 |
|---|---|---|---|
| **12.4%** FS, HR 1.50 | 0.977 | 0.932 | **0.901** |
| **12.4%** DINA, HR 1.50 | 0.9949 | 0.9887 | 0.9925 |
| **31%** FS, HR 1.50 | 0.9810 | 0.9710 | 0.9755 |
| **31%** DINA, HR 1.50 | 0.9905 | 0.9865 | 0.9840 |

**No — the decay does not appear on DINA, and at 31% it does not appear on FS either.** DINA is flat
and high at both prevalences (0.984–0.995). FS decays only at **12.4%** (0.977 → 0.932 → 0.901);
at 31% FS is flat (0.981 → 0.971 → 0.976). **The harm-block two-sided decay is a 12.4%, FS-specific
phenomenon in this grid**, not a property of the IJ two-term interval as such.

**The miss, split by side.** FS's 12.4% decay is **entirely an upper-side miss** — `miss_above`
0.016 → 0.057 → 0.093 at HR 1.50 with `miss_below` ≤ 0.010. DINA's misses stay **≤ 0.010 on both
sides at every cell of both blocks**, shifting from below-dominant at n 500 to above-dominant at
n 1500 without either side growing.

---

## Cells completed, deferred and dropped, with walls

| block | cell | status | wall |
|---|---|---|---|
| A (12.4%) | HR 1.50 n 500 | **complete** | 1840 s (30.7 min) |
| A | HR 1.50 n 1000 | **complete** | 1400 s (23.3 min) |
| A | HR 1.50 n 1500 | **complete** | 1186 s (19.8 min) |
| A | HR 1.75 n 500 | **complete** | 2061 s (34.4 min) |
| A | HR 1.75 n 1000 | **complete** | 1692 s (28.2 min) |
| A | HR 1.75 n 1500 | **complete** | 1543 s (25.7 min) |
| | **Block A total** | **6 / 6** | **9,722 s = 2.701 h** |
| B (31%) | HR 1.50 n 500 | **complete** | 4976 s (82.9 min) |
| B | HR 1.50 n 1000 | **complete** | 5569 s (92.8 min) |
| B | HR 1.50 n 1500 | **complete** | 5857 s (97.6 min) |
| B | HR 1.75 n 500 | **complete** | 5692 s (94.9 min) |
| B | HR 1.75 n 1000 | **complete** | 6610 s (110.2 min) |
| B | HR 1.75 n 1500 | **DEFERRED** (checkpoint, ratified) | — |
| | **Block B total** | **5 / 6** | **28,704 s = 7.973 h** |
| C (null) | HR 1.00, both prevalences, n 500/1000/1500 | **DEFERRED** (pinned) | — |
| | **CAMPAIGN TOTAL** | **11 / 18 cells** | **38,426 s = 10.674 h** |

**10.674 h against the 13 h ceiling — 18% headroom unused.** Nothing was dropped; the seven absent
cells are **deferred** to a follow-up session. Block B realized **7.973 h against 8.694 h
re-projected (8% under)**, against Block A's 27% over-run — so the by-n calibration the checkpoint
applied corrected in the right direction and slightly over-corrected. Every cell carries 2,000
replicates; none was thinned. No cell failed; no gate was stopped. The 16 h hard timeout was never
approached.

## Stage 3 artifact

`summary_dinamr.qmd` → `summary_dinamr.html`, **11 of 18 cells**, rendered clean. The absent-cell
guard works: the seven deferred cells are listed as deferred and **exactly one** chunk skips with a
named note — the null-cell table, which is the only one with nothing to show. The A3 behaviour does
not recur.
