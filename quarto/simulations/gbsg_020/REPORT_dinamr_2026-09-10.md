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
