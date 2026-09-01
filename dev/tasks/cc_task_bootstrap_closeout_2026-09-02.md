# CC task — bootstrap-efficiency close-out: two verifications, production-shape confirmation, and the register stamped

**File:** `dev/tasks/cc_task_bootstrap_closeout_2026-09-02.md` · **Issued:** 2026-09-02 by chat, commissioned by Larry
**Repo:** `~/Documents/GitHub/forestsearch` (pop-os) · **Branch:** `feature/glm-extension`
**Transport:** `~/Downloads` → `dev/tasks/`, committed first (§1). Hyphen-stripped stem `cc_task_bootstrap_closeout_20260902.md` is the same file.
**Reads (in-repo):** `dev/glm-continuous-sims/REPORT_assembly_skip_2026-09-02.md` (Settlements A and B — the seeding mechanics and the soundness verdict this task builds on), `assembly_settleA_2026-09-02.R` (the transplant source for §4), `medians_baseline_2026-09-02.R` / `assembly_battery_2026-09-02.R` (fixture blocks and the pruner, reused verbatim), `REPORT_bootstrap_reprofile_2026-09-01.md` and `REPORT_bootstrap_profile_2026-08-30.md` (the bucket lines cited in §5's closures).

**Why.** The implementation phase is complete (0.3.1 → 0.3.5, every change bit-identical). What remains before the workstream can close comprehensively: two verifications that gate decisions Larry has not yet made, one production-shape confirmation that turns Settlement A into campaign guidance, and formal closure lines for three memo items the profiles killed but no report ever stamped. **This task changes nothing** — no `R/` edit, no comment edit, no version bump, no defaults. Findings go in the record; the decisions they feed return to Larry.

---

## ⚠ CATEGORY — read-only

Writes: this task doc, the §3/§4 scripts (transplanted or fixture-block-verbatim), their outputs and logs, and the report — all under `dev/tasks/` and `dev/glm-continuous-sims/`. **Zero edits to `R/`, `DESCRIPTION`, `NEWS.md`, tests, drivers, applications, or documents.** The stale doc line at `bootstrap_analysis_dofuture.R:200` recorded at Settlement A stays untouched — it belongs to the queued pre-CRAN housekeeping pass; say so in the report's register.

**Compute:** §3 direct-call pair ≈15 s; §4 sequential nb = 100 continuous bootstrap ≈4 min + its 48-worker twin ≈1 min + one survival 48-worker nb = 100 wall ≈1 min; §2 and §5 are source reads and citations. **Estimate 10–15 minutes wall.**

**Unattended.** Gates stop the task; never ask, never work around. Provenance gate dirt-tolerant; stop only if dirt touches this task's own paths.

---

## 1. Provenance and first commit — GATE

Standard block; `packageVersion("forestsearch")` expect **0.3.5**. *GATE:* branch `feature/glm-extension`; the assembly-skip commits `4b826c98`, `934f31ed`, `f3975b99`, `e31e612d` in the log. Copy this document into `dev/tasks/`, commit **alone**, report filenames and hash. Record the `vi.grf.min` default in force.

## 2. Verification A — what does the replicate pipeline consume, and is `stop_threshold` forwardable?

This feeds Larry's early-stop decision. Source quotes only; no runs.

1. **Consumption.** From `bootstrap_analysis_dofuture.R` (and whatever it calls to build the per-replicate payload), quote every read of the inner `forestsearch()` result inside a replicate: which components of `grp.consistency` (or the full candidate/consistency table) reach the stored payload or any downstream computation — H/Hc estimation inputs, counters, diagnostic columns, anything. Then state the verdict in one of two forms: **(i)** the replicate consumes only winner-level content (selected subgroup definition/membership, its estimates, H/Hc fits) plus named bookkeeping fields — list those bookkeeping fields exactly, because they are what would legitimately differ under early stop; or **(ii)** something consumes the full evaluation (a per-candidate Pcons vector, the complete consistency table, a count used analytically) — quote it, because that would make early stop a content change, not a bookkeeping change.
2. **Forwardability.** Quote where `args_FS_boot` is assembled and state whether a caller-supplied `stop_threshold` rides through to the inner `forestsearch()` call unmodified. If it is forced, stripped, or defaulted by the wrapper, quote the line — enablement would then require a wrapper change (an `R/` edit for a future task), and the report must say so.
3. **Branch caveat, stated in the report:** this verdict holds for `feature/glm-extension`'s wrapper. The `feature/mr-in-replicates` branch adds in-replicate machinery with its own consumption; the early-stop decision must be re-checked at merge time. One recorded line; no work on that branch here.

## 3. Verification B — direct-call plan reproducibility (settles the seeding-unification item by measurement)

The recorded item asks whether the plain loop (sequential plan) and the batched path (parallel plan) produce identical results for the same seed on a **direct** `forestsearch()` call — replicates always take the plain loop, so this is direct-call reproducibility hygiene, not bootstrap speed. Settle it:

1. One script, fixture block **verbatim** from the committed battery (F1 continuous, `seedit = 8316952` as the drivers set), pruner sourced or copied verbatim from the battery: run F1 once under `parallel_args = list(plan = "sequential", workers = 1L)` and once under `list(plan = "multisession", workers = 6L)`; compare pruned results with `identical()`.
2. Verdict for the record: **(i)** identical — the unification item is moot, close it; or **(ii)** different — report the first differing component and values verbatim, quote the two RNG-consumption sites that diverge (the plain loop's draw pattern vs `future_lapply`'s `future.seed` streams in the consistency resampling), and the item remains an open decision for Larry with the mechanism named. Either way, one gbsg repeat (F2, both plans) so the verdict is stated per outcome type. No fix, no unification, no default change here.

## 4. Production-shape confirmation and campaign guidance

Settlement A proved reproducibility at nb = 40 / 20 workers. Turn it into guidance at production shape:

1. **Transplant** `assembly_settleA_2026-09-02.R` → `closeout_parallel_2026-09-02.R`, changing only the named lines: `nb_boots = 100`, the worker count (48), and the output filenames. *GATE:* show the diff; anything beyond those named lines is a stop.
2. Run the continuous F1 bootstrap at nb = 100: once sequential, once `multisession, workers = 48`. *GATE:* payloads `identical()` under the pruner (structurally guaranteed by the pre-generated index matrix and per-replicate streams — verify anyway; if they differ, that is a major finding: report first difference verbatim and stop). Report both walls, the per-replicate throughput, startup amortization (wall − nb × per-rep), and the implied B = 1000 wall at 48 workers.
3. One survival nb = 100 run at `multisession, workers = 48` — wall only (its reproducibility follows from the same mechanism plus Settlement A; say so rather than paying the ~7-minute sequential comparator). Implied survival B = 1000 wall at 48 workers.
4. In the report, a short guidance block for campaign launch blocks: the outer-plan invocation as the wrapper's source shows it should be written, the measured throughputs, the startup-amortization note, the reminder that inner searches stay sequential by construction, and the observed worker memory footprint (report `free`/RSS before and during the 48-worker run so campaigns on this machine can size worker counts). **Guidance in the record only — no driver or document is edited.**

## 5. The register, stamped — citations, no compute

One report section closing, with the exact committed lines cited (report + section, quoted figures):

1. **Memo item E** (combination-index hoisting): dead by measurement — enumeration/floors bucket ≤0.01 s per call across all profiled fixtures (08-30 and 09-01 bucket tables). Closed.
2. **Memo item F** (membership/floor checks): dead by measurement — its targets sit inside buckets totalling ≤0.2% per call at every profile through 0.3.5. Closed.
3. **Memo §4** (the `"split"` consistency fallback inside replicates): fired 0 of 1,616 (08-30) and 0 again (09-01). Closed.
4. **Memo item A, binary limb:** deferred by design to the OR workstream (IRLS bit-identity uncertifiable; glm.fit-level routing is that track's call). Recorded as handed off, not open here.
5. **Deep Cox surgery:** recommend-closed as mooted by reproducible outer parallelism; the residual per-call coxph internals are wall-clock-covered by §4's campaign configuration. (Larry may reopen; the record states the recommendation and the residual seconds.)
6. **Stale doc line** `bootstrap_analysis_dofuture.R:200`: assigned to the pre-CRAN housekeeping pass.

## 6. Report

`dev/glm-continuous-sims/REPORT_bootstrap_closeout_2026-09-02.md`:

1. Provenance.
2. §2's consumption quotes and verdict — including, under verdict (i), the exact named bookkeeping fields that would differ under early stop, and the forwardability finding; the mr-in-replicates caveat line.
3. §3's plan-reproducibility verdict per outcome type, with quotes or first-difference as the outcome dictates.
4. §4's walls, throughputs, amortization, memory footprint, implied B = 1000 lines, and the guidance block.
5. §5's register closures with citations.
6. **The decision hand-back:** a short section addressed to the workstream, stating the two decisions now sitting with Larry (early-stop enablement — contingent on §2's verdict and field list; medians removal — hygiene, ~0.22 s) and what each would require if commissioned (per §2.2: driver-arg-only or wrapper edit; per the medians record: NA-out + display recompute, separate task). No recommendation beyond what the findings state; the recommendation layer is chat's.
7. Ten-line verdict.

Commits, in order: (1) this task doc alone; (2) the §3 script and the §4 transplant (diff shown); (3) all outputs, logs, `.rds`, and the report. Explicit paths; tree clean at close. **No push. No render. No edit outside `dev/`.**
