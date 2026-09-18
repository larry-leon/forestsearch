# CC TASK — read-only audit: the admission criterion and the defaults, at HEAD

**Opened:** 2026-09-18
**Repository:** forestsearch (the package).
**Type:** read-only audit. **No edit to `R/`. No fit. No compute. No install. No push.**
**Purpose:** produce one reference picture of how the admission criterion and the defaults are currently set
up, at repository HEAD, serving both the threshold-specification workstream and the crash-diagnostics
workstream so they do not read the same code twice and reach two readings.

**Everything in this document is a question, not a statement.** Where it names an argument, a default, a line
or a behaviour, that name came from an audit of **installed 0.3.5** or from a handoff anchored at `35eca3a`,
and may be wrong or stale. **Answer from current source. Do not confirm anything this document asserts.**

---

## 1. Rules

- Read-only. The only commits are this task document into `dev/tasks/` at the start, and the report at the end.
- **Find things by search, never by line number.** Report the real file and line for everything.
- Record the HEAD SHA and whether the working tree is clean. A dirty tree is acceptable for a read; note it.
- Where a question has no answer in the source, say "not found" and say where you looked. Do not infer.
- **Findings only. Propose no fix, and attach no follow-up task.** If something looks blocking, say so in one
  line and stop there.

## 1a. Do not duplicate work already in the record

- `dev/tasks/TASK_binary_admission_check_2026-09-18.md` is already committed and executed. **Read it and its
  report first.** Where it has already established something this document asks for, cite it rather than
  re-deriving it, and say so. Where this document's question goes further, answer the further part.
- A partial read-only pass was run on 2026-09-18 at `82b421e9` before this document existed. Its findings are
  listed in §1b as claims to confirm or correct, not as established facts. **Extend that pass; do not restart
  it.**

## 1b. Claims from the partial pass, to confirm or correct from source

Each line below is a claim with a cited location. Confirm it, correct it, or report it as not found.

| Claim | Cited at |
|---|---|
| per-arm floor for binary counts events per arm and rejects below `d0.min` / `d1.min` | `R/subgroup_search.R:596-598` |
| per-arm floor for survival does the same via `calculate_event_counts()` / `meets_event_criteria()` | `R/subgroup_search.R:653-654`, `:727`, `:738-739` |
| **the per-arm floor is skipped entirely for continuous and count outcomes** | `R/subgroup_search.R:609` |
| size floor rejects when `nx <= n.min` — a candidate must exceed it, not merely meet it | `R/subgroup_search.R:613` (GLM), `:660` (survival) |
| **`n.min` defaults to 60 in `forestsearch()` but 30 in the `subgroup_search()` helper** | `R/forestsearch_main.R:1246`; `R/subgroup_search.R:83` |
| `d0.min` and `d1.min` both default to 10 | `R/forestsearch_main.R:1266-1267` |
| `sg_focus` defaults to `"hr"` | `R/forestsearch_main.R:1252` |
| `effect.threshold` and `consistency.threshold` both default to NULL and are resolved from `hr.threshold = 1.25` and `hr.consistency = 1.0` at a single site | `R/forestsearch_main.R:1248-1251`, resolution at `:1793-1794` |
| `adverse_outcome` defaults to NULL and resolves TRUE for binary and count only, at two sites | `:1287`, resolved `:1342-1343` and `:1718-1719` |
| ratio measures are converted to log internally; RD, IRD and MD are compared on the identity scale, with a silent remap of the survival default | `R/forestsearch_main.R:1781-1860` |

**Two of these are load-bearing and need more than confirmation:**

- The `n.min` 60-versus-30 split. Report which value a call actually receives at each entry point, and whether
  any committed caller reaches the helper's default rather than the front door's.
- The continuous-and-count skip at `:609`. Report exactly which outcome types bypass the per-arm floor, and
  what admits a candidate on those paths instead.

**Not yet covered by the partial pass, and still required:** `n.min.frac`, `m1.threshold`, `stop_threshold`,
the GRF and DINA argument surface, the estimand defaults per outcome type, and the reference table itself.

---

## 2. Part A — the admission criterion and its floors

The crash-diagnostics workstream reports that an admission floor was not applied, that the same selector code
serves survival, and that the committed DINA and GRF campaigns ran without the criterion. **That report is
not verified here and is not to be assumed.** Establish from source what exists.

1. **Inventory every floor that gates whether a candidate subgroup can be admitted or selected.** For each:
   the argument name, every alias, its default, its unit (a count or a proportion), whether it is applied
   per arm or pooled, whether it counts observations or events, and which `outcome_type` values it applies to.
2. **Distinguish the stages.** Is the same floor applied at candidate generation, at screening, and at the
   consistency stage, or are these separate floors with separate arguments? Report each stage and the floor
   that governs it.
3. **Where each floor is enforced**, and **what happens to a candidate that fails it** — silently dropped,
   warned about, errored on, or retained with a degenerate estimate. This distinction is the important one.
4. **Per-arm behaviour.** Is there any floor on the treatment and control arms separately, as opposed to the
   subgroup as a whole? If yes, on what quantity.
5. **Binary and GLM outcomes specifically.** Is there any guard against a zero cell — no events, or no
   non-events, in an arm? If a candidate with a zero cell reaches the fit, report what the fit returns or
   raises. Read the code path; do not run it.
6. **Identifier forwarding.** Are the floors passed to the `dina` and `grf` paths, do those paths carry their
   own floors under different argument names, or do they apply none? Report per identifier.
7. **Any path on which a floor is bypassed** — a branch that admits candidates without consulting it, or a
   wrapper that does not forward it. If the reported gap exists, this is where it will be.

## 3. Part B — the defaults and their resolution

1. **`effect_measure` resolution.** Every site where it is resolved from `outcome_type`, with the resolved
   value per outcome type. **Report how many sites exist** — a handoff cites two in `forestsearch_main.R`;
   confirm the count from search rather than trusting it.

   The partial pass reported `adverse_outcome` resolved at two sites adjacent to those cited line numbers but
   **never reported `effect_measure`'s own default per outcome type.** That default — specifically what the
   binary path resolves to when `effect_measure` is unset — is unconfirmed at HEAD and is the fact a pending
   behaviour change rests on. Answer it explicitly, with the file and line, and state whether the
   `effect_measure` and `adverse_outcome` resolutions live in the same duplicated blocks or in different ones.

1a. **Duplicated defaults across layers.** For every argument in Parts A and B, report whether the exported
   entry point and any internal helper declare **different** defaults for it, as the partial pass found for
   `n.min`. List every such divergence. This is the same structural defect as a duplicated resolution site and
   is the most likely place for further instances.
2. **The three threshold arguments** — the screening threshold on the subgroup's own effect, the per-split
   consistency threshold, and the proportion of splits that must clear it. For each: name, every alias,
   default, and the scale it is interpreted on, **per outcome type and per resolved estimand.**
3. **Any remapping** of threshold values when the estimand is on an identity scale rather than a ratio scale.
   Report the actual mapped values, from source.
4. **Where each threshold is compared**, on what scale, and whether the inequality is strict.
5. **The consistency-stage entry condition** — the exact condition under which that stage runs, and what is
   returned when it does not run. Report the downstream consequence, not just the branch.
6. **What `fpr_calibration()` and the operating-characteristics functions require** of the two effect
   thresholds, and on what scale they compare.
7. **Identifier invariance.** For a given outcome type, do `consistency`, `dina` and `grf` resolve the same
   estimand and the same screening threshold? Where does each identifier read that threshold from?
8. **Is "not supplied" detectable** for each threshold at the point where a derivation would have to act —
   for every alias spelling, and through any wrapper that passes arguments down explicitly? Report the
   mechanism, or report that it is lost.
9. `git log -S` on the consistency-threshold default and on any `0.8` derivation: did a derivation ever exist
   in this package's history, and on what scale?

## 4. Part C — what the campaign template actually passes

Locate the binary simulation campaign template in this repository and report, **verbatim from the file with
path and line numbers**, whether each of the following is supplied or inherited, and its value if supplied:

`effect_measure`; the screening threshold; the consistency threshold; the consistency proportion; and every
floor argument found in Part A.

The expectation on record is that it passes `OR` explicitly. **Report what the file says, not whether it
matches.**

---

## 5. Deliverable

`REPORT_criterion_and_defaults_audit_2026-09-18.md`, written **beside the repository's existing `REPORT_*`
files**, not in `dev/tasks/`. It opens with the HEAD SHA and the tree state, and contains:

1. **One reference table** — for each argument found in Parts A and B: name, aliases, default, unit, scale,
   which outcome types and identifiers it applies to, where it is resolved, and where it is enforced. This
   table is the deliverable both workstreams will work from; everything else in the report supports it.
2. Part A's answers, with the stage-by-stage floor map and the failed-candidate disposition for each.
3. Part B's answers, with the resolution-site count stated explicitly.
4. Part C's verbatim findings.
5. **A discrepancy list:** every place where current source differs from what this document asserted. That
   list is the most valuable output, because two workstreams are currently reasoning from the installed build.
6. Anything found that was not asked about, recorded as a finding with no task attached.

Commit the report. **Do not push.**
