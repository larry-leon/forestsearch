# CC TASK — documentation only: the threshold naming collision

**Opened:** 2026-09-18
**Repository:** forestsearch (the package).
**Authorized by:** Larry, 2026-09-18.
**Class:** touches `R/`, **documentation only.** No behaviour change, no method change, no code line altered.
**Compute:** none. **Install:** none. **Push:** never.

**Why this is first:** it is the only item on the threshold register that cannot alter a committed result
under any circumstance. There is therefore no baseline capture, no resolution probe, no exposure check and no
fits in this task.

**Evidence base:** `quarto/simulations/actg175/binary_020/REPORT_criterion_and_defaults_audit_2026-09-18.md`
(`ce5f1e84`). Every value that report states is a pointer — **re-read it from source before documenting it.**
Documentation that repeats a stale number is worse than no documentation.

---

## 1. Rules

- Copy this document into `dev/tasks/` and commit it before anything else.
- **Read the current roxygen from source first.** Do not import wording from this document, from any handoff,
  or from any report. Those are pointers; the file is the fact.
- **Every value written into the documentation must be traced to a source line**, and the report cites that
  line. If a value cannot be traced, do not document it — record it as not found.
- Gates are stop-on-failure, not stop-to-ask.

## 2. What must NOT appear in this task

These describe behaviour that does not exist at HEAD. Documenting them here would make the roxygen a
description of a future package.

- **No statement that c2 must be ≤ c1, and no statement that c2 > c1 is an error.** At HEAD it is permitted.
  That sentence belongs to the Directive A task, in the same commit that makes it true.
- **No mention of any derivation of c2 from c1.** The audit established (§3.10) that no such derivation has
  ever existed in the package's history, on any scale.
- **No mention of `OR` as the binary default.** At HEAD it is `RD`.
- **No code change of any kind.** If documenting something correctly seems to require a code change, stop and
  report; do not make it.

---

## 3. Part 1 — the threshold vocabulary (required)

### 3.1 The collision itself

`consistency.threshold` (alias `hr.consistency`) is **c2, an effect threshold** on the per-split subgroup
effect. `pconsistency.threshold` is **p\*, a proportion** — the fraction of splits that must clear c2. Two
different quantities; the word "consistency" appears in both names. This has already caused a
misunderstanding in applied work.

Find **every** roxygen `@param` block for these four arguments, across every exported function that carries
them — at minimum `forestsearch()`, `subgroup.search()`, `subgroup.consistency()` and
`consistency_resample()`, and any other the search turns up — and make each one state, in its own words:

- `effect.threshold` / `hr.threshold` — **c1**, the screening threshold on a candidate subgroup's own effect.
- `consistency.threshold` / `hr.consistency` — **c2, an effect threshold**, applied to the subgroup's effect
  within each split half. Explicitly: *not* a proportion, and *not* the consistency rate.
- `pconsistency.threshold` — **p\***, the proportion of splits that must clear c2. A rate in [0, 1].

Each of the two consistency arguments should name the other and say how they differ, so a reader landing on
either one cannot mistake it for the other.

### 3.2 The resolved defaults, per estimand (required)

This is required rather than optional, because documenting c2's default as `1.0` is **misleading on the GLM
paths** — the value a run actually uses is remapped per estimand. Documenting the argument honestly means
documenting what it resolves to.

Verify each row against source before writing it. The audit locates the resolution block at
`R/forestsearch_main.R:1796-1885`; confirm the range and the values.

| Resolved estimand | Screening (c1) | Per-split (c2) |
|---|---|---|
| `HR` (survival) | 1.25, compared on the natural HR scale by the search | 1.0 |
| `OR`, `RR`, `IRR` | `log(1.25)` | `log(1.0)` = 0 |
| `RD` | 0.05 | 0.0 |
| `IRD` | 0.01 | 0.0 |
| `MD` | 0.0 | 0.0 |

State that `pconsistency.threshold` is **not** remapped — it is a rate, identical for every outcome type.

### 3.3 The `fpr_calibration()` vocabulary bridge (required)

`fpr_calibration()` already names these quantities **`c1`** and **`c2`** (the audit cites
`R/fpr_calibration.R:49-50`; verify). Nothing in its documentation connects those names to
`effect.threshold` and `consistency.threshold`. Add that cross-reference in both directions: `c1` is the
screening threshold, `c2` the per-split threshold, passed through to `forestsearch()` on the natural scale.

This is the cheapest part of the whole task and it closes the gap between the two vocabularies the package
already uses for one pair of quantities.

### 3.4 The two meanings of "the screening threshold" (offered — strike at review if not wanted)

The audit's §3.3 records that on the survival path the search compares on the **natural HR scale**, while on
every GLM path it compares on the **link scale**, and that `threshold_config$screening` holds the *log* value
for survival — the MR and admission scale, not the search's. Both are correct for their consumer, but the
phrase names two different numbers depending on which object a reader inspects.

If kept, document it wherever `threshold_config` is described, as a note on which scale each consumer reads.
**This part is separable: striking it changes nothing in 3.1–3.3.**

---

## 4. Part 2 — the DINA caps roxygen (separable; strike at review if not wanted)

This is Directive C part 1 of `HANDOFF_forestsearch_package_items_2026-09-17.md`, which is documentation only
and independent of everything in Part 1. **Commit it separately from Part 1**, so either can be reverted alone.

Established by `quarto/gbsg/REPORT_dina_family_survey_2026-09-17.md` (`dba3073`) — verify from source before
documenting:

- Rewrite the two `dina_frontier()` caps' roxygen: they **trim the returned frontier table** — single cuts,
  never conjunctions, not the selection family — are **not consulted by `dina_subgroup()`**, and have **no
  effect under `subgroup_method = "dina"`**.
- Replace the phrase describing a cap as "the DINA analog of forestsearch's `max_subgroups_search`" with an
  explicit contrast between the two arguments: one truncates the evaluated pool and defaults `Inf`; the other
  trims a report and defaults finite.
- Extend `m_diff`'s "ignored on that path" sentence to **all seven** frontier keys, naming them.
- Retitle the `details`-time frontier print as a display of **proposed single cuts** shown beside the family
  counts — not as "candidates".

**Part 2 changes no behaviour.** The frontier-key warning and the `Inf` display defaults are Directive C
parts 2 and 3, and are **not** in this task.

---

## 5. Post-conditions (machine-checkable)

1. **The diff touches only documentation.** For every hunk under `R/`, every changed line is a roxygen
   comment line (begins with `#'` after optional whitespace). **No non-comment line in any `R/` file is
   added, removed or modified.** Assert this mechanically over the diff and state the assertion's output in
   the report. If it fails, the task fails.
2. `devtools::document()` has been run; the only other changed paths are under `man/`.
3. `NAMESPACE` is unchanged. If it moves, stop and report why.
4. `devtools::check()` shows **no new** ERROR, WARNING or NOTE relative to a check of the pre-change tree.
   Record the pre-change set; require only that it does not grow.
5. Every numeric value written into the documentation appears in the report with its source file and line.
6. A grep for the forbidden content of §2 over the diff returns nothing: no added text asserting c2 ≤ c1, no
   added text describing a derivation, no added text naming `OR` as the binary default.
7. `NEWS.md` carries one entry under documentation, saying that the threshold arguments' roles and their
   per-estimand resolved defaults are now documented, and that `consistency.threshold` is an effect threshold
   distinct from `pconsistency.threshold`.

---

## 6. Deliverable

`REPORT_threshold_docs_2026-09-18.md`, filed **beside the audit report it follows from**
(`quarto/simulations/actg175/binary_020/`, as the audit noted there is no repository-root `REPORT_*`
location). State the placement and the reason in the report, as the audit did.

The report contains:

1. Every roxygen block changed, with file and line.
2. Every value documented, with the source line it was verified against.
3. The output of post-condition 1's mechanical assertion.
4. The pre-change and post-change `check()` finding sets.
5. Whether §3.4 and Part 2 were kept or struck, and by whose decision.
6. Anything found while reading the roxygen that is **wrong at HEAD** rather than merely unclear — recorded as
   a finding, with no fix applied and no task attached.

Commits: this task document; Part 1; Part 2 (if kept). **Do not push.**

---

## 7. Out of scope

- Directive A's validation and derivation, and the c2 ≤ c1 sentence that goes with them.
- Directive B, the binary default estimand.
- Directive C parts 2 and 3 — the frontier-key warning and the `Inf` display defaults.
- The `args_call_all` sync defect (audit §3.9).
- Every other finding in the audit's §6. They are recorded there; none is documented here.
