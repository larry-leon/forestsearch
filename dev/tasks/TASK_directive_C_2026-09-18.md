# CC TASK — Directive C: the DINA guard, the frontier-key warning, honest displays

**Opened:** 2026-09-18 · **Repository:** forestsearch · **Authorized by:** Larry, 2026-09-18.
**Order fixed by Larry: the DINA guard lands first.**
**Evidence base:** the audit, sync and Directive A reports in `dev/reports/`; the docs task's report
(`REPORT_threshold_docs_2026-09-18.md`) for what is already documented. Pointers — re-verify from source.
**Compute:** at most two seeded micro-fits inside the acceptance tests (one `dina`, reused across parts),
hard abort 5 minutes. Everything else is resolution, search, or string assertions.
**Testing:** ONLY this task's own acceptance test files, hard abort 3 minutes. No R CMD check, no full suite.

## Parts, in order

| # | Part | Class |
|---|---|---|
| 1 | **The DINA guard** — error on identity-scale estimands where the floor derivation is meaningless | changes behaviour (new error on an opt-in combination) |
| 2 | Frontier-key warning under `subgroup_method = "dina"` | changes behaviour (messaging) |
| 3 | `dina_frontier()` display caps → `Inf` / `Inf`, warn when a finite value trims | changes behaviour (display object only) |
| 4 | Retitle the `details`-time frontier print as proposed single cuts (the docs task stopped on this — it is a string in code, legal here) | messaging |
| 5 | c2 / p\* echo annotation under `dina` / `grf` | display |

Each part is its own commit; any part can be reverted alone. **No selection logic is touched anywhere.**

## Part 1 — the DINA guard

**The defect (audit §3.8, F3):** DINA's floor derives as `m_diff <- if (family == "gaussian") hr.threshold
else log(hr.threshold)` — correct for ratio thresholds, meaningless for identity-scale ones. Binary with
explicit `effect_measure = "RD"` and a threshold of 0.07 yields a log-OR floor of `log(0.07) = −2.66`,
admitting nearly everything, silently.

**The fix is a refusal, not a derivation** — there is no RD-scale number DINA could correctly apply, because
it evaluates candidates on the family link scale. `stop()` when the resolved estimand is identity-scale
(`RD`, `IRD`) on a non-gaussian family and the DINA floor derivation would run, with a message of the shape:

```
subgroup_method = "dina" does not support identity-scale estimands (RD, IRD):
DINA's admission floor operates on the family link scale. Use a ratio estimand
(e.g. effect_measure = "OR") with dina, or subgroup_method = "consistency" or
"grf" for RD.
```

- **Place the guard at the derivation site**, not at the `subgroup_method` dispatch, so every route to the
  derivation is covered — including `use_dina` screening if it reaches the same derivation with an
  identity estimand. Establish the routes from source first and list them in the report.
- **Gaussian is untouched:** `MD` on the gaussian family is a correct identity-scale floor and `mddina` is a
  committed campaign — the guard must be provably unable to fire there (acceptance test).
- **Exposure gate:** confirm by search that no committed in-repo caller reaches DINA (any route) with an
  identity-scale estimand on a non-gaussian family. Binary and count DINA have zero committed campaign
  cells; verify rather than inherit. **If a committed caller is found, STOP and present pin-vs-proceed**,
  as the sync task did.

## Part 2 — the frontier-key warning

A `warning()` when `dina_args` carries any frontier key under `subgroup_method = "dina"`, naming the ignored
keys. The seven keys and their inertness are already established and documented (docs task: all seven inert
under `subgroup_method = "dina"`; under `use_dina` + `selected_only = TRUE` six inert, `digits` still acts) —
the warning is the runtime counterpart on the `subgroup_method = "dina"` path only. One warning per fit,
listing all offending keys at once, never one warning per key.

## Part 3 — honest display caps

`dina_frontier()`: `max_subgroups` and `max_per_covariate` default `Inf` / `Inf`; a user-supplied finite
value that actually trims warns, naming both counts (kept vs. available). Aligns with
`max_subgroups_search`'s own `Inf` doctrine; makes the frontier display truthful (a finite pooled cap drops
covariates whole and truncates the rank-2 pass — observed in applied vignettes).

- **Exposure gate:** inventory every in-repo caller of `dina_frontier()` that inherits the finite caps. The
  frontier table is a display object — `dina_subgroup()` never consults it — but if any committed rendered
  artifact's tables would change on re-render, STOP and present pin-vs-proceed with the list. Update the
  roxygen the docs task wrote for the caps (defaults change from `3L` / `10L`).

## Part 4 — the frontier print retitle

The `details`-time frontier print titles the table as candidates. Retitle as a display of **proposed single
cuts** shown beside the family counts. The docs task located the strings (`forestsearch_helpers.R` and
`forestsearch_main.R` per its report) — re-find by search. String change only; no structural change to what
is printed.

## Part 5 — the c2 / p\* echo annotation

Wherever c2 (`consistency.threshold` / `hr.consistency`) and p\* (`pconsistency.threshold`) are echoed in
output under `subgroup_method = "dina"` or `"grf"` — print/summary methods, fit banners, stored config
displays — annotate them as **not used on this path**. Find every echo site from source before annotating
any; list them in the report.

**Motivating case, cite it:** `quarto/simulations/gbsg_020/summary_grfmr.qmd` runs `hr.threshold = 0.90`
against the inert default c2 = 1.0 under `grf`, and nothing in its output says c2 was never consulted
(Directive A report's finding). This is display only: no warning, no error, no value changes — the sync/A
probe cells for dina/grf must remain byte-identical.

## Rules

Copy this document into `dev/tasks/` and commit it first. No install; `devtools::load_all()` only. Do not
push. Gates stop-on-failure (the two exposure gates above present options instead, as the sync did). Find
everything by search. Report ≤ ~120 lines, file + function citations, to
`dev/reports/REPORT_directive_C_2026-09-18.md`.

## Baseline and gates

Before any edit, extend the existing probe (one copy, `helper-threshold-sync.R`) or record equivalent
baseline evidence for: dina cells under each estimand (OR, RR, IRR, MD-gaussian pass; RD, IRD will error
after Part 1), grf cells unchanged everywhere, and the consistency-path cells byte-identical throughout —
Parts 1–5 must not move any consistency-path resolution.

After the edits:

- **Gate A:** the only new failures are dina + identity-scale on non-gaussian; the only new warnings are the
  frontier-key and cap-trim warnings, each firing exactly where specified and nowhere else.
- **Gate B:** `mddina`-shaped calls (gaussian, MD, dina) pass silently, byte-identical resolution.
- **Gate C:** every consistency-path and every grf resolution cell byte-identical to baseline.

## Acceptance tests (one file, 3-minute abort)

1. The guard fires for dina + RD and dina + IRD (non-gaussian), message as specified, **before any model
   runs**; it cannot fire for gaussian MD, or for OR/RR/IRR, or under `consistency` / `grf`.
2. Every route to the m_diff derivation found in Part 1 is covered by the guard (route list asserted).
3. Frontier keys under `subgroup_method = "dina"` produce one warning naming every offending key; no keys,
   no warning; the same args under `consistency` warn nothing.
4. `dina_frontier()` defaults are `Inf` / `Inf`; a finite trimming value warns with both counts; a finite
   non-trimming value is silent.
5. The retitled print string appears; the old title does not.
6. The echo annotation appears under dina and grf and does NOT appear under consistency.
7. The micro-fit(s) stay under the cap; wall clock recorded.

## NEWS and report

`NEWS.md`: entries for Parts 1–3 (who is affected; for Part 1: only explicit identity-scale DINA requests,
which nothing committed uses). Report: gates, the route list, both exposure inventories, echo-site list,
findings with no tasks attached. Commit. **Do not push.**

## Out of scope

The floors and the per-arm criterion (the diagnostics chat's D and E). The `use_dina` unoriented-floor
residual (md workstream's record) — report if observed, change nothing. RD/IRD/MD/IRR behaviour outside the
guard. fs-glms-interpretable.
