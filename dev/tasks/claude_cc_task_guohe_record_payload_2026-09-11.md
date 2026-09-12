# TASK — Guo & He record: factual corrections and the manuscript payload

**File:** `dev/tasks/claude_cc_task_guohe_record_payload_2026-09-11.md` · **Issued:** 2026-09-11 by chat, commissioned by Larry
**Machine:** either (Linux or Mac Studio) · **Repo:** `larry-leon/forestsearch` · **Branch:** `feature/glm-extension`, HEAD at or after `1d9401cb` (Larry pulls in GitHub Desktop before kickoff)
**Transport:** four files `~/Downloads` → `dev/tasks/`, verified and committed first (§1).
**Companion files (travel with this document):**

| file | SHA-256 |
|---|---|
| `claude_cc_task_guohe_record_payload_2026-09-11_edits.json` | `8c0c09b8f15c7b2f1a137d22f9374893eafb53598e0c09c64239f89c1c850b82` |
| `claude_cc_task_guohe_record_payload_2026-09-11_apply.py` | `3913c2017da21ab9e405d1a22e8d9b08d502c3832a7f4b58689fd2c2ad05edf2` |
| `claude_cc_task_guohe_record_payload_2026-09-11_check.R` | `0789fccb2046d2280235287f5aa6dc372b92728a037d8f0c83a21cf8716957a2` |

**Purpose.** `quarto/GuoHe/guohe_supp_section.qmd` was read against its sources on 2026-09-11 (Guo & He 2021; the committed records). It carries factual errors, and the manuscript section — now drafted separately, in the `fs-glms-interpretable` project — needs summaries no record holds yet. This task (1) corrects the record by exact-string edits without restructuring it, (2) adds one section that computes the manuscript payload from the committed bundles, and (3) renders, gates, commits and reports. The record stays the analysis record; the manuscript text is not written here.

---

## Standing conventions (govern this session)

1. **Unattended.** Gates are stop-on-failure: on a failed gate, commit what is green, write the record with the failure stated, stop. Never ask, never work around. STOP only for a result that cannot be trusted (input mismatch, a failed exact-edit or identity gate, anything requiring an `R/` change). Documentation gaps — including render failures and render warnings — go to OPEN ITEMS verbatim, and the run continues with every step that does not depend on the failed output.
2. **Exact single-string matching only** on the record: every edit is applied by the shipped applier, which checks the input hash, requires each old text to occur exactly once (or, for the two spans, both anchors exactly once and the span's SHA-256), and checks the output hash. Never hand-edit the record; never use line offsets.
3. **No `R/` change of any kind. No compute:** no simulation, no bootstrap, no re-run. Every number comes from committed bundles.
4. **No network git** (no fetch, pull or push). Every `git add` by explicit path. Larry pushes.
5. **Record numbers are computed, never typed:** paste script output verbatim.
6. **Protected text.** The two fixed blockquote sentences in the record's B3 callout are protected. Sentence 1 changes by one word under D1 ("hazard-ratio scale" → "log-hazard-ratio scale"); sentence 2 is untouched. The applier verifies both.
7. **Verify from source.** Stage 0 quotes by content from the current tree, not by line offset.
8. **Tracking conventions:** the rendered HTML stays untracked (`quarto/.gitignore`: `GuoHe/*.html`); the payload is tracked, following the `_payloads/<stem>/<stem>_payload.rds` pattern of `quarto/applications/actg175/analysis_actg175_continuous_oc.qmd`.

---

## ⚠ CATEGORY

**No `R/` change.** One analysis record edited (16 exact edits), one payload written, one report. **Compute: none.** Two renders of the record (baseline and edited), each reading committed bundles only; expected minutes, not hours. No compute go/no-go is required.

---

## Decisions (Larry, 2026-09-11)

- **D1** — fixed sentence 1 reads "log-hazard-ratio scale" (Guo & He's Table 8 is captioned on the log hazard ratio). Approved.
- **D2** — the E1 range citations leave the record (B6 limit 2 and the setup's `dev/notes` read); B6 names campaign `e1stud` only as where the screened identifier is evaluated, with its criterion outcome.
- **D3** — B6 limit 1 carries the theory-scope sentence (Guo & He's own continuum design; near-duplicate candidates strain condition (A5)).
- **Scope** — the manuscript section is drafted in `fs-glms-interpretable` from a brief pinned to this task's final commit; nothing here drafts it.

**What the 16 edits do** (full text in the edits file; `purpose` field per edit):

| id | where | change |
|---|---|---|
| E01 | YAML | title and subtitle: analysis record; no "data-built" for this family |
| E02 | setup comment | states the record's role and names this task |
| E03 | setup (span) | removes the E1 range-citation block and its `dev/notes` read |
| E04–E06 | opening | "data-built" removed; says what the order-statistic family is |
| E07 | B2 callout | Guo & He's conservativeness as they attribute it (their Section 2.5, small r), replacing an unsupported mechanism |
| E08 | B3, fixed sentence 1 | D1 |
| E09 | B3 | the complement-truth explanation now renders (it sat only in a never-taken branch); the structural argument is scoped to this nested family |
| E10–E11 | B5 table/caption | specificity and PPV shown as "—" at β₂ = 0 (no harm region); mean γ_ĉ column added |
| E12 | B5 callout | identification framing; the null cell as the argmax's reference; dilution of γ_ĉ |
| E13 | B6 (span) | limit 1 grounded on this design's records; limit 2 without E1 ranges; limit 3 acknowledges Guo & He Section 6; limit 5 added; the cert-note chunk removed |
| E14 | before Provenance | new section "Manuscript payload" (`export-payload` chunk) |
| E15–E16 | provenance table | the certified-range row becomes the payload row |

---

## 1. Provenance and first commit — GATE

```bash
cd <repo root>
hostname; pwd; git branch --show-current; git rev-parse HEAD; git status --porcelain; git log --oneline -3
git merge-base --is-ancestor 1d9401cb HEAD && echo "HEAD contains 1d9401cb"
Rscript -e 'cat(as.character(packageVersion("forestsearch")), R.version.string, "\n")'
```

*GATE:* branch `feature/glm-extension`; HEAD contains `1d9401cb` (if not: STOP — Larry pulls first); no dirt in this task's paths (`quarto/GuoHe/`, `dev/tasks/`). Pre-existing untracked files elsewhere are left alone.

Copy the four files from `~/Downloads` into `dev/tasks/`; verify each companion's SHA-256 against the table above (`shasum -a 256` on macOS, `sha256sum` on Linux) — any mismatch: STOP. Commit the four files alone, by explicit path. Record `HEAD` before this commit as **BASE** for §6.

## 2. Stage 0 — preconditions and source quotes — GATE 0

2.1 *Input hash:* `shasum -a 256 quarto/GuoHe/guohe_supp_section.qmd` must equal `0058800765baa846721e562b59756be2cf875fd9b2e23d319cd65e0b00f12a4a` (git blob `4afa7e42`, the text the edits were written against). Mismatch: STOP.

2.2 *Inputs present* (all under `quarto/GuoHe/`): `mr_field_vs_guohe_{t35_beta2_00..05, t6_k02, t6_k06, t6_k10, t6_k12, t7_beta2_00..05}.rds` (16), `guohe_repro_t7_beta2_00..05.rds` (6), `mr_field_complement_vs_guohe_t7_beta2_00..05.rds` (6), `guohe_adaptive_t7_beta2_0{0,3}_fixedr00833.rds` (2). Any absent: STOP.

2.3 *Quotes by content* (these define two payload rows; record each with its current line number):
- `R/fs_mr_inference.R`: `se_wald    <- sdv[sel]`; the harm naive interval `lower = to_eff(beta_naive - z975 * se_wald)`; the complement naive interval `lower = to_eff(bnc - z975 * sec)`; the complement debiased block `se = sec_used, se_ij = se_ijc_rep$se, se_wald = sec,`.
- `quarto/GuoHe/mr_field_complement_vs_guohe_run.R`: `cgate_naive_est = .lg(gc$naive$est)`, `cgate_est = .lg(gc$debiased$est)`, `cgate_se_ij = .nz(gc$debiased$se_ij)`, `cgate_se_wald = .nz(gc$debiased$se_wald)`.

If any quote is absent, record it in OPEN ITEMS and mark the payload's `comp_t7` rows `naive` and `ij` as **unconfirmed** in the report; continue.

2.4 *Quarto:* the record's setup comment says to render with RStudio's bundled Quarto binary on the executing machine; record its path and version. Render in a UTF-8 locale.

## 3. Stage 1 — baseline render — GATE 1

Render the unedited record: `quarto render quarto/GuoHe/guohe_supp_section.qmd` with the §2.4 binary. Keep the log and record the exit status and any `WARNING` lines. A failure or warning here concerns the environment, not the edits: OPEN ITEM, continue. The baseline HTML is untracked; nothing is committed in this stage.

## 4. Stage 2 — apply the edits — GATE 2

```bash
python3 dev/tasks/claude_cc_task_guohe_record_payload_2026-09-11_apply.py \
        dev/tasks/claude_cc_task_guohe_record_payload_2026-09-11_edits.json --check-only
python3 dev/tasks/claude_cc_task_guohe_record_payload_2026-09-11_apply.py \
        dev/tasks/claude_cc_task_guohe_record_payload_2026-09-11_edits.json
```

Paste both outputs in full. *GATE 2:* both exit 0. The applier has then verified: the input hash; all 16 edits applied by exact matching; each fixed sentence present exactly once (sentence 1 in its D1 form); `{#sec-b3}` once; zero occurrences of the retired terms and reads (`stable-pick`, `data-built`, `CERT_`, `cert_upper`, `cert_joint`, `The wider figures`, the `dev/notes` path, the old B5 headline, and case-insensitive `certif`); every chunk label once; the chunks `b2-build`, `b2-table`, `b2-adaptive-note`, `b2-consolidated`, `b3-table`, `b5-chat-plot`, `provenance-platform` byte-identical to the committed text (so every B2 and B3 table renders unchanged from unchanged inputs); and the edited file's SHA-256 equal to `15abd1209a4ea42432389a6e8d8c7ece1972c4a39f7f4632ac0d65a7b33c6866`. Any failure: STOP (the file is not written on failure).

## 5. Stage 3 — render and identity — GATE 3

5.1 Render the edited record as in §3. Record the exit status, any `WARNING` lines, whether the payload exists at `quarto/GuoHe/_payloads/guohe_supp_section/guohe_supp_section_payload.rds`, and whether the rendered B6 contains "Five limits belong on the record"; anything short of a clean render goes to OPEN ITEMS. The edited text is hash-verified by §4, so a render failure is an environment matter: the `.qmd` is committed regardless (§5.3). B6 carries two `stopifnot()` guards on the statements it makes (field bias positive in every cell; the complement's low point at β₂ = 0), so a false statement fails the render rather than printing; if either fires, quote the error in OPEN ITEMS. If the render stopped before writing the payload, skip §5.2 and the payload commit.

5.2 From the repo root: `Rscript dev/tasks/claude_cc_task_guohe_record_payload_2026-09-11_check.R > /tmp/guohe_payload_check.md`. *GATE 3b:* exit 0 and first line `GATE I: 186 comparisons, 0 outside tolerance`. The script compares the payload with the committed records wherever they overlap: `REPORT_mr_field_vs_guohe_2026-09-05.md` lines 12–27 (the 16-cell harm side at r = 1/30: field, IJ, G&H) and `REPORT_mr_complement_vs_guohe_2026-09-09.md` lines 167–172 and 183–188 (complement and joint). Tolerance is 0.6 units in the record's last printed digit: the 2026-09-05 record transcribed its headline table from 3–4-decimal displays and double-rounds six entries; every other value matches exactly. Any comparison outside tolerance: STOP, and the payload is not committed.

5.3 `shasum -a 256` the payload (record the value). Commit by explicit path: the `.qmd` (always, once §4 has passed) and the payload `.rds` (only if GATE 3b passed). Not the HTML.

## 6. Stage 4 — the record

`quarto/GuoHe/REPORT_guohe_record_payload_2026-09-11.md`, beside the other `REPORT_*` files, containing:

- **Provenance:** machine, BASE and final HEAD, R and forestsearch versions, the Quarto binary path and version, render walls and exit codes, `WARNING` counts.
- **Stage 0:** the §2.3 quotes with line numbers.
- **Stage 2:** both applier outputs, verbatim.
- **Stage 3:** the GATE I line, and then the markdown tables from `/tmp/guohe_payload_check.md` pasted verbatim (every payload element, the meta block, the input MD5s); the payload's SHA-256.
- **OPEN ITEMS**, then `git log --oneline BASE..HEAD`.

No interpretation: the tables are the record. Findings go in OPEN ITEMS; no task is proposed unless something blocks. Commit the report alone.

## 7. Out of scope

No `R/` change; no simulation; no restructuring of the record (block order, headings and every B2/B3 table unchanged); no edit to any other file beyond the four task files, the record, its payload and the report; no touch of `feature/glm-extension-mac`; no push; no manuscript text.
