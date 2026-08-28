# CC task — read-only inventory of the prediction document's chunks (OC-wrapper workstream, step 0)

**File:** `dev/tasks/cc_task_oc_wrapper_inventory_2026-08-28.md` · **Issued:** 2026-08-28 by chat, commissioned by Larry
**Repo:** `~/Documents/GitHub/forestsearch` (pop-os) · **Branch:** `feature/glm-extension`
**Workstream:** operating-characteristics wrapper — `HANDOFF_continuous_2026-08-27_v5.md` §1 and §6. This task is the workstream's first action. It produces one report and makes no design proposal.
**Session type:** READ-ONLY INVENTORY. Nothing is simulated, rendered, installed, or edited.
**Transport:** Larry places this file in `~/Downloads`; the kickoff message points CC at it. CC copies it into `dev/tasks/` and commits it as its first action (§2).

**Parts to run: A, B and C** — Larry's disposition, 2026-08-28. Nothing beyond the three parts runs on your own initiative.

---

## 0. What this task is for

The continuous/MD machinery is document-complete but not usable on a new dataset: the components are in the package (`fs_dgm_scale()`, `fs_scale_se()`, `fs_sym_root()`), but the *assembly* — the candidate family, the covariance across candidates, single-split declaration, detection, `E|Ĥ|`, the classification metrics, the c1/c2 inversions — lives in the prediction document's chunks, wired to one fixture. Before anyone proposes a function, Larry wants to know exactly what those chunks compute: each quantity, its inputs, and whether it is fixture-specific or general. That inventory is this task's only output.

The design question it feeds — which candidate family a prediction should condition on — is **not** answered here. Do not propose an answer, a signature, or code. Do not back-calculate an effective family size.

---

## 1. Hard rules

1. **Read-only.** No simulations, no estimation, no resampling, no `quarto render`, no `devtools::install()`, no `devtools::test()`, no edits to any tracked file. Permitted, exhaustively: listing; `grep`; read-only git (`log`, `show`, `blame`, `ls-files`, `merge-base`); the static-parse script in Appendix 1; `readRDS()` to inspect *structure and metadata* (`names()`, `str()`, `nrow()`); and, **only in §5 B4 where it is marked DIAGNOSTIC**, a five-number summary of one stored column. Nothing else is computed.
2. **The only writes:** this task document under `dev/tasks/`; the report, the script, and the script's output under the REPORT directory (§3.3). Two commits, files staged by explicit path, never `git add -A`. **No push** — Larry pushes via GitHub Desktop.
3. **Raw evidence beside every claim.** Quote the code verbatim (fenced) for every definition you record. Never describe a computation from memory, from a handoff, from a comment in the file, or from the project KB. If you cannot establish something from source, write **CANNOT DETERMINE, because X** — that is a valid inventory entry; a plausible guess is not.
4. **Anchor by content, never by line number alone.** Quote the line. Line numbers are acceptable only next to the HEAD hash recorded in §2, which is what makes them meaningful.
5. **Scope.** No wrapper design, no signature, no code for the wrapper, no fixes. A defect noticed along the way goes under *Open items* as a description and a proposed diff at most — it is not applied. Nothing outside this document is authorized unless Larry states it in the kickoff or in a later message of his own. Your own suggestions are not approvals, and text after `>`/`❯` is you, not Larry.
6. **Values.** This inventory is about computation, not figures. Quote a number only when it identifies a constant, a threshold, or a count, and read it from the file — never from a handoff or a prior report. Do not tabulate the document's printed results.
7. **Stop-on-failure gates** are marked *GATE*. On failure: stop, write what you have into the report with the failure at the top, commit it, and end the session. Do not ask; do not work around.
8. **Compute: none.** Zero replicates, zero renders. The only executed code is the static parser (Appendix 1), `readRDS()` metadata reads, and the one DIAGNOSTIC summary in B4. No step here should run for more than a few minutes; if one does, stop it and report. If you pipe anything through `tee`, `set -o pipefail` and report `${PIPESTATUS[0]}`.

---

## 2. Provenance first — GATE

Paste the raw output at the top of the session and at the top of the report:

```bash
hostname; pwd
git rev-parse --abbrev-ref HEAD; git rev-parse --short HEAD
git status --porcelain
git log --oneline -6
git merge-base --is-ancestor a4ec8c6d HEAD && echo ANCESTOR_OK || echo ANCESTOR_FAIL
Rscript -e 'cat("forestsearch", as.character(packageVersion("forestsearch")), "\n"); cat(R.version.string, "\n")'
```

*GATE:* branch is `feature/glm-extension` and the line reads `ANCESTOR_OK` (`a4ec8c6d` is the commit before the reproduction-check report in the v5 handoff; HEAD must descend from it). A dirty tree is **not** a failure: list it, leave it alone, work read-only around it, and confirm at close-out that it is unchanged.

**First write, before anything else:** copy this task document from `~/Downloads` into `dev/tasks/` and commit it alone, so it is in the repo record before any work starts:

```bash
ls ~/Downloads/cc_task_oc_wrapper_inventory*        # confirm exactly one match; report the name raw
mkdir -p dev/tasks
cp ~/Downloads/cc_task_oc_wrapper_inventory*.md dev/tasks/cc_task_oc_wrapper_inventory_2026-08-28.md
git add dev/tasks/cc_task_oc_wrapper_inventory_2026-08-28.md
git commit -m "task doc — OC-wrapper step 0: read-only inventory of the prediction document"
```

*GATE:* exactly one file matches. If the copy in `~/Downloads` has a different stem (past sessions have had hyphens stripped), copy it under the `dev/tasks/` name above regardless, and report both names.

---

## 3. Locate the inputs by content — GATE

Nothing below names a path from memory. Every path is discovered and recorded.

### 3.1 The prediction document

```bash
grep -rl --include='*.qmd' -e 'fs_sym_root' . | sort
grep -rl --include='*.qmd' -e 'warmup-floors' . | sort
```

The prediction document is the **single** `.qmd` that contains both `forestsearch::fs_sym_root` and `warmup-floors`. *GATE:* exactly one file qualifies (the frozen parent documents predate the `fs_sym_root()` promotion and should not match the namespaced call; if more than one file matches, stop and list them). Record: path; `wc -l`; `git log -1 --format='%h %ci %s' -- <path>`; `git status --porcelain -- <path>` (must be clean — the inventory describes the committed file). Call it `DOC` below; its stem is `STEM`. A second check closes in A1: the index must list each of the chunk labels `anchor`, `warmup-floors`, `worked-predictions`, `worked-sensitivity`, `worked-null`, `appendix-tier2`, `worked-prediction` exactly once — *GATE*.

Also list, without inventorying them, every other `.qmd` and `.R` file in DOC's directory — path and YAML title or first comment line — so the report shows what DOC sits beside (frozen parents, drivers, scratch).

### 3.2 The payload the document reads

From the script's §7 (file I/O) after A1 — or `grep -n 'readRDS' DOC` — record every `readRDS()`/`source()`/`file.path()` the document executes, resolved to an actual path on disk: exists? tracked (`git ls-files --error-unmatch`)? `git log -1 --format='%h %ci %s' -- <path>`. Then, metadata only:

```r
p <- readRDS("<path>")
names(p); names(p$scale); str(p$scale, max.level = 1); str(p$scale$regions)
```

If the document reads more than one payload, or a payload conditionally (the null branch), record each and the condition.

### 3.3 The REPORT directory

```bash
git log --diff-filter=A --name-only --format= -- '*REPORT_*' | head -3
```

The REPORT directory is the directory of the **first** line printed (the most recently added `REPORT_*` file). *GATE:* the command prints at least one line. Record it as `RDIR`. The report is `RDIR/REPORT_oc_wrapper_inventory_2026-08-28.md`.

### 3.4 The drivers and the measured payloads (for Part B)

```bash
git show --stat --format='%h %s' dafff647 d884adbf 2b180813
```

`dafff647` guarded the null cell's `dgm_scale` path, `d884adbf` made the save-site name vectors tolerate an absent `scale` element, `2b180813` committed the null-cell payload. The `.qmd`/`.R` files those commits touched that call `forestsearch(` — DOC excluded — are the drivers; the `.rds` that `2b180813` added is the null-cell measurement payload. Record paths. If the alternative-cell measurement payloads are not identified this way, find them from the document's `m500` block (A3 stage S14) — its comment says where its literals were tabulated from. If still unidentified: CANNOT DETERMINE, and Part B runs on whatever was found.

---

## 4. Part A — what the prediction document's chunks compute

### A1. Mechanical index (completeness layer)

Save Appendix 1 verbatim as `RDIR/qmd_chunk_index.R` and run it:

```bash
SHA=$(git rev-parse --short HEAD)
Rscript RDIR/qmd_chunk_index.R DOC RDIR/chunk_index_STEM_${SHA}.md
grep -cE '^\s*`{3,}\{r' DOC          # fenced R chunk openers, for the cross-check
```

The script is a static parser: it evaluates nothing. For every fenced chunk it records label, options, line range, objects assigned or modified in place, free symbols (inputs), `pkg::` calls, RNG/draw calls, analytic calls (`pnorm`, `pmvnorm`, `integrate`, `uniroot`, …), file I/O, numeric literals; then the chunk dependency graph (each free symbol attributed to the last earlier chunk that assigned it) and the objects the prose reads through inline `` `r …` ``.

Two known limits of a static parse, to read the index correctly: symbols inside data-masked calls (`with()`, `within()`, `subset()`, `transform()`, formula terms) appear as free symbols — they are column names, and their appearance is itself evidence for the F classification; and rules built as strings and evaluated (`eval(parse(text = …))`) are invisible to the parser, so S1 must read them from the code directly.

Record in the report: the R-chunk count from the index and from the `grep -c`, and whether they agree (if not, say which chunks differ and why — a fence variant the regex missed is the likely cause; report it, do not patch the script). Any chunk the script reports as a parse error: the document renders, so a parse error in an evaluated chunk means the script misread it — quote the chunk header and first three lines and continue. The index file is committed with the report; do not paste it into the report, link it.

*GATE (from §3.1):* the index lists each of the seven named chunk labels exactly once.

The index is the completeness check for A2: **every R chunk in the index appears in A2 at least once**, if only as "display / assertion / no quantity".

### A2. Quantity table

One row per **quantity**, not per chunk. A chunk that prints seven figures yields at least seven rows; intermediate objects that feed later quantities (a covariance matrix, a draw matrix, a family table) get rows too, because they are the inputs the later rows cite.

| id | quantity (symbol and, where printed, its printed label) | chunk | definition as implemented | inputs | constants and draw counts | kind | computation | inputs fixture-bound? | scale | evidence |
|---|---|---|---|---|---|---|---|---|---|---|

- **definition as implemented** — one line; a formula where it is short, otherwise the operation. Written from the code, not from the surrounding prose.
- **inputs** — the objects, payload fields, and package calls it consumes (from the index's free-symbol list, resolved to what they are).
- **constants and draw counts** — every numeric literal and every draw-count variable that enters, with its value.
- **kind** — one of: `payload-read` (read from the payload, not derived) · `analytic` (closed form: `pnorm`, `pmvnorm`, `integrate`, `uniroot` on a closed form) · `MC` (a Monte-Carlo proportion or expectation, with its draw count) · `typed` (a literal in code) · `assertion` · `display` (formatting, `kable`, MC-SE strings).
- **computation** — `G` or `F` by this criterion, and no other:
  - **G (general):** the operation is expressed only in terms of (i) the DGM's regions/scale table as `fs_dgm_scale()` returns it, or quantities the package derives from it; (ii) `n`, `c1`, `c2`; (iii) an abstract candidate family — membership indicators on the DGM population, or prevalences and pairwise overlaps derived from them; (iv) draw counts and seeds. Test: handed those inputs for a *different* continuous-outcome DGM, would the code run unchanged?
  - **F (fixture-specific):** the operation refers to `df_super` columns by name, to the specific Q rule, to ACTG175 variables, to `seed_base`, to payload elements other than the scale/regions table, or to typed numbers that belong to this fixture (the super-population size, the sample sizes, the family counts, |Q|).
- **inputs fixture-bound?** — separate from the computation column: `no`, or `yes: <which inputs>`. A general operation fed by a fixture-built family is `G` / `yes: family`. This split is the point of the table; do not collapse it.
- **scale** — `MD-only` if the row relies on the continuous-scale pathway (`fs_dgm_scale()`'s MD branch, the `C_g[mu_0, tau]` negligibility argument, or the `β(Ĥ)` constancy check — v5 §6 records why none of these carries to log-OR); `any` if it is scale-free once the scale table exists; `unknown` if you cannot tell from source.
- **evidence** — the verbatim code lines (fenced, below the table, keyed by id, if too long for a cell).

### A3. The assembly, stage by stage

For each stage below, answer every question from source with the code quoted. Where the document computes a stage in more than one place (a base cell, a sensitivity row, a tier-2 cell, a null cell), record each and say what differs. If a stage on this list does not exist in the document, say so. If the document computes something not on this list, add a stage — the list is a coverage check, not a bound.

**S1 — the candidate family.** Which object holds it. How rules are enumerated: covariates used; for continuous covariates, how cut-points are chosen (grid, quantiles, observed values — and how many); factor levels and both directions or not; conjunction depth (one-factor, two-factor); the count before filtering. The prevalence filter: the exact expression; the population it is evaluated on (`df_super` at its full size, or a sample of `n`); a lower bound only or also an upper bound; the count after filtering; whether Q is a member by construction, by forcing, or by chance; any further filter (duplicates, near-duplicates, containment). The per-member quantities computed (prevalence, effect, SE, overlap with Q, …) with their formulas. Finally, **the family's interface**: exactly which fields of the family object the downstream stages read.

**S2 — the covariance across candidates.** How the covariance matrix of the candidate statistics is formed — the formula in terms of overlaps, prevalences, `v_cond0`, `sigma`, `n`, or the package call that produces it; its dimension; where `fs_sym_root()` enters and what the root feeds; how a rank-deficient matrix is handled.

**S3 — single-split declaration.** The per-member declaration event: the exact conditions (v5 §2 speaks of a member's "three conditions"), the constants they compare against (`c1`, `c2`, others), and what "single split" is — how the consistency screen's split is represented. The family-level event (any member declares). Whether it is computed by MC on the joint-normal draws (which draw count), analytically (`pmvnorm`), or both — and which one the document reports.

**S4 — detection.** How the selected subgroup Ĥ is chosen among declaring members (the selection rule as coded). The detection event as coded (Ĥ identical to Q? an overlap criterion? something else?). Draw count.

**S5 — `E|Ĥ|`.** The definition (size of the selected subgroup — in subjects at `n`, or as a population prevalence), whether conditional on declaration, and how non-declaration enters.

**S6 — the classification metrics.** In `worked-predictions`, the figures printed after `E|Ĥ|`: name each **from the code** (its symbol and printed label — do not take the names from the handoff), its formula, and what it conditions on.

**S7 — the inversions.** `c1` for target power: the objective function (which probability it inverts), the root finder, bracket and tolerance, and what happens when the target is unattainable. The `c2` ceiling: formula and inputs. The screening-only clearance: formula and inputs.

**S8 — `warmup-floors`.** The four figures: the definition of each and whether each is MC or analytic.

**S9 — `worked-sensitivity`.** Which inputs are varied, over which values, what is held fixed, and whether the draws are shared across rows or re-drawn.

**S10 — `worked-null`.** How the null prediction is formed: what the scale table is under the null, any guard on the `dgm_scale` path under the null (cf. `dafff647` — locate it, whether in DOC or in a driver) and what it guards, the definition of false declaration as coded, the `P1` range, `L_eff`.

**S11 — `appendix-tier2`.** What tier 2 is (a different family, a different filter, a different `n`?), and how the top-7 selection frequencies are computed.

**S12 — `anchor` and the assertions.** Every `stopifnot()`/`if (…) stop()` in the document: what it asserts and against what. Which `.anchor` elements are read from the payload and which are derived (v5 §2 says `m_ITT`, `v_cond0` are read). Where `.mc_se_prop()` is defined and each site where it is applied.

**S13 — the draw structure.** Every draw-count variable (`Rg`, `Rdraw`, `R7`, and any others): where defined, its value, every quantity it feeds. Every `set.seed()` and what it seeds. Whether one set of draws serves several quantities (paired) or each quantity re-draws.

**S14 — measured comparators embedded in the document.** The `m500` block: quote it with its comment. Which payload(s) and which columns its literals were tabulated from, as far as the comment and `git log -L` / `git blame` on the block establish. Which other typed comparisons to measurement exist in the document (prediction-vs-measured lines), and where their measured side comes from. This is a read, not a conversion — the m500 residue item stays where it is.

### A4. Constants and knobs

From the index's §5 (numeric literals by chunk), produce one table covering **every** numeric literal in the R chunks and every YAML `params:` entry:

| chunk | literal | role (what it is, from the code) | class |
|---|---|---|---|

with `class` one of: `fixture` (a value belonging to this fixture — sizes, counts, the Q cut-points, seeds) · `knob` (a design choice a different analysis might set differently — the prevalence floor, `n`, `c1`, `c2`, the power targets, draw counts, rule depth, cut grid, tolerances that are choices) · `math` (0, 1, 2, 0.5, `1e-12`-type tolerances that are not choices) · `display` (digits, widths). The 08-27 literal sweep covered literals duplicating payload quantities; this one is the design-constant sweep, and it is defined by the index, so nothing can be missed by looking only at prose.

### A5. The package boundary and the in-document helpers

For each `forestsearch::` function the document calls (index §6), from the **current `R/` source** (never the installed help page or the KB): the file; the signature verbatim; the structure of what it returns (names); and the arguments the document passes at each call site, each marked fixture-bound or general per the A2 criterion.

Then every function the document defines itself (`name <- function(...)`, from index §4): one line on what it does, and whether it duplicates logic that exists in `R/` — name the package function only if you have read its body and it is the same operation; otherwise "no counterpart found in R/" or CANNOT DETERMINE.

Part C reads `R/` further; A5 reads only the functions the document calls.

---

## 5. Part B — what the measurement side records about the realized family

*Runs — named in the disposition above.* The prediction runs a fixed family; the search builds its family from the data, and its size can vary across replicates. Whether it does is a measurement, not a conjecture — if the existing payloads recorded it. This part finds out, and checks that prediction and measurement score the same events.

**B1.** Drivers and payloads from §3.4. Per driver: path, `git log -1`, the `forestsearch()` call verbatim (all arguments), the `saveRDS()` site(s) verbatim.

**B2.** From the driver source, the per-replicate record: every field written (the recorder block — quote it). State whether any field records the realized candidate-family size — candidates enumerated, candidates passing the consistency screen, rows of the returned candidate table — and which object it is taken from. If none does: say so, and identify by content the line where such a field would be written (a proposal in the report; no edit).

**B3.** Per measured payload (null cell from `2b180813`; the alternative-cell payloads if §3.4 identified them): `names()`, `str(max.level = 2)`, the per-replicate table's column names and `nrow()`. Whether the per-replicate record carries the selected subgroup's definition (so that detection could be re-scored under the document's S4 definition): yes/no, with the column named.

**B4 — DIAGNOSTIC (the only computation in this task).** If B2/B3 found a per-replicate family-size field: its five-number summary across replicates and the count of replicates at zero, per payload. This is a description of a stored column for this report only; it is not a result and must not be presented as one. If no such field exists, B4 is "not available" and nothing is computed.

**B5 — definition parity.** For each of: declaration / false declaration; detection; `E|Ĥ|`; each classification metric of S6: the document's definition (from A3) beside the driver's scoring expression (verbatim). Verdict per row: `same event` / `different event: <how>` / `CANNOT DETERMINE`. A prediction-vs-measurement gap is interpretable only where the row says `same event`.

---

## 6. Part C — how the search builds its family, from `R/` source

*Runs — named in the disposition above.* **This part reads `R/` and changes nothing; none of the three `R/` categories (moves code / changes behaviour / changes the method) applies.** The KB is stale for R code; read the files at HEAD.

**C1.** `forestsearch()`'s signature verbatim from its source file (locate with `grep -rn '^forestsearch *<- *function' R/`), every argument with its default.

**C2.** The candidate-family pipeline, from data in to Ĥ out, as a table: step | file and function | governing arguments (with defaults) | data-dependent at fixed `n`? (yes/no — does its outcome vary across replicates drawn from the same DGM) | evidence. Steps to locate, by function not by name: rule enumeration (which covariates, how cut-points are formed, conjunction depth); any pre-screen before the consistency step; subgroup-size constraints (minimum and maximum, and on what they are evaluated); the consistency screen (`consistency_method = "resample"` / `method = "closed"` and the `"split"` fallback, per v5 §4) and its threshold; near-duplicate removal (`remove_near_duplicate_subgroups()`; `.maxeff_membership_dedup` on the `maxeff` path, per v5 §5); ordering and selection (`setorder(-HR, K)`, `sg_focus`); any time or count cap. Names in this paragraph come from the v5 handoff and are to be verified, not assumed.

**C3.** What the returned object records about the realized family: the candidate table (`found.hrs` per v5 §4), the consistency table, and any count — with the field names from source. State plainly whether a caller can read the realized family size and the realized family's rule definitions from a normal `forestsearch()` return.

**C4.** The search's size constraints beside the document's prevalence filter (S1): the same quantity or not (population prevalence on `df_super` vs sample proportion at `n`; lower bound only vs two-sided), from source.

**C5.** Whether candidate enumeration is callable without running the search — an exported or internal function that turns a data frame and the `forestsearch()` arguments into the enumerated rule set before screening. Name it from source or write "none found".

---

## 7. The report

`RDIR/REPORT_oc_wrapper_inventory_2026-08-28.md`, in this order:

0. Provenance (raw, §2).
1. Inputs located (§3): paths, hashes, counts, gates passed.
2. Part A: A1 (counts cross-check, parse errors, link to the index file) · A2 table · A3 stages S1–S14 · A4 table · A5.
3. Part B (B1–B5).
4. Part C (C1–C5).
5. Open items: every CANNOT DETERMINE, every defect noticed (described, not fixed), everything skipped and why.
6. Verdict, ten lines or fewer: the facts most relevant to the family question, stated as facts. No recommendation, no design.

Evidence is fenced and verbatim. Nothing in the report is a figure copied from a handoff. Line numbers, where given, sit next to the HEAD hash.

---

## 8. Close-out

```bash
git status --porcelain            # only: the report, the script, the index file — plus any pre-existing dirt from §2, unchanged
git add RDIR/REPORT_oc_wrapper_inventory_2026-08-28.md RDIR/qmd_chunk_index.R RDIR/chunk_index_STEM_${SHA}.md
git commit -m "report — OC-wrapper step 0: read-only inventory of the prediction document"
git log --oneline -3
git status --porcelain
```

No push. End the session with the verdict block from §7 item 6 pasted into chat.

---

## 9. Out of scope — do not touch

- The prediction document, the drivers, anything in `R/`, `tests/`, `NEWS.md`, `DESCRIPTION`: no edits of any kind.
- The residue items in v5 §6 (m500 conversion, secondary prose literals, display precision, the three batch documents, the CI-length median flag, the roxygen `%` sweep): read where A3 requires it, change nothing.
- No render of any document; no `devtools::install()`; no simulation of any size, including smoke tests.
- No effective-M back-calculation, no design proposal, no signature, no wrapper code — in the report or in chat.

---

## Appendix 1 — `qmd_chunk_index.R` (save verbatim into `RDIR/`)

```r
#!/usr/bin/env Rscript
# qmd_chunk_index.R — read-only static index of the R code chunks in a .qmd
#
# Usage:  Rscript qmd_chunk_index.R <file.qmd> <out.md>
#
# Base R only; no package is loaded and nothing in the document is evaluated.
# For every code chunk the script records: label, header options and `#|`
# options, line range, top-level objects assigned (or modified in place),
# free symbols (used but not assigned earlier in the same chunk), functions
# called, namespace-qualified calls (pkg::fun), numeric literals, RNG /
# draw / root-finding / integration calls, and file I/O calls. It then
# derives the chunk dependency graph (which earlier chunk last assigned each
# free symbol) and lists the objects that the prose reads through inline
# `r ...` code. Anything the parser cannot handle is reported, never guessed.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) stop("usage: Rscript qmd_chunk_index.R <file.qmd> <out.md>")
qmd_path <- args[[1L]]
out_path <- args[[2L]]
stopifnot(file.exists(qmd_path))
lines <- readLines(qmd_path, warn = FALSE)

## ---- 1. locate chunks ------------------------------------------------------
open_re  <- "^\\s*(`{3,})\\{([A-Za-z=][A-Za-z0-9_.]*)(.*)\\}\\s*$"
close_re <- "^\\s*(`{3,})\\s*$"
chunks <- list()
i <- 1L
n <- length(lines)
while (i <= n) {
  if (grepl(open_re, lines[i])) {
    fence <- sub(open_re, "\\1", lines[i])
    lang  <- sub(open_re, "\\2", lines[i])
    hdr   <- trimws(sub(open_re, "\\3", lines[i]))
    j <- i + 1L
    while (j <= n && !(grepl(close_re, lines[j]) &&
                       nchar(sub(close_re, "\\1", lines[j])) >= nchar(fence))) j <- j + 1L
    body <- if (j > i + 1L) lines[(i + 1L):(j - 1L)] else character(0)
    chunks[[length(chunks) + 1L]] <- list(lang = lang, header = hdr, start = i,
                                          end = min(j, n), body = body,
                                          unterminated = j > n)
    i <- j + 1L
  } else i <- i + 1L
}

## ---- 2. per-chunk options and label --------------------------------------
parse_header <- function(hdr) {
  hdr <- sub("^,\\s*", "", hdr)
  if (!nzchar(hdr)) return(list(label = NA_character_, opts = character(0)))
  # split on commas outside quotes/parens (good enough for chunk headers)
  parts <- character(0); depth <- 0L; cur <- ""; inq <- ""
  for (ch in strsplit(hdr, "")[[1L]]) {
    if (nzchar(inq)) { cur <- paste0(cur, ch); if (ch == inq) inq <- ""; next }
    if (ch %in% c("'", '"')) { inq <- ch; cur <- paste0(cur, ch); next }
    if (ch %in% c("(", "[", "{")) depth <- depth + 1L
    if (ch %in% c(")", "]", "}")) depth <- depth - 1L
    if (ch == "," && depth == 0L) { parts <- c(parts, trimws(cur)); cur <- ""; next }
    cur <- paste0(cur, ch)
  }
  parts <- c(parts, trimws(cur))
  parts <- parts[nzchar(parts)]
  label <- NA_character_
  if (length(parts) && !grepl("=", parts[1L], fixed = TRUE)) {
    label <- parts[1L]; parts <- parts[-1L]
  }
  list(label = label, opts = parts)
}

for (k in seq_along(chunks)) {
  ch <- chunks[[k]]
  ph <- parse_header(ch$header)
  hash_opts <- grep("^\\s*#\\|", ch$body, value = TRUE)
  hash_opts <- trimws(sub("^\\s*#\\|\\s?", "", hash_opts))
  lab <- ph$label
  hl <- grep("^label:\\s*", hash_opts, value = TRUE)
  if (is.na(lab) && length(hl)) lab <- trimws(sub("^label:\\s*", "", hl[1L]))
  if (is.na(lab)) lab <- sprintf("unlabelled-%02d", k)
  chunks[[k]]$label <- lab
  chunks[[k]]$opts  <- c(ph$opts, hash_opts)
  chunks[[k]]$code  <- ch$body[!grepl("^\\s*#\\|", ch$body)]
}

## ---- 3. AST walkers ----------------------------------------------------------
ASSIGN_OPS <- c("<-", "<<-", "=", "->", "->>")
MODIFIER_FUNS <- c("$", "[[", "[", "@", "attr", "names", "dimnames", "levels",
                   "class", "dim", "colnames", "rownames", "body", "formals",
                   "environment", "attributes", "diag", "substr", "is.na")
RNG_FUNS <- c("set.seed", "rnorm", "runif", "rbinom", "rpois", "rexp", "rgamma",
              "rbeta", "rt", "rchisq", "sample", "sample.int", "rmvnorm", "mvrnorm",
              "RNGkind", "nextRNGStream", "nextRNGSubStream")
ANALYTIC_FUNS <- c("pnorm", "qnorm", "dnorm", "pmvnorm", "qmvnorm", "integrate",
                   "uniroot", "optimize", "optimise", "optim", "nlm", "nlminb",
                   "pt", "qt", "pchisq", "qchisq", "pbinom", "qbinom")
IO_FUNS <- c("readRDS", "saveRDS", "source", "load", "save", "read.csv",
             "read.table", "readr::read_csv", "fread", "file.path", "sys.source",
             "write.csv", "writeLines", "readLines", "here", "file.exists")

lhs_root <- function(e) {
  # returns list(root = symbol name or NA, text = deparsed LHS, modify = TRUE/FALSE)
  if (is.symbol(e)) return(list(root = as.character(e), text = as.character(e), modify = FALSE))
  if (is.character(e) && length(e) == 1L) return(list(root = e, text = e, modify = FALSE))
  if (is.call(e)) {
    f <- e[[1L]]
    fn <- if (is.symbol(f)) as.character(f) else paste(deparse(f), collapse = "")
    if (fn %in% MODIFIER_FUNS || grepl("<-$", fn)) {
      inner <- lhs_root(e[[2L]])
      return(list(root = inner$root, text = paste(deparse(e), collapse = ""), modify = TRUE))
    }
    if (fn %in% c("(", "{")) return(lhs_root(e[[2L]]))
  }
  list(root = NA_character_, text = paste(deparse(e), collapse = ""), modify = TRUE)
}

new_acc <- function() {
  a <- new.env()
  for (nm in c("assigned", "assigned_text", "uses", "funs", "ns", "nums", "strs",
               "rng", "analytic", "io", "io_text", "pkgs")) assign(nm, character(0), envir = a)
  a
}

is_empty <- function(x) identical(x, quote(expr = ))
walk_elems <- function(lst, acc, locals, from = 1L) {
  # iterate by index so the empty symbol (e.g. x[1, ]) is never bound to a variable
  if (length(lst) < from) return(invisible())
  for (i in from:length(lst)) if (!is_empty(lst[[i]]) && !is.null(lst[[i]])) walk(lst[[i]], acc, locals)
  invisible()
}

walk <- function(e, acc, locals = character(0)) {
  # acc: environment with character vectors: assigned, assigned_text, uses,
  #      funs, ns, nums, strs, rng, analytic, io, io_text, local_defs
  if (is.symbol(e)) {
    nm <- as.character(e)
    if (nzchar(nm) && !(nm %in% locals)) acc$uses <- c(acc$uses, nm)
    return(invisible())
  }
  if (is.numeric(e) && length(e) == 1L) { acc$nums <- c(acc$nums, format(e, digits = 15)); return(invisible()) }
  if (is.character(e) && length(e) == 1L) { acc$strs <- c(acc$strs, e); return(invisible()) }
  if (!is.call(e)) return(invisible())
  f <- e[[1L]]
  if (is.symbol(f)) {
    fn <- as.character(f)
    if (fn %in% ASSIGN_OPS && length(e) == 3L) {
      if (fn %in% c("->", "->>")) { lhs <- e[[3L]]; rhs <- e[[2L]] } else { lhs <- e[[2L]]; rhs <- e[[3L]] }
      # RHS first (uses), then LHS (assignment)
      walk(rhs, acc, locals)
      lr <- lhs_root(lhs)
      if (!is.na(lr$root)) {
        if (lr$modify) {
          # modifying an existing object reads it first
          if (!(lr$root %in% locals)) acc$uses <- c(acc$uses, lr$root)
          # index/attr arguments may use symbols too
          if (is.call(lhs) && !(as.character(lhs[[1L]])[1L] %in% c("$", "@")))
            walk_elems(as.list(lhs), acc, locals, from = 3L)
        }
        acc$assigned <- c(acc$assigned, lr$root)
        acc$assigned_text <- c(acc$assigned_text, lr$text)
      } else {
        walk(lhs, acc, locals)
      }
      return(invisible())
    }
    if (fn == "function") {
      fml <- e[[2L]]
      fnames <- if (is.null(fml)) character(0) else names(fml)
      # defaults may reference globals
      if (!is.null(fml)) walk_elems(as.list(fml), acc, c(locals, fnames))
      body_acc <- new_acc()
      walk(e[[3L]], body_acc, c(locals, fnames))
      inner_locals <- unique(c(fnames, body_acc$assigned))
      acc$uses <- c(acc$uses, setdiff(body_acc$uses, inner_locals))
      acc$funs <- c(acc$funs, setdiff(body_acc$funs, inner_locals))
      acc$ns <- c(acc$ns, body_acc$ns); acc$nums <- c(acc$nums, body_acc$nums)
      acc$strs <- c(acc$strs, body_acc$strs); acc$rng <- c(acc$rng, body_acc$rng)
      acc$analytic <- c(acc$analytic, body_acc$analytic); acc$io <- c(acc$io, body_acc$io)
      acc$io_text <- c(acc$io_text, body_acc$io_text); acc$pkgs <- c(acc$pkgs, body_acc$pkgs)
      return(invisible())
    }
    if (fn == "local" && length(e) >= 2L) {
      # local({ ... }) scopes its assignments; only its free symbols reach the chunk
      body_acc <- new_acc()
      walk(e[[2L]], body_acc, locals)
      acc$uses <- c(acc$uses, setdiff(body_acc$uses, body_acc$assigned))
      acc$funs <- c(acc$funs, "local", setdiff(body_acc$funs, body_acc$assigned))
      acc$ns <- c(acc$ns, body_acc$ns); acc$nums <- c(acc$nums, body_acc$nums)
      acc$strs <- c(acc$strs, body_acc$strs); acc$rng <- c(acc$rng, body_acc$rng)
      acc$analytic <- c(acc$analytic, body_acc$analytic); acc$io <- c(acc$io, body_acc$io)
      acc$io_text <- c(acc$io_text, body_acc$io_text); acc$pkgs <- c(acc$pkgs, body_acc$pkgs)
      return(invisible())
    }
    if (fn == "for" && length(e) == 4L) {
      acc$assigned <- c(acc$assigned, as.character(e[[2L]]))
      acc$assigned_text <- c(acc$assigned_text, paste0("for-var ", as.character(e[[2L]])))
      walk(e[[3L]], acc, locals); walk(e[[4L]], acc, c(locals, as.character(e[[2L]])))
      return(invisible())
    }
    if (fn == "assign" && length(e) >= 3L && is.character(e[[2L]])) {
      acc$assigned <- c(acc$assigned, e[[2L]])
      acc$assigned_text <- c(acc$assigned_text, paste0("assign(\"", e[[2L]], "\")"))
      walk_elems(as.list(e), acc, locals, from = 3L)
      return(invisible())
    }
    if (fn %in% c("library", "require", "requireNamespace", "loadNamespace") && length(e) >= 2L) {
      acc$funs <- c(acc$funs, fn)
      acc$pkgs <- c(acc$pkgs, paste(deparse(e[[2L]]), collapse = ""))
      return(invisible())
    }
    if (fn %in% c("::", ":::") && length(e) == 3L) {
      acc$ns <- c(acc$ns, paste0(deparse(e[[2L]]), fn, deparse(e[[3L]])))
      return(invisible())
    }
    if (fn %in% c("$", "@") && length(e) == 3L) {
      walk(e[[2L]], acc, locals)   # the name after $ is not a free symbol
      return(invisible())
    }
    if (!(fn %in% locals)) acc$funs <- c(acc$funs, fn)
    if (fn %in% RNG_FUNS) acc$rng <- c(acc$rng, fn)
    if (fn %in% ANALYTIC_FUNS) acc$analytic <- c(acc$analytic, fn)
    if (fn %in% IO_FUNS) {
      acc$io <- c(acc$io, fn)
      acc$io_text <- c(acc$io_text, paste(deparse(e, width.cutoff = 500L), collapse = " "))
    }
    walk_elems(as.list(e), acc, locals, from = 2L)
    return(invisible())
  }
  # call whose function is itself an expression, e.g. pkg::f(x) or obj$f(x)
  if (is.call(f) && is.symbol(f[[1L]]) && as.character(f[[1L]]) %in% c("::", ":::")) {
    q <- paste0(deparse(f[[2L]]), as.character(f[[1L]]), deparse(f[[3L]]))
    acc$ns <- c(acc$ns, q); acc$funs <- c(acc$funs, q)
    bare <- deparse(f[[3L]])
    if (bare %in% RNG_FUNS) acc$rng <- c(acc$rng, q)
    if (bare %in% ANALYTIC_FUNS) acc$analytic <- c(acc$analytic, q)
    if (bare %in% IO_FUNS) { acc$io <- c(acc$io, q); acc$io_text <- c(acc$io_text, paste(deparse(e, width.cutoff = 500L), collapse = " ")) }
    walk_elems(as.list(e), acc, locals, from = 2L)
    return(invisible())
  }
  walk_elems(as.list(e), acc, locals, from = 1L)
  invisible()
}


analyse_chunk <- function(code) {
  res <- list(parse_error = NA_character_, free = character(0), assigned = character(0),
              assigned_text = character(0), funs = character(0), ns = character(0),
              nums = character(0), rng = character(0), analytic = character(0),
              io_text = character(0), pkgs = character(0), n_expr = 0L)
  if (!length(code) || !any(nzchar(trimws(code)))) return(res)
  exprs <- tryCatch(parse(text = code, keep.source = FALSE),
                    error = function(e) { res$parse_error <<- conditionMessage(e); NULL })
  if (is.null(exprs)) return(res)
  res$n_expr <- length(exprs)
  defined <- character(0)
  acc_all <- new_acc()
  free <- character(0)
  for (ex in as.list(exprs)) {
    a <- new_acc()
    walk(ex, a)
    # free = used in this expression and not yet defined by an earlier expression
    # (a symbol both used and assigned inside ONE expression, e.g. x <- x + 1, counts as free)
    free <- c(free, setdiff(a$uses, defined))
    defined <- unique(c(defined, a$assigned))
    for (nm in ls(acc_all)) assign(nm, c(get(nm, envir = acc_all), get(nm, envir = a)), envir = acc_all)
  }
  res$free <- unique(free)
  res$assigned <- unique(acc_all$assigned)
  res$assigned_text <- unique(acc_all$assigned_text)
  fl <- unique(acc_all$funs)
  res$funs <- fl[grepl("^[A-Za-z.][A-Za-z0-9._]*$", fl) | grepl(":::?", fl)]
  res$funs <- setdiff(res$funs, c("if", "for", "while", "repeat", "break", "next", "function"))
  res$pkgs <- unique(acc_all$pkgs)
  res$ns <- unique(acc_all$ns)
  res$nums <- acc_all$nums
  res$rng <- unique(acc_all$rng)
  res$analytic <- unique(acc_all$analytic)
  res$io_text <- unique(acc_all$io_text)
  res
}

for (k in seq_along(chunks)) {
  chunks[[k]]$an <- if (chunks[[k]]$lang == "r") analyse_chunk(chunks[[k]]$code) else NULL
}

## ---- 4. dependency graph ---------------------------------------------------
# a symbol used freely in chunk k is attributed to the LAST earlier R chunk that assigned it
edges <- data.frame(chunk = character(0), depends_on = character(0), symbols = character(0),
                    stringsAsFactors = FALSE)
external <- list()
for (k in seq_along(chunks)) {
  an <- chunks[[k]]$an
  if (is.null(an)) next
  # function names that are not defined anywhere in the document are base/package calls;
  # only symbols used as *values* or as functions defined by an earlier chunk make edges
  cand <- unique(c(an$free, an$funs))
  ext <- character(0)
  by_src <- list()
  for (s in cand) {
    src <- NA_integer_
    if (k > 1L) for (j in (k - 1L):1L) {
      if (!is.null(chunks[[j]]$an) && s %in% chunks[[j]]$an$assigned) { src <- j; break }
    }
    if (is.na(src)) { if (s %in% an$free) ext <- c(ext, s) } else {
      key <- chunks[[src]]$label
      by_src[[key]] <- c(by_src[[key]], s)
    }
  }
  for (key in names(by_src))
    edges <- rbind(edges, data.frame(chunk = chunks[[k]]$label, depends_on = key,
                                     symbols = paste(sort(unique(by_src[[key]])), collapse = ", "),
                                     stringsAsFactors = FALSE))
  external[[chunks[[k]]$label]] <- sort(unique(ext))
}

## ---- 5. inline code in prose -------------------------------------------------
in_chunk <- rep(FALSE, n)
for (ch in chunks) in_chunk[ch$start:ch$end] <- TRUE
prose <- lines[!in_chunk]
inline <- unlist(regmatches(prose, gregexpr("`r [^`]+`", prose)))
inline_code <- sub("^`r ", "", sub("`$", "", inline))
inline_syms <- character(0)
inline_bad <- character(0)
for (ic in inline_code) {
  ex <- tryCatch(parse(text = ic, keep.source = FALSE), error = function(e) NULL)
  if (is.null(ex)) { inline_bad <- c(inline_bad, ic); next }
  a <- new_acc(); for (e1 in as.list(ex)) walk(e1, a)
  inline_syms <- c(inline_syms, a$uses)
}
inline_tab <- if (length(inline_syms)) sort(table(inline_syms), decreasing = TRUE) else NULL

## ---- 6. write the report -----------------------------------------------------
md <- character(0)
p <- function(...) md <<- c(md, paste0(...))
esc <- function(x) gsub("|", "\\|", x, fixed = TRUE)
p("# Static chunk index — `", basename(qmd_path), "`")
p("")
p("Generated by `qmd_chunk_index.R` (base R, static parse only; nothing evaluated).  ")
p("File: `", qmd_path, "` — ", n, " lines, ", length(chunks), " fenced code chunks (",
  sum(vapply(chunks, function(c) c$lang == "r", logical(1))), " R).  ")
p("Line numbers are valid only at the commit recorded in the provenance header of the report that embeds this index.")
p("")
p("## 1. Chunk index (document order)")
p("")
p("| # | label | lang | lines | n expr | options | assigns / modifies | free symbols (inputs) | pkg::calls | RNG / draw | analytic | parse |")
p("|---|---|---|---|---|---|---|---|---|---|---|---|")
for (k in seq_along(chunks)) {
  ch <- chunks[[k]]; an <- ch$an
  p("| ", k, " | `", esc(ch$label), "` | ", ch$lang, " | ", ch$start, "–", ch$end, " | ",
    if (is.null(an)) "" else an$n_expr, " | ", esc(paste(ch$opts, collapse = "; ")), " | ",
    if (is.null(an)) "" else esc(paste(sort(an$assigned), collapse = ", ")), " | ",
    if (is.null(an)) "" else esc(paste(sort(an$free), collapse = ", ")), " | ",
    if (is.null(an)) "" else esc(paste(sort(an$ns), collapse = ", ")), " | ",
    if (is.null(an)) "" else esc(paste(sort(an$rng), collapse = ", ")), " | ",
    if (is.null(an)) "" else esc(paste(sort(an$analytic), collapse = ", ")), " | ",
    if (is.null(an)) "n/a" else if (is.na(an$parse_error)) "ok" else paste0("ERROR: ", esc(gsub("\\s+", " ", an$parse_error))),
    if (isTRUE(ch$unterminated)) " (UNTERMINATED FENCE)" else "", " |")
}
p("")
p("## 2. Chunk dependency edges (free symbol → last earlier chunk that assigned it)")
p("")
if (nrow(edges)) {
  p("| chunk | depends on | via symbols |"); p("|---|---|---|")
  for (r in seq_len(nrow(edges))) p("| `", esc(edges$chunk[r]), "` | `", esc(edges$depends_on[r]), "` | ", esc(edges$symbols[r]), " |")
} else p("(no edges)")
p("")
p("```mermaid"); p("flowchart TD")
for (k in seq_along(chunks)) if (chunks[[k]]$lang == "r") p("  c", k, "[\"", gsub("\"", "'", chunks[[k]]$label), "\"]")
if (nrow(edges)) for (r in seq_len(nrow(edges))) {
  from <- which(vapply(chunks, function(c) c$label == edges$depends_on[r], logical(1)))[1L]
  to   <- which(vapply(chunks, function(c) c$label == edges$chunk[r], logical(1)))[1L]
  p("  c", from, " --> c", to)
}
p("```")
p("")
p("## 3. Symbols used but never assigned by any earlier chunk (external: package, base, or defined in prose/YAML)")
p("")
p("| chunk | external symbols |"); p("|---|---|")
for (key in names(external)) if (length(external[[key]])) p("| `", esc(key), "` | ", esc(paste(external[[key]], collapse = ", ")), " |")
p("")
p("## 4. Assignment targets in full (including in-place modifications such as `x$y <- …`)")
p("")
for (k in seq_along(chunks)) {
  an <- chunks[[k]]$an
  if (is.null(an) || !length(an$assigned_text)) next
  p("- `", esc(chunks[[k]]$label), "`: ", esc(paste(an$assigned_text, collapse = " ; ")))
}
p("")
p("## 5. Numeric literals by chunk (raw; classify their roles by hand)")
p("")
p("| chunk | literal | count |"); p("|---|---|---|")
for (k in seq_along(chunks)) {
  an <- chunks[[k]]$an
  if (is.null(an) || !length(an$nums)) next
  tb <- table(an$nums)
  for (v in names(tb)) p("| `", esc(chunks[[k]]$label), "` | ", v, " | ", tb[[v]], " |")
}
p("")
p("## 6. Namespace-qualified calls (all chunks)")
p("")
allns <- sort(unique(unlist(lapply(chunks, function(c) if (is.null(c$an)) NULL else c$an$ns))))
if (length(allns)) for (q in allns) {
  where <- vapply(chunks, function(c) !is.null(c$an) && q %in% c$an$ns, logical(1))
  p("- `", q, "` — chunks: ", paste(vapply(chunks[where], function(c) c$label, character(1)), collapse = ", "))
} else p("(none)")
p("")
p("## 6b. Packages attached / loaded (library, require, requireNamespace, loadNamespace)")
p("")
for (k in seq_along(chunks)) {
  an <- chunks[[k]]$an
  if (is.null(an) || !length(an$pkgs)) next
  p("- `", esc(chunks[[k]]$label), "`: ", esc(paste(an$pkgs, collapse = ", ")))
}
p("")
p("## 7. File I/O calls (verbatim)")
p("")
for (k in seq_along(chunks)) {
  an <- chunks[[k]]$an
  if (is.null(an) || !length(an$io_text)) next
  for (t in an$io_text) p("- `", esc(chunks[[k]]$label), "`: `", esc(t), "`")
}
p("")
p("## 8. Functions called, by chunk (bare names; base/package/document-defined not distinguished here)")
p("")
for (k in seq_along(chunks)) {
  an <- chunks[[k]]$an
  if (is.null(an) || !length(an$funs)) next
  p("- `", esc(chunks[[k]]$label), "`: ", esc(paste(sort(an$funs), collapse = ", ")))
}
p("")
p("## 9. Objects read by inline `r …` code in prose")
p("")
p(length(inline_code), " inline code spans; ", length(inline_bad), " failed to parse.")
p("")
if (!is.null(inline_tab)) {
  p("| symbol | inline uses |"); p("|---|---|")
  for (s in names(inline_tab)) p("| `", esc(s), "` | ", inline_tab[[s]], " |")
}
if (length(inline_bad)) { p(""); p("Unparseable inline spans:"); for (b in inline_bad) p("- `", esc(b), "`") }
p("")
writeLines(md, out_path)
cat("wrote", out_path, "-", length(chunks), "chunks;",
    sum(vapply(chunks, function(c) !is.null(c$an) && !is.na(c$an$parse_error), logical(1))),
    "parse errors\n")
```
