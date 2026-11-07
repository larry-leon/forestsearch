# TASK_fs-article-ms-alignment_2026-11-07

**Repo:** forestsearch. **Branch:** `fs-article-ms-alignment` from `master`. **Touches:**
`vignettes/articles/mrct-consistency.qmd` and its tracked render
`vignettes/articles/mrct-consistency.html` (re-render) only. **No `R/`.** **No simulation
compute:** render-only; expect roughly 10 minutes.

## Context (orientation, not evidence)

Five approved revisions (rendered preview, approved by Larry 2026-11-07), aligning the article with
the accompanying manuscript: (A) the trial-wide harmed subgroup leads §1.1's two explanations;
(B) S2's first occurrence gains a cross-reference; (C) a §1.4 callout points at the one-call
wrapper; (D) §2.6 gains a primary campaign-observed table in §1.3's unconditional row format, with
the existing conditional pipeline table kept beneath; (E) the manuscript's Table 1 (event notation)
is reproduced as a closing unnumbered appendix and referenced from §1.3. The appendix content
travels in the companion because the manuscript lives in the other repository.

**Companion file:** `EDITS_fs-article-ms-alignment_2026-11-07.md`.

## First action

Copy this task doc and the companion from `~/Downloads` into `dev/tasks/` and commit them together.
Run `export GIT_TERMINAL_PROMPT=0` at session start. `master` may be ahead of `origin/master`
(Larry pushes); not a stop.

## Legacy copies (standing instruction, before any edit)

Next free dated directory under `legacy_vignette/` (`06Nov2026/` through `06Nov2026_4/` and
`07Nov2026/` exist — never overwrite; next free suffix or the new date). Copy both files as they
stand at task start; `cmp`-verify. Untracked and local via `.git/info/exclude`; never touch
`.gitignore`. Record the copy paths in NOTES.

## Edit-application convention

Text blocks (A), (B), (E1): whitespace-normalized exact matching, one occurrence each. (C): insert
the callout block immediately before the paragraph beginning "$\hat H$ is learned on the non-region
patients only" (one occurrence). (D): a delimited replacement — keep the heading line matching
`^## The S2 campaign in brief`, replace everything after it up to (not including) the line matching
`^## The one-call form` with the companion's section content; both heading lines must occur exactly
once. (E2): append the companion's appendix section at the end of the file. New prose re-wrapped at
~100 characters; chunk code exactly as given.

## Row (2) of the new §2.6 table — specified either-way action

Inspect the column names of `core_full_numeric.csv` (the S2 row) and the other shipped core CSVs
for a count of region ITT > 1.0 **over all trials** (none is expected). If one exists, fill row (2)
from it and drop the caption's analytic-only clause; otherwise keep the dash and the clause as
given. Record the finding in NOTES.

## Gates (stop-on-failure)

- **G0 — state.** Working tree clean on `master`; record HEAD in NOTES.
- **G1 — anchors.** (A), (B), (E1) blocks and (C)'s insertion anchor each match exactly once; (D)'s
  two heading lines each occur exactly once. On any miss: STOP; NOTES pastes the passage.
- **G2 — apply.** `git diff` confined to the six regions ((A), (B), (C) insertion, (D) section,
  (E1), (E2) append). Byte-identical in particular: the (M1)–(M4) bullets, §1.3's named-rows
  bullets, all of §2.7, and the question/signal-generation sentences after (A).
- **G3 — chunks.** Every code chunk byte-identical except `case-campaign` (replaced) and
  `case-campaign-pipeline` (new); fenced-block count rises by exactly one.
- **G4 — render.** `quarto render` green; every `stopifnot` passes, the new identity assertions
  included.
- **G5 — render delta.** Against the tracked baseline HTML: outside the edited regions, timing
  positions and the render date, numeric tokens identical. Within the new §2.6 table: row (4) must
  render 0.240 (2,398) and row (8) 0.123 (1,228) — the baseline conditional table's own counts —
  and the conditional table beneath must be value-identical to the baseline's (caption aside). The
  appendix renders as a table; §1.3's appendix link resolves (anchor `sec-event-notation` present).
- **G6 — discipline.** No manuscript **values** enter (the appendix is definitions only); no
  provenance note beyond the approved row-(2) caption clause; inline-R additions are only the new
  §2.6 reading's reads of `obs`/`camp`.

## Close-out

- Commit the edits and the re-rendered tracked HTML; fast-forward merge to `master`; delete the
  branch (`git branch -r` checked first). Larry pushes. `git status --porcelain` empty at close.
- `NOTES_fs-article-ms-alignment_2026-11-07.md` to `~/Downloads`: HEADs, gate outcomes, the row-(2)
  finding, the G5 delta listing, the legacy-copy paths.
