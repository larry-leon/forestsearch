# TASK_fs-article-rationale-two-routes_2026-11-06

**Repo:** forestsearch. **Branch:** `fs-article-rationale-two-routes` from `master`. **Touches:**
`vignettes/articles/mrct-consistency.qmd` and its tracked render
`vignettes/articles/mrct-consistency.html` (re-render) only. **No `R/`.** **No simulation
compute:** render-only; expect roughly 10 minutes.

## Context (orientation, not evidence)

Two approved edits (rendered preview, approved by Larry 2026-11-06): the §1.2 hazard-scale bullet
gains its rationale and invariance details, and the §1.3 named-rows bullets spell out the two
routes to consistency (standard pathway vs the subgroup pathway in oracle form). Both anchors are
in text applied by TASK_fs-article-readability-bullets earlier today.

**Companion file:** `EDITS_fs-article-rationale-two-routes_2026-11-06.md` — the two old → new
blocks.

## First action

Copy this task doc and the companion from `~/Downloads` into `dev/tasks/` and commit them together.
Run `export GIT_TERMINAL_PROMPT=0` at session start. `master` may be ahead of `origin/master`
(Larry pushes); not a stop.

## Legacy copies (standing instruction, before any edit)

- Next free dated directory under `legacy_vignette/` (`06Nov2026/`, `06Nov2026_2/`, `06Nov2026_3/`
  exist — never overwrite; use the next free suffix or the new date).
- Copy both files as they stand at task start; verify each with `cmp`. Untracked and local via
  `.git/info/exclude`; never touch `.gitignore`. Record the copy paths in NOTES.

## Edit-application convention

Whitespace-normalized matching (any run of spaces/newlines ≡ one space); each `old` block occurs
exactly once; `new` written with its list structure, prose re-wrapped at ~100 characters.

## Gates (stop-on-failure)

- **G0 — state.** Working tree clean on `master`; record HEAD in NOTES.
- **G1 — anchors.** Both `old` blocks match exactly once. On a miss: STOP; NOTES pastes the passage
  as it now stands.
- **G2 — apply.** `git diff` confined to the two passages; the event bullets above the named rows
  ($A$, $B$, $C^*$) byte-identical; the other §1.2 bullets byte-identical.
- **G3 — code integrity.** Every code chunk byte-identical (prose-only task).
- **G4 — stated ranges vs the records (STOP for rewording if out).** One R session:
  (a) the seed's event rate — `a <- subset(speff2trial::ACTG175, arms %in% c(1L, 3L));
  mean(a$cens)` must lie in [0.15, 0.35] (the bullet's "roughly 20–30%");
  (b) the scaled design — from the committed build record
  (`inst/extdata/mrct-consistency/build_full_summary.csv`, row `S2_rmlog4`):
  `event_prob_40_nonregion` in [0.50, 0.70], and the prevalence-weighted overall
  `prevalence*event_prob_40_region + (1-prevalence)*event_prob_40_nonregion` in [0.45, 0.70]
  (the bullet's "toward a modern trial's 50–70%"). Record all values in NOTES.
- **G5 — render.** `quarto render` green; every `stopifnot` passes.
- **G6 — render delta.** Against the tracked baseline HTML: outside the two passages and timing
  positions, numeric tokens identical; inside them, additions are the companion's verbatim text
  only (expected constants: 3, 20, 30, 40, 50, 70, 1.0 — none computed). Both passages still
  render as lists.
- **G7 — discipline.** No inline R added or removed; no manuscript-value comparison; no
  provenance note.

## Close-out

- Commit the edits and the re-rendered tracked HTML; fast-forward merge to `master`; delete the
  branch (`git branch -r` checked first). Larry pushes. `git status --porcelain` empty at close.
- `NOTES_fs-article-rationale-two-routes_2026-11-06.md` to `~/Downloads`: HEADs, gate outcomes,
  the G4 values, the G6 delta listing, the legacy-copy paths.
