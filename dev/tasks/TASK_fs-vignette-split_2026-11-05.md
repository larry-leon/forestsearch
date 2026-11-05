# TASK_fs-vignette-split_2026-11-05

**Dispositions on record (Larry, 2026-11-05 13:50, clarified 13:57)**: split the mrct-consistency article into two independently executable pkgdown articles — `mrct-consistency` keeps **Sections 1–2** (Overview; The case study: the application example), and a new `mrct-consistency-simulations` carries **Sections 3–6** (Setting up simulations from your own trial; How replicate 5 was generated; Simulation results; Appendix). Because §§4–5 check themselves against the case-study objects, the simulations article opens with a condensed recomputation of the Section-2 analysis (`#sec-case-objects`).

**Repo**: `forestsearch` only. **No `R/` changes of any kind.** Files touched: `vignettes/articles/mrct-consistency.qmd`, new `vignettes/articles/mrct-consistency-simulations.qmd`, and `_pkgdown.yml` only if it lists articles explicitly. **STOP** if anything the spec requires is missing from the installed dev package.
**Branch**: short-lived from `master` (record HEAD; verify from source), fast-forward, delete branch (`git branch -r` first).
**Compute**: none beyond the articles' own renders. The split adds one recomputed `forestsearch()` + MR call in the simulations article (the bridge; seconds-to-a-minute) and removes §§3–6 execution from the application article, which becomes much faster.
**Session start**: `export GIT_TERMINAL_PROMPT=0`. Credentials-failed fetch recorded, not a stop.
**Files from ~/Downloads**: this task doc and `vignette-split-spec_2026-11-05.md`. **First action**: copy both into `dev/tasks/` and commit them together as the task commit.
**Gates** stop-on-failure. `NOTES_fs-vignette-split_2026-11-05.md` to `~/Downloads` at close.

## 1. Execute the split per the spec

The spec file carries the new verbatim content (sims front matter, opening paragraph, bridge section, app closing line and intro addition), the per-article `needed` lists, and the cross-reference rewrite tables with a bounded rule for any occurrence the tables miss. Mechanics:

1. Build `mrct-consistency-simulations.qmd` per spec §C (assembly order; moved sections otherwise byte-unchanged except the listed rewrites).
2. Reduce `mrct-consistency.qmd` per spec §B (remove from `#sec-own-sims` to end; drop `nsims_demo`; apply the rewrite table; App `needed` list; append the closing timing line; insert the intro addition).
3. The bridge chunk (spec §A.3) is **assembled from the app article's committed chunks** as itemized there — same calls, display lines dropped, every `stopifnot` kept. Introduce no new computation.

## 2. Gates (stop-on-failure)

1. **Independence**: each article renders green in a **fresh session** (`quarto render` one, then the other, separate R processes); every `stopifnot` in both passes — including the bridge's checks and `#sec-find`'s redraw/equality checks, now against the bridge objects.
2. **Cross-references**: in each rendered HTML, zero unresolved references (no "?@sec-", no "**?**"); grep the two `.qmd` sources — every `@sec-` id referenced is defined in the same file.
3. **Object provenance**: in the simulations article, the objects `trial`, `covs`, `seed_rep5`, `fs`, `mr`, `mv`, `fl`, `Hhat_def`, `in_Hhat`, `Hstar`, `reg`, `four`, `set_label` are defined in the bridge chunk and nowhere else before use (grep).
4. **No drift**: the bridge's `forestsearch()` call is byte-identical to the app article's `case-fs` call (diff the call text).
5. **Content accounting**: every line of the pre-split article is in exactly one of the two articles, modulo the listed rewrites, the removed `nsims_demo` lines, the per-article front matter/scaffolding, and the new content of spec §A — produce the accounting (e.g., a filtered diff) in NOTES.
6. `git diff` confined to the two article files (+ `_pkgdown.yml` if applicable) and the `dev/tasks/` copies; `R/`, `DESCRIPTION`, `NAMESPACE`, `inst/extdata/` untouched.
7. If `_pkgdown.yml` lists articles: new entry added; `pkgdown::check_pkgdown()` (or the lightest equivalent check) passes. Record which case held.
8. Wall-clock of each article's render recorded.

## 3. Close-out

Commit sequence: task commit (both files in `dev/tasks/`); the split commit(s); fast-forward `master`; delete branch. NOTES with: HEADs; per-article render minutes; the content accounting of gate 5; any rule-based rewrites beyond the spec's tables; the bridge chunk's printed line (patients, MR seconds, Ĥ); pkgdown outcome; flags.

## Out of scope

`R/`; records; `DESCRIPTION`; any rewording beyond the spec's rewrites and new content; vocabulary harmonization; the mrct-consistency repo.
