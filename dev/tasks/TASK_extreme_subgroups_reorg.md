# TASK — Reorganize `quarto/extreme_subgroups/`; move slides and v1 feasibility doc to `hta-subpopulations`

**Date:** 2026-10-19
**Basis:** `dev/reports/REPORT_extreme_subgroups_inventory.md` (commit `7aaee08f`). Larry approved proposals 1–8 on 2026-10-19.
**Repos:** `~/Documents/GitHub/forestsearch` (run CC here), `~/Documents/GitHub/hta-subpopulations`
**R/ impact:** none. The two comments in `R/` that mention the old path stay as they are, because they cite commit `002f4f37`.
**Compute:** none. Nothing is rendered. Validation is parse-only (`quarto inspect`).
**Push:** never.
**Gates:** stop-on-failure.

---

## Part 0 — Record this task

1. `cp ~/Downloads/TASK_extreme_subgroups_reorg.md dev/tasks/`
2. `git add dev/tasks/TASK_extreme_subgroups_reorg.md`
3. `git commit -m "tasks: extreme_subgroups reorganization"`

## Part A — Preflight and tag (gate)

- A1. In both repos, `git status --porcelain` is empty apart from `.claude/settings.local.json`. Leave that file untouched and out of every commit.
- A2. All 40 files from the inventory D1 table are present under `quarto/extreme_subgroups/`, and the folder contains no other files.
- A3. `which quarto` resolves.
- A4. Tag the current HEAD before making any change. Record the tag's commit hash.
  `git tag -a extreme-subgroups-pre-reorg -m "State of quarto/extreme_subgroups before the 2026-10-19 reorganization (see REPORT_extreme_subgroups_inventory.md)"`

## Part B — Copy material to `hta-subpopulations` (gate)

Let `H=~/Documents/GitHub/hta-subpopulations` and `ES=quarto/extreme_subgroups`.

B1. Copy these files:

```
mkdir -p $H/presentations/2026-extreme-subgroups $H/docs/feasibility/v1
cp $ES/fixed_random/extreme_subgroups_slides.qmd   $H/presentations/2026-extreme-subgroups/
cp $ES/fixed_random/extreme_subgroups_slides.pdf   $H/presentations/2026-extreme-subgroups/
cp $ES/fixed_random/sim_flow_part1.pdf             $H/presentations/2026-extreme-subgroups/
cp $ES/fixed_random/sim_flow_part2.pdf             $H/presentations/2026-extreme-subgroups/
cp $ES/forestsearch_fixed_baseline_slides.pptx     $H/presentations/2026-extreme-subgroups/
cp $ES/build_deck.js                               $H/presentations/2026-extreme-subgroups/
cp $ES/FEASIBILITY_extreme_subgroups_HTA.docx      $H/docs/feasibility/v1/
```

The `.tex` file is **not** copied. It is the `keep-tex` output and is regenerated on render.

B2. Gate: each of the 7 copied files has the same `shasum -a 256` hash in both repos.

B3. Update the `Layout` list in `$H/README.md`. Add `presentations/2026-extreme-subgroups/` with the description "extreme-subgroups talk decks (Beamer source + PDF; pptx + builder)". Add `docs/feasibility/v1/` with the description "v1 feasibility assessment (Sept 2026, superseded by v2)". Change no other README text.

B4. Commit in `hta-subpopulations`:

```
git add presentations/ docs/feasibility/v1/ README.md
git commit -m "Add extreme-subgroups talk decks and v1 feasibility assessment

Moved from forestsearch quarto/extreme_subgroups/ (tag
extreme-subgroups-pre-reorg): the Beamer deck source, its PDF and
sim_flow figures, the pptx deck and its builder, and the September v1
feasibility .docx."
```

## Part C — Reorganize `quarto/extreme_subgroups/` in forestsearch (gate)

C1. Rename the kept documents with `git mv`, from `$ES/fixed_random/` to `$ES/`. Do the same for their HTML so the rendered results stay on disk under the new names; C4 then untracks them.

| From (`fixed_random/`) | To (`quarto/extreme_subgroups/`) |
|---|---|
| `extreme_subgroups.qmd` / `.html` | `gbsg_survival_resample.qmd` / `.html` |
| `gbsg_extreme_subgroups_fixed.qmd` / `.html` | `gbsg_survival_fixed.qmd` / `.html` |
| `extreme_subgroups_design_comparison.qmd` / `.html` | `gbsg_design_comparison.qmd` / `.html` |
| `actg175_extreme_subgroups_continuous_fixed.qmd` / `.html` | `actg175_continuous_fixed.qmd` / `.html` |

C2. Delete everything else in the folder that was tracked at the tag, with `git rm`:

- Root:
  - `extreme_subgroups{,_short,_v1,_v2,_fixed,_design_comparison}.qmd`
  - `extreme_subgroups{,_fixed,_design_comparison}.html`
  - `extreme_subgroups_deck.{pdf,pptx}`
  - `extreme_subgroups_section_6_6_4.{pdf,pptx}`
  - `fixed_random.zip`
  - `forestsearch_fixed_baseline_slides.pptx`
  - `build_deck.js`
  - `FEASIBILITY_extreme_subgroups_HTA.docx`
- `fixed_random/`:
  - `extreme_subgroups_fixed.{qmd,html}` (duplicate of `gbsg_`)
  - `extreme_subgroups_fixed_wrapper.{qmd,html}`
  - `extreme_subgroups_slides.{qmd,pdf,tex}`
  - `sim_flow_part{1,2}.pdf`
  - all of `legacy/`

After C1 and C2, `fixed_random/` must be empty. Remove the directory.

C3. Update the old basenames in the **four kept `.qmd` files only**. Apply these replacements in this order, longest first. Leave everything else as is, including `/media/larryleon/bigstorage/...` paths and `extreme_sims_*` payload names.

1. `actg175_extreme_subgroups_continuous_fixed` → `actg175_continuous_fixed`
2. `gbsg_extreme_subgroups_fixed` → `gbsg_survival_fixed`
3. `extreme_subgroups_design_comparison` → `gbsg_design_comparison`
4. `extreme_subgroups_fixed_wrapper` → (should not occur; report if it does)
5. `extreme_subgroups_fixed.qmd` → `gbsg_survival_fixed.qmd`, and the same for `.html`
6. `extreme_subgroups.qmd` → `gbsg_survival_resample.qmd`, and the same for `.html`

Gate: in the four files, `grep -nE "extreme_subgroups(_fixed|_design_comparison|_fixed_wrapper)?\.(qmd|html)|gbsg_extreme_subgroups|actg175_extreme_subgroups"` returns nothing. In `git diff -M -U0` for the four files, every `+` content line contains one of the four new basenames.

C4. Stop tracking the rendered HTML but keep it on disk:
`git rm --cached $ES/*.html`

C5. Write `$ES/README.md`:

```markdown
# Extreme-subgroup simulation documents

Research documents for the behaviour of subgroup estimates under a uniform
treatment effect. They use `run_subgroup_sims()`, `subgroup_cox()` /
`subgroup_glm()` and `compare_subgroup_sims()`. The published pkgdown article
is a separate, adapted copy: `vignettes/articles/extreme_subgroups.qmd`.

| Document | Data / endpoint | Baseline design |
|---|---|---|
| `gbsg_survival_resample.qmd` | GBSG, recurrence-free survival | random-X (`baseline = "resample"`) |
| `gbsg_survival_fixed.qmd` | GBSG, recurrence-free survival | fixed-X (`baseline = "fixed"`) |
| `gbsg_design_comparison.qmd` | reads the two GBSG payloads | fixed vs random comparison (no simulation) |
| `actg175_continuous_fixed.qmd` | ACTG 175 (arms 1 vs 3), CD4 change | fixed-X |

## Rendering

The documents default to 10,000 replicates. Rendering is a compute step.

    mkdir -p results
    quarto render gbsg_survival_resample.qmd            # writes results/extreme_sims_resample_10000_payload.rds
    quarto render gbsg_survival_fixed.qmd               # writes results/extreme_sims_fixed_10000_payload.rds
    quarto render gbsg_design_comparison.qmd            # reads both payloads
    quarto render gbsg_survival_fixed.qmd -P n_sims:50  # quick check

`results/`, `*.rds` and rendered `*.html` are git-ignored. Earlier versions,
renders, slides and the August archive are in git history at tag
`extreme-subgroups-pre-reorg`. Slides and the v1 feasibility document now
live in the private `hta-subpopulations` repository.
```

C6. In `.gitignore`, replace `quarto/extreme_subgroups/fixed_random/results/` with these three lines:

```
quarto/extreme_subgroups/results/
quarto/extreme_subgroups/*.html
quarto/extreme_subgroups/*.rds
```

C7. In `dev/accept_phase2_payload_compare.R` (L21, L27), `dev/accept_phase3_memo_compare.R` (L12) and `dev/diag_phase2_timing.R` (L26), replace `quarto/extreme_subgroups/fixed_random/results` with `quarto/extreme_subgroups/results`. Change nothing else.

C8. Gates:

- `git ls-files quarto/extreme_subgroups` lists exactly five files: `README.md` and the four `.qmd`.
- The four `.html` files exist on disk and `git check-ignore` reports each one as ignored.
- `quarto inspect` succeeds on each of the four `.qmd`. It parses only and must not execute code.
- `git grep -n "extreme_subgroups/fixed_random" -- . ':!dev/tasks' ':!dev/reports'` returns only the two `R/` comment lines (`R/compare_subgroup_sims.R:4`, `R/run_subgroup_sims.R:6`).
- `git diff --cached --stat -- R/` is empty.
- `du -sh quarto/extreme_subgroups` counting tracked files only (`git ls-files … | xargs du -ch | tail -1`) is under 1 MB.

C9. Commit in forestsearch:

```
git add -A quarto/extreme_subgroups .gitignore dev/accept_phase2_payload_compare.R dev/accept_phase3_memo_compare.R dev/diag_phase2_timing.R
git commit -m "quarto/extreme_subgroups: reorganize to four current documents

Keeps the four current wrapper-based documents, flattened and renamed
dataset_endpoint_design: gbsg_survival_resample, gbsg_survival_fixed,
gbsg_design_comparison, actg175_continuous_fixed. Removes byte-identical
duplicates (fixed_random/extreme_subgroups_fixed, legacy/, deck copies of
quarto/guides/), superseded August versions and fixed_random.zip. Rendered
HTML is now ignored (74 MB -> <1 MB tracked). Slides and the v1 feasibility
document moved to the private hta-subpopulations repository. Prior state is
tagged extreme-subgroups-pre-reorg. dev/ payload paths updated; R/ untouched.
The pkgdown article (vignettes/articles/extreme_subgroups.qmd) is unaffected."
```

## Report back

- Tag hash (A4) and the commit hashes for Parts 0, B4 and C9, with the repo for each
- B2 hashes
- C3 replacement counts per file, plus any occurrence of rule 4
- C8 gate outputs
- One line per failed gate, if any
