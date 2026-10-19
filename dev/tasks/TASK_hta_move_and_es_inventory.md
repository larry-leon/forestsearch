# TASK — Move HTA manuscript work to `hta-subpopulations`; inventory `quarto/extreme_subgroups/`

**Date:** 2026-10-19
**Repos:** `~/Documents/GitHub/forestsearch` (source), `~/Documents/GitHub/hta-subpopulations` (new, private)
**R/ impact:** none. No package code is read for changes, moved or edited.
**Push:** never. Larry pushes both repos from GitHub Desktop.
**Gates:** stop-on-failure. A green report needs no chat round trip.

---

## Part 0 — Record this task

1. In `forestsearch`: `cp ~/Downloads/TASK_hta_move_and_es_inventory.md dev/tasks/`
2. Commit: `git add dev/tasks/TASK_hta_move_and_es_inventory.md && git commit -m "tasks: HTA move + extreme_subgroups inventory"`

## Part A — Preflight (gate)

Stop and report if any check fails.

- A1. `~/Documents/GitHub/hta-subpopulations` exists and is a git repo; `git -C … remote get-url origin` resolves to `larry-leon/hta-subpopulations`. If `gh` is available, `gh repo view larry-leon/hta-subpopulations --json isPrivate` returns `true`.
- A2. `hta-subpopulations` working tree is clean (`git status --porcelain` empty).
- A3. In `forestsearch`, `quarto/hta_manuscript/` is tracked at HEAD: `git ls-files quarto/hta_manuscript` lists the `.qmd`. `git status --porcelain quarto/hta_manuscript` is empty.
- A4. Record the forestsearch branch name and HEAD short hash. These go in the report and the commit messages.

## Part B — Seed `hta-subpopulations` (gate)

Target layout for this step only:

```
hta-subpopulations/
  README.md
  .gitignore
  docs/feasibility/FEASIBILITY_HTA_subpopulations_v2.qmd
  docs/feasibility/FEASIBILITY_HTA_subpopulations_v2.docx
```

B1. Copy the files. **Disposition: track the `.qmd` and `.docx` but not the `.html`.** Co-authors read Word, and the 2 MB self-contained HTML is regenerated on render.

```
mkdir -p ~/Documents/GitHub/hta-subpopulations/docs/feasibility
cp quarto/hta_manuscript/FEASIBILITY_HTA_subpopulations_v2.qmd  ~/Documents/GitHub/hta-subpopulations/docs/feasibility/
cp quarto/hta_manuscript/FEASIBILITY_HTA_subpopulations_v2.docx ~/Documents/GitHub/hta-subpopulations/docs/feasibility/
```

B2. Write `README.md`. If GitHub Desktop already created one, replace it.

```markdown
# hta-subpopulations

Manuscript and analyses on the informativeness of small PICO subpopulations
in health technology assessment (HTA): trial-specific simulation and
calibration for requested subpopulations under the EU Joint Clinical
Assessment and national HTA processes.

**Status:** private, pre-submission. Do not share outside the author group.

## Layout

- `docs/feasibility/` — feasibility assessment and venue plan (v2, Oct 2026)
- `manuscript/` — manuscript source (to come)
- `analyses/gbsg/`, `analyses/actg175/` — self-contained analysis documents (to come)
- `presentations/` — talks, including the PSI 2026 precursor (to come)

## Software

Analyses use the `forestsearch` R package, installed from GitHub:

    pak::pak("larry-leon/forestsearch@feature/glm-extension")

No package code lives in this repository.

## Working conventions

- Work on a branch and open a pull request into `main`.
- Commit `.qmd` sources. Rendered `.docx` may be committed for co-author
  reading. Rendered `.html` and `_files/` are not committed.
```

B3. Write `.gitignore`. Merge with any R template GitHub Desktop created, keeping both sets of lines.

```
# Quarto
*.html
*_files/
_freeze/
.quarto/
# R
.Rhistory
.RData
.Rproj.user/
# OS
.DS_Store
```

B4. Gate: run `sha256sum` on both copied files in both repos. The hashes must be identical.

B5. Commit in `hta-subpopulations`. Do not push.

```
git add README.md .gitignore docs/feasibility/
git commit -m "Seed repository: feasibility assessment v2

Moves the HTA manuscript feasibility assessment (v2, 18 Oct 2026) from
forestsearch quarto/hta_manuscript/ at <branch>@<short-hash>. Adds README
with planned layout and .gitignore (rendered HTML and Quarto caches
excluded; .docx retained for co-author reading)."
```

## Part C — Remove from `forestsearch` (gate)

C1. `git rm -r quarto/hta_manuscript`

C2. Gate: run `git grep -n "hta_manuscript" -- . ':!dev/tasks'`. It should return nothing. If it finds references, **list them in the report and do not edit them.**

C3. Commit in `forestsearch`. Do not push.

```
git commit -m "docs: move quarto/hta_manuscript to hta-subpopulations repo

The HTA manuscript feasibility assessment now lives in the private
repository larry-leon/hta-subpopulations (docs/feasibility/). No R/ changes."
```

## Part D — Inventory `quarto/extreme_subgroups/` (read-only)

Do **not** move, rename, delete or render anything in this part.

Write `REPORT_extreme_subgroups_inventory.md` in the same directory as the existing `REPORT_*` files. Find that directory with `git ls-files | grep REPORT_ | head`. If none exist, use `dev/reports/`. The report must contain the following sections.

D1. **File table.** List every file under `quarto/extreme_subgroups/`, tracked and untracked, recursively. Columns:
- relative path
- size
- tracked / untracked / ignored
- date and subject of the last commit
- category, one of: analysis .qmd · slides (revealjs/pptx/beamer format or "slides" in the name) · R script · data (.rds/.RData/.csv) · rendered output (.html/.pdf/.docx/`_files/`) · freeze/cache · other

D2. **Per-.qmd detail.** For each `.qmd`, give:
- `title` and `format:` from the YAML
- inputs it reads: `readRDS`, `load`, `source`, `read.*` paths, and any reads from `dev/`
- outputs it writes: `saveRDS`, `ggsave`, and similar
- whether it runs simulations, with the replicate counts (`n_sims*` values) and the `baseline` setting

D3. **External references.** List every place outside the folder that points into `quarto/extreme_subgroups/`. Check `_pkgdown.yml`, `vignettes/`, `README*`, `.Rbuildignore`, any `_quarto.yml`, `.github/workflows/`, `docs/` and `NEWS.md`. Use `git grep -n "extreme_subgroups"`. Show the matching lines.

D4. **Weight.** Give:
- the total size of the folder
- files over 5 MB, and whether each is tracked
- `.qmd` files that share a title or differ only by suffix (e.g. `_fixed`, `_v2`)

D5. **Vignette linkage.** State which document(s) the forestsearch site or vignette actually publishes, and by what mechanism (pkgdown article, copied HTML, link).

Commit the report in `forestsearch`. Do not push.

```
git commit -m "dev: inventory of quarto/extreme_subgroups (read-only)"
```

## Report back

Report the following:
- the commit hashes for Parts 0, B5, C3 and D, and the repo each one is in
- the B4 hashes
- the C2 grep result
- the path to the D report
- one line per failed gate, if any
