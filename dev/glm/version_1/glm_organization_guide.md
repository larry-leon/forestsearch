# GLM Extension — File Organization Guide

**Context:** forestsearch R package, GitHub Desktop, RStudio, `master` branch  
**Question:** Where does the GLM extension work live while in development?

---

## Recommendation: Git Feature Branch, Not a Gitignored `dev/` Directory

The short answer is: **do not use a gitignored `dev/` directory for this work**.
Here is why, and what to do instead.

---

## Why Gitignored `dev/` Is the Wrong Tool Here

A gitignored local directory is appropriate for:
- Scratch calculations you truly do not care about recovering
- Personal notes that should never be shared
- Large raw data files that cannot be committed

It is the wrong tool for the GLM extension because:

| Concern | Gitignored `dev/` | Feature branch |
|---|---|---|
| Backup if laptop fails | ❌ None | ✅ Pushed to GitHub |
| Full version history | ❌ None | ✅ Complete commit log |
| `master` stays CRAN-ready | ✅ | ✅ |
| Can share / review later | ❌ | ✅ |
| RStudio project stays clean | ✅ | ✅ |
| Can merge incrementally | ❌ | ✅ |

Losing a gitignored working directory to disk failure or an accidental `rm -rf`
means losing all the GLM work with no recovery path.

---

## The Correct Pattern: Feature Branch

### Step 1 — Create the branch in GitHub Desktop

1. Open GitHub Desktop with the `forestsearch` repo active
2. Click the **Current Branch** dropdown (top center) → **New Branch**
3. Name it: `feature/glm-extension`
4. Based on: `master` (current)
5. Click **Create Branch**

You are now on `feature/glm-extension`. The `master` branch is untouched.

### Step 2 — Directory layout on the branch

```
forestsearch/
├── R/                          ← CRAN source files (survival, unchanged on this branch initially)
│   ├── forestsearch_main.R
│   ├── subgroup_search.R
│   ├── subgroup_consistency.R
│   └── ... (46 existing files)
│
├── dev/                        ← NEW: tracked on this branch, gitignored on master
│   └── glm/
│       ├── README.md           ← design notes, decisions log
│       ├── glm_extension_design.md
│       ├── scratch/            ← throwaway exploration scripts (can add to .gitignore)
│       │   └── explore_convergence.R
│       └── tests/              ← informal test scripts (before formal testthat)
│           ├── test_binary_estimator.R
│           └── test_continuous_estimator.R
│
├── R/                          ← As GLM files graduate from dev/, they move here:
│   ├── glm_effect_estimators.R    (Phase 1 complete → move here)
│   ├── grf_subg_harm_glm.R        (Phase 2 complete → move here)
│   ├── glm_simulate.R             (Phase 3 complete → move here)
│   └── plot_glm_results.R         (Phase 4 complete → move here)
│
├── tests/testthat/             ← Formal tests move here when a file graduates
│   └── test-glm_effect_estimators.R
│
└── vignettes/                  ← GLM vignette added here in Phase 4
    └── glm_outcomes.Rmd
```

### Step 3 — The graduation workflow

Files move from `dev/glm/` to `R/` when they are:
- Fully roxygen-documented
- Passing `devtools::check()` locally
- Covered by at least a smoke test

This keeps `R/` clean and CRAN-check-passable at every stage of the branch.

### Step 4 — Keep `master` unaffected

While on `feature/glm-extension`, any commits to `R/` files that touch
survival infrastructure (e.g., threading `outcome_type` through `forestsearch_main.R`)
will only exist on this branch. The `master` branch sees none of these changes
until you deliberately merge.

To switch back to the CRAN-ready state at any time:

```
# GitHub Desktop: click Current Branch → master
# RStudio Build tab: devtools::load_all() refreshes to master state
```

### Step 5 — Pushing the branch to GitHub

In GitHub Desktop, after your first commit on the branch:
- Click **Publish branch** (or **Push origin** on subsequent commits)

This gives you:
- Full off-machine backup
- A browsable commit history at `github.com/larry-leon/forestsearch/tree/feature/glm-extension`
- The ability to open a Draft PR against `master` when ready for final review

---

## `.gitignore` Additions (Minimal)

You do not need to gitignore `dev/` on the feature branch — you *want* it tracked
there. The only addition useful on either branch:

```gitignore
# In dev/glm/scratch/ — truly throwaway exploration
dev/glm/scratch/
```

If you want `dev/` to be invisible on `master` (so CRAN check never sees it),
add a `.Rbuildignore` entry — not a `.gitignore` entry:

```
# .Rbuildignore (already in package, add this line)
^dev$
```

`R CMD check` and `devtools::check()` already ignore `dev/` by default if it is
listed in `.Rbuildignore`. This ensures the directory never appears in the built
tarball.

---

## `.Rbuildignore` vs `.gitignore`: The Distinction

| File | Purpose | What it controls |
|---|---|---|
| `.gitignore` | Hides files from Git tracking | What Git sees |
| `.Rbuildignore` | Excludes files from the R package tarball | What `R CMD check` / CRAN see |

`dev/` should be in `.Rbuildignore` (so CRAN never sees it) but **not** in
`.gitignore` on the feature branch (so Git tracks and backs it up).

---

## Handling the Four Output Files from the Design Session

The four files already created (`glm_extension_design.md`,
`glm_effect_estimators.R`, `grf_subg_harm_glm.R`, `glm_simulate.R`) should:

1. Be downloaded from the artifact panel
2. Placed in `dev/glm/` (design doc + stubs) after creating the branch
3. `glm_effect_estimators.R`, `grf_subg_harm_glm.R`, and `glm_simulate.R` are
   Phase 1–3 stubs — they live in `dev/glm/` until each passes a full
   `devtools::check()` pass, then move to `R/`

---

## Summary

```
feature/glm-extension branch
├── dev/glm/          ← tracked on branch, excluded from tarball via .Rbuildignore
│   ├── design docs
│   ├── stub R files (pre-graduation)
│   └── scratch/      ← optionally gitignored within dev/
└── R/                ← graduated, CRAN-compliant files move here
```

`master` branch stays: 0 errors, 0 warnings, 1 note — CRAN-ready at all times.
