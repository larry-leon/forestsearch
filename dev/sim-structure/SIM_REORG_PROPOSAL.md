# Structure proposal — `quarto/simulations/gbsg_redux`

Response to `dev/sim-structure/CC_BRIEF_sim_reorg.md`.

Nothing under `quarto/` was modified. No file was moved, renamed, or deleted.
Every number below is measured on the working tree at `e9eb514`
(branch `feature/mr-in-replicates`); the commands are reproducible from the
directory itself.

Two of the brief's premises did not survive measurement, and both change the
design. They are §2.7 and §2.8. Read those before the layout.

---

## 1. Inventory

### 1.1 What is actually tracked

The brief counts 204 top-level entries. That is the top level only; the
directory also contains `mr_sweep/`, which git tracks in full.

| | files | disk |
|---|---:|---:|
| `.qmd` (top level) | 73 | 4.6 M |
| `.html` (top level) | 69 | 182 M |
| `.rds` (top level) | 59 | 2.5 M |
| `.rds` (under `mr_sweep/`) | 198 | 31 M |
| `.zip` (under `mr_sweep/`) | 3 | 11 M |
| `betaHhat_truth.R` | 1 | 9 K |
| `_oracle_cells_kappa.md` | 1 | 4 K |
| **git-tracked total** | **404** | **232 M** |

Plus one untracked directory, `_sim_mr_coverage_h10_knoise_files/` (matched by
the `*_files/` rule in `.gitignore`).

Blob bytes already written to history for this path:

```
qmd     9.0 MB
html  200.7 MB
rds    32.6 MB
other  11.1 MB
```

### 1.2 The 73 `.qmd`, grouped by role

| role | count | what it is |
|---|---:|---|
| **`t1_t2` family** | **60** | `fs_` 49, `grf_` 6, `dina_` 4, `fsNEW_` 1. One 1231-line document, 60 times. |
| `mr_coverage_sweep_*` | 8 | A **different and better** machinery — see §1.4. Not superseded. |
| `sim_fs_maxcons_fb_mr_*` | 2 | The newest generation: new naming, `mr_draws`, `focus_tag` guard, `sg_focus = "maxcons"`. |
| `_sim_*` | 3 | Quarto **include-fragments** (leading `_` = excluded from project render), not standalone probes. See §2.6. |

### 1.3 The duplication, quantified

The brief measured siblings within a cell. The duplication is much wider than
that: it crosses cells *and* methods.

```
diff, differing lines (of 1231)

  batch siblings, same cell .......................  2–4
  fs_h15_knoise3_n1500  vs  fs_h15_knoise0_n500 ...  8      (different cells)
  fs_h10_knoise0_n500   vs  grf_h10_knoise0_n500 ...  2      (different methods)
  fs_h10_knoise0_n500   vs  dina_h10_knoise0_n500 ..  2
  fs_h10_knoise0_n500   vs  fsNEW_h10_knoise0_n500 .  4
```

The eight lines separating two different design cells are, in full:

```
95   n_sims       <- 50L        vs  100L
101  k_random_noise <- 3        vs  0
117  sim_id_start <- 151L       vs  101L
122  n_sample     <- 1500L      vs  500L
```

So **all 60 documents are one document.** They differ in five assignments —
`subgroup_method`, `target_hr_harm`, `k_random_noise`, `n_sample`, and the
`sim_id_start`/`n_sims` pair. Everything else, including `rds_stem`,
`rds_path`, `combine_glob` and `sim_id_end`, is already derived from those five.
That is ~74,000 lines of `.qmd` where ~1,240 plus a parameter table would do.

### 1.4 The batching machinery

```r
rds_stem      <- sprintf("%s_t1_t2_m1_h%02d_knoise%d_n%d",
                         method_tag, round(10*target_hr_harm), k_random_noise, n_sample)
sim_id_end    <- sim_id_start + n_sims - 1L
rds_path      <- sprintf("%s_res_%d_%d.rds", rds_stem, sim_id_start, sim_id_end)
combine_glob  <- sprintf("%s_res_*.rds", rds_stem)
combine_files <- NULL          # explicit override; NULL -> use combine_glob
combined_path <- NULL          # NULL -> "<stem>_combined_<min>_<max>.rds"
save_rds      <- identical(run_mode, "batch") && !is.null(rds_path)
```

`run_mode` takes exactly two values, `"batch"` and `"combine"`, both dispatched
off `identical()` at lines 248, 536 and 615. `run_smoke` (line 161) is a master
short-circuit. `combined_path` is deliberately chosen outside `combine_glob` so
a re-combine never re-ingests its own output — a good decision, and one the
proposal preserves.

**The important property: the `.rds` name is a pure function of the parameters.
The `.qmd` name is not.** It is assigned by hand and drifts (§2.3).

### 1.5 The agreement check

Combine mode gates on five fields and refuses to pool if any disagree:

```r
for (k in c("n_sample", "nb_boots", "mr_draws", "subgroup_method", "seed_base"))
```

plus a `sim_id` disjointness check that `stop()`s on any repeat, and a
`truth`-identity check that only `warning()`s. Provenance fields (`sg_focus`,
`focus_tag`, `consistency_method`, `forestsearch_version`, `n_workers`,
`parallel_mode`) are recorded and **not** gated, with the reason stated inline
at lines 678–686: a missing field resolves to `NA`, `NA` counts as a distinct
value, so gating would `stop()` on every pre-existing batch. That decision is
correct and this proposal does not touch it.

The newer `sim_fs_maxcons_*` documents add an NA-tolerant `focus_tag` guard
(unique over non-`NA` values only) — the right pattern, and the one any future
field should follow.

### 1.6 Result storage, and what is terminal

Read from every `.rds` in the directory:

| | files | terminal? |
|---|---:|---|
| `*_res_a_b.rds` | 50 | input to combine mode only |
| `*_combined_a_b.rds` | 8 | **terminal** — read by `dev/sg-focus-work/` (§2.8) |
| `mr_sweep/<run_tag>/*_res.rds` | 198 | input to the sweep's own grid accumulation |

`mr_sweep/` is a genuinely different design and a good one: **one** document
renders the whole `subgroup_method × n_sample` grid, is resumable (a cell is
skipped if its bundle exists), writes per-cell bundles into
`mr_sweep/<run_tag>/`, and accumulates `mr_coverage_grid_<run_tag>.rds`. The
pattern this proposal recommends for the `t1_t2` family already exists in the
same directory.

---

## 2. Findings

### 2.1 Drift is not a risk; it has already happened, in the files that produce the reported numbers

The κ / β(H) enhancement — a 25-line `beta-true-region` chunk, suppression of
the oracle's `b_betaHhat` / `Cov_betaHhat` cells, and an identification-factor
footnote — is present in exactly **4 of the 10 combine-mode documents**:

| cell | combine doc has κ? |
|---|---|
| `dina_h10_knoise0_n500` | yes |
| `fs_h10_knoise0_n500` | yes |
| `fs_h15_knoise0_n1000` | yes |
| `grf_h10_knoise0_n500` | yes |
| `fs_h10_knoise0_n1000` (×2 docs) | **no** |
| `fs_h15_knoise0_n500` | **no** |
| `fs_h15_knoise3_n1500` | **no** |
| `fs_h20_knoise3_n1500` (×2 docs) | **no** |

Six pooled cells therefore render a *different table for the same estimand*
than the other four. `_oracle_cells_kappa.md` records the intent verbatim:

> File modified: `fs_t1_t2_m1_h10_knoise0_n500_combine_1_500.qmd` — **that file
> only**. Not `betaHhat_truth.R`, not any `*_batch_*.qmd` […]

It then propagated to three of nine siblings and stopped. This is the brief's
argument made concrete, and it is the single strongest reason to act.

### 2.2 The frozen-cell problem is worse than "cannot pool"

Every one of the 58 legacy bundles carries `meta$gate_draws = 5000`. The one new
bundle (`fs_maxcons_fb_mr_..._res_1_100.rds`) carries `meta$mr_draws = 5000`.
All 73 source documents now write `mr_draws`; **zero** still write `gate_draws`.

Because the gate reads `mr_draws`, a legacy-only pool sees `NA` for every batch
— one unique value — and passes. Mix in one new batch and the set becomes
`{NA, 5000}` and it `stop()`s.

The consequence is not "some pooling is awkward". It is that **all eight
existing cells are frozen at their current replicate counts.** Each can be
re-combined from what exists; none can be extended by one more batch without
re-running the cell from scratch. Any layout that implies these are live,
appendable cells would be lying about them.

### 2.3 Filenames do not describe the files

| claim in name | measured content |
|---|---|
| `fs_..._h20_knoise3_n1500_combine_1_100.qmd` | `run_mode <- "batch"`, `sim_id_start = 1` |
| `..._combine_101_150`, `_151_200`, `_201_250`, `_251_350`, `_351_400` | all `run_mode <- "batch"` |
| `..._h20_knoise3_n1500_combine_1_250.qmd` | combine mode, but `combine_files <- NULL` → the glob takes all 7 batches → 500 rows, not 250 |
| `..._h10_knoise0_n1000_combine_1_400.qmd` | same: cannot produce a 1–400 pool |
| `..._h15_knoise3_n1500_combine_1_500.qmd` | combine mode with stale `sim_id_start = 351, n_sims = 50` (inert here, but unverifiable from the name) |

**Seven of the eight files named `..._h20_knoise3_n1500_combine_*` are batch
runs.** During migration, the filename must be treated as a label of unknown
accuracy; the `run_mode` / `sim_id_start` / `n_sims` triple in the file is the
only truth. §5 supplies the extracted table.

### 2.4 The 182 M of `.html` is unreferenced boilerplate

No `.qmd`, `.md`, `.yml` or `.R` anywhere in the repo links to any of them. Each
is ~2.8 M because `quarto/_quarto.yml` sets `embed-resources: true`, inlining
MathJax and Bootstrap — so the size is essentially independent of the content.
59 of the 69 are `t1_t2`-family renders whose unique content is one smoke-test
table.

### 2.5 Two artifacts have no source that reproduces them

* **`fs_t1_t2_m1_h20_knoise3_n1500_combined_1_250.rds`** (250 rows). Both h20
  combine documents have `combine_files <- NULL`; the glob now returns all seven
  batch files → 500 rows. Nothing on disk regenerates a 250-row pool.
* **`fs_maxcons_fb_mr_m1_h10_knoise0_n500_res_1_100.rds`** (100 rows,
  `mr_draws = 5000`, production `nb_boots`). Its only source,
  `sim_fs_maxcons_..._batch_1_100.qmd`, currently sits at `n_sims = 50L`,
  `nb_boots = 20L`, `mr_draws = 500L` — a smoke configuration that would write
  `..._res_1_50.rds`. This is also the file `dev/sim-check/fixcheck/` reads.

Both are **fixed points**: preserve, never regenerate. Conversely
`fs_t1_t2_m1_h15_knoise0_n1000_batch_1_20.qmd` has no `res_1_20.rds` — a smoke
stub that produced nothing kept.

### 2.6 The `_sim_*` files are include-fragments, not superseded probes

The brief groups them with `fsNEW_` as supersedes. They are not. The leading
underscore is Quarto's exclude-from-render marker; each is a manuscript figure
or table fragment that reads only the pooled `mr_coverage_grid_*.rds` and runs
in seconds. They are the *consumers* of the sweep, not competitors to it.

However — **no document in the main tree includes them.** A repo-wide search for
`{{< include >}}` referencing them returns nothing outside
`.claude/worktrees/`. Either the including manuscript was never committed, or
they are dead. That is a question of intent, not evidence (§6).

### 2.7 The brief's `.html` framing understates one thing and overstates another

Untracking the renders **stops the growth but does not shrink `.git`.** The
~201 M already in history stays there absent a `filter-repo` rewrite, which this
proposal does not recommend and which the CLAUDE.md history rule forbids without
explicit approval. Stated plainly so the saving is not oversold.

### 2.8 §2.5 of the brief does not hold — and this is what frees the design

The brief names two cross-directory consumers and calls them "the main
constraint on the design". I checked both:

* **`quarto/GuoHe/guohe_sec52.qmd`** reads only
  `guohe_repro_t7_beta2_*.rds`, `guohe_sec52_truth_beta2_*.rds` and
  `sens_B500/*` — all from its own directory. It never opens a `gbsg_redux`
  file.
* **`quarto/simulations/gbsg/fs_subgroup_id_sweep_n500.qmd:266`** reads
  `det$covs` from a `results` object built from `fs_idsweep_*.rds` under its own
  `cache_dir`.

Both are references to *recorder-schema column names* (`n_sel`, `covs`) — which
is why they matter for the schema, correctly flagged in §4 of the brief — not to
`gbsg_redux` artifacts. A repo-wide grep for `gbsg_redux` outside the directory
returns **no hit anywhere under `quarto/`**.

The real cross-directory consumers are three lines in `dev/`:

```
dev/sg-focus-work/R/phase0_cells.R:53      SIM_DIR <- "quarto/simulations/gbsg_redux"
dev/sg-focus-work/R/phase0_cells.R:56      sprintf("fs_t1_t2_m1_%s_combined_1_500.rds", cell)
dev/sg-focus-work/R/phase0_summarise.R:34  file.path(SIM_DIR, sprintf("fs_t1_t2_m1_%s_combined_1_500.rds", cell))
dev/sim-check/fixcheck/test_combine_guard.R:34  ../../../quarto/simulations/gbsg_redux/fs_maxcons_..._res_1_100.rds
```

One path pattern over six files, plus one single-file reference. That is a
three-line edit, not a design constraint — and it means the layout is free in
exactly the dimension the brief assumed it was pinned.

### 2.9 Side issue, outside this scope

`quarto/_quarto.yml` makes `quarto/` a Quarto **project** with `output-dir`
commented out. A bare `quarto render` in `quarto/` would attempt all 73 of these
documents at their committed `n_sims`. Nothing prevents that today. Worth a
`project: render:` allowlist, but it is a separate change and I have not made it.

---

## 3. Proposed layout

### 3.1 Target

```
quarto/simulations/gbsg_redux/
├── betaHhat_truth.R                       # unchanged, stays at this level
├── cells.yml                              # NEW — the cell registry
├── run_batch.sh                           # NEW — the driver
│
├── sim_fs_gate_fb_mr_m1.qmd               # the t1_t2 template, params-driven   (60 -> 1)
├── sim_fs_maxcons_fb_mr_m1.qmd            # the maxcons generation, params-driven ( 2 -> 1)
├── mr_coverage_sweep.qmd                  # already run_tag-driven               ( 8 -> 1)
├── _fig_mr_coverage.qmd                   # fragments, renamed for what they are
├── _fig_mr_selection_stability.qmd
│
├── results/                               # tracked — every .rds, names UNCHANGED
│   ├── fs_t1_t2_m1_h10_knoise0_n500_res_1_20.rds
│   ├── fs_t1_t2_m1_h10_knoise0_n500_combined_1_500.rds
│   └── …  (59 files, byte-identical, git mv only)
│
├── renders/                               # .gitignored
│   └── combined/                          # tracked exception — 10 combine renders
│
└── mr_sweep/<run_tag>/                    # untouched
```

Top level goes from 205 entries to 8 files and 3 directories.

**The documents stay at this level.** That is deliberate: `source("betaHhat_truth.R")`
appears in 70 of the 73 documents and resolves against the document's own
directory. Any layout that pushes the `.qmd` into subdirectories breaks all 70.
Only the `.rds` move, which costs the two lines in §3.4.

**`.rds` basenames do not change.** `results/` is a prefix, nothing more. The
existing artifacts keep the names that already encode their cell and range, so
the six-file pattern in `phase0_cells.R` needs one `file.path` component added,
not a rewrite.

### 3.2 How a batch is launched

No document here has a `params:` block today (measured: 0 of 73). Add one, and
replace the five literal assignments with reads from it. Nothing downstream
changes — `rds_stem`, `rds_path`, `combine_glob`, `sim_id_end` are already
derived.

```yaml
params:
  subgroup_method: "consistency"
  target_hr_harm:  1.5
  k_random_noise:  3
  n_sample:        1500
  run_mode:        "batch"
  sim_id_start:    1
  n_sims:          50
  nb_boots:        300
  mr_draws:        5000
```

```bash
quarto render sim_fs_gate_fb_mr_m1.qmd \
  -P target_hr_harm:1.5 -P k_random_noise:3 -P n_sample:1500 \
  -P sim_id_start:151 -P n_sims:50 \
  -o renders/fs_t1_t2_m1_h15_knoise3_n1500_batch_151_200.html
```

**This answers the constraint the brief identifies as the reason per-batch
copies were attractive.** Quarto's `X.qmd → X.html` default is overridable per
render: `quarto --version` reports 1.9.38, and `quarto render --help` lists both
`-P, --execute-param` and `-o, --output`. Per-batch copies were never required
for output naming.

One caveat, honestly flagged: I verified the flags **exist**; I did not run a
render, because rendering one of these documents executes the simulation. Before
step 4 of the migration, do one throwaway render at `n_sims: 2` and confirm that
(a) `-o` accepts a path with a directory component alongside `embed-resources`,
and (b) knitr's working directory is still the document's directory so
`source("betaHhat_truth.R")` resolves. If `-o` proves awkward with a directory,
the fallback is render-then-`mv` in `run_batch.sh` — one line, same result.

`cells.yml` + `run_batch.sh` exist so a launch is *recorded* rather than typed:

```yaml
# cells.yml
fs_t1_t2_m1_h15_knoise3_n1500:
  subgroup_method: consistency
  target_hr_harm: 1.5
  k_random_noise: 3
  n_sample: 1500
  batches: [[1,50],[51,50],[101,50],[151,50],[201,50],
            [251,50],[301,50],[351,50],[401,50],[451,50]]
```

```bash
./run_batch.sh fs_t1_t2_m1_h15_knoise3_n1500 151    # one batch
./run_batch.sh fs_t1_t2_m1_h15_knoise3_n1500 combine
```

This registry is the migration's real deliverable. The parameter table in §5 is
the only information the 60 retired documents hold, and `cells.yml` is where it
goes.

### 3.3 One document per cell, or one per batch — and why

**One parameterised document for the whole family**, not one per cell and not
one per batch. The measurement in §1.3 is decisive: cells differ from each other
by four parameter lines and methods by two, so "one per cell" would still be 11
copies of the same file with the same drift mechanism, just at lower multiplicity.
§2.1 shows the mechanism firing at 10 copies already.

This is a judgement call in one respect, and here is the alternative. Per-batch
files have a real property the fold gives up: **the source of a completed run is
frozen alongside its output.** You can `git log` a batch file and see exactly the
document that produced `res_151_200.rds`. Under the fold, that provenance moves
into `cells.yml` plus the template's history, which is a weaker link.

I still recommend the fold, for two reasons. First, §2.2: the cells are frozen
anyway — no existing `.rds` can be extended, so the freeze-the-source property
protects a re-run that cannot happen. Second, the property is not actually held
today: §2.3 shows seven files whose names misdescribe their own `run_mode`, so
the per-file provenance is already unreliable. If the maintainer weights
frozen provenance more highly than I do, the middle option is to fold to one
document per *method* (`fs` / `grf` / `dina`, 3 files) and keep per-cell
parameters in `cells.yml` — 60 → 3 rather than 60 → 1, capturing most of the
benefit.

### 3.4 What moving the `.rds` costs

Two lines in the template:

```r
rds_path     <- file.path("results", sprintf("%s_res_%d_%d.rds", rds_stem, sim_id_start, sim_id_end))
combine_glob <- file.path("results", sprintf("%s_res_*.rds", rds_stem))
```

and three lines in `dev/` (§2.8). `combined_path <- NULL` resolves relative to
the working directory, so it needs the same prefix or an explicit value.

That is the entire cost. Per §2.8 it is not the constraint the brief expected.

### 3.5 What git tracks

**Recommend: `.gitignore renders/`, with `renders/combined/` force-added.**

| option | tree | keeps the record? |
|---|---:|---|
| keep everything | 182 M, +2.8 M per render forever | yes |
| **combines only (recommended)** | **~28 M** | **yes, for the runs that report numbers** |
| ignore all `.html` | 0 M | no |

The tradeoff the brief names is real — a tracked render is a record of what was
actually produced, and this project has needed exactly that twice. The
recommendation keeps that property where it carries information. A combine
render contains the cell's reported estimation and coverage tables; a batch
render contains a smoke-test table plus 2.8 M of inlined MathJax and Bootstrap
(§2.4). Keeping 10 and ignoring 59 preserves the record and drops 85% of the
weight.

Per §2.7, this stops growth; it does not shrink `.git`.

---

## 4. Migration plan

Ordered. Each step leaves the repo working. No step re-runs a simulation, and
no step rewrites an existing `.rds`.

### Step 0 — baseline (nothing moves)

Record `git rev-parse HEAD`. Write
`dev/sim-structure/rds_manifest.csv`: for every `.rds` in the directory, its
path, `nrow(results)`, and every `meta` field. This is the invariant every later
step checks against.

### Step 1 — reconcile the drift, before restructuring anything

**This must come first, and it is the only step that changes a reported number.**

Decide whether the κ / β(H) treatment in `_oracle_cells_kappa.md` is correct.

* If yes: apply it to the six combine documents lacking it (§2.1) and re-render
  those six combines.
* If no: remove it from the four that carry it and re-render those four.

Do this while the per-cell files still exist, so the change is auditable one
cell at a time against a committed render. **Do not fold 60 files into one until
you know which version is the one to keep** — otherwise the fold silently elects
a winner and the six divergent cells are lost without a diff.

Working repo: yes. Only combine documents change; no batch artifact is touched.

### Step 2 — move the artifacts

`git mv` the 59 top-level `.rds` into `results/`. Update the two template lines
(§3.4) and the three `dev/` lines (§2.8). Re-run the step-0 manifest and require
it to match exactly.

Nothing renders in this step, and the only readers are the three `dev/` scripts,
both updated in the same commit. Working repo: yes.

### Step 3 — untrack the batch renders

Add `quarto/simulations/gbsg_redux/renders/` to `.gitignore`. `git mv` the 10
combine renders to `renders/combined/` and `git add -f` them. `git rm --cached`
the remaining 59; they stay on disk. Fully reversible; touches no artifact.

### Step 4 — introduce the template, and prove it

Create `sim_fs_gate_fb_mr_m1.qmd` from the **reconciled** combine document
(post-step 1) with the `params:` block from §3.2. Delete nothing yet.

Gate the fold on one check: render the template in `combine` mode for one cell
and diff the `combined_*.rds` it writes against the committed one, bitwise on
every non-timing column. `dev/sim-check/SIM_INTEGRITY_FINDINGS.md` establishes
the precedent and the exemptions (`fb_secs` / `fit_mr_secs`), and its
cross-architecture caveat applies: run the check on one machine.

Only that check licenses step 5.

### Step 5 — retire the copies

`git rm` the 60 `t1_t2` `.qmd`. They remain in history. Commit `cells.yml` and
`run_batch.sh` carrying the exact parameters each retired file held — the table
in §5.2, transcribed from the files rather than from their names (§2.3).

### Step 6 — the sweep family

Same fold for the eight `mr_coverage_sweep_*.qmd` → one document plus a
`run_tag` parameter. The design is already `run_tag`-driven, so this is
parameter extraction with no logic change. `mr_sweep/` is untouched.

### Step 7 — resolve the fragments

Per the maintainer's answer to §6.2.

### Fixed points — preserved, never migrated

* `fs_t1_t2_m1_h20_knoise3_n1500_combined_1_250.rds` — unreproducible (§2.5)
* `fs_maxcons_fb_mr_m1_h10_knoise0_n500_res_1_100.rds` — unreproducible (§2.5)
* all eight `*_combined_1_500.rds` — reproducible only while the glob set is
  unchanged, which step 2 preserves and no later step disturbs

Step 2 relocates these; nothing rewrites them.

---

## 5. Retirement candidates

Listed, not acted on.

### 5.1 Recommended

| candidate | evidence |
|---|---|
| `fsNEW_t1_t2_m1_h10_knoise0_n500_*` (1 qmd, 1 rds, 1 html) | differs from the `fs` template by 4 lines, one of which is `method_tag <- "fsNEW"`. Its `res_1_100.rds` was never combined — no `fsNEW_*_combined_*` exists. Repo-wide grep for `fsNEW` outside the directory: **0 hits** (only its own three files and a Quarto index cache). |
| `fs_t1_t2_m1_h10_knoise0_n1000_combine_1_400.qmd` | cannot produce a 1–400 pool: `combine_files <- NULL`, glob takes all 8 batches. No `combined_1_400.rds` exists. Functionally identical to its `combine_1_500` sibling. |
| `fs_t1_t2_m1_h20_knoise3_n1500_combine_1_250.qmd` | same, for 1–250. The `combined_1_250.rds` it is named for is a fixed point it can no longer produce (§2.5). |
| `fs_t1_t2_m1_h15_knoise0_n1000_batch_1_20.qmd` | smoke stub; no `res_1_20.rds` was ever kept. The four other `_batch_1_20` stubs did produce artifacts and should be kept as cells. |
| 7 of 8 `mr_coverage_sweep_*.qmd` | differ only in `run_tag`; the machinery is already `run_tag`-parameterised. All `mr_sweep/<run_tag>/` outputs are preserved and independent. |
| `mr_sweep/*.zip` (3 files, 11 M) | each sits beside a directory of the same name with matching size (3.4/3.9/4.1 M vs 3.5/4.0/4.2 M) and 22 files vs 23 zip entries. Almost certainly snapshots. **Verify contents match before removing** — I compared sizes and counts, not bytes. |

### 5.2 Not retirement — the parameters to transcribe

The extracted `run_mode` / `sim_id_start` / `n_sims` for all 73 documents is in
this repository's history of this analysis and reproduces with:

```bash
cd quarto/simulations/gbsg_redux
for f in *.qmd; do
  printf "%-58s %-8s %-6s %-6s\n" "$f" \
    "$(grep -m1 '^run_mode'     "$f" | sed -E 's/.*<- *"([a-z]+)".*/\1/')" \
    "$(grep -m1 '^sim_id_start' "$f" | sed -E 's/.*<- *([0-9]+)L.*/\1/')" \
    "$(grep -m1 '^n_sims'       "$f" | sed -E 's/.*<- *([0-9]+)L.*/\1/')"
done
```

Use this, not the filenames, to build `cells.yml` (§2.3).

### 5.3 Explicitly **not** retirement candidates

* **`mr_coverage_sweep_*` as a family.** The brief lists it under supersedes. It
  is the opposite — the newest and best-designed machinery here, and the model
  §3 recommends generalising. Retire seven copies; keep the design.
* **`_sim_*`.** Include-fragments, not probes (§2.6). Their status is §6.2.

---

## 6. Open questions for the maintainer

**6.1 Is the κ / β(H) treatment correct, and should it apply to every cell?**
This blocks step 1 and therefore everything after it. Four of ten combine
documents carry it and six do not (§2.1); the fold cannot proceed until one of
those is the answer. This is the only question that must be settled before work
starts.

**6.2 Are the three `_sim_*` fragments live?** Nothing in the main tree includes
them. Either their manuscript was never committed, or they are dead. Evidence
cannot distinguish these.

**6.3 Should the eight frozen cells ever be extended past their current
replicate counts?** Per §2.2, extending any of them means re-running the whole
cell — the `gate_draws` → `mr_draws` rename made the existing batches
unpoolable with anything new. If the answer is no, `cells.yml` should mark them
`frozen: true` so the layout does not imply an append that would fail at the
agreement check. If yes, budget for the full re-run and the layout should carry
`mr_draws` from the start.

**6.4 What did `fsNEW` mean?** The brief asks and I could not recover it. It is
one file, four lines from the `fs` template, with an orphan artifact and no
external reference.

**6.5 How much does frozen per-file provenance matter?** §3.3 sets out the one
real property the fold gives up and offers a 60 → 3 middle option. My
recommendation is 60 → 1, but the weighting is the maintainer's.

**6.6 Keep the three `mr_sweep/*.zip`?** They appear to duplicate the
directories beside them (§5.1).

---

## Summary

The measured duplication is wider than the brief estimated — 60 documents, not
22 stems, collapsing to one — and the drift it warns about has already produced
a live inconsistency across six of ten pooled cells (§2.1). The constraint the
brief expected to dominate the design, cross-directory `.rds` consumers, is
three lines in `dev/` rather than a hard dependency in `quarto/` (§2.8), and
Quarto's `-P` / `-o` flags dissolve the output-naming problem that made
per-batch copies attractive (§3.2).

The one thing that must happen before any restructuring is §2.1: decide which
version of the combine report is correct, while the per-cell files still exist
to be diffed.
