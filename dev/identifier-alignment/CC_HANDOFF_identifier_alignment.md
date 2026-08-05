---
bibliography: []
---

# Handoff: identifier alignment — implementation and re-analysis

Working directory: `forestsearch/dev/identifier-alignment`
Input files: `forestsearch/dev/identifier-alignment/sim_analyses`
Branch: `mr-in-replicates` (repo default branch is `master`, not `main`).
Do **not** merge to `feature/glm-extension` first: enforcement sites 2 and 3 are
`mr_in_replicates`, which only exists on this branch.

Everything in Phase A is settled and needs no further discussion. Phase C has
one open item flagged as **BLOCKER** that must be resolved with Larry before
the FS simulation is run.

---

## Background in one paragraph

The manuscript requires every identifier's selection map to rank on the
inferential coefficient $\hat\beta(g)$ — the same within-subgroup effect the
trial reports. This is the alignment condition (main §4.2), stated as an iff
for the multiplier-resampling (MR) correction reproducing the León bootstrap to
first order. DINA and GRF currently default to ranking on their own native
statistics (DINA's $\bar\tau$, GRF's doubly-robust score mean), so at package
defaults both are outside the manuscript's scope. Separately, MR requires a
fixed candidate family (§5); FS attains this by disabling its model-based front
ends, whereas DINA and GRF cannot attain it at all — which is itself the
manuscript's subject matter and must **not** be warned about.

---

## The six configurations

| # | `subgroup_method` | sub-mode | ranks on | condition 1 |
|---|---|---|---|---|
| 1 | `consistency` | any `sg_focus` | $\hat\beta$, raw or $\sigma_D$-standardized | pass |
| 2 | `dina` | `dina_select_statistic = "effect"` | $\hat\beta$ | pass |
| 3 | `dina` | `dina_select_statistic = "dina"` | $\bar\tau$ | fail |
| 4 | `grf` | `frontier` + `grf_select_statistic = "effect"` | $\hat\beta$ | pass |
| 5 | `grf` | `frontier` + `grf_select_statistic = "dr"` | mean DR score | fail |
| 6 | `grf` | `tree` (statistic ignored) | policy-tree objective | fail, structurally |

All six are **retained**. #3, #5, #6 become reachable only by explicit opt-in.

---

## Directory layout — verified

```
dev/identifier-alignment/
├── CC_HANDOFF_identifier_alignment.md
├── CC_PROMPT.md
└── sim_analyses/                    <- underscore, not hyphen
    ├── analysis_actg175_binary_multimethod.qmd
    ├── analysis_gbsg_survival_multimethod.qmd
    ├── sim_fs_maxeffCons_fb_mr_m1_h10_knoise0_n500_combine_1_500.qmd
    ├── actg175_table_payload.rds
    ├── gbsg_table_payload.rds
    ├── fs_maxcons_fb_mr_m1_h10_knoise0_n500_combined_1_500.rds
    ├── fs_t1_t2_m1_h10_knoise0_n500_combined_1_500.rds
    ├── dina_t1_t2_m1_h10_knoise0_n500_combined_1_500.rds
    └── grf_t1_t2_m1_h10_knoise0_n500_combined_1_500.rds
```

Four facts that follow, all of which supersede earlier drafts of this handoff:

1. **There is no `results/` directory.** Only the pooled `_combined_1_500.rds`
   bundles are on disk; the per-batch `_res_*.rds` files are not. So the FS
   batch-file collision described in C0 cannot occur on a first render. It
   becomes live the moment a render creates `results/`, so the rule in C0 still
   governs — it is simply not yet triggered.

2. **`actg175_table_payload.rds` is present and IS overwritten.** This collision
   is real and immediate. See B0.

3. **`betaHhat_truth.R` is not in this directory.** The three simulation
   documents `source()` it with a bare relative path. Locate it in the repo and
   either copy it beside them or repoint the `source()` call.

4. **There are two FS reference bundles.** `fs_maxcons_fb_mr_*` is the one to
   compare against — its `meta` records `sg_focus = "maxcons"`, matching the
   simulation documents. `fs_t1_t2_*` is the older run under the previous stem
   and is **not** the comparison target.

Nothing else is in scope. The coverage-sweep documents referenced in earlier
drafts (`_sim_mr_coverage_h10.qmd`, `mr_coverage_sweep_h10_knoise0.qmd`) are not
present here; do not go looking for them.

## Reporting standard

Comparisons are **not** expected to match exactly, and a difference is not a
failure. Resampling seeds, MR draw streams, package version drift, and the
deliberate configuration changes in this task all move numbers. Report:

- what matched, what did not, and by how much;
- for each difference, the most likely cause;
- anything that could not be run, and why.

Do not tune settings to force agreement with the reference. Do not assert
pass/fail thresholds. Where a difference follows from a decision recorded in
this handoff — the ACTG175 front-end change most of all — report it as the
expected consequence of that decision, not as a discrepancy.

## Phase A — package changes

### A1. Default changes

All six arguments are in the `forestsearch()` formals in `R/forestsearch_main.R`.
The last three resolve via `match.arg()` (at `:1126`, `:2067`, `:2068`), so
re-ordering the vector changes the default while leaving the permitted set
unchanged. Nothing is added or removed.

| line | argument | current | change to |
|---|---|---|---|
| 1023 | `use_lasso` | `TRUE` | `FALSE` |
| 1024 | `use_grf` | `TRUE` | `FALSE` |
| 1027 | `use_dina` | `FALSE` | *unchanged* |
| 1031 | `dina_select_statistic` | `c("dina", "effect")` | `c("effect", "dina")` |
| 1035 | `grf_selection` | `c("tree", "frontier")` | `c("frontier", "tree")` |
| 1036 | `grf_select_statistic` | `c("dr", "effect")` | `c("effect", "dr")` |

Update the corresponding `@param` roxygen at `:249` (`use_lasso`) and `:250`
(`use_grf`), both of which currently state "Default TRUE".

### A2. Shared validator

Add `.validate_mr_configuration()` to `R/forestsearch_helpers.R`. It raises a
hard `stop()` — never a warning, never a silent coercion — under either of two
conditions. Each message must name the offending argument and the exact setting
that resolves it.

**Condition 1 failures** (selection does not rank on $\hat\beta$):

- `subgroup_method == "dina"` and `dina_select_statistic == "dina"`
- `subgroup_method == "grf"` and `grf_select_statistic == "dr"`
- `subgroup_method == "grf"` and `grf_selection == "tree"`

**Condition 2 failure, FS only** (family not fixed, but trivially fixable):

- `subgroup_method == "consistency"` and any of `use_lasso`, `use_grf`,
  `use_dina` is `TRUE`

**Do not error** for #2 and #4 on condition 2. Their families regenerate under
resampling by construction; the manuscript studies that deliberately (§6 frames
those results as empirical evaluations rather than instances of Theorem 2).
Warning about it would fire on every DINA and GRF replicate across thousands of
simulations.

**Do not error** when `use_lasso`/`use_grf`/`use_dina` are `TRUE` alongside
`subgroup_method` of `"dina"` or `"grf"`. Verified inert on those paths: the
flags are consumed at `forestsearch_main.R:2199` (GRF cuts), `:2304` (DINA
cuts), and `:2456` (lasso), all downstream of the DINA branch's return at
`:2052` and the GRF branch's return at `:2195`. The only earlier reference,
`:1777`, feeds `interpret_search_config()`, which is diagnostic output only.

### A3. Enforcement sites

Four, of which three are code and one is documentation.

1. `forestsearch()` in `R/forestsearch_main.R`, when `mr_inference = TRUE`
   (formal at `:1103`, defaults `FALSE`).
2. `forestsearch_bootstrap_dofuture()` in `R/bootstrap_dofuture_main.R`, when
   `mr_in_replicates = TRUE`. This is a genuinely separate path: the function
   reconstructs the call from `fs.est$args_call_all` at `:242`, and under
   `mr_in_replicates = TRUE` propagates `mr_inference` unchanged into every
   replicate. A #3/#5/#6 fit would otherwise run a misaligned re-selection
   `nb_boots` times. Under the default `FALSE` the flag is stripped, so no
   error is needed there.
3. `forestsearch_Kfold()` in `R/forestsearch_cross_validation.R:287`, when
   `mr_in_replicates = TRUE` (formal at `:295`). Structurally identical to the
   bootstrap path — `mr_inference` is stripped at `:402–403` unless
   `mr_in_replicates`, and propagates into every fold when it is set. Same
   check, same message.
4. `fs_mr_inference()` in `R/fs_mr_inference.R:286` — **cannot** be guarded.
   Its signature is `(df, candidates, spec, selected_members, ...)` with no
   configuration visible. Add a roxygen note that it is an engine-level entry
   point assuming alignment has been established upstream.

### A3b. What is already implemented — do not duplicate

`.fs_mr_reselection_from_focus()` in `R/fs_mr_inference_methods.R:130` already
derives MR's re-selection rule per engine:

```r
hr_rule <- if (identical(engine, "consistency")) "maxcons" else "maxeff"
```

with `maxeffCons = "maxeff"` mapped explicitly at `:148`. So the
consistency-only asymmetry — FS ranks with a consistency term, DINA and GRF do
not — is **already correct inside the MR layer**. The gap is confined to the
identifier layer: `dina_subgroup()`'s whitelist and the GRF `frontier_rule`
`switch()`. Do not re-implement the mapping in MR.

Note also that MR's rule named `"maxeff"` is *gated* by `c_screen`,
`c_consistency`, and `p_star`, so it corresponds to `sg_focus = "maxeffCons"`,
not to `sg_focus = "maxeff"` (which is ungated at the identifier). The comment
at `:143–147` documents this. It implies a possible mismatch under
`sg_focus = "maxeff"` specifically — the identifier disables its effect floor
while MR's re-selection still applies `t_g`. Out of scope here; flag it if you
touch that path.

### A4. Family-status field

Record the fixed-family status on the object returned by `forestsearch()`, so a
sweep can tabulate it without parsing text. Suggested values: `"fixed"` (#1
with all three front ends `FALSE`), `"conditional-removable"` (#1 with a front
end on), `"conditional-inherent"` (#2–#6). Surface it in `print()`/`summary()`
alongside the estimate. It never raises a condition, so it has no volume cost.

### A5. Full bootstrap is unrestricted

`forestsearch_bootstrap_dofuture()` must remain valid for all six
configurations. It replays the entire pipeline per resample rather than
linearizing it, so §4.2's alignment condition does not bear on its validity —
the bootstrap is the reference standard in that iff, not a party to it. Note
in its roxygen that for #3, #5, #6 it correctly estimates the selection bias
that a native-statistic ranking induces on the reported $\hat\beta$, which is a
coherent quantity but outside the manuscript's framework.

### A6. NEWS.md

The `use_lasso`/`use_grf` default flip changes the candidate family on **every**
FS call, not only under MR. Flag as breaking. Note that reproducing León et al.
(2024), which used the prognostic lasso prefilter, now requires setting them
`TRUE` explicitly.

### A7. Verification

Run `R CMD check`. Do **not** build testthat scaffolding — testing convention
in this project is integration-style Quarto living documents.

---

## Phase B — application re-analyses

Revise, run, and compare against the manuscript payloads.

| analysis file | comparison target | built |
|---|---|---|
| `analysis_gbsg_survival_multimethod.qmd` | `gbsg_table_payload.rds` | 2026-06-21, fs 0.2.0 |
| `analysis_actg175_binary_multimethod.qmd` | `actg175_table_payload.rds` | 2026-06-21, fs 0.2.0 |

Each payload is a list whose `$table` is a 6-row data frame — `method` ∈
{FS, DINA, GRF} × `region` ∈ {H, Hᶜ} — with columns `n`, `pct`, `naive_est/lo/hi`,
`fb_est/lo/hi`, `mr_est/lo/hi`. `$labels` carries the selected subgroup
definitions, e.g. GBSG FS `{er <= 0} & {size <= 35}`, DINA
`{grade3 >= 1} & {pgr <= 10}`, GRF `{er <= 0}`.

Report as a side-by-side diff. Two tiers:

- **Labels and `n`** should match exactly. A changed subgroup definition means
  the identifier now selects differently and needs explaining, not tolerating.
- **Estimates and interval bounds** will differ slightly (resampling seeds,
  MR draws). Report absolute and relative differences; do not assert
  equivalence.

Confirm each analysis sets `dina_select_statistic = "effect"`,
`grf_select_statistic = "effect"`, `grf_selection = "frontier"` — under the
Phase A defaults these become redundant, but leaving them explicit documents
intent and guards against a future default change.

### B0. Both applications overwrite payloads on render

Each document ends by writing its own payload with `saveRDS()`:

- `analysis_gbsg_survival_multimethod.qmd:2630` → `gbsg_survival_multimethod_payload.rds`
- `analysis_actg175_binary_multimethod.qmd:2030` → `actg175_table_payload.rds`

`results_dir <- NULL` in both (GBSG `:189`, ACTG175 `:108`), so `.out_dir`
resolves to the directory the document renders from. **ACTG175's write filename
is identical to its reference bundle.** Rendering it in the same directory as
`actg175_table_payload.rds` destroys the comparison target.

Copy both reference payloads somewhere safe before rendering anything, or
render the documents outside `sim_analyses/`. Treat `sim_analyses/` as
read-only throughout.

**On `gbsg_table_payload.rds`.** The GBSG document writes
`gbsg_survival_multimethod_payload.rds`, not `gbsg_table_payload.rds` — so the
uploaded reference does not share its output filename. Its contents match what
that document's payload block constructs (same `table`/`labels`/`meta`/`extras`
structure, `built_at = 2026-06-21`, `forestsearch_version = "0.2.0"`), and the
document's own comment at `:2628` notes the payload filename has changed at
least once. Treat it as manuscript-era output of the same document, compare
against it, and say in the report that the filenames differ.

### B1. ACTG175 requires a configuration change — Larry has ruled

The two applications were **not** run under the same family regime.

- **GBSG** sets `use_lasso = FALSE`, `use_grf = FALSE`, `use_dina = FALSE` in
  its main fit (lines 700, 701, 710). Fixed family. Nothing to change.
- **ACTG175** sets `use_lasso = TRUE`, `use_grf = TRUE` (lines 541–542) and
  passes `mr_inference = run_mr` at `:564` with `run_mr <- TRUE` (`:100`).
  That is a condition-2 failure combined with MR — exactly what the Phase A
  validator now rejects. The published ACTG175 FS row is therefore a
  conditional-family result.

**Decision: follow the GBSG configuration.** Set `use_lasso <- FALSE` and
`use_grf <- FALSE` in the ACTG175 main fit. This aligns it with GBSG and with
what the supplement (S1, "A note on the lasso") calls the recommended
fully-enumerated default.

Consequences to expect and report, not to suppress:

1. **The FS subgroup may change.** Removing the prefilter changes the candidate
   family, so a different region may win. If the label differs from
   `actg175_table_payload.rds`, that is the **finding**, not a reproduction
   failure. Report the old and new labels side by side with the estimates for
   each.
2. **Only the FS row is affected.** The DINA and GRF calls (lines 1207 and
   1429) set `subgroup_method` explicitly, and the front-end flags are inert on
   those paths — the document's own comments say so. Those rows should
   reproduce normally.
3. **Watch for candidate-pool truncation.** Without the lasso prefilter the
   enumerated family is larger. `max_subgroups_search` defaults to `200`, and
   `subgroup_consistency_main.R:610–622` warns whenever the pool is actually
   truncated. If that warning fires, report it with both counts — it means the
   identifier ranked over fewer candidates than it enumerated.

Do not make the equivalent change to GBSG. It is already correct.

---

## Phase C — simulation re-runs, sim_id 1–20 only

**Do not run 1–500.** Run sim_id 1 through 20 and compare against the same
sim_ids extracted from the reference files.

### Comparison mechanics are verified and sound

- Every reference `.rds` is a list of `results` (data frame), `truth`, `meta`.
- `results` carries a `sim_id` column spanning 1–500, so
  `subset(results, sim_id <= 20)` is the comparison set.
- Per-replicate seeding is `seed_base + sim_id`, covering the DGM, the search,
  the bootstrap, and the noise covariates (documented in the sim `.qmd` at
  lines 108–109). `seed_base = 8316951` in all reference files. So sim_ids 1–20
  are exactly reproducible.
- The original FS run was batched with `res_1_20.rds` as its own unit, so
  1–20 is a natural boundary.

`results` columns of interest: `sim_id`, `detected`, `n_sel`, `label`,
`sg_def`, `covs`, the `nv_*` (naive), `t1_*` (FB), `t2_*` (MR) estimate blocks
for H and Hᶜ, `t2_harm_flag`, `ij_source`, and `sens`/`spec`/`ppv`/`npv`.

### C0. Where to render

**The constraint, not a layout.** No document may render in a directory that
holds a file it will overwrite. Two such files exist, both named in B0 and here:

| document | writes | collides with |
|---|---|---|
| `sim_fs_maxcons_*.qmd` | `results/fs_maxcons_fb_mr_m1_h10_knoise0_n500_res_1_20.rds` | the manuscript's first FS batch file |
| `analysis_actg175_binary_multimethod.qmd` | `actg175_table_payload.rds` | its own reference bundle |

The FS filename is byte-identical to what the manuscript run produced —
`meta$source_files` on the FS reference bundle lists
`fs_maxcons_fb_mr_m1_h10_knoise0_n500_res_1_20.rds`. That file is **not
currently on disk** (there is no `results/` directory), so this collision is
latent rather than live: it triggers on a second render, or if the batch files
are ever restored. The ACTG175 collision is live now. DINA and GRF do not
collide; their references use the older `_t1_t2_` stem.

Two ways to satisfy the constraint. Either is fine — pick one and say which in
the report:

- **Move the references.** Copy the reference `.rds` files out of the render
  directory first, then render in place. Fewest moving parts.
- **Move the renders.** Render in a separate directory. The three simulation
  documents also need `betaHhat_truth.R` reachable from there.

No directory needs creating by hand: the simulation documents create their own
`results/` via `dir.create(dirname(rds_path), recursive = TRUE)` (line 703), and
both application documents create `.out_dir`. Treat `sim_analyses/` as
read-only either way.

### C1. The three simulation documents

The original `sim_fs_maxeffCons_..._combine_1_500.qmd` set
`sg_focus <- "maxeffCons"`, which does **not** match any reference payload —
the FS reference records `sg_focus = "maxcons"`. Those are different selection
rules, not aliases: `maxcons` normalizes to canonical `hr` with FS key
`(-Pcons, -hr, K)`, consistency-led; `maxeffCons` has its own branch in
`sort_subgroups()`, `(-hr, K)` among qualifiers, effect-led. Larry has resolved
this in favour of `maxcons`.

**This is already fixed.** Three documents are supplied alongside this handoff,
one per engine, each with `subgroup_method` pinned:

| file | engine | reference |
|---|---|---|
| `sim_fs_maxcons_fb_mr_m1_h10_knoise0_n500.qmd` | `consistency` | `fs_maxcons_fb_mr_m1_h10_knoise0_n500_combined_1_500.rds` |
| `sim_dina_maxcons_fb_mr_m1_h10_knoise0_n500.qmd` | `dina` | `dina_t1_t2_m1_h10_knoise0_n500_combined_1_500.rds` |
| `sim_grf_maxcons_fb_mr_m1_h10_knoise0_n500.qmd` | `grf` | `grf_t1_t2_m1_h10_knoise0_n500_combined_1_500.rds` |

They are byte-identical apart from `subgroup_method`, the subtitle, and one
comment. Do not re-derive them from the original `sim_fs_maxeffCons_*` file.
Changes already applied and verified in all three:

- `sg_focus <- "maxcons"`, with the prose and comments made engine-aware.
- `sim_id_start <- 1L`, `n_sims <- 20L`, `run_mode <- "batch"`.
- A new **MR validity guard** chunk (`@sec-guard`) that hard-errors on any
  configuration failing condition 1 or 2, mirroring `.validate_mr_configuration()`.
  Verified to pass #1 (front ends off), #2, #4 and to reject #3, #5, #6 and #1
  with any front end on.
- Front-end flags re-annotated as the condition-2 requirement rather than a
  speed tradeoff; the stale "set `use_grf = TRUE` for FSlg" note replaced.
- All 12 R chunks parse under R 4.3.3.

Unchanged and matching the reference `meta`: `seed_base = 8316951`,
`nb_boots = 300`, `mr_draws = 5000`, `consistency_method = "resample"`,
`stop_threshold = NULL`.

### C2. Running all three

Render each document once. `subgroup_method` is pinned per file, so no editing
between runs. The `method_args` switch supplies each engine's own arguments and
the output stem renames itself from `method_tag`/`focus_tag`, giving
`fs_maxcons_*`, `dina_maxcons_*`, `grf_maxcons_*`.

`sg_focus = "maxcons"` is correct for all three. `.normalize_sg_focus()`
collapses `maxcons` and `eff` onto the same canonical focus, so it reproduces
the FS reference (recorded as `maxcons`) and the DINA and GRF references (run
under `eff`, which predates focus recording — their `meta` carries `gate_draws`
rather than `mr_draws` and no `sg_focus` field).

| run | `subgroup_method` | reference |
|---|---|---|
| 1 | `"consistency"` | `fs_maxcons_fb_mr_m1_h10_knoise0_n500_combined_1_500.rds` |
| 2 | `"dina"` | `dina_t1_t2_m1_h10_knoise0_n500_combined_1_500.rds` |
| 3 | `"grf"` | `grf_t1_t2_m1_h10_knoise0_n500_combined_1_500.rds` |

All three already set `dina_select_statistic = "effect"`,
`grf_selection = "frontier"`, `grf_select_statistic = "effect"`. Under the
Phase A defaults these become redundant; leave them explicit as documentation
of intent.


### C3. Add a pre-loop configuration check

The three documents deliberately contain **no** validation code — CC writes it,
reusing the Phase A validator rather than duplicating its logic. Add it to each.

It is needed because the run loop swallows errors. The `forestsearch()` call is
wrapped in `tryCatch(..., error = function(e) NULL)`, and on `NULL` the
replicate returns `.na_record()`, which sets `detected = 0L`. So a package-level
`stop()` raised inside the loop surfaces as a 0% detection rate rather than as a
failure — indistinguishable from an identifier that simply found nothing. That
is the same failure mode the original audit flagged for DINA under
`maxeffCons`.

Add the check **outside** the loop, after the setup chunk and before the run
chunk, so an invalid configuration halts the render immediately with a message
naming the offending argument. Keep it to a call into the Phase A validator
plus a one-line confirmation on success.

### C4. Reporting

Expect close but not identical results. Report per-replicate for sim_ids 1–20:
agreement on `detected`, on `sg_def`/`label`, on `n_sel`, and the differences
in the naive/FB/MR estimates. Summarize what changed and why, rather than
asserting a pass/fail threshold.

---

## Conventions

- Modify only what is explicitly requested; ask before touching other files.
- Read files immediately before editing; never infer unseen contents.
- Never convert `.qmd` ↔ `.Rmd`. Quarto is preferred.
- `Roxygen: list(markdown = TRUE)` is ON — write literal `% < > &`; never
  Rd-escape manually. `@section` titles are plain text only.
- Fixed abbreviations: FS, DINA, GRF, FB, MR. Never conflate β(Ĥ), the
  conditional estimand, with θ†(H), the marginal one.
- Report NPV alongside sensitivity, specificity, and PPV.
- `doFuture` workers need the installed package (`devtools::install()`, not
  `load_all()`).
- Larry's primary machine is Pop!_OS with ~127 cores.
