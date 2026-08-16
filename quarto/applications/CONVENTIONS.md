# Conventions for the analysis documents

Scope: the `.qmd` documents kept at the top level of `quarto/applications/gbsg/`
and `quarto/applications/actg175/`. Files under `_archive/` and `_broken/` are
superseded and follow none of this.

Everything below is transcribed from the documents as they currently stand. When
a convention and a document disagree, the document is the fact and this file is
the bug.

This directory is not shipped: `.Rbuildignore` contains `^quarto$`, which prunes
the whole `quarto/` tree, nested paths included.

---

## 1. Payload paths

Every document that exports results carries this block in its setup chunk, with
only the illustrative `dirout` string adapted to the document:

```r
# ---------------------------------------------------------------------------
# Output location for this document's results:  <base>/<dirout>/
#   results_dir = NULL -> base is <qmd dir>/_payloads
#   dirout      = NULL -> the .qmd stem
# Setting results_dir redirects the base only; each run still gets its own
# subdirectory named by dirout.
#
# Examples:
#   results_dir <- NULL;  dirout <- NULL
#     -> <qmd dir>/_payloads/<stem>/<stem>_payload.rds
#   results_dir <- NULL;  dirout <- "gbsg_survival_sgfocus"
#     -> <qmd dir>/_payloads/gbsg_survival_sgfocus/gbsg_survival_sgfocus_payload.rds
#   results_dir <- "~/GitHub/fs-post-selection-manuscript_jasa/_payloads"; dirout <- NULL
#     -> ~/GitHub/fs-post-selection-manuscript_jasa/_payloads/<stem>/<stem>_payload.rds
# ---------------------------------------------------------------------------
results_dir <- NULL
dirout      <- NULL
```

The export chunk at the end of the document resolves the path:

```r
.qmd_dir  <- tryCatch(dirname(knitr::current_input(dir = TRUE)),
                      error = function(e) getwd())
.qmd_stem <- tryCatch(tools::file_path_sans_ext(basename(knitr::current_input())),
                      error = function(e) "analysis")
.dirout   <- if (is.null(dirout)) .qmd_stem else dirout
.base     <- if (is.null(results_dir)) file.path(.qmd_dir, "_payloads") else results_dir
.out_dir  <- file.path(.base, .dirout)
dir.create(.out_dir, recursive = TRUE, showWarnings = FALSE)
.payload_file <- file.path(.out_dir, paste0(.dirout, "_payload.rds"))
```

Two properties follow, and both are load-bearing.

**`dirout` defaults to the `.qmd` stem, so payload paths cannot collide.** Two
files cannot share a name within one directory, so the stem is unique by
construction; the directory and the file inside it both take that name. Setting
`results_dir` redirects only the base, so each document still lands in its own
`<dirout>/` subdirectory underneath it.

**Payload filenames are never hardcoded literals.** The filename is always
`paste0(.dirout, "_payload.rds")` — derived, never typed. Building it from
`.qmd_stem` instead of `.dirout` is also wrong: the two coincide only while
`dirout` is `NULL`, and diverge the moment someone sets it, leaving the file
named differently from the directory holding it.

This is not hypothetical hygiene. Two documents once both wrote a literal
`actg175_table_payload.rds` and silently overwrote each other.

`analysis_actg175_continuous_compare_all.qmd` carries the same config block and
resolves `output_dir <- file.path(.base, .dirout)` identically, but writes
`selected_subgroups_continuous.rds` and `comparison_continuous.rds` rather than
a payload. Those two names are read by name downstream and must not be
renamed to fit the pattern.

---

## 2. Payload shape

Every payload writes the same seven top-level elements, in this order:

```
table  labels  meta  extras  est_scale  built_at  forestsearch_version
```

No document omits one and none adds an eighth. `est_scale` is `"hr"` in the
gbsg documents and `"or"` in the actg175 ones.

Divergence **below** the top level is deliberate. `extras$gh` exists only where
a Guo & He section runs; `extras$mr_confirm` only in `maxeff_mrconfirm`;
`extras$frozen_gate` and `extras$fb_mr_gap` only in `frozen_family`;
`extras$pareto` only where a banded focus can run the frontier. A reader can
still consume all payloads uniformly, because R returns `NULL` for a missing
list name — `payload$extras$gh` is `NULL` on a document that never had one, the
same value it takes on a document whose Guo & He section was switched off.

Every element is guarded so a section that did not run stores `NULL` rather
than erroring. The guards are `exists()`/`is.null()` based, via two helpers
defined at the top of each export chunk:

```r
.get <- function(nm) if (exists(nm, inherits = TRUE)) get(nm, inherits = TRUE) else NULL
.try <- function(expr) tryCatch(expr, error = function(e) NULL)
```

### The `n_events` guard

```r
n_total     = .try(nrow(.get("df.analysis"))),
# sum(NULL) is 0, so guard on the data being present: a render that never
# reached the data chunk must store NULL, not a spurious zero.
n_events    = .try(if (is.null(.get("df.analysis"))) NULL else
                   sum(.get("df.analysis")[[.get("event.name")]])),
```

`.try()` alone is not enough here, because nothing throws. When the data object
is absent, `.get()` returns `NULL`, `NULL[[NULL]]` returns `NULL`, and
`sum(NULL)` returns `0` — a valid number, so `tryCatch` never fires and the
payload records an event count of zero for a render that counted nothing. The
`is.null()` check is what turns that into `NULL`.

`n_total` needs no such guard: `nrow(NULL)` is already `NULL`, because
`dim(NULL)` is `NULL`. The asymmetry is the whole reason one line has the check
and the line above it does not — it is not an oversight.

The actg175 binary documents use their own object names (`actg_df`,
`outcome.name`) in the same construction. Object names are per-document; the
guard is not.

The five multimethod payloads express the equivalent guard differently, as
`n_total = if (!is.null(.fs)) nrow(.fs$df.est) else NA_integer_`, and carry no
`n_events`.

---

## 3. LOO cache keys

Four gbsg documents cache their leave-one-out result. The key is built
identically in all four:

```r
.cache_doc   <- tryCatch(tools::file_path_sans_ext(basename(knitr::current_input())),
                         error = function(e) "analysis")
.cache_focus <- if (is.null(fs$sg_focus)) "unknownfocus" else as.character(fs$sg_focus)[1]
.cache_rule  <- {
  .r <- fs$args_call_all$selection_rule
  if (is.null(.r) || !nzchar(as.character(.r)[1])) "neighborhood" else as.character(.r)[1]
}
loo_cache <- Sys.getenv(
  "LOO_CACHE",
  unset = file.path(gh_dir, sprintf("cv_out_%s_%s_%s.rds",
                                    .cache_doc, .cache_focus, .cache_rule)))
```

Resolving to:

| document | cache file |
|---|---|
| `analysis_gbsg_survival_effMaxSG` | `cv_out_analysis_gbsg_survival_effMaxSG_hrMaxSG_neighborhood.rds` |
| `analysis_gbsg_survival_sgfocus` | `cv_out_analysis_gbsg_survival_sgfocus_hrMaxSG_both.rds` |
| `analysis_gbsg_survival_frozen_family` | `cv_out_analysis_gbsg_survival_frozen_family_maxeff_neighborhood.rds` |
| `analysis_gbsg_survival_maxeff_mrconfirm` | `cv_out_analysis_gbsg_survival_maxeff_mrconfirm_maxeff_neighborhood.rds` |

### Which axis does the work

**The document stem is what guarantees uniqueness.** It is the only axis that
cannot collide, for the same reason `dirout` cannot: filenames are unique within
a directory. Every other axis can tie. `frozen_family` and `maxeff_mrconfirm`
both run `sg_focus = "maxeff"` under the default `selection_rule`, and are told
apart *only* by their candidate families — median-cut versus quantile-cut, a
difference that lives in `conf_force` and `conf.cont_jcuts` and that no scalar
setting in a filename can express. A focus-plus-rule key would map both to
`cv_out_maxeff_neighborhood.rds`.

**Focus and selection rule are diagnostic, not load-bearing.** They make the
filename say what it holds, so a directory listing is readable and a stale cache
is recognisable on sight. They are read from the *fitted object*
(`fs$sg_focus`, `fs$args_call_all$selection_rule`), not from the `fs_*` config
variables, so the name reports what the engine actually ran with rather than
what was typed above it.

### The bug this replaced

The key used to be `sprintf("cv_out_%s.rds", focus$canonical)` — focus only.
`analysis_gbsg_survival_effMaxSG.qmd` and `analysis_gbsg_survival_sgfocus.qmd`
are identical except for `fs_selection_rule` (`"neighborhood"` versus
`"both"`), and both canonicalize to `hrMaxSG`, so both resolved to
`cv_out_hrMaxSG.rds`. The read path is `if (file.exists(loo_cache)) cv_out <-
readRDS(loo_cache)`, so on a fresh run whichever document rendered second loaded
the first's leave-one-out result and reported it as its own — no recompute, no
warning. The comment above the old key claimed the focus key prevented silent
reuse; that was true for focus and silent about every other axis.

Fixed in `48c419eb`, which also deleted the three stale `cv_out_*.rds` then in
`quarto/GuoHe/`.

### Reading the filenames

A cache file named `cv_out_analysis_*` means the `knitr::current_input()`
fallback fired — the chunk ran outside a knitr render, `.cache_doc` fell back to
the literal `"analysis"`, and the document axis collapsed. Such a file carries
no document identity and is shared by anything else that hits the same fallback.
Treat it as untrusted and delete it.

`LOO_CACHE` overrides the whole construction, so an explicitly set env var
bypasses every guarantee above. That is intended, for pointing a render at a
known cache, but it is the caller's problem to keep unique.

---

## 4. `max_subgroups_search`

**Zero configuration uses remain in any kept document.** No variable assignment,
no focus-profile parameter, no argument in any `forestsearch()` call including
the DINA and GRF screening calls, no row in any config-record table. The
package default is `Inf`, and the documents now take it.

The string still appears, and every remaining occurrence is deliberate:

- **Provenance assertions.** `fs_param_provenance()` verifies the engine-forced
  values under `sg_focus = "maxeff"` against a literal expectation:
  ```r
  forced <- list(pconsistency.threshold = 0, stop_threshold = NULL,
                 use_twostage = FALSE, max_subgroups_search = Inf, minp = 0)
  ```
  Removing the entry would drop that row from the check. The `.ov` name vectors
  and `.fs_selection_knobs` list serve the same purpose.
- **Prose.** The Guo & He skip note explains that `maxeff` guarantees an
  untruncated family (`` `max_subgroups_search -> Inf` ``).
- **The `fixed_family` comment block**, which records that the published
  supplement Table S11 came from an earlier run capped at
  `max_subgroups_search = 50L`, that the cap is gone, and that results may
  therefore differ from the published table intentionally.

These describe the *package's* behaviour or the document's history. They are not
settings, and they should not be swept.

### It does not bear on the fixed-family condition

The fixed-family condition requires that the candidate family be **finite and
held fixed under resampling** — not that it be maximal. `max_subgroups_search`
controls how many candidates are evaluated; it does not control whether the
family regenerates. A truncated family is still a fixed family if its cuts do
not move, and an untruncated one still regenerates if they do.

What actually governs regeneration is the cut construction, which
`frozen_family` exists to demonstrate: under `conf.cont_jcuts` the *instruction*
"cut at J quantiles" is replayed on each bootstrap resample, so cut locations
move and the family regenerates; under a `conf_force` expression with a fixed
numeric threshold the cut is a literal string, replayed verbatim, and cannot
move. So removing the cap changed the family's *size*, and may change which
subgroup is selected, but it did not change whether any analysis satisfies the
fixed-family condition.

---

## 5. Open item

`analysis_actg175_crossanalysis_summary.qmd` reads two files that **no kept
document writes**:

- `selected_subgroups_binary.rds` (loaded at its `load_if_exists()` call)
- `comparison_binary.rds` (loaded for the per-subject agreement section)

Its own header table attributes both to `analysis_actg175_binary_compare_all.qmd`,
which now lives only under `actg175/_archive/`. The continuous halves
(`selected_subgroups_continuous.rds`, `comparison_continuous.rds`) are still
produced, by `analysis_actg175_continuous_compare_all.qmd`.

Both reads are wrapped in `file.exists()` guards, so the document degrades
quietly rather than erroring: the binary rows are absent and the per-subject
agreement section prints a skip message. That means the cross-analysis renders
clean while silently reporting only half of what it is named for. Unresolved —
either the writer is restored from the archive, or the binary half is removed
from the document.
