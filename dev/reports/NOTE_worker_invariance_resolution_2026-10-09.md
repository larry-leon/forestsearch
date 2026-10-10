# NOTE — Worker-count invariance: resolution (2026-10-09)

- Follow-up to `REPORT_worker_invariance_diagnosis_2026-10-09.md` (committed at `982bfe80`). That report classified the failure as **H0, a harness artifact**: `multisession` workers loaded a stale installed build (0.2.2), while the main session ran HEAD (0.4.0) through `load_all()`.

## Action taken (proposal 1 of the report, approved by Larry)

- I installed HEAD into the user library with `Rscript -e 'devtools::install(dependencies = FALSE, upgrade = FALSE)'`.
  - The command as specified, with `upgrade = "never"`, aborted before installing anything. This devtools version only accepts `TRUE`, `FALSE` or `NA` for `upgrade`. `upgrade = FALSE` has the same meaning: never upgrade dependencies.
- Installed state afterwards:
  - `packageVersion("forestsearch")`: **0.4.0** (previously 0.2.2).
  - Library path: `/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/forestsearch`.
  - Packaged: `2026-10-10 00:15:22 UTC; larryleon`, built from HEAD `982bfe80`.

## Verification

- I re-ran `devtools::test(filter = "search-reproducibility")`, the same route that failed in Stage 3.
- "search + consistency results are invariant to worker count": **PASS**. All 6 expectations passed, with 0 failed, 0 skipped and no errors.
- This confirms H0. With the main session and the workers on the same build, the invariance guarantee holds at HEAD.
- Only this test file was run, not the full suite.

## Open items (from the report; not applied)

- Proposal 2: add a test guard that skips the test with an explicit message when the installed version differs from the source `DESCRIPTION` version.
- Proposal 3: narrow the `suppressWarnings()` at `R/subgroup_consistency_main.R:986` so `future`'s package-version-mismatch warning is no longer hidden.
- Going forward: reinstall HEAD before any `devtools::test()` run that includes `multisession` tests. Otherwise the user library can fall behind again.
