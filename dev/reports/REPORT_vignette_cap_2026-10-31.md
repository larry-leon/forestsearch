# REPORT — the primary vignette: the declaration at the recommended cap, and the release re-check (2026-10-31)

Task: `dev/tasks/cc_task_vignette_cap_2026-10-31.md`. Branch `feature/vignette-cap`, from `master` at `748cec46`.
Record: `dev/reports/vignette_cap_2026-10-31/` (`logs/`, `results/`).
**Outcome: every gate passed. No §3.3 deviation.**

## §0 Provenance

- HEAD `748cec46ca357e229f5e9bb15b6fd7a38f8ea103` on `master`; no tracked changes; no fetch.
  **UNTRACKED_0: empty.**
- `vignettes/forestsearch.qmd` md5 at start: `037aad7ccfe9811bc59181356cc180a0` (matches).
- The zip's md5 was `0e4a71661e7b63f6209b397c1392343c` (matches). The zip has no top-level directory, so I
  unzipped it with `-d dev/tasks/vignette_cap_2026-10-31`. The filenames arrived unaltered.
  `MANIFEST.md5`: both files are `OK`.
- Commit 1: `35cb202d task: vignette — the declaration at the recommended cap — task document and inputs`.
- R 4.6.1 (2026-06-24). forestsearch 0.4.1, quarto (R) 1.5.1, knitr 1.52, rmarkdown 2.32, xml2 1.6.0.
  The Quarto CLI is `/usr/local/bin/quarto` 1.10.18 (on `PATH`). `pandoc` on `PATH` is `/usr/bin/pandoc` 3.1.3.
  `RSTUDIO_PANDOC` (pandoc 3.10) was exported for the build and the check. (`results/provenance.txt`)

## Gates

- §0.1 branch/HEAD/clean tree: PASS. §0.2 vignette md5: PASS. §0.4 zip md5 and manifest: PASS.
  §0.6 forestsearch 0.4.1 and xml2 installed: PASS.
- §1.1 placed md5 `f019e68cfa8d2538e889d834fc628d59`: PASS.
  §1.2 numstat `63	5	vignettes/forestsearch.qmd`: PASS (`results/gate1_numstat.txt`).
  §1.3 `parse OK`: PASS (`logs/parse.log`).
- §2.2 `R CMD build src`: exit 0, and the tarball exists: PASS.
  §2.3 the tarball has no path under `vignettes/articles/`, `dev/` or `docs/` (625 entries): PASS
  (`results/tarball_paths.txt`). §2.4 the tarball's `forestsearch.qmd` md5 is `f019e68c…`: PASS.
- §3.2 `vignette_cap_gate.R` on the built `inst/doc/forestsearch.html`: **`RESULT: OK`**, with all 20 checks passing
  (`results/gate_vignette.txt`).
  - Table 1 has 10 rows, identical. The declaration table has 5 rows, identical.
  - The new heading appears once. The prose counts are as specified. There is no `NaN`, `?@` or `#> Error`.
- §4.1 `R CMD check --as-cran --run-donttest`: **`Status: OK`**: PASS (`logs/00check.log`, `logs/check.log`).
- §5.2 post-conditions: all three hold (`results/postconditions.txt`).
  - Outside `dev/`, the diff from `748cec46` is only `vignettes/forestsearch.qmd`.
  - The committed vignette and the tarball's vignette both have md5 `f019e68cfa8d2538e889d834fc628d59`.
  - `git ls-files 'vignettes/articles/*.html'` prints nothing.

## Key values

- The cap table, from the built HTML:
  - cutoff of the screen at p\* = 0.99: **2.432**
  - T of the declared region: **2.438**
  - margin over the capped cutoff: **0.006**
  - uncapped calibrated cutoff (c0 = HR 0.75): **2.726**
  - diagnostic (8) at the capped cutoff: **0.186**
- Recorded diagnostic line: `recorded (not gated): cap diagnostic 0.186; the paper reports 0.186 (MATCH)`.

## Compute

- `R CMD build`: 97 s (cap 15 min; part B took 93 s). (`results/build_wall.txt`, `logs/build.log`)
- `R CMD check --as-cran --run-donttest`: 284 s (cap 45 min; part B took 294 s). (`results/gate4_check.txt`)
  - Vignette re-build: [73s/72s]. `cran-comments.md` says "about 75 s".
  - Examples: [41s/42s]. `cran-comments.md` says "about 40 s".
  - Both are within 50% of the stated figures, so there is no OPEN ITEM.
- Tarball `forestsearch_0.4.1.tar.gz` md5: `f8352f31acaf3b76d8ab939c84109c15` (`results/tarball_md5.txt`).

## Commits

1. `35cb202d` task: vignette — the declaration at the recommended cap — task document and inputs
2. `2b424163` docs: vignette — the declaration at the recommended cap (2026b §4.4, §6.2); field-s named as the studentized complement bound
3. dev: report — vignette cap reading (2026-10-31) (this report and the record)

## §3.3 deviation

- None. The `declaration-cap` chunk ran as written, and both of its `stopifnot()` checks passed.

## OPEN ITEMS

- None.

## Verdict

```
VERDICT — vignette cap (2026-10-31): PASS
  gates §0–§5: all PASS; no deviation
  cap table: 2.432 / 2.438 / 0.006 / 2.726 / 0.186
  cap diagnostic vs paper 0.186: MATCH
  R CMD check --as-cran --run-donttest: Status: OK (284 s); build 97 s
  tarball md5: f8352f31acaf3b76d8ab939c84109c15
  branch feature/vignette-cap, 3 commits on 748cec46; not pushed, not merged
  OPEN ITEMS: none
```
