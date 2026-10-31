#!/usr/bin/env python3
"""Edit set for articles-finish part B (2026-10-30): ?forestsearch examples, the vignette's
next steps, NEWS, README, cran-comments and _pkgdown.yml.  Same JSON shape as
alignment_edits_2026-10-30.json, applied by dev/tasks/apply_alignment_edits.py with
repo-relative file keys.

    python3 build_partB_edits.py <dir holding the six source files, repo layout> <out dir>
"""
import json, hashlib, os, sys

SRC, OUT = sys.argv[1], sys.argv[2]
E = []


def edit(f, i, sec, old, new, basis):
    E.append(dict(file=f, id=i, section=sec, old=old, new=new, basis=basis))


# --- R/forestsearch_main.R: roxygen @examples only (inside \dontrun) -----------------------
edit("R/forestsearch_main.R", "R1", "@examples, Example 1 header",
     "#' ## León and Anderson (2026a, Sec. 5) and of the survival-post-selection\n",
     "#' ## León and Anderson (2026a, Sec. 6) and of the survival-post-selection\n",
     "2026a Section 6 is Applications (Section 5 is Scope); Larry approved 2026-10-30.")

edit("R/forestsearch_main.R", "R2", "@examples, Example 1 — mr_estimates_table() needs a bootstrap result",
     "#' print(fs_surv)\n"
     "#' mr_estimates_table(fs_surv)\n",
     "#' print(fs_surv)\n"
     "#'\n"
     "#' ## mr_estimates_table() lays the multiplier-resampling estimates out in the\n"
     "#' ## format of the bootstrap tables, so it takes a bootstrap result as its\n"
     "#' ## template (a few replicates are enough for the layout; the\n"
     "#' ## survival-post-selection article compares the two at full size).\n"
     "#' fs_boot <- forestsearch_bootstrap_dofuture(\n"
     "#'   fs.est        = fs_surv,\n"
     "#'   nb_boots      = 20,\n"
     "#'   parallel_args = list(plan = \"multisession\", workers = 4)\n"
     "#' )\n"
     "#' mr_estimates_table(fs_surv, boot_results = fs_boot, est.scale = \"hr\")\n",
     "mr_estimates_table(fs, boot_results, est.scale, digits) has no default for boot_results "
     "(R/mr_estimates_table.R:43) and needs its FSsg_tab and SG_CIs; Larry approved 2026-10-30.")

edit("R/forestsearch_main.R", "R3", "@examples, Example 2 header — which paper's analysis it runs",
     "#' ## the binary application of León and Anderson (2026b). Adverse-scale\n"
     "#' ## outcome (1 = no CD4 improvement at week 20), OR thresholds at the null,\n"
     "#' ## same selection rule and multiplier-resampling gate as Example 1.\n",
     "#' ## the forest-search analysis of León and Anderson (2026a, Table S3); the\n"
     "#' ## binary analysis of León and Anderson (2026b) selects the same region\n"
     "#' ## with effect_neighborhood = 0.20. Adverse-scale outcome (1 = no CD4\n"
     "#' ## improvement at week 20), OR thresholds at the null, same selection rule\n"
     "#' ## and multiplier-resampling gate as Example 1.\n",
     "Example 2's call is 2026a's ACTG 175 forest-search call (effect_neighborhood 0.10, "
     "fs.splits 500: stage-0 settings.md §1, §3); 2026b's binary call uses 0.20. Comment only.")

# --- vignettes/forestsearch.qmd: Where to go next -----------------------------------------
edit("vignettes/forestsearch.qmd", "V1", "Where to go next — the Guo and He and ACTG 175 articles",
     "  agreement, and what each construction does and does not say.\n"
     "- The package website carries the original workflow",
     "  agreement, and what each construction does and does not say.\n"
     "- The [Guo and He article](https://larry-leon.github.io/forestsearch/articles/guohe-nested-continuum.html)\n"
     "  sets the bounds beside the de-biased estimate and one-sided bound of Guo and\n"
     "  He (2021) on a nested continuum of subgroups, Section S4 of León and\n"
     "  Anderson (2026b).\n"
     "- The ACTG 175 analyses are an article on the\n"
     "  [package website](https://larry-leon.github.io/forestsearch/articles/):\n"
     "  forest search on a binary and a continuous endpoint, with the bounds and\n"
     "  the declaration calibration (León and Anderson 2026b, Section S3), and the\n"
     "  three identifiers on the binary endpoint with the full bootstrap beside\n"
     "  multiplier resampling (León and Anderson 2026a, Table S3).\n"
     "- The package website carries the original workflow",
     "Larry approved 2026-10-30. The ACTG 175 article's own URL is not live until the merge, "
     "so the bullet links the articles index (URL checks).")

edit("vignettes/forestsearch.qmd", "V2", "Where to go next — ?forestsearch",
     "- `?forestsearch` documents every argument, with examples for this analysis and\n"
     "  for the binary and continuous ACTG 175 analyses of León and Anderson\n"
     "  (2026b), and a section on the options not tested against the manuscript\n"
     "  analyses.\n",
     "- `?forestsearch` documents every argument, with examples on this trial and on\n"
     "  the binary and continuous ACTG 175 endpoints, and a section on the options\n"
     "  not tested against the manuscript analyses.\n",
     "Example 1 is 2026a's GBSG call (epsilon 0.10) and Example 2 2026a's ACTG 175 call; "
     "neither is 2026b's (edit R3).")

# --- NEWS.md: 0.4.1 Documentation -----------------------------------------------------------
edit("NEWS.md", "N1", "0.4.1 Documentation — the ACTG 175 article and the examples",
     "* The README describes the current method and shows the analysis of `vignette(\"forestsearch\")`.\n",
     "* New article on the package website, \"Binary and continuous outcomes: the ACTG 175 analyses\": "
     "Section S3 of León and Anderson (2026b) — forest search on the binary and continuous endpoints, "
     "with the selection-adjusted bounds and the declaration calibration — and Table S3 of León and "
     "Anderson (2026a), the three identifiers on the binary endpoint with the full bootstrap, at a small "
     "number of replicates, beside multiplier resampling.\n"
     "* The README describes the current method and shows the analysis of `vignette(\"forestsearch\")`.\n"
     "* `?forestsearch`: Example 1 passes a bootstrap result to `mr_estimates_table()`, which needs "
     "one as its template (the example called it without one); Example 2 names the analysis it runs, "
     "León and Anderson (2026a, Table S3); a section reference is corrected.\n",
     "No URL to an article not yet on the site (URL checks).")

# --- README.md: Documentation ---------------------------------------------------------------
edit("README.md", "M1", "Documentation — the articles and the examples",
     "  articles on post-selection inference for the same trial (the full\n"
     "  bootstrap beside multiplier resampling), on the comparison with Guo and He\n"
     "  (2021), on treatment-effect definitions and their causal reading, on the\n"
     "  simulation framework, and on the original workflow of León et al. (2024).\n"
     "- `?forestsearch` documents every argument, with examples for this analysis\n"
     "  and for the binary and continuous ACTG 175 analyses of León and Anderson\n"
     "  (2026b).\n",
     "  articles on post-selection inference for the same trial (the full\n"
     "  bootstrap beside multiplier resampling), on the binary and continuous\n"
     "  analyses of the ACTG 175 trial, on the comparison with Guo and He (2021),\n"
     "  on treatment-effect definitions and their causal reading, on the\n"
     "  simulation framework, and on the original workflow of León et al. (2024).\n"
     "- `?forestsearch` documents every argument, with examples on this trial and\n"
     "  on the binary and continuous ACTG 175 endpoints.\n",
     "As V2.")

# --- cran-comments.md -----------------------------------------------------------------------
edit("cran-comments.md", "C1", "Submission — Documentation bullet",
     "* Documentation: a new getting-started vignette reproducing a published\n"
     "  analysis; the longer post-selection analysis is now an article on the\n"
     "  package website, so the package ships one vignette instead of two.\n",
     "* Documentation: a new getting-started vignette reproducing a published\n"
     "  analysis; the longer analyses are articles on the package website, so the\n"
     "  package ships one vignette instead of two; the `?forestsearch` examples\n"
     "  are corrected (one called `mr_estimates_table()` without the bootstrap\n"
     "  result it requires).\n",
     "The NEWS lines of N1.")

# --- _pkgdown.yml ---------------------------------------------------------------------------
edit("_pkgdown.yml", "P1", "navbar — Articles menu",
     "        - text: Comparison with Guo & He (2021)\n"
     "          href: articles/guohe-nested-continuum.html\n",
     "        - text: Comparison with Guo & He (2021)\n"
     "          href: articles/guohe-nested-continuum.html\n"
     "        - text: Binary and Continuous Outcomes (ACTG 175)\n"
     "          href: articles/actg175-glm-outcomes.html\n",
     "The new article.")

edit("_pkgdown.yml", "P2", "articles index — Post-Selection Inference",
     "      - articles/survival-post-selection\n"
     "      - articles/guohe-nested-continuum\n",
     "      - articles/survival-post-selection\n"
     "      - articles/guohe-nested-continuum\n"
     "      - articles/actg175-glm-outcomes\n",
     "pkgdown requires every article in the index (check_pkgdown).")

files = list(dict.fromkeys(e["file"] for e in E))
txt = {f: open(os.path.join(SRC, f), encoding="utf-8", newline="").read() for f in files}
md5 = lambda s: hashlib.md5(s.encode("utf-8")).hexdigest()
out = dict(txt)
for e in E:
    assert out[e["file"]].count(e["old"]) == 1, ("old", e["id"], out[e["file"]].count(e["old"]))
    assert out[e["file"]].count(e["new"]) == 0, ("new", e["id"])
    out[e["file"]] = out[e["file"]].replace(e["old"], e["new"], 1)
for e in E:
    assert out[e["file"]].count(e["new"]) == 1, ("post", e["id"])
for f in files:
    p = os.path.join(OUT, "edited", f)
    os.makedirs(os.path.dirname(p), exist_ok=True)
    open(p, "w", encoding="utf-8", newline="").write(out[f])
spec = dict(issued="2026-10-30", base_commit="da17be3c",
            files={f: dict(md5_before=md5(txt[f]), md5_after=md5(out[f]),
                           n_edits=sum(e["file"] == f for e in E)) for f in files},
            edits=E)
json.dump(spec, open(os.path.join(OUT, "articles_finish_partB_edits_2026-10-30.json"), "w",
                     encoding="utf-8"), ensure_ascii=False, indent=1)
for f in files:
    print(f, spec["files"][f])
