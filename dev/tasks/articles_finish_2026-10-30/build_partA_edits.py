#!/usr/bin/env python3
"""Edit set for NEWS.md and survival-post-selection.qmd (articles-finish Part A, 2026-10-30).
Same JSON shape as alignment_edits_2026-10-30.json, applied by apply_alignment_edits.py."""
import json, hashlib, os, sys
SRC = sys.argv[1]; OUT = sys.argv[2]
E = []
def edit(f, i, sec, old, new, basis): E.append(dict(file=f, id=i, section=sec, old=old, new=new, basis=basis))

edit("NEWS.md", "N1", "0.4.1 Documentation — the aligned articles and the new article",
"* The methodology and simulation articles on the package website carry a note placing them as the 2024 workflow, with the deprecated arguments removed from their listings.\n",
"* The methodology and simulation articles on the package website keep their text as the record of the 2024 workflow, with the deprecated arguments removed from their listings and a dated note at each point where the current method differs. The articles on treatment-effect definitions, biomarker effects and the causal interpretation of hazard ratios now state the estimand the package reports for an identified subgroup, as in León and Anderson (2026a, 2026b). All five articles are now Quarto documents.\n"
"* New article on the package website, \"Post-selection bounds on a nested continuum: the Guo and He (2021) design\": Section S4 of León and Anderson (2026b), with the complete code, a small default run, and the paper's results beside it.\n"
"* The README describes the current method and shows the analysis of `vignette(\"forestsearch\")`.\n",
"NEWS 0.4.1 Documentation section; no URL to an article not yet on the site (URL checks).")

edit("vignettes/articles/survival-post-selection.qmd", "P1", "The settings — DINA floor",
"- **DINA** is re-ranked on the Cox coefficient against a floor of 0.90, with a\n  two-fold cross-fit surface; **GRF**",
"- **DINA** is re-ranked on the Cox coefficient against a floor of 0.90, with a\n  two-fold cross-fit surface. The floor of 0.90 is the one in the document that\n  produced Table S2, which this call reproduces; Section 6 of the paper gives a\n  hazard ratio of 1.0. **GRF**",
"Open finding 9 of the 10-28 handoff; 2026a §6 'DR-score and mean-effect floors of 0 and log 1.0'; postsel stage-0 reproduced Table S2 at 2 dp with the document's calls.")

files = sorted({e["file"] for e in E})
txt = {f: open(os.path.join(SRC, os.path.basename(f)), encoding="utf-8").read() for f in files}
md5 = lambda s: hashlib.md5(s.encode("utf-8")).hexdigest()
out = dict(txt)
for e in E:
    assert out[e["file"]].count(e["old"]) == 1 and out[e["file"]].count(e["new"]) == 0, e["id"]
    out[e["file"]] = out[e["file"]].replace(e["old"], e["new"], 1)
os.makedirs(os.path.join(OUT, "partA"), exist_ok=True)
for f in files:
    open(os.path.join(OUT, "partA", os.path.basename(f)), "w", encoding="utf-8").write(out[f])
spec = dict(issued="2026-10-30", base_commit="0bdc2036",
            files={f: dict(md5_before=md5(txt[f]), md5_after=md5(out[f]),
                           n_edits=sum(e["file"] == f for e in E)) for f in files}, edits=E)
json.dump(spec, open(os.path.join(OUT, "articles_finish_edits_2026-10-30.json"), "w", encoding="utf-8"),
          ensure_ascii=False, indent=1)
for f in files: print(f, spec["files"][f])
