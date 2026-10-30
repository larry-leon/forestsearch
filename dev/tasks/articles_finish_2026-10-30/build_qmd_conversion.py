#!/usr/bin/env python3
"""Convert the five aligned R Markdown articles to Quarto (2026-10-31 task, Part A).

Input: the five .Rmd exactly as committed on master after the alignment merge
(md5 checked).  Output:
  qmd/<name>.qmd           -- the conversion (diagrams as native Quarto {dot} cells)
  qmd_grviz/<name>.qmd     -- fallback for the two articles with diagrams: identical
                              except that the DiagrammeR::grViz() chunks are kept
The body text is unchanged except for three repairs that the conversion needs:
  * methodology: the power-approximation display equation gets the label
    {#eq-power}, which the existing cross-reference "@eq-power" points to
    (R Markdown rendered it as "(eq-power?)");
  * treatment_effect_definitions: the bookdown reference
    "Section \\@ref(theoretical-relationships)" becomes a link;
  * the YAML header: R Markdown html_document options -> Quarto html options,
    in the form the package's other Quarto articles use.
"""
import hashlib, os, sys

SRC = sys.argv[1] if len(sys.argv) > 1 else "/home/claude/outputs/aligned"
OUT = sys.argv[2] if len(sys.argv) > 2 else "/home/claude/outputs"
EXPECT = {  # post-alignment md5 on master (a88f8ea4)
    "methodology.Rmd": "5d79611e4a538e2b54aadaa151fd27d1",
    "paper_simulations.Rmd": "6154d7ea3b2fd8df2b831b1ae5524443",
    "treatment_effect_definitions.Rmd": "2c46b87fccf0114413d602c4b74e618e",
    "biomarker_effects.Rmd": "42a79e8cbc708576fd5950929a95c4b3",
    "causal_effects_brief_review.Rmd": "94caca1f6166919499fce0f1d6b79aaa",
}

def md5(s): return hashlib.md5(s.encode("utf-8")).hexdigest()

def once(text, old, new, tag):
    n = text.count(old)
    if n != 1:
        sys.exit(f"{tag}: anchor occurs {n} times")
    return text.replace(old, new, 1)

def split_yaml(text):
    assert text.startswith("---\n")
    j = text.index("\n---\n", 4)
    return text[: j + 5], text[j + 5 :]

def vig(entry):
    return ("vignette: >\n"
            f"  %\\VignetteIndexEntry{{{entry}}}\n"
            "  %\\VignetteEngine{quarto::html}\n"
            "  %\\VignetteEncoding{UTF-8}\n")

def fmt(code_fold, fw, fh):
    s = ("format:\n  html:\n    embed-resources: true\n    toc: true\n    toc-depth: 3\n"
         "    number-sections: true\n")
    if code_fold:
        s += "    code-fold: true\n"
    s += f"    code-overflow: wrap\n    fig-width: {fw}\n    fig-height: {fh}\n"
    return s

YAML = {
"methodology.Rmd": "---\n"
 'title: "Methodology"\n'
 'subtitle: "Exploratory Subgroup Identification with Bootstrap Bias Correction"\n'
 'author: "Larry F. Leon, Thomas Jemielita, Zifang Guo, Rachel Marceau West, Keaven M. Anderson"\n'
 "date: today\n" + fmt(True, 8, 6) +
 "bibliography: references.bib\n"
 "nocite: |\n  @Leon_2024\n  @Athey2019\n  @Wager_2014\n  @Leon_Anderson_2026a\n  @Leon_Anderson_2026b\n"
 + vig("ForestSearch Methodology") + "---\n",
"paper_simulations.Rmd": "---\n"
 'title: "Simulation Studies for Evaluating ForestSearch Performance"\n'
 'subtitle: "Operating Characteristics and Power Analysis"\n'
 'author: "ForestSearch Package"\n'
 "date: today\n" + fmt(True, 9, 6)
 + vig("Simulation Studies for Evaluating ForestSearch Performance") + "---\n",
"treatment_effect_definitions.Rmd": "---\n"
 'title: "Treatment Effect Definitions in ForestSearch"\n'
 'subtitle: "Marginal Effects, Controlled Direct Effects, Average Hazard Ratios, and Simulation Truth"\n'
 'author: "ForestSearch Package Documentation"\n'
 "date: today\n" + fmt(True, 8, 5)
 + vig("Treatment Effect Definitions in ForestSearch") + "---\n",
"biomarker_effects.Rmd": "---\n"
 'title: "Biomarker Treatment Effects in the Weibull AFT/Cox Framework"\n'
 'subtitle: "Theory, Estimands, and Causal Foundations"\n'
 'author: "ForestSearch Package Documentation"\n'
 "date: today\n" + fmt(False, 9, 6)
 + vig("Biomarker Treatment Effects in the Weibull AFT/Cox Framework") + "---\n",
"causal_effects_brief_review.Rmd": "---\n"
 'title: "Causal Interpretation of Hazard Ratios in Survival Analysis"\n'
 'subtitle: "A Quick Review for ForestSearch Users"\n'
 'author: "ForestSearch Package Documentation"\n'
 "date: today\n" + fmt(False, 8, 5)
 + vig("Causal Interpretation of Hazard Ratios in Survival Analysis") + "---\n",
}

# --- the three diagrams as native Quarto {dot} cells (DOT text verbatim, with the
#     two R-string conventions undone: '\\n' -> '\n' and single -> double quotes) ---
def dot_cell(cap, dot):
    return "```{dot}\n" f'//| fig-cap: "{cap}"\n' + dot.strip("\n") + "\n```\n"

M_ALG_OLD_START = '```{r algorithm-diagram, echo=FALSE, fig.cap="ForestSearch Algorithm Overview"}'
M_BOOT_OLD_START = '```{r bootstrap-diagram4, echo=FALSE, fig.cap="Bootstrap Bias Correction Workflow"}'
P_WF_OLD_START = '```{r sim-workflow-diagram, echo=FALSE, fig.cap="Four-stage simulation workflow for evaluating subgroup identification performance.", fig.height=2.0, fig.width=9}'

def chunk_at(text, start):
    i = text.index(start)
    j = text.index("\n```\n", i) + len("\n```\n")
    return i, j, text[i:j]

def grviz_dot(chunk):
    a = chunk.index('grViz("') + len('grViz("')
    b = chunk.rindex('"')
    dot = chunk[a:b]
    dot = dot.replace("\\\\n", "\\n").replace("'", '"')
    # strip the uniform indentation of the paper_simulations string
    lines = dot.split("\n")
    while lines and not lines[0].strip(): lines.pop(0)
    while lines and not lines[-1].strip(): lines.pop()
    ind = min(len(l) - len(l.lstrip()) for l in lines if l.strip())
    return "\n".join(l[ind:] if len(l) >= ind else l for l in lines)

def convert(name, with_dot):
    text = open(os.path.join(SRC, name), encoding="utf-8").read()
    if md5(text) != EXPECT[name]:
        sys.exit(f"{name}: md5 {md5(text)} is not the post-alignment text")
    _, body = split_yaml(text)
    if name == "methodology.Rmd":
        body = once(body,
            "\\varphi(w_2; \\beta, 8/d) \\, dw_1 \\, dw_2,\n$$\n",
            "\\varphi(w_2; \\beta, 8/d) \\, dw_1 \\, dw_2,\n$$ {#eq-power}\n", "eq-power")
        if with_dot:
            i, j, ch = chunk_at(body, M_ALG_OLD_START)
            body = body[:i] + dot_cell("ForestSearch Algorithm Overview", grviz_dot(ch)) + body[j:]
            i, j, ch = chunk_at(body, M_BOOT_OLD_START)
            body = body[:i] + dot_cell("Bootstrap Bias Correction Workflow", grviz_dot(ch)) + body[j:]
    if name == "paper_simulations.Rmd" and with_dot:
        i, j, ch = chunk_at(body, P_WF_OLD_START)
        body = body[:i] + dot_cell(
            "Four-stage simulation workflow for evaluating subgroup identification performance.",
            grviz_dot(ch)) + body[j:]
    if name == "treatment_effect_definitions.Rmd":
        body = once(body, "as discussed in Section \\@ref(theoretical-relationships).",
                    "as discussed in [Theoretical Relationships](#theoretical-relationships).",
                    "bookdown ref")
    return YAML[name] + body

os.makedirs(os.path.join(OUT, "qmd"), exist_ok=True)
os.makedirs(os.path.join(OUT, "qmd_grviz"), exist_ok=True)
for name in EXPECT:
    q = name[:-4] + ".qmd"
    out = convert(name, with_dot=True)
    open(os.path.join(OUT, "qmd", q), "w", encoding="utf-8").write(out)
    print(f"qmd/{q}  md5 {md5(out)}  lines {out.count(chr(10))}")
    if name in ("methodology.Rmd", "paper_simulations.Rmd"):
        out2 = convert(name, with_dot=False)
        open(os.path.join(OUT, "qmd_grviz", q), "w", encoding="utf-8").write(out2)
        print(f"qmd_grviz/{q}  md5 {md5(out2)}  lines {out2.count(chr(10))}")
