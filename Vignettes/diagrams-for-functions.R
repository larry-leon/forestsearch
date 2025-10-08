---
  title: "Side by Side DiagrammeR Diagrams"
output: html_document
---

  ```{r, echo=FALSE, results='asis'}
library(DiagrammeR)

# Diagram 1
diagram1 <- grViz("
digraph get_FSdata_helpers {
  graph [rankdir = LR]
  node [shape=box, style=filled, fillcolor=\"#e6f2ff\", fontname=\"Arial\"]
  edge [fontname=\"Arial\"]
  get_FSdata [label=\"get_FSdata\"]
  get_FSdata -> lasso_selection
  get_FSdata -> get_conf_force
  get_FSdata -> filter_by_lassokeep
  get_FSdata -> is_continuous
  get_FSdata -> cut_var
  get_FSdata -> process_conf_force_expr
  get_FSdata -> dummy
  get_FSdata -> dummy2
}")

# Diagram 2
diagram2 <- grViz("
digraph grf_subg_helpers {
  graph [rankdir = LR]
  node [shape=box, style=filled, fillcolor=\"#e6f2ff\", fontname=\"Arial\"]
  edge [fontname=\"Arial\"]
  grf_subg [label=\"grf.subg.harm.survival\"]
  grf_subg -> policy_tree
  grf_subg -> causal_survival_forest
  grf_subg -> double_robust_scores
  grf_subg -> aggregate
}")

# Print HTML for side-by-side display
cat('<div style=\"display: flex; gap: 20px;\">')
cat('<div style=\"flex: 1;\">')
print(diagram1)
cat('</div>')
cat('<div style=\"flex: 1;\">')
print(diagram2)
cat('</div>')
cat('</div>')


library(DiagrammeR)
library(gridExtra)

# Save diagrams as SVG
export_graph(grViz("digraph {...}"), file_name = "diagram1.svg")
export_graph(grViz("digraph {...}"), file_name = "diagram2.svg")

# Read SVGs and display side by side
library(grid)
img1 <- grid::rasterGrob(png::readPNG("diagram1.png"), interpolate=TRUE)
img2 <- grid::rasterGrob(png::readPNG("diagram2.png"), interpolate=TRUE)
grid.arrange(img1, img2, ncol=2)

# Only run if not already installed
# install.packages("manipulateWidget")
library(DiagrammeR)
library(manipulateWidget)

diagram1 <- grViz("
digraph get_FSdata_helpers {
  graph [rankdir = LR]
  node [shape=box, style=filled, fillcolor=\"#e6f2ff\", fontname=\"Arial\"]
  edge [fontname=\"Arial\"]
  get_FSdata [label=\"get_FSdata\"]
  get_FSdata -> lasso_selection
  get_FSdata -> get_conf_force
  get_FSdata -> filter_by_lassokeep
  get_FSdata -> is_continuous
  get_FSdata -> cut_var
  get_FSdata -> process_conf_force_expr
  get_FSdata -> dummy
  get_FSdata -> dummy2
}")

diagram2 <- grViz("
digraph grf_subg_helpers {
  graph [rankdir = LR]
  node [shape=box, style=filled, fillcolor=\"#e6f2ff\", fontname=\"Arial\"]
  edge [fontname=\"Arial\"]
  grf_subg [label=\"grf.subg.harm.survival\"]
  grf_subg -> policy_tree
  grf_subg -> causal_survival_forest
  grf_subg -> double_robust_scores
  grf_subg -> aggregate
}")

splitLayout(
  cellWidths = c("50%", "50%"),
  diagram1,
  diagram2
)

```{r, results='asis', echo=FALSE}
library(DiagrammeR)
library(manipulateWidget)

diagram1 <- grViz("
digraph get_FSdata_helpers {
  graph [rankdir = LR]
  node [shape=box, style=filled, fillcolor=\"#e6f2ff\", fontname=\"Arial\"]
  edge [fontname=\"Arial\"]
  get_FSdata [label=\"get_FSdata\"]
  get_FSdata -> lasso_selection
  get_FSdata -> get_conf_force
  get_FSdata -> filter_by_lassokeep
  get_FSdata -> is_continuous
  get_FSdata -> cut_var
  get_FSdata -> process_conf_force_expr
  get_FSdata -> dummy
  get_FSdata -> dummy2
}")

diagram2 <- grViz("
digraph grf_subg_helpers {
  graph [rankdir = LR]
  node [shape=box, style=filled, fillcolor=\"#e6f2ff\", fontname=\"Arial\"]
  edge [fontname=\"Arial\"]
  grf_subg [label=\"grf.subg.harm.survival\"]
  grf_subg -> policy_tree
  grf_subg -> causal_survival_forest
  grf_subg -> double_robust_scores
  grf_subg -> aggregate
}")

splitLayout(
  cellWidths = c("50%", "50%"),
  diagram1,
  diagram2
)







