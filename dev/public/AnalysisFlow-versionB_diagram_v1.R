library(DiagrammeR)

library(DiagrammeRsvg)
library(rsvg)

graph_code <- 'digraph flowchart {
  node [shape = rectangle, style=filled, fillcolor=lightgrey]
  Start [label = "Non-AP data for subgroup identification"]
  
  Decision1 [label="Comprehensive evaluation of 
  biomarker(BM) and all disease/demographic factors;
  BM cuts e.g. {<1},{<2},...,{<10}", shape = box, fillcolor=lightblue]
  
  #Analysis1 [label = "H with largest hr, min(hr)>=0.9", shape = rectangle, fillcolor=lightyellow]
  
  #Analysis2 [label = "Smallest H, min(hr)>=0.9", shape = rectangle, fillcolor=lightyellow]
  
  #Analysis3 [label = "Largest B, max(hr)<= 0.6", shape = rectangle, fillcolor=lightyellow]
  
  Analysis [label = "Subgroup identification algorithm (specifications)", shape= rectangle, fillcolor=lightyellow]
  
  Decision2 [label="Benefitting subgroup G* biologically plausible?", shape=rectangle, fillcolor=lightgrey]
  
  Question1 [label = "Operating characteristics:
  false-discovery, estimation (bias,CI)", shape=diamond, fillcolor=lightgray]
  
  APanalysis [label="Primary analysis of AP
  within (benefitting) G* sub-population", shape=diamond, fillcolor=green]
  
  Start -> Decision1
  
  #Decision1 -> Analysis1 [label="option"]
  #Decision1 -> Analysis2 [label="option"]
  #Decision1 -> Analysis3 [label="option"]
  
  Decision1 -> Analysis [label="Team goal"]
  
  #Analysis2 -> Decision2 [label="Team conclusion"]
  Analysis -> Decision2 [label="Team conclusion"]
  
  Decision2 -> APanalysis[label="Apply to AP"]
  
  APanalysis -> Question1 [label = "Data-driven 
  simulation (DdS):
  SG sizes per data;
  Outcome/censoring models
  based on data (Weibull fits)"]
  
}'

mygraph <- grViz(graph_code)

# export graph
export_svg(mygraph) %>%
  charToRaw() %>%
  rsvg() %>%
  png::writePNG("MRCT-Flow-versionB.png")


