# Hand-off: forestsearch paper + simulation design strategy

## What I want from this conversation
Organize a paper and design the simulation study that supports it. I want to
start at the *claims* level --- what the paper argues and for whom --- and let
the simulation design follow from that, rather than designing simulations first.

## Where the software stands (context, not the topic)
`forestsearch` (R, CRAN v0.1.0; v0.2.0 in development) implements exploratory
harm-subgroup identification for clinical trials. It now has three selection
engines under one interface: **consistency** (the default, with optional GRF
screening), **DINA**, and **GRF** (policy-tree selection, plus an experimental
Pareto-frontier mode). The survival pipeline has been extended to **GLM
outcomes** (binary/OR, continuous/MD, count/IRR). All three engines share an
aligned `sg_focus` vocabulary (eff / effMaxSG / effMinSG / maxSG / minSG), and
embedded-vs-standalone equivalence is verified. Grounded in León et al. 2024
(*Stat Med*, DOI 10.1002/sim.10163).

## Paper-level questions I want to work through first
1. **Central contribution.** Is the paper primarily (a) the methodological
   extension (GLM outcomes, the unified multi-engine interface, the aligned
   selection vocabulary), (b) an operating-characteristics / benchmarking study
   across engines, (c) an applied demonstration, or some weighting of these?
2. **Target venue and audience** (methods journal vs. applied/clinical-trials
   journal), and how that shifts emphasis.
3. **Relationship to León et al. 2024** --- what is genuinely new here vs. what
   extends that work, and how to position it.
4. **Methods vs. application split** --- how much real-data demonstration, and
   which datasets.

## Simulation design questions that follow
5. What **claims** must the simulations substantiate, and what **DGMs** (signal
   regimes, null regimes, effect scales) and **comparators** are needed to
   support each --- no more, no less.
6. How to structure **operating characteristics** (detection, the
   sens/spec/PPV/NPV classification quartet, behavior under the null, any
   threshold calibration such as `dmin.grf` and `n.min`) into a coherent study
   rather than a scatter of runs.
7. **Scale and reproducibility** --- production sim counts, the Quarto-first
   runner pattern, parallel design.

Let's start with #1 and #2 --- the contribution and venue --- since those
determine almost everything downstream. I'll share my current thinking and would
like you to push on it.
