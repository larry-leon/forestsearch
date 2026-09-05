# TASK_mr_vs_guohe_2026-09-04 — Addendum A (post-Gate-0)

Date: 2026-09-04. Applies from Stage 1 onward. Copy to `dev/tasks/` beside the task document and commit before continuing.

Stage 0 report (46cb9602) is accepted as green. This addendum resolves the two open items it raised and adds the resulting Stage 2 record columns. Nothing here changes the gate's arithmetic.

## A1. Definitions (from the companion paper's supplement, S4.1 (A6) and S4.5)

Quantities the gate already forms. Let `B_eff` be the n × K influence matrix (column g = the dfbeta vector db_{g,·} of candidate g, zero outside g's members), Ĥ the observed winner, and Ĥ*_b the winner on multiplier draw b.

1. Cross-covariance matrix of the candidate effects (finite-sample; the supplement's n-scaling is a limit device and cancels below):
   `Σ̂ = B_effᵀ B_eff`, so `Σ̂[g,g] = σ̂²_{D,g}` and `Σ̂[Ĥ,g] = Σ_i db_{Ĥ,i} db_{g,i}`.
   On a disjoint family Σ̂ is exactly diagonal; on the nested §5.2 family it is not.

2. Re-selection frequencies: `p̂_g = (1/B) Σ_b 1{Ĥ*_b = g}`, Σ_g p̂_g = 1.

3. **(A6) diagnostic — "retained mass"** (supplement S4.1, S4.4, verbatim quantity):
   `A6_mass = 2 Σ_g p̂_g Σ̂[Ĥ,g] + σ̂²_{D,Ĥ}`; the assumption is `A6_mass ≥ 0`.
   Also record the standardized `A6_mass_std = A6_mass / σ̂²_{D,Ĥ}`, which equals 3 under separation (p̂_Ĥ = 1) and is ≥ 1 on any family whose off-diagonal Σ̂[Ĥ,g] are non-negative — i.e. (A6) holds by construction on both replication families (disjoint: off-diagonals zero; nested: overlapping members give positive cross-terms). Report it anyway; it is free, and the report should state that it is structurally satisfied here rather than empirically tested.

4. **M_eff — effective competition** (supplement S4.5). "M_eff" is the supplement's own symbol, defined by `c(M_eff) = m̂(0)`, where:
   - `c(M)` is the expected maximum of M independent standard normals, extended to real M ≥ 1 by
     `c(M) = M ∫ z φ(z) Φ(z)^{M−1} dz` (monotone, c(1) = 0, c(2) = 1/√π = 0.5642, c(10) = 1.5388);
   - `m̂(0) = E[max_g ζ_g]`, ζ ~ N(0, R̂), with `R̂ = D^{−1/2} Σ̂ D^{−1/2}` the correlation matrix of Σ̂ (unit diagonal, so M_eff is scale-free; note this normalization in the record). Compute by Monte Carlo with a fixed seed and ≥ 1e5 draws (K = 2 has the closed form `m̂(0) = √((1 − r)/π)` for correlation r); record the MC standard error.
   - Solve `c(M) − m̂(0) = 0` for M on [1, K + 1] (e.g. `uniroot`). M_eff is the number of independent candidates that would produce the same tie optimism.
   - Implied tie residual on the log-HR scale (supplement S4.5, Theorem 2(iii)): `tie_resid_implied = (1 − 2^{−1/2}) · m̂(0) · σ̂_{D,Ĥ}`. In tie scenarios this is the theory's prediction for the bias that MR leaves in β̃(Ĥ); the Stage 3 report compares it with the observed mean(MR − θ_{k̂}) on the log scale.

   Identity checks for the implementation (add to Stage 1 §1b):
   - exchangeable R̂ with M = 10 and ρ = 0.3, 0.6, 0.9 must give M_eff = 6.2, 3.7, 1.8 (supplement S4.5; continuous-M solve gives 6.22, 3.65, 1.80);
   - a disjoint family (R̂ = I) must give M_eff = K to within MC error;
   - flipping the treatment coding (the §5.1 orientation adapter) must negate every β̂(g) and every column of `B_eff` and leave Σ̂ unchanged, to relative tolerance 1e-10.

## A2. Exposure of the per-draw winners (Stage 0 finding 1)

Disposition requested from Larry as **D6**; recommendation in brackets.

- [Recommended] Add-only return: a new argument on `fs_mr_inference()` (working name `return_reselection = FALSE`) that, when TRUE, adds the per-draw winner index vector (length B) and the tabulated `p̂` over the family to the return list. Default FALSE reproduces the current return object exactly; the existing tests pass unchanged; the winner vector is already computed, so nothing in the arithmetic moves. R/ classification: adds code; moves nothing; changes no existing behaviour; changes no method.
- Alternative: the Stage-1 wrapper over the unchanged internals, if the return-list edit would touch anything beyond the return list. Either way, the Stage 1 record names the function and the lines changed.

## A3. Notes carried into Stages 2–3

- The comparison is against the *replicated* Guo–He numbers (same code, same replicates), not the paper's printed tables. The §5.2 record documents a ~0.010 conservatism deficit in the replicated proposed-method coverage in 4 of 30 cells; the Stage 3 report footnotes this so that MR-vs-G&H coverage differences of that order in §5.2 are not over-read.
- Stage 2 per-replicate record gains: `p_hat_H` (winner's own re-selection frequency), `p_hat_top3` (labels and values), `Sigma_HH`, `A6_mass`, `A6_mass_std`, `m0_hat`, `m0_mc_se`, `M_eff`, `tie_resid_implied`.
- Stage 3 adds, per scenario: mean M_eff; and in tie scenarios, `tie_resid_implied` (mean) against the observed log-scale retained bias of MR, with a Monte Carlo interval.
