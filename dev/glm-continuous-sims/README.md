# GLM / continuous simulations

Closed-form references for the GLM/continuous pathway of `forestsearch`.
Design stage. No harness has been built.

## Read order

1. `HANDOFF_glm_continuous_simulations.md` — the originating handoff, committed
   verbatim. **Three of its claims are wrong**; read the NOTE before relying on
   any number in it.
2. `NOTE_target_is_collapsibility.md` — what replaces them, with measurements.
3. `SPEC_betaHhat_md.md` — the one specification implementable now.
4. `design-checks/` — the scripts behind every number quoted above.

## Status

| item | state |
|---|---|
| conditional estimand $\beta(\widehat H)$ for MD | closed form established, exact |
| `betaHhat_truth_md.R` | specified, not built |
| simulation aim | **open** |
| harm vs benefit | **open** |

## Open decisions

Both block the harness design; neither blocks `SPEC_betaHhat_md.md`.

- **What the simulation is for.** Unbiasedness of the correction for
  $\beta(\widehat H)$; coverage of the IJ interval; or degradation as the
  candidate family grows. Different DGMs and replicate counts.
- **Harm or benefit.** `quarto/simulations/actg175/continuous/actg175_continuous_simulations.qmd`
  searches for *benefit* by treatment switching (target MD$(Q) = +40$), so its
  realized regions are $\widehat G/\widehat G^{c}$, not $\widehat H/\widehat H^{c}$.
  Benefit work is deferred elsewhere in the project. Staying on benefit inherits
  that deferral; moving to harm needs a new continuous DGM first.

## About `design-checks/`

Standalone scripts that check a *claim*, against the law, in a sandbox. Base R,
seeded, runnable with `Rscript`. Several minutes each on one core.

**None has ever called `forestsearch()`.** They do not verify the package. The
handoff's own closing sentence is that the gap between a correct closed form and
what the code returns is exactly where the errors live — that gap is still open.

| script | claim | verdict |
|---|---|---|
| `oracle_design.R` | four properties of the multiplier oracle | Q1 refuted, Q3 vacuous |
| `check_betaHhat_md.R` | closed-form $\beta(\widehat H)$ on a synthetic DGM | reproduces; premise too narrow |
| `check_oracle_precision_claim.R` | is the oracle 24–32x sharper? | **no** — artifact of draw budgets |
| `check_md_collapsibility.R` | why MD is special | collapsibility, not independence |
| `check_md_target_exact.R` | exact target vs the handoff's formula | handoff's form low by $\delta$ |

Two cautions carried from the preceding cycle.

- `check_betaHhat_md.R` returns $\tau$ on four of its seven rows, and its
  negation row is byte-identical to its exact row (both normalise to the same
  bounds before evaluation). It does not discriminate negation handling, and an
  implementation returning $\tau$ whenever $\widehat H \subseteq H_{\mathrm{true}}$
  passes four rows.
- `oracle_design.R`'s Q3 ($\Sigma = \mathtt{crossprod(B)}$ under overlap) cannot
  fail: for $\xi \sim \mathcal N(0, I_n)$ and $D = B^{\mathsf T}\xi$,
  $\mathrm{Cov}(D) = B^{\mathsf T}B$ identically. Overlap cannot perturb an
  algebraic identity in the multipliers.

## Layout note

Artifacts built from these specs live under `quarto/`, not here, following
`quarto/simulations/actg175/binary/betaHhat_truth_glm.R`.
