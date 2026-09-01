# Memo (rev. 2) — using the applied OC evaluation to interpret the ACTG175 ForestSearch analysis

**Supersedes** `memo_oc_applied_takeaways_2026-08-31.md`, which (a) attributed the 58.5% declaration rate to a no-subgroup truth — that figure is the **sub-threshold-subgroup** null at the fixed gate `(10, 5)` — and (b) quoted a single tail probability for the observed statistic, when the correct statement is a range across the null family.
**Primary artifact:** `quarto/applications/actg175/analysis_actg175_continuous_oc.qmd` — self-contained, `c2 = 0.8·c1` policy, both nulls, the Q knob. **Companion:** `analysis_actg175_continuous_oc_evaluation.qmd`, the deep fixed-`c2` run at 2×10⁵ draws.
**Precision:** the self-contained document runs 2×10⁴ draws — MC SE ≈ 0.004 on every rate. Re-render at higher `draws` via its YAML params for archival precision.

The analysis found **Ĥ = {age ≤ 37 & cd40 > 507}**: 66 of 1,083 subjects, fitted mean CD4-decline difference **T̂ = 87.9**, consistency 0.95, against an ITT of **−27.6** (the combination benefits on average). The evaluation asks what this procedure, exactly as configured, does under truths anchored to this data, and therefore what may be claimed from what it found.

---

## The finding has three layers, and they earn different confidence

**1. The declaration event is close to uninformative.** At the analyst's screening floor with the scaled consistency floor (`c1 = 10, c2 = 8`), the search declares *some* harm subgroup **42.3% of the time when no subgroup exists at all**, and **50.4%** when a true subgroup exists with zero effect against the beneficial background. At the original fixed gate `(10, 5)` the latter is ~0.59. Two things follow. First, "ForestSearch identified a subgroup" carries almost no evidential weight at this configuration — it is near a coin flip under nothing. Second, the *gap between the two nulls is only ~8 points*, so the false-declaration rate is driven by multiplicity over M = 4,508 correlated candidates and a liberal threshold, not by a lurking weak subgroup. The thresholds (10, 5) are clinical-relevance floors, not error controls.

**2. The observed magnitude is where the evidence lives — and it is marginal across every null considered.** The selection-honest tail probability of the observed statistic, computed against the distribution of the *search's own maximum* so multiplicity is already inside it, is

| null truth | P(T ≥ 87.9) |
|---|---|
| no subgroup (homogeneous) | ≤ 0.038 (payload; bounded by the row below) |
| Q = Ĥ, zero effect | **0.038** (deep run at 2×10⁵ draws: 0.041) |
| Q = wider (1.96× prevalence), zero effect | **0.048** |
| Q = broad (3.01× prevalence), zero effect | **0.057** |

So the answer is a **range, 0.04–0.06**, not a single p-value — and it depends on how broad one is willing to believe the truly affected population might be. Under the broad reading the observation is *not* surprising even with zero harm. Quote the range; the finding is suggestive, and it is not decisive under any null in the family.

**3. The magnitude of the harm is weakly identified, and its lower bound collapses if Q may be broader than Ĥ.** The consonance curve, read at the analyst's own gate, crosses 0.05 at **q ≈ 8.9** CD4 units under the primary Q — already an order of magnitude below the naive 87.9, and below the analysis's own clinical floor of 10. Under the supersets the 0.05 crossing falls to **0.54 (wider)** and **≤ 0.01 (broad)**: if the harmed population is not exactly the sharp rule the search returned, the data are consistent with *arbitrarily small* harm. **Do not quote 87.9 as an effect size.** The de-biased (MR/bootstrap) estimates are the estimation-side response; the consonance curve is the design-side statement of why they are needed.

## The threshold policy, and what error control costs

`c2` fixed at 5 makes the consistency screen decorative once `c1` grows: at `c1 = 50` a candidate clearing the effect screen passes a consistency floor a tenth its size automatically. Tying **`c2 = 0.8·c1`** makes the screen bind — for small, noisy candidates `c2 + z_p·se_g` overtakes `c1` — and the payoff is large: 5% family-wise control arrives at **`c1 ≈ 40`** (no-subgroup null) or **`c1 ≈ 50`** (sub-threshold null), against **86** under the fixed floor. The type-I surface over the `(c1, c2)` crossing is the picture of this.

The policy applies to type-I and power, which are threshold-*policy* quantities. It deliberately does **not** apply to the calibration read-off, which interrogates an analysis already conducted at the analyst's own consistency floor. (At `c1 = T̂_obs` the two agree to four decimals: the effect screen dominates eligibility there.)

Even so, no configuration on this ladder is a test. Strict control still costs most of the power at every planted effect examined, so the finding is hypothesis-generating and needs prespecified confirmation — which this same machinery can size.

## The true-subgroup knob: breadth versus severity

The truth is a modelling choice, and the document makes it a knob — a named list of cutpoint sets. The primary plants Ĥ itself; two supersets keep the same factors with more permissive thresholds drawn from the family's own population grid: **wider** = {age ≤ 39 & cd40 > 449}, prevalence 1.96× Ĥ; **broad** = {age ≤ 39 & cd40 > 412}, 3.01×. (Both buy their breadth almost entirely on the `cd40` axis — the age grid's next step from 37 is 39, which barely moves prevalence.)

Each variant carries its **own** sub-threshold null: at `q → 0` a candidate's true mean is `beta_treat·(1 − PQg)`, which depends on the planted Q, and the between-variant differences reach 16–20 CD4 units. Across primary → broad at the observed-statistic rung, **power rises 0.58 → 0.95 while sensitivity falls 0.78 → 0.34**, and the calibration curves shift left: the 0.5 crossing moves **76.3 → 60.2 → 54.0**, the 0.05 crossing **8.9 → 0.54 → ≤ 0.01**.

That is the reading the section exists for: **the observed statistic constrains a trade between breadth and severity.** A narrow true subgroup must be severely harmed to make T̂ = 87.9 typical; a broad one need only be mildly harmed — and at the low end a broad one need barely be harmed at all. "What harm is consistent with what we saw" has no single answer until you fix *who* is harmed.

## Concerns to state alongside any use of these numbers

1. **Fixed-family conditioning.** The evaluation attaches to the no-front-end `maxeffCons` anchor — the configuration the MR inference uses, and on this data it selected the same rule as the lasso/GRF arms. The front-ends-on arms search a data-dependent family; their operating characteristics are a different, conditional-on-family estimand.
2. **The two-point mixture stylization.** Truths are one background effect and one planted effect on the data's own covariate rows, with a GLM baseline surface. Real heterogeneity is graded. The consonance curve is a model-based consonance statement — not a posterior, not an all-cause interval.
3. **The null anchors at a point estimate.** The background is the fitted −27.0 (SE ≈ 7.9), and that uncertainty is not propagated. Treat the 0.04–0.06 range as an order of magnitude, not a knife-edge.
4. **Q is not exactly in the discovery family.** Ĥ's `cd40` cut (507) sits off the population grid (nearest 506); the nearest member matches at purity 0.9969. Immaterial here, but every "power" figure is power to declare *some* rule.
5. **A documented orientation limitation.** The wrapper orients candidate means by the harm direction, which has no meaning at a homogeneous truth — hence the no-subgroup null is built by substituting the common effect into the family rather than through the generator. Recorded in the document's prose; not patched.
6. **Machinery offsets.** The wrapper's expected declared-subgroup size carries a ~2-subject offset in its verification fixtures, and the prevalence-scaled SE approximation was stress-tested (the `se_direct` sensitivity moved every headline ≤ 1.4 points; not adopted).

## Three sentences for a results section

*ForestSearch as configured (c1 = 10, c2 = 5, consistency ≥ 0.90) declares some harm subgroup in 42% of trials under a no-subgroup truth anchored to this dataset and 50% when a true subgroup exists below the threshold, so the declaration alone is not evidence; the observed maximum statistic of 87.9 has selection-adjusted tail probability 0.04–0.06 across that family of nulls, depending on how broadly the affected population is defined. The magnitude is weakly identified: statistics of this size arise with ≥ 5% probability from true harms of ≈ 9 CD4 units if the affected group is exactly the rule identified, and from essentially any harm if it is threefold broader, so the naive estimate should be read through the de-biased estimates and this calibration rather than at face value. Holding the family-wise false-declaration rate to 5% requires a screening threshold near 40–50 CD4 units with the consistency floor scaled to it, at which power is substantially reduced — the finding is hypothesis-generating and requires prespecified confirmation.*
