---
bibliography: []
---

# Audit: does the code implement the manuscript's algorithm?

**For a fresh Claude Code session on `mr-in-replicates`, repo root
`~/Documents/GitHub/forestsearch`.**

This is a reading and tracing task. It produces findings, not commits. Do not
change any code until the findings are reviewed.

---

## Why this exists

Six defects were found in this codebase over one working session, and every one
had the same shape: **a quantity the manuscript defines once, reconstructed or
filtered at a second site.**

| defect | manuscript quantity | where it diverged |
|---|---|---|
| MR admission set | the identifier's selection domain | `t_g` rebuilt from raw thresholds at MR's call site |
| DINA/GRF consistency floor | no consistency screen exists for them | `c_consistency = 0` passed with `p_star` live, leaving `z·sigma_D` |
| GRF inclusion band | the `effMaxSG`/`effMinSG` band | hardcoded multiplicative band ignoring `selection_rule` |
| `maxeff`/`maxeffCons` frontier | the `sg_focus` dispatch | seven whitelisted values, five branches |
| `sigma_D` on the GLM path | `sqrt(sum(dfbeta^2))` | `stats::dfbeta()`, whose meaning changed underneath |
| MR family restriction | Algorithm Step 2's qualifying candidates | a second native-statistic band on top, discarding 99.84% of GRF's family |

Finding them one at a time cost a week. The purpose of this audit is to find the
remainder in one pass, or to establish that there are none.

**This is not a bug hunt for arbitrary defects.** It is a specific question:
for each quantity the manuscript defines, does the code compute it once, in one
place, as defined?

---

## Method — three rules, each learned by violating it

**1. Trace call sites, never formals or documentation.**

A declared default does not tell you the value used. `family_native_neighborhood`
defaulted to `NULL`, and its own documentation said `NULL` disables the
restriction — but both call sites substituted `effect_neighborhood` before
calling, so the restriction ran at 0.10 in every DINA and GRF analysis ever
performed. Auditing formals would have missed it entirely.

Likewise, a function's name and roxygen may describe behaviour it no longer has.
`.fs_mr_restrict_native()` bands on the *native* statistic per its name and
rationale, but under the current defaults it bands on `sel_effect`, the
inferential effect.

**2. A leaf check does not establish the wiring.**

30,000 numerical comparisons against `.grf_frontier_select()` passed while the
argument under test never reached it, because two intermediate argument-builders
lacked the formal. Verifying a function says nothing about whether its inputs
arrive. Where a claim depends on a value flowing from `forestsearch()` to a leaf,
trace every hop.

**3. Enumerate exhaustively; never `head`-limit a search you will draw a
conclusion from.**

A `grep | head -20` reported four MR enforcement sites when there were five —
the fifth appeared on the line after the cut. Use `grep -rn` unbounded and count.

---

## The reference definitions

Quoted from the manuscript. These are the standard; the code is the thing being
checked against them.

### ℱ — the candidate family

§2.1:

> the candidate family is every conjunction of at most two distinct conditions,
> ℱ = { φ : φ ∈ ℬ } ∪ { φ ∧ φ′ : φ, φ′ ∈ ℬ, φ ≠ φ′ }, |ℱ| ≤ L + (L choose 2).

> The bound in Equation 2 is not attained — empty intersections (the two sides
> of one cut) and candidates failing minimum subgroup-size and minimum-event
> guards are dropped

> Throughout, **fixed family** means precisely this: the list of subgroup
> definitions is one fixed object, and any resampling re-evaluates the same list
> rather than generating a new one (cut locations from observed quantiles are
> held fixed under all resampling).

Algorithm Step 2:

> For a model-based identifier, ℱ is its qualifying candidates and 𝒮 re-ranks
> them on β̂.

Supplement S2, DINA Step 4 (Qualify):

> Retain a candidate iff Σᵢ M_{s,i} ≥ n_min and τ̄_s ≥ m, where τ̄_s = (Σᵢ
> M_{s,i} τ̂ᵢ)/Σᵢ M_{s,i} is the region-mean effect.

Supplement S1, GRF:

> The remedy is identical: score the forest's **qualifying candidates** by the
> inferential effect β̂(g), demoting the forest to a candidate generator.

### β̂(g) and the influence db_{g,i}

§2.2:

> The estimate β̂(g) solves the within-subgroup score equation Σ_{i∈g} ψᵢ(β; g)
> = 0, and the per-subject empirical influence is the dfbeta
> db_{g,i} = ℐ_g⁻¹ ψᵢ(β̂(g); g) ≈ β̂(g) − β̂_(−i)(g), Σ_{i∈g} db_{g,i} = 0

### σ²_{D,g}

§2.2:

> The sum of squared influences, σ²_{D,g} = Σ_{i∈g} db²_{g,i}, is the robust
> (sandwich / infinitesimal-jackknife) variance of β̂(g)

### t_g — the admission threshold

§2.4:

> Collecting the screen and the calibrated consistency requirement into
> per-candidate thresholds t_g = max(c_screen, c_cons + z_{(1+p⋆)/2} σ_{D,g})

### 𝒮 — the selection map

§2.4:

> among the admitted {g : β̂(g) ≥ t_g}, take the largest effect, or the most
> consistent candidate, with subgroup size breaking ties — which we write as a
> selection map acting on the vector of candidate effects: Ĥ = 𝒮({β̂(g) : g ∈ ℱ})

> 𝒮 is deterministic given the candidate effects (with the σ_{D,g} entering only
> through the fixed thresholds t_g), and it ranks on the inferential effect — the
> same β̂(g) the report is about.

### The alignment condition

§4.2:

> the procedure's re-selection Ĥ*_b = 𝒮({β̂(g) + D_g^(b)}) (Equation 7) is the
> faithful linearization of the identifier's selection **exactly when S_g = β̂(g)**

> The de-biased estimate Equation 10 and its interval reproduce León's bootstrap
> to first order **if and only if** the identifier ranks on the inferential
> effect, S_g = β̂(g).

### The perturbation and re-selection

Algorithm Steps 3 and 4:

> 3. **Perturb all candidates jointly.** Draw centered, unit-variance multipliers
> Ξ (n × B) and form P = B_eff^⊤ Ξ, so one shared draw perturbs every candidate
> at once: β̂(g) + D_g^(b) with D_g^(b) = P[g, b].

> 4. **Re-select under each draw.** For b = 1, …, B, re-apply the rule to the
> perturbed effects, Ĥ*_b = 𝒮({β̂(g) + D_g^(b)}).

### The de-biasing terms

§3.1, Equations 9 and 10:

> On each draw the re-selected winner Ĥ*_b carries its own perturbation
> D_{Ĥ*_b}(b) = P[Ĥ*_b, b]. Because the rule rewards a large perturbed effect
> β̂(g) + D_g(b), conditioning on the winner draws its perturbation from the upper
> part of the candidates' distribution — a maximum-type order statistic over the
> competition — so it is positive on average. That average is the optimism of
> selection; alongside it sits a same-draws term evaluated at the fixed observed
> subgroup,
>
> bias_sel = (1/B) Σ_b D_{Ĥ*_b}(b),   bias_fix = (1/B) Σ_b D_{Ĥ}(b)      (9)
>
> and the de-biased estimate of the selected effect subtracts both,
>
> β̃(Ĥ) = β̂(Ĥ) − bias_sel − bias_fix                                    (10)

> The second is taken at the fixed Ĥ and so has expectation zero under centered
> multipliers (E D_Ĥ = 0); León et al. (2024) nonetheless keep it. [...] dropping
> it is harmless for the point estimate but not for its variance: mean-zero
> though it is, the same-draws discrepancy supplies the per-draw residual the
> interval below is built from. The complement subgroup Ĥ^c is de-biased
> identically.

### The interval and the IJ variance

§3.2, Equations 12, 13 and 14:

> Write the per-draw residual as the centered total bias,
>
> r_b = (bias_sel + bias_fix) − D_{Ĥ*_b}(b) − D_{Ĥ}(b)                     (12)
>
> — the total bias removed in Equation 10 minus the two perturbations that draw b
> contributed to it. The infinitesimal-jackknife variance and its finite-B
> correction are then
>
> Ṽ = Σ_i coṽ_i²,  coṽ_i = (1/B) Σ_b (K*_{bi} − K̄*_i) r_b,   V̂ = Ṽ − (n/B) r¯²   (13)
>
> where K̄*_i = B⁻¹ Σ_b K*_{bi} is subject i's mean multiplicity across draws and
> r¯² = B⁻¹ Σ_b r_b²; the subtracted term removes the upward Monte-Carlo bias of Ṽ
> at finite B.

> The same matrix delivers this with no second pass. The centered multiplicity
> K*_{bi} − K̄*_i is the row-centered multiplier, so stacking the residuals into
> r = (r_1, …, r_B)ᵀ and row-centering the draws, Ξ̃ = Ξ − rowmean Ξ,
>
> coṽ = (1/B) Ξ̃ r,   Ṽ = ‖coṽ‖²,   V̂ = Ṽ − (n/B) r¯²                        (14)

Algorithm Step 7:

> Report β̃(Ĥ) ± z_{1−α/2} √V̂, exponentiated to the hazard- or odds-ratio scale
> (Equation 15). The complement Ĥ^c is treated identically.

### The fixed-family requirement

§5:

> For the multiplier-resampling procedure to match it, the re-selection Equation
> 7 must range over the full enumerated family — every cut combination the search
> could form — re-applying the screen on each draw; restricting the perturbed
> competition to the handful of post-screen candidates would capture the rule's
> optimism but omit the screen's, under-correcting the bias.

---

## What to audit

For **each** quantity below: locate every site in `R/` that computes, filters,
re-derives, or overrides it. The count is the easy part and not the point — the
work is determining, where more than one site exists, whether they **agree with
each other** and whether each **matches the definition above**. Two sites that
agree are a consolidation opportunity; two that disagree are a defect; one that
departs from the definition is a defect regardless of how many there are.

Report a finding when any of these hold:

- a quantity is computed in more than one place
- a filter is applied that the manuscript does not describe
- a filter the manuscript describes is absent
- a value is substituted at a call site such that a declared default never applies
- a function's name or documentation describes behaviour it does not have
- a code path exists that no argument can reach, or an argument that no code path consults

### Phase 1 — identifier and MR paths

Report Phase 1 before starting Phase 2.

1. **ℱ, the candidate family.** Every site that builds, filters, truncates, or
   restricts the set of candidates — for each of the three engines, and
   separately for the identifier and for MR.

   One departure is already established, so do not spend time re-deriving it:
   `max_subgroups_search` caps how many candidates reach consistency
   evaluation, and the manuscript has no counterpart — §2.1's family is every
   conjunction of at most two conditions less empty intersections and
   size/event failures, and Algorithm Step 2 selects over that family. Record
   it, and report instead: at what point in the pipeline the cap is applied,
   whether anything downstream records that truncation occurred, and whether
   any other cap of the same kind exists elsewhere.

2. **β̂(g) and db_{g,i}.** Every site that fits a within-subgroup working model
   and extracts influences. Survival and GLM paths separately.

3. **σ²_{D,g}.** Every site that forms it. Confirm it is `sum(dfbeta^2)`
   everywhere it appears.

4. **t_g.** Every site that computes an admission threshold or filters candidates
   by one. This was fixed once, via `.fs_resolve_admission()` — verify no site
   still reconstructs it.

5. **𝒮, the selection map.** Every site that ranks candidates and picks a winner:
   identifier-side for all three engines, and MR's re-selection. For each
   `sg_focus` value and each engine, does the identifier's rule and MR's
   re-selection rule agree — ranking statistic **and** admissible domain?

6. **The inclusion band.** Every site that computes one. Should now be
   `.compute_inclusion_band()` only; verify.

7. **The alignment condition, S_g = β̂(g).** For every reachable combination of
   `subgroup_method`, `select_statistic`, `grf_selection` and `sg_focus`,
   determine the statistic the identifier ranks on and the statistic MR
   re-selects on. Report any mismatch.

### Phase 2 — the correction itself

Untouched so far, and worth the same scrutiny.

8. **P = B_eff^⊤ Ξ.** Verify the multipliers are centered and unit-variance for
   each `multiplier` option, that one shared draw perturbs every candidate, and
   that B_eff holds the influences of Equation 3.

9. **The de-biasing terms.** Equations 9 and 10, quoted above. Verify both
   terms are formed as stated, that `bias_fix` is retained rather than dropped,
   and that the complement is treated identically.

10. **The IJ variance.** Equations 12, 13 and 14, quoted above. Verify the
    per-draw residual is the centered total bias as written; that the finite-B
    correction subtracts (n/B) r¯²; and that the "same draws" property holds —
    the variance must use the identical P the point correction used, not a
    second draw.

### Phase 3 — replay paths

11. **The bootstrap.** `bootstrap_analysis_dofuture.R` reconstructs the
    `forestsearch()` call per replicate — the reconstruction pattern at its
    largest scale. Enumerate every argument it alters, drops, or nulls.

    The criterion is not "does the manuscript describe this alteration" — it
    describes almost none, because it does not specify bootstrap
    implementation. The criterion is: **does the alteration change what ℱ, 𝒮,
    t_g or σ_D are within a replicate, relative to the original fit?** An
    alteration that only affects parallelism, verbosity or plotting is not a
    finding. One that changes the family, the ranking, the admission set or
    the influence is.

12. **Cross-validation.** `forestsearch_Kfold()` and `forestsearch_tenfold()`,
    same criterion.

---

## Output

One Quarto document, `dev/identifier-alignment/code_theory_audit.qmd`, HTML.

Per quantity: the definition, a table of every site found with `file:line`, and
a verdict — **matches**, **duplicated**, **departs**, or **unreachable**. For
each departure, state what the code does, what the manuscript says, and the
direction of the consequence if determinable.

A summary table at the top: quantity, number of sites, verdict.

**Progress and stopping.** Report each phase when complete rather than holding
everything to the end. Within a phase, if one quantity is consuming
disproportionate effort, report what you have established, state precisely what
is unresolved, and move on — a partial finding with a clear boundary is worth
more than a complete one delivered late. Twelve quantities across the package is
a large surface; do not treat exhaustion of the list as the success criterion.

Do not fix anything. Do not commit code changes. If a departure looks trivial to
fix, still report it rather than fixing it — every defect in the table above
looked trivial in isolation, and several turned out to interact.

---

## Calibration

The six areas in the opening table were changed recently. **Audit them exactly
as you would any other — do not assume the changes were complete or correct.**
A fix that looked right and was not is the same failure this audit exists to
catch, and expecting a match is the prior most likely to let one through.

Use the six only to check your *method*: if your approach would not have
surfaced those six from a cold start, it will not surface the seventh.

---

## Ask rather than guess

- if a manuscript passage is ambiguous about what the code should do
- if a quantity appears in the code with no counterpart in the manuscript
- if the manuscript defines something with no counterpart in the code
- if a departure could be read as either a defect or a deliberate extension

The manuscript PDFs are at `dev/identifier-alignment/` — if not present, ask for
them rather than working from the definitions quoted here, which are extracts and
may lack context.
