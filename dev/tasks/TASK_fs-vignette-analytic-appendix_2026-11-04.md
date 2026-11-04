# TASK_fs-vignette-analytic-appendix_2026-11-04 — expand and align the vignette's analytic appendix

Repo: **forestsearch** (Pop!_OS box). File touched: `vignettes/articles/mrct-consistency.qmd` ONLY —
**R/ is not touched** (CC reads R/ and NAMESPACE read-only for Part C's table; no code changes; the
article is a pkgdown article, not a CRAN vignette, so the CRAN build is unaffected).
`export GIT_TERMINAL_PROMPT=0`; CC never pushes. Standing first action: copy this task document from
`~/Downloads/` into `dev/tasks/` (create if absent) and commit it alone. Short-lived branch off
`main` (check `git branch -r` first), fast-forward, delete. NOTES to `~/Downloads`.

Origin: Larry's read (2026-11-04) — Section 6 "Appendix: the analytic tools" (`{#sec-appendix}`,
currently the block from "# Appendix: the analytic tools" to end of file) is tersely written; expand
it and align it with the NEW manuscript §3.5. The alignment target below was transcribed by chat on
2026-11-04 from the current manuscript build (sim-manuscript-complete.pdf, the §3.5 produced by
TASK_mrct-sim-s35-si_2026-11-04); it is the spec of record for this task — CC does not need and
cannot read the mrct-consistency repo.

## Part A — structure and pointers

1. Replace the appendix's opening sentence (currently: condenses techreport Chapter 6 + "theory
   notes AHR_mrct_dgm.qmd and AHR_complete.qmd of the mrct-consistency repository") with: this
   appendix condenses §3.5 of the accompanying manuscript and its **Analytic Derivations**
   supporting document, where the derivations are (sections cited inline below as §B.x); keep the
   repository theory notes as a one-clause parenthetical (the development record). Exact citation
   wording for the manuscript is "the accompanying manuscript (in submission)" for now — Larry will
   set the final citation when the venue/preprint id exists.
2. Adopt the manuscript's role-tag convention for the appendix's displays — bold bracketed
   **[identity]**, **[computed]** (naming the forestsearch function), **[assumption]**,
   **[interpretation]** — exactly as in the transcription block below.
3. Name the modelling assumptions by the manuscript's labels (A1) 1:1 allocation, (A3) uniform
   staggered entry with calendar-time analysis and entry independent of event time, (A5) joint
   asymptotic normality — matching the uses below; the full list (A1)–(A6) lives in Analytic
   Derivations §B.1.3 (say so once).
4. Reorganize the existing bold paragraphs under run-in subheads matching the manuscript's: keep
   "The DGM in closed form" (vignette-specific, stays first), then "The clock and the arm
   functions", "The window function", "The expected-event measure of the window", "The
   events-weighted average over the window", "The least-false limit", "The events form", "The
   clock identity and the frailty reading", "The design layer", "Per-patient entry". ALL existing
   appendix content is KEPT (the Breslow-limit ratio, the first-order ω form with π(θ), the
   frailty reading, per-patient entry) — reorganized and expanded, nothing dropped.

## Part B — the aligned displays (transcription block = the spec)

Where the appendix already has one of these in inline or prose form, promote it to the displayed
form below with its tag and derivation pointer; where it is absent, add it. Symbols must match
exactly (set 𝒮; arms z ∈ {0,1}; u the operational clock; t_e, t_a; u_κ). Displays the appendix has
that match already (AHR/CDE, events form, design layer) get only the tag + pointer added.

1. **[identity]** Clock (derivation: Analytic Derivations, §B.3a):
   u = Λ0(s) = (s e^{-μ})^{1/τ};  P(Λ0(T) > u) = P(T > Λ0^{-1}(u)) = e^{-u e^{θ_z}} —
   every patient exponential with rate e^{θ_z} on the operational clock.
2. **[computed]** `ahr_limit()` — arm functions (derivation: §B.3a):
   A_z(u) = |𝒮|^{-1} Σ_{i∈𝒮} e^{-u e^{θ_{z,i}}};  f_z(u) = −A_z′(u) = |𝒮|^{-1} Σ_{i∈𝒮} e^{θ_{z,i}} e^{-u e^{θ_{z,i}}}
   — set-averaged arm survivor and event-density functions; arm z's marginal survivor in on-study
   time is A_z(Λ0(s)).
3. **[identity]** Window function (derivation: §B.3b). With entry E uniform on (0, t_e),
   independent of event time, analysis at calendar t_a: H(v) = min{1, max(v, 0)/t_e},
   Λ0^{-1}(u) = e^{μ} u^{τ}, and
   H_w(u) = H(t_a − Λ0^{-1}(u)) = 1 for 0 ≤ u ≤ u_κ;  (t_a − e^{μ}u^{τ})/t_e for u_κ < u < Λ0(t_a);
   0 for u ≥ Λ0(t_a),  with u_κ = Λ0(t_a − t_e) the entry kink — the last operational time at which
   every patient is still in their window; beyond it H_w falls (linearly in s) and is 0 from
   Λ0(t_a) on.
4. **[identity]** Expected-event measure (derivation: §B.3a for dN, §B.5 for its total):
   dN(u) = ½ H_w(u){f_1(u) + f_0(u)} du;  ∫dN = ½(D_1 + D_0);  D_z = ∫ H_w f_z du.
   Explanations to carry: the product H_w·f_z is exact because entry and event time are
   independent (A3); the ½ averages the arms under 1:1 allocation (A1); D_z is arm z's per-patient
   probability of an event in the window; ∫dN is the expected number of window events per patient
   of 𝒮; a trial of n patients with n_𝒮 ~ Bin(n, π_𝒮) in 𝒮 expects n π_𝒮 ∫dN window events from 𝒮
   — the expected-events number the layer reports (derivation: §B.1.2 and §B.8.2).
5. **[identity]** Events-weighted average (derivation: §B.5): ⟨g⟩_N = ∫ g dN / ∫ dN — and the
   sentence: every windowed summary of the layer is an integral or an average under dN (the score
   and information, the events form, the imbalance moments).
6. **[computed]** `ahr_population()`, `ahr_limit()` — least-false limit (derivation: §B.3a and
   §B.3c): U_𝒮(β) = ½ ∫_0^{Λ0(t_a)} H_w(u){f_1(u) − r(u, β)[f_1(u) + f_0(u)]} du,
   r(u, β) = e^β A_1(u) / {e^β A_1(u) + A_0(u)};  margHR(40) = e^{β*_𝒮} at t_a = 40; root unique
   because ∂U_𝒮/∂β = −∫ r(1 − r) dN < 0 (Struthers–Kalbfleisch written on the operational clock).
7. **[computed]** `ahr_limit()` — margHR(∞) (derivation: §B.3c):
   U_𝒮^∞(β) = ½ ∫_0^∞ {f_1(u) − r(u, β)[f_1(u) + f_0(u)]} du — the same root with no censoring,
   H_w ≡ 1 and the integral run to infinity; it depends on the set's patients alone, not on μ or τ.
8. **[computed]** `ahr_limit()` — AHR and CDE closed forms (derivation: §B.3c) — the appendix's
   existing displays, kept, with ψ = θ_1 − θ_0, p_𝒮 the set's H* patient share and
   ω⁰_𝒮 = Σ_{𝒮∩H*} e^{θ_0} / Σ_𝒮 e^{θ_0} the baseline-hazard share:
   AHR(𝒮) = exp{|𝒮|^{-1} Σ_𝒮 ψ} = b^{1−p_𝒮} h^{p_𝒮};  CDE(𝒮) = Σ_𝒮 e^{θ_1}/Σ_𝒮 e^{θ_0} = (1 − ω⁰_𝒮) b + ω⁰_𝒮 h.
9. Events form — existing content kept; tag the exact score **[identity]** and the expansion
   **[interpretation]**, with δ(u) = log{A_1(u)/A_0(u)}, P̄ = D_1/(D_1+D_0):
   P̄ = ⟨expit(β*_𝒮 + δ)⟩_N;  β*_𝒮 ≈ log(D_1/D_0) − δ̄ − ½(1 − 2P̄) Var_N(δ), δ̄ = ⟨δ⟩_N,
   Var_N(δ) = ⟨δ²⟩_N − δ̄²; the approximation ladder and its tabulation against the exact root are
   Analytic Derivations §B.6. The vignette's first-order ω form with π(θ) stays, cited to the
   Analytic Derivations (no § number supplied — cite the document).
10. **[identity]** Clock identity (display it; the appendix currently states it in prose):
    A_z^{(r)}(u) = A_z^{(0)}(e^r u);  f_z^{(r)}(u) = e^r f_z^{(0)}(e^r u) — a regional prognostic
    shift r is a clock on the region; margHR(∞), AHR, CDE do not depend on r; only the windowed
    margHR(40) moves; a region delayed by a is the undelayed region analysed at horizon t_a − a.
11. Design layer — existing content kept; tag the joint-normality display **[assumption]** (A5)
    and the probability displays **[computed]** `ahr_design()`; z_A = (−z_{0.975} σ̃_all − μ_all)/σ_all,
    z_B = −μ_R/σ_R, z_C = −μ_{G*}/σ_{G*}; P(trigger) = Φ₂(z_A, −z_B; −ρ_{all,R}); oracle joint
    Φ₃(z_A, −z_B, z_C; **R**). Cite the Analytic Derivations (no § number supplied — cite the
    document). Keep the vignette's variance/correlation material as is, same treatment.

Derivation-pointer rule: use ONLY the §B.x labels given above; for appendix material with no label
supplied here (Breslow limit, π(θ) form, per-patient entry, design-layer variances), cite
"the Analytic Derivations" without a section number — do not guess section numbers.

## Part C — function ↔ object table

Close the appendix with a small table mapping the analytic-layer exports to what they compute and
where the derivation lives: `ahr_population()`, `ahr_limit()`, `ahr_design()`,
`scale_dgm_hazard()`, plus any other analytic-layer export the appendix or manuscript §7 names —
verify the exact exported names and one-line purposes against the CURRENT R/ source and NAMESPACE
(read-only; verify from source, never from this description). Columns: function | computes
(object/display name from Part B) | derivation (§B.x or "Analytic Derivations"). If a function
named here is not an export in the current source, STOP and record (do not invent).

## Gates (stop-on-failure)

- G1 Article renders green by the repo's standard route for pkgdown articles; record the wall-clock
  in NOTES. The article's existing record-check and verbal-claim chunks must pass unchanged.
- G2 Scope: `git diff --name-only` confined to `vignettes/articles/mrct-consistency.qmd` (plus
  `dev/tasks/`); nothing under `R/`, `NAMESPACE`, `DESCRIPTION`, `man/`, `tests/`.
- G3 Spec fidelity: every display added or promoted matches Part B symbol-for-symbol; NOTES carries
  a checklist of the 11 items with done/already-present per item.
- G4 No invented pointers: grep the appendix for "§B." — every occurrence is one of
  {B.1.2, B.1.3, B.3a, B.3b, B.3c, B.5, B.6, B.8.2}.
- G5 Part C functions all exist as exports (checked against NAMESPACE).
- G6 Clean close: fast-forward `main`, delete branch, `git status --porcelain` empty.

NOT in this task: any edit in the mrct-consistency repo; the manuscript-citation wording beyond the
neutral placeholder (Larry sets it when the venue/preprint id exists — expected as a later one-line
amendment or a port-time edit).
