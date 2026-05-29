# DINA Inference Design Notes

Running notes on the theoretical decisions behind the
`dina_subgroup` / `dina_subgroup_bootstrap` inference machinery.
Maintained across the implementation work so a future methodology
Quarto can be written from these in one pass.

## 1. The forestsearch (FS) inference framework

Reference: León et al. (2024) *Statistics in Medicine* (DOI:
10.1002/sim.10163), Section 3.2.

FS uses a non-cross-fit Cox model fit within the selected subgroup
$\hat H$. The point estimate $\hat\beta(\hat H)$ is therefore subject
to *winner's curse* bias: the subgroup was selected because its Cox
log-HR exceeded a screening threshold, so the within-subgroup estimate
on the same data is upward-biased.

FS corrects this by a double-bootstrap construction (eq. 7):

$$\hat\beta^{*}(\hat H)
   = \hat\beta(\hat H)
     - \frac{1}{B}\sum_{b=1}^{B}\left[
        \eta^*_b(\hat H^*_b) + \eta^*_b(\hat H)
       \right],$$

with two bias terms:

  * $\eta^*_b(\hat H^*_b) = \hat\beta^*_b(\hat H^*_b) -
        \hat\beta(\hat H^*_b)$  — discrepancy between
    bootstrap-fitted-on-bootstrap-subgroup vs
    original-data-fitted-on-bootstrap-subgroup.
  * $\eta^*_b(\hat H) = \hat\beta^*_b(\hat H) - \hat\beta(\hat H)$
    — discrepancy on the *fixed original* subgroup.

Variance is via the infinitesimal jackknife (Wager–Hastie–Efron 2014),
treating $\hat\beta^{*}(\hat H)$ as a bagged estimator (eq. 9).


## 2. Why FS's bias correction is unnecessary for DINA

Two structural differences make DINA's inference setting
fundamentally different:

### 2a. Cross-fitting eliminates the first-order bias of $\hat\beta$

DINA estimates $\hat\beta$ by cross-fitting (Chernozhukov et al.
double-debiased ML): nuisance estimation (propensity, baseline) is
performed on training folds, and $\hat\beta$ is estimated on
held-out test folds. This breaks the dependency between the data
used to fit $\hat\beta$ and the data used to construct the score
function, giving asymptotic unbiasedness even for data-dependent
contrasts $a^\top \hat\beta$.

The two bias terms in eq. (7) of León et al. target exactly the
bias that cross-fitting eliminates. Applying both is correcting a
bias that is no longer there.

### 2b. dina_subgroup's target is at the boundary by construction

In FS, $\hat\beta(\hat H)$ is the *log-HR within* the subgroup —
the same quantity used to *select* $\hat H$ (above the screening
threshold of 1.25). Winner's curse is the dominant effect, and the
point estimate is upward-biased.

In `dina_subgroup`, $\bar{\hat\tau}_{S^*}$ is the *mean* of a
linear functional $a_S^\top \hat\beta$ over the selected subgroup,
where $a_S$ depends on the search-step choice $(j^*, q^*)$. The
search rule is "largest subgroup whose mean $\hat\tau$ exceeds
$m_{\text{diff}}$" — so the chosen subgroup is at the *boundary*:
adding one more patient would push the mean below $m_{\text{diff}}$.

Worked example. With DGM $\tau(x) = 0.4 + 1.2\, x_1$ on
$X \sim \mathrm{Unif}[-1,1]$ and $m_{\text{diff}} = 0.8$:
the population-optimal threshold solves
$0.4 + 1.2 \cdot \mathbb{E}[x_1 \mid x_1 \ge q] = 0.8$,
giving $q = -1/3$ and population mean $\bar\tau = 0.8$ exactly.

The population target *is* $m_{\text{diff}}$. There is no winner's
curse to correct because the estimand is at the boundary by
construction. The smoke-test result $\bar{\hat\tau} \approx 0.801$
across families is correct, not biased.


## 3. Why the "selection-adjusted CI" mechanically failed

The first implementation reported the percentile CI of the
bootstrap distribution of $\bar{\hat\tau}_{S^*_b}$ — each
iteration's mean over its bootstrap-selected subgroup. Per Section 2b,
each iteration's value is also pinned at the boundary, so the
distribution clusters tightly above $m_{\text{diff}}$ with widths
~0.005–0.015. The CI is statistically correct as a CI on
"the procedure's output value" but is **uninformative**: it
quantifies boundary slack (granularity of the threshold grid),
not effect uncertainty.

The hard floor $\bar{\hat\tau}_{S^*_b} \ge m_{\text{diff}}$ holds
by construction for every iteration, so the percentile lower
endpoint can never fall below $m_{\text{diff}}$. A "math-sanity"
assertion enforcing this is a tautology and was correctly noted as
hollow in retrospect.


## 4. The right inference target: fixed-subgroup effect

The clinically interpretable estimand is the linear functional
$a^{*\top} \beta$ where $a^* = (1, \bar{x}_{S^*})$ is the
covariate-mean vector of the *original*-data subgroup. This is the
expected treatment effect for patients whose covariates lie in the
discovered subgroup, treating the subgroup definition as fixed.

Bootstrap procedure: for each resample $b$, refit DINA → get
$\hat\beta_b$, and compute $T_b = a^{*\top}\hat\beta_b$ with $a^*$
held fixed from the original data. Percentile CI on
$\{T_b\}_{b=1}^B$.

Width of this CI reflects sampling variance of $\hat\beta$,
*not* boundary slack. Order of magnitude check on the smoke-test
DGM: with x1-coefficient SE $\approx 0.15$ and $\bar x_1 \approx 0.35$,
SE of $a^{*\top}\hat\beta \approx \sqrt{0.15^2 + 0.35^2 \cdot 0.15^2}
\approx 0.16$, giving a 95% percentile CI roughly $\pm 0.31$
around the point estimate — sensible.

The IJ-variance Wald CI returned by `dina_subgroup(dina_bagged_fit)`
targets the same quantity asymptotically. The percentile bootstrap
adds three things: (i) non-parametric construction (no asymptotic
normality assumption), (ii) the per-bootstrap selection-frequency
diagnostic on $(j_b, \text{dir}_b)$, and (iii) bootstrap distributions
on the structural quantities ($n_{\text{subgroup}}$ and threshold $q$).


## 5. Stability diagnostics

Three additional outputs complement the effect CI:

**selection_frequency** — counts of how often each
$(\text{covariate}, \text{direction})$ pair was selected across
the $B$ bootstraps. The diagnostic for "is the discovered subgroup's
defining covariate robust?" In the smoke test, $x_1$ was selected
in 96–100% of iterations across families, indicating the
*selection* itself is stable even though the *effect* CI is
non-degenerate.

**n_subgroup_ci** — percentile CI of $n_{\text{subgroup},b}$
across bootstrap iterations that selected the modal covariate.
Tells the user "how stable is the size of the discovered subgroup?"

**threshold_ci** — percentile CI of $q_b$ across bootstrap
iterations that selected the modal covariate. Tells the user
"how precisely is the boundary located?"

The restriction to the modal covariate is necessary because
thresholds on different covariates are not comparable
(different scales, different meanings).


## 6. When the full FS-style bias correction would be warranted

The argument in Section 2a depends on cross-fitting working
correctly. If propensity or baseline models are misspecified, or if
the search space is much larger than what we've tested, the
cross-fitting guarantee weakens and the FS-style two-bias-term
correction becomes more relevant. This would be a `bias_correct = TRUE`
option in a future version, not the default.

## 7. Related work — positioning in the HTE literature

The literature on data-driven subgroup identification, ITRs, and HTE
estimation is fragmented across three vocabularies (subgroup
identification in biostatistics; individualized treatment regimes /
rules in statistics and precision medicine; policy learning in
econometrics). The canonical taxonomy is Lipkovich, Dmitrienko &
D'Agostino (2017, *Stat Med* 36:136–196), substantially updated in
Lipkovich, Svensson, Ratitch & Dmitrienko (2024, *Stat Med*
43(22):4388–4436). They sort methods into four classes:

* **Global outcome modeling** — fit $Y \sim f(W, X)$, derive
  $\tau(x)$ as a derived quantity (S-, T-, X-learners).
* **Global treatment-effect modeling** — model $\tau(x)$ directly
  (causal forest, R-learner, DR-learner, **best_linear_projection**).
* **Optimal treatment regimes** — estimate the decision rule
  $d^*(x)$ directly (OWL, A-learning, Q-learning, **policy_tree**).
* **Local modeling** — recursive partitioning (interaction trees,
  SIDES, Virtual Twins, MOB, Causal Interaction Trees).

DINA sits in *global treatment-effect modeling*. The `dina_subgroup`
threshold-search step crosses into *local modeling*. forestsearch
sits in *local modeling* with an unusually exhaustive search.

### 7.1 GRF-native tools

The grf and policytree R packages provide three relevant pieces of
infrastructure:

**`causal_forest` / `causal_survival_forest`** — non-parametric
$\hat\tau(x)$ via forest weights with the Robinson decomposition
(causal_forest) or AIPW-style censoring corrections
(causal_survival_forest, Cui et al. 2023, *JRSS-B*).

**`best_linear_projection`** (Semenova & Chernozhukov 2021,
*Econometrics Journal*) — given a fitted causal forest and a
covariate vector $A_i$, returns a doubly-robust estimate of the
linear model $\tau(X_i) \approx \beta_0 + A_i^\top \beta$. The point
estimates come from a regression of the AIPW scores on $A_i$, with
HC3-robust SEs. **This is the direct GRF counterpart to DINA's
estimand**: both target a linear functional of the CATE.

**`rank_average_treatment_effect`** (RATE) and the **Targeting
Operator Characteristic (TOC)** curve — Yadlowsky, Fleming, Shah,
Brunskill & Wager 2025, *JRSS-B*. Given a prioritization rule
$S(X_i)$ (which could be a CATE estimate from any source), the TOC
plots $\mathbb{E}[\tau(X) \mid F_S(S(X)) \ge 1-q] - \mathrm{ATE}$
against $q$, and the RATE is a weighted integral. Provides a
hypothesis test for "does the prioritization rule effectively
stratify the sample". Doesn't pick a subgroup, but evaluates how
well any HTE estimator (including DINA) does so.

**`policytree::policy_tree`** (Athey & Wager 2021, *Econometrica*
89:133–161) — doubly-robust empirical welfare maximization over
shallow decision trees. Takes AIPW scores from a causal forest and
exhaustively searches for the depth-$d$ tree maximizing average
welfare. Produces an interpretable rule like
"treat if `age <= 50 AND er > 0`". The GRF ecosystem's answer to
"what's the best simple treatment policy given a fitted forest?".

### 7.2 Closely related non-GRF approaches

**CAPITAL** (Wang, Cai, Zeng, Song & Kosorok 2021, arXiv:2110.05636)
— constrained policy tree search. Explicitly targets *the largest
subgroup with average treatment effect at or above a clinically
meaningful threshold*. This is structurally **the same problem
`dina_subgroup` solves**, with a tree-search solver instead of a
univariate threshold search.

**SIDES / SIDEScreen** (Lipkovich et al., implemented in
`SIDES` / `subgroup.id` R packages and commercial tools) —
recursive partitioning generating multiple candidate subgroups,
each tested for treatment-by-subgroup interaction with multiplicity
control. The historical workhorse in pharma. Conceptually closest
to forestsearch in spirit (enumerate-and-screen) but uses recursive
partitioning rather than exhaustive combinatorial enumeration.

**Virtual Twins** (Foster, Taylor & Ruberg 2011, *Stat Med*
30:2867–2880; R package `aVirtualTwins`) — fits a T-learner CATE
estimator, then runs RPART on the predicted $\hat\tau(x)$ to find
subgroups. Two-step procedure with a tree as the second step.
forestsearch benchmarks against this in its 2024 SIM paper.

**Model-Based Recursive Partitioning (MOB)** (Seibold, Zeileis &
Hothorn 2016, *Biom J*; R package `partykit`) — partitions where
treatment-model *parameters* differ, not where the predicted
$\hat\tau(x)$ differs. A different formalization of "subgroup
identification" but related target.

**Causal Interaction Trees** (Su et al. 2009; Steingrimsson & Yang
2020, arXiv:2003.03042) — extends classical interaction trees with
three doubly-robust subgroup-specific estimators (IPW, g-formula,
DR). Provides theoretical guarantees for the observational case.

**`personalized` R package** (Huling & Yu 2021, *JSS* 98(5):1–60) —
unified framework with several subgroup-identification approaches
and built-in CV-based comparison across methods. Useful for
benchmarking.

**Best-selected subgroup inference** (Guo & He 2021, *JASA*
116:1498–1506) — bootstrap calibration for inference when the
subgroup is chosen by maximum effect estimate over a prespecified
set. Cited in León et al. 2024. Targets a different selection
rule than forestsearch (best-effect rather than maximally
consistent).

**R-learner** (Nie & Wager 2021, *Biometrika* 108:299–319) — meta-
learner for $\tau(x)$ using the Robinson decomposition with
arbitrary base learners. **DINA is essentially a parametric,
GLM-link-aware instantiation of the R-learner.**

**DR-learner** (Kennedy 2023, *Electronic J Stat* 17:3008–3049) —
doubly-robust meta-learner with sharper convergence guarantees than
the R-learner under weaker conditions. Same general structure.


## 8. Positioning verdict for DINA + `dina_subgroup` + forestsearch

The package combines three pieces with different lineages. The
honest positioning of each:

### 8.1 DINA itself

**Closest comparator: parametric R-learner / DR-learner with
exponential-family link functions; equivalently, a *parametric*
counterpart to `best_linear_projection`.**

Common features with the R/DR-learner family:
* Cross-fitting (Chernozhukov et al. 2018, *Econometrics J*) to
  break the bias from estimating nuisance functions on the same
  data used to fit the treatment-effect model.
* Robinson-style nuisance subtraction (in the Gaussian case).
* Asymptotic normality of $\hat\beta$ for inference on linear
  functionals.

Differences:
* DINA fixes a *parametric linear functional* $\tau(x) =
  \beta_0 + x^\top \beta$ on the natural-parameter scale. The
  R-learner and DR-learner allow flexible base learners (any ML
  method) for $\tau(\cdot)$.
* DINA extends naturally across the GLM family — log-HR for Cox,
  log-IRR for Poisson, log-OR for binomial, mean for Gaussian —
  in a unified framework. Most R-learner implementations are
  Gaussian-only or require manual adaptation per family.
* DINA carries sandwich variance and (via `dina_bagged`)
  infinitesimal-jackknife variance for free.

**The closest single function in existing software is
`best_linear_projection`**, but the targets differ in an important
way. `best_linear_projection` regresses non-parametric AIPW scores
from a fitted forest onto user-chosen covariates — it's a *post-hoc
linear summary* of a non-parametric estimate. DINA fits the linear
form *upfront*. The two should agree closely when the true CATE is
approximately linear on the natural-parameter scale and disagree
when it isn't.

Practical implication: an instructive head-to-head benchmark would
fit DINA and run `causal_forest()` followed by
`best_linear_projection()` on the same data with the same covariate
set, compare estimates and SEs. If they agree, DINA's parametric
assumption is defensible on that dataset; if they diverge, the true
CATE is non-linear and DINA's interpretation needs care.

### 8.2 `dina_subgroup`

**Closest comparator: CAPITAL (Wang et al. 2021) — restricted to
univariate threshold search.**

Both solve the *same constrained optimization*: among subgroups
with average treatment effect at or above a clinically meaningful
threshold, return the largest. The differences are in the search
space:

* CAPITAL: depth-$d$ tree over all candidate covariates.
  Multi-covariate splits allowed; subgroup definitions can be
  conjunctions like "$x_3 > c_1$ AND $x_5 \le c_2$".
* `dina_subgroup`: univariate threshold on a single covariate.
  Subgroup definitions are simple: "$x_j$ direction $q$".

The univariate restriction is a deliberate simplification, not a
limitation imposed by DINA. It buys interpretability (single-
biomarker subgroups, easier for clinical reporting) and inference
tractability (the bootstrap procedure naturally targets the linear
functional $a^*{}^\top \hat\beta$, and the stability CIs on
$n_{\text{subgroup}}$ and threshold are well-defined). CAPITAL's
tree search is more flexible but harder to characterize
inferentially.

Conceptually `dina_subgroup` is **closer to CAPITAL than to
Virtual Twins**: VT does a two-step procedure (CATE → tree) where
the tree-stage objective is variance reduction of $\hat\tau(x)$,
not satisfying a threshold constraint. CAPITAL and `dina_subgroup`
share the *constrained-search* framing.

### 8.3 forestsearch

**Closest comparator: SIDES with a bias-corrected bootstrap
inference layer drawing on Harrell et al. 1996 and the
Wager–Hastie–Efron 2014 infinitesimal jackknife.**

Shared structure with SIDES:
* Enumerate-and-screen rather than recursive partitioning.
* Test each candidate subgroup for differential treatment effect
  above a threshold.
* Multiplicity-aware selection.

Differences:
* forestsearch searches *exhaustively* over $\binom{L}{2} + L$
  two-factor combinations (with optional LASSO/GRF pre-selection
  to keep $L$ tractable). SIDES uses recursive partitioning,
  exploring a much smaller candidate set per pass.
* forestsearch's *splitting consistency* criterion is novel —
  random 50/50 within-subgroup splits, with both halves required
  to exceed a consistency threshold (HR $\ge 1$). This isn't the
  same as the interaction $p$-values used in SIDES.
* The two-bias-term bootstrap correction $\eta_b(\hat H^*_b) +
  \eta_b(\hat H)$ is more specific than Harrell's single-term
  correction and is custom to the FS selection-and-estimate
  structure.

What's distinctive in forestsearch is *the combination*:
exhaustive enumeration + splitting consistency + two-bias-term
bootstrap + IJ variance. No single existing method combines all
four.

### 8.4 What's distinctive about the package as a whole

The intellectual contribution of the package isn't the individual
pieces — DINA inherits from R-learner family, `dina_subgroup`
solves CAPITAL's problem with a simpler search, forestsearch
inherits from SIDES — but **the combination**:

1. A *parametric, GLM-family-aware* CATE estimator with built-in
   cross-fitting and sandwich/IJ variance (DINA), spanning four
   outcome families in a unified API.
2. An interpretable *univariate threshold search* against a
   clinically meaningful harm threshold, with multiple
   inferential outputs (Wald CI from `dina_subgroup`,
   selection-adjusted bootstrap percentile CI from
   `dina_subgroup_bootstrap`, plus stability CIs).
3. A separate *exhaustive combinatorial* discovery tool
   (forestsearch) with its own bias-corrected bootstrap inference,
   for users who want broader subgroups than single-biomarker
   thresholds.
4. The two discovery tools designed around the same
   threshold/consistency semantics (HR $\ge 1.25$ screening,
   $\ge 1.0$ consistency, $\ge 60$ patients minimum) so they
   can be applied side-by-side on the same data with comparable
   outputs.

The closest competitor to *the package* (rather than to any
individual function) is the `personalized` R package (Huling & Yu
2021), which also offers a unified framework. But `personalized`
spans different sub-methods (OWL, A-learning, etc.); this package
spans different *outcome families* within a coherent semi-parametric
framework, plus the dual discovery tools.

### 8.5 What would strengthen the positioning

Three concrete next moves if you want to push the comparative
case:

1. **Head-to-head with `best_linear_projection`** on the GBSG data,
   same covariate set. Establishes whether DINA's parametric form
   matches the non-parametric BLP target.
2. **Head-to-head with `policy_tree`** (depth = 1, which restricts
   the policy to a single-covariate threshold — the closest like-
   for-like comparison to `dina_subgroup`).
3. **Head-to-head with CAPITAL** on a simulation where the true
   harm subgroup is single-covariate vs. multi-covariate, to map
   where each method dominates.

These could be deferred until the methodology Quarto is closer to
publication-ready; they're benchmarks rather than core methodology.
