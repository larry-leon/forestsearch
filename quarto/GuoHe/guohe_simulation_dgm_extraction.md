---
title: "Guo & He (2021): Simulation Design Extraction"
subtitle: "Verbatim reading of Section 5, for confirmation before implementation"
bibliography: []
---

# Guo & He (2021) — Simulation Design Extraction

**Source.** Guo X, He X. Inference on Selected Subgroups in Clinical Trials.
*Journal of the American Statistical Association* 2021;116(535):1498–1506.
doi:10.1080/01621459.2020.1740096. Extracted from `GuoHe_2021.pdf` (text layer,
`pdftotext -layout`) and `GuoHe-Supplemental_2021.pdf`; implementation details
cross-read from the authors' `synthetic code.R`.

**Status.** This document records what the paper *says*. Nothing here has been
implemented. It exists to be confirmed or corrected before a single line of
simulation code is written, per the standing instruction that a mis-specified
DGM is the most likely failure mode.

---

## 1. Correction to the session brief: the simulations are Section 5

The brief directs reading of "Section 4". Section 4 is the **MONET1
synthetic-data application** (Tables 1–2). The **Monte Carlo simulation study is
Section 5**, comprising:

| Sub-section | Setting | Tables |
|---|---|---|
| §5.1 | Two predefined subgroups; then a $k$-sweep | 3, 4, 5, **6** |
| §5.2 | Post-hoc identified subgroups | 7 |
| §5.3 | Synthetic (MONET1) data-generating model | 8 |

Section 5 opens by stating the scope: censored outcomes throughout, treatment
effect measured by the log hazard ratio from the proportional hazards model,
with §5.1 and §5.2 covering predefined and post-hoc identified subgroups
respectively, and §5.3 comparing coverage under the synthetic data-generating
model of Section 4.

**The outcome type is survival (Cox) in every published simulation.** The
closed-form GLM fast path does not apply; the lean `coxph.fit` path is the
relevant one.

---

## 2. Correction to the target list: Tables 3, 4 and 5 are one simulation

Tables 3, 4 and 5 are three summaries of a **single** design — the two
predefined subgroups of §5.1 — reporting respectively:

- **Table 3** — empirical coverage of the 95% lower confidence bound of $\beta_s$
- **Table 4** — average distance between the 95% lower bound and $\beta_s$
- **Table 5** — empirical bias for $\beta_s$

They are not three targets requiring three runs. One replication loop, storing
the lower bound and the bias-reduced estimate per replication, yields all three.
This materially reduces the cost of the brief's targets 1 and 2, which are in
fact the same job.

---

## 3. §5.1 design — two predefined subgroups (Tables 3, 4, 5)

Stated in the text:

| Element | Specification |
|---|---|
| Sample size | $n = 400$ |
| Model | $\lambda(t) = \lambda_0(t)\,e^{\beta_i D}$ for subgroup $i = 1, 2$ |
| Baseline hazard | $\lambda_0(t)$ of Weibull$(1,1)$ |
| Subgroup membership | Each subject falls in one of the two subgroups w.p. $0.5$ |
| Treatment | $D$ random with equal probability |
| Censoring | Random right censoring by $C$ with $\log(C) \sim \mathrm{Unif}(-1.25,\,1.00)$ |
| Censoring rate | About 40% across the $\beta_i$ considered |
| Effect configuration | $\beta_1 = 0$ fixed; $\beta_2$ varied in $[0, 0.5]$ |
| Reported $\beta_2$ grid | $0,\ 1/10,\ 2/10,\ 3/10,\ 4/10,\ 5/10$ |
| Monte Carlo replications | 2000 |
| Tuning grid | $r = 1/3,\ 1/12,\ 1/21,\ 1/30$ |
| Adaptive | Data-dependent choice of $r$ with $v = 5$ folds |
| Comparator | "Naive" — select the better subgroup from $\hat\beta_i$, then proceed as if selection were independent of the data |
| Reported MC standard error | Around 0.005 for all Table 3 entries |

**Bootstrap replicate count $B$ is not stated in Section 5.** The authors' own
MONET1 script uses $B = 2000$ with $r = 0.03$, $n = 1090$, level $= 0.05$. This
is the only published value and is the natural default, but it is an inference,
not a quotation, and is flagged as such.

### Published values to be reproduced

**Table 3** — empirical coverage, 95% lower bound of $\beta_s$:

| $\beta_2$ | $r=1/3$ | $1/12$ | $1/21$ | $1/30$ | Naive | Adaptive |
|---|---|---|---|---|---|---|
| 0 | 0.933 | 0.950 | 0.952 | 0.952 | 0.896 | 0.943 |
| 1/10 | 0.926 | 0.945 | 0.947 | 0.947 | 0.912 | 0.936 |
| 2/10 | 0.928 | 0.949 | 0.951 | 0.951 | 0.910 | 0.939 |
| 3/10 | 0.941 | 0.957 | 0.959 | 0.959 | 0.919 | 0.947 |
| 4/10 | 0.939 | 0.955 | 0.956 | 0.957 | 0.927 | 0.945 |
| 5/10 | 0.952 | 0.965 | 0.965 | 0.966 | 0.934 | 0.953 |

**Table 4** — average distance between the 95% lower bound and $\beta_s$:

| $\beta_2$ | $r=1/3$ | $1/12$ | $1/21$ | $1/30$ | Naive | Adaptive |
|---|---|---|---|---|---|---|
| 0 | 0.248 | 0.265 | 0.266 | 0.266 | 0.213 | 0.258 |
| 1/10 | 0.252 | 0.269 | 0.270 | 0.270 | 0.218 | 0.262 |
| 2/10 | 0.267 | 0.285 | 0.288 | 0.287 | 0.233 | 0.277 |
| 3/10 | 0.290 | 0.311 | 0.313 | 0.314 | 0.258 | 0.302 |
| 4/10 | 0.301 | 0.326 | 0.328 | 0.329 | 0.273 | 0.313 |
| 5/10 | 0.310 | 0.339 | 0.342 | 0.343 | 0.286 | 0.323 |

**Table 5** — empirical bias for $\beta_s$:

| $\beta_2$ | $r=1/3$ | $1/12$ | $1/21$ | $1/30$ | Naive | Adaptive |
|---|---|---|---|---|---|---|
| 0 | 0.028 | 0.008 | 0.007 | 0.006 | 0.107 | 0.018 |
| 1/10 | 0.024 | 0.002 | 0.000 | −0.001 | 0.100 | 0.012 |
| 2/10 | 0.005 | −0.021 | −0.022 | −0.023 | 0.077 | −0.008 |
| 3/10 | −0.003 | −0.045 | −0.036 | −0.037 | 0.061 | −0.018 |
| 4/10 | −0.018 | −0.063 | −0.065 | −0.066 | 0.029 | −0.042 |
| 5/10 | −0.027 | −0.067 | −0.070 | −0.071 | 0.022 | −0.040 |

---

## 4. Table 6 — the $k$-sweep, absent from the session brief

This table was not among the brief's targets and is arguably the most
consequential one for `fs-glms-interpretable`. Design, continuing §5.1:

| Element | Specification |
|---|---|
| Model | $\lambda(t) = \lambda_0(t)\,e^{\beta_i D}$, $i = 1, \dots, k$ |
| Subgroup membership | Equal probability of falling into each of the $k$ subgroups |
| Sample size | $n = 200k$ (grows with $k$) |
| Treatment, censoring | As in the two-subgroup setting |
| Effect configuration | The most challenging scenario, $\beta_1 = \cdots = \beta_k = 0$ |
| $k$ grid | 2, 6, 10, 12 |
| Replications | 2000 |

**Table 6** — coverage of the 95% lower bound and empirical bias for $\beta_s$:

| $k$ | | $r=1/3$ | $1/12$ | $1/21$ | $1/30$ | Naive | Adaptive |
|---|---|---|---|---|---|---|---|
| 2 | Cover | 0.929 | 0.952 | 0.953 | 0.953 | 0.900 | 0.939 |
| | Bias | 0.029 | 0.006 | 0.004 | 0.004 | 0.105 | 0.014 |
| 6 | Cover | 0.891 | 0.941 | 0.943 | 0.945 | 0.739 | 0.930 |
| | Bias | 0.060 | 0.011 | 0.009 | 0.008 | 0.240 | 0.029 |
| 10 | Cover | 0.866 | 0.944 | 0.949 | 0.950 | 0.594 | 0.927 |
| | Bias | 0.066 | 0.009 | 0.006 | 0.005 | 0.290 | 0.031 |
| 12 | Cover | 0.860 | 0.946 | 0.950 | 0.950 | 0.543 | 0.925 |
| | Bias | 0.062 | 0.003 | 0.001 | 0.001 | 0.302 | 0.026 |

### Why this matters for the manuscript

The authors' own narrative around Table 6 states that results are somewhat more
sensitive to the choice of $r$ as the number of subgroups increases, that
smaller values of $r$ generally work better, and that the adaptive rule produced
more under-coverage at higher $k$ — closing with the observation that additional
research is needed to find a more reliable adaptive method across a wide range
of problem settings.

That is a **published, authors' own** statement of the direction the GBSG
$r$-sweep found at $k = 1744$. Reproducing Table 6 therefore does more than
validate the implementation: it establishes that the tuning sensitivity is a
documented property of the method that intensifies with $k$, and that our
large-family finding extrapolates a trend the authors themselves identified
rather than contradicting them. This is materially fairer to Guo & He than the
GBSG result standing alone, and correspondingly harder to dispute.

Note the confound to state plainly when reporting: in Table 6, $n = 200k$ grows
with $k$, so the degradation is attributable to $k$ despite a *growing* sample
size — which strengthens rather than weakens the reading.

---

## 5. §5.2 and §5.3 — secondary targets

**§5.2, post-hoc identified case (Table 7).** $n = 400$; hazard
$\lambda_0(t)e^{b(W)D}$ with $\lambda_0$ Weibull$(1,1)$; $D$ and $W$
independent, $D \sim \mathrm{Bernoulli}(1, 0.5)$, $W \sim \mathrm{Unif}[0,80]$;
censoring as in §5.1, again about 40%. Candidate subgroups $S(c) = \{W \le c\}$
for $c \in [30,60]$, with $\beta(c)$ the subgroup effect of $S(c)$. The effect
function is the step $b(w) = \beta_1$ for $w > 30$ and $\beta_2$ for
$w \le 30$, with $\beta_1 = 0$ fixed and $\beta_2$ varied in $[0,0.5]$;
homogeneity — and hence maximal selection bias — obtains at $\beta_2 = \beta_1$.
2000 Monte Carlo samples. Table 7 reports **no Adaptive column**, only
$r = 1/3, 1/12, 1/21, 1/30$ and Naive.

The paper notes that given $S(c)$ the event time need not follow a proportional
hazards model, but that $\beta(c)$ remains well defined by Lin and Wei (1989),
and that $\beta(c)$ is generally *not* $b(c)$ but rather a weighted average of
$b(\cdot)$ over $[0,c]$. **This makes the truth target non-trivial to compute**
and is the main implementation risk in §5.2.

**§5.3, synthetic data model (Table 8).** Coverage of the 95% **upper** bound of
the log hazard ratio of the best selected subgroup, at $r = 1/3, 1/9, 1/30$ plus
Naive and Adaptive: 0.917, 0.946, 0.950, 0.805, 0.935. Note the $r$ grid differs
from Tables 3–6. The supplement carries the design detail (Tables B.1–B.2:
subgroup proportions $p_i$, treatment proportion $p$, and the event- and
censoring-time distributions, all estimated from Kubota et al. 2014).

---

## 6. Orientation and estimand — the point most likely to be got wrong

From §2.1: without loss of generality a **larger** $\beta_i$ means a **better**
treatment effect. Two quantities are distinguished:

- $\beta_s$, the **best selected subgroup effect**, $s = \arg\max_{i \in [k]} \hat\beta_i$ — the true effect of the selected subgroup, a random target;
- $\beta_{\max} = \max_{i \in [k]} \beta_i$, the **best subgroup effect**, a fixed parameter.

Tables 3–6 target $\beta_s$. Corollary 1 is what licenses using the
$\beta_{\max}$ procedure for $\beta_s$: it establishes that the bootstrap
distribution of $\sqrt{n}(\beta^*_{\max,\text{modified}} - \hat\beta_{\max})$
and the sampling distribution of $\sqrt{n}(\hat\beta_{\max} - \beta_s)$ have
uniformly vanishing distance, so the same interval serves as an asymptotically
sharp prediction interval for $\beta_s$.

**Sign convention.** The DGM writes the hazard as $\lambda_0(t)e^{\beta_i D}$,
under which larger $\beta_i$ is larger hazard and therefore *worse*. The
inferential quantity is the **negated** log hazard ratio, exactly as in the
authors' code, which sets `X_bar[l] = -coxph(...)$coef` with the comment that
the sign is changed to make it consistent with the framework of the algorithm.
Selection is $\arg\max$ of the negated coefficient.

Under this reading, with $\beta_1 = 0$ and $\beta_2 \in [0, 0.5]$, subgroup 2
becomes progressively *worse* and subgroup 1 separates as the unique best. At
$\beta_2 = 0$ the two coincide, $|H| = 2$, and the $O(1/\sqrt n)$ bias of
§2.4 is present; as $\beta_2$ grows, $|H| = 1$ and the bias vanishes. This is
consistent with the naive bias in Table 5 declining monotonically from 0.107 to
0.022, and with naive coverage rising from 0.896 to 0.934.

**This is the single item most in need of explicit confirmation before coding.**

---

## 7. Algorithms, verbatim

**Algorithm 1 — lower confidence limit for $\beta_{\max}$.**

1. For $i = 1,\dots,k$, set $d_i = (1 - n^{r-0.5})(\hat\beta_{\max} - \hat\beta_i)$.
2. For $b = 1, \dots, B$: on bootstrap sample $b$ compute subgroup effects
   $\beta^*_{i,b}$, then $T^*_b = \sqrt{n}\{\max_{i \in [k]}(\beta^*_{i,b} + d_i) - \hat\beta_{\max}\}$.
3. Let $c_\alpha = \mathrm{quantile}(T^*_b, 1-\alpha)$. The level $1-\alpha$
   lower confidence limit is $\hat\beta_{\max} - c_\alpha/\sqrt{n}$.

The bootstrap is the **pair bootstrap**: draw $n$ subjects with replacement.

**Bias-reduced estimator (§2.4).**
$\hat\beta_{\max,\text{reduced}} = \hat\beta_{\max} - E^*[\beta^*_{\max,\text{modified}} - \hat\beta_{\max}]$,
where $\beta^*_{\max,\text{modified}} = \max_i(\beta^*_i + d_i)$.

**Algorithm 2 — cross-validated choice of $r$.** Let
$A = \{r_1,\dots,r_m\} \subset (0, 0.5)$.

1. Randomly partition the data into $v$ approximately equal subsamples.
2. For each $r_l$ and each fold $j$: use fold $j$ as the **reference** data and
   the remainder as **training**.
3. On the training data compute $\hat\beta_{\max,\text{reduced},j}(r_l)$.
4. On the reference data compute, for each $i \in [k]$, the effect estimate
   $\hat\beta_{i,j}$ and its standard error $\hat\sigma_{i,j}$.
5. Set $h_{i,j}(r_l) = (\hat\beta_{\max,\text{reduced},j}(r_l) - \hat\beta_{i,j})^2 - \hat\sigma_{i,j}^2$.
6. Choose $r = \arg\min_{r_l}\{\min_{i \in [k]}[\sum_{j=1}^{v} h_{i,j}(r_l)/v]\}$.

**The grid $A$ is not specified in the paper.** The tabulated
$\{1/3, 1/12, 1/21, 1/30\}$ is the natural candidate for Tables 3–6 and
$\{1/3, 1/9, 1/30\}$ for Table 8, but this is an inference. The adaptive column
is sensitive to this choice, and it is the one piece of `guohe_adaptive_r.R`
with no independent validation.

### Correspondence with the authors' code

The published $d_i$ and the code's correction term agree identically:
$c_n[l] = (n^{1/2} - n^{r})(MM - \bar X_l)/\sqrt{n} = (1 - n^{r-0.5})(\hat\beta_{\max} - \hat\beta_l)$.

The code also reveals the **symmetric candidate construction**: each indicator
column contributes both $\{$column $= 1\}$ and $\{$column $= 0\}$, so $i$
columns yield $2i$ candidates — the `symmetric = TRUE` convention.

---

## 8. Structural check already performed

The stated censoring rate is a free, executable check on the DGM reading, and it
passes. With $T \sim$ Weibull$(1,1) = \mathrm{Exp}(1)$ under $\beta = 0$ and
$C = e^U$, $U \sim \mathrm{Unif}(-1.25, 1.00)$, the implied censoring
probability averaged over $D \sim \mathrm{Bernoulli}(0.5)$ is

| $\beta_2$ | $P(\text{censored})$ |
|---|---|
| 0.00 | 0.410 |
| 0.25 | 0.372 |
| 0.50 | 0.336 |

against the paper's "about 40% across different choices of $\beta_i$
considered". The agreement at $\beta_2 = 0$ is close, and the mild decline with
$\beta_2$ is the expected consequence of a higher hazard in the treated half of
subgroup 2. This is consistent with the reading of Weibull$(1,1)$ as the unit
exponential and of the censoring distribution as stated.

This is a structural check on the DGM only. It says nothing about the estimator
and must not be allowed to stand in for one.

---

## 9. Materials still outstanding

The following were named in the session brief and are **not** present:

- `guohe_algorithm3.R` — the estimator under validation
- `guohe_adaptive_r.R` — Algorithm 2, the component with no independent validation
- `maxeff_tests.R` — the fail-loud style template

The uploaded archive supplied the papers, `synthetic code.R` and
`syntheticdata.csv`; the package archive supplied `R/` (90 files), in which Guo
& He appears only as a documentation reference to the `sg_focus = "maxeff"`
argmax primitive. No simulation can be wired around an estimator that is not in
hand.
