# REPORT — signed orientation for `fs_oc_family_enumerate()` (2026-08-31)

**Task:** `dev/tasks/cc_task_oc_signed_orientation_2026-08-31.md` (commit `bd233fd7`)
**Premise verdict:** stage 0's report (`REPORT_oc_applied_stage0_2026-08-31.md`,
commit `209a8e85`) §3.5: **the sign gate BINDS** at every rung of the applied
ACTG175 ladder (`m_tau[Qc] = −26.978725` at all q while every planted
`m_tau[Q] > 0`). The premise is met; the task ran.

*(Supersedes the placeholder `REPORT_oc_signed_orientation_2026-08-31.md`
committed at `6c071c74`, written when stage 0 had stopped at its provenance
gate and no verdict existed.)*

## §1 Provenance

- `pop-os` · branch `feature/glm-extension` · start HEAD `209a8e85` · tree
  clean · installed version at start **0.3.1** ✓ · stage-0 report commit in
  the log ✓.
- Task document already committed alone at `bd233fd7`
  (`~/Downloads` stem hyphen-stripped: `cc_task_oc_signed_orientation_20260831.md`
  → `dev/tasks/cc_task_oc_signed_orientation_2026-08-31.md`); per Larry's
  instruction the transport step was not repeated.
- **Parity baseline** (`devtools::test()` against pre-edit source):
  `[ FAIL 0 | WARN 31 | SKIP 3 | PASS 4844 ]` — WARN 31 as expected.

## §2 Source verification — GATE: PASS

From the pre-edit `R/fs_oc_family.R` (at `209a8e85`):

The same-sign `stop()` (:248–251):

```r
    if (length(unique(sign(reg$m_tau))) != 1L) {
      stop("the region effects m_tau do not share a sign; the oriented ",
           "(abs) reading of the scale table is not defined.", call. = FALSE)
    }
```

The `tauQc` / `bint` / `seQ1000` block (:252–257):

```r
    piQ     <- reg$P_g[iQ]
    m_Q     <- abs(reg$m_tau[iQ])
    tauQc   <- abs(reg$m_tau[iQc])
    bint    <- m_Q - tauQc
    seQ1000 <- sqrt(reg$V_eff[iQ] / (1000 * piQ))
    PQ      <- mean(inQ)
```

The `beta_g` construction at the interface-fields step (:401–406):

```r
  if (!is_null) {
    PQg    <- PgQ / Pg
    beta_g <- tauQc + bint * PQg
    se_g   <- seQ1000 * sqrt(1000 / n) * sqrt(piQ / Pg)
    sens_g <- PgQ / PQ
    spec_g <- 1 - PgQc / (1 - PQ)
```

**Downstream grep** (`grep -n "abs(\|positiv" R/fs_oc_predict.R R/fs_oc_grid.R`):
no `abs()` anywhere in either file; the only "positive" hits are argument
validations on `n`, `draws`, `block` (fs_oc_grid.R:104,115,118;
fs_oc_predict.R:184,197) — none on `beta_g`, `EbetaH`, or the gate. The
consumption is sign-agnostic throughout: draws are
`Bhat ~ N(beta_g, Sg)` (fs_oc_predict.R:251–252 resample; :257–261 split),
the gate compares `Bhat >= c1` and `Bhat - c2 >= z_p * se_g` directly
(:272–273, :276), and the functionals are plain weighted sums —
`EbetaH = sum(sel_c * beta_g)` (:305), `mass_below = sum(sel_c[beta_g < c1])`
(:308). `fs_oc_grid()` runs the same factored blocks (fs_oc_predict.R:238–242).

## §3 The edit

`R/fs_oc_family.R` alternative branch only: the same-sign `stop()` deleted;
`s <- sign(m_tau[Q])` with a `stop()` at `s == 0` directing to a nonzero
planted effect or the null path; `tauQc <- s * m_Qc`,
`bint <- s * (m_Q − m_Qc)` on the **signed** `m_tau` values, so
`beta_g = tauQc + bint * PQg` is the signed mixture
`s * (m_Qc + (m_Q − m_Qc) * PQg)` (bit-identical form: negation and
`abs()` are exact in IEEE arithmetic, so on same-sign families every
intermediate coincides with the former reading); `seQ1000` unchanged;
provenance stored in a new alt-only `orientation` element
(`s`, `m_tau_Q`, `m_tau_Qc`, `tauQc`, `bint`); roxygen updated (orientation
is the harm direction; opposite-sign families supported; `tauQc` may be
negative; benefit-direction candidates carry oriented-negative `beta_g`).
The null branch, `se_g`'s construction, and all downstream files untouched.

The full diff:

```diff
diff --git a/R/fs_oc_family.R b/R/fs_oc_family.R
index 911c4bff..9f8e563c 100644
--- a/R/fs_oc_family.R
+++ b/R/fs_oc_family.R
@@ -80,11 +80,21 @@
 #' records the branch taken.
 #'
 #' \strong{Scale.} Every mean and standard error is derived from
-#' \code{\link{fs_dgm_scale}(dgm)}: with \code{Q} the true harm region,
-#' \code{tauQc = |m_tau[Qc]|}, \code{bint = |m_tau[Q]| - |m_tau[Qc]|},
+#' \code{\link{fs_dgm_scale}(dgm)}.  The orientation is the harm direction
+#' \code{s = sign(m_tau[Q])}, so the planted effect is oriented positive:
+#' with \code{Q} the true harm region, \code{tauQc = s * m_tau[Qc]},
+#' \code{bint = s * (m_tau[Q] - m_tau[Qc])},
 #' \code{seQ1000 = sqrt(V_eff[Q] / (1000 * P(Q)))}, and for each candidate
-#' \code{beta_g = tauQc + bint * PQg} and
-#' \code{se_g = seQ1000 * sqrt(1000 / n) * sqrt(P(Q) / Pg)}.
+#' \code{beta_g = tauQc + bint * PQg} -- the signed mixture
+#' \code{s * (m_tau[Qc] + (m_tau[Q] - m_tau[Qc]) * PQg)} -- and
+#' \code{se_g = seQ1000 * sqrt(1000 / n) * sqrt(P(Q) / Pg)} (sign-free, from
+#' \code{V_eff[Q]}).  Opposite-sign families (\code{sign(m_tau[Qc]) != s})
+#' are supported: benefit-direction candidates carry oriented-negative
+#' \code{beta_g}, and \code{tauQc} may be negative for such families.  When
+#' both region effects share a sign the values coincide exactly with the
+#' former oriented-absolute reading.  A DGM with \code{m_tau[Q]} exactly zero
+#' is rejected (no harm direction to orient by): plant a nonzero Q effect or
+#' use the null path.
 #'
 #' @param dgm An object of class \code{"glm_dgm"} from
 #'   \code{\link{generate_glm_dgm}}, with the true-region indicator in
@@ -106,7 +116,9 @@
 #'     \item{\code{lab}}{Character, length M: the rule of each candidate.}
 #'     \item{\code{Pg}}{Population prevalence \code{P(g)}.}
 #'     \item{\code{PQg}}{Purity \code{P(g & Q) / P(g)}.}
-#'     \item{\code{beta_g}}{Mixture mean \code{tauQc + bint * PQg}.}
+#'     \item{\code{beta_g}}{Oriented mixture mean
+#'       \code{s * (m_tau[Qc] + (m_tau[Q] - m_tau[Qc]) * PQg)}, with
+#'       \code{s = sign(m_tau[Q])} the harm direction.}
 #'     \item{\code{se_g}}{Anchored standard error at \code{n}.}
 #'     \item{\code{sens_g}}{\code{P(g & Q) / P(Q)}.}
 #'     \item{\code{spec_g}}{\code{1 - P(g & Qc) / P(Qc)}.}
@@ -117,6 +129,12 @@
 #'       population membership of each candidate.}
 #'     \item{\code{null}}{Logical: \code{TRUE} when the null branch was
 #'       taken (Q empty).}
+#'     \item{\code{orientation}}{Alternative branch only: list with the
+#'       harm-direction sign \code{s}, the signed region effects
+#'       \code{m_tau_Q} and \code{m_tau_Qc}, and the oriented mixture
+#'       coefficients \code{tauQc = s * m_tau_Qc} (may be negative for
+#'       opposite-sign families) and \code{bint = s * (m_tau_Q - m_tau_Qc)}.
+#'       Absent on the null branch.}
 #'     \item{\code{scale}}{The \code{fs_dgm_scale} object used.}
 #'     \item{\code{n}}{The trial size.}
 #'     \item{\code{args_used}}{The \code{forestsearch_args} entries consumed,
@@ -245,14 +263,26 @@ fs_oc_family_enumerate <- function(dgm, forestsearch_args, n,
       stop("fs_dgm_scale(dgm) did not return both 'Q' and 'Qc' regions.",
            call. = FALSE)
     }
-    if (length(unique(sign(reg$m_tau))) != 1L) {
-      stop("the region effects m_tau do not share a sign; the oriented ",
-           "(abs) reading of the scale table is not defined.", call. = FALSE)
+    # Signed orientation: s is the harm direction sign(m_tau[Q]), so the
+    # planted effect is oriented positive and the general signed mixture
+    #   beta_g = s * (m_tau[Qc] + (m_tau[Q] - m_tau[Qc]) * PQg)
+    # applies.  When both region effects share the sign this is algebraically
+    # (and in IEEE arithmetic bit-) identical to the previous |.| reading
+    # (s*m_Qc = |m_Qc|, s*(m_Q - m_Qc) = |m_Q| - |m_Qc|); when they differ,
+    # benefit-direction candidates carry oriented-negative means, so tauQc
+    # may be negative.
+    m_Q_signed  <- reg$m_tau[iQ]
+    m_Qc_signed <- reg$m_tau[iQc]
+    s <- sign(m_Q_signed)
+    if (s == 0) {
+      stop("m_tau[Q] is exactly zero: there is no harm direction to orient ",
+           "by.  Plant a nonzero Q effect, or use the null path ",
+           "(generate_glm_dgm(model = \"null\")) for a homogeneous-effect ",
+           "family.", call. = FALSE)
     }
     piQ     <- reg$P_g[iQ]
-    m_Q     <- abs(reg$m_tau[iQ])
-    tauQc   <- abs(reg$m_tau[iQc])
-    bint    <- m_Q - tauQc
+    tauQc   <- s * m_Qc_signed
+    bint    <- s * (m_Q_signed - m_Qc_signed)
     seQ1000 <- sqrt(reg$V_eff[iQ] / (1000 * piQ))
     PQ      <- mean(inQ)
   } else {
@@ -429,6 +459,14 @@ fs_oc_family_enumerate <- function(dgm, forestsearch_args, n,
     cuts = cuts,
     counts = counts
   )
+  if (!is_null) {
+    # Orientation provenance (alternative branch only; the null branch is
+    # unchanged in structure): the harm-direction sign, the signed region
+    # effects, and the oriented mixture coefficients actually used.
+    out$orientation <- list(s = s, m_tau_Q = m_Q_signed,
+                            m_tau_Qc = m_Qc_signed,
+                            tauQc = tauQc, bint = bint)
+  }
   class(out) <- c("fs_oc_family", "list")
   out
 }
```

Note the deleted check also covered `m_tau[S]` (`unique(sign(reg$m_tau))`
runs over all three rows); `S` is a convex mixture of `Q` and `Qc` and never
enters the mixture formula, so on the same-sign domain the check's reach was
identical, and on the extended domain it has no role.

## §4 The identity gates — all PASS

### 4.1 MD40, bit-for-bit — GATE: PASS

Rebuilt exactly as the breadth ladder's §2 did: the driver frame (z1–z12,
switched treat coding, `str2` factors) and the direct
`generate_glm_dgm(k_inter = k40)` route with
`k40 = −13.744764123964362` from `oc_breadth_ladder_2026-08-30_gate.rds`;
enumerated at n = 500 with `CORR$fs_args`
(`oc_wrapper_grid_corrected_2026-08-30.rds`). Against the corrected run's
stored family (`CORR$alt$families[["500"]]`, M = 1696), `identical()` on
every stored field:

```
   lab     Pg    PQg beta_g   se_g sens_g spec_g      M
  TRUE   TRUE   TRUE   TRUE   TRUE   TRUE   TRUE   TRUE
```

Orientation on this family: `s = −1`, `tauQc = 26.2552358760`,
`bint = 13.7447641240` — the former absolute values exactly.

### 4.2 The standing guards — GATE: PASS

Re-run after the edit:

- `fidelity_fs_oc_predict_2026-08-28.R` → `FIDELITY GATE: PASS
  (bit-identical)`, exit 0 (all 13 quantities `identical()`, including the
  nine-quantity vector).
- `prerefactor_reference_2026-08-29.R check` → `REFACTOR GUARD: PASS
  (identical to the 0.2.4 reference)`, exit 0 (`hand_resample`,
  `hand_split`, `fam_resample`, `fam_split` all `identical: TRUE`).

### 4.3 Same-sign synthetic, pre- vs post-edit — GATE: PASS

Small deterministic DGM (seed 41, N = 2000, both region effects negative:
`m_tau[Q] = −66`, `m_tau[Qc] = −26`). **Pre-edit route: the installed 0.3.1
library** (the pre-edit build; source tree at `209a8e85`) — stated per the
task's "git stash or temp lib — say which". Post-edit via
`devtools::load_all()`. Every pre-edit field `identical()`:

```
  lab  Pg  PQg  beta_g  se_g  sens_g  spec_g  ovl  M  PQ  memb  null
  scale  n  args_used  cuts  counts        — all TRUE
fields added post-edit: orientation
```

(The task itself adds the `orientation` provenance element, so the gate is
field-wise `identical()` over the complete pre-edit field set; the single
added field is the documented extension.) For the record, the pre-edit
build's opposite-sign fixture stops with the old message verbatim: `the
region effects m_tau do not share a sign; the oriented (abs) reading of the
scale table is not defined.`

### 4.4 The new domain: the applied opposite-sign DGM at q = 20 — GATE: PASS

Stage 0's anchored ACTG175 DGM (Ĥ = `{age <= 37} & !{cd40 <= 507}`,
`k_inter = 20 − beta_treat`), enumerated at n = 1083, `max_M = 10000`:

- **No stop.** `M = 4508` — equal to stage 0's §4 M (same memberships).
- Orientation: `s = +1`, `m_tau_Q = 20`, `m_tau_Qc = −26.978725`,
  `tauQc = −26.978725` (negative, as the extended domain allows),
  `bint = 46.978725`.
- **`beta_g` range: [−26.978725, +19.852268].**
- **Independent recomputation** for five candidates spanning purity 0–1,
  directly from `df_super` memberships and the signed mixture formula:

```
  m=  13  PQg=0.000000  beta_direct=-26.97872450  beta_fam=-26.97872450  |diff|=0.000e+00  age > 37
  m= 967  PQg=0.008056  beta_direct=-26.60024535  beta_fam=-26.60024535  |diff|=0.000e+00  age > 35 & cd80 > 482
  m= 505  PQg=0.057325  beta_direct=-24.28567660  beta_fam=-24.28567660  |diff|=0.000e+00  age <= 30 & preanti > 672
  m=3497  PQg=0.088786  beta_direct=-22.80764881  beta_fam=-22.80764881  |diff|=0.000e+00  wtkg <= 92 & cd40 > 269
  m=1210  PQg=0.996855  beta_direct= 19.85226816  beta_fam= 19.85226816  |diff|=0.000e+00  age <= 37 & cd40 > 506
```

  max |diff| = 0 (well under 1e-12).
- **Coherence direction:** the candidate nearest Ĥ (member 1210,
  `age <= 37 & cd40 > 506`, Jaccard 0.9969) has the largest `beta_g`
  (19.8523 = the maximum), and all 1056 pure-benefit candidates
  (`PQg = 0`) have `beta_g = −26.978725 < 0`.

## §5 Tests

`tests/testthat/test-fs-oc-signed.R` — four tests, hand-built continuous
DGMs (seed 41, N = 800), no `forestsearch()` runs:

- (a) same-sign fixture: `beta_g` `identical()` to the expected values
  constructed in-test from the abs formula
  (`abs(m_Qc) + (abs(m_Q) − abs(m_Qc)) * PQg`); orientation provenance
  agrees with the abs reading.
- (b) opposite-sign fixture: enumeration succeeds; `beta_g` matches the
  signed mixture computed directly from the scale table (1e-12); monotone
  in `PQg`; `PQg = 0` candidates negative; `tauQc < 0`.
- (c) `m_tau[Q]` exactly zero stops with the new message
  (`"no harm direction to orient by"`).
- (d) the null path untouched: `orientation` absent, field names identical
  to the pre-edit structure, `beta_g` the oriented common effect.

Post-edit: `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 20 ]`.

**v5 §9 (test fails against the defect it guards):** with the edit
stashed (`git stash push -- R/fs_oc_family.R`), the file fails
`[ FAIL 7 | WARN 0 | SKIP 0 | PASS 9 ]` — tests (b)'s assertions fail on
`Error: the region effects m_tau do not share a sign; the oriented (abs)
reading of the scale table is not defined.` — then `git stash pop` restored
the edit.

## §6 Close-out

- `devtools::document()`: `man/fs_oc_family_enumerate.Rd` regenerated;
  `NAMESPACE` unchanged (no new exports).
- **`devtools::test()` — GATE: PASS.**
  Baseline (pre-edit): `[ FAIL 0 | WARN 31 | SKIP 3 | PASS 4844 ]`.
  Post-edit: `[ FAIL 0 | WARN 31 | SKIP 3 | PASS 4864 ]` —
  **FAIL 0; warning parity 31 = 31** (the new file contributes 0 warnings
  and +20 passes; no decrease, no unrelated increase).
- **`R CMD check`** (`devtools::check(document = FALSE, args =
  "--no-manual")`, `RSTUDIO_PANDOC =
  /usr/lib/rstudio/resources/app/bin/quarto/bin/tools/x86_64`) —
  **GATE: PASS**, status verbatim:

```
0 errors ✔ | 0 warnings ✔ | 0 notes ✔
```

  (A first invocation used a wrong `RSTUDIO_PANDOC` directory and failed in
  vignette re-building with "Pandoc is required ... but not available";
  environment fixed, re-run — the failing log is preserved in the session
  scratch, the passing run is the gate.)
- `DESCRIPTION` → **0.3.2**; `NEWS.md` 0.3.2 section: opposite-sign region
  effects supported via `s = sign(m_tau[Q])`; same-sign families
  bit-identical (the MD40 gate); `tauQc` may be negative; exact-zero
  planted effect rejected with guidance.
- `devtools::install(upgrade = FALSE)`; installed version confirmed
  **0.3.2**.
- Commits (explicit paths, no push): see below.

## Commits

- `bd233fd7` — task document (transport, from the first pass).
- `eb136a35` — the edit: `R/fs_oc_family.R`,
  `tests/testthat/test-fs-oc-signed.R`, `man/fs_oc_family_enumerate.Rd`,
  `DESCRIPTION`, `NEWS.md`.
- This report, committed separately by explicit path.

## Ten-line verdict

1. Stage 0 says the sign gate binds everywhere on the applied ladder; the
   premise is met and the edit ran.
2. Downstream verification confirmed `beta_g` is consumed sign-agnostically
   in `fs_oc_predict()`/`fs_oc_grid()` — the edit stayed within its
   authorized surface.
3. The oriented-abs reading was replaced by the signed mixture
   `beta_g = s·(m_Qc + (m_Q − m_Qc)·PQg)`, `s = sign(m_tau[Q])`.
4. On same-sign families the two readings coincide bit-for-bit: the MD40
   family gates `identical()` on every stored field (M = 1696).
5. Both standing guards pass bit-identically after the edit.
6. The synthetic same-sign family is field-wise `identical()` pre- vs
   post-edit; the only addition is the `orientation` provenance element.
7. The applied opposite-sign DGM now enumerates: M = 4508 (equal to stage
   0), `beta_g` ∈ [−26.98, +19.85], independent recomputation exact.
8. Coherence holds: the candidate nearest Ĥ carries the largest `beta_g`;
   pure-benefit candidates are oriented-negative.
9. Tests: 20 new passes, 0 new warnings, FAIL 0 overall, warning parity
   31 = 31; the guarded defect demonstrably fails the reverted code;
   check 0/0/0.
10. 0.3.2 installed; the applied OC evaluation task can now enumerate its
    ladder — its go/no-go remains Larry's, on stage 0's §5 projection.

## Deviations

- The first pass of this task (2026-08-31, HEAD `6c071c74`'s history)
  ended at the premise gate because stage 0 had stopped at its provenance
  gate; Larry's parking instruction unblocked stage 0, whose re-run
  produced the binding verdict this run keys on. The placeholder report
  from that pass is superseded by this file.
- `R CMD check` was invoked twice: the first invocation's `RSTUDIO_PANDOC`
  pointed at a non-existent directory (environment error, not a package
  failure); the gate is the corrected, passing run.
- None otherwise: no change to the null branch, `se_g`, downstream
  consumers, drivers or documents; no renders; no push.
