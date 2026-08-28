# REPORT — OC-wrapper step 0: read-only inventory of the prediction document

Task: `dev/tasks/cc_task_oc_wrapper_inventory_2026-08-28.md` (Parts A, B and C, per Larry's disposition of 2026-08-28).
Session type: read-only inventory. Nothing was simulated, rendered, installed, or edited; the only executed code was the static parser (`qmd_chunk_index.R`), `readRDS()` metadata reads, and no B4 diagnostic (B4 is "not available" — see Part B).

All line numbers in this report are valid at commit `86342569` (the task-doc commit; the inventoried files are identical at `6afb3df3`, the HEAD found at session start — the task-doc commit touched only `dev/tasks/`).

---

## 0. Provenance (raw, §2)

```
$ hostname; pwd
pop-os
/home/larryleon/Documents/GitHub/forestsearch
$ git rev-parse --abbrev-ref HEAD; git rev-parse --short HEAD
feature/glm-extension
6afb3df3
$ git status --porcelain
(empty — clean tree)
$ git log --oneline -6
6afb3df3 MRCT literature search
27dd09e2 docs(reports): record the 0.2.2-vs-0.2.0 application reproduction check
a4ec8c6d test(consistency): pin search reproducibility across worker counts
7c8a2e57 docs(news): record the four 0.2.2 additions missing from NEWS.md
cfd17663 docs(tasks): record the blueprint-gaps investigation task
dca84c93 refactor(analytic): read the symmetric root from the package
$ git merge-base --is-ancestor a4ec8c6d HEAD && echo ANCESTOR_OK || echo ANCESTOR_FAIL
ANCESTOR_OK
$ Rscript -e '...'
forestsearch 0.2.2
R version 4.6.1 (2026-06-24)
```

*GATE passed:* branch `feature/glm-extension`, `ANCESTOR_OK`, clean tree.

Task-doc transport note (§2 second GATE): exactly one file matched in `~/Downloads`; its stem there was `cc_task_oc_wrapper_inventory_20260828.md` (hyphens stripped by transport, as anticipated by the task). It was copied to `dev/tasks/cc_task_oc_wrapper_inventory_2026-08-28.md` and committed alone as `86342569` before any work.

---

## 1. Inputs located (§3)

### 1.1 The prediction document (DOC) — GATE passed

```
$ grep -rl --include='*.qmd' -e 'fs_sym_root' . | sort
quarto/simulations/actg175/continuous/analytic_verification_and_prediction_md_harm.qmd
$ grep -rl --include='*.qmd' -e 'warmup-floors' . | sort
quarto/simulations/actg175/continuous/analytic_verification_and_prediction_md_harm.qmd
```

Exactly one file qualifies on both patterns. `DOC = quarto/simulations/actg175/continuous/analytic_verification_and_prediction_md_harm.qmd`, `STEM = analytic_verification_and_prediction_md_harm`.

- `wc -l`: **2701** lines.
- `git log -1`: `dca84c93 2026-08-27 11:03:10 -0700 refactor(analytic): read the symmetric root from the package`
- `git status --porcelain -- DOC`: clean.
- Seven-label GATE (closed in A1): the index lists each of `anchor`, `warmup-floors`, `worked-predictions`, `worked-sensitivity`, `worked-null`, `appendix-tier2`, `worked-prediction` **exactly once** — passed.

DOC's directory neighbours (path :: YAML title or first comment line; not inventoried):

| file | title / first line |
|---|---|
| `actg175_continuous_simulations.qmd` | "Simulation Studies: Continuous Outcome (ACTG175)" |
| `analytic_derivations_mr_operating_characteristics.qmd` | "Analytic Derivations for Predicting MR Operating Characteristics" (frozen parent) |
| `analytic_verification_math_details.qmd` | "Analytic Verification: Mathematical Details" (frozen parent) |
| `md_harm_mr_simulation_summary.qmd` | "Simulation Summary: MD/Harm MR Coverage Grid (ACTG175 Continuous)" |
| `mr_coverage_sweep_md_harm.qmd` | "MR coverage sweep — continuous / MD, harm direction" |
| `sim_fs_maxeffCons_fb_mr_md40_knoise0_n500_batch_1_100.qmd` | "FS/maxeffCons on the MD/harm DGM — MR vs full-bootstrap, batchable" |
| `sim_fs_maxeffCons_fb_mr_md40_knoise0_n500_batch_1_20.qmd` | (same title) |
| `sim_fs_maxeffCons_fb_mr_md40_knoise0_n500_batch_1_50.qmd` | (same title) |
| `sim_fs_maxeffCons_mr_md40_knoise0_n500_batch_1_1000.qmd` | (same title) — **driver, Part B** |
| `sim_fs_maxeffCons_mr_md40_knoise0_n500_batch_1_100.qmd` | (same title) |
| `sim_fs_maxeffCons_mr_md40_knoise0_n700_batch_1_1000.qmd` | (same title) — **driver, Part B** |
| `betaHhat_truth_md.R` | `# betaHhat_truth_md.R -- SHIM` |
| `o1_maxeffCons_mr_grid.R` | `# O1 -- maxeffCons MR-only grid over n, sibling run of the committed FB/MR twin.` |
| `o1_readout.R` | `# O1 readout -- columns matched to REPORT_stopB_md_harm_grid.md section 1 so the` |
| `o2_fb_mr_batches.R` | `# O2 -- FB/MR batches 51 onward on the anchor cell.` |

### 1.2 The payloads the document reads

From `grep -n 'readRDS\|source(\|file.path' DOC` (index §7 agrees; the document contains no `source()` call):

**(a) Unconditional — the `anchor` chunk (DOC:35–39):**

```r
.pay_path <- file.path("mr_md_harm",
                       "fs_maxeffCons_mr_md40_knoise0_n500_s1000_d5000",
                       "fs_maxeffCons_mr_md40_knoise0_n500_res_1_1000.rds")
stopifnot(file.exists(.pay_path))
.scale <- readRDS(.pay_path)$scale
```

Resolved: `quarto/simulations/actg175/continuous/mr_md_harm/fs_maxeffCons_mr_md40_knoise0_n500_s1000_d5000/fs_maxeffCons_mr_md40_knoise0_n500_res_1_1000.rds` — exists, tracked, `git log -1`: `bb75cca6 2026-08-26 16:03:20 -0700 sims: re-run both production cells with computed scale and corrected oracle bounds`.

Metadata (readRDS, names/str only):

```
names(p): "results" "truth" "scale" "oc" "meta"
names(p$scale): "regions" "sigma" "outcome_type" "effect_measure" "rand_ratio" "p_treat" "n_super"
str(p$scale, max.level=1): regions df 3x16; sigma 128; outcome_type "continuous"; effect_measure "MD";
                           rand_ratio 1; p_treat 0.5; n_super 5000; class c("fs_dgm_scale","list")
str(p$scale$regions): region chr Q/Qc/S; n_g int 1723 3277 5000; P_g 0.345 0.655 1;
  m_mu0, m_mu1, m_tau (-40, -26.3, -31); V_mu0, V_mu1, V_tau, C_mu0_tau;
  v_cond0 16256 (all three); v_cond1; V_arm0; V_arm1; bracket 16863/17009/16973; V_eff 67454/68034/67891
```

**(b) Conditional — the `itt-live` chunk (DOC:749–751):**

```r
pfile <- "mr_sweep_md_harm/md_harm_s50_pilot/fs_md_harm_n1000_res.rds"
if (file.exists(pfile)) {
  b <- readRDS(pfile); tr <- b$truth; r <- b$results
```

Condition: `file.exists(pfile)`. Resolved path is **absent on disk and untracked** (a machine-local summary-era pilot; `git log -1 -- <path>` shows only a historical `ae8c5557 2026-08-10 continuous 1-100 md40`). The chunk skips silently, exactly as its guard note in prose states (DOC:762–764). No metadata read possible.

The document reads no other payload; in particular the **null branch reads nothing** — `worked-null` (DOC:2062–2077) is a pure computation on objects from earlier chunks (`grep -n 'mdnull' DOC` returns nothing).

### 1.3 The REPORT directory — GATE passed

```
$ git log --diff-filter=A --name-only --format= -- '*REPORT_*' | head -3
dev/glm-continuous-sims/REPORT_repro_check_vs_0.2.0_2026-08-27.md
dev/glm-continuous-sims/REPORT_overnight2.html
dev/glm-continuous-sims/REPORT_overnight2.md
```

`RDIR = dev/glm-continuous-sims`. This report is `RDIR/REPORT_oc_wrapper_inventory_2026-08-28.md`.

### 1.4 Drivers and measured payloads (for Part B)

```
$ git show --stat --format='%h %s' dafff647 d884adbf 2b180813
dafff647 fix(sims): guard the null cell's dgm_scale path so the branch can execute
  sim_fs_maxeffCons_mr_md40_knoise0_n500_batch_1_1000.qmd | 19 +++++++++++--------
  sim_fs_maxeffCons_mr_md40_knoise0_n700_batch_1_1000.qmd | 19 +++++++++++--------
d884adbf fix(sims): save-site name vectors tolerate an absent scale element
  (same two files, 6 ++/-- each)
2b180813 Create fs_maxeffCons_mr_mdnull_knoise0_n500_res_1_1000.rds
  quarto/simulations/actg175/continuous/mr_md_harm/fs_maxeffCons_mr_mdnull_knoise0_n500_s1000_d5000/
    fs_maxeffCons_mr_mdnull_knoise0_n500_res_1_1000.rds | Bin 0 -> 178654 bytes
```

Drivers (both call `forestsearch(` at their line 561):
- `quarto/simulations/actg175/continuous/sim_fs_maxeffCons_mr_md40_knoise0_n500_batch_1_1000.qmd`
- `quarto/simulations/actg175/continuous/sim_fs_maxeffCons_mr_md40_knoise0_n700_batch_1_1000.qmd`

Null-cell measurement payload (from `2b180813`):
`quarto/simulations/actg175/continuous/mr_md_harm/fs_maxeffCons_mr_mdnull_knoise0_n500_s1000_d5000/fs_maxeffCons_mr_mdnull_knoise0_n500_res_1_1000.rds` — exists, tracked.

Alternative-cell measurement payloads, identified from the `m500` block and the appendix prose (A3 S14): the `m500` comment reads "1000/1000, tabulated from the payload" and @sec-appendix-n500 names `sim_fs_maxeffCons_mr_md40_knoise0_n500_batch_1_1000`; the n=700 appendix names `sim_fs_maxeffCons_mr_md40_knoise0_n700_batch_1_1000`. Their tracked payloads:
- n=500: the same file the `anchor` chunk reads (1.2a above).
- n=700: `mr_md_harm/fs_maxeffCons_mr_md40_knoise0_n700_s1000_d5000/fs_maxeffCons_mr_md40_knoise0_n700_res_1_1000.rds` — exists, tracked, added in `bb75cca6`.

(One further tracked payload exists and is not read by DOC or the appendix: `mr_md_harm/fs_maxeffCons_mr_md40_knoise0_n500_s100_d5000/fs_maxeffCons_mr_md40_knoise0_n500_res_1_100.rds`.)

---

## 2. Part A — what the prediction document's chunks compute

### A1. Mechanical index

- Script saved verbatim as `RDIR/qmd_chunk_index.R` (extracted from the task doc's Appendix 1; 440 lines).
- Run: `Rscript RDIR/qmd_chunk_index.R DOC RDIR/chunk_index_analytic_verification_and_prediction_md_harm_86342569.md` → `wrote ... 15 chunks; 0 parse errors`.
- Cross-check: `grep -cE '^\s*`{3,}\{r' DOC` → **15**. Index reports **15 fenced chunks, all 15 R**. **The counts agree.** Zero parse errors; no unterminated fences.
- Seven-label GATE: passed (each named label exactly once; see the index).
- Index file: [`chunk_index_analytic_verification_and_prediction_md_harm_86342569.md`](chunk_index_analytic_verification_and_prediction_md_harm_86342569.md) (committed alongside this report; not pasted here).

The 15 chunks, for reference below (label — lines at `86342569`): `setup` 20–23 · `anchor` 25–103 · `itt-live` 748–759 · `clearance-example` 1155–1160 · `floors-single` 1362–1427 · `worked-prediction` 1510–1546 · `warmup-floors` 1571–1613 · `design-calculator` 1632–1638 · `worked-scenario` 1733–1790 · `worked-predictions` 1825–1875 · `worked-sensitivity` 1957–2018 · `worked-null` 2062–2077 · `appendix-identity-locks` 2245–2276 · `appendix-tier2` 2311–2371 · `appendix-n500` 2409–2441.

Every chunk in the index appears in the A2 table at least once (completeness check satisfied; `setup` appears as "display/no quantity").

### A2. Quantity table

Column key — kind: `payload-read` / `analytic` / `MC` (with draw count) / `typed` / `assertion` / `display`. computation: `G` = expressed only in the scale table, `n`/`c1`/`c2`, an abstract family (memberships / prevalences / overlaps), draw counts and seeds; `F` = refers to fixture-specific columns, the Q rule, ACTG175 variables, `seed_base`, non-scale payload elements, or fixture-typed numbers. scale: `MD-only` / `any` / `unknown` (v5 §6 criterion: relies on the MD branch of `fs_dgm_scale()`, the negligibility argument, or collapsibility of the mean difference → MD-only). Evidence: line ranges at `86342569`; verbatim code fenced below the table, keyed by chunk.

| id | quantity | chunk | definition as implemented | inputs | constants / draws | kind | comp. | inputs fixture-bound? | scale | evidence |
|---|---|---|---|---|---|---|---|---|---|---|
| A0 | (no quantity) | `setup` | `knitr::opts_chunk$set(comment = "#>")` | — | — | display | — | — | — | L20–23 |
| A1 | `.scale` | `anchor` | `readRDS(.pay_path)$scale`, asserted `fs_dgm_scale` | payload | — | payload-read | G (the `fs_dgm_scale` interface) | yes: the n=500 payload path is fixture-typed | any | L35–40 |
| A2 | `.reg`, `.iQ`, `.iQc`, `.iS` | `anchor` | `.scale$regions`; `match("Q"/"Qc"/"S", .reg$region)` | `.scale` | — | payload-read | G | yes: region names Q/Qc/S | any | L41–42, 60, 66 |
| A3 | `.anchor$sigma, piQ, bracket, V_eff, m_Q` | `anchor` | direct reads of `.scale$sigma`, `P_g[.iQ]`, `bracket[.iQ]`, `V_eff[.iQ]`, `abs(m_tau[.iQ])` | `.reg` | — | payload-read | G | yes: payload | MD-only (bracket/V_eff are the MD branch) | L45–51 |
| A4 | `.anchor$seQ1000` | `anchor` | `sqrt(V_eff / (1000 * piQ))` | A3 | 1000 | analytic | G | yes: A3 | MD-only | L54 |
| A5 | `.anchor$sd_eff` | `anchor` | `sqrt(bracket)` | A3 | — | analytic | G | yes: A3 | MD-only | L55 |
| A6 | `.anchor$c2_ceiling` | `anchor` | `pnorm((m_Q - 10) / (seQ1000*sqrt(2)))^2` | A3, A4 | 10 (c2 floor), 2 | analytic | G | yes: A3 | MD-only (inputs); form any | L57 |
| A7 | `.anchor$scr_clear_Q` | `anchor` | `pnorm((m_Q - 30) / seQ1000)` | A3, A4 | 30 (c1 floor) | analytic | G | yes: A3 | MD-only (inputs) | L58 |
| A8 | `.anchor$m_Qc, beta_int, V_mu0_Q, m_ITT, v_cond0` | `anchor` | `abs(m_tau[.iQc])`; `m_Q - m_Qc`; `V_mu0[.iQ]`; `abs(m_tau[.iS])`; `v_cond0[.iQ]` | `.reg` | — | payload-read (beta_int derived) | G | yes: payload | MD-only | L62–69 |
| A9 | regions-coherence assertions | `anchor` | `stopifnot(P_g[S]==1, n_g[S]==n_super, P_g[Q]+P_g[Qc]==1, n_g sums, one sign of m_tau, mixture m_ITT, v_cond0 region-invariant)` | `.reg`, `.scale` | — | assertion | G | yes: payload | any | L73–88 |
| A10 | `.fmt_tex`, `.mc_se_prop` | `anchor` | helpers: LaTeX thousands format; `sqrt(p*(1-p)/R)` | — | — | display / analytic helper | G | no | any | L92–101 |
| A11 | `pred_itt`, `obs_itt` + lock | `itt-live` | `tr$effect_Qc + tr$beta_inter*tr$prevalence_Q`; `unique(round(r$betaHhat_Hc[undetected],9))`; `stopifnot(len==1, |diff|<1e-8)` | pilot payload (`truth`, `results`) | 1e-8 | payload-read + assertion | F: reads `betaHhat_Hc`, `detected` columns of a specific bundle | yes: pilot payload | MD-only (mixture identity) | L748–759; chunk skips (file absent) |
| A12 | `z`, `pred` (marginal clearance of Q at n=1000) | `clearance-example` | `z <- (40-30)/.anchor$seQ1000; pred <- pnorm(z)` | `.anchor` | 40 (oriented truth, typed), 30 | analytic | G (Φ of standardized margin) | yes: 40 is the fixture truth typed; seQ1000 | MD-only (inputs) | L1156–1159 |
| A13 | `detect1(c1,c2,m,s)` | `floors-single` | exact single-candidate split-pair rate: `integrate(dnorm(w,m,s√2)*pnorm(pmax(c2, 2c1−w), m, s√2, lower=F), c2, Inf)` | args | rel.tol 1e-10 | analytic (helper feeding many rows) | G | no | any | L1366–1371 |
| A14 | `c1_for(rho,c2,m,s)` | `floors-single` | inversion: `uniroot(detect1−rho, c(m−8s, m+8s), tol=1e-8)`; `NA` if `rho >= pnorm(c2,m,s√2,lower=F)^2` (c2 ceiling) | A13 | 8, 1e-8 | analytic | G | no | any | L1372–1377 |
| A15 | `chan` (3 ideal channels) | `floors-single` | means `(tQ, ITT, tQc)`; prevalences `(piQ, 1, 1−piQ)`; `sigma_n1000 = sQ1*sqrt(piQ/prev)` | `.anchor` | — | analytic | G (prevalence-scaled scale) | yes: anchors | MD-only | L1378–1385 |
| A16 | `res` (channels × n table: `screen_only`, `P1_c2_10`, `P1_c2_25`) | `floors-single` | `pnorm(30, mean, σ·sc, lower=F)`; `detect1(30, {10,25}, mean_i, σ_i·sc)`, `sc = sqrt(1000/n)`, n ∈ {500,700,1000} | A13, A15 | 30, 10, 25, 500/700/1000 | analytic | G | yes: channel means | MD-only (inputs) | L1386–1399 |
| A17 | separation & `nstar(rho,alpha)` | `floors-single` | `bint*(1−piQ)`; `1000*((z_ρ·sQ1 + z_{1−α}·sQ1·√piQ)/(bint(1−piQ)))^2` at (.90,.10),(.80,.20),(.95,.05) | `.anchor` | 1000; the 3 (ρ,α) pairs | analytic | G | yes: anchors | MD-only | L1400–1407 |
| A18 | `tr` (c1 trade-off table) | `floors-single` | `detect1(c1, 10, m_channel, s_channel)` over `c1 = seq(28,42,2)` | A13, A15 | 28–42 by 2; c2=10 | analytic | G | yes | MD-only (inputs) | L1408–1415 |
| A19 | `inv` (exact c1 inversions) + ceiling readout | `floors-single` | `c1_for(c(.80,.90,.95), c2=10, tQ, sQ1)`; ceiling printed `pnorm(10, tQ, sQ1√2, lower=F)^2` | A14 | .80/.90/.95, 10 | analytic | G | yes | MD-only (inputs) | L1416–1418 |
| A20 | `Z` contour of `P1(c1,c2)` | `floors-single` | `outer(g1,g2, detect1(a,b,tQ,sQ1))`, `NA` above diagonal; contour at levels incl. .743 | A13 | grids 22–44/0–40 by .5 | analytic + display | G | yes | MD-only (inputs) | L1419–1426 |
| A21 | `cand` (M=3 warm-up family) | `worked-prediction` | typed: `pi = (piQ, .45, .31)`, `piQg = (1, .28/.45, 1)`; derived `beta = tauQc + bint*piQg`, `se = seQ*sqrt(piQ/pi)` | `.anchor`; typed constants | .45, .31, .28 | typed + analytic | mean/se construction G; the family itself F (typed) | yes: family | MD-only (mixture mean & se) | L1515–1520 |
| A22 | `ov`, `R`, `Sg` (M=3 covariance) | `worked-prediction` | typed overlaps (P∩D via `0.28*0.31/piQ` D⊂Q composition); `R = ov/sqrt(outer(pi,pi))`; `Sg = R*outer(se,se)` | A21 | — | typed + analytic | formula G; entries F | yes: overlaps | any (given se) | L1521–1525 |
| A23 | `p_marg` | `worked-prediction` | `pnorm(c1, beta, se, lower=F)` per member | A21 | c1=30 | analytic | G | yes: family | MD-only (inputs) | L1526 |
| A24 | `p_fam` (family declaration, screening-only) | `worked-prediction` | `1 − pmvnorm(upper=rep(c1,3), mean=beta, sigma=Sg, GenzBretz(abseps=1e-6))` | A21, A22 | 1e-6; seed 20260824 | analytic (quasi-MC quadrature, seeded) | G | yes: family | any (given beta,Sg) | L1511, 1527–1528 |
| A25 | `sel` (selection probabilities, screening-only argmax) | `worked-prediction` | per member m: `pmvnorm(lower=(c1,0,0))` of `A %*% beta`, `A Sg A'` with contrast rows `e_m`, `e_m − e_j` | A21, A22 | — | analytic | G | yes: family | any | L1529–1537 |
| A26 | `out` table + identity gap | `worked-prediction` | kable of beta/se/`p_marg`/`sel`; caption prints `p_fam` and `abs(sum(sel) − p_fam)` | A21–A25 | round 3/4 | display (gap = internal consistency readout) | G | yes | — | L1538–1545 |
| A27 | `P1_3` | `warmup-floors` | `detect1(c1, 10, cand$beta[m], cand$se[m])`, m=1..3 | A13, A21 | c2=10 | analytic | G | yes: family | MD-only (inputs) | L1573–1574 |
| A28 | `rt3`, `rt5` | `warmup-floors` | `forestsearch::fs_sym_root(2*Sg, scale=1)`; `fs_sym_root(2*S5, scale=1)` | A22, A33 | — | analytic | G | yes: Sg | any | L1575, 1597 |
| A29 | `fam_split(mu,root,M,c2v)` | `warmup-floors` | MC split-pair family rate: two draws `V1,V2 = Z %*% t(root) + mu`; `mean(rowSums((V1+V2 ≥ 2c1) & (V1 ≥ c2v) & (V2 ≥ c2v)) > 0)` | root, mu | `Rg = 2e5` | MC helper | G | no | any | L1576–1581 |
| A30 | `fam10`, `fam25` | `warmup-floors` | `fam_split(cand$beta, rt3, 3, {10, 25})` | A21, A28, A29 | Rg=2e5; c2 ∈ {10,25}; seed 20260826 | MC (2e5) | G | yes: family | any (given beta, Sg) | L1572, 1582–1583 |
| A31 | `beta5`, `prev5`, `se5` (5-member family) | `warmup-floors` | augment with S (ITT rule; prev 1, se `sQ1√piQ`) and Qc (prev 1−piQ, se `sQ1√(piQ/(1−piQ))`) | A21, `.anchor` | — | analytic | G construction; F members (m_ITT, m_Qc anchors) | yes | MD-only | L1584–1588 |
| A32 | `ov5` | `warmup-floors` | typed/derived overlap matrix: S row/col = prev5; Qc∩P = .45−.28; Qc∩{Q,D} = 0; diag = prev | A21 | .45, .28 | typed + analytic | F entries | yes | any | L1589–1594 |
| A33 | `R5`, `S5` | `warmup-floors` | `R5 = ov5/sqrt(outer(prev5,prev5))`; `S5 = R5*outer(se5,se5)` | A31, A32 | — | analytic | G | yes | any | L1595–1596 |
| A34 | `fam5_10` | `warmup-floors` | `fam_split(beta5, rt5, 5, 10)` | A29, A31, A28 | Rg=2e5 | MC (2e5) | G | yes: family | any | L1598 |
| A35 | `scr5`, `P1_5` + kable | `warmup-floors` | `pnorm(c1, beta5, se5, lower=F)`; `detect1(c1, 10, beta5[m], se5[m])` | A13, A31 | c1=30, c2=10 | analytic | G | yes | MD-only (inputs) | L1599–1607 |
| A36 | family-declaration readouts + MC SEs | `warmup-floors` | `cat` of fam10/fam25/fam5_10 with `.mc_se_prop(., Rg)` and `p_fam` | A30, A34, A24, A10 | — | display | — | — | — | L1608–1612 |
| A37 | `n_for(target, se_n1000, mu, c1)` | `design-calculator` | `1000*(se_n1000*qnorm(target)/(mu−c1))^2` at targets .80/.90/.95 | `.anchor$seQ1000` | 1000, mu=40, c1=30 | analytic | G | yes: defaults are fixture (mu=40) | MD-only (inputs) | L1633–1637 |
| A38 | `J3` (3×3 joint) + anchor assertions | `worked-scenario` | typed 9-cell table; `stopifnot` sums to 1, row marginals (.455,.045,.500), col marginals (.727,.023,.250) | typed | the 9 cells; 1e-12 | typed + assertion | F (scenario interior, hand-solved) | yes | any | L1741–1746 |
| A39 | `PQ`, `gstar_P`, `jac` | `worked-scenario` | `PQ = J3[3,1]+J3[3,2]` (Q = age>34 & pre≤744.5); `gstar_P = J3[2,1]+J3[3,1]`; `jac = J3[3,1]/(gstar_P+PQ−J3[3,1])`; asserted vs `piQ_fix` (1e-4) and 0.882 (1e-3) | A38, `.anchor` | .882 | analytic + assertion | F (Q rule encoded in cell indices) | yes | any | L1747–1753 |
| A40 | `cells` (18-cell law), `conds` (6 conditions) | `worked-scenario` | `expand.grid(a,p,v)`; `prob = J3 * ifelse(v, pV, 1−pV)`; `inQ = (a==3 & p<=2)`; 6 single conditions on representable cuts + V | A38 | pV = .42 | typed + analytic | F (cut/V structure fixture-styled) | yes | any | L1755–1762 |
| A41 | `fam`, `Pg`, `keep`, `M`, `lab` (candidate family) | `worked-scenario` | 6 singles + 12 cross-variable pairs = 18; `Pg = sum(cells$prob[m])`; **filter `Pg >= 0.12`** (n.min = 60 of 500, lower bound only); M = 16 after | A40 | 0.12; pair index list | analytic | family construction F; the filter form G | yes: family | any | L1763–1773 |
| A42 | `PQg`, `beta_g`, `se_g` | `worked-scenario` | `PQg = P(m ∩ Q)/Pg`; `beta_g = tauQc + bint*PQg` (mixture mean); `se_g = seQ1000*sqrt(2)*sqrt(piQ_fix/Pg)` (n=500 anchor scaling) | A40, A41, `.anchor` | 2 | analytic | G | yes: family + anchors | MD-only (mixture mean; anchor se) | L1774–1776 |
| A43 | `sens_g`, `spec_g` | `worked-scenario` | `P(m∩Q)/PQ`; `1 − P(m∩Qc)/(1−PQ)` | A40, A41 | — | analytic | G | yes: family | any | L1777–1779 |
| A44 | `ovl`, `Rho`, `Sg` (M=16 covariance) | `worked-scenario` | pairwise `P(m_i ∩ m_j)`; `Rho = ovl/sqrt(outer(Pg,Pg))`; `Sg = Rho*outer(se_g,se_g)` | A41, A42 | — | analytic | G | yes: family | any (given se) | L1780–1783 |
| A45 | `tab3` (family table) | `worked-scenario` | kable of prevalence/purity/beta/se/`clears_floor = beta_g >= c1` | A41–A42 | — | display | — | — | — | L1784–1789 |
| A46 | `Rdraw`, `root_half` | `worked-predictions` | `Rdraw = 2e5`; `root_half = forestsearch::fs_sym_root(Sg, scale = 2)` (cov = 2·Sg; rank deficiency handled by the root's zero-clamp) | A44 | 2e5 | analytic | G | yes: Sg | any | L1827, 1833 |
| A47 | `W1`, `W2`, `Bhat` | `worked-predictions` | half-sample draws `Z %*% t(root_half) + beta_g`; `Bhat = (W1+W2)/2` ~ N(beta, Sg) exactly | A42, A46 | seed 20260825 | MC (2e5) | G | yes: beta, Sg | any | L1826, 1834–1838 |
| A48 | `pass` (per-member single-split declaration) | `worked-predictions` | `(W1 + W2 >= 2*c1) & (W1 >= c2) & (W2 >= c2)` | A47 | c1=30, c2=10 | MC (2e5) | G | yes | any | L1839 |
| A49 | `det`, `det_rate` (family declaration / DETECTED share) | `worked-predictions` | `det = rowSums(pass) > 0`; `det_rate = mean(det)` | A48 | — | MC (2e5) | G | yes | any | L1840, 1845 |
| A50 | `winner` (selection: maxeffCons argmax) | `worked-predictions` | `Bmask = Bhat; Bmask[!pass] = -Inf; winner = max.col(Bmask, "first")`, NA when `!det` | A47, A48 | — | MC (2e5) | G | yes | any | L1841–1842 |
| A51 | `P1` (per-member gate rate) | `worked-predictions` | `colMeans(pass)` | A48 | — | MC (2e5) | G | yes | any | L1843 |
| A52 | `p_sel`, `sel_c` | `worked-predictions` | `tabulate(winner[det], M)/Rdraw`; `sel_c = p_sel/det_rate` (selection given detection) | A49, A50 | — | MC (2e5) | G | yes | any | L1844, 1846 |
| A53 | `EnH` = E\|Ĥ\| given detection | `worked-predictions` | `n_rep * sum(sel_c * Pg)` — selection-weighted **population prevalence × n**; conditional on detection by construction (`sel_c` sums to 1); non-declaration excluded by conditioning | A41, A52 | n_rep=500 | MC-derived | G | yes | any | L1847 |
| A54 | `Esens`, `Espec`, `Eppv`, `Enpv` | `worked-predictions` | selection-weighted `sens_g` / `spec_g` / `PQg` / NPV where `npv_m = (pHc − (PQ − PQg·Pg))/pHc` | A43, A42, A52 | — | MC-derived | G | yes | any | L1848–1851 |
| A55 | `EbetaH` | `worked-predictions` | `sum(sel_c * beta_g)` — selection-weighted mixture mean | A42, A52 | — | MC-derived | G | yes | MD-only (mixture identity) | L1852 |
| A56 | `Enaive_bias` | `worked-predictions` | `mean(noise[cbind(which(det), winner[det])])`, `noise = Bhat − beta_g` — mean winner's noise (winner's curse) | A47, A50 | — | MC (2e5) | G | yes | any | L1853–1854 |
| A57 | `mass_below` | `worked-predictions` | `sum(sel_c[beta_g < c1])` — selection mass on rules with true mean below floor | A42, A52 | — | MC-derived | G | yes | any | L1855 |
| A58 | top-6 selection table; predictions table (+MC SE) | `worked-predictions` | kables; captions map each row to its driver-report comparison column | A49–A57, A10 | top 6 | display | — | — | — | L1856–1874 |
| A59 | `.rng_state` save / restore | `worked-sensitivity` | `.rng_state <- .Random.seed` … `.Random.seed <- .rng_state` | RNG state | — | (hygiene) | G | no | any | L1960, 2017 |
| A60 | `build_J(x,y)` | `worked-sensitivity` | closed-form completion of the 3×3 joint from free `(x=J31, y=J22)` under the anchors; `stopifnot` non-negativity + marginals | A38 anchors | .882, marginals, 0.0946, 1e-9/1e-12 | analytic + assertion | F (anchor values) | yes | any | L1961–1971 |
| A61 | `predict_scn(x,y,pV,R)` | `worked-sensitivity` | re-runs the whole S1→S6 machinery per variant (family, floor 0.12, means, cov, `fs_sym_root(scale=2)`, split gate, argmax) → `(detection, EnH, sens, ppv, EbetaH, naive, below)`; common random numbers via `set.seed(20260827)` inside | A60; machinery of A40–A57 | R=2e5; seed 20260827 | MC (2e5 per variant) | G machinery; F interior | yes | MD-only (mixture mean/se inside) | L1972–2003 |
| A62 | `grid`, `sens_out` (7-variant sweep) | `worked-sensitivity` | variants: base; x=.325/.341; y=.0005/.0060; pV=.30/.55 | A61 | the 7 rows | MC | F (sweep values) | yes | — | L2004–2016 |
| A63 | `beta_null`, `W1n/W2n`, `pass_n` | `worked-null` | `beta_null = rep(tauQc, M)`; new draws through the **same `root_half`** (alternative-cell covariance reused); same gate | A42 (tauQc), A46 | seed 20260828; Rdraw=2e5 | MC (2e5) | G | yes: tauQc, root_half | any (given cov) — the *null scale table* question: see S10 | L2063–2069 |
| A64 | null family false declaration + MC SE | `worked-null` | `mean(rowSums(pass_n) > 0)` | A63 | — | MC (2e5) | G | yes | any | L2070–2072 |
| A65 | `p1n` range; implied `L_eff` | `worked-null` | `colMeans(pass_n)`; `L_eff = log(1 − fam)/log(1 − max(p1n))` | A63, A64 | — | MC-derived | G | yes | any | L2073–2076 |
| A66 | `meas` (n=700 measured record) | `appendix-identity-locks` | typed list of report values (sens .123, ppv .402, nH 74, n 700, est/bias per estimator, sd/width/cov) — cited in-code to report §6.3 / §6.7 | typed | 18 literals | typed | F (measured record) | yes | — | L2246–2252 |
| A67 | `backout`, `lock1` (Prop. 1) | `appendix-identity-locks` | `backout = est − bias` (3 rows); `lock1 = tauQc + bint * meas$ppv` | A66, A42 | — | analytic | G formula | yes | MD-only (mixture identity) | L2254–2255 |
| A68 | `lock2` (classification identity) | `appendix-identity-locks` | `meas$sens * (meas$n * piQ_fix) / meas$nH` | A66 | — | analytic | G | yes | any | L2258 |
| A69 | `cov_formula`, `lock3n`, `lock3m` | `appendix-identity-locks` | `pnorm((w/2−b)/s) − pnorm((−w/2−b)/s)` on naive and MR (bias, width, SD) | A66 | 2 | analytic | G | yes | any | L2260–2262 |
| A70 | `fam0…M7=18` (tier-2 family) | `appendix-tier2` | same `cells`/`conds` law; **filter `Pg7 >= 60/700`** (≈.0857) admits M7=18; se `seQ1000*sqrt(1000/700)*sqrt(piQ_fix/Pg7)`; same mixture means, overlap covariance, `fs_sym_root(scale=2)` | A40; `.anchor` | 60, 700, 1000; seed 20260825; R7=2e5 | analytic + MC (2e5) | G machinery; F family | yes | MD-only (means/ses) | L2312–2338 |
| A71 | tier-2 `pred` (9 quantities) | `appendix-tier2` | detection, EnH (×700), sens/spec/ppv/npv, EbH, naive, below — same functionals as A49–A57 at n=700 | A70 | — | MC (2e5) | G | yes | mixed as above | L2338–2354 |
| A72 | `measd` (typed measured column) | `appendix-tier2` | `c(0.999, 74, 0.123, 0.904, 0.402, 0.662, 31.76, 77.09)` + "not tabulated (bundle-computable)" for below-floor mass | typed | 8 literals | typed | F | yes | — | L2355–2358 |
| A73 | `top7` (top-3 predicted selections + MC SE) | `appendix-tier2` | `top7 <- order(-sel)[1:3]`; cat with `.mc_se_prop(sel[top7], R7)` — despite the name, **three** entries are printed | A71, A10 | 3 | display | — | — | — | L2367–2370 |
| A74 | `m500` (n=500 measured record) | `appendix-n500` | typed list: detection 1.000 ("1000/1000, tabulated from the payload"), nH 72, sens .17, spec .87, ppv .41, npv .66, oracle/naive/mr rows | typed | 18 literals | typed | F | yes | — | L2411–2417 |
| A75 | `bH500` (back-outs) + n=500 locks table | `appendix-n500` | `oracle − oracle_bias`, `naive − naive_bias`; locks re-applied: Prop. 1, classification identity `sens*(500*piQ_fix)/nH`, coverage formula | A74, A42 | 500 | analytic | G formulas | yes | MD-only (Prop. 1 rows) | L2418–2430 |
| A76 | Step-4-vs-n500 comparison table | `appendix-n500` | kable of `(det_rate, EnH, Esens, Espec, Eppv, Enpv, EbetaH, Enaive_bias)` (from `worked-predictions`, unchanged) against `m500` | A49–A56, A74 | — | display | — | — | — | L2431–2440 |

**Verbatim evidence (fenced, keyed by chunk; lines at `86342569`).** The `anchor`, `itt-live`, `clearance-example`, `floors-single`, `worked-prediction`, `warmup-floors`, `design-calculator`, `worked-scenario`, `worked-predictions`, `worked-sensitivity`, `worked-null`, `appendix-identity-locks`, `appendix-tier2` and `appendix-n500` definitions quoted in A3 below are the evidence for the table rows; each A3 stage quotes its chunk's defining lines in full. (The chunks total ~430 code lines; every definition row above is covered by a fenced quote in A3.)

### A3. The assembly, stage by stage

#### S1 — the candidate family

The document constructs **three distinct fixed families**, none by enumeration over data:

**(a) The M=3 warm-up family** (`worked-prediction`, L1515–1523) — typed outright:

```r
cand <- data.frame(
  name = c("Q", "P (overlap proxy)", "D (near-dup, subset of Q)"),
  pi   = c(piQ, 0.45, 0.31),
  piQg = c(1.00, 0.28/0.45, 1.00))          # P(Q | g)
cand$beta <- tauQc + bint * cand$piQg
cand$se   <- seQ * sqrt(piQ / cand$pi)
ov <- matrix(c(piQ, 0.28, 0.31,
               0.28, 0.45, 0.28*0.31/piQ,   # |P∩D| via D ⊂ Q composition
               0.31, 0.28*0.31/piQ, 0.31), 3, 3)
```

Q **is a member by construction** here. Extended to 5 members in `warmup-floors` (L1584–1594) by adding the ITT rule S (prevalence 1) and the complement Qc, with a typed overlap matrix.

**(b) The M=16 scenario family** (`worked-scenario`, L1756–1773) — enumerated over an 18-cell constructed law, not over data:

```r
cells <- expand.grid(a = 1:3, p = 1:3, v = 1:0)
cells$prob <- J3[cbind(cells$a, cells$p)] * ifelse(cells$v == 1, pV, 1 - pV)
cells$inQ  <- cells$a == 3 & cells$p <= 2
# Candidate conditions on representable cuts (both sides) + V
conds <- list(`age>33` = cells$a >= 2, `age<=33` = cells$a == 1,
              `pre<=677.9` = cells$p == 1, `pre>677.9` = cells$p >= 2,
              `V=1` = cells$v == 1, `V=0` = cells$v == 0)
fam <- c(lapply(names(conds), function(nm) list(nm = nm, m = conds[[nm]])),
         unlist(lapply(list(c(1,3), c(1,4), c(2,3), c(2,4),
                            c(1,5), c(1,6), c(2,5), c(2,6),
                            c(3,5), c(3,6), c(4,5), c(4,6)),
   function(ix) list(list(nm = paste(names(conds)[ix], collapse = " & "),
                          m = conds[[ix[1]]] & conds[[ix[2]]]))),
   recursive = FALSE))
Pg  <- vapply(fam, function(f) sum(cells$prob[f$m]), numeric(1))
keep <- Pg >= 0.12                                  # n.min = 60 of 500
fam <- fam[keep]; Pg <- Pg[keep]; M <- length(fam)
```

- Covariates: two 3-category axes standing for `age` (cuts ≤33 / (33,34] / >34) and `preanti` (≤677.9 / (677.9,744.5] / >744.5), plus a stylized binary distractor `V` (P(V=1)=0.42, independent).
- Cut-points: **one representable cut per continuous axis** (33 and 677.9, the J-grid positions the driver documents); the true boundaries (34, 744.5) are cell edges but **not cuttable** by any `conds` entry.
- Both directions: yes (`age>33`/`age<=33`, `pre<=677.9`/`pre>677.9`, `V=1`/`V=0`).
- Conjunction depth: singles and **cross-variable pairs only** (the 12 listed index pairs; no same-variable pairs) — 6 + 12 = **18 before filtering**.
- Prevalence filter: `Pg >= 0.12`, evaluated on the **population law** (`cells$prob`, i.e. df_super-style probabilities, not a sample of n); **lower bound only**; **M = 16 after** (count printed in the caption, L1788–1789). No upper bound, no duplicate/near-duplicate/containment filter.
- Q is **not a member** (not representable in `conds`); the best representable rule g\* = `age>33 & pre<=677.9` is a member.
- Per-member quantities (L1774–1783): `Pg`; purity `PQg = P(g∩Q)/Pg`; `beta_g = tauQc + bint*PQg`; `se_g = seQ1000*sqrt(2)*sqrt(piQ_fix/Pg)`; `sens_g = P(g∩Q)/PQ`; `spec_g = 1 − P(g∩Qc)/(1−PQ)`; pairwise overlaps `ovl`.

**(c) The M7=18 tier-2 family** (`appendix-tier2`, L2313–2323): the **same** 18 candidates from the same `cells`/`conds` law, re-filtered at `Pg7 >= 60/700` (all 18 pass).

**The family's interface** — exactly what downstream stages read: the label vector (`lab`), the prevalence vector (`Pg`), the purity vector (`PQg`), the mean vector (`beta_g`), the se vector (`se_g`), the classification vectors (`sens_g`, `spec_g`), the overlap matrix (`ovl` → `Rho` → `Sg`), and the count `M`. Nothing downstream reads the membership indicators again after `worked-scenario` computes these — except `appendix-tier2`, which re-reads `cells`/`conds` to rebuild them at the other floor.

#### S2 — the covariance across candidates

One formula at every site (`worked-prediction` L1524–1525, `warmup-floors` L1595–1596, `worked-scenario` L1780–1783, `predict_scn` L1988–1990, `appendix-tier2` L2333–2335):

```r
Rho <- ovl / sqrt(outer(Pg, Pg))
Sg  <- Rho * outer(se_g, se_g)
```

i.e. correlation = pairwise overlap over the geometric mean of prevalences, scaled by the per-member ses (which carry `v_cond0`/`V_eff`, `sigma`, `n` through `seQ1000`). Dimension: M×M (3, 5, 16, or 18). `fs_sym_root()` enters as the draw root:

- `warmup-floors`: `rt3 <- forestsearch::fs_sym_root(2 * Sg, scale = 1)` (L1575) and `rt5 <- fs_sym_root(2 * S5, scale = 1)` — root of the **half-sample** covariance 2·Sg, feeding `fam_split`'s two half-sample draws.
- `worked-predictions` L1833, `predict_scn` L1991, `appendix-tier2` L2336: `fs_sym_root(Sg, scale = 2)` — same matrix by the other route.

Rank deficiency: handled inside `fs_sym_root()` by the eigenvalue zero-clamp (`R/fs_sym_root.R:65–70`: `d <- sqrt(pmax(scale * eS$values, 0))`), with the in-document comment (L1828–1832) stating the population law is *exactly* rank-deficient (complement pairs) and that degenerate directions get zero variance, enforcing the linear constraints exactly. `worked-prediction`'s analytic path instead hands `Sg` to `pmvnorm` directly.

#### S3 — single-split declaration

Per-member declaration event — the three conditions, identical at every MC site (`worked-predictions` L1839, `fam_split` L1580, `predict_scn` L1994, `worked-null` L2069, `appendix-tier2` L2341):

```r
pass <- (W1 + W2 >= 2 * c1) & (W1 >= c2) & (W2 >= c2)
```

with `W1`, `W2` the two half-sample fitted effects (mean `beta_g`, covariance 2·Sg each, independent), `c1 = 30` (screening floor; the full-sample mean `(W1+W2)/2 ≥ c1`), `c2 = 10` (consistency floor; **each half** ≥ c2 — "single split" is one half-sample pair standing in for the driver's `p̂cons ≥ 0.90` over 400 splits, an idealization the document states with direction at L1806–1823). Family-level event: `rowSums(pass) > 0` (any member declares).

Computed **both** ways:
- **MC** on the joint-normal draws, 2×10⁵ draws, at every gate site (this is what the document reports for family declaration, detection, selection).
- **Analytically**: single-candidate exact rate `detect1()` (integrate, L1366–1371) and screening-only family exceedance `pmvnorm` (L1527–1528) and per-member marginal `pnorm`. The document reports the analytic figures for single-candidate/channel tables and the MC figures for family-with-consistency; `warmup-floors` prints both side by side (L1608–1610).

#### S4 — detection

Selection rule as coded (`worked-predictions` L1841–1842):

```r
Bmask <- Bhat; Bmask[!pass] <- -Inf
winner <- max.col(Bmask, ties.method = "first"); winner[!det] <- NA_integer_
```

— maxeffCons: the argmax of the **full-sample fitted effect** `Bhat` among gate-passers. The event the document calls "detection" (`det`, `det_rate`) is **family declaration** — any member passes — matching the driver's DETECTED status; there is **no Ĥ = Q or overlap-with-Q event coded anywhere in the document**. Recovery of Q is addressed only distributionally, through the selection distribution `p_sel` (compared rule-for-rule to the driver's `sg_def` frequency table) and the classification functionals. Draw count 2×10⁵.

#### S5 — E|Ĥ|

`worked-predictions` L1846–1847:

```r
sel_c <- p_sel / det_rate                       # selection | detection
EnH   <- n_rep * sum(sel_c * Pg)                # expected |Hhat| | detection
```

Definition: selection-weighted **population prevalence**, converted to **subjects at n** by the factor `n_rep` (500; 700 in tier-2). **Conditional on declaration** — `sel_c` renormalizes by `det_rate`, so non-declaration is excluded by conditioning (it enters only through the conditioning event, not as a zero).

#### S6 — the classification metrics

From the code of `worked-predictions` (L1848–1855) with their printed labels (L1864–1872):

| symbol | printed label | formula | conditions on |
|---|---|---|---|
| `det_rate` | "detection rate (DETECTED share)" | `mean(det)` | — |
| `EnH` | "E[\|Hhat\|] given detection (n_harm mean)" | `n_rep * sum(sel_c * Pg)` | declaration |
| `Esens` | "E[sens]" | `sum(sel_c * sens_g)` | declaration |
| `Espec` | "E[spec]" | `sum(sel_c * spec_g)` | declaration |
| `Eppv` | "E[PPV]" | `sum(sel_c * PQg)` | declaration |
| `Enpv` | "E[NPV]" | `sum(sel_c * npv_m)`, `npv_m = (pHc − (PQ − PQg[m]*Pg[m]))/pHc`, `pHc = 1 − Pg[m]` | declaration |
| `EbetaH` | "E[beta(Hhat)] oriented (betaHhat_H mean)" | `sum(sel_c * beta_g)` | declaration |
| `Enaive_bias` | "E[naive − beta(Hhat)] (naive bias vs beta(Hhat))" | `mean(noise[cbind(which(det), winner[det])])`, `noise = Bhat − beta_g` | declaration |
| `mass_below` | "selection mass on rules with true mean below the floor" | `sum(sel_c[beta_g < c1])` | declaration |

All are selection-weighted **population** functionals of the family (`sens_g` etc. are population overlaps on the cell law), not sample counts.

#### S7 — the inversions

- **c1 for target power** (`floors-single` L1372–1377): objective = the exact single-candidate rate `detect1(c1, c2, m, s)`; root finder `uniroot`, bracket `c(m − 8*s, m + 8*s)`, `tol = 1e-8`. Unattainable target: returns `NA_real_` when `rho >= ceiling` where `ceiling <- pnorm(c2, m, s*sqrt(2), lower.tail = FALSE)^2`; the printout labels 95% "infeasible (c2 ceiling %.3f)" (L1417). Family-level inversions are described in prose only ("the same one-dimensional search", L1640–1642, and Step 6 L2128–2138) — **not coded**.
- **The c2 ceiling** (`anchor` L57): `pnorm((m_Q − 10)/(seQ1000*sqrt(2)))^2` — inputs `m_Q`, the c2 floor 10, `seQ1000`. Same form re-derived locally in `c1_for` and in the L1418 printout.
- **Screening-only clearance** (`anchor` L58): `pnorm((m_Q − 30)/seQ1000)` — inputs `m_Q`, the c1 floor 30, `seQ1000`.
- **Sample size for target marginal clearance** (`design-calculator` L1633–1634): closed form `1000 * (se_n1000 * qnorm(target)/(mu − c1))^2`.
- **n\* for joint (ρ, α) feasibility** (`floors-single` L1401–1402): closed form `1000 * ((z_ρ·sQ1 + z_{1−α}·sQ1·√piQ)/(bint·(1−piQ)))^2`.

#### S8 — `warmup-floors`

The chunk prints four family-declaration figures plus a table:

1. `fam10` — 3-candidate single-split family declaration at (c1=30, c2=10): **MC**, 2×10⁵ draws (`fam_split`).
2. `fam25` — same at c2=25: **MC**, 2×10⁵.
3. screening-only `p_fam` (reprinted from `worked-prediction`): **analytic** (`pmvnorm`).
4. `fam5_10` — 5-member family (ITT rule S and Qc added) at c2=10: **MC**, 2×10⁵.

The kable's per-member columns are **analytic**: `scr5 = pnorm(c1, beta5, se5, lower=F)` and `P1_5 = detect1(c1, 10, beta5[m], se5[m])`.

#### S9 — `worked-sensitivity`

- **Varied**: the two free interior parameters of the 3×3 joint — `x = J31` ∈ {0.3357 (base), 0.325, 0.341} and `y = J22` ∈ {0.0020 (base), 0.0005, 0.0060} — and the distractor prevalence `pV` ∈ {0.42 (base), 0.30, 0.55}; seven variants, one factor at a time (grid, L2004–2009).
- **Held fixed**: all documented anchors (four marginals, P(Q), Jaccard 0.882 — enforced inside `build_J` by `stopifnot`), the floors, n_rep, the anchor scales, the family construction and the 0.12 floor.
- **Draws**: re-drawn per variant but under **common random numbers** — `predict_scn` calls `set.seed(20260827)` before drawing (L1973), so all seven variants share the same normal draws. The chunk saves and restores `.Random.seed` (L1960, 2017) so downstream chunks are unaffected.

#### S10 — `worked-null`

```r
set.seed(20260828)   # pinned: was inheriting the advanced L1653 chain
beta_null <- rep(tauQc, M)
W1n <- matrix(rnorm(Rdraw * M), Rdraw, M) %*% t(root_half) +
       rep(beta_null, each = Rdraw)
W2n <- ...
pass_n <- (W1n + W2n >= 2 * c1) & (W1n >= c2) & (W2n >= c2)
```

- **The scale table under the null**: there is none — the document's null prediction **reuses the alternative cell's covariance root** (`root_half`, built from the M=16 scenario `Sg`) and only collapses every mean to `tauQc`. On the measurement side, the driver likewise has no null scale: `dgm_scale <- if (isTRUE(null_cell)) NULL else fs_dgm_scale(dgm)`.
- **The `dafff647` guard**: located **in the drivers, not in DOC** — five edits (driver L348 `dgm_scale <- if (isTRUE(null_cell)) NULL else fs_dgm_scale(dgm)`; guarded `.payload$scale` / `.cpayload$scale` writes; guarded SCALE-CHECK condition; the null-branch assertion replaced by `stopifnot(sum(inQ) == 0L, is.na(truth$effect_Q), isTRUE(all.equal(truth$beta_inter, 0)))`). What it guards: under `model = "null"` Q is empty, `fs_dgm_scale(dgm)` has no Q row, and the old assertion compared against `effect_Q = NA` — the branch could not execute.
- **False declaration as coded**: `mean(rowSums(pass_n) > 0)` — family-level single-split declaration when every candidate mean equals the common effect (a "near-threshold null", the floor 0.13–0.28 SEs above the means, per the prose at L2056–2060).
- **`P1` range**: `p1n <- colMeans(pass_n)`; printed `min(p1n)`–`max(p1n)` with `max(.mc_se_prop(p1n, Rdraw))`.
- **`L_eff`**: `log(1 - mean(rowSums(pass_n) > 0)) / log(1 - max(p1n))` — the implied effective family size against the largest per-candidate rate.

#### S11 — `appendix-tier2`

Tier 2 is **not** a different family or filter rule — it is the **same 18-candidate scenario law re-evaluated at n = 700**: the same `cells`/`conds`, the prevalence floor re-expressed at the new n (`keep7 <- Pg7 >= 60/700`, i.e. ≈8.6%, admitting all M7 = 18 candidates instead of 16), ses rescaled by `sqrt(1000/700)`, same gate (30, 10), same maxeffCons argmax, R7 = 2×10⁵, seed 20260825. Its predictions table carries a typed measured column (A72). The "top-7 selection frequencies" of the task list do not exist as such: the object is named `top7` but holds the **top 3** — `top7 <- order(-sel)[1:3]` (L2367), printed with per-share MC SEs. Selection frequencies are `sel <- tabulate(win[det], M7)/R7` (L2345), same construction as S4.

#### S12 — `anchor` and the assertions

Every `stopifnot()` in the document (there are no `if (…) stop()` calls in DOC):

| site | asserts | against |
|---|---|---|
| L38 | payload file exists | `.pay_path` |
| L40 | `.scale` inherits `"fs_dgm_scale"` | class |
| L43/61/67 | Q / Qc / S rows exist | `match()` non-NA |
| L73–88 | regions coherence: `P_g[S]==1`; `n_g[S]==n_super`; `P_g[Q]+P_g[Qc]==1`; `n_g[Q]+n_g[Qc]==n_g[S]`; one common sign of `m_tau`; the mixture identity `P_Q·m_Q + P_Qc·m_Qc == m_S`; `v_cond0` region-invariant | the payload's own regions table |
| L755 (`itt-live`) | undetected-branch `betaHhat_Hc` is a single value equal to `tauQc + beta_int·π` within 1e-8 | pilot bundle |
| L1744–1746 | `J3` sums to 1; row marginals (.455,.045,.500); col marginals (.727,.023,.250) | typed joint |
| L1748 | `PQ` matches the payload `piQ_fix` within 1e-4 | anchor |
| L1751 | Jaccard(g\*,Q) = 0.882 within 1e-3 | driver-documented anchor |
| L1967–1969 (`build_J`) | completed joint non-negative; sums to 1 (1e-9); both marginal sets (1e-9) | construction |

`.anchor` elements **read from the payload**: `sigma`, `piQ`, `bracket`, `V_eff`, `m_Q`, `m_Qc`, `V_mu0_Q`, `m_ITT`, `v_cond0` (the last two confirming v5 §2's "read" claim, from `.reg$m_tau[.iS]` and `.reg$v_cond0[.iQ]`). **Derived**: `seQ1000`, `sd_eff`, `c2_ceiling`, `scr_clear_Q`, `beta_int`.

`.mc_se_prop` is defined at `anchor` L101 (`function(p, R) sqrt(p * (1 - p) / R)`) and applied at exactly six sites: `warmup-floors` L1609 (fam10, fam25), L1612 (fam5_10); `worked-predictions` L1863 (det_rate); `worked-null` L2072 (family rate), L2075 (max over p1n); `appendix-tier2` L2370 (top-3 selection shares).

#### S13 — the draw structure

| variable | defined | value | feeds |
|---|---|---|---|
| `Rg` | `warmup-floors` L1576 | 2e5 | `fam_split` → fam10, fam25, fam5_10 (and their MC SEs) |
| `Rdraw` | `worked-predictions` L1827 | 2e5 | W1/W2 → every Step-4 quantity; **reused by `worked-null`** for W1n/W2n and its MC SEs |
| `R` (formal of `predict_scn`) | `worked-sensitivity` L1972 | default 2e5 | the 7-variant sweep |
| `R7` | `appendix-tier2` L2337 | 2e5 | tier-2 draws and MC SEs |

`set.seed()` sites and what each seeds:

| seed | site | seeds |
|---|---|---|
| 20260824 | `worked-prediction` L1511 | `pmvnorm`'s quasi-MC randomization (the gap readout) |
| 20260826 | `warmup-floors` L1572 | the `fam_split` draws (all three family rates) |
| 20260825 | `worked-scenario` L1735 | (chunk draws nothing itself; pins the state formerly inherited downstream) |
| 20260825 | `worked-predictions` L1826 | W1/W2 — "pinned: identical to the state formerly inherited from L1653" |
| 20260827 | inside `predict_scn` L1973 | common random numbers across the 7 sweep variants |
| 20260828 | `worked-null` L2063 | W1n/W2n — "pinned: was inheriting the advanced L1653 chain" |
| 20260825 | `appendix-tier2` L2312 | tier-2 W1/W2 |

Pairing: within each chunk one set of draws serves all its quantities (P1, selection, EnH, classification, bias are all functionals of the same W1/W2). Across chunks, each re-draws under its own pinned seed — except `worked-null`, which re-draws but **shares `root_half` and `Rdraw`** with `worked-predictions`, and the sensitivity sweep, which shares draws across variants by re-seeding. `worked-sensitivity` additionally saves/restores `.Random.seed` so its insertion leaves downstream chunks bit-identical.

#### S14 — measured comparators embedded in the document

**The `m500` block** (`appendix-n500` L2410–2417), with its comment:

```r
# Deterministic: no random draws, so downstream chunks are unaffected.
m500 <- list(detection = 1.000,      # 1000/1000, tabulated from the payload
             nH = 72, sens = 0.17, spec = 0.87, ppv = 0.41, npv = 0.66,
             oracle = 39.43, oracle_bias = 7.67, oracle_sd = 20.21,
             oracle_se = 19.81,
             naive = 106.04, naive_bias = 74.27, naive_sd = 18.04,
             naive_w = 121.66, naive_cov = 0.23,
             mr = 52.01, mr_bias = 20.25, mr_cov = 0.98)
```

Provenance as far as the comment and history establish: the in-code comment attributes only the detection line to "the payload"; the section prose (L2404–2407) attributes the rest to "that report's estimation table (harm block) and its identification summary", i.e. the rendered `sim_fs_maxeffCons_mr_md40_knoise0_n500_batch_1_1000` report, whose payload is the tracked n=500 bundle of §1.4. `git log -L 2409,2441:DOC` shows the block entered at `a36a8d13` (2026-08-25, "…the n=500 companion") and was touched by `3bedf5e0` (2026-08-25, "Reconcile the 997/1000 detection citation against its investigated source"). Which payload *columns* each literal was tabulated from is **not stated in the file** beyond the detection line — CANNOT DETERMINE from source alone, because the numbers are transcribed from a rendered report, not read from the payload (this is the standing m500-conversion residue item of v5 §6, which stays where it is).

**Other typed comparisons to measurement**:
- `meas` in `appendix-identity-locks` (L2246–2252): the n=700 record, cited in-code to "report sections 6.3 (estimation table, harm block) and 6.7 (classification kable); sizes from the 6.3 block header" of the rendered n=700 report.
- `measd` in `appendix-tier2` (L2355–2358): `c(0.999, 74, 0.123, 0.904, 0.402, 0.662, 31.76, 77.09)` — same n=700 record (31.76 is the locks' back-out, not a report cell).
- `clearance-example` (L1158–1159): "measured (overnight O1): 36/50 = 0.720, MC SE 0.064" — typed from the O1 overnight anchor run inside a `cat` string.
- The `itt-live` chunk's comparison is live (reads the pilot bundle) but the surrounding prose (L735–746) types the same identities' observed values.

### A4. Constants and knobs

Source: index §5 (every numeric literal, by chunk, from the parsed AST) plus YAML. **The YAML has no `params:` block** (DOC L1–18: title/subtitle/author/date/format/execute only), so the table is exactly the index's literal inventory, classified. Classes: `fixture` / `knob` / `math` / `display`.

| chunk | literal | role (from the code) | class |
|---|---|---|---|
| anchor | 0 | `.fmt_tex(digits = 0)` default | display |
| anchor | 1 | `p*(1-p)`; `P_g` sums to 1; `length(unique(sign)) == 1L` | math |
| anchor | 10 | the c2 floor inside `c2_ceiling` | knob |
| anchor | 30 | the c1 floor inside `scr_clear_Q` | knob |
| anchor | 1000 | reference n in `seQ1000` | knob |
| anchor | 2 | `sqrt(2)`, `^2` | math |
| itt-live | 0 | `detected %in% 0L` (undetected branch) | math |
| itt-live | 1 | `length(obs_itt) == 1L` | math |
| itt-live | 1e-08 | lock tolerance | math |
| itt-live | 9 | `round(., 9)` | display |
| clearance-example | 30 | c1 floor | knob |
| clearance-example | 40 | oriented truth on Q, typed | fixture |
| floors-single | 0 | integrate/contour bounds; `abline(0,1)` | math |
| floors-single | 0.05 / 0.1 / 0.2 | the α of the (ρ, α) pairs in `nstar` | knob |
| floors-single | 0.5 | contour level; grid step `by = 0.5` | display |
| floors-single | 0.6 / 0.7 / 0.8 / 0.9 / 0.95 | contour levels; ρ targets (.80/.90/.95 in `inv` and `nstar`) | knob |
| floors-single | 0.743 | contour level marking the fixture's own P1 | fixture |
| floors-single | 1 | identities; `abline(0, 1)` | math |
| floors-single | 10 | c2 floor (columns, `tr`, `inv`, fixture point) | knob |
| floors-single | 1000 | reference n | knob |
| floors-single | 19 | `pch = 19` | display |
| floors-single | 1e-08 | `uniroot` tol | knob |
| floors-single | 1e-10 | `integrate` rel.tol | knob |
| floors-single | 2 | `sqrt(2)`, `2*c1`, `seq(…, 2)` | math |
| floors-single | 22 / 44 | contour grid `g1 = seq(22, 44, .5)` | knob |
| floors-single | 25 | alternative c2 | knob |
| floors-single | 28 / 42 | trade-off grid `c1 = seq(28, 42, 2)` | knob |
| floors-single | 3 | `seq_len(3)` channels; `lty = 3`; round digits | math |
| floors-single | 30 | c1 floor | knob |
| floors-single | 4 | `tr[, 2:4]` | display |
| floors-single | 40 | contour grid `g2 = seq(0, 40, .5)` upper | knob |
| floors-single | 5 | `res[, 3:5]` | display |
| floors-single | 500 / 700 | the two design sample sizes | knob |
| floors-single | 8 | `uniroot` bracket half-width `8*s` | knob |
| floors-single | Inf | integrate upper limit | math |
| floors-single | NA | `NA_real_` returns | math |
| worked-prediction | 0 | contrast lower bounds | math |
| worked-prediction | 0.28 | \|P∩Q\| overlap | fixture |
| worked-prediction | 0.31 | prevalence of D | fixture |
| worked-prediction | 0.45 | prevalence of P | fixture |
| worked-prediction | 1 | purity 1.00; diag | math |
| worked-prediction | 1e-06 | `GenzBretz(abseps)` | knob |
| worked-prediction | 20260824 | seed | knob |
| worked-prediction | 3 | M = 3 family size (diag(3), 1:3, matrix dims) | fixture |
| worked-prediction | 30 | c1 | knob |
| worked-prediction | 4 | round digits | display |
| worked-prediction | Inf | `pmvnorm` upper | math |
| warmup-floors | 0 | ov5 zero overlaps (Qc with Q, D) | fixture |
| warmup-floors | 0.28 / 0.45 | Qc∩P = 0.45 − 0.28 | fixture |
| warmup-floors | 1 | prevalence of S; identities | math |
| warmup-floors | 10 | c2 | knob |
| warmup-floors | 2 | `2*Sg`, `2*c1` | math |
| warmup-floors | 20260826 | seed | knob |
| warmup-floors | 25 | alternative c2 | knob |
| warmup-floors | 2e+05 | `Rg` draw count | knob |
| warmup-floors | 3 / 5 | family sizes M | fixture |
| warmup-floors | 4 | round/format digits; ov5 index 4 (S row) | display / fixture (index) |
| design-calculator | 0.8 / 0.9 / 0.95 | clearance targets | knob |
| design-calculator | 1000 | reference n | knob |
| design-calculator | 2 | `^2` | math |
| design-calculator | 30 | default c1 | knob |
| design-calculator | 40 | default mu (fixture truth) | fixture |
| design-calculator | 80 / 90 / 95 | printed percent labels | display |
| worked-scenario | 0.3553 / 0.0121 / 0.0876 / 0.0360 / 0.0020 / 0.0070 / 0.3357 / 0.0089 / 0.1554 | the nine J3 cells (hand-solved interior) | fixture |
| worked-scenario | 0.455 / 0.045 / 0.5 | age row marginals | fixture |
| worked-scenario | 0.727 / 0.023 / 0.25 | preanti column marginals | fixture |
| worked-scenario | 0.42 | stylized P(V=1) | fixture (stylized) |
| worked-scenario | 0.882 | Jaccard(g*, Q) anchor | fixture |
| worked-scenario | 0.12 | prevalence floor = 60/500 | knob |
| worked-scenario | 0.9 | `pcons_gate` (declared, not used downstream in-chunk) | knob |
| worked-scenario | 10 | c2 | knob |
| worked-scenario | 30 | c1 | knob |
| worked-scenario | 500 | n_rep | knob |
| worked-scenario | 1e-12 / 1e-04 / 0.001 | assertion tolerances (0.001 = 1e-3 on the Jaccard) | math |
| worked-scenario | 20260825 | seed | knob |
| worked-scenario | 1 / 2 / 3 / 4 / 5 / 6 | cell/condition indices, pair list, `J3[3,1]` etc. | math (indices) |
| worked-predictions | 2e+05 | `Rdraw` | knob |
| worked-predictions | 20260825 | seed | knob |
| worked-predictions | 2 | `scale = 2`, `(W1+W2)/2`, `2*c1` | math |
| worked-predictions | 6 | top-6 display cut | display |
| worked-predictions | 0 / 1 / 3 | `rowSums > 0`; indices; round digits | math / display |
| worked-predictions | Inf / NA | mask value; NA winner | math |
| worked-sensitivity | 0.3357 / 0.0020 / 0.42 | base interior (repeated per grid row) | fixture |
| worked-sensitivity | 0.325 / 0.341 / 0.0005 (5e-04) / 0.006 / 0.3 / 0.55 | sweep variant values | knob |
| worked-sensitivity | 0.882 | Jaccard anchor in `build_J` | fixture |
| worked-sensitivity | 0.455 / 0.045 / 0.5 / 0.727 / 0.023 / 0.25 / 0.1554 / 0.0946 | marginal and derived-cell constants in `build_J` | fixture |
| worked-sensitivity | 0.12 | prevalence floor | knob |
| worked-sensitivity | 2e+05 | R default | knob |
| worked-sensitivity | 20260827 | common-random-numbers seed | knob |
| worked-sensitivity | 1e-09 / 1e-12 | `build_J` tolerances | math |
| worked-sensitivity | 0 / 1 / 2 / 3 / 4 / 5 / 6 | indices, `2*c1`, `(A+B)/2` | math |
| worked-null | 20260828 | seed | knob |
| worked-null | 0 / 1 / 2 | `rowSums > 0`; `1 −`; `2*c1` | math |
| appendix-identity-locks | 0.123 / 0.402 / 74 / 700 / 40.66 / 108.85 / 53.94 / 8.9 / 77.09 / 22.17 / 17.57 / 121.4 / 0.18 / 22.37 / 144.49 / 0.97 | the typed n=700 measured record | fixture |
| appendix-identity-locks | 2 | `w/2` | math |
| appendix-identity-locks | 3 | round digits | display |
| appendix-tier2 | 60 | n.min | knob |
| appendix-tier2 | 700 | n | knob |
| appendix-tier2 | 1000 | se rescale reference | knob |
| appendix-tier2 | 2e+05 | R7 | knob |
| appendix-tier2 | 20260825 | seed | knob |
| appendix-tier2 | 0.999 / 74 / 0.123 / 0.904 / 0.402 / 0.662 / 31.76 / 77.09 | typed measured column | fixture |
| appendix-tier2 | 0 / 1 / 2 / 3 / 4 / 5 / 6 / Inf | gate/mask/index arithmetic; top `[1:3]` | math / display |
| appendix-n500 | 1 / 72 / 0.17 / 0.87 / 0.41 / 0.66 / 39.43 / 7.67 / 20.21 / 19.81 / 106.04 / 74.27 / 18.04 / 121.66 / 0.23 / 52.01 / 20.25 / 0.98 / 74.27 | the typed `m500` record (1 = detection 1.000) | fixture |
| appendix-n500 | 500 | n in the classification identity | knob |
| appendix-n500 | 2 | `w/2`, round(., 2) | math / display |
| appendix-n500 | 3 | round digits | display |

(The 08-27 literal sweep covered literals duplicating payload quantities; the residual fixture-typed values above — 40, 0.28/0.31/0.45, the J3 interior, 0.882, and the three measured-record blocks — are the fixture constants this design-sweep classification adds.)

### A5. The package boundary and the in-document helpers

**`forestsearch::` functions the document calls** (index §6): exactly one — `forestsearch::fs_sym_root` (plus `knitr::kable`/`knitr::opts_chunk`, and `library(mvtnorm)` in two chunks).

`fs_sym_root` — file `R/fs_sym_root.R`; signature and body verbatim (L65–70):

```r
fs_sym_root <- function(S, scale = 2) {
  eS <- eigen((S + t(S)) / 2, symmetric = TRUE)
  V  <- eS$vectors
  d  <- sqrt(pmax(scale * eS$values, 0))
  V %*% diag(d, nrow = length(d)) %*% t(V)
}
```

Returns: an unnamed symmetric numeric matrix R with `R %*% t(R) = scale * sym(S)` (no list structure). Call sites and the A2-criterion marking of each argument:

| site | call | `S` | `scale` |
|---|---|---|---|
| warmup-floors L1575 | `fs_sym_root(2 * Sg, scale = 1)` | fixture-bound (M=3 typed-family covariance) | general |
| warmup-floors L1597 | `fs_sym_root(2 * S5, scale = 1)` | fixture-bound (5-member covariance) | general |
| worked-predictions L1833 | `fs_sym_root(Sg, scale = 2)` | fixture-bound (M=16 scenario covariance) | general |
| worked-sensitivity L1991 | `fs_sym_root(Sg, scale = 2)` | fixture-bound (variant covariance) | general |
| appendix-tier2 L2336 | `fs_sym_root(Sg7, scale = 2)` | fixture-bound (M7=18 covariance) | general |

The operation itself is G — it consumes only a covariance matrix.

**Functions the document defines itself** (index §4), one line each, with the R/ counterpart check (searched `R/` for `mc_se`, `sqrt(p * (1 - p)`, `detect1`, `integrate(`, thousands formatting):

| helper | defined | what it does | counterpart in R/ |
|---|---|---|---|
| `.fmt_tex(x, digits)` | anchor L92–94 | `formatC` with `big.mark = "{,}"` (LaTeX thousands) | no counterpart found in R/ |
| `.mc_se_prop(p, R)` | anchor L101 | `sqrt(p*(1-p)/R)` | no counterpart found in R/ |
| `detect1(c1, c2, m, s)` | floors-single L1366–1371 | exact single-candidate split-pair rate by 1-D quadrature | no counterpart found in R/ |
| `c1_for(rho, c2, m, s)` | floors-single L1372–1377 | uniroot inversion of `detect1` with c2-ceiling NA | no counterpart found in R/ |
| `zr(x)` | floors-single L1400 | alias for `qnorm` | (trivially `stats::qnorm`) |
| `nstar(rho, alpha)` | floors-single L1401–1402 | closed-form n for joint (ρ, α) separation | no counterpart found in R/ |
| `n_for(target, se_n1000, mu, c1)` | design-calculator L1633–1634 | closed-form n for marginal clearance | no counterpart found in R/ |
| `fam_split(mu, root, Mloc, c2v)` | warmup-floors L1577–1581 | MC family single-split declaration rate | no counterpart found in R/ |
| `build_J(x, y)` | worked-sensitivity L1961–1971 | closed-form 3×3 joint from the two free interior cells | no counterpart found in R/ (scenario-specific) |
| `predict_scn(x, y, pV, R)` | worked-sensitivity L1972–2003 | the full S1→S6 pipeline for one interior variant | no counterpart found in R/ — this is the closest thing to the future wrapper that exists anywhere, and it lives in a chunk |

None of these duplicates logic present in `R/` (I read `fs_sym_root.R`, `fs_dgm_scale.R` headers, `fs_mr_oc_summary.R`, and searched the quoted patterns; no package function computes a split-pair rate, a family declaration, or an MC SE of a proportion).

---

## 3. Part B — what the measurement side records about the realized family

### B1. Drivers and payloads

Both drivers are **byte-identical except one line** (`diff`: only L133, `n_sample <- 500L` vs `700L`). Shared `git log -1`: `d884adbf 2026-08-27 10:45:33 -0700 fix(sims): save-site name vectors tolerate an absent scale element`.

The `forestsearch()` call, verbatim (n500 driver L561–591; comments elided at `# ALIGNED:` block L578–587 for space — the argument is `parallel_args = inner_parallel`):

```r
fs.est <- suppressWarnings(forestsearch(
  df.analysis = df, confounders.name = confs,
  outcome.name = outcome_name, treat.name = treat_name, id.name = id_name,
  outcome_type = "continuous", effect_measure = "MD",
  effect.threshold = md_threshold, consistency.threshold = md_consistency,
  pconsistency.threshold = pconsistency, fs.splits = fs_splits,
  n.min = n_min, d0.min = d0_min, d1.min = d1_min, maxk = maxk,
  vi.grf.min = vi_grf_min, sg_focus = sg_focus,
  selection_rule = selection_rule,
  effect_neighborhood = effect_neighborhood,
  stop_threshold = stop_threshold,
  consistency_method = consistency_method,
  conf.cont_jcuts = fs_conf.cont_jcuts,
  use_lasso = use_lasso, use_dina = use_dina, use_grf = use_grf,
  use_twostage = use_twostage, is.RCT = is_rct,
  adverse_outcome = adverse_outcome,
  details = FALSE, quiet = FALSE, seedit = sd_i,
  parallel_args = inner_parallel,
  mr_inference = TRUE,
  mr_inference_args = list(ci_method = "ij", draws = mr_draws,
                           include_complement = TRUE)))
```

with the argument values from the setup chunk (L106–216): `md_threshold = 30`; `md_consistency = 10`; `pconsistency = 0.90`; `fs_splits = 400L`; `maxk = 2L`; `n_min = 60L`; `d0_min = 12L`; `d1_min = 12L`; `consistency_method = "resample"`; `sg_focus = "maxeffCons"`; `subgroup_method = "consistency"`; `selection_rule = "neighborhood"`; `effect_neighborhood = 0.10`; `stop_threshold = NULL`; `fs_conf.cont_jcuts = list(age = 10, preanti = 10)`; `use_lasso/use_dina/use_grf = FALSE`; `use_twostage = TRUE`; `is_rct = TRUE`; `vi_grf_min = -0.2`; `adverse_outcome = FALSE`; `mr_draws = 5000L`; `seed_base = 8316951L`, `sd_i = seed_base + sim_id`.

`saveRDS()` sites, verbatim:

Batch mode (L874–877):

```r
if (!is.null(dgm_scale)) .payload$scale <- dgm_scale
.payload$oc    <- fs_mr_oc_summary(.payload)
saveRDS(.payload[intersect(c("results", "truth", "scale", "oc", "meta"),
                           names(.payload))], out)
```

Combine mode (L967–970):

```r
if (!is.null(dgm_scale)) .cpayload$scale <- dgm_scale
.cpayload$oc    <- fs_mr_oc_summary(.cpayload)
saveRDS(.cpayload[intersect(c("results", "truth", "scale", "oc", "meta"),
                          names(.cpayload))], cpath)
```

### B2. The per-replicate record

The recorder template, verbatim (n500 driver L463–488):

```r
.na_record <- function(sim_id) data.frame(
  sim_id = sim_id, status = NA_character_, detected = NA_integer_,
  mr_ok = NA_integer_, err_msg = NA_character_, mr_msg = NA_character_,
  n_sel = NA_integer_, n_harm = NA_integer_, n_true = NA_integer_,
  sg_def = NA_character_, covs = NA_character_,
  betaHhat_H = NA_real_, betaHhat_Hc = NA_real_,
  fb_secs = NA_real_, fit_mr_secs = NA_real_, fb_err = NA_character_,
  fb_src1 = NA_real_, fb_src2 = NA_real_, fb_nres = NA_integer_,
  or_H_est = NA_real_, or_H_lo = NA_real_, or_H_hi = NA_real_, or_H_se = NA_real_,
  or_Hc_est = NA_real_, or_Hc_lo = NA_real_, or_Hc_hi = NA_real_, or_Hc_se = NA_real_,
  nv_H_est = NA_real_, nv_H_lo = NA_real_, nv_H_hi = NA_real_, nv_H_se = NA_real_,
  nv_Hc_est = NA_real_, nv_Hc_lo = NA_real_, nv_Hc_hi = NA_real_, nv_Hc_se = NA_real_,
  fb_H_est = NA_real_, fb_H_lo = NA_real_, fb_H_hi = NA_real_, fb_H_se = NA_real_,
  fb_Hc_est = NA_real_, fb_Hc_lo = NA_real_, fb_Hc_hi = NA_real_, fb_Hc_se = NA_real_,
  mr_H_est = NA_real_, mr_H_lo = NA_real_, mr_H_hi = NA_real_, mr_H_se_ij = NA_real_,
  mr_Hc_est = NA_real_, mr_Hc_lo = NA_real_, mr_Hc_hi = NA_real_, mr_Hc_se_ij = NA_real_,
  ij_source = NA_character_,
  sens = NA_real_, spec = NA_real_, ppv = NA_real_, npv = NA_real_,
  stringsAsFactors = FALSE)
```

(the schema-change comment block at L471–477 elided; `fs_attach_betaHhat` later appends `betaHhat_status`, `nH_eval`, `nHc_eval`). The fields written per replicate beyond the template seeds: `n_true` (L552), oracle rows (L554–557), timing/messages, `status`/`detected` (L601–604), `sg_def <- paste(fs.est$sg.harm, collapse = " & ")` (L612), `n_harm <- sum(sgv == 1L)` with `sgv <- fs.est$grp.consistency$sg.harm.id` (L626–628), the `.classify` outputs (L629), and the MR block reads (L642–667), including:

```r
rec$n_sel     <- g$n_selected %||% NA_integer_
```

**Does any field record the realized candidate-family size?** **No.** Every count field is a subject count: `n_sel = g$n_selected = sz[sel]` — the number of *subjects* in the MR-selected subgroup (`R/fs_mr_inference.R:603`; empirically `n_sel == n_harm` in all three payloads); `n_harm` = subjects flagged by `sg.harm.id`; `n_true` = subjects with the true harm flag; `nH_eval`/`nHc_eval` = df_super evaluation counts; `fb_nres` = full-bootstrap resample count (NA here). Nothing records candidates enumerated, candidates passing the screen, or rows of the candidate table.

**Where such a field would be written** (proposal only; no edit): in `record_replicate()` immediately after the detection branch — the natural line is beside L626 (`sgv <- fs.est$grp.consistency$sg.harm.id`), reading any of the three counts the return object already carries (see C3): e.g. `rec$n_cand <- nrow(fs.est$find.grps$out.found$hr.subgroups)` (post-screen family), `rec$n_cand_eval <- fs.est$grp.consistency$n_candidates_evaluated` / `...$n_passed` (consistency stage), or `rec$n_family_mr <- fs.est$mr_inference$n_family` (MR's re-enumerated superset — note this last is a *different* family: all ≤maxk combinations at the `n.min` floor only, without d0/d1 minima, dedup, or the effect screen, per `R/forestsearch_main.R:3260–3272`). A NO-DETECTION replicate would need the first two read before the early return at L602.

### B3. The measured payloads (metadata only)

| payload | `names()` | results nrow × ncol | notes |
|---|---|---|---|
| null cell `…mdnull_…n500_res_1_1000.rds` | `"results" "truth" NA "oc" "meta"` | 1000 × 59 | the third element has **name `NA` and value `NULL`** — written before `d884adbf` replaced the fixed name vector with `intersect()` at the save site (see Open items); no `scale` element, as designed for the null cell |
| alt n=500 `…md40_…n500_res_1_1000.rds` | `"results" "truth" "scale" "oc" "meta"` | 1000 × 59 | `scale` = the `fs_dgm_scale` object the anchor reads |
| alt n=700 `…md40_…n700_res_1_1000.rds` | `"results" "truth" "scale" "oc" "meta"` | 1000 × 59 | same schema |

Common structure (str, max.level 2): `results` the 59-column per-replicate frame (identical column set in all three, listed in B2 plus `betaHhat_status`, `nH_eval`, `nHc_eval`); `truth` a 9-element list (`effect_Q`, `effect_Qc`, `beta_inter`, `prevalence_Q`, `effect_ITT`, `marg_H`, `marg_H_se`, `marg_Hc`, `marg_Hc_se`); `oc` an `fs_mr_oc` list (`estimation` 6×15, `identification` 2×8, `targets`, `meta`); `meta` an 18-element list (n_sample, n_sims, nb_boots = NULL, mr_draws 5000, sg_focus "maxeffCons", selection_rule, consistency_method "resample", subgroup_method, target_md_harm −40, effect_threshold 30, consistency_threshold 10, adverse_outcome FALSE, seed_base 8316951, sim_id_start 1, parallel_mode, n_workers, pkg_version "0.2.2", built_at). Null-cell truth: `effect_Q = NA`, `beta_inter = 0`, `prevalence_Q = 0`, `effect_Qc = effect_ITT = −26.3`.

**Does the per-replicate record carry the selected subgroup's definition** (so detection could be re-scored under the document's S4 selection law)? **Yes** — the `sg_def` column holds the realized rule verbatim (e.g. `"!{cd40 <= 415} & !{cd80 <= 1022}"`), in all three payloads, alongside the subject-level consequences (`n_harm`, sens/spec/ppv/npv, `betaHhat_H`). What it does **not** carry is the rest of the realized family the winner was chosen from.

### B4 — DIAGNOSTIC

**Not available.** B2/B3 found no per-replicate family-size field in any payload, so nothing was computed (per the task's rule, the diagnostic runs only on a stored family-size column).

### B5 — definition parity

Document definition (from A3) beside the driver's scoring expression, verbatim; verdict per row.

**(1) Declaration / false declaration.**
- DOC (S3): `pass <- (W1 + W2 >= 2*c1) & (W1 >= c2) & (W2 >= c2)`; family event `rowSums(pass) > 0`, on a fixed family. Null cell: same event at `beta_null = rep(tauQc, M)`.
- Driver (L601–604): `if (is.null(fs.est) || is.null(fs.est$sg.harm)) { rec$status <- "NO-DETECTION"; rec$detected <- 0L; ... }` else `rec$status <- "DETECTED"; rec$detected <- 1L` — i.e. the search returned a subgroup, which requires a candidate passing `p̂cons ≥ 0.90` over `fs.splits = 400` (resample closed-form rate, `R/consistency_resample.R`: `rate ≈ 2Φ((β̂ − c)/σ_D) − 1` per candidate) after the full enumeration/screen.
- Verdict: **different event** — single-split gate on a fixed family vs the 0.90-over-400-splits gate on the regenerated realized family. The document states this idealization and its direction explicitly (L1806–1823); the comparison is benchmark-grade, not identity-grade.

**(2) Detection.**
- DOC (S4): "detection rate (DETECTED share)" = the same family-declaration event; no Ĥ-vs-Q event exists.
- Driver: `det_rate <- mean(results$detected %in% 1L)` (L1505) — DETECTED share.
- Verdict: **same event** at the level both sides define it (share of replicates on which the procedure declares *some* subgroup), inheriting row (1)'s gate difference. Neither side scores "found Q" as an event; both compare selection *distributions* (DOC `p_sel` vs the driver's `sg_def` frequency table).

**(3) E|Ĥ|.**
- DOC: `EnH <- n_rep * sum(sel_c * Pg)` — selection-weighted population prevalence × n, conditional on declaration.
- Driver: `rec$n_harm <- sum(sgv == 1L, na.rm = TRUE)` — realized subject count on the analysis sample; summarized over detected replicates (report's `n_harm` mean).
- Verdict: **same event** (subgroup size in subjects at n, conditional on detection), evaluated as population expectation vs sample realization.

**(4) Classification metrics.**
- DOC: `Esens/Espec/Eppv/Enpv` = selection-weighted **population** overlaps on the cell law (A2 rows A54).
- Driver (`.classify`, L449–460): `tp/(tp+fn)`, `tn/(tn+fp)`, `tp/(tp+fp)`, `tn/(tn+fn)` from `sg_hat` vs the true harm flag **on the analysis sample**; means over detected reps.
- Verdict: **different event in the population-vs-sample sense** — per-replicate sample proportions vs population functionals. The document itself records that for a *selected* rule sample PPV is systematically biased upward relative to population purity (L2283–2292), so the gap is signed, not noise.

**(5) E[β(Ĥ)].**
- DOC: `EbetaH <- sum(sel_c * beta_g)` — selection-weighted population mixture mean (oriented).
- Driver: `betaHhat_H` via `fs_attach_betaHhat(results, eval_df, focus = "harm", outcome_type = "continuous", effect_measure = "MD")` (L849–851) — exact finite-population mean of τ over Ĥ on `df_super` (raw scale; oriented downstream by `rd$bH_or <- .orient * rd$betaHhat_H`, L998); `oc$estimation` uses `beta = orient * rd$betaHhat_H` (`R/fs_mr_oc_summary.R:130`).
- Verdict: **same event** — both are the population-exact conditional estimand at the realized/selected rule, oriented; the only difference is which rules the selection law ranges over.

**(6) Naive bias.**
- DOC: `Enaive_bias <- mean(noise[cbind(which(det), winner[det])])` — mean winner's noise = E[naive − β(Ĥ)] over declaring draws.
- Driver/OC: `bias_beta = .fs_oc_mean(e - b$beta)` (`R/fs_mr_oc_summary.R:159`) with `e` the oriented naive estimate and `b$beta` the oriented `betaHhat_H`, per replicate, averaged over detected reps.
- Verdict: **same event** (mean of naive minus exact conditional estimand over detected replicates) — with the caveat that DOC's "naive" is the winner's fitted Gaussian draw while the driver's is the refitted `lm` treatment coefficient on the identified subgroup; both are the selected subgroup's unpenalized fit, so I record **same event**, evaluation differing only through the model-vs-Gaussian idealization already counted in row (1).

**(7) Below-floor selection mass** (S6's last row): DOC `mass_below = sum(sel_c[beta_g < c1])`; driver: **no counterpart column** — the document itself notes it is "not tabulated (bundle-computable)" (L2357) and proposes the per-replicate analogue share of `betaHhat_H < c1`. Verdict: **CANNOT DETERMINE** (nothing to compare on the measured side).

---

## 4. Part C — how the search builds its family, from `R/` source at HEAD

### C1. `forestsearch()` signature

`R/forestsearch_main.R:1141–1238`, verbatim (every formal with its default):

```r
forestsearch <- function(df.analysis,
                         outcome.name = "tte",
                         event.name = "event",
                         treat.name = "treat",
                         id.name = "id",
                         potentialOutcome.name = NULL,
                         flag_harm.name = NULL,
                         confounders.name = NULL,
                         parallel_args = list(plan = "multisession",
                                              workers = .default_parallel_workers(),
                                              show_message = TRUE),
                         df.predict = NULL,
                         df.test = NULL,
                         is.RCT = TRUE,
                         seedit = 8316951,
                         est.scale = "hr",
                         use_lasso = FALSE,
                         use_grf = FALSE,
                         grf_res = NULL,
                         grf_cuts = NULL,
                         use_dina = FALSE,
                         dina_res = NULL,
                         dina_cuts = NULL,
                         dina_args = list(),
                         dina_select_statistic = c("effect", "dina"),
                         subgroup_method = c("consistency", "dina", "grf"),
                         max_n_confounders = 1000,
                         grf_depth = 2,
                         grf_selection = c("frontier", "tree"),
                         grf_select_statistic = c("effect", "dr"),
                         dmin.grf = 0.0,
                         frac.tau = 0.8,
                         return_selected_cuts_only = TRUE,
                         conf_force = NULL,
                         defaultcut_names = NULL,
                         cut_type = "default",
                         exclude_cuts = NULL,
                         replace_med_grf = FALSE,
                         cont.cutoff = 4,
                         conf.cont_medians = NULL,
                         conf.cont_medians_force = NULL,
                         conf.cont_jcuts = NULL,
                         collapse_cuts = TRUE,
                         collapse_cuts_args = list(),
                         n.min = 60,
                         n.min.frac = 0.10,
                         effect.threshold = NULL,
                         consistency.threshold = NULL,
                         hr.threshold = 1.25,
                         hr.consistency = 1.0,
                         sg_focus = "hr",
                         selection_rule = "neighborhood",
                         effect_neighborhood = 0.10,
                         fs.splits = 1000,
                         m1.threshold = Inf,
                         pconsistency.threshold = 0.90,
                         stop_threshold = pconsistency.threshold,
                         show_candidate_summary = FALSE,
                         max_print = 10,
                         d0.min = 10,
                         d1.min = 10,
                         max.minutes = 3,
                         minp = 0.025,
                         details = FALSE,
                         quiet = FALSE,
                         maxk = 2,
                         by.risk = 12,
                         plot.sg = FALSE,
                         plot.grf = FALSE,
                         max_subgroups_search = Inf,
                         vi.grf.min = -0.2,
                         use_twostage = TRUE,
                         twostage_args = list(),
                         consistency_method = c("resample", "split"),
                         outcome_type = c("survival", "binary", "continuous",
                                          "count"),
                         effect_measure = NULL,
                         offset.name = NULL,
                         adverse_outcome = NULL,
                         overdispersion = c("none", "quasi", "negbin"),
                         grf_count_transform = c("log", "identity"),
                         tune_grf = FALSE,
                         adjust_covariates = NULL,
                         ps_method = NULL,
                         ps_adjust_method = c("none", "iptw", "dr_gcomp"),
                         ps_hat = NULL,
                         mr_inference = FALSE,
                         mr_inference_args = list()) {
```

(Comment lines inside the formals at L1197–1201 elided; `stop_threshold = pconsistency.threshold` is written as a promise, per the in-source comment.)

### C2. The candidate-family pipeline, data in → Ĥ out

"Data-dependent at fixed n" = its outcome varies across replicates drawn from the same DGM.

| step | file / function | governing arguments (defaults) | data-dependent? | evidence |
|---|---|---|---|---|
| 1. Cut construction (continuous → binary factors, both directions) | `R/get_fsdata.R` `get_FSdata()`; `R/get_FSdata_helpers.R` `cut_var_jq()` | `cut_type = "default"` (q25/median/q75); `conf.cont_jcuts` (J-quantile grids **replace** defaults for named vars); `cont.cutoff = 4`; `collapse_cuts = TRUE`; `conf.cont_medians*`, `conf_force`, `exclude_cuts` | **yes** — cuts are placed at *empirical* quantiles of the replicate's sample: `cut_var_jq()` emits `"x <= qj(x, k, J+1)"` expressions | `R/get_FSdata_helpers.R:107–115`: `paste0(x, " <= qj(", x, ", ", ks, ", ", J + 1L, ")")`; driver comment L197–199 |
| 2. GRF variable-importance pre-screen | `forestsearch_main.R` §5 (L2848–2833 region) | `vi.grf.min = -0.2` (ratio floor; `NULL` disables), `max_n_confounders = 1000`, `tune_grf` | **yes** — a causal forest fit on the replicate ranks factors; with `vi.grf.min = -0.2` the ratio test `vi_ratio > vi.grf.min` keeps everything with positive max VI, but ordering/cap are sample-fit-dependent | `R/forestsearch_main.R:2808–2830`: `selected.vars <- which(vi_ratio > vi.grf.min)` |
| 3. Dummy expansion — both directions of each cut | `forestsearch_main.R` §6 | — | no (mechanical) | L2841–2846: "dummy() expands each 2-level factor into two indicator columns … the search explores BOTH directions" |
| 4. Rule enumeration: all combinations of ≤ maxk factor columns | `R/subgroup_search.R` `generate_combination_indices()` inside `subgroup.search()` | `maxk = 2` | no (combinatorial in L; L itself depends on step 1–2 outputs) | `subgroup_search.R:339–370` |
| 5. Per-combination filter chain | `subgroup_search.R` `evaluate_combination_with_status()` | `minp = 0.025`; `rmin = 5`; `d0.min`/`d1.min` (**skipped for continuous**); `n.min = 60` | **yes** — every check is on the replicate sample | L553–673, quoted below |
| 6. Effect screen | same | `hr.threshold` ← `effect.threshold` (driver: 30); `disable_effect_floor` under maxeff | **yes** — per-subgroup `lm` fit on the sample: keep if `glm_result$hr > threshold` | L626–629 |
| 7. Candidate table assembly + preview sort | `subgroup_search.R` `format_search_results()` | — | yes (contents) | L876–881: `names(hr.out) <- c("grp","K","n","E","d1","m1","m0","HR","L(HR)","U(HR)", colnames(Z))`; `setorder(hr.out, -HR, K)` |
| 8. Consistency-stage re-filter | `R/subgroup_consistency_main.R` L511–534 | `hr.threshold` (= `effect_threshold` on GLM), `m1.threshold = Inf`; **maxeff keeps everything** | yes | `found.hrs <- hr.subgroups[hr.subgroups$HR >= hr.threshold, ]` |
| 9. Near-duplicate removal | `subgroup_consistency_main.R` L559–585 | focus-dependent: maxeff → `.maxeff_membership_dedup()` (exact-membership); **all other foci incl. maxeffCons → `remove_near_duplicate_subgroups()`** (statistics-keyed) | yes | quoted at A3-source read; branch on `identical(sg_focus, "maxeff")` |
| 10. Preview ordering + pool cap | `R/subgroup_consistency_helpers.R` `sort_subgroups_preview()`; cap L610–626 | maxeffCons: `setorder(result_new, -HR, K)`; `stop_Kgroups = max_subgroups_search = Inf` (driver default) | yes (order is sample HR) | helpers L641–664; main L625 `found.hrs <- found.hrs[seq_len(maxsgs), ]` |
| 11. Consistency screen | `subgroup_consistency_main.R` loop → `evaluate_subgroup_consistency()`; `R/consistency_resample.R` (`consistency_method = "resample"`); `"split"` = repeated 50/50 refits | `fs.splits` (`n.splits = 400`), `pconsistency.threshold = 0.90`, `consistency.threshold` (c2 = 10), `use_twostage = TRUE` (`twostage_args`: screen at `n.splits.screen = 30`, `screen.threshold = max(0.5, pcons − 2.5·SE)`), `stop_threshold` | **yes** — `rate ≈ 2Φ((β̂ − c)/σ_D) − 1` from the replicate's fit and dfbeta | `consistency_resample.R:1–29` (header formula, quoted below); `subgroup_consistency_main.R:480–491` |
| 12. Winner selection | `subgroup_consistency_helpers.R` `sort_subgroups()` | `sg_focus = "maxeffCons"`: `setorder(result_new, -hr, K)` **among consistency-qualifying rows only** ("every row reaching here already cleared pconsistency.threshold, because evaluate_subgroup_consistency() returns NULL below it") | yes | helpers L570–574 (comment verbatim), quoted below |
| 13. Time cap | `subgroup.search` | `max.minutes = 3` | potentially (a slow replicate truncates enumeration) | `subgroup_search.R` formal `max.minutes = 30` (function default; forestsearch passes 3) |

Key filter evidence (`R/subgroup_search.R:553–673`, GLM/continuous branch):

```r
  # Status 1: Check prevalence
  if (!meets_prevalence_threshold(x, minp)) ...        # colMeans(x) >= minp, per FACTOR column
  # Status 2: Check redundancy
  redundancy_check <- extract_idx_flagredundancy(x, rmin)   # each added factor must shrink
                                                            # membership by > rmin subjects
  ...
  # Continuous: skip Status 3 entirely — d0.min/d1.min do not apply
  # Status 4: Check total subgroup size
  nx <- sum(id.x)
  if (nx <= n.min) return(list(status = 4L, result = NULL))
  # Status 5: Fit GLM via estimator closure
  # Status 6: Check effect threshold. ... if (glm_result$hr <= hr.threshold) ...
```

`meets_prevalence_threshold` (L688–691): `all(colMeans(x) >= minp)` — a **per-constituent-factor** sample-prevalence floor, not a subgroup floor. The subgroup-size floor is `nx > n.min` in **subjects of the analysis sample**; there is **no maximum-size constraint anywhere in the chain** (no upper bound on `nx` or `Pg`).

Consistency-screen formula (`R/consistency_resample.R:10–21`, header, verbatim):

```
# replaced by a SINGLE subgroup fit plus a multiplier representation of the
# split pair. ... the random-split half effects are
#       beta_hat +/- D,   D = sum_i G_i * dfbeta[i, treat],
# ... A split is consistent iff both
# halves clear the threshold on the COMPARISON scale c:
#       rate ~ 2 * Phi((beta_hat - c) / sigma_D) - 1,
```

Winner rule (`R/subgroup_consistency_helpers.R:570–574`):

```r
  # maxeffCons: the effect maximiser among CONSISTENCY-QUALIFYING candidates.
  ...
  if (sg_focus == "maxeffCons") {
    data.table::setorder(result_new, -hr, K)
    return(result_new)
  }
```

The v5-named pieces verified: `consistency_method = c("resample", "split")` (formal, C1); `remove_near_duplicate_subgroups()` and `.maxeff_membership_dedup()` (step 9 — note the maxeff dedup applies to `"maxeff"` only, **not** to `"maxeffCons"`, which takes the statistics-keyed near-duplicate path); `setorder(-HR, K)` (steps 7/10/12); `sg_focus` (step 12). No `method = "closed"` value exists in source — the closed-form computation is what `consistency_method = "resample"` *does* (the header above); the literal argument values are `"resample"` and `"split"`.

### C3. What the returned object records about the realized family

`forestsearch()` return (`R/forestsearch_main.R:3348–3387`), relevant fields:

```r
  out <- list(
    grp.consistency = grp.consistency,
    find.grps = find.grps,
    confounders.candidate = FSconfounders.name,
    confounders.evaluated = confs_labels,
    ...
    sg.harm = sg.harm,
    mr_inference = mr_out,
    ...)
```

- **The candidate table**: `find.grps$out.found$hr.subgroups` — one row per screened candidate with columns `grp, K, n, E, d1, m1, m0, HR, L(HR), U(HR)` plus one 0/1 indicator column per cut factor (`colnames(Z)`), i.e. **the realized rule definitions in indicator form** (`format_search_results`, `subgroup_search.R:876–881`). This is the pre-dedup, post-screen family; its `nrow` is the realized screened-family size.
- **The consistency table**: `grp.consistency$out_sg` (built by `sg_consistency_out()`; only candidates that passed consistency have rows in `out_sg$result` — in-source comment at `subgroup_consistency_main.R:1114–1116`), plus the explicit counts `grp.consistency$n_candidates_total`, `n_candidates_evaluated`, `n_passed`, `early_stop_triggered`, `early_stop_candidate` (return list, L1135–1148).
- **A third count**: `mr_inference$n_family = length(asm$names)` and `n_selected = sz[sel]` (`R/fs_mr_inference.R:603`) — the size of MR's **re-enumerated** family (all ≤ maxk combinations at the `n.min` floor only; a stated superset of the search's family, `forestsearch_main.R:3260–3272`).

**Plainly**: yes — a caller holding a normal `forestsearch()` return can read the realized family size (three variants of it, each with a stated meaning) and the realized rules' definitions (indicator columns of `hr.subgroups`, plus `sg.harm` as cut-expression strings for the winner). What is **not** kept anywhere is the per-replicate persistence of any of this: the drivers discard `fs.est` after reading `sg.harm`, `sg.harm.id`, and `mr_inference` scalars (B2).

### C4. The search's size constraints beside the document's prevalence filter (S1)

| | document (S1) | search (C2 step 5) |
|---|---|---|
| quantity | `Pg >= 0.12` — **population** probability on the constructed cell law (df_super-style) | `sum(id.x) > n.min` — **subject count of the replicate's analysis sample** at n (plus `colMeans(x) >= minp` per constituent factor, plus `rmin` marginal-shrink redundancy) |
| bound type | lower only | lower only (no upper bound exists in either) |
| same quantity? | **No** — population prevalence vs realized sample proportion; the document's `0.12 = 60/500` reproduces `n.min/n` in expectation but the search's floor binds on the sampled count, so membership of near-floor rules varies across replicates (data-dependent, C2), while the document's family is fixed. The document's S1 has no analogue of `minp` (per-factor) or `rmin` (redundancy). |

### C5. Is candidate enumeration callable without running the search?

Partially, from three exported pieces — but **no single function** turns a data frame plus the `forestsearch()` arguments into the enumerated rule set before screening:

- `get_FSdata()` (exported, `NAMESPACE:165`) — data frame + cut arguments → the dichotomized candidate-factor frame (the cut set).
- `get_combinations_info()` (exported, `NAMESPACE:166`; `subgroup_search.R:394`) — L, maxk → the combination index structure.
- `subgroup.search()` (exported, `NAMESPACE:235`) — but it **fuses enumeration with the screen and the per-subgroup model fits**; there is no argument that stops it after enumeration (`disable_effect_floor` disables only the effect screen, not the fits or the other filters).

The internal pieces that would compose an enumeration-only path (`generate_combination_indices()`, `get_covs_in()`, `get_subgroup_membership()`) are `@keywords internal`, not exported — though `forestsearch_main.R:3264–3277` (the MR family re-enumeration) demonstrates exactly that composition in eleven lines. Verdict: **none found** as a single callable; the composition exists in-source at the MR call site.

---

## 5. Open items

**CANNOT DETERMINE:**
1. (A3 S14) Which payload *columns* the non-detection `m500` literals and the `meas`/`measd` literals were tabulated from: the in-code comments cite the rendered reports' sections, not payload columns; only `m500$detection` is attributed to the payload directly. (This is the standing m500-conversion residue item of v5 §6 — read, not converted, per task scope.)
2. (B5 row 7) Below-floor selection mass has no measured counterpart column; parity untestable as stored.

**Defects noticed (described, not fixed):**
1. **Null-cell payload carries a stray `NA`-named `NULL` element** (`names(p)` = `"results" "truth" NA "oc" "meta"`): it was written by the pre-`d884adbf` save site `saveRDS(.payload[c("results","truth","scale","oc","meta")], out)` while `.payload` had no `scale`; indexing by an absent name yields the `<NA> = NULL` slot. The committed fix (`d884adbf`, `intersect(...)`) post-dates the payload commit (`2b180813`, 2026-08-26 23:44 vs 2026-08-27 10:45). Consequence: cosmetic for readers using `$`-access; a consumer iterating `names()` hits an `NA` name. Possible dispositions (no action taken): re-save the bundle through the fixed save site, or leave as a recorded artifact.
2. **Task-list naming vs code**: the "top-7 selection frequencies" of A3 S11 are top **3** in code (`top7 <- order(-sel)[1:3]`, DOC L2367) — the variable name overstates the printout. No document change made or proposed beyond noting it.
3. **`n_sel` vs `n_harm` redundancy**: the driver records both, and in all three payloads they are identical by construction (both count the selected subgroup's subjects, one via MR's `n_selected`, one via `sg.harm.id`) except where MR is absent. Not a defect, but a column pair a schema consumer could mistake for family-vs-subgroup counts; recorded here because B2's question turns exactly on that misreading.
4. **`pcons_gate` is assigned but unused in DOC** (`worked-scenario` L1737 defines `pcons_gate <- 0.90`; no later chunk reads it — confirmed by the index's dependency table, which lists no consumer). Dead assignment only; the prose carries the 0.90 discussion.

**Skipped, and why:**
- The B4 diagnostic: no per-replicate family-size field exists (B2/B3), so per the task nothing was computed.
- The `itt-live` payload metadata read (§3.2): the file is absent on this machine (machine-local, untracked); recorded as such.
- The `…n500_s100_d5000` payload: not read beyond existence/tracking — no part of the task consumes it.

---

## 6. Verdict

1. The document computes its operating characteristics on three **fixed, hand-constructed families** (M = 3/5 typed; M = 16 and M7 = 18 enumerated over a constructed 18-cell law), never on data; Q is a member only of the typed warm-up family.
2. Every downstream stage consumes the family only through one interface: labels, prevalences, purities, means, ses, overlap matrix, M.
3. The gate machinery (single-split declaration, maxeffCons argmax, selection-weighted functionals) is expressed entirely in that interface plus (n, c1, c2, draws) — it runs unchanged on any family; the family construction and every measured comparator are fixture-typed.
4. The means and ses are MD-only (mixture identity, `fs_dgm_scale` MD branch); the covariance, gate, and selection layers are scale-free given them.
5. The search builds its family per replicate from empirical quantile cuts, sample-count floors (`n > n.min`, `minp` per factor, `rmin`), a sample-fit effect screen, statistics-keyed dedup, and a closed-form resample consistency screen at 0.90/400 — every stage data-dependent at fixed n.
6. The realized family size is **computed and returned** by `forestsearch()` in three variants (`nrow(find.grps$out.found$hr.subgroups)`, `grp.consistency$n_candidates_*`, `mr_inference$n_family`) but **no driver records any of them per replicate**, and no payload carries one.
7. The payloads do carry the winner's rule verbatim (`sg_def`) per replicate, so selection distributions are comparable rule-for-rule; the rest of the realized family is discarded at the recorder.
8. Definition parity: detection share, E|Ĥ|, E[β(Ĥ)], and naive bias score the same events (population vs sample evaluation); declaration differs by the stated single-split-vs-400-split idealization; classification differs population-vs-sample with a signed, document-acknowledged gap; below-floor mass has no measured counterpart.
9. The document's prevalence floor and the search's `n.min` are different quantities (population probability vs sampled subject count), both lower-only, with no upper bound on either side.
10. Enumeration without the search is not callable as one function; its pieces exist (two exported, three internal) and the composition is demonstrated inside the MR block of `forestsearch_main.R`.
