---
bibliography: []
---

# Handoff: forestsearch code/theory alignment

For a fresh chat. Written 2026-08-05, after a long session that fixed seven
defects and produced a systematic audit.

Repo `~/Documents/GitHub/forestsearch`, branch `mr-in-replicates`. The branch
was pushed on 2026-08-05; see `HANDOFF_2026-08-05_session2.md` for what has
landed since this document was written.

---

## The through-line

Every defect found in this project has had one shape: **a quantity the
manuscript defines once, reconstructed or filtered at a second site.** That is
the lens. It is not a general bug hunt.

The manuscript is `fs-post-selection` — "Post-Selection Inference for the
Standard Trial Analysis of Data-Driven Subgroups". Its definitions, quoted
verbatim, are in `dev/identifier-alignment/CODE_THEORY_AUDIT_SPEC.md` — use
that rather than paraphrasing from memory. The PDFs are **not** in the repo.

---

## State

**Committed** (thirteen commits, none pushed): identifier alignment defaults and
`.validate_mr_configuration()`; the two application re-analyses; three
simulation documents; four dev notes; output bundles; the `maxeff`/`maxeffCons`
frontier fix; the MR admission-set resolver; the GRF inclusion band; the MR
`.inband()` consolidation; the `dfbeta_glm()` swap; the MR family-restriction
removal; the F1 `lm()` fix.

**Uncommitted, deliberately** — `max_subgroups_search` default set to `Inf` in
`R/forestsearch_main.R`, its `NEWS.md` entry, the regenerated
`man/forestsearch.Rd`, and both analysis documents with their four
`max_subgroups_search` caps removed. These ride along with the next commit that
runs a check.

**Also in the tree**: `CODE_THEORY_AUDIT_SPEC.md`, `code_theory_audit.qmd` and
its rendered HTML, and `actg175_realignment_report.qmd` (Larry's).

**Held**: regeneration of all five analysis documents. Do not run it until the
audit findings are settled — it has already gone stale twice.

---

## Audit findings, and where each stands

From `dev/identifier-alignment/code_theory_audit.qmd`. Fourteen findings, eight
questions.

| id | status |
|---|---|
| **F1** — `.dfbeta_glm()` errored on every unweighted `lm()` fit | **fixed**, committed `411a448` |
| **F2** — GRF identifier admits `HR >= 1` (hardcoded `dmin = 1`) while MR admits `log(hr.threshold)` | **in flight** |
| **F3** — DINA identifier applies no `beta-hat` floor while MR applies `log(hr.threshold)` | **in flight** |
| **F4** — under `ps_method != "none"`, identifier ranks on the IPTW effect while `.consistency_glm_pieces()` fits the plain conditional model | open |
| **F9** — two `sg_focus` dispatch sites the completeness assert does not read | open, low |
| **F10** — `.mr_cscr`, `c_screen_mr`, `c_consistency_mr` assigned and never consulted | open, low |
| **F11** — `max.minutes` is a formal no path consults | open, low |
| **F12** — cross-validation omits the `grf_res`/`grf_cuts`/`dina_res`/`dina_cuts` nulling the bootstrap performs, so a cached fit leaks the held-out fold into the training family | open, **medium** |
| **F13** — `bias_sel` averages over draws with a winner, `bias_fix` over all draws; Eq. 9 divides both by B | open, low |
| **F14** — both replay paths re-derive quantile cuts per replicate, so ℱ is not resample-invariant | noted, not a defect |

**The F2/F3 answer has been sent to CC** and is the manuscript's, not a
judgement call. §2.4 admits *"among the admitted {g : β̂(g) ≥ t_g}"*; §1.1's
common form for all three identifiers *"retains those clearing a harm
threshold"*; §6's DR-score and mean-effect floors are the identifiers'
**native-scale** qualification, applied before re-ranking, not a substitute for
$t_g$. So 𝒮 admits on $\hat\beta \ge t_g$ for every identifier, and for DINA and
GRF there is no consistency term, so $t_g = c_{\mathrm{screen}}$. **MR is
correct; both identifiers are wrong; the fix moves the identifiers to MR.**

---

## Settled — do not relitigate

- All six `sg_focus` × engine configurations are **retained**; the three that
  fail the alignment condition are opt-in only and hard-error under MR.
- **No warnings** for DINA/GRF on the fixed-family condition. Their families
  regenerate by construction and that is the manuscript's subject matter, not a
  defect. At sweep volume a warning would drown everything else.
- **MR is not required to mimic the full bootstrap.** This killed
  `.fs_mr_restrict_native()`, whose whole rationale was mirroring FB's
  neighbourhood. MR stands on its own and must not be limited to make the two
  agree.
- **Freezing the cut grid is not to be implemented.** FB refits everything per
  resample by design and that is correct. Whether a run has a fixed family is a
  property of the configuration; `conf_force` with literal cut expressions is
  the route if anyone wants it.
- `family_status = "no-front-end"`, not `"fixed"` — the manuscript's §2.1 fixed
  family additionally requires resample-invariant cut locations.
- **No epsilon guard** in `.dfbeta_glm()` at leverage 1. The LOO estimate does
  not exist there; `NaN` is honest and the caller drops the candidate.
- Manuscript revision is **out of scope**. Everything gets regenerated; do not
  track corrections against current tables.

---

## Queued, in order

1. **F2/F3** — in flight with CC.
2. **F12** — the CV nulling gap. Needs a decision on whether cached fits should
   be reusable across folds at all.
3. **The MD closed-form fixture** — see below. Not yet specified.
4. **Regenerate all five documents.** Expect large movement: GRF's MR family
   goes from 2 candidates to 1289, and the FS labels may shift in both
   applications from the uncapping.
5. **`selection_criteria.qmd` additions** — Larry's spec document. Needs: a
   scope statement that its key table is `subgroup_method = "consistency"` only;
   the DINA and GRF keys, which contain no `Pcons` term; that `maxeff` and
   `maxeffCons` are synonyms outside FS; that MR reuses `maxeff`/`maxcons` for
   different rules; the band-per-engine position after the GRF fix; and the
   `conf_force` freeze recipe.
6. **Process-discipline drafts** — two documents exist and are unreviewed, with
   five open questions in `review_note_process_discipline.md`.

---

## The MD closed-form fixture (not yet written up)

Larry's derivation document `postselection_theory_derivations.qmd` already gives
`db` and $\sigma^2_D$ in closed form for MD and log-OR (§sec-worked), explicitly
as implementation fixtures. Three further results were derived and numerically
verified in the last session and are **not** in that document:

**With Gaussian multipliers**, $D_g^{(b)} \sim N(0, \sigma^2_{D,g})$ exactly, and
$\mathrm{Cov}(D_g, D_h) = \sum_i \mathrm{db}_{g,i}\mathrm{db}_{h,i}$.

**For two candidates with no admission floor**, with $d = \hat\beta_1 -
\hat\beta_2$ and $\tau^2 = \sigma_1^2 + \sigma_2^2 - 2\sigma_{12}$:

$$\mathrm{bias}_{\mathrm{sel}} = \tau\,\varphi(d/\tau)$$

Verified against 400,000 Monte Carlo draws across five configurations including
correlated influence vectors. At $d = 0$ with independent unit variances it
returns $1/\sqrt{\pi} = 0.5642 = E[\max$ of two standard normals$]$.

**$\mathrm{bias}_{\mathrm{fix}} = 0$** exactly under centered multipliers.

This matters because the selection-bias term is the only quantity in the
pipeline with no independent check — everything else can be verified against a
refit. Under MD the score is linear, so Lemma 1's curvature term vanishes
identically and the whole pipeline is analytic.

An MD fixture would have caught F1 immediately: the MD path errors *entirely*,
so asserting $\sigma^2_D = \hat s_1^2/n_1 + \hat s_0^2/n_0$ fails on contact.

Open: the $M = 1$ degenerate IJ variance. A sketch suggested $\tilde V \to
4\sigma^2_D$, which was surprising enough to want derived properly rather than
asserted.

---

## Method — earned, not theorised

Each of these was learned by violating it.

**Trace call sites, never formals or documentation.**
`family_native_neighborhood` defaulted to `NULL`, and its own roxygen said
`NULL` disables the restriction — but both call sites substituted
`effect_neighborhood` first, so it ran at 0.10 in every DINA and GRF analysis
ever performed. A function's name and rationale may describe behaviour it no
longer has.

**A leaf check does not establish the wiring.** 30,000 numerical comparisons
against `.grf_frontier_select()` passed while the argument under test never
reached it, because two intermediate argument-builders lacked the formal.

**Enumerate exhaustively.** A `grep | head -20` reported four MR enforcement
sites when there were five.

**A control that passes because the quantity is easy is not evidence the code
path is reachable.** F1's gaussian control was `glm(family = gaussian())`, where
`fit$weights` is a vector of ones; the MD path uses `lm()`, where it is `NULL`.
The control exercised the gaussian *family* and never the `lm` *path*. This is
CC's formulation and it is the best statement of the failure anyone produced.

**Before any re-run, enumerate what every artifact writes.** Four separate times
a document's output filename collided with the reference it was being compared
against. Any task shaped as "re-run X and compare to what X produced before"
carries this by construction.

**Run FB only where the comparison requires it.** MR-versus-MR checks need no
bootstrap. Where FB is needed to detect whether a gap closed rather than to
estimate it, 100 replicates suffices.

---

## Corrections made in the last session — do not repeat

I asserted several things that turned out to be wrong. Recording them so the
pattern is visible rather than the specifics.

- Called `effect_neighborhood` a "tie-break" when it is a deliberate
  effect-versus-size trade — a candidate with a slightly weaker effect but
  substantially larger *n* can win.
- Said winner's curse "should be minimal" with a small family. For $M = 2$ it is
  $0.56\sigma$, which is not minimal.
- Reported a duplicated `X <- stats::model.matrix(fit)` line that did not exist
   — read a diff as a state.
- Framed the FB/MR gap as something to close, when MR is not required to mimic
  FB.
- Supplied the `dfbeta_glm()` implementation that caused F1, and validated its
  gaussian control on the wrong object.

The common failure is asserting from a fragment rather than reading the source.
Larry has corrected this repeatedly and is right to.

---

## Working conventions

- Every response creating or modifying a file ends with `present_files`.
- Recommend-then-implement: present options, "proceed" triggers application.
- CC writes the code; the chat writes specifications. Code in a handoff
  duplicates what the executor is about to build.
- Never infer file contents. Ask for the file.
- Deliverables `.qmd` with PDF target, `scrartcl`; `\mathrm`, LaTeX in `$...$`.
- Fixed abbreviations: FS, DINA, GRF, FB, MR. β(Ĥ) conditional versus θ†(H)
  marginal — never conflate.
- Testing means integration-style Quarto living documents, not testthat
  scaffolding — with the standing exception of contract tests for static
  properties, where an assertion nothing invokes would be decorative.
- `R CMD check` clean is the bar for a commit; `devtools::test()` alone is not.
- Larry's machine is Pop!_OS, ~127 cores, R 4.6.1. The chat sandbox has R 4.3.3
  and no Quarto.
