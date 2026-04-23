# Next session: Investigate bootstrap 0/20 identification rate

**Purpose of this document:** A self-contained starter for a fresh
Claude session to investigate why two smoke-test Quarto documents show
0 / 20 bootstrap identification rates despite primary `forestsearch()`
and cross-validation both succeeding.

**Context:** `forestsearch` R package, v0.2.0-dev branch,
`feature/glm-extension`. CRAN-published v0.1.0 in March 2026.


## How to use this document

1. Open a fresh chat with Claude in your `forestsearch` project.
2. Paste the "Opening message" block below as your first message.
3. Have the "Artifacts to attach" ready in case Claude asks for them.
4. Follow the suggested workflow. Redirect Claude if it drifts into
   testthat territory or re-audits files already fixed.


## Opening message (copy-paste this)

> I want to investigate a bootstrap identification issue in
> `forestsearch`. Two smoke-test Quarto documents
> (`smoke_count_glm.qmd`, `smoke_continuous_glm.qmd`) show the same
> unexpected pattern:
>
> - Primary `forestsearch()` call: identifies subgroup
>   `{age = 62}` ✓
> - CV via `forestsearch_tenfold()`: 50 / 50 identifying folds ✓
> - Bootstrap via `forestsearch_bootstrap_dofuture()`:
>   **0 / 20 iterations identify a subgroup** ✗
>
> Same synthetic DGM (harm signal at
> `age <= 50 & biomarker_hi == 1`), same primary config. The
> survival bootstrap smoke test (`test_bootstrap_grf_cuts.qmd`,
> GBSG) works fine at 31 / 31, so this appears specific to the
> GLM outcome paths in the count and continuous smoke tests.
>
> Before jumping to fixes, I want a structured diagnostic pass.
> Please start by **auditing without modifying**:
>
> 1. Read `R/bootstrap_analysis_dofuture.R`,
>    `R/bootstrap_dofuture_main.R`, and
>    `R/bootstrap_calculations_helpers.R` to understand how
>    bootstrap iterations call `forestsearch()` and how GLM
>    parameters are propagated.
> 2. Identify what could cause a config that works on CV (which
>    uses held-out folds) to fail on bootstrap resamples. Likely
>    suspects: parameter propagation through `fs.est$args_call_all`,
>    scale / threshold mismatches specific to bootstrap-scale noise
>    on small synthetic DGMs, seed handling.
> 3. Give me a severity-organized finding list
>    (Critical → High → Moderate → Minor) before proposing any
>    fix.
>
> Do not write testthat tests. Smoke-test Quarto docs are the
> testing mechanism.


## Why this framing

Three things in the opening message are deliberate:

1. **"Auditing without modifying"** invokes the preferred discipline
   (severity-organized cross-file audit before any fixes). Without
   that instruction, Claude can be tempted to jump to fixes after
   spotting the first plausible cause.

2. **"Do not write testthat tests"** is an explicit guardrail
   against a known rabbit hole. The userMemories memory entry
   should trigger this anyway, but stating it upfront in session 1
   reinforces it.

3. **The specific symptom AND non-symptom together** — telling
   Claude that CV works (50 / 50) and that survival bootstrap
   works (31 / 31) constrains the hypothesis space to "something
   specific to the GLM bootstrap path." Without that context,
   time could be wasted on hypotheses that CV would already have
   exposed.


## Artifacts to have ready

You may be asked for any of these:

- **HTML renders from this session** (evidence of the symptom):
  - `smoke_count_glm.html`
  - `smoke_continuous_glm.html`
- **Qmd sources** (the exact configs that reproduce):
  - `smoke_count_glm.qmd`
  - `smoke_continuous_glm.qmd`
- **Per-iteration results data.frame** (if Claude asks).
  Run the smoke test, then `head(fs_bc$results, 20)` —
  the `sg1_b`, `sg2_b`, `grf_cuts_b` columns will tell
  you whether bootstrap iterations are returning "subgroup
  found but rejected" or "no subgroup candidate at all."


## What NOT to share at session start

- The full transcript summary from this session. The new session
  should approach the bootstrap issue fresh, not carry forward
  Phase A / B / C context that is not directly relevant. The
  userMemories will give it enough project background; the
  transcript would add noise.
- Memory edits from the last session unless Claude specifically
  asks for them.


## What to redirect if it happens

- If Claude tries to re-audit `R/get_fsdata.R`, redirect: that
  work is done and committed.
- If Claude proposes rebuilding smoke tests, redirect: smoke
  tests exist and render cleanly. The issue is a finding
  revealed by them, not a test infrastructure problem.
- If Claude proposes testthat tests, redirect: that is the
  rabbit hole memory #4 flags.


## Likely hypotheses (what to watch for)

These are guesses, not conclusions. The new session should
confirm or reject each through code audit:

1. **`fs.est$args_call_all` dropping an outcome-type-relevant
   parameter** (`effect_measure`, `offset.name`, or
   `outcome_type` itself) when bootstrap copies the config.
   Bootstrap workers might be getting a survival-defaulted
   config on a count / continuous DGM.

2. **Bootstrap seed handling producing degenerate resamples**
   on small-N synthetic DGMs where the harm signal is carried
   by a small fraction of observations. With N = 800 and ~23%
   harm prevalence in the smoke test DGMs, a bootstrap resample
   might routinely dilute the signal below threshold.

3. **GRF on bootstrap resamples returning no cuts** for GLM
   outcome types specifically, leading to `use_grf <- FALSE`
   being silently triggered and the pipeline falling through
   to no candidates.

4. **Parallel worker scoping** — bootstrap uses a different
   parallel path than CV (`foreach` + `doFuture` vs
   `future_lapply`). Config might not survive the handoff to
   workers cleanly.

A careful audit should pick which of these is actually
happening and whether it is a genuine bug or a smoke-test-DGM
issue. Either outcome is useful to know.


## Suggested sequence within the session

1. Audit the three bootstrap source files and give the
   severity-organized finding list.
2. Agree on the top one or two hypotheses to test.
3. Add minimal instrumentation to `R/bootstrap_analysis_dofuture.R`
   (probably a few lines of `message()` calls that print the
   config the worker received and whether GRF returned any
   cuts) — scoped to a debug flag that can be off by default.
4. Re-render `smoke_count_glm.qmd` with the instrumentation
   on. Review the output together.
5. Decide whether the finding is a real bug (fix it) or a
   smoke-test DGM issue (adjust the DGM or document the
   limitation).
6. If it was a real bug: remove instrumentation, commit
   the fix, re-render to confirm.


## Definition of done

- Root cause of the 0 / 20 symptom identified and documented
  (even if the conclusion is "DGM is under-powered for
  bootstrap resampling; increase N in the smoke test").
- If it is a code bug: fix applied, smoke tests re-rendered
  green, commit staged.
- If it is a DGM issue: smoke test DGM adjusted OR a
  comment added to the qmd documenting the expected behavior.


## Known good baseline

As of the end of the previous session (April 23, 2026):

- `R/get_fsdata.R` robustness hardening committed.
- All smoke tests render. 164 / 164 assertions pass across
  the Phase A / B / C test docs. No regressions from the
  robustness work.
- Canonical simulation documents
  (`gbsg_poisson_simulations.qmd`,
  `actg175_continuous_simulations.qmd`) in the main suite
  do NOT have the 0 / 20 bootstrap symptom, so it is specific
  to the smoke-test DGMs' small-N configurations, not a
  pipeline regression.
