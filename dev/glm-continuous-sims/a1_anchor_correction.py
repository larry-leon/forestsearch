#!/usr/bin/env python3
"""
a1_anchor_correction.py -- correct the analytic verification document's scale.

WHAT IS WRONG
    The document anchors its scale on a MEASURED sigma[betahat_1000(Q)] =
    13.6786, implying an effective bracket of 16,119.  The fixture's residual
    variance is 16,256.28, and on the true region the bracket is
    sigma^2 + V_Q[mu_0] with V_Q[mu_0] >= 0, so the anchored value sat BELOW
    THE NOISE FLOOR -- structurally impossible, and checkable from sigma alone.

    Direct enumeration gives bracket(Q) = 16,863.44, hence V_eff = 67,454 and
    sigma[betahat_1000(Q)] = 13.9909.  That value reproduces the measured mean
    fitted standard error to 0.15% at n = 500 and 0.03% at n = 700.

    Two further prose claims are wrong and are rewritten, not renumbered:
      - the Jensen magnitudes (1.15% / 0.29%; exact values 0.295% / 0.192%),
        and the claim that they cover "roughly half the discrepancy";
      - the CATE negligibility argument, which bounds V_g[tau]/2 <=
        beta_int^2/8 = 23.6 but never mentions C_g[mu_0,tau], which is
        UNBOUNDED and equals -33.59 on S -- larger in magnitude than the bound
        on the term that is discussed.

    The document's FORMULAS are already correct (L396 carries the full
    sigma^2 + V_g[mu_0] + V_g[tau]/2 + C_g[mu_0,tau]).  Only the calibration
    and these claims change.

STRUCTURE
    Part 1: six chunk knobs, 13.6786 -> 13.9909.  Pure substitution; the chunk
            output then matches /tmp/a1/new.out exactly.
    Part 2: nine anchored prose edits.  Four are renumberings, four reframe the
            anchor from measured to computed, one is a full paragraph rewrite.

SAFETY
    Dry run by default.  Every edit gated on its anchor matching byte for byte
    exactly once; a single mismatch aborts before anything is written.  The
    write is the last operation.

USAGE
    python3 a1_anchor_correction.py <document.qmd> [--apply]
"""

import sys

OLD_KNOB, NEW_KNOB = "13.6786", "13.9909"
N_KNOBS = 6

EDITS = []


def edit(name, old, new):
    EDITS.append((name, old, new))


# -- 1. anchor definition: measured -> computed --------------------------------
edit("anchor-definition",
"""by the law above, would pin the same constant. Concretely, the committed
measurement is $\\sigma[\\hat\\beta_{1000}(Q)] = 13.679$, which fixes
$V_{\\mathrm{eff}} = 1000\\,\\pi_Q\\,\\sigma[\\hat\\beta_{1000}(Q)]^2
= 64{,}476$ and hence
$\\sigma[\\hat\\beta_{1000}(\\mathcal S)] =
\\sqrt{V_{\\mathrm{eff}}/1000} = 8.030$, agreeing with the prevalence
rule's $13.679\\sqrt{0.3446}$. **No relationship between $g$ and $Q$ is
being asserted**: the region at which the measurement happens to sit is
""",
"""by the law above, would pin the same constant. The constant is not
measured but *computed*: enumerating $\\sigma^2 + V_Q[\\mu_0]$ on the stored
super-population gives an effective bracket of $16{,}863.44$, hence
$V_{\\mathrm{eff}} = 67{,}454$ and
$\\sigma[\\hat\\beta_{1000}(Q)] =
\\sqrt{V_{\\mathrm{eff}}/(1000\\,\\pi_Q)} = 13.991$, and hence
$\\sigma[\\hat\\beta_{1000}(\\mathcal S)] =
\\sqrt{V_{\\mathrm{eff}}/1000} = 8.213$, agreeing with the prevalence
rule's $13.991\\sqrt{0.3446}$. **No relationship between $g$ and $Q$ is
being asserted**: the region at which the constant is evaluated is
""")

# -- 2. symbol table -----------------------------------------------------------
edit("symbol-table",
"| $\\sigma[\\hat\\beta(g)]$ | the estimator's SD; anchor $\\sigma[\\hat\\beta_{1000}(Q)] = 13.679$ |",
"| $\\sigma[\\hat\\beta(g)]$ | the estimator's SD; computed $\\sigma[\\hat\\beta_{1000}(Q)] = 13.991$ |")

# -- 3. the check has been run -------------------------------------------------
edit("check-run",
"""This yields a falsifiable check on the anchor. The committed
$\\sigma[\\hat\\beta_{1000}(Q)] = 13.679$ implies
$V_{\\mathrm{eff}} = 1000 \\cdot \\pi_Q \\cdot \\operatorname{Var}[\\hat\\beta_{1000}(Q)] = 64{,}476$ and hence an
effective bracket $\\sigma^2 + V_Q[\\mu_0] =
16{,}119$, an effective per-subject standard deviation of about $127$
cells/$\\mu$L.
""",
"""This yields a falsifiable check on the scale, and the check has been run.
Direct enumeration on the stored super-population gives $\\sigma = 127.500$
and $V_Q[\\mu_0] = 607.158$, hence an effective bracket
$\\sigma^2 + V_Q[\\mu_0] = 16{,}863.44$,
$V_{\\mathrm{eff}} = 4 \\times 16{,}863.44 = 67{,}454$, and an effective
per-subject standard deviation of $129.9$ cells/$\\mu$L.
""")

# -- 4. the discrepancy passage: full rewrite ----------------------------------
edit("discrepancy-passage",
"""of $20.21$ and a mean fitted standard error of $19.81$, against the
anchor-implied $13.679\\sqrt2 = 19.344$: the anchor is low by $2.4\\%$
against the fitted standard error and $4.5\\%$ against the empirical
spread — the latter $1.9$ Monte-Carlo standard errors, so the two agree
at the edge of resolution. In bracket terms the effective per-subject
standard deviation is $127.0$ from the anchor against $130.0$ and $132.6$
from the two measurements. Part of the gap is identified and expected:
the balanced-arms and fixed-region-size idealizations understate the
scale, since $\\mathbb E[1/n_1 + 1/n_0]$ and $\\mathbb E[1/n_g]$ exceed
their fixed-count counterparts by Jensen's inequality, worth about
$1.15\\%$ and $0.29\\%$ in the standard deviation at this design — roughly
half the discrepancy. The consequences are inside the tolerance already
declared: adopting the fitted value would move the $Q$-channel clearance
from $0.768$ to $0.762$ and the sample-size threshold $n^\\ast$ of
@sec-floors from $9{,}538$ to about $10{,}000$. What remains open is the
direct evaluation: computing $\\sigma$ and
$V_Q[\\mu_0]$ on the stored super-population would
separate anchor measurement error from the idealization terms, which the
oracle comparison alone cannot.
""",
"""of $20.21$ and a mean fitted standard error of $19.81$, against the
computed $13.991\\sqrt2 = 19.786$. The agreement is close and, more
tellingly, scale-free: the same computation at $n = 700$ gives $16.722$
against a measured mean fitted standard error of $16.752$. Applying the
exact Jensen factor and the two concavity corrections a fitted standard
error carries, the computed bracket reproduces the measured mean fitted
standard error to $0.15\\%$ at $n = 500$ and $0.03\\%$ at $n = 700$, with
both empirical standard deviations inside one Monte-Carlo standard error.

An earlier version of this document anchored the scale on a *measured*
$\\sigma[\\hat\\beta_{1000}(Q)] = 13.679$, implying a bracket of
$16{,}119$. That value lies below the fixture's residual variance of
$16{,}256$; since the bracket on the true region is
$\\sigma^2 + V_Q[\\mu_0]$ with $V_Q[\\mu_0] \\ge 0$, it was not merely
imprecise but structurally impossible — a check available from $\\sigma$
alone. The idealizations are real but small: $\\mathbb E[1/n_1 + 1/n_0]$
and $\\mathbb E[1/n_g]$ exceed their fixed-count counterparts by Jensen's
inequality, worth $0.295\\%$ and $0.192\\%$ in the standard deviation at
this design, $0.490\\%$ jointly at $n = 500$ and $0.348\\%$ at $n = 700$.
Every rate below uses the computed constant.
""")

# -- 5. marginal recipe at n = 500 ---------------------------------------------
edit("marginal-recipe",
"""13.679\\sqrt 2 = 19.34$,
giving marginal clearance $\\Phi(10/19.34) = 0.697$ for $Q$ — against
""",
"""13.991\\sqrt 2 = 19.79$,
giving marginal clearance $\\Phi(10/19.79) = 0.693$ for $Q$ — against
""")

# -- 6. c2 cost at the fixture -------------------------------------------------
edit("c2-cost",
"""$s = \\sigma[\\hat\\beta_{1000}(Q)] = 13.679$, floors $30/10$) this costs
about two
points: screening-only $0.768$ versus $P_1 = 0.745$ — the
""",
"""$s = \\sigma[\\hat\\beta_{1000}(Q)] = 13.991$, floors $30/10$) this costs
about two
points: screening-only $0.763$ versus $P_1 = 0.738$ — the
""")

# -- 7. three-candidate family SEs ---------------------------------------------
edit("family-SEs",
"""(prevalence $0.31$, purity $1$). SEs are anchored to the measured
$\\sigma[\\hat\\beta_{1000}(Q)] = 13.679$ and scaled by the balanced-arms
""",
"""(prevalence $0.31$, purity $1$). SEs use the computed
$\\sigma[\\hat\\beta_{1000}(Q)] = 13.991$ and are scaled by the balanced-arms
""")

# -- 8. CATE bound: both terms, not one ----------------------------------------
edit("cate-bound",
"""settles the matter here: back-solving the committed anchor gives
$V_{\\mathrm{eff}} = 1000\\cdot\\pi_Q\\cdot
\\operatorname{Var}[\\hat\\beta_{1000}(Q)] = 64{,}476$, an effective
bracket of $16{,}119$ (per-subject effective standard deviation
$\\approx 127$), against a maximal CATE term of
$\\beta_{\\mathrm{int}}^2/8 = 23.6$ — $0.15\\%$ of the variance and
$0.07\\%$ of the standard deviation. The CATE contribution to the bracket
is therefore negligible on this fixture, which means the measured
""",
"""settles the matter here: direct enumeration gives
$V_{\\mathrm{eff}} = 67{,}454$, an effective
bracket of $16{,}863.44$ (per-subject effective standard deviation
$129.9$). Two CATE-dependent terms enter that bracket, and only one of
them is bounded. The variance term satisfies
$\\tfrac12 V_g[\\tau] \\le \\beta_{\\mathrm{int}}^2/8 = 23.6$, attained at
$p = 1/2$; on $\\mathcal S$ it equals $21.33$. The cross term
$C_g[\\mu_0,\\tau]$ admits no such bound — it scales with
$\\beta_{\\mathrm{int}}$ times the prognostic contrast between $Q$ and its
complement — and on $\\mathcal S$ it equals $-33.59$, larger in magnitude
than the bound on the first and opposite in sign. Both are negligible
against a bracket of $16{,}863$ ($0.13\\%$ and $0.20\\%$), and on $Q$ and
$Q^c$ both vanish identically because $\\tau$ is constant there. The CATE
contribution is therefore negligible on this fixture, which means the measured
""")

# -- 9. scenario recap ---------------------------------------------------------
edit("scenario-recap",
"""the selection rule (maxeffCons); and the scale anchor (the committed
measured $\\sigma[\\hat\\beta_{1000}(Q)] = 13.679$, which the fitted-SE
observation
of Level 3 licenses as a population-effective scale, rescaled to $n = 500$
""",
"""the selection rule (maxeffCons); and the scale (the computed
$\\sigma[\\hat\\beta_{1000}(Q)] = 13.991$, which the Level 3 fitted-SE
observation
confirms as a population-effective scale, rescaled to $n = 500$
""")


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        return 1
    path = sys.argv[1]
    apply_ = "--apply" in sys.argv
    text = open(path, encoding="utf-8").read()

    print(f"document: {path}  ({text.count(chr(10)) + 1} lines)\n")

    nk = text.count(OLD_KNOB)
    print(f"Part 1  chunk knobs {OLD_KNOB} -> {NEW_KNOB}: found {nk}, expected {N_KNOBS}"
          f"   {'OK' if nk == N_KNOBS else 'MISMATCH'}\n")

    print("Part 2  prose edits")
    print(f"  {'edit':<22} {'count':>5}  status")
    print("  " + "-" * 44)
    ok = (nk == N_KNOBS)
    for name, old, _ in EDITS:
        c = text.count(old)
        good = (c == 1)
        ok = ok and good
        print(f"  {name:<22} {c:>5}  {'READY' if good else 'ANCHOR-MISMATCH'}")

    print(f"\n{'ALL ANCHORS MATCHED' if ok else '*** ANCHORS DID NOT MATCH -- nothing will be written ***'}")
    if not apply_:
        print("\n(dry run -- nothing written; re-run with --apply)")
        return 0 if ok else 1
    if not ok:
        print("\n!! refusing to write.")
        return 1

    new = text.replace(OLD_KNOB, NEW_KNOB)
    for name, old, rep in EDITS:
        assert new.count(old) == 1, name
        new = new.replace(old, rep)

    assert OLD_KNOB not in new, "knob survived"
    assert new.count(NEW_KNOB) == N_KNOBS, "knob count wrong"
    # One 13.679 is intentional: the sentence recording what the old anchor
    # was.  Whitelisted exactly as 16{,}119 is on the following line.
    assert new.count("13.679") == 1, "unexpected 13.679 count"
    assert "64{,}476" not in new, "a prose 64,476 survived"
    assert "16{,}119" not in new.replace("$16{,}119$. That value", "X"), "unexpected 16,119"

    open(path, "w", encoding="utf-8").write(new)
    print(f"\nwritten: {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
