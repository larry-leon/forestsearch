---
title: "F13 — two open questions, recorded not chased"
date: 2026-08-05
bibliography: []
---

# Context

F13 fixed an incoherence between the two multiplier-resampling bias terms'
denominators: `selection_bias` averaged over the draws that produced a
re-selected winner, `fixed_bias` over all `B`. Both are now taken over the
draws on which a selection occurred — **convention (iv), conditional on
identification** — rather than both over `B` with a zero substituted on
no-winner draws (**(iii), Eq. 9 read literally**).

Both conventions center the residual to machine precision, so centering does
not choose between them; they differ in **estimand**. The reasoning for (iv) is
recorded in `SPEC_f13_reverse_to_conditional.md` and is not re-opened here.

Two questions surfaced during that work. Both are recorded so they are not
lost, and neither was chased.

---

# 1. What does the full bootstrap record for a replicate that identifies nothing?

**Why it matters.** MR is a linearization of the full bootstrap (FB). If FB
**drops** a replicate that identifies no subgroup, then FB is itself conditional
on identification, and (iv) aligns MR with it — which would be an independent
argument for the convention rather than a restatement of it. If FB
**substitutes** something instead, that is worth knowing, because MR and FB
would then condition differently and the two would diverge exactly on the
no-winner draws.

**Status: unverified.** `bootstrap_analysis_dofuture.R` was not read for this
question. The relevant behaviour is what happens when a replicate's
`forestsearch()` call returns `sg.harm = NULL` — whether the replicate
contributes to the bias average, is dropped from it, or contributes a
substituted value.

**Cheap to settle**, and worth settling before the manuscript's Theorem 1
discussion is finalised, since that theorem has MR reproducing FB on a fixed
family.

---

# 2. Is the manuscript silent on the no-winner case, or contrary to (iv)?

Eq. 9 divides both bias terms by `B`. That form **appears to assume every draw
admits a candidate** — under which assumption (iii) and (iv) coincide exactly
and the question does not arise.

If that reading is right, the manuscript is **silent** on the case where a draw
admits nobody, not contrary to (iv), and no amendment is implied by this
change.

**Status: unverified, and flagged as such deliberately.** Settling it needs a
read of §3.1 to see whether the no-winner case is addressed anywhere, which was
out of scope. It should not be asserted either way in the manuscript file until
that read happens.

A related point that **is** settled: Eq. 13 needs no amendment. Its uncentered
second moment is Wager–Hastie–Efron's centered quantity once the two bias terms
share a denominator, under either convention. The apparent gap was the
implementation, not the theory.

---

# Supporting evidence produced for this change

In `dev/identifier-alignment/verification/`:

- `verify_eq9_alignment.R` (V8) — the three-arm comparison establishing that
  (i) is not centered, and that (iii) is. Written when (iii) was the intended
  target; retained because the diagnosis it records is unchanged.
- `verify_conditional_convention.R` (V9) — adds (iv) and shows it is equally
  centered (`mean(r)` = 9.7e−18), leaves `tilde_V` unchanged (+0.0%), and leaves
  `selection_bias` unchanged, where (iii) multiplies it by `selection_rate`.

## The conditional `fixed_bias` behaves like a selection effect

Checked because under (iv) `fixed_bias` is a conditional expectation and is no
longer mean-zero — at `selection_rate` 0.865 it was 0.0909 against 0.0053
unconditional, and it was worth confirming that this is selection rather than
an artefact. Synthetic reference, `B` = 40,000, sd(D_Ĥ) = 0.992:

| `selection_rate` | conditional `fixed_bias` | ÷ sd(D_Ĥ) |
|---:|---:|---:|
| 0.756 | 0.2039 | 0.205 |
| 0.796 | 0.1739 | 0.175 |
| 0.822 | 0.1532 | 0.154 |
| 0.850 | 0.1311 | 0.132 |
| 0.881 | 0.1064 | 0.107 |
| 0.919 | 0.0754 | 0.076 |
| **1.000** | **−0.0098** | −0.010 |

Strictly monotone in the selection rate, and at rate 1 it equals the
unconditional value **exactly**. The unconditional value, −0.0098, is about two
Monte-Carlo standard errors from zero at this `B` and is consistent with the
mean-zero multipliers. The conditional value is roughly proportional to the
excluded fraction, ≈ 0.85–0.93 × (1 − rate).

That is the signature of conditioning on an event correlated with the quantity
averaged: draws clearing the admission floor are exactly those where the
perturbation pushed effects up. It vanishes as the conditioning event becomes
vacuous, which an artefact would not do.
