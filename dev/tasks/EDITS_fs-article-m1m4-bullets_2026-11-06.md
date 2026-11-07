# EDITS_fs-article-m1m4-bullets_2026-11-06 — companion to the task doc

One old → new replacement in `vignettes/articles/mrct-consistency.qmd`, §1.2
"The data-generating mechanism". Match `old` with whitespace-normalized exact matching (one
occurrence); write `new` with the list structure as given, prose re-wrapped at ~100 characters.
The paragraph that follows ("Writing the effect with treatment-by-region interactions …") is
untouched.

**old**

The mechanism therefore assumes, inside and outside the region alike: **(M1)** proportional hazards within each patient, with patient-specific ratio $e^{\psi(X)}$; **(M2)** a common effect law across the region boundary — $\psi$ does not involve $R$, so $(b, h)$ are the same inside and outside the region; **(M3)** a purely prognostic region — $r$ shifts $\theta_0$ for all region patients and neither arm is favored; and **(M4)** free regional composition — the distribution of $X$, and with it the region's $H^*$ share, may differ between the region and the non-region population. Under (M1)–(M4) any divergence of the region's windowed Cox limit from the trial's is produced by composition and by the finite analysis window (the clock identity of the simulations article), not by any region-dependent treatment effect.

**new**

The mechanism therefore assumes, inside and outside the region alike:

- **(M1) Proportional hazards within each patient**, with patient-specific ratio $e^{\psi(X)}$.
- **(M2) A common effect law across the region boundary** — $\psi$ does not involve $R$, so $(b, h)$ are the same inside and outside the region. The extension below relaxes this: with treatment-by-region interactions $\gamma_R$ and $\gamma_{H^*\!R}$, $\psi$ may depend on $R$, the effects inside the region become $(b_R, h_R)$, and (M2) is then the testable null $\gamma_R = \gamma_{H^*\!R} = 0$.
- **(M3) A purely prognostic region** — $r$ shifts $\theta_0$ for all region patients and neither arm is favored.
- **(M4) Free regional composition** — the distribution of $X$, and with it the region's $H^*$ share, may differ between the region and the non-region population.

Under (M1)–(M4) any divergence of the region's windowed Cox limit from the trial's is produced by composition and by the finite analysis window (the clock identity of the simulations article), not by any region-dependent treatment effect.

(The γ notation in the new (M2) sentence — `$\gamma_{H^*\!R}$` — matches the following paragraph's
existing notation exactly.)
