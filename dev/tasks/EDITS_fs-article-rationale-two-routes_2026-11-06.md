# EDITS_fs-article-rationale-two-routes_2026-11-06 — companion to the task doc

Two old → new replacements in `vignettes/articles/mrct-consistency.qmd`. Match each `old` with
whitespace-normalized exact matching (one occurrence); write `new` with its list structure, prose
re-wrapped at ~100 characters.

## E-A — §1.2, the hazard-scale bullet

**old**

- **hazard scale** — all hazards multiplied by 3;

**new**

- **hazard scale** — all hazards multiplied by 3: the Weibull fit to ACTG 175 carries that trial's own event rate, roughly 20–30% of patients within the 40-month window, and the multiplier raises the design toward a modern trial's 50–70%. It rescales event frequency only: every patient's hazard ratio — and with it $b$, $h$, and each set's AHR and CDE — is exactly as set (the windowed Cox limits do move, since the window sits later on the faster clock);

(The surrounding bullets — effects, region, the trial — are untouched.)

## E-B — §1.3, the named rows

**old**

- **consistency** = $A \setminus B$;
- **non-consistency (the trigger)** = $A \cap B$;
- **oracle consistency (the oracle joint)** = $A \cap B \cap C^*$;
- **row (7)** adds the two routes to consistency.

**new**

- **consistency** = $A \setminus B$ — the standard pathway: the trial is positive and the region is not above 1.0;
- **non-consistency (the trigger)** = $A \cap B$ — the standard pathway fails, and the investigation begins;
- **oracle consistency (the oracle joint)** = $A \cap B \cap C^*$ — the subgroup pathway, in its oracle form: the trigger occurs and removing $H^*$ exactly leaves the region below 1.0;
- **row (7)** — the two routes together: consistency and oracle consistency are mutually exclusive (the region is above 1.0 in one and not the other), so their probabilities add.

(The event bullets above these — $A$, $B$, $C^*$ — and the preceding sentence are untouched.)
