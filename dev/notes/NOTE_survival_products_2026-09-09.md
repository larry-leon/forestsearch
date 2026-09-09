# NOTE — Survival post-selection products (certified 2026-09-09)

Date: 2026-09-09. Source: `dev/tasks/TASK_cimethod_note_2026-09-09.md`, Part N2 (text as specified there).
Supersedes `dev/notes/NOTE_complement_product_2026-09-08.md`. Evidence: `REPORT_cert20_2026-09-08.md`,
`REPORT_tier2_2026-09-08.md`, `REPORT_fixedphat_ij2s_2026-09-09.md` (C-4),
`REPORT_field_studentize_e1_2026-09-08.md`, `REPORT_cimethod_flip_2026-09-09.md` (the defaults).

**Survival post-selection products (certified 2026-09-09).** Defaults are the recommendations: `ci_method = "field"`, `field_complement = TRUE`, `field_scale_complement = "selected"`, `return_reselection = TRUE`.

**Harm subgroup Ĥ — one-sided lower bound on β(Ĥ), the field.** Coverage 0.944–0.974 across ten harm cells at 12.4% and 31% prevalence, n = 500 / 1000 / 1500; flat in n.

**Complement Ĥᶜ — one-sided upper bound on β(Ĥᶜ), field-s** (the studentized complement field, R1). **The small-sample shortfall closes with n:** 0.912 → 0.942 → 0.947 (31%, HR 1.50) and 0.919 → 0.942 → 0.946 (HR 1.75); ≥ 0.94 at every cell of the 12.4% dominated regime (0.941 / 0.956 / 0.961). field-s is at or above the unscaled field in all thirteen cells.

**Joint two-subgroup claim:** Bonferroni, γ = 0.025 each side; 0.939–0.963 at every cell.

**Two-sided intervals are not certified.** The IJ two-term two-sided interval is retained and reported as the secondary, conservative option and is the only two-sided construction offered, but its harm-block coverage falls to 0.913–0.917 at 12.4% prevalence with n ≥ 1000 (0.971–0.981 at 31%). The field's own two-sided is lower still (0.878–0.930). Read two-sided statements at low prevalence and large n with that caveat.

**Analysis-time diagnostic.** p̂(Ĥ), the field's re-selection frequency, is recorded and reported; **no construction reads it**. Harm-block bias is a monotone increasing function of p̂, crossing zero near p̂ ≈ 0.5: **over-correction at low p̂** (bias −0.11 to −0.28 log units at 12.4%, n = 1500) and **under-correction at high p̂** (the stable-pick regime, +0.02). The one-sided products are each exposed to one pole only, and in the conservative direction, which is why they certify while the two-sided interval — exposed to both — does not. Complement coverage is a stable function of p̂ (flat-to-improving at fixed band in all 18 usable band × series combinations); the harm block's is not common across prevalence, so its p̂ flag is directional, not calibrated.

**Caveats on record.** ε > 0.25 is not adoptable (Larry, 2026-09-08); under a pure-size pick (`maxSG`) the naive SE mis-calibrates (0.876 of the error SD) and field-s inherits it; fixed p̂ bands are fixed in value but not in meaning, since the p̂ distribution shifts with n (replicates below p̂ = 0.20 are 52% of the 12.4% HR 1.75 cell at n = 500 and 8% at n = 1500).
