# REPORT — reconciling the FS one-sided product ranges

- **Date:** 2026-09-12
- **Scope:** reading committed bundles only. No compute, no re-run.
- **Outcome:** the certification records are **correct as written** and were **not edited**. Two of the
  six "corrections" made in the previous pass of `current_status.md` were **wrong and are withdrawn**.

## The three disagreements, and what each actually was

| product | certification record | recomputed here | what differs |
|---|---|---|---|
| field lower on β(Ĥ) | 0.944–0.974 | 0.941–0.975 (12 cells) | **cell set** — the record predates `p12ext` |
| field-s upper on β(Ĥᶜ) | 0.912–0.960 | 0.9125–0.9605 | **nothing** — rounding convention |
| Bonferroni joint | 0.939–0.963 | 0.932–0.964 | **construction** — record is `joint_s`, mine was `joint` |

## The per-cell values, so the comparison is checkable

Detected replicates with all four field quantities finite. Designated comparator per cell.

| prevalence | HR | n | campaign | n_eval | field lower β(Ĥ) | field-s upper β(Ĥᶜ) | joint | **joint_s** |
|---|---|---|---|---|---|---|---|---|
| 12.4% | 1.50 | 500 | p12ext | 1822 | 0.9654 | 0.9374 | 0.9528 | 0.9555 |
| 12.4% | 1.50 | 1000 | p12ext | 1948 | **0.9410** | 0.9533 | 0.9425 | 0.9435 |
| 12.4% | 1.50 | 1500 | p12ext | 1976 | 0.9565 | 0.9585 | 0.9570 | 0.9580 |
| 12.4% | 1.75 | 500 | tier2 | 1900 | 0.9663 | 0.9405 | 0.9558 | 0.9568 |
| 12.4% | 1.75 | 1000 | tier2 | 1990 | 0.9437 | 0.9563 | 0.9412 | **0.9412** |
| 12.4% | 1.75 | 1500 | tier2 | 1998 | 0.9600 | 0.9605 | 0.9635 | **0.9640** |
| 31% | 1.50 | 500 | e1stud | 1999 | **0.9745** | **0.9125** | **0.9320** | 0.9420 |
| 31% | 1.50 | 1000 | cert20 | 2000 | 0.9585 | 0.9420 | 0.9460 | 0.9500 |
| 31% | 1.50 | 1500 | cert20 | 2000 | 0.9615 | 0.9465 | 0.9440 | 0.9505 |
| 31% | 1.75 | 500 | e1stud | 1999 | 0.9700 | 0.9195 | 0.9335 | **0.9395** |
| 31% | 1.75 | 1000 | cert20 | 2000 | 0.9525 | 0.9420 | 0.9445 | 0.9485 |
| 31% | 1.75 | 1500 | cert20 | 2000 | 0.9610 | 0.9460 | 0.9425 | 0.9490 |

Range endpoints in bold.

| set | field lower | field-s upper | joint | joint_s |
|---|---|---|---|---|
| all 12 harm cells | 0.9410–0.9745 | 0.9125–0.9605 | 0.9320–0.9635 | 0.9395–0.9640 |
| 9 cells, excluding 12.4% HR 1.50 | **0.9437–0.9745** | 0.9125–0.9605 | 0.9320–0.9635 | **0.9395–0.9640** |

## Where each older figure came from

### 1. field lower — a cell-set difference, both correct

`NOTE_survival_products_2026-09-09.md` reports **0.944–0.974**. Its evidence list is `REPORT_cert20`,
`REPORT_tier2`, `REPORT_fixedphat_ij2s`, `REPORT_field_studentize_e1`, `REPORT_cimethod_flip` — it does
**not** include `REPORT_p12ext_2026-09-09`, and `p12ext` is exactly the campaign supplying the three
12.4% HR 1.50 cells. Drop those three and the committed bundles give **0.9437–0.9745 → 0.944–0.974**,
reproducing the record exactly. Keep them and the range is 0.941–0.975, the low end being 12.4% HR 1.50
n 1000 at 0.9410 — the only harm cell below 0.944.

**Unreconciled detail, flagged not resolved:** the note says "ten harm cells"; **nine** is what reproduces
the range, and no tenth bundle on disk carries the field columns (the `map1` / `map1c` / `map1w` and
`p30` bundles at that coordinate return `n_eval` 0 with `effect_neighborhood` NULL — they predate the
field construction). This does not move the range either way.

### 2. field-s upper — no difference at all

The value is **0.9125–0.9605**. "0.912–0.960" and "0.913–0.961" are the same two endpoints under
different rounding. The record's per-cell quotes reproduce **exactly**:

| record says | recomputed |
|---|---|
| 31% HR 1.50: 0.912 → 0.942 → 0.947 | 0.9125, 0.9420, 0.9465 |
| 31% HR 1.75: 0.919 → 0.942 → 0.946 | 0.9195, 0.9420, 0.9460 |
| 12.4% dominated regime: 0.941 / 0.956 / 0.961 | 0.9405, 0.9563, 0.9605 |

The previous pass recorded a "correction" here. **There was nothing to correct; it is withdrawn.**

### 3. Bonferroni joint — a different construction

The record's **0.939–0.963** is the **studentized** pair `fld_joint_s_bonf_loH` / `fld_joint_s_bonf_upHc`,
whose 9-cell range is **0.9395–0.9640**. The 0.932–0.964 quoted against it was the **unscaled** pair
`fld_joint_bonf_loH` / `fld_joint_bonf_upHc`. These are different quantities — the studentized pair runs
0.4 to 1.0 points higher at every 31% cell, and the gap is widest exactly where the unscaled minimum sits
(31% HR 1.50 n 500: `joint` 0.9320 against `joint_s` 0.9420). `REVIEW_E1_fields_2026-09-08.md` shows the
same split independently (`joint_s` 0.942 / 0.939 / 0.949 / 0.951–0.952 against `joint` 0.932 / 0.933 /
0.935 / 0.935).

The previous pass recorded a "correction" here too. **It compared two different constructions; it is
withdrawn.** Name the construction whenever this number is quoted.

## What could not be checked

`SUMMARY_survival_properties_2026-09-10.md` and `REVIEW_certification_2026-09-09.md` are **not in this
repository**; only `dev/notes/NOTE_survival_products_2026-09-09.md` is. The figure **0.941–0.980**
attributed to them therefore could not be verified against its source. **0.980 does not reproduce from any
field-lower computation on the committed bundles** — the maximum over all twelve harm cells is 0.9745. The
nearest 0.98 anywhere in the certification record is the IJ two-sided interval at 31% prevalence,
0.971–0.981, a different product.

## Standing lesson

Quote a coverage range with **the cell set** and **the construction** attached. Two of the three
disagreements above were not disagreements at all, and the third was a cell-set difference in which both
figures are right.
