# What the sign of `effect.threshold` does on the MD path

Prompted by an observation with no obvious explanation: on one draw from a
calibrated harm DGM, `effect.threshold = -30` and `effect.threshold = +30`
returned the **identical** subgroup under `adverse_outcome = FALSE`.

Reading only. Nothing is changed, and nothing here proposes a change.

## Answer

**The sign of `effect.threshold` is not a direction switch.** It is a one-sided
floor on an *orientation-corrected* scale where positive always means harm.
Orientation is set by `adverse_outcome`, never by the threshold.

So `-30` is not "the harm-direction version of `+30`". It is a strictly
**looser** floor than `+30`: it admits everything `+30` admits, plus every
candidate whose oriented effect falls in $[-30, +30)$ — which includes
candidates where treatment *helps*. The two runs agreed because the winner's
oriented effect was far above both floors.

## The executing path

`effect.threshold` is an alias, resolved immediately:

```
R/forestsearch_main.R:1228   user_set_threshold <- !is.null(effect.threshold) || !missing(hr.threshold)
R/forestsearch_main.R:1230   if (!is.null(effect.threshold)) hr.threshold <- effect.threshold
```

For `effect_measure = "MD"` it is then carried through **untransformed**:

```
R/forestsearch_main.R:1642   effect_threshold <- hr.threshold
R/forestsearch_main.R:1734   effect_threshold <- log(effect_threshold)   # ratio measures ONLY
```

`:1734` sits inside the branch for `OR`/`RR`/`IRR`. `MD` is an identity-scale
measure (`is_identity <- effect_measure %in% c("RD","IRD","MD")`, `:1739`), so
it takes no log and no remap. The value reaching the engines is the number the
caller typed, sign included.

It is then handed to two places:

```
R/forestsearch_main.R:2850   search_overrides$effect_threshold <- effect_threshold
R/forestsearch_main.R:2978   consistency_overrides$hr.threshold <- effect_threshold
```

### Every site that compares it

Three, all one-sided, all on a signed quantity:

| site | comparison |
|---|---|
| `R/subgroup_search.R:626` | `if (!disable_effect_floor && glm_result$hr <= hr.threshold) return(status 6)` |
| `R/subgroup_search.R:667` | `if (!disable_effect_floor && cox_result$hr <= hr.threshold) return(status 6)` |
| `R/subgroup_consistency_main.R:528, 531` | `found.hrs <- hr.subgroups[hr.subgroups$HR >= hr.threshold, ]` |

`:626` is the GLM screening reject inside `search_combinations_parallel()`;
`:531` is the re-filter inside `subgroup.consistency()`. Both are "keep
candidates whose effect is **at or above** the floor".

**The threshold is never `abs()`'d.** `grep` for `abs(hr.threshold)`,
`abs(effect_threshold)`, `abs(effect.threshold)`, `abs(screen_threshold)` over
all of `R/` returns nothing. It is never negated, and its sign is never
inspected to choose a comparison direction.

The only site that reads the sign at all is cosmetic:

```
R/subgroup_search.R:194   if (is_ratio && hr.threshold > 0) { ...exp() for display... }
```

That chooses whether to exponentiate for a printed label. `is_ratio` is FALSE
for MD, so on this path it is not even reached.

## Where the direction actually comes from

`R/glm_effect_estimators.R:814-821`, inside the continuous estimator closure:

```r
function(data_slice) {
  # When higher Y = better (adverse_outcome = FALSE), negate Y
  # so that MD > 0 consistently indicates treatment HARM,
  # aligning with the survival convention where HR > 1 = harm
  # and the binary convention where OR > 1 on adverse Y = harm.
  if (!adverse_outcome) {
    data_slice[[outcome.name]] <- -data_slice[[outcome.name]]
  }
```

Every candidate effect entering the three comparisons above has already passed
through this. So the `HR` column is not the mean difference on the analyst's
scale — it is the mean difference on the scale where **positive means harm**.

`adverse_outcome` for continuous outcomes defaults to `FALSE`, resolved at
`R/forestsearch_main.R:1220-1221` from `NULL` via
`outcome_type %in% c("binary","count")`.

## The measurement that closes it

Harm DGM calibrated to MD$(Q) = -40$ on the raw `cd4_change` scale, one draw at
$n = 4000$, `adverse_outcome = FALSE`:

| `effect.threshold` | realized `sg.harm` | in-trial MD, raw scale |
|---|---|---|
| $-30$ | `!{wtkg <= 83} & {symptom}` | $-81.8697$ |
| $+30$ | `!{wtkg <= 83} & {symptom}` | $-81.8697$ |

The selected region's raw MD is $-81.87$ — treatment reduces CD4 change by 82
units, i.e. strong harm. `adverse_outcome = FALSE` negates $Y$, so the value
actually compared against the floor is $+81.87$:

$$+81.87 \ge -30 \quad\text{TRUE}, \qquad +81.87 \ge +30 \quad\text{TRUE}.$$

Both floors admit it, the same winner is selected from the admitted set, and
the runs are bit-identical. Nothing about the sign of the threshold entered the
outcome.

## Consequence worth knowing

The two runs agreeing is not evidence that the sign is ignored in general —
it is evidence that *this* winner cleared both floors. The floors are genuinely
different sets:

- `effect.threshold = +30` keeps candidates with oriented effect $\ge +30$:
  subgroups where treatment is harmful by at least 30 units.
- `effect.threshold = -30` keeps candidates with oriented effect $\ge -30$:
  the above, **plus** every subgroup where treatment is harmful by less than 30,
  neutral, or beneficial by up to 30.

So a negative `effect.threshold` on the MD path is not a harm-direction
specifier. It is a near-vacuous floor. An analyst writing `effect.threshold =
-30` to mean "find harm of at least 30" gets close to the opposite: a filter
that admits almost the whole family, leaving the selection rule to pick
unconstrained.

The correct way to search for harm on a beneficial continuous outcome is
`adverse_outcome = FALSE` (the default) with a **positive** threshold on the
oriented scale. That is what the ACTG175 continuous vignette's
`fs_md_threshold = 30` already is, and it is why the benefit search there needs
`adverse_outcome = TRUE` (vignette lines 138-162) rather than a sign change on
the threshold.

## Not established here

Whether any committed sweep passed a negative `effect.threshold` on an
identity-scale measure. This note traces the mechanism only; it does not audit
call sites. If such a run exists, its candidate family was far larger than the
threshold suggests, and the effect on the selection rate would be substantial.
