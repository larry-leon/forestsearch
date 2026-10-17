# Map `sg_focus` to MR's re-selection rule, faithfully per engine.

MR must re-select under the same rule the search used. The only rule
that differs by engine is the `"hr"`/`"eff"` focus: the consistency
search ranks by consistency rate, whose MR analog is `"maxcons"`;
DINA/GRF rank by effect, whose analog is `"maxeff"`. Every other focus
maps identically (`maxSG`/`minSG` pass through; `hrMaxSG`/`hrMinSG` -\>
the size-within-effect- neighborhood rules `effMaxSG`/`effMinSG`). This
makes `sg_focus` the single source of truth, so callers never need to
specify `reselection` by hand.

## Usage

``` r
.fs_mr_reselection_from_focus(sg_focus, engine = c("effect", "consistency"))
```

## Arguments

- sg_focus:

  The search focus (`"eff"`, `"maxSG"`, `"effMaxSG"`, ...).

- engine:

  `"consistency"` (so `"hr"` -\> `"maxcons"`) or `"effect"` (DINA/GRF;
  `"hr"` -\> `"maxeff"`).
