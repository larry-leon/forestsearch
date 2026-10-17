# Extract DINA per-covariate Pareto frontiers as forestsearch cuts

Runs the same univariate (covariate, direction, threshold) search as
[`dina_subgroup()`](https://larry-leon.github.io/forestsearch/reference/dina_subgroup.md),
but instead of selecting one subgroup it returns a pool of
**per-covariate Pareto-non-dominated** cuts. Each cut is rendered as a
forestsearch cut expression (`"x1 <= 0.5"`), so the pool can be fed to
[`forestsearch()`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
as screening-stage candidates in the same way GRF tree cuts are. This is
the DINA analog of GRF candidate discovery.

## Usage

``` r
dina_frontier(
  fit,
  df,
  covariates,
  scope = c("wide", "harm"),
  m_diff = NULL,
  n_min = 60L,
  direction = c("both", "left", "right"),
  max_per_covariate = Inf,
  max_subgroups = Inf,
  digits = 3L
)
```

## Arguments

- fit:

  a fitted DINA object (class `"dina"` or `"dina_bagged"`).

- df:

  data frame containing the covariate columns.

- covariates:

  character vector of numeric covariate columns to search over.

- scope:

  one of `"wide"` (default) or `"harm"`. `"wide"` casts a broad net:
  every candidate with size `>= n_min` is eligible, and the
  per-covariate frontiers span the whole effect-size boundary (mirroring
  GRF's wide-net role, leaving forestsearch's consistency stage as the
  gatekeeper). `"harm"` restricts the eligible set to candidates with
  mean tau-hat `>= m_diff` before taking the frontiers.

- m_diff:

  scalar harm floor on the link (natural-parameter) scale. Required when
  `scope = "harm"`; ignored when `scope = "wide"`. Default `NULL`.

- n_min:

  positive integer minimum subgroup size. Default `60L`.

- direction:

  one of `"both"` (default), `"left"`, `"right"`.

- max_per_covariate:

  positive integer, or `Inf` for no limit: the most cuts any single
  covariate may contribute to the **returned table** (its top cuts by
  effect). A report-trimming cap, not a search control. Default `Inf` –
  show everything. A finite value that actually trims warns, naming the
  kept and available counts.

- max_subgroups:

  positive integer, or `Inf` for no limit: the cap on the total pooled
  cuts in the **returned table**. When the pool is larger it is trimmed
  round-robin by within-covariate effect rank. Default `Inf`; a finite
  value that actually trims warns, naming the kept and available counts.

  Both caps default to `Inf` for the same reason
  [`forestsearch()`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)'s
  `max_subgroups_search` does: a finite pooled cap drops covariates
  whole and truncates the rank-2 pass, so the table stops being a
  faithful picture of the frontier. The two are still not the same knob
  – `max_subgroups_search` truncates the pool of candidate subgroups
  that forestsearch actually **evaluates** at the consistency stage,
  while `max_subgroups` trims the rows of the **table this function
  returns**. One changes what is computed; the other changes what is
  reported.

- digits:

  significant figures for rounding the emitted threshold in `cut_expr`.
  Default `3L`.

## Value

A data frame (one row per cut, after the caps and threshold-rounding
deduplication) with columns:

- covariate, direction, threshold:

  the subgroup definition (raw threshold).

- n_subgroup:

  subgroup size.

- mean_tau_hat:

  subgroup-mean tau-hat on the link scale.

- effect:

  the same on the natural scale (exponentiated for ratio families
  binomial/poisson/cox), used for the frontiers and the effect ranking.

- cut_expr:

  the canonical forestsearch cut expression.

The result carries attributes `n_frontier` (total per-covariate frontier
size, deduped, before the caps), `n_qualifying`, and `scope`. Returns a
0-row data frame (same columns) when no candidate qualifies.

## Details

**Why per-covariate and not a single global frontier.** The global
Pareto frontier in (effect, N) space collapses onto whichever single
covariate dominates the trade-off boundary, returning many micro-stepped
cuts on that one covariate. That is faithful to the definition but the
opposite of what forestsearch needs: its value is in *composing* cuts
drawn from several covariates. So this function computes a separate
frontier for each covariate, then pools them, mirroring how a GRF tree
splits on multiple variables.

**Two caps, each optional.** Both caps act on **this function's return
value** and nothing else: they trim the frontier table of **single
cuts** that is returned here. They are not a cap on a selection family –
the rows they trim are univariate cuts, never conjunctions – and
[`dina_subgroup()`](https://larry-leon.github.io/forestsearch/reference/dina_subgroup.md)
does not take them, so nothing they do reaches the subgroup selector.
Either may be `Inf` to disable it:

- `max_per_covariate` bounds how many cuts any single covariate
  contributes (its top cuts by effect), limiting within-covariate
  redundancy.

- `max_subgroups` bounds the total pool size across covariates.

Under `forestsearch(subgroup_method = "dina")` neither cap has any
effect: that path selects with
[`dina_subgroup()`](https://larry-leon.github.io/forestsearch/reference/dina_subgroup.md)
and never consults a frontier table. They bite only where this
function's output is used directly – including
`forestsearch(use_dina = TRUE, dina_args = list(selected_only = FALSE))`,
where `cut_expr` becomes the screening-stage candidate pool. When the
pool exceeds `max_subgroups`, the global trim is round-robin on
within-covariate effect rank – every covariate's best cut first (ordered
by effect), then every covariate's second-best, and so on – so a finite
budget preserves cross-covariate spread instead of refilling with one
covariate's micro-steps. Setting `max_per_covariate = Inf` recovers a
pure global budget; setting `max_subgroups = Inf` recovers
per-covariate-only limits. Both default to `Inf`, and a finite cap that
actually removes rows warns rather than trimming the display silently.

**Cut form.** The expression is always the canonical
`"<covariate> <= <threshold>"`. forestsearch's factor machinery exposes
both the cut and its complement, so one expression covers both of DINA's
left/right directions – the `direction` column is retained for reference
only. Thresholds are rounded to `digits` significant figures so they
deduplicate cleanly against existing GRF / median cuts and read well;
raw values remain in `threshold`.

These cuts are **candidates**, not forced selections: forestsearch
composes (up to `maxk`) and consistency-gatekeeps them, exactly as it
does GRF cuts. To use them without disturbing any user-supplied forced
cuts, append rather than overwrite, e.g.
`forestsearch(..., conf_force = c(my_forced_cuts, fr$cut_expr))`.

## See also

[`dina_subgroup()`](https://larry-leon.github.io/forestsearch/reference/dina_subgroup.md)
for the single-subgroup selector;
[`forestsearch()`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
for the screening stage that consumes the cuts;
[`compute_pareto_frontier()`](https://larry-leon.github.io/forestsearch/reference/compute_pareto_frontier.md)
for the post-hoc frontier reporting on a forestsearch result.

## Examples

``` r
set.seed(1)
n <- 400
df_demo <- data.frame(
  w  = stats::rbinom(n, 1, 0.5),
  x1 = stats::runif(n, -1, 1),
  x2 = stats::runif(n, -1, 1),
  x3 = stats::runif(n, -1, 1)
)
tau_x      <- 0.3 + 1.2 * df_demo$x1 - 0.4 * df_demo$x2
df_demo$y  <- 0.5 * df_demo$x1 + df_demo$w * tau_x + stats::rnorm(n)

fit <- dina(df_demo, outcome = "y", treatment = "w",
            covariates = c("x1", "x2", "x3"),
            family = "gaussian", seed = 1L)

fr <- dina_frontier(fit, df_demo, covariates = c("x1", "x2", "x3"),
                    n_min = 60L)
fr
#>     covariate direction    threshold n_subgroup mean_tau_hat    effect
#> 1          x1     right  0.703335280         60    1.2688848 1.2688848
#> 2          x2      left -0.680430355         71    0.6245141 0.6245141
#> 3          x3      left -0.721954309         70    0.5733911 0.5733911
#> 4          x1     right  0.702152439         61    1.2598885 1.2598885
#> 5          x2      left -0.678410738         73    0.6035410 0.6035410
#> 6          x3      left -0.718981634         71    0.5660133 0.5660133
#> 7          x1     right  0.699104208         62    1.2598806 1.2598806
#> 8          x2      left -0.676807748         74    0.5910245 0.5910245
#> 9          x3      left -0.713463377         72    0.5544840 0.5544840
#> 10         x1     right  0.692491505         64    1.2570180 1.2570180
#> 11         x2      left -0.673584746         76    0.5865057 0.5865057
#> 12         x3      left -0.602254387         93    0.5517411 0.5517411
#> 13         x1     right  0.686240297         65    1.2564136 1.2564136
#> 14         x2      left -0.671508636         78    0.5638887 0.5638887
#> 15         x3      left -0.600961623         94    0.5458879 0.5458879
#> 16         x1     right  0.685179807         66    1.2553659 1.2553659
#> 17         x2      left -0.666575667         79    0.5578463 0.5578463
#> 18         x3      left -0.552653160        100    0.5448111 0.5448111
#> 19         x1     right  0.676768133         67    1.2498341 1.2498341
#> 20         x2      left -0.651964913         80    0.5457674 0.5457674
#> 21         x3      left -0.530570640        104    0.5407226 0.5407226
#> 22         x1     right  0.671035900         68    1.2489137 1.2489137
#> 23         x2      left -0.650676642         81    0.5427627 0.5427627
#> 24         x3      left -0.520786579        107    0.5398300 0.5398300
#> 25         x1     right  0.660491408         69    1.2468367 1.2468367
#> 26         x2      left -0.644187837         83    0.5413039 0.5413039
#> 27         x3      left -0.520241661        108    0.5291471 0.5291471
#> 28         x1     right  0.658648609         71    1.2394233 1.2394233
#> 29         x2      left -0.585591721         91    0.5393589 0.5393589
#> 30         x3      left -0.518151912        109    0.5277379 0.5277379
#> 31         x1     right  0.658402128         72    1.2369583 1.2369583
#> 32         x2      left -0.583921973         92    0.5378181 0.5378181
#> 33         x3      left -0.505056653        110    0.5216154 0.5216154
#> 34         x1     right  0.657068720         73    1.2322015 1.2322015
#> 35         x2      left -0.554588046         96    0.5276573 0.5276573
#> 36         x3      left -0.493302190        111    0.5199359 0.5199359
#> 37         x1     right  0.656150036         74    1.2264731 1.2264731
#> 38         x2      left -0.478592269        113    0.5219530 0.5219530
#> 39         x3      left -0.492405869        112    0.5099996 0.5099996
#> 40         x1     right  0.655435539         75    1.2238430 1.2238430
#> 41         x2      left -0.475135616        114    0.5181378 0.5181378
#> 42         x3      left -0.490128237        113    0.5015058 0.5015058
#> 43         x1     right  0.642324636         77    1.2236166 1.2236166
#> 44         x2      left -0.470758745        115    0.5153154 0.5153154
#> 45         x3      left -0.488291352        114    0.4959922 0.4959922
#> 46         x1     right  0.640559987         78    1.2189771 1.2189771
#> 47         x2      left -0.468926560        116    0.5122064 0.5122064
#> 48         x3      left -0.404048530        123    0.4930466 0.4930466
#> 49         x1     right  0.640343234         79    1.2156520 1.2156520
#> 50         x2      left -0.348568012        141    0.5120922 0.5120922
#> 51         x3      left -0.402062191        124    0.4867269 0.4867269
#> 52         x1     right  0.628503501         80    1.2148243 1.2148243
#> 53         x2      left -0.343684039        142    0.5100082 0.5100082
#> 54         x3      left -0.398158496        125    0.4859818 0.4859818
#> 55         x1     right  0.613083655         81    1.2111468 1.2111468
#> 56         x2      left -0.342915395        143    0.5059999 0.5059999
#> 57         x3      left -0.336409070        134    0.4839188 0.4839188
#> 58         x1     right  0.608781543         82    1.2094876 1.2094876
#> 59         x2      left -0.324358549        150    0.4968166 0.4968166
#> 60         x3      left -0.266508638        147    0.4830362 0.4830362
#> 61         x1     right  0.603375411         83    1.2037742 1.2037742
#> 62         x2      left -0.311815692        153    0.4887928 0.4887928
#> 63         x3      left -0.265900185        148    0.4809056 0.4809056
#> 64         x1     right  0.597629099         84    1.1990347 1.1990347
#> 65         x2      left -0.286105680        159    0.4868585 0.4868585
#> 66         x3      left -0.256946704        149    0.4803424 0.4803424
#> 67         x1     right  0.585818048         85    1.1951582 1.1951582
#> 68         x2      left -0.266053581        164    0.4860626 0.4860626
#> 69         x3      left -0.250869562        152    0.4797803 0.4797803
#> 70         x1     right  0.584757512         86    1.1912002 1.1912002
#> 71         x2      left -0.245392972        165    0.4792209 0.4792209
#> 72         x3      left -0.244002384        154    0.4781767 0.4781767
#> 73         x1     right  0.578188944         87    1.1881781 1.1881781
#> 74         x2      left -0.230703273        170    0.4766202 0.4766202
#> 75         x3      left -0.242830940        155    0.4723653 0.4723653
#> 76         x1     right  0.566195342         88    1.1821718 1.1821718
#> 77         x2      left -0.229521150        171    0.4742572 0.4742572
#> 78         x3      left -0.239259750        156    0.4707339 0.4707339
#> 79         x1     right  0.554289971         92    1.1766400 1.1766400
#> 80         x2      left -0.226510676        172    0.4697066 0.4697066
#> 81         x3      left -0.195871952        169    0.4677359 0.4677359
#> 82         x1     right  0.540450711         94    1.1759758 1.1759758
#> 83         x2      left -0.204587816        175    0.4695183 0.4695183
#> 84         x3      left -0.195326980        170    0.4654322 0.4654322
#> 85         x1     right  0.536289101         95    1.1756181 1.1756181
#> 86         x2      left -0.194025596        181    0.4658506 0.4658506
#> 87         x3      left -0.192771824        172    0.4565489 0.4565489
#> 88         x1     right  0.517393614         98    1.1733108 1.1733108
#> 89         x2      left -0.173322094        186    0.4622235 0.4622235
#> 90         x3      left -0.187805741        173    0.4544102 0.4544102
#> 91         x1     right  0.513289815         99    1.1718685 1.1718685
#> 92         x2      left -0.168279505        187    0.4576871 0.4576871
#> 93         x3      left -0.174825494        176    0.4507896 0.4507896
#> 94         x1     right  0.506665401        100    1.1695748 1.1695748
#> 95         x2      left -0.143305997        191    0.4542592 0.4542592
#> 96         x3      left -0.173938198        177    0.4484864 0.4484864
#> 97         x1     right  0.503008336        101    1.1671145 1.1671145
#> 98         x2      left -0.132396476        192    0.4512316 0.4512316
#> 99         x3      left -0.161392022        179    0.4462556 0.4462556
#> 100        x1     right  0.502443971        102    1.1633215 1.1633215
#> 101        x2      left -0.127041490        194    0.4490260 0.4490260
#> 102        x3      left -0.157440621        180    0.4445442 0.4445442
#> 103        x1     right  0.495270127        105    1.1584027 1.1584027
#> 104        x2      left -0.122919864        195    0.4459759 0.4459759
#> 105        x3      left -0.125951818        182    0.4435965 0.4435965
#> 106        x1     right  0.479283492        106    1.1578676 1.1578676
#> 107        x2      left -0.115683153        196    0.4407365 0.4407365
#> 108        x3      left -0.112592765        184    0.4400100 0.4400100
#> 109        x1     right  0.474610936        107    1.1530293 1.1530293
#> 110        x3      left -0.110194743        185    0.4385777 0.4385777
#> 111        x2      left -0.075285903        205    0.4278624 0.4278624
#> 112        x1     right  0.465436035        108    1.1478119 1.1478119
#> 113        x3      left -0.069317629        194    0.4372533 0.4372533
#> 114        x2      left -0.074539450        206    0.4242609 0.4242609
#> 115        x1     right  0.462743158        109    1.1466991 1.1466991
#> 116        x3      left -0.062480016        195    0.4307211 0.4307211
#> 117        x2      left -0.061241925        207    0.4212646 0.4212646
#> 118        x1     right  0.461865768        110    1.1418303 1.1418303
#> 119        x3      left -0.061237852        196    0.4237950 0.4237950
#> 120        x2      left  0.011062160        221    0.4204542 0.4204542
#> 121        x1     right  0.459606576        111    1.1350202 1.1350202
#> 122        x3      left -0.060140077        197    0.4237074 0.4237074
#> 123        x2      left  0.011763969        222    0.4191157 0.4191157
#> 124        x1     right  0.455534009        112    1.1331609 1.1331609
#> 125        x3      left -0.029038675        202    0.4198615 0.4198615
#> 126        x2      left  0.248729095        254    0.4186161 0.4186161
#> 127        x1     right  0.448981151        114    1.1301391 1.1301391
#> 128        x3      left -0.011183774        203    0.4174679 0.4174679
#> 129        x2      left  0.262097495        255    0.4137661 0.4137661
#> 130        x1     right  0.447381502        115    1.1248530 1.1248530
#> 131        x3      left -0.008717267        205    0.4171907 0.4171907
#> 132        x2      left  0.271899697        256    0.4122713 0.4122713
#> 133        x1     right  0.433178958        116    1.1229985 1.1229985
#> 134        x3      left  0.069179517        221    0.4169936 0.4169936
#> 135        x2      left  0.277511939        257    0.4091690 0.4091690
#> 136        x1     right  0.431666567        117    1.1201249 1.1201249
#> 137        x3      left  0.081284562        222    0.4144320 0.4144320
#> 138        x2      left  0.280236016        258    0.4048044 0.4048044
#> 139        x1     right  0.429932425        118    1.1156146 1.1156146
#> 140        x3      left  0.104187589        223    0.4097992 0.4097992
#> 141        x2      left  0.380783410        275    0.4029959 0.4029959
#> 142        x1     right  0.429435899        119    1.1140052 1.1140052
#> 143        x3      left  0.110808264        224    0.4074548 0.4074548
#> 144        x2      left  0.464842395        294    0.4029163 0.4029163
#> 145        x1     right  0.425080245        120    1.1087085 1.1087085
#> 146        x3      left  0.295529748        251    0.4058188 0.4058188
#> 147        x2      left  0.470007293        295    0.3986011 0.3986011
#> 148        x1     right  0.424449119        122    1.1019124 1.1019124
#> 149        x3      left  0.299424498        252    0.4043498 0.4043498
#> 150        x2      left  0.490220272        301    0.3956395 0.3956395
#> 151        x1     right  0.417950397        123    1.0991991 1.0991991
#> 152        x3      left  0.304508257        253    0.4032650 0.4032650
#> 153        x2      left  0.503091746        302    0.3911285 0.3911285
#> 154        x1     right  0.417480102        124    1.0963168 1.0963168
#> 155        x3      left  0.355308401        263    0.4018807 0.4018807
#> 156        x2      left  0.557089280        312    0.3904338 0.3904338
#> 157        x1     right  0.413323064        126    1.0896928 1.0896928
#> 158        x3      left  0.357681620        264    0.3963081 0.3963081
#> 159        x2      left  0.557678083        313    0.3878836 0.3878836
#> 160        x1     right  0.411718013        127    1.0879574 1.0879574
#> 161        x3      left  0.464560100        286    0.3934001 0.3934001
#> 162        x2      left  0.569870956        314    0.3875589 0.3875589
#> 163        x1     right  0.392393999        128    1.0870470 1.0870470
#> 164        x3      left  0.479212996        288    0.3929228 0.3929228
#> 165        x2      left  0.578462979        315    0.3866540 0.3866540
#> 166        x1     right  0.390024493        129    1.0849139 1.0849139
#> 167        x3      left  0.483934311        290    0.3923850 0.3923850
#> 168        x2      left  0.579828223        317    0.3848849 0.3848849
#> 169        x1     right  0.383855344        130    1.0784992 1.0784992
#> 170        x3      left  0.487565712        291    0.3911434 0.3911434
#> 171        x2      left  0.602283774        320    0.3832145 0.3832145
#> 172        x1     right  0.380513446        131    1.0782152 1.0782152
#> 173        x3      left  0.495411849        292    0.3899298 0.3899298
#> 174        x2      left  0.604258017        322    0.3805709 0.3805709
#> 175        x1     right  0.378914733        132    1.0725781 1.0725781
#> 176        x3      left  0.500796756        293    0.3850739 0.3850739
#> 177        x2      left  0.627611567        325    0.3781067 0.3781067
#> 178        x1     right  0.376550475        133    1.0704258 1.0704258
#> 179        x3      left  0.508248102        294    0.3840109 0.3840109
#> 180        x2      left  0.631128349        326    0.3774838 0.3774838
#> 181        x1     right  0.370215010        134    1.0692598 1.0692598
#> 182        x3      left  0.508534445        295    0.3833653 0.3833653
#> 183        x2      left  0.632111010        327    0.3744077 0.3744077
#> 184        x1     right  0.353137545        135    1.0676849 1.0676849
#> 185        x3      left  0.509554386        297    0.3750066 0.3750066
#> 186        x2      left  0.634378314        328    0.3726455 0.3726455
#> 187        x1     right  0.351707521        136    1.0664197 1.0664197
#> 188        x3      left  0.520140645        300    0.3744625 0.3744625
#> 189        x2      left  0.678607816        336    0.3707058 0.3707058
#> 190        x1     right  0.342816656        138    1.0632814 1.0632814
#> 191        x3      left  0.524174781        301    0.3743594 0.3743594
#> 192        x2      left  0.686161104        337    0.3678463 0.3678463
#> 193        x1     right  0.336101208        139    1.0626701 1.0626701
#> 194        x3      left  0.525966711        303    0.3736196 0.3736196
#> 195        x2      left  0.717374906        339    0.3670193 0.3670193
#> 196        x1     right  0.335495255        140    1.0597462 1.0597462
#> 197        x3      left  0.539632931        304    0.3705315 0.3705315
#> 198        x2      left  0.727201308        341    0.3632711 0.3632711
#> 199        x1     right  0.326685451        141    1.0560898 1.0560898
#> 200        x3      left  0.585190728        307    0.3698059 0.3698059
#> 201        x2      left  0.729235653        342    0.3604138 0.3604138
#> 202        x1     right  0.326021935        142    1.0519770 1.0519770
#> 203        x3      left  0.589600580        309    0.3663399 0.3663399
#> 204        x2      left  0.739314530        343    0.3588550 0.3588550
#> 205        x1     right  0.317755218        143    1.0485114 1.0485114
#> 206        x3      left  0.611831109        312    0.3655173 0.3655173
#> 207        x2      left  0.748565988        345    0.3573137 0.3573137
#> 208        x1     right  0.315656139        145    1.0477226 1.0477226
#> 209        x3      left  0.627273175        316    0.3649294 0.3649294
#> 210        x2      left  0.753716824        346    0.3535486 0.3535486
#> 211        x1     right  0.309465752        146    1.0440179 1.0440179
#> 212        x3      left  0.645078209        320    0.3644048 0.3644048
#> 213        x2      left  0.757688481        347    0.3493098 0.3493098
#> 214        x1     right  0.306698886        147    1.0412650 1.0412650
#> 215        x3      left  0.757645758        354    0.3618841 0.3618841
#> 216        x2      left  0.762084763        348    0.3485445 0.3485445
#> 217        x1     right  0.303881887        148    1.0377794 1.0377794
#> 218        x3      left  0.762670510        356    0.3581664 0.3581664
#> 219        x2      left  0.765993481        349    0.3485031 0.3485031
#> 220        x1     right  0.297101895        149    1.0349834 1.0349834
#> 221        x3      left  0.800239094        361    0.3545041 0.3545041
#> 222        x2      left  0.785820825        352    0.3457153 0.3457153
#> 223        x1     right  0.296026892        150    1.0349439 1.0349439
#> 224        x3      left  0.805284707        363    0.3539009 0.3539009
#> 225        x2      left  0.804293122        357    0.3456488 0.3456488
#> 226        x1     right  0.293358380        151    1.0307823 1.0307823
#> 227        x3      left  0.806931432        364    0.3503217 0.3503217
#> 228        x2      left  0.812646602        358    0.3422658 0.3422658
#> 229        x1     right  0.292134554        152    1.0305349 1.0305349
#> 230        x3      left  0.814397507        365    0.3472268 0.3472268
#> 231        x2      left  0.841654934        364    0.3397702 0.3397702
#> 232        x1     right  0.288214457        153    1.0271473 1.0271473
#> 233        x3      left  0.867012500        374    0.3458660 0.3458660
#> 234        x2      left  0.850385122        366    0.3375560 0.3375560
#> 235        x1     right  0.275817446        154    1.0234049 1.0234049
#> 236        x3      left  0.869432650        375    0.3423095 0.3423095
#> 237        x2      left  0.850646952        367    0.3360744 0.3360744
#> 238        x1     right  0.269789672        155    1.0202580 1.0202580
#> 239        x3      left  0.870193809        376    0.3410977 0.3410977
#> 240        x2      left  0.947453372        391    0.3356597 0.3356597
#> 241        x1     right  0.253309695        157    1.0149015 1.0149015
#> 242        x3      left  0.920534328        385    0.3408662 0.3408662
#> 243        x2      left  0.957116416        392    0.3345091 0.3345091
#> 244        x1     right  0.250490671        158    1.0121261 1.0121261
#> 245        x3      left  0.922198355        386    0.3398027 0.3398027
#> 246        x2      left  0.960729562        393    0.3312434 0.3312434
#> 247        x1     right  0.246871419        159    1.0109691 1.0109691
#> 248        x3      left  0.923529374        388    0.3383089 0.3383089
#> 249        x2      left  0.977667186        395    0.3293403 0.3293403
#> 250        x1     right  0.243546559        160    1.0095391 1.0095391
#> 251        x3      left  0.928411382        390    0.3376231 0.3376231
#> 252        x2     right -0.995961535        399    0.3280153 0.3280153
#> 253        x1     right  0.242152588        161    1.0077010 1.0077010
#> 254        x3      left  0.950453431        393    0.3374991 0.3374991
#> 255        x2      left  0.999861187        400    0.3263429 0.3263429
#> 256        x1     right  0.240855851        162    1.0049668 1.0049668
#> 257        x3      left  0.951341521        394    0.3347248 0.3347248
#> 258        x2     right -0.997370687        400    0.3263429 0.3263429
#> 259        x1     right  0.219497870        163    1.0028384 1.0028384
#> 260        x3      left  0.954685066        396    0.3312586 0.3312586
#> 261        x1     right  0.212881989        165    0.9979469 0.9979469
#> 262        x3      left  0.956756141        397    0.3288099 0.3288099
#> 263        x1     right  0.207584803        166    0.9944128 0.9944128
#> 264        x3      left  0.968657797        398    0.3280966 0.3280966
#> 265        x1     right  0.207191383        167    0.9895208 0.9895208
#> 266        x3      left  0.992543569        400    0.3263429 0.3263429
#> 267        x1     right  0.196790836        168    0.9868198 0.9868198
#> 268        x3     right -0.996715478        400    0.3263429 0.3263429
#> 269        x1     right  0.194480971        169    0.9820911 0.9820911
#> 270        x1     right  0.192896911        170    0.9817705 0.9817705
#> 271        x1     right  0.176838381        171    0.9780376 0.9780376
#> 272        x1     right  0.176132786        172    0.9770374 0.9770374
#> 273        x1     right  0.156568535        174    0.9697487 0.9697487
#> 274        x1     right  0.156323833        175    0.9679280 0.9679280
#> 275        x1     right  0.144340519        176    0.9669261 0.9669261
#> 276        x1     right  0.135374795        177    0.9619321 0.9619321
#> 277        x1     right  0.130537955        178    0.9572577 0.9572577
#> 278        x1     right  0.123976185        179    0.9525006 0.9525006
#> 279        x1     right  0.118620541        180    0.9482741 0.9482741
#> 280        x1     right  0.108354122        181    0.9443337 0.9443337
#> 281        x1     right  0.106951407        182    0.9421706 0.9421706
#> 282        x1     right  0.102527934        183    0.9375243 0.9375243
#> 283        x1     right  0.096340829        184    0.9362617 0.9362617
#> 284        x1     right  0.084672020        185    0.9353344 0.9353344
#> 285        x1     right  0.083958560        186    0.9328107 0.9328107
#> 286        x1     right  0.078924056        187    0.9294394 0.9294394
#> 287        x1     right  0.077588599        188    0.9250237 0.9250237
#> 288        x1     right  0.075729846        189    0.9228262 0.9228262
#> 289        x1     right  0.066121212        190    0.9222153 0.9222153
#> 290        x1     right  0.065202439        191    0.9206021 0.9206021
#> 291        x1     right  0.061648728        192    0.9172664 0.9172664
#> 292        x1     right  0.056097758        193    0.9136337 0.9136337
#> 293        x1     right  0.054192389        194    0.9119724 0.9119724
#> 294        x1     right  0.050611007        195    0.9092913 0.9092913
#> 295        x1     right  0.044496552        196    0.9066200 0.9066200
#> 296        x1     right  0.035901141        197    0.9051976 0.9051976
#> 297        x1     right  0.034019653        198    0.9033804 0.9033804
#> 298        x1     right  0.032247829        199    0.9006123 0.9006123
#> 299        x1     right  0.031052301        200    0.8991126 0.8991126
#> 300        x1     right  0.029465306        201    0.8972185 0.8972185
#> 301        x1     right  0.014748510        202    0.8946959 0.8946959
#> 302        x1     right  0.011686048        203    0.8911508 0.8911508
#> 303        x1     right  0.004501631        204    0.8881840 0.8881840
#> 304        x1     right  0.003509650        205    0.8845257 0.8845257
#> 305        x1     right  0.001700975        206    0.8834972 0.8834972
#> 306        x1     right -0.004263896        207    0.8794257 0.8794257
#> 307        x1     right -0.006629560        208    0.8766525 0.8766525
#> 308        x1     right -0.009871789        209    0.8732474 0.8732474
#> 309        x1     right -0.010690948        210    0.8681073 0.8681073
#> 310        x1     right -0.028655786        211    0.8658745 0.8658745
#> 311        x1     right -0.035660854        212    0.8621492 0.8621492
#> 312        x1     right -0.040328411        213    0.8582475 0.8582475
#> 313        x1     right -0.045073578        214    0.8556960 0.8556960
#> 314        x1     right -0.047623963        215    0.8534360 0.8534360
#> 315        x1     right -0.048349923        216    0.8492606 0.8492606
#> 316        x1     right -0.049204159        217    0.8462390 0.8462390
#> 317        x1     right -0.052993803        218    0.8425054 0.8425054
#> 318        x1     right -0.054883769        219    0.8389561 0.8389561
#> 319        x1     right -0.055024158        220    0.8355549 0.8355549
#> 320        x1     right -0.055531595        221    0.8321564 0.8321564
#> 321        x1     right -0.059742765        222    0.8285685 0.8285685
#> 322        x1     right -0.065769769        223    0.8253127 0.8253127
#> 323        x1     right -0.078209641        224    0.8233048 0.8233048
#> 324        x1     right -0.080849176        225    0.8212639 0.8212639
#> 325        x1     right -0.087716532        226    0.8196725 0.8196725
#> 326        x1     right -0.096316974        227    0.8153699 0.8153699
#> 327        x1     right -0.099029126        228    0.8124721 0.8124721
#> 328        x1     right -0.103145070        229    0.8092410 0.8092410
#> 329        x1     right -0.114979784        230    0.8067847 0.8067847
#> 330        x1     right -0.115980429        231    0.8027928 0.8027928
#> 331        x1     right -0.117081714        232    0.8003767 0.8003767
#> 332        x1     right -0.123685206        233    0.7960940 0.7960940
#> 333        x1     right -0.128047197        234    0.7926166 0.7926166
#> 334        x1     right -0.130198756        235    0.7911333 0.7911333
#> 335        x1     right -0.130564719        236    0.7884809 0.7884809
#> 336        x1     right -0.131938303        237    0.7858073 0.7858073
#> 337        x1     right -0.134453740        238    0.7833255 0.7833255
#> 338        x1     right -0.137088040        239    0.7808101 0.7808101
#> 339        x1     right -0.143049238        240    0.7777438 0.7777438
#> 340        x1     right -0.147803798        241    0.7752464 0.7752464
#> 341        x1     right -0.151692948        242    0.7727529 0.7727529
#> 342        x1     right -0.152883163        243    0.7697241 0.7697241
#> 343        x1     right -0.153802781        244    0.7661050 0.7661050
#> 344        x1     right -0.154809684        245    0.7643461 0.7643461
#> 345        x1     right -0.160053360        246    0.7617403 0.7617403
#> 346        x1     right -0.163322017        247    0.7598204 0.7598204
#> 347        x1     right -0.169363997        248    0.7575310 0.7575310
#> 348        x1     right -0.183444480        249    0.7537699 0.7537699
#> 349        x1     right -0.188419000        250    0.7501480 0.7501480
#> 350        x1     right -0.195605875        252    0.7481584 0.7481584
#> 351        x1     right -0.201794888        253    0.7467464 0.7467464
#> 352        x1     right -0.204312745        254    0.7453718 0.7453718
#> 353        x1     right -0.205704133        255    0.7423222 0.7423222
#> 354        x1     right -0.212101886        256    0.7404056 0.7404056
#> 355        x1     right -0.215268230        257    0.7364938 0.7364938
#> 356        x1     right -0.223843022        258    0.7336451 0.7336451
#> 357        x1     right -0.227800734        259    0.7303008 0.7303008
#> 358        x1     right -0.234990753        260    0.7278266 0.7278266
#> 359        x1     right -0.236735703        261    0.7252033 0.7252033
#> 360        x1     right -0.240254817        262    0.7232669 0.7232669
#> 361        x1     right -0.241446550        263    0.7206107 0.7206107
#> 362        x1     right -0.242481119        264    0.7174983 0.7174983
#> 363        x1     right -0.256408526        265    0.7158225 0.7158225
#> 364        x1     right -0.258232925        266    0.7123948 0.7123948
#> 365        x1     right -0.259285868        267    0.7094842 0.7094842
#> 366        x1     right -0.261832752        270    0.6996305 0.6996305
#> 367        x1     right -0.266872767        271    0.6967918 0.6967918
#> 368        x1     right -0.275403433        272    0.6928818 0.6928818
#> 369        x1     right -0.276211630        273    0.6912014 0.6912014
#> 370        x1     right -0.285188510        274    0.6882010 0.6882010
#> 371        x1     right -0.286717187        275    0.6861853 0.6861853
#> 372        x1     right -0.298748860        276    0.6834599 0.6834599
#> 373        x1     right -0.302518382        277    0.6815979 0.6815979
#> 374        x1     right -0.303551198        280    0.6736103 0.6736103
#> 375        x1     right -0.307675800        281    0.6711122 0.6711122
#> 376        x1     right -0.311025213        282    0.6681144 0.6681144
#> 377        x1     right -0.314513051        283    0.6659287 0.6659287
#> 378        x1     right -0.316380198        284    0.6646788 0.6646788
#> 379        x1     right -0.319534477        285    0.6625953 0.6625953
#> 380        x1     right -0.326243382        286    0.6591804 0.6591804
#> 381        x1     right -0.327109284        287    0.6577808 0.6577808
#> 382        x1     right -0.329358480        288    0.6542874 0.6542874
#> 383        x1     right -0.332788970        289    0.6521201 0.6521201
#> 384        x1     right -0.342512567        290    0.6503461 0.6503461
#> 385        x1     right -0.350269719        291    0.6484086 0.6484086
#> 386        x1     right -0.368243176        292    0.6462284 0.6462284
#> 387        x1     right -0.371596340        293    0.6440110 0.6440110
#> 388        x1     right -0.383895102        294    0.6414617 0.6414617
#> 389        x1     right -0.396613459        295    0.6393226 0.6393226
#> 390        x1     right -0.403676941        296    0.6363541 0.6363541
#> 391        x1     right -0.404611086        297    0.6347935 0.6347935
#> 392        x1     right -0.412016402        298    0.6325810 0.6325810
#> 393        x1     right -0.414102554        299    0.6298844 0.6298844
#> 394        x1     right -0.415337210        300    0.6270396 0.6270396
#> 395        x1     right -0.433045581        301    0.6244817 0.6244817
#> 396        x1     right -0.438499573        302    0.6215394 0.6215394
#> 397        x1     right -0.444683827        303    0.6190742 0.6190742
#> 398        x1     right -0.446647021        304    0.6160231 0.6160231
#> 399        x1     right -0.457255803        305    0.6129954 0.6129954
#> 400        x1     right -0.471501901        306    0.6092956 0.6092956
#> 401        x1     right -0.474925681        307    0.6067562 0.6067562
#> 402        x1     right -0.478140526        308    0.6050647 0.6050647
#> 403        x1     right -0.488035503        309    0.6028599 0.6028599
#> 404        x1     right -0.495980405        310    0.6005054 0.6005054
#> 405        x1     right -0.514086782        311    0.5981059 0.5981059
#> 406        x1     right -0.516841555        312    0.5952719 0.5952719
#> 407        x1     right -0.519767454        313    0.5927619 0.5927619
#> 408        x1     right -0.525470293        314    0.5902259 0.5902259
#> 409        x1     right -0.549733194        315    0.5867745 0.5867745
#> 410        x1     right -0.550738982        316    0.5843464 0.5843464
#> 411        x1     right -0.574709308        317    0.5804654 0.5804654
#> 412        x1     right -0.582732049        318    0.5783131 0.5783131
#> 413        x1     right -0.587708788        319    0.5763021 0.5763021
#> 414        x1     right -0.589538136        320    0.5722199 0.5722199
#> 415        x1     right -0.596033965        321    0.5688145 0.5688145
#> 416        x1     right -0.601115634        322    0.5656514 0.5656514
#> 417        x1     right -0.605706573        323    0.5630004 0.5630004
#> 418        x1     right -0.612431857        324    0.5592134 0.5592134
#> 419        x1     right -0.615404607        325    0.5559094 0.5559094
#> 420        x1     right -0.620440165        326    0.5528244 0.5528244
#> 421        x1     right -0.625147258        327    0.5493572 0.5493572
#> 422        x1     right -0.629860071        328    0.5465641 0.5465641
#> 423        x1     right -0.644100634        329    0.5433414 0.5433414
#> 424        x1     right -0.644676496        330    0.5404592 0.5404592
#> 425        x1     right -0.652694893        331    0.5360167 0.5360167
#> 426        x1     right -0.659512840        332    0.5315440 0.5315440
#> 427        x1     right -0.664569484        333    0.5277379 0.5277379
#> 428        x1     right -0.669169755        334    0.5249349 0.5249349
#> 429        x1     right -0.670217707        335    0.5214491 0.5214491
#> 430        x1     right -0.671346851        336    0.5191316 0.5191316
#> 431        x1     right -0.673829043        337    0.5151635 0.5151635
#> 432        x1     right -0.675781836        338    0.5134753 0.5134753
#> 433        x1     right -0.683208087        339    0.5100155 0.5100155
#> 434        x1     right -0.683560451        340    0.5068532 0.5068532
#> 435        x1     right -0.689197555        341    0.5030439 0.5030439
#> 436        x1     right -0.694423717        342    0.4988086 0.4988086
#> 437        x1     right -0.705037910        343    0.4951515 0.4951515
#> 438        x1     right -0.706876314        344    0.4919536 0.4919536
#> 439        x1     right -0.712425439        345    0.4894117 0.4894117
#> 440        x1     right -0.713666827        346    0.4854348 0.4854348
#> 441        x1     right -0.727024828        347    0.4823271 0.4823271
#> 442        x1     right -0.731777324        348    0.4797268 0.4797268
#> 443        x1     right -0.734344374        349    0.4761192 0.4761192
#> 444        x1     right -0.742870647        350    0.4730048 0.4730048
#> 445        x1     right -0.752552870        352    0.4671464 0.4671464
#> 446        x1     right -0.757057007        354    0.4615241 0.4615241
#> 447        x1     right -0.759435488        355    0.4584188 0.4584188
#> 448        x1     right -0.764558776        356    0.4551421 0.4551421
#> 449        x1     right -0.765604900        357    0.4514300 0.4514300
#> 450        x1     right -0.769315350        358    0.4490906 0.4490906
#> 451        x1     right -0.786703361        359    0.4463928 0.4463928
#> 452        x1     right -0.789899723        360    0.4442357 0.4442357
#> 453        x1     right -0.804382925        361    0.4412914 0.4412914
#> 454        x1     right -0.809695733        362    0.4392887 0.4392887
#> 455        x1     right -0.811219465        363    0.4358571 0.4358571
#> 456        x1     right -0.811908222        364    0.4321903 0.4321903
#> 457        x1     right -0.827209250        365    0.4292846 0.4292846
#> 458        x1     right -0.829377988        366    0.4272526 0.4272526
#> 459        x1     right -0.841120614        367    0.4248487 0.4248487
#> 460        x1     right -0.846205924        368    0.4221402 0.4221402
#> 461        x1     right -0.853506261        369    0.4190534 0.4190534
#> 462        x1     right -0.855078733        370    0.4157340 0.4157340
#> 463        x1     right -0.859934642        371    0.4135217 0.4135217
#> 464        x1     right -0.870100162        372    0.4114243 0.4114243
#> 465        x1     right -0.883371215        374    0.4035955 0.4035955
#> 466        x1     right -0.884135864        375    0.3994916 0.3994916
#> 467        x1     right -0.885254643        376    0.3969020 0.3969020
#> 468        x1     right -0.889671427        377    0.3949933 0.3949933
#> 469        x1     right -0.891619141        378    0.3924438 0.3924438
#> 470        x1     right -0.895701973        379    0.3889630 0.3889630
#> 471        x1     right -0.896809354        380    0.3854500 0.3854500
#> 472        x1     right -0.898457671        381    0.3831813 0.3831813
#> 473        x1     right -0.900692832        382    0.3811513 0.3811513
#> 474        x1     right -0.908721582        383    0.3785685 0.3785685
#> 475        x1     right -0.913563482        384    0.3750653 0.3750653
#> 476        x1     right -0.914685762        385    0.3721323 0.3721323
#> 477        x1     right -0.916943142        386    0.3687421 0.3687421
#> 478        x1     right -0.929444466        387    0.3661917 0.3661917
#> 479        x1     right -0.938486858        388    0.3630919 0.3630919
#> 480        x1     right -0.941245669        389    0.3589732 0.3589732
#> 481        x1     right -0.943214134        390    0.3554045 0.3554045
#> 482        x1     right -0.945004242        391    0.3522036 0.3522036
#> 483        x1     right -0.954517552        392    0.3504354 0.3504354
#> 484        x1     right -0.966775578        393    0.3472340 0.3472340
#> 485        x1     right -0.976480531        394    0.3441561 0.3441561
#> 486        x1     right -0.977700983        395    0.3413814 0.3413814
#> 487        x1     right -0.983743089        397    0.3354051 0.3354051
#> 488        x1     right -0.990467109        398    0.3328388 0.3328388
#> 489        x1     right -0.996134331        399    0.3298407 0.3298407
#> 490        x1      left  0.998910810        400    0.3263429 0.3263429
#>           cut_expr
#> 1      x1 >= 0.703
#> 2      x2 <= -0.68
#> 3     x3 <= -0.722
#> 4      x1 >= 0.702
#> 5     x2 <= -0.678
#> 6     x3 <= -0.719
#> 7      x1 >= 0.699
#> 8     x2 <= -0.677
#> 9     x3 <= -0.713
#> 10     x1 >= 0.692
#> 11    x2 <= -0.674
#> 12    x3 <= -0.602
#> 13     x1 >= 0.686
#> 14    x2 <= -0.672
#> 15    x3 <= -0.601
#> 16     x1 >= 0.685
#> 17    x2 <= -0.667
#> 18    x3 <= -0.553
#> 19     x1 >= 0.677
#> 20    x2 <= -0.652
#> 21    x3 <= -0.531
#> 22     x1 >= 0.671
#> 23    x2 <= -0.651
#> 24    x3 <= -0.521
#> 25      x1 >= 0.66
#> 26    x2 <= -0.644
#> 27     x3 <= -0.52
#> 28     x1 >= 0.659
#> 29    x2 <= -0.586
#> 30    x3 <= -0.518
#> 31     x1 >= 0.658
#> 32    x2 <= -0.584
#> 33    x3 <= -0.505
#> 34     x1 >= 0.657
#> 35    x2 <= -0.555
#> 36    x3 <= -0.493
#> 37     x1 >= 0.656
#> 38    x2 <= -0.479
#> 39    x3 <= -0.492
#> 40     x1 >= 0.655
#> 41    x2 <= -0.475
#> 42     x3 <= -0.49
#> 43     x1 >= 0.642
#> 44    x2 <= -0.471
#> 45    x3 <= -0.488
#> 46     x1 >= 0.641
#> 47    x2 <= -0.469
#> 48    x3 <= -0.404
#> 49      x1 >= 0.64
#> 50    x2 <= -0.349
#> 51    x3 <= -0.402
#> 52     x1 >= 0.629
#> 53    x2 <= -0.344
#> 54    x3 <= -0.398
#> 55     x1 >= 0.613
#> 56    x2 <= -0.343
#> 57    x3 <= -0.336
#> 58     x1 >= 0.609
#> 59    x2 <= -0.324
#> 60    x3 <= -0.267
#> 61     x1 >= 0.603
#> 62    x2 <= -0.312
#> 63    x3 <= -0.266
#> 64     x1 >= 0.598
#> 65    x2 <= -0.286
#> 66    x3 <= -0.257
#> 67     x1 >= 0.586
#> 68    x2 <= -0.266
#> 69    x3 <= -0.251
#> 70     x1 >= 0.585
#> 71    x2 <= -0.245
#> 72    x3 <= -0.244
#> 73     x1 >= 0.578
#> 74    x2 <= -0.231
#> 75    x3 <= -0.243
#> 76     x1 >= 0.566
#> 77     x2 <= -0.23
#> 78    x3 <= -0.239
#> 79     x1 >= 0.554
#> 80    x2 <= -0.227
#> 81    x3 <= -0.196
#> 82      x1 >= 0.54
#> 83    x2 <= -0.205
#> 84    x3 <= -0.195
#> 85     x1 >= 0.536
#> 86    x2 <= -0.194
#> 87    x3 <= -0.193
#> 88     x1 >= 0.517
#> 89    x2 <= -0.173
#> 90    x3 <= -0.188
#> 91     x1 >= 0.513
#> 92    x2 <= -0.168
#> 93    x3 <= -0.175
#> 94     x1 >= 0.507
#> 95    x2 <= -0.143
#> 96    x3 <= -0.174
#> 97     x1 >= 0.503
#> 98    x2 <= -0.132
#> 99    x3 <= -0.161
#> 100    x1 >= 0.502
#> 101   x2 <= -0.127
#> 102   x3 <= -0.157
#> 103    x1 >= 0.495
#> 104   x2 <= -0.123
#> 105   x3 <= -0.126
#> 106    x1 >= 0.479
#> 107   x2 <= -0.116
#> 108   x3 <= -0.113
#> 109    x1 >= 0.475
#> 110    x3 <= -0.11
#> 111  x2 <= -0.0753
#> 112    x1 >= 0.465
#> 113  x3 <= -0.0693
#> 114  x2 <= -0.0745
#> 115    x1 >= 0.463
#> 116  x3 <= -0.0625
#> 117  x2 <= -0.0612
#> 118    x1 >= 0.462
#> 119  x3 <= -0.0612
#> 120   x2 <= 0.0111
#> 121     x1 >= 0.46
#> 122  x3 <= -0.0601
#> 123   x2 <= 0.0118
#> 124    x1 >= 0.456
#> 125   x3 <= -0.029
#> 126    x2 <= 0.249
#> 127    x1 >= 0.449
#> 128  x3 <= -0.0112
#> 129    x2 <= 0.262
#> 130    x1 >= 0.447
#> 131 x3 <= -0.00872
#> 132    x2 <= 0.272
#> 133    x1 >= 0.433
#> 134   x3 <= 0.0692
#> 135    x2 <= 0.278
#> 136    x1 >= 0.432
#> 137   x3 <= 0.0813
#> 138     x2 <= 0.28
#> 139     x1 >= 0.43
#> 140    x3 <= 0.104
#> 141    x2 <= 0.381
#> 142    x1 >= 0.429
#> 143    x3 <= 0.111
#> 144    x2 <= 0.465
#> 145    x1 >= 0.425
#> 146    x3 <= 0.296
#> 147     x2 <= 0.47
#> 148    x1 >= 0.424
#> 149    x3 <= 0.299
#> 150     x2 <= 0.49
#> 151    x1 >= 0.418
#> 152    x3 <= 0.305
#> 153    x2 <= 0.503
#> 154    x1 >= 0.417
#> 155    x3 <= 0.355
#> 156    x2 <= 0.557
#> 157    x1 >= 0.413
#> 158    x3 <= 0.358
#> 159    x2 <= 0.558
#> 160    x1 >= 0.412
#> 161    x3 <= 0.465
#> 162     x2 <= 0.57
#> 163    x1 >= 0.392
#> 164    x3 <= 0.479
#> 165    x2 <= 0.578
#> 166     x1 >= 0.39
#> 167    x3 <= 0.484
#> 168     x2 <= 0.58
#> 169    x1 >= 0.384
#> 170    x3 <= 0.488
#> 171    x2 <= 0.602
#> 172    x1 >= 0.381
#> 173    x3 <= 0.495
#> 174    x2 <= 0.604
#> 175    x1 >= 0.379
#> 176    x3 <= 0.501
#> 177    x2 <= 0.628
#> 178    x1 >= 0.377
#> 179    x3 <= 0.508
#> 180    x2 <= 0.631
#> 181     x1 >= 0.37
#> 182    x3 <= 0.509
#> 183    x2 <= 0.632
#> 184    x1 >= 0.353
#> 185     x3 <= 0.51
#> 186    x2 <= 0.634
#> 187    x1 >= 0.352
#> 188     x3 <= 0.52
#> 189    x2 <= 0.679
#> 190    x1 >= 0.343
#> 191    x3 <= 0.524
#> 192    x2 <= 0.686
#> 193    x1 >= 0.336
#> 194    x3 <= 0.526
#> 195    x2 <= 0.717
#> 196    x1 >= 0.335
#> 197     x3 <= 0.54
#> 198    x2 <= 0.727
#> 199    x1 >= 0.327
#> 200    x3 <= 0.585
#> 201    x2 <= 0.729
#> 202    x1 >= 0.326
#> 203     x3 <= 0.59
#> 204    x2 <= 0.739
#> 205    x1 >= 0.318
#> 206    x3 <= 0.612
#> 207    x2 <= 0.749
#> 208    x1 >= 0.316
#> 209    x3 <= 0.627
#> 210    x2 <= 0.754
#> 211    x1 >= 0.309
#> 212    x3 <= 0.645
#> 213    x2 <= 0.758
#> 214    x1 >= 0.307
#> 215    x3 <= 0.758
#> 216    x2 <= 0.762
#> 217    x1 >= 0.304
#> 218    x3 <= 0.763
#> 219    x2 <= 0.766
#> 220    x1 >= 0.297
#> 221      x3 <= 0.8
#> 222    x2 <= 0.786
#> 223    x1 >= 0.296
#> 224    x3 <= 0.805
#> 225    x2 <= 0.804
#> 226    x1 >= 0.293
#> 227    x3 <= 0.807
#> 228    x2 <= 0.813
#> 229    x1 >= 0.292
#> 230    x3 <= 0.814
#> 231    x2 <= 0.842
#> 232    x1 >= 0.288
#> 233    x3 <= 0.867
#> 234     x2 <= 0.85
#> 235    x1 >= 0.276
#> 236    x3 <= 0.869
#> 237    x2 <= 0.851
#> 238     x1 >= 0.27
#> 239     x3 <= 0.87
#> 240    x2 <= 0.947
#> 241    x1 >= 0.253
#> 242    x3 <= 0.921
#> 243    x2 <= 0.957
#> 244     x1 >= 0.25
#> 245    x3 <= 0.922
#> 246    x2 <= 0.961
#> 247    x1 >= 0.247
#> 248    x3 <= 0.924
#> 249    x2 <= 0.978
#> 250    x1 >= 0.244
#> 251    x3 <= 0.928
#> 252   x2 >= -0.996
#> 253    x1 >= 0.242
#> 254     x3 <= 0.95
#> 255        x2 <= 1
#> 256    x1 >= 0.241
#> 257    x3 <= 0.951
#> 258   x2 >= -0.997
#> 259    x1 >= 0.219
#> 260    x3 <= 0.955
#> 261    x1 >= 0.213
#> 262    x3 <= 0.957
#> 263    x1 >= 0.208
#> 264    x3 <= 0.969
#> 265    x1 >= 0.207
#> 266    x3 <= 0.993
#> 267    x1 >= 0.197
#> 268   x3 >= -0.997
#> 269    x1 >= 0.194
#> 270    x1 >= 0.193
#> 271    x1 >= 0.177
#> 272    x1 >= 0.176
#> 273    x1 >= 0.157
#> 274    x1 >= 0.156
#> 275    x1 >= 0.144
#> 276    x1 >= 0.135
#> 277    x1 >= 0.131
#> 278    x1 >= 0.124
#> 279    x1 >= 0.119
#> 280    x1 >= 0.108
#> 281    x1 >= 0.107
#> 282    x1 >= 0.103
#> 283   x1 >= 0.0963
#> 284   x1 >= 0.0847
#> 285    x1 >= 0.084
#> 286   x1 >= 0.0789
#> 287   x1 >= 0.0776
#> 288   x1 >= 0.0757
#> 289   x1 >= 0.0661
#> 290   x1 >= 0.0652
#> 291   x1 >= 0.0616
#> 292   x1 >= 0.0561
#> 293   x1 >= 0.0542
#> 294   x1 >= 0.0506
#> 295   x1 >= 0.0445
#> 296   x1 >= 0.0359
#> 297    x1 >= 0.034
#> 298   x1 >= 0.0322
#> 299   x1 >= 0.0311
#> 300   x1 >= 0.0295
#> 301   x1 >= 0.0147
#> 302   x1 >= 0.0117
#> 303   x1 >= 0.0045
#> 304  x1 >= 0.00351
#> 305   x1 >= 0.0017
#> 306 x1 >= -0.00426
#> 307 x1 >= -0.00663
#> 308 x1 >= -0.00987
#> 309  x1 >= -0.0107
#> 310  x1 >= -0.0287
#> 311  x1 >= -0.0357
#> 312  x1 >= -0.0403
#> 313  x1 >= -0.0451
#> 314  x1 >= -0.0476
#> 315  x1 >= -0.0483
#> 316  x1 >= -0.0492
#> 317   x1 >= -0.053
#> 318  x1 >= -0.0549
#> 319   x1 >= -0.055
#> 320  x1 >= -0.0555
#> 321  x1 >= -0.0597
#> 322  x1 >= -0.0658
#> 323  x1 >= -0.0782
#> 324  x1 >= -0.0808
#> 325  x1 >= -0.0877
#> 326  x1 >= -0.0963
#> 327   x1 >= -0.099
#> 328   x1 >= -0.103
#> 329   x1 >= -0.115
#> 330   x1 >= -0.116
#> 331   x1 >= -0.117
#> 332   x1 >= -0.124
#> 333   x1 >= -0.128
#> 334    x1 >= -0.13
#> 335   x1 >= -0.131
#> 336   x1 >= -0.132
#> 337   x1 >= -0.134
#> 338   x1 >= -0.137
#> 339   x1 >= -0.143
#> 340   x1 >= -0.148
#> 341   x1 >= -0.152
#> 342   x1 >= -0.153
#> 343   x1 >= -0.154
#> 344   x1 >= -0.155
#> 345    x1 >= -0.16
#> 346   x1 >= -0.163
#> 347   x1 >= -0.169
#> 348   x1 >= -0.183
#> 349   x1 >= -0.188
#> 350   x1 >= -0.196
#> 351   x1 >= -0.202
#> 352   x1 >= -0.204
#> 353   x1 >= -0.206
#> 354   x1 >= -0.212
#> 355   x1 >= -0.215
#> 356   x1 >= -0.224
#> 357   x1 >= -0.228
#> 358   x1 >= -0.235
#> 359   x1 >= -0.237
#> 360    x1 >= -0.24
#> 361   x1 >= -0.241
#> 362   x1 >= -0.242
#> 363   x1 >= -0.256
#> 364   x1 >= -0.258
#> 365   x1 >= -0.259
#> 366   x1 >= -0.262
#> 367   x1 >= -0.267
#> 368   x1 >= -0.275
#> 369   x1 >= -0.276
#> 370   x1 >= -0.285
#> 371   x1 >= -0.287
#> 372   x1 >= -0.299
#> 373   x1 >= -0.303
#> 374   x1 >= -0.304
#> 375   x1 >= -0.308
#> 376   x1 >= -0.311
#> 377   x1 >= -0.315
#> 378   x1 >= -0.316
#> 379    x1 >= -0.32
#> 380   x1 >= -0.326
#> 381   x1 >= -0.327
#> 382   x1 >= -0.329
#> 383   x1 >= -0.333
#> 384   x1 >= -0.343
#> 385    x1 >= -0.35
#> 386   x1 >= -0.368
#> 387   x1 >= -0.372
#> 388   x1 >= -0.384
#> 389   x1 >= -0.397
#> 390   x1 >= -0.404
#> 391   x1 >= -0.405
#> 392   x1 >= -0.412
#> 393   x1 >= -0.414
#> 394   x1 >= -0.415
#> 395   x1 >= -0.433
#> 396   x1 >= -0.438
#> 397   x1 >= -0.445
#> 398   x1 >= -0.447
#> 399   x1 >= -0.457
#> 400   x1 >= -0.472
#> 401   x1 >= -0.475
#> 402   x1 >= -0.478
#> 403   x1 >= -0.488
#> 404   x1 >= -0.496
#> 405   x1 >= -0.514
#> 406   x1 >= -0.517
#> 407    x1 >= -0.52
#> 408   x1 >= -0.525
#> 409    x1 >= -0.55
#> 410   x1 >= -0.551
#> 411   x1 >= -0.575
#> 412   x1 >= -0.583
#> 413   x1 >= -0.588
#> 414    x1 >= -0.59
#> 415   x1 >= -0.596
#> 416   x1 >= -0.601
#> 417   x1 >= -0.606
#> 418   x1 >= -0.612
#> 419   x1 >= -0.615
#> 420    x1 >= -0.62
#> 421   x1 >= -0.625
#> 422    x1 >= -0.63
#> 423   x1 >= -0.644
#> 424   x1 >= -0.645
#> 425   x1 >= -0.653
#> 426    x1 >= -0.66
#> 427   x1 >= -0.665
#> 428   x1 >= -0.669
#> 429    x1 >= -0.67
#> 430   x1 >= -0.671
#> 431   x1 >= -0.674
#> 432   x1 >= -0.676
#> 433   x1 >= -0.683
#> 434   x1 >= -0.684
#> 435   x1 >= -0.689
#> 436   x1 >= -0.694
#> 437   x1 >= -0.705
#> 438   x1 >= -0.707
#> 439   x1 >= -0.712
#> 440   x1 >= -0.714
#> 441   x1 >= -0.727
#> 442   x1 >= -0.732
#> 443   x1 >= -0.734
#> 444   x1 >= -0.743
#> 445   x1 >= -0.753
#> 446   x1 >= -0.757
#> 447   x1 >= -0.759
#> 448   x1 >= -0.765
#> 449   x1 >= -0.766
#> 450   x1 >= -0.769
#> 451   x1 >= -0.787
#> 452    x1 >= -0.79
#> 453   x1 >= -0.804
#> 454    x1 >= -0.81
#> 455   x1 >= -0.811
#> 456   x1 >= -0.812
#> 457   x1 >= -0.827
#> 458   x1 >= -0.829
#> 459   x1 >= -0.841
#> 460   x1 >= -0.846
#> 461   x1 >= -0.854
#> 462   x1 >= -0.855
#> 463    x1 >= -0.86
#> 464    x1 >= -0.87
#> 465   x1 >= -0.883
#> 466   x1 >= -0.884
#> 467   x1 >= -0.885
#> 468    x1 >= -0.89
#> 469   x1 >= -0.892
#> 470   x1 >= -0.896
#> 471   x1 >= -0.897
#> 472   x1 >= -0.898
#> 473   x1 >= -0.901
#> 474   x1 >= -0.909
#> 475   x1 >= -0.914
#> 476   x1 >= -0.915
#> 477   x1 >= -0.917
#> 478   x1 >= -0.929
#> 479   x1 >= -0.938
#> 480   x1 >= -0.941
#> 481   x1 >= -0.943
#> 482   x1 >= -0.945
#> 483   x1 >= -0.955
#> 484   x1 >= -0.967
#> 485   x1 >= -0.976
#> 486   x1 >= -0.978
#> 487   x1 >= -0.984
#> 488    x1 >= -0.99
#> 489   x1 >= -0.996
#> 490    x1 <= 0.999
fr$cut_expr   # ready to append to forestsearch(conf_force = ...)
#>   [1] "x1 >= 0.703"    "x2 <= -0.68"    "x3 <= -0.722"   "x1 >= 0.702"   
#>   [5] "x2 <= -0.678"   "x3 <= -0.719"   "x1 >= 0.699"    "x2 <= -0.677"  
#>   [9] "x3 <= -0.713"   "x1 >= 0.692"    "x2 <= -0.674"   "x3 <= -0.602"  
#>  [13] "x1 >= 0.686"    "x2 <= -0.672"   "x3 <= -0.601"   "x1 >= 0.685"   
#>  [17] "x2 <= -0.667"   "x3 <= -0.553"   "x1 >= 0.677"    "x2 <= -0.652"  
#>  [21] "x3 <= -0.531"   "x1 >= 0.671"    "x2 <= -0.651"   "x3 <= -0.521"  
#>  [25] "x1 >= 0.66"     "x2 <= -0.644"   "x3 <= -0.52"    "x1 >= 0.659"   
#>  [29] "x2 <= -0.586"   "x3 <= -0.518"   "x1 >= 0.658"    "x2 <= -0.584"  
#>  [33] "x3 <= -0.505"   "x1 >= 0.657"    "x2 <= -0.555"   "x3 <= -0.493"  
#>  [37] "x1 >= 0.656"    "x2 <= -0.479"   "x3 <= -0.492"   "x1 >= 0.655"   
#>  [41] "x2 <= -0.475"   "x3 <= -0.49"    "x1 >= 0.642"    "x2 <= -0.471"  
#>  [45] "x3 <= -0.488"   "x1 >= 0.641"    "x2 <= -0.469"   "x3 <= -0.404"  
#>  [49] "x1 >= 0.64"     "x2 <= -0.349"   "x3 <= -0.402"   "x1 >= 0.629"   
#>  [53] "x2 <= -0.344"   "x3 <= -0.398"   "x1 >= 0.613"    "x2 <= -0.343"  
#>  [57] "x3 <= -0.336"   "x1 >= 0.609"    "x2 <= -0.324"   "x3 <= -0.267"  
#>  [61] "x1 >= 0.603"    "x2 <= -0.312"   "x3 <= -0.266"   "x1 >= 0.598"   
#>  [65] "x2 <= -0.286"   "x3 <= -0.257"   "x1 >= 0.586"    "x2 <= -0.266"  
#>  [69] "x3 <= -0.251"   "x1 >= 0.585"    "x2 <= -0.245"   "x3 <= -0.244"  
#>  [73] "x1 >= 0.578"    "x2 <= -0.231"   "x3 <= -0.243"   "x1 >= 0.566"   
#>  [77] "x2 <= -0.23"    "x3 <= -0.239"   "x1 >= 0.554"    "x2 <= -0.227"  
#>  [81] "x3 <= -0.196"   "x1 >= 0.54"     "x2 <= -0.205"   "x3 <= -0.195"  
#>  [85] "x1 >= 0.536"    "x2 <= -0.194"   "x3 <= -0.193"   "x1 >= 0.517"   
#>  [89] "x2 <= -0.173"   "x3 <= -0.188"   "x1 >= 0.513"    "x2 <= -0.168"  
#>  [93] "x3 <= -0.175"   "x1 >= 0.507"    "x2 <= -0.143"   "x3 <= -0.174"  
#>  [97] "x1 >= 0.503"    "x2 <= -0.132"   "x3 <= -0.161"   "x1 >= 0.502"   
#> [101] "x2 <= -0.127"   "x3 <= -0.157"   "x1 >= 0.495"    "x2 <= -0.123"  
#> [105] "x3 <= -0.126"   "x1 >= 0.479"    "x2 <= -0.116"   "x3 <= -0.113"  
#> [109] "x1 >= 0.475"    "x3 <= -0.11"    "x2 <= -0.0753"  "x1 >= 0.465"   
#> [113] "x3 <= -0.0693"  "x2 <= -0.0745"  "x1 >= 0.463"    "x3 <= -0.0625" 
#> [117] "x2 <= -0.0612"  "x1 >= 0.462"    "x3 <= -0.0612"  "x2 <= 0.0111"  
#> [121] "x1 >= 0.46"     "x3 <= -0.0601"  "x2 <= 0.0118"   "x1 >= 0.456"   
#> [125] "x3 <= -0.029"   "x2 <= 0.249"    "x1 >= 0.449"    "x3 <= -0.0112" 
#> [129] "x2 <= 0.262"    "x1 >= 0.447"    "x3 <= -0.00872" "x2 <= 0.272"   
#> [133] "x1 >= 0.433"    "x3 <= 0.0692"   "x2 <= 0.278"    "x1 >= 0.432"   
#> [137] "x3 <= 0.0813"   "x2 <= 0.28"     "x1 >= 0.43"     "x3 <= 0.104"   
#> [141] "x2 <= 0.381"    "x1 >= 0.429"    "x3 <= 0.111"    "x2 <= 0.465"   
#> [145] "x1 >= 0.425"    "x3 <= 0.296"    "x2 <= 0.47"     "x1 >= 0.424"   
#> [149] "x3 <= 0.299"    "x2 <= 0.49"     "x1 >= 0.418"    "x3 <= 0.305"   
#> [153] "x2 <= 0.503"    "x1 >= 0.417"    "x3 <= 0.355"    "x2 <= 0.557"   
#> [157] "x1 >= 0.413"    "x3 <= 0.358"    "x2 <= 0.558"    "x1 >= 0.412"   
#> [161] "x3 <= 0.465"    "x2 <= 0.57"     "x1 >= 0.392"    "x3 <= 0.479"   
#> [165] "x2 <= 0.578"    "x1 >= 0.39"     "x3 <= 0.484"    "x2 <= 0.58"    
#> [169] "x1 >= 0.384"    "x3 <= 0.488"    "x2 <= 0.602"    "x1 >= 0.381"   
#> [173] "x3 <= 0.495"    "x2 <= 0.604"    "x1 >= 0.379"    "x3 <= 0.501"   
#> [177] "x2 <= 0.628"    "x1 >= 0.377"    "x3 <= 0.508"    "x2 <= 0.631"   
#> [181] "x1 >= 0.37"     "x3 <= 0.509"    "x2 <= 0.632"    "x1 >= 0.353"   
#> [185] "x3 <= 0.51"     "x2 <= 0.634"    "x1 >= 0.352"    "x3 <= 0.52"    
#> [189] "x2 <= 0.679"    "x1 >= 0.343"    "x3 <= 0.524"    "x2 <= 0.686"   
#> [193] "x1 >= 0.336"    "x3 <= 0.526"    "x2 <= 0.717"    "x1 >= 0.335"   
#> [197] "x3 <= 0.54"     "x2 <= 0.727"    "x1 >= 0.327"    "x3 <= 0.585"   
#> [201] "x2 <= 0.729"    "x1 >= 0.326"    "x3 <= 0.59"     "x2 <= 0.739"   
#> [205] "x1 >= 0.318"    "x3 <= 0.612"    "x2 <= 0.749"    "x1 >= 0.316"   
#> [209] "x3 <= 0.627"    "x2 <= 0.754"    "x1 >= 0.309"    "x3 <= 0.645"   
#> [213] "x2 <= 0.758"    "x1 >= 0.307"    "x3 <= 0.758"    "x2 <= 0.762"   
#> [217] "x1 >= 0.304"    "x3 <= 0.763"    "x2 <= 0.766"    "x1 >= 0.297"   
#> [221] "x3 <= 0.8"      "x2 <= 0.786"    "x1 >= 0.296"    "x3 <= 0.805"   
#> [225] "x2 <= 0.804"    "x1 >= 0.293"    "x3 <= 0.807"    "x2 <= 0.813"   
#> [229] "x1 >= 0.292"    "x3 <= 0.814"    "x2 <= 0.842"    "x1 >= 0.288"   
#> [233] "x3 <= 0.867"    "x2 <= 0.85"     "x1 >= 0.276"    "x3 <= 0.869"   
#> [237] "x2 <= 0.851"    "x1 >= 0.27"     "x3 <= 0.87"     "x2 <= 0.947"   
#> [241] "x1 >= 0.253"    "x3 <= 0.921"    "x2 <= 0.957"    "x1 >= 0.25"    
#> [245] "x3 <= 0.922"    "x2 <= 0.961"    "x1 >= 0.247"    "x3 <= 0.924"   
#> [249] "x2 <= 0.978"    "x1 >= 0.244"    "x3 <= 0.928"    "x2 >= -0.996"  
#> [253] "x1 >= 0.242"    "x3 <= 0.95"     "x2 <= 1"        "x1 >= 0.241"   
#> [257] "x3 <= 0.951"    "x2 >= -0.997"   "x1 >= 0.219"    "x3 <= 0.955"   
#> [261] "x1 >= 0.213"    "x3 <= 0.957"    "x1 >= 0.208"    "x3 <= 0.969"   
#> [265] "x1 >= 0.207"    "x3 <= 0.993"    "x1 >= 0.197"    "x3 >= -0.997"  
#> [269] "x1 >= 0.194"    "x1 >= 0.193"    "x1 >= 0.177"    "x1 >= 0.176"   
#> [273] "x1 >= 0.157"    "x1 >= 0.156"    "x1 >= 0.144"    "x1 >= 0.135"   
#> [277] "x1 >= 0.131"    "x1 >= 0.124"    "x1 >= 0.119"    "x1 >= 0.108"   
#> [281] "x1 >= 0.107"    "x1 >= 0.103"    "x1 >= 0.0963"   "x1 >= 0.0847"  
#> [285] "x1 >= 0.084"    "x1 >= 0.0789"   "x1 >= 0.0776"   "x1 >= 0.0757"  
#> [289] "x1 >= 0.0661"   "x1 >= 0.0652"   "x1 >= 0.0616"   "x1 >= 0.0561"  
#> [293] "x1 >= 0.0542"   "x1 >= 0.0506"   "x1 >= 0.0445"   "x1 >= 0.0359"  
#> [297] "x1 >= 0.034"    "x1 >= 0.0322"   "x1 >= 0.0311"   "x1 >= 0.0295"  
#> [301] "x1 >= 0.0147"   "x1 >= 0.0117"   "x1 >= 0.0045"   "x1 >= 0.00351" 
#> [305] "x1 >= 0.0017"   "x1 >= -0.00426" "x1 >= -0.00663" "x1 >= -0.00987"
#> [309] "x1 >= -0.0107"  "x1 >= -0.0287"  "x1 >= -0.0357"  "x1 >= -0.0403" 
#> [313] "x1 >= -0.0451"  "x1 >= -0.0476"  "x1 >= -0.0483"  "x1 >= -0.0492" 
#> [317] "x1 >= -0.053"   "x1 >= -0.0549"  "x1 >= -0.055"   "x1 >= -0.0555" 
#> [321] "x1 >= -0.0597"  "x1 >= -0.0658"  "x1 >= -0.0782"  "x1 >= -0.0808" 
#> [325] "x1 >= -0.0877"  "x1 >= -0.0963"  "x1 >= -0.099"   "x1 >= -0.103"  
#> [329] "x1 >= -0.115"   "x1 >= -0.116"   "x1 >= -0.117"   "x1 >= -0.124"  
#> [333] "x1 >= -0.128"   "x1 >= -0.13"    "x1 >= -0.131"   "x1 >= -0.132"  
#> [337] "x1 >= -0.134"   "x1 >= -0.137"   "x1 >= -0.143"   "x1 >= -0.148"  
#> [341] "x1 >= -0.152"   "x1 >= -0.153"   "x1 >= -0.154"   "x1 >= -0.155"  
#> [345] "x1 >= -0.16"    "x1 >= -0.163"   "x1 >= -0.169"   "x1 >= -0.183"  
#> [349] "x1 >= -0.188"   "x1 >= -0.196"   "x1 >= -0.202"   "x1 >= -0.204"  
#> [353] "x1 >= -0.206"   "x1 >= -0.212"   "x1 >= -0.215"   "x1 >= -0.224"  
#> [357] "x1 >= -0.228"   "x1 >= -0.235"   "x1 >= -0.237"   "x1 >= -0.24"   
#> [361] "x1 >= -0.241"   "x1 >= -0.242"   "x1 >= -0.256"   "x1 >= -0.258"  
#> [365] "x1 >= -0.259"   "x1 >= -0.262"   "x1 >= -0.267"   "x1 >= -0.275"  
#> [369] "x1 >= -0.276"   "x1 >= -0.285"   "x1 >= -0.287"   "x1 >= -0.299"  
#> [373] "x1 >= -0.303"   "x1 >= -0.304"   "x1 >= -0.308"   "x1 >= -0.311"  
#> [377] "x1 >= -0.315"   "x1 >= -0.316"   "x1 >= -0.32"    "x1 >= -0.326"  
#> [381] "x1 >= -0.327"   "x1 >= -0.329"   "x1 >= -0.333"   "x1 >= -0.343"  
#> [385] "x1 >= -0.35"    "x1 >= -0.368"   "x1 >= -0.372"   "x1 >= -0.384"  
#> [389] "x1 >= -0.397"   "x1 >= -0.404"   "x1 >= -0.405"   "x1 >= -0.412"  
#> [393] "x1 >= -0.414"   "x1 >= -0.415"   "x1 >= -0.433"   "x1 >= -0.438"  
#> [397] "x1 >= -0.445"   "x1 >= -0.447"   "x1 >= -0.457"   "x1 >= -0.472"  
#> [401] "x1 >= -0.475"   "x1 >= -0.478"   "x1 >= -0.488"   "x1 >= -0.496"  
#> [405] "x1 >= -0.514"   "x1 >= -0.517"   "x1 >= -0.52"    "x1 >= -0.525"  
#> [409] "x1 >= -0.55"    "x1 >= -0.551"   "x1 >= -0.575"   "x1 >= -0.583"  
#> [413] "x1 >= -0.588"   "x1 >= -0.59"    "x1 >= -0.596"   "x1 >= -0.601"  
#> [417] "x1 >= -0.606"   "x1 >= -0.612"   "x1 >= -0.615"   "x1 >= -0.62"   
#> [421] "x1 >= -0.625"   "x1 >= -0.63"    "x1 >= -0.644"   "x1 >= -0.645"  
#> [425] "x1 >= -0.653"   "x1 >= -0.66"    "x1 >= -0.665"   "x1 >= -0.669"  
#> [429] "x1 >= -0.67"    "x1 >= -0.671"   "x1 >= -0.674"   "x1 >= -0.676"  
#> [433] "x1 >= -0.683"   "x1 >= -0.684"   "x1 >= -0.689"   "x1 >= -0.694"  
#> [437] "x1 >= -0.705"   "x1 >= -0.707"   "x1 >= -0.712"   "x1 >= -0.714"  
#> [441] "x1 >= -0.727"   "x1 >= -0.732"   "x1 >= -0.734"   "x1 >= -0.743"  
#> [445] "x1 >= -0.753"   "x1 >= -0.757"   "x1 >= -0.759"   "x1 >= -0.765"  
#> [449] "x1 >= -0.766"   "x1 >= -0.769"   "x1 >= -0.787"   "x1 >= -0.79"   
#> [453] "x1 >= -0.804"   "x1 >= -0.81"    "x1 >= -0.811"   "x1 >= -0.812"  
#> [457] "x1 >= -0.827"   "x1 >= -0.829"   "x1 >= -0.841"   "x1 >= -0.846"  
#> [461] "x1 >= -0.854"   "x1 >= -0.855"   "x1 >= -0.86"    "x1 >= -0.87"   
#> [465] "x1 >= -0.883"   "x1 >= -0.884"   "x1 >= -0.885"   "x1 >= -0.89"   
#> [469] "x1 >= -0.892"   "x1 >= -0.896"   "x1 >= -0.897"   "x1 >= -0.898"  
#> [473] "x1 >= -0.901"   "x1 >= -0.909"   "x1 >= -0.914"   "x1 >= -0.915"  
#> [477] "x1 >= -0.917"   "x1 >= -0.929"   "x1 >= -0.938"   "x1 >= -0.941"  
#> [481] "x1 >= -0.943"   "x1 >= -0.945"   "x1 >= -0.955"   "x1 >= -0.967"  
#> [485] "x1 >= -0.976"   "x1 >= -0.978"   "x1 >= -0.984"   "x1 >= -0.99"   
#> [489] "x1 >= -0.996"   "x1 <= 0.999"   
```
