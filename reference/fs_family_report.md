# Report which stages of the candidate family are data-dependent

A user who sets `use_lasso = FALSE`, `use_grf = FALSE`,
`use_dina = FALSE` and `vi.grf.min = NULL` may reasonably believe the
candidate family
[`forestsearch`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
searches is now fixed. It is not. Continuous cuts are placed at *sample*
quantiles
([`get_FSdata()`](https://larry-leon.github.io/forestsearch/reference/get_FSdata.md));
the prevalence, redundancy and size floors (`minp`, `rmin`, `n.min`) are
applied to *sample* counts
([`subgroup.search()`](https://larry-leon.github.io/forestsearch/reference/subgroup.search.md));
and the consistency stage removes near-duplicate candidates keyed on
their fitted *statistics*, with no argument that turns it off
([`remove_near_duplicate_subgroups()`](https://larry-leon.github.io/forestsearch/reference/remove_near_duplicate_subgroups.md)).
Those stages are the method, not screening. This function says so, for a
given argument set, in a form that cannot overpromise.

## Usage

``` r
fs_family_report(x, data = NULL, outcome_type = NULL)

# S3 method for class 'fs_family_report'
print(x, ...)
```

## Arguments

- x:

  Either a named list of
  [`forestsearch`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
  arguments (as one would pass to `do.call(forestsearch, x)`;
  unspecified arguments take
  [`forestsearch()`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)'s
  formal defaults), or a fitted `forestsearch` object, in which case its
  `args_call_all` is read.

- data:

  Optional data frame. When supplied, the report is grounded in counts:
  the number of cut columns
  [`get_FSdata()`](https://larry-leon.github.io/forestsearch/reference/get_FSdata.md)
  produces on this data under these arguments, and the number of
  up-to-`maxk` combinations that implies. Nothing is fitted; the search
  is not run.

- outcome_type:

  One of `"survival"`, `"binary"`, `"continuous"`, `"count"`. Taken from
  `x` when present there; required otherwise, because the per-arm floor
  row depends on it.

- ...:

  Ignored.

## Value

A data frame of class `c("fs_family_report", "data.frame")` with one row
per stage and columns `stage`, `arguments` (the governing
[`forestsearch()`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
formals), `values` (their resolved values, formatted), `status` (one of
`"deterministic"`, `"disabled"`, `"inert"`, `"data-dependent"`,
`"data-dependent (not disableable)"`) and `note` (one line: why, and
what would change it; `NA` where nothing would). Attributes: `verdict`
(one sentence), `status_counts` (a named integer vector),
`data_supplied` (logical), and – when `data` is supplied –
`n_cut_columns` and `n_combinations`.

## Details

**It reports and changes nothing.** It never fits a model, never runs
the search, and is called by nothing else in the package. **No
combination of arguments makes the family deterministic while cuts are
placed at sample quantiles**; the best a caller can do is switch off the
disableable stages, and the printed footer lists what remains.

The report mirrors
[`forestsearch()`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)'s
own argument resolution rather than reading arguments at face value:
`sg_focus` aliases are normalised; `sg_focus = "maxeff"` zeroes `minp`,
`rmin` and `pconsistency.threshold`, sets `stop_threshold` to `NULL`,
`use_twostage` to `FALSE` and `max_subgroups_search` to `Inf`, and
disables the effect floor; `stop_threshold` is reset to `NULL` for every
focus other than `"maxeffCons"` (it is meaningful only there);
`max_n_confounders` is applied only inside the GRF variable-importance
block, so it is inert when `vi.grf.min` is `NULL`; and the per-arm event
floors `d0.min` / `d1.min` are skipped for continuous and count
outcomes.

## See also

[`forestsearch`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md),
[`fs_oc_family_enumerate`](https://larry-leon.github.io/forestsearch/reference/fs_oc_family_enumerate.md)
(the population-frame enumeration used by the operating-characteristic
wrapper, which is deterministic in the DGM precisely because it does not
cut at sample quantiles).

## Examples

``` r
rep <- fs_family_report(
  list(confounders.name = c("age", "preanti", "hemo"),
       conf.cont_jcuts = list(age = 10, preanti = 10),
       n.min = 60, maxk = 2, sg_focus = "maxeffCons",
       use_lasso = FALSE, use_grf = FALSE, use_dina = FALSE,
       vi.grf.min = NULL),
  outcome_type = "continuous")
print(rep)
#> Candidate family is data-dependent: 8 of 20 stages vary with the sample; 4 are disableable, 4 are intrinsic to the method.
#> 
#>  stage                            status                          
#>  cut construction                 data-dependent (not disableable)
#>  LASSO screen                     disabled                        
#>  DINA / GRF subgroup paths        disabled                        
#>  GRF variable-importance ordering disabled                        
#>  confounder cap                   inert                           
#>  dummy expansion                  deterministic                   
#>  combination enumeration          deterministic                   
#>  per-factor prevalence floor      data-dependent                  
#>  redundancy (rmin)                data-dependent (not disableable)
#>  subgroup size floor              data-dependent (not disableable)
#>  per-arm floors                   inert                           
#>  estimability boundary            inert                           
#>  effect screen                    data-dependent                  
#>  near-duplicate removal           data-dependent (not disableable)
#>  candidate cap                    inert                           
#>  m1 filter                        inert                           
#>  consistency screen               data-dependent                  
#>  early stopping                   data-dependent                  
#>  winner selection                 deterministic                   
#>  time cap                         inert                           
#>  values                                                      
#>  conf.cont_jcuts = list(age = 10, preanti = 10); cut_type ...
#>  use_lasso = FALSE                                           
#>  use_dina = FALSE; use_grf = FALSE; subgroup_method = "con...
#>  vi.grf.min = NULL                                           
#>  max_n_confounders = 1000                                    
#>  each 2-level factor -> two indicator columns (both direct...
#>  maxk = 2                                                    
#>  minp = 0.025                                                
#>  rmin = 5 (subgroup.search() default)                        
#>  n.min = 60; n.min.frac = 0.1                                
#>  d0.min = 10; d1.min = 10                                    
#>  effect_measure = "MD" -> no existence condition             
#>  effect.threshold = NULL                                     
#>  sg_focus = "maxeffCons" -> statistics-keyed dedup           
#>  max_subgroups_search = Inf                                  
#>  m1.threshold = Inf                                          
#>  consistency_method = "resample"; pconsistency.threshold =...
#>  stop_threshold = 0.9                                        
#>  sg_focus = "maxeffCons" (canonical "maxeffCons"); selecti...
#>  max.minutes = 3                                             
#> 
#> Intrinsic to the method -- no argument switches these off:
#>   * cut construction: continuous cuts are placed at SAMPLE quantiles (get_FSdata(): get_conf_force() / cut_var_jq(), medians for cut_type = "median"); only a full set of cut VALUES supplied through conf_force with the defaults excluded would fix them, and the floors and near-duplicate removal below still vary
#>   * redundancy (rmin): each added factor must shrink the SAMPLE membership by more than rmin subjects, walking the columns in their current order (extract_idx_flagredundancy()); rmin = 0 (maxeff) still drops exact-membership prefixes, so it cannot be switched off, and rmin is not a forestsearch() formal
#>   * subgroup size floor: sum(id.x) > n.min on the SAMPLE (subgroup.search() status 4); a floor on a sample count with no off switch
#>   * near-duplicate removal: candidates whose (K, n, E, d1, m1, m0, HR, L, U) agree to 0.001 collapse to the first (remove_near_duplicate_subgroups(), subgroup_consistency_main.R L579); statistics are sample-fitted; NO argument disables it
attr(rep, "verdict")
#> [1] "Candidate family is data-dependent: 8 of 20 stages vary with the sample; 4 are disableable, 4 are intrinsic to the method."

# grounded in counts: cut columns and combinations on a data frame
set.seed(1)
df <- data.frame(age = round(rnorm(200, 40, 8)),
                 preanti = round(rexp(200, 1 / 500)),
                 hemo = rbinom(200, 1L, 0.1))
rep2 <- fs_family_report(
  list(confounders.name = c("age", "preanti", "hemo"),
       conf.cont_jcuts = list(age = 10, preanti = 10), maxk = 2,
       use_lasso = FALSE, use_grf = FALSE, vi.grf.min = NULL),
  data = df, outcome_type = "continuous")
attr(rep2, "n_cut_columns"); attr(rep2, "n_combinations")
#> [1] 42
#> [1] 903
```
