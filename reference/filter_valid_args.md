# Filter a Parameter List to a Function's Formals

Given a named list of user-supplied parameters and a target function,
returns only the parameter entries whose names appear in the target's
[`formals`](https://rdrr.io/r/base/formals.html). Entries with names
that do not match any formal (typos, renamed parameters, parameters from
a different function) are either silently dropped (default) or reported
via [`warning()`](https://rdrr.io/r/base/warning.html).

## Usage

``` r
filter_valid_args(
  params,
  target_fun,
  exclude = character(0),
  warn_unknown = FALSE
)
```

## Arguments

- params:

  Named list of user-supplied parameter values.

- target_fun:

  A function (not its name). Its
  [`formals`](https://rdrr.io/r/base/formals.html) determine the set of
  valid parameter names.

- exclude:

  Character vector of formal-argument names to explicitly *exclude* from
  the allowlist, even though they appear in `formals(target_fun)`.
  Useful for arguments the caller manages internally (e.g.,
  `df.analysis`, `confounders.name`, `details`, `plot.sg`, `quiet`) and
  does not want users to override via the params list.

- warn_unknown:

  Logical. If `TRUE`, emit a
  [`warning()`](https://rdrr.io/r/base/warning.html) listing any
  `params` names that were dropped because they do not appear in
  `formals(target_fun)` (typo-detection aid). Default `FALSE` (silent
  drop, matching legacy behavior).

## Value

A named list, the subset of `params` whose names are in
`formals(target_fun)` and not in `exclude`.

## Details

This replaces the older "hand-maintained `valid_pnames` whitelist"
pattern in `.run_fs_analysis_gen()` and `.run_grf_analysis_gen()`. New
[`forestsearch()`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
or `grf.subg.harm.*()` arguments flow through automatically; authors of
new parameters no longer need to update a whitelist in a separate file
to make the new parameter reachable from user-level `fs_params` /
`grf_params` lists.
