# Recover the GRF arguments used inside a forestsearch run

Returns the exact list of arguments that
[`forestsearch()`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
passed to its inner GRF subgroup-identification call
([`grf.subg.harm.survival()`](https://larry-leon.github.io/forestsearch/reference/grf.subg.harm.survival.md)
for survival outcomes,
[`grf.subg.harm.glm()`](https://larry-leon.github.io/forestsearch/reference/grf.subg.harm.glm.md)
for binary, continuous, and count outcomes). The result is suitable for
`do.call(grf.subg.harm.glm, args)` or
`do.call(grf.subg.harm.survival, args)` to reproduce the standalone GRF
run with identical settings.

## Usage

``` r
extract_grf_args(fs)
```

## Arguments

- fs:

  A fitted `forestsearch` object.

## Value

A named list of arguments. For survival outcomes the list targets
[`grf.subg.harm.survival`](https://larry-leon.github.io/forestsearch/reference/grf.subg.harm.survival.md);
for binary, continuous, or count outcomes it targets
[`grf.subg.harm.glm`](https://larry-leon.github.io/forestsearch/reference/grf.subg.harm.glm.md).
The list carries an attribute `"grf_function"` holding the target
function name as a character string.

## Details

Typical uses:

- Sensitivity analyses: modify one or two arguments and re-fit GRF
  without re-running the full forestsearch pipeline, e.g.
  `do.call(grf.subg.harm.glm, modifyList(extract_grf_args(fs), list(maxdepth = 3L)))`.

- Reproducibility audits: confirm what the inner GRF saw.

- Independent variable-importance or policy-tree extraction from a
  forest fit with the same configuration.

This function reads from `fs$args_call_all` (populated by
[`forestsearch()`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
via [`mget()`](https://rdrr.io/r/base/get.html) of its formals) and from
`fs$outcome_type`. It is the user-facing complement to the internal
builders that the
[`forestsearch()`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
call sites use.

## See also

[`forestsearch`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md),
[`grf.subg.harm.glm`](https://larry-leon.github.io/forestsearch/reference/grf.subg.harm.glm.md),
[`grf.subg.harm.survival`](https://larry-leon.github.io/forestsearch/reference/grf.subg.harm.survival.md).

## Examples

``` r
if (FALSE) { # \dontrun{
# Reproduce the inner GRF run
args <- extract_grf_args(fs)
fn   <- match.fun(attr(args, "grf_function"))
grf_standalone <- do.call(fn, args)

# Sensitivity: deeper trees, everything else unchanged
args_deep <- modifyList(args, list(maxdepth = 3L))
grf_deep  <- do.call(fn, args_deep)
} # }
```
