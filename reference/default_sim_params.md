# Default ForestSearch parameters (general)

Returns the default `fs_params` list used by
[`run_simulation_analysis()`](https://larry-leon.github.io/forestsearch/reference/run_simulation_analysis.md).
The `parallel_args` entry defaults to `list(plan = "sequential")`
because
[`run_simulation_analysis()`](https://larry-leon.github.io/forestsearch/reference/run_simulation_analysis.md)
is designed to be called inside a `foreach() %dofuture%` loop (one
replicate per worker). Running the inner
[`forestsearch()`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
multisession in that context produces nested parallelism: each outer
worker tries to spawn its own pool of inner workers, which `parallelly`
rejects with a 300\\ limit error. Users who call
[`run_simulation_analysis()`](https://larry-leon.github.io/forestsearch/reference/run_simulation_analysis.md)
once at the top level (i.e., not inside `%dofuture%`) can opt back into
multisession by passing
`fs_params = list(parallel_args = list(plan = "multisession", workers = N))`.

## Usage

``` r
default_sim_params()
```
