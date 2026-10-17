# Compute node metrics for a policy tree

Aggregates scores by leaf node and calculates treatment effect
differences

## Usage

``` r
compute_node_metrics(data, dr.scores, tree, X, n.min)
```

## Arguments

- data:

  Data frame. Original data

- dr.scores:

  Matrix. Doubly robust scores

- tree:

  Policy tree object

- X:

  Matrix. Covariate matrix

- n.min:

  Integer. Minimum subgroup size

## Value

Data frame with node metrics

## Examples

``` r
if (FALSE) { # \dontrun{
# compute_node_metrics() is called internally by grf.subg.harm.survival().
# See grf.subg.harm.survival() for the standard entry point.
} # }
```
