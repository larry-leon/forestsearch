# Summarize Bootstrap Event Counts

Provides summary statistics for event counts across bootstrap
iterations, helping assess the reliability of HR estimates when events
are sparse.

## Usage

``` r
summarize_bootstrap_events(boot_results, threshold = 5)
```

## Arguments

- boot_results:

  List. Output from forestsearch_bootstrap_dofuture()

- threshold:

  Integer. Minimum event threshold for flagging low counts (default: 5)

## Value

Invisibly returns a list with summary statistics:

- threshold:

  The event threshold used

- nb_boots:

  Total number of bootstrap iterations

- n_successful:

  Number of iterations that found a new subgroup

- original_H:

  List with low event counts for original H on bootstrap samples

- original_Hc:

  List with low event counts for original Hc on bootstrap samples

- new_Hstar:

  List with low event counts for new H\* on original data

- new_Hcstar:

  List with low event counts for new Hc\* on original data

## Details

This function summarizes event counts in four scenarios:

1.  ORIGINAL subgroup H evaluated on BOOTSTRAP samples

2.  ORIGINAL subgroup Hc evaluated on BOOTSTRAP samples

3.  NEW subgroup H\* (found in bootstrap) evaluated on ORIGINAL data

4.  NEW subgroup Hc\* (found in bootstrap) evaluated on ORIGINAL data

Low event counts (below threshold) can lead to unstable HR estimates.
This summary helps identify potential issues with sparse events.
