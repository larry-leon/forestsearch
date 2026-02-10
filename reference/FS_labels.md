# Convert Factor Code to Label

Converts q-indexed codes to human-readable labels using the confs_labels
mapping. Supports both full format ("q1.1", "q3.0") and short format
("q1", "q3"). Handles vector input via recursion.

## Usage

``` r
FS_labels(Qsg, confs_labels)
```

## Arguments

- Qsg:

  Character. Factor code in format "q." or "q". For the full format,
  action 0 = NOT (complement), action 1 = IN (member). Short format
  defaults to action 1. Can also be a character vector for vectorized
  input.

- confs_labels:

  Character vector. Labels for each factor, indexed by factor number.

## Value

Character. Human-readable label wrapped in braces, e.g., "{age \<= 50}"
or "!{age \<= 50}" for complement. Returns the original code if no match
is found.

## Examples

``` r
if (FALSE) { # \dontrun{
labels <- c("age <= 50", "tumor size <= 20", "nodes <= 3")
FS_labels("q1.1", labels)
FS_labels("q1.0", labels)
FS_labels("q2", labels)
FS_labels(c("q1", "q3"), labels)
} # }
```
