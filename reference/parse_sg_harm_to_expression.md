# Parse sg.harm Factor Names to Expression

Converts ForestSearch factor names (e.g., "er.0", "grade3.1") into
human-readable R expressions (e.g., "er \<= 0", "grade3 == 1").

## Usage

``` r
parse_sg_harm_to_expression(sg_harm, fs.est = NULL)
```

## Arguments

- sg_harm:

  Character vector of factor names from fs.est\$sg.harm.

- fs.est:

  ForestSearch object (for accessing confs_labels if available).

## Value

Character string expression or NULL if parsing fails.
