# Add ID Column to Data Frame

Ensures that a data frame has a unique ID column. If `id.name` is not
provided, a column named `"id"` is added. If `id.name` is provided but
does not exist in the data frame, it is created with unique integer
values.

## Usage

``` r
add_id_column(df.analysis, id.name = NULL)
```

## Arguments

- df.analysis:

  Data frame to which the ID column will be added.

- id.name:

  Character. Name of the ID column to add (default is `NULL`, which uses
  `"id"`).

## Value

Data frame with the ID column added if necessary.
