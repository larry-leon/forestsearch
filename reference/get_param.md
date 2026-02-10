# Get Parameter with Default Fallback

Safely retrieves a named element from a list, returning a default value
if the element is missing or `NULL`.

## Usage

``` r
get_param(args_list, param_name, default_value)
```

## Arguments

- args_list:

  List to extract from.

- param_name:

  Character. Name of the element to retrieve.

- default_value:

  Default value to return if element is missing or `NULL`.

## Value

The value of `args_list[[param_name]]` if present and non-`NULL`,
otherwise `default_value`.
