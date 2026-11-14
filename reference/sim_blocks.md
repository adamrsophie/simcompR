# Inspect the blocks and connections of a model

Inspect the blocks and connections of a model

## Usage

``` r
sim_blocks(model)

sim_connections(model)
```

## Arguments

- model:

  A
  [`sim_model()`](https://adamrsophie.github.io/simcompR/reference/sim_model.md).

## Value

`sim_blocks()` returns a data frame with one row per block;
`sim_connections()` returns the connection table.

## Examples

``` r
m <- sim_model()
sim_add_block(m, "constant", "a", value = 2)
sim_add_block(m, "gain", "k", gain = 3)
sim_connect(m, "a", "k")
sim_blocks(m)
#>   name     type category
#> 1    a constant   source
#> 2    k     gain     math
sim_connections(m)
#>   from from_port to to_port
#> 1    a         1  k       1
```
