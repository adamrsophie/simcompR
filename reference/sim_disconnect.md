# Remove a connection or a block

Remove a connection or a block

## Usage

``` r
sim_disconnect(model, to, input = 1)

sim_remove_block(model, block)
```

## Arguments

- model:

  A
  [`sim_model()`](https://adamrsophie.github.io/simcompR/reference/sim_model.md).

- to, input:

  Destination block and input port of the connection to remove.

- block:

  Block to remove, together with all its connections.

## Value

The model, invisibly.

## Examples

``` r
m <- sim_model()
sim_add_block(m, "constant", "a", value = 2)
sim_add_block(m, "gain", "k", gain = 3)
sim_connect(m, "a", "k")
sim_disconnect(m, "k", 1)
sim_remove_block(m, "a")
m
#> <sim_model 'model'> 1 block(s), 0 connection(s), 0 event(s)
#>   k [gain]
```
