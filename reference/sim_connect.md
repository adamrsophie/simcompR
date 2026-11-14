# Connect two blocks

Connects output port `output` of block `from` to input port `input` of
block `to`. An output may feed several inputs; an input accepts exactly
one connection.

## Usage

``` r
sim_connect(model, from, to, output = 1, input = 1)
```

## Arguments

- model:

  A
  [`sim_model()`](https://adamrsophie.github.io/simcompR/reference/sim_model.md).

- from, to:

  Blocks, given as names or `sim_block` handles.

- output, input:

  Port indices or port names.

## Value

The model, invisibly.

## Examples

``` r
m <- sim_model()
sim_add_block(m, "constant", "a", value = 2)
sim_add_block(m, "gain", "k", gain = 3)
sim_connect(m, "a", "k")
m
#> <sim_model 'model'> 2 block(s), 1 connection(s), 0 event(s)
#>   a [constant]
#>   k [gain]
#> Connections:
#>   a:1 -> k:1
```
