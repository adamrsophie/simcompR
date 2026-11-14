# Change block parameters in a model

Change block parameters in a model

## Usage

``` r
sim_set_block_params(model, block, ...)
```

## Arguments

- model:

  A
  [`sim_model()`](https://adamrsophie.github.io/simcompR/reference/sim_model.md).

- block:

  Block name or handle.

- ...:

  New parameter values.

## Value

The model, invisibly.

## Examples

``` r
m <- sim_model()
sim_add_block(m, "gain", "k", gain = 3)
sim_set_block_params(m, "k", gain = 5)
```
