# Add a block to a model

Add a block to a model

## Usage

``` r
sim_add_block(model, type, name, ...)
```

## Arguments

- model:

  A
  [`sim_model()`](https://adamrsophie.github.io/simcompR/reference/sim_model.md).

- type:

  Block type, see
  [`sim_block_types()`](https://adamrsophie.github.io/simcompR/reference/sim_block_types.md).

- name:

  Unique block name. Names may not contain `/`, which separates
  subsystem levels in compiled names.

- ...:

  Block parameters. Any parameter may be a
  [`sim_param()`](https://adamrsophie.github.io/simcompR/reference/sim_param.md)
  reference.

## Value

A `sim_block` handle that can be passed to
[`sim_connect()`](https://adamrsophie.github.io/simcompR/reference/sim_connect.md).
The model is modified in place.

## Examples

``` r
m <- sim_model()
b <- sim_add_block(m, "sine", "u", amplitude = 2, frequency = 3)
b
#> <sim_block 'u' of type sine in model 'model'>
```
