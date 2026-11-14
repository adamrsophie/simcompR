# Create a block-diagram model

A `sim_model` is a mutable container of blocks, connections and events.
Functions such as
[`sim_add_block()`](https://adamrsophie.github.io/simcompR/reference/sim_add_block.md)
and
[`sim_connect()`](https://adamrsophie.github.io/simcompR/reference/sim_connect.md)
modify it in place. A model is never simulated directly:
[`sim_compile()`](https://adamrsophie.github.io/simcompR/reference/sim_compile.md)
takes a snapshot and turns it into an immutable [execution
plan](https://adamrsophie.github.io/simcompR/reference/sim_compile.md),
so later edits to the model never change an existing plan.

## Usage

``` r
sim_model(name = "model", params = list())
```

## Arguments

- name:

  Model name.

- params:

  Named list of model parameters, referenced from blocks with
  [`sim_param()`](https://adamrsophie.github.io/simcompR/reference/sim_param.md).

## Value

An object of class `sim_model`.

## See also

[`sim_add_block()`](https://adamrsophie.github.io/simcompR/reference/sim_add_block.md),
[`sim_connect()`](https://adamrsophie.github.io/simcompR/reference/sim_connect.md),
[`sim_compile()`](https://adamrsophie.github.io/simcompR/reference/sim_compile.md).

## Examples

``` r
m <- sim_model("feedback_system")
ref <- sim_add_block(m, "constant", "reference", value = 1)
err <- sim_add_block(m, "sum", "error", signs = "+-")
plant <- sim_add_block(m, "integrator", "plant", initial = 0)
sim_connect(m, ref, err, input = 1)
sim_connect(m, plant, err, input = 2)
sim_connect(m, err, plant)
m
#> <sim_model 'feedback_system'> 3 block(s), 3 connection(s), 0 event(s)
#>   reference [constant]
#>   error [sum]
#>   plant [integrator]
#> Connections:
#>   reference:1 -> error:1
#>   plant:1 -> error:2
#>   error:1 -> plant:1
```
