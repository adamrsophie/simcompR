# Refer to a model or subsystem parameter

A block parameter given as `sim_param("name")` is looked up when the
model is compiled: first in the parameters of the enclosing subsystem
instance, then in the parameters of the parent, up to the top-level
model. Changing a referenced parameter with
[`sim_set_params()`](https://adamrsophie.github.io/simcompR/reference/sim_set_params.md)
or
[`sim_sweep()`](https://adamrsophie.github.io/simcompR/reference/sim_sweep.md)
does not require editing the block.

## Usage

``` r
sim_param(name)
```

## Arguments

- name:

  Parameter name.

## Value

An object of class `sim_param`.

## Examples

``` r
m <- sim_model("decay", params = list(a = 2))
sim_add_block(m, "integrator", "x", initial = 1)
sim_add_block(m, "gain", "k", gain = sim_param("a"))
sim_connect(m, "x", "k")
sim_add_block(m, "gain", "neg", gain = -1)
sim_connect(m, "k", "neg")
sim_connect(m, "neg", "x")
plan <- sim_compile(m)
sim_param_values(plan)
#> $a
#> [1] 2
#> 
```
