# Hierarchical subsystems

A subsystem is an ordinary
[`sim_model()`](https://adamrsophie.github.io/simcompR/reference/sim_model.md)
whose external interface is declared with `"inport"` and `"outport"`
blocks. Port numbers come from their `index` parameter. The subsystem is
placed into a parent model with
`sim_add_block(parent, "subsystem", name, system = sub, params = ...)`.

## Usage

``` r
sim_subsystem(name, params = list())
```

## Arguments

- name:

  Subsystem name.

- params:

  Default parameters of the subsystem. Instance `params` override them.

## Value

A `sim_model`.

## Details

At compile time every subsystem instance is flattened into the parent:
inner block `b` of instance `s` becomes block `s/b`. Each instance gets
its own copy of the inner blocks, so instances never share state. A
[`sim_param()`](https://adamrsophie.github.io/simcompR/reference/sim_param.md)
reference inside the subsystem is resolved first in the instance's
`params`, then in the parent's scope. The flattened model is
mathematically identical to building the same blocks by hand.

`sim_subsystem()` is a convenience wrapper around
[`sim_model()`](https://adamrsophie.github.io/simcompR/reference/sim_model.md).

## Examples

``` r
# A first-order lag y' = (u - y) / tau, with inv_tau = 1 / tau.
lag <- sim_subsystem("lag", params = list(inv_tau = 1))
sim_add_block(lag, "inport", "u", index = 1)
sim_add_block(lag, "subtract", "e")
sim_add_block(lag, "gain", "k", gain = sim_param("inv_tau"))
sim_add_block(lag, "integrator", "y", initial = 0)
sim_add_block(lag, "outport", "out", index = 1)
sim_connect(lag, "u", "e", input = 1)
sim_connect(lag, "y", "e", input = 2)
sim_connect(lag, "e", "k")
sim_connect(lag, "k", "y")
sim_connect(lag, "y", "out")

m <- sim_model("two_lags")
sim_add_block(m, "step", "u", time = 0)
sim_add_block(m, "subsystem", "fast", system = lag, params = list(inv_tau = 10))
sim_add_block(m, "subsystem", "slow", system = lag, params = list(inv_tau = 1))
sim_connect(m, "u", "fast")
sim_connect(m, "u", "slow")
r <- sim_run(sim_compile(m), tspan = c(0, 3))
tail(sim_signal(r, c("fast", "slow")), 1)
#>        fast      slow
#> [201,]    1 0.9502129
```
