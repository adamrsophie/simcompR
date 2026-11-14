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

sim_param_expr(expr)
```

## Arguments

- name:

  Parameter name.

- expr:

  A portable expression (string) in parameter names.

## Value

An object of class `sim_param`.

## Details

References can be combined with numbers and with each other by
arithmetic (`-sim_param("a")`, `1 / sim_param("tau")`,
`2 * sim_param("a") + sim_param("b")`) and by the functions `abs`,
`sqrt`, `exp`, `log`, `sin`, `cos`, `tan`, `floor`, `ceiling` and
`sign`. The result is a parameter expression, kept as a portable string
(so it can be saved with
[`sim_write_json()`](https://adamrsophie.github.io/simcompR/reference/sim_write_json.md))
and evaluated at compile time. `sim_param_expr()` builds one directly
from a string. References may also appear inside list parameters. A
vector cannot hold references (R would drop their class); use an
expression such as `sim_param_expr("c(a, 2 * a)")` instead.

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

# The same model with one gain: k = -a
m2 <- sim_model("decay", params = list(a = 2))
sim_add_block(m2, "integrator", "x", initial = 1)
sim_add_block(m2, "gain", "k", gain = -sim_param("a"))
sim_connect(m2, "x", "k")
sim_connect(m2, "k", "x")
sim_param_expr("-a / 2")
#> <sim_param expression -a / 2>
```
