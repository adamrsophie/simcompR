# Change parameters of a compiled plan

Returns a new plan with updated model parameters. Blocks are rebuilt
with the new values. When block sizes, feedthrough, state dimensions,
sample times and connections are unchanged the dependency analysis is
reused (`plan$structure_reused` is `TRUE`); otherwise the graph is
analyzed again, so a change that alters the structure (for example
making the `D` matrix of a state-space block nonzero) is handled
correctly.

## Usage

``` r
sim_set_params(plan, ...)

sim_param_values(plan)
```

## Arguments

- plan:

  A `sim_plan`.

- ...:

  Named parameter values, or a single named list.

## Value

A new `sim_plan`.

## Examples

``` r
m <- sim_model("decay", params = list(a = 1))
sim_add_block(m, "integrator", "x", initial = 1)
sim_add_block(m, "gain", "k", gain = sim_param("a"))
sim_connect(m, "x", "k")
sim_add_block(m, "gain", "neg", gain = -1)
sim_connect(m, "k", "neg")
sim_connect(m, "neg", "x")
p1 <- sim_compile(m)
p2 <- sim_set_params(p1, a = 3)
p2$structure_reused
#> [1] TRUE
sim_states(sim_run(p2, c(0, 1), times = 1))
#>               x
#> [1,] 0.04978721
```
