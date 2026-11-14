# Compile a model into an execution plan

Compilation turns a mutable
[`sim_model()`](https://adamrsophie.github.io/simcompR/reference/sim_model.md)
into an immutable `sim_plan` in these stages:

## Usage

``` r
sim_compile(
  model,
  algebraic_loops = c("error", "solve"),
  loop_tol = 1e-10,
  loop_maxit = 50
)
```

## Arguments

- model:

  A
  [`sim_model()`](https://adamrsophie.github.io/simcompR/reference/sim_model.md).

- algebraic_loops:

  `"error"` (default) or `"solve"`.

- loop_tol:

  Absolute residual tolerance for loop solving.

- loop_maxit:

  Maximum Newton iterations per loop solve.

## Value

A `sim_plan` object. Use
[`sim_run()`](https://adamrsophie.github.io/simcompR/reference/sim_run.md)
to simulate it,
[`sim_set_params()`](https://adamrsophie.github.io/simcompR/reference/sim_set_params.md)
to change parameters and
[`sim_diagnostics()`](https://adamrsophie.github.io/simcompR/reference/sim_diagnostics.md)
to inspect the analysis.

## Details

1.  Flatten subsystems and resolve
    [`sim_param()`](https://adamrsophie.github.io/simcompR/reference/sim_param.md)
    references.

2.  Build every block and validate its parameters.

3.  Validate connections: blocks and ports exist and every input has
    exactly one source.

4.  Propagate signal sizes along connections and check compatibility.

5.  Build the direct-feedthrough dependency graph. There is an edge from
    block A to block B when an output of A feeds an input of B that has
    direct feedthrough. Inputs of integrators, unit delays and other
    state-output blocks create no edge, so feedback through them is
    valid.

6.  Find strongly connected components of that graph. A component with
    more than one block, or a block feeding itself, is an algebraic
    loop.

7.  Order the blocks topologically. Ties are broken by block name, so
    the order, the state vector layout and the results do not depend on
    the order in which blocks or connections were added.

8.  Lay out the signal, continuous state and discrete state vectors and
    collect sample times, zero-crossing functions and events.

Algebraic loops are rejected by default. With
`algebraic_loops = "solve"` each loop is solved at every evaluation by
Newton iterations on the loop outputs, with a finite-difference Jacobian
and backtracking; the run fails with an error when the residual does not
fall below `loop_tol`.

## Examples

``` r
m <- sim_model("feedback")
sim_add_block(m, "constant", "r", value = 1)
sim_add_block(m, "subtract", "e")
sim_add_block(m, "integrator", "x")
sim_connect(m, "r", "e", input = 1)
sim_connect(m, "x", "e", input = 2)
sim_connect(m, "e", "x")
plan <- sim_compile(m)
plan
#> <sim_plan 'feedback'> 3 block(s), 1 continuous state(s), 0 discrete state(s), 3 signal element(s)
#> Execution order: r > x > e
sim_diagnostics(plan)$order
#> [1] "r" "x" "e"

# A loop of direct-feedthrough blocks is an algebraic loop.
m2 <- sim_model("loop")
sim_add_block(m2, "constant", "u", value = 1)
sim_add_block(m2, "sum", "s", signs = "+-")
sim_add_block(m2, "gain", "k", gain = 2)
sim_connect(m2, "u", "s", input = 1)
sim_connect(m2, "k", "s", input = 2)
sim_connect(m2, "s", "k")
try(sim_compile(m2))
#> Error : Algebraic loop(s) detected: every block on the cycle has direct feedthrough, so the outputs depend on themselves at the same time.
#>   k -> s -> k
#> Break the loop with a stateful block (integrator, unit_delay, ...) or compile with algebraic_loops = "solve".
plan2 <- sim_compile(m2, algebraic_loops = "solve")
sim_signal(sim_run(plan2, c(0, 1), times = 1), "s") # y = 1 - 2y, so 1/3
#>         s 
#> 0.3333333 
```
