# Algebraic loops

``` r

library(simcompR)
```

A feedback cycle made only of direct-feedthrough blocks defines a signal
implicitly, y = f(y, u), at every instant. Such a cycle is an *algebraic
loop*. A cycle that passes through an integrator, a delay or any other
block without direct feedthrough is ordinary dynamic feedback and is not
a loop.

## Detection

Compilation builds the graph of direct-feedthrough dependencies and
finds its strongly connected components. A component with more than one
block, or a block wired to itself, is a loop. By default it is an error,
with one cycle printed:

``` r

m <- sim_model("loop")
sim_add_block(m, "constant", "u", value = 1)
sim_add_block(m, "sum", "s", signs = "+-")
sim_add_block(m, "gain", "k", gain = 2)
sim_connect(m, "u", "s", input = 1)
sim_connect(m, "k", "s", input = 2)
sim_connect(m, "s", "k")
sim_compile(m)
#> Error:
#> ! Algebraic loop(s) detected: every block on the cycle has direct feedthrough, so the outputs depend on themselves at the same time.
#>   k -> s -> k
#> Break the loop with a stateful block (integrator, unit_delay, ...) or compile with algebraic_loops = "solve".
```

`simcompR` never breaks a loop by inserting a hidden delay: that would
change the mathematics silently.

## Solving loops

With `algebraic_loops = "solve"` each loop is solved at every
evaluation. The unknowns are the outputs z of the blocks in the loop,
and the residual is F(z) - z, where F evaluates the loop blocks once
from the guessed outputs. Newton’s method with a finite-difference
Jacobian and a backtracking line search drives the residual below
`loop_tol`.

``` r

plan <- sim_compile(m, algebraic_loops = "solve")
plan
#> <sim_plan 'loop'> 3 block(s), 0 continuous state(s), 0 discrete state(s), 3 signal element(s)
#> Execution order: u > k > s
#> Algebraic loops (solved by Newton iteration):
#>   {k, s}
r <- sim_run(plan, c(0, 1), times = c(0, 1))
sim_signal(r, "s")   # s = 1 - 2 s, so 1/3
#> [1] 0.3333333 0.3333333
sim_diagnostics(r)$loop_residuals
#>      time     residual iterations
#> [1,]    0 5.551115e-17          2
#> [2,]    1 5.551115e-17          2
```

The residual and iteration count at every recorded time are part of the
result diagnostics. A loop that does not converge is an error, never a
silently inconsistent result.

### A nonlinear loop

y = \cos(y) has the unique solution y \approx 0.739085:

``` r

nl <- sim_model("cos_loop")
sim_add_block(nl, "math", "c", fun = "cos")
sim_add_block(nl, "gain", "y", gain = 1)
sim_connect(nl, "c", "y")
sim_connect(nl, "y", "c")
r <- sim_run(sim_compile(nl, algebraic_loops = "solve"), c(0, 1), times = 1)
sim_signal(r, "y")
#>         y 
#> 0.7390851
```

### Loops without a solution

y = y + 1 has no solution. The Jacobian of the residual is singular and
the run fails with an explanation:

``` r

bad <- sim_model("inconsistent")
sim_add_block(bad, "constant", "one")
sim_add_block(bad, "sum", "s")
sim_connect(bad, "one", "s", input = 1)
sim_connect(bad, "s", "s", input = 2)
sim_run(sim_compile(bad, algebraic_loops = "solve"), c(0, 1))
#> Error:
#> ! Simulation failed at t = 0: Algebraic loop {s} at t = 0: singular Jacobian, the loop has no unique solution..
```

## Loops inside dynamics

Loops can sit inside a continuous model. Here \dot x = s with s = -x -
s, so \dot x = -x/2:

``` r

dyn <- sim_model("dynamic_loop")
sim_add_block(dyn, "integrator", "x", initial = 1)
sim_add_block(dyn, "sum", "s", signs = "--")
sim_add_block(dyn, "gain", "k", gain = 1)
sim_connect(dyn, "x", "s", input = 1)
sim_connect(dyn, "k", "s", input = 2)
sim_connect(dyn, "s", "k")
sim_connect(dyn, "s", "x")
r <- sim_run(sim_compile(dyn, algebraic_loops = "solve"), c(0, 4),
             rtol = 1e-10, atol = 1e-12)
max(abs(sim_signal(r, "x") - exp(-sim_times(r) / 2)))
#> [1] 4.410894e-11
```

Every derivative evaluation solves the loop, which costs extra function
evaluations. When possible, solve the loop by hand (here, replace the
loop by a gain of 1/2) or break it with a physically meaningful state.

## Feedthrough that depends on parameters

Whether a state-space block has direct feedthrough depends on its `D`
matrix. A plan whose `D` changes from zero to nonzero through
[`sim_set_params()`](https://adamrsophie.github.io/simcompR/reference/sim_set_params.md)
is analyzed again, so a new loop is never missed. See
[`vignette("compilation")`](https://adamrsophie.github.io/simcompR/articles/compilation.md).

## Limits

- Loop solving uses a dense finite-difference Jacobian; it is meant for
  small loops.
- Switching blocks inside a loop keep the mode they had at the start of
  the segment.
- Consistent initialization of differential-algebraic systems is not
  covered here; see the roadmap.
