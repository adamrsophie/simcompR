# Building your first block diagram

`simcompR` builds block diagrams in R code and simulates them. A diagram
is a set of blocks (sources, gains, sums, integrators, delays, …) whose
output ports are wired to input ports of other blocks. This vignette
builds the classic negative-feedback loop, explains each step and checks
the result against the exact solution.

``` r

library(simcompR)
```

## The model

We want to simulate

\frac{dx}{dt} = r - x, \qquad r = 1, \qquad x(0) = 0,

whose exact solution is x(t) = 1 - e^{-t}. As a block diagram this is a
constant reference r, a sum block computing the error r - x, and an
integrator whose input is the error and whose output is x. The output of
the integrator is fed back into the sum.

## Step 1: create a model and add blocks

A model is a container.
[`sim_add_block()`](https://adamrsophie.github.io/simcompR/reference/sim_add_block.md)
adds a block of a registered type with a unique name and its parameters,
and returns a handle.

``` r

model <- sim_model("feedback_system")

reference <- sim_add_block(model, type = "constant", name = "reference", value = 1)
error <- sim_add_block(model, type = "sum", name = "error", signs = c("+", "-"))
plant <- sim_add_block(model, type = "integrator", name = "plant", initial = 0)
```

The sum block has one input port per sign, so `signs = c("+", "-")` (or
the string `"+-"`) gives two inputs and computes `in1 - in2`.

Models are modified in place: there is no need to reassign `model`.

## Step 2: connect ports

[`sim_connect()`](https://adamrsophie.github.io/simcompR/reference/sim_connect.md)
wires an output port of one block to an input port of another. Ports are
numbered from 1; the default is port 1.

``` r

sim_connect(model, from = reference, to = error, input = 1)
sim_connect(model, from = plant, to = error, input = 2)
sim_connect(model, from = error, to = plant)
model
#> <sim_model 'feedback_system'> 3 block(s), 3 connection(s), 0 event(s)
#>   reference [constant]
#>   error [sum]
#>   plant [integrator]
#> Connections:
#>   reference:1 -> error:1
#>   plant:1 -> error:2
#>   error:1 -> plant:1
```

Connections are checked immediately where possible: connecting to a port
that does not exist, or connecting two sources to the same input, is an
error.

``` r

sim_connect(model, from = reference, to = error, input = 3)
#> Error:
#> ! Block 'error' has 2 input port(s); port 3 does not exist.
sim_connect(model, from = reference, to = plant)
#> Error:
#> ! Input 1 of block 'plant' is already connected to 'error' output 1. Disconnect it first with sim_disconnect().
```

## Step 3: compile

Compilation validates the whole diagram and produces an immutable
execution plan.

``` r

plan <- sim_compile(model)
plan
#> <sim_plan 'feedback_system'> 3 block(s), 1 continuous state(s), 0 discrete state(s), 3 signal element(s)
#> Execution order: plant > reference > error
```

Notice the execution order: `plant` comes before `error`, even though it
was added later. An integrator’s output depends only on its state, not
on its input at the same instant (it has no *direct feedthrough*), so it
can be evaluated first. The sum does have direct feedthrough, so it must
wait for both of its inputs. The feedback loop is therefore valid: it
passes through a stateful block. See
[`vignette("compilation")`](https://adamrsophie.github.io/simcompR/articles/compilation.md)
for the details.

## Step 4: run

``` r

result <- sim_run(plan, tspan = c(0, 10), times = seq(0, 10, by = 0.01))
result
#> <sim_result 'feedback_system'> status: success
#>   1001 recorded time(s) in [0, 10]
#>   3 signal element(s), 1 continuous and 0 discrete state(s), 0 event(s)
```

`times` are the output times. They are independent of the internal steps
the adaptive solver takes.

## Step 5: inspect and check

``` r

x <- sim_signal(result, "plant")
t <- sim_times(result)
max(abs(x - (1 - exp(-t))))
#> [1] 1.046722e-08
```

The error is far below the default solver tolerance (`rtol = 1e-6`,
`atol = 1e-8`). Tighter tolerances give a smaller error:

``` r

tight <- sim_run(plan, tspan = c(0, 10), times = t, rtol = 1e-10, atol = 1e-12)
max(abs(sim_signal(tight, "plant") - (1 - exp(-t))))
#> [1] 8.902844e-13
```

Every block output is recorded, so the error signal is also available:

``` r

head(sim_signal(result, c("plant", "error")))
#>            plant     error
#> [1,] 0.000000000 1.0000000
#> [2,] 0.009950177 0.9900498
#> [3,] 0.019801323 0.9801987
#> [4,] 0.029554462 0.9704455
#> [5,] 0.039210556 0.9607894
#> [6,] 0.048770570 0.9512294
plot(result, c("plant", "error"))
```

![](simcompR_files/figure-html/unnamed-chunk-10-1.png)

[`summary()`](https://rdrr.io/r/base/summary.html) reports the solver
statistics and final states:

``` r

summary(result)
#> Simulation status: success
#> Time range: [0, 10], 1001 row(s)
#> Solver: lsoda, rtol = 1e-06, atol = 1e-08
#> Segments: 1, steps: 1004, RHS evaluations: 1467, roots located: 0
#> Final continuous states:
#>     plant 
#> 0.9999546
```

## Where to go next

- [`vignette("ports-and-signals")`](https://adamrsophie.github.io/simcompR/articles/ports-and-signals.md):
  port numbering, names, vector signals and the sizing rules.
- [`vignette("continuous-simulation")`](https://adamrsophie.github.io/simcompR/articles/continuous-simulation.md):
  solvers, tolerances and output times.
- [`vignette("discrete-time")`](https://adamrsophie.github.io/simcompR/articles/discrete-time.md)
  and
  [`vignette("multirate")`](https://adamrsophie.github.io/simcompR/articles/multirate.md):
  sampled blocks.
- [`vignette("hybrid-events")`](https://adamrsophie.github.io/simcompR/articles/hybrid-events.md):
  events, zero crossings and resets.
- [`vignette("subsystems")`](https://adamrsophie.github.io/simcompR/articles/subsystems.md),
  [`vignette("algebraic-loops")`](https://adamrsophie.github.io/simcompR/articles/algebraic-loops.md),
  [`vignette("state-machines")`](https://adamrsophie.github.io/simcompR/articles/state-machines.md),
  [`vignette("experiments")`](https://adamrsophie.github.io/simcompR/articles/experiments.md).
