# Extract simulation output

`sim_signal()` returns recorded signals by block name. A block with one
output port is named by the block name; other ports are `block.port`.
Subsystem instances are named by the instance name, inner blocks by
`instance/block`, and a `scope` block by its own name (its input
signal). `sim_states()` and `sim_discrete_states()` return the state
trajectories; `sim_events()` the event log; `sim_times()` the recorded
times.

## Usage

``` r
sim_signal(result, names = NULL, drop = TRUE)

sim_states(result)

sim_discrete_states(result)

sim_events(result)

sim_times(result)
```

## Arguments

- result:

  A `sim_result`.

- names:

  Signal names. `NULL` returns all signals.

- drop:

  Drop dimensions of extent one, so a single scalar signal is returned
  as a vector.

## Value

A numeric vector or matrix (signals and states) or a data frame
(events).

## Details

Rows follow `result$time`. At an event the time appears twice: the row
before and the row after the reset.

## Examples

``` r
m <- sim_model()
sim_add_block(m, "sine", "u")
sim_add_block(m, "integrator", "x", initial = -1)
sim_connect(m, "u", "x")
r <- sim_run(sim_compile(m), c(0, pi), times = c(0, pi / 2, pi))
sim_times(r)
#> [1] 0.000000 1.570796 3.141593
sim_signal(r, "x")      # -cos(t)
#> [1] -1.000000e+00  1.705167e-06  9.999982e-01
sim_states(r)
#>                  x
#> [1,] -1.000000e+00
#> [2,]  1.705167e-06
#> [3,]  9.999982e-01
```
