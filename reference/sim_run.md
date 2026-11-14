# Simulate a compiled plan

Runs a
[`sim_compile()`](https://adamrsophie.github.io/simcompR/reference/sim_compile.md)
plan over `tspan` and returns a `sim_result`.

## Usage

``` r
sim_run(
  plan,
  tspan,
  times = NULL,
  method = "lsoda",
  rtol = 1e-06,
  atol = 1e-08,
  hmax = NULL,
  maxsteps = 5000,
  initial = NULL,
  include_ticks = TRUE,
  max_events = 1000,
  zc_tol = 1e-09,
  on_failure = c("error", "warn"),
  ...
)
```

## Arguments

- plan:

  A compiled `sim_plan`.

- tspan:

  Numeric `c(start, end)`.

- times:

  Output times within `tspan`. Defaults to 201 equally spaced points.
  Solver output times are independent of the internal steps taken by the
  integrator.

- method:

  Integration method passed to
  [`deSolve::ode()`](https://rdrr.io/pkg/deSolve/man/ode.html). Models
  with zero crossings need one of `"lsoda"`, `"lsode"`, `"lsodes"` or
  `"lsodar"`.

- rtol, atol:

  Relative and absolute tolerances.

- hmax:

  Optional maximum step size. When the model has zero-crossing functions
  and `hmax` is not given, it is `(tf - t0) / 100`, because the root
  finder cannot see two crossings inside one step; give a smaller value
  for signals that cross faster.

- maxsteps:

  Maximum number of internal steps per segment.

- initial:

  Optional named list overriding initial states, keyed by block name.

- include_ticks:

  Also record a row at every sample hit, so held signals are captured at
  their changes.

- max_events:

  Maximum number of zero-crossing events before the run is stopped with
  status `"event_limit"`. Guards against chattering and Zeno behavior.

- zc_tol:

  A zero-crossing function whose value at the start of a segment is
  within `zc_tol` of zero is treated as exactly zero there, which
  prevents retriggering the same root.

- on_failure:

  `"error"` raises an error when the simulation fails; `"warn"` returns
  the partial result with a warning.

- ...:

  Further arguments for
  [`deSolve::ode()`](https://rdrr.io/pkg/deSolve/man/ode.html).

## Value

A `sim_result`. See
[`sim_signal()`](https://adamrsophie.github.io/simcompR/reference/sim_signal.md),
[`sim_states()`](https://adamrsophie.github.io/simcompR/reference/sim_signal.md),
[`sim_events()`](https://adamrsophie.github.io/simcompR/reference/sim_signal.md),
[`sim_diagnostics()`](https://adamrsophie.github.io/simcompR/reference/sim_diagnostics.md)
and [`plot()`](https://rdrr.io/r/graphics/plot.default.html).

## Details

**Timeline.** The interval is cut into segments at every sample hit of a
discrete block, every scheduled event time and every breakpoint of a
source (for example the time of a step). Sample hits are computed as
`offset + k * period` with integer `k`, never by accumulating `period`,
so they do not drift. Hits closer than `1e-12 * max(1, |t|)` are treated
as simultaneous.

**At a segment boundary** the simulator, in this order:

1.  applies scheduled events (in priority order);

2.  evaluates the outputs of all blocks with a sample hit, in execution
    order, and holds them;

3.  updates the discrete states of those blocks, using the inputs from
    step 2. Outputs are computed before states are updated, so a unit
    delay produces `y[k] = u[k-1]`;

4.  refreshes the modes of switching blocks from their current inputs.

**Inside a segment** discrete states and held outputs are constant and
the continuous states are integrated with
[`deSolve::ode()`](https://rdrr.io/pkg/deSolve/man/ode.html). The
derivative function has no side effects. When zero-crossing functions
exist (from
[`sim_add_event()`](https://adamrsophie.github.io/simcompR/reference/sim_add_event.md)
or from `switch`, `relay` and `stop` blocks) the root finder stops the
solver at the crossing; the event is applied and integration restarts
there. A crossing that ends exactly at zero is not detected again at the
restart.

**Recording.** Signals are computed after the solve from the accepted
states, never inside the derivative function. A value recorded at a
sample hit is the value after the hit (right-continuous). At an event
the result has two rows with the same time: the state before and after
the reset.

**Failure.** When the solver fails, a block signals an error or the
event limit is reached, the result is marked with a status other than
`"success"` and contains only the times actually reached; no values are
extrapolated. With `on_failure = "error"` (the default for solver and
block failures) an error is raised instead.

## Examples

``` r
m <- sim_model("feedback")
sim_add_block(m, "constant", "r", value = 1)
sim_add_block(m, "subtract", "e")
sim_add_block(m, "integrator", "x")
sim_connect(m, "r", "e", input = 1)
sim_connect(m, "x", "e", input = 2)
sim_connect(m, "e", "x")
r <- sim_run(sim_compile(m), tspan = c(0, 5), times = 0:5)
cbind(numeric = sim_signal(r, "x"), exact = 1 - exp(-(0:5)))
#>        numeric     exact
#> [1,] 0.0000000 0.0000000
#> [2,] 0.6321204 0.6321206
#> [3,] 0.8646649 0.8646647
#> [4,] 0.9502131 0.9502129
#> [5,] 0.9816845 0.9816844
#> [6,] 0.9932621 0.9932621
```
