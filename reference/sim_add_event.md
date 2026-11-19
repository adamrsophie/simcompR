# Model events: zero crossings, scheduled times and state resets

`sim_add_event()` attaches an event to a model. When the trigger fires
the event can reset the state of any block and can stop the simulation.

## Usage

``` r
sim_add_event(
  model,
  name,
  trigger,
  reset = list(),
  terminal = FALSE,
  priority = 0
)

sim_remove_event(model, name)

sim_zero_crossing(
  signal,
  direction = "either",
  level = 0,
  port = 1,
  element = 1
)

sim_at(times)
```

## Arguments

- model:

  A
  [`sim_model()`](https://adamrsophie.github.io/simcompR/reference/sim_model.md).

- name:

  Unique event name; for `sim_remove_event()`, the event to remove.

- trigger:

  A trigger from `sim_zero_crossing()` or `sim_at()`.

- reset:

  Named list of state resets (see Details).

- terminal:

  If `TRUE`, the simulation stops after the instant at which the event
  fires.

- priority:

  Ordering among simultaneous events (smaller first).

- signal:

  Block (name or handle) whose output is watched.

- direction:

  `"rising"`, `"falling"` or `"either"`.

- level:

  Crossing level.

- port:

  Output port of `signal`.

- element:

  Element of a vector signal.

- times:

  Event times.

## Value

`sim_add_event()` and `sim_remove_event()` return the model invisibly;
the trigger constructors return trigger objects.

## Details

A zero-crossing trigger watches `signal - level`, where `signal` is an
output of a block. The crossing time is located by the root finder of
the integrator (methods `"lsoda"`, `"lsode"`, `"lsodes"` and
`"lsodar"`), the solver is stopped there, the resets are applied and
integration restarts from the new state. `direction` selects rising
crossings (negative to positive), falling crossings or both. A crossing
in the other direction is ignored. A scheduled trigger fires at fixed
times; the solver always stops exactly at those times.

A reset is a named list. Each name is a block that has state; the value
is either a function `function(x, s, t)` or a string holding a portable
expression. `x` is the current state of that block (continuous states
for continuous blocks, discrete states otherwise), `s` is a named list
of all block outputs at the event time (before the reset), and `t` is
the time. Model parameters are also visible in string expressions. Only
string resets can be saved with
[`sim_write_json()`](https://adamrsophie.github.io/simcompR/reference/sim_write_json.md).

A `signal` or reset target inside a subsystem instance is written as a
path, `"instance/block"` (deeper levels add more `/` parts). A crossing
function that sits exactly on `level` fires when the signal moves away
in the trigger direction, including at the initial time.

When several events fire at the same instant (scheduled events, located
crossings and crossings caused by a jump at that instant) they are
applied in increasing `priority`, then by name, and every reset sees the
state left by the previous one. A reset that makes another event's
signal cross its level fires that event at the same instant. A terminal
event stops the run once the instant is settled. See
[`sim_run()`](https://adamrsophie.github.io/simcompR/reference/sim_run.md)
for the full order of what happens at an instant.

## Examples

``` r
# Bouncing ball: h'' = -g, reverse and damp the velocity at impact.
m <- sim_model("ball", params = list(e = 0.8))
sim_add_block(m, "constant", "gravity", value = -9.81)
sim_add_block(m, "integrator", "v", initial = 0)
sim_add_block(m, "integrator", "h", initial = 1)
sim_connect(m, "gravity", "v")
sim_connect(m, "v", "h")
sim_add_event(m, "impact",
  trigger = sim_zero_crossing("h", direction = "falling"),
  reset = list(v = "-e * x"))
r <- sim_run(sim_compile(m), tspan = c(0, 1.5))
sim_events(r)
#>        time  event          type source direction
#> 1 0.4515236 impact zero_crossing      h   falling
#> 2 1.1739615 impact zero_crossing      h   falling
```
