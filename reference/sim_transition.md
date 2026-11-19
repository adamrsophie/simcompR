# State machine transitions

Helper to declare one transition of a `"state_machine"` block.

## Usage

``` r
sim_transition(
  from,
  to,
  guard = "TRUE",
  after = 0,
  action = NULL,
  priority = NA
)
```

## Arguments

- from, to:

  State names.

- guard:

  Guard expression (string) or function; `"TRUE"` by default.

- after:

  Minimum time in `from` before the transition may fire.

- action:

  Optional transition action.

- priority:

  Optional priority; smaller values are tried first.

## Value

A `sim_transition` object.

## Details

A state machine block is a discrete block that executes at every sample
hit of its clock (`period`, `offset`). At a hit it:

1.  collects the transitions leaving the current state, sorted by
    `priority` (smaller first) and then by their position in the list;

2.  takes the first transition whose guard holds and whose `after` time
    (time spent in the current state) has elapsed;

3.  runs the exit action of the old state, the transition action and the
    entry action of the new state, in that order, then records the entry
    time;

4.  when no transition fires, runs the `during` action of the current
    state.

The entry action of the initial state runs at the first sample hit,
before transitions are considered, with the inputs and time of that hit.

At most one transition fires per sample hit, so evaluation never loops.
Outputs are computed after the transition, which means the outputs at a
hit depend on the inputs at the same hit (direct feedthrough).

Guards and actions are portable strings or R functions. A guard string
is an expression in `u` (the input vector), `t`, `elapsed` (time in the
current state), the data variables and the parameters of the enclosing
model or subsystem. An action string assigns data variables, for example
`"count <- count + 1"`; each data variable must still hold one finite
number afterwards. A guard function has signature
`function(u, t, elapsed, data)` and returns `TRUE` or `FALSE`; an action
function has signature `function(data, u, t)` and returns the new data
vector: either unnamed with one value per variable, in the order of
`data`, or named with exactly the data variable names, in any order
(values are matched by name). Guards and actions must be pure functions
of their arguments, because the block evaluates them to compute its
outputs and again to compute its next state.

`after` compares the time spent in the state, `t - entry time`, with the
requested duration, allowing only for the rounding of that subtraction,
so the behavior does not depend on where the time axis starts.

The block's parameters are `states` (character), `initial` (state name),
`transitions` (a list of `sim_transition()`), `entry`, `exit`, `during`
(named lists of actions by state), `data` (named numeric initial values
of data variables), `outputs` (optional named list of numeric output
vectors by state), `n_inputs`, `period` and `offset`. Output ports are
`state` (index of the active state in `states`), then `y` when `outputs`
is given, then `data` when there are data variables.

## See also

[`sim_check_state_machine()`](https://adamrsophie.github.io/simcompR/reference/sim_check_state_machine.md).

## Examples

``` r
# A thermostat that heats below 19 degrees and stops above 21.
m <- sim_model("thermostat")
sim_add_block(m, "integrator", "temp", initial = 18)
sim_add_block(m, "state_machine", "ctrl", period = 0.1,
  states = c("heating", "idle"),
  transitions = list(
    sim_transition("heating", "idle", "u[1] > 21"),
    sim_transition("idle", "heating", "u[1] < 19")),
  outputs = list(heating = 1, idle = 0))
sim_add_block(m, "expression", "rate", expr = "2 * u[1] - 0.1 * (u[2] - 10)")
sim_add_block(m, "mux", "mx")
sim_connect(m, "temp", "ctrl")
sim_connect(m, "ctrl", "mx", output = "y", input = 1)
sim_connect(m, "temp", "mx", input = 2)
sim_connect(m, "mx", "rate")
sim_connect(m, "rate", "temp")
r <- sim_run(sim_compile(m), c(0, 10))
range(sim_signal(r, "temp"))
#> [1] 18.00000 21.02948
```
