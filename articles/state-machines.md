# State machines

``` r

library(simcompR)
```

The `state_machine` block is a finite state machine inspired by
Stateflow charts. It is a discrete block: it runs at the sample hits of
its clock and holds its outputs in between.

## Semantics

At every sample hit the block:

1.  collects the transitions leaving the active state, sorted by
    `priority` (smaller first) and then by their position in the list;
2.  fires the first one whose guard is true and whose `after` time (time
    spent in the active state) has elapsed;
3.  runs, in order, the exit action of the old state, the transition
    action and the entry action of the new state, and records the entry
    time;
4.  if no transition fires, runs the `during` action of the active
    state.

At most one transition fires per hit, so evaluation always terminates,
and the ordering never depends on R list iteration details. Outputs are
computed after the transition: the outputs at a hit reflect the inputs
at that hit.

Guards and actions are strings (portable, checked against a whitelist)
or R functions. Guards see `u` (the input vector), `t`, `elapsed` and
the data variables; actions assign data variables.

Guards and actions must be pure: their result may depend only on their
arguments. The block evaluates them once to compute its outputs at a hit
and again to compute its next state, and the event engine may evaluate
them at more points of one instant, so a function with side effects
(printing, drawing random numbers, counting calls) would see more calls
than hits.

### Elapsed time

`after` compares the true time in the state, `t` minus the entry time,
with the requested duration. The only allowance is for the rounding of
the two times, at most 2\varepsilon\max(\|t\|, \|t\_\text{entry}\|), so
a wait of one period far from the origin fires on time. A state can be
entered between sample hits, through an initial override or an event
reset; it then waits its full duration and fires at the first hit after
it has elapsed. Here the machine enters `a` at t = 0.4 and its one
second wait fires at t = 2, not at the hit at t = 1:

``` r

fsm <- sim_model("late_entry")
sim_add_block(fsm, "constant", "c")
sim_add_block(fsm, "state_machine", "sm", period = 1, states = c("a", "b"),
              transitions = list(sim_transition("a", "b", after = 1)))
sim_connect(fsm, "c", "sm")
r <- sim_run(sim_compile(fsm), c(0.4, 3), times = c(0.4, 1, 2, 3),
             include_ticks = FALSE, initial = list(sm = c(1, 0.4)))
sim_signal(r, "sm")
#> [1] 1 1 2 2
```

The discrete states of the block are the active state index, the entry
time and the data variables (see
[`sim_discrete_states()`](https://adamrsophie.github.io/simcompR/reference/sim_signal.md)).
Before the first hit the entry time is `-Inf`, meaning “not entered
yet”; this is the one discrete state allowed to be non-finite.

## A timed traffic light

``` r

m <- sim_model("traffic")
sim_add_block(m, "constant", "enable", value = 1)
sim_add_block(m, "state_machine", "light", period = 1,
  states = c("red", "green", "yellow"),
  transitions = list(
    sim_transition("red", "green", after = 3),
    sim_transition("green", "yellow", after = 2),
    sim_transition("yellow", "red", after = 1)),
  outputs = list(red = 1, green = 2, yellow = 3))
sim_connect(m, "enable", "light")
r <- sim_run(sim_compile(m), c(0, 12), times = 0:12)
data.frame(t = sim_times(r), state = c("red", "green", "yellow")[sim_signal(r, "light.state")])
#>     t  state
#> 1   0    red
#> 2   1    red
#> 3   2    red
#> 4   3  green
#> 5   4  green
#> 6   5 yellow
#> 7   6    red
#> 8   7    red
#> 9   8    red
#> 10  9  green
#> 11 10  green
#> 12 11 yellow
#> 13 12    red
```

The output ports are `state` (the index of the active state), `y` (the
per-state `outputs`, when given) and `data` (when there are data
variables).

## Guards, priorities and data

A counter of how many times the machine re-entered state `a`:

``` r

m <- sim_model("counter")
sim_add_block(m, "time", "t")
sim_add_block(m, "state_machine", "sm", period = 1,
  states = c("a", "b"),
  transitions = list(
    sim_transition("a", "b", "u[1] >= 2", action = "n <- n + 10"),
    sim_transition("b", "a", "TRUE")),
  entry = list(a = "n <- n + 1"),
  data = c(n = 0))
sim_connect(m, "t", "sm")
r <- sim_run(sim_compile(m), c(0, 6), times = 0:6)
cbind(t = sim_times(r), state = sim_signal(r, "sm.state"), n = sim_signal(r, "sm.data"))
#>      t state  n
#> [1,] 0     1  1
#> [2,] 1     1  1
#> [3,] 2     2 11
#> [4,] 3     1 12
#> [5,] 4     2 22
#> [6,] 5     1 23
#> [7,] 6     2 33
```

The entry action of the initial state runs once, at the first hit, with
the inputs and time of that hit, before any transition is tried.

## A thermostat in closed loop

A heater controlled by a two-state machine with hysteresis. The room
cools toward 10 degrees and the heater adds 2 degrees per unit time.

``` r

th <- sim_model("thermostat")
sim_add_block(th, "integrator", "temp", initial = 18)
sim_add_block(th, "state_machine", "ctrl", period = 0.1,
  states = c("heating", "idle"),
  transitions = list(
    sim_transition("heating", "idle", "u[1] > 21"),
    sim_transition("idle", "heating", "u[1] < 19")),
  outputs = list(heating = 1, idle = 0))
sim_add_block(th, "mux", "mx")
sim_add_block(th, "expression", "rate", expr = "2 * u[1] - 0.1 * (u[2] - 10)")
sim_connect(th, "temp", "ctrl")
sim_connect(th, "ctrl", "mx", output = "y", input = 1)
sim_connect(th, "temp", "mx", input = 2)
sim_connect(th, "mx", "rate")
sim_connect(th, "rate", "temp")
r <- sim_run(sim_compile(th), c(0, 30))
plot(r, c("temp", "ctrl.y"))
```

![](state-machines_files/figure-html/unnamed-chunk-6-1.png)

The machine samples every 0.1 time units, so the temperature overshoots
the thresholds slightly; this is the expected behavior of a sampled
controller.

## Static checks

[`sim_check_state_machine()`](https://adamrsophie.github.io/simcompR/reference/sim_check_state_machine.md)
reports unreachable states (ignoring guards) and states without outgoing
transitions, and errors on unknown states or ambiguous priorities. The
block runs the structural checks when it is built.

``` r

sim_check_state_machine(c("off", "on", "fault"),
  list(sim_transition("off", "on", "u[1] > 0"),
       sim_transition("on", "off", "u[1] <= 0")),
  initial = "off")
#> $unreachable
#> [1] "fault"
#> 
#> $absorbing
#> [1] "fault"
sim_check_state_machine(c("a", "b"),
  list(sim_transition("a", "b", priority = 1),
       sim_transition("a", "a", priority = 1)))
#> Error:
#> ! State 'a' has transitions with the same priority; the order would be ambiguous.
```

Formal verification of state machine properties (for example with the
sibling package ‘z3r’) is planned, not part of this release.
