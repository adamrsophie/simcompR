# Check a state machine definition

Static checks of a state machine: unknown states, ambiguous priorities,
states that cannot be reached from the initial state and states without
outgoing transitions. Reachability ignores guards, so a state reported
as reachable may still never be entered.

## Usage

``` r
sim_check_state_machine(states, transitions, initial = states[1])
```

## Arguments

- states:

  Character vector of state names.

- transitions:

  List of
  [`sim_transition()`](https://adamrsophie.github.io/simcompR/reference/sim_transition.md)
  objects.

- initial:

  Initial state.

## Value

A list with `unreachable` and `absorbing` state names. Errors are raised
for structural problems.

## Examples

``` r
sim_check_state_machine(c("a", "b", "c"),
  list(sim_transition("a", "b", "u[1] > 0")), initial = "a")
#> $unreachable
#> [1] "c"
#> 
#> $absorbing
#> [1] "b" "c"
#> 
```
