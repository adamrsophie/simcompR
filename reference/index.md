# Package index

## Building models

- [`sim_model()`](https://adamrsophie.github.io/simcompR/reference/sim_model.md)
  : Create a block-diagram model
- [`sim_add_block()`](https://adamrsophie.github.io/simcompR/reference/sim_add_block.md)
  : Add a block to a model
- [`sim_connect()`](https://adamrsophie.github.io/simcompR/reference/sim_connect.md)
  : Connect two blocks
- [`sim_disconnect()`](https://adamrsophie.github.io/simcompR/reference/sim_disconnect.md)
  [`sim_remove_block()`](https://adamrsophie.github.io/simcompR/reference/sim_disconnect.md)
  : Remove a connection or a block
- [`sim_set_block_params()`](https://adamrsophie.github.io/simcompR/reference/sim_set_block_params.md)
  : Change block parameters in a model
- [`sim_blocks()`](https://adamrsophie.github.io/simcompR/reference/sim_blocks.md)
  [`sim_connections()`](https://adamrsophie.github.io/simcompR/reference/sim_blocks.md)
  : Inspect the blocks and connections of a model
- [`sim_copy()`](https://adamrsophie.github.io/simcompR/reference/sim_copy.md)
  : Copy a model
- [`sim_param()`](https://adamrsophie.github.io/simcompR/reference/sim_param.md)
  [`sim_param_expr()`](https://adamrsophie.github.io/simcompR/reference/sim_param.md)
  : Refer to a model or subsystem parameter
- [`sim_block_types()`](https://adamrsophie.github.io/simcompR/reference/sim_block_types.md)
  : List registered block types

## Subsystems, events and state machines

- [`sim_subsystem()`](https://adamrsophie.github.io/simcompR/reference/sim_subsystem.md)
  : Hierarchical subsystems
- [`sim_add_event()`](https://adamrsophie.github.io/simcompR/reference/sim_add_event.md)
  [`sim_remove_event()`](https://adamrsophie.github.io/simcompR/reference/sim_add_event.md)
  [`sim_zero_crossing()`](https://adamrsophie.github.io/simcompR/reference/sim_add_event.md)
  [`sim_at()`](https://adamrsophie.github.io/simcompR/reference/sim_add_event.md)
  : Model events: zero crossings, scheduled times and state resets
- [`sim_transition()`](https://adamrsophie.github.io/simcompR/reference/sim_transition.md)
  : State machine transitions
- [`sim_check_state_machine()`](https://adamrsophie.github.io/simcompR/reference/sim_check_state_machine.md)
  : Check a state machine definition

## Compiling and running

- [`sim_compile()`](https://adamrsophie.github.io/simcompR/reference/sim_compile.md)
  : Compile a model into an execution plan
- [`sim_set_params()`](https://adamrsophie.github.io/simcompR/reference/sim_set_params.md)
  [`sim_param_values()`](https://adamrsophie.github.io/simcompR/reference/sim_set_params.md)
  : Change parameters of a compiled plan
- [`sim_run()`](https://adamrsophie.github.io/simcompR/reference/sim_run.md)
  : Simulate a compiled plan
- [`sim_diagnostics()`](https://adamrsophie.github.io/simcompR/reference/sim_diagnostics.md)
  : Compilation diagnostics

## Results

- [`sim_signal()`](https://adamrsophie.github.io/simcompR/reference/sim_signal.md)
  [`sim_states()`](https://adamrsophie.github.io/simcompR/reference/sim_signal.md)
  [`sim_discrete_states()`](https://adamrsophie.github.io/simcompR/reference/sim_signal.md)
  [`sim_events()`](https://adamrsophie.github.io/simcompR/reference/sim_signal.md)
  [`sim_times()`](https://adamrsophie.github.io/simcompR/reference/sim_signal.md)
  : Extract simulation output
- [`plot(`*`<sim_result>`*`)`](https://adamrsophie.github.io/simcompR/reference/plot.sim_result.md)
  : Plot simulation results

## Experiments

- [`sim_sweep()`](https://adamrsophie.github.io/simcompR/reference/sim_sweep.md)
  [`sim_monte_carlo()`](https://adamrsophie.github.io/simcompR/reference/sim_sweep.md)
  : Parameter sweeps and Monte Carlo experiments

## Files, diagrams and extension

- [`sim_write_json()`](https://adamrsophie.github.io/simcompR/reference/sim_write_json.md)
  [`sim_read_json()`](https://adamrsophie.github.io/simcompR/reference/sim_write_json.md)
  [`sim_model_to_list()`](https://adamrsophie.github.io/simcompR/reference/sim_write_json.md)
  [`sim_model_from_list()`](https://adamrsophie.github.io/simcompR/reference/sim_write_json.md)
  : Save and load models as JSON
- [`sim_to_dot()`](https://adamrsophie.github.io/simcompR/reference/sim_to_dot.md)
  [`plot(`*`<sim_model>`*`)`](https://adamrsophie.github.io/simcompR/reference/sim_to_dot.md)
  : Draw or export a block diagram
- [`sim_register_block_type()`](https://adamrsophie.github.io/simcompR/reference/sim_register_block_type.md)
  : Register a block type

## Package

- [`simcompR`](https://adamrsophie.github.io/simcompR/reference/simcompR-package.md)
  [`simcompR-package`](https://adamrsophie.github.io/simcompR/reference/simcompR-package.md)
  : simcompR: block-diagram modeling and dynamic system simulation
