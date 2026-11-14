# Changelog

## simcompR 0.1.0

First release.

- Block-diagram models with named and numbered ports, model parameters
  and editing functions.
- Compilation into immutable plans: size propagation, direct-feedthrough
  graph, algebraic loop detection, deterministic execution order and
  reuse of the structure when only parameters change.
- Continuous simulation with deSolve, stopping exactly at source
  breakpoints, sample hits and events.
- Discrete and multirate execution with exact sample times.
- Scheduled and zero-crossing events, switch, relay and stop blocks.
- Optional Newton solution of algebraic loops with residual diagnostics.
- Hierarchical subsystems with parameter namespaces.
- State machine blocks.
- Parameter sweeps and Monte Carlo experiments with fixed seeds.
- Versioned JSON model files.
- Fifteen vignettes.

### Known limitations

- Kinks inside memoryless blocks (`saturation`, `dead_zone`, `abs`,
  `compare`, `lookup`) have no crossing function; put an event on the
  input when the exact time matters.
- Top-level events can watch a subsystem output but cannot reset a state
  inside a subsystem; define such events inside the subsystem.
- Algebraic loops are solved by Newton iteration with a
  finite-difference Jacobian, started from zero at each evaluation;
  there is no initial guess argument yet.
- Signals are real vectors of fixed size; no buses, matrices or
  variable-size signals. State machines are flat.
- [`sim_param()`](https://adamrsophie.github.io/simcompR/reference/sim_param.md)
  references cannot be combined in arithmetic (write `-a` with a gain
  block or a separate parameter).
- A zero-crossing event whose function starts a run exactly on its level
  does not fire at the start, even when it then moves in the trigger
  direction.
- Block parameters whose names are a prefix of an argument of
  [`sim_add_block()`](https://adamrsophie.github.io/simcompR/reference/sim_add_block.md)
  (`model`, `type`, `name`) or
  [`sim_set_params()`](https://adamrsophie.github.io/simcompR/reference/sim_set_params.md)
  (`plan`) are captured by R’s partial matching; avoid names such as `n`
  or `p`.
- One [`deSolve::ode()`](https://rdrr.io/pkg/deSolve/man/ode.html) call
  per segment: fast sample rates over long spans are slow.
- `method = "radau"` and `method = "iteration"` are not supported; zero
  crossings need `lsoda`, `lsode`, `lsodes` or `lsodar`.

### Roadmap

- Buses, hierarchical state machines, crossing functions for saturation
  and similar blocks, loop initial guesses, a benchmarks directory.

- Graphical diagram editing.

- Physical (acausal) component libraries and DAE initialization.

- Integration with modelcompR, controlcompR and z3r.
