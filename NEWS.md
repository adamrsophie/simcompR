# simcompR 0.1.0

First release.

* Block-diagram models with named and numbered ports, model parameters and
  editing functions, including `sim_remove_event()`.
* Compilation into immutable plans: size propagation, direct-feedthrough
  graph, algebraic loop detection, deterministic execution order and reuse
  of the structure when only parameters change. A plan keeps the block
  types it was compiled with.
* Continuous simulation with deSolve, stopping exactly at source
  breakpoints, sample hits, events and the kinks of memoryless blocks
  (saturation, dead zone, abs, compare, lookup, minmax and the jumps of
  `sign`, `floor`, `ceiling` and `round`).
* Discrete and multirate execution with exact sample times.
* Scheduled and zero-crossing events, including crossings caused by jumps
  at segment boundaries and crossings that start exactly on their level;
  switch, relay and stop blocks. Top-level events can watch and reset
  blocks inside subsystem instances (`"instance/block"`).
* Optional Newton solution of algebraic loops with residual diagnostics
  and starting values (`sim_compile(loop_guess = )`).
* Hierarchical subsystems with parameter namespaces; misspelled instance
  parameters are errors.
* Parameter expressions: `-sim_param("a") / sim_param("tau")` or
  `sim_param_expr("c(1, 2 * zeta * w, w^2)")`, references inside list
  parameters, and model parameters visible in `expression` blocks and
  state machine guards and actions.
* State machine blocks.
* Autonomous `ode` blocks (`n_inputs = 0`).
* Parameter sweeps and Monte Carlo experiments with fixed seeds.
* Versioned JSON model files.
* Fifteen vignettes and a `benchmarks/` directory.

## Known limitations

* Signals are real vectors of fixed size; no buses, matrices or
  variable-size signals. State machines are flat.
* Each segment between sample hits, breakpoints, events or kinks is one
  `deSolve::ode()` call. On a typical machine a 1 kHz sampled loop costs
  about 0.75 s of CPU time per simulated second (see `benchmarks/`), so
  very fast sample rates over long spans are slow.
* `method = "radau"` and `method = "iteration"` are not supported. Zero
  crossings need `lsoda`, `lsode`, `lsodes` or `lsodar`; other methods
  skip kinks and refuse models with events.

## Roadmap

* Buses and hierarchical state machines.
* Graphical diagram editing.
* Physical (acausal) component libraries and DAE initialization.
* Integration with modelcompR, controlcompR and z3r.
