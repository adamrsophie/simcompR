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

### Roadmap

- Graphical diagram editing.
- Physical (acausal) component libraries and DAE initialization.
- Integration with modelcompR, controlcompR and z3r.
