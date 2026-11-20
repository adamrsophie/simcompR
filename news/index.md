# Changelog

## simcompR 0.1.0

First release.

- Block-diagram models with named and numbered ports, model parameters
  and editing functions, including
  [`sim_remove_event()`](https://adamrsophie.github.io/simcompR/reference/sim_add_event.md).
- Compilation into immutable plans: size propagation, direct-feedthrough
  graph, algebraic loop detection, deterministic execution order and
  reuse of the structure when only parameters change. A plan keeps the
  block types it was compiled with.
- Continuous simulation with deSolve, stopping exactly at source
  breakpoints, sample hits, events and the kinks of memoryless blocks
  (saturation, dead zone, abs, compare, lookup, minmax and the jumps of
  `sign`, `floor`, `ceiling` and `round`). The outputs of `compare`,
  `sign`, `floor`, `ceiling` and `round` jump at the located instant, so
  switches, stop blocks and events downstream see the jump; the value
  there is the exact one-sided limit, so results do not change when
  inputs and thresholds are far from zero. The run’s `zc_tol` is the one
  band for crossings, kinks and switch modes. A source that jumps at the
  final time is seen from the left by the last step. The fixed-step
  methods `"euler"` and `"rk4"` step by `hmax` (by default
  `(tf - t0) / 1000`), not by the output grid. A block failure inside a
  segment keeps the output times reached before it. A `jacfunc` passed
  to deSolve receives the absolute time.
- Discrete and multirate execution with exact sample times that do not
  depend on where the time axis starts.
- Scheduled and zero-crossing events, including crossings caused by
  jumps at segment boundaries and departures from a level (also with
  zero first derivative); switch, relay and stop blocks. All events of
  an instant run in one priority order; every handler, reset or sample
  update that carries a function across its level fires the events
  watching it at the same instant, and terminal events stop once the
  instant is settled. Every handler, stop blocks included, counts toward
  `max_events`, and an exhausted budget is reported as `"event_limit"`.
  Events can watch and reset blocks inside subsystem instances at any
  depth (`"outer/inner/block"`); resets inside a subsystem read its
  signals by their local names.
- Time instants are grouped once, so a sample hit, a breakpoint and a
  scheduled event computed a few units in the last place apart meet at
  one boundary at any time origin, and only blocks that jump at a
  boundary are evaluated at their left limit there.
- Block outputs and states must stay finite; a non-finite value fails
  the run with a message naming its source.
- Optional Newton solution of algebraic loops with residual diagnostics
  and starting values (`sim_compile(loop_guess = )`).
- Hierarchical subsystems with parameter namespaces; misspelled instance
  parameters are errors.
- Parameter expressions: `-sim_param("a") / sim_param("tau")` or
  `sim_param_expr("c(1, 2 * zeta * w, w^2)")`, references inside list
  parameters, and model parameters visible in `expression` blocks and
  state machine guards and actions.
- State machine blocks; `after` measures the true time in the state,
  also for states entered between sample hits.
- Autonomous `ode` blocks (`n_inputs = 0`).
- Parameter sweeps and Monte Carlo experiments with fixed seeds;
  summaries are matched by name, and `white_noise` uses its own
  counter-based generator.
- Versioned JSON model files that keep every double bit for bit (short
  decimals where exact, C99 hexadecimal otherwise), integer, missing and
  named values and integer matrices, and refuse values and files they
  cannot represent, including numbers out of the range of a double and
  encoded values with more than one type tag. Member order does not
  matter. A parameter given twice to a block is refused everywhere.
- Sweep tables keep the class of scalar parameters (for example `Date`).
- Fifteen vignettes and a `benchmarks/` directory.

### Known limitations

- Signals are real vectors of fixed size; no buses, matrices or
  variable-size signals. State machines are flat.
- Each segment between sample hits, breakpoints, events or kinks is one
  [`deSolve::ode()`](https://rdrr.io/pkg/deSolve/man/ode.html) call. In
  the measurements of `benchmarks/run.R` a 1 kHz sampled loop cost about
  0.7 to 0.8 s of CPU time per simulated second, so very fast sample
  rates over long spans are slow. Locating kinks adds one segment per
  kink (about twice the cost of the saturated sine benchmark, for the
  accuracy gained).
- With `RNGkind(normal.kind = "Box-Muller")`,
  [`sim_monte_carlo()`](https://adamrsophie.github.io/simcompR/reference/sim_sweep.md)
  cannot restore the cached normal deviate of the caller’s stream and
  warns.
- `method = "radau"` and `method = "iteration"` are not supported, and
  neither is `method = "ode23"`: its implementation in deSolve (1.40 and
  1.42) can take one uncontrolled step across a whole segment and return
  a wrong result. Zero crossings need `lsoda`, `lsode`, `lsodes` or
  `lsodar`; other methods skip kinks and refuse models with events.
  deSolve arguments that set times, events or root functions are
  refused.
- A function that comes within `zc_tol` of its level counts as reaching
  it: `exp(-t)` “crosses” zero once it is below `zc_tol`, at the next
  boundary. Set `zc_tol = 0` for exact sign changes only.

### Roadmap

- Buses and hierarchical state machines.
- Graphical diagram editing.
- Physical (acausal) component libraries and DAE initialization.
- Integration with modelcompR, controlcompR and z3r.
