# simcompR: block-diagram modeling and dynamic system simulation

`simcompR` builds block diagrams in code, compiles them into a
dependency-aware execution plan and simulates continuous, discrete,
multirate and hybrid dynamics.

## Details

The workflow has three stages:

1.  Build a model with
    [`sim_model()`](https://adamrsophie.github.io/simcompR/reference/sim_model.md),
    [`sim_add_block()`](https://adamrsophie.github.io/simcompR/reference/sim_add_block.md)
    and
    [`sim_connect()`](https://adamrsophie.github.io/simcompR/reference/sim_connect.md).

2.  Compile it with
    [`sim_compile()`](https://adamrsophie.github.io/simcompR/reference/sim_compile.md).
    Compilation validates ports and signal sizes, detects algebraic
    loops and fixes the execution order.

3.  Run the plan with
    [`sim_run()`](https://adamrsophie.github.io/simcompR/reference/sim_run.md)
    and inspect the
    [sim_result](https://adamrsophie.github.io/simcompR/reference/sim_run.md)
    with
    [`sim_signal()`](https://adamrsophie.github.io/simcompR/reference/sim_signal.md),
    [`sim_states()`](https://adamrsophie.github.io/simcompR/reference/sim_signal.md),
    [`sim_events()`](https://adamrsophie.github.io/simcompR/reference/sim_signal.md)
    and [`plot()`](https://rdrr.io/r/graphics/plot.default.html).

See
[`vignette("simcompR")`](https://adamrsophie.github.io/simcompR/articles/simcompR.md)
for a first example and
[`sim_block_types()`](https://adamrsophie.github.io/simcompR/reference/sim_block_types.md)
for the block library.

## See also

Useful links:

- <https://github.com/adamrsophie/simcompR>

- <https://adamrsophie.github.io/simcompR/>

- Report bugs at <https://github.com/adamrsophie/simcompR/issues>

## Author

**Maintainer**: Adam R Sophie <adamrsophie@gmail.com> \[copyright
holder\]

Authors:

- Adam R Sophie <adamrsophie@gmail.com> \[copyright holder\]
