# Register a block type

Every block in a model is an instance of a registered block type. A type
declares its parameters and a `build` function that turns resolved
parameter values into a block definition: ports, sizes, direct
feedthrough, states, sample time and the functions that compute outputs,
derivatives and discrete updates. The built-in library is registered
when the package loads; `sim_register_block_type()` adds user types.

## Usage

``` r
sim_register_block_type(
  type,
  build,
  params = list(),
  required = names(params)[vapply(params, is.null, TRUE)],
  category = "user",
  description = "",
  portable = TRUE,
  overwrite = FALSE
)
```

## Arguments

- type:

  Type name (a single string).

- build:

  Function of one argument, the named list of resolved parameters,
  returning a block definition.

- params:

  Named list of parameter defaults.

- required:

  Names of parameters the user must supply. By default, every parameter
  whose default is `NULL`.

- category:

  Category used by
  [`sim_block_types()`](https://adamrsophie.github.io/simcompR/reference/sim_block_types.md).

- description:

  One line description.

- portable:

  `TRUE` when every parameter can be written to JSON, so the block can
  be saved with
  [`sim_write_json()`](https://adamrsophie.github.io/simcompR/reference/sim_write_json.md).

- overwrite:

  Allow replacing an existing type.

## Value

The type name, invisibly.

## Details

The list returned by `build(p)` may contain:

- `inputs`, `outputs`: integer vectors of port sizes. `NA` means the
  size is inherited from the connected signal (inputs) or computed by
  `size_out` (outputs). Names, when given, become port names.

- `size_out(in_sizes)`: returns output sizes once input sizes are known.

- `check_sizes(in_sizes, out_sizes)`: optional; signals an error when
  the resolved sizes are inconsistent.

- `feedthrough`: logical, one per input port. `TRUE` when the outputs at
  time `t` depend on that input at time `t`.

- `nx`, `x0`: number and initial value of continuous states.

- `nd`, `d0`: number and initial value of discrete states.

- `sample`: `NULL` for continuous (or inherited) blocks, otherwise
  `c(period, offset)`. Discrete blocks compute outputs only at their
  sample hits and hold them in between.

- `output(t, x, d, u)`: returns a list with one numeric vector per
  output port (a plain vector is accepted for one port). `u` is a list
  of input vectors.

- `deriv(t, x, d, u)`: continuous state derivative.

- `update(t, x, d, u)`: next discrete state, called at sample hits after
  all outputs at that instant have been computed.

- `held0(size, d0)`: initial held output of a discrete block before its
  first sample hit. `size` holds the output port sizes and `d0`, when
  the function accepts it, the initial discrete state (including
  overrides given to
  [`sim_run()`](https://adamrsophie.github.io/simcompR/reference/sim_run.md)).

- `zc(t, x, d, u)` and `on_zc(i, direction, t, x, d, u)`: zero-crossing
  functions and their handler, which returns
  `list(d = , x = , terminal = )`. `zc_direction` gives one of
  `"rising"`, `"falling"`, `"either"` per component.

- `breakpoints(t0, tf)`: times where the block's output is
  discontinuous; the integrator restarts there.

- `modes(t, x, d, u)`: for switching blocks whose output branch is
  stored in the discrete state `d` (like `switch` and `relay`), returns
  the mode for the current inputs. It is called at the start and at
  every segment boundary, repeatedly until all modes settle, never
  inside a segment; the crossing functions in `zc` move the mode inside
  a segment.

- `stop_if(t, x, d, u)`: returns `TRUE` to end the simulation; checked
  at every boundary and after every event (see the `stop` block).

- `kinks(t, x, d, u)` and `kink_size(in_sizes)`: crossing functions
  whose zeros are points where a memoryless output or its slope jumps
  (the corners of a saturation), and their number for the given input
  sizes. The integrator restarts exactly at each kink; nothing is
  logged. They are ignored by methods without root finding.

All functions must be free of side effects: the solver calls them at
trial points that it may later reject.

## Examples

``` r
sim_register_block_type(
  "cube", overwrite = TRUE,
  build = function(p) list(
    inputs = NA, outputs = NA, feedthrough = TRUE,
    size_out = function(n) n,
    output = function(t, x, d, u) u[[1]]^3
  ),
  category = "math", description = "Elementwise cube"
)
m <- sim_model("cube_demo")
sim_add_block(m, "time", "t")
sim_add_block(m, "cube", "c")
sim_connect(m, "t", "c")
r <- sim_run(sim_compile(m), tspan = c(0, 2), times = c(0, 1, 2))
sim_signal(r, "c")
#> [1] 0 1 8
```
