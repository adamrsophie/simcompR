# Save and load models as JSON

Models are written as versioned, declarative JSON: block types, names,
parameters, connections, events and nested subsystems. Numbers are
stored as strings with 17 significant digits, so they read back
bit-for-bit, and matrices keep their dimensions.

## Usage

``` r
sim_write_json(model, path, pretty = TRUE)

sim_read_json(path)

sim_model_to_list(model)

sim_model_from_list(x)
```

## Arguments

- model:

  A
  [`sim_model()`](https://adamrsophie.github.io/simcompR/reference/sim_model.md).

- path:

  File path.

- pretty:

  Indent the JSON.

- x:

  A list produced by `sim_model_to_list()`.

## Value

`sim_write_json()` returns `path` invisibly; `sim_read_json()` and
`sim_model_from_list()` return a `sim_model`; `sim_model_to_list()`
returns a list.

## Details

Only portable content can be written. Block types defined by R functions
(`fcn`, `ode`, `signal_function`), user types registered with
`portable = FALSE`, parameters holding functions and event resets given
as functions are rejected with an error naming the offending block,
rather than silently dropped. Use the string forms (`expression` blocks,
string resets and guards), which are parsed against a whitelist of
functions when the model is read, so loading a file never runs arbitrary
R code.

## Examples

``` r
m <- sim_model("feedback", params = list(k = 2))
sim_add_block(m, "step", "r", time = 1)
sim_add_block(m, "subtract", "e")
sim_add_block(m, "gain", "g", gain = sim_param("k"))
sim_add_block(m, "integrator", "x")
sim_connect(m, "r", "e", input = 1)
sim_connect(m, "x", "e", input = 2)
sim_connect(m, "e", "g")
sim_connect(m, "g", "x")
f <- tempfile(fileext = ".json")
sim_write_json(m, f)
m2 <- sim_read_json(f)
r1 <- sim_run(sim_compile(m), c(0, 3), times = 3)
r2 <- sim_run(sim_compile(m2), c(0, 3), times = 3)
identical(sim_signal(r1, "x"), sim_signal(r2, "x"))
#> [1] TRUE
```
