# Writing custom blocks

``` r

library(simcompR)
```

Every built-in block is registered with the same public function that is
available to users:
[`sim_register_block_type()`](https://adamrsophie.github.io/simcompR/reference/sim_register_block_type.md).
A block type is a `build` function that turns resolved parameters into a
block definition.

## The block definition

`build(p)` receives the named list of parameters and returns a list
with:

| Field | Meaning |
|----|----|
| `inputs` | input port sizes (`NA` for “inferred”); names give port names |
| `outputs` | output port sizes, likewise |
| `size_out` | optional `function(n_in)` giving output sizes from input sizes |
| `feedthrough` | `TRUE`/`FALSE`, or one value per input port |
| `nx`, `x0` | number of continuous states and their initial values |
| `nd`, `d0` | number of discrete states and their initial values |
| `sample` | `c(period, offset)` for discrete blocks |
| `output` | `function(t, x, d, u)` returning the outputs |
| `deriv` | `function(t, x, d, u)` returning \dot x |
| `update` | `function(t, x, d, u)` returning the next discrete state |
| `held0` | `function(size, d0)` giving held outputs before the first hit |
| `zc`, `on_zc`, `zc_direction` | zero-crossing functions and handlers |
| `breakpoints` | times where the block’s output is discontinuous |

`u` is a list with one numeric vector per input port. `output` must be
free of side effects: it is called at trial points of the solver.

The most important field to get right is `feedthrough`. It must be
`TRUE` if `output` reads `u`. Declaring it `FALSE` when the output does
depend on the input would hide algebraic loops and produce wrong
results.

## A stateless block

A dead-band rectifier, y = \max(\|u\| - w, 0):

``` r

sim_register_block_type("deadband_abs", category = "user",
  description = "max(|u| - width, 0)",
  params = list(width = 0.5),
  build = function(p) {
    w <- p$width
    list(inputs = NA, outputs = NA, feedthrough = TRUE,
         size_out = function(n) n,
         output = function(t, x, d, u) pmax(abs(u[[1]]) - w, 0))
  })

m <- sim_model("custom")
sim_add_block(m, "sine", "s", amplitude = 2)
sim_add_block(m, "deadband_abs", "db", width = 1)
sim_connect(m, "s", "db")
plot(sim_run(sim_compile(m), c(0, 10)), c("s", "db"))
```

![](custom-blocks_files/figure-html/unnamed-chunk-3-1.png)

## A continuous block

A first-order low-pass filter with time constant \tau, \tau\dot x = u -
x, y = x. The output does not read the input, so there is no
feedthrough:

``` r

sim_register_block_type("lowpass", category = "user",
  params = list(tau = 1, initial = 0),
  build = function(p) {
    tau <- p$tau
    if (!is.numeric(tau) || length(tau) != 1L || tau <= 0) stop("`tau` must be positive.")
    list(inputs = 1, outputs = 1, feedthrough = FALSE,
         nx = 1, x0 = p$initial,
         output = function(t, x, d, u) x,
         deriv = function(t, x, d, u) (u[[1]] - x) / tau)
  })

f <- sim_model("filter")
sim_add_block(f, "step", "u", time = 0)
sim_add_block(f, "lowpass", "lp", tau = 0.5)
sim_connect(f, "u", "lp")
r <- sim_run(sim_compile(f), c(0, 3), times = 3, rtol = 1e-10, atol = 1e-12)
c(simulated = sim_signal(r, "lp"), exact = 1 - exp(-3 / 0.5))
#> simulated.lp        exact 
#>    0.9975212    0.9975212
```

Because `lowpass` has no feedthrough, it can close a feedback loop on
its own.

## A discrete block

A moving average of the last three samples. Its discrete state holds the
previous two inputs:

``` r

sim_register_block_type("moving_average3", category = "user",
  params = list(period = 1, offset = 0),
  build = function(p) {
    list(inputs = 1, outputs = 1, feedthrough = TRUE,
         nd = 2, d0 = c(0, 0), sample = c(p$period, p$offset),
         output = function(t, x, d, u) (u[[1]] + d[1] + d[2]) / 3,
         update = function(t, x, d, u) c(u[[1]], d[1]),
         held0 = function(size, d0) 0)
  })

ma <- sim_model("ma")
sim_add_block(ma, "time", "t")
sim_add_block(ma, "moving_average3", "avg", period = 1)
sim_connect(ma, "t", "avg")
sim_signal(sim_run(sim_compile(ma), c(0, 5), times = 0:5), "avg")
#> [1] 0.0000000 0.3333333 1.0000000 2.0000000 3.0000000 4.0000000
```

## Portability

Block types whose parameters are plain values can be saved to JSON (see
[`vignette("model-files")`](https://adamrsophie.github.io/simcompR/articles/model-files.md)).
Register a type with `portable = FALSE` when its parameters include R
functions; such models are refused by
[`sim_write_json()`](https://adamrsophie.github.io/simcompR/reference/sim_write_json.md)
instead of being saved incompletely. The type itself must be registered
again, with the same name, in any session that reads the file.

## Listing types

``` r

head(sim_block_types(), 10)
#>                          type category
#> 1                         abs     math
#> 2                     compare    logic
#> 3                    constant   source
#> 4                   dead_zone     math
#> 5                deadband_abs     user
#> 6                       demux  routing
#> 7         discrete_integrator discrete
#> 8        discrete_state_space discrete
#> 9  discrete_transfer_function discrete
#> 10                     divide     math
#>                                               description
#> 1                                          Absolute value
#> 2     1 when `u op threshold` holds, else 0 (elementwise)
#> 3                                         Constant signal
#> 4       Zero inside [lower, upper], shifted input outside
#> 5                                     max(|u| - width, 0)
#> 6               Split a vector into pieces of given sizes
#> 7  Accumulator with forward, backward or trapezoidal rule
#> 8        x[k+1] = A x[k] + B u[k], y[k] = C x[k] + D u[k]
#> 9              SISO num(z)/den(z), descending powers of z
#> 10                                              in1 / in2
#>                               parameters
#> 1                                       
#> 2                          op, threshold
#> 3                                  value
#> 4                           lower, upper
#> 5                                  width
#> 6                                  sizes
#> 7  gain, initial, period, offset, method
#> 8    A, B, C, D, initial, period, offset
#> 9      num, den, initial, period, offset
#> 10
subset(sim_block_types(), category == "user")
#>               type category         description     parameters
#> 5     deadband_abs     user max(|u| - width, 0)          width
#> 19         lowpass     user                       tau, initial
#> 22 moving_average3     user                     period, offset
```
