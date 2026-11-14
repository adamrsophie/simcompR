# Plot simulation results

Draws one panel per signal against time with base graphics. Signals
produced by discrete blocks are drawn as steps (held between sample
hits); continuous signals as lines. Event times are marked with dotted
vertical lines. By default the inputs of `scope` blocks are shown, or,
without scopes, the outputs of continuous and discrete dynamic blocks.

## Usage

``` r
# S3 method for class 'sim_result'
plot(x, signals = NULL, events = TRUE, ...)
```

## Arguments

- x:

  A `sim_result`.

- signals:

  Signal names to plot.

- events:

  Mark event times.

- ...:

  Passed to
  [`graphics::plot()`](https://rdrr.io/r/graphics/plot.default.html).

## Value

`x`, invisibly.

## Examples

``` r
m <- sim_model()
sim_add_block(m, "sine", "u")
sim_add_block(m, "zero_order_hold", "zoh", period = 0.5)
sim_connect(m, "u", "zoh")
r <- sim_run(sim_compile(m), c(0, 6))
plot(r, c("u", "zoh"))
```
