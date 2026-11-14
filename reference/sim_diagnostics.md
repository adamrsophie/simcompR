# Compilation diagnostics

Compilation diagnostics

## Usage

``` r
sim_diagnostics(x)
```

## Arguments

- x:

  A `sim_plan` or a `sim_result`.

## Value

For a plan, a list with the execution `order`, algebraic `loops`,
direct-feedthrough `edges`, sample `rates` and `notes` about rate
transitions. For a result, the solver status, statistics, warnings and
algebraic loop residuals.

## Examples

``` r
m <- sim_model()
sim_add_block(m, "sine", "u")
sim_add_block(m, "zero_order_hold", "zoh", period = 0.1)
sim_add_block(m, "unit_delay", "z", period = 0.2)
sim_connect(m, "u", "zoh")
sim_connect(m, "zoh", "z")
sim_diagnostics(sim_compile(m))
#> $order
#> [1] "u"   "z"   "zoh"
#> 
#> $loops
#> list()
#> 
#> $edges
#>   from  to
#> 1    u zoh
#> 
#> $rates
#>     period offset blocks
#> zoh    0.1      0    zoh
#> z      0.2      0      z
#> 
#> $notes
#> [1] "'z' (period 0.2) reads 'zoh' (period 0.1): the value is the one held at the reader's sample hit (zero-order hold across rates)."
#> 
```
