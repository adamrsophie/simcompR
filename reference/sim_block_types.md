# List registered block types

List registered block types

## Usage

``` r
sim_block_types(category = NULL)
```

## Arguments

- category:

  Optional category filter.

## Value

A data frame with columns `type`, `category`, `description` and
`parameters`.

## Examples

``` r
head(sim_block_types())
#>                  type category
#> 1                 abs     math
#> 2             compare    logic
#> 3            constant   source
#> 4           dead_zone     math
#> 5               demux  routing
#> 6 discrete_integrator discrete
#>                                              description
#> 1                                         Absolute value
#> 2    1 when `u op threshold` holds, else 0 (elementwise)
#> 3                                        Constant signal
#> 4      Zero inside [lower, upper], shifted input outside
#> 5              Split a vector into pieces of given sizes
#> 6 Accumulator with forward, backward or trapezoidal rule
#>                              parameters
#> 1                                      
#> 2                         op, threshold
#> 3                                 value
#> 4                          lower, upper
#> 5                                 sizes
#> 6 gain, initial, period, offset, method
sim_block_types("discrete")
#>                         type category
#> 1        discrete_integrator discrete
#> 2       discrete_state_space discrete
#> 3 discrete_transfer_function discrete
#> 4                 unit_delay discrete
#> 5            zero_order_hold discrete
#>                                              description
#> 1 Accumulator with forward, backward or trapezoidal rule
#> 2       x[k+1] = A x[k] + B u[k], y[k] = C x[k] + D u[k]
#> 3             SISO num(z)/den(z), descending powers of z
#> 4                                          y[k] = u[k-1]
#> 5               Sample the input at each hit and hold it
#>                              parameters
#> 1 gain, initial, period, offset, method
#> 2   A, B, C, D, initial, period, offset
#> 3     num, den, initial, period, offset
#> 4               initial, period, offset
#> 5               period, offset, initial
```
