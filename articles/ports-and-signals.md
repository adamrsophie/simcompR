# Block ports and signals

``` r

library(simcompR)
```

Every block declares its input and output ports explicitly. This
vignette explains how ports are numbered and named, how signal sizes are
determined and which connections are rejected.

## Ports, numbers and names

Ports are numbered from 1. Most blocks have one output; some have
several (for example `demux`, or `state_machine`). Ports also have
names, which can be used instead of numbers once the block’s parameters
are known.

``` r

m <- sim_model("ports")
sim_add_block(m, "constant", "v", value = c(1, 2, 3))
sim_add_block(m, "demux", "split", sizes = c(1, 2))
sim_add_block(m, "gain", "double", gain = 2)
sim_connect(m, "v", "split")
sim_connect(m, "split", "double", output = "out2")
r <- sim_run(sim_compile(m), c(0, 1), times = 1)
sim_signal(r, c("split.out1", "split.out2", "double"))
#>    split.out1 split.out2[1] split.out2[2]     double[1]     double[2] 
#>             1             2             3             4             6
```

Signals are named after the block that produces them. A block with one
output port is named by the block name (`"double"`); otherwise the port
name is appended (`"split.out2"`). Vector signals have one column per
element.

## The block library

[`sim_block_types()`](https://adamrsophie.github.io/simcompR/reference/sim_block_types.md)
lists every registered type, its category and its parameters.

``` r

types <- sim_block_types()
table(types$category)
#> 
#> continuous   discrete      logic       math    routing       sink     source 
#>          4          5          5         13          3          3         10 
#>  subsystem 
#>          3
types[types$category == "math", c("type", "parameters")]
#>          type        parameters
#> 1         abs                  
#> 4   dead_zone      lower, upper
#> 9      divide                  
#> 10 expression        expr, size
#> 11        fcn         fun, size
#> 13       gain              gain
#> 17     lookup x, y, extrapolate
#> 18       math               fun
#> 19     minmax     fun, n_inputs
#> 24    product               ops
#> 28 saturation      lower, upper
#> 38   subtract                  
#> 39        sum             signs
```

## Signal types

All signals are real numeric vectors of a fixed size. Logical values are
carried as numbers: 1 for true and 0 for false, and logic blocks treat
any nonzero input as true. There is no silent conversion of R logical
vectors or matrices at ports: a block that produces the wrong length is
an error.

## Signal sizes

Sizes are either fixed by a block (a `constant` with three values has an
output of size 3, an integrator has the size of its initial state) or
inherited (a scalar `gain` has the size of its input). Compilation
propagates sizes along connections until every port has a size:

``` r

m <- sim_model("sizes")
sim_add_block(m, "constant", "r", value = c(1, -1))
sim_add_block(m, "sum", "e", signs = "+-")
sim_add_block(m, "integrator", "x", initial = c(0, 0))
sim_connect(m, "r", "e", input = 1)
sim_connect(m, "x", "e", input = 2)
sim_connect(m, "e", "x")
plan <- sim_compile(m)
plan$sizes$out
#> $e
#> [1] 2
#> 
#> $r
#> [1] 2
#> 
#> $x
#> [1] 2
```

### Scalar expansion, and nothing else

Elementwise blocks (`sum`, `product`, `minmax`, `logic`, `switch`)
accept inputs of equal size or scalars. A scalar input is expanded to
the common size. This is the only implicit size rule. R’s usual
recycling of a length 2 vector against a length 3 vector is refused:

``` r

m <- sim_model("recycling")
sim_add_block(m, "constant", "a", value = c(1, 2))
sim_add_block(m, "constant", "b", value = c(1, 2, 3))
sim_add_block(m, "sum", "s")
sim_connect(m, "a", "s", input = 1)
sim_connect(m, "b", "s", input = 2)
sim_compile(m)
#> Error:
#> ! Block 's' (sum): input sizes 2, 3 are incompatible; inputs must have equal sizes or be scalars.
```

A fixed-size input fed by a signal of another size is also refused:

``` r

m <- sim_model("mismatch")
sim_add_block(m, "constant", "a", value = c(1, 2))
sim_add_block(m, "gain", "K", gain = matrix(1, 2, 3))
sim_connect(m, "a", "K")
sim_compile(m)
#> Error:
#> ! Signal size mismatch: 'a' output 1 has size 2 but input 'in1' of block 'K' expects size 3.
```

## Connection rules

The compiler rejects:

- inputs that are not connected;
- an input with more than one source;
- ports or blocks that do not exist;
- incompatible sizes.

``` r

m <- sim_model("open")
sim_add_block(m, "constant", "a")
sim_add_block(m, "product", "p")
sim_connect(m, "a", "p", input = 1)
sim_compile(m)
#> Error:
#> ! Unconnected input port(s): p:in2.
```

An output may feed any number of inputs, and an output that feeds
nothing is fine (its value is still recorded).

## Routing blocks

`mux` concatenates inputs, `demux` splits a vector, and `selector` picks
elements by index:

``` r

m <- sim_model("routing")
sim_add_block(m, "time", "t")
sim_add_block(m, "sine", "s")
sim_add_block(m, "mux", "both", n_inputs = 2)
sim_add_block(m, "selector", "second", index = 2)
sim_connect(m, "t", "both", input = 1)
sim_connect(m, "s", "both", input = 2)
sim_connect(m, "both", "second")
r <- sim_run(sim_compile(m), c(0, 1), times = c(0, 0.5, 1))
cbind(sim_signal(r, "both"), second = sim_signal(r, "second"))
#>      both[1]   both[2]    second
#> [1,]     0.0 0.0000000 0.0000000
#> [2,]     0.5 0.4794255 0.4794255
#> [3,]     1.0 0.8414710 0.8414710
```

## Recording with scopes

Every block output is recorded. A `scope` block marks a signal of
interest: its input is available under the scope’s name and
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) shows scopes by
default.

``` r

sim_add_block(m, "scope", "view")
sim_connect(m, "second", "view")
r <- sim_run(sim_compile(m), c(0, 6))
plot(r)
```

![](ports-and-signals_files/figure-html/unnamed-chunk-10-1.png)

## Units

Signals carry no units in this release. Keep units consistent within a
model; unit checking is on the roadmap.
