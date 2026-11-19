# Discrete-time simulation

``` r

library(simcompR)
```

Discrete blocks compute

y_k = g(x_k, u_k, k), \qquad x\_{k+1} = f(x_k, u_k, k)

at the sample hits t_k = \text{offset} + k \cdot \text{period} of their
clock, and **hold** their outputs constant between hits.

## Semantics at a sample hit

At each hit time the simulator

1.  evaluates the outputs of all blocks with a hit, in execution order,
    and holds them;
2.  then updates the discrete states of those blocks, using the inputs
    read in step 1.

Outputs before updates is what makes a unit delay a delay:

y\[k\] = x\[k\], \qquad x\[k+1\] = u\[k\] \quad\Longrightarrow\quad
y\[k\] = u\[k-1\].

``` r

m <- sim_model("delay")
sim_add_block(m, "time", "t")
sim_add_block(m, "zero_order_hold", "hold", period = 1)
sim_add_block(m, "unit_delay", "delay", period = 1, initial = -1)
sim_connect(m, "t", "hold")
sim_connect(m, "t", "delay")
r <- sim_run(sim_compile(m), c(0, 5), times = 0:5)
cbind(t = sim_times(r), hold = sim_signal(r, "hold"), delay = sim_signal(r, "delay"))
#>      t hold delay
#> [1,] 0    0    -1
#> [2,] 1    1     0
#> [3,] 2    2     1
#> [4,] 3    3     2
#> [5,] 4    4     3
#> [6,] 5    5     4
```

The zero-order hold has direct feedthrough: `hold[k] = u[k]`. The delay
does not: `delay[k] = u[k-1]`, with the initial value before the first
input.

Values recorded at a sample hit are the values **after** the hit
(right-continuous). Between hits, held values do not change:

``` r

r <- sim_run(sim_compile(m), c(0, 3), times = seq(0, 3, 0.25))
plot(r, c("t", "hold", "delay"))
```

![](discrete-time_files/figure-html/unnamed-chunk-4-1.png)

Signals of discrete blocks are drawn as steps.

## Exact sample times

Sample hits are computed as `offset + k * period` with an integer `k`,
never by repeatedly adding `period`. A period of 0.1, which has no exact
binary representation, produces exactly 1000 hits in \[0, 100):

``` r

acc <- sim_model("count")
sim_add_block(acc, "constant", "one")
sim_add_block(acc, "discrete_integrator", "n", period = 0.1, gain = 10)
sim_connect(acc, "one", "n")
r <- sim_run(sim_compile(acc), c(0, 100), times = 100, include_ticks = FALSE)
sim_signal(r, "n")
#>    n 
#> 1000
```

## Discrete blocks

| Block | Output | Update | Feedthrough |
|----|----|----|----|
| `unit_delay` | y = x | x^+ = u | no |
| `zero_order_hold` | y = u |  | yes |
| `discrete_integrator` (forward) | y = x | x^+ = x + KTu | no |
| `discrete_integrator` (backward) | y = x + KTu | x^+ = y | yes |
| `discrete_integrator` (trapezoidal) | y = x + \frac{KT}{2}u | x^+ = y + \frac{KT}{2}u | yes |
| `discrete_state_space` | y = Cx + Du | x^+ = Ax + Bu | if D \neq 0 |
| `discrete_transfer_function` | realization of N(z)/D(z) |  | if proper |
| `white_noise` | seeded Gaussian draw | sample counter | no |
| `state_machine` | see [`vignette("state-machines")`](https://adamrsophie.github.io/simcompR/articles/state-machines.md) |  | yes |

## A first-order recurrence

x\_{k+1} = 0.5\\x_k + u_k, checked against a plain R loop:

``` r

rec <- sim_model("recurrence")
sim_add_block(rec, "constant", "u", value = 1)
sim_add_block(rec, "discrete_state_space", "x", A = 0.5, B = 1, C = 1,
              initial = 4, period = 0.1)
sim_connect(rec, "u", "x")
r <- sim_run(sim_compile(rec), c(0, 1), times = seq(0, 1, 0.1))
ref <- Reduce(function(x, k) 0.5 * x + 1, 1:10, accumulate = TRUE, 4)
all.equal(unname(sim_signal(r, "x")), ref)
#> [1] TRUE
```

The same system as a transfer function in z, \frac{1}{z - 0.5}:

``` r

tfz <- sim_model("tfz")
sim_add_block(tfz, "constant", "u", value = 1)
sim_add_block(tfz, "discrete_transfer_function", "x", num = 1, den = c(1, -0.5),
              initial = 4, period = 0.1)
sim_connect(tfz, "u", "x")
all.equal(unname(sim_signal(sim_run(sim_compile(tfz), c(0, 1), times = seq(0, 1, 0.1)), "x")), ref)
#> [1] TRUE
```

## Sampled control of a continuous plant

A plant \dot x = -x + u is controlled by u_k = K(r - x(t_k)), held
between samples. The exact zero-order-hold discretization gives the
reference x\_{k+1} = \phi x_k + \gamma u_k with \phi = e^{-T}, \gamma =
1 - e^{-T}.

``` r

K <- 2; Ts <- 0.2
ctl <- sim_model("sampled_control")
sim_add_block(ctl, "constant", "ref", value = 1)
sim_add_block(ctl, "zero_order_hold", "sensor", period = Ts)
sim_add_block(ctl, "subtract", "error")
sim_add_block(ctl, "gain", "K", gain = K)
sim_add_block(ctl, "zero_order_hold", "u", period = Ts)
sim_add_block(ctl, "subtract", "rhs")
sim_add_block(ctl, "integrator", "x")
sim_connect(ctl, "x", "sensor")
sim_connect(ctl, "ref", "error", input = 1)
sim_connect(ctl, "sensor", "error", input = 2)
sim_connect(ctl, "error", "K")
sim_connect(ctl, "K", "u")
sim_connect(ctl, "u", "rhs", input = 1)
sim_connect(ctl, "x", "rhs", input = 2)
sim_connect(ctl, "rhs", "x")
r <- sim_run(sim_compile(ctl), c(0, 4), times = Ts * 0:20, rtol = 1e-10, atol = 1e-12)

phi <- exp(-Ts); gam <- 1 - phi
xk <- Reduce(function(x, k) phi * x + gam * K * (1 - x), 1:20, accumulate = TRUE, 0)
max(abs(sim_signal(r, "x") - xk))
#> [1] 1.681196e-10
plot(sim_run(sim_compile(ctl), c(0, 4)), c("x", "u"))
```

![](discrete-time_files/figure-html/unnamed-chunk-8-1.png)

Between hits the plant evolves continuously while `u` is held: the
continuous state is integrated over each sample interval with the held
input, and the integrator restarts at every hit.

## Discrete-only models

Models with no continuous states never call the ODE solver:

``` r

sim_diagnostics(sim_run(sim_compile(rec), c(0, 1)))$stats
#>        segments           steps rhs_evaluations       jacobians           roots 
#>               0               0               0              NA               0
```
