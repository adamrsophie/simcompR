# Continuous-time simulation

``` r

library(simcompR)
```

A diagram with integrators, transfer functions, state-space blocks or
general `ode` blocks describes a system

\dot x(t) = f\bigl(t, x(t), u(t), \theta\bigr),

where x stacks the continuous states of all blocks, u are the source
signals and \theta the parameters.
[`sim_run()`](https://adamrsophie.github.io/simcompR/reference/sim_run.md)
integrates it with
[deSolve](https://cran.r-project.org/package=deSolve).

## How the derivative is evaluated

At every call from the solver, at a trial time t and trial state x:

1.  each block, in execution order, reads its inputs from the signal
    vector and computes its outputs from (t, x, \text{inputs});
2.  each block with states computes its state derivative;
3.  the derivatives are returned to the solver.

This function has **no side effects**. Adaptive solvers evaluate the
derivative at trial points that are later rejected, and they do not call
it once per output time. Nothing is appended, counted or committed
during these calls. Recorded signals are recomputed after the solve,
from the accepted states only.

## Output times versus internal steps

`times` controls where results are reported, not how the solver steps:

``` r

m <- sim_model("decay")
sim_add_block(m, "integrator", "x", initial = 1)
sim_add_block(m, "gain", "k", gain = -2)
sim_connect(m, "x", "k")
sim_connect(m, "k", "x")
plan <- sim_compile(m)

coarse <- sim_run(plan, c(0, 3), times = c(0, 3))
fine <- sim_run(plan, c(0, 3), times = seq(0, 3, by = 0.001))
c(coarse = sim_signal(coarse, "x")[2], fine = tail(sim_signal(fine, "x"), 1),
  exact = exp(-6))
#>      coarse        fine       exact 
#> 0.002478752 0.002478752 0.002478752
sim_diagnostics(coarse)$stats
#>        segments           steps rhs_evaluations       jacobians           roots 
#>               1              52             111               0               0
```

Both runs take a similar number of internal steps; the fine run only
interpolates more output points.

## Solvers and tolerances

Any method of
[`deSolve::ode()`](https://rdrr.io/pkg/deSolve/man/ode.html) that can be
told to stop at given times can be used, which is every method except
`"radau"`. Models with zero crossings need a root-finding method:
`"lsoda"`, `"lsode"`, `"lsodes"` or `"lsodar"`. The default is
`"lsoda"`, which switches automatically between stiff and nonstiff
methods. Tolerances are `rtol = 1e-6` and `atol = 1e-8` by default.

``` r

err <- sapply(c(1e-3, 1e-6, 1e-9, 1e-12), function(tol) {
  r <- sim_run(plan, c(0, 3), times = seq(0, 3, 0.1), rtol = tol, atol = tol)
  max(abs(sim_signal(r, "x") - exp(-2 * sim_times(r))))
})
data.frame(tolerance = c(1e-3, 1e-6, 1e-9, 1e-12), max_error = signif(err, 3))
#>   tolerance max_error
#> 1     1e-03  9.65e-04
#> 2     1e-06  1.54e-06
#> 3     1e-09  1.77e-09
#> 4     1e-12  3.91e-12
```

Fixed-step methods are available for comparison:

``` r

r <- sim_run(plan, c(0, 3), times = seq(0, 3, 0.01), method = "rk4")
abs(tail(sim_signal(r, "x"), 1) - exp(-6))
#> [1] 2.016337e-11
```

Pass `hmax` to cap the step size, `maxsteps` to allow more steps per
segment, and any other `deSolve` argument through `...`.

## Continuous blocks

| Block | Equations | Feedthrough |
|----|----|----|
| `integrator` | \dot x = u, y = x | no |
| `transfer_function` | Y(s) = \frac{N(s)}{D(s)} U(s) | only if proper but not strictly proper |
| `state_space` | \dot x = Ax + Bu, y = Cx + Du | only if D \neq 0 |
| `ode` | \dot x = f(t,x,u), y = g(t,x,u) | declared by the user |

### Harmonic oscillator

\ddot q = -\omega^2 q written as two integrators:

``` r

w <- 2
osc <- sim_model("oscillator")
sim_add_block(osc, "integrator", "q", initial = 1)
sim_add_block(osc, "integrator", "qdot", initial = 0)
sim_add_block(osc, "gain", "spring", gain = -w^2)
sim_connect(osc, "qdot", "q")
sim_connect(osc, "q", "spring")
sim_connect(osc, "spring", "qdot")
r <- sim_run(sim_compile(osc), c(0, 10), rtol = 1e-10, atol = 1e-12)
max(abs(sim_signal(r, "q") - cos(w * sim_times(r))))
#> [1] 5.529261e-10
plot(r, c("q", "qdot"))
```

![](continuous-simulation_files/figure-html/unnamed-chunk-6-1.png)

### Transfer functions

A transfer function is realized in controllable canonical form. Its
state vector has `length(den) - 1` entries.

``` r

tf <- sim_model("second_order")
sim_add_block(tf, "step", "u", time = 1)
sim_add_block(tf, "transfer_function", "G", num = 4, den = c(1, 1.2, 4))
sim_connect(tf, "u", "G")
r <- sim_run(sim_compile(tf), c(0, 10))
plot(r, c("u", "G"))
```

![](continuous-simulation_files/figure-html/unnamed-chunk-7-1.png)

### State-space systems and coupled compartments

A two-compartment model with exchange rates, checked against the matrix
exponential:

``` r

A <- matrix(c(-0.7, 0.3, 0.2, -0.4), 2, 2)
cm <- sim_model("compartments")
sim_add_block(cm, "constant", "dose", value = 0)
sim_add_block(cm, "state_space", "body", A = A, B = matrix(0, 2, 1),
              C = diag(2), initial = c(10, 0))
sim_connect(cm, "dose", "body")
r <- sim_run(sim_compile(cm), c(0, 5), times = c(0, 1, 5), rtol = 1e-10, atol = 1e-12)
e <- eigen(A)
exact <- t(sapply(c(0, 1, 5), function(t) {
  Re(e$vectors %*% diag(exp(e$values * t)) %*% solve(e$vectors) %*% c(10, 0))
}))
max(abs(sim_signal(r, "body") - exact))
#> [1] 9.923984e-11
```

### General nonlinear blocks

An `ode` block takes R functions. Here a logistic growth law with the
rate supplied as an input signal:

``` r

lg <- sim_model("logistic")
sim_add_block(lg, "constant", "rate", value = 1.5)
sim_add_block(lg, "ode", "N", initial = 0.05,
              deriv = function(t, x, u) u * x * (1 - x))
sim_connect(lg, "rate", "N")
r <- sim_run(sim_compile(lg), c(0, 8))
plot(r, "N")
```

![](continuous-simulation_files/figure-html/unnamed-chunk-9-1.png)

With `n_inputs = 0` an `ode` block is autonomous, so a system such as
the Lorenz equations needs no dummy input:

``` r

lz <- sim_model("lorenz")
sim_add_block(lz, "ode", "xyz", n_inputs = 0, initial = c(1, 1, 1),
              deriv = function(t, x, u) c(10 * (x[2] - x[1]),
                                          x[1] * (28 - x[3]) - x[2],
                                          x[1] * x[2] - 8 / 3 * x[3]))
r <- sim_run(sim_compile(lz), c(0, 20), times = seq(0, 20, by = 0.01))
xyz <- sim_signal(r, "xyz")
plot(xyz[, 1], xyz[, 3], type = "l", xlab = "x", ylab = "z")
```

![](continuous-simulation_files/figure-html/unnamed-chunk-10-1.png)

`ode` and `fcn` blocks hold R closures; they cannot be saved to portable
JSON files (see
[`vignette("model-files")`](https://adamrsophie.github.io/simcompR/articles/model-files.md)).

## Discontinuous inputs

Sources such as `step`, `ramp`, `piecewise` and `pulse` declare the
times where they jump. The simulator stops the integrator exactly there
and restarts it, rather than integrating across the jump as if the
right-hand side were smooth:

``` r

pw <- sim_model("piecewise")
sim_add_block(pw, "piecewise", "u", times = c(0, 1, 2.5), values = c(1, -2, 0.5))
sim_add_block(pw, "integrator", "x")
sim_connect(pw, "u", "x")
r <- sim_run(sim_compile(pw), c(0, 4), times = c(1, 2.5, 4))
sim_signal(r, "x")  # exactly 1, -2 and -1.25
#> [1]  1.00 -2.00 -1.25
```

## Initial conditions

Initial states come from block parameters and can be overridden per run:

``` r

r <- sim_run(plan, c(0, 1), times = 1, initial = list(x = 5))
sim_signal(r, "x") / exp(-2)
#>        x 
#> 5.000017
```

## Solver failures

When the solver cannot reach the end of the interval, or a signal or
state stops being finite,
[`sim_run()`](https://adamrsophie.github.io/simcompR/reference/sim_run.md)
raises an error naming the cause. With `on_failure = "warn"` it returns
the times actually reached, marked as failed, and never extrapolates:

``` r

blow <- sim_model("blowup")
sim_add_block(blow, "integrator", "x", initial = 1)
sim_add_block(blow, "math", "sq", fun = "square")
sim_connect(blow, "x", "sq")
sim_connect(blow, "sq", "x")          # x' = x^2 explodes at t = 1
r <- suppressWarnings(sim_run(sim_compile(blow), c(0, 2), on_failure = "warn"))
#> DLSODA-  Warning..Internal T (=R1) and H (=R2) are
#>       such that in the machine, T + H = T on the next step  
#>      (H = step size). Solver will continue anyway.
#> In above message, R1 = 0.999998, R2 = 4.61465e-17
#>  
#> DLSODA-  Warning..Internal T (=R1) and H (=R2) are
#>       such that in the machine, T + H = T on the next step  
#>      (H = step size). Solver will continue anyway.
#> In above message, R1 = 0.999998, R2 = 4.61465e-17
#>  
#> DLSODA-  Warning..Internal T (=R1) and H (=R2) are
#>       such that in the machine, T + H = T on the next step  
#>      (H = step size). Solver will continue anyway.
#> In above message, R1 = 0.999998, R2 = 4.61465e-17
#>  
#> DLSODA-  Warning..Internal T (=R1) and H (=R2) are
#>       such that in the machine, T + H = T on the next step  
#>      (H = step size). Solver will continue anyway.
#> In above message, R1 = 0.999998, R2 = 3.68895e-17
#>  
#> DLSODA-  Warning..Internal T (=R1) and H (=R2) are
#>       such that in the machine, T + H = T on the next step  
#>      (H = step size). Solver will continue anyway.
#> In above message, R1 = 0.999998, R2 = 3.68895e-17
#>  
#> DLSODA-  Warning..Internal T (=R1) and H (=R2) are
#>       such that in the machine, T + H = T on the next step  
#>      (H = step size). Solver will continue anyway.
#> In above message, R1 = 0.999998, R2 = 3.0557e-17
#>  
#> DLSODA-  Warning..Internal T (=R1) and H (=R2) are
#>       such that in the machine, T + H = T on the next step  
#>      (H = step size). Solver will continue anyway.
#> In above message, R1 = 0.999998, R2 = 3.0557e-17
#>  
#> DLSODA-  Warning..Internal T (=R1) and H (=R2) are
#>       such that in the machine, T + H = T on the next step  
#>      (H = step size). Solver will continue anyway.
#> In above message, R1 = 0.999998, R2 = 3.0557e-17
#>  
#> DLSODA-  Warning..Internal T (=R1) and H (=R2) are
#>       such that in the machine, T + H = T on the next step  
#>      (H = step size). Solver will continue anyway.
#> In above message, R1 = 0.999998, R2 = 2.44272e-17
#>  
#> DLSODA-  Warning..Internal T (=R1) and H (=R2) are
#>       such that in the machine, T + H = T on the next step  
#>      (H = step size). Solver will continue anyway.
#> In above message, R1 = 0.999998, R2 = 2.44272e-17
#>  
#> DLSODA-  Above warning has been issued I1 times.  
#>      It will not be issued again for this problem.
#> In above message, I1 = 10
#> 
r$status
#> [1] "failed"
r$message
#> [1] "in [0, 2]: Block 'sq' returned Inf on port 1; signals must be finite"
max(sim_times(r))
#> [1] 0
```
