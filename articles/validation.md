# Validation against exact solutions

``` r

library(simcompR)
tight <- function(plan, tspan, ...) sim_run(plan, tspan, rtol = 1e-10, atol = 1e-12, ...)
```

This vignette collects reference problems with known answers. Each one
is also part of the test suite, which runs on every change on Linux,
macOS and Windows. The table at the end is computed when this vignette
is built; it is not copied from elsewhere.

``` r

checks <- list()
```

## Exponential decay

\dot x = -a x, x(0) = 1, so x(t) = e^{-at}.

``` r

m <- sim_model("decay")
sim_add_block(m, "integrator", "x", initial = 1)
sim_add_block(m, "gain", "k", gain = -1.5)
sim_connect(m, "x", "k")
sim_connect(m, "k", "x")
r <- tight(sim_compile(m), c(0, 4), times = seq(0, 4, 0.1))
checks$decay <- max(abs(sim_signal(r, "x") - exp(-1.5 * sim_times(r))))
```

## Harmonic oscillator

\ddot q = -\omega^2 q, q(0) = 1, so q(t) = \cos \omega t.

``` r

o <- sim_model("osc")
sim_add_block(o, "integrator", "q", initial = 1)
sim_add_block(o, "integrator", "v", initial = 0)
sim_add_block(o, "gain", "k", gain = -9)
sim_connect(o, "v", "q")
sim_connect(o, "q", "k")
sim_connect(o, "k", "v")
r <- tight(sim_compile(o), c(0, 10), times = seq(0, 10, 0.05))
checks$oscillator <- max(abs(sim_signal(r, "q") - cos(3 * sim_times(r))))
```

## Linear system and the matrix exponential

``` r

A <- matrix(c(-1, 2, -2, -1), 2, 2)
s <- sim_model("linear")
sim_add_block(s, "constant", "u", value = 0)
sim_add_block(s, "state_space", "g", A = A, B = matrix(0, 2, 1), C = diag(2),
              initial = c(1, 0))
sim_connect(s, "u", "g")
r <- tight(sim_compile(s), c(0, 3), times = c(1, 2, 3))
# exp(A t) for A = [[-1, -2], [2, -1]] is exp(-t) * rotation(2 t)
ex <- t(sapply(c(1, 2, 3), function(t) exp(-t) * c(cos(2 * t), sin(2 * t))))
checks$state_space <- max(abs(sim_signal(r, "g") - ex))
```

## Second-order step response

G(s) = \frac{\omega^2}{s^2 + 2\zeta\omega s + \omega^2} with \zeta =
0.3, \omega = 2.

``` r

z <- 0.3; w <- 2; wd <- w * sqrt(1 - z^2)
tf <- sim_model("tf")
sim_add_block(tf, "step", "u", time = 0)
sim_add_block(tf, "transfer_function", "G", num = w^2, den = c(1, 2 * z * w, w^2))
sim_connect(tf, "u", "G")
r <- tight(sim_compile(tf), c(0, 10), times = seq(0, 10, 0.05))
tt <- sim_times(r)
exact <- 1 - exp(-z * w * tt) * (cos(wd * tt) + z * w / wd * sin(wd * tt))
checks$second_order <- max(abs(sim_signal(r, "G") - exact))
```

## Discrete recurrence

``` r

d <- sim_model("recurrence")
sim_add_block(d, "constant", "u", value = 1)
sim_add_block(d, "discrete_state_space", "x", A = 0.9, B = 0.5, C = 1,
              initial = 0, period = 0.5)
sim_connect(d, "u", "x")
r <- sim_run(sim_compile(d), c(0, 10), times = seq(0, 10, 0.5))
ref <- Reduce(function(x, k) 0.9 * x + 0.5, 1:20, accumulate = TRUE, 0)
checks$recurrence <- max(abs(sim_signal(r, "x") - ref))
```

## Sampled control and exact discretization

``` r

Ts <- 0.25; K <- 1.5
sc <- sim_model("sampled")
sim_add_block(sc, "constant", "ref", value = 1)
sim_add_block(sc, "zero_order_hold", "y_k", period = Ts)
sim_add_block(sc, "subtract", "e")
sim_add_block(sc, "gain", "K", gain = K)
sim_add_block(sc, "zero_order_hold", "u_k", period = Ts)
sim_add_block(sc, "subtract", "rhs")
sim_add_block(sc, "integrator", "x")
sim_connect(sc, "x", "y_k")
sim_connect(sc, "ref", "e", input = 1)
sim_connect(sc, "y_k", "e", input = 2)
sim_connect(sc, "e", "K")
sim_connect(sc, "K", "u_k")
sim_connect(sc, "u_k", "rhs", input = 1)
sim_connect(sc, "x", "rhs", input = 2)
sim_connect(sc, "rhs", "x")
r <- tight(sim_compile(sc), c(0, 5), times = Ts * 0:20)
phi <- exp(-Ts)
ref <- Reduce(function(x, k) phi * x + (1 - phi) * K * (1 - x), 1:20, accumulate = TRUE, 0)
checks$sampled_control <- max(abs(sim_signal(r, "x") - ref))
```

## Bouncing ball impact times

``` r

b <- sim_model("ball")
sim_add_block(b, "constant", "g", value = -9.81)
sim_add_block(b, "integrator", "v")
sim_add_block(b, "integrator", "h", initial = 1)
sim_connect(b, "g", "v")
sim_connect(b, "v", "h")
sim_add_event(b, "impact", sim_zero_crossing("h", direction = "falling"),
              reset = list(v = "-0.7 * x"))
r <- tight(sim_compile(b), c(0, 2))
t1 <- sqrt(2 / 9.81); v <- 0.7 * 9.81 * t1; ref <- t1
for (i in 1:3) { ref <- c(ref, tail(ref, 1) + 2 * v / 9.81); v <- 0.7 * v }
checks$impacts <- max(abs(sim_events(r)$time[1:4] - ref))
```

## Algebraic loop

y = \cos y:

``` r

al <- sim_model("loop")
sim_add_block(al, "math", "c", fun = "cos")
sim_add_block(al, "gain", "y", gain = 1)
sim_connect(al, "c", "y")
sim_connect(al, "y", "c")
r <- sim_run(sim_compile(al, algebraic_loops = "solve"), c(0, 1), times = 1)
checks$algebraic_loop <- abs(sim_signal(r, "y") - 0.7390851332151607)
```

## Summary

``` r

data.frame(problem = names(checks), max_abs_error = signif(unlist(checks), 3),
           row.names = NULL)
#>           problem max_abs_error
#> 1           decay      1.06e-10
#> 2      oscillator      8.30e-10
#> 3     state_space      2.87e-11
#> 4    second_order      3.65e-10
#> 5      recurrence      0.00e+00
#> 6 sampled_control      1.08e-10
#> 7         impacts      4.37e-14
#> 8  algebraic_loop      0.00e+00
```

Continuous errors are controlled by the tolerances (here 10^{-10}
relative, 10^{-12} absolute). Discrete results are exact up to
floating-point round-off. Event times are located by the root finder to
within the same tolerances.
