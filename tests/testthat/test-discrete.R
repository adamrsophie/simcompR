test_that("unit delay outputs the previous input", {
  m <- sim_model()
  sim_add_block(m, "time", "t")
  sim_add_block(m, "unit_delay", "z", period = 1, initial = -1)
  sim_add_block(m, "zero_order_hold", "h", period = 1)
  sim_connect(m, "t", "z")
  sim_connect(m, "t", "h")
  r <- sim_run(sim_compile(m), c(0, 5), times = 0:5)
  expect_equal(unname(sim_signal(r, "h")), 0:5)               # y[k] = u[k]
  expect_equal(unname(sim_signal(r, "z")), c(-1, 0:4))        # y[k] = u[k-1]
  # between hits the held value does not change
  r2 <- sim_run(sim_compile(m), c(0, 2), times = c(0.5, 1.5), include_ticks = FALSE)
  expect_equal(unname(sim_signal(r2, "h")), c(0, 1))
  expect_equal(unname(sim_signal(r2, "z")), c(-1, 0))
})

test_that("first-order recurrence matches an independent loop", {
  # x[k+1] = 0.5 x[k] + u[k], u = 1, y = x
  m <- sim_model()
  sim_add_block(m, "constant", "u", value = 1)
  sim_add_block(m, "discrete_state_space", "g", A = 0.5, B = 1, C = 1, D = 0,
                initial = 4, period = 0.1)
  sim_connect(m, "u", "g")
  r <- sim_run(sim_compile(m), c(0, 1), times = seq(0, 1, 0.1))
  x <- numeric(11); x[1] <- 4
  for (k in 1:10) x[k + 1] <- 0.5 * x[k] + 1
  expect_equal(unname(sim_signal(r, "g")), x)

  # the same system as a z-domain transfer function 1 / (z - 0.5)
  m2 <- sim_model()
  sim_add_block(m2, "constant", "u", value = 1)
  sim_add_block(m2, "discrete_transfer_function", "g", num = 1, den = c(1, -0.5),
                initial = 4, period = 0.1)
  sim_connect(m2, "u", "g")
  r2 <- sim_run(sim_compile(m2), c(0, 1), times = seq(0, 1, 0.1))
  expect_equal(unname(sim_signal(r2, "g")), x)
})

test_that("discrete state-space with feedthrough and vector signals", {
  A <- matrix(c(0.9, 0.1, 0, 0.8), 2, 2)
  B <- matrix(c(1, 0, 0, 1), 2, 2)
  C <- matrix(c(1, 0, 1, 1), 2, 2)
  D <- diag(0.5, 2)
  m <- sim_model()
  sim_add_block(m, "time", "t")
  sim_add_block(m, "mux", "u")
  sim_add_block(m, "math", "s", fun = "sin")
  sim_connect(m, "t", "u", input = 1)
  sim_connect(m, "t", "s")
  sim_connect(m, "s", "u", input = 2)
  sim_add_block(m, "discrete_state_space", "g", A = A, B = B, C = C, D = D,
                initial = c(1, -1), period = 0.25)
  sim_connect(m, "u", "g")
  r <- sim_run(sim_compile(m), c(0, 2), times = seq(0, 2, 0.25))
  x <- c(1, -1)
  y <- matrix(0, 9, 2)
  for (k in 0:8) {
    tk <- 0.25 * k
    uk <- c(tk, sin(tk))
    y[k + 1, ] <- C %*% x + D %*% uk
    x <- A %*% x + B %*% uk
  }
  expect_equal(unname(sim_signal(r, "g")), y, tolerance = 1e-12)
})

test_that("discrete integrator rules follow their recurrences", {
  run <- function(method) {
    m <- sim_model()
    sim_add_block(m, "time", "t")
    sim_add_block(m, "discrete_integrator", "di", gain = 2, period = 0.5,
                  method = method, initial = 1)
    sim_connect(m, "t", "di")
    unname(sim_signal(sim_run(sim_compile(m), c(0, 2), times = seq(0, 2, 0.5)), "di"))
  }
  u <- seq(0, 2, 0.5); kt <- 2 * 0.5
  fwd <- 1 + kt * c(0, cumsum(u)[-5])
  expect_equal(run("forward"), fwd)
  expect_equal(run("backward"), 1 + kt * cumsum(u))
  trap <- numeric(5); x <- 1
  for (k in 1:5) { trap[k] <- x + kt / 2 * u[k]; x <- x + kt * u[k] }
  expect_equal(run("trapezoidal"), trap)
})

test_that("sample hits do not drift for non-representable periods", {
  m <- sim_model()
  sim_add_block(m, "constant", "one")
  sim_add_block(m, "discrete_integrator", "count", period = 0.1, method = "forward",
                gain = 10)
  sim_connect(m, "one", "count")
  r <- sim_run(sim_compile(m), c(0, 100), times = 100, include_ticks = FALSE)
  # 1000 hits in [0, 100): forward rule gives count[k] = k at t = k * 0.1
  expect_equal(unname(sim_signal(r, "count")), 1000)
})

test_that("multirate models follow the reference timeline", {
  # fast sampler at 0.1 (offset 0), slow at 0.25 with offset 0.05
  m <- sim_model()
  sim_add_block(m, "time", "t")
  sim_add_block(m, "zero_order_hold", "fast", period = 0.1)
  sim_add_block(m, "zero_order_hold", "slow", period = 0.25, offset = 0.05)
  sim_add_block(m, "unit_delay", "lag", period = 0.25, offset = 0.05)
  sim_connect(m, "t", "fast")
  sim_connect(m, "fast", "slow")
  sim_connect(m, "fast", "lag")
  times <- seq(0, 1, 0.01)
  r <- sim_run(sim_compile(m), c(0, 1), times = times, include_ticks = FALSE)
  fast_ref <- floor(times / 0.1 + 1e-9) * 0.1
  slow_hit <- floor((times - 0.05) / 0.25 + 1e-9)
  slow_ref <- ifelse(slow_hit < 0, 0, floor((0.05 + slow_hit * 0.25) / 0.1 + 1e-9) * 0.1)
  lag_ref <- ifelse(slow_hit < 1, 0,
                    floor((0.05 + (slow_hit - 1) * 0.25) / 0.1 + 1e-9) * 0.1)
  expect_equal(unname(sim_signal(r, "fast")), fast_ref, tolerance = 1e-12)
  expect_equal(unname(sim_signal(r, "slow")), slow_ref, tolerance = 1e-12)
  expect_equal(unname(sim_signal(r, "lag")), lag_ref, tolerance = 1e-12)
  plan <- sim_compile(m)
  expect_equal(nrow(plan$rates$table), 2L)
  expect_length(plan$diagnostics$notes, 2L)
})

test_that("coincident hits evaluate outputs before updates", {
  # a delay chain at the same rate: each stage lags one more sample
  m <- sim_model()
  sim_add_block(m, "time", "t")
  sim_add_block(m, "unit_delay", "z1", period = 1)
  sim_add_block(m, "unit_delay", "z2", period = 1)
  sim_connect(m, "t", "z1")
  sim_connect(m, "z1", "z2")
  r <- sim_run(sim_compile(m), c(0, 4), times = 0:4)
  expect_equal(unname(sim_signal(r, "z1")), c(0, 0:3))
  expect_equal(unname(sim_signal(r, "z2")), c(0, 0, 0:2))
})

test_that("sampled control of a continuous plant matches exact discretization", {
  # plant x' = -a x + u, controller u[k] = -K x(t_k) + K r, zero-order hold
  a <- 1; K <- 2; Ts <- 0.2; rr <- 1
  m <- sim_model()
  sim_add_block(m, "constant", "ref", value = rr)
  sim_add_block(m, "zero_order_hold", "sample", period = Ts)
  sim_add_block(m, "subtract", "err")
  sim_add_block(m, "gain", "K", gain = K)
  sim_add_block(m, "zero_order_hold", "u", period = Ts)
  sim_add_block(m, "sum", "rhs", signs = "+-")
  sim_add_block(m, "gain", "a", gain = a)
  sim_add_block(m, "integrator", "x", initial = 0)
  sim_connect(m, "x", "sample")
  sim_connect(m, "ref", "err", input = 1)
  sim_connect(m, "sample", "err", input = 2)
  sim_connect(m, "err", "K")
  sim_connect(m, "K", "u")
  sim_connect(m, "u", "rhs", input = 1)
  sim_connect(m, "x", "a")
  sim_connect(m, "a", "rhs", input = 2)
  sim_connect(m, "rhs", "x")
  times <- seq(0, 4, Ts)
  r <- tight(sim_compile(m), c(0, 4), times = times)
  # exact ZOH discretization: x[k+1] = phi x[k] + gam u[k]
  phi <- exp(-a * Ts); gam <- (1 - phi) / a
  x <- numeric(length(times))
  for (k in seq_len(length(times) - 1L)) {
    u <- K * (rr - x[k])
    x[k + 1] <- phi * x[k] + gam * u
  }
  expect_lt(max(abs(sim_signal(r, "x") - x)), 1e-8)
})

test_that("invalid sample times are rejected", {
  m <- sim_model()
  expect_error(sim_add_block(m, "unit_delay", "z", period = 0), "positive")
  expect_error(sim_add_block(m, "unit_delay", "z", period = 1, offset = 1), "offset")
})

test_that("discrete-only models run without an ODE solver", {
  m <- sim_model()
  sim_add_block(m, "constant", "u", value = 2)
  sim_add_block(m, "discrete_integrator", "acc", period = 1)
  sim_connect(m, "u", "acc")
  r <- sim_run(sim_compile(m), c(0, 3), times = 0:3)
  expect_equal(unname(sim_signal(r, "acc")), c(0, 2, 4, 6))
  expect_equal(sim_diagnostics(r)$stats[["segments"]], 0)
})
