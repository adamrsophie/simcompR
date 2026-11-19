test_that("negative feedback integrator matches 1 - exp(-t)", {
  m <- connect_feedback(feedback_model(1:3))
  times <- seq(0, 10, by = 0.01)
  r <- tight(sim_compile(m), c(0, 10), times = times)
  expect_true(r$success)
  expect_equal(r$time, times)
  expect_lt(max(abs(sim_signal(r, "plant") - (1 - exp(-times)))), 1e-8)
  # the error signal is recorded consistently with the state
  expect_lt(max(abs(sim_signal(r, "error") - exp(-times))), 1e-8)
})

test_that("first-order decay is exact for several rates", {
  for (a in c(0.1, 1, 5)) {
    r <- tight(sim_compile(decay_model(a, x0 = 2)), c(0, 3), times = seq(0, 3, 0.1))
    expect_lt(max(abs(sim_signal(r, "x") - 2 * exp(-a * r$time))), 1e-8)
  }
})

test_that("harmonic oscillator follows cos and sin", {
  w <- 2.5
  m <- sim_model("osc")
  sim_add_block(m, "integrator", "x1", initial = 1)
  sim_add_block(m, "integrator", "x2", initial = 0)
  sim_add_block(m, "gain", "k", gain = -w^2)
  sim_connect(m, "x2", "x1")
  sim_connect(m, "x1", "k")
  sim_connect(m, "k", "x2")
  r <- tight(sim_compile(m), c(0, 10), times = seq(0, 10, 0.05))
  expect_lt(max(abs(sim_signal(r, "x1") - cos(w * r$time))), 1e-7)
  expect_lt(max(abs(sim_signal(r, "x2") + w * sin(w * r$time))), 1e-6)
})

test_that("two-compartment system agrees with the matrix exponential", {
  A <- matrix(c(-0.7, 0.3, 0.2, -0.4), 2, 2)
  x0 <- c(10, 0)
  m <- sim_model("compartments")
  sim_add_block(m, "state_space", "sys", A = A, B = matrix(0, 2, 1), C = diag(2),
                initial = x0)
  sim_add_block(m, "constant", "zero", value = 0)
  sim_connect(m, "zero", "sys")
  times <- c(0, 0.5, 1, 2, 5)
  r <- tight(sim_compile(m), c(0, 5), times = times)
  # independent reference: eigen decomposition of A
  e <- eigen(A)
  ref <- t(vapply(times, function(t) {
    Re(e$vectors %*% diag(exp(e$values * t)) %*% solve(e$vectors) %*% x0)
  }, numeric(2)))
  expect_lt(max(abs(sim_signal(r, "sys") - ref)), 1e-7)
  expect_equal(colnames(sim_states(r)), c("sys[1]", "sys[2]"))
})

test_that("results do not depend on block or connection insertion order", {
  ref <- tight(sim_compile(connect_feedback(feedback_model(1:3), 1:3)), c(0, 4))
  for (ord in list(c(3, 1, 2), c(2, 3, 1))) {
    m <- connect_feedback(feedback_model(ord), rev(ord))
    p <- sim_compile(m)
    expect_equal(p$diagnostics$order, sim_compile(connect_feedback(feedback_model(1:3)))$diagnostics$order)
    r <- tight(p, c(0, 4))
    expect_identical(r$signals, ref$signals)
    expect_identical(r$states, ref$states)
  }
})

test_that("output grids are independent of internal steps", {
  plan <- sim_compile(decay_model(1))
  r1 <- tight(plan, c(0, 2), times = c(0, 2))
  r2 <- tight(plan, c(0, 2), times = seq(0, 2, 0.001))
  expect_equal(r1$time, c(0, 2))
  expect_equal(unname(sim_signal(r1, "x")[2]), exp(-2), tolerance = 1e-9)
  expect_equal(unname(sim_signal(r2, "x")[2001]), exp(-2), tolerance = 1e-9)
  expect_error(sim_run(plan, c(0, 2), times = 3), "within")
})

test_that("tighter tolerances reduce the error", {
  plan <- sim_compile(decay_model(3))
  err <- vapply(c(1e-4, 1e-7, 1e-10), function(tol) {
    r <- sim_run(plan, c(0, 2), times = seq(0, 2, 0.1), rtol = tol, atol = tol * 1e-2)
    max(abs(sim_signal(r, "x") - exp(-3 * r$time)))
  }, 0)
  expect_true(err[3] < err[2] && err[2] < err[1])
  expect_lt(err[3], 1e-8)
})

test_that("transfer functions match analytic step responses", {
  # 1 / (s^2 + 2 z w s + w^2) * w^2, underdamped second order
  w <- 2; z <- 0.3
  m <- sim_model("tf")
  sim_add_block(m, "step", "u", time = 0)
  sim_add_block(m, "transfer_function", "g", num = w^2, den = c(1, 2 * z * w, w^2))
  sim_connect(m, "u", "g")
  times <- seq(0, 6, 0.1)
  r <- tight(sim_compile(m), c(0, 6), times = times)
  wd <- w * sqrt(1 - z^2)
  exact <- 1 - exp(-z * w * times) * (cos(wd * times) + z / sqrt(1 - z^2) * sin(wd * times))
  expect_lt(max(abs(sim_signal(r, "g") - exact)), 1e-7)

  # proper (not strictly proper) transfer function has direct feedthrough
  m2 <- sim_model("lead")
  sim_add_block(m2, "constant", "u", value = 1)
  sim_add_block(m2, "transfer_function", "g", num = c(2, 1), den = c(1, 1))
  sim_connect(m2, "u", "g")
  r2 <- tight(sim_compile(m2), c(0, 3), times = c(0, 1, 3))
  # (2s + 1)/(s + 1) = 2 - 1/(s + 1): step response 2 - (1 - exp(-t))
  expect_lt(max(abs(sim_signal(r2, "g") - (1 + exp(-c(0, 1, 3))))), 1e-8)
  expect_error(sim_add_block(m2, "transfer_function", "bad", num = c(1, 0, 0), den = c(1, 1)),
               "improper")
})

test_that("general ode blocks and sources integrate correctly", {
  m <- sim_model("logistic")
  sim_add_block(m, "ode", "p", initial = 0.1,
                deriv = function(t, x, u) u * x * (1 - x))
  sim_add_block(m, "constant", "r", value = 2)
  sim_connect(m, "r", "p")
  times <- seq(0, 4, 0.5)
  r <- tight(sim_compile(m), c(0, 4), times = times)
  exact <- 1 / (1 + 9 * exp(-2 * times))
  expect_lt(max(abs(sim_signal(r, "p") - exact)), 1e-8)

  m2 <- sim_model("ramp_int")
  sim_add_block(m2, "ramp", "u", slope = 2, start = 1)
  sim_add_block(m2, "integrator", "x")
  sim_connect(m2, "u", "x")
  r2 <- tight(sim_compile(m2), c(0, 3), times = c(0, 1, 2, 3))
  expect_lt(max(abs(sim_signal(r2, "x") - pmax(0, c(0, 1, 2, 3) - 1)^2)), 1e-9)
})

test_that("sine, product and math blocks evaluate as documented", {
  m <- sim_model()
  sim_add_block(m, "sine", "s", amplitude = 2, frequency = 3, phase = 0.5, bias = 1)
  sim_add_block(m, "time", "t")
  sim_add_block(m, "product", "p", ops = "*/")
  sim_add_block(m, "math", "e", fun = "exp")
  sim_add_block(m, "saturation", "sat", lower = 0, upper = 2)
  sim_add_block(m, "expression", "ex", expr = "u[1]^2 + t")
  sim_connect(m, "s", "p", input = 1)
  sim_connect(m, "t", "p", input = 2)
  sim_connect(m, "t", "e")
  sim_connect(m, "s", "sat")
  sim_connect(m, "s", "ex")
  times <- c(0.5, 1, 2)
  # s / t is infinite at t = 0, so the run starts later
  r <- sim_run(sim_compile(m), c(0.25, 2), times = times)
  s <- 1 + 2 * sin(3 * times + 0.5)
  expect_equal(unname(sim_signal(r, "p")), s / times)
  expect_equal(unname(sim_signal(r, "e")), exp(times))
  expect_equal(unname(sim_signal(r, "sat")), pmin(pmax(s, 0), 2))
  expect_equal(unname(sim_signal(r, "ex")), s^2 + times)
})

test_that("initial state overrides are applied", {
  plan <- sim_compile(decay_model(1))
  r <- tight(plan, c(0, 1), times = 1, initial = list(x = 5))
  expect_equal(unname(sim_signal(r, "x")), 5 * exp(-1), tolerance = 1e-9)
  expect_error(sim_run(plan, c(0, 1), initial = list(k = 1)), "no state")
})

test_that("fixed-step methods from deSolve can be used", {
  r <- sim_run(sim_compile(decay_model(1)), c(0, 1), times = seq(0, 1, 0.01), method = "rk4")
  expect_lt(abs(tail(sim_signal(r, "x"), 1) - exp(-1)), 1e-8)
})
