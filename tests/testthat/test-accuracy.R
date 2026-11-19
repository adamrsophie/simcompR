# Regression tests with independently derived answers for timing, forcing,
# event and clock semantics.

integrate_source <- function(add_source, tspan = c(0, 1), ...) {
  m <- sim_model()
  add_source(m)
  sim_add_block(m, "integrator", "x", initial = rep(0, 1))
  sim_connect(m, "u", "x")
  sim_run(sim_compile(m), tspan, times = tspan, ...)
}

test_that("narrow linear features of from_data are integrated with a sparse output grid", {
  for (case in list(c(0.5, 0.01), c(0.123, 0.002), c(0.9, 0.05))) {
    mid <- case[1]; hw <- case[2]
    r <- integrate_source(function(m) {
      sim_add_block(m, "from_data", "u", times = c(0, mid - hw, mid, mid + hw, 1),
                    values = c(0, 0, 1, 0, 0), method = "linear")
    })
    expect_equal(unname(tail(sim_signal(r, "x"), 1)), hw, tolerance = 1e-6)
  }
})

test_that("vector from_data tables integrate each column, holding the end values", {
  m <- sim_model()
  sim_add_block(m, "from_data", "u", times = c(0.2, 0.3, 0.4),
                values = cbind(c(0, 2, 0), c(1, 1, 3)))
  sim_add_block(m, "integrator", "x", initial = c(0, 0))
  sim_connect(m, "u", "x")
  r <- sim_run(sim_compile(m), c(0, 1), times = c(0, 1))
  # column 1: triangle of area 0.2; column 2: 1 * 0.2 + 1 * 0.1 + 2 * 0.1 + 3 * 0.6
  expect_equal(unname(sim_signal(r, "x")[2, ]), c(0.2, 0.2 + 0.1 + 0.2 + 1.8),
               tolerance = 1e-6)
})

test_that("a single-point from_data table is a constant signal", {
  m <- sim_model()
  sim_add_block(m, "from_data", "f", times = 0, values = 4)
  r <- sim_run(sim_compile(m), c(0, 1), times = c(0, 1))
  expect_equal(unname(sim_signal(r, "f")), c(4, 4))
})

test_that("a pulse that starts after the run creates no breakpoints", {
  m <- sim_model()
  sim_add_block(m, "pulse", "p", period = 0.001, delay = 1e12)
  p <- sim_compile(m)
  expect_length(p$blocks$p$def$breakpoints(0, 1), 0L)
  expect_equal(unname(sim_signal(sim_run(p, c(0, 1), times = c(0, 1)), "p")), c(0, 0))
})

test_that("sample clocks keep every hit when the time axis is shifted", {
  m <- sim_model()
  sim_add_block(m, "constant", "one")
  sim_add_block(m, "discrete_integrator", "n", gain = 1000, period = 0.001)
  sim_add_block(m, "discrete_integrator", "q", gain = 1, period = 0.01, offset = 0.005)
  sim_connect(m, "one", "n")
  sim_connect(m, "one", "q")
  p <- sim_compile(m)
  run <- function(t0) sim_run(p, t0 + c(0, 1), times = t0 + c(0, 1), include_ticks = FALSE)
  ref <- run(0)
  expect_equal(unname(sim_signal(ref, "n"))[2], 1000)
  for (t0 in c(1e6, 1e9)) {
    r <- run(t0)
    expect_equal(max(r$time) - t0, 1, info = format(t0))
    expect_equal(sim_signal(r, c("n", "q")), sim_signal(ref, c("n", "q")), info = format(t0))
    expect_equal(sim_discrete_states(r), sim_discrete_states(ref), info = format(t0))
  }
  expect_error(sim_run(p, 1e13 + c(0, 1)), "too short to be resolved")
})

test_that("short periods and continuous states stay exact far from zero", {
  m <- sim_model()
  sim_add_block(m, "constant", "one")
  sim_add_block(m, "discrete_integrator", "n", gain = 1, period = 1e-4)
  sim_add_block(m, "integrator", "x")
  sim_connect(m, "one", "n")
  sim_connect(m, "one", "x")
  p <- sim_compile(m)
  for (t0 in c(0, 1e9, 1.7e9)) {
    tf <- t0 + 0.1
    r <- sim_run(p, c(t0, tf), times = c(t0, tf), include_ticks = FALSE)
    # 1001 hits in [t0, tf], each adding one period
    expect_equal(unname(sim_discrete_states(r)[2, 1]), 0.1001, tolerance = 1e-12, info = format(t0))
    # the integrated length is the represented length of the interval
    expect_equal(unname(tail(sim_signal(r, "x"), 1)), tf - t0, tolerance = 1e-12, info = format(t0))
  }
})
