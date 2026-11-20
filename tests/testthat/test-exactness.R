# Results that must not depend on the time origin, the magnitude of signals
# or the tolerance chosen, and contracts at the edges of the interface.

add <- sim_add_block
wire <- sim_connect
run1 <- function(m, span = c(0, 1), times = span, ...) {
  sim_run(sim_compile(m), span, times = times, ...)
}
last <- function(r, name) unname(tail(sim_signal(r, name), 1))

offset_model <- function(B, fun) {
  m <- sim_model()
  add(m, "ramp", "u", initial = B)
  add(m, "math", "f", fun = fun)
  add(m, "constant", "base", value = B)
  add(m, "subtract", "y")
  add(m, "integrator", "x")
  wire(m, "u", "f")
  wire(m, "f", "y", input = 1)
  wire(m, "base", "y", input = 2)
  wire(m, "y", "x")
  m
}

test_that("jump functions are exact far from the origin", {
  exact <- c(floor = 1, ceiling = 3, round = 2)
  for (fun in names(exact)) for (B in c(0, 1e8, 1e9, -1e9)) {
    r <- run1(offset_model(B, fun), c(0, 2), times = seq(0, 2, 0.25), hmax = 0.01,
              rtol = 1e-10, atol = 1e-12)
    # the kink time can only be resolved to one unit in the last place of B
    expect_equal(last(r, "x"), exact[[fun]], tolerance = 4 * .Machine$double.eps * max(1, abs(B)) + 1e-8,
                 info = paste(fun, B))
    # the recorded output is the function itself, never a shifted input
    y <- unname(sim_signal(r, "y"))
    t <- sim_times(r)
    expect_equal(y[t == 0.25], switch(fun, floor = 0, ceiling = 1, round = 0), info = paste(fun, B))
  }
})

test_that("a comparison keeps its pattern when input and threshold shift together", {
  compare_model <- function(B) {
    m <- sim_model()
    add(m, "constant", "base", value = B - 0.5)
    add(m, "time", "t")
    add(m, "step", "drop", time = 1, after = -2)
    add(m, "sum", "u", signs = "+++")
    wire(m, "base", "u", input = 1)
    wire(m, "t", "u", input = 2)
    wire(m, "drop", "u", input = 3)
    add(m, "compare", "cmp", threshold = B)
    wire(m, "u", "cmp")
    add(m, "integrator", "x")
    wire(m, "cmp", "x")
    m
  }
  for (B in c(0, 1e12, -1e12)) {
    r <- run1(compare_model(B), c(0, 2), times = seq(0, 2, 0.25), hmax = 0.01)
    # at 1e12 one unit in the last place of the input is 1.2e-4
    expect_equal(last(r, "x"), 0.5, tolerance = 1e-4, info = format(B))
  }
})

test_that("a switch uses the run's zc_tol, including zero", {
  switch_model <- function(amp) {
    m <- sim_model()
    add(m, "step", "u", time = 1, before = -amp, after = amp)
    add(m, "constant", "one")
    add(m, "constant", "zero", value = 0)
    add(m, "switch", "sw")
    add(m, "integrator", "x")
    wire(m, "one", "sw", input = 1)
    wire(m, "u", "sw", input = 2)
    wire(m, "zero", "sw", input = 3)
    wire(m, "sw", "x")
    m
  }
  for (tol in c(0, 1e-12)) {
    r <- run1(switch_model(5e-10), c(0, 2), times = c(0, 1, 2), zc_tol = tol)
    expect_equal(last(r, "x"), 1, info = tol)
  }
  # inside the default band the control rests on the threshold: ">=" passes
  r <- run1(switch_model(5e-10), c(0, 2), times = c(0, 1, 2))
  expect_equal(last(r, "x"), 2)
})

test_that("a source jumping at the final time is seen from the left", {
  for (method in c("lsoda", "rk4", "euler")) for (tf in c(1, 2)) {
    m <- sim_model()
    add(m, "step", "u", time = 1)
    add(m, "integrator", "x")
    wire(m, "u", "x")
    r <- run1(m, c(0, tf), method = method, hmax = 0.1)
    expect_equal(last(r, "x"), tf - 1, tolerance = 1e-10, info = paste(method, tf))
    expect_equal(last(r, "u"), 1, info = paste(method, tf))
  }
  for (b in c(0, 2^20)) {
    m <- sim_model()
    add(m, "pulse", "u", period = 1, width = 2^-20, delay = b, amplitude = 2^20)
    add(m, "integrator", "x")
    wire(m, "u", "x")
    r <- run1(m, b + c(0, 1), rtol = 1e-10, atol = 1e-12)
    expect_equal(last(r, "x"), 1, tolerance = 1e-8, info = b)
  }
})

test_that("a state machine entered between hits waits its full duration", {
  m <- sim_model()
  add(m, "constant", "c")
  add(m, "state_machine", "sm", period = 1, states = c("a", "b"),
      transitions = list(sim_transition("a", "b", after = 1)))
  wire(m, "c", "sm")
  # entered at 0.4 through an initial override: 0.6 has elapsed at t = 1
  r <- run1(m, c(0.4, 3), times = c(0.4, 1, 2, 3), include_ticks = FALSE,
            initial = list(sm = c(1, 0.4)))
  expect_equal(unname(sim_signal(r, "sm")), c(1, 1, 2, 2))
  # the same entry made by an event reset at 0.4
  sim_add_event(m, "restart", sim_at(0.4), reset = list(sm = "c(1, t)"))
  r <- run1(m, c(0, 3), times = 0:3, include_ticks = FALSE)
  expect_equal(unname(sim_signal(r, "sm")[match(0:3, sim_times(r))]), c(1, 1, 2, 2))
})

test_that("a Jacobian callback sees the absolute time", {
  m <- sim_model()
  add(m, "ode", "x", initial = 1, n_inputs = 0, deriv = function(t, x, u) -1000 * t * x)
  m2 <- sim_copy(m)
  sim_add_event(m2, "kick", sim_at(5.5), reset = list(x = "2*x"))
  for (mod in list(m, m2)) {
    seen <- numeric()
    r <- run1(mod, c(5, 6), method = "lsode", jactype = "fullusr",
              jacfunc = function(t, x, p) {
                seen <<- c(seen, t)
                matrix(-1000 * t, 1, 1)
              })
    expect_identical(r$status, "success")
    expect_true(length(seen) > 0 && all(seen >= 5 & seen <= 6))
    expect_lt(abs(last(r, "x")), 1e-8)
  }
})

test_that("unsupported solver settings are refused early", {
  m <- sim_model()
  add(m, "integrator", "x", initial = 1)
  add(m, "gain", "k", gain = -1)
  wire(m, "x", "k")
  wire(m, "k", "x")
  p <- sim_compile(m)
  expect_error(sim_run(p, c(0, 1), method = "ode23"), "ode23")
  expect_error(sim_run(p, c(0, 1), events = list(data = data.frame())), "events")
  expect_equal(last(sim_run(p, c(0, 1), method = "ode45", rtol = 1e-10, atol = 1e-12), "x"),
               exp(-1), tolerance = 1e-8)
})

test_that("an infinite output of a loop cut by a held block is refused", {
  m <- sim_model()
  add(m, "zero_order_hold", "h", period = 1)
  add(m, "expression", "f", expr = "1 / (t - 0.5) + u / 2")
  wire(m, "h", "f")
  wire(m, "f", "h")
  p <- sim_compile(m, algebraic_loops = "solve")
  expect_error(sim_run(p, c(0, 1), times = c(0, 0.5, 1)), "finite")
})

test_that("every stop counts against max_events", {
  m <- sim_model()
  add(m, "constant", "one")
  for (nm in c("a", "b", "c")) {
    add(m, "stop", nm)
    wire(m, "one", nm)
  }
  r <- suppressWarnings(run1(m, max_events = 1))
  expect_identical(r$status, "event_limit")
  expect_equal(nrow(sim_events(r)), 1)
  expect_equal(nrow(sim_events(run1(m))), 3)
})

test_that("block sizes are validated before integer conversion", {
  sim_register_block_type("test_fraction_size", overwrite = TRUE,
    build = function(p) list(inputs = 1, outputs = NA, feedthrough = TRUE,
                             size_out = function(n) 2.5,
                             output = function(t, x, d, u) c(u[[1]], u[[1]])))
  sim_register_block_type("test_huge_size", overwrite = TRUE,
    build = function(p) list(outputs = 2^31, output = function(t, x, d, u) 42))
  on.exit({
    rm("test_fraction_size", envir = simcompR:::.registry)
    rm("test_huge_size", envir = simcompR:::.registry)
  })
  m <- sim_model()
  add(m, "constant", "u")
  add(m, "test_fraction_size", "bad")
  wire(m, "u", "bad")
  expect_error(sim_compile(m), "size_out")
  expect_error(add(sim_model(), "test_huge_size", "bad"), "2147483647")
})

test_that("a parameter named pi takes precedence in every form", {
  value_of <- function(v) {
    m <- sim_model(params = list(pi = 10))
    add(m, "constant", "c", value = v)
    last(run1(m), "c")
  }
  expect_equal(value_of(sim_param("pi")), 10)
  expect_equal(value_of(sim_param("pi") + 1), 11)
  expect_equal(value_of(sim_param_expr("2 * pi")), 20)
  m <- sim_model(params = list(a = 1))
  add(m, "constant", "c", value = sim_param("a") + sim_param_expr("pi"))
  expect_equal(last(run1(m), "c"), 1 + pi)
})

test_that("parameter arithmetic keeps subnormal literals", {
  for (v in c(5e-324, -5e-324, .Machine$double.xmin / 3)) {
    m <- sim_model(params = list(a = 0))
    add(m, "constant", "c", value = sim_param("a") + v)
    expect_identical(last(run1(m), "c"), v)
  }
})

test_that("a parameter given twice is refused", {
  m <- sim_model()
  expect_error(add(m, "constant", "c", value = 1, value = 2), "more than once")
  add(m, "constant", "c")
  expect_error(sim_set_block_params(m, "c", value = 1, value = 2), "more than once")
  # names inside one vector value are not parameter keys
  add(m, "constant", "v", value = c(a = 1, a = 2))
  expect_equal(unname(sim_signal(run1(m), "v")[1, ]), c(1, 2))
})

test_that("sweep tables keep the class of scalar parameters", {
  m <- sim_model(params = list(a = 1, d = as.Date("2020-01-01")))
  add(m, "constant", "c", value = sim_param("a"))
  s <- sim_sweep(m, list(list(a = 1, d = as.Date("2021-01-01")), list(a = 2)), c(0, 1), times = 1)
  expect_s3_class(s$parameters$d, "Date")
  expect_equal(s$parameters$d, as.Date(c("2021-01-01", "2020-01-01")))
})
