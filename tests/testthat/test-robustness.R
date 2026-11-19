test_that("portable expressions cannot reach arbitrary functions", {
  expect_error(simcompR:::parse_safe("(stop)('x')", "e"), "computed function")
  expect_error(simcompR:::parse_safe("{ sum <- c; 1 }", "e"), "must not assign")
  expect_error(simcompR:::parse_safe("system('ls')", "e"), "not allowed")
  expect_error(simcompR:::eval_safe(quote(get("system")), list()))
  # a model file with a bypass in a state machine action is refused on load
  m <- sim_model()
  sim_add_block(m, "constant", "c")
  sim_add_block(m, "state_machine", "sm", states = "a", data = c(n = 0),
                entry = list(a = "n <- 1"))
  sim_connect(m, "c", "sm")
  x <- sim_model_to_list(m)
  x$blocks[[2]]$params$entry$`$list`[[1]]$`$chr`[[1]] <- "(stop)('bypass')"
  expect_error(sim_model_from_list(x), "computed function")
})

test_that("state machine entry runs at the first hit with real inputs and time", {
  m <- sim_model()
  sim_add_block(m, "constant", "c", value = 7)
  sim_add_block(m, "state_machine", "sm", states = "a", data = c(n = 0),
                entry = list(a = "n <- u[1] + t"))
  sim_connect(m, "c", "sm")
  r <- sim_run(sim_compile(m), c(5, 6))
  expect_equal(unname(tail(sim_signal(r, "sm.data"), 1)), 12)
})

test_that("algebraic loops through sampled blocks hold the sampled output", {
  m <- sim_model()
  sim_add_block(m, "constant", "u", value = 1)
  sim_add_block(m, "subtract", "s")
  sim_add_block(m, "zero_order_hold", "h", period = 1)
  sim_connect(m, "u", "s", input = 1)
  sim_connect(m, "h", "s", input = 2)
  sim_connect(m, "s", "h")
  r <- sim_run(sim_compile(m, algebraic_loops = "solve"), c(0, 2), times = c(0, 0.5, 1.5),
               include_ticks = FALSE)
  expect_equal(unname(sim_signal(r, "s")), rep(0.5, 3), tolerance = 1e-9)
})

test_that("mode initialization does not depend on block names", {
  for (nm in c("a", "zz")) {
    m <- sim_model()
    sim_add_block(m, "constant", "m", value = 1)
    sim_add_block(m, "relay", nm)
    sim_connect(m, "m", nm)
    expect_equal(unname(sim_signal(sim_run(sim_compile(m), c(0, 1), times = 1), nm)), 1)
  }
  # a strict switch at equality takes the third input
  m <- sim_model()
  for (v in c(10, 0, 20)) sim_add_block(m, "constant", paste0("c", v), value = v)
  sim_add_block(m, "switch", "sw", threshold = 0, criterion = ">")
  sim_connect(m, "c10", "sw", input = 1)
  sim_connect(m, "c0", "sw", input = 2)
  sim_connect(m, "c20", "sw", input = 3)
  expect_equal(unname(sim_signal(sim_run(sim_compile(m), c(0, 1), times = 1), "sw")), 20)
})

test_that("crossings caused by jumps or reaching a boundary are not lost", {
  counter <- function(m, signal, level = 0, direction = "either") {
    sim_add_block(m, "integrator", "count")
    sim_add_block(m, "constant", "zero", value = 0)
    sim_connect(m, "zero", "count")
    sim_add_event(m, "e", sim_zero_crossing(signal, level = level, direction = direction),
                  reset = list(count = "x + 1"))
    m
  }
  m <- sim_model()
  sim_add_block(m, "step", "u", time = 1, before = -1, after = 1)
  r <- sim_run(sim_compile(counter(m, "u", direction = "rising")), c(0, 2))
  expect_equal(unname(tail(sim_signal(r, "count"), 1)), 1)
  expect_equal(sim_events(r)$time, 1)

  m <- sim_model()
  sim_add_block(m, "time", "t")
  sim_add_block(m, "zero_order_hold", "h", period = 1)
  sim_connect(m, "t", "h")
  r <- sim_run(sim_compile(counter(m, "t", level = 1)), c(0, 2))
  expect_equal(unname(tail(sim_signal(r, "count"), 1)), 1)
})

test_that("simultaneous events see the effects of earlier handlers", {
  m <- sim_model()
  sim_add_block(m, "time", "t")
  sim_add_block(m, "constant", "z", value = 0)
  sim_add_block(m, "integrator", "a", initial = 1)
  sim_add_block(m, "integrator", "b")
  sim_connect(m, "z", "a")
  sim_connect(m, "z", "b")
  sim_add_event(m, "e1", sim_zero_crossing("t", level = 0.5), reset = list(a = "2"), priority = 0)
  sim_add_event(m, "e2", sim_zero_crossing("t", level = 0.5), reset = list(b = "s[['a']]"),
                priority = 1)
  expect_equal(unname(tail(sim_signal(sim_run(sim_compile(m), c(0, 1)), "b"), 1)), 2)
})

test_that("scheduled terminal events record the terminal state", {
  m <- sim_model()
  sim_add_block(m, "constant", "one")
  sim_add_block(m, "integrator", "x")
  sim_connect(m, "one", "x")
  sim_add_event(m, "end", sim_at(1), reset = list(x = "10"), terminal = TRUE)
  r <- sim_run(sim_compile(m), c(0, 2), times = c(0, 2))
  expect_equal(max(r$time), 1)
  expect_equal(unname(tail(sim_signal(r, "x"), 1)), 10)
  expect_equal(r$kind[r$time == 1], c("pre-event", "post-event"))
})

test_that("clocks far from zero keep every sample hit", {
  m <- sim_model()
  sim_add_block(m, "time", "t")
  sim_add_block(m, "zero_order_hold", "h", period = 0.01)
  sim_connect(m, "t", "h")
  r <- sim_run(sim_compile(m), c(1e6, 1e6 + 1), times = 1e6 + 1, include_ticks = FALSE)
  expect_equal(unname(sim_signal(r, "h")), 1e6 + 1, tolerance = 1e-12)
})

test_that("events can watch subsystem outputs", {
  s <- sim_subsystem("pass")
  sim_add_block(s, "inport", "i", index = 1)
  sim_add_block(s, "outport", "o", index = 1)
  sim_connect(s, "i", "o")
  m <- sim_model()
  sim_add_block(m, "time", "t")
  sim_add_block(m, "subsystem", "ss", system = s)
  sim_connect(m, "t", "ss")
  sim_add_event(m, "end", sim_zero_crossing("ss", level = 1), terminal = TRUE)
  expect_equal(max(sim_run(sim_compile(m), c(0, 2))$time), 1, tolerance = 1e-7)
})

test_that("discrete initial overrides also set held outputs", {
  m <- sim_model()
  sim_add_block(m, "constant", "z0", value = 0)
  sim_add_block(m, "unit_delay", "z", period = 1, offset = 0.5)
  sim_connect(m, "z0", "z")
  r <- sim_run(sim_compile(m), c(0, 0.25), initial = list(z = 5))
  expect_true(all(sim_signal(r, "z") == 5))
})

test_that("copies do not share nested subsystems", {
  s <- sim_subsystem("g")
  sim_add_block(s, "inport", "i", index = 1)
  sim_add_block(s, "gain", "k", gain = 2)
  sim_add_block(s, "outport", "o", index = 1)
  sim_connect(s, "i", "k")
  sim_connect(s, "k", "o")
  m <- sim_model()
  sim_add_block(m, "constant", "c")
  sim_add_block(m, "subsystem", "g", system = s)
  sim_connect(m, "c", "g")
  m2 <- sim_copy(m)
  sim_set_block_params(m2$blocks$g$params$system, "k", gain = 9)
  out <- function(model) unname(sim_signal(sim_run(sim_compile(model), c(0, 1), times = 1), "g"))
  expect_equal(out(m), 2)
  expect_equal(out(m2), 9)
})

test_that("invalid outputs, sizes and names are reported, not hidden", {
  m <- sim_model()
  sim_add_block(m, "constant", "c", value = -1)
  sim_add_block(m, "math", "q", fun = "sqrt")
  sim_connect(m, "c", "q")
  expect_error(suppressWarnings(sim_run(sim_compile(m), c(0, 1))), "NaN")

  m <- sim_model()
  sim_add_block(m, "constant", "one", value = 1)
  sim_add_block(m, "sum", "s")
  sim_add_block(m, "expression", "y", expr = "c(0, 0)")
  sim_connect(m, "one", "s", input = 1)
  sim_connect(m, "y", "s", input = 2)
  sim_connect(m, "s", "y")
  expect_error(sim_compile(m, algebraic_loops = "solve"), "size")

  m <- sim_model()
  sim_add_block(m, "constant", "v", value = c(1, 2))
  sim_add_block(m, "demux", "d")
  sim_connect(m, "v", "d")
  sim_add_block(m, "constant", "d.out1", value = 99)
  expect_error(sim_compile(m), "used twice")
})

test_that("static transfer functions are allowed", {
  m <- sim_model()
  sim_add_block(m, "constant", "c", value = 3)
  sim_add_block(m, "transfer_function", "g", num = 2, den = 1)
  sim_connect(m, "c", "g")
  expect_equal(unname(sim_signal(sim_run(sim_compile(m), c(0, 1), times = 1), "g")), 6)
})

test_that("a failing run does not abort a sweep", {
  m <- sim_model(params = list(k = 1))
  sim_add_block(m, "constant", "c", value = sim_param("k"))
  sim_add_block(m, "integrator", "x")
  sim_connect(m, "c", "x")
  s <- sim_sweep(m, data.frame(k = c(1, NA, 2)), tspan = c(0, 1),
                 summary = function(r) c(end = tail(sim_signal(r, "x"), 1)))
  expect_equal(s$status, c("success", "failed", "success"))
  expect_equal(s$summary$end, c(1, NA, 2))
})

test_that("noise does not depend on the caller's generator settings", {
  m <- sim_model()
  sim_add_block(m, "white_noise", "w", seed = 1, period = 0.1)
  plan <- sim_compile(m)
  w1 <- sim_signal(sim_run(plan, c(0, 1)), "w")
  old <- RNGkind()
  on.exit(RNGkind(old[1], old[2], old[3]))
  RNGkind(normal.kind = "Box-Muller")
  w2 <- sim_signal(sim_run(plan, c(0, 1)), "w")
  expect_identical(w1, w2)
  expect_equal(RNGkind()[2], "Box-Muller")
})

test_that("non-finite states end the run with a failure", {
  m <- sim_model()
  sim_add_block(m, "integrator", "x", initial = 1)
  sim_add_block(m, "expression", "f", expr = "log(u[1] - 2)")
  sim_connect(m, "x", "f")
  sim_connect(m, "f", "x")
  expect_error(suppressWarnings(sim_run(sim_compile(m), c(0, 1))), "NaN")
})

test_that("a failure at the first row still returns a partial result", {
  m <- sim_model()
  sim_add_block(m, "constant", "c")
  sim_add_block(m, "fcn", "f", fun = function(u, t) if (t > 0.5) stop("boom") else u)
  sim_connect(m, "c", "f")
  r <- suppressWarnings(sim_run(sim_compile(m), c(1, 2), on_failure = "warn"))
  expect_false(r$success)
  expect_match(r$message, "boom")
})

test_that("crossings are found whatever the output grid", {
  m <- sim_model()
  sim_add_block(m, "sine", "s", frequency = 2 * pi)
  sim_add_block(m, "constant", "zero", value = 0)
  sim_add_block(m, "discrete_integrator", "n", period = 100)
  sim_connect(m, "zero", "n")
  sim_add_event(m, "up", sim_zero_crossing("s", "rising"), reset = list(n = "x + 1"))
  r <- sim_run(sim_compile(m), c(0.1, 10.1), times = c(0.1, 10.1))
  expect_equal(nrow(sim_events(r)), 10)
})

test_that("a breakpoint on a sample hit is seen at that hit", {
  m <- sim_model()
  sim_add_block(m, "step", "u", time = 0.9)
  sim_add_block(m, "zero_order_hold", "h", period = 0.3)
  sim_connect(m, "u", "h")
  r <- sim_run(sim_compile(m), c(0, 1.2), times = c(0.9, 1.2), include_ticks = FALSE)
  expect_equal(unname(sim_signal(r, "h")), c(1, 1))
})

test_that("a lookup table needs two points", {
  m <- sim_model()
  sim_add_block(m, "constant", "c")
  expect_error(sim_add_block(m, "lookup", "l", x = 1, y = 1), "two points")
})

test_that("multi-port dynamic blocks plot through their port signals", {
  m <- sim_model()
  sim_add_block(m, "constant", "c")
  sim_add_block(m, "state_machine", "sm", states = c("a", "b"), data = c(n = 0),
                transitions = list(sim_transition("a", "b", "t > 0.5")))
  sim_connect(m, "c", "sm")
  r <- sim_run(sim_compile(m), c(0, 1))
  expect_true(all(simcompR:::default_plot_signals(r) %in% c("sm.state", "sm.data")))
  pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_invisible(plot(r))
})

test_that("a switch inside a solved loop starts in the right mode", {
  m <- sim_model()
  sim_add_block(m, "constant", "u", value = 1)
  sim_add_block(m, "constant", "ctl", value = -1)
  sim_add_block(m, "constant", "zero", value = 0)
  sim_add_block(m, "sum", "s", signs = "+-")
  sim_add_block(m, "switch", "sw")
  sim_connect(m, "u", "s", input = 1)
  sim_connect(m, "sw", "s", input = 2)
  sim_connect(m, "s", "sw", input = 1)
  sim_connect(m, "ctl", "sw", input = 2)
  sim_connect(m, "zero", "sw", input = 3)
  sim_add_block(m, "integrator", "x")
  sim_connect(m, "s", "x")
  r <- sim_run(sim_compile(m, algebraic_loops = "solve"), c(0, 1), times = 1)
  expect_equal(unname(sim_signal(r, "s")), 1)
  expect_equal(unname(sim_signal(r, "sw")), 0)
})

test_that("direction filters hold for jumps of discontinuous signals", {
  m <- sim_model()
  sim_add_block(m, "time", "t")
  sim_add_block(m, "compare", "c", threshold = 0.5, op = ">=")
  sim_connect(m, "t", "c")
  sim_add_block(m, "integrator", "n")
  sim_add_block(m, "constant", "z", value = 0)
  sim_connect(m, "z", "n")
  sim_add_event(m, "dn", sim_zero_crossing("c", "falling", level = 0.5),
                reset = list(n = "x + 1"))
  r <- sim_run(sim_compile(m), c(0, 1))
  expect_equal(nrow(sim_events(r)), 0)
})

test_that("smooth problems far from t = 0 are integrated at the true time", {
  m <- sim_model()
  sim_add_block(m, "ramp", "u", start = 1e9)
  sim_add_block(m, "integrator", "x")
  sim_connect(m, "u", "x")
  r <- sim_run(sim_compile(m), c(1e9, 1e9 + 1), times = 1e9 + 1, rtol = 1e-10, atol = 1e-12)
  expect_equal(unname(tail(sim_signal(r, "x"), 1)), 0.5, tolerance = 1e-5)
})

test_that("solver warnings are kept in the diagnostics", {
  m <- sim_model()
  sim_add_block(m, "integrator", "x", initial = 1)
  sim_add_block(m, "gain", "k", gain = -1)
  sim_connect(m, "x", "k")
  sim_connect(m, "k", "x")
  r <- suppressWarnings(sim_run(sim_compile(m), c(0, 100), times = c(0, 100),
                                maxsteps = 3, on_failure = "warn"))
  expect_equal(r$status, "failed")
  expect_type(sim_diagnostics(r)$warnings, "character")
  expect_gt(length(sim_diagnostics(r)$warnings), 0)
})

test_that("infinite signals and states are rejected", {
  m <- sim_model()
  sim_add_block(m, "constant", "one")
  sim_add_block(m, "constant", "zero", value = 0)
  sim_add_block(m, "divide", "q")
  sim_connect(m, "one", "q", input = 1)
  sim_connect(m, "zero", "q", input = 2)
  expect_error(sim_run(sim_compile(m), c(0, 1)), "signals must be finite")

  m <- sim_model()
  sim_add_block(m, "constant", "c", value = 0)
  sim_add_block(m, "integrator", "x")
  sim_connect(m, "c", "x")
  sim_add_event(m, "bad", sim_at(1), reset = list(x = "1 / 0"))
  expect_error(sim_run(sim_compile(m), c(0, 1)), "states must be finite")
  sim_remove_event(m, "bad")
  sim_add_event(m, "bad", sim_at(0.5), reset = list(x = "1 / 0"), terminal = TRUE)
  expect_error(sim_run(sim_compile(m), c(0, 1)), "states must be finite")

  m <- sim_model()
  sim_add_block(m, "constant", "big", value = 1e308)
  sim_add_block(m, "discrete_integrator", "d", gain = 10, period = 1)
  sim_connect(m, "big", "d")
  expect_error(sim_run(sim_compile(m), c(0, 2)), "states must be finite")
})

test_that("solver statistics report Jacobians for ODEPACK solvers only", {
  m <- sim_model()
  sim_add_block(m, "gain", "k", gain = -1000)
  sim_add_block(m, "integrator", "x", initial = 1)
  sim_connect(m, "x", "k")
  sim_connect(m, "k", "x")
  p <- sim_compile(m)
  r <- sim_run(p, c(0, 1), times = c(0, 1))
  # the same solve called directly, with the arguments sim_run() passes
  direct <- deSolve::ode(1, c(0, 1), function(t, y, parms) list(-1000 * y), NULL,
                         method = "lsoda", rtol = 1e-6, atol = 1e-8, maxsteps = 5000,
                         tcrit = 1)
  ist <- attr(direct, "istate")
  expect_equal(r$diagnostics$stats[["steps"]], ist[2])
  expect_equal(r$diagnostics$stats[["jacobians"]], ist[14])
  expect_gt(r$diagnostics$stats[["jacobians"]], 0)
  expect_true(is.na(sim_run(p, c(0, 1), times = c(0, 1), method = "rk4")$diagnostics$stats[["jacobians"]]))
})
