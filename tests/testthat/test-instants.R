# Event semantics at a single instant and the timeline of boundaries.

add <- sim_add_block
wire <- sim_connect
run1 <- function(m, span = c(0, 1), times = span, ...) {
  sim_run(sim_compile(m), span, times = times, ...)
}
last <- function(r, name) unname(tail(sim_signal(r, name), 1))
flat_state <- function(m, name = "x", initial = 0) {
  add(m, "constant", paste0(name, "_rate"), value = 0)
  add(m, "integrator", name, initial = initial)
  wire(m, paste0(name, "_rate"), name)
}

test_that("a jump made by a block handler is seen by model events", {
  m <- sim_model()
  add(m, "time", "t")
  add(m, "relay", "q", on_point = 0.5)
  wire(m, "t", "q")
  sim_add_event(m, "halt", sim_zero_crossing("q", "rising", level = 0.5), terminal = TRUE)
  r <- run1(m)
  expect_identical(r$status, "terminated")
  expect_equal(sim_events(r)$event, c("q", "halt"))
  expect_equal(tail(r$time, 1), 0.5)
})

test_that("jumps onto, off and across a level all fire", {
  for (v in list(c(-1, 0), c(0, 1), c(-1, 1))) {
    m <- sim_model()
    add(m, "step", "u", time = 0.5, before = v[1], after = v[2])
    sim_add_event(m, "halt", sim_zero_crossing("u", "rising"), terminal = TRUE)
    r <- run1(m)
    expect_identical(r$status, "terminated")
    expect_equal(tail(r$time, 1), 0.5)
  }
})

test_that("a terminal event does not skip a departure due at the same instant", {
  m <- sim_model()
  flat_state(m)
  add(m, "time", "t")
  sim_add_event(m, "finish", sim_at(0), terminal = TRUE)
  sim_add_event(m, "leave", sim_zero_crossing("t", "rising"), reset = list(x = "1"))
  r <- run1(m)
  expect_equal(sim_events(r)$event, c("finish", "leave"))
  expect_equal(last(r, "x"), 1)
})

test_that("the probe sees a departure that reverses later and a quartic plateau", {
  m <- sim_model()
  add(m, "signal_function", "u", fun = function(t) t * (1 - 2000 * t))
  sim_add_event(m, "fall", sim_zero_crossing("u", "falling"), terminal = TRUE)
  expect_equal(sim_events(run1(m, hmax = 1e-5))$time, 5e-4, tolerance = 1e-8)

  m <- sim_model()
  add(m, "signal_function", "u", fun = function(t) pmax(0, t - 0.5))
  sim_add_event(m, "rise", sim_zero_crossing("u", "rising"), terminal = TRUE)
  r <- run1(m, times = c(0, 0.5, 1), hmax = 0.01)
  expect_identical(r$status, "terminated")
  expect_equal(sim_events(r)$time, 0.5, tolerance = 1e-8)
})

test_that("the event budget wins over a terminal event and counts handlers", {
  m <- sim_model()
  add(m, "time", "t")
  for (q in c("a", "b", "c")) {
    add(m, "relay", q, on_point = 0.5)
    wire(m, "t", q)
  }
  r <- suppressWarnings(run1(m, max_events = 1))
  expect_identical(r$status, "event_limit")
  expect_equal(nrow(sim_events(r)), 1)

  m <- sim_model()
  flat_state(m)
  sim_add_event(m, "a", sim_at(0.5), terminal = TRUE)
  sim_add_event(m, "b", sim_at(0.5), reset = list(x = "7"))
  r <- suppressWarnings(run1(m, max_events = 1))
  expect_identical(r$status, "event_limit")
  expect_equal(last(r, "x"), 0)
})

test_that("each reset that crosses a level is seen, even if a later one undoes it", {
  m <- sim_model()
  flat_state(m, "x", -1)
  flat_state(m, "n")
  sim_add_event(m, "a", sim_at(0.5), reset = list(x = "1"))
  sim_add_event(m, "b", sim_at(0.5), reset = list(x = "-1"))
  sim_add_event(m, "up", sim_zero_crossing("x", "rising"), reset = list(n = "x + 1"))
  r <- run1(m)
  expect_equal(sim_events(r)$event, c("a", "b", "up"))
  expect_equal(last(r, "n"), 1)
})

test_that("a stop block logs once and records the row before it", {
  m <- sim_model()
  add(m, "step", "u", time = 0.5)
  add(m, "stop", "end", threshold = 0.5)
  wire(m, "u", "end")
  r <- run1(m)
  expect_equal(nrow(sim_events(r)), 1)
  expect_equal(r$kind, c("sample", "pre-event", "post-event"))
})

test_that("an output time at a root is recorded once, as the event rows", {
  m <- sim_model()
  add(m, "time", "t")
  sim_add_event(m, "halt", sim_zero_crossing("t", level = 0.5 + 2 * .Machine$double.eps),
                terminal = TRUE)
  r <- run1(m, c(0, 10), times = c(0, 0.5, 10), hmax = 0.01)
  expect_equal(r$kind, c("sample", "pre-event", "post-event"))
})

test_that("switch criteria hold on a plateau and after a departure", {
  for (op in c("!=", ">=")) {
    m <- sim_model()
    add(m, "constant", "one")
    add(m, "constant", "zero", value = 0)
    if (op == "!=") add(m, "time", "u") else add(m, "step", "u", time = 0.5, before = -1, after = 0)
    add(m, "switch", "sw", criterion = op)
    for (i in 1:3) wire(m, c("one", "u", "zero")[i], "sw", input = i)
    add(m, "integrator", "x")
    wire(m, "sw", "x")
    expect_equal(last(run1(m), "x"), if (op == "!=") 1 else 0.5, tolerance = 1e-8)
  }
  m <- sim_model()
  add(m, "step", "ctl", time = 1, before = 1e12 - 1, after = 1e12 + 1)
  add(m, "constant", "hi", value = 2)
  add(m, "constant", "lo", value = 1)
  add(m, "switch", "sw", threshold = 1e12)
  wire(m, "hi", "sw", input = 1)
  wire(m, "ctl", "sw", input = 2)
  wire(m, "lo", "sw", input = 3)
  expect_equal(unname(sim_signal(run1(m, c(0, 2)), "sw")), c(1, 2))
})

test_that("times far from zero keep their own boundaries and left limits", {
  T0 <- 2^40
  u <- 2^-12
  m <- sim_model()
  add(m, "piecewise", "v", times = T0 + 1 + c(0, 14, 20, 33) * u, values = c(0, 0, 1, 1))
  add(m, "zero_order_hold", "z", period = 2, offset = 1 + 20 * u)
  wire(m, "v", "z")
  r <- run1(m, T0 + c(0, 2))
  expect_equal(last(r, "z"), 1)

  # a ramp is integrated exactly by rk4; no block sees a left limit
  b <- 1e9
  m <- sim_model()
  add(m, "ramp", "u", start = b, slope = 2^24)
  add(m, "integrator", "x")
  wire(m, "u", "x")
  expect_equal(last(run1(m, b + c(0, 2^-10), method = "rk4"), "x"), 8, tolerance = 1e-12)

  # a clock hit exactly at a step reads the stepped value
  m <- sim_model()
  add(m, "step", "u", time = b)
  add(m, "zero_order_hold", "h", period = 1)
  wire(m, "u", "h")
  suppressWarnings(sim_add_event(m, "anchor", sim_at(b + 2^-23 * c(-30, -1, 20))))
  expect_equal(last(run1(m, b + c(-1, 1) / 64), "h"), 1)
})

test_that("short spans are integrated to their end", {
  m <- sim_model()
  add(m, "constant", "u", value = 1e15)
  add(m, "integrator", "x")
  wire(m, "u", "x")
  expect_equal(last(run1(m, c(0, 1e-15)), "x"), 1)

  tt <- 1e9 + c(0, 2^-20)
  m <- sim_model()
  add(m, "ode", "x", initial = 0, n_inputs = 0, deriv = function(t, x, u) 2^20)
  r <- run1(m, tt, rtol = 1e-11, atol = 1e-13)
  expect_identical(r$status, "success")
  expect_equal(last(r, "x"), 1, tolerance = 1e-9)
})

test_that("a root snapped to a boundary carries the states to the boundary", {
  tt <- 1e9 + c(0, 1 / 64)
  m <- sim_model()
  add(m, "ode", "x", initial = 0, n_inputs = 0, deriv = function(t, x, u) 1)
  sim_add_event(m, "near", sim_zero_crossing("x", "rising", level = 1 / 64 - 2^-19),
                reset = list(x = "x"))
  r <- run1(m, tt, rtol = 1e-11, atol = 1e-13)
  expect_equal(last(r, "x"), 1 / 64, tolerance = 1e-12)
})

test_that("a sample reads the switch mode left by a coincident event", {
  m <- sim_model()
  add(m, "ode", "x", initial = -1, n_inputs = 0, deriv = function(t, x, u) 0)
  add(m, "constant", "hi", value = 1)
  add(m, "constant", "lo", value = -1)
  add(m, "switch", "sw")
  add(m, "zero_order_hold", "h", period = 1)
  wire(m, "hi", "sw", input = 1)
  wire(m, "x", "sw", input = 2)
  wire(m, "lo", "sw", input = 3)
  wire(m, "sw", "h")
  sim_add_event(m, "flip", sim_at(1), reset = list(x = "1"))
  r <- run1(m, c(0, 2))
  post <- r$signals[r$time == 1 & r$kind == "post-event", c("sw", "h")]
  expect_equal(unname(post), c(1, 1))
})

test_that("unresolvable pulse edges are reported promptly", {
  m <- sim_model()
  add(m, "pulse", "p", period = 1)
  expect_error(sim_run(sim_compile(m), c(2^53, 2^53 + 64), times = 2^53),
               "resolution of double precision")
})

test_that("outputs that jump at a kink switch downstream blocks at the instant", {
  gate <- function(op) {
    m <- sim_model()
    add(m, "time", "t")
    add(m, "compare", "c", op = op, threshold = 3)
    add(m, "constant", "one")
    add(m, "constant", "zero", value = 0)
    add(m, "switch", "sw", threshold = 0.5)
    add(m, "integrator", "x")
    wire(m, "t", "c")
    wire(m, "one", "sw", input = 1)
    wire(m, "c", "sw", input = 2)
    wire(m, "zero", "sw", input = 3)
    wire(m, "sw", "x")
    last(run1(m, c(0, 5)), "x")
  }
  expect_equal(gate(">"), 2, tolerance = 1e-8)
  expect_equal(gate(">="), 2, tolerance = 1e-8)

  # starts on the kink: sin(0) = 0
  m <- sim_model()
  add(m, "sine", "s")
  add(m, "compare", "c")
  add(m, "constant", "hi", value = 1)
  add(m, "constant", "lo", value = -1)
  add(m, "switch", "sw", threshold = 0.5)
  add(m, "integrator", "x")
  wire(m, "s", "c")
  wire(m, "hi", "sw", input = 1)
  wire(m, "c", "sw", input = 2)
  wire(m, "lo", "sw", input = 3)
  wire(m, "sw", "x")
  expect_equal(last(run1(m, c(0, 9)), "x"), 9 - 2 * pi, tolerance = 1e-6)

  m <- sim_model()
  add(m, "time", "t")
  add(m, "compare", "c", threshold = 3)
  add(m, "stop", "st", threshold = 0.5)
  wire(m, "t", "c")
  wire(m, "c", "st")
  expect_equal(tail(run1(m, c(0, 5))$time, 1), 3, tolerance = 1e-8)

  m <- sim_model()
  add(m, "time", "t")
  add(m, "math", "f", fun = "ceiling")
  wire(m, "t", "f")
  add(m, "constant", "one")
  add(m, "constant", "zero", value = 0)
  add(m, "switch", "sw", threshold = 0.5)
  add(m, "integrator", "x")
  wire(m, "one", "sw", input = 1)
  wire(m, "f", "sw", input = 2)
  wire(m, "zero", "sw", input = 3)
  wire(m, "sw", "x")
  expect_equal(last(run1(m, c(0, 5)), "x"), 5, tolerance = 1e-8)
})

test_that("events see the output of a switch or relay flip", {
  m <- sim_model()
  add(m, "integrator", "x")
  add(m, "relay", "q", on_point = 0.5, off_point = -0.5, on_value = 1, off_value = -1,
      initial = "on")
  add(m, "gain", "neg", gain = -1)
  wire(m, "x", "q")
  wire(m, "q", "neg")
  wire(m, "neg", "x")
  flat_state(m, "n")
  sim_add_event(m, "off", sim_zero_crossing("q", "falling"), reset = list(n = "x + 1"))
  sim_add_event(m, "on", sim_zero_crossing("q", "rising"), reset = list(n = "x + 100"))
  expect_equal(last(run1(m, c(0, 5)), "n"), 203)
})

test_that("approach and departure inside the band are told apart", {
  m <- sim_model()
  add(m, "integrator", "x", initial = 1)
  add(m, "gain", "k", gain = -1)
  wire(m, "x", "k")
  wire(m, "k", "x")
  add(m, "zero_order_hold", "z", period = 0.01)
  wire(m, "x", "z")
  suppressWarnings(sim_add_event(m, "empty", sim_zero_crossing("x", "falling")))
  r <- run1(m, c(0, 60))
  expect_identical(r$status, "success")
  expect_lte(nrow(sim_events(r)), 1)

  m <- sim_model()
  add(m, "constant", "g", value = -1)
  add(m, "integrator", "v")
  add(m, "integrator", "h")
  wire(m, "g", "v")
  wire(m, "v", "h")
  suppressWarnings(sim_add_event(m, "below", sim_zero_crossing("h", "falling")))
  expect_equal(sim_events(run1(m, c(0, 2)))$time, 0)

  m <- sim_model()
  add(m, "constant", "one")
  add(m, "discrete_integrator", "c", period = 1)
  wire(m, "one", "c")
  suppressWarnings(sim_add_event(m, "e", sim_zero_crossing("c", "rising")))
  expect_equal(sim_events(run1(m, c(0, 3)))$time, 1)
})

test_that("fixed-step methods do not take the output grid as their step", {
  m <- sim_model()
  add(m, "integrator", "x", initial = 1)
  add(m, "gain", "k", gain = -1)
  wire(m, "x", "k")
  wire(m, "k", "x")
  expect_equal(last(run1(m, c(0, 5), method = "rk4"), "x"), exp(-5), tolerance = 1e-6)
})

test_that("chattering at a kink and block failures end with what was reached", {
  m <- sim_model()
  add(m, "integrator", "x", initial = 1)
  add(m, "math", "s", fun = "sign")
  add(m, "gain", "k", gain = -1)
  wire(m, "x", "s")
  wire(m, "s", "k")
  wire(m, "k", "x")
  r <- suppressWarnings(run1(m, c(0, 1.2), times = c(0, 1, 1.2), on_failure = "warn"))
  expect_identical(r$status, "failed")

  m <- sim_model()
  add(m, "constant", "c", value = -1)
  add(m, "integrator", "x", initial = 1)
  wire(m, "c", "x")
  add(m, "math", "lg", fun = "log")
  wire(m, "x", "lg")
  r <- suppressWarnings(run1(m, c(0, 2), times = seq(0, 2, 0.25), on_failure = "warn"))
  expect_true(all(c(0.25, 0.5, 0.75) %in% sim_times(r)))
  expect_match(r$message, "at t = 1")
})

test_that("run arguments are validated", {
  m <- sim_model()
  flat_state(m)
  p <- sim_compile(m)
  expect_error(sim_run(p, c(0, 1), zc_tol = "a"), "zc_tol")
  expect_error(sim_run(p, c(0, 1), maxsteps = 0), "maxsteps")
  expect_error(sim_run(p, c(0, 1), initial = list(x = 1, x = 2)), "more than once")
  expect_error(sim_monte_carlo(p, 2, function(i) list(), seed = 1.5, tspan = c(0, 1)), "seed")
  expect_error(sim_register_block_type("two_abbrev", params = list(n = 1, na = 2),
                                       build = function(p) list(outputs = 1L,
                                                                output = function(t, x, d, u) 1)),
               "abbreviate")
})
