ball_reference <- function(h0, e, n, g = 9.81) {
  # impact times and velocities of an ideal bouncing ball
  t <- sqrt(2 * h0 / g)
  times <- t
  v <- g * t * e
  for (i in seq_len(n - 1L)) {
    t <- t + 2 * v / g
    times <- c(times, t)
    v <- v * e
  }
  times
}

test_that("bouncing ball impacts and restitution are exact", {
  plan <- sim_compile(ball_model(e = 0.8, h0 = 2))
  r <- tight(plan, c(0, 4), times = seq(0, 4, 0.01))
  ev <- sim_events(r)
  ref <- ball_reference(2, 0.8, 5)
  expect_gte(nrow(ev), 5)
  expect_equal(ev$time[1:5], ref, tolerance = 1e-7)
  expect_true(all(ev$direction == "falling"))
  # pre-event and post-event rows share the impact time
  rows <- which(abs(r$time - ev$time[1]) < 1e-12)
  expect_equal(r$kind[rows], c("pre-event", "post-event"))
  v <- sim_signal(r, "v")[rows]
  expect_equal(unname(v[2]), unname(-0.8 * v[1]), tolerance = 1e-12)
  # free fall between impacts: h(t) = h0 - g t^2 / 2 before the first impact
  pre <- r$time < ref[1] - 1e-9
  expect_lt(max(abs(sim_signal(r, "h")[pre] - (2 - 9.81 * r$time[pre]^2 / 2))), 1e-8)
  # height never goes meaningfully negative
  expect_gt(min(sim_signal(r, "h")), -1e-6)
})

test_that("a chattering ball is stopped by the event limit, not by a wrong answer", {
  plan <- sim_compile(ball_model(e = 0.5, h0 = 1))
  expect_warning(r <- sim_run(plan, c(0, 5), max_events = 10), "max_events")
  expect_equal(r$status, "event_limit")
  expect_false(r$success)
  expect_lt(max(r$time), 5)
})

test_that("crossings in the other direction are ignored", {
  m <- sim_model()
  sim_add_block(m, "sine", "s")
  sim_add_block(m, "integrator", "count")
  sim_add_block(m, "constant", "zero", value = 0)
  sim_connect(m, "zero", "count")
  sim_add_event(m, "up", trigger = sim_zero_crossing("s", direction = "rising"),
                reset = list(count = "x + 1"))
  r <- tight(sim_compile(m), c(0.1, 20))
  ev <- sim_events(r)
  expect_equal(ev$time, 2 * pi * 1:3, tolerance = 1e-7)
  expect_equal(unname(tail(sim_signal(r, "count"), 1)), 3)
})

test_that("scheduled events reset states at exact times", {
  m <- decay_model(a = 1, x0 = 1)
  sim_add_event(m, "kick", trigger = sim_at(c(1, 2)), reset = list(x = "x + 1"))
  r <- tight(sim_compile(m), c(0, 3), times = c(0, 0.5, 1, 1.5, 2, 2.5, 3))
  # piecewise exponential with jumps of +1 at t = 1 and t = 2
  x1 <- exp(-1) + 1
  x2 <- x1 * exp(-1) + 1
  ref <- c(1, exp(-0.5), x1, x1 * exp(-0.5), x2, x2 * exp(-0.5), x2 * exp(-1))
  # at an event time the pre-event row comes first, then the post-event row
  keep <- !duplicated(r$time, fromLast = TRUE)
  expect_equal(r$time[keep], c(0, 0.5, 1, 1.5, 2, 2.5, 3))
  expect_lt(max(abs(sim_signal(r, "x")[keep] - ref)), 1e-8)
  expect_equal(r$kind[r$time == 1], c("pre-event", "post-event"))
  expect_equal(sim_events(r)$type, c("scheduled", "scheduled"))
})

test_that("steps restart the integrator at the discontinuity", {
  m <- sim_model()
  sim_add_block(m, "step", "u", time = 0.7, before = 0, after = 3)
  sim_add_block(m, "integrator", "x")
  sim_connect(m, "u", "x")
  r <- sim_run(sim_compile(m), c(0, 2), times = c(0, 0.7, 1, 2))
  expect_equal(unname(sim_signal(r, "x")), c(0, 0, 0.9, 3.9), tolerance = 1e-10)
})

test_that("piecewise and pulse sources are integrated exactly", {
  m <- sim_model()
  sim_add_block(m, "piecewise", "u", times = c(0, 1, 2.5), values = c(1, -2, 0.5))
  sim_add_block(m, "integrator", "x")
  sim_connect(m, "u", "x")
  r <- sim_run(sim_compile(m), c(0, 4), times = c(1, 2.5, 4))
  expect_equal(unname(sim_signal(r, "x")), c(1, 1 - 3, -2 + 0.75), tolerance = 1e-10)

  m2 <- sim_model()
  sim_add_block(m2, "pulse", "p", amplitude = 2, period = 1, width = 0.25)
  sim_add_block(m2, "integrator", "x")
  sim_connect(m2, "p", "x")
  r2 <- sim_run(sim_compile(m2), c(0, 3), times = 3)
  expect_equal(unname(sim_signal(r2, "x")), 3 * 2 * 0.25, tolerance = 1e-10)
})

test_that("a threshold-triggered switch changes dynamics at the located time", {
  # x' = 1 while x < 1, then x' = 2 (switch on control x >= 1)
  m <- sim_model()
  sim_add_block(m, "integrator", "x")
  sim_add_block(m, "constant", "fast", value = 2)
  sim_add_block(m, "constant", "slow", value = 1)
  sim_add_block(m, "switch", "sw", threshold = 1, criterion = ">=")
  sim_connect(m, "fast", "sw", input = 1)
  sim_connect(m, "x", "sw", input = 2)
  sim_connect(m, "slow", "sw", input = 3)
  sim_connect(m, "sw", "x")
  r <- tight(sim_compile(m), c(0, 1.5), times = c(0.5, 1.5))
  ev <- sim_events(r)
  expect_equal(nrow(ev), 1L)
  expect_equal(ev$time, 1, tolerance = 1e-8)
  expect_equal(ev$direction, "rising")
  expect_equal(unname(sim_signal(r, "x")[r$kind == "sample"]), c(0.5, 2), tolerance = 1e-8)
  expect_equal(unname(sim_discrete_states(r)[nrow(r$discrete), ]), 1)
})

test_that("relay hysteresis produces a limit cycle with the predicted period", {
  # x' = 0.5 - relay(x): the relay switches on above 0.5 and off below -0.5
  m <- sim_model()
  sim_add_block(m, "integrator", "x")
  sim_add_block(m, "relay", "r", on_point = 0.5, off_point = -0.5, initial = "off")
  sim_add_block(m, "constant", "bias", value = 0.5)
  sim_add_block(m, "subtract", "rhs")
  sim_connect(m, "x", "r")
  sim_connect(m, "bias", "rhs", input = 1)
  sim_connect(m, "r", "rhs", input = 2)
  sim_connect(m, "rhs", "x")
  r <- tight(sim_compile(m), c(0, 8.5))
  ev <- sim_events(r)
  # off from x = 0 up to 0.5 (t = 1), on down to -0.5 (t = 3), and so on
  expect_equal(ev$time, c(1, 3, 5, 7), tolerance = 1e-7)
  expect_equal(ev$direction, c("rising", "falling", "rising", "falling"))
  expect_lt(max(sim_signal(r, "x")), 0.5 + 1e-6)
  expect_gt(min(sim_signal(r, "x")), -0.5 - 1e-6)
})

test_that("stop blocks and terminal events end the run", {
  m <- decay_model(a = 1, x0 = 1)
  sim_add_block(m, "constant", "lim", value = 0.5)
  sim_add_block(m, "subtract", "below")
  sim_add_block(m, "stop", "halt")
  sim_connect(m, "lim", "below", input = 1)
  sim_connect(m, "x", "below", input = 2)
  sim_connect(m, "below", "halt")
  r <- tight(sim_compile(m), c(0, 5))
  expect_equal(r$status, "terminated")
  expect_true(r$success)
  expect_equal(max(r$time), log(2), tolerance = 1e-7)

  m2 <- decay_model(a = 1, x0 = 1)
  sim_add_event(m2, "end", trigger = sim_zero_crossing("x", level = 0.25), terminal = TRUE)
  r2 <- tight(sim_compile(m2), c(0, 5))
  expect_equal(max(r2$time), log(4), tolerance = 1e-7)
})

test_that("events need a root-finding method", {
  expect_error(sim_run(sim_compile(ball_model()), c(0, 1), method = "rk4"), "zero-crossing")
})

test_that("events on models without continuous states still locate crossings", {
  m <- sim_model()
  sim_add_block(m, "time", "t")
  sim_add_block(m, "stop", "halt", threshold = 1.234)
  sim_connect(m, "t", "halt")
  r <- sim_run(sim_compile(m), c(0, 3))
  expect_equal(max(r$time), 1.234, tolerance = 1e-7)
})

test_that("event resets may use functions, signals and parameters", {
  m <- ball_model(e = 0.9, h0 = 1)
  sim_add_block(m, "gain", "vv", gain = 1)
  sim_connect(m, "v", "vv")
  m$events$impact$reset$v <- function(x, s, t) -0.9 * s$vv
  r <- tight(sim_compile(m), c(0, 1))
  t1 <- sqrt(2 / 9.81)
  rows <- which(abs(r$time - t1) < 1e-6)
  v <- sim_signal(r, "v")[rows]
  expect_equal(unname(v[2]), unname(-0.9 * v[1]), tolerance = 1e-10)
  expect_error(sim_add_event(m, "bad", sim_at(1), reset = list(v = "system('ls')")),
               "not allowed")
})

held_states <- function(init) {
  m <- sim_model()
  for (nm in names(init)) {
    sim_add_block(m, "constant", paste0(nm, "_in"), value = 0)
    sim_add_block(m, "integrator", nm, initial = init[[nm]])
    sim_connect(m, paste0(nm, "_in"), nm)
  }
  m
}

test_that("blocks and events leaving their level at the start fire there", {
  # switch control t departs 0 upward at once: output 2 for all t > 0
  m <- sim_model()
  sim_add_block(m, "time", "t")
  sim_add_block(m, "constant", "hi", value = 2)
  sim_add_block(m, "constant", "lo", value = -1)
  sim_add_block(m, "switch", "sw", criterion = ">", threshold = 0)
  sim_add_block(m, "integrator", "x")
  sim_connect(m, "hi", "sw", input = 1)
  sim_connect(m, "t", "sw", input = 2)
  sim_connect(m, "lo", "sw", input = 3)
  sim_connect(m, "sw", "x")
  r <- sim_run(sim_compile(m), c(0, 2), times = c(0, 2))
  expect_equal(unname(tail(sim_signal(r, "x"), 1)), 4, tolerance = 1e-7)

  # relay at its on point turns on as the input rises
  m <- sim_model()
  sim_add_block(m, "time", "t")
  sim_add_block(m, "relay", "rl", on_point = 0, off_point = -0.5, on_value = 3, off_value = 0)
  sim_add_block(m, "integrator", "x")
  sim_connect(m, "t", "rl")
  sim_connect(m, "rl", "x")
  r <- sim_run(sim_compile(m), c(0, 2), times = c(0, 2))
  expect_equal(unname(tail(sim_signal(r, "x"), 1)), 6, tolerance = 1e-7)

  # a stop block whose input rises from its threshold stops at once
  m <- sim_model()
  sim_add_block(m, "time", "t")
  sim_add_block(m, "stop", "halt", threshold = 0)
  sim_connect(m, "t", "halt")
  r <- sim_run(sim_compile(m), c(0, 2), times = c(0, 2))
  expect_equal(r$status, "terminated")
  expect_equal(max(sim_times(r)), 0)

  # x = t^3 / 3 has zero slope at 0 but rises: the event fires at 0
  m <- sim_model()
  sim_add_block(m, "ode", "x", n_inputs = 0, initial = 0, deriv = function(t, x, u) t^2)
  sim_add_event(m, "up", sim_zero_crossing("x", direction = "rising"), terminal = TRUE)
  sim_add_event(m, "down", sim_zero_crossing("x", direction = "falling"), terminal = TRUE)
  r <- sim_run(sim_compile(m), c(0, 2), times = c(0, 2))
  expect_equal(r$status, "terminated")
  expect_equal(sim_events(r)$event, "up")
  expect_equal(max(sim_times(r)), 0)
})

test_that("resets that cross other events cascade at the same instant", {
  m <- held_states(c(p = -2, q = -2, k = 0))
  sim_add_block(m, "time", "t")
  sim_add_event(m, "first", sim_zero_crossing("t", level = 0.3), reset = list(p = "5"))
  sim_add_event(m, "second", sim_zero_crossing("p"), reset = list(q = "7"))
  sim_add_event(m, "third", sim_zero_crossing("q"), reset = list(k = "x + 10"))
  r <- sim_run(sim_compile(m), c(0, 1), times = c(0, 1))
  ev <- sim_events(r)
  expect_equal(ev$event, c("first", "second", "third"))
  expect_equal(ev$time, rep(0.3, 3), tolerance = 1e-8)
  expect_equal(unname(tail(sim_signal(r, "k"), 1)), 10)
  # one row before the first event and one after the last
  expect_equal(sum(abs(sim_times(r) - 0.3) < 1e-8), 2L)

  m <- held_states(c(p = -2, q = -2, k = 0))
  sim_add_event(m, "first", sim_at(0.6), reset = list(p = "5"))
  sim_add_event(m, "second", sim_zero_crossing("p"), reset = list(q = "7"))
  sim_add_event(m, "third", sim_zero_crossing("q"), reset = list(k = "x + 10"))
  r <- sim_run(sim_compile(m), c(0, 1), times = c(0, 1))
  expect_equal(sim_events(r)$event, c("first", "second", "third"))
  expect_equal(unname(tail(sim_signal(r, "k"), 1)), 10)
})

test_that("events of one instant follow priority across trigger types", {
  for (term in c(FALSE, TRUE)) {
    m <- held_states(c(v = 1))
    sim_add_block(m, "time", "t")
    sim_add_event(m, "late", sim_zero_crossing("t", level = 0.4),
                  reset = list(v = "x * 10"), priority = 5, terminal = term)
    sim_add_event(m, "early", sim_at(0.4), reset = list(v = "x + 4"), priority = 1)
    r <- sim_run(sim_compile(m), c(0, 1), times = c(0, 1))
    expect_equal(sim_events(r)$event, c("early", "late"), info = term)
    expect_equal(unname(tail(sim_signal(r, "v"), 1)), 50, info = term)
  }
  # a terminal event still lets later events of the same instant run
  m <- held_states(c(v = 1))
  sim_add_event(m, "stop_first", sim_at(0.4), terminal = TRUE, priority = 1)
  sim_add_event(m, "after", sim_at(0.4), reset = list(v = "x + 1"), priority = 2)
  r <- sim_run(sim_compile(m), c(0, 1), times = c(0, 1))
  expect_equal(r$status, "terminated")
  expect_equal(unname(tail(sim_signal(r, "v"), 1)), 2)
})

test_that("cascades that never settle end at the event limit", {
  # each event flips the signal the other one watches
  m <- held_states(c(a = -1, b = -1))
  sim_add_event(m, "kick", sim_at(0.5), reset = list(a = "1"))
  sim_add_event(m, "ping", sim_zero_crossing("a"), reset = list(b = "-x"))
  sim_add_event(m, "pong", sim_zero_crossing("b"), reset = list(a = "-x"))
  r <- suppressWarnings(sim_run(sim_compile(m), c(0, 1), max_events = 20, on_failure = "warn"))
  expect_equal(r$status, "event_limit")
})

test_that("a root located on a boundary does not hide jump crossings there", {
  m <- sim_model()
  sim_add_block(m, "constant", "one")
  sim_add_block(m, "integrator", "y")
  sim_connect(m, "one", "y")
  sim_add_block(m, "step", "s", time = 0.5)
  sim_add_block(m, "integrator", "x")
  sim_connect(m, "s", "x")
  sim_add_event(m, "on_step", sim_zero_crossing("s", level = 0.5), reset = list(x = "x + 10"))
  sim_add_event(m, "on_y", sim_zero_crossing("y", level = 0.5), reset = list(x = "x + 100"))
  r <- sim_run(sim_compile(m), c(0, 1), times = c(0, 1))
  expect_setequal(sim_events(r)$event, c("on_step", "on_y"))
  expect_equal(sim_events(r)$time, c(0.5, 0.5), tolerance = 1e-8)
  expect_equal(unname(tail(sim_signal(r, "x"), 1)), 110.5, tolerance = 1e-7)
})

test_that("a slow signal far from zero that returns to its level keeps its direction", {
  # y = 100 + 0.01 (t - t^2) leaves 100 rising at 0 and returns falling at 1
  m <- sim_model()
  sim_add_block(m, "ode", "y", n_inputs = 0, initial = 100,
                deriv = function(t, x, u) 0.01 * (1 - 2 * t))
  sim_add_event(m, "up", sim_zero_crossing("y", "rising", level = 100), terminal = FALSE,
                reset = list(y = "x"))
  sim_add_event(m, "down", sim_zero_crossing("y", "falling", level = 100), terminal = TRUE)
  r <- sim_run(sim_compile(m), c(0, 1.5), times = c(0, 1.5), rtol = 1e-10, atol = 1e-10)
  ev <- sim_events(r)
  expect_equal(ev$event, c("up", "down"))
  expect_equal(ev$direction, c("rising", "falling"))
  expect_equal(ev$time[2], 1, tolerance = 1e-6)
})
