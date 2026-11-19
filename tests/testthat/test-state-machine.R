test_that("timed transitions cycle a traffic light deterministically", {
  m <- sim_model()
  sim_add_block(m, "constant", "go", value = 1)
  sim_add_block(m, "state_machine", "light", period = 1,
    states = c("red", "green", "yellow"),
    transitions = list(
      sim_transition("red", "green", after = 3),
      sim_transition("green", "yellow", after = 2),
      sim_transition("yellow", "red", after = 1)),
    outputs = list(red = 1, green = 2, yellow = 3))
  sim_connect(m, "go", "light")
  r <- sim_run(sim_compile(m), c(0, 12), times = 0:12)
  expect_equal(unname(sim_signal(r, "light.y")),
               c(1, 1, 1, 2, 2, 3, 1, 1, 1, 2, 2, 3, 1))
  expect_equal(unname(sim_signal(r, "light.state")),
               unname(sim_signal(r, "light.y")))
})

test_that("guards are tried by priority and actions update data", {
  m <- sim_model()
  sim_add_block(m, "time", "t")
  sim_add_block(m, "state_machine", "sm", period = 1,
    states = c("a", "b", "c"),
    transitions = list(
      sim_transition("a", "c", "u[1] >= 2", priority = 2),
      sim_transition("a", "b", "u[1] >= 2", priority = 1, action = "n <- n + 10"),
      sim_transition("b", "a", "TRUE")),
    entry = list(a = "n <- n + 1"),
    data = c(n = 0))
  sim_connect(m, "t", "sm")
  r <- sim_run(sim_compile(m), c(0, 5), times = 0:5)
  # t=2: a -> b (priority 1, n = 11); t=3: b -> a (entry, n = 12);
  # t=4: a -> b (n = 22); t=5: b -> a (n = 23). Initial entry gives n = 1.
  expect_equal(unname(sim_signal(r, "sm.state")), c(1, 1, 2, 1, 2, 1))
  expect_equal(unname(sim_signal(r, "sm.data")), c(1, 1, 11, 12, 22, 23))
})

test_that("static checks catch structural problems", {
  chk <- sim_check_state_machine(c("a", "b", "c"),
    list(sim_transition("a", "b")), initial = "a")
  expect_equal(chk$unreachable, "c")
  expect_equal(chk$absorbing, c("b", "c"))
  expect_error(sim_check_state_machine(c("a", "b"), list(sim_transition("a", "z"))),
               "unknown state")
  expect_error(sim_check_state_machine(c("a", "b"),
    list(sim_transition("a", "b", priority = 1), sim_transition("a", "a", priority = 1))),
    "same priority")
  expect_error(sim_transition("a", "b", guard = "Sys.time() > 0"), "not allowed")
})

fsm_model <- function(...) {
  m <- sim_model()
  sim_add_block(m, "constant", "c")
  sim_add_block(m, "state_machine", "sm", ...)
  sim_connect(m, "c", "sm")
  m
}

test_that("timed transitions do not depend on where the time axis starts", {
  p <- sim_compile(fsm_model(states = c("a", "b"), period = 0.1,
                             transitions = list(sim_transition("a", "b", after = 1))))
  for (t0 in c(0, 1e6, 1e9)) {
    r <- sim_run(p, t0 + c(0, 2), times = t0 + c(0, 0.5, 0.95, 1, 2), include_ticks = FALSE)
    expect_equal(unname(sim_signal(r, "sm")), c(1, 1, 1, 2, 2), info = format(t0))
  }
})

test_that("held outputs before the first hit follow an overridden initial state", {
  m <- fsm_model(states = c("a", "b"), outputs = list(a = 10, b = 20), offset = 0.5)
  sim_add_block(m, "integrator", "x")
  sim_connect(m, "sm", "x", output = "y")
  r <- sim_run(sim_compile(m), c(0, 0.5), times = c(0, 0.25, 0.5), include_ticks = FALSE,
               initial = list(sm = c(2, 0)), rtol = 1e-10, atol = 1e-12)
  expect_equal(unname(sim_signal(r, "sm.y")), c(20, 20, 20))
  expect_equal(unname(sim_signal(r, "x"))[3], 10, tolerance = 1e-8)
  expect_error(sim_run(sim_compile(m), c(0, 1), initial = list(sm = c(3, 0))), "state index")
})

test_that("action data are matched by name and checked", {
  run_data <- function(...) {
    r <- sim_run(sim_compile(fsm_model(states = "a", ...)), c(0, 1), times = 0,
                 include_ticks = FALSE)
    unname(sim_signal(r, "sm.data"))
  }
  expect_equal(run_data(data = c(a = 1, b = 2), during = list(
    a = function(data, u, t) c(b = unname(data["b"]) + 1, a = unname(data["a"]) + 10))), c(11, 3))
  expect_equal(run_data(data = c(a = 1, b = 2), during = list(
    a = function(data, u, t) data + 1)), c(2, 3))
  expect_error(run_data(data = c(a = 1, b = 2), during = list(
    a = function(data, u, t) c(a = 1, c = 2))), "exactly the data variables")
  expect_error(run_data(data = c(a = 1, b = 2), during = list(
    a = function(data, u, t) 1)), "1 value")
  expect_error(run_data(data = c(a = 0, a = 10), during = list(a = "a <- a + 1")), "duplicate")
  expect_error(run_data(data = c(a = 0, b = 0), during = list(a = "a <- c(3, 4)")),
               "one finite number")
  expect_error(run_data(data = c(a = 0, t = 1)), "reserved")
  expect_error(run_data(data = c(a = Inf)), "finite")
})
