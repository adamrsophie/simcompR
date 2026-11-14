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
