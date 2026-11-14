test_that("kinks of memoryless blocks are located exactly", {
  run <- function(type, ..., connect = NULL) {
    m <- sim_model()
    sim_add_block(m, "sine", "u")
    sim_add_block(m, type, "k", ...)
    if (is.null(connect)) sim_connect(m, "u", "k") else connect(m)
    sim_add_block(m, "integrator", "x")
    sim_connect(m, "k", "x")
    r <- sim_run(sim_compile(m), c(0, 2 * pi), times = 2 * pi, rtol = 1e-10, atol = 1e-12,
                 include_ticks = FALSE)
    unname(tail(sim_signal(r, "x"), 1))
  }
  expect_equal(run("saturation", lower = -0.5, upper = 0.5), 0, tolerance = 1e-8)
  expect_equal(run("abs"), 4, tolerance = 1e-8)
  expect_equal(run("dead_zone", lower = -0.5, upper = 0.5), 0, tolerance = 1e-8)
  expect_equal(run("compare", op = ">", threshold = 0), pi, tolerance = 1e-8)
  expect_equal(run("lookup", x = c(-1, 0, 1), y = c(0, 1, 0)), 2 * pi - 4, tolerance = 1e-8)
  expect_equal(run("math", fun = "floor"), -pi, tolerance = 1e-8)
  expect_equal(run("minmax", fun = "max", n_inputs = 2, connect = function(m) {
    sim_add_block(m, "constant", "z", value = 0)
    sim_connect(m, "u", "k", input = 1)
    sim_connect(m, "z", "k", input = 2)
  }), 2, tolerance = 1e-8)
})

test_that("kinks are silent and optional", {
  m <- sim_model()
  sim_add_block(m, "sine", "u")
  sim_add_block(m, "saturation", "s", lower = -0.5, upper = 0.5)
  sim_add_block(m, "integrator", "x")
  sim_connect(m, "u", "s")
  sim_connect(m, "s", "x")
  p <- sim_compile(m)
  r <- sim_run(p, c(0, 2 * pi))
  expect_equal(nrow(sim_events(r)), 0)
  expect_gt(sim_diagnostics(r)$stats[["roots"]], 0)
  r4 <- sim_run(p, c(0, 1), method = "rk4", times = c(0, 1))
  expect_equal(r4$status, "success")
})

test_that("an event that starts on its level fires when it leaves in its direction", {
  m <- sim_model()
  sim_add_block(m, "constant", "c", value = -1)
  sim_add_block(m, "integrator", "x")
  sim_connect(m, "c", "x")
  sim_add_event(m, "e", sim_zero_crossing("x", "falling", level = 0), reset = list(x = "1"))
  r <- sim_run(sim_compile(m), c(0, 2.5))
  expect_equal(sim_events(r)$time, c(0, 1, 2), tolerance = 1e-7)
  expect_equal(r$kind[r$time == 0], c("pre-event", "post-event"))
  # resting on the level fires nothing
  m2 <- sim_model()
  sim_add_block(m2, "constant", "z", value = 0)
  sim_add_block(m2, "integrator", "x")
  sim_connect(m2, "z", "x")
  sim_add_event(m2, "e", sim_zero_crossing("x", "rising"), terminal = TRUE)
  expect_equal(sim_run(sim_compile(m2), c(0, 1))$status, "success")
})

test_that("top-level events reach blocks inside subsystems", {
  lag <- sim_subsystem("lag")
  sim_add_block(lag, "inport", "u")
  sim_add_block(lag, "integrator", "y")
  sim_add_block(lag, "outport", "out")
  sim_connect(lag, "u", "y")
  sim_connect(lag, "y", "out")
  top <- sim_model()
  sim_add_block(top, "constant", "c")
  sim_add_block(top, "subsystem", "s1", system = lag)
  sim_connect(top, "c", "s1")
  sim_add_event(top, "zero", sim_at(1), reset = list("s1/y" = "0"))
  sim_add_event(top, "half", sim_zero_crossing("s1/y", "rising", level = 0.5),
                reset = list("s1/y" = "x + 10"))
  r <- sim_run(sim_compile(top), c(0, 2), times = 2)
  expect_equal(sim_events(r)$time, c(0.5, 1, 1.5), tolerance = 1e-7)
  expect_equal(unname(tail(sim_signal(r, "s1"), 1)), 11, tolerance = 1e-7)
  expect_error(sim_add_event(top, "bad", sim_at(1), reset = list("s1/nope" = "0")), "nope")
  expect_error(sim_add_event(top, "bad", sim_at(1), reset = list("c/y" = "0")),
               "not a subsystem")
  sim_remove_event(top, "zero")
  expect_false("zero" %in% names(top$events))
  expect_error(sim_remove_event(top, "zero"), "no event")
})

test_that("loop starting values select a solution and NaN is reported clearly", {
  m <- sim_model()
  sim_add_block(m, "constant", "four", value = 4)
  sim_add_block(m, "divide", "y")
  sim_connect(m, "four", "y", input = 1)
  sim_connect(m, "y", "y", input = 2)
  expect_error(sim_run(sim_compile(m, algebraic_loops = "solve"), c(0, 1), times = 1),
               "loop_guess")
  for (g in c(1, -1)) {
    p <- sim_compile(m, algebraic_loops = "solve", loop_guess = list(y = g))
    expect_equal(unname(sim_signal(sim_run(p, c(0, 1), times = 1), "y")), 2 * g)
    expect_equal(unname(sim_signal(sim_run(sim_set_params(p), c(0, 1), times = 1), "y")), 2 * g)
  }
  expect_error(sim_compile(m, algebraic_loops = "solve", loop_guess = list(zz = 1)), "zz")
  m2 <- sim_model()
  sim_add_block(m2, "expression", "y", expr = "log(u[1] - 5)")
  sim_connect(m2, "y", "y")
  expect_error(suppressWarnings(sim_run(sim_compile(m2, algebraic_loops = "solve"), c(0, 1),
                                        times = 1)), "not finite")
})

test_that("short parameter names are not captured by partial matching", {
  sim_register_block_type("rep_n", params = list(n = 1), overwrite = TRUE,
    build = function(p) list(inputs = 1L, outputs = p$n, feedthrough = TRUE,
                             output = function(t, x, d, u) rep(u[[1]], p$n)))
  m <- sim_model(params = list(p = 1))
  sim_add_block(m, "constant", "c", value = sim_param("p"))
  sim_add_block(m, "rep_n", "r", n = 3)
  sim_connect(m, "c", "r")
  expect_length(sim_signal(sim_run(sim_compile(m), c(0, 1), times = 1), "r"), 3)
  sim_set_block_params(m, "r", n = 2)
  plan <- sim_set_params(sim_compile(m), p = 5)
  expect_equal(unname(sim_signal(sim_run(plan, c(0, 1), times = 1), "r")), c(5, 5))
})

test_that("parameter expressions, nested references and scoped expressions", {
  m <- sim_model("decay", params = list(a = 2, tau = 4))
  sim_add_block(m, "integrator", "x", initial = 1)
  sim_add_block(m, "gain", "k", gain = -sim_param("a") / sim_param("tau"))
  sim_connect(m, "x", "k")
  sim_connect(m, "k", "x")
  sim_add_block(m, "expression", "e", expr = "a * u[1] + tau")
  sim_connect(m, "x", "e")
  sim_add_block(m, "lookup", "lk", x = sim_param_expr("c(0, tau)"), y = c(0, 1))
  sim_connect(m, "x", "lk")
  p <- sim_compile(m)
  out <- function(plan) {
    r <- sim_run(plan, c(0, 1), times = 1, rtol = 1e-10, atol = 1e-12)
    unname(sim_signal(r, c("x", "e", "lk")))
  }
  expect_equal(out(p), c(exp(-0.5), 2 * exp(-0.5) + 4, exp(-0.5) / 4), tolerance = 1e-8)
  expect_equal(out(sim_set_params(p, a = 4))[1:2], c(exp(-1), 4 * exp(-1) + 4), tolerance = 1e-8)
  f <- withr::local_tempfile(fileext = ".json")
  sim_write_json(m, f)
  expect_equal(out(sim_compile(sim_read_json(f))), out(p))
  expect_s3_class(sqrt(2 * sim_param("a")), "sim_param")
  expect_error(sim_param("a") %% 2, "not available")
  expect_error(sim_param("a") + 1:2, "single finite")
  # a nested reference inside a list parameter
  sm <- sim_model(params = list(hi = 3))
  sim_add_block(sm, "constant", "c")
  sim_add_block(sm, "state_machine", "s", states = c("a", "b"),
                outputs = list(a = 0, b = sim_param("hi")),
                transitions = list(sim_transition("a", "b", "t >= hi - 2.5")))
  sim_connect(sm, "c", "s")
  r <- sim_run(sim_compile(sm), c(0, 1), times = 1)
  expect_equal(unname(tail(sim_signal(r, "s.y"), 1)), 3)
})

test_that("plans keep the block types they were compiled with", {
  sim_register_block_type("dbl_t", overwrite = TRUE, build = function(p) list(
    inputs = 1L, outputs = 1L, feedthrough = TRUE, output = function(t, x, d, u) 2 * u[[1]]))
  m <- sim_model(params = list(q = 1))
  sim_add_block(m, "constant", "c", value = sim_param("q"))
  sim_add_block(m, "dbl_t", "d")
  sim_connect(m, "c", "d")
  p <- sim_compile(m)
  sim_register_block_type("dbl_t", overwrite = TRUE, build = function(p) list(
    inputs = 1L, outputs = 1L, feedthrough = TRUE, output = function(t, x, d, u) 3 * u[[1]]))
  expect_equal(unname(sim_signal(sim_run(sim_set_params(p, q = 1), c(0, 1), times = 1), "d")), 2)
})

test_that("subsystem instance parameters are checked", {
  lag <- sim_subsystem("lag", params = list(rate = 1))
  sim_add_block(lag, "inport", "u")
  sim_add_block(lag, "gain", "k", gain = sim_param("rate"))
  sim_add_block(lag, "outport", "o")
  sim_connect(lag, "u", "k")
  sim_connect(lag, "k", "o")
  m <- sim_model()
  sim_add_block(m, "constant", "c")
  sim_add_block(m, "subsystem", "s", system = lag, params = list(rat = 10))
  sim_connect(m, "c", "s")
  expect_error(sim_compile(m), "rat")
  sim_set_block_params(m, "s", params = list(rate = 10))
  expect_equal(unname(sim_signal(sim_run(sim_compile(m), c(0, 1), times = 1), "s")), 10)
})

test_that("result data frames have unique columns and include every state", {
  m <- sim_model()
  sim_add_block(m, "constant", "one")
  sim_add_block(m, "integrator", "x")
  sim_add_block(m, "unit_delay", "z", period = 0.5)
  sim_connect(m, "one", "x")
  sim_connect(m, "x", "z")
  df <- as.data.frame(sim_run(sim_compile(m), c(0, 1), times = c(0, 1)))
  expect_false(anyDuplicated(names(df)) > 0)
  expect_true(all(c("time", "kind", "x", "state:x", "dstate:z") %in% names(df)))
})

test_that("scheduled times outside tspan are reported", {
  m <- sim_model()
  sim_add_block(m, "constant", "one")
  sim_add_block(m, "integrator", "x")
  sim_connect(m, "one", "x")
  sim_add_event(m, "late", sim_at(5), reset = list(x = "0"))
  expect_warning(sim_run(sim_compile(m), c(0, 1)), "outside tspan")
})

test_that("autonomous ode blocks need no input", {
  m <- sim_model()
  sim_add_block(m, "ode", "y", n_inputs = 0, initial = 1,
                deriv = function(t, x, u) -x)
  r <- sim_run(sim_compile(m), c(0, 1), times = 1, rtol = 1e-10, atol = 1e-12)
  expect_equal(unname(sim_signal(r, "y")), exp(-1), tolerance = 1e-8)
  expect_error(sim_add_block(m, "ode", "w", n_inputs = 2, initial = 1,
                             deriv = function(t, x, u) -x), "0 or 1")
})
