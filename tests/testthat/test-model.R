test_that("models are built in place and blocks get handles", {
  m <- sim_model("demo")
  b <- sim_add_block(m, "gain", "k", gain = 2)
  expect_s3_class(b, "sim_block")
  expect_equal(b$name, "k")
  expect_equal(nrow(sim_blocks(m)), 1L)
  expect_output(print(m), "demo")
})

test_that("invalid model definitions are rejected", {
  m <- sim_model("bad")
  expect_error(sim_add_block(m, "nope", "a"), "Unknown block type")
  sim_add_block(m, "constant", "a")
  expect_error(sim_add_block(m, "constant", "a"), "already has a block")
  expect_error(sim_add_block(m, "constant", "b", valu = 1), "unknown parameter")
  expect_error(sim_add_block(m, "constant", "x/y"), "without '/'")
  expect_error(sim_add_block(m, "constant", "c", value = c(1, NA)), "finite numeric")
  expect_error(sim_add_block(m, "lookup", "lk"), "required parameter")
  sim_add_block(m, "gain", "k")
  expect_error(sim_connect(m, "a", "missing"), "no block named")
  expect_error(sim_connect(m, "a", "k", input = 2), "does not exist")
  expect_error(sim_connect(m, "a", "k", output = 0), "positive integer")
  expect_error(sim_connect(m, "a", "k", input = "zz"), "no input port")
  sim_connect(m, "a", "k")
  expect_error(sim_connect(m, "a", "k"), "already connected")
})

test_that("unconnected inputs and size mismatches fail at compile time", {
  m <- sim_model()
  sim_add_block(m, "constant", "a", value = c(1, 2))
  sim_add_block(m, "sum", "s")
  sim_connect(m, "a", "s", input = 1)
  expect_error(sim_compile(m), "Unconnected input port\\(s\\): s:in2")
  sim_add_block(m, "constant", "b", value = c(1, 2, 3))
  sim_connect(m, "b", "s", input = 2)
  expect_error(sim_compile(m), "incompatible")

  m2 <- sim_model()
  sim_add_block(m2, "constant", "a", value = c(1, 2))
  sim_add_block(m2, "integrator", "x", initial = c(0, 0, 0))
  sim_connect(m2, "a", "x")
  expect_error(sim_compile(m2), "size mismatch")

  m3 <- sim_model()
  sim_add_block(m3, "constant", "a", value = c(1, 2))
  sim_add_block(m3, "gain", "g", gain = matrix(1, 2, 3))
  sim_connect(m3, "a", "g")
  expect_error(sim_compile(m3), "size mismatch")
})

test_that("scalar expansion is explicit and other recycling is refused", {
  m <- sim_model()
  sim_add_block(m, "constant", "v", value = c(1, 2, 3))
  sim_add_block(m, "constant", "s", value = 10)
  sim_add_block(m, "sum", "add")
  sim_connect(m, "v", "add", input = 1)
  sim_connect(m, "s", "add", input = 2)
  r <- sim_run(sim_compile(m), c(0, 1), times = 1)
  expect_equal(unname(sim_signal(r, "add")), c(11, 12, 13))
  expect_error(sim_add_block(m, "step", "st", before = c(0, 0), after = 1), "equal lengths")
})

test_that("disconnect, remove and parameter edits keep the model consistent", {
  m <- sim_model()
  sim_add_block(m, "constant", "a", value = 2)
  sim_add_block(m, "gain", "k", gain = 3)
  sim_connect(m, "a", "k")
  r <- sim_run(sim_compile(m), c(0, 1), times = 1)
  expect_equal(unname(sim_signal(r, "k")), 6)
  sim_set_block_params(m, "k", gain = 5)
  r <- sim_run(sim_compile(m), c(0, 1), times = 1)
  expect_equal(unname(sim_signal(r, "k")), 10)
  sim_disconnect(m, "k")
  expect_error(sim_compile(m), "Unconnected")
  expect_error(sim_disconnect(m, "k"), "not connected")
  sim_remove_block(m, "a")
  expect_equal(sim_blocks(m)$name, "k")
})

test_that("compiled plans do not change when the model is edited", {
  m <- sim_model()
  sim_add_block(m, "constant", "a", value = 2)
  plan <- sim_compile(m)
  sim_set_block_params(m, "a", value = 7)
  r <- sim_run(plan, c(0, 1), times = 1)
  expect_equal(unname(sim_signal(r, "a")), 2)
  m2 <- sim_copy(m, "copy")
  sim_add_block(m2, "gain", "g")
  expect_equal(nrow(sim_blocks(m)), 1L)
})

test_that("ports can be addressed by name", {
  m <- sim_model()
  sim_add_block(m, "constant", "a", value = c(1, 2, 3))
  sim_add_block(m, "demux", "d", sizes = c(1, 2))
  sim_add_block(m, "gain", "g", gain = 2)
  sim_connect(m, "a", "d")
  sim_connect(m, "d", "g", output = "out2")
  r <- sim_run(sim_compile(m), c(0, 1), times = 1)
  expect_equal(unname(sim_signal(r, "g")), c(4, 6))
  expect_equal(unname(sim_signal(r, "d.out1")), 1)
})

test_that("the block library covers the documented categories", {
  types <- sim_block_types()
  expect_true(all(c("constant", "step", "ramp", "sine", "time", "piecewise",
                    "gain", "sum", "subtract", "product", "divide", "saturation",
                    "abs", "math", "integrator", "transfer_function", "state_space",
                    "ode", "unit_delay", "discrete_integrator",
                    "discrete_state_space", "zero_order_hold", "scope",
                    "compare", "switch", "relay", "logic", "stop",
                    "state_machine", "subsystem") %in% types$type))
  expect_true(all(sim_block_types("discrete")$category == "discrete"))
})

test_that("user block types can be registered", {
  sim_register_block_type(
    "test_square", overwrite = TRUE,
    build = function(p) list(inputs = NA, outputs = NA, feedthrough = TRUE,
                             size_out = function(n) n,
                             output = function(t, x, d, u) u[[1]]^2))
  expect_error(sim_register_block_type("test_square", function(p) list()), "already registered")
  m <- sim_model()
  sim_add_block(m, "time", "t")
  sim_add_block(m, "test_square", "sq")
  sim_connect(m, "t", "sq")
  r <- sim_run(sim_compile(m), c(0, 3), times = 0:3)
  expect_equal(unname(sim_signal(r, "sq")), (0:3)^2)
})

test_that("blocks returning wrong output lengths are reported", {
  m <- sim_model()
  sim_add_block(m, "time", "t")
  sim_add_block(m, "fcn", "f", fun = function(u, t) if (t > 0.5) c(u, u) else u, size = 1)
  sim_connect(m, "t", "f")
  expect_error(sim_run(sim_compile(m), c(0, 1), times = c(0, 1)), "length 2")
})

test_that("signal names are unique across blocks, ports and subsystem outputs", {
  sub <- sim_subsystem("two_out")
  sim_add_block(sub, "inport", "i")
  sim_add_block(sub, "outport", "o1", index = 1)
  sim_add_block(sub, "outport", "o2", index = 2)
  sim_connect(sub, "i", "o1")
  sim_connect(sub, "i", "o2")
  m <- sim_model()
  sim_add_block(m, "constant", "u")
  sim_add_block(m, "subsystem", "s", system = sub)
  sim_add_block(m, "constant", "s.o2", value = 5)
  sim_connect(m, "u", "s")
  expect_error(sim_compile(m), "'s.o2' is used twice")

  m <- sim_model()
  for (nm in c("time", "kind")) expect_error(sim_add_block(m, "constant", nm), "reserved")
  for (nm in c("a[1]", "state:x", "a/b")) expect_error(sim_add_block(m, "constant", nm), "without")
})

test_that("forwarded arguments bind by exact name", {
  sim_register_block_type("repeat_one", params = list(n = 1), overwrite = TRUE,
    build = function(p) list(outputs = p$n, output = function(t, x, d, u) rep(1, p$n)))
  m <- sim_model()
  add <- function(...) sim_add_block(m, "repeat_one", "r", ...)
  add(n = 4)
  expect_length(sim_signal(sim_run(sim_compile(m), c(0, 1), times = 1), "r"), 4)
  wrap <- function(...) sim_add_block(...)
  wrap(m, "repeat_one", "q", n = 2)
  expect_length(sim_signal(sim_run(sim_compile(m), c(0, 1), times = 1), "q"), 2)
})

test_that("block handles only work on their own model and its copies", {
  a <- sim_model("a")
  h <- sim_add_block(a, "constant", "c", value = 9)
  b <- sim_model("b")
  sim_add_block(b, "constant", "c")
  sim_add_block(b, "gain", "g")
  expect_error(sim_connect(b, h, "g"), "belongs to model 'a'")
  expect_error(sim_add_event(b, "e", sim_zero_crossing(h), terminal = TRUE), "belongs to model 'a'")
  a2 <- sim_copy(a)
  sim_add_block(a2, "gain", "g")
  expect_silent(sim_connect(a2, h, "g"))
})

test_that("unnamed parameter edits and fractional port sizes are errors", {
  m <- sim_model()
  sim_add_block(m, "constant", "c", value = 3)
  expect_error(sim_set_block_params(m, "c", 100), "must be named")
  sim_register_block_type("half_port", overwrite = TRUE,
    build = function(p) list(outputs = 2.5, output = function(t, x, d, u) c(1, 2)))
  expect_error(sim_add_block(m, "half_port", "h"), "whole numbers")
})
