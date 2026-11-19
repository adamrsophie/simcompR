roundtrip <- function(m) {
  f <- withr::local_tempfile(fileext = ".json")
  sim_write_json(m, f)
  sim_read_json(f)
}

test_that("portable models round-trip exactly through JSON", {
  m <- sim_model("rt", params = list(k = 1 / 3, e = 0.8))
  sim_add_block(m, "step", "u", time = 0.3, after = 2)
  sim_add_block(m, "gain", "g", gain = sim_param("k"))
  sim_add_block(m, "state_space", "ss", A = matrix(c(-1, 0.5, 0, -2), 2), B = matrix(c(1, 0), 2),
                C = matrix(c(1, 1), 1), initial = c(0.1, 0.2))
  sim_add_block(m, "discrete_transfer_function", "dz", num = c(0.2, 0.1), den = c(1, -0.7),
                period = 0.1)
  sim_add_block(m, "expression", "ex", expr = "u[1]^2")
  sim_add_block(m, "scope", "sc")
  sim_connect(m, "u", "g"); sim_connect(m, "g", "ss"); sim_connect(m, "ss", "dz")
  sim_connect(m, "dz", "ex"); sim_connect(m, "ex", "sc")
  sim_add_event(m, "bump", sim_at(1), reset = list(ss = "x * e"))
  m2 <- roundtrip(m)
  expect_identical(sim_model_to_list(m2), sim_model_to_list(m))
  p1 <- sim_compile(m); p2 <- sim_compile(m2)
  r1 <- sim_run(p1, c(0, 2)); r2 <- sim_run(p2, c(0, 2))
  expect_identical(r1$signals, r2$signals)
  expect_identical(r1$events, r2$events)
})

test_that("subsystems and state machines round-trip", {
  lag <- sim_subsystem("lag", params = list(rate = 2))
  sim_add_block(lag, "inport", "u", index = 1)
  sim_add_block(lag, "transfer_function", "tf", num = 1, den = c(1, 1))
  sim_add_block(lag, "outport", "y", index = 1)
  sim_connect(lag, "u", "tf"); sim_connect(lag, "tf", "y")
  m <- sim_model("outer")
  sim_add_block(m, "sine", "s")
  sim_add_block(m, "subsystem", "a", system = lag, params = list(rate = 3))
  sim_add_block(m, "state_machine", "sm", period = 0.5, states = c("lo", "hi"),
    transitions = list(sim_transition("lo", "hi", "u[1] > 0.5"),
                       sim_transition("hi", "lo", "u[1] < -0.5", after = 0.5)),
    data = c(n = 0), entry = list(hi = "n <- n + 1"))
  sim_connect(m, "s", "a"); sim_connect(m, "a", "sm")
  m2 <- roundtrip(m)
  r1 <- sim_run(sim_compile(m), c(0, 10))
  r2 <- sim_run(sim_compile(m2), c(0, 10))
  expect_identical(r1$signals, r2$signals)
})

test_that("non-portable content is refused instead of dropped", {
  m <- sim_model()
  sim_add_block(m, "fcn", "f", fun = function(u, t) u)
  expect_error(sim_model_to_list(m), "cannot be written")
  m2 <- ball_model()
  m2$events$impact$reset$v <- function(x, s, t) -x
  expect_error(sim_model_to_list(m2), "string expression")
  m3 <- sim_model()
  sim_add_block(m3, "state_machine", "sm", states = "a",
                transitions = list(sim_transition("a", "a", guard = function(...) TRUE)))
  expect_error(sim_model_to_list(m3), "function guards")
})

test_that("reading rejects foreign files, future schemas and unsafe expressions", {
  f <- withr::local_tempfile(fileext = ".json")
  writeLines('{"schema": "other"}', f)
  expect_error(sim_read_json(f), "Not a simcompR model")
  writeLines('{"schema": "simcompR.model", "schema_version": 99}', f)
  expect_error(sim_read_json(f), "schema version 99")
  m <- sim_model()
  sim_add_block(m, "expression", "ex", expr = "u + 1")
  sim_add_block(m, "time", "t")
  sim_connect(m, "t", "ex")
  x <- sim_model_to_list(m)
  x$blocks[[1]]$params$expr$`$chr`[[1]] <- "system('echo hacked')"
  expect_error(sim_model_from_list(x), "not allowed")
})

roundtrip <- function(m) {
  f <- withr::local_tempfile(fileext = ".json")
  sim_write_json(m, f)
  sim_read_json(f)
}

test_that("values round-trip with their type, missing values and names", {
  vals <- list(chr = c("a", NA, "b"), chr0 = character(), chr_na_first = c(NA, "x"),
               int = c(1L, NA, 3L), dbl = c(a = 1.5, b = NA, c = NaN, d = Inf),
               lgl = c(TRUE, NA), named_na = stats::setNames(1:2, c("p", NA)),
               mat = matrix(c(0.1, 2, 3, 4), 2), lst = list(a = "x", b = list(NA_character_)))
  m <- sim_model(params = vals)
  back <- roundtrip(m)$params
  for (nm in names(vals)) expect_identical(back[[nm]], vals[[nm]], info = nm)
})

test_that("a missing character value keeps model behavior after a round trip", {
  m <- sim_model(params = list(labels = c("a", NA_character_, "b")))
  sim_add_block(m, "constant", "c")
  sim_add_block(m, "expression", "e", expr = "length(labels)")
  sim_connect(m, "c", "e")
  run <- function(mm) unname(sim_signal(sim_run(sim_compile(mm), c(0, 1), times = 1), "e"))
  expect_equal(run(roundtrip(m)), run(m))
})

test_that("values that would change on reload are rejected when writing", {
  bad <- list(arr = array(1:8, c(2, 2, 2)), chr_mat = matrix("a", 1, 1),
              dimnamed = matrix(1, 1, 1, dimnames = list("r", "c")),
              fac = factor("a"), df = data.frame(a = 1), date = Sys.Date())
  for (nm in names(bad)) {
    m <- sim_model(params = stats::setNames(list(bad[[nm]]), nm))
    expect_error(sim_model_to_list(m), "cannot be written", info = nm)
  }
})

test_that("schema versions must be whole numbers in range", {
  x <- sim_model_to_list(sim_model())
  for (v in list(1.9, 0, -1, "1", NA)) {
    x$schema_version <- v
    expect_error(sim_model_from_list(x), "schema version", info = format(v))
  }
  x$schema_version <- 2
  expect_error(sim_model_from_list(x), "reads version")
})
