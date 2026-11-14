test_that("tarjan_scc and condensed_order find components deterministically", {
  # 1 -> 2 -> 3 -> 1, 3 -> 4, 5 isolated, 4 -> 4
  adj <- list(2L, 3L, c(1L, 4L), 4L, integer())
  comps <- simcompR:::tarjan_scc(5L, adj)
  sets <- lapply(comps, sort)
  expect_true(list(1:3) %in% sets)
  expect_true(list(4L) %in% sets)
  ord <- simcompR:::condensed_order(5L, adj, comps)
  expect_equal(ord[[1]], 1:3)
  pos <- vapply(ord, min, 1L)
  expect_lt(which(pos == 1L), which(pos == 4L))
  expect_equal(simcompR:::find_cycle(1:3, adj), c(1L, 2L, 3L, 1L))
})

test_that("feedback through an integrator is not an algebraic loop", {
  plan <- sim_compile(connect_feedback(feedback_model(1:3)))
  expect_length(plan$diagnostics$loops, 0)
  # plant output does not depend on its input, so it comes first
  expect_lt(match("plant", plan$diagnostics$order), match("error", plan$diagnostics$order))
  expect_false(any(plan$diagnostics$edges$to == "plant"))
})

test_that("direct-feedthrough cycles are detected and reported", {
  m <- sim_model()
  sim_add_block(m, "constant", "u")
  sim_add_block(m, "sum", "s", signs = "+-")
  sim_add_block(m, "gain", "a", gain = 2)
  sim_add_block(m, "gain", "b", gain = 0.5)
  sim_connect(m, "u", "s", input = 1)
  sim_connect(m, "b", "s", input = 2)
  sim_connect(m, "s", "a")
  sim_connect(m, "a", "b")
  expect_error(sim_compile(m), "Algebraic loop.*a -> b -> s -> a")

  self <- sim_model()
  sim_add_block(self, "sum", "s")
  sim_add_block(self, "constant", "c")
  sim_connect(self, "c", "s", input = 1)
  sim_connect(self, "s", "s", input = 2)
  expect_error(sim_compile(self), "s -> s")
})

test_that("unit delays and state-space blocks with D = 0 break loops", {
  m <- sim_model()
  sim_add_block(m, "constant", "u")
  sim_add_block(m, "sum", "s")
  sim_add_block(m, "unit_delay", "z", period = 1)
  sim_connect(m, "u", "s", input = 1)
  sim_connect(m, "z", "s", input = 2)
  sim_connect(m, "s", "z")
  expect_s3_class(sim_compile(m), "sim_plan")

  ss <- function(D) {
    m <- sim_model(params = list(D = D))
    sim_add_block(m, "constant", "u")
    sim_add_block(m, "sum", "s", signs = "+-")
    sim_add_block(m, "state_space", "g", A = -1, B = 1, C = 1, D = sim_param("D"))
    sim_connect(m, "u", "s", input = 1)
    sim_connect(m, "g", "s", input = 2)
    sim_connect(m, "s", "g")
    m
  }
  expect_s3_class(sim_compile(ss(0)), "sim_plan")
  expect_error(sim_compile(ss(1)), "Algebraic loop")
  # changing D through parameters re-analyzes the graph
  p <- sim_compile(ss(0), algebraic_loops = "solve")
  expect_true(sim_set_params(p, D = 0)$structure_reused)
  p1 <- sim_set_params(p, D = 1)
  expect_false(p1$structure_reused)
  expect_length(p1$graph$loops, 1)
})

test_that("linear and nonlinear algebraic loops are solved with diagnostics", {
  m <- sim_model()
  sim_add_block(m, "constant", "u", value = 1)
  sim_add_block(m, "sum", "s", signs = "+-")
  sim_add_block(m, "gain", "k", gain = 2)
  sim_connect(m, "u", "s", input = 1)
  sim_connect(m, "k", "s", input = 2)
  sim_connect(m, "s", "k")
  r <- sim_run(sim_compile(m, algebraic_loops = "solve"), c(0, 1), times = c(0, 1))
  expect_equal(unname(sim_signal(r, "s")), rep(1 / 3, 2), tolerance = 1e-10)
  lr <- sim_diagnostics(r)$loop_residuals
  expect_true(all(lr[, "residual"] <= 1e-10))

  # y = cos(y): fixed point 0.739085...
  m2 <- sim_model()
  sim_add_block(m2, "math", "c", fun = "cos")
  sim_add_block(m2, "gain", "y", gain = 1)
  sim_connect(m2, "c", "y")
  sim_connect(m2, "y", "c")
  r2 <- sim_run(sim_compile(m2, algebraic_loops = "solve"), c(0, 1), times = 1)
  expect_equal(unname(sim_signal(r2, "y")), 0.7390851332151607, tolerance = 1e-9)

  # an inconsistent loop y = y + 1 has no solution and must fail clearly
  m3 <- sim_model()
  sim_add_block(m3, "constant", "one")
  sim_add_block(m3, "sum", "s")
  sim_connect(m3, "one", "s", input = 1)
  sim_connect(m3, "s", "s", input = 2)
  expect_error(sim_run(sim_compile(m3, algebraic_loops = "solve"), c(0, 1)),
               "singular Jacobian|did not converge")
})

test_that("algebraic loops inside continuous dynamics integrate correctly", {
  # x' = y, y = -x - y  (loop through s and k), so x' = -x / 2
  m <- sim_model()
  sim_add_block(m, "integrator", "x", initial = 1)
  sim_add_block(m, "sum", "s", signs = "--")
  sim_add_block(m, "gain", "k", gain = 1)
  sim_connect(m, "x", "s", input = 1)
  sim_connect(m, "k", "s", input = 2)
  sim_connect(m, "s", "k")
  sim_connect(m, "s", "x")
  r <- tight(sim_compile(m, algebraic_loops = "solve"), c(0, 2), times = c(0, 1, 2))
  expect_lt(max(abs(sim_signal(r, "x") - exp(-c(0, 1, 2) / 2))), 1e-7)
})

test_that("sizes propagate through feedback loops", {
  m <- sim_model()
  sim_add_block(m, "constant", "r", value = c(1, 2))
  sim_add_block(m, "sum", "e", signs = "+-")
  sim_add_block(m, "integrator", "x", initial = c(0, 0))
  sim_connect(m, "r", "e", input = 1)
  sim_connect(m, "x", "e", input = 2)
  sim_connect(m, "e", "x")
  plan <- sim_compile(m)
  expect_equal(plan$sizes$out$e, 2L)
  r <- tight(plan, c(0, 1), times = 1)
  expect_equal(unname(sim_signal(r, "x")), c(1, 2) * (1 - exp(-1)), tolerance = 1e-9)
})

test_that("DOT export marks non-feedthrough edges", {
  dot <- sim_to_dot(connect_feedback(feedback_model(1:3)))
  expect_match(dot, '"error" -> "plant" \\[label="1:1", style=dashed\\]')
  expect_match(dot, '"plant" -> "error" \\[label="1:2"\\]')
})
