# Shared model builders for the tests.

feedback_model <- function(order = 1) {
  m <- sim_model("feedback")
  add <- list(
    function() sim_add_block(m, "constant", "reference", value = 1),
    function() sim_add_block(m, "sum", "error", signs = "+-"),
    function() sim_add_block(m, "integrator", "plant", initial = 0)
  )
  for (i in order) add[[i]]()
  m
}

connect_feedback <- function(m, order = 1:3) {
  cons <- list(
    function() sim_connect(m, "reference", "error", input = 1),
    function() sim_connect(m, "plant", "error", input = 2),
    function() sim_connect(m, "error", "plant")
  )
  for (i in order) cons[[i]]()
  m
}

decay_model <- function(a = 1, x0 = 1) {
  m <- sim_model("decay", params = list(a = a))
  sim_add_block(m, "integrator", "x", initial = x0)
  sim_add_block(m, "gain", "k", gain = sim_param("a"))
  sim_add_block(m, "gain", "neg", gain = -1)
  sim_connect(m, "x", "k")
  sim_connect(m, "k", "neg")
  sim_connect(m, "neg", "x")
  m
}

ball_model <- function(e = 0.8, h0 = 1) {
  m <- sim_model("ball", params = list(e = e))
  sim_add_block(m, "constant", "gravity", value = -9.81)
  sim_add_block(m, "integrator", "v", initial = 0)
  sim_add_block(m, "integrator", "h", initial = h0)
  sim_connect(m, "gravity", "v")
  sim_connect(m, "v", "h")
  sim_add_event(m, "impact", trigger = sim_zero_crossing("h", direction = "falling"),
                reset = list(v = "-e * x"))
  m
}

tight <- function(plan, tspan, ...) {
  sim_run(plan, tspan, rtol = 1e-10, atol = 1e-12, ...)
}
