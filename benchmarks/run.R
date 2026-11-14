# Timed scenarios. Run with an installed simcompR: Rscript benchmarks/run.R
library(simcompR)

time_it <- function(label, plan, tspan, ...) {
  tm <- system.time(r <- sim_run(plan, tspan, ...)); el <- tm[["user.self"]] + tm[["sys.self"]]
  st <- sim_diagnostics(r)$stats
  cat(sprintf("%-34s %8.3f s CPU segments %6d  rhs %7d\n", label, el,
              as.integer(st[["segments"]]), as.integer(st[["rhs_evaluations"]])))
  invisible(r)
}

# 1. Sampled feedback loop at 1 kHz for one second (1000 segments)
m <- sim_model("sampled_loop")
sim_add_block(m, "constant", "r", value = 1)
sim_add_block(m, "subtract", "e")
sim_add_block(m, "zero_order_hold", "h", period = 0.001)
sim_add_block(m, "gain", "k", gain = 2)
sim_add_block(m, "integrator", "x")
sim_connect(m, "r", "e", input = 1)
sim_connect(m, "x", "e", input = 2)
sim_connect(m, "e", "h")
sim_connect(m, "h", "k")
sim_connect(m, "k", "x")
time_it("1 kHz sampled loop, 1 s", sim_compile(m), c(0, 1), include_ticks = FALSE)

# 2. Continuous mass-spring-damper, 100 s
msd <- sim_model("msd")
sim_add_block(msd, "step", "f", time = 1)
sim_add_block(msd, "transfer_function", "G", num = 1, den = c(1, 0.4, 4))
sim_connect(msd, "f", "G")
time_it("continuous second order, 100 s", sim_compile(msd), c(0, 100))

# 3. Bouncing ball until the event limit
ball <- sim_model("ball", params = list(e = 0.8))
sim_add_block(ball, "constant", "g", value = -9.81)
sim_add_block(ball, "integrator", "v")
sim_add_block(ball, "integrator", "h", initial = 1)
sim_connect(ball, "g", "v")
sim_connect(ball, "v", "h")
sim_add_event(ball, "impact", sim_zero_crossing("h", "falling"), reset = list(v = "-e * x"))
time_it("bouncing ball, 3 s", sim_compile(ball), c(0, 3), max_events = 200)

# 4. Saturated sine with located kinks, 100 periods
sat <- sim_model("sat")
sim_add_block(sat, "sine", "u")
sim_add_block(sat, "saturation", "s", lower = -0.5, upper = 0.5)
sim_add_block(sat, "integrator", "x")
sim_connect(sat, "u", "s")
sim_connect(sat, "s", "x")
time_it("saturation kinks, 100 periods", sim_compile(sat), c(0, 200 * pi))

# 5. Discrete-only loop, with and without a saturation (no continuous state,
#    so kinks must cost nothing)
dl <- function(with_sat) {
  m <- sim_model("discrete_loop")
  sim_add_block(m, "sine", "u")
  sim_add_block(m, "unit_delay", "z", period = 0.01)
  if (with_sat) {
    sim_add_block(m, "saturation", "s", lower = -0.5, upper = 0.5)
    sim_connect(m, "u", "s")
    sim_connect(m, "s", "z")
  } else {
    sim_connect(m, "u", "z")
  }
  sim_compile(m)
}
time_it("discrete loop, 10 s", dl(FALSE), c(0, 10))
time_it("discrete loop with saturation, 10 s", dl(TRUE), c(0, 10))
