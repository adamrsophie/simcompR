test_that("sweeps reuse the plan and match individual runs", {
  m <- decay_model(1)
  s <- sim_sweep(m, data.frame(a = c(0.5, 2)), tspan = c(0, 1), times = c(0, 1),
                 rtol = 1e-10, atol = 1e-12,
                 summary = function(r) c(end = unname(tail(sim_signal(r, "x"), 1))))
  expect_s3_class(s, "sim_sweep")
  expect_equal(s$summary$end, exp(-c(0.5, 2)), tolerance = 1e-8)
  expect_equal(s$status, c("success", "success"))
  expect_length(s$results, 2)
  expect_output(print(s), "2 run")
  expect_error(sim_sweep(m, data.frame(b = 1), tspan = c(0, 1)), "Unknown model parameter")
})

test_that("Monte Carlo runs are reproducible and keep the global stream", {
  m <- decay_model(1)
  set.seed(123)
  before <- .Random.seed
  run <- function() sim_monte_carlo(m, n = 4, sampler = function(i) list(a = runif(1)),
                                    seed = 7, tspan = c(0, 1), times = 1, keep = FALSE,
                                    summary = function(r) c(end = unname(sim_signal(r, "x"))))
  a <- run(); b <- run()
  expect_identical(a$summary, b$summary)
  expect_identical(.Random.seed, before)
  expect_null(a$results)
  expect_equal(a$summary$end, exp(-a$summary$a), tolerance = 1e-6)
})

test_that("white noise is reproducible, seeded and held between hits", {
  m <- sim_model()
  sim_add_block(m, "white_noise", "w", sd = 2, period = 0.5, seed = 11)
  r1 <- sim_run(sim_compile(m), c(0, 5), times = seq(0, 5, 0.25), include_ticks = FALSE)
  set.seed(1); x <- runif(1)
  r2 <- sim_run(sim_compile(m), c(0, 5), times = seq(0, 5, 0.25), include_ticks = FALSE)
  expect_identical(r1$signals, r2$signals)
  w <- sim_signal(r1, "w")
  expect_equal(w[c(1, 3, 5)], w[c(2, 4, 6)])
  expect_false(isTRUE(all.equal(w[1], w[3])))
  sim_set_block_params(m, "w", seed = 12)
  r3 <- sim_run(sim_compile(m), c(0, 5), times = seq(0, 5, 0.25))
  expect_false(identical(r1$signals[1, ], r3$signals[1, ]))
})

test_that("failed runs inside a sweep are reported, not hidden", {
  m <- sim_model(params = list(c = 1))
  sim_add_block(m, "time", "t")
  sim_add_block(m, "fcn", "f", fun = function(u, t) if (t > 0.5 && u > 0) stop("boom") else u, size = 1)
  sim_add_block(m, "gain", "g", gain = sim_param("c"))
  sim_connect(m, "t", "g"); sim_connect(m, "g", "f")
  expect_warning(s <- sim_sweep(m, data.frame(c = c(-1, 1)), tspan = c(0, 1),
                                summary = function(r) c(n = length(r$time))), "boom")
  expect_equal(s$status, c("success", "failed"))
  expect_true(is.na(s$summary$n[2]))
})
