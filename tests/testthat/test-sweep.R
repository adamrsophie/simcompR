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

test_that("summary values are matched by name, never by position", {
  m <- sim_model(params = list(a = 1, b = 2))
  sim_add_block(m, "constant", "c", value = sim_param("a"))
  sw <- function(f) sim_sweep(m, data.frame(a = c(1, 2, 3)), c(0, 1), times = 1, summary = f)
  s <- sw(function(r) switch(r$metadata$params$a, c(mean = 10, sd = 20), c(sd = 30, mean = 40),
                             c(mean = 50)))
  expect_equal(s$summary$mean, c(10, 40, 50))
  expect_equal(s$summary$sd, c(20, 30, NA))
  s <- sw(function(r) if (r$metadata$params$a == 2) c(extra = 1) else c(y = 1))
  expect_equal(s$summary$y, c(1, NA, 1))
  expect_equal(s$summary$extra, c(NA, 1, NA))
  s <- sw(function(r) if (r$metadata$params$a == 2) c(1, 2) else c(y = 1))
  expect_equal(s$status, c("success", "failed", "success"))
  expect_match(s$messages[2], "name")
  s <- sw(function(r) if (r$metadata$params$a == 1) c(y = 1, y = 2) else c(y = 1))
  expect_match(s$messages[1], "duplicate")
  s <- sw(function(r) stop("no"))
  expect_equal(s$status, rep("failed", 3))
  expect_null(s$summary)
  s <- sw(function(r) if (r$metadata$params$a == 1) stop("no") else c(y = 5))
  expect_equal(s$summary$y, c(NA, 5, 5))
})

test_that("partial and empty parameter sets keep one row per run", {
  m <- sim_model(params = list(a = 1, b = 2))
  sim_add_block(m, "constant", "c", value = sim_param("a"))
  sim_add_block(m, "constant", "d", value = sim_param("b"))
  s <- sim_sweep(m, list(list(a = 3), list(b = 4)), c(0, 1), times = 1,
                 summary = function(r) c(y = unname(sim_signal(r, "c")),
                                         z = unname(sim_signal(r, "d"))))
  expect_equal(s$summary$y, c(3, 1))
  expect_equal(s$summary$z, c(2, 4))
  expect_equal(s$parameters$a, c(3, 1))
  expect_equal(s$parameters$b, c(2, 4))
  mc <- sim_monte_carlo(m, 3, function(i) list(), 7, c(0, 1), times = 1,
                        summary = function(r) c(y = unname(sim_signal(r, "c"))))
  expect_equal(nrow(mc$summary), 3L)
  expect_equal(mc$summary$run, 1:3)
  v <- data.frame(a = 1:2)
  v$b <- list(c(1, 2), c(3, 4))
  s <- sim_sweep(m, v, c(0, 1), times = 1)
  expect_equal(s$status, c("success", "success"))
  expect_equal(s$parameters$b[[2]], c(3, 4))
})

test_that("noise blocks never touch the caller's generator", {
  old <- RNGkind()
  on.exit(RNGkind(old[1], old[2], old[3]), add = TRUE)
  RNGkind(normal.kind = "Box-Muller")
  set.seed(44); invisible(rnorm(1)); expected <- rnorm(3)
  m <- sim_model()
  sim_add_block(m, "white_noise", "w", period = 0.1)
  p <- sim_compile(m)
  set.seed(44); invisible(rnorm(1))
  invisible(sim_run(p, c(0, 1)))
  expect_identical(rnorm(3), expected)
})

test_that("large noise seeds give distinct draws and streams", {
  for (seed in c(1, 1e9, 1e12, 2^53 - 1)) {
    m <- sim_model()
    sim_add_block(m, "white_noise", "w", seed = seed, period = 0.1)
    w <- sim_signal(sim_run(sim_compile(m), c(0, 1), times = seq(0, 1, 0.1)), "w")
    expect_length(unique(w), 11L)
  }
  m <- sim_model()
  expect_error(sim_add_block(m, "white_noise", "w", seed = 2^53), "2\\^53")
  x <- noise_normal(3, 0:4999, 1)
  expect_gt(stats::ks.test(x, "pnorm")$p.value, 1e-4)
  expect_lt(abs(cor(x[-1], x[-5000])), 0.06)
})

test_that("Monte Carlo warns that a Box-Muller cache cannot be restored", {
  old <- RNGkind()
  on.exit(RNGkind(old[1], old[2], old[3]), add = TRUE)
  RNGkind(normal.kind = "Box-Muller")
  m <- decay_model(1)
  expect_warning(sim_monte_carlo(m, 2, function(i) list(a = runif(1)), 7, c(0, 1), times = 1),
                 "Box-Muller")
  expect_equal(RNGkind()[2], "Box-Muller")
})

test_that("sweeps forward on_failure and keep columns unambiguous", {
  m <- sim_model(params = list(k = 1, run = 0))
  sim_add_block(m, "constant", "u", value = sim_param("k"))
  s <- sim_sweep(m, data.frame(k = 2), c(0, 1), times = 1, on_failure = "warn")
  expect_identical(s$status, "success")
  expect_error(sim_sweep(m, data.frame(k = 2, run = 4), c(0, 1), times = 1), "cannot be varied")
  expect_error(sim_sweep(m, data.frame(k = c(2, 3)), c(0, 1), times = 1,
                         summary = function(r) c(k = 1)), "also parameter")
})
