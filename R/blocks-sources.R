# Source blocks -------------------------------------------------------------

check_same_length <- function(..., what) {
  lens <- lengths(list(...))
  if (length(unique(lens)) != 1L) {
    sim_abort(what, " must have equal lengths (no recycling).")
  }
  lens[[1]]
}

register_source_blocks <- function() {
  reg <- function(...) sim_register_block_type(..., overwrite = TRUE)

  reg("constant", category = "source",
      description = "Constant signal",
      params = list(value = 1),
      build = function(p) {
        v <- check_numeric(p$value, "`value`")
        list(outputs = length(v), output = function(t, x, d, u) v)
      })

  reg("step", category = "source",
      description = "Step from `before` to `after` at `time`",
      params = list(time = 1, before = 0, after = 1),
      build = function(p) {
        t1 <- check_number(p$time, "`time`")
        b <- check_numeric(p$before, "`before`")
        a <- check_numeric(p$after, "`after`")
        n <- check_same_length(b, a, what = "`before` and `after`")
        list(outputs = n,
             output = function(t, x, d, u) if (t < t1) b else a,
             breakpoints = function(t0, tf) t1)
      })

  reg("ramp", category = "source",
      description = "Ramp with `slope` starting at `start`",
      params = list(slope = 1, start = 0, initial = 0),
      build = function(p) {
        s <- check_numeric(p$slope, "`slope`")
        t1 <- check_number(p$start, "`start`")
        y0 <- check_numeric(p$initial, "`initial`")
        n <- check_same_length(s, y0, what = "`slope` and `initial`")
        list(outputs = n,
             output = function(t, x, d, u) y0 + s * max(0, t - t1),
             breakpoints = function(t0, tf) t1)
      })

  reg("sine", category = "source",
      description = "bias + amplitude * sin(frequency * t + phase), frequency in rad/s",
      params = list(amplitude = 1, frequency = 1, phase = 0, bias = 0),
      build = function(p) {
        a <- check_numeric(p$amplitude, "`amplitude`")
        w <- check_numeric(p$frequency, "`frequency`")
        ph <- check_numeric(p$phase, "`phase`")
        b <- check_numeric(p$bias, "`bias`")
        n <- check_same_length(a, w, ph, b,
                               what = "`amplitude`, `frequency`, `phase` and `bias`")
        list(outputs = n,
             output = function(t, x, d, u) b + a * sin(w * t + ph))
      })

  reg("time", category = "source",
      description = "Simulation time",
      build = function(p) list(outputs = 1L, output = function(t, x, d, u) t))

  reg("piecewise", category = "source",
      description = "Piecewise-constant signal: values[i] from times[i] on",
      params = list(times = NULL, values = NULL),
      build = function(p) {
        tt <- check_numeric(p$times, "`times`")
        if (is.unsorted(tt, strictly = TRUE)) {
          sim_abort("`times` must be strictly increasing.")
        }
        vals <- p$values
        if (is.matrix(vals)) {
          vals <- check_matrix(vals, "`values`", nrow = length(tt))
        } else {
          vals <- matrix(check_numeric(vals, "`values`"), ncol = 1L)
          if (nrow(vals) != length(tt)) {
            sim_abort("`values` must have one entry (or matrix row) per time.")
          }
        }
        list(outputs = ncol(vals),
             output = function(t, x, d, u) {
               i <- findInterval(t, tt)
               vals[max(1L, i), ]
             },
             breakpoints = function(t0, tf) tt)
      })

  reg("pulse", category = "source",
      description = "Periodic rectangular pulse train",
      params = list(amplitude = 1, period = 1, width = 0.5, delay = 0),
      build = function(p) {
        a <- check_number(p$amplitude, "`amplitude`")
        per <- check_number(p$period, "`period`")
        w <- check_number(p$width, "`width`")
        del <- check_number(p$delay, "`delay`")
        if (per <= 0) sim_abort("`period` must be positive.")
        if (w <= 0 || w >= 1) sim_abort("`width` must be in (0, 1), a fraction of the period.")
        list(outputs = 1L,
             output = function(t, x, d, u) {
               if (t < del) return(0)
               # compare with the breakpoint times exactly as they are computed
               k <- floor((t - del) / per)
               while (del + (k + 1) * per <= t) k <- k + 1
               while (k > 0 && del + k * per > t) k <- k - 1
               if (t < del + (k + w) * per) a else 0
             },
             breakpoints = function(t0, tf) {
               k <- seq.int(max(0, floor((t0 - del) / per)),
                            ceiling((tf - del) / per))
               c(del + k * per, del + (k + w) * per)
             })
      })

  reg("from_data", category = "source",
      description = "Recorded signal interpolated in time",
      params = list(times = NULL, values = NULL, method = "linear"),
      build = function(p) {
        tt <- check_numeric(p$times, "`times`")
        if (is.unsorted(tt, strictly = TRUE)) {
          sim_abort("`times` must be strictly increasing.")
        }
        vals <- if (is.matrix(p$values)) {
          check_matrix(p$values, "`values`", nrow = length(tt))
        } else {
          matrix(check_numeric(p$values, "`values`"), ncol = 1L)
        }
        if (nrow(vals) != length(tt)) {
          sim_abort("`values` must have one entry (or matrix row) per time.")
        }
        method <- match.arg(p$method, c("linear", "constant"))
        list(outputs = ncol(vals),
             output = function(t, x, d, u) {
               vapply(seq_len(ncol(vals)), function(j) {
                 stats::approx(tt, vals[, j], xout = t, method = method,
                               rule = 2, f = 0, ties = "ordered")$y
               }, numeric(1))
             },
             breakpoints = function(t0, tf) if (method == "constant") tt else numeric())
      })

  reg("signal_function", category = "source", portable = FALSE,
      description = "Signal given by an R function of time (not portable)",
      params = list(fun = NULL, size = 1),
      build = function(p) {
        f <- p$fun
        if (!is.function(f)) sim_abort("`fun` must be a function of time.")
        n <- p$size
        if (!is_count(n)) sim_abort("`size` must be a positive integer.")
        list(outputs = as.integer(n),
             output = function(t, x, d, u) {
               y <- as.numeric(f(t))
               if (length(y) != n) {
                 sim_abort("`fun` returned ", length(y), " values, expected ", n, ".")
               }
               y
             })
      })

  reg("white_noise", category = "source",
      description = "Band-limited Gaussian noise, held between sample hits",
      params = list(mean = 0, sd = 1, period = 0.1, offset = 0, seed = 1),
      build = function(p) {
        mu <- check_numeric(p$mean, "`mean`")
        s <- check_numeric(p$sd, "`sd`")
        n <- check_same_length(mu, s, what = "`mean` and `sd`")
        if (any(s < 0)) sim_abort("`sd` must be non-negative.")
        seed <- p$seed
        if (!is_scalar_number(seed) || seed != round(seed)) {
          sim_abort("`seed` must be an integer.")
        }
        # The discrete state is the sample counter. Each draw depends only on
        # (seed, counter), so results are reproducible and independent of the
        # caller's random stream.
        # ponytail: one set.seed() per sample hit; precompute blocks of draws
        # if noise-heavy models become slow.
        draw <- function(k) {
          with_seed(seed * 1000003 + k, stats::rnorm(n, mu, s))
        }
        list(outputs = n, nd = 1L, d0 = 0, sample = c(p$period, p$offset),
             output = function(t, x, d, u) draw(d[[1]]),
             update = function(t, x, d, u) d + 1,
             held0 = function(size) mu)
      })
}

# Run expr with a temporary seed and restore the caller's random stream.
# Evaluate expr with a fixed generator configuration (Mersenne-Twister,
# Inversion, Rejection), so that results do not depend on the caller's
# RNGkind(), then restore the caller's kinds and stream.
with_seed <- function(seed, expr) {
  genv <- globalenv()
  had <- exists(".Random.seed", envir = genv, inherits = FALSE)
  if (had) old <- get(".Random.seed", envir = genv, inherits = FALSE)
  kinds <- RNGkind()
  on.exit({
    suppressWarnings(RNGkind(kinds[1], kinds[2], kinds[3]))
    if (had) assign(".Random.seed", old, envir = genv)
    else if (exists(".Random.seed", envir = genv, inherits = FALSE)) {
      rm(".Random.seed", envir = genv)
    }
  })
  set.seed(seed %% .Machine$integer.max, kind = "Mersenne-Twister",
           normal.kind = "Inversion", sample.kind = "Rejection")
  expr
}
