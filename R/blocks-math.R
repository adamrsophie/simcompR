# Mathematical blocks ------------------------------------------------------

# Combine sizes of inputs that are combined elementwise. Scalars expand to
# the common size (an explicit, documented rule); any other mismatch is an
# error, so R recycling never happens silently.
broadcast_size <- function(sizes) {
  if (all(is.na(sizes))) return(NA_integer_)
  sizes <- sizes[!is.na(sizes)]
  big <- unique(sizes[sizes != 1L])
  if (length(big) > 1L) {
    sim_abort("input sizes ", paste(sizes, collapse = ", "),
              " are incompatible; inputs must have equal sizes or be scalars.")
  }
  if (length(big)) big else 1L
}

expand <- function(v, n) if (length(v) == n) v else rep(v, n)

elementwise_block <- function(fun) {
  list(inputs = NA, outputs = NA, feedthrough = TRUE,
       size_out = function(n) n,
       output = function(t, x, d, u) fun(u[[1]]))
}

math_functions <- list(
  exp = exp, log = log, log10 = log10, sqrt = sqrt,
  square = function(v) v^2, reciprocal = function(v) 1 / v,
  sin = sin, cos = cos, tan = tan, asin = asin, acos = acos, atan = atan,
  sinh = sinh, cosh = cosh, tanh = tanh, abs = abs, sign = sign,
  floor = floor, ceiling = ceiling, round = round
)

parse_signs <- function(signs, allowed, what) {
  if (is.character(signs) && length(signs) == 1L) {
    signs <- strsplit(signs, "", fixed = TRUE)[[1]]
  }
  if (!is.character(signs) || !length(signs) || !all(signs %in% allowed)) {
    sim_abort(what, " must use only ", paste0("'", allowed, "'", collapse = " and "), ".")
  }
  signs
}

register_math_blocks <- function() {
  reg <- function(...) sim_register_block_type(..., overwrite = TRUE)

  reg("gain", category = "math",
      description = "Scalar, elementwise or matrix gain",
      params = list(gain = 1),
      build = function(p) {
        k <- p$gain
        if (is.matrix(k)) {
          k <- check_matrix(k, "`gain`")
          return(list(inputs = ncol(k), outputs = nrow(k), feedthrough = TRUE,
                      output = function(t, x, d, u) as.numeric(k %*% u[[1]])))
        }
        k <- check_numeric(k, "`gain`")
        if (length(k) == 1L) {
          return(list(inputs = NA, outputs = NA, feedthrough = TRUE,
                      size_out = function(n) n,
                      output = function(t, x, d, u) k * u[[1]]))
        }
        list(inputs = length(k), outputs = length(k), feedthrough = TRUE,
             output = function(t, x, d, u) k * u[[1]])
      })

  sum_build <- function(p) {
    s <- parse_signs(p$signs, c("+", "-"), "`signs`")
    w <- ifelse(s == "+", 1, -1)
    list(inputs = rep(NA, length(s)), outputs = NA, feedthrough = TRUE,
         size_out = function(n) broadcast_size(n),
         output = function(t, x, d, u) {
           n <- max(lengths(u))
           y <- numeric(n)
           for (i in seq_along(u)) y <- y + w[[i]] * expand(u[[i]], n)
           y
         })
  }
  reg("sum", category = "math", description = "Signed sum of inputs",
      params = list(signs = "++"), build = sum_build)
  reg("subtract", category = "math", description = "in1 - in2",
      build = function(p) sum_build(list(signs = "+-")))

  prod_build <- function(p) {
    ops <- parse_signs(p$ops, c("*", "/"), "`ops`")
    list(inputs = rep(NA, length(ops)), outputs = NA, feedthrough = TRUE,
         size_out = function(n) broadcast_size(n),
         output = function(t, x, d, u) {
           n <- max(lengths(u))
           y <- rep(1, n)
           for (i in seq_along(u)) {
             v <- expand(u[[i]], n)
             y <- if (ops[[i]] == "*") y * v else y / v
           }
           y
         })
  }
  reg("product", category = "math", description = "Elementwise product and quotient",
      params = list(ops = "**"), build = prod_build)
  reg("divide", category = "math", description = "in1 / in2",
      build = function(p) prod_build(list(ops = "*/")))

  reg("saturation", category = "math",
      description = "Clamp the input to [lower, upper]",
      params = list(lower = -1, upper = 1),
      build = function(p) {
        lo <- check_number(p$lower, "`lower`", finite = FALSE)
        hi <- check_number(p$upper, "`upper`", finite = FALSE)
        if (lo > hi) sim_abort("`lower` must not exceed `upper`.")
        elementwise_block(function(v) pmin(pmax(v, lo), hi))
      })

  reg("dead_zone", category = "math",
      description = "Zero inside [lower, upper], shifted input outside",
      params = list(lower = -0.5, upper = 0.5),
      build = function(p) {
        lo <- check_number(p$lower, "`lower`")
        hi <- check_number(p$upper, "`upper`")
        if (lo > hi) sim_abort("`lower` must not exceed `upper`.")
        elementwise_block(function(v) ifelse(v > hi, v - hi, ifelse(v < lo, v - lo, 0)))
      })

  reg("abs", category = "math", description = "Absolute value",
      build = function(p) elementwise_block(abs))

  reg("math", category = "math",
      description = paste("Elementwise function:", paste(names(math_functions), collapse = ", ")),
      params = list(fun = "exp"),
      build = function(p) {
        if (!is.character(p$fun) || length(p$fun) != 1L ||
            !p$fun %in% names(math_functions)) {
          sim_abort("`fun` must be one of ", paste(names(math_functions), collapse = ", "), ".")
        }
        elementwise_block(math_functions[[p$fun]])
      })

  reg("minmax", category = "math",
      description = "Minimum or maximum across inputs (or of one vector input)",
      params = list(fun = "min", n_inputs = 1),
      build = function(p) {
        fun <- match.arg(p$fun, c("min", "max"))
        if (!is_count(p$n_inputs)) sim_abort("`n_inputs` must be a positive integer.")
        n <- as.integer(p$n_inputs)
        if (n == 1L) {
          f <- if (fun == "min") min else max
          return(list(inputs = NA, outputs = 1L, feedthrough = TRUE,
                      output = function(t, x, d, u) f(u[[1]])))
        }
        f <- if (fun == "min") pmin else pmax
        list(inputs = rep(NA, n), outputs = NA, feedthrough = TRUE,
             size_out = function(s) broadcast_size(s),
             output = function(t, x, d, u) do.call(f, u))
      })

  reg("mux", category = "routing", description = "Concatenate inputs into one vector",
      params = list(n_inputs = 2),
      build = function(p) {
        if (!is_count(p$n_inputs)) sim_abort("`n_inputs` must be a positive integer.")
        list(inputs = rep(NA, p$n_inputs), outputs = NA, feedthrough = TRUE,
             size_out = function(s) sum(s),
             output = function(t, x, d, u) unlist(u, use.names = FALSE))
      })

  reg("demux", category = "routing", description = "Split a vector into pieces of given sizes",
      params = list(sizes = c(1, 1)),
      build = function(p) {
        s <- p$sizes
        if (!is.numeric(s) || !length(s) || any(s < 1) || any(s != round(s))) {
          sim_abort("`sizes` must be positive integers.")
        }
        s <- as.integer(s)
        ends <- cumsum(s)
        starts <- ends - s + 1L
        list(inputs = sum(s), outputs = s, feedthrough = TRUE,
             output = function(t, x, d, u) {
               lapply(seq_along(s), function(i) u[[1]][starts[i]:ends[i]])
             })
      })

  reg("selector", category = "routing", description = "Select elements by index",
      params = list(index = 1),
      build = function(p) {
        idx <- p$index
        if (!is.numeric(idx) || !length(idx) || any(idx < 1) || any(idx != round(idx))) {
          sim_abort("`index` must be positive integers.")
        }
        idx <- as.integer(idx)
        list(inputs = NA, outputs = length(idx), feedthrough = TRUE,
             check_sizes = function(n_in, n_out) {
               if (max(idx) > n_in[[1]]) {
                 sim_abort("`index` refers to element ", max(idx),
                           " of a signal of size ", n_in[[1]], ".")
               }
             },
             output = function(t, x, d, u) u[[1]][idx])
      })

  reg("lookup", category = "math",
      description = "One-dimensional table lookup with linear interpolation",
      params = list(x = NULL, y = NULL, extrapolate = "hold"),
      build = function(p) {
        xs <- check_numeric(p$x, "`x`")
        ys <- check_numeric(p$y, "`y`")
        check_same_length(xs, ys, what = "`x` and `y`")
        if (length(xs) < 2L) sim_abort("`x` and `y` need at least two points.")
        if (is.unsorted(xs, strictly = TRUE)) sim_abort("`x` must be strictly increasing.")
        ext <- match.arg(p$extrapolate, c("hold", "linear"))
        n <- length(xs)
        f <- function(v) {
          y <- stats::approx(xs, ys, xout = v, rule = 2, ties = "ordered")$y
          if (ext == "linear" && n >= 2L) {
            lo <- v < xs[1]
            hi <- v > xs[n]
            y[lo] <- ys[1] + (v[lo] - xs[1]) * (ys[2] - ys[1]) / (xs[2] - xs[1])
            y[hi] <- ys[n] + (v[hi] - xs[n]) * (ys[n] - ys[n - 1]) / (xs[n] - xs[n - 1])
          }
          y
        }
        elementwise_block(f)
      })

  reg("expression", category = "math",
      description = "Portable expression in u (input vector) and t",
      params = list(expr = NULL, size = NA),
      build = function(p) {
        e <- parse_safe(p$expr, "`expr`")
        f <- function(t, v) as.numeric(eval_safe(e, list(u = v, t = t)))
        size <- p$size
        list(inputs = NA, outputs = if (is.na(size)) NA else as.integer(size),
             feedthrough = TRUE,
             size_out = function(n) length(f(0, numeric(n))),
             output = function(t, x, d, u) f(t, u[[1]]))
      })

  reg("fcn", category = "math", portable = FALSE,
      description = "User R function of the input vector and time (not portable)",
      params = list(fun = NULL, size = NA),
      build = function(p) {
        f <- p$fun
        if (!is.function(f)) sim_abort("`fun` must be a function(u, t).")
        size <- p$size
        list(inputs = NA, outputs = if (is.na(size)) NA else as.integer(size),
             feedthrough = TRUE,
             size_out = function(n) length(f(numeric(n), 0)),
             output = function(t, x, d, u) as.numeric(f(u[[1]], t)))
      })
}
