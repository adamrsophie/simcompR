# Continuous and discrete dynamic blocks -------------------------------------

# Controllable canonical realization of num(s) / den(s) (or of z for discrete
# systems), with coefficients in descending powers.
tf_to_ss <- function(num, den) {
  num <- check_numeric(num, "`num`")
  den <- check_numeric(den, "`den`")
  while (length(den) > 1L && den[1] == 0) den <- den[-1]
  while (length(num) > 1L && num[1] == 0) num <- num[-1]
  if (den[1] == 0) sim_abort("`den` must not be identically zero.")
  n <- length(den) - 1L
  if (length(num) - 1L > n) {
    sim_abort("the transfer function is improper: deg(num) > deg(den).")
  }
  a <- den / den[1]
  b <- c(rep(0, n + 1L - length(num)), num) / den[1]
  D <- matrix(b[1], 1L, 1L)
  if (n == 0L) {
    return(list(A = matrix(0, 0, 0), B = matrix(0, 0, 1), C = matrix(0, 1, 0), D = D))
  }
  A <- matrix(0, n, n)
  A[1, ] <- -a[-1]
  if (n > 1L) A[cbind(2:n, 1:(n - 1L))] <- 1
  B <- matrix(c(1, rep(0, n - 1L)), n, 1L)
  C <- matrix(b[-1] - a[-1] * b[1], 1L, n)
  list(A = A, B = B, C = C, D = D)
}

check_ss <- function(p) {
  A <- check_matrix(p$A, "`A`")
  n <- nrow(A)
  if (ncol(A) != n) sim_abort("`A` must be square.")
  B <- check_matrix(p$B, "`B`", nrow = n)
  C <- check_matrix(p$C, "`C`", ncol = n)
  D <- if (is.null(p$D)) matrix(0, nrow(C), ncol(B)) else
    check_matrix(p$D, "`D`", nrow = nrow(C), ncol = ncol(B))
  list(A = A, B = B, C = C, D = D)
}

initial_state <- function(initial, n) {
  x0 <- check_numeric(initial, "`initial`")
  if (length(x0) == 1L && n > 1L && x0 == 0) return(rep(0, n))
  if (length(x0) != n) sim_abort("`initial` must have length ", n, ".")
  x0
}

ss_continuous <- function(s, x0) {
  n <- nrow(s$A)
  ft <- any(s$D != 0)
  out <- list(inputs = ncol(s$B), outputs = nrow(s$C), feedthrough = ft,
              nx = n, x0 = x0)
  if (n == 0L) {
    out$output <- function(t, x, d, u) as.numeric(s$D %*% u[[1]])
    return(out)
  }
  out$output <- if (ft) {
    function(t, x, d, u) as.numeric(s$C %*% x + s$D %*% u[[1]])
  } else {
    function(t, x, d, u) as.numeric(s$C %*% x)
  }
  out$deriv <- function(t, x, d, u) as.numeric(s$A %*% x + s$B %*% u[[1]])
  out
}

ss_discrete <- function(s, x0, period, offset) {
  ft <- any(s$D != 0)
  list(inputs = ncol(s$B), outputs = nrow(s$C), feedthrough = ft,
       nd = nrow(s$A), d0 = x0, sample = c(period, offset),
       output = function(t, x, d, u) as.numeric(s$C %*% d + s$D %*% u[[1]]),
       update = function(t, x, d, u) as.numeric(s$A %*% d + s$B %*% u[[1]]),
       held0 = function(size) {
         # before the first hit the block holds C x0 (input taken as zero)
         as.numeric(s$C %*% x0)
       })
}

register_dynamic_blocks <- function() {
  reg <- function(...) sim_register_block_type(..., overwrite = TRUE)

  reg("integrator", category = "continuous",
      description = "Continuous integrator, dx/dt = u, y = x",
      params = list(initial = 0),
      build = function(p) {
        x0 <- check_numeric(p$initial, "`initial`")
        list(inputs = length(x0), outputs = length(x0), feedthrough = FALSE,
             nx = length(x0), x0 = x0,
             output = function(t, x, d, u) x,
             deriv = function(t, x, d, u) u[[1]])
      })

  reg("transfer_function", category = "continuous",
      description = "SISO transfer function num(s)/den(s)",
      params = list(num = 1, den = c(1, 1), initial = 0),
      build = function(p) {
        s <- tf_to_ss(p$num, p$den)
        ss_continuous(s, initial_state(p$initial, nrow(s$A)))
      })

  reg("state_space", category = "continuous",
      description = "dx/dt = A x + B u, y = C x + D u",
      params = list(A = NULL, B = NULL, C = NULL, D = NULL, initial = 0),
      build = function(p) {
        s <- check_ss(p)
        ss_continuous(s, initial_state(p$initial, nrow(s$A)))
      })

  reg("ode", category = "continuous", portable = FALSE,
      description = "General dx/dt = f(t, x, u), y = g(t, x, u) from R functions (not portable)",
      params = list(deriv = NULL, output = NULL, initial = NULL,
                    feedthrough = FALSE, size = NA),
      build = function(p) {
        f <- p$deriv
        g <- p$output %||% function(t, x, u) x
        if (!is.function(f) || !is.function(g)) {
          sim_abort("`deriv` and `output` must be functions(t, x, u).")
        }
        x0 <- check_numeric(p$initial, "`initial`")
        ft <- isTRUE(p$feedthrough)
        size <- p$size
        list(inputs = NA, outputs = if (is.na(size)) NA else as.integer(size),
             feedthrough = ft, nx = length(x0), x0 = x0,
             size_out = function(n) length(g(0, x0, numeric(n))),
             output = function(t, x, d, u) as.numeric(g(t, x, u[[1]])),
             deriv = function(t, x, d, u) {
               dx <- as.numeric(f(t, x, u[[1]]))
               if (length(dx) != length(x0)) {
                 sim_abort("`deriv` returned ", length(dx), " values, expected ",
                           length(x0), ".")
               }
               dx
             })
      })

  reg("unit_delay", category = "discrete",
      description = "y[k] = u[k-1]",
      params = list(initial = 0, period = 1, offset = 0),
      build = function(p) {
        x0 <- check_numeric(p$initial, "`initial`")
        list(inputs = length(x0), outputs = length(x0), feedthrough = FALSE,
             nd = length(x0), d0 = x0, sample = c(p$period, p$offset),
             output = function(t, x, d, u) d,
             update = function(t, x, d, u) u[[1]],
             held0 = function(size) x0)
      })

  reg("zero_order_hold", category = "discrete",
      description = "Sample the input at each hit and hold it",
      params = list(period = 1, offset = 0, initial = 0),
      build = function(p) {
        y0 <- check_numeric(p$initial, "`initial`")
        list(inputs = NA, outputs = NA, feedthrough = TRUE,
             sample = c(p$period, p$offset),
             size_out = function(n) n,
             check_sizes = function(n_in, n_out) {
               if (length(y0) != 1L && length(y0) != n_in[[1]]) {
                 sim_abort("`initial` must be a scalar or match the input size.")
               }
             },
             output = function(t, x, d, u) u[[1]],
             held0 = function(size) expand(y0, size[[1]]))
      })

  reg("discrete_integrator", category = "discrete",
      description = "Accumulator with forward, backward or trapezoidal rule",
      params = list(gain = 1, initial = 0, period = 1, offset = 0,
                    method = "forward"),
      build = function(p) {
        k <- check_number(p$gain, "`gain`")
        x0 <- check_numeric(p$initial, "`initial`")
        method <- match.arg(p$method, c("forward", "backward", "trapezoidal"))
        h <- check_number(p$period, "`period`")
        kt <- k * h
        n <- length(x0)
        # Discrete state: x[k] (accumulator). Simulink-compatible forms:
        # forward:     y = x,               x+ = x + K T u
        # backward:    y = x + K T u,       x+ = y
        # trapezoidal: y = x + K T/2 u,     x+ = y + K T/2 u
        out <- switch(method,
          forward = list(ft = FALSE, y = function(d, v) d,
                         up = function(d, v) d + kt * v),
          backward = list(ft = TRUE, y = function(d, v) d + kt * v,
                          up = function(d, v) d + kt * v),
          trapezoidal = list(ft = TRUE, y = function(d, v) d + kt / 2 * v,
                             up = function(d, v) d + kt * v))
        list(inputs = n, outputs = n, feedthrough = out$ft, nd = n, d0 = x0,
             sample = c(p$period, p$offset),
             output = function(t, x, d, u) out$y(d, u[[1]]),
             update = function(t, x, d, u) out$up(d, u[[1]]),
             held0 = function(size) x0)
      })

  reg("discrete_state_space", category = "discrete",
      description = "x[k+1] = A x[k] + B u[k], y[k] = C x[k] + D u[k]",
      params = list(A = NULL, B = NULL, C = NULL, D = NULL, initial = 0,
                    period = 1, offset = 0),
      build = function(p) {
        s <- check_ss(p)
        ss_discrete(s, initial_state(p$initial, nrow(s$A)), p$period, p$offset)
      })

  reg("discrete_transfer_function", category = "discrete",
      description = "SISO num(z)/den(z), descending powers of z",
      params = list(num = 1, den = c(1, -0.5), initial = 0, period = 1, offset = 0),
      build = function(p) {
        s <- tf_to_ss(p$num, p$den)
        n <- nrow(s$A)
        if (n == 0L) {
          s$A <- matrix(0, 1, 1); s$B <- matrix(0, 1, 1); s$C <- matrix(0, 1, 1)
          n <- 1L
        }
        ss_discrete(s, initial_state(p$initial, n), p$period, p$offset)
      })
}
