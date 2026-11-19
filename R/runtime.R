# Pure evaluation functions used by the simulator. None of them modifies
# anything outside its own frame: the solver may call them at trial points
# that are later rejected.

block_inputs <- function(b, sig) {
  n <- length(b$in_idx)
  if (n == 1L) return(list(sig[b$in_idx[[1L]]]))
  if (!n) return(list())
  lapply(b$in_idx, function(i) sig[i])
}

write_outputs <- function(b, y, sig, strict = FALSE) {
  # fast path: one port, a plain numeric vector of the right length
  if (!is.list(y) && length(b$out_idx) == 1L && is.double(y) && is.null(dim(y))) {
    idx <- b$out_idx[[1L]]
    if (length(y) == length(idx) && !(strict && !all(is.finite(y)))) {
      sig[idx] <- y
      return(sig)
    }
  }
  if (!is.list(y)) y <- list(y)
  if (length(y) != length(b$out_idx)) {
    sim_abort("Block '", b$name, "' returned ", length(y), " output(s), expected ",
              length(b$out_idx), ".")
  }
  for (j in seq_along(y)) {
    v <- y[[j]]
    if (!is.numeric(v) || (length(dim(v)) > 1L && sum(dim(v) > 1L) > 1L)) {
      sim_abort("Block '", b$name, "' returned a ", class(v)[1], " on port ", j,
                "; outputs must be numeric vectors.")
    }
    v <- as.numeric(v)
    if (strict && !all(is.finite(v))) {
      sim_abort("Block '", b$name, "' returned ", format(v[!is.finite(v)][1]), " on port ", j,
                "; signals must be finite.")
    }
    if (length(v) != length(b$out_idx[[j]])) {
      sim_abort("Block '", b$name, "' returned an output of length ", length(v),
                " on port ", j, "; expected ", length(b$out_idx[[j]]), ".")
    }
    sig[b$out_idx[[j]]] <- v
  }
  sig
}

eval_block <- function(b, t, x, d, sig) {
  write_outputs(b, b$def$output(t, x[b$x_idx], d[b$d_idx], block_inputs(b, sig)), sig,
                strict = TRUE)
}

# Evaluate all block outputs at time t.
#
# `ticking` is a logical vector over blocks: discrete blocks with a sample
# hit now compute their outputs; all other discrete blocks keep the held
# value. With refresh = TRUE, blocks with modes (switch, relay) recompute
# their mode from their current inputs before computing outputs.
#
# Returns list(sig, d, loop) where loop records algebraic loop residuals.
evaluate <- function(plan, t, x, d, held, ticking = NULL, refresh = FALSE,
                     t_left = Inf) {
  if (refresh && any(vapply(plan$blocks, `[[`, TRUE, "has_modes"))) {
    # Modes are recomputed from fully evaluated signals, never from a
    # partially filled signal vector, and repeated until they settle, so the
    # result does not depend on the execution order of mode blocks.
    for (iter in seq_len(50L)) {
      ev <- evaluate(plan, t, x, d, held, ticking, refresh = FALSE)
      dn <- d
      for (b in plan$blocks) {
        if (b$has_modes) {
          dn[b$d_idx] <- b$def$modes(t, x[b$x_idx], d[b$d_idx], block_inputs(b, ev$sig))
        }
      }
      if (identical(dn, d)) return(ev)
      d <- dn
    }
    sim_abort("Block modes did not settle at t = ", format(t),
              " (switching blocks keep changing each other).")
  }
  sig <- held
  blocks <- plan$blocks
  loop_info <- NULL
  for (item in plan$schedule) {
    if (!item$loop) {
      if (item$skip) next
      if (item$discrete && (is.null(ticking) || !ticking[item$blocks])) next
      b <- blocks[[item$blocks]]
      # sources with breakpoints show their left limit at a segment's end
      tb <- if (item$has_bp && t > t_left) t_left else t
      sig <- eval_block(b, tb, x, d, sig)
      next
    }
    res <- solve_loop(plan, item, t, x, d, sig, ticking)
    sig <- res$sig
    loop_info <- rbind(loop_info, c(residual = res$residual, iterations = res$iterations))
  }
  list(sig = sig, d = d, loop = loop_info)
}

# Newton iteration on the outputs z of one algebraic loop: z = F(z).
solve_loop <- function(plan, item, t, x, d, sig, ticking) {
  blocks <- plan$blocks[item$blocks]
  zi <- item$z_idx
  active <- vapply(item$blocks, function(i) {
    b <- plan$blocks[[i]]
    !b$discrete || (!is.null(ticking) && ticking[i])
  }, TRUE)
  # outputs of sampled blocks between their hits are held, not unknowns
  zi <- unlist(lapply(blocks[active], function(b) unlist(b$out_idx)))
  if (!length(zi)) return(list(sig = sig, residual = 0, iterations = 0L))
  F <- function(z) {
    s <- sig
    s[zi] <- z
    out <- s
    for (k in seq_along(blocks)) {
      if (active[k]) out <- write_outputs(blocks[[k]], blocks[[k]]$def$output(
        t, x[blocks[[k]]$x_idx], d[blocks[[k]]$d_idx], block_inputs(blocks[[k]], s)), out)
    }
    out[zi] - z
  }
  tol <- plan$opts$loop_tol
  loop_name <- paste0("Algebraic loop {", paste(names(blocks), collapse = ", "), "} at t = ",
                      format(t))
  z <- sig[zi]
  g <- plan$loop_guess[zi]
  if (length(g)) z[!is.na(g)] <- g[!is.na(g)]
  r <- F(z)
  nr <- max(abs(r))
  if (!is.finite(nr)) {
    sim_abort(loop_name, ": the residual is not finite (NaN or infinite) at the ",
              "starting point; give a starting value with sim_compile(loop_guess = ).")
  }
  it <- 0L
  while (nr > tol && it < plan$opts$loop_maxit) {
    it <- it + 1L
    n <- length(z)
    J <- matrix(0, n, n)
    for (j in seq_len(n)) {
      h <- 1e-7 * max(1, abs(z[j]))
      zp <- z; zp[j] <- zp[j] + h
      J[, j] <- (F(zp) - r) / h
    }
    if (!all(is.finite(J))) {
      sim_abort(loop_name, ": the residual is not finite near the current iterate; ",
                "give a starting value with sim_compile(loop_guess = ).")
    }
    step <- tryCatch(qr.solve(J, -r), error = function(e) NULL)
    if (is.null(step)) {
      sim_abort(loop_name, ": singular Jacobian at the current iterate. The loop may ",
                "have no unique solution, or the start is at a stationary point; ",
                "give a starting value with sim_compile(loop_guess = ).")
    }
    lambda <- 1
    repeat {
      zn <- z + lambda * step
      rn <- F(zn)
      if (all(is.finite(rn)) && max(abs(rn)) < nr) break
      lambda <- lambda / 2
      if (lambda < 1e-8) break
    }
    z <- zn; r <- rn; nr <- max(abs(r))
  }
  if (!is.finite(nr) || nr > tol) {
    sim_abort(loop_name, " did not converge: residual ",
              format(nr), " after ", it, " iteration(s).")
  }
  s <- sig
  s[zi] <- z
  list(sig = s, residual = nr, iterations = it)
}

derivatives <- function(plan, t, x, d, sig) {
  dx <- numeric(plan$nx)
  for (i in plan$continuous) {
    b <- plan$blocks[[i]]
    v <- as.numeric(b$def$deriv(t, x[b$x_idx], d[b$d_idx], block_inputs(b, sig)))
    if (length(v) != length(b$x_idx)) {
      sim_abort("Block '", b$name, "' returned ", length(v), " derivative(s), expected ",
                length(b$x_idx), ".")
    }
    dx[b$x_idx] <- v
  }
  dx
}

# Zero-crossing functions of blocks followed by those of events.
zero_crossings <- function(plan, t, x, d, sig) {
  g <- numeric(plan$zc$n + sum(vapply(plan$events, function(e) {
    e$trigger$kind == "zero_crossing"
  }, TRUE)))
  for (i in plan$zc$blocks) {
    b <- plan$blocks[[i]]
    g[plan$zc$block_idx[[i]]] <- b$def$zc(t, x[b$x_idx], d[b$d_idx], block_inputs(b, sig))
  }
  for (e in plan$events) {
    if (e$trigger$kind == "zero_crossing") {
      g[e$zc_index] <- sig[e$trigger$sig_index] - e$trigger$level
    }
  }
  g
}

signal_list <- function(plan, sig) {
  lapply(plan$signal_of, function(i) sig[i])
}
