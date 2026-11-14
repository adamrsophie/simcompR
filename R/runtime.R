# Pure evaluation functions used by the simulator. None of them modifies
# anything outside its own frame: the solver may call them at trial points
# that are later rejected.

block_inputs <- function(b, sig) {
  lapply(b$in_idx, function(i) sig[i])
}

write_outputs <- function(b, y, sig) {
  if (!is.list(y)) y <- list(y)
  if (length(y) != length(b$out_idx)) {
    sim_abort("Block '", b$name, "' returned ", length(y), " output(s), expected ",
              length(b$out_idx), ".")
  }
  for (j in seq_along(y)) {
    v <- as.numeric(y[[j]])
    if (length(v) != length(b$out_idx[[j]])) {
      sim_abort("Block '", b$name, "' returned an output of length ", length(v),
                " on port ", j, "; expected ", length(b$out_idx[[j]]), ".")
    }
    sig[b$out_idx[[j]]] <- v
  }
  sig
}

eval_block <- function(b, t, x, d, sig) {
  write_outputs(b, b$def$output(t, x[b$x_idx], d[b$d_idx], block_inputs(b, sig)), sig)
}

# Evaluate all block outputs at time t.
#
# `ticking` is a logical vector over blocks: discrete blocks with a sample
# hit now compute their outputs; all other discrete blocks keep the held
# value. With refresh = TRUE, blocks with modes (switch, relay) recompute
# their mode from their current inputs before computing outputs.
#
# Returns list(sig, d, loop) where loop records algebraic loop residuals.
evaluate <- function(plan, t, x, d, held, ticking = NULL, refresh = FALSE) {
  sig <- held
  blocks <- plan$blocks
  loop_info <- NULL
  for (item in plan$schedule) {
    if (!item$loop) {
      b <- blocks[[item$blocks]]
      if (b$discrete && (is.null(ticking) || !ticking[item$blocks])) next
      if (refresh && b$has_modes) {
        d[b$d_idx] <- b$def$modes(t, x[b$x_idx], d[b$d_idx], block_inputs(b, sig))
      }
      if (length(b$out_idx)) sig <- eval_block(b, t, x, d, sig)
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
  z <- sig[zi]
  r <- F(z)
  nr <- max(abs(r))
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
    step <- tryCatch(qr.solve(J, -r), error = function(e) NULL)
    if (is.null(step)) {
      sim_abort("Algebraic loop {", paste(names(blocks), collapse = ", "),
                "} at t = ", format(t), ": singular Jacobian, the loop has no ",
                "unique solution.")
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
    sim_abort("Algebraic loop {", paste(names(blocks), collapse = ", "),
              "} at t = ", format(t), " did not converge: residual ",
              format(nr), " after ", it, " iteration(s).")
  }
  s <- sig
  s[zi] <- z
  list(sig = s, residual = nr, iterations = it)
}

derivatives <- function(plan, t, x, d, sig) {
  dx <- numeric(plan$nx)
  for (b in plan$blocks) {
    if (!length(b$x_idx)) next
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
