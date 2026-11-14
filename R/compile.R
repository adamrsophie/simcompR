#' Compile a model into an execution plan
#'
#' Compilation turns a mutable [sim_model()] into an immutable
#' `sim_plan` in these stages:
#'
#' 1. Flatten subsystems and resolve [sim_param()] references.
#' 2. Build every block and validate its parameters.
#' 3. Validate connections: blocks and ports exist and every input has
#'    exactly one source.
#' 4. Propagate signal sizes along connections and check compatibility.
#' 5. Build the direct-feedthrough dependency graph. There is an edge from
#'    block A to block B when an output of A feeds an input of B that has
#'    direct feedthrough. Inputs of integrators, unit delays and other
#'    state-output blocks create no edge, so feedback through them is valid.
#' 6. Find strongly connected components of that graph. A component with
#'    more than one block, or a block feeding itself, is an algebraic loop.
#' 7. Order the blocks topologically. Ties are broken by block name, so the
#'    order, the state vector layout and the results do not depend on the
#'    order in which blocks or connections were added.
#' 8. Lay out the signal, continuous state and discrete state vectors and
#'    collect sample times, zero-crossing functions and events.
#'
#' Algebraic loops are rejected by default. With `algebraic_loops = "solve"`
#' each loop is solved at every evaluation by Newton iterations on the loop
#' outputs, with a finite-difference Jacobian and backtracking; the run
#' fails with an error when the residual does not fall below `loop_tol`.
#'
#' @param model A [sim_model()].
#' @param algebraic_loops `"error"` (default) or `"solve"`.
#' @param loop_tol Absolute residual tolerance for loop solving.
#' @param loop_maxit Maximum Newton iterations per loop solve.
#' @return A `sim_plan` object. Use [sim_run()] to simulate it,
#'   [sim_set_params()] to change parameters and [sim_diagnostics()] to
#'   inspect the analysis.
#' @examples
#' m <- sim_model("feedback")
#' sim_add_block(m, "constant", "r", value = 1)
#' sim_add_block(m, "subtract", "e")
#' sim_add_block(m, "integrator", "x")
#' sim_connect(m, "r", "e", input = 1)
#' sim_connect(m, "x", "e", input = 2)
#' sim_connect(m, "e", "x")
#' plan <- sim_compile(m)
#' plan
#' sim_diagnostics(plan)$order
#'
#' # A loop of direct-feedthrough blocks is an algebraic loop.
#' m2 <- sim_model("loop")
#' sim_add_block(m2, "constant", "u", value = 1)
#' sim_add_block(m2, "sum", "s", signs = "+-")
#' sim_add_block(m2, "gain", "k", gain = 2)
#' sim_connect(m2, "u", "s", input = 1)
#' sim_connect(m2, "k", "s", input = 2)
#' sim_connect(m2, "s", "k")
#' try(sim_compile(m2))
#' plan2 <- sim_compile(m2, algebraic_loops = "solve")
#' sim_signal(sim_run(plan2, c(0, 1), times = 1), "s") # y = 1 - 2y, so 1/3
#' @export
sim_compile <- function(model, algebraic_loops = c("error", "solve"),
                        loop_tol = 1e-10, loop_maxit = 50) {
  check_model(model)
  algebraic_loops <- match.arg(algebraic_loops)
  check_number(loop_tol, "`loop_tol`")
  if (!is_count(loop_maxit)) sim_abort("`loop_maxit` must be a positive integer.")
  flat <- flatten_model(model)
  opts <- list(name = model$name, algebraic_loops = algebraic_loops,
               loop_tol = loop_tol, loop_maxit = as.integer(loop_maxit),
               model_version = model$version)
  compile_flat(flat, model$params, opts)
}

compile_flat <- function(flat, params, opts, previous = NULL) {
  bnames <- sort(names(flat$blocks), method = "radix")
  nb <- length(bnames)
  if (!nb) sim_abort("Model '", opts$name, "' has no blocks.")
  fb <- flat$blocks[bnames]

  # 1-2: resolve parameters and build definitions -----------------------
  defs <- lapply(fb, function(b) {
    sc <- b$scopes
    sc[[length(sc)]] <- params
    rp <- resolve_params(b$params, sc, b$name)
    def <- build_block(get_block_type(b$type), rp, b$name)
    def$type <- b$type
    def$resolved <- rp
    def
  })
  names(defs) <- bnames

  # 3: connections ---------------------------------------------------------
  cn <- flat$connections
  fi <- match(cn$from, bnames)
  ti <- match(cn$to, bnames)
  bad <- is.na(fi) | is.na(ti)
  if (any(bad)) {
    sim_abort("Connection refers to unknown block(s): ",
              paste(unique(c(cn$from[is.na(fi)], cn$to[is.na(ti)])), collapse = ", "), ".")
  }
  for (k in seq_len(nrow(cn))) {
    nout <- length(defs[[fi[k]]]$outputs)
    nin <- length(defs[[ti[k]]]$inputs)
    if (cn$from_port[k] > nout) {
      sim_abort("Block '", cn$from[k], "' has ", nout, " output port(s); ",
                "connection uses output ", cn$from_port[k], ".")
    }
    if (cn$to_port[k] > nin) {
      sim_abort("Block '", cn$to[k], "' has ", nin, " input port(s); ",
                "connection uses input ", cn$to_port[k], ".")
    }
  }
  key <- paste(ti, cn$to_port)
  if (anyDuplicated(key)) {
    k <- which(duplicated(key))[1]
    sim_abort("Input ", cn$to_port[k], " of block '", cn$to[k], "' has more ",
              "than one source.")
  }
  src <- lapply(defs, function(d) matrix(NA_integer_, length(d$inputs), 2L))
  for (k in seq_len(nrow(cn))) src[[ti[k]]][cn$to_port[k], ] <- c(fi[k], cn$from_port[k])
  unconnected <- unlist(lapply(seq_len(nb), function(i) {
    miss <- which(is.na(src[[i]][, 1]))
    if (length(miss)) paste0(bnames[i], ":", defs[[i]]$input_names[miss])
  }))
  if (length(unconnected)) {
    sim_abort("Unconnected input port(s): ", paste(unconnected, collapse = ", "), ".")
  }

  # 4: sizes -----------------------------------------------------------------
  sizes <- propagate_sizes(defs, src, bnames)

  # 5-7: dependency graph and schedule ---------------------------------------
  adj <- vector("list", nb)
  for (i in seq_len(nb)) adj[[i]] <- integer()
  for (k in seq_len(nrow(cn))) {
    if (defs[[ti[k]]]$feedthrough[cn$to_port[k]]) adj[[fi[k]]] <- c(adj[[fi[k]]], ti[k])
  }
  adj <- lapply(adj, function(a) sort(unique(a)))
  signature <- list(sizes = sizes,
                    ft = lapply(defs, `[[`, "feedthrough"),
                    nx = vapply(defs, `[[`, 0L, "nx"),
                    nd = vapply(defs, `[[`, 0L, "nd"),
                    sample = lapply(defs, `[[`, "sample"),
                    nzc = vapply(defs, `[[`, 0L, "nzc"),
                    conns = cn)
  if (!is.null(previous) && identical(previous$signature, signature)) {
    comps <- previous$graph$components
    loops <- previous$graph$loops
    structure_reused <- TRUE
  } else {
    comps <- condensed_order(nb, adj, tarjan_scc(nb, adj))
    loops <- Filter(function(cmp) length(cmp) > 1L || cmp[1] %in% adj[[cmp[1]]], comps)
    structure_reused <- FALSE
  }
  if (length(loops) && opts$algebraic_loops == "error") {
    msgs <- vapply(loops, function(cmp) {
      cyc <- find_cycle(cmp, adj)
      paste0("  ", paste(bnames[cyc], collapse = " -> "))
    }, "")
    sim_abort("Algebraic loop(s) detected: every block on the cycle has direct ",
              "feedthrough, so the outputs depend on themselves at the same time.\n",
              paste(msgs, collapse = "\n"),
              "\nBreak the loop with a stateful block (integrator, unit_delay, ...) ",
              "or compile with algebraic_loops = \"solve\".")
  }

  # 8: layouts -----------------------------------------------------------------
  layout <- layout_vectors(defs, sizes, bnames)
  blocks <- lapply(seq_len(nb), function(i) {
    d <- defs[[i]]
    in_idx <- lapply(seq_along(d$inputs), function(j) {
      layout$out_idx[[src[[i]][j, 1]]][[src[[i]][j, 2]]]
    })
    list(name = bnames[i], type = d$type, def = d,
         in_idx = in_idx, out_idx = layout$out_idx[[i]],
         x_idx = layout$x_idx[[i]], d_idx = layout$d_idx[[i]],
         discrete = !is.null(d$sample), has_modes = is.function(d$modes),
         nzc = d$nzc)
  })
  names(blocks) <- bnames
  signal_of <- layout$signal_of
  for (a in names(flat$aliases)) {
    if (is.null(signal_of[[a]])) signal_of[[a]] <- signal_of[[flat$aliases[[a]]]]
  }
  for (b in blocks) {
    if (b$type == "scope" && is.null(signal_of[[b$name]])) {
      signal_of[[b$name]] <- unlist(b$in_idx)
    }
  }

  schedule <- lapply(comps, function(cmp) {
    if (length(cmp) == 1L && !(cmp[1] %in% adj[[cmp[1]]])) {
      return(list(loop = FALSE, blocks = cmp))
    }
    list(loop = TRUE, blocks = cmp,
         z_idx = unlist(lapply(cmp, function(i) unlist(blocks[[i]]$out_idx))))
  })

  held0 <- numeric(layout$nsig)
  for (b in blocks) {
    if (!b$discrete || !length(b$out_idx)) next
    h <- if (is.function(b$def$held0)) b$def$held0(sizes$out[[b$name]]) else 0
    if (!is.list(h)) h <- list(h)
    for (j in seq_along(b$out_idx)) {
      v <- as.numeric(h[[min(j, length(h))]])
      n <- length(b$out_idx[[j]])
      if (length(v) == 1L) v <- rep(v, n)
      if (length(v) != n) {
        sim_abort("Block '", b$name, "': initial held output has length ",
                  length(v), ", expected ", n, ".")
      }
      held0[b$out_idx[[j]]] <- v
    }
  }

  rates <- sample_rates(blocks)
  zc <- layout_zero_crossings(blocks)
  events <- compile_events(flat$events, blocks, layout, params, zc$n)

  diag <- list(
    order = bnames[unlist(lapply(schedule, `[[`, "blocks"))],
    loops = lapply(loops, function(cmp) bnames[cmp]),
    edges = data.frame(
      from = rep(bnames, lengths(adj)),
      to = bnames[unlist(adj)],
      stringsAsFactors = FALSE),
    rates = rates$table,
    notes = c(sizes$notes, rate_notes(blocks, src, bnames))
  )

  structure(list(
    name = opts$name, opts = opts, flat = flat, params = params,
    blocks = blocks, schedule = schedule,
    signal_names = layout$signal_names, signal_of = signal_of,
    nsig = layout$nsig, nx = layout$nx, nd = layout$nd,
    x0 = layout$x0, d0 = layout$d0, held0 = held0,
    state_names = layout$state_names, dstate_names = layout$dstate_names,
    rates = rates, zc = zc, events = events,
    sizes = sizes, signature = signature,
    graph = list(components = comps, loops = loops),
    structure_reused = structure_reused,
    diagnostics = diag
  ), class = "sim_plan")
}

propagate_sizes <- function(defs, src, bnames) {
  nb <- length(defs)
  ins <- lapply(defs, `[[`, "inputs")
  outs <- lapply(defs, `[[`, "outputs")
  partial <- FALSE
  notes <- character()
  repeat {
    changed <- FALSE
    for (i in seq_len(nb)) {
      for (j in seq_along(ins[[i]])) {
        s <- outs[[src[[i]][j, 1]]][src[[i]][j, 2]]
        if (is.na(s)) next
        if (is.na(ins[[i]][j])) {
          ins[[i]][j] <- s
          changed <- TRUE
        } else if (ins[[i]][j] != s) {
          sim_abort("Signal size mismatch: '", bnames[src[[i]][j, 1]], "' output ",
                    src[[i]][j, 2], " has size ", s, " but input '",
                    defs[[i]]$input_names[j], "' of block '", bnames[i],
                    "' expects size ", ins[[i]][j], ".")
        }
      }
      if (anyNA(outs[[i]]) && (!anyNA(ins[[i]]) || (partial && any(!is.na(ins[[i]]))))) {
        f <- defs[[i]]$size_out
        if (!is.function(f)) {
          sim_abort("Block '", bnames[i], "' has outputs of unknown size and no `size_out`.")
        }
        if (anyNA(ins[[i]])) {
          # partial pass inside a feedback loop: accept only a definite answer
          new <- tryCatch(suppressWarnings(as.integer(f(ins[[i]]))), error = function(e) NA)
          if (anyNA(new) || length(new) != length(outs[[i]])) next
        } else {
          new <- withCallingHandlers(as.integer(f(ins[[i]])), error = function(e) {
            sim_abort("Block '", bnames[i], "' (", defs[[i]]$type, "): ", conditionMessage(e))
          })
        }
        if (length(new) != length(outs[[i]]) || anyNA(new) || any(new < 1L)) {
          sim_abort("Block '", bnames[i], "': could not determine output sizes.")
        }
        known <- !is.na(outs[[i]])
        if (any(outs[[i]][known] != new[known])) {
          sim_abort("Block '", bnames[i], "': declared and computed output sizes differ.")
        }
        outs[[i]] <- new
        changed <- TRUE
      }
    }
    if (!changed) {
      if (!partial) {
        partial <- TRUE
        next
      }
      # Last resort inside a feedback loop with no size information at all
      # (for example y = cos(y)): take one unknown input as scalar and say so.
      stuck <- which(vapply(ins, anyNA, TRUE))
      if (!length(stuck)) break
      i <- stuck[1]
      ins[[i]][is.na(ins[[i]])] <- 1L
      notes <- c(notes, paste0("Input size of block '", bnames[i], "' could not be ",
                               "inferred from the diagram and was taken as 1 (scalar)."))
      partial <- FALSE
    } else {
      partial <- FALSE
    }
  }
  unresolved <- bnames[vapply(seq_len(nb), function(i) anyNA(ins[[i]]) || anyNA(outs[[i]]), TRUE)]
  if (length(unresolved)) {
    sim_abort("Could not determine signal sizes for block(s): ",
              paste(unresolved, collapse = ", "),
              ". A feedback loop needs at least one block with a known output size.")
  }
  for (i in seq_len(nb)) {
    f <- defs[[i]]$check_sizes
    if (is.function(f)) {
      withCallingHandlers(f(ins[[i]], outs[[i]]), error = function(e) {
        sim_abort("Block '", bnames[i], "' (", defs[[i]]$type, "): ", conditionMessage(e))
      })
    }
  }
  names(ins) <- names(outs) <- bnames
  list(inputs = ins, out = outs, notes = notes)
}

layout_vectors <- function(defs, sizes, bnames) {
  out_idx <- x_idx <- d_idx <- vector("list", length(defs))
  sig_names <- character(); signal_of <- list()
  state_names <- dstate_names <- character()
  x0 <- d0 <- numeric()
  pos <- 0L
  for (i in seq_along(defs)) {
    d <- defs[[i]]
    so <- sizes$out[[i]]
    out_idx[[i]] <- lapply(seq_along(so), function(j) pos + sum(so[seq_len(j - 1L)]) + seq_len(so[j]))
    for (j in seq_along(so)) {
      nm <- if (length(so) == 1L) bnames[i] else paste0(bnames[i], ".", d$output_names[j])
      sig_names <- c(sig_names, element_names(nm, so[j]))
      signal_of[[nm]] <- out_idx[[i]][[j]]
    }
    pos <- pos + sum(so)
    x_idx[[i]] <- length(x0) + seq_len(d$nx)
    if (d$nx) state_names <- c(state_names, element_names(bnames[i], d$nx))
    x0 <- c(x0, d$x0)
    d_idx[[i]] <- length(d0) + seq_len(d$nd)
    if (d$nd) dstate_names <- c(dstate_names, element_names(bnames[i], d$nd))
    d0 <- c(d0, d$d0)
  }
  list(out_idx = out_idx, x_idx = x_idx, d_idx = d_idx, nsig = pos,
       signal_names = sig_names, signal_of = signal_of,
       state_names = state_names, dstate_names = dstate_names,
       x0 = x0, d0 = d0, nx = length(x0), nd = length(d0))
}

sample_rates <- function(blocks) {
  disc <- Filter(function(b) b$discrete, blocks)
  if (!length(disc)) {
    return(list(table = data.frame(period = numeric(), offset = numeric(),
                                   blocks = character()),
                period = numeric(), offset = numeric(), block_rate = integer()))
  }
  st <- t(vapply(disc, function(b) b$def$sample, numeric(2)))
  key <- unique(st)
  key <- key[order(key[, 1], key[, 2]), , drop = FALSE]
  rate_of <- vapply(disc, function(b) {
    which(key[, 1] == b$def$sample[1] & key[, 2] == b$def$sample[2])
  }, integer(1))
  block_rate <- integer(length(blocks))
  block_rate[match(names(disc), names(blocks))] <- rate_of
  tab <- data.frame(period = key[, 1], offset = key[, 2],
                    blocks = vapply(seq_len(nrow(key)), function(r) {
                      paste(names(disc)[rate_of == r], collapse = ", ")
                    }, ""), stringsAsFactors = FALSE)
  list(table = tab, period = key[, 1], offset = key[, 2], block_rate = block_rate)
}

layout_zero_crossings <- function(blocks) {
  pos <- 0L
  idx <- vector("list", length(blocks))
  dirs <- character()
  for (i in seq_along(blocks)) {
    n <- blocks[[i]]$nzc
    idx[[i]] <- pos + seq_len(n)
    pos <- pos + n
    if (n) dirs <- c(dirs, blocks[[i]]$def$zc_direction)
  }
  list(block_idx = idx, n = pos, directions = dirs,
       blocks = which(vapply(blocks, function(b) b$nzc > 0L, TRUE)))
}

compile_events <- function(events, blocks, layout, params, zc_offset) {
  if (!length(events)) return(list())
  bnames <- names(blocks)
  out <- lapply(events, function(e) {
    tr <- e$trigger
    if (tr$kind == "zero_crossing") {
      b <- blocks[[tr$signal]]
      if (is.null(b)) sim_abort("Event '", e$name, "': unknown block '", tr$signal, "'.")
      if (tr$port > length(b$out_idx)) {
        sim_abort("Event '", e$name, "': block '", tr$signal, "' has no output ", tr$port, ".")
      }
      idx <- b$out_idx[[tr$port]]
      if (tr$element > length(idx)) {
        sim_abort("Event '", e$name, "': output of '", tr$signal, "' has only ",
                  length(idx), " element(s).")
      }
      tr$sig_index <- idx[tr$element]
    }
    vars <- scope_values(e$scopes, params)
    resets <- lapply(names(e$reset), function(target) {
      b <- blocks[[target]]
      if (is.null(b)) sim_abort("Event '", e$name, "': unknown reset target '", target, "'.")
      kind <- if (length(b$x_idx)) "x" else if (length(b$d_idx)) "d" else
        sim_abort("Event '", e$name, "': block '", target, "' has no state to reset.")
      r <- e$reset[[target]]
      fun <- if (is.function(r)) r else {
        expr <- parse_safe(r, paste0("Event '", e$name, "' reset of '", target, "'"))
        function(x, s, t) eval_safe(expr, c(vars, list(x = x, s = s, t = t)))
      }
      list(target = target, kind = kind,
           idx = if (kind == "x") b$x_idx else b$d_idx, fun = fun)
    })
    list(name = e$name, trigger = tr, resets = resets, terminal = e$terminal,
         priority = e$priority)
  })
  ord <- order(vapply(out, `[[`, 0, "priority"), vapply(out, `[[`, "", "name"),
               method = "radix")
  out <- out[ord]
  # zero-crossing components of events follow those of blocks
  k <- zc_offset
  for (i in seq_along(out)) {
    if (out[[i]]$trigger$kind == "zero_crossing") {
      k <- k + 1L
      out[[i]]$zc_index <- k
    }
  }
  names(out) <- vapply(out, `[[`, "", "name")
  out
}

# Parameter values visible to string expressions of an event: every name of
# every enclosing scope, resolved through parameter references.
scope_values <- function(scopes, params) {
  if (is.null(scopes) || !length(scopes)) scopes <- list(params)
  scopes[[1]] <- params
  inner_first <- rev(scopes)
  nms <- unique(unlist(lapply(scopes, names)))
  vals <- lapply(nms, function(n) {
    tryCatch(resolve_value(sim_param(n), inner_first, 1L, "event"),
             error = function(e) NULL)
  })
  names(vals) <- nms
  Filter(Negate(is.null), vals)
}

rate_notes <- function(blocks, src, bnames) {
  notes <- character()
  for (i in seq_along(blocks)) {
    b <- blocks[[i]]
    if (!b$discrete) next
    for (j in seq_len(nrow(src[[i]]))) {
      s <- blocks[[src[[i]][j, 1]]]
      if (s$discrete && !identical(s$def$sample, b$def$sample)) {
        notes <- c(notes, paste0(
          "'", b$name, "' (period ", b$def$sample[1], ") reads '", s$name,
          "' (period ", s$def$sample[1], "): the value is the one held at the ",
          "reader's sample hit (zero-order hold across rates)."))
      }
    }
  }
  notes
}

#' Change parameters of a compiled plan
#'
#' Returns a new plan with updated model parameters. Blocks are rebuilt
#' with the new values. When block sizes, feedthrough, state dimensions,
#' sample times and connections are unchanged the dependency analysis is
#' reused (`plan$structure_reused` is `TRUE`); otherwise the graph is
#' analyzed again, so a change that alters the structure (for example making
#' the `D` matrix of a state-space block nonzero) is handled correctly.
#'
#' @param plan A `sim_plan`.
#' @param ... Named parameter values, or a single named list.
#' @return A new `sim_plan`.
#' @examples
#' m <- sim_model("decay", params = list(a = 1))
#' sim_add_block(m, "integrator", "x", initial = 1)
#' sim_add_block(m, "gain", "k", gain = sim_param("a"))
#' sim_connect(m, "x", "k")
#' sim_add_block(m, "gain", "neg", gain = -1)
#' sim_connect(m, "k", "neg")
#' sim_connect(m, "neg", "x")
#' p1 <- sim_compile(m)
#' p2 <- sim_set_params(p1, a = 3)
#' p2$structure_reused
#' sim_states(sim_run(p2, c(0, 1), times = 1))
#' @export
sim_set_params <- function(plan, ...) {
  check_plan(plan)
  new <- list(...)
  if (length(new) == 1L && is.null(names(new)) && is.list(new[[1]])) new <- new[[1]]
  check_param_list(new, "Parameters")
  unknown <- setdiff(names(new), names(plan$params))
  if (length(unknown)) {
    sim_abort("Unknown model parameter(s): ", paste(unknown, collapse = ", "),
              ". Model parameters: ", paste(names(plan$params), collapse = ", "), ".")
  }
  params <- plan$params
  params[names(new)] <- new
  compile_flat(plan$flat, params, plan$opts, previous = plan)
}

#' @rdname sim_set_params
#' @export
sim_param_values <- function(plan) {
  check_plan(plan)
  plan$params
}

check_plan <- function(plan) {
  if (!inherits(plan, "sim_plan")) {
    sim_abort("Expected a compiled plan; call sim_compile() first.")
  }
  invisible(plan)
}

#' Compilation diagnostics
#'
#' @param x A `sim_plan` or a `sim_result`.
#' @return For a plan, a list with the execution `order`, algebraic `loops`,
#'   direct-feedthrough `edges`, sample `rates` and `notes` about rate
#'   transitions. For a result, the solver status, statistics, warnings and
#'   algebraic loop residuals.
#' @examples
#' m <- sim_model()
#' sim_add_block(m, "sine", "u")
#' sim_add_block(m, "zero_order_hold", "zoh", period = 0.1)
#' sim_add_block(m, "unit_delay", "z", period = 0.2)
#' sim_connect(m, "u", "zoh")
#' sim_connect(m, "zoh", "z")
#' sim_diagnostics(sim_compile(m))
#' @export
sim_diagnostics <- function(x) {
  if (inherits(x, "sim_plan")) return(x$diagnostics)
  if (inherits(x, "sim_result")) return(x$diagnostics)
  sim_abort("`x` must be a sim_plan or a sim_result.")
}

#' @export
print.sim_plan <- function(x, ...) {
  cat("<sim_plan '", x$name, "'> ", length(x$blocks), " block(s), ",
      x$nx, " continuous state(s), ", x$nd, " discrete state(s), ",
      x$nsig, " signal element(s)\n", sep = "")
  cat("Execution order: ", paste(x$diagnostics$order, collapse = " > "), "\n", sep = "")
  if (nrow(x$rates$table)) {
    cat("Sample rates:\n")
    for (r in seq_len(nrow(x$rates$table))) {
      cat("  period ", x$rates$table$period[r], ", offset ", x$rates$table$offset[r],
          ": ", x$rates$table$blocks[r], "\n", sep = "")
    }
  }
  if (length(x$graph$loops)) {
    cat("Algebraic loops (solved by Newton iteration):\n")
    for (l in x$diagnostics$loops) cat("  {", paste(l, collapse = ", "), "}\n", sep = "")
  }
  if (x$zc$n || length(x$events)) {
    cat("Zero-crossing functions: ", x$zc$n, "; events: ", length(x$events), "\n", sep = "")
  }
  invisible(x)
}
