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
#' Every solve starts from `loop_guess` where given and from zero
#' otherwise, so results never depend on the history of solver trial
#' points. A loop with several solutions (for example `y = 4 / y`) is
#' steered to one of them by its starting value.
#'
#' @param model A [sim_model()].
#' @param algebraic_loops `"error"` (default) or `"solve"`.
#' @param loop_tol Absolute residual tolerance for loop solving.
#' @param loop_maxit Maximum Newton iterations per loop solve.
#' @param loop_guess Optional named list of starting values for loop
#'   solving, keyed by signal name (as in [sim_signal()]); each value has
#'   the length of that signal.
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
                        loop_tol = 1e-10, loop_maxit = 50, loop_guess = NULL) {
  check_model(model)
  check_acyclic(model)
  algebraic_loops <- match.arg(algebraic_loops)
  check_number(loop_tol, "`loop_tol`")
  if (!is_count(loop_maxit)) sim_abort("`loop_maxit` must be a positive integer.")
  if (!is.null(loop_guess) && (!is.list(loop_guess) || is.null(names(loop_guess)))) {
    sim_abort("`loop_guess` must be a named list of numeric values.")
  }
  flat <- flatten_model(model)
  opts <- list(name = model$name, algebraic_loops = algebraic_loops,
               loop_tol = loop_tol, loop_maxit = as.integer(loop_maxit),
               loop_guess = loop_guess, model_version = model$version)
  compile_flat(flat, model$params, opts)
}

compile_flat <- function(flat, params, opts, previous = NULL) {
  bnames <- sort(names(flat$blocks), method = "radix")
  nb <- length(bnames)
  if (!nb) sim_abort("Model '", opts$name, "' has no blocks.")
  fb <- flat$blocks[bnames]

  # 1-2: resolve parameters and build definitions -----------------------
  # A plan keeps the block types it was compiled with, so registering or
  # replacing a type later never changes what sim_set_params() rebuilds.
  if (is.null(flat$type_defs)) {
    types <- unique(vapply(fb, `[[`, "", "type"))
    flat$type_defs <- stats::setNames(lapply(types, get_block_type), types)
  }
  defs <- lapply(fb, function(b) {
    sc <- b$scopes
    sc[[length(sc)]] <- params
    rp <- resolve_params(b$params, sc, b$name)
    tdef <- flat$type_defs[[b$type]]
    if (isTRUE(tdef$scoped)) rp$.scope <- scope_values(rev(sc), params)
    def <- build_block(tdef, rp, b$name)
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
  defs <- lapply(seq_len(nb), function(i) add_kinks(defs[[i]], sizes$inputs[[i]], bnames[i]))

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
         sizes_out = sizes$out[[i]],
         discrete = !is.null(d$sample), has_modes = is.function(d$modes),
         nzc = d$nzc)
  })
  names(blocks) <- bnames
  signal_of <- layout$signal_of
  for (a in names(flat$aliases)) {
    if (!is.null(signal_of[[a]])) {
      sim_abort("Signal name '", a, "' is used twice: it names an output of subsystem ",
                "instance '", sub("[.][^.]*$", "", a), "' and another signal; rename one of them.")
    }
    signal_of[[a]] <- signal_of[[flat$aliases[[a]]]]
  }
  for (b in blocks) {
    if (b$type != "scope") next
    if (!is.null(signal_of[[b$name]])) {
      sim_abort("Signal name '", b$name, "' is used twice: it names scope '", b$name,
                "' and another signal; rename one of them.")
    }
    signal_of[[b$name]] <- unlist(b$in_idx)
  }

  schedule <- lapply(comps, function(cmp) {
    if (length(cmp) == 1L && !(cmp[1] %in% adj[[cmp[1]]])) {
      b <- blocks[[cmp]]
      return(list(loop = FALSE, blocks = cmp, discrete = b$discrete,
                  has_bp = is.function(b$def$breakpoints), skip = !length(b$out_idx)))
    }
    list(loop = TRUE, blocks = cmp,
         z_idx = unlist(lapply(cmp, function(i) unlist(blocks[[i]]$out_idx))))
  })

  guess <- rep(NA_real_, layout$nsig)
  for (nm in names(opts$loop_guess)) {
    idx <- signal_of[[nm]]
    if (is.null(idx)) sim_abort("`loop_guess`: unknown signal '", nm, "'.")
    v <- check_numeric(opts$loop_guess[[nm]], paste0("`loop_guess$", nm, "`"))
    if (length(v) != length(idx)) {
      sim_abort("`loop_guess$", nm, "` must have length ", length(idx), ".")
    }
    guess[idx] <- v
  }

  held0 <- numeric(layout$nsig)
  for (b in blocks) held0 <- set_held0(b, layout$d0[b$d_idx], held0)

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
    blocks = blocks, schedule = schedule, loop_guess = guess,
    continuous = which(vapply(blocks, function(b) length(b$x_idx) > 0L, TRUE)),
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
  declared <- lapply(outs, function(o) !is.na(o))
  partial <- FALSE
  notes <- character()
  # visit blocks in the order sizes flow along connections, so a chain is
  # resolved in one pass whatever the block names
  dep <- vector("list", nb)
  for (i in seq_len(nb)) for (j in seq_along(ins[[i]])) {
    k <- src[[i]][j, 1]
    dep[[k]] <- c(dep[[k]], i)
  }
  dep <- lapply(dep, function(v) sort(unique(as.integer(v))))
  visit <- unlist(condensed_order(nb, dep, tarjan_scc(nb, dep)), use.names = FALSE)
  repeat {
    changed <- FALSE
    for (i in visit) {
      for (j in seq_along(ins[[i]])) {
        s <- outs[[src[[i]][j, 1]]][src[[i]][j, 2]]
        if (is.na(s)) {
          # a consumer that needs a known size fixes its producer's port
          if (!is.na(ins[[i]][j])) {
            outs[[src[[i]][j, 1]]][src[[i]][j, 2]] <- ins[[i]][j]
            changed <- TRUE
          }
          next
        }
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
  # A partial pass inside a loop may have guessed an output size before all
  # input sizes were known. Check every inferred size against the final
  # input sizes; report a clear error instead of failing at run time.
  for (i in seq_len(nb)) {
    f <- defs[[i]]$size_out
    # declared sizes are not inferred and are checked at run time instead
    if (!is.function(f) || all(declared[[i]])) next
    final <- tryCatch(suppressWarnings(as.integer(f(ins[[i]]))), error = function(e) NULL)
    inf <- !declared[[i]]
    if (!is.null(final) && length(final) == length(outs[[i]]) &&
        any(final[inf] != outs[[i]][inf])) {
      sim_abort("Block '", bnames[i], "': output size inferred inside a feedback loop (",
                paste(outs[[i]], collapse = ", "), ") does not match its final input sizes (",
                paste(final, collapse = ", "), "). Declare the sizes, for example with a ",
                "constant of the right length or a block parameter.")
    }
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
      if (!is.null(signal_of[[nm]])) {
        sim_abort("Signal name '", nm, "' is used twice (a block name collides with ",
                  "the name of another block's output port); rename one of them.")
      }
      en <- element_names(nm, so[j])
      dup <- intersect(en, sig_names)
      if (length(dup)) {
        sim_abort("Signal element name '", dup[1], "' is used twice; rename one of the blocks.")
      }
      sig_names <- c(sig_names, en)
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

# Kinks are points where a memoryless block's output, or its slope, jumps
# (the corners of a saturation, the switching point of a compare). A block
# declares them with `kinks(t, x, d, u)` and `kink_size(in_sizes)`; they
# become extra crossing functions whose only effect is to restart the
# integrator exactly at the kink, so the right-hand side is smooth inside
# every step. They are optional: methods without root finding ignore them.
add_kinks <- function(d, in_sizes, name) {
  d$zc_optional <- rep(FALSE, d$nzc)
  if (!is.function(d$kinks)) return(d)
  nk <- as.integer(d$kink_size(in_sizes))
  if (length(nk) != 1L || is.na(nk) || nk < 0L) {
    sim_abort("Block '", name, "': `kink_size` must return one count.")
  }
  if (!nk) return(d)
  nz <- d$nzc
  zc0 <- d$zc
  on0 <- d$on_zc
  kf <- d$kinks
  d$zc <- function(t, x, dd, u) {
    k <- as.numeric(kf(t, x, dd, u))
    if (length(k) != nk) {
      sim_abort("Block '", name, "': `kinks` returned ", length(k), " value(s), expected ", nk, ".")
    }
    c(if (nz) zc0(t, x, dd, u), k)
  }
  d$on_zc <- function(i, direction, t, x, dd, u) {
    if (i > nz) list(silent = TRUE) else on0(i, direction, t, x, dd, u)
  }
  d$zc_direction <- c(if (nz) d$zc_direction, rep("either", nk))
  d$zc_optional <- c(rep(FALSE, nz), rep(TRUE, nk))
  d$nzc <- nz + nk
  if (isTRUE(d$kink_side) && !nz && nk == in_sizes[[1]]) d <- add_kink_sides(d, nk, kf)
  d
}

# A memoryless block whose output jumps at its kinks (compare, sign, floor,
# ceiling, round) keeps, per element, the side of the kink the input is
# heading to: set by the direction of a located kink (or of a departure
# from it) and cleared once the input is off the kink. On the kink the
# output is taken just past it on that side, so the jump happens at the
# instant itself, where switches, stop blocks and events see it.
add_kink_sides <- function(d, nk, kf) {
  nd0 <- d$nd
  out0 <- d$output
  band <- function(v) 1e-8 * (1 + abs(v))
  sides <- function(dd) dd[nd0 + seq_len(nk)]
  d$nd <- nd0 + nk
  d$d0 <- c(d$d0, numeric(nk))
  d$output <- function(t, x, dd, u) {
    s <- sides(dd)
    if (any(s != 0)) {
      v <- u[[1]]
      on <- s != 0 & abs(kf(t, x, dd, u)) <= band(v)
      v[on] <- v[on] + s[on] * 10 * band(v[on])
      u[[1]] <- v
    }
    out0(t, x, dd[seq_len(nd0)], u)
  }
  d$on_zc <- function(i, direction, t, x, dd, u) {
    v <- u[[1]]
    h <- 10 * band(v[i])
    up <- u; up[[1]][i] <- v[i] + h
    dn <- u; dn[[1]][i] <- v[i] - h
    # direction of the input from the direction of the kink function
    slope <- sign(kf(t, x, dd, up)[i] - kf(t, x, dd, dn)[i])
    dirs <- c(rising = 1, falling = -1, either = 0)
    out <- dd
    out[nd0 + i] <- dirs[[direction]] * slope
    list(silent = TRUE, d = out)
  }
  d$modes <- function(t, x, dd, u) {
    s <- sides(dd)
    s[abs(kf(t, x, dd, u)) > band(u[[1]])] <- 0
    c(dd[seq_len(nd0)], s)
  }
  d
}

layout_zero_crossings <- function(blocks) {
  pos <- 0L
  idx <- vector("list", length(blocks))
  dirs <- character()
  optional <- logical()
  for (i in seq_along(blocks)) {
    n <- blocks[[i]]$nzc
    idx[[i]] <- pos + seq_len(n)
    pos <- pos + n
    if (n) {
      dirs <- c(dirs, blocks[[i]]$def$zc_direction)
      optional <- c(optional, blocks[[i]]$def$zc_optional)
    }
  }
  list(block_idx = idx, n = pos, directions = dirs, optional = optional,
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
         priority = e$priority, prefix = e$prefix %||% "")
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
  fx <- exact_match_args(sys.call(), "plan", parent.frame())
  if (is.null(fx)) {
    new <- list(...)
  } else {
    if (is.null(fx$fixed$plan)) sim_abort("sim_set_params() needs `plan`.")
    plan <- fx$fixed$plan
    new <- fx$dots
  }
  check_plan(plan)
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

# Held outputs of a discrete block before its first sample hit. `held0` may
# take the initial discrete state as a second argument, so that overriding
# the state in sim_run(initial = ) also changes what the block shows.
set_held0 <- function(b, d0, held) {
  if (!b$discrete || !length(b$out_idx)) return(held)
  h <- if (!is.function(b$def$held0)) 0 else
    if (length(formals(b$def$held0)) >= 2L) b$def$held0(b$sizes_out, d0) else
      b$def$held0(b$sizes_out)
  if (!is.list(h)) h <- list(h)
  for (j in seq_along(b$out_idx)) {
    v <- as.numeric(h[[min(j, length(h))]])
    n <- length(b$out_idx[[j]])
    if (length(v) == 1L) v <- rep(v, n)
    if (length(v) != n) {
      sim_abort("Block '", b$name, "': initial held output has length ",
                length(v), ", expected ", n, ".")
    }
    if (!all(is.finite(v))) {
      sim_abort("Block '", b$name, "': initial held output ", j,
                " must be finite, not ", deparse_short(v), ".")
    }
    held[b$out_idx[[j]]] <- v
  }
  held
}
