#' Hierarchical subsystems
#'
#' A subsystem is an ordinary [sim_model()] whose external interface is
#' declared with `"inport"` and `"outport"` blocks. Port numbers come from
#' their `index` parameter. The subsystem is placed into a parent model
#' with `sim_add_block(parent, "subsystem", name, system = sub, params = ...)`.
#'
#' At compile time every subsystem instance is flattened into the parent:
#' inner block `b` of instance `s` becomes block `s/b`. Each instance gets
#' its own copy of the inner blocks, so instances never share state. A
#' [sim_param()] reference inside the subsystem is resolved first in the
#' instance's `params`, then in the parent's scope. The flattened model is
#' mathematically identical to building the same blocks by hand.
#'
#' `sim_subsystem()` is a convenience wrapper around [sim_model()].
#'
#' @param name Subsystem name.
#' @param params Default parameters of the subsystem. Instance `params`
#'   override them.
#' @return A `sim_model`.
#' @examples
#' # A first-order lag y' = (u - y) / tau, with inv_tau = 1 / tau.
#' lag <- sim_subsystem("lag", params = list(inv_tau = 1))
#' sim_add_block(lag, "inport", "u", index = 1)
#' sim_add_block(lag, "subtract", "e")
#' sim_add_block(lag, "gain", "k", gain = sim_param("inv_tau"))
#' sim_add_block(lag, "integrator", "y", initial = 0)
#' sim_add_block(lag, "outport", "out", index = 1)
#' sim_connect(lag, "u", "e", input = 1)
#' sim_connect(lag, "y", "e", input = 2)
#' sim_connect(lag, "e", "k")
#' sim_connect(lag, "k", "y")
#' sim_connect(lag, "y", "out")
#'
#' m <- sim_model("two_lags")
#' sim_add_block(m, "step", "u", time = 0)
#' sim_add_block(m, "subsystem", "fast", system = lag, params = list(inv_tau = 10))
#' sim_add_block(m, "subsystem", "slow", system = lag, params = list(inv_tau = 1))
#' sim_connect(m, "u", "fast")
#' sim_connect(m, "u", "slow")
#' r <- sim_run(sim_compile(m), tspan = c(0, 3))
#' tail(sim_signal(r, c("fast", "slow")), 1)
#' @export
sim_subsystem <- function(name, params = list()) {
  sim_model(name, params)
}

subsystem_ports <- function(system, kind) {
  is_kind <- vapply(system$blocks, function(b) b$type == kind, logical(1))
  blocks <- system$blocks[is_kind]
  if (!length(blocks)) return(character())
  idx <- vapply(blocks, function(b) {
    i <- b$params$index %||% 1
    if (!is_count(i)) sim_abort("Subsystem '", system$name, "': ", kind, " '",
                                b$name, "' needs a positive integer `index`.")
    as.integer(i)
  }, integer(1))
  if (anyDuplicated(idx) || !setequal(idx, seq_along(idx))) {
    sim_abort("Subsystem '", system$name, "': ", kind, " indices must be 1..",
              length(idx), " without gaps or duplicates.")
  }
  names(blocks)[order(idx)]
}

register_subsystem_blocks <- function() {
  reg <- function(...) sim_register_block_type(..., overwrite = TRUE)
  # a declared size seeds size inference, also inside feedback loops
  pass <- function(p) list(inputs = p$size, outputs = p$size, feedthrough = TRUE,
                           size_out = function(n) n,
                           check_sizes = function(n_in, n_out) {
                             if (!is.na(p$size) && n_in[[1]] != p$size) {
                               sim_abort("signal size ", n_in[[1]],
                                         " does not match the declared size ", p$size, ".")
                             }
                           },
                           output = function(t, x, d, u) u[[1]])
  reg("inport", category = "subsystem", description = "Subsystem input port",
      params = list(index = 1, size = NA), build = pass)
  reg("outport", category = "subsystem", description = "Subsystem output port",
      params = list(index = 1, size = NA), build = pass)
  reg("subsystem", category = "subsystem",
      description = "Instance of a hierarchical subsystem model",
      params = list(system = NULL, params = list()),
      build = function(p) {
        sys <- p$system
        if (!inherits(sys, "sim_model")) sim_abort("`system` must be a sim_model.")
        nin <- length(subsystem_ports(sys, "inport"))
        nout <- length(subsystem_ports(sys, "outport"))
        list(inputs = rep(NA, nin), outputs = rep(NA, nout),
             input_names = subsystem_ports(sys, "inport"),
             output_names = subsystem_ports(sys, "outport"),
             feedthrough = rep(TRUE, nin), subsystem = TRUE,
             output = function(t, x, d, u) {
               sim_abort("subsystem blocks are flattened before simulation.")
             })
      })
}

# Names a subsystem may read from its scope: parameter references and the
# free variables of string expressions (block parameters, event resets,
# state machine guards and actions), at any depth.
# Stop with a clear error when a subsystem contains, at any depth, a model
# that also contains it: such a hierarchy has no finite flattening. `path`
# holds the models above `model`, `trail` their names and block names.
check_acyclic <- function(model, path = list(), trail = model$name) {
  for (i in seq_along(path)) {
    if (identical(path[[i]], model)) {
      sim_abort("Subsystem cycle: ", paste(trail, collapse = " > "),
                " leads back to model '", model$name, "', which would contain itself.")
    }
  }
  path <- c(path, list(model))
  for (b in model$blocks) {
    sys <- if (identical(b$type, "subsystem")) b$params$system
    if (inherits(sys, "sim_model")) {
      check_acyclic(sys, path, c(trail, paste0(b$name, " (", sys$name, ")")))
    }
  }
  invisible(model)
}

referenced_params <- function(model) {
  out <- character()
  walk <- function(v) {
    if (inherits(v, "sim_param")) {
      out <<- c(out, if (is.null(v$expr)) v$name else v$vars)
    } else if (inherits(v, "sim_model")) {
      out <<- c(out, names(v$params), referenced_params(v))
    } else if (inherits(v, "sim_transition")) {
      walk(v$guard); walk(v$action)
    } else if (is.character(v)) {
      for (s in v) {
        e <- tryCatch(str2lang(s), error = function(e) NULL)
        if (!is.null(e)) out <<- c(out, all.vars(e))
      }
    } else if (is.list(v)) {
      for (el in v) walk(el)
    }
  }
  for (b in model$blocks) walk(b$params)
  for (e in model$events) walk(e$reset)
  unique(out)
}

# Flatten a model into a list of leaf blocks and connections. Each flat block
# keeps its raw parameters and the chain of parameter scopes it lives in.
flatten_model <- function(model, prefix = "", scopes = list(), depth = 0L) {
  if (depth > 50L) sim_abort("Subsystems are nested more than 50 levels deep ",
                             "(a subsystem may contain itself).")
  scopes <- c(scopes, list(model$params))
  blocks <- list()
  aliases <- list()
  out_map <- list()
  conns <- model$connections
  conns$from <- paste0(prefix, conns$from, recycle0 = TRUE)
  conns$to <- paste0(prefix, conns$to, recycle0 = TRUE)
  events <- lapply(model$events, function(e) {
    e$name <- paste0(prefix, e$name)
    if (e$trigger$kind == "zero_crossing") e$trigger$signal <- paste0(prefix, e$trigger$signal)
    if (length(e$reset)) names(e$reset) <- paste0(prefix, names(e$reset))
    e$scopes <- scopes
    e
  })
  names(events) <- vapply(events, `[[`, "", "name")
  for (b in model$blocks) {
    full <- paste0(prefix, b$name)
    if (b$type != "subsystem") {
      blocks[[full]] <- list(name = full, type = b$type, params = b$params,
                             scopes = rev(scopes))
      next
    }
    sys <- b$params$system
    if (!inherits(sys, "sim_model")) {
      sim_abort("Block '", full, "': `system` must be a sim_model.")
    }
    inst <- b$params$params %||% list()
    check_param_list(inst, paste0("Block '", full, "': `params`"))
    unused <- setdiff(names(inst), c(names(sys$params), referenced_params(sys)))
    if (length(unused)) {
      sim_abort("Block '", full, "': instance parameter(s) ", paste(unused, collapse = ", "),
                " are neither parameters of subsystem '", sys$name,
                "' nor used inside it. Subsystem parameters: ",
                if (length(sys$params)) paste(names(sys$params), collapse = ", ") else "none", ".")
    }
    # instance parameters override the subsystem defaults; both form the
    # innermost scope of the subsystem's blocks
    sub_params <- sys$params
    sub_params[names(inst)] <- inst
    sys_view <- list(name = sys$name, params = sub_params, blocks = sys$blocks,
                     connections = sys$connections, events = sys$events)
    inner <- flatten_model(sys_view, paste0(full, "/"), scopes, depth + 1L)
    blocks <- c(blocks, inner$blocks)
    events <- c(events, inner$events)
    aliases <- c(aliases, inner$aliases)
    inports <- paste0(full, "/", subsystem_ports(sys, "inport"))
    outports <- paste0(full, "/", subsystem_ports(sys, "outport"))
    # the instance name refers to its outputs in results
    alias_names <- if (length(outports) == 1L) full else
      paste0(full, ".", subsystem_ports(sys, "outport"))
    aliases[alias_names] <- as.list(outports)
    out_map[[full]] <- outports
    into <- conns$to == full
    if (any(conns$to_port[into] > length(inports))) {
      sim_abort("Block '", full, "' has ", length(inports), " input(s).")
    }
    conns$to[into] <- inports[conns$to_port[into]]
    conns$to_port[into] <- 1L
    outof <- conns$from == full
    if (any(conns$from_port[outof] > length(outports))) {
      sim_abort("Block '", full, "' has ", length(outports), " output(s).")
    }
    conns$from[outof] <- outports[conns$from_port[outof]]
    conns$from_port[outof] <- 1L
    conns <- rbind(conns, inner$connections)
  }
  # events of this level that watch a subsystem instance watch its outport
  for (k in seq_along(events)) {
    tr <- events[[k]]$trigger
    if (tr$kind != "zero_crossing" || is.null(out_map[[tr$signal]])) next
    ports <- out_map[[tr$signal]]
    if (tr$port > length(ports)) {
      sim_abort("Event '", events[[k]]$name, "': subsystem '", tr$signal, "' has ",
                length(ports), " output(s).")
    }
    tr$signal <- ports[tr$port]
    tr$port <- 1L
    events[[k]]$trigger <- tr
  }
  list(blocks = blocks, connections = conns, events = events, aliases = aliases)
}
