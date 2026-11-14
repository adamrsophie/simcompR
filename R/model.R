#' Create a block-diagram model
#'
#' A `sim_model` is a mutable container of blocks, connections and events.
#' Functions such as [sim_add_block()] and [sim_connect()] modify it in
#' place. A model is never simulated directly: [sim_compile()] takes a
#' snapshot and turns it into an immutable [execution plan][sim_compile()],
#' so later edits to the model never change an existing plan.
#'
#' @param name Model name.
#' @param params Named list of model parameters, referenced from blocks with
#'   [sim_param()].
#' @return An object of class `sim_model`.
#' @seealso [sim_add_block()], [sim_connect()], [sim_compile()].
#' @examples
#' m <- sim_model("feedback_system")
#' ref <- sim_add_block(m, "constant", "reference", value = 1)
#' err <- sim_add_block(m, "sum", "error", signs = "+-")
#' plant <- sim_add_block(m, "integrator", "plant", initial = 0)
#' sim_connect(m, ref, err, input = 1)
#' sim_connect(m, plant, err, input = 2)
#' sim_connect(m, err, plant)
#' m
#' @export
sim_model <- function(name = "model", params = list()) {
  if (!is.character(name) || length(name) != 1L || !nzchar(name)) {
    sim_abort("`name` must be a non-empty string.")
  }
  check_param_list(params, "`params`")
  m <- new.env(parent = emptyenv())
  m$name <- name
  m$params <- params
  m$blocks <- list()
  m$connections <- data.frame(from = character(), from_port = integer(),
                              to = character(), to_port = integer(),
                              stringsAsFactors = FALSE)
  m$events <- list()
  m$version <- 0L
  class(m) <- "sim_model"
  m
}

check_param_list <- function(params, what) {
  if (!is.list(params)) sim_abort(what, " must be a named list.")
  if (length(params) &&
      (is.null(names(params)) || any(!nzchar(names(params))) ||
       anyDuplicated(names(params)))) {
    sim_abort(what, " must be a list with unique, non-empty names.")
  }
  invisible(params)
}

check_model <- function(model) {
  if (!inherits(model, "sim_model")) sim_abort("`model` must be a sim_model.")
  invisible(model)
}

touch <- function(model) {
  model$version <- model$version + 1L
  invisible(model)
}

#' Add a block to a model
#'
#' @param model A [sim_model()].
#' @param type Block type, see [sim_block_types()].
#' @param name Unique block name. Names may not contain `/`, which separates
#'   subsystem levels in compiled names.
#' @param ... Block parameters. Any parameter may be a [sim_param()]
#'   reference.
#' @return A `sim_block` handle that can be passed to [sim_connect()].
#'   The model is modified in place.
#' @examples
#' m <- sim_model()
#' b <- sim_add_block(m, "sine", "u", amplitude = 2, frequency = 3)
#' b
#' @export
sim_add_block <- function(model, type, name, ...) {
  check_model(model)
  if (!is.character(name) || length(name) != 1L || !nzchar(name) ||
      grepl("/", name, fixed = TRUE)) {
    sim_abort("Block `name` must be a non-empty string without '/'.")
  }
  if (name %in% names(model$blocks)) {
    sim_abort("Model '", model$name, "' already has a block named '", name, "'.")
  }
  type_def <- get_block_type(type)
  params <- list(...)
  if (length(params) && (is.null(names(params)) || any(!nzchar(names(params))))) {
    sim_abort("Block parameters must be named.")
  }
  unknown <- setdiff(names(params), names(type_def$params))
  if (length(unknown)) {
    sim_abort("Block '", name, "' (", type, "): unknown parameter(s) ",
              paste(unknown, collapse = ", "), ". Known: ",
              paste(names(type_def$params), collapse = ", "), ".")
  }
  ports <- probe_ports(type_def, params, model$params, name)
  model$blocks[[name]] <- list(name = name, type = type, params = params,
                               ports = ports)
  touch(model)
  invisible(block_handle(model, name))
}

# Try to build the block to learn its port layout at connection time. The
# build can fail when parameters refer to subsystem parameters that are not
# known yet; port checks are then left to sim_compile().
probe_ports <- function(type_def, params, model_params, name) {
  resolved <- tryCatch(resolve_params(params, list(model_params), name),
                       error = function(e) NULL)
  if (is.null(resolved)) return(NULL)
  def <- build_block(type_def, resolved, name)
  list(inputs = def$input_names, outputs = def$output_names)
}

block_handle <- function(model, name) {
  b <- model$blocks[[name]]
  structure(list(name = name, type = b$type, model = model$name),
            class = "sim_block")
}

#' @export
print.sim_block <- function(x, ...) {
  cat("<sim_block '", x$name, "' of type ", x$type, " in model '", x$model,
      "'>\n", sep = "")
  invisible(x)
}

block_name_of <- function(model, block, what) {
  if (inherits(block, "sim_block")) block <- block$name
  if (!is.character(block) || length(block) != 1L) {
    sim_abort("`", what, "` must be a block name or a sim_block handle.")
  }
  if (!block %in% names(model$blocks)) {
    sim_abort("Model '", model$name, "' has no block named '", block, "'.")
  }
  block
}

port_index <- function(model, block, port, direction) {
  ports <- model$blocks[[block]]$ports
  names_ <- if (is.null(ports)) NULL else ports[[paste0(direction, "s")]]
  if (is.character(port) && length(port) == 1L) {
    if (is.null(names_)) {
      sim_abort("Block '", block, "': ports can be named only when the ",
                "block's parameters are known; use an index.")
    }
    idx <- match(port, names_)
    if (is.na(idx)) {
      sim_abort("Block '", block, "' has no ", direction, " port '", port,
                "'. Ports: ", paste(names_, collapse = ", "), ".")
    }
    return(idx)
  }
  if (!is_count(port)) {
    sim_abort("Block '", block, "': ", direction, " port must be a positive ",
              "integer index or a port name.")
  }
  port <- as.integer(port)
  if (!is.null(names_) && port > length(names_)) {
    sim_abort("Block '", block, "' has ", length(names_), " ", direction,
              " port(s); port ", port, " does not exist.")
  }
  port
}

#' Connect two blocks
#'
#' Connects output port `output` of block `from` to input port `input` of
#' block `to`. An output may feed several inputs; an input accepts exactly
#' one connection.
#'
#' @param model A [sim_model()].
#' @param from,to Blocks, given as names or `sim_block` handles.
#' @param output,input Port indices or port names.
#' @return The model, invisibly.
#' @examples
#' m <- sim_model()
#' sim_add_block(m, "constant", "a", value = 2)
#' sim_add_block(m, "gain", "k", gain = 3)
#' sim_connect(m, "a", "k")
#' m
#' @export
sim_connect <- function(model, from, to, output = 1, input = 1) {
  check_model(model)
  from <- block_name_of(model, from, "from")
  to <- block_name_of(model, to, "to")
  output <- port_index(model, from, output, "output")
  input <- port_index(model, to, input, "input")
  cn <- model$connections
  dup <- cn$to == to & cn$to_port == input
  if (any(dup)) {
    sim_abort("Input ", input, " of block '", to, "' is already connected ",
              "to '", cn$from[dup], "' output ", cn$from_port[dup],
              ". Disconnect it first with sim_disconnect().")
  }
  model$connections <- rbind(cn, data.frame(from = from, from_port = output,
                                            to = to, to_port = input,
                                            stringsAsFactors = FALSE))
  touch(model)
}

#' Remove a connection or a block
#'
#' @param model A [sim_model()].
#' @param to,input Destination block and input port of the connection to
#'   remove.
#' @param block Block to remove, together with all its connections.
#' @return The model, invisibly.
#' @examples
#' m <- sim_model()
#' sim_add_block(m, "constant", "a", value = 2)
#' sim_add_block(m, "gain", "k", gain = 3)
#' sim_connect(m, "a", "k")
#' sim_disconnect(m, "k", 1)
#' sim_remove_block(m, "a")
#' m
#' @export
sim_disconnect <- function(model, to, input = 1) {
  check_model(model)
  to <- block_name_of(model, to, "to")
  input <- port_index(model, to, input, "input")
  cn <- model$connections
  hit <- cn$to == to & cn$to_port == input
  if (!any(hit)) sim_abort("Input ", input, " of block '", to, "' is not connected.")
  model$connections <- cn[!hit, , drop = FALSE]
  touch(model)
}

#' @rdname sim_disconnect
#' @export
sim_remove_block <- function(model, block) {
  check_model(model)
  block <- block_name_of(model, block, "block")
  cn <- model$connections
  model$connections <- cn[cn$from != block & cn$to != block, , drop = FALSE]
  model$blocks[[block]] <- NULL
  touch(model)
}

#' Change block parameters in a model
#'
#' @param model A [sim_model()].
#' @param block Block name or handle.
#' @param ... New parameter values.
#' @return The model, invisibly.
#' @examples
#' m <- sim_model()
#' sim_add_block(m, "gain", "k", gain = 3)
#' sim_set_block_params(m, "k", gain = 5)
#' @export
sim_set_block_params <- function(model, block, ...) {
  check_model(model)
  block <- block_name_of(model, block, "block")
  b <- model$blocks[[block]]
  new <- list(...)
  type_def <- get_block_type(b$type)
  unknown <- setdiff(names(new), names(type_def$params))
  if (length(unknown)) {
    sim_abort("Block '", block, "': unknown parameter(s) ",
              paste(unknown, collapse = ", "), ".")
  }
  b$params[names(new)] <- new
  b$ports <- probe_ports(type_def, b$params, model$params, block)
  model$blocks[[block]] <- b
  touch(model)
}

#' Inspect the blocks and connections of a model
#'
#' @param model A [sim_model()].
#' @return `sim_blocks()` returns a data frame with one row per block;
#'   `sim_connections()` returns the connection table.
#' @examples
#' m <- sim_model()
#' sim_add_block(m, "constant", "a", value = 2)
#' sim_add_block(m, "gain", "k", gain = 3)
#' sim_connect(m, "a", "k")
#' sim_blocks(m)
#' sim_connections(m)
#' @export
sim_blocks <- function(model) {
  check_model(model)
  if (!length(model$blocks)) {
    return(data.frame(name = character(), type = character(),
                      category = character(), stringsAsFactors = FALSE))
  }
  data.frame(
    name = names(model$blocks),
    type = vapply(model$blocks, `[[`, "", "type"),
    category = vapply(model$blocks, function(b) get_block_type(b$type)$category, ""),
    stringsAsFactors = FALSE, row.names = NULL
  )
}

#' @rdname sim_blocks
#' @export
sim_connections <- function(model) {
  check_model(model)
  cn <- model$connections
  rownames(cn) <- NULL
  cn
}

#' @export
print.sim_model <- function(x, ...) {
  cat("<sim_model '", x$name, "'> ", length(x$blocks), " block(s), ",
      nrow(x$connections), " connection(s), ", length(x$events),
      " event(s)\n", sep = "")
  if (length(x$params)) {
    cat("Parameters: ", paste(names(x$params), collapse = ", "), "\n", sep = "")
  }
  if (length(x$blocks)) {
    b <- sim_blocks(x)
    for (i in seq_len(nrow(b))) cat("  ", b$name[i], " [", b$type[i], "]\n", sep = "")
  }
  if (nrow(x$connections)) {
    cn <- x$connections
    cat("Connections:\n")
    cat(paste0("  ", cn$from, ":", cn$from_port, " -> ", cn$to, ":", cn$to_port),
        sep = "\n")
  }
  invisible(x)
}

#' Copy a model
#'
#' Models are modified in place, so assigning a model to a new variable
#' does not copy it. `sim_copy()` returns an independent copy.
#'
#' @param model A [sim_model()].
#' @param name Optional new name.
#' @return A new `sim_model`.
#' @examples
#' m <- sim_model("a")
#' m2 <- sim_copy(m, "b")
#' sim_add_block(m2, "constant", "c", value = 1)
#' length(m$blocks)
#' length(m2$blocks)
#' @export
sim_copy <- function(model, name = model$name) {
  check_model(model)
  out <- sim_model(name, model$params)
  for (f in c("blocks", "connections", "events")) out[[f]] <- model[[f]]
  out
}
