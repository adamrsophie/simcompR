# Block type registry ------------------------------------------------------

.registry <- new.env(parent = emptyenv())

#' Register a block type
#'
#' Every block in a model is an instance of a registered block type. A type
#' declares its parameters and a `build` function that turns resolved
#' parameter values into a block definition: ports, sizes, direct
#' feedthrough, states, sample time and the functions that compute outputs,
#' derivatives and discrete updates. The built-in library is registered when
#' the package loads; `sim_register_block_type()` adds user types.
#'
#' The list returned by `build(p)` may contain:
#'
#' * `inputs`, `outputs`: integer vectors of port sizes. `NA` means the size
#'   is inherited from the connected signal (inputs) or computed by
#'   `size_out` (outputs). Names, when given, become port names.
#' * `size_out(in_sizes)`: returns output sizes once input sizes are known.
#' * `check_sizes(in_sizes, out_sizes)`: optional; signals an error when the
#'   resolved sizes are inconsistent.
#' * `feedthrough`: logical, one per input port. `TRUE` when the outputs at
#'   time `t` depend on that input at time `t`.
#' * `nx`, `x0`: number and initial value of continuous states.
#' * `nd`, `d0`: number and initial value of discrete states.
#' * `sample`: `NULL` for continuous (or inherited) blocks, otherwise
#'   `c(period, offset)`. Discrete blocks compute outputs only at their
#'   sample hits and hold them in between.
#' * `output(t, x, d, u)`: returns a list with one numeric vector per output
#'   port (a plain vector is accepted for one port). `u` is a list of input
#'   vectors.
#' * `deriv(t, x, d, u)`: continuous state derivative.
#' * `update(t, x, d, u)`: next discrete state, called at sample hits after
#'   all outputs at that instant have been computed.
#' * `held0(size, d0)`: initial held output of a discrete block before its
#'   first sample hit. `size` holds the output port sizes and `d0`, when the
#'   function accepts it, the initial discrete state (including overrides
#'   given to [sim_run()]).
#' * `zc(t, x, d, u)` and `on_zc(i, direction, t, x, d, u)`: zero-crossing
#'   functions and their handler, which returns `list(d = , x = , terminal = )`.
#'   `zc_direction` gives one of `"rising"`, `"falling"`, `"either"` per
#'   component.
#' * `breakpoints(t0, tf)`: times where the block's output is
#'   discontinuous; the integrator restarts there.
#' * `modes(t, x, d, u)`: for switching blocks whose output branch is
#'   stored in the discrete state `d` (like `switch` and `relay`), returns
#'   the mode for the current inputs. It is called at the start and at
#'   every segment boundary, repeatedly until all modes settle, never
#'   inside a segment; the crossing functions in `zc` move the mode inside a
#'   segment.
#' * `stop_if(t, x, d, u)`: returns `TRUE` to end the simulation; checked
#'   at every boundary and after every event (see the `stop` block).
#' * `kinks(t, x, d, u)` and `kink_size(in_sizes)`: crossing functions whose
#'   zeros are points where a memoryless output or its slope jumps (the
#'   corners of a saturation), and their number for the given input sizes.
#'   The integrator restarts exactly at each kink; nothing is logged. They
#'   are ignored by methods without root finding.
#'
#' All functions must be free of side effects: the solver calls them at
#' trial points that it may later reject.
#'
#' @param type Type name (a single string).
#' @param build Function of one argument, the named list of resolved
#'   parameters, returning a block definition.
#' @param params Named list of parameter defaults.
#' @param required Names of parameters the user must supply. By default,
#'   every parameter whose default is `NULL`.
#' @param category Category used by [sim_block_types()].
#' @param description One line description.
#' @param portable `TRUE` when every parameter can be written to JSON, so
#'   the block can be saved with [sim_write_json()].
#' @param overwrite Allow replacing an existing type.
#' @return The type name, invisibly.
#' @examples
#' sim_register_block_type(
#'   "cube", overwrite = TRUE,
#'   build = function(p) list(
#'     inputs = NA, outputs = NA, feedthrough = TRUE,
#'     size_out = function(n) n,
#'     output = function(t, x, d, u) u[[1]]^3
#'   ),
#'   category = "math", description = "Elementwise cube"
#' )
#' m <- sim_model("cube_demo")
#' sim_add_block(m, "time", "t")
#' sim_add_block(m, "cube", "c")
#' sim_connect(m, "t", "c")
#' r <- sim_run(sim_compile(m), tspan = c(0, 2), times = c(0, 1, 2))
#' sim_signal(r, "c")
#' @export
sim_register_block_type <- function(type, build, params = list(),
                                    required = names(params)[vapply(params, is.null, TRUE)],
                                    category = "user", description = "",
                                    portable = TRUE, overwrite = FALSE) {
  if (!is.character(type) || length(type) != 1L || !nzchar(type)) {
    sim_abort("`type` must be a non-empty string.")
  }
  if (!is.function(build)) sim_abort("`build` must be a function.")
  if (!is.list(params) || (length(params) && is.null(names(params)))) {
    sim_abort("`params` must be a named list.")
  }
  if (exists(type, envir = .registry, inherits = FALSE) && !overwrite) {
    sim_abort("Block type '", type, "' is already registered. ",
              "Use `overwrite = TRUE` to replace it.")
  }
  force(required)
  if (length(setdiff(required, names(params)))) {
    sim_abort("`required` must name parameters listed in `params`.")
  }
  assign(type, list(type = type, build = build, params = params, required = required,
                    category = category, description = description,
                    portable = portable), envir = .registry)
  invisible(type)
}

#' List registered block types
#'
#' @param category Optional category filter.
#' @return A data frame with columns `type`, `category`, `description` and
#'   `parameters`.
#' @examples
#' head(sim_block_types())
#' sim_block_types("discrete")
#' @export
sim_block_types <- function(category = NULL) {
  types <- sort(ls(.registry))
  rows <- lapply(types, function(t) {
    def <- get(t, envir = .registry)
    data.frame(type = t, category = def$category,
               description = def$description,
               parameters = paste(names(def$params), collapse = ", "),
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  if (!is.null(category)) out <- out[out$category %in% category, , drop = FALSE]
  rownames(out) <- NULL
  out
}

get_block_type <- function(type) {
  if (!is.character(type) || length(type) != 1L) {
    sim_abort("Block `type` must be a single string.")
  }
  if (!exists(type, envir = .registry, inherits = FALSE)) {
    sim_abort("Unknown block type '", type, "'. See sim_block_types().")
  }
  get(type, envir = .registry)
}

# Parameter references -------------------------------------------------------

#' Refer to a model or subsystem parameter
#'
#' A block parameter given as `sim_param("name")` is looked up when the
#' model is compiled: first in the parameters of the enclosing subsystem
#' instance, then in the parameters of the parent, up to the top-level
#' model. Changing a referenced parameter with [sim_set_params()] or
#' [sim_sweep()] does not require editing the block.
#'
#' References can be combined with numbers and with each other by
#' arithmetic (`-sim_param("a")`, `1 / sim_param("tau")`,
#' `2 * sim_param("a") + sim_param("b")`) and by the functions `abs`,
#' `sqrt`, `exp`, `log`, `sin`, `cos`, `tan`, `floor`, `ceiling` and
#' `sign`. The result is a parameter expression, kept as a portable string
#' (so it can be saved with [sim_write_json()]) and evaluated at compile
#' time. `sim_param_expr()` builds one directly from a string. References
#' may also appear inside list parameters. A vector cannot hold
#' references (R would drop their class); use an expression such as
#' `sim_param_expr("c(a, 2 * a)")` instead.
#'
#' @param name Parameter name.
#' @param expr A portable expression (string) in parameter names.
#' @return An object of class `sim_param`.
#' @examples
#' m <- sim_model("decay", params = list(a = 2))
#' sim_add_block(m, "integrator", "x", initial = 1)
#' sim_add_block(m, "gain", "k", gain = sim_param("a"))
#' sim_connect(m, "x", "k")
#' sim_add_block(m, "gain", "neg", gain = -1)
#' sim_connect(m, "k", "neg")
#' sim_connect(m, "neg", "x")
#' plan <- sim_compile(m)
#' sim_param_values(plan)
#'
#' # The same model with one gain: k = -a
#' m2 <- sim_model("decay", params = list(a = 2))
#' sim_add_block(m2, "integrator", "x", initial = 1)
#' sim_add_block(m2, "gain", "k", gain = -sim_param("a"))
#' sim_connect(m2, "x", "k")
#' sim_connect(m2, "k", "x")
#' sim_param_expr("-a / 2")
#' @export
sim_param <- function(name) {
  if (!is.character(name) || length(name) != 1L || !nzchar(name)) {
    sim_abort("`name` must be a non-empty string.")
  }
  structure(list(name = name), class = "sim_param")
}

#' @rdname sim_param
#' @export
sim_param_expr <- function(expr) {
  if (!is.character(expr) || length(expr) != 1L || !nzchar(expr)) {
    sim_abort("`expr` must be a non-empty string.")
  }
  e <- parse_safe(expr, "Parameter expression")
  vars <- setdiff(all.vars(e), c("pi", "TRUE", "FALSE"))
  structure(list(name = NULL, expr = expr, vars = vars), class = "sim_param")
}

param_text <- function(v) if (is.null(v$expr)) v$name else paste0("(", v$expr, ")")

#' @export
Ops.sim_param <- function(e1, e2) {
  if (!.Generic %in% c("+", "-", "*", "/", "^")) {
    sim_abort("Operator `", .Generic, "` is not available for sim_param(); use +, -, *, / or ^.")
  }
  txt <- function(v) {
    if (is_param_ref(v)) return(param_text(v))
    if (is.numeric(v) && length(v) == 1L && is.finite(v)) return(sprintf("%.17g", v))
    sim_abort("Arithmetic on sim_param() takes parameters and single finite numbers.")
  }
  sim_param_expr(if (missing(e2)) paste0(.Generic, txt(e1)) else
    paste(txt(e1), .Generic, txt(e2)))
}

#' @export
Math.sim_param <- function(x, ...) {
  ok <- c("abs", "sqrt", "exp", "log", "sin", "cos", "tan", "floor", "ceiling", "sign")
  if (!.Generic %in% ok || length(list(...))) {
    sim_abort("Function `", .Generic, "` is not available for sim_param(); use one of ",
              paste(ok, collapse = ", "), ".")
  }
  sim_param_expr(paste0(.Generic, "(", param_text(x), ")"))
}

#' @export
print.sim_param <- function(x, ...) {
  if (is.null(x$expr)) cat("<sim_param ", x$name, ">\n", sep = "")
  else cat("<sim_param expression ", x$expr, ">\n", sep = "")
  invisible(x)
}

is_param_ref <- function(x) inherits(x, "sim_param")

# Resolve a block parameter list against a chain of scopes.
# scopes[[1]] is the innermost scope, the last element is the model scope.
resolve_params <- function(params, scopes, block_name) {
  lapply(params, resolve_value, scopes = scopes, level = 1L,
         block_name = block_name)
}

resolve_value <- function(value, scopes, level, block_name) {
  if (is.list(value) && !is.object(value)) {
    return(lapply(value, resolve_value, scopes = scopes, level = level,
                  block_name = block_name))
  }
  if (!is_param_ref(value)) return(value)
  if (!is.null(value$expr)) {
    vals <- lapply(value$vars, function(n) resolve_value(sim_param(n), scopes, level, block_name))
    names(vals) <- value$vars
    return(tryCatch(eval_safe(parse_safe(value$expr, "Parameter expression"), vals),
                    error = function(e) {
                      sim_abort("Block '", block_name, "': parameter expression '",
                                value$expr, "' failed: ", conditionMessage(e))
                    }))
  }
  name <- value$name
  for (i in seq.int(level, length.out = max(0L, length(scopes) - level + 1L))) {
    scope <- scopes[[i]]
    if (name %in% names(scope)) {
      return(resolve_value(scope[[name]], scopes, i + 1L, block_name))
    }
  }
  sim_abort("Block '", block_name, "': parameter reference '", name,
            "' is not defined in any enclosing scope.")
}

# Build a block definition from a type and resolved parameters, normalizing
# every optional field so the compiler can rely on them.
build_block <- function(type_def, params, block_name) {
  scope <- params$.scope
  params$.scope <- NULL
  defaults <- type_def$params
  unknown <- setdiff(names(params), names(defaults))
  if (length(unknown)) {
    sim_abort("Block '", block_name, "' (", type_def$type, "): unknown ",
              "parameter(s) ", paste(unknown, collapse = ", "), ". Known: ",
              paste(names(defaults), collapse = ", "), ".")
  }
  p <- defaults
  p[names(params)] <- params
  if (isTRUE(type_def$scoped)) p$.scope <- scope
  missing <- setdiff(type_def$required, names(params))
  if (length(missing)) {
    sim_abort("Block '", block_name, "' (", type_def$type, "): required ",
              "parameter(s) missing: ", paste(missing, collapse = ", "), ".")
  }
  def <- withCallingHandlers(
    type_def$build(p),
    error = function(e) {
      msg <- conditionMessage(e)
      if (!startsWith(msg, "Block '")) {
        sim_abort("Block '", block_name, "' (", type_def$type, "): ", msg)
      }
    }
  )
  normalize_definition(def, block_name, type_def$type)
}

normalize_definition <- function(def, block_name, type) {
  where <- paste0("Block '", block_name, "' (", type, ")")
  in_nm <- names(def$inputs)
  out_nm <- names(def$outputs)
  def$inputs <- as.integer(def$inputs %||% integer())
  def$outputs <- as.integer(def$outputs %||% integer())
  ni <- length(def$inputs)
  no <- length(def$outputs)
  def$input_names <- def$input_names %||% in_nm %||%
    if (ni) paste0("in", seq_len(ni)) else character()
  def$output_names <- def$output_names %||% out_nm %||%
    if (no) paste0("out", seq_len(no)) else character()
  def$feedthrough <- as.logical(def$feedthrough %||% rep(FALSE, ni))
  if (length(def$feedthrough) == 1L && ni > 1L) {
    def$feedthrough <- rep(def$feedthrough, ni)
  }
  if (length(def$feedthrough) != ni) {
    sim_abort(where, ": `feedthrough` must have one value per input port.")
  }
  def$x0 <- as.numeric(def$x0 %||% numeric())
  def$nx <- as.integer(def$nx %||% length(def$x0))
  if (length(def$x0) != def$nx) {
    sim_abort(where, ": `x0` must have length `nx`.")
  }
  def$d0 <- as.numeric(def$d0 %||% numeric())
  def$nd <- as.integer(def$nd %||% length(def$d0))
  if (length(def$d0) != def$nd) {
    sim_abort(where, ": `d0` must have length `nd`.")
  }
  if (def$nx > 0L && !is.function(def$deriv)) {
    sim_abort(where, ": blocks with continuous states need `deriv`.")
  }
  if (no > 0L && !is.function(def$output)) {
    sim_abort(where, ": blocks with outputs need `output`.")
  }
  if (!is.null(def$sample)) {
    def$sample <- check_sample_time(def$sample[[1]], def$sample[[2]], where)
    if (def$nx > 0L) {
      sim_abort(where, ": a discrete block cannot have continuous states.")
    }
  }
  if (!is.null(def$update) && is.null(def$sample)) {
    sim_abort(where, ": `update` requires a sample time.")
  }
  def$nzc <- 0L
  if (is.function(def$zc)) {
    def$zc_direction <- def$zc_direction %||% "either"
    def$nzc <- as.integer(length(def$zc_direction))
    bad <- setdiff(def$zc_direction, c("rising", "falling", "either"))
    if (length(bad)) {
      sim_abort(where, ": `zc_direction` must be 'rising', 'falling' or 'either'.")
    }
    if (!is.function(def$on_zc)) sim_abort(where, ": `zc` requires `on_zc`.")
  }
  def
}
