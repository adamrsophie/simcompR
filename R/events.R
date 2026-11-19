#' Model events: zero crossings, scheduled times and state resets
#'
#' `sim_add_event()` attaches an event to a model. When the trigger fires
#' the event can reset the state of any block and can stop the simulation.
#'
#' A zero-crossing trigger watches `signal - level`, where `signal` is an
#' output of a block. The crossing time is located by the root finder of the
#' integrator (methods `"lsoda"`, `"lsode"`, `"lsodes"` and `"lsodar"`), the
#' solver is stopped there, the resets are applied and integration restarts
#' from the new state. `direction` selects rising crossings (negative to
#' positive), falling crossings or both. A crossing in the other direction
#' is ignored. A scheduled trigger fires at fixed times; the solver always
#' stops exactly at those times.
#'
#' A reset is a named list. Each name is a block that has state; the value
#' is either a function `function(x, s, t)` or a string holding a portable
#' expression. `x` is the current state of that block (continuous states
#' for continuous blocks, discrete states otherwise), `s` is a named list of
#' all block outputs at the event time (before the reset), and `t` is the
#' time. Model parameters are also visible in string expressions. Only
#' string resets can be saved with [sim_write_json()].
#'
#' A `signal` or reset target inside a subsystem instance is written as a
#' path, `"instance/block"` (deeper levels add more `/` parts). A crossing
#' function that sits exactly on `level` fires when the signal moves away
#' in the trigger direction, including at the initial time.
#'
#' When several events fire at the same instant (scheduled events, located
#' crossings and crossings caused by a jump at that instant) they are
#' applied in increasing `priority`, then by name, and every reset sees the
#' state left by the previous one. A reset that makes another event's
#' signal cross its level fires that event at the same instant. A terminal
#' event stops the run once the instant is settled. See [sim_run()] for the
#' full order of what happens at an instant.
#'
#' @param model A [sim_model()].
#' @param name Unique event name; for `sim_remove_event()`, the event to
#'   remove.
#' @param trigger A trigger from `sim_zero_crossing()` or `sim_at()`.
#' @param reset Named list of state resets (see Details).
#' @param terminal If `TRUE`, the simulation stops after the instant at
#'   which the event fires.
#' @param priority Ordering among simultaneous events (smaller first).
#' @param signal Block (name or handle) whose output is watched.
#' @param direction `"rising"`, `"falling"` or `"either"`.
#' @param level Crossing level.
#' @param port Output port of `signal`.
#' @param element Element of a vector signal.
#' @param times Event times.
#' @return `sim_add_event()` and `sim_remove_event()` return the model
#'   invisibly; the trigger
#'   constructors return trigger objects.
#' @examples
#' # Bouncing ball: h'' = -g, reverse and damp the velocity at impact.
#' m <- sim_model("ball", params = list(e = 0.8))
#' sim_add_block(m, "constant", "gravity", value = -9.81)
#' sim_add_block(m, "integrator", "v", initial = 0)
#' sim_add_block(m, "integrator", "h", initial = 1)
#' sim_connect(m, "gravity", "v")
#' sim_connect(m, "v", "h")
#' sim_add_event(m, "impact",
#'   trigger = sim_zero_crossing("h", direction = "falling"),
#'   reset = list(v = "-e * x"))
#' r <- sim_run(sim_compile(m), tspan = c(0, 1.5))
#' sim_events(r)
#' @export
sim_add_event <- function(model, name, trigger, reset = list(),
                          terminal = FALSE, priority = 0) {
  check_model(model)
  if (!is.character(name) || length(name) != 1L || is.na(name) || !nzchar(name)) {
    sim_abort("Event `name` must be a non-empty string.")
  }
  if (grepl("/", name, fixed = TRUE)) {
    sim_abort("Event name '", name, "' must not contain '/': it separates the names ",
              "of subsystem instances, which qualify the events inside them.")
  }
  if (name %in% names(model$events)) {
    sim_abort("Model '", model$name, "' already has an event named '", name, "'.")
  }
  if (!inherits(trigger, "sim_trigger")) {
    sim_abort("`trigger` must come from sim_zero_crossing() or sim_at().")
  }
  if (trigger$kind == "zero_crossing") {
    trigger$signal <- block_path_of(model, trigger$owner %||% trigger$signal, "signal")
    trigger$owner <- NULL
  }
  if (!is.list(reset) || (length(reset) && is.null(names(reset)))) {
    sim_abort("`reset` must be a named list of block resets.")
  }
  rn <- names(reset)
  if (anyNA(rn) || any(!nzchar(rn)) || anyDuplicated(rn)) {
    sim_abort("`reset` must name each block once; ",
              if (anyDuplicated(rn)) paste0("'", rn[anyDuplicated(rn)], "' is repeated.")
              else "a name is missing.")
  }
  for (b in names(reset)) {
    block_path_of(model, b, "reset target")
    r <- reset[[b]]
    if (is.character(r)) parse_safe(r, paste0("Reset of block '", b, "'"))
    else if (!is.function(r)) {
      sim_abort("Reset of block '", b, "' must be a function or a string.")
    }
  }
  check_number(priority, "`priority`")
  if (!is.logical(terminal) || length(terminal) != 1L || is.na(terminal)) {
    sim_abort("`terminal` must be TRUE or FALSE.")
  }
  if (!length(reset) && !terminal) {
    sim_warn("Event '", name, "' has no reset and is not terminal; it will only be logged.")
  }
  model$events[[name]] <- list(name = name, trigger = trigger, reset = reset,
                               terminal = terminal, priority = priority)
  touch(model)
}

#' @rdname sim_add_event
#' @export
sim_remove_event <- function(model, name) {
  check_model(model)
  if (!is.character(name) || length(name) != 1L || !name %in% names(model$events)) {
    sim_abort("Model '", model$name, "' has no event named '", paste(name, collapse = ", "), "'.")
  }
  model$events[[name]] <- NULL
  touch(model)
}

#' @rdname sim_add_event
#' @export
sim_zero_crossing <- function(signal, direction = "either", level = 0,
                              port = 1, element = 1) {
  owner <- NULL
  if (inherits(signal, "sim_block")) {
    owner <- signal
    signal <- signal$name
  }
  if (!is.character(signal) || length(signal) != 1L) {
    sim_abort("`signal` must be a block name or handle.")
  }
  direction <- match.arg(direction, c("either", "rising", "falling"))
  check_number(level, "`level`")
  if (!is_count(port) || !is_count(element)) {
    sim_abort("`port` and `element` must be positive integers.")
  }
  structure(list(kind = "zero_crossing", signal = signal, direction = direction,
                 level = as.numeric(level), port = as.integer(port),
                 element = as.integer(element), owner = owner),
            class = "sim_trigger")
}

#' @rdname sim_add_event
#' @export
sim_at <- function(times) {
  times <- check_numeric(times, "`times`")
  structure(list(kind = "scheduled", times = sort(unique(times))),
            class = "sim_trigger")
}

#' @export
print.sim_trigger <- function(x, ...) {
  if (x$kind == "scheduled") {
    cat("<sim_trigger at ", paste(format(x$times), collapse = ", "), ">\n", sep = "")
  } else {
    cat("<sim_trigger ", x$direction, " crossing of ", x$signal, ":", x$port,
        "[", x$element, "] through ", x$level, ">\n", sep = "")
  }
  invisible(x)
}
