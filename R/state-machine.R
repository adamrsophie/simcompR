#' State machine transitions
#'
#' Helper to declare one transition of a `"state_machine"` block.
#'
#' A state machine block is a discrete block that executes at every sample
#' hit of its clock (`period`, `offset`). At a hit it:
#'
#' 1. collects the transitions leaving the current state, sorted by
#'    `priority` (smaller first) and then by their position in the list;
#' 2. takes the first transition whose guard holds and whose `after` time
#'    (time spent in the current state) has elapsed;
#' 3. runs the exit action of the old state, the transition action and the
#'    entry action of the new state, in that order, then records the entry
#'    time;
#' 4. when no transition fires, runs the `during` action of the current
#'    state.
#'
#' The entry action of the initial state runs at the first sample hit,
#' before transitions are considered, with the inputs and time of that hit.
#'
#' At most one transition fires per sample hit, so evaluation never loops.
#' Outputs are computed after the transition, which means the outputs at a
#' hit depend on the inputs at the same hit (direct feedthrough).
#'
#' Guards and actions are portable strings or R functions. A guard string
#' is an expression in `u` (the input vector), `t`, `elapsed` (time in the
#' current state), the data variables and the parameters of the enclosing
#' model or subsystem. An action string assigns data
#' variables, for example `"count <- count + 1"`; each data variable must
#' still hold one finite number afterwards. A guard function has
#' signature `function(u, t, elapsed, data)` and returns `TRUE` or `FALSE`;
#' an action function has signature `function(data, u, t)` and returns the
#' new data vector: either unnamed with one value per variable, in the
#' order of `data`, or named with exactly the data variable names, in any
#' order (values are matched by name). Guards and actions must be pure
#' functions of their arguments, because the block evaluates them to
#' compute its outputs and again to compute its next state.
#'
#' `after` compares the time spent in the state, `t - entry time`, with
#' the requested duration, allowing only for the rounding of that
#' subtraction, so the behavior does not depend on where the time axis
#' starts.
#'
#' The block's parameters are `states` (character), `initial` (state name),
#' `transitions` (a list of `sim_transition()`), `entry`, `exit`, `during`
#' (named lists of actions by state), `data` (named numeric initial values
#' of data variables), `outputs` (optional named list of numeric output
#' vectors by state), `n_inputs`, `period` and `offset`. Output ports are
#' `state` (index of the active state in `states`), then `y` when `outputs`
#' is given, then `data` when there are data variables.
#'
#' @param from,to State names.
#' @param guard Guard expression (string) or function; `"TRUE"` by default.
#' @param after Minimum time in `from` before the transition may fire.
#' @param action Optional transition action.
#' @param priority Optional priority; smaller values are tried first.
#' @return A `sim_transition` object.
#' @seealso [sim_check_state_machine()].
#' @examples
#' # A thermostat that heats below 19 degrees and stops above 21.
#' m <- sim_model("thermostat")
#' sim_add_block(m, "integrator", "temp", initial = 18)
#' sim_add_block(m, "state_machine", "ctrl", period = 0.1,
#'   states = c("heating", "idle"),
#'   transitions = list(
#'     sim_transition("heating", "idle", "u[1] > 21"),
#'     sim_transition("idle", "heating", "u[1] < 19")),
#'   outputs = list(heating = 1, idle = 0))
#' sim_add_block(m, "expression", "rate", expr = "2 * u[1] - 0.1 * (u[2] - 10)")
#' sim_add_block(m, "mux", "mx")
#' sim_connect(m, "temp", "ctrl")
#' sim_connect(m, "ctrl", "mx", output = "y", input = 1)
#' sim_connect(m, "temp", "mx", input = 2)
#' sim_connect(m, "mx", "rate")
#' sim_connect(m, "rate", "temp")
#' r <- sim_run(sim_compile(m), c(0, 10))
#' range(sim_signal(r, "temp"))
#' @export
sim_transition <- function(from, to, guard = "TRUE", after = 0, action = NULL,
                           priority = NA) {
  for (v in list(from, to)) {
    if (!is.character(v) || length(v) != 1L) sim_abort("`from` and `to` must be state names.")
  }
  check_number(after, "`after`")
  if (after < 0) sim_abort("`after` must be non-negative.")
  if (!is.function(guard)) parse_safe(guard, "Transition guard")
  if (!is.null(action) && !is.function(action)) parse_safe(action, "Transition action")
  structure(list(from = from, to = to, guard = guard, after = as.numeric(after),
                 action = action, priority = as.numeric(priority)),
            class = "sim_transition")
}

#' Check a state machine definition
#'
#' Static checks of a state machine: unknown states, ambiguous priorities,
#' states that cannot be reached from the initial state and states without
#' outgoing transitions. Reachability ignores guards, so a state reported as
#' reachable may still never be entered.
#'
#' @param states Character vector of state names.
#' @param transitions List of [sim_transition()] objects.
#' @param initial Initial state.
#' @return A list with `unreachable` and `absorbing` state names. Errors are
#'   raised for structural problems.
#' @examples
#' sim_check_state_machine(c("a", "b", "c"),
#'   list(sim_transition("a", "b", "u[1] > 0")), initial = "a")
#' @export
sim_check_state_machine <- function(states, transitions, initial = states[1]) {
  if (!is.character(states) || !length(states) || anyDuplicated(states) || anyNA(states)) {
    sim_abort("`states` must be unique state names.")
  }
  if (!is.character(initial) || length(initial) != 1L || !initial %in% states) {
    sim_abort("`initial` must be one of the states.")
  }
  if (!is.list(transitions)) sim_abort("`transitions` must be a list of sim_transition().")
  for (tr in transitions) {
    if (!inherits(tr, "sim_transition")) sim_abort("`transitions` must contain sim_transition() objects.")
    bad <- setdiff(c(tr$from, tr$to), states)
    if (length(bad)) sim_abort("Transition refers to unknown state(s): ", paste(bad, collapse = ", "), ".")
  }
  from <- vapply(transitions, `[[`, "", "from")
  pri <- vapply(transitions, `[[`, 0, "priority")
  for (s in unique(from)) {
    p <- pri[from == s & !is.na(pri)]
    if (anyDuplicated(p)) {
      sim_abort("State '", s, "' has transitions with the same priority; ",
                "the order would be ambiguous.")
    }
  }
  to <- vapply(transitions, `[[`, "", "to")
  seen <- initial
  frontier <- initial
  while (length(frontier)) {
    nxt <- setdiff(to[from %in% frontier], seen)
    seen <- c(seen, nxt)
    frontier <- nxt
  }
  list(unreachable = setdiff(states, seen), absorbing = setdiff(states, from))
}

compile_action <- function(a, what, scope = list()) {
  if (is.null(a)) return(NULL)
  if (is.function(a)) {
    return(function(data, u, t) align_data(a(data, u, t), data, what))
  }
  e <- parse_safe(a, what)
  function(data, u, t) {
    env <- safe_frame(c(scope, as.list(data), list(u = u, t = t)))
    eval(e, env)
    out <- vapply(names(data), function(n) {
      v <- get(n, envir = env)
      if (!is.numeric(v) && !is.logical(v) || length(v) != 1L || !is.finite(v)) {
        sim_abort(what, " left data variable '", n, "' as ", deparse_short(v),
                  "; data variables hold one finite number.")
      }
      as.numeric(v)
    }, 0)
    names(out) <- names(data)
    out
  }
}

# The data vector returned by an action function: unnamed values in the
# order of the data variables, or named values matched by name.
align_data <- function(v, data, what) {
  if (!is.numeric(v) && !is.logical(v)) {
    sim_abort(what, " must return a numeric data vector, not a ", class(v)[1], ".")
  }
  nm <- names(v)
  if (!is.null(nm)) {
    if (anyDuplicated(nm) || !setequal(nm, names(data)) || length(v) != length(data)) {
      sim_abort(what, " returned names {", paste(nm, collapse = ", "), "}; expected ",
                "exactly the data variables {", paste(names(data), collapse = ", "), "}.")
    }
    v <- v[names(data)]
  } else if (length(v) != length(data)) {
    sim_abort(what, " returned ", length(v), " value(s); there are ", length(data),
              " data variable(s).")
  }
  v <- as.numeric(v)
  if (any(!is.finite(v))) sim_abort(what, " returned a non-finite data value.")
  names(v) <- names(data)
  v
}

deparse_short <- function(v) {
  s <- paste(deparse(v, width.cutoff = 60L), collapse = " ")
  if (nchar(s) > 60L) paste0(substr(s, 1L, 57L), "...") else s
}

compile_guard <- function(g, what, scope = list()) {
  if (is.function(g)) return(g)
  e <- parse_safe(g, what)
  function(u, t, elapsed, data) {
    v <- eval_safe(e, c(scope, as.list(data), list(u = u, t = t, elapsed = elapsed)))
    if (length(v) != 1L || is.na(v)) sim_abort(what, " must evaluate to TRUE or FALSE.")
    isTRUE(as.logical(v))
  }
}

register_state_machine_block <- function() {
  sim_register_block_type(
    "state_machine", category = "logic", overwrite = TRUE,
    description = "Finite state machine with guarded and timed transitions",
    params = list(states = NULL, initial = NA, transitions = list(),
                  entry = list(), exit = list(), during = list(),
                  data = numeric(), outputs = NULL, n_inputs = 1,
                  period = 1, offset = 0),
    required = "states",
    build = function(p) {
      scope <- p$.scope %||% list()
      states <- p$states
      initial <- if (is.na(p$initial[1])) states[1] else p$initial
      sim_check_state_machine(states, p$transitions, initial)
      data0 <- p$data
      if (length(data0)) {
        dn <- names(data0)
        if (!is.numeric(data0) || is.null(dn) || anyNA(dn) || any(!nzchar(dn))) {
          sim_abort("`data` must be a named numeric vector.")
        }
        if (anyDuplicated(dn)) {
          sim_abort("`data` has duplicate names: ", paste(unique(dn[duplicated(dn)]), collapse = ", "), ".")
        }
        reserved <- intersect(dn, c("u", "t", "elapsed", "data"))
        if (length(reserved)) {
          sim_abort("`data` must not use the reserved name(s) ", paste(reserved, collapse = ", "), ".")
        }
        if (any(!is.finite(data0))) sim_abort("`data` values must be finite.")
      }
      data0 <- as.numeric(data0)
      names(data0) <- names(p$data)
      for (kind in c("entry", "exit", "during")) {
        bad <- setdiff(names(p[[kind]]), states)
        if (length(bad)) sim_abort("`", kind, "` refers to unknown state(s): ", paste(bad, collapse = ", "), ".")
      }
      acts <- lapply(c(entry = "entry", exit = "exit", during = "during"), function(kind) {
        lapply(stats::setNames(states, states), function(s) {
          compile_action(p[[kind]][[s]], paste0("`", kind, "` action of '", s, "'"), scope)
        })
      })
      trs <- lapply(seq_along(p$transitions), function(i) {
        tr <- p$transitions[[i]]
        list(from = match(tr$from, states), to = match(tr$to, states), after = tr$after,
             guard = compile_guard(tr$guard, paste0("Guard of transition ", i), scope),
             action = compile_action(tr$action, paste0("Action of transition ", i), scope),
             key = c(if (is.na(tr$priority)) Inf else tr$priority, i))
      })
      ord <- order(vapply(trs, function(z) z$key[1], 0), vapply(trs, function(z) z$key[2], 0))
      trs <- trs[ord]
      outs <- p$outputs
      if (!is.null(outs)) {
        if (!is.list(outs) || !setequal(names(outs), states)) {
          sim_abort("`outputs` must be a named list with one entry per state.")
        }
        outs <- lapply(outs[states], check_numeric, what = "`outputs`")
        if (length(unique(lengths(outs))) != 1L) sim_abort("`outputs` entries must have equal lengths.")
      }
      nin <- p$n_inputs
      if (!is_count(nin)) sim_abort("`n_inputs` must be a positive integer.")
      nvar <- length(data0)
      # outgoing transitions of each state, already in priority order
      by_state <- lapply(seq_along(states), function(s) Filter(function(tr) tr$from == s, trs))
      run <- function(s, entered, data, u, t) {
        v <- unlist(u, use.names = FALSE)
        elapsed <- t - entered
        # only the rounding of t - entered is forgiven, never a fraction of t
        slack <- 8 * .Machine$double.eps * max(1, abs(t), abs(entered))
        for (tr in by_state[[s]]) {
          if (elapsed < tr$after - slack) next
          if (!tr$guard(v, t, t - entered, data)) next
          if (!is.null(acts$exit[[s]])) data <- acts$exit[[s]](data, v, t)
          if (!is.null(tr$action)) data <- tr$action(data, v, t)
          if (!is.null(acts$entry[[tr$to]])) data <- acts$entry[[tr$to]](data, v, t)
          return(list(s = tr$to, entered = t, data = data))
        }
        if (!is.null(acts$during[[s]])) data <- acts$during[[s]](data, v, t)
        list(s = s, entered = entered, data = data)
      }
      unpack <- function(d) {
        s <- d[[1]]
        if (!is.finite(s) || s != round(s) || s < 1 || s > length(states)) {
          sim_abort("state index ", format(s), " is not one of 1..", length(states), ".")
        }
        data <- d[seq_len(nvar) + 2L]
        names(data) <- names(data0)
        list(s = s, entered = d[[2]], data = data)
      }
      make_out <- function(s, data) {
        o <- list(state = s)
        if (!is.null(outs)) o$y <- outs[[s]]
        if (nvar) o$data <- unname(data)
        o
      }
      s0 <- match(initial, states)
      # The entry action of the initial state runs at the first sample hit,
      # with the actual inputs and time (never while the model is built).
      start <- function(st, u, t) {
        if (is.finite(st$entered)) return(st)
        st$entered <- t
        if (!is.null(acts$entry[[st$s]])) {
          st$data <- acts$entry[[st$s]](st$data, unlist(u, use.names = FALSE), t)
        }
        st
      }
      osz <- c(state = 1L)
      if (!is.null(outs)) osz <- c(osz, y = length(outs[[1]]))
      if (nvar) osz <- c(osz, data = nvar)
      list(inputs = rep(NA, nin), outputs = osz, feedthrough = TRUE,
           nd = 2L + nvar, d0 = c(s0, -Inf, data0),
           sample = c(p$period, p$offset),
           output = function(t, x, d, u) {
             st <- start(unpack(d), u, t)
             r <- run(st$s, st$entered, st$data, u, t)
             make_out(r$s, r$data)
           },
           update = function(t, x, d, u) {
             st <- start(unpack(d), u, t)
             r <- run(st$s, st$entered, st$data, u, t)
             c(r$s, r$entered, r$data)
           },
           # entry time -Inf means "not entered yet": the one non-finite
           # value this block's state may hold
           nonfinite_ok = 2L,
           held0 = function(size, d0) {
             st <- unpack(d0)
             make_out(st$s, st$data)
           })
    })
}
