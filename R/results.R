# Recompute every recorded signal from accepted states. Stops at the first
# row that cannot be evaluated and reports it, so no value is ever invented.
record_signals <- function(plan, rows) {
  n <- length(rows)
  S <- matrix(NA_real_, n, plan$nsig, dimnames = list(NULL, plan$signal_names))
  loop_res <- NULL
  for (k in seq_len(n)) {
    r <- rows[[k]]
    ev <- tryCatch(evaluate(plan, r$t, r$x, r$d, r$held), error = function(e) e)
    if (inherits(ev, "error")) {
      return(list(S = S[seq_len(k - 1L), , drop = FALSE], loop = loop_res, ok = k - 1L,
                  error = paste0("at t = ", format(r$t), ": ", conditionMessage(ev))))
    }
    S[k, ] <- ev$sig
    if (!is.null(ev$loop)) loop_res <- rbind(loop_res, cbind(time = r$t, ev$loop))
  }
  list(S = S, loop = loop_res, ok = n, error = NULL)
}

build_result <- function(plan, rows, rec, ev_log, status, message, stats, solver,
                         x0, d0, initial, warnings = character()) {
  n <- length(rows)
  time <- vapply(rows, `[[`, 0, "t")
  X <- matrix(as.numeric(unlist(lapply(rows, `[[`, "x"))), n, plan$nx, byrow = TRUE,
              dimnames = list(NULL, plan$state_names))
  D <- matrix(as.numeric(unlist(lapply(rows, `[[`, "d"))), n, plan$nd, byrow = TRUE,
              dimnames = list(NULL, plan$dstate_names))
  S <- rec$S
  loop_res <- rec$loop
  events <- if (length(ev_log)) do.call(rbind, ev_log) else
    data.frame(time = numeric(), event = character(), type = character(),
               source = character(), direction = character(),
               stringsAsFactors = FALSE)
  sampled <- unlist(lapply(plan$blocks, function(b) {
    if (b$discrete) plan$signal_names[unlist(b$out_idx)]
  }))
  structure(list(
    time = time,
    kind = vapply(rows, `[[`, "", "kind"),
    states = X, discrete = D, signals = S,
    signal_of = plan$signal_of, sampled = sampled,
    scopes = names(Filter(function(b) b$type == "scope", plan$blocks)),
    dynamic = names(Filter(function(b) {
      (length(b$x_idx) || b$discrete) && length(b$out_idx)
    }, plan$blocks)),
    events = events,
    status = status, success = status %in% c("success", "terminated"),
    message = message,
    diagnostics = list(status = status, message = message, solver = solver,
                       stats = stats, loop_residuals = loop_res,
                       n_events = nrow(events), warnings = warnings),
    metadata = list(model = plan$name, params = plan$params,
                    initial_continuous = stats::setNames(x0, plan$state_names),
                    initial_discrete = stats::setNames(d0, plan$dstate_names),
                    initial_override = initial,
                    simcompR = pkg_version("simcompR"),
                    deSolve = pkg_version("deSolve"),
                    R = R.version.string, created = Sys.time())
  ), class = "sim_result")
}

check_result <- function(x) {
  if (!inherits(x, "sim_result")) sim_abort("`result` must be a sim_result.")
  invisible(x)
}

#' Extract simulation output
#'
#' `sim_signal()` returns recorded signals by block name. A block with one
#' output port is named by the block name; other ports are `block.port`.
#' Subsystem instances are named by the instance name, inner blocks by
#' `instance/block`, and a `scope` block by its own name (its input
#' signal). `sim_states()` and `sim_discrete_states()` return the state
#' trajectories; `sim_events()` the event log; `sim_times()` the recorded
#' times.
#'
#' Rows follow `result$time`. At an event the time appears twice: the row
#' before and the row after the reset.
#'
#' @param result A `sim_result`.
#' @param names Signal names. `NULL` returns all signals.
#' @param drop Drop dimensions of extent one, so a single scalar signal is
#'   returned as a vector.
#' @return A numeric vector or matrix (signals and states) or a data frame
#'   (events).
#' @examples
#' m <- sim_model()
#' sim_add_block(m, "sine", "u")
#' sim_add_block(m, "integrator", "x", initial = -1)
#' sim_connect(m, "u", "x")
#' r <- sim_run(sim_compile(m), c(0, pi), times = c(0, pi / 2, pi))
#' sim_times(r)
#' sim_signal(r, "x")      # -cos(t)
#' sim_states(r)
#' @export
sim_signal <- function(result, names = NULL, drop = TRUE) {
  check_result(result)
  if (is.null(names)) return(result$signals)
  miss <- setdiff(names, base::names(result$signal_of))
  if (length(miss)) {
    sim_abort("Unknown signal(s): ", paste(miss, collapse = ", "), ". Available: ",
              paste(base::names(result$signal_of), collapse = ", "), ".")
  }
  cols <- lapply(names, function(nm) {
    idx <- result$signal_of[[nm]]
    m <- result$signals[, idx, drop = FALSE]
    colnames(m) <- element_names(nm, length(idx))
    m
  })
  out <- do.call(cbind, cols)
  if (drop) drop(out) else out
}

#' @rdname sim_signal
#' @export
sim_states <- function(result) {
  check_result(result)
  result$states
}

#' @rdname sim_signal
#' @export
sim_discrete_states <- function(result) {
  check_result(result)
  result$discrete
}

#' @rdname sim_signal
#' @export
sim_events <- function(result) {
  check_result(result)
  result$events
}

#' @rdname sim_signal
#' @export
sim_times <- function(result) {
  check_result(result)
  result$time
}

#' @export
as.data.frame.sim_result <- function(x, row.names = NULL, optional = FALSE, ...) {
  data.frame(time = x$time, x$signals, x$states, check.names = FALSE,
             row.names = row.names)
}

#' @export
print.sim_result <- function(x, ...) {
  cat("<sim_result '", x$metadata$model, "'> status: ", x$status, "\n", sep = "")
  if (nzchar(x$message)) cat("  ", x$message, "\n", sep = "")
  if (length(x$time)) {
    cat("  ", length(x$time), " recorded time(s) in [", format(min(x$time)), ", ",
        format(max(x$time)), "]\n", sep = "")
  }
  cat("  ", ncol(x$signals), " signal element(s), ", ncol(x$states),
      " continuous and ", ncol(x$discrete), " discrete state(s), ",
      nrow(x$events), " event(s)\n", sep = "")
  invisible(x)
}

#' @export
summary.sim_result <- function(object, ...) {
  out <- list(
    status = object$status, message = object$message,
    time = if (length(object$time)) range(object$time) else c(NA, NA),
    rows = length(object$time), events = table(object$events$type),
    solver = object$diagnostics$solver, stats = object$diagnostics$stats,
    final_states = if (length(object$time) && ncol(object$states)) object$states[length(object$time), ] else numeric(),
    signal_ranges = if (length(object$time)) apply(object$signals, 2, range, na.rm = TRUE) else NULL
  )
  class(out) <- "summary.sim_result"
  out
}

#' @export
print.summary.sim_result <- function(x, ...) {
  cat("Simulation status: ", x$status, if (nzchar(x$message)) paste0(" (", x$message, ")"), "\n", sep = "")
  cat("Time range: [", format(x$time[1]), ", ", format(x$time[2]), "], ", x$rows, " row(s)\n", sep = "")
  cat("Solver: ", x$solver$method, ", rtol = ", x$solver$rtol, ", atol = ", x$solver$atol, "\n", sep = "")
  st <- x$stats
  cat("Segments: ", st[["segments"]], ", steps: ", st[["steps"]], ", RHS evaluations: ",
      st[["rhs_evaluations"]], ", roots located: ", st[["roots"]], "\n", sep = "")
  if (length(x$events)) {
    cat("Events: ", paste(names(x$events), x$events, sep = " = ", collapse = ", "), "\n", sep = "")
  }
  if (length(x$final_states)) {
    cat("Final continuous states:\n")
    print(x$final_states)
  }
  invisible(x)
}

#' Plot simulation results
#'
#' Draws one panel per signal against time with base graphics. Signals
#' produced by discrete blocks are drawn as steps (held between sample
#' hits); continuous signals as lines. Event times are marked with dotted
#' vertical lines. By default the inputs of `scope` blocks are shown, or,
#' without scopes, the outputs of continuous and discrete dynamic blocks.
#'
#' @param x A `sim_result`.
#' @param signals Signal names to plot.
#' @param events Mark event times.
#' @param ... Passed to [graphics::plot()].
#' @return `x`, invisibly.
#' @examples
#' m <- sim_model()
#' sim_add_block(m, "sine", "u")
#' sim_add_block(m, "zero_order_hold", "zoh", period = 0.5)
#' sim_connect(m, "u", "zoh")
#' r <- sim_run(sim_compile(m), c(0, 6))
#' plot(r, c("u", "zoh"))
#' @export
plot.sim_result <- function(x, signals = NULL, events = TRUE, ...) {
  if (is.null(signals)) signals <- default_plot_signals(x)
  if (!length(signals)) sim_abort("Nothing to plot; give `signals`.")
  n <- length(signals)
  op <- graphics::par(mfrow = c(n, 1L), mar = c(3, 4, 1.5, 1))
  on.exit(graphics::par(op))
  for (nm in signals) {
    y <- sim_signal(x, nm, drop = FALSE)
    held <- all(colnames(x$signals)[x$signal_of[[nm]]] %in% x$sampled)
    graphics::matplot(x$time, y, type = if (held) "s" else "l", lty = 1,
                      col = seq_len(ncol(y)), xlab = "", ylab = nm, main = "", ...)
    if (events && nrow(x$events)) {
      graphics::abline(v = unique(x$events$time), lty = 3, col = "grey50")
    }
  }
  invisible(x)
}

default_plot_signals <- function(x) {
  if (length(x$scopes)) return(x$scopes)
  if (length(x$dynamic)) {
    # a block with several output ports is plotted through its port signals
    sigs <- unlist(lapply(x$dynamic, function(b) {
      nms <- names(x$signal_of)
      if (!is.null(x$signal_of[[b]])) b else nms[startsWith(nms, paste0(b, "."))]
    }))
    return(utils::head(sigs, 6L))
  }
  utils::head(names(x$signal_of), 6L)
}
