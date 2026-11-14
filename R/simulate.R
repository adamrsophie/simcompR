root_methods <- c("lsoda", "lsode", "lsodes", "lsodar")

#' Simulate a compiled plan
#'
#' Runs a [sim_compile()] plan over `tspan` and returns a `sim_result`.
#'
#' **Timeline.** The interval is cut into segments at every sample hit of a
#' discrete block, every scheduled event time and every breakpoint of a
#' source (for example the time of a step). Sample hits are computed as
#' `offset + k * period` with integer `k`, never by accumulating `period`,
#' so they do not drift. Hits closer than `1e-12 * max(1, |t|)` are treated
#' as simultaneous.
#'
#' **At a segment boundary** the simulator, in this order:
#' 1. applies scheduled events (in priority order);
#' 2. evaluates the outputs of all blocks with a sample hit, in execution
#'    order, and holds them;
#' 3. updates the discrete states of those blocks, using the inputs from
#'    step 2. Outputs are computed before states are updated, so a unit
#'    delay produces `y[k] = u[k-1]`;
#' 4. refreshes the modes of switching blocks from their current inputs.
#'
#' **Inside a segment** discrete states and held outputs are constant and
#' the continuous states are integrated with [deSolve::ode()]. The
#' derivative function has no side effects. When zero-crossing functions
#' exist (from [sim_add_event()] or from `switch`, `relay` and `stop`
#' blocks) the root finder stops the solver at the crossing; the event is
#' applied and integration restarts there. A crossing that ends exactly at
#' zero is not detected again at the restart.
#'
#' **Recording.** Signals are computed after the solve from the accepted
#' states, never inside the derivative function. A value recorded at a
#' sample hit is the value after the hit (right-continuous). At an event
#' the result has two rows with the same time: the state before and after
#' the reset.
#'
#' **Failure.** When the solver fails, a block signals an error or the
#' event limit is reached, the result is marked with a status other than
#' `"success"` and contains only the times actually reached; no values are
#' extrapolated. With `on_failure = "error"` (the default for solver and
#' block failures) an error is raised instead.
#'
#' @param plan A compiled `sim_plan`.
#' @param tspan Numeric `c(start, end)`.
#' @param times Output times within `tspan`. Defaults to 201 equally spaced
#'   points. Solver output times are independent of the internal steps
#'   taken by the integrator.
#' @param method Integration method passed to [deSolve::ode()]. Models with
#'   zero crossings need one of `"lsoda"`, `"lsode"`, `"lsodes"` or
#'   `"lsodar"`.
#' @param rtol,atol Relative and absolute tolerances.
#' @param hmax Optional maximum step size. When the model has zero-crossing
#'   functions and `hmax` is not given, it is `(tf - t0) / 100`, because the
#'   root finder cannot see two crossings inside one step; give a smaller
#'   value for signals that cross faster.
#' @param maxsteps Maximum number of internal steps per segment.
#' @param initial Optional named list overriding initial states, keyed by
#'   block name.
#' @param include_ticks Also record a row at every sample hit, so held
#'   signals are captured at their changes.
#' @param max_events Maximum number of zero-crossing events before the run
#'   is stopped with status `"event_limit"`. Guards against chattering and
#'   Zeno behavior.
#' @param zc_tol A zero-crossing function whose value at the start of a
#'   segment is within `zc_tol` of zero is treated as exactly zero there,
#'   which prevents retriggering the same root.
#' @param on_failure `"error"` raises an error when the simulation fails;
#'   `"warn"` returns the partial result with a warning.
#' @param ... Further arguments for [deSolve::ode()].
#' @return A `sim_result`. See [sim_signal()], [sim_states()],
#'   [sim_events()], [sim_diagnostics()] and `plot()`.
#' @examples
#' m <- sim_model("feedback")
#' sim_add_block(m, "constant", "r", value = 1)
#' sim_add_block(m, "subtract", "e")
#' sim_add_block(m, "integrator", "x")
#' sim_connect(m, "r", "e", input = 1)
#' sim_connect(m, "x", "e", input = 2)
#' sim_connect(m, "e", "x")
#' r <- sim_run(sim_compile(m), tspan = c(0, 5), times = 0:5)
#' cbind(numeric = sim_signal(r, "x"), exact = 1 - exp(-(0:5)))
#' @export
sim_run <- function(plan, tspan, times = NULL, method = "lsoda",
                    rtol = 1e-6, atol = 1e-8, hmax = NULL, maxsteps = 5000,
                    initial = NULL, include_ticks = TRUE, max_events = 1000,
                    zc_tol = 1e-9, on_failure = c("error", "warn"), ...) {
  check_plan(plan)
  on_failure <- match.arg(on_failure)
  tspan <- check_numeric(tspan, "`tspan`")
  if (length(tspan) != 2L || tspan[2] <= tspan[1]) {
    sim_abort("`tspan` must be c(start, end) with end > start.")
  }
  t0 <- tspan[1]; tf <- tspan[2]
  tol <- time_tol(c(t0, tf))
  if (is.null(times)) times <- seq(t0, tf, length.out = 201L)
  times <- check_numeric(times, "`times`")
  if (any(times < t0 - tol | times > tf + tol)) sim_abort("`times` must lie within `tspan`.")
  times <- sort(unique(pmin(pmax(times, t0), tf)))
  check_number(rtol, "`rtol`"); check_number(atol, "`atol`")
  if (!is_count(max_events)) sim_abort("`max_events` must be a positive integer.")
  n_event_zc <- sum(vapply(plan$events, function(e) e$trigger$kind == "zero_crossing", TRUE))
  nzc <- plan$zc$n + n_event_zc
  if (!is.character(method) || length(method) != 1L) sim_abort("`method` must be a string.")
  if (method %in% c("radau", "iteration")) {
    sim_abort("method = \"", method, "\" is not supported: it cannot be told to stop at ",
              "segment boundaries.")
  }
  # kinks of memoryless blocks are optional crossing functions
  optional <- c(plan$zc$optional %||% rep(FALSE, plan$zc$n), rep(FALSE, n_event_zc))
  n_strict <- sum(!optional)
  if (nzc && !method %in% root_methods) {
    if (n_strict) {
      sim_abort("This model has zero-crossing functions; use method = one of ",
                paste0('"', root_methods, '"', collapse = ", "), ".")
    }
    nzc <- 0L
  }

  x <- plan$x0
  d <- plan$d0
  held <- plan$held0
  if (!is.null(initial)) {
    st <- apply_initial(plan, initial, x, d)
    x <- st$x; d <- st$d
    for (nm in names(initial)) held <- set_held0(plan$blocks[[nm]], d[plan$blocks[[nm]]$d_idx], held)
  }
  x_init <- x
  d_init <- d

  # Timeline: boundaries, sample hits and scheduled events ---------------
  rates <- plan$rates
  ticks <- lapply(seq_along(rates$period), function(r) {
    P <- rates$period[r]; o <- rates$offset[r]
    k <- seq.int(ceiling((t0 - o) / P - 1e-9), floor((tf - o) / P + 1e-9))
    o + k * P
  })
  sched <- Filter(function(e) e$trigger$kind == "scheduled", plan$events)
  outside <- vapply(sched, function(e) any(e$trigger$times < t0 - tol | e$trigger$times > tf + tol), TRUE)
  if (any(outside)) {
    sim_warn("Scheduled event(s) ", paste(names(sched)[outside], collapse = ", "),
             " have times outside tspan [", format(t0), ", ", format(tf), "]; those times are skipped.")
  }
  brk <- unlist(lapply(plan$blocks, function(b) {
    if (is.function(b$def$breakpoints)) b$def$breakpoints(t0, tf)
  }))
  brk <- brk[brk > t0 & brk < tf]
  bounds <- merge_times(c(t0, tf, unlist(ticks),
                          unlist(lapply(sched, function(e) e$trigger$times)), brk), tol)
  bounds <- bounds[bounds >= t0 - tol & bounds <= tf + tol]
  bounds[1] <- t0
  bounds[length(bounds)] <- tf
  near <- function(v) match_tol(v, bounds, tol)
  tick_at <- vector("list", length(bounds))
  for (r in seq_along(ticks)) {
    for (i in stats::na.omit(near(ticks[[r]]))) tick_at[[i]] <- c(tick_at[[i]], r)
  }
  sched_at <- vector("list", length(bounds))
  for (e in names(sched)) {
    for (i in stats::na.omit(near(sched[[e]]$trigger$times))) sched_at[[i]] <- c(sched_at[[i]], e)
  }
  rec_times <- times
  if (include_ticks && length(ticks)) {
    tk <- unlist(ticks)
    tk <- tk[tk >= t0 - tol & tk <= tf + tol]
    extra <- unique(tk[is.na(match_tol(tk, rec_times, tol))])
    rec_times <- sort(c(rec_times, extra))
  }

  ode_args <- list(...)
  # With crossing functions, cap the step so that two crossings inside one
  # step (which the root finder cannot see) need a fast signal.
  if (is.null(hmax) && nzc && n_strict) hmax <- (tf - t0) / 100
  if (!is.null(hmax)) ode_args$hmax <- hmax
  dummy <- plan$nx == 0L && nzc > 0L

  rows <- list()
  add_row <- function(t, x, d, held, kind) {
    rows[[length(rows) + 1L]] <<- list(t = t, x = x, d = d, held = held, kind = kind)
  }
  ev_log <- list()
  log_event <- function(t, name, type, source, direction = NA_character_) {
    ev_log[[length(ev_log) + 1L]] <<- data.frame(
      time = t, event = name, type = type, source = source,
      direction = direction, stringsAsFactors = FALSE)
  }
  stats_acc <- c(segments = 0, steps = 0, rhs_evaluations = 0, jacobians = 0,
                 roots = 0)
  status <- "success"
  message <- ""
  warnings <- character()
  t <- t0
  bi <- 1L
  at_boundary <- TRUE
  post_event <- FALSE
  stop_now <- FALSE
  nev <- 0L
  restarts <- 0L
  stuck <- 0L
  fired_t <- NA_real_  # time of the last fired model crossing events
  fired_set <- character()
  g_prev <- NULL      # crossing functions at the end of the last segment
  g_seg0 <- NULL      # and at its start

  failed <- function(msg) {
    status <<- "failed"
    message <<- msg
  }

  repeat {
    res <- tryCatch({
      if (at_boundary) {
        if (length(sched_at[[bi]])) add_row(t, x, d, held, "pre-event")
        for (en in sched_at[[bi]]) {
          e <- plan$events[[en]]
          st <- apply_event_resets(plan, e, t, x, d, held)
          x <- st$x; d <- st$d
          log_event(t, en, "scheduled", "model")
          nev <- nev + 1L
          post_event <- TRUE
          if (e$terminal) stop_now <- TRUE
        }
        if (length(tick_at[[bi]])) {
          ticking <- plan$rates$block_rate %in% tick_at[[bi]] &
            vapply(plan$blocks, `[[`, TRUE, "discrete")
          ev <- evaluate(plan, t, x, d, held, ticking, refresh = TRUE)
          d <- ev$d
          for (i in which(ticking)) {
            b <- plan$blocks[[i]]
            idx <- unlist(b$out_idx)
            held[idx] <- ev$sig[idx]
          }
          dnew <- d
          for (i in which(ticking)) {
            b <- plan$blocks[[i]]
            if (is.function(b$def$update)) {
              v <- as.numeric(b$def$update(t, x[b$x_idx], d[b$d_idx], block_inputs(b, ev$sig)))
              if (length(v) != length(b$d_idx)) {
                sim_abort("Block '", b$name, "' update returned ", length(v),
                          " value(s), expected ", length(b$d_idx), ".")
              }
              dnew[b$d_idx] <- v
            }
          }
          d <- dnew
        }
      }
      ev <- evaluate(plan, t, x, d, held, NULL, refresh = TRUE)
      d <- ev$d
      if (at_boundary && !is.null(g_prev) && !stop_now) {
        bx <- boundary_crossings(plan, t, x, d, held, ev$sig, g_prev, g_seg0, zc_tol)
        if (length(bx$log)) {
          if (!post_event) add_row(t, x, d, held, "pre-event")
          for (h in bx$log) log_event(t, h$name, h$type, h$source, h$direction)
          fired_t <- t
          fired_set <- vapply(bx$log, `[[`, "", "name")
          x <- bx$x; d <- bx$d
          nev <- nev + 1L
          post_event <- TRUE
          if (bx$terminal) stop_now <- TRUE
          ev <- evaluate(plan, t, x, d, held, NULL, refresh = TRUE)
          d <- ev$d
        }
      }
      g_prev <- NULL
      for (b in plan$blocks) {
        if (is.function(b$def$stop_if) &&
            isTRUE(b$def$stop_if(t, x[b$x_idx], d[b$d_idx], block_inputs(b, ev$sig)))) {
          log_event(t, b$name, "stop", b$name)
          stop_now <- TRUE
        }
      }
      if (post_event || stop_now || any(abs(rec_times - t) <= tol)) {
        add_row(t, x, d, held, if (post_event) "post-event" else "sample")
      }
      post_event <- FALSE
      list(x = x, d = d, held = held)
    }, error = function(e) {
      failed(paste0("at t = ", format(t), ": ", conditionMessage(e)))
      NULL
    })
    if (is.null(res)) break
    x <- res$x; d <- res$d; held <- res$held
    if (stop_now) {
      status <- "terminated"
      message <- paste0("stopped by a terminal event at t = ", format(t))
      break
    }
    if (nev >= max_events) {
      status <- "event_limit"
      message <- paste0("reached max_events = ", max_events, " at t = ", format(t),
                        "; the model may be chattering or Zeno")
      break
    }
    if (t >= tf - tol) break

    tn <- bounds[bi + 1L]

    # A watched signal that sits exactly on its level and leaves it in the
    # trigger direction crosses here; the root finder cannot see that.
    if (n_event_zc && nzc) {
      skip <- if (identical(fired_t, t)) fired_set else character()
      dep <- tryCatch(departure_crossings(plan, t, tn, x, d, held, skip, zc_tol),
                      error = function(e) e)
      if (inherits(dep, "error")) {
        failed(paste0("at t = ", format(t), ": ", conditionMessage(dep)))
        break
      }
      if (length(dep$log)) {
        if (!length(rows) || rows[[length(rows)]]$t != t) add_row(t, x, d, held, "pre-event")
        else rows[[length(rows)]]$kind <- if (rows[[length(rows)]]$kind == "post-event")
          "post-event" else "pre-event"
        for (h in dep$log) log_event(t, h$name, h$type, h$source, h$direction)
        fired_t <- t
        fired_set <- c(skip, vapply(dep$log, `[[`, "", "name"))
        x <- dep$x; d <- dep$d
        nev <- nev + 1L
        post_event <- TRUE
        if (dep$terminal) stop_now <- TRUE
        at_boundary <- FALSE
        next
      }
    }
    inner <- rec_times[rec_times > t + tol & rec_times < tn - tol]
    if (plan$nx == 0L && !dummy) {
      for (tt in inner) add_row(tt, x, d, held, "sample")
      t <- tn; bi <- bi + 1L; at_boundary <- TRUE
      next
    }

    seg <- integrate_segment(plan, t, tn, inner, x, d, held, dummy, method,
                             rtol, atol, maxsteps, zc_tol, nzc, ode_args)
    stats_acc <- stats_acc + seg$stats
    warnings <- unique(c(warnings, seg$warnings))
    for (k in seq_along(seg$times)) add_row(seg$times[k], seg$x[k, ], d, held, "sample")
    if (!is.null(seg$error)) {
      failed(seg$error)
      break
    }
    if (is.null(seg$root)) {
      x <- seg$x_end
      g_prev <- seg$g_end
      g_seg0 <- seg$g_start
      t <- tn; bi <- bi + 1L; at_boundary <- TRUE
      next
    }

    # Zero crossing ---------------------------------------------------------
    tr <- seg$root$time
    xr <- seg$root$x
    hits <- handle_crossings(plan, tr, xr, d, held, seg$root$which, dummy,
                             seg$g_start, t, zc_tol)
    if (!is.null(hits$error)) {
      failed(hits$error)
      break
    }
    if (hits$fired) {
      add_row(tr, xr, d, held, "pre-event")
      for (h in hits$log) log_event(tr, h$name, h$type, h$source, h$direction)
      fired_t <- tr
      fired_set <- vapply(hits$log, `[[`, "", "name")
      nev <- nev + 1L
      post_event <- TRUE
    } else if (hits$silent) {
      # a restart at a kink; guard against restarts that do not advance
      stuck <- if (abs(tr - t) <= tol) stuck + 1L else 0L
      if (stuck > 100L) {
        failed(paste0("restarts at kinks do not advance at t = ", format(tr)))
        break
      }
    } else {
      restarts <- restarts + 1L
      if (restarts > 100L * max_events) {
        failed("too many restarts at crossings in the ignored direction")
        break
      }
    }
    x <- hits$x; d <- hits$d
    if (hits$terminal) {
      add_row(tr, x, d, held, "post-event")
      status <- "terminated"
      message <- paste0("stopped by a terminal event at t = ", format(tr))
      break
    }
    if (nev >= max_events) {
      add_row(tr, x, d, held, "post-event")
      status <- "event_limit"
      message <- paste0("reached max_events = ", max_events, " at t = ", format(tr),
                        "; the model may be chattering or Zeno")
      break
    }
    if (abs(tr - tn) <= tol) {
      t <- tn; bi <- bi + 1L; at_boundary <- TRUE
    } else {
      t <- tr; at_boundary <- FALSE
    }
  }

  rec <- record_signals(plan, rows)
  if (!is.null(rec$error)) {
    rows <- rows[seq_len(rec$ok)]
    if (status != "failed") {
      status <- "failed"
      message <- rec$error
    }
  }
  if (status == "failed") {
    message <- sub("[.]+$", "", message)
    if (on_failure == "error") sim_abort("Simulation failed ", message, ".")
    sim_warn("Simulation failed ", message, ". Returning the partial result.")
  }
  if (status == "event_limit") sim_warn("Simulation stopped: ", message, ".")

  build_result(plan, rows, rec, ev_log, status, message, stats_acc,
               list(tspan = tspan, method = method, rtol = rtol, atol = atol,
                    hmax = hmax, maxsteps = maxsteps, max_events = max_events),
               x0 = x_init, d0 = d_init, initial = initial, warnings = warnings)
}

apply_initial <- function(plan, initial, x, d) {
  if (!is.list(initial) || is.null(names(initial))) {
    sim_abort("`initial` must be a named list keyed by block name.")
  }
  for (nm in names(initial)) {
    b <- plan$blocks[[nm]]
    if (is.null(b)) sim_abort("`initial`: unknown block '", nm, "'.")
    v <- check_numeric(initial[[nm]], paste0("`initial$", nm, "`"))
    if (length(b$x_idx)) {
      if (length(v) != length(b$x_idx)) sim_abort("`initial$", nm, "` must have length ", length(b$x_idx), ".")
      x[b$x_idx] <- v
    } else if (length(b$d_idx)) {
      if (length(v) != length(b$d_idx)) sim_abort("`initial$", nm, "` must have length ", length(b$d_idx), ".")
      d[b$d_idx] <- v
    } else {
      sim_abort("`initial`: block '", nm, "' has no state.")
    }
  }
  list(x = x, d = d)
}

integrate_segment <- function(plan, t, tn, inner, x, d, held, dummy, method,
                              rtol, atol, maxsteps, zc_tol, nzc, ode_args) {
  y0 <- if (dummy) 0 else x
  xs <- function(y) if (dummy) numeric() else y
  # Inside a segment, sources with breakpoints take their left limit at the
  # end point: a source that jumps at tn must not leak its next value into
  # this segment when the solver evaluates the right-hand side exactly at
  # tn. All other blocks see the true solver time.
  t_left <- tn - 4 * .Machine$double.eps * max(1, abs(tn))
  # The root finder usually asks for the point the derivative was just
  # computed at; remember the last evaluation so it is not repeated. The
  # cache lives in this segment's frame only and holds a pure function's
  # value, so it changes no result.
  last_t <- NA_real_; last_y <- NULL; last_sig <- NULL
  signals_at <- function(tt, y) {
    if (identical(tt, last_t) && identical(y, last_y)) return(last_sig)
    s <- evaluate(plan, tt, xs(y), d, held, t_left = t_left)$sig
    last_t <<- tt; last_y <<- y; last_sig <<- s
    s
  }
  func <- function(tt, y, parms) {
    sig <- signals_at(tt, y)
    list(if (dummy) 0 else derivatives(plan, tt, y, d, sig))
  }
  args <- c(list(y = y0, times = c(t, inner, tn), func = func, parms = NULL,
                 method = method, rtol = rtol, atol = atol, maxsteps = maxsteps,
                 tcrit = tn),
            ode_args)
  g0 <- NULL
  if (nzc) {
    g0 <- zero_crossings(plan, t, x, d, evaluate(plan, t, x, d, held)$sig)
    at0 <- abs(g0) <= zc_tol
    off <- ifelse(at0, g0, 0)
    args$rootfunc <- function(tt, y, parms) {
      xx <- xs(y)
      v <- zero_crossings(plan, tt, xx, d, signals_at(tt, y)) - off
      # A function that starts at zero and stays exactly there (a constant
      # signal sitting on its threshold) would make the root finder stop with
      # "root too near the initial point". Past the start, report an exact
      # zero as positive; at the start itself keep the zero, so the root
      # finder takes its reference sign from the first step.
      if (tt > t) v[at0 & v == 0] <- 1
      v
    }
  }
  empty <- list(times = numeric(), x = matrix(numeric(), 0L, length(x)),
                stats = c(segments = 1, steps = 0, rhs_evaluations = 0, jacobians = 0, roots = 0))
  warns <- character()
  out <- tryCatch(withCallingHandlers(do.call(deSolve::ode, args), warning = function(w) {
    warns <<- c(warns, conditionMessage(w))
    invokeRestart("muffleWarning")
  }), error = function(e) e)
  if (inherits(out, "error")) {
    empty$error <- paste0("in [", format(t), ", ", format(tn), "]: ",
                          sub("[.]+$", "", conditionMessage(out)))
    empty$warnings <- warns
    return(empty)
  }
  ist <- attr(out, "istate")
  stats <- c(segments = 1, steps = ist[2] %||% NA, rhs_evaluations = ist[3] %||% NA,
             jacobians = ist[4] %||% NA, roots = 0)
  stats[is.na(stats)] <- 0
  tt <- out[, 1]
  ymat <- out[, -1, drop = FALSE]
  if (dummy) ymat <- matrix(numeric(), nrow(out), 0L)
  troot <- attr(out, "troot")
  is_root <- length(troot) > 0L
  last <- nrow(out)
  keep <- which(tt %in% inner & (if (is_root) tt < troot[1] else TRUE))
  res <- list(times = tt[keep], x = ymat[keep, , drop = FALSE], stats = stats,
              warnings = warns)
  if (!all(is.finite(ymat))) {
    bad <- which(!apply(is.finite(ymat), 1L, all))[1]
    res$times <- res$times[res$times < tt[bad]]
    res$x <- res$x[seq_along(res$times), , drop = FALSE]
    res$error <- paste0("the states became NaN or infinite at t = ", format(tt[bad]))
    return(res)
  }
  flag <- ist[1] %||% 2
  if (!is_root && (flag < 0 || abs(tt[last] - tn) > 1e-9 * max(1, abs(tn)))) {
    res$error <- paste0("the integrator stopped at t = ", format(tt[last]),
                        " before reaching ", format(tn), " (deSolve istate ", flag, ")")
    return(res)
  }
  if (is_root) {
    stats["roots"] <- 1
    res$stats <- stats
    iroot <- attr(out, "iroot")
    which_root <- if (length(iroot) == nzc) which(iroot != 0) else as.integer(iroot)
    res$root <- list(time = tt[last], x = ymat[last, ], which = which_root)
  } else {
    res$x_end <- ymat[last, ]
    if (nzc) {
      xe <- xs(res$x_end)
      res$g_end <- zero_crossings(plan, tn, xe, d,
                                  evaluate(plan, tn, xe, d, held, t_left = t_left)$sig)
    }
  }
  res$g_start <- g0
  res
}

# Decide which located crossings are real events (direction check) and
# apply their handlers. Returns the new state.
handle_crossings <- function(plan, t, x, d, held, which, dummy, g_start, t_start,
                             zc_tol) {
  out <- tryCatch({
    sig <- evaluate(plan, t, x, d, held)$sig
    g <- zero_crossings(plan, t, x, d, sig)
    # The solver stops at the first sign change in the segment, so the sign
    # at the segment start gives the direction. When that value was itself
    # zero, look a short step back along the trajectory, inside the segment.
    back <- NULL
    dir_of <- function(k) {
      if (abs(g_start[k]) > zc_tol) return(if (g_start[k] < 0) "rising" else "falling")
      if (is.null(back)) {
        h <- min(1e-7 * max(1, abs(t)), (t - t_start) / 2)
        f <- if (plan$nx) derivatives(plan, t, x, d, sig) else numeric()
        xb <- x - h * f
        back <<- zero_crossings(plan, t - h, xb, d, evaluate(plan, t - h, xb, d, held)$sig)
      }
      if (back[k] < g[k]) "rising" else if (back[k] > g[k]) "falling" else "either"
    }
    log <- list()
    fired <- FALSE
    silent <- FALSE
    terminal <- FALSE
    dnew <- d; xnew <- x
    # block crossings, in execution order
    for (i in plan$zc$blocks) {
      b <- plan$blocks[[i]]
      comps <- plan$zc$block_idx[[i]]
      for (j in seq_along(comps)) {
        k <- comps[j]
        if (!k %in% which) next
        want <- b$def$zc_direction[j]
        dir <- dir_of(k)
        if (want != "either" && dir != want && dir != "either") next
        r <- b$def$on_zc(j, dir, t, xnew[b$x_idx], dnew[b$d_idx], block_inputs(b, sig))
        if (isTRUE(r$silent)) {
          silent <- TRUE
          next
        }
        if (!is.null(r$d)) dnew[b$d_idx] <- r$d
        if (!is.null(r$x)) xnew[b$x_idx] <- r$x
        if (isTRUE(r$terminal)) terminal <- TRUE
        fired <- TRUE
        log[[length(log) + 1L]] <- list(name = b$name, type = if (isTRUE(r$terminal)) "stop" else "block",
                                        source = b$name, direction = dir)
      }
    }
    for (e in plan$events) {
      if (e$trigger$kind != "zero_crossing" || !e$zc_index %in% which) next
      dir <- dir_of(e$zc_index)
      want <- e$trigger$direction
      if (want != "either" && dir != want && dir != "either") next
      # each handler sees the signals left by the handlers before it
      if (fired) sig <- evaluate(plan, t, xnew, dnew, held)$sig
      st <- apply_event_resets(plan, e, t, xnew, dnew, held, sig)
      xnew <- st$x; dnew <- st$d
      if (e$terminal) terminal <- TRUE
      fired <- TRUE
      log[[length(log) + 1L]] <- list(name = e$name, type = "zero_crossing",
                                      source = e$trigger$signal, direction = dir)
    }
    list(fired = fired, silent = silent, terminal = terminal, x = xnew, d = dnew, log = log)
  }, error = function(e) list(error = paste0("at t = ", format(t), ": ", conditionMessage(e))))
  out
}

# Model events whose crossing function changed sign across a segment
# boundary: a source jump, a sample hit or a reset moved the watched signal
# past its level without any continuous crossing for the solver to locate.
# A function that reached zero exactly at the boundary counts as crossing
# there, with the direction it approached from.
boundary_crossings <- function(plan, t, x, d, held, sig, g_end, g_start, zc_tol) {
  g <- zero_crossings(plan, t, x, d, sig)
  log <- list(); terminal <- FALSE; fired <- FALSE
  for (e in plan$events) {
    if (e$trigger$kind != "zero_crossing") next
    k <- e$zc_index
    a <- if (abs(g_end[k]) > zc_tol) sign(g_end[k]) else
      if (abs(g_start[k]) > zc_tol && abs(g[k]) <= zc_tol) sign(g_start[k]) else 0
    b <- if (abs(g[k]) > zc_tol) sign(g[k]) else if (a != 0 && abs(g_end[k]) <= zc_tol) -a else 0
    if (a == 0 || b == 0 || a == b) next
    dir <- if (b > 0) "rising" else "falling"
    want <- e$trigger$direction
    if (want != "either" && dir != want) next
    if (fired) sig <- evaluate(plan, t, x, d, held)$sig
    st <- apply_event_resets(plan, e, t, x, d, held, sig)
    x <- st$x; d <- st$d
    fired <- TRUE
    if (e$terminal) terminal <- TRUE
    log[[length(log) + 1L]] <- list(name = e$name, type = "zero_crossing",
                                    source = e$trigger$signal, direction = dir)
  }
  list(x = x, d = d, log = log, terminal = terminal)
}

# Model events whose function is zero at t (within zc_tol) and moves away
# from zero in the trigger direction. The direction comes from a short step
# along the trajectory inside the coming segment. Events in `skip` fired at
# this instant already and are not fired again.
departure_crossings <- function(plan, t, tn, x, d, held, skip, zc_tol) {
  evs <- Filter(function(e) e$trigger$kind == "zero_crossing" && !e$name %in% skip,
                plan$events)
  if (!length(evs)) return(NULL)
  sig <- evaluate(plan, t, x, d, held)$sig
  g <- zero_crossings(plan, t, x, d, sig)
  at <- vapply(evs, function(e) abs(g[e$zc_index]) <= zc_tol, TRUE)
  if (!any(at)) return(NULL)
  h <- min(1e-7 * max(1, abs(t)), (tn - t) / 4)
  f <- if (plan$nx) derivatives(plan, t, x, d, sig) else numeric()
  xf <- x + h * f
  t_left <- tn - 4 * .Machine$double.eps * max(1, abs(tn))
  gf <- zero_crossings(plan, t + h, xf, d,
                       evaluate(plan, t + h, xf, d, held, t_left = t_left)$sig)
  log <- list(); terminal <- FALSE; fired <- FALSE
  for (e in evs[at]) {
    k <- e$zc_index
    dg <- gf[k] - g[k]
    if (dg == 0) next
    dir <- if (dg > 0) "rising" else "falling"
    want <- e$trigger$direction
    if (want != "either" && dir != want) next
    if (fired) sig <- evaluate(plan, t, x, d, held)$sig
    st <- apply_event_resets(plan, e, t, x, d, held, sig)
    x <- st$x; d <- st$d
    fired <- TRUE
    if (e$terminal) terminal <- TRUE
    log[[length(log) + 1L]] <- list(name = e$name, type = "zero_crossing",
                                    source = e$trigger$signal, direction = dir)
  }
  list(x = x, d = d, log = log, terminal = terminal)
}

apply_event_resets <- function(plan, e, t, x, d, held, sig = NULL) {
  if (!length(e$resets)) return(list(x = x, d = d))
  if (is.null(sig)) sig <- evaluate(plan, t, x, d, held)$sig
  s <- signal_list(plan, sig)
  for (r in e$resets) {
    cur <- if (r$kind == "x") x[r$idx] else d[r$idx]
    v <- as.numeric(r$fun(cur, s, t))
    if (length(v) != length(r$idx) || anyNA(v)) {
      sim_abort("Event '", e$name, "': reset of '", r$target, "' returned ",
                length(v), " value(s), expected ", length(r$idx), ".")
    }
    if (r$kind == "x") x[r$idx] <- v else d[r$idx] <- v
  }
  list(x = x, d = d)
}
