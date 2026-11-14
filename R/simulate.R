root_methods <- c("lsoda", "lsode", "lsodes", "lsodar")

#' Simulate a compiled plan
#'
#' Runs a [sim_compile()] plan over `tspan` and returns a `sim_result`.
#'
#' **Timeline.** The interval is cut into segments at every sample hit of a
#' discrete block, every scheduled event time and every breakpoint of a
#' source (for example the time of a step). Sample hits are computed as
#' `offset + k * period` with integer `k`, never by accumulating `period`,
#' so they do not drift. Hits closer than `1e-10 * max(1, |t|)` are treated
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
#' @param hmax Optional maximum step size.
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
  if (method == "radau") {
    sim_abort("method = \"radau\" is not supported: it cannot be told to stop at ",
              "segment boundaries, so it would evaluate inputs past discontinuities.")
  }
  if (nzc && !method %in% root_methods) {
    sim_abort("This model has zero-crossing functions; use method = one of ",
              paste0('"', root_methods, '"', collapse = ", "), ".")
  }

  x <- plan$x0
  d <- plan$d0
  held <- plan$held0
  if (!is.null(initial)) {
    st <- apply_initial(plan, initial, x, d)
    x <- st$x; d <- st$d
  }

  # Timeline: boundaries, sample hits and scheduled events ---------------
  rates <- plan$rates
  ticks <- lapply(seq_along(rates$period), function(r) {
    P <- rates$period[r]; o <- rates$offset[r]
    k <- seq.int(ceiling((t0 - o) / P - 1e-9), floor((tf - o) / P + 1e-9))
    o + k * P
  })
  sched <- Filter(function(e) e$trigger$kind == "scheduled", plan$events)
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

  failed <- function(msg) {
    status <<- "failed"
    message <<- msg
  }

  repeat {
    res <- tryCatch({
      if (at_boundary) {
        for (en in sched_at[[bi]]) {
          e <- plan$events[[en]]
          st <- apply_event_resets(plan, e, t, x, d, held)
          x <- st$x; d <- st$d
          log_event(t, en, "scheduled", "model")
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
      for (b in plan$blocks) {
        if (is.function(b$def$stop_if) &&
            isTRUE(b$def$stop_if(t, x[b$x_idx], d[b$d_idx], block_inputs(b, ev$sig)))) {
          log_event(t, b$name, "stop", b$name)
          stop_now <- TRUE
        }
      }
      if (post_event || any(abs(rec_times - t) <= tol)) add_row(t, x, d, held, if (post_event) "post-event" else "sample")
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
    if (t >= tf - tol) break

    tn <- bounds[bi + 1L]
    inner <- rec_times[rec_times > t + tol & rec_times < tn - tol]
    if (plan$nx == 0L && !dummy) {
      for (tt in inner) add_row(tt, x, d, held, "sample")
      t <- tn; bi <- bi + 1L; at_boundary <- TRUE
      next
    }

    seg <- integrate_segment(plan, t, tn, inner, x, d, held, dummy, method,
                             rtol, atol, maxsteps, zc_tol, nzc, ode_args)
    stats_acc <- stats_acc + seg$stats
    for (k in seq_along(seg$times)) add_row(seg$times[k], seg$x[k, ], d, held, "sample")
    if (!is.null(seg$error)) {
      failed(seg$error)
      break
    }
    if (is.null(seg$root)) {
      x <- seg$x_end
      t <- tn; bi <- bi + 1L; at_boundary <- TRUE
      next
    }

    # Zero crossing ---------------------------------------------------------
    tr <- seg$root$time
    xr <- seg$root$x
    hits <- handle_crossings(plan, tr, xr, d, held, seg$root$which, dummy)
    if (!is.null(hits$error)) {
      failed(hits$error)
      break
    }
    if (hits$fired) {
      add_row(tr, xr, d, held, "pre-event")
      for (h in hits$log) log_event(tr, h$name, h$type, h$source, h$direction)
      nev <- nev + 1L
      post_event <- TRUE
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
    if (on_failure == "error") sim_abort("Simulation failed ", message, ".")
    sim_warn("Simulation failed ", message, ". Returning the partial result.")
  }
  if (status == "event_limit") sim_warn("Simulation stopped: ", message, ".")

  build_result(plan, rows, rec, ev_log, status, message, stats_acc,
               list(tspan = tspan, method = method, rtol = rtol, atol = atol,
                    hmax = hmax, maxsteps = maxsteps, max_events = max_events),
               x0 = plan$x0, d0 = plan$d0, initial = initial)
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
  func <- function(tt, y, parms) {
    sig <- evaluate(plan, tt, xs(y), d, held)$sig
    list(if (dummy) 0 else derivatives(plan, tt, y, d, sig))
  }
  args <- c(list(y = y0, times = c(t, inner, tn), func = func, parms = NULL,
                 method = method, rtol = rtol, atol = atol, maxsteps = maxsteps,
                 tcrit = tn),
            ode_args)
  if (nzc) {
    g0 <- zero_crossings(plan, t, x, d, evaluate(plan, t, x, d, held)$sig)
    off <- ifelse(abs(g0) <= zc_tol, g0, 0)
    args$rootfunc <- function(tt, y, parms) {
      xx <- xs(y)
      zero_crossings(plan, tt, xx, d, evaluate(plan, tt, xx, d, held)$sig) - off
    }
  }
  empty <- list(times = numeric(), x = matrix(numeric(), 0L, length(x)),
                stats = c(segments = 1, steps = 0, rhs_evaluations = 0, jacobians = 0, roots = 0))
  out <- tryCatch(suppressWarnings(do.call(deSolve::ode, args)),
                  error = function(e) e)
  if (inherits(out, "error")) {
    empty$error <- paste0("in [", format(t), ", ", format(tn), "]: ", conditionMessage(out))
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
  res <- list(times = tt[keep], x = ymat[keep, , drop = FALSE], stats = stats)
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
  }
  res
}

# Decide which located crossings are real events (direction check) and
# apply their handlers. Returns the new state.
handle_crossings <- function(plan, t, x, d, held, which, dummy) {
  out <- tryCatch({
    sig <- evaluate(plan, t, x, d, held)$sig
    g <- zero_crossings(plan, t, x, d, sig)
    # time derivative of g along the trajectory, by a forward difference
    h <- 1e-7 * max(1, abs(t))
    f <- if (plan$nx) derivatives(plan, t, x, d, sig) else numeric()
    xh <- x + h * f
    gh <- zero_crossings(plan, t + h, xh, d, evaluate(plan, t + h, xh, d, held)$sig)
    slope <- (gh - g) / h
    dir_of <- function(k) if (slope[k] > 0) "rising" else if (slope[k] < 0) "falling" else "either"
    log <- list()
    fired <- FALSE
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
      st <- apply_event_resets(plan, e, t, xnew, dnew, held, sig)
      xnew <- st$x; dnew <- st$d
      if (e$terminal) terminal <- TRUE
      fired <- TRUE
      log[[length(log) + 1L]] <- list(name = e$name, type = "zero_crossing",
                                      source = e$trigger$signal, direction = dir)
    }
    list(fired = fired, terminal = terminal, x = xnew, d = dnew, log = log)
  }, error = function(e) list(error = paste0("at t = ", format(t), ": ", conditionMessage(e))))
  out
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
