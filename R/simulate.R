root_methods <- c("lsoda", "lsode", "lsodes", "lsodar")

#' Simulate a compiled plan
#'
#' Runs a [sim_compile()] plan over `tspan` and returns a `sim_result`.
#'
#' **Timeline.** The interval is cut into segments at every sample hit of a
#' discrete block, every scheduled event time and every breakpoint of a
#' source (for example the time of a step). Sample hits are computed as
#' `offset + k * period` with integer `k`, never by accumulating `period`,
#' so they do not drift and do not depend on where the time axis starts.
#' Two computed times within 16 units in the last place of the largest
#' time in `tspan` are the same instant; a sample period that is not well
#' above that resolution is an error. The solver integrates each segment
#' on time measured from the segment start, so continuous states keep
#' their accuracy when the time axis starts far from zero.
#'
#' **At an instant** (a segment boundary or a located crossing) the
#' simulator, in this order:
#' 1. runs the handlers of `switch`, `relay` and `stop` blocks whose
#'    crossing was located by the solver;
#' 2. runs the model events due now in increasing `priority`, then by name:
#'    scheduled events, located crossings and crossings caused by a jump
#'    of a signal at a segment boundary (a step, a sample hit, a reset).
#'    Each event sees the states left by the one before;
#' 3. evaluates the outputs of all blocks with a sample hit, in execution
#'    order, holds them, and then updates the discrete states of those
#'    blocks from the same inputs, so a unit delay produces
#'    `y[k] = u[k-1]`; then refreshes the modes of switching blocks;
#' 4. fires, in further layers, the events whose crossing functions were
#'    crossed by steps 2 and 3, comparing each layer with the values
#'    before it, until no new crossing occurs. An event does not fire
#'    twice in a row from the same layer;
#' 5. fires block handlers and events whose function sits on its level
#'    and moves off it in the trigger direction, found by integrating a
#'    short probe (1/1000 of the coming segment), then repeats step 4.
#'    So an event can fire at the initial time, and a function whose first
#'    derivative is zero there is still seen to move.
#'
#' A terminal event ends the run after the instant is settled, so other
#' events of the same instant still run. Each handler execution counts
#' toward `max_events`; a cascade that does not settle ends the run with
#' status `"event_limit"`.
#'
#' **Inside a segment** discrete states and held outputs are constant and
#' the continuous states are integrated with [deSolve::ode()]. The
#' derivative function has no side effects. When zero-crossing functions
#' exist (from [sim_add_event()], from `switch`, `relay` and `stop` blocks,
#' or from the kinks of memoryless blocks such as `saturation`) the root
#' finder stops the solver at the crossing and the instant is handled as
#' above (a kink only restarts the solver). A crossing that ends exactly
#' at zero is not detected again at the restart.
#'
#' **Recording.** Signals are computed after the solve from the accepted
#' states, never inside the derivative function. A value recorded at a
#' sample hit is the value after the hit (right-continuous). At an instant
#' where events fire the result has exactly two rows with the same time:
#' the values before the first event and after the last.
#'
#' **Finite values.** Block outputs, continuous states and discrete states
#' must stay finite. A block output, reset, sample update or block handler
#' that produces `NA`, `NaN` or an infinite value fails the run with a
#' message naming it.
#'
#' **Failure.** When the solver fails, a block signals an error or the
#' event limit is reached, the result is marked with a status other than
#' `"success"` and contains only the times actually reached; no values are
#' extrapolated. With `on_failure = "error"` (the default for solver and
#' block failures) an error is raised instead.
#'
#' **Statistics.** `result$diagnostics$stats` sums, over all segments, the solver
#' steps, right-hand side evaluations and located roots as reported by
#' deSolve. `jacobians` counts Jacobian evaluations and LU decompositions
#' for the ODEPACK solvers (`lsoda`, `lsode`, `lsodes`, `lsodar`, `vode`,
#' `bdf`, `adams`) and is `NA` for solvers that do not report it.
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
#' @param hmax Optional maximum step size. When the model has events or
#'   switching blocks and `hmax` is not given, it is `(tf - t0) / 100`
#'   (kinks alone do not set it), because the
#'   root finder cannot see two crossings inside one step; give a smaller
#'   value for signals that cross faster.
#' @param maxsteps Maximum number of internal steps per segment.
#' @param initial Optional named list overriding initial states, keyed by
#'   block name.
#' @param include_ticks Also record a row at every sample hit, so held
#'   signals are captured at their changes.
#' @param max_events Maximum number of event and block handler executions
#'   before the run is stopped with status `"event_limit"`. Guards against
#'   chattering and Zeno behavior.
#' @param zc_tol A zero-crossing function within `zc_tol` of zero is on its
#'   level: arriving there counts as a crossing, the same root is not
#'   located again at a restart, and leaving the level fires the event
#'   only in its trigger direction (step 5 above).
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
  small <- rates$period <= 4 * tol
  if (any(small)) {
    sim_abort("Sample period ", format(rates$period[small][1]), " is too short to be ",
              "resolved at times near ", format(max(abs(tspan))), " in double precision; ",
              "shift the time axis closer to zero or use a longer period.")
  }
  # Hit k of a rate is offset + k * period. The first and last indices are
  # adjusted by one where rounding of (t - offset) / period would drop a hit
  # at t0 or tf, or add one outside [t0, tf].
  ticks <- lapply(seq_along(rates$period), function(r) {
    P <- rates$period[r]; o <- rates$offset[r]
    lo <- ceiling((t0 - o) / P); hi <- floor((tf - o) / P)
    if (o + (lo - 1) * P >= t0 - tol) lo <- lo - 1
    if (o + lo * P < t0 - tol) lo <- lo + 1
    if (o + (hi + 1) * P <= tf + tol) hi <- hi + 1
    if (o + hi * P > tf + tol) hi <- hi - 1
    if (hi < lo) return(numeric())
    o + seq.int(lo, hi) * P
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
  # without continuous states nothing is integrated, so kinks need no stop
  if (plan$nx == 0L && !n_strict) nzc <- 0L
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
  stats_acc <- c(segments = 0, steps = 0, rhs_evaluations = 0, jacobians = NA_real_,
                 roots = 0)
  status <- "success"
  message <- ""
  warnings <- character()
  t <- t0
  bi <- 1L
  at_boundary <- TRUE
  nev <- 0L
  restarts <- 0L
  stuck <- 0L
  root <- NULL     # crossings located by the solver at t
  g_end <- NULL    # crossing functions at the left limit of t
  g_seg0 <- NULL   # and at the start of the segment that ended at t
  zc_events <- Filter(function(e) e$trigger$kind == "zero_crossing", plan$events)
  strict_zc <- !optional[seq_len(plan$zc$n)]
  d_inf_ok <- nonfinite_allowed(plan)
  ode_probe <- c(list(method = method, rtol = rtol, atol = atol, maxsteps = maxsteps), ode_args)
  ode_probe$hmax <- NULL

  failed <- function(msg) {
    status <<- "failed"
    message <<- msg
  }

  # Everything that happens at the instant t, in this order:
  #   1. block handlers of crossings located by the solver;
  #   2. model events due now, in increasing priority then name: scheduled
  #      events, located crossings and crossings caused by a jump of a
  #      signal at t; each sees the state left by the one before;
  #   3. sample hits (outputs, then updates), then the modes of switching
  #      blocks;
  #   4. crossings caused by steps 2 and 3, in further layers, each layer
  #      compared with the values before it, until no new crossing occurs;
  #   5. crossing functions that sit on their level and move off it in the
  #      trigger direction (found with a short probe integration), then
  #      step 4 again.
  # A terminal event ends the run after the instant is settled.
  instant <- function() {
    fired_ev <- character()
    fired_blk <- character()
    out <- list(terminal = FALSE, limit = FALSE, fired = FALSE, silent = FALSE)
    pre_done <- FALSE
    note_pre <- function() {
      if (!pre_done) add_row(t, x, d, held, "pre-event")
      pre_done <<- TRUE
      out$fired <<- TRUE
    }
    count <- function() {
      nev <<- nev + 1L
      if (nev >= max_events) out$limit <<- TRUE
    }
    run_events <- function(cand) {
      # cand: named character vector of directions (NA for scheduled)
      if (!length(cand)) return(character())
      cand <- cand[order(match(names(cand), names(plan$events)))]
      done <- character()
      for (en in names(cand)) {
        if (out$limit) break
        e <- plan$events[[en]]
        note_pre()
        st <- apply_event_resets(plan, e, t, x, d, held)
        check_states(plan, st$x, st$d, d_inf_ok, paste0("Event '", en, "'"))
        x <<- st$x; d <<- st$d
        sched <- e$trigger$kind == "scheduled"
        log_event(t, en, if (sched) "scheduled" else "zero_crossing",
                  if (sched) "model" else e$trigger$signal, cand[[en]])
        if (e$terminal) out$terminal <<- TRUE
        done <- c(done, en)
        count()
      }
      fired_ev <<- c(fired_ev, done)
      done
    }
    run_block <- function(b, j, dir, sig) {
      r <- b$def$on_zc(j, dir, t, x[b$x_idx], d[b$d_idx], block_inputs(b, sig))
      if (isTRUE(r$silent)) {
        out$silent <<- TRUE
        return(FALSE)
      }
      xn <- x; dn <- d
      if (!is.null(r$d)) dn[b$d_idx] <- r$d
      if (!is.null(r$x)) xn[b$x_idx] <- r$x
      if (identical(xn, x) && identical(dn, d) && !isTRUE(r$terminal)) return(FALSE)
      check_states(plan, xn, dn, d_inf_ok, paste0("Block '", b$name, "'"))
      note_pre()
      x <<- xn; d <<- dn
      if (isTRUE(r$terminal)) out$terminal <<- TRUE
      log_event(t, b$name, if (isTRUE(r$terminal)) "stop" else "block", b$name, dir)
      fired_blk <<- c(fired_blk, paste(b$name, j))
      count()
      TRUE
    }
    want_ok <- function(want, dir) want == "either" || dir == want || dir == "either"

    # 1. located crossings
    cand <- character()
    if (!is.null(root)) {
      sig <- evaluate(plan, t, x, d, held)$sig
      g <- zero_crossings(plan, t, x, d, sig)
      back <- NULL
      dir_of <- function(k) {
        if (abs(root$g_start[k]) > zc_tol) return(if (root$g_start[k] < 0) "rising" else "falling")
        # the function started on its level: look a short way back along
        # the trajectory, lengthening the step until the value differs (a
        # slow signal far from zero needs more than a few units of rounding)
        span <- (t - root$t_start) / 2
        h <- min(64 * tol, span)
        f <- if (plan$nx) derivatives(plan, t, x, d, sig) else numeric()
        repeat {
          if (is.null(back) || back$h != h) {
            xb <- x - h * f
            back <<- list(h = h, g = zero_crossings(plan, t - h, xb, d,
                                                    evaluate(plan, t - h, xb, d, held)$sig))
          }
          if (back$g[k] != g[k] || h >= span) break
          h <- min(h * 100, span)
        }
        if (back$g[k] < g[k]) "rising" else if (back$g[k] > g[k]) "falling" else "either"
      }
      for (i in plan$zc$blocks) {
        b <- plan$blocks[[i]]
        comps <- plan$zc$block_idx[[i]]
        for (j in seq_along(comps)) {
          k <- comps[j]
          if (!k %in% root$which) next
          dir <- dir_of(k)
          if (!want_ok(b$def$zc_direction[j], dir)) next
          if (run_block(b, j, dir, sig)) sig <- evaluate(plan, t, x, d, held)$sig
        }
      }
      for (e in zc_events) {
        if (!e$zc_index %in% root$which) next
        dir <- dir_of(e$zc_index)
        if (want_ok(e$trigger$direction, dir)) cand[e$name] <- dir
      }
    }
    # 2. model events due now
    if (at_boundary) for (en in sched_at[[bi]]) cand[en] <- NA_character_
    sig <- evaluate(plan, t, x, d, held)$sig
    g_ref <- zero_crossings(plan, t, x, d, sig)
    if (at_boundary && !is.null(g_end)) {
      for (e in zc_events) {
        if (e$name %in% names(cand)) next
        k <- e$zc_index
        dir <- jump_direction(g_end[k], g_seg0[k], g_ref[k], zc_tol)
        if (!is.na(dir) && want_ok(e$trigger$direction, dir)) cand[e$name] <- dir
      }
    }
    excl <- run_events(cand)
    # 3. sample hits and modes
    if (at_boundary && length(tick_at[[bi]])) {
      ticking <- plan$rates$block_rate %in% tick_at[[bi]] &
        vapply(plan$blocks, `[[`, TRUE, "discrete")
      ev <- evaluate(plan, t, x, d, held, ticking, refresh = TRUE)
      d <<- ev$d
      for (i in which(ticking)) {
        idx <- unlist(plan$blocks[[i]]$out_idx)
        held[idx] <<- ev$sig[idx]
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
      check_states(plan, x, dnew, d_inf_ok, "A sample update")
      d <<- dnew
    }
    ev <- evaluate(plan, t, x, d, held, NULL, refresh = TRUE)
    d <<- ev$d
    # 4. and 5. cascades and departures
    probes <- 0L
    repeat {
      layers <- 0L
      repeat {
        if (out$limit) break
        g_now <- zero_crossings(plan, t, x, d, ev$sig)
        cand <- character()
        for (e in zc_events) {
          if (e$name %in% excl) next
          dir <- cross_direction(g_ref[e$zc_index], g_now[e$zc_index], zc_tol)
          if (!is.na(dir) && want_ok(e$trigger$direction, dir)) cand[e$name] <- dir
        }
        if (!length(cand)) break
        layers <- layers + 1L
        if (layers > max_events) {
          out$limit <- TRUE
          break
        }
        excl <- run_events(cand)
        g_ref <- g_now
        ev <- evaluate(plan, t, x, d, held, NULL, refresh = TRUE)
        d <<- ev$d
      }
      if (out$limit || out$terminal || t >= tf - tol || !nzc || probes >= 10L) break
      g <- zero_crossings(plan, t, x, d, ev$sig)
      on_level <- abs(g) <= zc_tol & c(strict_zc, rep(TRUE, length(zc_events)))
      for (e in zc_events) if (e$name %in% fired_ev) on_level[e$zc_index] <- FALSE
      for (i in plan$zc$blocks) {
        comps <- plan$zc$block_idx[[i]]
        hit <- paste(plan$blocks[[i]]$name, seq_along(comps)) %in% fired_blk
        on_level[comps[hit]] <- FALSE
      }
      if (!any(on_level)) break
      probes <- probes + 1L
      gf <- departure_probe(plan, t, bounds[bi + 1L], x, d, held, ode_probe)
      dg <- gf - g
      dir_of <- function(k) if (dg[k] > 0) "rising" else if (dg[k] < 0) "falling" else NA
      moved <- FALSE
      for (i in plan$zc$blocks) {
        b <- plan$blocks[[i]]
        comps <- plan$zc$block_idx[[i]]
        for (j in seq_along(comps)) {
          k <- comps[j]
          if (!on_level[k] || is.na(dir_of(k))) next
          if (!want_ok(b$def$zc_direction[j], dir_of(k))) next
          if (run_block(b, j, dir_of(k), ev$sig)) {
            moved <- TRUE
            ev <- evaluate(plan, t, x, d, held, NULL)
          }
        }
      }
      cand <- character()
      for (e in zc_events) {
        k <- e$zc_index
        if (on_level[k] && !is.na(dir_of(k)) && want_ok(e$trigger$direction, dir_of(k))) {
          cand[e$name] <- dir_of(k)
        }
      }
      excl <- run_events(cand)
      if (!moved && !length(excl)) break
      g_ref <- g
      ev <- evaluate(plan, t, x, d, held, NULL, refresh = TRUE)
      d <<- ev$d
    }
    for (b in plan$blocks) {
      if (is.function(b$def$stop_if) &&
          isTRUE(b$def$stop_if(t, x[b$x_idx], d[b$d_idx], block_inputs(b, ev$sig)))) {
        log_event(t, b$name, "stop", b$name)
        out$terminal <- TRUE
      }
    }
    if (out$fired || out$terminal || any(abs(rec_times - t) <= tol)) {
      add_row(t, x, d, held, if (out$fired) "post-event" else "sample")
    }
    out
  }

  repeat {
    inst <- tryCatch(instant(), error = function(e) {
      failed(paste0("at t = ", format(t), ": ", conditionMessage(e)))
      NULL
    })
    if (is.null(inst)) break
    if (!is.null(root) && !inst$fired) {
      if (inst$silent) {
        # a restart at a kink; guard against restarts that do not advance
        stuck <- stuck + 1L
        if (stuck > 100L) {
          failed(paste0("restarts at kinks do not advance at t = ", format(t)))
          break
        }
      } else {
        restarts <- restarts + 1L
        if (restarts > 100L * max_events) {
          failed("too many restarts at crossings in the ignored direction")
          break
        }
      }
    }
    root <- NULL; g_end <- NULL; g_seg0 <- NULL
    if (inst$terminal) {
      status <- "terminated"
      message <- paste0("stopped by a terminal event at t = ", format(t))
      break
    }
    if (inst$limit) {
      status <- "event_limit"
      message <- paste0("reached max_events = ", max_events, " at t = ", format(t),
                        "; the model may be chattering or Zeno")
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
                             rtol, atol, maxsteps, zc_tol, nzc, ode_args, tol)
    stats_acc <- add_stats(stats_acc, seg$stats)
    warnings <- unique(c(warnings, seg$warnings))
    for (k in seq_along(seg$times)) add_row(seg$times[k], seg$x[k, ], d, held, "sample")
    if (!is.null(seg$error)) {
      failed(seg$error)
      break
    }
    g_end <- seg$g_end
    g_seg0 <- seg$g_start
    if (is.null(seg$root)) {
      x <- seg$x_end
      stuck <- 0L
      t <- tn; bi <- bi + 1L; at_boundary <- TRUE
      next
    }
    tr <- seg$root$time
    if (abs(tr - t) > tol) stuck <- 0L
    x <- seg$root$x
    root <- list(which = seg$root$which, g_start = seg$g_start, t_start = t)
    if (abs(tr - tn) <= tol) {
      t <- tn; bi <- bi + 1L; at_boundary <- TRUE
    } else {
      t <- tr; at_boundary <- FALSE; g_end <- NULL
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

odepack_methods <- c("lsoda", "lsode", "lsodes", "lsodar", "vode", "bdf", "bdf_d",
                     "adams", "impAdams", "impAdams_d")

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
                              rtol, atol, maxsteps, zc_tol, nzc, ode_args, tol) {
  y0 <- if (dummy) 0 else x
  xs <- function(y) if (dummy) numeric() else y
  # Inside a segment, sources with breakpoints take their left limit at the
  # end point: a source that jumps at tn must not leak its next value into
  # this segment when the solver evaluates the right-hand side exactly at
  # tn. All other blocks see the true solver time.
  t_left <- tn - 4 * .Machine$double.eps * max(1, abs(tn))
  # The solver runs on time since the segment start, tau in [0, L], so the
  # length it integrates is exact even when t is far from zero (at t = 1e9
  # one unit in the last place is 1.2e-7). Blocks see the absolute time.
  L <- tn - t
  at <- function(tau) if (tau >= L) tn else t + tau
  signals_at <- function(tau, y) evaluate(plan, at(tau), xs(y), d, held, t_left = t_left)$sig
  func <- function(tau, y, parms) {
    sig <- signals_at(tau, y)
    list(if (dummy) 0 else derivatives(plan, at(tau), y, d, sig))
  }
  inner_local <- inner - t
  args <- c(list(y = y0, times = c(0, inner_local, L), func = func, parms = NULL,
                 method = method, rtol = rtol, atol = atol, maxsteps = maxsteps,
                 tcrit = L),
            ode_args)
  g0 <- NULL
  if (nzc) {
    g0 <- zero_crossings(plan, t, x, d, evaluate(plan, t, x, d, held)$sig)
    at0 <- abs(g0) <= zc_tol
    off <- ifelse(at0, g0, 0)
    args$rootfunc <- function(tau, y, parms) {
      xx <- xs(y)
      v <- zero_crossings(plan, at(tau), xx, d, signals_at(tau, y)) - off
      # A function that starts at zero and stays exactly there (a constant
      # signal sitting on its threshold) would make the root finder stop with
      # "root too near the initial point". Past the start, report an exact
      # zero as positive; at the start itself keep the zero, so the root
      # finder takes its reference sign from the first step.
      if (tau > 0) v[at0 & v == 0] <- 1
      v
    }
  }
  empty <- list(times = numeric(), x = matrix(numeric(), 0L, length(x)),
                stats = c(segments = 1, steps = NA, rhs_evaluations = NA, jacobians = NA, roots = 0))
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
  # ODEPACK solvers report Jacobian evaluations and LU decompositions in
  # istate[14]; the Runge-Kutta solvers have no Jacobian (NA)
  jac <- if (method %in% odepack_methods && length(ist) >= 14L) ist[14] else NA_real_
  stats <- c(segments = 1, steps = ist[2] %||% NA, rhs_evaluations = ist[3] %||% NA,
             jacobians = jac, roots = 0)
  tau <- out[, 1]
  ymat <- out[, -1, drop = FALSE]
  if (dummy) ymat <- matrix(numeric(), nrow(out), 0L)
  troot <- attr(out, "troot")
  is_root <- length(troot) > 0L
  last <- nrow(out)
  # back to absolute time; requested output times keep their exact values
  tt <- vapply(tau, at, 0)
  hit <- match(tau, inner_local)
  tt[!is.na(hit)] <- inner[hit[!is.na(hit)]]
  keep <- which(!is.na(hit) & (if (is_root) tau < troot[1] else TRUE))
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
  if (!is_root && (flag < 0 || abs(tt[last] - tn) > tol)) {
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
    # a root on the segment end is also a boundary: keep the left limit, so
    # crossings caused by jumps there are still seen
    if (abs(tt[last] - tn) <= tol) {
      xe <- xs(ymat[last, ])
      res$g_end <- zero_crossings(plan, tn, xe, d,
                                  evaluate(plan, tn, xe, d, held, t_left = t_left)$sig)
    }
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

# Direction of a crossing between the value a just before an instant and
# the value b after it, or NA. A value within zc_tol of zero is on the
# level: reaching the level counts as crossing it, leaving it does not (a
# departure from the level is found by departure_probe()).
cross_direction <- function(a, b, zc_tol) {
  sa <- if (abs(a) > zc_tol) sign(a) else 0
  sb <- if (abs(b) > zc_tol) sign(b) else 0
  if (sa == 0 || sb == sa) return(NA_character_)
  if (sb == 0) return(if (sa < 0) "rising" else "falling")
  if (sb > 0) "rising" else "falling"
}

# Crossing across a segment boundary: g_end is the left limit, g_start the
# value at the start of the segment, g the value at the boundary. A
# function that arrived exactly at its level at the boundary crosses there
# in the direction it came from.
jump_direction <- function(g_end, g_start, g, zc_tol) {
  a <- if (abs(g_end) > zc_tol) sign(g_end) else
    if (abs(g_start) > zc_tol && abs(g) <= zc_tol) sign(g_start) else 0
  b <- if (abs(g) > zc_tol) sign(g) else if (a != 0 && abs(g_end) <= zc_tol) -a else 0
  if (a == 0 || b == 0 || a == b) return(NA_character_)
  if (b > 0) "rising" else "falling"
}

# Crossing functions a short way into the coming segment, for functions that
# sit on their level now. The states are integrated (not extrapolated), so
# a function whose first derivative vanishes, such as x = t^2, is still seen
# to move. The probe covers 1/1000 of the segment: a function that leaves
# its level and returns within that span is not detected.
departure_probe <- function(plan, t, tn, x, d, held, ode_args) {
  h <- (tn - t) / 1000
  th <- t + h
  t_left <- tn - 4 * .Machine$double.eps * max(1, abs(tn))
  xf <- x
  if (plan$nx) {
    # local time, as in integrate_segment()
    func <- function(tau, y, parms) {
      tt <- t + tau
      list(derivatives(plan, tt, y, d, evaluate(plan, tt, y, d, held, t_left = t_left)$sig))
    }
    out <- suppressWarnings(do.call(deSolve::ode, c(list(y = x, times = c(0, h), func = func,
                                                          parms = NULL), ode_args)))
    xf <- out[nrow(out), -1L]
    if (nrow(out) < 2L || !all(is.finite(xf))) sim_abort("the departure probe failed.")
  }
  zero_crossings(plan, th, xf, d, evaluate(plan, th, xf, d, held, t_left = t_left)$sig)
}

# Discrete states that may hold -Inf: those a block declares with
# `nonfinite_ok` (the state machine's "not entered yet" time).
nonfinite_allowed <- function(plan) {
  ok <- logical(plan$nd)
  for (b in plan$blocks) {
    k <- b$def$nonfinite_ok
    if (length(k)) ok[b$d_idx[k]] <- TRUE
  }
  ok
}

check_states <- function(plan, x, d, d_inf_ok, who) {
  if (length(x) && !all(is.finite(x))) {
    bad <- which(!is.finite(x))[1]
    sim_abort(who, " made continuous state '", plan$state_names[bad], "' ",
              format(x[bad]), "; states must be finite.")
  }
  badd <- !is.finite(d) & !(d_inf_ok & d == -Inf)
  if (any(badd)) {
    bad <- which(badd)[1]
    sim_abort(who, " made discrete state '", plan$dstate_names[bad], "' ",
              format(d[bad]), "; states must be finite.")
  }
}

add_stats <- function(acc, s) {
  for (nm in names(s)) {
    acc[[nm]] <- if (is.na(acc[[nm]]) && is.na(s[[nm]])) NA_real_ else
      sum(acc[[nm]], s[[nm]], na.rm = TRUE)
  }
  acc
}

apply_event_resets <- function(plan, e, t, x, d, held, sig = NULL) {
  if (!length(e$resets)) return(list(x = x, d = d))
  if (is.null(sig)) sig <- evaluate(plan, t, x, d, held)$sig
  s <- signal_list(plan, sig)
  for (r in e$resets) {
    cur <- if (r$kind == "x") x[r$idx] else d[r$idx]
    v <- as.numeric(r$fun(cur, s, t))
    if (length(v) != length(r$idx)) {
      sim_abort("Event '", e$name, "': reset of '", r$target, "' returned ",
                length(v), " value(s), expected ", length(r$idx), ".")
    }
    if (r$kind == "x") x[r$idx] <- v else d[r$idx] <- v
  }
  list(x = x, d = d)
}
