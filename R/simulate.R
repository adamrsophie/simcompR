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
#'   `"lsodar"`. The fixed-step methods `"euler"` and `"rk4"` take steps of
#'   at most `hmax`, by default `(tf - t0) / 1000`, shortened to end at
#'   every boundary and output time.
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
#'   only in its trigger direction (step 5 above). Leaving means moving
#'   away from zero; a function that approaches zero inside the band, or
#'   stays on it, is watched until it leaves the band, and the event fires
#'   then in the direction of the side it leaves on. A function that enters
#'   the band without reaching zero (an asymptotic approach) is reported
#'   once, at the next segment boundary.
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
  # never so coarse that the two ends of a short span would merge
  tol <- min(time_tol(c(t0, tf)), (tf - t0) / 4)
  if (is.null(times)) times <- seq(t0, tf, length.out = 201L)
  times <- check_numeric(times, "`times`")
  if (any(times < t0 - tol | times > tf + tol)) sim_abort("`times` must lie within `tspan`.")
  times <- sort(unique(pmin(pmax(times, t0), tf)))
  check_number(rtol, "`rtol`"); check_number(atol, "`atol`")
  if (!is_count(max_events)) sim_abort("`max_events` must be a positive integer.")
  if (!is_count(maxsteps)) sim_abort("`maxsteps` must be a positive integer.")
  if (!is_scalar_number(zc_tol) || zc_tol < 0) {
    sim_abort("`zc_tol` must be a single finite number of at least 0.")
  }
  if (!is.logical(include_ticks) || length(include_ticks) != 1L || is.na(include_ticks)) {
    sim_abort("`include_ticks` must be TRUE or FALSE.")
  }
  if (!is.null(hmax) && (!is_scalar_number(hmax) || hmax <= 0)) {
    sim_abort("`hmax` must be a single positive number.")
  }
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
  # every time that cuts the interval, tagged with what it is
  bp_list <- lapply(seq_along(plan$blocks), function(i) {
    b <- plan$blocks[[i]]
    v <- if (is.function(b$def$breakpoints)) b$def$breakpoints(t0, tf)
    v[v > t0 & v < tf]
  })
  sched_times <- lapply(sched, function(e) {
    v <- e$trigger$times
    v[v >= t0 - tol & v <= tf + tol]
  })
  all_t <- c(t0, tf, unlist(ticks), unlist(sched_times), unlist(bp_list))
  kind <- c(0L, 0L, rep(1L, length(unlist(ticks))), rep(2L, length(unlist(sched_times))),
            rep(3L, length(unlist(bp_list))))
  who <- c(0L, 0L, rep(seq_along(ticks), lengths(ticks)),
           rep(seq_along(sched_times), lengths(sched_times)),
           rep(seq_along(bp_list), lengths(bp_list)))
  mg <- merge_groups(all_t, tol)
  bounds <- mg$rep
  # the groups of t0 and tf are the end points themselves
  bounds[mg$group[1]] <- t0
  bounds[mg$group[2]] <- tf
  if (mg$group[1] == mg$group[2]) {
    sim_abort("`tspan` is too short to be resolved in double precision.")
  }
  tick_at <- sched_at <- bp_at <- vector("list", length(bounds))
  for (k in which(kind == 1L)) tick_at[[mg$group[k]]] <- c(tick_at[[mg$group[k]]], who[k])
  for (k in which(kind == 2L)) {
    sched_at[[mg$group[k]]] <- c(sched_at[[mg$group[k]]], names(sched)[who[k]])
  }
  for (k in which(kind == 3L)) bp_at[[mg$group[k]]] <- c(bp_at[[mg$group[k]]], who[k])
  tick_at <- lapply(tick_at, unique)
  sched_at <- lapply(sched_at, unique)
  rec_times <- times
  if (include_ticks && length(ticks)) {
    # a hit is recorded at its boundary, the instant where it takes effect
    tk <- bounds[vapply(tick_at, length, 0L) > 0L]
    extra <- unique(tk[is.na(match_tol(tk, rec_times, tol))])
    rec_times <- sort(c(rec_times, extra))
  }

  ode_args <- list(...)
  # With crossing functions, cap the step so that two crossings inside one
  # step (which the root finder cannot see) need a fast signal.
  if (is.null(hmax) && nzc && n_strict) hmax <- (tf - t0) / 100
  if (!is.null(hmax)) ode_args$hmax <- hmax
  # fixed-step methods step from one requested time to the next; give them
  # a grid of their own so results do not depend on `times`
  fixed_step <- if (method %in% c("euler", "rk4")) hmax %||% ((tf - t0) / 1000)
  # without continuous states nothing is integrated, so kinks need no stop
  if (plan$nx == 0L && !n_strict) nzc <- 0L
  dummy <- plan$nx == 0L && nzc > 0L

  rows <- list()
  add_row <- function(t, x, d, held, kind, loop = NULL) {
    rows[[length(rows) + 1L]] <<- list(t = t, x = x, d = d, held = held, kind = kind,
                                       loop = loop)
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
  resting <- logical(nzc)  # functions on their level that the probe saw stay
  depart <- numeric(nzc)   # and the direction of those it saw leave
  root <- NULL     # crossings located by the solver at t
  g_end <- NULL    # crossing functions at the left limit of t
  g_seg0 <- NULL   # and at the start of the segment that ended at t
  zc_events <- Filter(function(e) e$trigger$kind == "zero_crossing", plan$events)
  strict_zc <- !optional[seq_len(plan$zc$n)]
  # kinks of outputs that jump are probed when the input sits on them
  for (i in plan$zc$blocks) {
    if (isTRUE(plan$blocks[[i]]$def$kink_side)) strict_zc[plan$zc$block_idx[[i]]] <- TRUE
  }
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
  #      signal at t (including jumps made by step 1); each sees the state
  #      left by the one before;
  #   3. sample hits: modes brought up to date with the state after step 2,
  #      outputs held, updates, then the modes again;
  #   4. events whose function was carried across its level by any single
  #      handler, sample hit or mode change above, in further layers, until
  #      no new crossing occurs;
  #   5. block handlers and events whose function sits on its level and
  #      moves off it in the trigger direction (found with short probe
  #      integrations), then step 4 again.
  # A terminal event ends the run after the instant is settled; every
  # handler execution counts toward max_events.
  instant <- function() {
    fired_ev <- character()
    fired_blk <- character()
    out <- list(terminal = FALSE, limit = FALSE, fired = FALSE, silent = FALSE)
    pre_done <- FALSE
    pending <- character()
    tick_loop <- NULL
    note_pre <- function() {
      if (!pre_done) add_row(t, x, d, held, "pre-event")
      pre_done <<- TRUE
      out$fired <<- TRUE
    }
    g_here <- function() zero_crossings(plan, t, x, d, evaluate(plan, t, x, d, held)$sig)
    # events carried across their level between crossing values a and b
    add_crossings <- function(a, b, excl = character()) {
      for (e in zc_events) {
        if (e$name %in% excl) next
        dir <- cross_direction(a[e$zc_index], b[e$zc_index], zc_tol)
        if (!is.na(dir) && want_ok(e$trigger$direction, dir)) pending[e$name] <<- dir
      }
    }
    budget_left <- function() {
      if (nev >= max_events) {
        out$limit <<- TRUE
        return(FALSE)
      }
      TRUE
    }
    run_events <- function(cand) {
      # cand: named character vector of directions (NA for scheduled)
      if (!length(cand)) return(character())
      cand <- cand[order(match(names(cand), names(plan$events)))]
      done <- character()
      g_cur <- g_here()
      for (en in names(cand)) {
        if (!budget_left()) break
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
        nev <<- nev + 1L
        g_new <- g_here()
        # an event does not retrigger itself from the layer it fired in
        add_crossings(g_cur, g_new, excl = names(cand))
        g_cur <- g_new
      }
      fired_ev <<- c(fired_ev, done)
      done
    }
    run_block <- function(b, j, dir, sig) {
      if (out$limit) return(FALSE)
      r <- b$def$on_zc(j, dir, t, x[b$x_idx], d[b$d_idx], block_inputs(b, sig))
      if (isTRUE(r$silent)) {
        out$silent <<- TRUE
        # a kink may move the side a jumping output is on: not an event,
        # but what it carries across a level is
        if (!is.null(r$d) && !identical(as.numeric(r$d), d[b$d_idx])) {
          g_a <- g_here()
          d[b$d_idx] <<- r$d
          add_crossings(g_a, g_here())
          return(TRUE)
        }
        return(FALSE)
      }
      xn <- x; dn <- d
      for (f in c("d", "x")) {
        n_want <- length(b[[paste0(f, "_idx")]])
        if (!is.null(r[[f]]) && (!is.numeric(r[[f]]) || length(r[[f]]) != n_want)) {
          sim_abort("Block '", b$name, "' on_zc returned `", f, "` of length ",
                    length(r[[f]]), ", expected ", n_want, " number(s).")
        }
      }
      if (!is.null(r$d)) dn[b$d_idx] <- r$d
      if (!is.null(r$x)) xn[b$x_idx] <- r$x
      if (identical(xn, x) && identical(dn, d) && !isTRUE(r$terminal)) return(FALSE)
      if (!budget_left()) return(FALSE)
      check_states(plan, xn, dn, d_inf_ok, paste0("Block '", b$name, "'"))
      note_pre()
      g_a <- g_here()
      x <<- xn; d <<- dn
      add_crossings(g_a, g_here())
      if (isTRUE(r$terminal)) out$terminal <<- TRUE
      log_event(t, b$name, if (isTRUE(r$terminal)) "stop" else "block", b$name, dir)
      fired_blk <<- c(fired_blk, paste(b$name, j))
      nev <<- nev + 1L
      TRUE
    }
    want_ok <- function(want, dir) want == "either" || dir == want || dir == "either"
    refresh_modes <- function() {
      g_a <- g_here()
      ev <- evaluate(plan, t, x, d, held, NULL, refresh = TRUE)
      d <<- ev$d
      add_crossings(g_a, g_here())
      ev
    }

    # 1. located crossings
    cand <- character()
    if (!is.null(root)) {
      sig <- evaluate(plan, t, x, d, held)$sig
      g <- zero_crossings(plan, t, x, d, sig)
      back <- NULL
      dir_of <- function(k) {
        # a resting function left its band on the side it is now on
        if (isTRUE(root$resting[k])) return(if (g[k] > 0) "rising" else "falling")
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
      for (e in zc_events) {
        if (!e$zc_index %in% root$which) next
        dir <- dir_of(e$zc_index)
        if (want_ok(e$trigger$direction, dir)) cand[e$name] <- dir
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
    }
    # 2. model events due now
    if (at_boundary) for (en in sched_at[[bi]]) cand[en] <- NA_character_
    if (at_boundary && !is.null(g_end)) {
      g_ref <- g_here()
      for (e in zc_events) {
        k <- e$zc_index
        dir <- jump_direction(g_end[k], g_seg0[k], g_ref[k], zc_tol)
        if (!is.na(dir) && want_ok(e$trigger$direction, dir)) pending[e$name] <- dir
      }
    }
    for (en in setdiff(names(pending), names(cand))) cand[en] <- pending[[en]]
    pending <- character()
    run_events(cand)
    # 3. sample hits and modes
    if (at_boundary && length(tick_at[[bi]])) {
      g_a <- g_here()
      ticking <- plan$rates$block_rate %in% tick_at[[bi]] &
        vapply(plan$blocks, `[[`, TRUE, "discrete")
      ev <- evaluate(plan, t, x, d, held, ticking, refresh = TRUE)
      # the loop statistics of the solve at the hit, kept for the record
      tick_loop <- ev$loop
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
      add_crossings(g_a, g_here())
    }
    ev <- refresh_modes()
    # 4. and 5. cascades and departures
    resting <<- logical(nzc)
    depart <<- numeric(nzc)
    probes <- 0L
    repeat {
      layers <- 0L
      while (length(pending) && !out$limit) {
        layers <- layers + 1L
        if (layers > max_events) {
          out$limit <- TRUE
          break
        }
        cand <- pending
        pending <- character()
        run_events(cand)
        ev <- refresh_modes()
      }
      if (out$limit || t >= tf - tol || !nzc) break
      g <- zero_crossings(plan, t, x, d, ev$sig)
      level_any <- abs(g) <= zc_tol & c(strict_zc, rep(TRUE, length(zc_events)))
      on_level <- level_any
      for (e in zc_events) if (e$name %in% fired_ev) on_level[e$zc_index] <- FALSE
      for (i in plan$zc$blocks) {
        comps <- plan$zc$block_idx[[i]]
        hit <- paste(plan$blocks[[i]]$name, seq_along(comps)) %in% fired_blk
        on_level[comps[hit]] <- FALSE
      }
      resting <<- logical(nzc)
      depart <<- numeric(nzc)
      if (!any(level_any)) break
      probes <- probes + 1L
      if (probes > max_events) {
        out$limit <- TRUE
        break
      }
      dg <- departure_probe(plan, t, bounds[bi + 1L], x, d, held, ode_probe, hmax,
                            bp_at[[bi + 1L]], g)
      # a departure moves away from the level; a function that stays on it
      # or approaches it inside the band rests and is watched for leaving
      away <- abs(g + dg) > abs(g)
      resting <<- level_any & !away
      depart <<- structure(ifelse(level_any & away, sign(dg), 0), tau = attr(dg, "h"),
                           amp = abs(g) + abs(dg))
      dir_of <- function(k) {
        if (!away[k]) NA else if (dg[k] > 0) "rising" else if (dg[k] < 0) "falling" else NA
      }
      moved <- FALSE
      # blocks whose input rests on a level take their mode for the plateau
      for (i in plan$zc$blocks) {
        b <- plan$blocks[[i]]
        if (!is.function(b$def$on_rest)) next
        comps <- plan$zc$block_idx[[i]]
        for (j in seq_along(comps)) {
          if (!resting[comps[j]]) next
          r <- b$def$on_rest(j, t, x[b$x_idx], d[b$d_idx], block_inputs(b, ev$sig))
          if (!is.null(r$d) && !identical(as.numeric(r$d), d[b$d_idx])) {
            g_a <- g_here()
            d[b$d_idx] <<- r$d
            add_crossings(g_a, g_here())
            ev <- evaluate(plan, t, x, d, held, NULL)
            moved <- TRUE
          }
        }
      }
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
      done <- run_events(cand)
      if (!moved && !length(done)) break
      resting <<- logical(nzc)
      depart <<- numeric(nzc)
      ev <- refresh_modes()
    }
    for (b in plan$blocks) {
      if (!is.function(b$def$stop_if)) next
      if (any(startsWith(fired_blk, paste0(b$name, " ")))) next
      if (isTRUE(b$def$stop_if(t, x[b$x_idx], d[b$d_idx], block_inputs(b, ev$sig)))) {
        note_pre()
        log_event(t, b$name, "stop", b$name)
        out$terminal <- TRUE
      }
    }
    if (out$fired || out$terminal || !is.na(match_tol(t, rec_times, tol))) {
      add_row(t, x, d, held, if (out$fired) "post-event" else "sample", tick_loop)
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
          failed(paste0("restarts at kinks do not advance at t = ", format(t),
                        "; the model may be chattering at a kink (a sliding mode)"))
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
    if (inst$limit) {
      status <- "event_limit"
      message <- paste0("reached max_events = ", max_events, " at t = ", format(t),
                        "; the model may be chattering or Zeno")
      break
    }
    if (inst$terminal) {
      status <- "terminated"
      message <- paste0("stopped by a terminal event at t = ", format(t))
      break
    }
    if (t >= tf - tol) break

    tn <- bounds[bi + 1L]
    # requested times strictly inside the segment, found by bisection
    lo <- findInterval(t + tol, rec_times) + 1L
    hi <- findInterval(tn - tol, rec_times, left.open = TRUE)
    inner <- if (hi >= lo) rec_times[lo:hi] else numeric()
    inner <- inner[inner > t + tol & inner < tn - tol]
    if (plan$nx == 0L && !dummy) {
      for (tt in inner) add_row(tt, x, d, held, "sample")
      t <- tn; bi <- bi + 1L; at_boundary <- TRUE
      next
    }

    seg <- integrate_segment(plan, t, tn, inner, x, d, held, dummy, method,
                             rtol, atol, maxsteps, zc_tol, nzc, ode_args, tol,
                             left = bp_at[[bi + 1L]], resting = resting, depart = depart,
                             fixed_step = fixed_step)
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
    # restarts at kinks that advance by less than about 1e-8 of the time
    # scale are chattering (a sliding mode); counted below
    if (tr - t > sqrt(.Machine$double.eps) * max(1, abs(t))) stuck <- 0L
    x <- seg$root$x
    root <- list(which = seg$root$which, g_start = seg$g_start, t_start = t,
                 resting = resting)
    if (abs(tr - tn) <= tol) {
      # the root and the boundary are one instant; carry the states over
      # the (at most tol) gap so the boundary row is consistent with its time
      if (tn != tr && plan$nx) {
        sig_r <- evaluate(plan, tr, x, d, held, t_left = left_of(t, tn),
                          left = bp_at[[bi + 1L]])$sig
        x <- x + (tn - tr) * derivatives(plan, tr, x, d, sig_r)
      }
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
  if (anyDuplicated(names(initial))) {
    sim_abort("`initial` names block '", names(initial)[anyDuplicated(names(initial))],
              "' more than once.")
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
                              rtol, atol, maxsteps, zc_tol, nzc, ode_args, tol,
                              left = integer(), resting = logical(), depart = numeric(),
                              fixed_step = NULL) {
  y0 <- if (dummy) 0 else x
  xs <- function(y) if (dummy) numeric() else y
  # Inside a segment, sources with breakpoints take their left limit at the
  # end point: a source that jumps at tn must not leak its next value into
  # this segment when the solver evaluates the right-hand side exactly at
  # tn. All other blocks see the true solver time.
  t_left <- left_of(t, tn)
  # The solver runs on time since the segment start, tau in [0, L], so the
  # length it integrates is exact even when t is far from zero (at t = 1e9
  # one unit in the last place is 1.2e-7). Blocks see the absolute time.
  L <- tn - t
  at <- function(tau) if (tau >= L) tn else t + tau
  signals_at <- function(tau, y) evaluate(plan, at(tau), xs(y), d, held, t_left = t_left, left = left)$sig
  tau_fail <- NA_real_
  func <- function(tau, y, parms) {
    tryCatch({
      sig <- signals_at(tau, y)
      list(if (dummy) 0 else derivatives(plan, at(tau), y, d, sig))
    }, error = function(e) {
      # remember the earliest time a block failed at
      if (is.na(tau_fail) || tau < tau_fail) tau_fail <<- tau
      stop(e)
    })
  }
  inner_local <- inner - t
  grid <- numeric()
  if (!is.null(fixed_step)) {
    n <- ceiling(L / fixed_step * (1 - 1e-12))
    if (n > 1L) grid <- L * seq_len(n - 1L) / n
  }
  args <- c(list(y = y0, times = sort(unique(c(0, inner_local, grid, L))), func = func, parms = NULL,
                 method = method, rtol = rtol, atol = atol, maxsteps = maxsteps,
                 tcrit = L),
            ode_args)
  g0 <- NULL
  if (nzc) {
    g0 <- zero_crossings(plan, t, x, d, evaluate(plan, t, x, d, held)$sig)
    if (length(resting) != length(g0)) resting <- logical(length(g0))
    if (length(depart) != length(g0)) depart <- numeric(length(g0))
    dep <- depart != 0 & !resting
    at0 <- abs(g0) <= zc_tol & !resting
    off <- ifelse(at0, g0, 0)
    args$rootfunc <- function(tau, y, parms) {
      xx <- xs(y)
      g <- zero_crossings(plan, at(tau), xx, d, signals_at(tau, y))
      v <- g - off
      # A function that starts at zero and stays exactly there (a constant
      # signal sitting on its threshold) would make the root finder stop with
      # "root too near the initial point". Past the start, report an exact
      # zero as positive; at the start itself keep the zero, so the root
      # finder takes its reference sign from the first step.
      if (tau > 0) v[at0 & !dep & v == 0] <- 1
      # a function seen leaving its level keeps the sign it leaves with
      # until it comes back across the band
      # (an offset that fades out over the time the probe needed to see the
      # move, so the root finder starts from a definite sign and the level
      # itself is located exactly later on)
      if (any(dep)) {
        w <- pmax(0, 1 - tau / attr(depart, "tau")[dep])
        v[dep] <- g[dep] + depart[dep] * attr(depart, "amp")[dep] * w
      }
      # a function resting on its level is watched for leaving it in
      # either direction: |g| grows past zc_tol
      if (any(resting)) v[resting] <- abs(v[resting]) - zc_tol
      # a function seen leaving its level keeps the sign it leaves with
      # until it comes back across the band
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
    msg <- sub("[.]+$", "", conditionMessage(out))
    ok <- inner[inner - t < tau_fail]
    if (!is.na(tau_fail) && length(ok)) {
      # keep the output times reached before the failure: integrate again
      # up to the last of them (which may also find an earlier event)
      te <- max(ok)
      res <- integrate_segment(plan, t, te, ok[ok < te], x, d, held, dummy, method,
                               rtol, atol, maxsteps, zc_tol, nzc, ode_args, tol,
                               left = integer(), resting = resting, depart = depart,
                               fixed_step = fixed_step)
      if (!is.null(res$root) || !is.null(res$error)) return(res)
      res$times <- c(res$times, te)
      res$x <- rbind(res$x, res$x_end)
      res$error <- paste0("at t = ", format(at(tau_fail)), ": ", msg)
      res$warnings <- unique(c(res$warnings, warns))
      return(res)
    }
    empty$error <- paste0(if (is.na(tau_fail)) paste0("in [", format(t), ", ", format(tn), "]")
                          else paste0("at t = ", format(at(tau_fail))), ": ", msg)
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
  # output times at the root instant itself belong to the instant's rows
  keep <- which(!is.na(hit) & (if (is_root) tt < tt[last] - tol else TRUE))
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
                                  evaluate(plan, tn, xe, d, held, t_left = t_left, left = left)$sig)
    }
  } else {
    res$x_end <- ymat[last, ]
    if (nzc) {
      xe <- xs(res$x_end)
      res$g_end <- zero_crossings(plan, tn, xe, d,
                                  evaluate(plan, tn, xe, d, held, t_left = t_left, left = left)$sig)
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
  if (sb == sa) return(NA_character_)
  # leaving the level by a jump departs in the direction of the jump
  if (sa == 0) return(if (sb > 0) "rising" else "falling")
  if (sb == 0) return(if (sa < 0) "rising" else "falling")
  if (sb > 0) "rising" else "falling"
}

# Crossing across a segment boundary: g_end is the left limit, g_start the
# value at the start of the segment, g the value at the boundary. Reaching
# the level, leaving it and jumping across it all count; a function that
# rested on its level through the segment and jumps off it departs in the
# direction of the jump. A function that arrived at its level at the end of
# the segment crosses in the direction it came from.
jump_direction <- function(g_end, g_start, g, zc_tol) {
  side <- function(v) if (abs(v) > zc_tol) sign(v) else 0
  a <- side(g_end); b <- side(g)
  dir <- function(up) if (up) "rising" else "falling"
  if (a != 0) {
    if (b == a) return(NA_character_)
    return(dir(if (b == 0) a < 0 else b > 0))
  }
  s0 <- side(g_start)
  if (s0 != 0) {
    # arrived at the level during the segment
    if (b == s0) return(NA_character_)
    return(dir(s0 < 0))
  }
  if (b == 0) NA_character_ else dir(b > 0)
}

# Change of the crossing functions just after t, for functions that sit on
# their level now. The states are integrated (not extrapolated), so a
# function whose first derivative vanishes, such as x = t^3, is still seen
# to move. The probe looks at a ladder of short times, from 1e-8 of the
# largest step up to 1/1000 of the coming segment (or hmax), plus the
# shortest step that changes t at all; for each function it keeps the
# change at the shortest time where the function has moved, so a later
# reversal does not hide the direction it leaves in. Returns 0 where the
# function did not move.
departure_probe <- function(plan, t, tn, x, d, held, ode_args, hmax, left, g) {
  H <- (tn - t) / 1000
  if (!is.null(hmax)) H <- min(H, hmax)
  s_min <- 2 * .Machine$double.eps * abs(t)
  hs <- H * 10^(-8:0)
  hs <- hs[hs > s_min]
  if (s_min > 0 && s_min <= (tn - t) / 2) hs <- c(hs, s_min)
  if (!length(hs)) hs <- (tn - t) / 2
  hs <- sort(unique(hs))
  t_left <- left_of(t, tn)
  xs <- matrix(x, length(hs), length(x), byrow = TRUE)
  if (plan$nx) {
    # local time, as in integrate_segment()
    func <- function(tau, y, parms) {
      tt <- t + tau
      list(derivatives(plan, tt, y, d, evaluate(plan, tt, y, d, held, t_left = t_left,
                                                left = left)$sig))
    }
    out <- suppressWarnings(do.call(deSolve::ode, c(list(y = x, times = c(0, hs), func = func,
                                                          parms = NULL), ode_args)))
    if (nrow(out) < length(hs) + 1L || !all(is.finite(out[, -1L]))) {
      sim_abort("the integration just after this instant failed; the model may be ",
                "chattering at a discontinuity (a sliding mode).")
    }
    xs <- out[-1L, -1L, drop = FALSE]
  }
  dg <- numeric(length(g))
  h_seen <- numeric(length(g))
  for (i in seq_along(hs)) {
    th <- t + hs[i]
    gi <- zero_crossings(plan, th, xs[i, ], d,
                         evaluate(plan, th, xs[i, ], d, held, t_left = t_left, left = left)$sig)
    todo <- dg == 0 & gi != g
    dg[todo] <- (gi - g)[todo]
    h_seen[todo] <- hs[i]
    if (all(dg != 0)) break
  }
  structure(dg, h = h_seen)
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
  if (nzchar(e$prefix %||% "")) {
    # inside a subsystem instance its own signals go by their local names
    loc <- names(s)[startsWith(names(s), e$prefix)]
    s_loc <- s[loc]
    names(s_loc) <- substring(loc, nchar(e$prefix) + 1L)
    s <- c(s_loc, s[setdiff(names(s), names(s_loc))])
  }
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
