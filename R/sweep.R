#' Parameter sweeps and Monte Carlo experiments
#'
#' `sim_sweep()` runs one simulation per row of `parameters`, a data frame
#' (or list of named lists) of model parameter values. The model is
#' compiled once; each run uses [sim_set_params()], which reuses the
#' dependency analysis whenever the structure is unchanged.
#'
#' `sim_monte_carlo()` draws `n` parameter sets with `sampler(i)` under a
#' fixed `seed` and runs them like a sweep. The caller's random number
#' stream is saved and restored, so the experiment is reproducible and does
#' not disturb the global random state.
#'
#' With `cores > 1` runs use [parallel::mclapply()] (not available on
#' Windows, where runs are sequential). Results do not depend on `cores`.
#'
#' @param x A [sim_model()] or a compiled `sim_plan`.
#' @param parameters Data frame of parameter values, one run per row, or a
#'   list of named lists.
#' @param n Number of Monte Carlo runs.
#' @param sampler Function of the run index returning a named list of
#'   parameter values.
#' @param seed Integer seed for `sampler`.
#' @param tspan,times Passed to [sim_run()].
#' @param summary Optional function of a `sim_result` returning a named
#'   numeric vector. Its values are collected in `$summary`.
#' @param keep Keep the full results (`TRUE`) or only the summaries.
#' @param cores Number of parallel workers.
#' @param ... Further arguments for [sim_run()].
#' @return A `sim_sweep` object with elements `parameters`, `results` (when
#'   `keep = TRUE`), `status` and `summary` (a data frame of parameters and
#'   summary values).
#' @examples
#' m <- sim_model("decay", params = list(a = 1))
#' sim_add_block(m, "integrator", "x", initial = 1)
#' sim_add_block(m, "gain", "k", gain = sim_param("a"))
#' sim_add_block(m, "gain", "neg", gain = -1)
#' sim_connect(m, "x", "k")
#' sim_connect(m, "k", "neg")
#' sim_connect(m, "neg", "x")
#' s <- sim_sweep(m, data.frame(a = c(0.5, 1, 2)), tspan = c(0, 1),
#'   summary = function(r) c(x_end = tail(sim_signal(r, "x"), 1)))
#' s$summary
#' mc <- sim_monte_carlo(m, n = 5, sampler = function(i) list(a = runif(1, 0, 2)),
#'   seed = 42, tspan = c(0, 1), keep = FALSE,
#'   summary = function(r) c(x_end = tail(sim_signal(r, "x"), 1)))
#' mc$summary
#' @export
sim_sweep <- function(x, parameters, tspan, times = NULL, summary = NULL,
                      keep = TRUE, cores = 1, ...) {
  plan <- if (inherits(x, "sim_model")) sim_compile(x) else check_plan(x)
  sets <- parameter_sets(parameters)
  run_sets(plan, sets, tspan, times, summary, keep, cores, ...)
}

#' @rdname sim_sweep
#' @export
sim_monte_carlo <- function(x, n, sampler, seed, tspan, times = NULL,
                            summary = NULL, keep = TRUE, cores = 1, ...) {
  plan <- if (inherits(x, "sim_model")) sim_compile(x) else check_plan(x)
  if (!is_count(n)) sim_abort("`n` must be a positive integer.")
  if (!is.function(sampler)) sim_abort("`sampler` must be a function of the run index.")
  if (!is_scalar_number(seed)) sim_abort("`seed` must be a number.")
  sets <- with_seed(seed, lapply(seq_len(n), function(i) {
    s <- sampler(i)
    check_param_list(s, "`sampler()` result")
    s
  }))
  out <- run_sets(plan, sets, tspan, times, summary, keep, cores, ...)
  out$seed <- seed
  out
}

parameter_sets <- function(parameters) {
  if (is.data.frame(parameters)) {
    if (!nrow(parameters)) sim_abort("`parameters` has no rows.")
    return(lapply(seq_len(nrow(parameters)), function(i) {
      as.list(parameters[i, , drop = FALSE])
    }))
  }
  if (is.list(parameters) && length(parameters) && all(vapply(parameters, is.list, TRUE))) {
    return(parameters)
  }
  sim_abort("`parameters` must be a data frame or a list of named lists.")
}

run_sets <- function(plan, sets, tspan, times, summary, keep, cores, ...) {
  if (!is.null(summary) && !is.function(summary)) sim_abort("`summary` must be a function.")
  if (!is_count(cores)) sim_abort("`cores` must be a positive integer.")
  one <- function(i) {
    p <- sim_set_params(plan, sets[[i]])
    r <- tryCatch(sim_run(p, tspan, times = times, on_failure = "warn", ...),
                  error = function(e) e)
    if (inherits(r, "error")) {
      return(list(result = NULL, status = "failed", summary = NULL,
                  message = conditionMessage(r)))
    }
    s <- if (!is.null(summary) && r$success) summary(r) else NULL
    list(result = if (keep) r, status = r$status, summary = s, message = r$message)
  }
  idx <- seq_along(sets)
  runs <- if (cores > 1L && .Platform$OS.type != "windows") {
    parallel::mclapply(idx, one, mc.cores = cores)
  } else {
    lapply(idx, one)
  }
  pars <- do.call(rbind, lapply(sets, function(s) {
    as.data.frame(lapply(s, function(v) if (length(v) == 1L && is.atomic(v)) v else I(list(v))),
                  stringsAsFactors = FALSE)
  }))
  status <- vapply(runs, `[[`, "", "status")
  summ <- NULL
  if (!is.null(summary)) {
    vals <- lapply(runs, `[[`, "summary")
    template <- Find(Negate(is.null), vals)
    if (!is.null(template)) {
      vals <- lapply(vals, function(v) if (is.null(v)) template * NA else v)
      summ <- cbind(pars, as.data.frame(do.call(rbind, vals)), status = status)
    }
  }
  structure(list(parameters = pars,
                 results = if (keep) lapply(runs, `[[`, "result"),
                 status = status,
                 messages = vapply(runs, function(r) r$message %||% "", ""),
                 summary = summ),
            class = "sim_sweep")
}

#' @export
print.sim_sweep <- function(x, ...) {
  cat("<sim_sweep> ", length(x$status), " run(s): ",
      paste(names(table(x$status)), table(x$status), sep = " = ", collapse = ", "), "\n", sep = "")
  if (!is.null(x$summary)) print(utils::head(x$summary, 10L))
  invisible(x)
}
