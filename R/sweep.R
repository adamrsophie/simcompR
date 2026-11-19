#' Parameter sweeps and Monte Carlo experiments
#'
#' `sim_sweep()` runs one simulation per row of `parameters`, a data frame
#' (or list of named lists) of model parameter values. The model is
#' compiled once; each run uses [sim_set_params()], which reuses the
#' dependency analysis whenever the structure is unchanged.
#'
#' `sim_monte_carlo()` draws `n` parameter sets with `sampler(i)` under a
#' fixed `seed` and runs them like a sweep. The sampler runs with the
#' Mersenne-Twister, Inversion and Rejection generators whatever the
#' caller's [RNGkind()], so the draws are reproducible. Afterwards the
#' caller's `.Random.seed` and generator kinds are restored. One exception
#' cannot be repaired from R: with `normal.kind = "Box-Muller"`, R keeps the
#' second deviate of each pair outside `.Random.seed`, and that cached value
#' is lost, so the caller's next normal draw differs; a warning says so.
#' `white_noise` blocks do not use R's generator at all.
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
#'   numeric vector with unique, nonempty names. Values are matched to
#'   columns by name, never by position: a run that returns the names in
#'   another order keeps their meaning, a run that omits a name gets `NA`
#'   there, and a new name adds a column (`NA` for the other runs). A run
#'   whose summary is not a named numeric vector is marked `"failed"` with
#'   the reason in `$messages`.
#' @param keep Keep the full results (`TRUE`) or only the summaries.
#' @param cores Number of parallel workers.
#' @param ... Further arguments for [sim_run()].
#' @return A `sim_sweep` object with elements `parameters`, `results` (when
#'   `keep = TRUE`), `status`, `messages` and `summary`. `parameters` has one
#'   row per run, a `run` column and one column per parameter that any run
#'   sets; a run that does not set a parameter shows the model's value.
#'   Vector-valued parameters are list columns. `summary` is `parameters`
#'   joined with the summary values and the run status.
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
  if (!is_scalar_number(seed) || seed != round(seed) || abs(seed) > .Machine$integer.max) {
    sim_abort("`seed` must be a whole number in the integer range; set.seed() would ",
              "otherwise give distinct seeds the same stream.")
  }
  if (RNGkind()[2] == "Box-Muller") {
    sim_warn("normal.kind is \"Box-Muller\": its cached deviate cannot be restored, ",
             "so the next normal draw after sim_monte_carlo() differs from an ",
             "uninterrupted stream. .Random.seed and RNGkind() are restored.")
  }
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
    # list columns hold one value per row
    return(lapply(seq_len(nrow(parameters)), function(i) {
      lapply(parameters, function(col) if (is.list(col)) col[[i]] else col[i])
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
  # a misspelled parameter is a mistake in the design, not a failed run
  unknown <- setdiff(unlist(lapply(sets, names)), names(plan$params))
  if (length(unknown)) {
    sim_abort("Unknown model parameter(s): ", paste(unique(unknown), collapse = ", "),
              ". Model parameters: ", paste(names(plan$params), collapse = ", "), ".")
  }
  clash <- intersect(unlist(lapply(sets, names)), c("run", "status"))
  if (length(clash)) {
    sim_abort("Parameter(s) ", paste(unique(clash), collapse = ", "), " cannot be varied in ",
              "a sweep: `run` and `status` name the columns of the result.")
  }
  run_args <- list(...)
  if (is.null(run_args$on_failure)) run_args$on_failure <- "warn"
  one <- function(i) {
    fail <- function(e) list(result = NULL, status = "failed", summary = NULL,
                             message = conditionMessage(e))
    r <- tryCatch(do.call(sim_run, c(list(sim_set_params(plan, sets[[i]]), tspan,
                                          times = times), run_args)),
                  error = function(e) e)
    if (inherits(r, "error")) return(fail(r))
    s <- NULL
    if (!is.null(summary) && r$success) {
      s <- tryCatch(check_summary(summary(r)), error = function(e) e)
      if (inherits(s, "error")) {
        return(list(result = if (keep) r, status = "failed", summary = NULL,
                    message = paste("summary failed:", conditionMessage(s))))
      }
    }
    list(result = if (keep) r, status = r$status, summary = s, message = r$message)
  }
  idx <- seq_along(sets)
  runs <- if (cores > 1L && .Platform$OS.type != "windows") {
    parallel::mclapply(idx, one, mc.cores = cores)
  } else {
    lapply(idx, one)
  }
  pars <- effective_parameters(plan, sets)
  status <- vapply(runs, `[[`, "", "status")
  summ <- NULL
  if (!is.null(summary)) {
    vals <- lapply(runs, `[[`, "summary")
    cols <- unique(unlist(lapply(vals, names)))
    if (length(cols)) {
      both <- intersect(cols, names(pars))
      if (length(both)) {
        sim_abort("Summary name(s) ", paste(both, collapse = ", "), " are also parameter ",
                  "names; rename them in `summary` so each column means one thing.")
      }
      tab <- matrix(NA_real_, length(vals), length(cols), dimnames = list(NULL, cols))
      for (i in seq_along(vals)) if (length(vals[[i]])) tab[i, names(vals[[i]])] <- vals[[i]]
      summ <- cbind(pars, as.data.frame(tab, optional = TRUE), status = status)
    }
  }
  structure(list(parameters = pars,
                 results = if (keep) lapply(runs, `[[`, "result"),
                 status = status,
                 messages = vapply(runs, function(r) r$message %||% "", ""),
                 summary = summ),
            class = "sim_sweep")
}

# A summary is a named numeric vector; its names identify the columns.
check_summary <- function(s) {
  if (!(is.numeric(s) || is.logical(s)) || !is.null(dim(s))) {
    stop("it must return a named numeric vector, not a ", class(s)[1], ".", call. = FALSE)
  }
  nm <- names(s)
  if (length(s) && (is.null(nm) || anyNA(nm) || any(!nzchar(nm)))) {
    stop("every value it returns must have a name.", call. = FALSE)
  }
  if (anyDuplicated(nm)) {
    stop("it returned duplicate names: ", paste(unique(nm[duplicated(nm)]), collapse = ", "),
         ".", call. = FALSE)
  }
  if (any(nm %in% c("run", "status"))) {
    stop("the names `run` and `status` are reserved.", call. = FALSE)
  }
  stats::setNames(as.numeric(s), nm)
}

# One row per run with the value each parameter actually had: the run's own
# value, or the model's when the run does not set it.
effective_parameters <- function(plan, sets) {
  cols <- unique(unlist(lapply(sets, names)))
  out <- data.frame(run = seq_along(sets))
  for (nm in cols) {
    v <- lapply(sets, function(s) if (nm %in% names(s)) s[[nm]] else plan$params[[nm]])
    scalar <- all(vapply(v, function(z) is.atomic(z) && length(z) == 1L, TRUE))
    out[[nm]] <- if (scalar) unlist(v) else I(v)
  }
  out
}

#' @export
print.sim_sweep <- function(x, ...) {
  cat("<sim_sweep> ", length(x$status), " run(s): ",
      paste(names(table(x$status)), table(x$status), sep = " = ", collapse = ", "), "\n", sep = "")
  if (!is.null(x$summary)) print(utils::head(x$summary, 10L))
  invisible(x)
}
