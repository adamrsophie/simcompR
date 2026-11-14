`%||%` <- function(a, b) if (is.null(a)) b else a

sim_abort <- function(..., call = NULL) {
  stop(simpleError(paste0(...), call = call))
}

sim_warn <- function(...) {
  warning(paste0(...), call. = FALSE)
}

is_scalar_number <- function(x) {
  is.numeric(x) && length(x) == 1L && !is.na(x) && is.finite(x)
}

is_count <- function(x) {
  is_scalar_number(x) && x >= 1 && x == round(x)
}

check_number <- function(x, what, finite = TRUE) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || (finite && !is.finite(x))) {
    sim_abort(what, " must be a single ", if (finite) "finite ", "number.")
  }
  invisible(as.numeric(x))
}

check_numeric <- function(x, what) {
  if (!is.numeric(x) || length(x) < 1L || anyNA(x) || any(!is.finite(x))) {
    sim_abort(what, " must be a non-empty finite numeric vector.")
  }
  if (is.matrix(x)) {
    sim_abort(what, " must be a vector, not a matrix.")
  }
  invisible(as.numeric(x))
}

check_matrix <- function(x, what, nrow = NULL, ncol = NULL) {
  if (is.numeric(x) && !is.matrix(x) && length(x) == 1L) x <- matrix(x, 1L, 1L)
  if (!is.matrix(x) || !is.numeric(x) || anyNA(x) || any(!is.finite(x))) {
    sim_abort(what, " must be a finite numeric matrix.")
  }
  if (!is.null(nrow) && nrow(x) != nrow) {
    sim_abort(what, " must have ", nrow, " rows, not ", nrow(x), ".")
  }
  if (!is.null(ncol) && ncol(x) != ncol) {
    sim_abort(what, " must have ", ncol, " columns, not ", ncol(x), ".")
  }
  storage.mode(x) <- "double"
  x
}

check_sample_time <- function(period, offset, what) {
  check_number(period, paste0(what, ": `period`"))
  check_number(offset, paste0(what, ": `offset`"))
  if (period <= 0) sim_abort(what, ": `period` must be positive.")
  if (offset < 0 || offset >= period) {
    sim_abort(what, ": `offset` must satisfy 0 <= offset < period.")
  }
  c(period = as.numeric(period), offset = as.numeric(offset))
}

# Absolute time tolerance used to merge clock hits and boundaries.
# About 4500 units in the last place: far above the rounding error of
# offset + k * period, far below any meaningful time step.
time_tol <- function(tspan) {
  1e-12 * max(1, abs(tspan))
}

# R binds a named argument to a formal by partial matching when no formal
# matches exactly, so a block parameter `n` would bind to `name`. This
# re-matches the call with exact names only. It returns NULL when nothing
# was matched partially, otherwise list(fixed = formals, dots = the rest).
exact_match_args <- function(call, formal_names, env) {
  args <- as.list(call)[-1L]
  if (!length(args)) return(NULL)
  nms <- names(args) %||% character(length(args))
  exact <- nzchar(nms) & nms %in% formal_names
  free <- setdiff(formal_names, nms[exact])
  partial <- nzchar(nms) & !exact &
    vapply(nms, function(n) any(startsWith(free, n)), TRUE)
  if (!any(partial)) return(NULL)
  if (any(vapply(args, identical, TRUE, quote(...)))) return(NULL)
  vals <- lapply(args, eval, envir = env)
  fixed <- list()
  dots <- list()
  for (i in seq_along(vals)) {
    if (exact[i]) {
      fixed[nms[i]] <- list(vals[[i]])
    } else if (!nzchar(nms[i]) && length(free)) {
      fixed[free[1]] <- list(vals[[i]])
      free <- free[-1L]
    } else {
      dots <- c(dots, vals[i])
    }
  }
  list(fixed = fixed, dots = dots)
}

# Merge sorted times that agree within tol. A cluster starts at its first
# time and never extends beyond tol from it, so close times cannot chain
# together. The cluster is represented by its largest time: a sample hit
# computed as 0.8999999999999999 and a step at 0.9 meet at 0.9, where the
# step has already switched.
merge_times <- function(times, tol) {
  times <- sort(unique(times))
  if (length(times) < 2L) return(times)
  out <- numeric()
  start <- -Inf
  for (v in times) {
    if (v - start > tol) {
      out <- c(out, v)
      start <- v
    } else {
      out[length(out)] <- v
    }
  }
  out
}

# Index of the element of sorted `ref` within tol of each value, else NA.
match_tol <- function(v, ref, tol) {
  if (!length(v)) return(integer())
  if (!length(ref)) return(rep(NA_integer_, length(v)))
  i <- pmax(findInterval(v, ref), 1L)
  j <- pmin(i + 1L, length(ref))
  pick <- ifelse(abs(ref[j] - v) < abs(ref[i] - v), j, i)
  ifelse(abs(ref[pick] - v) <= tol, pick, NA_integer_)
}

pkg_version <- function(pkg) {
  tryCatch(as.character(utils::packageVersion(pkg)), error = function(e) NA_character_)
}

# Names for a vector signal: "name" when scalar, "name[i]" otherwise.
element_names <- function(name, n) {
  if (n == 1L) name else paste0(name, "[", seq_len(n), "]")
}

# ---------------------------------------------------------------------------
# Portable expressions
#
# Expressions written as strings (event resets, state machine guards and
# actions, math functions) are parsed once and checked against a whitelist
# of functions. Reading a model file therefore never runs arbitrary R code.

safe_functions <- c(
  "+", "-", "*", "/", "^", "%%", "%/%", "(", "{", "[", "[[", "<-", "=",
  "==", "!=", "<", ">", "<=", ">=", "&", "|", "!", "&&", "||",
  "abs", "sqrt", "exp", "log", "log10", "log2", "log1p", "expm1",
  "sin", "cos", "tan", "asin", "acos", "atan", "atan2", "sinh", "cosh",
  "tanh", "asinh", "acosh", "atanh", "floor", "ceiling", "round", "trunc",
  "sign", "min", "max", "pmin", "pmax", "sum", "prod", "mean", "length",
  "c", "ifelse", "if", "rep", "seq_len", "numeric", "pi", "TRUE", "FALSE",
  "is.na", "as.numeric", "cumsum", "rev", "norm", "crossprod", "%*%", "t"
)

#' Parse a portable expression
#'
#' Besides checking called functions against the whitelist, the expression
#' must call functions only by name (no computed call heads such as
#' `(f)(x)`) and may assign only to plain variables or their elements,
#' never to a name of a permitted function.
#'
#' @param text A single string.
#' @param what Description used in error messages.
#' @return A parsed expression.
#' @noRd
parse_safe <- function(text, what) {
  if (!is.character(text) || length(text) != 1L || is.na(text)) {
    sim_abort(what, " must be a single string.")
  }
  expr <- tryCatch(str2lang(text), error = function(e) {
    sim_abort(what, " could not be parsed: ", conditionMessage(e))
  })
  check_safe_expr(expr, what)
  expr
}

check_safe_expr <- function(expr, what) {
  if (!is.call(expr)) return(invisible())
  fn <- expr[[1L]]
  if (!is.symbol(fn)) {
    sim_abort(what, " calls a computed function; only named functions are allowed ",
              "in portable expressions.")
  }
  name <- as.character(fn)
  if (!name %in% safe_functions) {
    sim_abort(what, " uses functions that are not allowed in portable ",
              "expressions: ", name, ".")
  }
  if (name %in% c("<-", "=")) {
    target <- expr[[2L]]
    while (is.call(target) && is.symbol(target[[1L]]) &&
           as.character(target[[1L]]) %in% c("[", "[[")) {
      lapply(as.list(target)[-c(1L, 2L)], check_safe_expr, what = what)
      target <- target[[2L]]
    }
    if (!is.symbol(target)) {
      sim_abort(what, " assigns to something other than a variable.")
    }
    if (as.character(target) %in% c(safe_functions, "u", "t", "x", "s", "elapsed")) {
      sim_abort(what, " must not assign to '", as.character(target), "'.")
    }
    check_safe_expr(expr[[3L]], what)
    return(invisible())
  }
  lapply(as.list(expr)[-1L], check_safe_expr, what = what)
  invisible()
}

# Environment holding only the permitted functions and constants. Its parent
# is the empty environment, so nothing else in R is reachable by name.
safe_env <- local({
  env <- NULL
  function() {
    if (is.null(env)) {
      env <<- new.env(parent = emptyenv())
      for (f in setdiff(safe_functions, c("TRUE", "FALSE"))) {
        assign(f, get(f, envir = baseenv()), envir = env)
      }
      lockEnvironment(env, bindings = TRUE)
    }
    env
  }
})

safe_frame <- function(vars) {
  list2env(vars, parent = safe_env())
}

eval_safe <- function(expr, vars) {
  eval(expr, safe_frame(vars))
}
