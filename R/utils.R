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
time_tol <- function(tspan) {
  1e-10 * max(1, abs(tspan))
}

# Merge sorted times that agree within tol. Keeps the first of each cluster.
merge_times <- function(times, tol) {
  times <- sort(unique(times))
  if (length(times) < 2L) return(times)
  keep <- c(TRUE, diff(times) > tol)
  # cluster representatives: first element of each run
  times[keep]
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
  bad <- setdiff(called_functions(expr), safe_functions)
  if (length(bad)) {
    sim_abort(what, " uses functions that are not allowed in portable ",
              "expressions: ", paste(bad, collapse = ", "), ".")
  }
  expr
}

called_functions <- function(expr) {
  if (is.call(expr)) {
    fn <- expr[[1L]]
    head <- if (is.symbol(fn)) as.character(fn) else called_functions(fn)
    c(head, unlist(lapply(as.list(expr)[-1L], called_functions)))
  } else {
    character()
  }
}

eval_safe <- function(expr, vars) {
  env <- list2env(vars, parent = baseenv())
  eval(expr, env)
}
