schema_name <- "simcompR.model"
schema_version <- 1L

#' Save and load models as JSON
#'
#' Models are written as versioned, declarative JSON: block types, names,
#' parameters, connections, events and nested subsystems. Numbers are
#' stored as strings with 17 significant digits, so they read back
#' bit-for-bit; integer, logical and character vectors keep their type,
#' missing values and names, and numeric matrices keep their dimensions.
#' A value that could not be read back unchanged (an array with more than
#' two dimensions, a matrix that is not numeric or has dimension names, a
#' factor, a data frame, or any other value with a class or attributes
#' other than names) is rejected when the model is written, instead of
#' being saved in a different form.
#'
#' Only portable content can be written. Block types defined by R functions
#' (`fcn`, `ode`, `signal_function`), user types registered with
#' `portable = FALSE`, parameters holding functions and event resets given
#' as functions are rejected with an error naming the offending block,
#' rather than silently dropped. Use the string forms (`expression` blocks,
#' string resets and guards), which are parsed against a whitelist of
#' functions when the model is read, so loading a file never runs arbitrary
#' R code.
#'
#' @param model A [sim_model()].
#' @param path File path.
#' @param x A list produced by `sim_model_to_list()`.
#' @param pretty Indent the JSON.
#' @return `sim_write_json()` returns `path` invisibly; `sim_read_json()`
#'   and `sim_model_from_list()` return a `sim_model`; `sim_model_to_list()`
#'   returns a list.
#' @examples
#' m <- sim_model("feedback", params = list(k = 2))
#' sim_add_block(m, "step", "r", time = 1)
#' sim_add_block(m, "subtract", "e")
#' sim_add_block(m, "gain", "g", gain = sim_param("k"))
#' sim_add_block(m, "integrator", "x")
#' sim_connect(m, "r", "e", input = 1)
#' sim_connect(m, "x", "e", input = 2)
#' sim_connect(m, "e", "g")
#' sim_connect(m, "g", "x")
#' f <- tempfile(fileext = ".json")
#' sim_write_json(m, f)
#' m2 <- sim_read_json(f)
#' r1 <- sim_run(sim_compile(m), c(0, 3), times = 3)
#' r2 <- sim_run(sim_compile(m2), c(0, 3), times = 3)
#' identical(sim_signal(r1, "x"), sim_signal(r2, "x"))
#' @export
sim_write_json <- function(model, path, pretty = TRUE) {
  x <- sim_model_to_list(model)
  writeLines(jsonlite::toJSON(x, auto_unbox = TRUE, pretty = pretty, null = "null"), path)
  invisible(path)
}

#' @rdname sim_write_json
#' @export
sim_read_json <- function(path) {
  x <- jsonlite::fromJSON(path, simplifyVector = FALSE)
  sim_model_from_list(x)
}

#' @rdname sim_write_json
#' @export
sim_model_to_list <- function(model) {
  check_model(model)
  check_acyclic(model)
  blocks <- lapply(model$blocks, function(b) {
    type_def <- get_block_type(b$type)
    if (!isTRUE(type_def$portable)) {
      sim_abort("Block '", b$name, "' has type '", b$type, "', which is defined by R ",
                "functions and cannot be written to a portable model file.")
    }
    list(name = b$name, type = b$type,
         params = encode_params(b$params, paste0("Block '", b$name, "'")))
  })
  cn <- model$connections
  conns <- lapply(seq_len(nrow(cn)), function(k) {
    list(from = cn$from[k], output = cn$from_port[k], to = cn$to[k], input = cn$to_port[k])
  })
  events <- lapply(model$events, function(e) {
    resets <- lapply(names(e$reset), function(b) {
      r <- e$reset[[b]]
      if (!is.character(r)) {
        sim_abort("Event '", e$name, "': the reset of '", b, "' is an R function and ",
                  "cannot be written to a portable model file; use a string expression.")
      }
      list(block = b, expr = r)
    })
    tr <- e$trigger
    trig <- if (tr$kind == "scheduled") list(kind = "scheduled", times = enc(tr$times, "times"))
      else list(kind = "zero_crossing", signal = tr$signal, direction = tr$direction,
                level = enc(tr$level, "level"), port = tr$port, element = tr$element)
    list(name = e$name, trigger = trig, reset = resets, terminal = e$terminal,
         priority = enc(e$priority, "priority"))
  })
  list(schema = schema_name, schema_version = schema_version, name = model$name,
       params = encode_params(model$params, paste0("Model '", model$name, "'")),
       blocks = unname(blocks), connections = conns, events = unname(events))
}

#' @rdname sim_write_json
#' @export
sim_model_from_list <- function(x) {
  if (!is.list(x) || !identical(x$schema, schema_name)) {
    sim_abort("Not a simcompR model file (missing schema '", schema_name, "').")
  }
  v <- x$schema_version
  if (!is.numeric(v) || length(v) != 1L || is.na(v) || v != round(v) || v < 1) {
    sim_abort("Model file has an invalid schema version; expected a whole number ",
              "from 1 to ", schema_version, ".")
  }
  if (v > schema_version) {
    sim_abort("Model file uses schema version ", x$schema_version,
              "; this version of simcompR reads version ", schema_version, " and older.")
  }
  m <- sim_model(x$name, decode_params(x$params))
  for (b in x$blocks) {
    do.call(sim_add_block, c(list(model = m, type = b$type, name = b$name), decode_params(b$params)))
  }
  for (cn in x$connections) {
    sim_connect(m, cn$from, cn$to, output = cn$output, input = cn$input)
  }
  for (e in x$events) {
    tr <- e$trigger
    trig <- if (tr$kind == "scheduled") sim_at(dec(tr$times)) else
      sim_zero_crossing(tr$signal, tr$direction, dec(tr$level), tr$port, tr$element)
    reset <- list()
    for (r in e$reset) reset[[r$block]] <- r$expr
    sim_add_event(m, e$name, trig, reset = reset, terminal = isTRUE(e$terminal),
                  priority = dec(e$priority))
  }
  m
}

encode_params <- function(params, where) {
  out <- lapply(names(params), function(n) enc(params[[n]], paste0(where, " parameter '", n, "'")))
  names(out) <- names(params)
  if (!length(out)) return(structure(list(), names = character()))
  out
}

decode_params <- function(x) {
  if (!length(x)) return(list())
  out <- lapply(x, dec)
  names(out) <- names(x)
  out
}

num_str <- function(v) {
  out <- sprintf("%.17g", v)
  out[is.na(v)] <- "NA"
  out[is.nan(v)] <- "NaN"
  out
}

enc <- function(v, where) {
  if (is.null(v)) return(list(`$null` = TRUE))
  if (is_param_ref(v)) {
    return(if (is.null(v$expr)) list(`$param` = v$name) else list(`$param_expr` = v$expr))
  }
  if (inherits(v, "sim_model")) return(list(`$model` = sim_model_to_list(v)))
  if (inherits(v, "sim_transition")) {
    if (is.function(v$guard) || is.function(v$action)) {
      sim_abort(where, ": transitions with R function guards or actions cannot be ",
                "written to a portable model file; use string expressions.")
    }
    return(list(`$transition` = list(from = v$from, to = v$to, guard = v$guard,
                                     after = enc(v$after, where),
                                     action = enc(v$action, where),
                                     priority = enc(v$priority, where))))
  }
  if (is.function(v)) {
    sim_abort(where, " is an R function and cannot be written to a portable model file.")
  }
  extra <- setdiff(names(attributes(v)), "names")
  if (is.matrix(v) && is.double(v) && identical(extra, "dim")) {
    return(list(`$matrix` = list(nrow = nrow(v), ncol = ncol(v), data = as.list(num_str(as.vector(v))))))
  }
  if (length(extra)) {
    what <- if (length(dim(v)) > 2L) "an array with more than two dimensions"
      else if (is.matrix(v)) "a matrix that is not a plain numeric matrix"
      else paste0("a value with attribute(s) ", paste(extra, collapse = ", "))
    sim_abort(where, " is ", what, " (class '", class(v)[1], "'), which cannot be ",
              "written to a portable model file without changing it.")
  }
  nm <- names(v)
  if (is.integer(v)) {
    out <- list(`$int` = as.list(num_str(v)))
  } else if (is.double(v)) {
    out <- list(`$num` = as.list(num_str(v)))
  } else if (is.character(v)) {
    # missing elements are JSON nulls
    out <- list(`$chr` = lapply(v, function(s) if (is.na(s)) NULL else s))
  } else if (is.logical(v)) {
    out <- list(`$lgl` = as.list(ifelse(is.na(v), "NA", ifelse(v, "TRUE", "FALSE"))))
  } else if (is.list(v)) {
    items <- lapply(seq_along(v), function(i) enc(v[[i]], where))
    out <- list(`$list` = items)
  } else {
    sim_abort(where, " has class '", class(v)[1], "', which cannot be written to JSON.")
  }
  if (!is.null(nm)) out$names <- lapply(nm, function(s) if (is.na(s)) NULL else s)
  out
}

# JSON array of strings and nulls to a character vector, by position
dec_chr <- function(items) {
  vapply(items, function(s) if (is.null(s)) NA_character_ else as.character(s), "")
}

dec <- function(x) {
  if (!is.list(x) || !length(x)) return(x)
  key <- names(x)[1]
  nm <- if (!is.null(x$names)) dec_chr(x$names) else NULL
  v <- switch(key,
    `$null` = NULL,
    `$param` = sim_param(x[[1]]),
    `$param_expr` = sim_param_expr(x[[1]]),
    `$model` = sim_model_from_list(x[[1]]),
    `$transition` = {
      t <- x[[1]]
      sim_transition(t$from, t$to, guard = t$guard, after = dec(t$after),
                     action = dec(t$action), priority = dec(t$priority))
    },
    `$matrix` = matrix(as.numeric(unlist(x[[1]]$data)), x[[1]]$nrow, x[[1]]$ncol),
    `$num` = suppressWarnings(as.numeric(unlist(x[[1]]))),
    `$int` = suppressWarnings(as.integer(unlist(x[[1]]))),
    `$chr` = dec_chr(x[[1]]),
    `$lgl` = as.logical(unlist(x[[1]])),
    `$list` = lapply(x[[1]], dec),
    sim_abort("Unknown encoded value '", key, "' in model file.")
  )
  if (key == "$num" && length(v) == 0L) v <- numeric()
  if (key == "$int" && length(v) == 0L) v <- integer()
  if (key == "$chr" && length(v) == 0L) v <- character()
  if (key == "$list" && length(v) == 0L) v <- list()
  if (!is.null(nm) && !is.null(v)) names(v) <- nm
  v
}
