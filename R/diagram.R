#' Draw or export a block diagram
#'
#' `sim_to_dot()` writes the model as a Graphviz DOT graph. Edges into
#' inputs with direct feedthrough are solid; edges into inputs without
#' direct feedthrough (integrators, delays) are dashed, which makes valid
#' dynamic feedback visible. `plot()` draws the same graph with 'igraph'
#' when it is installed.
#'
#' @param x,model A [sim_model()].
#' @param ... Passed to `igraph::plot.igraph()`.
#' @return `sim_to_dot()` returns a string; `plot()` returns `x` invisibly.
#' @examples
#' m <- sim_model("feedback")
#' sim_add_block(m, "constant", "r")
#' sim_add_block(m, "subtract", "e")
#' sim_add_block(m, "integrator", "x")
#' sim_connect(m, "r", "e", input = 1)
#' sim_connect(m, "x", "e", input = 2)
#' sim_connect(m, "e", "x")
#' cat(sim_to_dot(m))
#' if (requireNamespace("igraph", quietly = TRUE)) plot(m)
#' @export
sim_to_dot <- function(model) {
  check_model(model)
  ft <- edge_feedthrough(model)
  cn <- model$connections
  esc <- function(s) gsub('"', '\\\\"', s)
  nodes <- vapply(names(model$blocks), function(b) {
    sprintf('  "%s" [shape=box, label="%s\\n(%s)"];', esc(b), esc(b), model$blocks[[b]]$type)
  }, "")
  edges <- vapply(seq_len(nrow(cn)), function(k) {
    sprintf('  "%s" -> "%s" [label="%d:%d"%s];', esc(cn$from[k]), esc(cn$to[k]),
            cn$from_port[k], cn$to_port[k], if (isFALSE(ft[k])) ", style=dashed" else "")
  }, "")
  paste(c(sprintf('digraph "%s" {', esc(model$name)), "  rankdir=LR;", nodes, edges, "}"),
        collapse = "\n")
}

edge_feedthrough <- function(model) {
  cn <- model$connections
  vapply(seq_len(nrow(cn)), function(k) {
    b <- model$blocks[[cn$to[k]]]
    def <- tryCatch(build_block(get_block_type(b$type),
                                resolve_params(b$params, list(model$params), b$name), b$name),
                    error = function(e) NULL)
    if (is.null(def)) NA else def$feedthrough[cn$to_port[k]]
  }, TRUE)
}

#' @rdname sim_to_dot
#' @export
plot.sim_model <- function(x, ...) {
  if (!requireNamespace("igraph", quietly = TRUE)) {
    sim_abort("Plotting a model needs the 'igraph' package; use sim_to_dot() otherwise.")
  }
  cn <- x$connections
  g <- igraph::graph_from_data_frame(
    data.frame(from = cn$from, to = cn$to, stringsAsFactors = FALSE),
    vertices = data.frame(name = names(x$blocks), stringsAsFactors = FALSE))
  ft <- edge_feedthrough(x)
  igraph::plot.igraph(g, layout = igraph::layout_with_sugiyama(g)$layout,
                      vertex.shape = "rectangle", vertex.size = 30, vertex.size2 = 12,
                      vertex.color = "#F4EFE6", vertex.label.color = "#13213A",
                      edge.lty = ifelse(ft %in% FALSE, 2, 1),
                      edge.arrow.size = 0.4, main = x$name, ...)
  invisible(x)
}
