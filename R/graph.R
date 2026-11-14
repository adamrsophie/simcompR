# Graph algorithms on integer vertex ids 1..n with adjacency lists.

# Tarjan's strongly connected components. Returns a list of integer vectors.
# Vertices are visited in increasing id order, so the result is deterministic.
tarjan_scc <- function(n, adj) {
  index <- integer(n); index[] <- NA_integer_
  low <- integer(n)
  on_stack <- logical(n)
  stack <- integer()
  counter <- 0L
  comps <- list()
  # iterative depth-first search to avoid deep recursion on long chains
  for (root in seq_len(n)) {
    if (!is.na(index[root])) next
    work <- list(c(root, 1L))
    while (length(work)) {
      top <- work[[length(work)]]
      v <- top[1]; i <- top[2]
      if (i == 1L) {
        counter <- counter + 1L
        index[v] <- counter; low[v] <- counter
        stack <- c(stack, v); on_stack[v] <- TRUE
      }
      nb <- adj[[v]]
      if (i <= length(nb)) {
        w <- nb[i]
        work[[length(work)]] <- c(v, i + 1L)
        if (is.na(index[w])) {
          work[[length(work) + 1L]] <- c(w, 1L)
        } else if (on_stack[w]) {
          low[v] <- min(low[v], index[w])
        }
        next
      }
      work[[length(work)]] <- NULL
      if (length(work)) {
        parent <- work[[length(work)]][1]
        low[parent] <- min(low[parent], low[v])
      }
      if (low[v] == index[v]) {
        pos <- match(v, stack)
        comp <- stack[pos:length(stack)]
        stack <- stack[seq_len(pos - 1L)]
        on_stack[comp] <- FALSE
        comps[[length(comps) + 1L]] <- sort(comp)
      }
    }
  }
  comps
}

# Kahn topological sort of the condensation. Ties are broken by the smallest
# vertex id, which makes the order a function of the graph only.
condensed_order <- function(n, adj, comps) {
  comp_of <- integer(n)
  for (k in seq_along(comps)) comp_of[comps[[k]]] <- k
  nc <- length(comps)
  cadj <- vector("list", nc)
  indeg <- integer(nc)
  for (v in seq_len(n)) {
    for (w in adj[[v]]) {
      a <- comp_of[v]; b <- comp_of[w]
      if (a != b && !(b %in% cadj[[a]])) {
        cadj[[a]] <- c(cadj[[a]], b)
        indeg[b] <- indeg[b] + 1L
      }
    }
  }
  key <- vapply(comps, min, integer(1))
  ready <- which(indeg == 0L)
  out <- integer()
  while (length(ready)) {
    k <- ready[which.min(key[ready])]
    ready <- ready[ready != k]
    out <- c(out, k)
    for (b in cadj[[k]]) {
      indeg[b] <- indeg[b] - 1L
      if (indeg[b] == 0L) ready <- c(ready, b)
    }
  }
  comps[out]
}

# Shortest cycle through the first vertex of a strongly connected
# component (breadth-first search), for error messages.
find_cycle <- function(comp, adj) {
  start <- comp[1]
  parent <- integer(max(comp)); parent[] <- NA_integer_
  queue <- start
  while (length(queue)) {
    v <- queue[1]; queue <- queue[-1]
    for (w in intersect(adj[[v]], comp)) {
      if (w == start) {
        path <- v
        while (path[1] != start) path <- c(parent[path[1]], path)
        return(c(path, start))
      }
      if (is.na(parent[w]) && w != start) {
        parent[w] <- v
        queue <- c(queue, w)
      }
    }
  }
  c(start, start)
}
