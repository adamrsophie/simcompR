# Graph algorithms on integer vertex ids 1..n with adjacency lists.

# Tarjan's strongly connected components. Returns a list of integer vectors.
# Vertices are visited in increasing id order, so the result is deterministic.
# The depth-first search is iterative (no deep recursion on long chains) and
# all stacks are preallocated, so the cost is linear in vertices and edges.
tarjan_scc <- function(n, adj) {
  index <- rep(NA_integer_, n)
  low <- integer(n)
  on_stack <- logical(n)
  stack <- integer(n); sp <- 0L
  work_v <- integer(n); work_i <- integer(n); wp <- 0L
  counter <- 0L
  comps <- vector("list", n); nc <- 0L
  for (root in seq_len(n)) {
    if (!is.na(index[root])) next
    wp <- 1L; work_v[1L] <- root; work_i[1L] <- 1L
    while (wp > 0L) {
      v <- work_v[wp]; i <- work_i[wp]
      if (i == 1L) {
        counter <- counter + 1L
        index[v] <- counter; low[v] <- counter
        sp <- sp + 1L; stack[sp] <- v; on_stack[v] <- TRUE
      }
      nb <- adj[[v]]
      if (i <= length(nb)) {
        w <- nb[i]
        work_i[wp] <- i + 1L
        if (is.na(index[w])) {
          wp <- wp + 1L; work_v[wp] <- w; work_i[wp] <- 1L
        } else if (on_stack[w] && index[w] < low[v]) {
          low[v] <- index[w]
        }
        next
      }
      wp <- wp - 1L
      if (wp > 0L) {
        parent <- work_v[wp]
        if (low[v] < low[parent]) low[parent] <- low[v]
      }
      if (low[v] == index[v]) {
        pos <- sp
        while (stack[pos] != v) pos <- pos - 1L
        comp <- stack[pos:sp]
        sp <- pos - 1L
        on_stack[comp] <- FALSE
        nc <- nc + 1L
        comps[[nc]] <- sort(comp)
      }
    }
  }
  comps[seq_len(nc)]
}

# Kahn topological sort of the condensation. Ties are broken by the smallest
# vertex id, which makes the order a function of the graph only. The ready
# components are kept in a binary heap keyed by that id.
condensed_order <- function(n, adj, comps) {
  comp_of <- integer(n)
  for (k in seq_along(comps)) comp_of[comps[[k]]] <- k
  nc <- length(comps)
  from <- comp_of[rep.int(seq_len(n), lengths(adj))]
  to <- comp_of[unlist(adj, use.names = FALSE)]
  keep <- from != to
  from <- from[keep]; to <- to[keep]
  dup <- duplicated(from * (nc + 1) + to)
  from <- from[!dup]; to <- to[!dup]
  cadj <- split(to, factor(from, levels = seq_len(nc)))
  indeg <- tabulate(to, nbins = nc)
  key <- vapply(comps, min, integer(1))
  heap <- integer(nc); hn <- 0L
  push <- function(k) {
    hn <<- hn + 1L; i <- hn; heap[i] <<- k
    while (i > 1L) {
      p <- i %/% 2L
      if (key[heap[p]] <= key[heap[i]]) break
      tmp <- heap[p]; heap[p] <<- heap[i]; heap[i] <<- tmp; i <- p
    }
  }
  pop <- function() {
    top <- heap[1L]
    heap[1L] <<- heap[hn]; hn <<- hn - 1L
    i <- 1L
    repeat {
      l <- 2L * i; r <- l + 1L; m <- i
      if (l <= hn && key[heap[l]] < key[heap[m]]) m <- l
      if (r <= hn && key[heap[r]] < key[heap[m]]) m <- r
      if (m == i) break
      tmp <- heap[m]; heap[m] <<- heap[i]; heap[i] <<- tmp; i <- m
    }
    top
  }
  # a sorted vector is already a heap
  start <- which(indeg == 0L)
  hn <- length(start)
  heap[seq_len(hn)] <- start[order(key[start])]
  out <- integer(nc); no <- 0L
  while (hn > 0L) {
    k <- pop()
    no <- no + 1L; out[no] <- k
    for (b in cadj[[k]]) {
      indeg[b] <- indeg[b] - 1L
      if (indeg[b] == 0L) push(b)
    }
  }
  comps[out[seq_len(no)]]
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
