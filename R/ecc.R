#' Finds an edge clique cover of a graph G.
#'
#' @param G igraph graph
#' @returns cliques A list of cliques. Each clique is
#' a vector integers where the integers correspond to
#' the nodes of G.
edge_clique_cover <- function(G) {
  igraph_available <- requireNamespace("igraph", quietly=T)
  if (! igraph_available) {
    stop("edge_clique_cover functionality requires 'igraph' to be installed.")
  }
  G <- igraph::simplify(G) # remove repeat edges and self loops
  # TODO: if implementing "count_signals" functionality,
  # it's crucial to add all self loops. See the python implementation.
  nodes <- as.numeric(igraph::V(G))
  n <- length(nodes)
  degrees <- igraph::degree(G)
  cliques <- list()
  covered <- as.integer(rep(0, length(igraph::E(G))))
  # iterate through edges
  for (edgenum in 1:length(covered)) {
    # Check if edge hasbeen covered
    if (covered[edgenum]) {
      next
    }
    # Signal that edge will be covered
    covered[edgenum] <- T
    e <- igraph::ends(G, edgenum)
    v1 <- e[1]
    v2 <- e[2]
    degrees[c(v1,v2)] <- degrees[c(v1,v2)] - 1

    # else, create a new clique
    clique <- c(v1, v2)
    neighbors <- igraph::neighbors(G, v1)
    neighbors <- intersect(neighbors, igraph::neighbors(G, v2))
    #neighbors <- setdiff(neighbors, c(v1, v2))
    # if (any(degrees[neighbors] > 0)) {
    while (length(neighbors) > 0) {
      vnew <- neighbors[which.max(degrees[neighbors])]
      # Update covered and degrees
      for (v in clique) {
        enum <- igraph::get.edge.ids(G, c(v, vnew), error=T)
        if (! covered[enum]) {
          degrees[c(v, vnew)] <- degrees[c(v, vnew)] - 1
          covered[enum] <- T
        }
      }
      # Add vnew to clique and update neighbors
      neighbors <- intersect(neighbors, igraph::neighbors(G, vnew))
      #neighbors <- setdiff(neighbors, c(vnew))
      clique <- c(clique, c(vnew))
    }
    #}
    cliques[[length(cliques)+1]] <- unique(clique)
  }
  return(cliques)
}
