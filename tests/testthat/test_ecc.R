#' Checks that the ECC output is a valid clique cover.
#' @param am Adjacency matrix of graph
#' @param cliques Edge clique cover
check_edge_clique_cover <- function(am, cliques) {
  # am2 is a copy to see if we can uncover the edges
  am2 <- am
  for (j in 1:dim(am2)[1]) {
    am[j, j] = T
  }
  for (cn in 1:length(cliques)) {
    lc <- cliques[[cn]]
    # Mark these edges as covered
    for (j in lc) {
      for (j2 in lc) {
        am2[j, j2] <- F
      }
    }
    # check that lc is a valid clique
    testthat::expect_true(
      all(am[lc,lc] == 1)
    )
  }
  testthat::expect_true(
    all(am2 == 0)
  )
}

test_that('edge clique cover yields valid results', {
  # Generate fake data
  set.seed(123)
  thresh <- 0.5; p <- 500
  data <- generate_regression_data(p=p, n=2*p)
  V <- cor(data$X)
  flags <- V >= thresh
  for (j in 1:p) {
    flags[j,j] = F
  }
  G <- igraph::graph_from_adjacency_matrix(flags, mode='undirected')

  # test that ecc gives the correct outcome
  cliques <- edge_clique_cover(G)
  #print(lapply(cliques, function(x) {sort(as.numeric(x))}))
  check_edge_clique_cover(flags, cliques)

  # Check that this is a good solution
  max_cliques <- igraph::max_cliques(G)
  #print("Max cliques are:")
  #print(lapply(max_cliques, function(x) {sort(as.numeric(x))}))
  nmc <- length(max_cliques)
  testthat::expect_true(
    length(cliques) <= nmc
  )
})

