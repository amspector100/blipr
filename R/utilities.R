#' Generate synthetic sparse regression data.
#'
#'@param n Number of data points
#'@param p Dimensionality
#'@param a Correlations are sampled as Beta(a,b) random variables.
#'@param b Correlations are sampled as Beta(a,b) random variables.
#'@param sparsity Proportion of nonzero coefficients.
#'@param coeff_size Size of coefficients
#'
#'@export
generate_regression_data <- function(
  n=100, p=500, a=5, b=1, sparsity=0.05, coeff_size=1
) {
  # Create design matrix
  X <- matrix(stats::rnorm(n*p), n, p)
  rhos <- stats::rbeta(p-1, a, b)
  for (j in 2:p) {
    X[,j] <- sqrt(1 - rhos[j-1]**2) * X[,j] + rhos[j-1] * X[,j-1]
  }

  # Create coefficients
  beta <- rep(0, p)
  sp <- ceiling(sparsity * p)
  inds <- sample(1:p, size=sp, replace=F)
  beta[inds] <-sqrt(coeff_size) * stats::rnorm(sp)

  # Create y
  y <- X %*% beta + stats::rnorm(n)
  return(list(X=X, y=y, beta=beta))
}

#' Generate synthetic change point data.
#'
#'@param n Number of data points
#'@param sparsity Proportion of change points
#'@param tau2 Size of change points
#'
#'@export
generate_changepoint_data <- function(
  n=100, sparsity=0.01, tau2=1
) {
  # Create change points
  beta <- rep(0, n)
  sp <- ceiling(sparsity * n)
  inds <- sample(1:n, size=sp, replace=F)
  beta[inds] <-sqrt(tau2) * stats::rnorm(sp)

  # Create data
  mu <- cumsum(beta)
  y <- mu + stats::rnorm(n)
  return(list(y=y, mu=mu, beta=beta))
}
