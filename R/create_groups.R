#' Calculates all sequential candidate groups below max_size.
#'@param inclusions (N,p)-shaped matrix of posterior samples
#'where a nonzero value indicates the presence of a signal.
#'@param q The nominal level at which to control the error rate
#' (optional)
#'@param max_pep The maximum posterior error probability
#' (PEP) allowed in a candidate group. Default is 1.
#'@param max_size maximum allowable size for each group.
#'@param prenarrow If true, prenarrows the candidate groups
#' as described in the paper. Defaults to TRUE.
#'
#'@export
sequential_groups <- function(
  inclusions,
  q=0,
  max_pep=1,
  max_size=25,
  prenarrow=TRUE
) {
  # preprocessing
  inclusions <- inclusions != 0
  N <- dim(inclusions)[1]
  p <- dim(inclusions)[2]
  max_size <- min(max_size, p)
  # Precompute cumulative sums for speed
  cum_incs = cbind(
    matrix(0, N, 1),
    t(apply(inclusions, 1, cumsum))
  )

  # Compute successive groups of size m
  all_peps <- list()
  for (m in 1:max_size) {
    cum_diffs <- cum_incs[, (m+1):(p+1)] - cum_incs[, 1:(p-m+1)]
    if (m == p) {
        all_peps[[m]] <- mean(cum_diffs == 0)
      } else {
        all_peps[[m]] <- colMeans(cum_diffs == 0)
      }
  }

  # each index is the first (smallest) variable in the group
  # which has size m
  active_inds = list()
  for (m in 1:max_size) {
    active_inds[[m]] <- which(all_peps[[m]] < max_pep)
  }

  # Prenarrowing: this iteratively updates elim_inds
  # so that when considering the set of groups of size m,
  # elim_inds are all the indices that are redundant
  if (prenarrow & q > 0) {
    elim_inds = which(all_peps[[1]] < q)
    for (m in 2:max_size) {
      # if index j is eliminated for level m-1,
      # indexes j, j-1 are eliminated for level m
      elim_inds = union(elim_inds, elim_inds - 1)
      # update active_inds[m]
      active_inds[[m]] = setdiff(active_inds[[m]], elim_inds)
      # Update elim_inds (groups with PEP < q)
      elim_inds = union(
        elim_inds, which(all_peps[[m]] < q)
      )
    }
  }

  # Create candidate groups
  cand_groups = list()
  ncg = 1
  for (m in 1:max_size) {
    for (ind in active_inds[[m]]) {
      group <- c(ind:(ind+m-1))
      pep <- all_peps[[m]][ind]
      cand_groups[[ncg]] = list(
        group=group,
        pep=pep
      )
      ncg <- ncg + 1
    }
  }
  return(cand_groups)
}

#' After prefiltering groups, some features/locations
#' may not appear in any candidate groups. When this happens,
#' this function reindexes the locations to improve efficiency.
elim_redundant_features <- function(cand_groups) {
  # Step 1: find relevant features
  active_features <- Reduce(
    union, lapply(cand_groups, function(x) {x$group})
  )
  # Step 2: change indexes to save computation
  nrel = length(active_features)
  orig2new = rep(0, nrel)
  for (i in 1:nrel) {
    orig2new[active_features[i]] <- i
  }
  # Step 3: reindex
  for (i in 1:length(cand_groups)) {
    blip_group <- sapply(
      cand_groups[[i]]$group, function(x) {orig2new[x]}
    )
    cand_groups[[i]]$blip_group <- blip_group
  }
  return(cand_groups)
}

#' Creates groups based on dist_matrix using hierarchical clustering
dist_matrix_to_groups <- function(dist_matrix) {
  groups <- list()
  for (method in c('single', 'average', 'complete')) {
    # Create clustering
    tree <- stats::hclust(dist_matrix, method=method)
    p <- length(tree$height) + 1
    # TODO: may be a more efficient version of this somewhere
    all_groupings <- sapply(1:p, function(k) {stats::cutree(tree, k=k)})
    for (k in 1:p) {
      # Turn the hclust style vector specifyer into a list of the group indices
      new_groups <- lapply(
        unique(all_groupings[,k]), 
        function(j) {which(all_groupings[,k] == j)}
      ) 
      groups <- c(groups, new_groups)
    }
  }
  return(unique(groups))
}

#' Creates hierarchically structured candidate groups
#' based on dist_matrix.
#' @param inclusions (N,p)-shaped array of posterior samples
#' where a nonzero value indicates the presence of a signal.
#' @param dist_matrix A distance matrix corresponding to
#' distances between locations, used for hierarchical clustering.
#' @param X The design matrix in regression problems, which will
#' be used to create dist_matrix if dist_matrix is not provided.
#'@param max_pep The maximum posterior error probability
#' (PEP) allowed in a candidate group. Default is 1.
#'@param max_size maximum allowable size for each group.
#'@param filter_sequential If TRUE, ignore sequential groups
#' of variables to avoid duplication.
#' @export
hierarchical_groups <- function(
  inclusions,
  dist_matrix=NULL,
  X=NULL,
  max_pep=1,
  max_size=25,
  filter_sequential=FALSE
) {
  # Preprocessing
  N <- dim(inclusions)[1]
  p <- dim(inclusions)[2]
  inclusions <- inclusions != 0
  # Trivial cases with zero/one feature
  if (p == 0) {return(list())}
  if (p == 1) {
    pep <- 1 - mean(inclusions)
    return(list(list(pep=pep, group=c(1))))
  }
  # Estimate cov matrix from inclusions/X
  if (is.null(dist_matrix)) {
    if (! is.null(X)) {
      dist_matrix <- as.dist(1 - abs(stats::cor(X)))
    } else {
      # Ensure standard deviation is not zero for any cols
      nvals <- apply(inclusions, 2, function(x) length(unique(x)))
      if (any(nvals <= 1)) {
        precorr <- rbind(rep(1, N), rep(0, N), inclusions)
      } else {
        precorr <- inclusions
      }
      # negative correlations = super close together
      dist_matrix <- as.dist(1 + stats::cor(precorr))
    }
  }
  # Initialize output, create groups
  cand_groups <- list()
  ncg <- 1
  groups <- dist_matrix_to_groups(dist_matrix)
  # Create PEPs
  for (group in groups) {
    gsize <- length(group)
    if (gsize > max_size) {next}
    # possibly filter out sequential groups
    if (filter_sequential) {
      if (max(group) - min(group) == gsize - 1) {next}
    }
    # Calculate PEP
    if (gsize > 1) {
      pep = 1 - mean(apply(inclusions[,group], 1, any))
    } else {
      pep = 1 - mean(inclusions[,group])
    }
    if (pep < max_pep) {
      cand_groups[[ncg]] <- list(
        pep=pep,
        group=group
      )
      ncg <- ncg + 1
    }
  }
  return(cand_groups)
}