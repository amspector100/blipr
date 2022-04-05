#' #' @keywords internal
# normalize_locs <- function(locs) {
#   flags <- !is.na(locs)
#   shift <- min(locs[flags])
#   if (shift == -Inf) {
#     stop("infinite values are not permitted in locs")
#   }
#   locs <- locs - shift
#   scale <- max(locs[flags])
#   if (scale == Inf) {
#     stop("infinite values are not permitted in locs")
#   }
#   if (scale == 0) {
#     stop("scale == 0: locs appear to have only one unique value.")
#   }
#   locs <- locs / scale
#   return(list(locs=locs, shift=shift, scale=scale))
# }

#' @param M1 a k1 x d matrix
#' @param M2 a k2 x d matrix
#' @returns a k1 x k2 x d array of differences
#' @keywords internal
broadcast_diffs <- function(M1, M2) {
  k1 <- dim(M1)[1]
  k2 <- dim(M2)[1]
  d <- dim(M1)[2]
  diffs <- array(NaN, c(k1, k2, d))
  for (j in 1:k1) {
    diffs[j,,] <- sweep(M2, 2, M1[j,])
  }
  return(-1*diffs)
}

#' Maps (center, rad) pairs to string keys.
#' @param centers 1 x k x d array of centers
#' @param rad radius
#' @keywords internal
centerrad2key <- function(centers, rad, dec=NULL) {
  # perhaps rounding, the degree of rounding depends on rad
  if (is.null(dec)) {dec <- max(ceiling(-1*log(rad, base=10) + 3))}
  centers <- round(centers, digits=dec)
  rad <- round(rad, digits=dec)
  keys <- apply(centers, 2, function(x) {paste(x, collapse=',')})
  if (length(rad) == 1) {
    keys <- lapply(keys, function(x) {paste(x, rad, sep=',')})
  } else {
    keys <- lapply(1:length(keys), function(i) {paste(keys[i], rad[i], sep=',')})
  }
  return(keys)
}

#' Inverse function of centerrad2key
#' @param x list of keys
#' @keywords internal
key2centerrad <- function(keys) {
  # Split by separator
  nkeys <- length(keys)
  splits <- sapply(
    keys, function(x) {base::strsplit(x, split=',')[[1]]}
  )
  d <- dim(splits)[1] - 1
  # Recreate matrix of centers / vector of radii
  radii <- as.numeric(splits[(d+1),])
  centers <- t(matrix(
    as.numeric(splits[-(d+1),]), d, nkeys
  ))
  return(list(centers=centers, radii=radii))
}

#' Finds additional centers for 2d circles
#' @param samples 1 x k x d array of post. samples
#' @param centers 1 x k x d array of centers
#' @param radius radius of circles
#' @param gsize grid size
#' @keywords internal
additional_circle_centers <- function(
  samples, centers, radius, gsize, TOL=1e-5
) {
  k <- dim(samples)[2]
  d <- dim(samples)[3]
  radius2 <- radius * radius
  centers_vec <- c()
  # Check upper/lower/left/right
  for (j in 1:2) {
    for (offset in c(-1, 1)) {
      centers_new <- centers # automatically creates a copy of inds we modify
      centers_new[,,j] <- centers_new[,,j] + offset / gsize
      # Check if points are in offset centers
      included <- base::rowSums((samples - centers_new)**2, dims=2) <= radius2 # 1 x k
      included <- which(included)
      # Add to list
      if (length(included) > 0) {
        centers_vec <- unique(c(
          centers_vec,
          centerrad2key(centers_new[,included,,drop=F], radius)
        ))
      }
    }
  }
  return(centers_vec)
}

#' @keywords internal
find_centers <- function(samples, gsize, shape) {
  # samples is 1 x k x d
  k = dim(samples)[2]
  if (k == 0) {
    return(as.character(c()))
  }
  d = dim(samples)[3]
  # lower-left corner
  corners <- floor(samples * gsize) / gsize

  # Adjust to make centers
  centers <- corners + 1 / (2 * gsize)
  if (shape == 'circle') {
    if (d != 2) {
      stop(paste("Shape=circle only implemented for 2-d data, detected", d, "dimensions"))
    }
    radius <- sqrt(d) / (2 * gsize)
  } else if (shape == 'square') {
    radius <- 1 / (2*gsize)
  } else {
    stop(paste("Unrecognized shape=", shape, ", must be one of square/circle", sep=''))
  }

  # Add extra centers for circles only
  centers_keys <- unique(centerrad2key(centers, radius))
  if (shape == 'circle') {
    additional_centers <- additional_circle_centers(
      samples, centers, radius, gsize
    )
    centers_keys <- unique(c(centers_keys, additional_centers))
  }
  return(centers_keys)
}

#' Computes contiguous candidate groups on a lattice R^d
#'
#' @param locs A (N, num_disc, d)-dimensional array. Here, N is the
#' number of samples from the posterior, d is the number
#' of dimensions of the space, and each point corresponds
#' to a signal in a particular posterior sample.
#' @param grid_sizes List of grid sizes to split up the locations.
#' The grid size is inversely proportional to the distance between
#' lattice points.
#' @param extra_centers An (ncenters, d)-dimensional matrix. At
#' each resolution, candidate groups will be computed with centers
#' at this location.
#' @param max_pep The maximum allowable PEP for output candidate groups.
#' Defaults to 0.25.
#' @param shape One of ('circle', 'square').
#' @param log_interval if not NULL, will log progress every log_interval
#' updates.
#'
#' @return list of candidate groups
#' @export
lattice_peps <- function(
  locs,
  grid_sizes,
  extra_centers=NULL,
  max_pep=0.25,
  shape='square',
  log_interval=NULL
) {
  # Normalize locs to be in [0,1]^d
  #norm <- normalize_locs(locs)
  #locs <- norm$locs
  N <- dim(locs)[1]
  k <- dim(locs)[2]
  d <- dim(locs)[3]
  # preprocess extra centers
  if (! is.null(extra_centers)) {
    nextra <- dim(extra_centers)[1]
    norm_ecent <- extra_centers #(extra_centers - norm$shift) / norm$scale
  }
  # Create PIPs
  pips <- list()
  for (j in 1:N) {
    if (! is.null(log_interval)) {
      if (j %% log_interval == 0) {
        cat("Starting sample j=", j, " out of ", N, ".\n", sep='')
      }
    }
    # Ignore NaNs
    samples <- locs[j,,]
    active <- which(apply(
      ! is.nan(abind::adrop(locs[j,,,drop=F], drop=1)),
      1, all
    ))
    samples <- locs[j,active,,drop=F] # 1 x k x d
    # Loop through grid sizes and find centers
    all_centers <- c()
    for (gsize in grid_sizes) {
      all_centers <- c(all_centers, find_centers(samples, gsize, shape))
    }
    # Repeat for manually added centers
    all_extra_centers <- c()
    if (! is.null(extra_centers)) {
      dists <- broadcast_diffs(
        norm_ecent,
        abind::adrop(samples, drop=1)
      ) # nextra x k x d
      if (shape == 'square') {
        dists <- apply(abs(dists), c(1,2), max)
      } else {
        dists <- sqrt(apply(dists**2, c(1,2), sum))
      }
      min_dists <- apply(dists, 1, min)
      for (gsize in grid_sizes) {
        radius <- 1 / (2 * gsize)
        if (shape == 'circle') {
          radius <- sqrt(d) * radius
        }
        ncs <- which(min_dists <= radius)
        if (length(ncs) > 0) {
          new_keys <- centerrad2key(
            array(norm_ecent[ncs,], c(1,length(ncs), d)),
            radius
          )
          all_extra_centers <- c(all_extra_centers, new_keys)
        }
      }
    }
    # Update PIPs
    all_centers <- unique(c(all_centers, all_extra_centers))
    for (key in all_centers) {
      if (is.null(pips[[key]])) {
        pips[[key]] <- 1 / N
      } else {
        pips[[key]] <- pips[[key]] + 1 / N
      }
    }
  }
  # Filter pips
  filtered_peps <- list()
  for (key in names(pips)) {
    pep <- 1 - pips[[key]]
    if (pep <= max_pep) {
      filtered_peps[[key]] <- pep
    }
  }
  return(filtered_peps)
}

#' Turns the output of the lattice_peps function into
#' a list of list of candidate groups. Each sub-list
#' corresponds to a list of completely disconnected
#' candidate groups which can be fed to BLiP separately
#' (this saves computation).
#'
#' @param filtered_peps An output of the lattice_peps function
#' @param min_blip_size Combines connected components so all
#' subproblems are at least this size.
#' @param verbose If true, will give progress reports over time.
#' @param shape One of 'square' or 'circle'.
#' @param max_pep The maximum pep for candidate groups.
#'
#' @return A list of list of candidate groups.
#'
#' @export
lattice_peps_to_cand_groups <- function(
  filtered_peps,
  min_blip_size=5000,
  verbose=F,
  shape='square',
  max_pep=1
) {
  if (! requireNamespace("igraph", quietly=F)) {
    stop("lattice_peps_to_cand_groups requires 'igraph' to be installed.")
  }

  # Step 0: create PEPs, centers, radii
  peps <- as.numeric(unlist(filtered_peps))
  keys <- names(filtered_peps)
  out <- key2centerrad(keys)
  centers <- out$centers # n x d
  d <- dim(centers)[2]
  radii <- out$radii

  # Step 1: create adjacency matrix
  ngroups <- length(radii)
  if (ngroups > 50000) {
    msg <- paste("Computing adjancency matrix may be expensive for ngroups=", ngroups, sep='')
    msg <- paste(msg, ". Try decreasing max_pep.", sep='')
    warning(msg)
  }

  # Step 2: Construct constraints
  if (verbose) {
    cat("Constructing constraint matrix with ngroups=", ngroups, ".\n", sep='')
  }
  deltas <- broadcast_diffs(
    matrix(radii, ngroups, 1), matrix(-1*radii, ngroups, 1)
  )[,,1]
  diffs <- broadcast_diffs(
    centers, centers
  )
  if (shape == 'square'){
    dists <- apply(abs(diffs), c(1,2), max)
  } else if (shape == 'circle') {
    dists <- sqrt(apply(diffs**2, c(1,2), sum))
  } else {
    stop(paste("Unrecognized shape=", shape, ", must be one of square/circle", sep=''))
  }
  constraints <- dists < deltas
  # save memory
  rm(dists)
  rm(deltas)

  # Step 3: Split problem into connected components
  if (verbose) {
    cat("Constructing graph and isolating connected components.\n")
  }
  G <- igraph::graph_from_adjacency_matrix(constraints, mode='undirected')
  components <- igraph::components(G)
  nc <- length(unique(components$membership))
  merged_components <- list(c())
  nmc <- 1
  for (j in 1:nc) {
    comp <- which(components$membership == j)
    if (length(merged_components[[nmc]]) > min_blip_size) {
      nmc <- nmc + 1
      merged_components[[nmc]] <- comp
    } else {
      merged_components[[nmc]] <- c(merged_components[[nmc]], comp)
    }
  }
  rm(G) # save memory
  if (verbose) {
    cat("Decomposed into ", nmc, " components. Constructing cand groups.\n", sep='')
  }

  # Step 4: construct cand groups for BLiP
  all_cand_groups <- list()
  for (compnum in 1:nmc) {
    # Create subgraph
    component <- merged_components[[compnum]]
    component_cand_groups <- list()
    component_groups <- list()
    subG <- igraph::graph_from_adjacency_matrix(
      constraints[component, component], mode='undirected'
    )
    # Run edge clique cover
    if (verbose) {
      cat("Finding edge clique cover for component ", compnum, " of ", nmc, ".\n", sep='')
    }
    cliques <- edge_clique_cover(subG)
    maxj <- max(sapply(cliques, max))
    # initialize
    component_groups[[maxj+1]] <- c(-1)
    for (cn in 1:length(cliques)) {
      for (j in cliques[[cn]]) {
        if (length(component_groups[[j]]) == 0) {
          component_groups[[j]] <- c(cn)
        } else {
          component_groups[[j]] <- c(component_groups[[j]], cn)
        }
      }
    }
    if (verbose) {
      cat("Found clique cover of length ", length(cliques), " for component ", compnum, ", now constructing groups.\n", sep='')
    }
    for (ii in 1:length(component)) {
      group <- as.vector(unlist(component_groups[ii]))
      j <- component[ii]
      pep <- peps[j]
      component_cand_groups <- c(
        component_cand_groups,
        list(list(
          pep=pep,
          group=group,
          center=centers[j,],
          radius=radii[j]
        ))
      )
    }
    all_cand_groups[[compnum]] <- component_cand_groups
  }
  return(list(
    all_cand_groups=all_cand_groups,
    merged_components=merged_components
  ))
}
