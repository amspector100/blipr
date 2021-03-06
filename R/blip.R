default_solver <- function() {
	solvers <- CVXR::installed_solvers()
	if ('GUROBI' %in% solvers) { return('GUROBI') }
	if ('CBC' %in% solvers) {return('CBC')}
	warning("Using ECOS solver, which will be slightly slower. Consider installing CBC.")
	return("ECOS")
}

BINARY_TOL <- 1e-3
ERROR_OPTIONS <- c("fdr", "local_fdr", "fwer", "pfer")
WEIGHT_FN_STRS <- c("inverse_size", "log_inverse_size")
DEFAULT_GRID_SIZES <- exp(log(10)*seq(log10(4),4, length.out=25))

#' Given samples from a posterior or a list of candidate groups, BLiP performs
#' resolution-adaptive signal detection to maximize power while controlling
#' (e.g.) the FDR.
#'
#' Note: when working with image data or a continuous set
#' of locations, consider using BLiP_cts.
#'
#' @param samples (N,p)-shaped matrix of posterior samples where a nonzero
#' value indicates the presence of a signal.
#' @param cand_groups A list of lists, where the inner lists must have a
#' "group" attribute, corresponding to the features in the group
#' and a "pep" attribute, corresponding to a posterior error probability.
#' @param error Bayesian error rate to control: one of "fwer", "pfer", "fdr", "local_fdr".
#' @param q The level at which to control the Bayesian FWER/PFER/FDR/local FDR.
#' @param max_pep Never select any group with a pep greater than max_pep.
#' @param deterministic Whether or not BLiP should return a deterministic solution.
#' @param weight_fn How to weight discoveries. Can be one of 'inverse_size'
#' or 'log_inverse_size' or a function which takes a candidate group as an input
#' and returns a weight.
#' @param verbose If TRUE, gives occasional progress reports.
#' @param perturb If TRUE, adds a tiny (random) perturbation to the weights to ensure
#' the existence of a unique optimal solution.
#' @param max_iters Maximum number of binary-search iterations for FWER when.
#' @param search_method For FWER control, how to find the optimal parameter for
#' the LP. Either 'binary' (defalt) or 'none'.
#' @param solver The solver to use within CVXR. By default, will use Gurobi, CBC,
#' or ECOS (in that order), depending on whether they are installed.
#' @return A list of candidate groups, which asserts that each group contains
#' a signal.
#'
#' @examples
#' # Example 1: sparse linear regression
#' set.seed(123); n <- 100; p <- 200
#' data <- blipr::generate_regression_data(n=n, p=p)
#' # sample from the posterior, e.g., using NPrior
#' nprior <- NPrior::NPrior_run(
#'   X=data$X, y=data$y, N=5000, prior='SpSL-G'
#' )
#' # run blip on posterior samples
#' detections <- blipr::BLiP(samples=t(nprior$ThetaSamples), q=0.1, error='fdr')
#'
#' # Example 2: Running BLiP directly on candidate groups
#' cand_groups <- list(
#'   list(group=c(1), pep=0.1),
#'   list(group=c(2), pep=0.5),
#'   list(group=c(1,2), pep=0.01)
#'  )
#' detections <- blipr::BLiP(cand_groups=cand_groups, q=0.1, error='fdr')
#' @export
BLiP <- function(
  samples=NULL,
  cand_groups=NULL,
  weight_fn='inverse_size',
  error='fdr',
  q=0.1,
  max_pep=1,
  deterministic=T,
  verbose=F,
  perturb=T,
  max_iters=100,
  search_method='binary',
  solver=NULL
) {
  # Preprocessing
  error <- tolower(error)
  if (! error %in% ERROR_OPTIONS) {
    stop(paste("error (", error, ") must be one of", paste(ERROR_OPTIONS, collapse=', ')))
  }
  if (is.null(cand_groups) & is.null(samples)) {
    stop("At least one of cand_groups and samples must be provided.")
  }
  if (error %in% c("fwer", "pfer", "local_fdr")) {max_pep <- min(max_pep, q)}
  if (is.null(solver)) {solver <- default_solver()}
  logLevel <- ifelse(verbose, 1, 0)

  # Create cand_groups if necessary
  if (is.null(cand_groups)) {
    seq_groups <- sequential_groups(
      samples=samples, q=q, max_pep=max_pep, prenarrow=T
    )
    hier_groups <- hierarchical_groups(
      samples=samples, max_pep=max_pep, filter_sequential=T
    )
    cand_groups <- c(seq_groups, hier_groups)
  }

  # Prefilter and prune
  cand_groups <- prefilter(cand_groups, max_pep=max_pep)
  if (length(cand_groups) == 0) {return(list())} # edge case with no rejections
  cand_groups <- elim_redundant_features(cand_groups)
  nrel <- length(Reduce(union, lapply(cand_groups, function(x) x$group)))

  # Weights for each candidate group
  ngroups <- length(cand_groups)
  if (methods::is(weight_fn, 'character')) {
    if (weight_fn == 'prespecified') {
      weights <- sapply(cand_groups, function(x) x$weight)
    } else {
      if (weight_fn == 'inverse_size') {
        weight_fn <- inverse_size
      } else if (weight_fn == 'log_inverse_size') {
        weight_fn <- log_inverse_size
      } else {
        stop(
          paste("Unrecognized weight_fn (", weight_fn, ") must be one of", paste(WEIGHT_FN_STRS, collapse=', '))
        )
      }
      weights <- sapply(cand_groups, weight_fn)
    }
  } else {
    weights <- sapply(cand_groups, weight_fn)
  }

  # Perturb to ensure unique solution
  orig_weights <- weights
  if (perturb) {
    weights <- weights * (1 + 0.0001 * stats::runif(ngroups))
  }

  # Extract peps
  peps <- sapply(cand_groups, function(x) x$pep)

  # Constraints to ensure selected groups are disjoint
  if (verbose) {
    cat("BLiP problem has", ngroups, "groups in contention, with", nrel, "active features/locations.\n")
  }
  A <- matrix(0, ngroups, nrel)
  for (gj in 1:length(cand_groups)) {
    for (loc in cand_groups[[gj]]$blip_group) {
      A[gj, loc] = 1
    }
  }

  # Determine if a binary search for FWER is necessary
  binary_search <- (! is.null(samples)) & (error == 'fwer') & (search_method == 'binary')

  # Assemble variables and constraints
  x <- CVXR::Variable(ngroups)
  v_param <- CVXR::Parameter()
  v_var <- CVXR::Variable(pos=TRUE) # for FDR only
  objective <- CVXR::Maximize(sum(weights * x * (1 - peps)))
  b <- rep(1, nrel)
  constraints <- list(
    x >= 0,
    x <= 1,
    t(A) %*% x <= b
  )
  if (error %in% c("pfer", "fwer")) {
    CVXR::value(v_param) <- q
    constraints <- c(constraints, list(sum(x * peps) <= v_param))
  }
  if (error == 'fdr') {
    constraints <- c(constraints, list(
      sum(x * peps) <= v_var,
      sum(x) >= v_var / q,
      v_var >= 0,
      v_var <= q * nrel + 1
    ))
  }

  # Create problem
  problem <- CVXR::Problem(
    objective=objective, constraints=constraints
  )
  if (! binary_search) {
    result <- CVXR::solve(
      problem, solver=solver, verbose=verbose, logLevel=logLevel
    )
    selections <- as.numeric(result$getValue(x))
  } else {
    N <- dim(samples)[1]
    # min val not controlling FWER
    v_upper <- q * nrel + 1
    # max val controlling FWER
    v_lower <- 0
    for (niter in 1:max_iters) {
      # Set current value
      v <- (v_upper + v_lower) / 2
      CVXR::value(v_param) <- v
      # Solve
      result <- CVXR::solve(
        problem, solver=solver, warm_start=TRUE, verbose=verbose, logLevel=logLevel
      )
      # TODO could be more clever than this for FWER
      selections <- round(as.numeric(result$getValue(x)))
      # Calculate FWER for these selections
      false_disc <- rep(0, N)
      for (gj in which(selections > BINARY_TOL)) {
        group <- cand_groups[[gj]]$group
        if (length(group) == 1) {
          false_disc <- false_disc | samples[,group] == 0
        } else {
          false_disc <- false_disc | apply(samples[,group] == 0, 1, all)
        }
      }
      fwer <- mean(false_disc)
      if (fwer > q) {
        v_upper <- v
      } else {
        v_lower <- v
      }
      if (v_upper - v_lower < 1e-4) {
        break
      }
    }
    # Solve with v_lower for final solution
    CVXR::value(v_param) <- v_lower
    result <- CVXR::solve(
      problem, solver=solver, warm_start=TRUE, verbose=verbose, logLevel=logLevel
    )
    selections <- round(as.numeric(result$getValue(x)))
  }
  # Save information
  for (gj in 1:ngroups) {
    cand_groups[[gj]]$sprob <- selections[gj]
    cand_groups[[gj]]$weight <- orig_weights[gj]
  }
  if (error == 'fdr') {
    v_opt <- as.numeric(result$getValue(v_var))
  } else {
    v_opt <- as.numeric(CVXR::value(v_param))
  }

  return(binarize_selections(
    cand_groups=cand_groups,
    q=q,
    error=error,
    deterministic=deterministic
  ))
}

#' BLiP when the set of locations is continuous, e.g., when
#' working with image data.
#'
#' @param locs A (N, num_disc, d)-dimensional array. Here, N is the
#' number of samples from the posterior, d is the number
#' of dimensions of the space, and each point corresponds
#' to a signal in a particular posterior sample.
#' @param grid_sizes List of grid sizes to split up the locations.
#' The grid size is inversely proportional to the distance between
#' lattice points.
#' @param weight_fn Weight function which maps candidate groups
#' to weights. Defaults to inverse radius.
#' @param extra_centers An (ncenters, d)-dimensional matrix. At
#' each resolution, candidate groups will be computed with centers
#' at this location.
#' @param max_pep The maximum allowable PEP for output candidate groups.
#' Defaults to 0.25.
#' @param shape One of ('circle', 'square').
#' @param min_blip_size Combines connected components so all subproblems
#' are at least this size.
#' @param verbose If T, gives occasional status updates. Defaults to T.
#' @param ... Other arguments to the underlying BliP call,
#' such as the error rate or nominal level.
#'
#' @export
BLiP_cts <- function(
  locs,
  grid_sizes=DEFAULT_GRID_SIZES,
  weight_fn=inverse_radius,
  extra_centers=NULL,
  max_pep=0.25,
  shape='circle',
  min_blip_size=5000,
  verbose=T,
  ...
) {
  # 1. Calculate filtered PEPs
  N <- dim(locs)[1]
  d <- dim(locs)[2]
  if (verbose) {
    log_interval <- max(1, round(N / 10))
  } else {
    log_interval <- NULL
  }
  peps <- lattice_peps(
    locs=locs,
    grid_sizes=grid_sizes,
    extra_centers=extra_centers,
    max_pep=max_pep,
    shape=shape,
    log_interval=log_interval
  )

  # 2. Calculate candidate groups
  out <- lattice_peps_to_cand_groups(
    filtered_peps=peps,
    min_blip_size=min_blip_size,
    verbose=verbose,
    shape=shape,
    max_pep=max_pep
  )
  all_cand_groups <- out$all_cand_groups

  # 3. Run BLiP
  all_rej <- list()
  for (i in 1:length(all_cand_groups)) {
    rej <- BLiP(
      cand_groups=all_cand_groups[[i]],
      weight_fn=weight_fn,
      verbose=verbose,
      max_pep=max_pep,
      ...
    )
    all_rej <- c(all_rej, rej)
  }
  return(all_rej)
}

binarize_selections <- function(
  cand_groups,
  q,
  error,
  deterministic,
  tol=1e-3,
  nsample=10,
  verbose=F
) {
  output <- list()
  nontriv_cand_groups <- list() # non-integer solutions

  # Account for integer solutions
  for (cg in cand_groups) {
    if (cg$sprob < tol) {
      next
    } else if (cg$sprob > 1 - tol) {
      output <- c(output, list(cg))
    } else {
      nontriv_cand_groups <- c(nontriv_cand_groups, list(cg))
    }
  }

  # The easy cases: 0 or 1 non-integer values
  ngroups <- length(nontriv_cand_groups)
  if (ngroups == 0 | ngroups == 1) { return(output) }
  # if (ngroups == 1) {
  # 	if ((! deterministic) & stats::runif(1) < nontriv_cand_groups[[1]]$sprob) {
  # 		output <- c(output, nontriv_cand_groups)
  # 	}
  # 	return(output)
  # }

  # Preprocessing for hard cases (reindex to make the problem smaller)
  nontriv_cand_groups <- elim_redundant_features(nontriv_cand_groups)
  nrel <- length(Reduce(union, lapply(cand_groups, function(x) x$group)))
  # disjointness constraints
  A <- matrix(0, ngroups, nrel)
  for (gj in 1:ngroups) {
    for (loc in nontriv_cand_groups[[gj]]$blip_group) {
      A[gj, loc] = 1
    }
  }
  # The hard cases. Method 1: integer LP
  if (deterministic) {
    # construct integer linear program
    peps <- sapply(nontriv_cand_groups, function(x) x$pep)
    weights <- sapply(nontriv_cand_groups, function(x) x$weight)
    # Assemble variables and construct problem
    x <- CVXR::Variable(ngroups, integer=TRUE)
    objective <- CVXR::Maximize(sum(weights * x * (1 - peps)))
    constraints <- list(
      x >= 0,
      x <= 1,
      t(A) %*% x <= rep(1, nrel)
    )
    ndisc_output <- length(output)
    v_output <- ifelse(
      ndisc_output == 0,
      0,
      sum(sapply(output, function(x) x$pep))
    )
    # pfer / fwer specific constraints
    if (error %in% c("pfer", "fwer")) {
      if (error == 'pfer') {
        v_opt = q
      } else {
        v_opt = sum(sapply(cand_groups, function(x) {x$pep * x$sprob}))
      }
      v_new <- v_opt - v_output
      constraints <- c(constraints, list(sum(x * peps) <= v_new))
    }
    # No backtracking required for PFER, FWER, or local FDR
    logLevel = ifelse(verbose, 1, 0)
    if (error %in% c("pfer", "fwer", "local_fdr")) {
      problem <- CVXR::Problem(objective=objective, constraints=constraints)
      # We use CBC by default instead of GLPK:
      # See https://github.com/cvxpy/cvxpy/issues/1112
      if ('CBC' %in% CVXR::installed_solvers()) {
        result <- CVXR::solve(
          problem, solver='CBC', verbose=verbose, logLevel=logLevel
        )
      } else {
        result <- CVXR::solve(problem, verbose=verbose)
      }
    }
    if (error == 'fdr') {
      # Sort output in order of peps
      output <- output[order(sapply(output, function(x) {x$pep}))]
      # Iteratively try to solve problem and then backtrack if infeasible
      # (backtracking is extremely rare)
      while (length(output) >= 0) {
        v_var <- CVXR::Variable(pos=TRUE)
        constraints_fdr <- c(constraints, list(
          sum(x * peps) <= v_var,
          sum(x) >= (v_var + v_output) / q - ndisc_output
        ))
        problem <- CVXR::Problem(objective=objective, constraints=constraints_fdr)
        # Try to solve this problem with the default solver. This can fail in
        # small problem settings, in which case we try CBC if it is available.
        result <- tryCatch(
          {
            result <- CVXR::solve(
              problem, solver='GLPK', verbose=verbose, logLevel=logLevel
            )
          },
          error=function(cond) {
            if ('CBC' %in% CVXR::installed_solvers()) {
              return(CVXR::solve(
                problem, solver='CBC', verbose=verbose, logLevel=logLevel
              ))
            } else {
              stop(cond)
            }
          }
        )
        # Break if the solution is feasible, else backtrack.
        if (result$status != "infeasible") {
          break
        } else {
          # This should never be triggered (it is mathematically impossible)
          if (length(output) == 0) {
            warning("Backtracking for FDR control failed. Try installing more MILP solvers or setting deterministic=False.")
            output <- c()
            break
          }
          v_output <- v_output - output[[length(output)]]$pep
          ndisc_output <- ndisc_output - 1
          if (length(output) == 1) {
            output <- c()
          } else {
            output <- output[1:length(output)-1]
          }
        }
      }
    }
    # Add new solutions to output
    binarized_probs <- as.numeric(result$getValue(x))
    to_add <- nontriv_cand_groups[binarized_probs > 1 - tol]
    if (is.null(output)) { output <- to_add } else {output <- c(output, to_add)}
    # Method 2: sampling
  } else {
    # Add selection probs to A
    sprobs <- sapply(nontriv_cand_groups, function(x) x$sprob)
    A <- cbind(A, sprobs) # A[,nrel+1] corresponds to sprobs

    # Sort features in terms of marg. prob of selection
    marg_probs <- sapply(1:nrel, function(j) {sum(A[A[,j] == 1, nrel+1])})
    inds <- order(-1*marg_probs)

    # Initialize
    expected_powers <- rep(0, nsample)
    all_outputs <- list()
    for (ii in 1:nsample) {
      eliminated_groups <- rep(F, ngroups)
      selected_groups <- c()
      # Loop through by feature and sample
      for (feature in inds) {
        if (all(eliminated_groups)) { break }
        # Subset of available groups which contain the feature
        available.flags <- (A[,feature] == 1) & (! eliminated_groups)
        if (sum(available.flags) == 0) {
          subset <- NULL
        } else if (sum(available.flags) == 1) {
          subset <- matrix(A[available.flags,], nrow=1, ncol=nrel+1)
        } else {
          subset = A[available.flags,]
        }
        if (! is.null(subset)) {
          # Scale up cond probabilities
          prev_elim = A[(A[,feature] == 1) & (eliminated_groups),]
          if (is.null(dim(prev_elim))) {
            prev_elim <- matrix(prev_elim, ncol=length(prev_elim))
          }
          scale <- 1 - sum(prev_elim[,nrel+1])
          new_probs = subset[,nrel+1] / scale

          # select nothing with some probability
          if (stats::runif(1) <= 1 - sum(new_probs)) {
            eliminated_groups[A[,feature] == 1] = 1
		      } else {
		        # Else continue and select one of the groups
		        selected_group <- which(stats::rmultinom(1, 1, new_probs) == 1)
		        selected_group <- which(available.flags)[selected_group]
		        selected_groups <- c(selected_groups, selected_group)
		        # Eliminate groups mutually exclusive with the one we just selected
		        group_features <- which(A[selected_group,1:nrel] == 1)
		        if (length(group_features) == 1) {
		          new_elim_groups <- which(A[,group_features] != 0)
		        } else {
		          new_elim_groups <- which(rowSums(A[,group_features]) != 0)
		        }
		        eliminated_groups[new_elim_groups] = T
		      }
		    }
		    eliminated_groups <- as.integer(eliminated_groups)
		  }
		  output_ii <- c(output, nontriv_cand_groups[selected_groups])
		  # For local FDR, no backtracking required
		  if (error != "local_fdr" & length(output_ii) > 0) {
		    haterror <- extract_error(output_ii, error)
		    output_ii <- output_ii[order(sapply(output_ii, function(x) {x$pep}))]
		    while (haterror > q) {
		      if (length(output_ii) == 1) {
		        output_ii <- c()
		        break
		      } else {
		        output_ii <- output_ii[1:length(output_ii)-1]
		      }
		      haterror <- extract_error(output_ii, error)
		    }
		  }
		  # Calculate expected power
		  expected_powers[ii] <- ifelse(
		    length(output_ii) > 0,
		    sum(sapply(output_ii, function(x) {x$pep * x$weight})),
		    0
		  )
		  all_outputs[[ii]] <- output_ii
		}
		# Pick output with largest expected power
		output <- all_outputs[[which.max(expected_powers)]]
	}
	return(output)
}

#' @keywords internal
#' This is only for use within binarize_selections.
extract_error <- function(discoveries, error) {
  if (length(discoveries) == 0) {
    return(0)
  }
  peps <- sapply(discoveries, function(x) {x$pep})
  if (error == 'fdr') {return(mean(peps))}
  # For FWER, whenever this function is called,
  # we make the PFER approximation.
  return(sum(peps))
}
