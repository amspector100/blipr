check_unique <- function(cand_groups) {
  groups <- lapply(
    cand_groups,
    function(x) {x$group}
  )
  testthat::expect_equal(
    length(groups), length(unique(groups))
  )
}

g_in_groups <- function(groups, g) {
  return(any(sapply(groups, function(x) setequal(x, g))))
}

calc_expected_pep <- function(x, y, radius, locs, shape) {
  N <- dim(locs)[1]
  flags <- ! is.na(locs)
  centers <- matrix(c(x, y), 1, 2)
  pip <- 0
  for (j in 1:N) {
    active <- which(apply(
      ! is.nan(abind::adrop(locs[j,,,drop=F], drop=1)),
      1, all
    ))
    if (length(active) > 0) {
      dists <- broadcast_diffs(
        centers,
        abind::adrop(locs[j,active,,drop=F], drop=1)
      ) # 1 x k x 2
      if (shape == 'circle') {
        dists <- sqrt(apply(dists**2, c(2), sum)) # k-length
      } else {
        dists <- apply(abs(dists), c(2), max) # k-length
      }
      if (min(dists) <= radius) {
        pip <- pip + 1 / N
      }
    }
  }
  expected_pep <- 1 - pip
  return(expected_pep)
}

check_cand_groups_correct <- function(cand_groups, shape) {
  # check that overlap calculations are correct
  for (j1 in 2:length(cand_groups)) {
    cg1 <- cand_groups[[j1]]
    for (j2 in 1:(j1-1)) {
      cg2 <- cand_groups[[j2]]
      g1 <- cg1$group
      g2 <- cg2$group
      r1 <- cg1$radius
      r2 <- cg2$radius
      if (shape == 'circle') {
        dist <- sqrt(sum((cg1$center-cg2$center)**2))
      } else {
        dist <- max(abs(cg1$center - cg2$center))
      }
      if (length(intersect(g1, g2)) > 0) {
        testthat::expect_true(
          dist < r1 + r2
        )
      } else {
        if (dist < r1 + r2) {
          cat("About to fail test: note that ")
          cat("g1=", g1, ' ')
          cat("g2=", g2, ' ')
          # print(cg1$group)
          # cat("cg2=")
          # print(cg2$group)
          cat("j1=", j1, "j2=", j2, "\n")
          cat("center1=", cg1$center, "center2=", cg2$center, "\n")
          cat("dist =", dist, "shape =", shape, ' ')
          cat("r1 =", r1, "r2 =", r2, '\n')
        }
        testthat::expect_true(
          dist >= r1 + r2
        )
      }
    }
  }
}

test_that('key functions work properly', {
  # Generate fake data
  for (d in c(1,2,3)) {
    for (k in c(1,2,5)) {
      centers_orig <- round(array(runif(d*k), c(1,k,d)), 3)
      radius <- round(stats::runif(1), 3)
      # make keys
      keys <- centerrad2key(centers_orig, radius)
      # Recover centers/radii
      output <- key2centerrad(keys)
      # Check the results are equal
      testthat::expect_equal(
        output$centers[,,drop=T],
        centers_orig[1,,,drop=T]
      )
      testthat::expect_equal(
        min(output$radii),
        radius
      )
    }
  }
})

test_that('broadcast_diffs works properly', {
  # Generate fake data
  k1 <- 10
  k2 <- 25
  d <- 5
  M1 <- matrix(rnorm(k1*d), k1, d)
  M2 <- matrix(rnorm(k2*d), k2, d)
  diffs <- broadcast_diffs(M1, M2)
  for (j1 in 1:k1) {
    for (j2 in 1:k2) {
      testthat::expect_equal(
        diffs[j1,j2,],
        M1[j1,] - M2[j2,]
      )
    }
  }
})

test_that('test square groups cts', {
  N <- 3
  k <- 2
  d <- 2
  locs <- array(
    c(0.5502, 0.549, 0.553,
      0.300001, 0.305, NaN,
      0.4502, 0.451, 0.456,
      0.200001, 0.201, NaN),
    c(N, k, 2)
  )

  # Run with one manual center added
  xc1 <- 0.55012
  yc1 <- 0.45012
  ecents <- array(c(xc1, yc1), c(1, 2))
  peps <- lattice_peps(
    locs,
    grid_sizes=c(10,100,1000),
    extra_centers=ecents,
    max_pep=1,
    shape='square'
  )
  # Check that PEPs are right
  centers <- array(
    c(0.55, 0.555, 0.3005, xc1,
      0.45, 0.455, 0.2005, yc1),
    c(1, 4, 2)
  )
  radii <-c(1/20, 1/200, 1/2000, 1/20)
  keys <- centerrad2key(centers, radii)
  expected <- c(0, 1/3, 2/3, 0)
  for (j in 1:length(keys)) {
    testthat::expect_equal(
      peps[[keys[[j]]]],
      expected[j]
    )
  }
})

test_that('test circular groups cts', {
  N <- 6
  k <- 2
  d <- 2
  locs <- array(
    c(0.13054, 0.12958, 0.46001, NaN, NaN, NaN,
      NaN, 0.46009, 0.1302, 0.95, NaN, NaN,
      0.410234, 0.406639, 0.45119, NaN, NaN, NaN,
      NaN, 0.459, 0.4126, 0.899, NaN, NaN),
    c(N, k, 2)
  )

  # Run with one manual center added
  peps <- lattice_peps(
    locs,
    grid_sizes=c(10,100,1000),
    max_pep=1,
    shape='circle'
  )
  # Check that PEPs are right
  r1 <- sqrt(2) / 20; r2 <- sqrt(2) / 200
  xs <- c(0.45, 0.135, 0.135, 0.95, 0.95, 0.455)
  ys <- c(0.45, 0.405, 0.415, 0.95, 0.85, 0.465)
  centers <- array(
    c(xs, ys),
    c(1, 6, 2)
  )
  radii <- c(r1, r2, r2, r1, r1, r2)
  for (j in 1:length(radii)) {
    expected_pep <- calc_expected_pep(
      locs=locs,
      x=xs[j],
      y=ys[j],
      radius=radii[j],
      shape='circle'
    )
    centers <- matrix(c(xs[j], ys[j]), 1, 2)
    key <- centerrad2key(
      array(centers, c(1, 1, 2)), radii[j]
    )[[1]]
    if (expected_pep != 1) {
      testthat::expect_equal(
        peps[[key]],
        expected_pep
      )
    }
  }
  pepkeys <- names(peps)
  splits <- sapply(
    pepkeys, function(x) {base::strsplit(x, split=',')[[1]]}
  )
  d <- dim(splits)[1] - 1
  # Recreate matrix of centers / vector of radii
  pepradii <- as.numeric(splits[(d+1),])
  #print(pepkeys[order(pepradii)])
})

test_that('test lattice_peps_to_cand_groups', {
  # Helper function
  create_fpeps <- function(keys, peps) {
    out <- list()
    for (j in 1:length(keys)) {out[[keys[j]]] <- peps[j]}
    return(out)
  }
  # Notes:
  # (1) The first two overlap if they are rectangles, but not if
  # they are circles.
  # (2) The second two do overlap
  # (3) The last two do not overlap

  # (center_x, center_y, key) --> pep
  keys1 <- c(
    paste(0.12, 0.15, 0.01, sep=','),
    paste(0.13, 0.16, 0.001, sep=','),
    paste(0.25, 0.23, 0.011, sep=','),
    paste(0.25, 0.23, 0.053, sep=','),
    paste(0.421, 0.493, 0.01414, sep=','),
    paste(0.419, 0.522, 0.01, sep=','),
    paste(0.419, 0.522, 0.05, sep=','),
    paste(0.419, 0.522, 0.10, sep=',')
  )
  peps1 <- c(
    0.113, 0.514, 0.153, 0.0513,
    0.0999, 0.01, 0.001, 0.0
  )

  # random example
  set.seed(123)
  ncg <- 50
  peps2 <- rep(0, ncg)
  keys2 <- rep("0", ncg)
  for (j in 1:ncg) {
    keys2[j] <- paste(round(stats::runif(3),4), collapse=",")
    peps2[j] <- stats::runif(1)
  }

  # Examples where we know the optimal solution
  # to the ECC problem
  keys3 <- c(
    paste(0.12, 0.15, 0.01, sep=','),
    paste(0.12, 0.15, 0.02, sep=','),
    paste(0.12, 0.15, 0.03, sep=','),
    paste(0.12, 0.15, 0.04, sep=','),
    paste(0.12, 0.15, 0.05, sep=','),
    paste(0.12, 0.15, 0.06, sep=',')
  )
  peps3 <- rep(0, length(keys3))

  keys4 <- c(
    paste(0.12, 0.15, 0.01, sep=','),
    paste(0.12, 0.15, 0.02, sep=','),
    paste(0.12, 0.15, 0.05, sep=','),
    paste(0.20, 0.20, 0.05, sep=','),
    paste(0.20, 0.20, 0.001, sep=','),
    paste(0.20, 0.20, 0.002, sep=',')
  )
  peps4 <- rep(0, length(keys4))

  all_peps <- list(peps1, peps2, peps3, peps4)
  all_keys <- list(keys1, keys2, keys3, keys4)
  opt_nloc <- c(NaN, NaN, 1, 3)
  # Test both examples
  for (test in 1:length(all_peps)) {
    for (shape in c("circle", "square")) {
      out <- lattice_peps_to_cand_groups(
        create_fpeps(all_keys[[test]], all_peps[[test]]),
        verbose=F,
        shape=shape
      )
      # In simple example, there should be only one merged connected component
      cand_groups <- out$all_cand_groups
      testthat::expect_true(
        length(cand_groups)==1
      )
      cand_groups <- cand_groups[[1]] # extract first component

      # Check correctness
      check_cand_groups_correct(cand_groups, shape)

      # Check number of locations
      if (! is.nan(opt_nloc[test])) {
        nlocs <- length(Reduce(union, lapply(cand_groups, function(x) x$group)))
        testthat::expect_true(
          nlocs <= opt_nloc[test]
        )
      }
    }
  }
})

test_that('test BLiP_cts', {
  # Create locs
  set.seed(123)
  N <- 30
  k <- 3
  for (d in c(1:3)) {
    mu <- matrix(runif(k*d), k, d)
    sigma <- runif(k) / 10
    locs <- array(runif(N*k*d), c(N, k, d))
    for (j in 1:k) {
      locs[,j,] <- locs[,j,] * sigma[j]
    }
    for (j in 1:N) {
      locs[j,,] <- locs[j,,] + mu
    }
    if (d == 2) {
      shapes <- c('square','circle')
    } else {
      shapes <- c("square")
    }
    for (shape in shapes) {
      rej <- BLiP_cts(
        locs=locs,
        extra_centers=mu,
        max_pep=0.25,
        shape=shape,
        min_blip_size=50001,
        verbose=F
      )
      # check disjointness is correct
      check_cand_groups_correct(rej, shape=shape)
      check_unique(rej)
      # Check peps are correct
      if (d == 2) {
        for (rj in 1:length(rej)) {
          expected_pep <- calc_expected_pep(
            x=rej[[rj]]$center[1],
            y=rej[[rj]]$center[2],
            radius=rej[[rj]]$radius,
            locs=locs,
            shape=shape
          )
          testthat::expect_equal(
            expected_pep,
            rej[[rj]]$pep
          )
        }
      }
    }
  }
})
