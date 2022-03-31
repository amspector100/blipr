check_lp_output <- function(selections, expected) {
  # Check the peps are the same
  peps_out <- sort(as.numeric(lapply(selections, function(x) {x$pep})))
  peps_expected <- sort(as.numeric(lapply(expected, function(x) {x$pep})))
  testthat::expect_equal(peps_out, peps_expected)

  # Check the groups are the same (sorting to avoid spurrious errors)
  groups_out <- sort(as.character(lapply(selections, function(x) {paste0(x$group, collapse="-")})))
  groups_expected <- sort(as.character(lapply(expected, function(x) {paste0(x$group, collapse="-")})))
  testthat::expect_equal(groups_out, groups_expected)
}

check_disjoint <- function(selections) {
  # Iteratively test that the groups are disjoint
  groups <- lapply(selections, function(x) x$group)
  all_locs <- c()
  for (g in groups) {
    testthat::expect_equal(
      numeric(0),
      intersect(all_locs, g)
    )
    all_locs <- union(all_locs, g)
  }
}

check_error <- function(selections, error, q, info=NULL) {
  # Check that the selections control the error rate
  if (length(selections) == 0) {return(0)}
  peps <- sapply(selections, function(x) {x$pep})
  if (error == 'fdr') {achieved_error <- mean(peps)}
  else if (error == 'pfer' | error == "fwer") {achieved_error <- sum(peps)}
  else if (error == 'local_fdr') {achieved_error <- max(peps)}
  else {stop("Unrecognized error rate")}
  if (achieved_error > q) {
    print("About to fail a check_error test. Printing output:")
    cat("error=", error, "q=", q, "selections=\n")
    print(selections)
    cat("info=")
    print(info)
  }
  testthat::expect_true(
    achieved_error <= q
  )
}

test_that('blip (fdr) makes correct detections in tricky examples', {
   # Cand groups where FDR solution is not an optimal PFER solution
   cand_groups <- list(
     list(group=c(1), pep=0.05, weight=1),
     list(group=c(2), pep=0.1, weight=1),
     list(group=c(3), pep=0.05, weight=1/100),
     list(group=c(4), pep=0.05, weight=1/200)
   )
   detections <- BLiP(
     cand_groups=cand_groups,
     error='fdr',
     weight_fn='prespecified',
     q=0.05
   )
   check_lp_output(
     detections,
     cand_groups[c(1, 3, 4)]
   )
   # cand groups where higher PFER leads to lower power
   cand_groups <- list(
     list(group=c(1,2), pep=0.001, weight=0.05),
     list(group=c(2,3), pep=0.05, weight=1),
     list(group=c(3), pep=0.0001, weight=0.05)
   )
   detections <- BLiP(
     cand_groups=cand_groups,
     error='fdr',
     weight_fn='prespecified',
     q=0.05
   )
   check_lp_output(
     detections,
     cand_groups[c(2)]
   )
 })

test_that('blip (fdr) can backtrack correctly when necessary', {
  # Make sure BLiP returns a valid solution
  q <- 0.1
  cand_groups <- list(
    list(group=c(1,2), pep=0.0, weight=1.01),
    list(group=c(2,3), pep=0.0, weight=1),
    list(group=c(1,3), pep=0.0, weight=1),
    list(group=c(4), pep=0.25, weight=1)
  )
  detections <- BLiP(
    cand_groups=cand_groups,
    error='fdr',
    weight_fn='prespecified',
    q=q,
  )
  check_lp_output(
    detections,
    cand_groups[c(1)]
  )

  # Make sure BLiP returns the optimal solution
  cand_groups <- list(
    list(group=c(1,2), pep=0.0, weight=1.01),
    list(group=c(2,3), pep=0.0, weight=1),
    list(group=c(1,3), pep=0.0, weight=1),
    list(group=c(4), pep=0.1, weight=1),
    list(group=c(5), pep=0.21, weight=1)
  )
  detections <- BLiP(
    cand_groups=cand_groups,
    error='fdr',
    weight_fn='prespecified',
    q=q,
  )
  check_lp_output(
    detections,
    cand_groups[c(1,4)]
  )
})



 test_that('binarize_selections makes correct choices in simple example', {
   # Fake problem data
   cand_groups <- list(
     list(group=c(1,2,3), pep=0.01),
     list(group=c(1,2), pep=0.05),
     list(group=c(1), pep=0.1),
     list(group=c(2), pep=0.08),
     list(group=c(3), pep=0.5)
   )

   # Derandomized version for FDR, FWER, PFER
   selections <- BLiP(cand_groups=cand_groups, error="fdr", q=0.1)
   check_lp_output(
     selections,
     cand_groups[c(3, 4)]
   )
   selections <- BLiP(cand_groups=cand_groups, error="fwer", q=0.1, search_method='none')
   check_lp_output(
     selections,
     cand_groups[c(4)]
   )
   selections <- BLiP(cand_groups=cand_groups, error="PFER", q=0.679)
   check_lp_output(
     selections,
     cand_groups[c(3, 4)]
   )
})

test_that('binarize_selections and BLiP deterministically control the error rate', {
  # Fake problem data: two very simple settings
  cand_groups1 <- list(
    list(group=c(1,2,3), pep=0.05),
    list(group=c(2,3), pep=0.15)
  )
  cand_groups2 <- list(
    list(group=c(1), pep=0.5, weight=0.5),
    list(group=c(2), pep=0.15, weight=0.5)
  )
  all_cand_groups <- list(cand_groups1, cand_groups2)
  # ten additional more complex settings
  set.seed(123)
  num_cand_groups <- 30
  p <- 200
  m <- 5
  nsetting <- 10
  for (j in 1:nsetting) {
    new_cand_groups <- list()
    for (i in 1:num_cand_groups) {
      new_cand_groups[[i]] <- list(
        group=sample(1:p, size=m, replace=F),
        pep=stats::runif(1)
      )
    }
    all_cand_groups[[2 + j]] = new_cand_groups
  }
  # Check error control and disjointness
  reps <- 3 # number of repitions for randomized solution
  for (error in c("fdr")) {#c("local_fdr", "fwer", "pfer", "fdr")) {
    for (j in 1:length(all_cand_groups)) {
      q <- 0.2 * runif(1)
      # run deterministic BLiP
      print(all_cand_groups[[j]])
      detections <- BLiP(
        cand_groups=all_cand_groups[[j]], error=error, q=q, deterministic=T
      )
      check_disjoint(detections)
      check_error(detections, error=error, q=q, info=cg_input)
      # Run non-deterministic version
      print(all_cand_groups[[j]])
      for (r in reps) {
        detections <- BLiP(
          cand_groups=all_cand_groups[[j]], error=error, q=q, deterministic=F
        )
        check_disjoint(detections)
        check_error(detections, error=error, q=q, info=cg_input)
      }
    }
  }

  # Repeat with specific input sprobs
  cand_groups <- list(
    list(group=c(1), pep=0.2, weight=0.5, sprob=0.5),
    list(group=c(1,2), pep=0.0, weight=0.25, sprob=0.5),
    list(group=c(3), pep=0.2, weight=1, sprob=0.5),
    list(group=c(3,4), pep=0.0, weight=0.5, sprob=0.5)
  )
  # Check randomized result
  reps <- 10
  q <- 0.2
  for (error in c("fdr", "pfer", "local_fdr")) {
    # Check deterministic result
    signals <- binarize_selections(
      cand_groups, q=q, error=error, deterministic=T
    )
    check_disjoint(signals)
    check_error(signals, error=error, q=q)
    # if (error == 'pfer') {
    #   check_lp_output(detections, cand_groups[c(2,3)])
    # } else if (error == 'fdr') {
    #   check_lp_output(detections, cand_groups[c(2,4)])
    # }
    # Repeat many times for randomized solution
    for (j in 1:reps) {
      signals <- binarize_selections(
        cand_groups, q=q, error=error, deterministic=F
      )
      check_disjoint(signals)
      check_error(signals, error=error, q=q)
    }
  }
})
