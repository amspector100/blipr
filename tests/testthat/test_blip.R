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

# test_that('blip (fdr) makes correct detections in tricky examples', {
#    # Cand groups where FDR solution is not an optimal PFER solution
#    cand_groups <- list(
#      list(group=c(1), pep=0.05, weight=1),
#      list(group=c(2), pep=0.1, weight=1),
#      list(group=c(3), pep=0.05, weight=1/100),
#      list(group=c(4), pep=0.05, weight=1/200)
#    )
#    detections <- BLiP(
#      cand_groups=cand_groups,
#      error='fdr',
#      weight_fn='prespecified',
#      q=0.05
#    )
#    check_lp_output(
#      detections,
#      cand_groups[c(1, 3, 4)]
#    )
#    # cand groups where higher PFER leads to lower power
#    cand_groups <- list(
#      list(group=c(1,2), pep=0.001, weight=0.05),
#      list(group=c(2,3), pep=0.05, weight=1),
#      list(group=c(3), pep=0.0001, weight=0.05)
#    )
#    detections <- BLiP(
#      cand_groups=cand_groups,
#      error='fdr',
#      weight_fn='prespecified',
#      q=0.05
#    )
#    check_lp_output(
#      detections,
#      cand_groups[c(2)]
#    )
#  })

test_that('blip (fdr) can backtrack correctly when necessary', {
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
})



#  test_that('binarize_selections makes correct choices in simple example', {
#    # Fake problem data
#    cand_groups <- list(
#      list(group=c(1,2,3), pep=0.01),
#      list(group=c(1,2), pep=0.05),
#      list(group=c(1), pep=0.1),
#      list(group=c(2), pep=0.08),
#      list(group=c(3), pep=0.5)
#    )
#
#    # Derandomized version for FDR, FWER, PFER
#    selections <- BLiP(cand_groups=cand_groups, error="fdr", q=0.1)
#    check_lp_output(
#      selections,
#      cand_groups[c(3, 4)]
#    )
#    selections <- BLiP(cand_groups=cand_groups, error="fwer", q=0.1, search_method='none')
#    check_lp_output(
#      selections,
#      cand_groups[c(4)]
#    )
#    selections <- BLiP(cand_groups=cand_groups, error="PFER", q=0.679)
#    check_lp_output(
#      selections,
#      cand_groups[c(3, 4)]
#    )
#
#    # Check randomization for pfer
#    set.seed(123)
#    selections <- BLiP(cand_groups=cand_groups, error="PFER", q=0.679, deterministic=F)
#    check_lp_output(
#      selections,
#      cand_groups[c(3, 4, 5)]
#    )
# })
#
# test_that('binarize_selections correctly rounds/samples in explicit tradeoff', {
#   # Fake problem data
#   cand_groups <- list(
#     list(group=c(1,2,3), pep=0.05),
#     list(group=c(2,3), pep=0.15)
#   )
#
#   # Randomized FDR is approximately 50-50 between the groups
#   # make sure it chooses each group at least once in 5 runs
#   set.seed(110)
#   group1_flag <- F
#   group2_flag <- F
#   for (j in 1:5) {
#     selections <- BLiP(cand_groups=cand_groups, error="fdr", q=0.1, deterministic=F)
#     testthat::expect_equal(length(selections), 1)
#     if (length(selections[[1]]$group) == 3) {
#       group1_flag <- T
#     } else {
#       group2_flag <- T
#     }
#   }
#   testthat::expect_true(group1_flag & group2_flag)
#
#   # Derandomized for FDR makes correct decision
#   selections <- BLiP(cand_groups=cand_groups, error="fdr", q=0.1)
#   check_lp_output(
#     selections,
#     cand_groups[1]
#   )
#
#   # Repeat for case with randomized pairs
#   cand_groups <- list(
#     list(group=c(1), pep=0.2, weight=0.5, sprob=0.5),
#     list(group=c(1,2), pep=0.0, weight=0.25, sprob=0.5),
#     list(group=c(3), pep=0.2, weight=1, sprob=0.5),
#     list(group=c(3,4), pep=0.0, weight=0.5, sprob=0.5)
#   )
#   # Check deterministic result
#   detections <- binarize_selections(
#     cand_groups, q=0.2, error='pfer', deterministic=T
#   )
#   check_lp_output(
#     detections,
#     cand_groups[c(2,3)]
#   )
#   # Check sampled result
#   reps <- 128
#   powers <- rep(0, reps)
#   for (j in 1:reps) {
#     signals <- binarize_selections(
#       cand_groups, q=0.2, error='pfer', deterministic=F
#     )
#     check_disjoint(signals)
#     powers[j] <- ifelse(length(signals) == 0, 0, sum(sapply(signals, function(x) x$weight)))
#   }
#   expected <- sum(sapply(cand_groups, function(x) x$weight)) / 2
#   testthat::expect_true(
#     abs(mean(powers) - expected) < 4*(sd(powers) / sqrt(reps))
#   )
#
# })
