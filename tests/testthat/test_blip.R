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


test_that('binarize.selections makes correct choices in simple example', {
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

  # Check randomization for pfer
  set.seed(123)
  selections <- BLiP(cand_groups=cand_groups, error="PFER", q=0.679, deterministic=F)
  check_lp_output(
    selections,
    cand_groups[c(3, 4, 5)]
  )
})

test_that('binarize.selections correctly rounds/samples in explicit tradeoff', {
  # Fake problem data
  cand_groups <- list(
    list(group=c(1,2,3), pep=0.05),
    list(group=c(2,3), pep=0.15)
  )

  # Randomized FDR is approximately 50-50 between the groups
  # make sure it chooses each group at least once in 5 runs
  set.seed(110)
  group1_flag <- F
  group2_flag <- F
  for (j in 1:5) {
    selections <- BLiP(cand_groups=cand_groups, error="fdr", q=0.1, deterministic=F)
    testthat::expect_equal(length(selections), 1)
    if (length(selections[[1]]$group) == 3) {
      group1_flag <- T
    } else {
      group2_flag <- T
    }
  }
  testthat::expect_true(group1_flag & group2_flag)

  # Derandomized for FDR makes correct decision
  selections <- BLiP(cand_groups=cand_groups, error="fdr", q=0.1)
  check_lp_output(
    selections,
    cand_groups[1]
  )
})
