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

test_that('sequential_groups calculates peps and groups correctly', {
  # Generate fake data
  inclusions = rbind(
    c(0, 1, 0, 0, 1),
    c(1, 1, 0, 1, 1),
    c(1, 0, 0, 0, 1),
    c(0, 1, 0, 1, 1),
    c(1, 0, 0, 1, 1)
  )
  # Create groups
  cand_groups <- sequential_groups(inclusions, max_size=4)
  check_unique(cand_groups)
  expected_length <- 5 + 4 + 3 + 2 - 1 # ignore column with all 0s
  testthat::expect_equal(
    length(cand_groups), expected_length
  )
  # Check that max_peps works
  for (max_pep in c(0.2, 0.4, 0.6, 0.8, 1)) {
    cand_groups <- sequential_groups(
      inclusions, max_size=4, max_pep=max_pep
    )
    max_pep_obs = max(sapply(cand_groups, function(x) {x$pep}))
    testthat::expect_true(
      max_pep_obs < max_pep,
    )
    # Check for duplicates for good measure
    check_unique(cand_groups)
  }
  # Check that prenarrowing works
  q = 0.5
  max_pep = 0.5
  cand_groups <- sequential_groups(
    inclusions, prenarrow=TRUE, q=q, max_pep=max_pep
  )
  check_unique(cand_groups)
  groups = lapply(cand_groups, function(x) {x$group})
  testthat::expect_true(
    g_in_groups(groups, c(1))
  )
  testthat::expect_true(
    g_in_groups(groups, c(1,2)) == 0
  )
  testthat::expect_true(
    g_in_groups(groups, c(1,2,3)) == 0
  )
  # Finally, check that we can successfully eliminate the last redundant feature
  cand_groups = elim_redundant_features(cand_groups)
  nrel = length(Reduce(union, lapply(cand_groups, function(x) {x$blip_group})))
  testthat::expect_equal(
    nrel, 4
  )
})


test_that('hierarchical_groups calculates peps and groups correctly', {
  # Generate fake data
  inclusions = rbind(
    c(0, 1, 0, 0, 1),
    c(1, 1, 0, 1, 1),
    c(1, 0, 0, 0, 1),
    c(0, 1, 0, 1, 1),
    c(1, 0, 0, 1, 1)
  )
  # Create groups
  cand_groups <- hierarchical_groups(inclusions, max_size=4, max_pep=1)
  groups <- lapply(cand_groups, function(x) {x$group})
  expected_groups <- list(
    c(1), c(2), c(4), c(5), c(1,2)
  )
  for (g in expected_groups) {
    testthat::expect_true(
      g_in_groups(groups, g)
    )
  }

  # Check that max_peps, max_size works
  max_size = 4
  for (max_pep in c(0.2, 0.4, 0.6, 0.8, 1)) {
    cand_groups <- hierarchical_groups(
      inclusions, max_size=max_size, max_pep=max_pep
    )
    max_pep_obs = max(sapply(cand_groups, function(x) x$pep))
    testthat::expect_true(
      max_pep_obs < max_pep,
    )
    max_size_obs = max(sapply(cand_groups, function(x) length(x$group)))
    testthat::expect_true(
      max_size_obs <= max_size,
    )
    # Check for duplicates for good measure
    check_unique(cand_groups)
  }

})

