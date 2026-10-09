test_that("joint Tucker3 group selection evaluates every G", {
  set.seed(2601)

  n <- 60L
  X <- matrix(rnorm(n * 4), nrow = n)

  result <- select_scr_tucker3_groups(
    X = X,
    groups = c(2, 3),
    variable_covariance = diag(2),
    occasion_covariance = diag(2),
    n_starts = 1,
    centroid_ranks = 1,
    variable_ranks = 1,
    occasion_ranks = 1,
    max_iter = 5,
    inner_max_iter = 5,
    seed = 11
  )

  expect_s3_class(result, "scr_tucker3_group_selection")
  expect_equal(result$comparison$G, c(2L, 3L))
  expect_equal(nrow(result$comparison), 2L)
  expect_equal(sum(result$comparison$selected), 1L)
  expect_true(result$selected_structure[["G"]] %in% c(2L, 3L))
  expect_true(all(result$comparison$P == 1L))
  expect_true(all(result$comparison$Q == 1L))
  expect_true(all(result$comparison$R == 1L))
})


test_that("joint group selection is reproducible under a fixed seed", {
  set.seed(2602)
  X <- matrix(rnorm(48 * 4), nrow = 48)

  first <- select_scr_tucker3_groups(
    X = X,
    groups = c(2, 3),
    variable_covariance = diag(2),
    occasion_covariance = diag(2),
    n_starts = 2,
    centroid_ranks = 1,
    variable_ranks = 1,
    occasion_ranks = 1,
    max_iter = 4,
    inner_max_iter = 4,
    seed = 123
  )

  second <- select_scr_tucker3_groups(
    X = X,
    groups = c(2, 3),
    variable_covariance = diag(2),
    occasion_covariance = diag(2),
    n_starts = 2,
    centroid_ranks = 1,
    variable_ranks = 1,
    occasion_ranks = 1,
    max_iter = 4,
    inner_max_iter = 4,
    seed = 123
  )

  expect_equal(first$comparison, second$comparison)
  expect_equal(first$selected_structure, second$selected_structure)
  expect_equal(first$criterion_value, second$criterion_value)
})


test_that("outer selector supports ICL", {
  comparison <- data.frame(
    G = c(2L, 3L),
    P = c(1L, 1L),
    Q = c(1L, 1L),
    R = c(1L, 1L),
    parameters = c(20L, 24L),
    bic = c(100, 101),
    icl = c(98, 95)
  )

  expect_equal(
    scr3way:::select_best_group_candidate(
      comparison,
      criterion = "BIC"
    ),
    2L
  )
  expect_equal(
    scr3way:::select_best_group_candidate(
      comparison,
      criterion = "ICL"
    ),
    1L
  )
})


test_that("outer criterion ties prefer the lower-complexity structure", {
  comparison <- data.frame(
    G = c(3L, 2L, 4L),
    P = c(1L, 1L, 1L),
    Q = c(1L, 1L, 1L),
    R = c(1L, 1L, 1L),
    parameters = c(22L, 20L, 20L),
    bic = c(100, 100 + 5e-9, 100 + 4e-9),
    icl = c(90, 90, 90)
  )

  expect_equal(
    scr3way:::select_best_group_candidate(
      comparison,
      criterion = "BIC",
      criterion_tolerance = 1e-8
    ),
    2L
  )
})


test_that("custom group membership initializers are supported", {
  set.seed(2603)
  labels <- rep(1:2, each = 20)
  means <- rbind(c(-2, 0, -1, 0), c(2, 0, 1, 0))
  X <- means[labels, , drop = FALSE] +
    matrix(rnorm(40 * 4, sd = 0.35), nrow = 40)

  initializer <- function(X, groups, start, seed) {
    ordered <- rank(X[, 1L], ties.method = "first")
    labels <- pmin(
      groups,
      ceiling(ordered * groups / nrow(X))
    )
    membership <- matrix(1, nrow = nrow(X), ncol = groups)
    membership[cbind(seq_len(nrow(X)), labels)] <- 3 + start
    membership
  }

  result <- select_scr_tucker3_groups(
    X = X,
    groups = 2,
    variable_covariance = diag(2),
    occasion_covariance = diag(2),
    n_starts = 1,
    centroid_ranks = 1,
    variable_ranks = 1,
    occasion_ranks = 1,
    membership_initializer = initializer,
    max_iter = 3,
    inner_max_iter = 3,
    seed = 5
  )

  expect_equal(result$comparison$G, 2L)
  expect_equal(result$comparison$start, 1L)
})


test_that("group selector rejects impossible group and centroid settings", {
  X <- matrix(rnorm(40), nrow = 10, ncol = 4)

  expect_error(
    select_scr_tucker3_groups(
      X = X,
      groups = 1,
      variable_covariance = diag(2),
      occasion_covariance = diag(2)
    ),
    "greater than or equal to two"
  )

  expect_error(
    select_scr_tucker3_groups(
      X = X,
      groups = 2,
      variable_covariance = diag(2),
      occasion_covariance = diag(2),
      centroid_ranks = 2,
      variable_ranks = 1,
      occasion_ranks = 1,
      max_iter = 2,
      inner_max_iter = 2
    ),
    "no admissible centroid_ranks remain"
  )
})
