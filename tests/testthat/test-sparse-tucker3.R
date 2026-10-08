test_that("zero penalty preserves an already metric-normalized basis", {
  covariance <- matrix(c(2, 0.2, 0.2, 1.5), 2, 2)
  raw <- matrix(c(1, 0, 0, 1), 2, 2)
  basis <- threeway:::covariance_normalize_basis(raw, covariance)

  result <- scr_sparse_metric_basis(
    basis = basis,
    covariance = covariance,
    penalty = 0
  )

  expect_equal(result$basis, basis, tolerance = 1e-10)
  expect_equal(result$active, 1:2)
  expect_equal(
    result$metric_gram,
    diag(2),
    tolerance = 1e-9
  )
})


test_that("weak loading rows are thresholded exactly to zero", {
  covariance <- diag(3)
  basis <- matrix(
    c(
      1, 0,
      0.05, 0,
      0, 1
    ),
    nrow = 3,
    byrow = TRUE
  )

  result <- scr_sparse_metric_basis(
    basis = basis,
    covariance = covariance,
    penalty = 0.1
  )

  expect_equal(result$basis[2, ], c(0, 0))
  expect_false(2L %in% result$active)
  expect_equal(
    result$metric_gram,
    diag(2),
    tolerance = 1e-9
  )
})


test_that("metric normalization preserves thresholded zero rows", {
  covariance <- matrix(
    c(
      2, 0.2, 0,
      0.2, 1.5, 0.1,
      0, 0.1, 1
    ),
    3,
    3
  )
  basis <- threeway:::covariance_normalize_basis(
    matrix(
      c(
        1, 0,
        0.02, 0.01,
        0, 1
      ),
      nrow = 3,
      byrow = TRUE
    ),
    covariance
  )

  penalty <- min(
    sqrt(sum(basis[2, ]^2)) * 1.1,
    min(
      sqrt(sum(basis[1, ]^2)),
      sqrt(sum(basis[3, ]^2))
    ) * 0.9
  )

  result <- scr_sparse_metric_basis(
    basis = basis,
    covariance = covariance,
    penalty = penalty
  )

  expect_equal(
    result$basis[2, ],
    c(0, 0),
    tolerance = 1e-12
  )
  expect_equal(
    crossprod(
      result$basis,
      solve(covariance, result$basis)
    ),
    diag(2),
    tolerance = 1e-8
  )
})


test_that("larger penalties produce nested row supports", {
  covariance <- diag(4)
  basis <- matrix(
    c(
      1, 0,
      0.5, 0,
      0, 1,
      0, 0.25
    ),
    nrow = 4,
    byrow = TRUE
  )

  low <- scr_sparse_metric_basis(
    basis,
    covariance,
    penalty = 0.1
  )
  high <- scr_sparse_metric_basis(
    basis,
    covariance,
    penalty = 0.3
  )

  expect_true(
    all(high$active %in% low$active)
  )
})


test_that("zero-penalty sparse Tucker3 projection recovers exact means", {
  probabilities <- c(0.2, 0.3, 0.5)
  group_mass <- 100 * probabilities

  raw_A <- matrix(c(1, -1, 0), ncol = 1)
  A <- scr_centroid_basis(raw_A, probabilities)
  A <- A / sqrt(
    as.numeric(crossprod(A, group_mass * A))
  )

  B <- matrix(c(1, 0, 0), ncol = 1)
  C <- matrix(c(1, 0), ncol = 1)
  core <- array(2, dim = c(1, 1, 1))
  grand <- matrix(0, nrow = 3, ncol = 2)

  means <- scr_tucker3_means(
    grand,
    A,
    B,
    C,
    core,
    probabilities
  )

  result <- scr_sparse_tucker3_projection(
    group_centroids = means,
    group_mass = group_mass,
    centroid_basis = A,
    variable_basis = B,
    occasion_basis = C,
    variable_covariance = diag(3),
    occasion_covariance = diag(2),
    variable_penalty = 0,
    occasion_penalty = 0
  )

  expect_equal(result$means, means, tolerance = 1e-9)
  expect_equal(result$core, core, tolerance = 1e-9)
  expect_lt(result$weighted_residual, 1e-9)
})


test_that("sparse Tucker3 path retains infeasible penalties explicitly", {
  probabilities <- c(0.4, 0.6)
  group_mass <- c(40, 60)

  raw_A <- matrix(c(1, -1), ncol = 1)
  A <- scr_centroid_basis(raw_A, probabilities)
  A <- A / sqrt(
    as.numeric(crossprod(A, group_mass * A))
  )

  B <- matrix(c(1, 0.1, 0), ncol = 1)
  C <- matrix(c(1, 0.1), ncol = 1)
  core <- array(1.5, dim = c(1, 1, 1))
  grand <- matrix(0, nrow = 3, ncol = 2)

  means <- scr_tucker3_means(
    grand,
    A,
    B,
    C,
    core,
    probabilities
  )

  path <- scr_sparse_tucker3_path(
    group_centroids = means,
    group_mass = group_mass,
    centroid_basis = A,
    variable_basis = B,
    occasion_basis = C,
    variable_covariance = diag(3),
    occasion_covariance = diag(2),
    variable_penalties = c(0, 2),
    occasion_penalties = c(0, 2)
  )

  expect_s3_class(path, "scr_sparse_tucker3_path")
  expect_equal(nrow(path$summary), 4L)
  expect_true(any(path$summary$feasible))
  expect_true(any(!path$summary$feasible))
})


test_that("rank collapse is rejected clearly", {
  expect_error(
    scr_sparse_metric_basis(
      basis = diag(2),
      covariance = diag(2),
      penalty = 2
    ),
    "fewer active rows"
  )
})
