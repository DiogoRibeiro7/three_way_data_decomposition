test_that("Tucker3 mean update recovers an exact low-rank mean tensor", {
  probabilities <- c(0.2, 0.3, 0.5)
  group_mass <- 100 * probabilities

  raw_A <- matrix(c(1, -1, 0), ncol = 1)
  A <- scr_centroid_basis(raw_A, probabilities)
  A <- A / sqrt(
    as.numeric(crossprod(A, group_mass * A))
  )

  B <- qr.Q(qr(matrix(c(1, 0, 1, 1, 1, 0), nrow = 3, ncol = 2)))
  C <- matrix(c(1, 0), ncol = 1)
  core <- array(c(2, -1), dim = c(1, 2, 1))
  grand <- matrix(
    c(1, 2, 3, 4, 5, 6),
    nrow = 3,
    ncol = 2
  )

  means <- scr_tucker3_means(
    grand,
    A,
    B,
    C,
    core,
    probabilities
  )

  update <- scr_tucker3_mean_update(
    group_centroids = means,
    group_mass = group_mass,
    variable_covariance = diag(3),
    occasion_covariance = diag(2),
    centroid_rank = 1,
    variable_rank = 2,
    occasion_rank = 1,
    inner_max_iter = 30,
    inner_tolerance = 1e-10
  )

  expect_equal(update$means, means, tolerance = 1e-7)
  expect_lt(update$weighted_residual, 1e-7)
  expect_equal(
    as.numeric(crossprod(probabilities, update$centroid_basis)),
    0,
    tolerance = 1e-8
  )
  expect_equal(
    crossprod(
      update$centroid_basis,
      group_mass * update$centroid_basis
    ),
    matrix(1),
    tolerance = 1e-8
  )
})

test_that("Tucker3 mean update respects covariance metrics", {
  set.seed(123)

  centroids <- matrix(rnorm(4 * 6), nrow = 4)
  mass <- c(20, 30, 25, 25)
  SV <- matrix(c(2, 0.3, 0.3, 1.5), 2, 2)
  SO <- matrix(c(1.5, 0.2, 0.2, 1.2), 2, 2)

  # Use three variables? Here SV is 2x2 and SO is 2x2, so JK = 4.
  centroids <- centroids[, 1:4, drop = FALSE]

  update <- scr_tucker3_mean_update(
    group_centroids = centroids,
    group_mass = mass,
    variable_covariance = SV,
    occasion_covariance = SO,
    centroid_rank = 2,
    variable_rank = 1,
    occasion_rank = 1,
    inner_max_iter = 20
  )

  expect_equal(
    crossprod(
      update$variable_basis,
      solve(SV, update$variable_basis)
    ),
    matrix(1),
    tolerance = 1e-7
  )
  expect_equal(
    crossprod(
      update$occasion_basis,
      solve(SO, update$occasion_basis)
    ),
    matrix(1),
    tolerance = 1e-7
  )
})

test_that("Tucker3 SCR fit returns a coherent fitted model", {
  set.seed(456)

  n <- 120
  probabilities <- c(0.25, 0.35, 0.40)
  group_mass <- n * probabilities

  raw_A <- matrix(c(1, -1, 0), ncol = 1)
  A <- scr_centroid_basis(raw_A, probabilities)
  A <- A / sqrt(
    as.numeric(crossprod(A, group_mass * A))
  )

  B <- matrix(c(1, 0, 0), ncol = 1)
  C <- matrix(c(1, 0), ncol = 1)
  core <- array(4, dim = c(1, 1, 1))
  grand <- matrix(0, nrow = 3, ncol = 2)

  means <- scr_tucker3_means(
    grand,
    A,
    B,
    C,
    core,
    probabilities
  )
  covariance <- diag(6)

  generated <- generate_gaussian_mixture(
    n = n,
    probabilities = probabilities,
    means = means,
    covariances = replicate(3, covariance, simplify = FALSE)
  )

  U0 <- matrix(runif(n * 3), nrow = n, ncol = 3)
  U0 <- U0 / rowSums(U0)

  fit <- fit_scr_s3_tucker3(
    X = generated$X,
    membership = U0,
    centroid_rank = 1,
    variable_rank = 1,
    occasion_rank = 1,
    variable_covariance = diag(3),
    occasion_covariance = diag(2),
    tolerance = 1e-5,
    max_iter = 30,
    inner_max_iter = 20,
    inner_tolerance = 1e-7
  )

  expect_equal(dim(fit$U), c(n, 3L))
  expect_equal(dim(fit$A), c(3L, 1L))
  expect_equal(dim(fit$TB), c(3L, 1L))
  expect_equal(dim(fit$TC), c(2L, 1L))
  expect_equal(dim(fit$core), c(1L, 1L, 1L))
  expect_equal(dim(fit$M), c(3L, 6L))
  expect_equal(rowSums(fit$U), rep(1, n), tolerance = 1e-10)
  expect_equal(sum(fit$probabilities), 1, tolerance = 1e-10)
  expect_true(is.finite(fit$like))
  expect_true(is.finite(fit$bic))

  expected_bic <- 2 * fit$like -
    log(n) * scr_tucker3_parameter_count(
      groups = 3,
      variables = 3,
      occasions = 2,
      centroid_rank = 1,
      variable_rank = 1,
      occasion_rank = 1
    )
  expect_equal(fit$bic, expected_bic)
})

test_that("Tucker3 fitter rejects impossible centroid ranks", {
  X <- matrix(rnorm(40), nrow = 10, ncol = 4)
  U <- matrix(0.5, nrow = 10, ncol = 2)

  expect_error(
    fit_scr_s3_tucker3(
      X = X,
      membership = U,
      centroid_rank = 2,
      variable_rank = 1,
      occasion_rank = 1,
      variable_covariance = diag(2),
      occasion_covariance = diag(2)
    ),
    "centroid_rank cannot exceed groups - 1"
  )
})
