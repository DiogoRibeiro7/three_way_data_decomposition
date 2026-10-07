test_that("nearest Kronecker approximation recovers exact separability", {
  SV <- matrix(c(2, 0.4, 0.4, 1.2), 2, 2)
  SO <- matrix(c(1.5, 0.2, 0.2, 0.9), 2, 2)
  covariance <- kronecker(SO, SV)

  result <- scr_nearest_kronecker_covariance(
    covariance = covariance,
    variables = 2,
    occasions = 2
  )

  expect_lt(result$relative_error, 1e-10)
  expect_equal(result$rank_one_fraction, 1, tolerance = 1e-10)
  expect_equal(
    mean(diag(result$variable_covariance)),
    1,
    tolerance = 1e-10
  )
  expect_equal(
    result$covariance,
    covariance,
    tolerance = 1e-9
  )
})


test_that("nonseparable perturbation produces positive separability gap", {
  SV <- matrix(c(1.5, 0.2, 0.2, 1), 2, 2)
  SO <- matrix(c(1.2, 0.1, 0.1, 0.8), 2, 2)
  covariance <- kronecker(SO, SV)

  perturbation <- matrix(
    c(
      0.15, 0, 0, 0.05,
      0, 0.05, 0.02, 0,
      0, 0.02, 0.08, 0,
      0.05, 0, 0, 0.12
    ),
    nrow = 4,
    byrow = TRUE
  )
  covariance <- covariance + perturbation

  result <- scr_nearest_kronecker_covariance(
    covariance = covariance,
    variables = 2,
    occasions = 2
  )

  expect_gt(result$relative_error, 0)
  expect_gte(result$rank_one_fraction, 0)
  expect_lte(result$rank_one_fraction, 1)
  expect_silent(chol(result$variable_covariance))
  expect_silent(chol(result$occasion_covariance))
})


test_that("zero nugget reproduces the exact Kronecker covariance", {
  SV <- matrix(c(2, 0.25, 0.25, 1), 2, 2)
  SO <- matrix(c(1.4, 0.1, 0.1, 0.8), 2, 2)

  covariance <- scr_kronecker_nugget_covariance(
    SV,
    SO,
    nugget = 0
  )

  expect_equal(covariance, kronecker(SO, SV))
})


test_that("positive nugget covariance has consistent spectral diagnostics", {
  SV <- matrix(c(2, 0.25, 0.25, 1), 2, 2)
  SO <- matrix(c(1.4, 0.1, 0.1, 0.8), 2, 2)

  result <- scr_kronecker_nugget_covariance(
    SV,
    SO,
    nugget = 0.3,
    details = TRUE
  )

  expect_silent(chol(result$covariance))
  direct_logdet <- as.numeric(
    determinant(result$covariance, logarithm = TRUE)$modulus
  )
  expect_equal(result$log_determinant, direct_logdet, tolerance = 1e-10)
  expect_gt(result$condition_number, 1)
  expect_true(all(result$eigenvalues > 0))
})


test_that("covariance parameter counts are consistent", {
  separable <- scr_covariance_parameter_count(
    variables = 3,
    occasions = 2,
    model = "separable"
  )
  nugget <- scr_covariance_parameter_count(
    variables = 3,
    occasions = 2,
    model = "nugget"
  )
  unrestricted <- scr_covariance_parameter_count(
    variables = 3,
    occasions = 2,
    model = "unrestricted"
  )

  expect_equal(separable, 8L)
  expect_equal(nugget, 9L)
  expect_equal(unrestricted, 21L)
  expect_lt(separable, nugget)
  expect_lt(nugget, unrestricted)
})


test_that("covariance utilities validate invalid controls", {
  expect_error(
    scr_kronecker_nugget_covariance(
      diag(2),
      diag(2),
      nugget = -0.1
    ),
    "non-negative"
  )

  expect_error(
    scr_nearest_kronecker_covariance(
      covariance = diag(4),
      variables = 3,
      occasions = 2
    ),
    "dimensions"
  )
})
