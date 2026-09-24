test_that("generate_gaussian_mixture returns valid SCR simulation objects", {
  set.seed(123)

  means <- matrix(
    c(
      0, 0,
      3, 3
    ),
    nrow = 2,
    byrow = TRUE
  )
  covariances <- list(diag(2), diag(2))

  result <- generate_gaussian_mixture(
    n = 50,
    probabilities = c(0.4, 0.6),
    means = means,
    covariances = covariances
  )

  expect_equal(dim(result$X), c(50L, 2L))
  expect_equal(dim(result$U), c(50L, 2L))
  expect_length(result$z, 50L)
  expect_true(all(result$z %in% 1:2))
  expect_equal(rowSums(result$U), rep(1, 50), tolerance = 1e-12)
  expect_true(all(result$U >= 0))
  expect_true(all(result$U <= 1))
})

test_that("posterior memberships match the Gaussian reference formula", {
  set.seed(456)

  probabilities <- c(0.25, 0.75)
  means <- matrix(
    c(
      0, 0,
      2, 2
    ),
    nrow = 2,
    byrow = TRUE
  )

  result <- generate_gaussian_mixture(
    n = 20,
    probabilities = probabilities,
    means = means,
    covariances = list(diag(2), diag(2))
  )

  log_kernel <- cbind(
    -0.5 * rowSums(
      sweep(result$X, 2, means[1, ], "-")^2
    ),
    -0.5 * rowSums(
      sweep(result$X, 2, means[2, ], "-")^2
    )
  )
  expected <- exp(pmax(log_kernel, -700))
  expected <- sweep(expected, 2, probabilities, "*")
  expected <- expected / rowSums(expected)

  expect_equal(result$U, expected, tolerance = 1e-12)
})

test_that("MATLAB-style stacked covariance input is supported", {
  means <- matrix(
    c(
      0, 0,
      1, 1
    ),
    nrow = 2,
    byrow = TRUE
  )
  covariance_list <- list(
    matrix(c(1.0, 0.2, 0.2, 2.0), 2, 2),
    matrix(c(2.0, 0.1, 0.1, 1.5), 2, 2)
  )
  stacked <- do.call(rbind, covariance_list)

  set.seed(789)
  from_list <- generate_gaussian_mixture(
    30,
    c(0.5, 0.5),
    means,
    covariance_list
  )

  set.seed(789)
  from_stack <- generate_gaussian_mixture(
    30,
    c(0.5, 0.5),
    means,
    stacked
  )

  expect_equal(from_stack, from_list)
})

test_that("genmixhet is a deterministic compatibility wrapper", {
  means <- matrix(c(0, 1), nrow = 2)
  covariance <- list(matrix(1), matrix(1))

  set.seed(101)
  legacy <- genmixhet(
    n = 10,
    p = c(0.5, 0.5),
    M = means,
    Sig = covariance
  )

  set.seed(101)
  modern <- generate_gaussian_mixture(
    n = 10,
    probabilities = c(0.5, 0.5),
    means = means,
    covariances = covariance
  )

  expect_equal(legacy, modern)
})

test_that("Gaussian mixture inputs are validated", {
  means <- matrix(c(0, 0, 1, 1), nrow = 2, byrow = TRUE)

  expect_error(
    generate_gaussian_mixture(
      10,
      c(0.2, 0.2),
      means,
      list(diag(2), diag(2))
    ),
    "sum to one"
  )

  expect_error(
    generate_gaussian_mixture(
      10,
      c(0.5, 0.5),
      means,
      list(diag(2), matrix(c(1, 2, 2, 1), 2, 2))
    ),
    "positive definite"
  )

  expect_error(
    generate_gaussian_mixture(
      0,
      c(0.5, 0.5),
      means,
      list(diag(2), diag(2))
    ),
    "positive integer"
  )
})
