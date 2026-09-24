test_that("S3 fits a small three-way Gaussian mixture", {
  set.seed(123)

  n <- 80
  groups <- 2
  variables <- 3
  occasions <- 2
  variable_rank <- 2
  occasion_rank <- 1

  B <- qr.Q(qr(matrix(rnorm(variables * variable_rank), variables)))
  C <- qr.Q(qr(matrix(rnorm(occasions * occasion_rank), occasions)))

  latent_means <- matrix(
    c(
      -2, -1,
       2,  1
    ),
    nrow = groups,
    byrow = TRUE
  )
  means <- latent_means %*% kronecker(t(C), t(B))

  covariance <- kronecker(diag(occasions), diag(variables))
  generated <- generate_gaussian_mixture(
    n = n,
    probabilities = c(0.5, 0.5),
    means = means,
    covariances = list(covariance, covariance)
  )

  U <- matrix(runif(n * groups), nrow = n, ncol = groups)
  U <- U / rowSums(U)

  TB <- qr.Q(qr(matrix(rnorm(variables * variable_rank), variables)))
  TC <- qr.Q(qr(matrix(rnorm(occasions * occasion_rank), occasions)))

  fit <- fit_scr_s3(
    X = generated$X,
    membership = U,
    variable_basis = TB,
    occasion_basis = TC,
    variable_covariance = diag(variables),
    occasion_covariance = diag(occasions),
    tolerance = 1e-6,
    max_iter = 100
  )

  expect_equal(dim(fit$U), c(n, groups))
  expect_equal(dim(fit$TB), c(variables, variable_rank))
  expect_equal(dim(fit$TC), c(occasions, occasion_rank))
  expect_equal(dim(fit$SV), c(variables, variables))
  expect_equal(dim(fit$SO), c(occasions, occasions))
  expect_equal(
    dim(fit$Y),
    c(groups, variable_rank * occasion_rank)
  )
  expect_equal(dim(fit$M), c(groups, variables * occasions))
  expect_equal(rowSums(fit$U), rep(1, n), tolerance = 1e-12)
  expect_equal(sum(fit$probabilities), 1, tolerance = 1e-12)
  expect_true(all(eigen(fit$SV, symmetric = TRUE)$values > 0))
  expect_true(all(eigen(fit$SO, symmetric = TRUE)$values > 0))
  expect_true(is.finite(fit$like))
  expect_true(is.finite(fit$bic))
})

test_that("S3 bases satisfy their covariance-metric normalizations", {
  set.seed(456)

  n <- 100
  variables <- 3
  occasions <- 3

  X <- matrix(rnorm(n * variables * occasions), nrow = n)
  U <- matrix(runif(n * 2), nrow = n, ncol = 2)
  U <- U / rowSums(U)

  TB <- qr.Q(qr(matrix(rnorm(variables * 2), variables, 2)))
  TC <- qr.Q(qr(matrix(rnorm(occasions * 2), occasions, 2)))

  fit <- fit_scr_s3(
    X = X,
    membership = U,
    variable_basis = TB,
    occasion_basis = TC,
    variable_covariance = diag(variables),
    occasion_covariance = diag(occasions),
    tolerance = 1e-6,
    max_iter = 100
  )

  expect_equal(
    crossprod(fit$TB, solve(fit$SV, fit$TB)),
    diag(2),
    tolerance = 1e-6
  )
  expect_equal(
    crossprod(fit$TC, solve(fit$SO, fit$TC)),
    diag(2),
    tolerance = 1e-6
  )
})

test_that("S3 reconstructed means have Tucker2 form", {
  set.seed(789)

  variables <- 3
  occasions <- 2
  X <- matrix(rnorm(60 * variables * occasions), nrow = 60)
  U <- matrix(runif(60 * 2), nrow = 60, ncol = 2)
  U <- U / rowSums(U)

  TB <- qr.Q(qr(matrix(rnorm(variables * 2), variables, 2)))
  TC <- matrix(c(1, 0), ncol = 1)

  fit <- fit_scr_s3(
    X,
    U,
    TB,
    TC,
    diag(variables),
    diag(occasions),
    tolerance = 1e-6,
    max_iter = 100
  )

  expected <- fit$Y %*% kronecker(t(fit$TC), t(fit$TB))
  expect_equal(fit$M, expected, tolerance = 1e-10)
})

test_that("S3 BIC matches the MATLAB parameter-count formula", {
  set.seed(101)

  n <- 80
  groups <- 3
  variables <- 4
  occasions <- 3
  q <- 2
  r <- 2

  X <- matrix(rnorm(n * variables * occasions), nrow = n)
  U <- matrix(runif(n * groups), nrow = n, ncol = groups)
  U <- U / rowSums(U)
  TB <- qr.Q(qr(matrix(rnorm(variables * q), variables, q)))
  TC <- qr.Q(qr(matrix(rnorm(occasions * r), occasions, r)))

  fit <- fit_scr_s3(
    X,
    U,
    TB,
    TC,
    diag(variables),
    diag(occasions),
    tolerance = 1e-6,
    max_iter = 100
  )

  n_parameters <- groups - 1 +
    variables * occasions +
    (groups - 1) * q * r +
    (variables - q) * q +
    (occasions - r) * r +
    (variables^2 + variables) / 2 +
    (occasions^2 + occasions) / 2 -
    1

  expect_equal(
    fit$bic,
    2 * fit$like - log(n) * n_parameters
  )
})

test_that("t3mixs preserves the historical eight-field interface", {
  set.seed(202)

  variables <- 3
  occasions <- 2
  X <- matrix(rnorm(60 * variables * occasions), nrow = 60)
  U <- matrix(runif(60 * 2), nrow = 60, ncol = 2)
  U <- U / rowSums(U)
  TB <- qr.Q(qr(matrix(rnorm(variables * 2), variables, 2)))
  TC <- matrix(c(1, 0), ncol = 1)

  result <- t3mixs(
    X = X,
    U = U,
    TB = TB,
    TC = TC,
    SV = diag(variables),
    SO = diag(occasions),
    eps = 1e-6,
    dis = 0,
    max_iter = 100
  )

  expect_named(
    result,
    c("U", "TB", "TC", "SO", "SV", "Y", "like", "bic")
  )
  expect_equal(rowSums(result$U), rep(1, nrow(X)), tolerance = 1e-12)
})

test_that("S3 validates dimensions and covariance inputs", {
  X <- matrix(rnorm(60), nrow = 10, ncol = 6)
  U <- matrix(0.5, nrow = 10, ncol = 2)
  TB <- matrix(rnorm(3 * 2), nrow = 3, ncol = 2)
  TC <- matrix(c(1, 0), nrow = 2, ncol = 1)

  expect_error(
    fit_scr_s3(
      X[, 1:5],
      U,
      TB,
      TC,
      diag(3),
      diag(2)
    ),
    "ncol\(X\)"
  )

  expect_error(
    fit_scr_s3(
      X,
      U,
      TB,
      TC,
      matrix(c(1, 2, 2, 1), 2, 2),
      diag(2)
    ),
    "dimensions must match"
  )

  bad_SV <- diag(3)
  bad_SV[3, 3] <- 0
  expect_error(
    fit_scr_s3(
      X,
      U,
      TB,
      TC,
      bad_SV,
      diag(2)
    ),
    "positive definite"
  )
})
