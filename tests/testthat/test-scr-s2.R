test_that("S2 fits a reduced Gaussian mixture", {
  set.seed(123)

  generated <- generate_gaussian_mixture(
    n = 100,
    probabilities = c(0.5, 0.5),
    means = matrix(
      c(
        -2, -2, 0, 0,
         2,  2, 0, 0
      ),
      nrow = 2,
      byrow = TRUE
    ),
    covariances = list(diag(4), diag(4))
  )

  initial_membership <- matrix(
    runif(100 * 2),
    nrow = 100,
    ncol = 2
  )
  initial_membership <- initial_membership /
    rowSums(initial_membership)

  initial_basis <- qr.Q(qr(matrix(rnorm(4 * 2), 4, 2)))

  fit <- fit_scr_s2(
    X = generated$X,
    membership = initial_membership,
    basis = initial_basis,
    tolerance = 1e-7,
    max_iter = 200
  )

  expect_equal(dim(fit$U), c(100L, 2L))
  expect_equal(dim(fit$TB), c(4L, 2L))
  expect_equal(dim(fit$SV), c(4L, 4L))
  expect_equal(dim(fit$Y), c(2L, 2L))
  expect_equal(dim(fit$M), c(2L, 4L))
  expect_equal(rowSums(fit$U), rep(1, 100), tolerance = 1e-12)
  expect_equal(sum(fit$probabilities), 1, tolerance = 1e-12)
  expect_true(all(eigen(fit$SV, symmetric = TRUE)$values > 0))
  expect_true(is.finite(fit$like))
  expect_true(is.finite(fit$bic))
  expect_true(fit$it <= 200)
})

test_that("S2 basis satisfies the covariance-metric normalization", {
  set.seed(456)

  X <- matrix(rnorm(120 * 4), nrow = 120, ncol = 4)
  U <- matrix(runif(120 * 3), nrow = 120, ncol = 3)
  U <- U / rowSums(U)
  TB <- qr.Q(qr(matrix(rnorm(4 * 2), 4, 2)))

  fit <- fit_scr_s2(
    X = X,
    membership = U,
    basis = TB,
    tolerance = 1e-6,
    max_iter = 100
  )

  metric_gram <- crossprod(fit$TB, solve(fit$SV, fit$TB))
  expect_equal(metric_gram, diag(2), tolerance = 1e-6)
})

test_that("S2 BIC matches the MATLAB parameter-count formula", {
  set.seed(789)

  X <- matrix(rnorm(80 * 5), nrow = 80, ncol = 5)
  U <- matrix(runif(80 * 3), nrow = 80, ncol = 3)
  U <- U / rowSums(U)
  TB <- qr.Q(qr(matrix(rnorm(5 * 2), 5, 2)))

  fit <- fit_scr_s2(
    X,
    U,
    TB,
    tolerance = 1e-6,
    max_iter = 100
  )

  groups <- 3
  dimension <- 5
  rank <- 2
  n_parameters <- groups - 1 +
    dimension +
    (groups - 1) * rank +
    (dimension - rank) * rank +
    (dimension^2 + dimension) / 2 -
    1

  expect_equal(
    fit$bic,
    2 * fit$like - log(nrow(X)) * n_parameters
  )
})

test_that("t2mixt preserves the historical six-field interface", {
  set.seed(101)

  X <- matrix(rnorm(60 * 4), nrow = 60, ncol = 4)
  U <- matrix(runif(60 * 2), nrow = 60, ncol = 2)
  U <- U / rowSums(U)
  TB <- qr.Q(qr(matrix(rnorm(4 * 2), 4, 2)))

  result <- t2mixt(
    X = X,
    U = U,
    TB = TB,
    eps = 1e-6,
    dis = 0,
    max_iter = 100
  )

  expect_named(result, c("U", "TB", "SV", "Y", "like", "bic"))
  expect_equal(rowSums(result$U), rep(1, nrow(X)), tolerance = 1e-12)
})

test_that("S2 validates membership and basis inputs", {
  X <- matrix(rnorm(40), nrow = 10, ncol = 4)
  U <- matrix(0.5, nrow = 10, ncol = 2)
  TB <- qr.Q(qr(matrix(rnorm(4 * 2), 4, 2)))

  expect_error(
    fit_scr_s2(X, U, matrix(1, nrow = 3, ncol = 1)),
    "ncol\(X\) rows"
  )

  rank_deficient <- cbind(rep(1, 4), rep(1, 4))
  expect_error(
    fit_scr_s2(X, U, rank_deficient),
    "full column rank"
  )

  bad_membership <- U
  bad_membership[1, ] <- c(0.8, 0.8)
  expect_error(
    fit_scr_s2(X, bad_membership, TB),
    "sum to one"
  )
})
