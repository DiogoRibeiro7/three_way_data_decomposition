test_that("H fits a simple two-component Gaussian mixture", {
  set.seed(123)

  generated <- generate_gaussian_mixture(
    n = 80,
    probabilities = c(0.5, 0.5),
    means = matrix(c(-2, -2, 2, 2), nrow = 2, byrow = TRUE),
    covariances = list(diag(2), diag(2))
  )

  initial <- matrix(runif(80 * 2), nrow = 80, ncol = 2)
  initial <- initial / rowSums(initial)

  fit <- fit_homoscedastic_gaussian_mixture(
    X = generated$X,
    membership = initial,
    tolerance = 1e-7,
    max_iter = 200
  )

  expect_equal(dim(fit$U), c(80L, 2L))
  expect_equal(dim(fit$Mmu), c(2L, 2L))
  expect_equal(dim(fit$Sig), c(2L, 2L))
  expect_equal(rowSums(fit$U), rep(1, 80), tolerance = 1e-12)
  expect_equal(sum(fit$probabilities), 1, tolerance = 1e-12)
  expect_true(all(eigen(fit$Sig, symmetric = TRUE)$values > 0))
  expect_true(is.finite(fit$like))
  expect_true(is.finite(fit$bic))
  expect_true(fit$it <= 200)
})

test_that("H likelihood trajectory is non-decreasing up to numerical tolerance", {
  set.seed(456)

  X <- rbind(
    matrix(rnorm(60, mean = -1), ncol = 2),
    matrix(rnorm(60, mean = 1), ncol = 2)
  )
  initial <- matrix(runif(nrow(X) * 2), nrow = nrow(X), ncol = 2)
  initial <- initial / rowSums(initial)

  fit <- fit_homoscedastic_gaussian_mixture(
    X,
    initial,
    tolerance = 1e-8,
    max_iter = 100
  )

  increments <- diff(fit$log_likelihood_trace)
  expect_true(all(increments >= -1e-7))
})

test_that("mixhom preserves the historical seven-field interface", {
  set.seed(789)

  X <- rbind(
    matrix(rnorm(40, mean = -1), ncol = 2),
    matrix(rnorm(40, mean = 1), ncol = 2)
  )
  U <- matrix(runif(nrow(X) * 2), nrow = nrow(X), ncol = 2)
  U <- U / rowSums(U)

  result <- mixhom(
    X = X,
    U = U,
    eps = 1e-6,
    dis = 0,
    max_iter = 100
  )

  expect_named(
    result,
    c("U", "Mmu", "Sig", "dif", "like", "bic", "it")
  )
  expect_equal(rowSums(result$U), rep(1, nrow(X)), tolerance = 1e-12)
})

test_that("H parameter count matches the MATLAB BIC formula", {
  set.seed(101)

  X <- matrix(rnorm(60), nrow = 20, ncol = 3)
  U <- matrix(runif(20 * 2), nrow = 20, ncol = 2)
  U <- U / rowSums(U)

  fit <- fit_homoscedastic_gaussian_mixture(
    X,
    U,
    tolerance = 1e-6,
    max_iter = 50
  )

  groups <- 2
  dimension <- 3
  n_parameters <- groups - 1 +
    groups * dimension +
    (dimension^2 + dimension) / 2

  expected_bic <- 2 * fit$like - log(nrow(X)) * n_parameters
  expect_equal(fit$bic, expected_bic)
})

test_that("H validates memberships and convergence controls", {
  X <- matrix(rnorm(20), nrow = 10, ncol = 2)
  U <- matrix(0.5, nrow = 10, ncol = 2)

  bad_rows <- U
  bad_rows[1, ] <- c(0.8, 0.8)

  expect_error(
    fit_homoscedastic_gaussian_mixture(X, bad_rows),
    "sum to one"
  )
  expect_error(
    fit_homoscedastic_gaussian_mixture(X, U, tolerance = 0),
    "positive finite"
  )
  expect_error(
    fit_homoscedastic_gaussian_mixture(X, U, max_iter = 0),
    "positive integer"
  )
})
