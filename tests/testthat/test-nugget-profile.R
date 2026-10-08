test_that("tau zero nugget likelihood matches direct Kronecker evaluation", {
  set.seed(3201)

  SV <- matrix(c(1.5, 0.2, 0.2, 1), 2, 2)
  SO <- matrix(c(1.2, 0.1, 0.1, 0.9), 2, 2)
  means <- rbind(c(-1, 0, 0, 0), c(1, 0, 0, 0))
  probabilities <- c(0.4, 0.6)
  X <- matrix(rnorm(40 * 4), nrow = 40)

  result <- scr_tucker3_nugget_loglik(
    X = X,
    means = means,
    probabilities = probabilities,
    variable_covariance = SV,
    occasion_covariance = SO,
    nugget = 0
  )

  covariance <- kronecker(SO, SV)
  inverse <- solve(covariance)
  logdet <- as.numeric(
    determinant(covariance, logarithm = TRUE)$modulus
  )

  log_kernel <- matrix(0, nrow = nrow(X), ncol = 2)
  for (g in 1:2) {
    centered <- sweep(X, 2, means[g, ], "-")
    quadratic <- rowSums((centered %*% inverse) * centered)
    log_kernel[, g] <- log(probabilities[g]) -
      0.5 * logdet -
      0.5 * quadratic
  }

  maxima <- apply(log_kernel, 1, max)
  direct <- sum(
    maxima +
      log(rowSums(exp(sweep(log_kernel, 1, maxima, "-"))))
  )

  expect_equal(result$log_likelihood, direct, tolerance = 1e-10)
  expect_equal(rowSums(result$membership), rep(1, nrow(X)), tolerance = 1e-12)
})


test_that("profiled nugget uses exactly one additional parameter", {
  set.seed(3202)

  fit <- list(
    M = rbind(c(-1, 0, 0, 0), c(1, 0, 0, 0)),
    SV = diag(2),
    SO = diag(2),
    probabilities = c(0.5, 0.5),
    ranks = c(P = 1L, Q = 1L, R = 1L)
  )
  X <- matrix(rnorm(60 * 4), nrow = 60)

  profile <- profile_scr_tucker3_nugget(
    X = X,
    fit = fit,
    upper = 2
  )

  expect_equal(
    unname(profile$parameters[["nugget"]] -
      profile$parameters[["separable"]]),
    1L
  )
  expect_equal(
    profile$bic,
    2 * profile$log_likelihood -
      log(nrow(X)) * profile$parameters[["nugget"]],
    tolerance = 1e-10
  )
})


test_that("nugget-contaminated synthetic mixture profiles to positive tau", {
  set.seed(3203)

  probabilities <- c(0.5, 0.5)
  means <- rbind(c(-2, 0, 0, 0), c(2, 0, 0, 0))
  base <- kronecker(diag(2), diag(2))
  covariance <- base + 0.8 * diag(4)

  generated <- generate_gaussian_mixture(
    n = 400,
    probabilities = probabilities,
    means = means,
    covariances = list(covariance, covariance)
  )

  fit <- list(
    M = means,
    SV = diag(2),
    SO = diag(2),
    probabilities = probabilities,
    ranks = c(P = 1L, Q = 1L, R = 1L)
  )

  profile <- profile_scr_tucker3_nugget(
    X = generated$X,
    fit = fit,
    upper = 3,
    boundary_tolerance = 1e-6
  )

  expect_gt(profile$nugget, 0)
  expect_gt(profile$likelihood_improvement, 0)
})


test_that("profile posterior memberships are normalized and ICL is finite", {
  set.seed(3204)

  fit <- list(
    M = rbind(c(-1, 0, 0, 0), c(1, 0, 0, 0)),
    SV = matrix(c(1.3, 0.1, 0.1, 1), 2, 2),
    SO = matrix(c(1.2, 0.05, 0.05, 0.9), 2, 2),
    probabilities = c(0.45, 0.55),
    ranks = c(P = 1L, Q = 1L, R = 1L)
  )
  X <- matrix(rnorm(80 * 4), nrow = 80)

  profile <- profile_scr_tucker3_nugget(
    X = X,
    fit = fit
  )

  expect_equal(
    rowSums(profile$membership),
    rep(1, nrow(X)),
    tolerance = 1e-10
  )
  expect_true(is.finite(profile$icl))
  expect_true(is.finite(profile$bic))
  expect_gte(profile$nugget, 0)
})


test_that("automatic nugget search is deterministic", {
  set.seed(3205)

  fit <- list(
    M = rbind(c(-1, 0, 0, 0), c(1, 0, 0, 0)),
    SV = diag(2),
    SO = diag(2),
    probabilities = c(0.5, 0.5),
    ranks = c(P = 1L, Q = 1L, R = 1L)
  )
  X <- matrix(rnorm(70 * 4), nrow = 70)

  first <- profile_scr_tucker3_nugget(X = X, fit = fit)
  second <- profile_scr_tucker3_nugget(X = X, fit = fit)

  expect_equal(first$nugget, second$nugget)
  expect_equal(first$search, second$search)
  expect_equal(first$log_likelihood, second$log_likelihood)
})


test_that("nugget profile validates fit objects and bounds", {
  X <- matrix(rnorm(40), nrow = 10, ncol = 4)

  expect_error(
    profile_scr_tucker3_nugget(
      X = X,
      fit = list()
    ),
    "fit must be a Tucker3 fit"
  )

  fit <- list(
    M = rbind(rep(0, 4), rep(1, 4)),
    SV = diag(2),
    SO = diag(2),
    probabilities = c(0.5, 0.5),
    ranks = c(P = 1L, Q = 1L, R = 1L)
  )

  expect_error(
    profile_scr_tucker3_nugget(
      X = X,
      fit = fit,
      upper = 0
    ),
    "upper must be NULL or a positive"
  )
})
