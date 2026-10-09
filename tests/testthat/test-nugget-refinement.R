test_that("zero refinement rounds returns one-shot fit and profile", {
  set.seed(3401)

  n <- 80L
  groups <- 2L
  labels <- rep(1:2, each = n / 2)
  means <- rbind(c(-2, 0, 0, 0), c(2, 0, 0, 0))
  X <- means[labels, , drop = FALSE] +
    matrix(rnorm(n * 4, sd = 0.6), nrow = n)

  membership <- matrix(0.1, nrow = n, ncol = groups)
  membership[cbind(seq_len(n), labels)] <- 0.9
  membership <- membership / rowSums(membership)

  result <- refine_scr_tucker3_nugget(
    X = X,
    membership = membership,
    centroid_rank = 1,
    variable_rank = 1,
    occasion_rank = 1,
    variable_covariance = diag(2),
    occasion_covariance = diag(2),
    refinement_rounds = 0,
    max_iter = 6,
    inner_max_iter = 6
  )

  expect_s3_class(result, "scr_tucker3_nugget_refinement")
  expect_equal(result$accepted_rounds, 0L)
  expect_equal(length(result$tau_trace), 1L)
  expect_equal(length(result$log_likelihood_trace), 1L)
  expect_true(result$converged)
})


test_that("accepted nugget refinement objective is non-decreasing", {
  set.seed(3402)

  n <- 100L
  labels <- rep(1:2, each = n / 2)
  means <- rbind(c(-2, 0, 0, 0), c(2, 0, 0, 0))
  X <- means[labels, , drop = FALSE] +
    matrix(rnorm(n * 4, sd = 1.2), nrow = n)

  membership <- matrix(runif(n * 2), nrow = n)
  membership <- membership / rowSums(membership)

  result <- refine_scr_tucker3_nugget(
    X = X,
    membership = membership,
    centroid_rank = 1,
    variable_rank = 1,
    occasion_rank = 1,
    variable_covariance = diag(2),
    occasion_covariance = diag(2),
    refinement_rounds = 3,
    objective_tolerance = 1e-5,
    max_iter = 8,
    inner_max_iter = 8
  )

  differences <- diff(result$log_likelihood_trace)
  expect_true(all(differences >= -1e-5))
  expect_equal(
    length(result$tau_trace),
    result$accepted_rounds + 1L
  )
  expect_equal(length(result$bic_trace), length(result$tau_trace))
  expect_equal(length(result$icl_trace), length(result$tau_trace))
})


test_that("refined memberships remain normalized", {
  set.seed(3403)

  X <- matrix(rnorm(60 * 4), nrow = 60)
  membership <- matrix(runif(60 * 2), nrow = 60)
  membership <- membership / rowSums(membership)

  result <- refine_scr_tucker3_nugget(
    X = X,
    membership = membership,
    centroid_rank = 1,
    variable_rank = 1,
    occasion_rank = 1,
    variable_covariance = diag(2),
    occasion_covariance = diag(2),
    refinement_rounds = 1,
    max_iter = 5,
    inner_max_iter = 5
  )

  expect_equal(
    rowSums(result$membership),
    rep(1, nrow(X)),
    tolerance = 1e-10
  )
})


test_that("nugget refinement is deterministic for fixed inputs", {
  set.seed(3404)

  X <- matrix(rnorm(70 * 4), nrow = 70)
  membership <- matrix(runif(70 * 2), nrow = 70)
  membership <- membership / rowSums(membership)

  first <- refine_scr_tucker3_nugget(
    X = X,
    membership = membership,
    centroid_rank = 1,
    variable_rank = 1,
    occasion_rank = 1,
    variable_covariance = diag(2),
    occasion_covariance = diag(2),
    refinement_rounds = 2,
    max_iter = 5,
    inner_max_iter = 5
  )

  second <- refine_scr_tucker3_nugget(
    X = X,
    membership = membership,
    centroid_rank = 1,
    variable_rank = 1,
    occasion_rank = 1,
    variable_covariance = diag(2),
    occasion_covariance = diag(2),
    refinement_rounds = 2,
    max_iter = 5,
    inner_max_iter = 5
  )

  expect_equal(first$tau_trace, second$tau_trace)
  expect_equal(first$log_likelihood_trace, second$log_likelihood_trace)
  expect_equal(first$membership, second$membership)
})


test_that("nugget-contaminated data can retain positive tau after refinement", {
  set.seed(3405)

  probabilities <- c(0.5, 0.5)
  means <- rbind(c(-2, 0, 0, 0), c(2, 0, 0, 0))
  covariance <- diag(4) + 0.7 * diag(4)

  generated <- generate_gaussian_mixture(
    n = 250,
    probabilities = probabilities,
    means = means,
    covariances = list(covariance, covariance)
  )

  membership <- generated$U

  result <- refine_scr_tucker3_nugget(
    X = generated$X,
    membership = membership,
    centroid_rank = 1,
    variable_rank = 1,
    occasion_rank = 1,
    variable_covariance = diag(2),
    occasion_covariance = diag(2),
    refinement_rounds = 2,
    max_iter = 8,
    inner_max_iter = 8
  )

  expect_gte(result$profile$nugget, 0)
  expect_true(any(result$tau_trace > 0))
})


test_that("refinement step helper rejects worsening objective", {
  decision <- scr3way:::evaluate_nugget_refinement_step(
    previous_loglik = 10,
    proposed_loglik = 9.5,
    objective_tolerance = 0.1
  )

  expect_false(decision$accept)
  expect_equal(decision$decrease, 0.5)
})


test_that("refinement controls validate invalid rounds", {
  X <- matrix(rnorm(40), nrow = 10, ncol = 4)
  membership <- matrix(0.5, nrow = 10, ncol = 2)

  expect_error(
    refine_scr_tucker3_nugget(
      X = X,
      membership = membership,
      centroid_rank = 1,
      variable_rank = 1,
      occasion_rank = 1,
      variable_covariance = diag(2),
      occasion_covariance = diag(2),
      refinement_rounds = -1
    ),
    "refinement_rounds"
  )

  expect_error(
    refine_scr_tucker3_nugget(
      X = X,
      membership = membership,
      centroid_rank = 1,
      variable_rank = 1,
      occasion_rank = 1,
      variable_covariance = diag(2),
      occasion_covariance = diag(2),
      objective_tolerance = -1
    ),
    "objective_tolerance"
  )
})
