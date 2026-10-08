test_that("Student-t posterior memberships are normalized", {
  set.seed(4001)

  X <- matrix(rnorm(50 * 4), nrow = 50)
  means <- rbind(c(-1, 0, 0, 0), c(1, 0, 0, 0))

  result <- scr_tucker3_student_loglik(
    X = X,
    means = means,
    probabilities = c(0.4, 0.6),
    variable_scale = diag(2),
    occasion_scale = diag(2),
    degrees_freedom = 5
  )

  expect_equal(
    rowSums(result$membership),
    rep(1, nrow(X)),
    tolerance = 1e-12
  )
  expect_true(is.finite(result$log_likelihood))
})


test_that("Student-t robustness weights decrease with Mahalanobis distance", {
  X <- rbind(
    c(0, 0, 0, 0),
    c(4, 0, 0, 0)
  )
  means <- rbind(
    c(0, 0, 0, 0),
    c(10, 0, 0, 0)
  )

  result <- scr_tucker3_student_loglik(
    X = X,
    means = means,
    probabilities = c(0.5, 0.5),
    variable_scale = diag(2),
    occasion_scale = diag(2),
    degrees_freedom = 6
  )

  expect_gt(
    result$component_weights[1, 1],
    result$component_weights[2, 1]
  )
})


test_that("large degrees of freedom approaches Gaussian package likelihood", {
  set.seed(4002)

  X <- matrix(rnorm(80 * 4), nrow = 80)
  means <- rbind(c(-1, 0, 0, 0), c(1, 0, 0, 0))
  probabilities <- c(0.45, 0.55)
  SV <- matrix(c(1.3, 0.1, 0.1, 1), 2, 2)
  SO <- matrix(c(1.2, 0.05, 0.05, 0.9), 2, 2)

  student <- scr_tucker3_student_loglik(
    X = X,
    means = means,
    probabilities = probabilities,
    variable_scale = SV,
    occasion_scale = SO,
    degrees_freedom = 1e7
  )

  gaussian <- scr_tucker3_nugget_loglik(
    X = X,
    means = means,
    probabilities = probabilities,
    variable_covariance = SV,
    occasion_covariance = SO,
    nugget = 0,
    return_membership = TRUE
  )

  expect_equal(
    student$log_likelihood,
    gaussian$log_likelihood,
    tolerance = 1e-5
  )
  expect_equal(
    student$membership,
    gaussian$membership,
    tolerance = 1e-6
  )
})


test_that("remote contamination receives a smaller observation weight", {
  X <- rbind(
    c(0, 0, 0, 0),
    c(0.2, 0, 0, 0),
    c(15, 15, 15, 15)
  )
  means <- rbind(
    c(0, 0, 0, 0),
    c(5, 0, 0, 0)
  )

  result <- scr_tucker3_student_loglik(
    X = X,
    means = means,
    probabilities = c(0.7, 0.3),
    variable_scale = diag(2),
    occasion_scale = diag(2),
    degrees_freedom = 4
  )

  expect_lt(
    result$observation_weights[3],
    min(result$observation_weights[1:2])
  )
})


test_that("Student-t parameter count matches Gaussian count for fixed nu", {
  gaussian <- scr_tucker3_parameter_count(
    groups = 3,
    variables = 3,
    occasions = 2,
    centroid_rank = 1,
    variable_rank = 1,
    occasion_rank = 1
  )

  fixed <- scr_tucker3_student_parameter_count(
    groups = 3,
    variables = 3,
    occasions = 2,
    centroid_rank = 1,
    variable_rank = 1,
    occasion_rank = 1,
    estimate_degrees_freedom = FALSE
  )

  estimated <- scr_tucker3_student_parameter_count(
    groups = 3,
    variables = 3,
    occasions = 2,
    centroid_rank = 1,
    variable_rank = 1,
    occasion_rank = 1,
    estimate_degrees_freedom = TRUE
  )

  expect_equal(fixed, gaussian)
  expect_equal(estimated, gaussian + 1L)
})


test_that("Student-t likelihood rejects invalid degrees of freedom", {
  X <- matrix(rnorm(40), nrow = 10, ncol = 4)
  means <- rbind(rep(0, 4), rep(1, 4))

  expect_error(
    scr_tucker3_student_loglik(
      X = X,
      means = means,
      probabilities = c(0.5, 0.5),
      variable_scale = diag(2),
      occasion_scale = diag(2),
      degrees_freedom = 2
    ),
    "greater than two"
  )
})
