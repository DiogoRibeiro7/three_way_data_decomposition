test_that("mpca_decomposition preserves the measurement mode", {
  set.seed(123)
  x <- array(rnorm(3 * 4 * 5), dim = c(3, 4, 5))

  result <- suppressMessages(
    mpca_decomposition(
      x,
      ranks = c(2, 2),
      max_iter = 4,
      tol = 1e-4
    )
  )

  expect_named(
    result,
    c("Z_ext", "U", "conv", "est", "norm_percent", "fnorm_resid", "all_resids")
  )
  expect_equal(result$Z_ext@modes, c(2L, 2L, 5L))
  expect_length(result$U, 3L)
  expect_equal(dim(result$U[[1L]]), c(3L, 2L))
  expect_equal(dim(result$U[[2L]]), c(4L, 2L))
  expect_true(is.logical(result$conv))
  expect_true(is.finite(result$fnorm_resid))
})

test_that("mpca_decomposition accepts an rTensor Tensor", {
  set.seed(321)
  x <- rTensor::as.tensor(
    array(rnorm(2 * 3 * 4), dim = c(2, 3, 4))
  )

  result <- suppressMessages(
    mpca_decomposition(
      x,
      ranks = c(1, 2),
      max_iter = 4,
      tol = 1e-4
    )
  )

  expect_equal(result$Z_ext@modes, c(1L, 2L, 4L))
})

test_that("mpca_decomposition validates ranks and tuning parameters", {
  x <- array(seq_len(60), dim = c(3, 4, 5))

  expect_error(
    mpca_decomposition(x, ranks = 2),
    "one positive integer for each compressed tensor mode"
  )
  expect_error(
    mpca_decomposition(x, ranks = c(4, 2)),
    "cannot exceed"
  )
  expect_error(
    mpca_decomposition(x, ranks = c(2, 2), max_iter = 1),
    "greater than or equal to 2"
  )
  expect_error(
    mpca_decomposition(x, ranks = c(2, 2), tol = 0),
    "positive finite number"
  )
})
