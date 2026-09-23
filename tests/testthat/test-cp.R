test_that("CPfunc adapts the rTensor CP result", {
  set.seed(123)
  x <- array(rnorm(2 * 3 * 4), dim = c(2, 3, 4))

  result <- CPfunc(
    x,
    R = 1,
    max_iter = 4,
    conv_eps = 1e-4
  )

  expect_named(result, c("A", "B", "C", "lambda"))
  expect_equal(dim(result$A), c(2L, 1L))
  expect_equal(dim(result$B), c(3L, 1L))
  expect_equal(dim(result$C), c(4L, 1L))
  expect_length(result$lambda, 1L)
  expect_true(all(is.finite(result$lambda)))
})

test_that("CPfunc accepts an rTensor Tensor", {
  set.seed(321)
  x <- rTensor::as.tensor(array(rnorm(3 * 3 * 3), dim = c(3, 3, 3)))

  result <- CPfunc(
    x,
    R = 1,
    max_iter = 4,
    conv_eps = 1e-4
  )

  expect_equal(dim(result$A), c(3L, 1L))
  expect_equal(dim(result$B), c(3L, 1L))
  expect_equal(dim(result$C), c(3L, 1L))
})

test_that("CPfunc validates tensor order and tuning parameters", {
  expect_error(
    CPfunc(matrix(1:6, nrow = 2), R = 1),
    "only three-mode tensors"
  )

  x <- array(seq_len(24), dim = c(2, 3, 4))

  expect_error(CPfunc(x, R = 0), "positive integer")
  expect_error(CPfunc(x, R = 1.5), "positive integer")
  expect_error(CPfunc(x, R = 1, max_iter = 1), "greater than or equal to 2")
  expect_error(CPfunc(x, R = 1, conv_eps = 0), "positive finite number")
})

test_that("cp_decomposition_wrapper delegates to CPfunc", {
  set.seed(456)
  x <- array(rnorm(2 * 3 * 4), dim = c(2, 3, 4))

  result <- cp_decomposition_wrapper(
    x,
    dims = 1,
    max_iter = 4,
    conv_eps = 1e-4
  )

  expect_named(result, c("A", "B", "C", "lambda"))
  expect_equal(dim(result$A), c(2L, 1L))
  expect_equal(dim(result$B), c(3L, 1L))
  expect_equal(dim(result$C), c(4L, 1L))
})
