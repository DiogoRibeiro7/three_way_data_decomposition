test_that("perform_matrix_operations returns consistent Kronecker identities", {
  a <- matrix(c(1, 2), ncol = 1)
  b <- matrix(c(3, 4), ncol = 1)

  result <- perform_matrix_operations(a, b)

  expect_equal(result$result_ab_t, a %*% t(b))
  expect_equal(as.vector(result$result_diff1), c(0, 0, 0, 0))
  expect_equal(as.vector(result$result_diff2), c(0, 0, 0, 0))
})

test_that("perform_matrix_operations rejects non-matrix inputs", {
  expect_error(
    perform_matrix_operations(c(1, 2), matrix(c(3, 4), ncol = 1)),
    "Both arguments must be matrices"
  )
})

test_that("tensor_products returns all ordered third-order products", {
  i1 <- matrix(c(1, 0), ncol = 1)
  i2 <- matrix(c(0, 1), ncol = 1)

  result <- tensor_products(i1, i2)

  expect_named(
    result,
    c(
      "i1xi1xi1",
      "i1xi1xi2",
      "i1xi2xi1",
      "i1xi2xi2",
      "i2xi1xi1",
      "i2xi1xi2",
      "i2xi2xi1",
      "i2xi2xi2"
    )
  )
  expect_true(all(vapply(result, is.matrix, logical(1))))
  expect_true(all(vapply(result, nrow, integer(1)) == 8L))
  expect_true(all(vapply(result, ncol, integer(1)) == 1L))

  expect_equal(as.vector(result$i1xi1xi1), c(1, 0, 0, 0, 0, 0, 0, 0))
  expect_equal(as.vector(result$i2xi2xi2), c(0, 0, 0, 0, 0, 0, 0, 1))
})

test_that("tensor_products rejects non-numeric inputs", {
  expect_error(
    tensor_products(c("a", "b"), c(1, 2)),
    "must be numeric"
  )
})
