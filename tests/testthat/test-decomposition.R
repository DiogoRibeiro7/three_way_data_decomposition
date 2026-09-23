test_that("hosvd_decomposition reconstructs a full-rank tensor", {
  x <- array(seq_len(24), dim = c(2, 3, 4))

  result <- suppressMessages(
    hosvd_decomposition(x, ranks = c(2, 3, 4))
  )

  expect_named(result, c("Z", "U", "est", "fnorm_resid"))
  expect_equal(result$Z@modes, c(2L, 3L, 4L))
  expect_length(result$U, 3L)
  expect_lt(result$fnorm_resid, 1e-10)
})

test_that("hosvd_decomposition respects truncated ranks", {
  x <- array(seq_len(24), dim = c(2, 3, 4))

  result <- suppressMessages(
    hosvd_decomposition(x, ranks = c(1, 2, 2))
  )

  expect_equal(result$Z@modes, c(1L, 2L, 2L))
  expect_equal(dim(result$U[[1L]]), c(2L, 1L))
  expect_equal(dim(result$U[[2L]]), c(3L, 2L))
  expect_equal(dim(result$U[[3L]]), c(4L, 2L))
})

test_that("tucker_decomposition returns the expected structure", {
  set.seed(123)
  x <- array(rnorm(3 * 4 * 5), dim = c(3, 4, 5))

  result <- suppressMessages(
    tucker_decomposition(
      x,
      ranks = c(2, 2, 2),
      max_iter = 4,
      tol = 1e-4
    )
  )

  expect_named(
    result,
    c("Z", "U", "conv", "est", "norm_percent", "fnorm_resid", "all_resids")
  )
  expect_equal(result$Z@modes, c(2L, 2L, 2L))
  expect_length(result$U, 3L)
  expect_true(is.logical(result$conv))
  expect_true(is.finite(result$fnorm_resid))
})

test_that("decomposition adapters validate ranks and tuning parameters", {
  x <- array(seq_len(24), dim = c(2, 3, 4))

  expect_error(
    hosvd_decomposition(x, ranks = c(1, 2)),
    "one positive integer for each tensor mode"
  )
  expect_error(
    hosvd_decomposition(x, ranks = c(3, 2, 2)),
    "cannot exceed"
  )
  expect_error(
    tucker_decomposition(x, ranks = c(1, 2, 2), max_iter = 1),
    "greater than or equal to 2"
  )
  expect_error(
    tucker_decomposition(x, ranks = c(1, 2, 2), tol = 0),
    "positive finite number"
  )
})
