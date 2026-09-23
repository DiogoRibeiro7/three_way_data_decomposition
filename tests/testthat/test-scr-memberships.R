test_that("hard_partition matches the MATLAB ftoh convention", {
  membership <- matrix(
    c(
      0.1, 0.6, 0.3,
      0.7, 0.2, 0.1,
      0.5, 0.5, 0.0
    ),
    nrow = 3,
    byrow = TRUE
  )

  result <- hard_partition(membership)

  expect_equal(
    result,
    matrix(
      c(
        0L, 1L, 0L,
        1L, 0L, 0L,
        1L, 0L, 0L
      ),
      nrow = 3,
      byrow = TRUE
    )
  )
  expect_equal(ftoh(membership), result)
})

test_that("adjusted_rand_index is one for identical hard partitions", {
  membership <- matrix(
    c(
      1, 0, 0,
      1, 0, 0,
      0, 1, 0,
      0, 1, 0,
      0, 0, 1,
      0, 0, 1
    ),
    nrow = 6,
    byrow = TRUE
  )

  contingency <- crossprod(membership, membership)

  expect_equal(adjusted_rand_index(contingency), 1)
  expect_equal(mrand(contingency), 1)
})

test_that("SCR contingency orientation uses hard-membership crossproducts", {
  truth <- matrix(
    c(
      1, 0, 0,
      1, 0, 0,
      0, 1, 0,
      0, 1, 0,
      0, 0, 1,
      0, 0, 1
    ),
    nrow = 6,
    byrow = TRUE
  )
  estimate <- matrix(
    c(
      1, 0, 0,
      0, 1, 0,
      0, 1, 0,
      1, 0, 0,
      0, 0, 1,
      0, 0, 1
    ),
    nrow = 6,
    byrow = TRUE
  )

  contingency <- crossprod(truth, estimate)

  n <- sum(contingency)
  row_pairs <- sum(choose(rowSums(contingency), 2))
  col_pairs <- sum(choose(colSums(contingency), 2))
  cell_pairs <- sum(choose(as.vector(contingency), 2))
  total_pairs <- choose(n, 2)
  expected <- row_pairs * col_pairs / total_pairs
  reference <- (cell_pairs - expected) /
    (0.5 * (row_pairs + col_pairs) - expected)

  expect_equal(adjusted_rand_index(contingency), reference)
  expect_equal(
    adjusted_rand_index(
      crossprod(hard_partition(truth), hard_partition(estimate))
    ),
    reference
  )
})

test_that("membership metrics validate malformed inputs", {
  expect_error(hard_partition(c(0.2, 0.8)), "numeric matrix")
  expect_error(
    hard_partition(matrix(c(0.2, NA_real_), nrow = 1)),
    "finite values"
  )
  expect_error(
    adjusted_rand_index(matrix(c(1, -1, 0, 1), nrow = 2)),
    "non-negative"
  )
  expect_error(
    adjusted_rand_index(matrix(1, nrow = 1)),
    "at least two observations"
  )
})
