# Load the necessary libraries
library(testthat)
source("tuckal.R")

test_that("tuckal function works correctly with valid input", {
    # Mock matrices to simulate input
    matrixA <- matrix(rnorm(20), nrow = 5, ncol = 4)
    matrixB <- matrix(rnorm(15), nrow = 5, ncol = 3)
    matrixX1 <- matrix(rnorm(50), nrow = 5, ncol = 10)
    matrixX2 <- matrix(rnorm(50), nrow = 5, ncol = 10)
    diagonalElements <- 3
    nComp1 <- 2
    nComp2 <- 2
    nComp3 <- 3
    totalSumSquares <- 1000
    maxIterations <- 10

    # Run the function
    result <- tuckal(matrixA, matrixB, diagonalElements, nComp1, nComp2, nComp3, matrixX1, matrixX2, totalSumSquares, maxIterations)

    # Check if results are as expected
    expect_is(result, "list")
    expect_true(length(result$matrixA) == nComp1 * nrow(matrixA))
    expect_true(length(result$matrixB) == nComp2 * nrow(matrixB))
    expect_true(result$iterations <= maxIterations)
    expect_true(result$explainedVarianceRatio <= 100 && result$explainedVarianceRatio >= 0)
})

test_that("tuckal handles zero iterations correctly", {
    matrixA <- matrix(rnorm(20), nrow = 5, ncol = 4)
    matrixB <- matrix(rnorm(15), nrow = 5, ncol = 3)
    matrixX1 <- matrix(rnorm(50), nrow = 5, ncol = 10)
    matrixX2 <- matrix(rnorm(50), nrow = 5, ncol = 10)
    diagonalElements <- 3
    nComp1 <- 2
    nComp2 <- 2
    nComp3 <- 3
    totalSumSquares <- 1000
    maxIterations <- 0  # Zero iterations

    result <- tuckal(matrixA, matrixB, diagonalElements, nComp1, nComp2, nComp3, matrixX1, matrixX2, totalSumSquares, maxIterations)
    
    # With zero iterations, we expect no change in matrices
    expect_equal(result$iterations, 0)
    expect_true(all(result$matrixA == matrixA[, 1:nComp1]))
    expect_true(all(result$matrixB == matrixB[, 1:nComp2]))
})