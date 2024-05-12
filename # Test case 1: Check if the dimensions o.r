# Test case 1: Check if the dimensions of the Kronecker product are correct
expected_dimensions <- paste(I, "x", J * K)
actual_dimensions <- paste(nrow(test_matrix_A), "x", ncol(test_matrix_A))
assertthat::assert_that(actual_dimensions == expected_dimensions)

# Test case 2: Check if the dimensions of the result after mode-1 multiplication are correct
expected_dimensions <- paste(I, J, K)
actual_dimensions <- paste(dim(test_result_A)[1], dim(test_result_A)[2], dim(test_result_A)[3])
assertthat::assert_that(actual_dimensions == expected_dimensions)

# Test case 3: Check if the result tensor has the correct values
expected_values <- X %*% test_matrix_A
actual_values <- test_result_A
assertthat::assert_that(all.equal(actual_values, expected_values))