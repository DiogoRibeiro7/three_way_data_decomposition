# Tucker Decomposition Algorithm (Tuckal Version 2)
# Fits a Tucker decomposition model to the provided matrices using an iterative approach.
#
# Args:
#   matrixA: The first input matrix.
#   matrixB: The second input matrix.
#   diagonalElements: The diagonal elements of the third input matrix.
#   nComp1: The number of components for matrix A.
#   nComp2: The number of components for matrix B.
#   nComp3: The number of components for matrix C.
#   matrixX1: The first data matrix.
#   matrixX2: The second data matrix.
#   totalSumSquares: The total sum of squares.
#   maxIterations: The maximum number of iterations.
#
# Returns:
#   A list containing the following elements:
#   - matrixA: The matrix A.
#   - matrixB: The matrix B.
#   - matrixC: The matrix C.
#   - coreMatrix: The core matrix G.
#   - explainedVarianceRatio: The explained variance ratio.
#   - iterations: The number of iterations performed.
tuckal <- function(matrixA, matrixB, diagonalElements, nComp1, nComp2, nComp3, matrixX1, matrixX2, totalSumSquares, maxIterations) {
    # Initial setup
    matrixA <- matrixA[, 1:nComp1]
    matrixB <- matrixB[, 1:nComp2]
    matrixC <- diag(nComp3)
    
    threshold <- 0.05
    previousLoss <- 0
    iter <- 0
    convergenceDelta <- Inf
    
    # Iterative process to fit the model
    while (abs(convergenceDelta) >= threshold && iter < maxIterations) {
        # Update step for matrix G
        kroneckerProduct <- kronecker(matrixC, matrixB)
        matrixG <- crossprod(matrixA, matrixX1) %*% kroneckerProduct
        reconstruct <- matrixA %*% matrixG %*% kronecker(t(matrixC), t(matrixB))
        residuals <- matrixX1 - reconstruct
        loss <- sum(diag(tcrossprod(residuals)))
        
        # Update matrix A
        descompX1 <- svd(matrixX1 %*% kroneckerProduct %*% t(matrixX1 %*% kroneckerProduct))
        matrixA <- descompX1$u[, 1:nComp1]
        
        # Update matrix B
        kroneckerProduct <- kronecker(matrixC, matrixA)
        descompX2 <- svd(matrixX2 %*% kroneckerProduct %*% t(matrixX2 %*% kroneckerProduct))
        matrixB <- descompX2$u[, 1:nComp2]
        
        # Compute changes and update control variables
        convergenceDelta <- loss - previousLoss
        previousLoss <- loss
        iter <- iter + 1
    }
    
    # Final calculations
    explainedVarianceRatio <- (totalSumSquares - loss) / totalSumSquares * 100
    
    # Construct result object
    results <- list(
        matrixA = matrixA,
        matrixB = matrixB,
        matrixC = matrixC,
        coreMatrix = matrixG,
        explainedVarianceRatio = explainedVarianceRatio,
        iterations = iter
    )
    
    return(results)
}

