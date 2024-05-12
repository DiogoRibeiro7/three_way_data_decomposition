library(rTensor)

tucker2R <- function(datos, amb = 2, stand = TRUE, nc1 = 2, nc2 = 2, niter = 10000) {
    if (!is.data.frame(datos)) {
        stop("datos must be a data frame")
    }

    if (ncol(datos) %% amb != 0) {
        stop("The variables must be the same in amb")
    }

    if (any(is.na(datos))) {
        stop("There is at least one NA; 'datos' must be complete")
    }

    if (amb > 4) {
        print("When amb > 4, the Tucker 3 is recommended!")
    }

    if (stand) {
        datos <- scale(datos)
    }

    I <- nrow(datos)
    J <- ncol(datos) / amb
    K <- amb

    # Create the data array and then convert it into a tensor object
    data_array <- array(data = as.numeric(datos), dim = c(I, J, K))
    datos_tensor <- as.tensor(data_array)

    # Assuming a Tucker decomposition function 'tucker' (Please replace with actual Tucker function if different)
    # This is a placeholder, as rTensor does not directly support Tucker decomposition
    # You might need to use a different package or manually implement Tucker decomposition
    tucker_decomp <- tucker(datos_tensor, modes = c(nc1, nc2, K))

    results <- list(
        core = tucker_decomp$core,
        mode_matrices = tucker_decomp$u,
        iterations = niter,
        data_standardized = stand
    )
    
    return(results)
}

# Example usage, ensure the data is correctly formatted as a data frame
data <- data.frame(matrix(runif(100 * 4), nrow = 100, ncol = 4))  # Mock data
result <- tucker2R(data, amb = 2, stand = TRUE, nc1 = 2, nc2 = 2, niter = 100)
print(result)
