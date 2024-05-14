# simulaG5ari function: Simulate ARI results and save to a file
#
# This function simulates ARI (Adjusted Rand Index) results for a specified
# number of scenarios and saves the results to a file named 'G5ari.RData'.
#
# The function iterates over combinations of N (number of observations),
# G (number of groups), nrep (number of repetitions), and dgp (data generating
# processes) to compute and store ARI values.
#
# Args:
#   None
#
# Returns:
#   None. The results are saved to 'G5ari.RData'.

simulaG5ari <- function() {
  ns <- 250
  ARI <- array(0, dim = c(ns, 3, 16))
  idx <- matrix(c(1:4, 5:8, 9:12, 13:16), nrow = 4, byrow = TRUE)
  
  G <- 5
  
  for (N in c(300, 500)) {
    for (nrep in c(1, 3)) {
      r <- determine_r(N, nrep)
      
      for (dgp in 1:4) {
        message(sprintf("N=%d, G=%d, nrep=%d, dgp=%d", N, G, nrep, dgp))
        
        # Scenario 2 (commented out; uncomment if needed)
        # ari <- tryCatch(
        #   {
        #     simula2(N, G, nrep, dgp, ns)
        #   },
        #   error = function(e) {
        #     warning(sprintf("Error in simula2(N=%d, G=%d, nrep=%d, dgp=%d): %s", N, G, nrep, dgp, e$message))
        #     matrix(NA, nrow = ns, ncol = 3)  # Return a matrix of NAs in case of error
        #   }
        # )
        
        # Scenario 1
        ari <- tryCatch(
          {
            simula1(N, G, nrep, dgp, ns)
          },
          error = function(e) {
            warning(sprintf("Error in simula1(N=%d, G=%d, nrep=%d, dgp=%d): %s", N, G, nrep, dgp, e$message))
            matrix(NA, nrow = ns, ncol = 3)  # Return a matrix of NAs in case of error
          }
        )
        
        idd <- idx[r, dgp]
        ARI[,,idd] <- ari
      }
    }
  }
  
  # Save the ARI results to a .RData file
  save(ARI, file = "G5ari.RData")
}

# Determine the index r based on the values of N and nrep
determine_r <- function(N, nrep) {
  if (N == 300 && nrep == 1) {
    return(1)
  } else if (N == 300 && nrep == 3) {
    return(2)
  } else if (N == 500 && nrep == 1) {
    return(3)
  } else if (N == 500 && nrep == 3) {
    return(4)
  } else {
    stop("Invalid combination of N and nrep")
  }
}


# Call the function to run the simulation
simulaG5ari()
