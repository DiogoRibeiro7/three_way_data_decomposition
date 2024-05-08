# Load the library
library(kernlab)

# Generate some sample data
set.seed(123)
data <- matrix(rnorm(100*2), ncol=2)

# Normalize the data
data <- scale(data)

# Perform Kernel K-means Clustering
# Using the Gaussian Radial Basis Function (RBF) kernel
clustering_result <- kkmeans(x = data, centers = 5, kernel = "rbfdot")

# Extract cluster membership
clusters <- as.integer(clustering_result)

# Plot the results
plot(data, col = clusters, pch = 19, xlab = "Feature 1", ylab = "Feature 2")
legend("topright", legend = unique(clusters), col = unique(clusters), pch = 19, title = "Clusters")


rbf_kernel <- rbfdot(sigma = 0.2)
clustering_result <- kkmeans(x = data, centers = 5, kernel = rbf_kernel)
print(clustering_result)
