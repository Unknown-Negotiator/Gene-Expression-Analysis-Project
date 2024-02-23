#=======================Up to Estimate the observed correlation matrix==========
# Load necessary libraries
library(MASS)  # For mvrnorm function
library(huge)  # For huge.generator function

# Set parameters
p <- 11  # Number of genes
n_values <- c(10, 100, 1000)  # Different sample sizes
prob <- 0.2  # Probability for random network
epsilon <- 1e-6  # Small value for regularization

simulated_datasets <- list()
true_networks <- list()
frobenius_norms <- numeric(length(n_values))

# Loop through different sample sizes
for (i in seq_along(n_values)) {
  n <- n_values[i]
  
  # Generate random network
  true_network <- huge.generator(n = n, d = p, graph = "random", prob = prob)
  
  # Generate covariance matrix
  cov_matrix <- cov(true_network$data)
  
  # Generate data from multivariate normal distribution with the random graph structure
  X <- true_network$data
  # Calculate correlation matrix
  cor_matrix <- cor(X)
  # Threshold correlation matrix to obtain adjacency matrix
  adjacency_matrix <- ifelse(abs(cor_matrix) > 0.1, 1, 0)  # Adjust threshold as needed
  # Save true network adjacency matrix as CSV
  write.csv(adjacency_matrix, file = paste0("true_network_n", n, ".csv"), row.names = FALSE)
  # Store true network in list
  true_networks[[paste0("true_network_n", n)]] <- adjacency_matrix
  
  
  # Regularize covariance matrix
  cov_matrix_reg <- cov_matrix + epsilon * diag(p)
  
  # Generate precision matrix
  precision_matrix <- solve(cov_matrix_reg)
  
  # Generate data from multivariate normal distribution with precision matrix reflecting the network structure
  simulated_data <- mvrnorm(n = n, mu = rep(0, p), Sigma = precision_matrix)
  
  # Store simulated dataset in list
  simulated_datasets[[paste0("simulated_data_n", n)]] <- simulated_data
  
  # Calculate observed correlation matrix
  obs_cor_matrix <- cor(simulated_data)
  
  # True correlation matrix
  true_cor_matrix <- cor(true_network$data)
  
  # Calculate Frobenius norm
  frobenius_norm <- norm(obs_cor_matrix - true_cor_matrix, "F")
  
  # Store Frobenius norm value
  frobenius_norms[i] <- frobenius_norm
}

# Accessing the simulated datasets
simulated_datasets[["simulated_data_n10"]] # Access dataset with n = 10
simulated_datasets[["simulated_data_n100"]] # Access dataset with n = 100
simulated_datasets[["simulated_data_n1000"]] # Access dataset with n = 1000

# Print Frobenius norm values
frobenius_norms

# Accessing the true networks
true_networks[["true_network_n10"]] # True network with n = 10
true_networks[["true_network_n100"]] # True network with n = 100
true_networks[["true_network_n1000"]] # True network with n = 1000


#=======================Use Meinhausen-Bühlman version of glass=================
# Load necessary libraries
library(huge)

# Set penalty parameter
penalty <- 0.3

# Initialize list to store estimated precision matrices
estimated_precision_matrices <- list()

# Loop through different sample sizes
for (n in n_values) {
  # Load simulated data
  simulated_data <- read.csv(paste0("simulated_data_n", n, ".csv"))
  
  # Estimate precision matrix using glasso
  estimated_precision_matrix <- huge.glasso(x = as.matrix(simulated_data), lambda = penalty)
  
  # Store estimated precision matrix in list
  estimated_precision_matrices[[paste0("estimated_precision_matrix_n", n)]] <- estimated_precision_matrix$path[[1]]
  write.csv(estimated_precision_matrix$path[[1]], file = paste0("estimated_precision_matrix_n", n, ".csv"), row.names = FALSE)
}

# Accessing the estimated precision matrices
estimated_precision_matrices[["estimated_precision_matrix_n10"]] # Access precision matrix with n = 10
estimated_precision_matrices[["estimated_precision_matrix_n100"]] # Access precision matrix with n = 100
estimated_precision_matrices[["estimated_precision_matrix_n1000"]] # Access precision matrix with n = 1000

#=======================TP, TN, FP, FN estimations=================
# Initialize vectors to store evaluation metrics
true_positives <- numeric(length(n_values))
true_negatives <- numeric(length(n_values))
false_positives <- numeric(length(n_values))
false_negatives <- numeric(length(n_values))

# Loop through different sample sizes
for (i in seq_along(n_values)) {
  n <- n_values[i]
  
  # Load true network and estimated precision matrix
  true_network_file <- paste0("true_network_n", n, ".csv")
  estimated_precision_file <- paste0("estimated_precision_matrix_n", n, ".csv")
  
  # Check if files exist
  if (!file.exists(true_network_file)) {
    cat("True network file does not exist:", true_network_file, "\n")
    next
  }
  if (!file.exists(estimated_precision_file)) {
    cat("Estimated precision matrix file does not exist:", estimated_precision_file, "\n")
    next
  }
  
  # Read true network and estimated precision matrix
  true_network <- as.matrix(read.csv(true_network_file, header = FALSE))
  estimated_precision_matrix <- tryCatch(
    as.matrix(read.csv(estimated_precision_file, header = FALSE)),
    error = function(e) {
      cat("Error reading estimated precision matrix for sample size", n, ":", conditionMessage(e), "\n")
      return(NULL)
    }
  )
  
  # Check if the estimated precision matrix is NULL
  if (is.null(estimated_precision_matrix)) {
    cat("Estimated precision matrix is NULL for sample size", n, "\n")
    next
  }
  
  # Check if the estimated precision matrix is numeric
  if (!is.numeric(estimated_precision_matrix)) {
    cat("Estimated precision matrix is not numeric for sample size", n, "\n")
    next
  }
  
  # Convert any non-numeric values to NA
  estimated_precision_matrix[!is.numeric(estimated_precision_matrix)] <- NA
  
  # Convert estimated precision matrix to binary adjacency matrix
  estimated_adjacency_matrix <- ifelse(abs(estimated_precision_matrix) > 0, 1, 0)
  
  # Calculate evaluation metrics
  true_positives[i] <- sum(true_network == 1 & estimated_adjacency_matrix == 1, na.rm = TRUE)
  true_negatives[i] <- sum(true_network == 0 & estimated_adjacency_matrix == 0, na.rm = TRUE)
  false_positives[i] <- sum(true_network == 0 & estimated_adjacency_matrix == 1, na.rm = TRUE)
  false_negatives[i] <- sum(true_network == 1 & estimated_adjacency_matrix == 0, na.rm = TRUE)
}

# Display evaluation metrics
evaluation_metrics <- data.frame(
  Sample_Size = n_values,
  True_Positives = true_positives,
  True_Negatives = true_negatives,
  False_Positives = false_positives,
  False_Negatives = false_negatives
)

print(evaluation_metrics)

#+++++++++++++Debugging+++++++++++++++++++++++++++++++++++

# Loop through different sample sizes
for (n in n_values) {
  # Load simulated data
  simulated_data <- read.csv(paste0("simulated_data_n", n, ".csv"))
  
  # Convert to numeric matrix if necessary
  if (!is.numeric(simulated_data)) {
    simulated_data <- as.matrix(simulated_data)
  }
  
  # Estimate precision matrix using glasso
  estimated_precision_matrix <- huge.glasso(x = simulated_data, lambda = penalty)
  
  # Extract the estimated network structure from the object
  estimated_network <- estimated_precision_matrix$path[[1]]
  
  # Save estimated network structure as CSV file
  write.csv(estimated_network, file = paste0("estimated_network_n", n, ".csv"), row.names = FALSE)
}

# Read CSV file into a matrix without headers
estimated_precision_matrix <- as.matrix(read.csv("estimated_precision_matrix_n10.csv", header = FALSE))

# Print the estimated precision matrix
print(estimated_precision_matrix)


# Read CSV file into a matrix without headers
estimated_precision_matrix <- as.matrix(read.csv("estimated_precision_matrix_n10.csv", header = FALSE))

# Remove the first row (headers) if needed
estimated_precision_matrix <- estimated_precision_matrix[-1, ]

# Print the estimated precision matrix
print(estimated_precision_matrix)

estimated_precision_matrix <- apply(estimated_precision_matrix, 2, as.numeric)

# Print the estimated precision matrix
print(estimated_precision_matrix)






