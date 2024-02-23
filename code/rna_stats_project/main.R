
library(huge)
L<-huge.generator(n=20, d=11, graph="random", prob =0.2)

X<-L$data
Rho<-cor(X) # correlation matrix of the data
plot(L) # True graph


#=======================Simulate 3 Datasets=====================================
##======My way
L10<-huge.generator(n=10, d=11, graph="random", prob =0.2)
L100<-huge.generator(n=100, d=11, graph="random", prob =0.2)
L1000<-huge.generator(n=1000, d=11, graph="random", prob =0.2)

X10<-L10$data
Rho<-cor(X10) # correlation matrix of the data
plot(L10) # True graph

X100<-L100$data
Rho<-cor(X100) # correlation matrix of the data
plot(L100) # True graph

X1000<-L1000$data
Rho<-cor(X1000) # correlation matrix of the data
plot(L1000) # True graph

# Load necessary libraries
library(huge)

##======Their Way 1
# Set parameters
p <- 11 # Number of genes
n_values <- c(10, 100, 1000) # Different sample sizes
prob <- 0.2 # Probability for random network
epsilon <- 1e-6 # Small value for regularization

# Loop through different sample sizes
for (n in n_values) {
  # Generate random network
  L <- huge.generator(n = n, d = p, graph = "random", prob = prob)
  
  # Generate covariance matrix
  cov_matrix <- cov(L$data)
  
  # Regularize covariance matrix
  cov_matrix_reg <- cov_matrix + epsilon * diag(p)
  
  # Generate precision matrix
  precision_matrix <- solve(cov_matrix_reg)
  
  # Generate data from multivariate normal distribution with precision matrix reflecting the network structure
  simulated_data <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = precision_matrix)
  
  # Save or use the simulated data as required
  # For example, you can save the simulated data to a file
  write.csv(simulated_data, file = paste0("simulated_data_n", n, ".csv"), row.names = FALSE)
}

##======Their Way 2
# Load necessary libraries
library(huge)

# Set parameters
p <- 11 # Number of genes
n_values <- c(10, 100, 1000) # Different sample sizes
prob <- 0.2 # Probability for random network
epsilon <- 1e-6 # Small value for regularization

# Initialize list to store simulated datasets
simulated_datasets <- list()

# Loop through different sample sizes
for (n in n_values) {
  # Generate random network
  L <- huge.generator(n = n, d = p, graph = "random", prob = prob)
  
  # Generate covariance matrix
  cov_matrix <- cov(L$data)
  
  # Regularize covariance matrix
  cov_matrix_reg <- cov_matrix + epsilon * diag(p)
  
  # Generate precision matrix
  precision_matrix <- solve(cov_matrix_reg)
  
  # Generate data from multivariate normal distribution with precision matrix reflecting the network structure
  simulated_data <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = precision_matrix)
  
  # Store simulated dataset in list
  simulated_datasets[[paste0("simulated_data_n", n)]] <- simulated_data
}

# Accessing the simulated datasets
simulated_datasets[["simulated_data_n10"]] # Access dataset with n = 10
simulated_datasets[["simulated_data_n100"]] # Access dataset with n = 100
simulated_datasets[["simulated_data_n1000"]] # Access dataset with n = 1000

#=======================Estimate the observed correlation matrix================
# Load necessary libraries
library(MASS)  # For mvrnorm function
library(huge)  # For huge.generator function

# Set parameters
p <- 11  # Number of genes
n_values <- c(10, 100, 1000)  # Different sample sizes
prob <- 0.2  # Probability for random network
epsilon <- 1e-6  # Small value for regularization

# Initialize list to store Frobenius norm values
frobenius_norms <- numeric(length(n_values))

# Loop through different sample sizes
for (i in seq_along(n_values)) {
  n <- n_values[i]
  
  # Generate random network
  L <- huge.generator(n = n, d = p, graph = "random", prob = prob)
  
  # Generate covariance matrix
  cov_matrix <- cov(L$data)
  
  # Regularize covariance matrix
  cov_matrix_reg <- cov_matrix + epsilon * diag(p)
  
  # Generate precision matrix
  precision_matrix <- solve(cov_matrix_reg)
  
  # Generate data from multivariate normal distribution with precision matrix reflecting the network structure
  simulated_data <- mvrnorm(n = n, mu = rep(0, p), Sigma = precision_matrix)
  
  # Calculate observed correlation matrix
  obs_cor_matrix <- cor(simulated_data)
  
  # True correlation matrix
  true_cor_matrix <- cor(L$data)
  
  # Calculate Frobenius norm
  frobenius_norm <- norm(obs_cor_matrix - true_cor_matrix, "F")
  
  # Store Frobenius norm value
  frobenius_norms[i] <- frobenius_norm
}

# Print Frobenius norm values
frobenius_norms

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
  estimated_precision_matrices[[paste0("estimated_precision_matrix_n", n)]] <- estimated_precision_matrix
}

# Accessing the estimated precision matrices
estimated_precision_matrices[["estimated_precision_matrix_n10"]] # Access precision matrix with n = 10
estimated_precision_matrices[["estimated_precision_matrix_n100"]] # Access precision matrix with n = 100
estimated_precision_matrices[["estimated_precision_matrix_n1000"]] # Access precision matrix with n = 1000


# Define true network structure (example, replace with your actual true network)
true_network <- matrix(c(
  0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0,
  1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
  1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
  1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
  0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0,
  0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0,
  0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0,
  0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1,
  0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0,
  0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1
), nrow = 11, byrow = TRUE)

# Initialize counts
true_positives <- 0
true_negatives <- 0
false_positives <- 0
false_negatives <- 0

# Loop through different sample sizes
for (n in n_values) {
  # Get estimated precision matrix
  estimated_precision_matrix <- estimated_precision_matrices[[paste0("estimated_precision_matrix_n", n)]]
  
  # Compare with true network structure
  for (i in 1:nrow(true_network)) {
    for (j in 1:ncol(true_network)) {
      if (true_network[i, j] == 1) {
        # True edge
        if (estimated_precision_matrix$path[[1]][i, j] == 1) {
          true_positives <- true_positives + 1
        } else {
          false_negatives <- false_negatives + 1
        }
      } else {
        # True non-edge
        if (estimated_precision_matrix$path[[1]][i, j] == 0) {
          true_negatives <- true_negatives + 1
        } else {
          false_positives <- false_positives + 1
        }
      }
    }
  }
}

# Output counts
cat("True positives:", true_positives, "\n")
cat("True negatives:", true_negatives, "\n")
cat("False positives:", false_positives, "\n")
cat("False negatives:", false_negatives, "\n")

























