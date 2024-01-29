
library(huge)
L<-huge.generator(n=20, d=11, graph="random", prob =0.2)

X<-L$data
Rho<-cor(X) # correlation matrix of the data
plot(L) # True graph

# Load the 'MASS' package for mvrnorm function
library(MASS)

# Function to simulate data given a precision matrix
simulate_data <- function(n, precision_matrix) {
  # Generate multivariate normal data with the given precision matrix
  data <- mvrnorm(n = n, mu = rep(0, ncol(precision_matrix)), Sigma = solve(precision_matrix))
  return(data)
}

# Simulate data for n = 10, 100, 1000
n_values <- c(10, 100, 1000)

for (n in n_values) {
  # Simulate Erdös-Renyi random network
  random_network <- huge.generator(n = n, d = 11, graph = "random", prob = 0.2)
  
  # Get precision matrix from the random network
  precision_matrix <- solve(cor(random_network$data))
  
  # Simulate data based on the precision matrix
  simulated_data <- simulate_data(n, precision_matrix)
  
  # You can do further analysis or visualization with the simulated data as needed
  # ...
  
  # Print summary or plot if desired
  cat("Summary for n =", n, ":\n")
  print(summary(simulated_data))
  cat("\n")
}


# Load the 'MASS' package for mvrnorm function
library(MASS)

# Function to simulate data given a precision matrix
simulate_data <- function(n, precision_matrix) {
  # Generate multivariate normal data with the given precision matrix
  data <- mvrnorm(n = n, mu = rep(0, ncol(precision_matrix)), Sigma = solve(precision_matrix))
  return(data)
}

# Simulate data for n = 10, 100, 1000
n_values <- c(10, 100, 1000)

for (n in n_values) {
  # Simulate Erdös-Renyi random network
  random_network <- huge.generator(n = n, d = 11, graph = "random", prob = 0.2)
  
  # Get correlation matrix from the random network data
  correlation_matrix <- cor(random_network$data)
  
  # Add a small positive constant to the diagonal for regularization
  correlation_matrix <- correlation_matrix + diag(length(correlation_matrix)) * 1e-6
  
  # Get precision matrix from the regularized correlation matrix
  precision_matrix <- solve(correlation_matrix)
  
  # Simulate data based on the precision matrix
  simulated_data <- simulate_data(n, precision_matrix)
  
  # You can do further analysis or visualization with the simulated data as needed
  # ...
  
  # Print summary or plot if desired
  cat("Summary for n =", n, ":\n")
  print(summary(simulated_data))
  cat("\n")
}


