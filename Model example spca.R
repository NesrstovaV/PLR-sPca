###################################################################
## Model example on how to use sparse PCA with pairwise logratios
###################################################################

### Clean workspace
rm(list = ls())

# 1) Load necessary functions and libraries ------------------------------------
library(sparsepca)
library(MASS)

source("Auxiliary_functions_spca.R")
source("spca_function.R")

# 2) Generate sample data ------------------------------------------------------
n <- 100              # observations
D <- 10               # parts/variables
Sig <- diag(D-1)      # positive-definite symmetric matrix -> covariance matrix
mu <- c(rep(0, D-1))  # means of variables

set.seed(1234)
# ilr coordinates
Z <- mvrnorm(n,mu,Sigma = Sig) 

# Z -> CoDa X
V <- ilrBase(D = D)  # ilrBase() in library(compositions)

X <- as.matrix(as.data.frame(acomp(exp(Z%*%t(V)))))




# 3) Apply sPCA to pairwise logratios ------------------------------------------

alpha_max <- 1 # specify max value of tuning parameter
alpha_nbr <- 50 # specify number of tuning parameters
alpha_ratio <- 1000 # specify ratio of largest to smallest tuning parameter
alpha_grid <- c(exp(seq(from = log(alpha_max), to = log(alpha_max/alpha_ratio), length = alpha_nbr)), 0)

a <- sort(alpha_grid,decreasing=F)        # zero included

# calculating PC1 and PC2
models <- list()

for(i in 1:length(a)){
  models[[i]] <- spca_funkce(X=X, alpha = a[i], k = 2, draw = F) 
}



