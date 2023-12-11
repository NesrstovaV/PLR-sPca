##################################################
### Empirical data set: Archaeometric data -> easyCODA: cups
##################################################

### Clean workspace
rm(list = ls())

# R version: R version 4.0.5 (2021-03-31) Shake and Throw (DELL)
# R vevrsion: R version 4.2.1 (2022-06-23) Funny-Looking Kid (Mac)

library(compositions)   # to handle compositional data; Aar massif data set; version 2.0-4
library(sparsepca)      # function spca(); version 0.1.2
library(easyCODA)       # function STEP(); data set "cups";version 0.35.2
library(tictoc)         # simmulation time measurement; version 1.0.1
library(ggplot2)
library(stringr)
library(tidyr)
library(tibble)


### Source functions -----------------------------------------------------------

# Key function: a function making a matrix of D(D-1) logratios and calculating sparse PCA
source("spca_function.R")


# Auxiliary functions for calculations
source("Auxiliary_functions_spca.R")

### Data -----------------------------------------------------------------------

data(cups)

X <- cups

D <- ncol(X)
n <- nrow(X)

## How many zeros are in the original data:

# number of zeros in each component (column)
sapply(X, function(x){ length(which(x==0))}) # no zeros


### SPCA models ----------------------------------------------------------------

alpha_max <- 0.281176870 # specify max value of tuning parameter
alpha_nbr <- 50 # specify number of tuning parameters
alpha_ratio <- 1000 # specify ratio of largest to smallest tuning parameter
alpha_grid <- c(exp(seq(from = log(alpha_max), to = log(alpha_max/alpha_ratio), length = alpha_nbr)), 0)

# a vector of sparsity parameters alpha
a <- sort(alpha_grid,decreasing=F)        # zero (-> fully dense model) included

models <- list()
for(i in 1:length(a)){
  models[[i]] <- spca_funkce(X, alpha = a[i], k = 2, draw = F)  # first two variables are not comp. parts
  print(i)
}

variances <- c()  # stores explained variance by PC1 and PC2
for(i in 1:length(a)){
  variances[i] <- sum(models[[i]]$expl.var)
}

zeros <- c()
for(i in 1:length(a)){
  zeros[i] <- models[[i]]$`number of zero logratios`   # "full-zero" logratios
}

X.pairs <- models[[1]]$X.pairwise    # matrix of all pairwise logratios is the same in all models


## Explained variability by each pairwise logratio 
logratio.vars <- apply(X.pairs,2,var)  
total.variance <- sum(logratio.vars)/(2*D)   # 2*11 = 2*D -> 11 components 
percent <- (logratio.vars/total.variance)*100

# pick just odd pairwise log-ratios (var(ln(a/b)) = var(ln(b/a)))
sort(percent,decreasing=T)[c(TRUE,FALSE)]
exvar <- round(sort(percent,decreasing=T)[c(TRUE,FALSE)],3)  # round sorted odd pairwise logratios

# percentage of zero logratios
nlogr <- D*(D-1)
prop <- (zeros/nlogr)*100  

tab <- as.data.frame(cbind(a,variances,zeros,prop))
colnames(tab) <- c("alpha", "var", "ZeroLogr","Percentage")


### Visualisation --------------------------------------------------------------

# a) Figure 9: development of sparsity - percentage of zero PLRs and expl. variance

plot_line(input = tab, a=a)

# b) Figure 10: stability of pairwise logratios

plot_heat_PLRs(input = models, a = a, X = X,exvar = exvar)

# c) Figure 11: stability of parts: 

plot_heat_parts(input = models, a = a, X = X)
