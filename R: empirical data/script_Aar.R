##################################################
### Empirical data set: Aar massif
##################################################

### Clean workspace
rm(list = ls())


# R version: R version 4.0.5 (2021-03-31) Shake and Throw (DELL)
# R vevrsion: R version 4.2.1 (2022-06-23) Funny-Looking Kid (Mac)

library(compositions)   # to handle compositional data; Aar massif data set; version 2.0-4
library(sparsepca)      # function spca(); version 0.1.2
library(easyCODA)       # function STEP(); version 0.35.2
library(tictoc)         # simulation time measurement; version 1.0.1
library(ggplot2)
library(stringr)
library(tidyr)
library(tibble)


### Source functions -----------------------------------------------------------

# Key function: a function making a matrix of D(D-1) logratios and calculating sparse PCA
source("spca_function.R")


# Auxiliary functions for calculations
source("Auxiliary_functions_spca.R")

### Data ------------------------------------------------------------------------------------------------------------

data(Aar)

data <- Aar

data2 <- (data[,14:29]/1000000)*100     # express trace elements in %

X <- cbind(data[,1:12],data2)  # components + trace elements, only %, not columns with sum
X1 <- X[,c(1:12)]              # without trace elements



D <- ncol(X)
D1 <- ncol(X1)

## How many zeros are in the data without trace elements:
sapply(X1[,-c(1,2)], function(x){ length(which(x==0))})  # there are no zeros


### SPCA models ----------------------------------------------------------------

alpha_max <- 0.244205309 # specify max value of tuning parameter
alpha_nbr <- 50 # specify number of tuning parameters
alpha_ratio <- 1000 # specify ratio of largest to smallest tuning parameter
alpha_grid <- c(exp(seq(from = log(alpha_max), to = log(alpha_max/alpha_ratio), length = alpha_nbr)), 0)

# a vector of sparsity parameters alpha
a <- sort(alpha_grid,decreasing=F)        # zero (-> fully dense model) included


models <- list()
for(i in 1:length(a)){
  models[[i]] <- spca_funkce(X1[,-c(1,2)],alpha = a[i], k = 2, draw = F)  # first two variables are not comp. parts
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
total.variance <- sum(logratio.vars)/(2*10)   # totvar(X); 2*10 = 2*D -> 10 components 
percent <- (logratio.vars/total.variance)*100

# pick just odd pairwise log-ratios (var(ln(a/b)) = var(ln(b/a)))
exvar <- round(sort(percent,decreasing=T)[c(TRUE,FALSE)],3)  # round sorted odd pairwise logratios


# percentage of zero logratios
D2 <- ncol(X1[,-c(1,2)])
nlogr <- D2*(D2-1)
prop <- (zeros/nlogr)*100  

tab1 <- as.data.frame(cbind(a,variances,zeros,prop))
colnames(tab1) <- c("alpha", "var", "ZeroLogr","Percentage")


### Visualisation --------------------------------------------------------------

# a) Figure 6: development of sparsity - percentage of zero PLRs and expl. variance

plot_line(input = tab1, a=a)

# b) Figure 7: stability of pairwise logratios

plot_heat_PLRs(input = models, a = a, X = X1[,-c(1,2)], exvar = exvar)

# c) Figure 8: stability of parts: 

plot_heat_parts(input = models, a = a, X = X1[,-c(1,2)])














