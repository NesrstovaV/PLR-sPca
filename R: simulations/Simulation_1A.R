##################################################
### Scenario A: D = 10
##################################################

### Clean workspace
rm(list = ls())


# R version:

library(dplyr)          # version 1.1.0 -> needs to be loader before "plyr" (loaded within other package loading)!
library(tidyr)          # version 1.2.0
library(ggnewscale)     # adding multiple colour scales in ggplot graphs; version: 0.4.9
library(MASS)           # mvrnorm(); version 7.3-57
library(compositions)   # to handle compositional data; version 2.0-4
library(sparsepca)      # function spca(); version 0.1.2
library(easyCODA)       # function STEP(); version 0.35.2
library(tictoc)         # simmulation time measurement; version 1.0.1
library(ggplot2)


### Source functions -----------------------------------------------------------

# Key function: a function making a matrix of D(D-1) logratios and calculating sparse PCA
source("spca_function.R")


# Auxiliary functions for calculations
source("Auxiliary_functions_spca.R")

### Intro ----------------------------------------------------------------------

alpha_max <- 0.281176870      # specify max value of tuning parameter
alpha_nbr <- 50               # specify number of tuning parameters
alpha_ratio <- 1000           # specify ratio of largest to smallest tuning parameter
alpha_grid <- c(exp(seq(from = log(alpha_max), to = log(alpha_max/alpha_ratio), length = alpha_nbr)), 0)

# a vector of sparsity parameters alpha
a <- sort(alpha_grid,decreasing=F)    # zero (-> fully dense model) included


R <- 100  # number of simulation runs
n <- 100  # number of observations

D1 <- 10  # number of components (parts)

### Data simulation parameters -------------------------------------------------
# covariance matrix: diagonal elements = a, nondiagonal elements = b
# 5 balances relevant, 4 balances noise

# SBP codes:
codes1 <- rbind(c(rep(1,5),rep(-1,5)), c(1,1,1,1,-1,rep(0,5)), c(1,1,1,-1,rep(0,6)), 
                c(1,1,-1,rep(0,7)), c(1,-1,rep(0,8)),c(rep(0,5),1,-1,-1,-1,-1),c(rep(0,6),1,-1,-1,-1),
                c(rep(0,7),1,-1,-1),c(rep(0,8),1,-1))
codes1

# Relevant balances: 1-5
mu1 <- rep(0,5)
C1 <- diag(1,ncol = D1/2,nrow = D1/2)   # first version of C1 was: 4 = diag, 1.5 = nondiag
for(i in 1:nrow(C1)){
  for(j in 1:ncol(C1)){
    if(i!=j){
      C1[i,j] <- 0.7
    }
  }
}
C1
eigen(C1)$values   # check for positive-definite matrix C1 (eigen values > 0 -> p.d. matrix)

### Simulation runs ------------------------------------------------------------

# MacBook Pro
# Procesor: M1
# RAM: 16GB

# Windows laptop
# Processor: Intel(R) Core(TM) i7-4810MQ
# Cores: 4
# Logical processors: 8
# Freq: 2.80GHz

# Approximate time (MacBook Pro): cca 2.5min for spca_models(); step_models() runs faster 
# Approximate time (Windows laptop): cca 4min for spca_models(); step_models() runs faster 

# Sparse PCA models for different values of sparsity parameter alpha
res1A <- list()
set.seed(1234)
tic()
for(i in 1:R){
  res1A[[i]] <- spca_models(n=n,D=D1,mu=mu1,C=C1,codes=codes1,a=a,nrel=5,nnoise=4)
  print(i)
}
toc() 

# Collect explained variances (in each model) by the first "k" principal components
variances1A <- matrix(NA,nrow=R,ncol=length(a))
for(i in 1:R){
  for(j in 1:length(a)){
    variances1A[i,j] <- sum(res1A[[i]][[j]]$expl.var)
  }
}
means1A <- colMeans(variances1A)

# Collect number of zero logratios in each model
zeros1A <-  matrix(NA,nrow=R,ncol=length(a))
for(i in 1:R){
  for(j in 1:length(a)){
    zeros1A[i,j] <- res1A[[i]][[j]]$`number of zero logratios` 
  }
}
mean1Azeros <- colMeans(zeros1A)

# (average) percentage of zero logratios (counted from mean1Azeros)
nlogr1A <- D1*(D1-1) 
prop1A <- (mean1Azeros/nlogr1A)*100  

tab1A <- as.data.frame(cbind(a,means1A,mean1Azeros,prop1A))
colnames(tab1A) <- c("alpha", "var", "ZeroLogr","Percentage")

### FPR and FNR ----------------------------------------------------------------
# True Negative = these logratios should have zero in corresponding loading
logr.names1A <- rownames(res1A[[1]][[1]]$loadings)

vars1A <- c("x6","x7","x8","x9","x10")
combs1A<-as.data.frame(t(combn(vars1A, 2)))

vec1 <- c()
for(i in 1:nrow(combs1A)){
  vec1[i] <- paste0("ln(",combs1A[i,1],"/",combs1A[i,2],")")
}

vec2<- c()
for(i in 1:nrow(combs1A)){
  vec2[i] <- paste0("ln(",combs1A[i,2],"/",combs1A[i,1],")")
}

tn1A <- c(vec1,vec2)

tf1A <- rep(T,length(logr.names1A))  # T/F vector; logratios for which F holds will be added
names(tf1A) <- logr.names1A
tf1A[tn1A] <- F

rates1A <- matrix(NA,nrow=length(a),ncol = 2)   # average over simulaton runs; final matrix of rates
colnames(rates1A) <- c("FNR","FPR")

mat <- matrix(NA,nrow=R,ncol=2)  # auxiliary matrix
for(j in 1:length(a)){
  for(i in 1:R){
    v <- ifelse(res1A[[i]][[j]]$loadings[,1]==0 & res1A[[i]][[j]]$loadings[,2]==0,0,1)
    r1 <- myFNR(as.numeric(tf1A),v)   # newly programmed functions for FPR and FNR
    r2 <- myFPR(as.numeric(tf1A),v)
    mat[i,] <- c(r1,r2)
  }
  rates1A[j,] <- colMeans(mat)
  mat <- matrix(NA,nrow=R,ncol=2)  # auxiliary matrix
}

### Results for STEP -----------------------------------------------------------

res1A.STEP <- list()
set.seed(1234)
tic()
for(i in 1:R){
  res1A.STEP[[i]] <- step_models(n=n,D=D1,mu=mu1,C=C1,codes=codes1,nrel=5,nnoise=4)
  print(i)
}
toc()  

# FPR and FNR

tn1A;tf1A  # these vectors hold -> results of STEP will be compared with these

rates1A.STEP <- c()  # only two values as it did not run for each alpha 

mat.STEP <- matrix(NA,nrow=R,ncol=2)  # auxiliary matrix
for(i in 1:R){
  ratios <- res1A.STEP[[i]]$names   # save pairwise logratios chosen with STEP
  
  # because STEP choses just A/B (not B/A then), we need to "create" these logratios
  # in our analysis, we consider both A/B and B/A, so to make FNR and FPR comparable
  spl <- strsplit(res1A.STEP[[i]]$names, split = "/")  # split parts of logratios
  df.aux <- data.frame(spl)
  
  ratios2 <- c()
  for(j in 1:ncol(df.aux)){
    ratios2[j] <- paste0("ln(",df.aux[2,j],"/",df.aux[1,j],")")
  }
  v <- c(paste0("ln(",res1A.STEP[[i]]$names,")"), ratios2)  # combine ratios and ratios2
  
  tf.vector <- rep(0,length(tf1A))
  names(tf.vector) <- names(tf1A)
  tf.vector[v] <- 1    # 1 = logratios chosen by STEP, 0 = not chosen
  
  r1 <- myFNR(as.numeric(tf1A),tf.vector)   # newly programmed functions for FPR and FNR
  r2 <- myFPR(as.numeric(tf1A),tf.vector)
  mat.STEP[i,] <- c(r1,r2)
}
rates1A.STEP <- colMeans(mat.STEP)
names(rates1A.STEP) <- c("FNR","FPR") 


### Comparison of sPCA and STEP based on ranks ---------------------------------

# for search, we use complete vector of important PLR, but to compare with STEP, we use the one without duplicities as well

loads1 <- rownames(res1A[[1]][[1]]$loadings)   # in all models, D = 10 (and scenarios A,B,C), loadings have the same names

# remove "duplicite" PLR in loads1
names1 <- loads1[c(1:9,11:18,21:27,31:36,41:45,51:54,61:63,71,72,81)] 




# Choose only names of important pairwise logratios:
important1A <- tf1A[tf1A==T]                                      # contains duplicities (i.e. both ln(A/B) and ln(B/A))
important1A.half <- important1A[names(important1A) %in% names1]   # IMPORTANT PLR without duplicities


resulting.list1A <- list()
for(i in 1:R){
  resulting.list1A[[i]] <- RateImportant(models = res1A[[i]], loads = loads1, names = names1, important = important1A.half, D = D1, a = a)
}


resulting.list1A.STEP <- list()
for(i in 1:R){
  resulting.list1A.STEP[[i]] <- RateImportant.G(models = res1A.STEP[[i]], modelsSPCA = res1A[[i]], important = important1A)   # tot.var = totvar1A[i]
}

cumsum.ideal.1A <- cumsum(1:(D1-1))  # to calculate rank.method minus ideal.rank


tabCumSums.1A <- matrix(0, nrow = R, ncol = D1-1)   # cumulative sum of rank
tabSums.1A <- matrix(0, nrow = R, ncol = D1-1)      # cumulative sum of % of important PLR
tab.Importance <- matrix(0, nrow = R, ncol = D1-1)  # 1 = PLR was important in the run, 0 = PLR was unimportant in the run


for(i in 1:R){
  tabCumSums.1A[i,] <- resulting.list1A[[i]]$cumsum_rank_both
}

for(i in 1:R){
  tabSums.1A[i,] <- resulting.list1A[[i]]$cummulative_sum
}

for(i in 1:R){
  tab.Importance[i,] <- resulting.list1A[[i]]$Is_important
}



tabCumSums.1A.STEP <- matrix(0, nrow = R, ncol = D1-1)
tabSums.1A.STEP <- matrix(0, nrow = R, ncol = D1-1)      
tab.Importance.STEP <- matrix(0, nrow = R, ncol = D1-1)  


for(i in 1:R){
  tabCumSums.1A.STEP[i,] <- resulting.list1A.STEP[[i]]$cumsum_rank_both
}

for(i in 1:R){
  tabSums.1A.STEP[i,] <- resulting.list1A.STEP[[i]]$`cummulative sum`
}

for(i in 1:R){
  tab.Importance.STEP[i,] <- resulting.list1A.STEP[[i]]$Is_important
}

# subtract cumsum.ideal.1A from each row of a matrix tabCumSums.1A
difference.1A <- tabCumSums.1A - rep(cumsum.ideal.1A, each = nrow(tabCumSums.1A))
difference.1A.STEP <- tabCumSums.1A.STEP - rep(cumsum.ideal.1A, each = nrow(tabCumSums.1A))


### Visualisation --------------------------------------------------------------

# a) Fig1: illustrative stability plot: 
plot_heat_PLRs_illustration(input.sPCA = res1A, input.STEP = res1A.STEP, input.STEP2 = resulting.list1A.STEP,a=a)

# b) Fig2 a,b: comparison of sPCA and STEP

plot_line_step(input.sPCA = tabSums.1A, input.STEP = tabSums.1A.STEP, input.importantPLR = important1A.half, D = D1)

plot_boxes(input.sPCA = difference.1A, input.STEP = difference.1A.STEP, D = D1)


# c) Fig3 a,b: development of sparsity - percentage of zero PLRs and expl. variance, FNR and FPR

plot_line(input=tab1A,a=a)
plot_rates(input = rates1A,a=a)


