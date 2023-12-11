##################################################
### Scenario B: D = 10
##################################################

### Clean workspace
rm(list = ls())


# R version:

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

alpha_max <- 0.281176870 # specify max value of tuning parameter
alpha_nbr <- 50 # specify number of tuning parameters
alpha_ratio <- 1000 # specify ratio of largest to smallest tuning parameter
alpha_grid <- c(exp(seq(from = log(alpha_max), to = log(alpha_max/alpha_ratio), length = alpha_nbr)), 0)

# a vector of sparsity parameters alpha
a <- sort(alpha_grid,decreasing=F)        # zero (-> fully dense model) included


R <- 100                       # number of simulation runs
n <- 100

D1 <- 10

### Data simulation parameters -------------------------------------------------
# covariance matrix: diagonal elements = a, nondiagonal elements = b
# 7 balances relevant, 2 balances noise

# SBP codes:
codes1B <- rbind(c(rep(1,7),rep(-1,3)), c(rep(1,6),-1,rep(0,3)), c(rep(1,5),-1,rep(0,4)), 
                 c(rep(1,4),-1,rep(0,5)), c(1,1,1,-1,rep(0,6)),c(1,1,-1,rep(0,7)),c(1,-1,rep(0,8)),
                 c(rep(0,7),1,-1,-1),c(rep(0,8),1,-1))
codes1B

# Relevant balances: 1-7
mu1B <- rep(0,7)
C1B <- diag(1,ncol = D1-3,nrow = D1-3)   # first version of C1 was: 4 = diag, 1.5 = nondiag
for(i in 1:nrow(C1B)){
  for(j in 1:ncol(C1B)){
    if(i!=j){
      C1B[i,j] <- 0.7
    }
  }
}
C1B
eigen(C1B)$values  # check for positive definite matrix (all eigen values must be positive)

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
# Approximate time (Windows laptop): cca 6min for spca_models(); step_models() runs faster 

res1B <- list()
set.seed(1234)
tic()
for(i in 1:R){
  res1B[[i]] <- spca_models(n=n,D=D1,mu=mu1B,C=C1B,codes=codes1B,a=a,nrel=7,nnoise=2)
  print(i)
}
toc()  

variances1B <- matrix(NA,nrow=R,ncol=length(a))
for(i in 1:R){
  for(j in 1:length(a)){
    variances1B[i,j] <- sum(res1B[[i]][[j]]$expl.var)
  }
}
means1B <- colMeans(variances1B)

zeros1B <-  matrix(NA,nrow=R,ncol=length(a))
for(i in 1:R){
  for(j in 1:length(a)){
    zeros1B[i,j] <- res1B[[i]][[j]]$`number of zero logratios` 
  }
}
mean1Bzeros <- colMeans(zeros1B)


# (average) percentage of zero logratios (counted from mean1Azeros)
nlogr1B <- D1*(D1-1)
prop1B <- (mean1Bzeros/nlogr1B)*100  

tab1B <- as.data.frame(cbind(a,means1B,mean1Bzeros,prop1B))
colnames(tab1B) <- c("alpha", "var", "ZeroLogr","Percentage")

### FPR and FNR ----------------------------------------------------------------
# True Negative = these logratios should have zero in corresponding loading

logr.names1B <- rownames(res1B[[1]][[1]]$loadings)

vars1B <- c("x8","x9","x10")
combs1B<-as.data.frame(t(combn(vars1B, 2)))

vec1 <- c()
for(i in 1:nrow(combs1B)){
  vec1[i] <- paste0("ln(",combs1B[i,1],"/",combs1B[i,2],")")
}

vec2<- c()
for(i in 1:nrow(combs1B)){
  vec2[i] <- paste0("ln(",combs1B[i,2],"/",combs1B[i,1],")")
}

tn1B <- c(vec1,vec2)

tf1B <- rep(T,length(logr.names1B))  # T/F vector; logratios for which F holds will be added
names(tf1B) <- logr.names1B
tf1B[tn1B] <- F

rates1B <- matrix(NA,nrow=length(a),ncol = 2)   # average over simulaton runs; final matrix of rates
colnames(rates1B) <- c("FNR","FPR")

mat <- matrix(NA,nrow=R,ncol=2)  # auxiliary matrix
for(j in 1:length(a)){
  for(i in 1:R){
    v <- ifelse(res1B[[i]][[j]]$loadings[,1]==0 & res1B[[i]][[j]]$loadings[,2]==0,0,1)
    r1 <- myFNR(as.numeric(tf1B),v)   # newly programmed functions for FPR and FNR
    r2 <- myFPR(as.numeric(tf1B),v)
    mat[i,] <- c(r1,r2)
  }
  rates1B[j,] <- colMeans(mat)
  mat <- matrix(NA,nrow=R,ncol=2)  # auxiliary matrix
}



### Results for STEP -----------------------------------------------------------

res1B.STEP <- list()
set.seed(1234)
tic()
for(i in 1:R){
  res1B.STEP[[i]] <- step_models(n=n,D=D1,mu=mu1B,C=C1B,codes=codes1B,nrel=7,nnoise=2)
  print(i)
}
toc()  

# FPR and FNR

tn1B;tf1B  # these vectors hold -> results of STEP will be compared with these

rates1B.STEP <- c()  # only two values as it did not run for each alpha 

mat.STEP <- matrix(NA,nrow=R,ncol=2)  # auxiliary matrix
for(i in 1:R){
  ratios <- res1B.STEP[[i]]$names   # save pairwise logratios chosen with STEP
  
  # because STEP choses just A/B (not B/A then), we need to "create" these logratios
  # in our analysis, we consider both A/B and B/A, so to make FNR and FPR comparable
  spl <- strsplit(res1B.STEP[[i]]$names, split = "/")  # split parts of logratios
  df.aux <- data.frame(spl)
  
  ratios2 <- c()
  for(j in 1:ncol(df.aux)){
    ratios2[j] <- paste0("ln(",df.aux[2,j],"/",df.aux[1,j],")")
  }
  v <- c(paste0("ln(",res1B.STEP[[i]]$names,")"), ratios2)  # combine ratios and ratios2
  
  tf.vector <- rep(0,length(tf1B))
  names(tf.vector) <- names(tf1B)
  tf.vector[v] <- 1    # 1 = logratios chosen by STEP, 0 = not chosen
  
  r1 <- myFNR(as.numeric(tf1B),tf.vector)   # newly programmed functions for FPR and FNR
  r2 <- myFPR(as.numeric(tf1B),tf.vector)
  mat.STEP[i,] <- c(r1,r2)
}
rates1B.STEP <- colMeans(mat.STEP)
names(rates1B.STEP) <- c("FNR","FPR")   # resulting values to compare with our results


### Comparison of sPCA and STEP based on ranks ---------------------------------

# for search, we use complete vector of important PLR, but to compare with STEP, we use the one without duplicities as well

loads1 <- rownames(res1B[[1]][[1]]$loadings)   # in all models, D = 10 (and scenarios A,B,C), loadings have the same names

# remove "duplicite" PLR in loads1 -> zatim takto nesikovne, nez me napadne sofistikovanejsi cesta
names1 <- loads1[c(1:9,11:18,21:27,31:36,41:45,51:54,61:63,71,72,81)] 




# Choose only names of important pairwise logratios:
important1B <- tf1B[tf1B==T]
important1B.half <- important1B[names(important1B) %in% names1]  # important PLR without duplicities

resulting.list1B <- list()
for(i in 1:R){
  resulting.list1B[[i]] <- RateImportant(models = res1B[[i]], loads = loads1, names = names1, important = important1B.half, D = D1, a = a)
}


resulting.list1B.STEP <- list()
for(i in 1:R){
  resulting.list1B.STEP[[i]] <- RateImportant.G(models = res1B.STEP[[i]], modelsSPCA = res1B[[i]], important = important1B)  # , tot.var = totvar1B[i]
}


cumsum.ideal.1B <- cumsum(1:(D1-1))  # to calculate rank.method minus ideal.rank

tabCumSums.1B <- matrix(0, nrow = R, ncol = D1-1)   # cumulative sum of rank
tabSums.1B <- matrix(0, nrow = R, ncol = D1-1)      # cumulative sum of % of important PLR
tab.Importance.1B <- matrix(0, nrow = R, ncol = D1-1)  # 1 = PLR was important in the run, 0 = PLR was unimportant in the run

for(i in 1:R){
  tabCumSums.1B[i,] <- resulting.list1B[[i]]$cumsum_rank_both
}

for(i in 1:R){
  tabSums.1B[i,] <- resulting.list1B[[i]]$cummulative_sum
}

for(i in 1:R){
  tab.Importance.1B[i,] <- resulting.list1B[[i]]$Is_important
}




tabCumSums.1B.STEP <- matrix(0, nrow = R, ncol = D1-1)
tabSums.1B.STEP <- matrix(0, nrow = R, ncol = D1-1)      
tab.Importance.1B.STEP <- matrix(0, nrow = R, ncol = D1-1) 


for(i in 1:R){
  tabCumSums.1B.STEP[i,] <- resulting.list1B.STEP[[i]]$cumsum_rank_both
}

for(i in 1:R){
  tabSums.1B.STEP[i,] <- resulting.list1B.STEP[[i]]$`cummulative sum`
}

for(i in 1:R){
  tab.Importance.1B.STEP[i,] <- resulting.list1B.STEP[[i]]$Is_important
}

# subtract cumsum.ideal.1A from each row of a matrix tabCumSums.1A
difference.1B <- tabCumSums.1B - rep(cumsum.ideal.1B, each = nrow(tabCumSums.1B))
difference.1B.STEP <- tabCumSums.1B.STEP - rep(cumsum.ideal.1B, each = nrow(tabCumSums.1B))


### Visualisation --------------------------------------------------------------


# a) Fig2 c,d: comparison of sPCA and STEP
plot_line_step(input.sPCA = tabSums.1B, input.STEP = tabSums.1B.STEP, input.importantPLR = important1B.half,D = D1)

plot_boxes(input.sPCA = difference.1B, input.STEP = difference.1B.STEP,D = D1)


# b) Fig4 a,b: development of sparsity - percentage of zero PLRs and expl. variance, FNR and FPR

plot_line(input=tab1B,a=a)
plot_rates(input = rates1B,a=a)
























