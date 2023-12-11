##################################################
### Scenario C: D = 10
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

alpha_max <- 0.323745754 # specify max value of tuning parameter
alpha_nbr <- 50 # specify number of tuning parameters
alpha_ratio <- 1000 # specify ratio of largest to smallest tuning parameter
alpha_grid <- c(exp(seq(from = log(alpha_max), to = log(alpha_max/alpha_ratio), length = alpha_nbr)), 0)

# a vector of sparsity parameters alpha
a <- sort(alpha_grid,decreasing=F)        # zero (-> fully dense model) included


R <- 100                       # number of simulation runs
n <- 100

D1 <- 10

### Data simulation parameters--------------------------------------------------
# covariance matrix: diagonal elements = a, nondiagonal elements = b
# 2 balances relevant, 7 balances noise

# SBP codes:
codes1C <- rbind(c(1,1,rep(-1,8)), c(1,-1,rep(0,8)), c(0,0,1,rep(-1,7)), 
                 c(rep(0,3),1,rep(-1,6)), c(rep(0,4),1,rep(-1,5)),c(rep(0,5),1,rep(-1,4)),c(rep(0,6),1,rep(-1,3)),
                 c(rep(0,7),1,-1,-1),c(rep(0,8),1,-1))
codes1C

# Relevant balances: 1-2
mu1C <- rep(0,2)
C1C <- diag(1,ncol = D1-8,nrow = D1-8)  # first version of C1 was: 4 = diag, 1.5 = nondiag
C1C[1,2] <- 0.7
C1C[2,1] <- 0.7
C1C
eigen(C1C)$values

### Simulation runs ------------------------------------------------------------

# MacBook Pro
# Procesor: M1
# RAM: 16GB

# Windows laptop
# Processor: Intel(R) Core(TM) i7-4810MQ
# Cores: 4
# Logical processors: 8
# Freq: 2.80GHz

# Approximate time (MacBook Pro): cca min for spca_models(); step_models() runs faster 
# Approximate time (Windows laptop): cca 2.5min for spca_models(); step_models() runs faster 

res1C <- list()
set.seed(1234)
tic()
for(i in 1:R){
  res1C[[i]] <- spca_models(n=n,D=D1,mu=mu1C,C=C1C,codes=codes1C,a=a,nrel=2,nnoise=7)
  print(i)
}
toc() 


variances1C <- matrix(NA,nrow=R,ncol=length(a))
for(i in 1:R){
  for(j in 1:length(a)){
    variances1C[i,j] <- sum(res1C[[i]][[j]]$expl.var)
  }
}
means1C <- colMeans(variances1C)

zeros1C <-  matrix(NA,nrow=R,ncol=length(a))
for(i in 1:R){
  for(j in 1:length(a)){
    zeros1C[i,j] <- res1C[[i]][[j]]$`number of zero logratios` 
  }
}
mean1Czeros <- colMeans(zeros1C)



# (average) percentage of zero logratios (counted from mean1Azeros)
nlogr1C <- D1*(D1-1)   # same as nlogr1A, nlogr1B
prop1C <- (mean1Czeros/nlogr1C)*100  

tab1C <- as.data.frame(cbind(a,means1C,mean1Czeros,prop1C))
colnames(tab1C) <- c("alpha", "var", "ZeroLogr","Percentage")


### FPR and FNR ----------------------------------------------------------------

# True Negative = these logratios should have zero in corresponding loading

logr.names1C <- rownames(res1C[[1]][[1]]$loadings)

vars1C <- c("x3","x4","x5","x6","x7","x8","x9","x10")
combs1C<-as.data.frame(t(combn(vars1C, 2)))

vec1 <- c()
for(i in 1:nrow(combs1C)){
  vec1[i] <- paste0("ln(",combs1C[i,1],"/",combs1C[i,2],")")
}

vec2<- c()
for(i in 1:nrow(combs1C)){
  vec2[i] <- paste0("ln(",combs1C[i,2],"/",combs1C[i,1],")")
}

tn1C <- c(vec1,vec2)

tf1C <- rep(T,length(logr.names1C))  # T/F vector; logratios for which F holds will be added
names(tf1C) <- logr.names1C
tf1C[tn1C] <- F

rates1C <- matrix(NA,nrow=length(a),ncol = 2)   # average over simulaton runs; final matrix of rates
colnames(rates1C) <- c("FNR","FPR")

mat <- matrix(NA,nrow=R,ncol=2)  # auxiliary matrix
for(j in 1:length(a)){
  for(i in 1:R){
    v <- ifelse(res1C[[i]][[j]]$loadings[,1]==0 & res1C[[i]][[j]]$loadings[,2]==0,0,1)
    r1 <- myFNR(as.numeric(tf1C),v)   # newly programmed functions for FPR and FNR
    r2 <- myFPR(as.numeric(tf1C),v)
    mat[i,] <- c(r1,r2)
  }
  rates1C[j,] <- colMeans(mat)
  mat <- matrix(NA,nrow=R,ncol=2)  # auxiliary matrix
}

### Results for STEP -----------------------------------------------------------

res1C.STEP <- list()
set.seed(1234)
tic()
for(i in 1:R){
  res1C.STEP[[i]] <- step_models(n=n,D=D1,mu=mu1C,C=C1C,codes=codes1C,nrel=2,nnoise=7)
  print(i)
}
toc()

tn1C;tf1C  # these vectors hold -> results of STEP will be compared with these

rates1C.STEP <- c()  # only two values as it did not run for each alpha 

mat.STEP <- matrix(NA,nrow=R,ncol=2)  # auxiliary matrix
for(i in 1:R){
  ratios <- res1C.STEP[[i]]$names   # save pairwise logratios chosen with STEP
  
  # because STEP choses just A/B (not B/A then), we need to "create" these logratios
  # in our analysis, we consider both A/B and B/A, so to make FNR and FPR comparable
  spl <- strsplit(res1C.STEP[[i]]$names, split = "/")  # split parts of logratios
  df.aux <- data.frame(spl)
  
  ratios2 <- c()
  for(j in 1:ncol(df.aux)){
    ratios2[j] <- paste0("ln(",df.aux[2,j],"/",df.aux[1,j],")")
  }
  v <- c(paste0("ln(",res1C.STEP[[i]]$names,")"), ratios2)  # combine ratios and ratios2
  
  tf.vector <- rep(0,length(tf1C))
  names(tf.vector) <- names(tf1C)
  tf.vector[v] <- 1    # 1 = logratios chosen by STEP, 0 = not chosen
  
  r1 <- myFNR(as.numeric(tf1C),tf.vector)   # newly programmed functions for FPR and FNR
  r2 <- myFPR(as.numeric(tf1C),tf.vector)
  mat.STEP[i,] <- c(r1,r2)
}
rates1C.STEP <- colMeans(mat.STEP)
names(rates1C.STEP) <- c("FNR","FPR")   # resulting values to compare with our results



### Comparison of sPCA and STEP based on ranks ---------------------------------

# for search, we use complete vector of important PLR, but to compare with STEP, we use the one without duplicities as well

loads1 <- rownames(res1C[[1]][[1]]$loadings)   # in all models, D = 10 (and scenarios A,B,C), loadings have the same names

# remove "duplicite" PLR in loads1 -> zatim takto nesikovne, nez me napadne sofistikovanejsi cesta
names1 <- loads1[c(1:9,11:18,21:27,31:36,41:45,51:54,61:63,71,72,81)] 


important1C <- tf1C[tf1C==T]
important1C.half <- important1C[names(important1C) %in% names1]  # important PLR without duplicities


resulting.list1C <- list()
for(i in 1:R){
  resulting.list1C[[i]] <- RateImportant(models = res1C[[i]], loads = loads1, names = names1, important = important1C.half, D = D1, a = a)
}


resulting.list1C.STEP <- list()
for(i in 1:R){
  resulting.list1C.STEP[[i]] <- RateImportant.G(models = res1C.STEP[[i]], modelsSPCA = res1C[[i]], important = important1C)  # ,  tot.var = totvar1C[i]
}


cumsum.ideal.1C <- cumsum(1:(D1-1))


tabCumSums.1C <- matrix(0, nrow = R, ncol = D1-1)   # cumulative sum of rank
tabSums.1C <- matrix(0, nrow = R, ncol = D1-1)      # cumulative sum of % of important PLR
tab.Importance.1C <- matrix(0, nrow = R, ncol = D1-1)  # 1 = PLR was important in the run, 0 = PLR was unimportant in the run

for(i in 1:R){
  tabCumSums.1C[i,] <- resulting.list1C[[i]]$cumsum_rank_both
}

for(i in 1:R){
  tabSums.1C[i,] <- resulting.list1C[[i]]$cummulative_sum
}

for(i in 1:R){
  tab.Importance.1C[i,] <- resulting.list1C[[i]]$Is_important
}


tabCumSums.1C.STEP <- matrix(0, nrow = R, ncol = D1-1)
tabSums.1C.STEP <- matrix(0, nrow = R, ncol = D1-1)      
tab.Importance.1C.STEP <- matrix(0, nrow = R, ncol = D1-1) 

for(i in 1:R){
  tabCumSums.1C.STEP[i,] <- resulting.list1C.STEP[[i]]$cumsum_rank_both
}

for(i in 1:R){
  tabSums.1C.STEP[i,] <- resulting.list1C.STEP[[i]]$`cummulative sum`
}

for(i in 1:R){
  tab.Importance.1C.STEP[i,] <- resulting.list1C.STEP[[i]]$Is_important
}


# subtract cumsum.ideal.1A from each row of a matrix tabCumSums.1A
difference.1C <- tabCumSums.1C - rep(cumsum.ideal.1C, each = nrow(tabCumSums.1C))
difference.1C.STEP <- tabCumSums.1C.STEP - rep(cumsum.ideal.1C, each = nrow(tabCumSums.1C))



### Visualisation --------------------------------------------------------------

# a) Fig2 e,f: comparison of sPCA and STEP
plot_line_step(input.sPCA = tabSums.1C, input.STEP = tabSums.1C.STEP, input.importantPLR = important1C.half,D = D1)

plot_boxes(input.sPCA = difference.1C, input.STEP = difference.1C.STEP, D = D1)


# b) Fig5 a,b: development of sparsity - percentage of zero PLRs and expl. variance, FNR and FPR

plot_line(input=tab1C,a=a)
plot_rates(input = rates1C,a=a)
