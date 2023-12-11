##################################################
### Scenario A: D = 20
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
library(stringr)


### Source functions -----------------------------------------------------------

# Key function: a function making a matrix of D(D-1) logratios and calculating sparse PCA
source("spca_function.R")


# Auxiliary functions for calculations
source("Auxiliary_functions_spca.R")

### Intro ----------------------------------------------------------------------

alpha_max <- 0.184206997 # specify max value of tuning parameter
alpha_nbr <- 50 # specify number of tuning parameters
alpha_ratio <- 1000 # specify ratio of largest to smallest tuning parameter
alpha_grid <- c(exp(seq(from = log(alpha_max), to = log(alpha_max/alpha_ratio), length = alpha_nbr)), 0)

# a vector of sparsity parameters alpha
a <- sort(alpha_grid,decreasing=F)        # zero (-> fully dense model) included


R <- 100                       # number of simulation runs
n <- 100
D2 <- 20


### Data simulation parameters -------------------------------------------------

# covariance matrix: diagonal elements = a, nondiagonal elements = b
# 10 balances relevant,  9 balances noise

# SBP codes:
codes2A <- rbind(c(rep(1,10),rep(-1,10)), c(rep(1,9),-1,rep(0,10)), c(rep(1,8),-1,rep(0,11)), 
                 c(rep(1,7),-1,rep(0,12)), c(rep(1,6),-1,rep(0,13)), c(rep(1,5),-1,rep(0,14)),
                 c(rep(1,4),-1,rep(0,15)), c(rep(1,3),-1,rep(0,16)), c(rep(1,2),-1,rep(0,17)),
                 c(1,-1,rep(0,18)), c(rep(0,10),1,rep(-1,9)), c(rep(0,11),1,rep(-1,8)),
                 c(rep(0,12),1,rep(-1,7)), c(rep(0,13),1,rep(-1,6)), c(rep(0,14),1,rep(-1,5)),
                 c(rep(0,15),1,rep(-1,4)), c(rep(0,16),1,rep(-1,3)), c(rep(0,17),1,rep(-1,2)),
                 c(rep(0,18),1,-1))


# Relevant balances: 1-10
mu2A <- rep(0,10)
C2A <- diag(1,ncol = D2/2,nrow = D2/2)   
for(i in 1:nrow(C2A)){
  for(j in 1:ncol(C2A)){
    if(i!=j){
      C2A[i,j] <- 0.7
    }
  }
}
eigen(C2A)$values


### Simulation runs ------------------------------------------------------------

# MacBook Pro
# Procesor: M1
# RAM: 16GB

# Windows laptop
# Processor: Intel(R) Core(TM) i7-4810MQ
# Cores: 4
# Logical processors: 8
# Freq: 2.80GHz


# Approximate time (MacBook Pro): cca 21 min for spca_models(); step_models() runs faster 
# Approximate time (Windows laptop): cca 30 min for spca_models(); step_models() cca 10min

res2A <- list()
set.seed(1234)
tic()
for(i in 1:R){
  res2A[[i]] <- spca_models(n=n,D=D2,mu=mu2A,C=C2A,codes=codes2A,a=a,nrel=10,nnoise=9)
  print(i)
}
toc()  

variances2A <- matrix(NA,nrow=R,ncol=length(a))
for(i in 1:R){
  for(j in 1:length(a)){
    variances2A[i,j] <- sum(res2A[[i]][[j]]$expl.var)
  }
}
means2A <- colMeans(variances2A)

zeros2A <-  matrix(NA,nrow=R,ncol=length(a))
for(i in 1:R){
  for(j in 1:length(a)){
    zeros2A[i,j] <- res2A[[i]][[j]]$`number of zero logratios` 
  }
}
mean2Azeros <- colMeans(zeros2A)


# (average) percentage of zero logratios (counted from mean1Azeros)
nlogr2A <- D2*(D2-1)
prop2A <- (mean2Azeros/nlogr2A)*100  


tab2A <- as.data.frame(cbind(a,means2A,mean2Azeros,prop2A))
colnames(tab2A) <- c("alpha", "var", "ZeroLogr","Percentage")



### FPR and FNR ----------------------------------------------------------------

logr.names2A <- rownames(res2A[[1]][[1]]$loadings)

# True Negative = these logratios should have zero in corresponding loading
vars2A <- c("x11","x12","x13","x14","x15","x16","x17","x18","x19","x20")
combs2A<-as.data.frame(t(combn(vars2A, 2)))

vec1 <- c()
for(i in 1:nrow(combs2A)){
  vec1[i] <- paste0("ln(",combs2A[i,1],"/",combs2A[i,2],")")
}

vec2<- c()
for(i in 1:nrow(combs2A)){
  vec2[i] <- paste0("ln(",combs2A[i,2],"/",combs2A[i,1],")")
}

tn2A <- c(vec1,vec2)

tf2A <- rep(T,length(logr.names2A))  # T/F vector; logratios for which F holds will be added
names(tf2A) <- logr.names2A
tf2A[tn2A] <- F

rates2A <- matrix(NA,nrow=length(a),ncol = 2)   # average over simulaton runs; final matrix of rates
colnames(rates2A) <- c("FNR","FPR")

mat <- matrix(NA,nrow=R,ncol=2)  # auxiliary matrix
for(j in 1:length(a)){
  for(i in 1:R){
    v <- ifelse(res2A[[i]][[j]]$loadings[,1]==0 & res2A[[i]][[j]]$loadings[,2]==0,0,1)
    r1 <- myFNR(as.numeric(tf2A),v)   # newly programmed functions for FPR and FNR
    r2 <- myFPR(as.numeric(tf2A),v)
    mat[i,] <- c(r1,r2)
  }
  rates2A[j,] <- colMeans(mat)
  mat <- matrix(NA,nrow=R,ncol=2)  # auxiliary matrix
}

### Results for STEP -----------------------------------------------------------

res2A.STEP <- list()
set.seed(1234)
tic()
for(i in 1:R){
  res2A.STEP[[i]] <- step_models(n=n,D=D2,mu=mu2A,C=C2A,codes=codes2A,nrel=10,nnoise=9)
  print(i)
}
toc()  # +- 13 min on my DELL laptop

# FPR and FNR

tn2A;tf2A  # these vectors hold -> results of STEP will be compared with these

rates2A.STEP <- c()  # only two values as it did not run for each alpha 

mat.STEP <- matrix(NA,nrow=R,ncol=2)  # auxiliary matrix
for(i in 1:R){
  ratios <- res2A.STEP[[i]]$names   # save pairwise logratios chosen with STEP
  
  # because STEP choses just A/B (not B/A then), we need to "create" these logratios
  # in our analysis, we consider both A/B and B/A, so to make FNR and FPR comparable
  spl <- strsplit(res2A.STEP[[i]]$names, split = "/")  # split parts of logratios
  df.aux <- data.frame(spl)
  
  ratios2 <- c()
  for(j in 1:ncol(df.aux)){
    ratios2[j] <- paste0("ln(",df.aux[2,j],"/",df.aux[1,j],")")
  }
  v <- c(paste0("ln(",res2A.STEP[[i]]$names,")"), ratios2)  # combine ratios and ratios2
  
  tf.vector <- rep(0,length(tf2A))
  names(tf.vector) <- names(tf2A)
  tf.vector[v] <- 1    # 1 = logratios chosen by STEP, 0 = not chosen
  
  r1 <- myFNR(as.numeric(tf2A),tf.vector)   # newly programmed functions for FPR and FNR
  r2 <- myFPR(as.numeric(tf2A),tf.vector)
  mat.STEP[i,] <- c(r1,r2)
}
rates2A.STEP <- colMeans(mat.STEP)
names(rates2A.STEP) <- c("FNR","FPR")   # resulting values to compare with our results



### Comparison of sPCA and STEP based on ranks ---------------------------------

loads2 <- rownames(res2A[[1]][[1]]$loadings)   # the same holds for D = 20

# remove "duplicite" PLR in loads2
list.logrs <- str_split(loads2, "[ln(/)]")
df.logrs <- matrix(NA, ncol = 2, nrow = length(list.logrs))
for(i in 1:nrow(df.logrs)){
  df.logrs[i,] <- list.logrs[[i]][c(4,5)]
}

df.logrs.half <- as.data.frame(df.logrs[!duplicated(lapply(as.data.frame(t(df.logrs), stringsAsFactors=FALSE), sort)),])

names2 <- paste0("ln(",df.logrs.half[,1],"/",df.logrs.half[,2],")")  # PLR withour duplicites

loads2;tf2A

important2A <- tf2A[tf2A==T]
important2A.half <- important2A[names(important2A) %in% names2]  # IMPORTANT PLR without duplicities

resulting.list2A <- list()
for(i in 1:R){
  resulting.list2A[[i]] <- RateImportant(models = res2A[[i]], loads = loads2, names = names2, important = important2A.half, D = D2, a = a)
}


resulting.list2A.STEP <- list()
for(i in 1:R){
  resulting.list2A.STEP[[i]] <- RateImportant.G(models = res2A.STEP[[i]], modelsSPCA = res2A[[i]], important = important2A) # , tot.var = totvar2A[i]
}


cumsum.ideal.2A <- cumsum(1:(D2-1)) # to calculate rank.method minus ideal.rank



tabCumSums.2A <- matrix(0, nrow = R, ncol = D2-1)   # cumulative sum of rank
tabSums.2A <- matrix(0, nrow = R, ncol = D2-1)      # cumulative sum of % of important PLR
tab.Importance.2A <- matrix(0, nrow = R, ncol = D2-1)  # 1 = PLR was important in the run, 0 = PLR was unimportant in the run

for(i in 1:R){
  tabCumSums.2A[i,] <- resulting.list2A[[i]]$cumsum_rank_both
}

for(i in 1:R){
  tabSums.2A[i,] <- resulting.list2A[[i]]$cummulative_sum
}

for(i in 1:R){
  tab.Importance.2A[i,] <- resulting.list2A[[i]]$Is_important
}


tabCumSums.2A.STEP <- matrix(0, nrow = R, ncol = D2-1)
tabSums.2A.STEP <- matrix(0, nrow = R, ncol = D2-1)      
tab.Importance.2A.STEP <- matrix(0, nrow = R, ncol = D2-1)  

for(i in 1:R){
  tabCumSums.2A.STEP[i,] <- resulting.list2A.STEP[[i]]$cumsum_rank_both
}

for(i in 1:R){
  tabSums.2A.STEP[i,] <- resulting.list2A.STEP[[i]]$`cummulative sum`
}


for(i in 1:R){
  tab.Importance.2A.STEP[i,] <- resulting.list2A.STEP[[i]]$Is_important
}

# subtract cumsum.ideal.1A from each row of a matrix tabCumSums.1A
difference.2A <- tabCumSums.2A - rep(cumsum.ideal.2A, each = nrow(tabCumSums.2A))
difference.2A.STEP <- tabCumSums.2A.STEP - rep(cumsum.ideal.2A, each = nrow(tabCumSums.2A))

### Visualisation --------------------------------------------------------------

# a) Fig12 a,b: comparison of sPCA and STEP

plot_line_step(input.sPCA = tabSums.2A, input.STEP = tabSums.2A.STEP, input.importantPLR = important2A.half, D = D2)

plot_boxes(input.sPCA = difference.2A, input.STEP = difference.2A.STEP,D = D2)


# b) Fig13 a,b: development of sparsity - percentage of zero PLRs and expl. variance, FNR and FPR

plot_line(input=tab2A,a=a)
plot_rates(input = rates2A,a=a)




