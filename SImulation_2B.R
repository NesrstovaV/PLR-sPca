##################################################
### Scenario B: D = 20
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


alpha_max <- 0.159985872 # specify max value of tuning parameter
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
# 15 balances relevant, 4 balances noise

# SBP codes:
codes2B <- rbind(c(rep(1,15),rep(-1,5)), c(rep(1,14),-1,rep(0,5)), c(rep(1,13),-1,rep(0,6)), 
                 c(rep(1,12),-1,rep(0,7)), c(rep(1,11),-1,rep(0,8)), c(rep(1,10),-1,rep(0,9)),
                 c(rep(1,9),-1,rep(0,10)), c(rep(1,8),-1,rep(0,11)), c(rep(1,7),-1,rep(0,12)),
                 c(rep(1,6),-1,rep(0,13)), c(rep(1,5),-1,rep(0,14)), c(rep(1,4),-1,rep(0,15)),
                 c(rep(1,3),-1,rep(0,16)), c(rep(1,2),-1,rep(0,17)), c(1,-1,rep(0,18)),
                 c(rep(0,15),1,rep(-1,4)), c(rep(0,16),1,rep(-1,3)), c(rep(0,17),1,rep(-1,2)),
                 c(rep(0,18),1,-1))

# Relevant balances: 1-15
mu2B <- rep(0,15)
C2B <- diag(1,ncol = D2-5,nrow = D2-5)   
for(i in 1:nrow(C2B)){
  for(j in 1:ncol(C2B)){
    if(i!=j){
      C2B[i,j] <- 0.7
    }
  }
}
eigen(C2B)$values 




### Simulation runs ------------------------------------------------------------

# MacBook Pro
# Procesor: M1
# RAM: 16GB

# Windows laptop
# Processor: Intel(R) Core(TM) i7-4810MQ
# Cores: 4
# Logical processors: 8
# Freq: 2.80GHz


# Approximate time (MacBook Pro): cca X min for spca_models(); step_models() runs faster 
# Approximate time (Windows laptop): cca 30 min for spca_models(); step_models() cca 10min

res2B <- list()
set.seed(1234)
tic()
for(i in 1:R){
  res2B[[i]] <- spca_models(n=n,D=D2,mu=mu2B,C=C2B,codes=codes2B,a=a,nrel=15,nnoise=4)
  print(i)
}
toc()  


variances2B <- matrix(NA,nrow=R,ncol=length(a))
for(i in 1:R){
  for(j in 1:length(a)){
    variances2B[i,j] <- sum(res2B[[i]][[j]]$expl.var[1:2])
  }
}

means2B <- colMeans(variances2B)

zeros2B <-  matrix(NA,nrow=R,ncol=length(a))
for(i in 1:R){
  for(j in 1:length(a)){
    zeros2B[i,j] <- res2B[[i]][[j]]$`number of zero logratios` 
  }
}

mean2Bzeros <- colMeans(zeros2B)


# (average) percentage of zero logratios (counted from mean1Azeros)
nlogr2B <- D2*(D2-1)
prop2B <- (mean2Bzeros/nlogr2B)*100  

tab2B <- as.data.frame(cbind(a,means2B,mean2Bzeros,prop2B))
colnames(tab2B) <- c("alpha", "var", "ZeroLogr","Percentage")


### FPR and FNR ----------------------------------------------------------------

logr.names2B <- rownames(res2B[[1]][[1]]$loadings)

# True Negative = these logratios should have zero in corresponding loading
vars2B <- c("x16","x17","x18","x19","x20")
combs2B<-as.data.frame(t(combn(vars2B, 2)))

vec1 <- c()
for(i in 1:nrow(combs2B)){
  vec1[i] <- paste0("ln(",combs2B[i,1],"/",combs2B[i,2],")")
}

vec2<- c()
for(i in 1:nrow(combs2B)){
  vec2[i] <- paste0("ln(",combs2B[i,2],"/",combs2B[i,1],")")
}

tn2B <- c(vec1,vec2)

tf2B <- rep(T,length(logr.names2B))  # T/F vector; logratios for which F holds will be added
names(tf2B) <- logr.names2B
tf2B[tn2B] <- F

rates2B <- matrix(NA,nrow=length(a),ncol = 2)   # average over simulaton runs; final matrix of rates
colnames(rates2B) <- c("FNR","FPR")

mat <- matrix(NA,nrow=R,ncol=2)  # auxiliary matrix
for(j in 1:length(a)){
  for(i in 1:R){
    #v <- ifelse(res2B[[i]][[j]]$loadings[,1]!=0 | res2B[[i]][[j]]$loadings[,2]!=0,1,0)
    v <- ifelse(res2B[[i]][[j]]$loadings[,1]==0 & res2B[[i]][[j]]$loadings[,2]==0,0,1)
    r1 <- myFNR(as.numeric(tf2B),v)  # newly programmed functions for FPR and FNR
    r2 <- myFPR(as.numeric(tf2B),v)  
    mat[i,] <- c(r1,r2)
  }
  rates2B[j,] <- colMeans(mat)
  mat <- matrix(NA,nrow=R,ncol=2)  # auxiliary matrix
}




### Results for STEP -----------------------------------------------------------

res2B.STEP <- list()
set.seed(1234)
tic()
for(i in 1:R){
  res2B.STEP[[i]] <- step_models(n=n,D=D2,mu=mu2B,C=C2B,codes=codes2B,nrel=15,nnoise=4)
  print(i)
}
toc()  

# FPR and FNR

tn2B;tf2B  # these vectors hold -> results of STEP will be compared with these

rates2B.STEP <- c()  # only two values as it did not run for each alpha 

mat.STEP <- matrix(NA,nrow=R,ncol=2)  # auxiliary matrix
for(i in 1:R){
  ratios <- res2B.STEP[[i]]$names   # save pairwise logratios chosen with STEP
  
  # because STEP choses just A/B (not B/A then), we need to "create" these logratios
  # in our analysis, we consider both A/B and B/A, so to make FNR and FPR comparable
  spl <- strsplit(res2B.STEP[[i]]$names, split = "/")  # split parts of logratios
  df.aux <- data.frame(spl)
  
  ratios2 <- c()
  for(j in 1:ncol(df.aux)){
    ratios2[j] <- paste0("ln(",df.aux[2,j],"/",df.aux[1,j],")")
  }
  v <- c(paste0("ln(",res2B.STEP[[i]]$names,")"), ratios2)  # combine ratios and ratios2
  
  tf.vector <- rep(0,length(tf2B))
  names(tf.vector) <- names(tf2B)
  tf.vector[v] <- 1    # 1 = logratios chosen by STEP, 0 = not chosen
  
  r1 <- myFNR(as.numeric(tf2B),tf.vector)   # newly programmed functions for FPR and FNR
  r2 <- myFPR(as.numeric(tf2B),tf.vector)
  mat.STEP[i,] <- c(r1,r2)
}
rates2B.STEP <- colMeans(mat.STEP)
names(rates2B.STEP) <- c("FNR","FPR")   # resulting values to compare with our results

### Comparison of sPCA and STEP based on ranks ---------------------------------

loads2 <- rownames(res2B[[1]][[1]]$loadings)   # the same holds for D = 20

# remove "duplicite" PLR in loads2
list.logrs <- str_split(loads2, "[ln(/)]")
df.logrs <- matrix(NA, ncol = 2, nrow = length(list.logrs))
for(i in 1:nrow(df.logrs)){
  df.logrs[i,] <- list.logrs[[i]][c(4,5)]
}

df.logrs.half <- as.data.frame(df.logrs[!duplicated(lapply(as.data.frame(t(df.logrs), stringsAsFactors=FALSE), sort)),])

names2 <- paste0("ln(",df.logrs.half[,1],"/",df.logrs.half[,2],")")  # PLR withour duplicites


important2B <- tf2B[tf2B==T]
important2B.half <- important2B[names(important2B) %in% names2]  # IMPORTANT PLR without duplicities


resulting.list2B <- list()
for(i in 1:R){
  resulting.list2B[[i]] <- RateImportant(models = res2B[[i]], loads = loads2, names = names2, important = important2B.half, D = D2, a = a)
}



resulting.list2B.STEP <- list()
for(i in 1:R){
  resulting.list2B.STEP[[i]] <- RateImportant.G(models = res2B.STEP[[i]], modelsSPCA = res2B[[i]], important = important2B)  # , tot.var = totvar2B[i]
}

cumsum.ideal.2B <- cumsum(1:(D2-1)) # to calculate rank.method minus ideal.rank


tabCumSums.2B <- matrix(0, nrow = R, ncol = D2-1)   # cumulative sum of rank
tabSums.2B <- matrix(0, nrow = R, ncol = D2-1)      # cumulative sum of % of important PLR
tab.Importance.2B <- matrix(0, nrow = R, ncol = D2-1)  # 1 = PLR was important in the run, 0 = PLR was unimportant in the run

for(i in 1:R){
  tabCumSums.2B[i,] <- resulting.list2B[[i]]$cumsum_rank_both
}

for(i in 1:R){
  tabSums.2B[i,] <- resulting.list2B[[i]]$cummulative_sum
}

for(i in 1:R){
  tab.Importance.2B[i,] <- resulting.list2B[[i]]$Is_important
}


tabCumSums.2B.STEP <- matrix(0, nrow = R, ncol = D2-1)
tabSums.2B.STEP <- matrix(0, nrow = R, ncol = D2-1)      
tab.Importance.2B.STEP <- matrix(0, nrow = R, ncol = D2-1) 

for(i in 1:R){
  tabCumSums.2B.STEP[i,] <- resulting.list2B.STEP[[i]]$cumsum_rank_both
}

for(i in 1:R){
  tabSums.2B.STEP[i,] <- resulting.list2B.STEP[[i]]$`cummulative sum`
}

for(i in 1:R){
  tab.Importance.2B.STEP[i,] <- resulting.list2B.STEP[[i]]$Is_important
}


# subtract cumsum.ideal.2B from each row of a matrix tabCumSums.2B
difference.2B <- tabCumSums.2B - rep(cumsum.ideal.2B, each = nrow(tabCumSums.2B))
difference.2B.STEP <- tabCumSums.2B.STEP - rep(cumsum.ideal.2B, each = nrow(tabCumSums.2B))



### Visualisation --------------------------------------------------------------


# a) Fig12 c,d: comparison of sPCA and STEP

plot_line_step(input.sPCA = tabSums.2B, input.STEP = tabSums.2B.STEP, input.importantPLR = important2B.half, D = D2)

plot_boxes(input.sPCA = difference.2B, input.STEP = difference.2B.STEP,D = D2)


# b) Fig13 c,d: development of sparsity - percentage of zero PLRs and expl. variance, FNR and FPR

plot_line(input=tab2B,a=a)
plot_rates(input = rates2B,a=a)



