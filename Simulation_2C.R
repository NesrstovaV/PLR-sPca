##################################################
### Scenario C: D = 20
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
alpha_max <- 0.212095089 # specify max value of tuning parameter
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
# 4 balances relevant, 15 balances noise

# SBP codes:
codes2C <- rbind(c(rep(1,4),rep(-1,16)), c(rep(1,3),-1,rep(0,16)), c(rep(1,2),-1,rep(0,17)), 
                 c(1,-1,rep(0,18)), c(rep(0,4),1,rep(-1,15)),c(rep(0,5),1,rep(-1,14)),
                 c(rep(0,6),1,rep(-1,13)), c(rep(0,7),1,rep(-1,12)), c(rep(0,8),1,rep(-1,11)),
                 c(rep(0,9),1,rep(-1,10)), c(rep(0,10),1,rep(-1,9)), c(rep(0,11),1,rep(-1,8)),
                 c(rep(0,12),1,rep(-1,7)), c(rep(0,13),1,rep(-1,6)), c(rep(0,14),1,rep(-1,5)),
                 c(rep(0,15),1,rep(-1,4)), c(rep(0,16),1,rep(-1,3)), c(rep(0,17),1,rep(-1,2)),
                 c(rep(0,18),1,-1)
)
codes2C

# Relevant balances: 1-4
mu2C <- rep(0,4)
C2C <- diag(1,ncol = 4,nrow = 4)   
for(i in 1:nrow(C2C)){
  for(j in 1:ncol(C2C)){
    if(i!=j){
      C2C[i,j] <- 0.7
    }
  }
}
eigen(C2C)$values



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
# Approximate time (Windows laptop): cca 12 min for spca_models(); step_models() cca 10min

res2C <- list()
set.seed(1234)
tic()
for(i in 1:R){
  res2C[[i]] <- spca_models(n=n,D=D2,mu=mu2C,C=C2C,codes=codes2C,a=a,nrel=4,nnoise=15)
  print(i)
}
toc() 


variances2C <- matrix(NA,nrow=R,ncol=length(a))
for(i in 1:R){
  for(j in 1:length(a)){
    variances2C[i,j] <- sum(res2C[[i]][[j]]$expl.var[1:2])
  }
} 

means2C <- colMeans(variances2C)

zeros2C <-  matrix(NA,nrow=R,ncol=length(a))
for(i in 1:R){
  for(j in 1:length(a)){
    zeros2C[i,j] <- res2C[[i]][[j]]$`number of zero logratios` 
  }
}

mean2Czeros <- colMeans(zeros2C)



# (average) percentage of zero logratios (counted from mean1Azeros)
nlogr2C <- D2*(D2-1)
prop2C <- (mean2Czeros/nlogr2C)*100  

tab2C <- as.data.frame(cbind(a,means2C,mean2Czeros,prop2C))
colnames(tab2C) <- c("alpha", "var", "ZeroLogr","Percentage")


### FPR and FNR ----------------------------------------------------------------
logr.names2C <- rownames(res2C[[1]][[1]]$loadings)


# True Negative = these logratios should have zero in corresponding loading
vars2C <- c("x5","x6","x7","x8","x9","x10","x11","x12","x13","x14","x15","x16","x17","x18","x19","x20")
combs2C<-as.data.frame(t(combn(vars2C, 2)))

vec1 <- c()
for(i in 1:nrow(combs2C)){
  vec1[i] <- paste0("ln(",combs2C[i,1],"/",combs2C[i,2],")")
}

vec2<- c()
for(i in 1:nrow(combs2C)){
  vec2[i] <- paste0("ln(",combs2C[i,2],"/",combs2C[i,1],")")
}

tn2C <- c(vec1,vec2)


tf2C <- rep(T,length(logr.names2C))  # T/F vector; logratios for which F holds will be added
names(tf2C) <- logr.names2C
tf2C[tn2C] <- F

rates2C <- matrix(NA,nrow=length(a),ncol = 2)   # average over simulaton runs; final matrix of rates
colnames(rates2C) <- c("FNR","FPR")

mat <- matrix(NA,nrow=R,ncol=2)  # auxiliary matrix
for(j in 1:length(a)){
  for(i in 1:R){
    #v <- ifelse(res2C[[i]][[j]]$loadings[,1]!=0 | res2C[[i]][[j]]$loadings[,2]!=0,1,0)
    v <- ifelse(res2C[[i]][[j]]$loadings[,1]==0 & res2C[[i]][[j]]$loadings[,2]==0,0,1)
    r1 <- myFNR(as.numeric(tf2C),v)  # newly programmed functions for FPR and FNR
    r2 <- myFPR(as.numeric(tf2C),v)  
    mat[i,] <- c(r1,r2)
  }
  rates2C[j,] <- colMeans(mat)
  mat <- matrix(NA,nrow=R,ncol=2)  # auxiliary matrix
}



### Results for STEP -----------------------------------------------------------

res2C.STEP <- list()
set.seed(1234)
tic()
for(i in 1:R){
  res2C.STEP[[i]] <- step_models(n=n,D=D2,mu=mu2C,C=C2C,codes=codes2C,nrel=4,nnoise=15)
  print(i)
}
toc()

# FPR and FNR

tn2C;tf2C  # these vectors hold -> results of STEP will be compared with these

rates2C.STEP <- c()  # only two values as it did not run for each alpha 

mat.STEP <- matrix(NA,nrow=R,ncol=2)  # auxiliary matrix
for(i in 1:R){
  ratios <- res2C.STEP[[i]]$names   # save pairwise logratios chosen with STEP
  
  # because STEP choses just A/B (not B/A then), we need to "create" these logratios
  # in our analysis, we consider both A/B and B/A, so to make FNR and FPR comparable
  spl <- strsplit(res2C.STEP[[i]]$names, split = "/")  # split parts of logratios
  df.aux <- data.frame(spl)
  
  ratios2 <- c()
  for(j in 1:ncol(df.aux)){
    ratios2[j] <- paste0("ln(",df.aux[2,j],"/",df.aux[1,j],")")
  }
  v <- c(paste0("ln(",res2C.STEP[[i]]$names,")"), ratios2)  # combine ratios and ratios2
  
  tf.vector <- rep(0,length(tf2C))
  names(tf.vector) <- names(tf2C)
  tf.vector[v] <- 1    # 1 = logratios chosen by STEP, 0 = not chosen
  
  r1 <- myFNR(as.numeric(tf2C),tf.vector)   # newly programmed functions for FPR and FNR
  r2 <- myFPR(as.numeric(tf2C),tf.vector)
  mat.STEP[i,] <- c(r1,r2)
}
rates2C.STEP <- colMeans(mat.STEP)
names(rates2C.STEP) <- c("FNR","FPR")   # resulting values to compare with our results






### Comparison of sPCA and STEP based on ranks ---------------------------------

loads2 <- rownames(res2C[[1]][[1]]$loadings)   # the same holds for D = 20

# remove "duplicite" PLR in loads2
list.logrs <- str_split(loads2, "[ln(/)]")
df.logrs <- matrix(NA, ncol = 2, nrow = length(list.logrs))
for(i in 1:nrow(df.logrs)){
  df.logrs[i,] <- list.logrs[[i]][c(4,5)]
}

df.logrs.half <- as.data.frame(df.logrs[!duplicated(lapply(as.data.frame(t(df.logrs), stringsAsFactors=FALSE), sort)),])

names2 <- paste0("ln(",df.logrs.half[,1],"/",df.logrs.half[,2],")")  # PLR withour duplicites


important2C <- tf2C[tf2C==T]
important2C.half <- important2C[names(important2C) %in% names2]  # IMPORTANT PLR without duplicities

resulting.list2C <- list()
for(i in 1:R){
  resulting.list2C[[i]] <- RateImportant(models = res2C[[i]], loads = loads2, names = names2, important = important2C.half, D = D2, a = a)
}

resulting.list2C.STEP <- list()
for(i in 1:R){
  resulting.list2C.STEP[[i]] <- RateImportant.G(models = res2C.STEP[[i]], modelsSPCA = res2C[[i]], important = important2C)  # , tot.var = totvar2C[i]
}

cumsum.ideal.2C <- cumsum(1:(D2-1)) # to calculate rank.method minus ideal.rank


tabCumSums.2C <- matrix(0, nrow = R, ncol = D2-1)   # cumulative sum of rank
tabSums.2C <- matrix(0, nrow = R, ncol = D2-1)      # cumulative sum of % of important PLR
tab.Importance.2C <- matrix(0, nrow = R, ncol = D2-1)  # 1 = PLR was important in the run, 0 = PLR was unimportant in the run

for(i in 1:R){
  tabCumSums.2C[i,] <- resulting.list2C[[i]]$cumsum_rank_both
}

for(i in 1:R){
  tabSums.2C[i,] <- resulting.list2C[[i]]$cummulative_sum
}

for(i in 1:R){
  tab.Importance.2C[i,] <- resulting.list2C[[i]]$Is_important
}

tabCumSums.2C.STEP <- matrix(0, nrow = R, ncol = D2-1)
tabSums.2C.STEP <- matrix(0, nrow = R, ncol = D2-1)      
tab.Importance.2C.STEP <- matrix(0, nrow = R, ncol = D2-1)  

for(i in 1:R){
  tabCumSums.2C.STEP[i,] <- resulting.list2C.STEP[[i]]$cumsum_rank_both
}

for(i in 1:R){
  tabSums.2C.STEP[i,] <- resulting.list2C.STEP[[i]]$`cummulative sum`
}

for(i in 1:R){
  tab.Importance.2C.STEP[i,] <- resulting.list2C.STEP[[i]]$Is_important
}


# subtract cumsum.ideal.1A from each row of a matrix tabCumSums.1A
difference.2C <- tabCumSums.2C - rep(cumsum.ideal.2C, each = nrow(tabCumSums.2C))
difference.2C.STEP <- tabCumSums.2C.STEP - rep(cumsum.ideal.2C, each = nrow(tabCumSums.2C))



### Visualisation --------------------------------------------------------------

# a) Fig12 e,f: comparison of sPCA and STEP

plot_line_step(input.sPCA = tabSums.2C, input.STEP = tabSums.2C.STEP, input.importantPLR = important2C.half, D = D2)

plot_boxes(input.sPCA = difference.2C, input.STEP = difference.2C.STEP,D = D2)


# b) Fig13 e,f: development of sparsity - percentage of zero PLRs and expl. variance, FNR and FPR

plot_line(input=tab2C,a=a)
plot_rates(input = rates2C,a=a)










