################################################################################
### AUXILIARY FUNCTIONS FOR SPCA SIMULATIONS
################################################################################

# True positive (TP): number of truly non-sparse logratios
# True negative (TN): number of truly sparse logratios
# False positive (FP): number of logratios identified as non-sparse, which are in fact SPARSE
# False negative (FN): number of logratios identified as sparse, which are in fact NON-SPARSE


###########
## FPR
###########
myFPR <- function(actual,predicted){
  TN_logical <- (actual==0)+(predicted==0)
  TN <- sum(TN_logical==2)
  
  FP_logical <- (actual==0)+(predicted==1)
  FP <- sum(FP_logical==2)
  
  FPR <- FP/(FP+TN)
  return(FPR)
}

###########
## FNR
###########
myFNR <- function(actual,predicted){
  FN_logical <- (actual==1)+(predicted==0)
  FN <- sum(FN_logical==2)
  
  TP_logical <- (actual==1)+(predicted==1)
  TP <- sum(TP_logical==2)
  
  FNR <- FN/(FN+TP)
  return(FNR)
}



#################################################################
## Function (for simulations) to calculate and store spca models
#################################################################
spca_models <- function(n = n, D = D, mu = mu, C = C, codes = codes, a = a, nrel = nrel, nnoise = nnoise){
  # n, D: dimensions of data
  # mu: vector mu
  # C: covariance matrix
  # codes: given SBP
  # a: vector of sparsity parameters alpha
  # nrel = number of relevant balances
  # nnoise = number of noise balances
  
  
  ## STEP 1: data
  
  # Relevant balances
  rel <- mvrnorm(n, mu, Sigma = C)
  
  # Noise balances: 6-9
  noise <- matrix(NA,ncol = nnoise, nrow = n)
  for(j in 1:nnoise){
    noise[,j] <- runif(n,-2,2)
  }
  
  # all balances into one data frame:
  bal <- cbind(rel,noise)
  
  # from ilr to clr
  V <- gsi.buildilrBase(t(codes)) 
  
  Bal_clr <- bal %*% t(V)  
  
  # from clr to composition
  X <-  as.matrix(as.data.frame(acomp(exp(Bal_clr)))) 
  colnames(X) <- paste0(rep("x",D),1:D)  # generic colnames of matrix X
  
  ## STEP 2: sparse PCA models
  models <- list()
  for(k in 1:length(a)){
    models[[k]] <- spca_funkce(X,alpha = a[k], k = 2, draw = F) # calculates PC1 and PC2
  }
  return(models)
}


#######################################################################################
## Function to simulate data and do stepwise procedure using STEP function 
#######################################################################################
step_models <- function(n = n, D = D, mu = mu, C = C, codes = codes, nrel = nrel, nnoise = nnoise){
  # n, D: dimensions of data
  # mu: vector mu
  # C: covariance matrix
  # codes: given SBP
  # a: vector of sparsity parameters alpha
  # nrel = number of relevant balances
  # nnoise = number of noise balances
  
  ## STEP 1: data
  
  # Relevant balances
  rel <- mvrnorm(n, mu, Sigma = C)
  
  # Noise balances: 6-9
  noise <- matrix(NA,ncol = nnoise, nrow = n)
  for(j in 1:nnoise){
    noise[,j] <- runif(n,-2,2)
  }
  
  # all balances into one data frame:
  bal <- cbind(rel,noise)
  
  # from ilr to clr
  V <- gsi.buildilrBase(t(codes)) 
  
  Bal_clr <- bal %*% t(V)  
  
  # from clr to composition
  X <-  as.matrix(as.data.frame(acomp(exp(Bal_clr)))) 
  colnames(X) <- paste0(rep("x",D),1:D)  # generic colnames of matrix X
  
  ## STEP 2: performing selection using easyCODA::STEP
  model <- easyCODA::STEP(X, weight = F)  # chooses pairwise logratios that explain maximum variance
  
  
}



#######################################################################################
## Functions to give ranks to logratios (sPCA and STEP)
#######################################################################################

RateImportant <- function(models = models, loads = NULL, names = NULL, important = NULL, D = NULL, a = NULL){
  # THIS FUNCTION IS WRITTEN FOR RESULTS OF SPCA
  # models = results of simulations (either spca or Greenacre)
  # loads = names of loadings = all pairwise logratios
  # names = vector of pairwise logratios without duplicite PLR (contains important AND unimportant PLR)
  # important = vector of (names of) important pairwise logratios (without duplicities)
  # D = number of components/parts
  # a = vector of sparsity parameters
  
  # tab.aux, tab.aux2  = auxiliary tables to store results for further calculations
  
  tab.aux <- matrix(NA, nrow = D*(D-1), ncol = length(a))  # contains values of loadings
  rownames(tab.aux) <- loads
  colnames(tab.aux) <- a
  
  for(j in 1:ncol(tab.aux)){
    SumSq <- apply((models[[j]]$loadings)^2,1,sum)
    tab.aux[,j] <- SumSq
  }
  
  tab.aux2 <- tab.aux  # contains just 1 (for non-zero loadings) or 0
  for(k in 1:nrow(tab.aux2)){
    for(l in 1:ncol(tab.aux2)){
      if(tab.aux2[k,l]>0){
        tab.aux2[k,l] <- 1
      } 
    }
  }
  
  tab.aux2 <- as.data.frame(tab.aux2)
  tab.aux2$total <- rowSums(tab.aux2)  # for each PLR calculates how many times it was chosen to a sparse model (sums all ones in tab.aux2)
  sorted.tab.aux2 <- tab.aux2[order(tab.aux2$total, decreasing=T),]  # sorts in descending order (top PLR -> most frequently chosen ones)
  
  half.tab.aux2 <- sorted.tab.aux2[rownames(sorted.tab.aux2) %in% names,]  # Lose duplicities -> choose rows according to vector "names" (without duplicities)
  pick <- rownames(half.tab.aux2)[1:(D-1)]   # choose first D-1 PLR -> "pick" is arranged in DECREASING ORDER
  
  IsImp <- pick %in% names(important)   # are these PLR in the set of important ones? We search in the vector of important PLR without duplicities
  
  IsImpNum <- as.numeric(IsImp)  # this vector will be used to calculate how many times (over sim runs) was each pairwise logratio important
  
  resulting.vector <- cumsum(IsImp)  # cummulative sum -> tells us which PLR was within the important ones and which was not (vector without duplicities should be used)
  
  ### explained variance
  
  # a) Calculate variance of all D(D-1) PLR -> lose duplicities to have only D(D-1)/2 
  
  logrs.var <- c(apply(models[[1]]$X.pairwise,2,var))   # variances of all D*(D-1) pairwise logratios
  logrs.var.half <- logrs.var[names(logrs.var) %in% names]   # Lose duplicities -> choose values according to vector "names" (without duplicities)
  
  # b) Sort the variance in DECREASING ORDER -> add ranks to them
  logrs.var.sort <- sort(logrs.var.half,decreasing = T)
  logrs.var.rank <- 1:length(logrs.var.sort)
  names(logrs.var.rank) <- names(logrs.var.sort)  # this vector gives info about the rank of each pairwise logratio; 1 = the highest variance
  
  # c) From the vector of ranks, chose important and unimportant pairwise logratios
  logrs.var.imp <- names(logrs.var.rank) %in% names(important)      # which logratios are important (TRUE) and unimportant (FALSE)
  
  rank.important <- logrs.var.rank[logrs.var.imp==T]
  rank.unimportant <- logrs.var.rank[logrs.var.imp==F]
  
  rank.both <- logrs.var.rank   # rank. both should be defined in this way -> order of PLR must be kept (not like on previous line - first important, then unimportant)
  
  # d) Cumulative sum of ranks -> "pick" contains the most stable PLR, so we chose these and then cumulative sum
  
  aux <- ifelse(is.na(rank.important[pick]),0,rank.important[pick])    # the rank of important PLR in the order as PLR are in "pick"; 0 for NA
  
  sum.rank.important <- cumsum(aux)
  sum.rank.both <- cumsum(rank.both[pick])
  
  # e) Variances of picked pairwise logratios 
  var.important <- logrs.var.half[names(aux)]
  
  var.both <- logrs.var.half[pick]
  
  ###
  
  # Which output will be plotted: cummulative_sum, Is_important, cumsum_rank_important, cumsum_rank_both
  
  result <- list("cummulative_sum" = resulting.vector, "Is_important" = IsImpNum, "rank_important" = rank.important, "rank_unimportant" = rank.unimportant, "rank_both" = rank.both,
                 "cumsum_rank_important" = sum.rank.important, "cumsum_rank_both" = sum.rank.both, "var_important" = var.important,
                 "var_both" = var.both)  # "explvar" = prop2; "explvar" = prop, "totvar" = tot.var
  
}

RateImportant.G <- function(models = models, modelsSPCA = modelsSPCA,important = NULL){   # , tot.var = tot.var
  # THIS FUNCTION IS WRITTEN FOR RESULTS OF STEP (GREENACRE)
  # modelsSPCA = loads models from spca -> obtain a matrix of all pairwise logratios
  # important = has to be a vector of all important PLR (with duplicities) -> we use the vector with duplicities because STEP chooses e.g. ln(A/B) but not ln(B/A) and we could miss it
  # totvar = total variance of the matrix of D(D-1)/2 pairwise logratios -> calculated in RateImportant function
  
  pick1.Green <- paste0("ln(",models$names,")")       # match the names with ours -> add ln(/)
  IsImp1.Green <- pick1.Green %in% names(important)   # important1A  -> important PLR with duplicities (both ln(A/B) and ln(B/A)) 
  
  IsImpNum <- as.numeric(IsImp1.Green)  # this vector will be used to calculate how many times (over sim runs) was each pairwise logratio important
  
  resulting.vector.Green <- cumsum(IsImp1.Green)
  
  ### explained variance
  
  # a) Calculate variance of D-1 PLR calculated by STEP
  
  logrs.var <- c(apply(modelsSPCA[[1]]$X.pairwise,2,var))   # variances of all D*(D-1) pairwise logratios -> from the matrix of all pairwise logratios
  
  # b) Sort the variance in DECREASING ORDER -> add ranks to them
  logrs.var.sort <- sort(logrs.var, decreasing = T)
  logrs.var.rank <- rep(1:(length(logrs.var.sort)/2), each = 2) # why each = 2 -> now we have all D(D-1) logratios and var(A/B) = var(B/A); (length(logrs.var.sort)/2) -> D(D-1)/2 
  names(logrs.var.rank) <- names(logrs.var.sort)  # this vector gives info about the rank of each pairwise logratio; 1 = the highest variance
  
  # as names of Greenacre's PLR do not always match names of spca results (A/B vs B/A), all D(D-1) pairwise logratios are considered. However, rank of A/B is the same as of B/A and it's taken into account.
  
  # c) From the vector of ranks, chose important and unimportant pairwise logratios
  logrs.var.imp <- names(logrs.var.rank) %in% names(important)      # which logratios are important (TRUE) and unimportant (FALSE)
  
  rank.important <- logrs.var.rank[logrs.var.imp==T]
  rank.unimportant <- logrs.var.rank[logrs.var.imp==F]
  
  rank.both <- logrs.var.rank   # rank. both should be defined in this way -> order of PLR must be kept 
  
  # d) Cumulative sum of ranks -> "pick" contains the most stable PLR, so we chose these and then cumulative sum
  
  aux <- ifelse(is.na(rank.important[pick1.Green]),0,rank.important[pick1.Green])    # the rank of important PLR in the order as PLR are in "pick"; 0 for NA
  
  sum.rank.important <- cumsum(aux)
  sum.rank.both <- cumsum(rank.both[pick1.Green])
  
  # e) Variances of picked pairwise logratios
  var.important <- logrs.var[names(aux)]
  
  var.both <- logrs.var[pick1.Green]
  
  ###
  
  # Which output will be plotted: cummulative_sum, Is_important, cumsum_rank_important, cumsum_rank_both
  
  result <- list("cummulative sum" = resulting.vector.Green, "Is_important" = IsImpNum, "rank_important" = rank.important, "rank_unimportant" = rank.unimportant, "rank_both" = rank.both,
                 "cumsum_rank_important" = sum.rank.important, "cumsum_rank_both" = sum.rank.both, "var_important" = var.important,
                 "var_both" = var.both)     # "explvar" = prop2; "explvar" = prop
}


#######################################################################################
## Functions for plots
#######################################################################################

## a) Line plot: sparsity and explained variance

plot_line <- function(input = input, a = a){
  
  tab1A <- input
  
  laby <- seq(0,100,10)
  labx <- round(a,4)
  
  T1A <- xtabs(`Percentage`~`alpha`, data = tab1A)
  T1A
  
  expl.vars.1A <- (tab1A$`var`)*100
  expl.vars.1A
  
  T1A.df <- as.data.frame(cbind(T1A,expl.vars.1A))
  head(T1A.df)
  T1A.df$alpha <- rownames(T1A.df)
  
  colnames(T1A.df) <- c("Percentage","Var","alpha")
  T1A.df$alpha <- as.numeric(T1A.df$alpha)  
  T1A.df$ind <- 1:51     # just for the case of plot -> it has problems plotting log(alpha); ticks will be alpha
  head(T1A.df)
  
  T1A.df$alpha <- as.numeric(T1A.df$alpha)   # log: values on the x axis will be spread evenly
  
  name1A <- rep("Percentage",51)
  part1A <- cbind(T1A.df$ind,T1A.df$Percentage,T1A.df$alpha,name1A)
  name2A <- rep("Var",51)
  part2A <- cbind(T1A.df$ind,T1A.df$Var,T1A.df$alpha,name2A)
  partsA <- as.data.frame(rbind(part1A,part2A))
  colnames(partsA) <- c("ind","value","alpha","Key")   # "Key" instead of "name" because of the plot legend
  partsA$ind <- as.numeric(partsA$ind)
  partsA$value <- as.numeric(partsA$value)
  partsA$alpha <- as.numeric(partsA$alpha)
  partsA$Key <- as.factor(partsA$Key)
  head(partsA)
  
  labx.pick <- rep("",51)
  labx.pick[c(T,F)] <- as.character(labx)[c(T,F)]  # [c(T,F,F,F,F)]
  labx.pick
  
  g1 <- ggplot2::ggplot(data=partsA, aes(x=ind, y=value, group=Key, shape = Key)) +
    geom_line(aes(color=Key))+
    geom_point(aes(color=Key),size = 4)+
    scale_color_manual(values=c('Percentage'="black",'Var'="firebrick3"),
                       breaks=c('Percentage','Var'),
                       labels=c('Zero logratios', 'Explained variability'))+
    scale_shape_manual(values=c(16, 17),labels=c('Zero logratios', 'Explained variability')) + # change shape of point manually
    theme_bw()+
    theme(axis.title=element_text(size=30,face="bold"))+
    theme(axis.text.x = element_text(size=15,angle = 90))+
    theme(axis.text.y = element_text(size=19))+
    scale_x_continuous(breaks=1:51, labels=as.character(labx.pick),expand = c(0, 1.5))+  # plot every fifth label
    scale_y_continuous(breaks = seq(0,105,10))+
    #ggtitle("Same-sized marker blocks")+
    labs(x ="Sparsity parameter", y = "Range in %")+  # title="Explained variability based on sparsity",
    #theme(plot.title = element_text(size=18))+
    theme(legend.key.size = unit(1.5, 'cm'))+   # size of a legend
    theme(legend.position="top")+               # legend position
    theme(legend.title = element_text(size=28))+ # legend title size
    theme(legend.text = element_text(size=28))   # legend text size
  return(g1)
}


## b) Line plot: FNR and FPR

plot_rates <- function(input,a=a){
  
  rates1A <- input
  
  labx <- round(a,4)
  labx.pick <- rep("",51)
  labx.pick[c(T,F)] <- as.character(labx)[c(T,F)]  # [c(T,F,F,F,F)]
  labx.pick
  
  typ <- factor(c(rep("FNR",51),rep("FPR",51)))
  v1 <- c(rates1A[,"FNR"],rates1A[,"FPR"])
  a.vec <- log(c(a,a))   # log scale of a sparsity parameter to space it evenly
  df1 <- data.frame(v1,typ,a.vec)
  colnames(df1) <- c("value","Rate","alpha")
  
  p1 <- ggplot2::ggplot(data=df1, aes(x=alpha, y=value, group=Rate, shape = Rate, color = Rate)) +
    geom_line(aes(group=Rate))+
    geom_point(aes(shape=Rate),size = 4)+
    scale_color_manual(values=c('FPR'="black",'FNR'="dodgerblue3"), 
                       breaks=c('FNR','FPR'),
                       labels=c('FNR', 'FPR'))+
    theme_bw()+
    labs(x="Sparsity parameter", y = "Values of FNR and FPR")+
    theme(axis.title=element_text(size=30,face="bold"))+
    theme(axis.text.x = element_text(size=15, angle = 90))+
    theme(axis.text.y = element_text(size=19))+
    scale_x_continuous(breaks = round(log(a),4), labels=as.character(labx.pick))+   # breaks = c(1,seq(5,D-1,5)), labels = round(a,4)
    theme(plot.title = element_text(size=18))+
    theme(legend.key.size = unit(1.5, 'cm'))+   # size of a legend
    theme(legend.position="top")+               # legend position
    theme(legend.title = element_text(size=28))+ # legend title size
    theme(legend.text = element_text(size=28)) #+  # legend text size
  return(p1)
}


## c) Line-step plot: comparison of sPCA and easyCODA::STEP function

plot_line_step <- function(input.sPCA = input.sPCA, input.STEP = input.STEP, input.importantPLR = input.importantPLR, D = D){
  
  tabSums.1A <- input.sPCA
  tabSums.1A.STEP <- input.STEP
  important1A.half <- input.importantPLR
  
  pl1A.df <- data.frame(number = rep(1:(D-1),2),
                        value = c(100*colMeans(tabSums.1A)/length(important1A.half),100*colMeans(tabSums.1A.STEP)/length(important1A.half)),
                        Method = c(rep("sPCA",(D-1)),rep("STEP",(D-1))))
  
  pl1A <- ggplot2::ggplot(pl1A.df, aes(x = number, y = value, group = Method))+
    geom_step(aes(color = Method), size = 0.8)+
    ylim(0,25)+
    scale_x_discrete(name ="Order of PLR", 
                     limits=paste0(1:(D-1)))+
    scale_color_manual(values=c('sPCA'="firebrick3",'STEP'="dodgerblue3"),
                       breaks=c('sPCA','STEP'),
                       labels=c('sPCA', 'STEP'))+
    theme_bw()+
    theme(axis.title=element_text(size=25,face="bold"))+
    theme(axis.text.x = element_text(size=25))+
    theme(axis.text.y = element_text(size=25))+
    labs(y = "Average % of important PLRs")+
    theme(legend.key.size = unit(1.5, 'cm'))+   # size of a legend
    theme(legend.title = element_text(size=20))+ # legend title size
    theme(legend.text = element_text(size=20))   # legend text size
  
  return(pl1A)
  
}

## d) Boxplots: comparison of sPCA and easyCODA::STEP function

plot_boxes <- function(input.sPCA = input.sPCA, input.STEP = input.STEP, D = D){
  
  difference.1A <- input.sPCA
  difference.1A.STEP <- input.STEP
  
  
  melted1A <- rbind(reshape2::melt(difference.1A),reshape2::melt(difference.1A.STEP))
  box1A <- data.frame(order = factor(rep(rep(1:(D-1),each=100),2)),
                      group = factor(rep(c("sPCA","STEP"),each = 100*(D-1))),
                      value = melted1A$value)
  
  # a) boxplot
  
  b1A <- ggplot2::ggplot(box1A, aes(x=order, y=value, fill=group)) + 
    geom_boxplot()+
    scale_fill_manual(values=c('sPCA'="firebrick3",'STEP'="dodgerblue3"),
                      breaks=c('sPCA','STEP'),
                      labels=c('sPCA', 'STEP'))+
    theme_bw()+
    theme(axis.title=element_text(size=25,face="bold"))+
    theme(axis.text.x = element_text(size=25))+
    theme(axis.text.y = element_text(size=25))+
    #scale_y_continuous(breaks = seq(0,150,25))+
    labs(x = "Order of PLR", y = "Difference of a rank.method and rank.ideal", fill = "Method")+
    theme(legend.key.size = unit(1.5, 'cm'))+   # size of a legend
    theme(legend.title = element_text(size=20))+ # legend title size
    theme(legend.text = element_text(size=20))   # legend text size
  return(b1A)
}

## e) Stability of PLRs: heatmap

plot_heat_PLRs <- function(input = input, a = a, X = X,exvar = exvar){
  
  #X1 <- X
  D2 <- ncol(X)
  models <- input
 
  load.names <- rownames(models[[1]]$loadings)
  
  tabulka <- matrix(NA, nrow = D2*(D2-1), ncol = length(a))  # 10*9 = D*(D-1) -> without trace elements
  rownames(tabulka) <- load.names
  colnames(tabulka) <- a
  
  for(i in 1:ncol(tabulka)){
    SumSq <- apply((models[[i]]$loadings)^2,1,sum)
    tabulka[,i] <- SumSq
  }
  
  tabulka2 <- tabulka
  for(i in 1:nrow(tabulka2)){
    for(j in 1:ncol(tabulka2)){
      if(tabulka2[i,j]>0){
        tabulka2[i,j] <- 1
      } 
    }
  }
  
  trial.table <- as.data.frame(tabulka2)
  trial.table$total <- rowSums(trial.table)
  sorted.trial.table <- trial.table[order(trial.table$total, decreasing=T),]
  
  half.table <- sorted.trial.table[rownames(sorted.trial.table) %in% names(exvar),]  # so that names in half.table coincide with those in exvar; values of exvar will be added to the table
  
  # reorder exvar according to rownames of a half.table
  reordered.exvar <- exvar[match(rownames(half.table),names(exvar))]
  
  half.table$exvar <- reordered.exvar      #exvar
  half.table$logrs <- rownames(half.table)
  colnames(half.table) <- c(round(as.numeric(colnames(half.table)[-c(52,53,54)]),5),"total","exvar","logrs")

  heatmapa <- half.table %>%
    pivot_longer(-logrs, names_to = "alpha", values_to = "Values")
  
  heat.complex <-ggplot(mapping = aes(x = alpha, y = logrs)) +
    geom_tile(data = filter(heatmapa, logrs != "exvar", alpha != "exvar"), aes(fill = Values),color = "#00000022") +
    geom_point(data = filter(heatmapa, logrs == "exvar" | alpha == "exvar"), aes(color = Values), size = 13, shape = 15) + #shape = 15 -> filled square
    geom_text(data = heatmapa, aes(label = round(Values, 1)), size = 5) + 
    scale_color_gradient2(name = "Explained variability in %",low = "red", mid = "white", high = "grey64", midpoint = 0) +
    scale_fill_gradientn(colours = c("white", "dodgerblue3", "yellowgreen","yellow", "orange","firebrick1" ,"firebrick4"),
                         values = c(0, 1, 10, 20, 30, 40, 51)/51, limits = c(0, 51),
                         name = "Count") +
    scale_x_discrete(limits = unique(heatmapa$alpha),position = "top") +
    scale_y_discrete(limits = rev(unique(heatmapa$logrs))) +
    theme_minimal() +
    labs(x ="Sparsity parameter", y = "Pairwise logratios")+
    theme(
      legend.key.size = unit(1.5, 'cm'), #change legend key size
      legend.text = element_text(size=28),   # legend text size
      text = element_text(size = 27),
      axis.title=element_text(size=30,face="bold"),
      axis.text.y = element_text(size = 20, colour = "black"),
      strip.text.y = element_text(angle = 0),
      axis.text.x = element_text(size = 18, colour = "black", angle = 90),
      legend.direction = "horizontal",
      legend.position = "bottom",
      panel.grid.major = element_blank()
    )
  heat.complex
  
  
  
}



## f) Stability of parts: heatmap 

plot_heat_parts <- function(input = input, a = a, X = X){
  
  X1 <- X
  models <- input
  D2 <- ncol(X1)
  
  load.names <- rownames(models[[1]]$loadings)
  
  tabulka <- matrix(NA, nrow = D2*(D2-1), ncol = length(a)) 
  rownames(tabulka) <- load.names
  colnames(tabulka) <- a
  
  for(i in 1:ncol(tabulka)){
    SumSq <- apply((models[[i]]$loadings)^2,1,sum)
    tabulka[,i] <- SumSq
  }
  
  tabulka2 <- tabulka
  for(i in 1:nrow(tabulka2)){
    for(j in 1:ncol(tabulka2)){
      if(tabulka2[i,j]>0){
        tabulka2[i,j] <- 1
      } 
    }
  }
  
  tab.counts <- matrix(NA,nrow=ncol(X1),ncol=length(a))  # to collect how many times each part appears in a logratio (throughout different sparsity settings)
  rownames(tab.counts) <- sort(colnames(X1))
  for(i in 1:ncol(tabulka2)){
    nonzero.names <-names(tabulka2[,i][tabulka2[,i]==1])  # which logratios are non-zero
    if(length(nonzero.names)==0){
      tab.counts[,i] <- rep(0,nrow(tab.counts))
    } else {
      inside.brackets <-gsub("[\\(\\)]", "", regmatches(nonzero.names, gregexpr("\\(.*?\\)", nonzero.names)))  # separate "the bracket" of the name, i.e. from ln(a/b) get a/b
      M <- as.data.frame(str_split(inside.brackets, "/", n = Inf, simplify = T))  # separate a/b and put into matrix
      M[,2] <- factor(M[,2], levels = rownames(tab.counts))   # to have all components as factor levels; with higher sparsity, some parts will have 0 in table()
      counts <- c(table(M$V2))   # How many times each component appears in a non-zero logratio
      tab.counts[,i] <- counts
    }
  }
  
  df.tab.counts <- as.data.frame(tab.counts)
  colnames(df.tab.counts) <- as.character(a)  # colnames = values of sparsity parameters
  
  melted <- df.tab.counts %>%
    rownames_to_column() %>%
    gather(colname, value, -rowname)
  colnames(melted) <- c("Parts","alpha","Count")
  head(melted)
  
  melted$Parts <- factor(melted$Parts)
  melted$alpha <- factor(melted$alpha)
  melted$categ <- factor(melted$Count)
  kat <- melted$categ
  
  laby <- seq(0,100,10)
  labx <- round(a,5)
  
  
  plot1 <- ggplot(melted,aes(alpha, Parts)) +
    geom_tile(aes(fill = kat, width=1), color = "black") +    
    coord_fixed(ratio=1)+
    scale_fill_manual(breaks = levels(melted$categ),values=colorRampPalette(c("white","darkblue"))(11))+   # values = pal
    theme(panel.background = element_blank(),   # remove grey background
          plot.background = element_blank())+
    theme(axis.title=element_text(size=26,face="bold"))+
    theme(axis.text.x = element_text(size=19, angle = 90))+
    theme(axis.text.y = element_text(size=16))+
    scale_x_discrete(expand = c(0,0),labels=as.character(labx))+
    scale_y_discrete(limits = rev(levels(melted$Parts)),expand = c(0, 0))+   # ,labels = as.character(levels(melt$Parts))
    theme(legend.key.size = unit(0.6, 'cm'))+   # size of a legend
    theme(legend.title = element_text(size=22))+ # legend title size
    theme(legend.text = element_text(size=22))+   # legend text size
    #theme(legend.position="bottom")+            # position of a legend
    labs(fill = "Frequency")+
    xlab("Sparsity parameter alpha") + 
    ylab("Compositional parts") 
  return(plot1)
}



## g) Stability of PLRs: heatmap with ranks given by STEP (Illustrative example in Fig.1)

plot_heat_PLRs_illustration <- function(input.sPCA = input.sPCA, input.STEP = input.STEP, input.STEP2 = input.STEP2,a=a){
  
  res1A <- input.sPCA
  res1A.STEP <- input.STEP
  resulting.list1A.STEP <- input.STEP2
  
  X.pairs <- res1A[[1]][[1]]$X.pairwise   # matrix of all pairwise logratios is the same in all models
  
  logratio.vars <- apply(X.pairs,2,var)  
  
  sort(logratio.vars, decreasing=T) 
  # pick just odd pairwise log-ratios (var(ln(a/b)) = var(ln(b/a)))
  sort(logratio.vars, decreasing=T)[c(TRUE,FALSE)] 
  exvar <- round(sort(logratio.vars, decreasing=T)[c(TRUE,FALSE)] ,3)
  exvar
  
  # exvar2 -> really variance explained by each logratio
  total.variance <- sum(logratio.vars)/(2*D1) 
  percent <- (logratio.vars/total.variance)*100
  sort(percent,decreasing=T)    
  # pick just odd pairwise log-ratios (var(ln(a/b)) = var(ln(b/a)))
  sort(percent,decreasing=T)[c(TRUE,FALSE)]
  exvar2 <- round(sort(percent,decreasing=T)[c(TRUE,FALSE)],3)  # round sorted odd pairwise logratio
  
  logr.names <- rownames(res1A[[1]][[1]]$loadings) # in each model, rownames are the same
  
  tabulka <- matrix(NA, nrow = D1*(D1-1), ncol = length(a))  
  rownames(tabulka) <- logr.names
  colnames(tabulka) <- a
  
  # table collecting non-zero and zero values of both loadings (hence the sum of loadings^2 -> to get one value); logratios in rows, in columns results for each choice of alpha (a)
  for(i in 1:ncol(tabulka)){
    SumSq <- apply((res1A[[1]][[i]]$loadings)^2,1,sum)
    tabulka[,i] <- SumSq
  }
  
  # table of 0,1 values: 0 = logratio is sparse/zero, 1 = logratio is non-zero; heatmap is based on this table
  tabulka2 <- tabulka
  for(i in 1:nrow(tabulka2)){
    for(j in 1:ncol(tabulka2)){
      if(tabulka2[i,j]>0){
        tabulka2[i,j] <- 1
      } 
    }
  }
  
  trial.table <- as.data.frame(tabulka2)
  trial.table$total <- rowSums(trial.table)
  sorted.trial.table <- trial.table[order(trial.table$total, decreasing=T),]  # trial table is sorted 
  
  # reorder exvar2 according to rownames of a half.table
  
  half.table <- sorted.trial.table[rownames(sorted.trial.table) %in% names(exvar2),]  # so that names in half.table coincide with those in exvar; values of exvar will be added to the table
  
  #reorder exvar2 according to rownames of a half.table
  reordered.exvar <- exvar2[match(rownames(half.table),names(exvar2))]
  
  half.table$exvar <- reordered.exvar    #exvar
  half.table$logrs <- rownames(half.table)
  colnames(half.table) <- c(round(as.numeric(colnames(half.table)[-c(52,53,54)]),5),"total","exvar","logrs")
  
  # -> new column: showing the rank 1-9 of STEP-chosen pairwise logratios (which PLR was the first by STE? which was the second etc.)
  res1A.STEP[[1]]$names  # that's how STEP chose pairwise logratios
  resulting.list1A.STEP[[1]]$var_both
  
  # this is the rank of the chosen pairwise logratios (by STEP) given by the variance of each pairwise logratio
  resulting.list1A.STEP[[1]]$rank_both[names(resulting.list1A.STEP[[1]]$var_both)]
  
  names(resulting.list1A.STEP[[1]]$var_both) %in% half.table$logrs
  
  STEP.order <- rep(NA,length(half.table$logrs))  # empty vector of NA -> to be filled by 1-9 for STEP-chosen PLRs, NA (blank) otherwise
  names(STEP.order) <- half.table$logrs
  STEP.order
  
  STEP.order.chosen <- 1:9
  names(STEP.order.chosen) <- names(resulting.list1A.STEP[[1]]$var_both)  # giving ranks 1-9 to PLRs chosen by STEP; var_both has names in the format ln(a/b)
  STEP.order.chosen
  
  STEP.order[names(STEP.order.chosen)] <- STEP.order.chosen  # combining information together -> ranks 1-9 for STEP-chosen PLRs and 0 otherwise
  STEP.order
  
  # add this new vector to half.table
  half.table$STEP <- STEP.order

  tab1 <- half.table[,-c(52,53,55)]
  
  tab2 <- tab1 %>%
    pivot_longer(-logrs, names_to = "alpha", values_to = "Values")
  
  list_correct_order <- rownames(tab1)  # for arrangement of rows in a heatmap
  
  df.diff <- ggplot(tab2, aes(x = alpha, y = logrs)) +
    geom_tile(aes(fill = Values),color = "#00000022") +
    geom_text(size = 3, aes(label = Values)) + #displays cell values
    scale_fill_gradient2(low = "white", #colors
                         high = "dodgerblue3") +
    guides(fill = "none") +   # turns off legend for fill with blue and white (fill of the heatmap)
    scale_y_discrete(limits=rev(list_correct_order))  # arrange rows in a heatmap
  
  c_total <- tab2 %>% 
    group_by(logrs) %>% 
    summarise(Values = sum(Values)) %>% 
    mutate(alpha = 'Total')
  c_total <- c_total[order(match(c_total$logrs,list_correct_order)),] # order rows of c_total
  
  
  c_exvar <- data.frame(logrs = list_correct_order,
                        Values = reordered.exvar,
                        alpha = rep("exvar",length(list_correct_order)))
  
  STEP.order[is.na(STEP.order)] <- ""
  c_step <- data.frame(logrs = list_correct_order,
                       Values = STEP.order,
                       alpha = rep("STEP",length(list_correct_order)))
  
  p <- df.diff + 
    # Total
    new_scale_color()+
    geom_point(data = c_total, 
               aes(color = Values), 
               size = 13.5,
               shape = 15) +
    scale_color_gradientn(name = "Count",colors = c("dodgerblue4", "dodgerblue3", "yellowgreen","yellow", "orange","firebrick1" ,"firebrick4")) +
    geom_text(data = c_total, size = 4.5, aes(label = round(Values,2))) +
    # exvar
    new_scale_color()+
    geom_point(data = c_exvar, 
               aes(color = Values), 
               size = 13.5,
               shape = 15)+
    scale_color_gradient2(name = "Explained variance",low = "white", high = "grey64")+
    geom_text(data = c_exvar, size = 4.5, aes(label = round(Values,1)))+
    # STEP
    new_scale_color() +
    geom_point(data = c_step, 
               aes(color = Values), 
               size = 12,
               shape = 15) +
    scale_color_manual(breaks = unique(STEP.order), values = rep("white",length(unique(STEP.order)))) +
    geom_text(data = c_step, size = 5, aes(label = Values))+
    guides(color = "none")+
    theme_minimal() +
    labs(x ="Sparsity parameter", y = "Pairwise logratios")
  
  # Final plot adjustments:
  p1 <- p+scale_x_discrete(limits = c(unique(tab2$alpha),"Total","exvar","STEP"),position = "top")+
    theme(
      legend.key.size = unit(1.5, 'cm'), #change legend key size
      legend.text = element_text(size=28),   # legend text size
      text = element_text(size = 27),
      axis.title=element_text(size=30,face="bold"),
      axis.text.y = element_text(size = 20, colour = "black"),
      strip.text.y = element_text(angle = 0),
      axis.text.x = element_text(size = 18, colour = "black", angle = 90),
      legend.direction = "horizontal",
      legend.position = "bottom",
      panel.grid.major = element_blank()
    )
  
  p1
}



