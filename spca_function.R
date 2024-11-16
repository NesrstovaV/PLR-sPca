#############################################################################
## Function making a matrix of D(D-1) logratios and calculating sparse PCA
#############################################################################

spca_funkce <- function(X,alpha = 0.01, beta = 1e-04, k = (D-1), draw = T){
  # X = input CoDa matrix
  # alpha = sparsity parameter
  # beta = tuning parameter; beta = 0 indicates lasso penalty
  # k = number of PCs to be calculated
  # draw = draw a biplot (T) or not
  
  D <- ncol(X)
  n <- nrow(X)
  
  X.pairs <- matrix(NA,nrow = n, ncol = D*(D-1))  # matrix of pairwise logratios
  col.names <- c()  # a vector to store colnames of a matrix X.pairs
  ind <- 1
  for(i in 1:D){
    for(j in 1:D){
      if(i!=j){
        logr <- log(X[,i]/X[,j])
        name.logr <- c(paste0("ln(",colnames(X)[i],"/",colnames(X)[j],")"))
        X.pairs[,ind] <- logr
        col.names[ind] <- name.logr
        ind <- ind+1
      } else {
        next
      }
    }
  }
  colnames(X.pairs) <- col.names
  
  model <- sparsepca::spca(X.pairs, alpha=alpha, beta = beta, k = k, verbose = F)
  
  # Model loadings:
  P <- model$loadings
  rownames(P) <- colnames(X.pairs)
  
  # explained variance:
  summary.model <- summary(model)
  expl.var <- summary(model)[3,1:k]  # 3rd line = proportion of variance, cols 1:k
  
  # How many logratios have zero loadings ("full zero" loadings; all loadings are 0 for a certain logratio)
  isZero <- sum(rowSums(P)==0)
  zeros_sub <- apply(P, 1, function(row) all(row == 0))
  allZero <- P[zeros_sub,]
  
  # Which variables form non-zero logratios
  row_sub <- apply(P, 1, function(row) any(row !=0 ))  # "all" instead of "any" chooses only those rows without zero (any)
  nonZero <- P[row_sub,]
  
  # Final output: store a matrix of pairwise logratios and a spca model
  result <- list("X.pairwise" = X.pairs, "model" = model, "loadings" = P, "model summary" = summary.model, "expl.var" = expl.var,
                 "number of zero logratios" = isZero, "table of all zero" = allZero) 
  
  
}