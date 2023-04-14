
##################### divide by its sum each row (to obtain probabilities)
normByRow <- function(Mat){
  U <- t(scale(t(Mat), center = FALSE, scale = colSums(t(Mat))))
  attr(U,"scaled:scale")<- NULL
  return(U)
}

##################### IN VE step, to obtain the tau  
fromBtoTau <- function(B,eps = 10^-5){
  B <- B - matrix(apply(B,1,max),nrow = nrow(B),ncol = ncol(B),byrow = FALSE)
  temp <- exp(B)
  temp2 <- normByRow(temp)
  temp2[temp2 < eps] <- eps
  temp2[temp2 > (1 - eps)] <- 1 - eps
  Tau <- normByRow(temp2)
  attr(Tau,"scaled:scale")<- NULL
  return(Tau)
}

##################### Difference between tau and tauOld
distTau  <- function(collecTau,collecTauOld)
{
  M <- length(collecTau)
  vdis <- sapply(1:M,function(m){
    drow <- sqrt(sum(as.vector(collecTau[[m]]$row - collecTauOld[[m]]$row)^2))
    dcol <- sqrt(sum(as.vector(collecTau[[m]]$col - collecTauOld[[m]]$col)^2))
    return(drow + dcol)
  })
  return(sum(vdis))
}

##################### Difference between paramPostOld and paramPost

disthyperparamPost <- function(hyperparamPost,hyperparamPostOld){
  
  d1 <- sqrt(sum(as.vector((hyperparamPost$connectParam$alpha - hyperparamPostOld$connectParam$alpha)^2)))
  d2 <- sqrt(sum(as.vector((hyperparamPost$connectParam$beta -  hyperparamPostOld$connectParam$beta )^2)))
  d3 <- sqrt(sum(as.vector((hyperparamPost$blockProp$row -  hyperparamPostOld$blockProp$row )^2)))
  d4 <- sqrt(sum(as.vector((hyperparamPost$blockProp$col -  hyperparamPostOld$blockProp$col )^2)))
  return(d1 + d2 + d3 +d4)
}


####################### Reorder groups row / cols. 
reorderBlocks  <- function(hyperparamPost, collecTau ,model){
  
  meanPost <- hyperparamPost$connectParam$alpha / (hyperparamPost$connectParam$alpha + hyperparamPost$connectParam$beta)
  colMean <- apply(meanPost,2,mean); 
  ordCol <- order(colMean,decreasing = TRUE)
  rowMean <- apply(meanPost,1,mean)
  ordRow <- order(rowMean,decreasing = TRUE)
  
  hyperparamPost$connectParam <- lapply(hyperparamPost$connectParam,function(u){u[ordRow,ordCol]})
  if(model == 'iidColBipartiteSBM'){  
    hyperparamPost$blockProp$row <- hyperparamPost$blockProp$row[ordRow]
    hyperparamPost$blockProp$col <- hyperparamPost$blockProp$col[ordCol]
  }
  if(model == 'piColBipartiteSBM'){
    hyperparamPost$blockProp$row <- hyperparamPost$blockProp$row[,ordRow]
    hyperparamPost$blockProp$col <- hyperparamPost$blockProp$col[,ordCol]
  }
  
  
  collecTau <- lapply(collecTau ,function(tau){
    tauNew <- tau
    tauNew$row <-tau$row[,ordRow]
    tauNew$col <-tau$col[,ordCol]
    return(tauNew)}
    )
  return(list(hyperparamPost = hyperparamPost,collecTau  = collecTau))
  
  
}


#######################""
mylBetaFunction = function(alpha){
   sum(lgamma(alpha))-lgamma(sum(alpha))
}


################################## set prior param
setHyperparamPrior <- function(M,KRow,KCol,emissionDist, model){
  # M  : number of bipartite networks
  # Krow  : number of blocks in row
  # Kcol  : number of blocks in col
  # emissionDistr : poisson or bernoulli
  
  hyperparamPrior <- list()
  
  
  
  hyperparamPrior$blockProp <- list()
  if (model == 'iidColBipartiteSBM'){
    hyperparamPrior$blockProp$row = rep(1/KRow,KRow) ### Jeffreys dirichlet prior the pi and rho
    hyperparamPrior$blockProp$col = rep(1/KCol,KCol)
  }
  if (model == 'piColBipartiteSBM'){
    hyperparamPrior$blockProp$row = matrix(1/KRow,M,KRow) ### Jeffreys dirichlet prior the pi^m and rho^m
    hyperparamPrior$blockProp$col = matrix(1/KCol,M,KCol)
  }
  
  hyperparamPrior$connectParam<- list()
  if (emissionDist =='bernoulli'){
    hyperparamPrior$connectParam$alpha <- matrix(1,KRow,KCol) ### Uniform prior on the alpha_{kl}
    hyperparamPrior$connectParam$beta <- matrix(1,KRow,KCol)
  }
  if (emissionDist =='poisson'){
    hyperparamPrior$connectParam$alpha <- matrix(1,KRow,KCol) ### Exp de param 1/100
    hyperparamPrior$connectParam$beta <- matrix(1/100,KRow,KCol)
  }
  return(hyperparamPrior)
}





