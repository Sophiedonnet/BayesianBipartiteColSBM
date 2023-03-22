
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

distPostParam <- function(postParam,postParamOld){
  
  d1 <- sqrt(sum(as.vector((postParam$connectParam$alpha - postParamOld$connectParam$alpha)^2)))
  d2 <- sqrt(sum(as.vector((postParam$connectParam$beta -  postParamOld$connectParam$beta )^2)))
  d3 <- sqrt(sum(as.vector((postParam$blockProp$row -  postParamOld$blockProp$row )^2)))
  d4 <- sqrt(sum(as.vector((postParam$blockProp$col -  postParamOld$blockProp$col )^2)))
  return(d1 + d2 + d3 +d4)
}


####################### Reorder groups row / cols. 
reorderBlocks  <- function(postParam, collecTau){
  
  meanPost <- postParam$connectParam$alpha / (postParam$connectParam$alpha + postParam$connectParam$beta)
  colMean <- apply(meanPost,2,mean); 
  ordCol <- order(colMean,decreasing = TRUE)
  rowMean <- apply(meanPost,1,mean)
  ordRow <- order(rowMean,decreasing = TRUE)
  
  postParam$connectParam <- lapply(postParam$connectParam,function(u){u[ordRow,ordCol]})
  postParam$blockProp$row <- postParam$blockProp$row[,ordRow]
  postParam$blockProp$col <- postParam$blockProp$col[,ordCol]
  collecTau <- lapply(collecTau ,function(tau){
    tauNew <- tau
    tauNew$row <-tau$row[,ordRow]
    tauNew$col <-tau$col[,ordCol]
    return(tauNew)}
    )
  return(list(postParam = postParam,collecTau  = collecTau))
  
  
}

############################## log marg likelihood
computeLogLikMarg <- function(collecNetwork, collecTau,priorParam,model){

  
  
  collecZMAP <- lapply(collecTau,function(tau){
    indZ <- tau
    indZ$row <- t(sapply(1:nrow(tau$row),function(i){u <- 0*tau$row[i,]; u[which.max(tau$row[i,])]=1; return(u)}))
    if(min(dim(indZ$row))==1){indZ$row = matrix(indZ$row,ncol=1)}
    indZ$col <- t(sapply(1:nrow(tau$col),function(i){u <- 0*tau$col[i,]; u[which.max(tau$col[i,])]=1; return(u)}))
    if(min(dim(indZ$col))==1){indZ$col = matrix(indZ$col,ncol=1)}
    return(indZ)
  })

  #browser()
  if(model=='bernoulli'){
    S1 <- lapply(1:length(collecNetwork),function(m){t(collecZMAP[[m]]$row)%*% collecNetwork[[m]] %*% collecZMAP[[m]]$col})
    S0 <- lapply(1:length(collecNetwork),function(m){t(collecZMAP[[m]]$row)%*% (1-collecNetwork[[m]]) %*% collecZMAP[[m]]$col})
    
    ahat <- Reduce('+',S1)
    bhat <- Reduce('+',S0)
    lprobYZ <- sum(lbeta(ahat+priorParam$connectParam$alpha,bhat+priorParam$connectParam$beta) -lbeta(priorParam$connectParam$alpha,priorParam$connectParam$beta))
   # lprobYZ2 <- sum(lbeta(as.vector(S1)+priorParam$connectParam$alpha,bhat+priorParam$connectParam$beta)) - length(collecNetwork)*lbeta(priorParam$connectParam$alpha,priorParam$connectParam$beta)
    
    
    
  }
  
  if(model=='poisson'){
    lprobYZ  = NA
    }
  
  
  lprobZ <- Reduce("+", 
                   lapply(1:length(collecNetwork),
                          function(m){
                            dPostRow <- apply(collecZMAP[[m]]$row,2,sum) 
                            LRow  <- mylBetaFunction(dPostRow + priorParam$blockProp$row[m,]) - mylBetaFunction(priorParam$blockProp$row[m,])
                            dPostCol <- apply(collecZMAP[[m]]$col,2,sum) + priorParam$blockProp$col[m,]
                            LCol  <- mylBetaFunction(dPostCol + priorParam$blockProp$col[m,]) - mylBetaFunction(priorParam$blockProp$col[m,])
                            return(LRow  +  LCol)
                            }
                          )
                   ) 

  res <- c(lprobYZ,lprobZ)
  names(res) =c('lprobYZ','lprobZ')
  return(res)
  
}

#######################""
mylBetaFunction = function(alpha){
   sum(lgamma(alpha))-lgamma(sum(alpha))
}


################################## set prior param
setPriorParam <- function(M,KRow,KCol,model){
  # M  : number of bipartite networks
  # Krow  : number of blocks in row
  # Kcol  : number of blocks in col
  # model : poisson or bernoulli
  
  priorParam <- list()
  
  priorParam$blockProp <- list()
  priorParam$blockProp$row = matrix(1/KRow,M,KRow) ### Jeffreys dirichlet prior the pi^m and rho^m
  priorParam$blockProp$col = matrix(1/KCol,M,KCol)
  
  priorParam$connectParam<- list()
  if (model=='bernoulli'){
    priorParam$connectParam$alpha <- matrix(1,KRow,KCol) ### Uniform prior on the alpha_{kl}
    priorParam$connectParam$beta <- matrix(1,KRow,KCol)
  }
  if (model=='poisson'){
    priorParam$connectParam$alpha <- matrix(1,KRow,KCol) ### Exp de param 1/100
    priorParam$connectParam$beta <- matrix(1/100,KRow,KCol)
  }
  return(priorParam)
}





