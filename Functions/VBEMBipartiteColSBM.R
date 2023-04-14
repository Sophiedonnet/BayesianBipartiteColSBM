VBEMBipartiteColSBM = function(collecNetworks,hyperparamPrior,collecTau,estimOptions, emissionDist, model ){

  
  if(is.null(estimOptions)){
    estimOptions <- list(maxIterVB = 100,
                        maxIterVE = 100,
                        valStopCritVE = 10^-5,
                        valStopCritVB = 10^-5)
  } 

  
  
  
  M <- length(collecNetworks)
  nRow <- sapply(collecNetworks,function(X)(nrow(X)))
  nCol <- sapply(collecNetworks,function(X)(ncol(X)))

  KRow <-nrow(hyperparamPrior$connectParam$alpha)
  KCol <-ncol(hyperparamPrior$connectParam$alpha)
  
  #-------------------- initialisation
  hyperparamPost <- Mstep(collecNetworks,M, nRow,nCol,collecTau,hyperparamPrior, emissionDist, model)
  noConvergence <- 0;
  iterVB <- 0
  stopVB <- 0 
  
  #-------------- RUN VBEM
  while (iterVB < estimOptions$maxIterVB & stopVB == 0){
  
    iterVB  <- iterVB  + 1
    hyperparamPostOld <- hyperparamPost
  
    ##--------------- VE step
    res_Estep <- Estep(collecNetworks,M, nRow,nCol,KRow,KCol,collecTau,hyperparamPost,estimOptions,emissionDist, model)
    collecTau <- res_Estep$collecTau
    #----- end  VB M step
    
    ##--------------- VB step
    hyperparamPost <- Mstep(collecNetworks,M, nRow,nCol,collecTau,hyperparamPrior, emissionDist, model)
    #----- end  VB M step
    
    ##-------------- stop criteria
    if (disthyperparamPost(hyperparamPost,hyperparamPostOld) < estimOptions$valStopCritVB) {stopVB <- 1}
    print(iterVB)
  }
  
  return(reorderBlocks(hyperparamPost, collecTau, model))
    

}


############################################################
####################" MSTEP from CollecTau to hyperparamPost
###############################################################
Mstep <- function(collecNetworks,M, nRow,nCol,collecTau,hyperparamPrior,emissionDist, model){
  
  
  hyperparamPost <- hyperparamPrior
  S <- lapply(1:M,function(m){t(collecTau[[m]]$row) %*% matrix(1,nRow[m],nCol[m]) %*% collecTau[[m]]$col})
  SY <- lapply(1:M,function(m){t(collecTau[[m]]$row) %*% collecNetworks[[m]] %*% collecTau[[m]]$col})
  
  
  hyperparamPost$connectParam$alpha  <- hyperparamPrior$connectParam$alpha +  Reduce(`+`, SY)
  hyperparamPost$connectParam$beta  <- hyperparamPrior$connectParam$beta +   Reduce(`+`, S) - (emissionDist == 'bernoulli') * Reduce(`+`, SY)
  
  Rtau <-  t(sapply(1:M,function(m){colSums(collecTau[[m]]$row)}))
  Ctau <-  t(sapply(1:M,function(m){colSums(collecTau[[m]]$col)}))
  
  if(model == 'iidColBipartiteSBM'){  
    hyperparamPost$blockProp$row  <- hyperparamPrior$blockProp$row + apply(Rtau,2,sum)
    hyperparamPost$blockProp$col  <- hyperparamPrior$blockProp$col + apply(Ctau,2,sum)
  }
  if(model == 'piColBipartiteSBM'){
    hyperparamPost$blockProp$row  <- hyperparamPrior$blockProp$row + Rtau
    hyperparamPost$blockProp$col  <- hyperparamPrior$blockProp$col + Ctau
  }
  
  
  return(hyperparamPost)
}
############################################################
####################" MSTEP from CollecTau to hyperparamPost
###############################################################
Estep <- function(collecNetworks,M, nRow,nCol,KRow,KCol,collecTau,hyperparamPost,estimOptions,emissionDist, model){
  
  iterVE <- 0
  stopVE <- 0
  
  while ((iterVE < estimOptions$maxIterVE) & (stopVE == 0)){
    collecTauOld <- collecTau #useful ?
    if(KRow == 1){for (m in 1:M){collecTau[[m]]$row <-  matrix(1,ncol  = 1,nrow = nRow[m])}}
    if(KCol == 1){for (m in 1:M){collecTau[[m]]$col <-  matrix(1,ncol  = 1,nrow = nCol[m])}}
    if((KCol>0) | (KRow>1)){
    ## useful quantities
      DiGAlpha <- digamma(hyperparamPost$connectParam$alpha) ### useful for Poisson and Bernoulli
      if(emissionDist == 'bernoulli'){
        DiGBeta <- digamma(hyperparamPost$connectParam$beta) 
        DiGAlphaBeta  <- digamma(hyperparamPost$connectParam$beta + hyperparamPost$connectParam$alpha) 
      }
      
      for(m in 1:M){
        if(model == 'iidColBipartiteSBM'){
          DiblockProp_m_col <- digamma(hyperparamPost$blockProp$col) - digamma(sum(hyperparamPost$blockProp$col))
          DiblockProp_m_row <- digamma(hyperparamPost$blockProp$row) - digamma(sum(hyperparamPost$blockProp$row))
        }
        if(model == 'piColBipartiteSBM'){
          DiblockProp_m_col <- digamma(hyperparamPost$blockProp$col[m,]) - digamma(sum(hyperparamPost$blockProp$col[m,]))
          DiblockProp_m_row <- digamma(hyperparamPost$blockProp$row[m,]) - digamma(sum(hyperparamPost$blockProp$row[m,]))
        }
        

        # row
        if(emissionDist == 'poisson'){
          lY_m  <- collecNetworks[[m]] %*% tcrossprod(collecTau[[m]]$col,log(hyperparamPost$connectParam$beta) + DiGAlpha)  
          lY_m <- lY_m - matrix(1,nRow[m], nCol[m])  %*% tcrossprod(collecTau[[m]]$col,hyperparamPost$connectParam$alpha/hyperparamPost$connectParam$beta)
        }
        if(emissionDist == 'bernoulli'){
          lY_m  <- collecNetworks[[m]] %*% tcrossprod(collecTau[[m]]$col,DiGAlpha- DiGAlphaBeta)  +  (1-collecNetworks[[m]]) %*% tcrossprod(collecTau[[m]]$col,DiGBeta- DiGAlphaBeta)
        }
        l3_m  <- matrix(DiblockProp_m_row,nrow = nRow[m],ncol = KRow,byrow = TRUE)
        collecTau[[m]]$row <- fromBtoTau(lY_m + l3_m) 
        
        
        # col
        if(emissionDist == 'poisson'){
          lY_m  <- t(collecNetworks[[m]]) %*% tcrossprod(collecTau[[m]]$row,t(log(hyperparamPost$connectParam$beta) + DiGAlpha))  
          lY_m <- lY_m - matrix(1, nCol[m],nCol[m])  %*% tcrossprod(collecTau[[m]]$row,t(hyperparamPost$connectParam$alpha/hyperparamPost$connectParam$beta))
        }
        if(emissionDist == 'bernoulli'){
          lY_m <-  t(collecNetworks[[m]])  %*% tcrossprod(collecTau[[m]]$row,t(DiGAlpha- DiGAlphaBeta))  + t(1-collecNetworks[[m]]) %*% tcrossprod(collecTau[[m]]$row,t(DiGBeta- DiGAlphaBeta))
        }
        l3_m  <- matrix(DiblockProp_m_col,nrow = nCol[m],ncol = KCol,byrow = TRUE)
        collecTau[[m]]$col <- fromBtoTau(lY_m + l3_m) 
      }
    }

    deltaTau <- distTau(collecTau,collecTauOld)
    if (deltaTau < estimOptions$valStopCritVE) {stopVE <- 1}
    iterVE <- iterVE + 1
    noConvergence <- 1*(iterVE  == estimOptions$maxIterVE)
  }
  
  return(list(collecTau = collecTau,noConvergence  = noConvergence))
}


############################## log marg likelihood
computeLogLikMarg_VB <- function(collecNetworks, collecTau, hyperparamPrior, emissionDist, model){
  
  collecZMAP <- lapply(collecTau,function(tau){
    indZ <- tau
    indZ$row <- t(sapply(1:nrow(tau$row),function(i){u <- 0*tau$row[i,]; u[which.max(tau$row[i,])]=1; return(u)}))
    if(min(dim(indZ$row))==1){indZ$row = matrix(indZ$row,ncol=1)}
    indZ$col <- t(sapply(1:nrow(tau$col),function(i){u <- 0*tau$col[i,]; u[which.max(tau$col[i,])]=1; return(u)}))
    if(min(dim(indZ$col))==1){indZ$col = matrix(indZ$col,ncol=1)}
    return(indZ)
  })
  
  if(emissionDist =='bernoulli'){
    S1 <- lapply(1:length(collecNetworks),function(m){t(collecZMAP[[m]]$row)%*% collecNetworks[[m]] %*% collecZMAP[[m]]$col})
    S0 <- lapply(1:length(collecNetworks),function(m){t(collecZMAP[[m]]$row)%*% (1-collecNetworks[[m]]) %*% collecZMAP[[m]]$col})
    
    ahat <- Reduce('+',S1)
    bhat <- Reduce('+',S0)
    lprobYZ <- sum(lbeta(ahat+hyperparamPrior$connectParam$alpha,bhat+hyperparamPrior$connectParam$beta) -lbeta(hyperparamPrior$connectParam$alpha,hyperparamPrior$connectParam$beta))
  }
  
  if(model=='poisson'){
    lprobYZ  = NA
  }
  
  if(model == 'piColBipartiteSBM'){
    lprobZ <- Reduce("+", 
                     lapply(1:length(collecNetworks),
                            function(m){
                              dPostRow <- apply(collecZMAP[[m]]$row,2,sum) 
                              LRow  <- mylBetaFunction(dPostRow + hyperparamPrior$blockProp$row[m,]) - mylBetaFunction(hyperparamPrior$blockProp$row[m,])
                              dPostCol <- apply(collecZMAP[[m]]$col,2,sum) 
                              LCol  <- mylBetaFunction(dPostCol + hyperparamPrior$blockProp$col[m,]) - mylBetaFunction(hyperparamPrior$blockProp$col[m,])
                              return(LRow  +  LCol)
                            }
                     )
    )
  }
  if(model == 'iidColBipartiteSBM'){
    dPostRow <- rowSums(sapply(1:length(collecNetworks),function(m){apply(collecZMAP[[m]]$row,2,sum)}))
    dPostCol <- rowSums(sapply(1:length(collecNetworks),function(m){apply(collecZMAP[[m]]$col,2,sum)}))
    LRow  <- mylBetaFunction(dPostRow + hyperparamPrior$blockProp$row) - mylBetaFunction(hyperparamPrior$blockProp$row)
    LCol  <- mylBetaFunction(dPostCol + hyperparamPrior$blockProp$col) - mylBetaFunction(hyperparamPrior$blockProp$col)
    lprobZ <-  LRow  +  LCol
  }
  
  
  res <- c(lprobYZ,lprobZ)
  names(res) =c('lprobYZ','lprobZ')
  return(res)
  
}

