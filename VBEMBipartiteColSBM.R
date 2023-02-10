VBEMBipartiteColSBM = function(collecNetwork,priorParam,collecTau,estimOptions)

  
  estimOptions <- list(maxIterVB = 100,
                     maxIterVE = 100,
                     valStopCritVE = 10^-5,
                     valStopCritVB = 10^-5)

  
  
  M <- length(collecNetwork)
  nRow <- sapply(collecNetwork,function(X)(nrow(X)))
  nCol <- sapply(collecNetwork,function(X)(ncol(X)))
  KRow <-ncol(priorParam$blockProp$row)
  KCol <-ncol(priorParam$blockProp$col)
  
  
  #-------------------- initialisation
  postParam <- Mstep(collecNetwork,M, nRow,nCol,collecTau,priorParam)
  noConvergence <- 0;
  iterVB <- 0
  stopVB <- 0 
  
  #-------------- RUN VBEM
  while (iterVB < estimOptions$maxIterVB & stopVB == 0){
  
    iterVB  <- iterVB  + 1
    postParamOld <- postParam
  
    ##--------------- VE step
    res_Estp <- EStep(collecNetwork,M, nRow,nCol,KRow,KCol,collecTau,postParam,estimOptions,model)
    collecTau <- res_Estp$collecTau
    #----- end  VB M step
    
    ##--------------- VB step
    postParam <- MStep(collecNetwork,M, nRow,nCol,collecTau,priorParam, distri = 'Bernoulli')
    #----- end  VB M step
    
    ##-------------- stop criteria
    if (distPostParam(postParam,postParamOld) < valStopCritVB) {stopVB <- 1}
    print(c(iterVB,deltaPostParam))
}


############################################################
####################" MSTEP from CollecTau to postParam
###############################################################
Mstep <- function(collecNetwork,M, nRow,nCol,collecTau,priorParam,model){
  
  
  postParam <- priorParam
  S <- lapply(1:M,function(m){t(collecTau[[m]]$row) %*% matrix(1,nRow[m],nCol[m]) %*% collecTau[[m]]$col})
  SY <- lapply(1:M,function(m){t(collecTau[[m]]$row) %*% collecNetwork[[m]] %*% collecTau[[m]]$col})
  
  postParam$connectParam$alpha  <- priorParam$connectParam$alpha +  Reduce(`+`, SY)
  postParam$connectParam$beta  <- priorParam$connectParam$beta +   Reduce(`+`, S) - Reduce(`+`, SY)
  
  Rtau <-  t(sapply(1:M,function(m){colSums(collecTau[[m]]$row)}))
  Ctau <-  t(sapply(1:M,function(m){colSums(collecTau[[m]]$col)}))
  
  postParam$blockProp$row  <- priorParam$blockProp$row + Rtau
  postParam$blockProp$col  <- priorParam$blockProp$col + Ctau
  
  return(postParam)
}
############################################################
####################" MSTEP from CollecTau to postParam
###############################################################
Estep <- function(collecNetwork,M, nRow,nCol,KRow,KCol,collecTau,postParam,estimOptions,model){
  
  iterVE <- 0
  stopVE <- 0
  
  while ((iterVE < estimOptions$maxIterVE) & (stopVE == 0)){
    collecTauOld <- collecTau #useful ?
    if(KRow == 1){for (m in 1:M){collecTau[[m]]$row <-  matrix(1,ncol  = 1,nrow = nRow[m])}}
    if(KCol == 1){for (m in 1:M){collecTau[[m]]$col <-  matrix(1,ncol  = 1,nrow = nCol[m])}}
    if((KCol>0) | (KRow>1)){
    ## useful quantities
      DiGAlpha <- digamma(postParam$connectParam$alpha) ### useful for Poisson and Bernoulli
      if(model == 'bernoulli'){
        DiGBeta <- digamma(postParam$connectParam$beta) 
        DiGAlphaBeta  <- digamma(postParam$connectParam$beta + postParam$connectParam$alpha) 
      }
      
      for(m in 1:M){
        DiblockProp_m_col <- digamma(postParam$blockProp$col[m,]) - digamma(sum(postParam$blockProp$col[m,]))
        DiblockProp_m_row <- digamma(postParam$blockProp$row[m,]) - digamma(sum(postParam$blockProp$row[m,]))
        # row
        if(model=='poisson'){
          lY_m  <- collecNetwork[[m]] %*% tcrossprod(collecTau[[m]]$col,log(postParam$connectParam$beta) + DiGAlpha)  
          lY_m <- lY_m - matrix(1,nRow[m], nCol[m])  %*% tcrossprod(collecTau[[m]]$col,postParam$connectParam$alpha/postParam$connectParam$beta)
        }
        if(model=='bernoulli'){
          lY_m  <- collecNetwork[[m]] %*% tcrossprod(collecTau[[m]]$col,DiGAlpha- DiGAlphaBeta)  +  (1-collecNetwork[[m]]) %*% tcrossprod(collecTau[[m]]$col,DiGBeta- DiGAlphaBeta)
        }
        l3_m  <- matrix(DiblockProp_m_row,nrow = nRow[m],ncol = KRow,byrow = TRUE)
        collecTau[[m]]$row <- fromBtoTau(lY_m + l3_m) 
        
        
        # col
        if(model=='poisson'){
          lY_m  <- t(collecNetwork[[m]]) %*% tcrossprod(collecTau[[m]]$row,t(log(postParam$connectParam$beta) + DiGAlpha))  
          lY_m <- lY_m - matrix(1, nCol[m],nCol[m])  %*% tcrossprod(collecTau[[m]]$row,t(postParam$connectParam$alpha/postParam$connectParam$beta))
        }
        if(model=='bernoulli'){
          lY_m <-  t(collecNetwork[[m]])  %*% tcrossprod(collecTau[[m]]$row,t(DiGAlpha- DiGAlphaBeta))  + t(1-collecNetwork[[m]]) %*% tcrossprod(collecTau[[m]]$row,t(DiGBeta- DiGAlphaBeta))
        }
        l3_m  <- matrix(DiblockProp_m_row,nrow = nCol[m],ncol = KCol,byrow = TRUE)
        collecTau[[m]]$col <- fromBtoTau(lY_m + l3_m) 
      }
    }

    deltaTau <- distTau(collecTau,collecTauOld)
    if (deltaTau < estimOptions$valStopCritVE) {stopVE <- 1}
    iterVE <- iterVE + 1
    if(iterVE  == maxIterVE){noConvergence = 1}
  }
  
  return(list(collecTau = collecTau,noConvergence  = noConvergence))
}

