
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
  temp2 <- scaleByRow(temp)
  temp2[temp2 < eps] <- eps
  temp2[temp2 > (1 - eps)] <- 1 - eps
  Tau <- scaleByRow(temp2)
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

distPostParam = function(PostParam,PostParamOld){
  
  d1 <- sqrt(sum(as.vector((postParam$connectParam$alpha - postParamOld$connectParam$alpha)^2)))
  d2 <- sqrt(sum(as.vector((postParam$connectParam$beta -  postParamOld$connectParam$beta )^2)))
  d3 <- sqrt(sum(as.vector((postParam$blockProp$row -  postParamOld$blockProp$row )^2)))
  d4 <- sqrt(sum(as.vector((postParam$blockProp$col -  postParamOld$blockProp$col )^2)))
  return(d1 + d2 + d3 +d4)
}


####################" MSTEP from CollecTau to postParam

Mstep <- function(collecNetwork,M, nRow,nCol,collecTau,priorParam){
  
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
