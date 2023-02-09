library(sbm)
source('utils.R')
nbNodes <- 
blockProp <- list(c(.5, .5), c(1/3, 1/3, 1/3)) # group proportions
blockProp2 <- list(c(.3, .7), c(0.5, 0.5,0.000001))
seed  <- 14579
set.seed(seed)
means <- matrix(runif(6), 2, 3)  # connectivity matrix
connectParam <- list(mean = round(means,2))
mySampler1 <- sampleBipartiteSBM(nbNodes = c(60, 80),  blockProp,  connectParam)
mySampler2 <- sampleBipartiteSBM(nbNodes = c(50, 100), blockProp2, connectParam)






##################"DATA
collecNetwork <- list()
collecNetwork[[1]] <- mySampler1$networkData
collecNetwork[[2]] <- mySampler2$networkData

# from data
M <- length(collecNetwork)
nRow <- sapply(collecNetwork,function(m){nrow(m)})
nCol <- sapply(collecNetwork,function(m){ncol(m)})

########## MODEL and prior
KRow <- 2 
KCol <- 3
priorParam <- list()
priorParam$connectParam<- list()
priorParam$connectParam$alpha <- matrix(1,KRow,KCol) ### Uniform prior on the alpha_{kl}
priorParam$connectParam$beta <- matrix(1,KRow,KCol)
priorParam$blockProp <- list()
priorParam$blockProp$row = matrix(1/KRow,M,KRow) ### Jeffreys dirichlet prior the pi^m and rho^m
priorParam$blockProp$col = matrix(1/KCol,M,KCol)

############ Post  
postParam <- priorParam 
postParam$blockProp$row = matrix(NA,M,KRow)
postParam$blockProp$col = matrix(NA,M,KCol)
postParam$connectParam$alpha <- matrix(NA,KRow,KCol)
postParam$connectParam$beta <- matrix(NA,KRow,KCol)


############### init Tau

initBipartite <- lapply(collecNetwork,estimateBipartiteSBM)
# Init tau 
collecTau <- lapply(initBipartite,function(m){
  eps  <-  10^-5
  tau <- m$probMemberships
  kr <- ncol(tau$row)
  if(kr < KRow){tau$row = cbind(tau$row,matrix(eps,nrow(tau$row),KRow-kr)); tau$row =normByRow(tau$row) }
  kc <- ncol(tau$col)
  if(kc < KCol){tau$col = cbind(tau$col,matrix(eps,nrow(tau$col),KCol-kc));tau$col =normByRow(tau$col) }
  return(tau)
  })
 
postParamInit <- Mstep(collecNetwork,M, nRow,nCol,collecTau,priorParam)


noConvergence <- 0;
maxIterVB <- 100
maxIterVE <- 100
iterVB <- 0
stopVB <- 0 
valStopCritVE <- 10^-5
valStopCritVB <- 10^-5
postParam<- postParamInit

while (iterVB < maxIterVB & stopVB == 0){
  
  iterVB  <- iterVB  + 1
  postParamOld <- postParam
  
  ############### #VE step

  iterVE <- 0
  stopVE <- 0
  while ((iterVE < maxIterVE) & (stopVE == 0)){
    collecTauOld <- collecTau #useful ?
    if(KRow == 1){for (m in 1:M){collecTau[[m]]$row <-  matrix(1,ncol  = 1,nrow = nRow[m])}}
    if(KCol == 1){for (m in 1:M){collecTau[[m]]$col <-  matrix(1,ncol  = 1,nrow = nCol[m])}}
    if((KCol>0) | (KRow>1)){
      ## useful quantities
      DiGAlpha <- digamma(postParam$connectParam$alpha) 
      DiGBeta <- digamma(postParam$connectParam$beta) 
      DiGAlphaBeta  <- digamma(postParam$connectParam$beta + postParam$connectParam$alpha) 
      for(m in 1:M){
        DiblockProp_m_col <- digamma(postParam$blockProp$col[m,]) - digamma(sum(postParam$blockProp$col[m,]))
        DiblockProp_m_row <- digamma(postParam$blockProp$row[m,]) - digamma(sum(postParam$blockProp$row[m,]))
        # row
        l1_m  <- collecNetwork[[m]] %*% tcrossprod(collecTau[[m]]$col,DiGAlpha- DiGAlphaBeta) 
        l2_m  <- (1-collecNetwork[[m]]) %*% tcrossprod(collecTau[[m]]$col,DiGBeta- DiGAlphaBeta)
        l3_m  <- matrix(DiblockProp_m_row,nrow = nRow[m],ncol = KRow,byrow = TRUE)
        collecTau[[m]]$row <- fromBtoTau(l1_m + l2_m + l3_m) 
        # col 
        l1_m  <-  t(collecNetwork[[m]])  %*% tcrossprod(collecTau[[m]]$row,t(DiGAlpha- DiGAlphaBeta))  
        l2_m  <- t(1-collecNetwork[[m]]) %*% tcrossprod(collecTau[[m]]$row,t(DiGBeta- DiGAlphaBeta))
        l3_m  <- matrix(DiblockProp_m_row,nrow = nCol[m],ncol = KCol,byrow = TRUE)
        collecTau[[m]]$col <- fromBtoTau(l1_m + l2_m + l3_m) 
        }
    deltaTau <- distTau(collecTau,collecTauOld)
    if (deltaTau < valStopCritVE) {stopVE <- 1}
    iterVE <- iterVE + 1
    if(iterVE  == maxIterVE){noConvergence = 1}
    }
  }### end VE
  
  ###### VB M step
  postParam<- Mstep(collecNetwork,M, nRow,nCol,collecTau,priorParam)
   
  #end  VB M step
  #-------stop criteria
  deltaPostParam <- distPostParam(postParam,postParamOld)
  if (deltaPostParam < valStopCritVB) {stopVB <- 1}
  print(c(iterVB,deltaPostParam))
}
    
  