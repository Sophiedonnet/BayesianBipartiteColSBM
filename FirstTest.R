rm(list=ls())
library(sbm)
source('utils.R')
source('VBEMBipartiteColSBM.R')
library(gtools)



seed  <- 14579
set.seed(seed)


M = 6
KRow = 4
KCol = 3


########### block proportions
blockProp <- list()
blockProp$row <-  rdirichlet(M,rep(1/(KRow-1),KRow)) 
blockProp$col <- rdirichlet(M,rep(1/(KCol-1),KCol))

blockProp$row[1,] <- rep(1/KRow,KRow) 
blockProp$col[1,] <- rep(1/KCol,KCol)

print(blockProp)

############## connectivity matrix
means <- matrix(rbeta(KRow*KCol,1/1.1,1/1.1), KRow, KCol)  
connectParam <- list(mean = round(means,2))
oCol <- order(colSums(connectParam$mean[,]),decreasing = TRUE)
oRow <- order(rowSums(connectParam$mean[,]),decreasing = TRUE)
connectParam$mean <- connectParam$mean[oRow,oCol]


########### sizes of networks
nbNodes <- matrix(sample(10*c(6:10),M*2,replace = TRUE),M,2)


############################################################
##################################"""""""" SIMULATION
############################################################


mySampler <- lapply(1:M, function(m){sampleBipartiteSBM(nbNodes = nbNodes[m,],  list(row=blockProp$row[m,],col=blockProp$col[m,]),  connectParam)})
collecNetwork <- lapply(mySampler,function(l){l$networkData})





############### init Tau

initBipartite <- lapply(collecNetwork,estimateBipartiteSBM)









initBipartite[[1]]$nbBlocks

###################### Estim avec Kcol Krow true


estimOptions <- list(maxIterVB = 1000,
                     maxIterVE = 100,
                     valStopCritVE = 10^-5,
                     valStopCritVB = 10^-5)


KRowEstim <- initBipartite[[1]]$nbBlocks[1]
KColEstim <- initBipartite[[1]]$nbBlocks[2]

#   Making the groups correspond between networks

collecTau <- lapply(1:M,function(m){
    eps  <-  10^-5
    alpha_1 <- initBipartite[[1]]$connectParam$mean
    KRowEstim_1 <- nrow(alpha_1)
    KColEstim_1 <- ncol(alpha_1)
    
    alpha_m <- initBipartite[[m]]$connectParam$mean
    Perm_row <- permutations(KRowEstim_1,nrow(alpha_m))
    Perm_col <- permutations(KColEstim_1,ncol(alpha_m))
    score  = Inf
    for (i in 1:nrow(Perm_row)){
      for (j in 1:nrow(Perm_col)){
        score_new = sum((alpha_1[Perm_row[i,],Perm_col[j,]] - alpha_m)^2)
        if (score_new < score){
          best_perm_row <- Perm_row[i,]
          best_perm_col <- Perm_col[j,]
          score <- score_new
        }
      }
    }
  tau_m <- list()
  tau_m$row <- matrix(eps,nbNodes[m,1],KRowEstim_1)
  tau_m$row[ , best_perm_row] = initBipartite[[m]]$probMemberships$row
  tau_m$row <- normByRow(tau_m$row)
  tau_m$col <- matrix(eps,nbNodes[m,2],KColEstim_1)
  tau_m$col[ , best_perm_col] = initBipartite[[m]]$probMemberships$col
  tau_m$col <- normByRow(tau_m$col)
  return(tau_m)}
)






priorParam <- setPriorParam(M,KRowEstim,KColEstim,model ='bernoulli')
resEstim  <- VBEMBipartiteColSBM(collecNetwork,priorParam,collecTau,estimOptions,model  ='bernoulli')
logLikMarg <- computeLogLikMarg(collecNetwork, resEstim$collecTau,priorParam,model = 'bernoulli')


postMeanEstim <- list()
postMeanEstim$connectParam <- resEstim$postParam$connectParam$alpha/(resEstim$postParam$connectParam$alpha + resEstim$postParam$connectParam$beta)
postMeanEstim$blockProp <- lapply(resEstim$postParam$blockProp,function(l){normByRow(l)}) 
postMemberships <- lapply(collecTau,function(tau){
  Z <- list()
  Z$row <- apply(tau$row,1,which.max)
  Z$col <- apply(tau$col,1,which.max)
  return(Z)
  })


lapply(1:M,function(m){table(postMemberships[[m]]$row,mySampler[[m]]$memberships$row)})

lapply(1:M,function(m){table(postMemberships[[m]]$col,mySampler[[m]]$memberships$col)})

 
################ sepSBM 


logLikMarg_sep <- matrix(0,M,2)
for (m in 1:M){
  print(paste0('Network ',m))
  KRow_m <- initBipartite[[m]]$nbBlocks[1]
  KCol_m <- initBipartite[[m]]$nbBlocks[2]
  tau_m <- initBipartite[[m]]$probMemberships
  priorParam_m <- setPriorParam(1, KRow_m,KCol_m,model ='bernoulli')
  resEstim_m  <- VBEMBipartiteColSBM(list(collecNetwork[[m]]),priorParam_m,list(tau_m),estimOptions,model  ='bernoulli')
  logLikMarg_sep[m,] <-  computeLogLikMarg(list(collecNetwork[[m]]), list(tau_m),priorParam_m,model = 'bernoulli')
}



cbind(sum(logLikMarg),sum(logLikMarg_sep))
    

muPost <-   postParam$connectParam$alpha/(postParam$connectParam$beta + postParam$connectParam$alpha)
oCol <- order(colSums(muPost),decreasing = TRUE)
oRow <- order(rowSums(muPost),decreasing = TRUE)
postParam$connectParam$alpha <-postParam$connectParam$alpha[oRow,oCol]
postParam$connectParam$beta <-postParam$connectParam$beta[oRow,oCol]

postParam$blockProp$row <- postParam$blockProp$row[,oRow]
postParam$blockProp$cow <- postParam$blockProp$cow[,oCol]

muPost <-   postParam$connectParam$alpha/(postParam$connectParam$beta + postParam$connectParam$alpha)

############## Plot post 
par(mfrow = c(KRow,KCol))
for (k in 1:KRow){
  for (l in 1:KCol){
  curve(dbeta(x,postParam$connectParam$alpha[k,l],postParam$connectParam$beta[k,l]),ylab = 'post')
  abline(v=connectParam$mean[k,l],col='red')
  abline(v=muPost[k,l],col='green')
  }
}
  