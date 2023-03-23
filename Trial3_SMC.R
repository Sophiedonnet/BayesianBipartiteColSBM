rm(list=ls())
library(sbm)
library(gtools)

source('Functions/utils.R')
source('Functions/VBEMBipartiteColSBM.R')
source('Functions/initializationCollecTau.R')
source('Functions/Functions-SMC.R')
source('Functions/Functions-colBipartiteSBM.R')


emissionDist = 'bernoulli'
model = 'piColBipartiteSBM'

###########################################################################################
############# simulation avec iid colSBM mais certains très petits réseaux. 
#############################################################################################


M = 6
KRow = 4
KCol = 3
#---  block proportions simul iidColSBM 
blockProp <- list()
blockProp$row <-  matrix(rdirichlet(1,rep(KRow,KRow)),byrow = TRUE,nrow = M, ncol = KRow) #### emptying some blocks in certain netwokrs 
blockProp$col <- matrix(rdirichlet(1,rep(KCol,KCol)),byrow = TRUE,nrow = M, ncol = KCol)   #### emptying some blocks in certain netwokros
print(blockProp)

#-----  connectivity matrix
means <- matrix(rbeta(KRow*KCol,1/1.1,1/1.1), KRow, KCol)  
connectParam <- list(mean = round(means,2))
oCol <- order(colSums(connectParam$mean[,]),decreasing = TRUE)
oRow <- order(rowSums(connectParam$mean[,]),decreasing = TRUE)
connectParam$mean <- connectParam$mean[oRow,oCol]


#-------  sizes of networks
nbNodes <- matrix(sample(10*c(6:10),M*2,replace = TRUE),M,2)
nbNodes[5,] = c(30,20)
nbNodes[6,] = c(20,30)

#---------  SIMULATION
mySampler <- lapply(1:M, function(m){sampleBipartiteSBM(nbNodes = nbNodes[m,],  list(row=blockProp$row[m,],col=blockProp$col[m,]),  connectParam)})
collecNetworks <- lapply(mySampler,function(l){l$networkData})



############### init CollecTau

initBipartite <- lapply(collecNetworks,estimateBipartiteSBM)
myRef = which.max(rowSums(t(sapply(initBipartite,function(m){m$nbBlocks}))))
collecTau <- initCollecTau(initBipartite,ref = myRef)
KRow <- initBipartite[[myRef]]$nbBlocks[1]
KCol <- initBipartite[[myRef]]$nbBlocks[2]


################################ Prior distribution

hyperparamPrior <- setHyperparamPrior(M,KRow,KCol, emissionDist, model)
hyperparam <- hyperparamPrior
MC = 3000
Hsample <- rParamZ(MC, hyperparam, emissionDist, model ,nbNodes, collecTau = NULL)


# 
checkSample = function(Hsample,MC,M,KRow,KCol,hyperparam,collecTau,emissionDist,model){
  
  blockPropSample <- Hsample$blockPropSample
  ZSample <- Hsample$ZSample
  connectParamSample <- Hsample$connectParamSample  
  
  #---------------  connectParamSample
  print('dim of connectParamSample')
  print(sum(dim(connectParamSample)==c(KRow,KCol,MC))==3)
  
  expected_mean <- hyperparam$connectParam$alpha / (hyperparam$connectParam$beta + (emissionDist == 'bernoulli')*hyperparam$connectParam$alpha)
  print(cbind(expected_mean,apply(connectParamSample,c(1,2),mean)))
  
  #--------------- 
  print('dim of blockPropSample')
  if(model == 'iidColBipartiteSBM'){
    print(dim(blockPropSample$row)==c(KRow,MC))
    print(dim(blockPropSample$col)==c(KCol,MC))
  }
  if(model == 'piColBipartiteSBM'){
    print(dim(blockPropSample$row)==c(M,KRow,MC))
    print(dim(blockPropSample$col)==c(M,KCol,MC))
  }
  #------------------------------- 
  
  
  
  
  
  
}
##################################


resEstim_iid  <- VBEMBipartiteColSBM(collecNetworks,hyperparamPrior_iid,collecTau,estimOptions, emissionDist  = 'bernoulli' , model ='iidColBipartiteSBM')
logLikMarg_iid <- computeLogLikMarg(collecNetworks, resEstim_iid$collecTau,hyperparamPrior_iid, emissionDist  = 'bernoulli' , model ='iidColBipartiteSBM')






###################### Estim avec Kcol Krow true


estimOptions <- list(maxIterVB = 1000,
                     maxIterVE = 100,
                     valStopCritVE = 10^-5,
                     valStopCritVB = 10^-5)







hyperparamPrior_picol <- setHyperparamPrior(M,KRowEstim,KColEstim, emissionDist  = 'bernoulli' , model ='piColBipartiteSBM')
resEstim_picol  <- VBEMBipartiteColSBM(collecNetworks,hyperparamPrior_picol,collecTau,estimOptions, emissionDist  = 'bernoulli' , model ='piColBipartiteSBM')
logLikMarg_picol <- computeLogLikMarg(collecNetworks, resEstim_picol$collecTau,hyperparamPrior_picol, emissionDist  = 'bernoulli' , model ='piColBipartiteSBM')
