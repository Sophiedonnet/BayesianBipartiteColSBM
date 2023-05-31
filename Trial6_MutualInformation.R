rm(list=ls())
#setwd('~/WORK_ALL/RECHERCHE/TRAVAUX_RECHERCHE/Barbillon-Chabert/BayesianBipartiteColSBM/')
if(Sys.info()[[4]]=='sophie-Latitude-5310'){
  setwd('/home/sophie/WORK_LOCAL/RECHERCHE/TRAVAUX_DE_RECHERCHE/Barbillon-Chabert/BayesianBipartiteColSBM')
}

library(sbm)
library(gtools)
sapply(list.files(paste0(getwd(),'/Functions'),full.names = TRUE), source)





###############################################################
######## MODELS 
###############################################################

emissionDist = 'bernoulli'
model = 'piColBipartiteSBM'


###########################################################################################
############# simulation 
#############################################################################################
M = 3
KRow = 4
KCol = 3


#-----  block proportions simul piColSBM avec classes vides
mySeed = sample(1:10000,1)
whereSaveData <- paste0(getwd(),'/Simu/Data/',model,'/data_simu_',mySeed,'.Rdata')
print(mySeed)
set.seed(mySeed)
blockPropTrue <- list()
blockPropTrue$row <-  rdirichlet(M,rep(1/1.2,KRow))  #### emptying some blocks in certain netwokrs
blockPropTrue$col <- rdirichlet(M,rep(1/1.2,KCol))   #### emptying some blocks in certain netwokros
blockPropTrue$row[1,] <- rep(1/KRow,KRow)
blockPropTrue$col[1,] <- rep(1/KCol,KCol)

print(blockPropTrue)

#-----  connectivity matrix
connectParamTrue <- list(mean = round( matrix(rbeta(KRow*KCol,1/1.01,1/1.01), KRow, KCol)  ,9))
oCol <- order(colSums(connectParamTrue$mean[,]),decreasing = TRUE)
oRow <- order(rowSums(connectParamTrue$mean[,]),decreasing = TRUE)
connectParamTrue$mean <- connectParamTrue$mean[oRow,oCol]

#-------  sizes of networks
nbNodes <- matrix(sample(20*c(5:10),M*2,replace = TRUE),M,2)

#nbNodes[1,] = c(150,200)
#nbNodes[5,] = c(30,20)

#---------  SIMULATION
mySampler <- lapply(1:M, function(m){sampleBipartiteSBM(nbNodes = nbNodes[m,],  list(row=blockPropTrue$row[m,],col=blockPropTrue$col[m,]),  connectParamTrue)})
collecNetworks <- lapply(mySampler,function(l){l$networkData})
mydata <- list(collecNetworks = collecNetworks, M= M, nbNodes = nbNodes)


##############################################################
#---------- load data
##############################################################
 


#--------------- init CollecTau
initBipartite <- lapply(collecNetworks,estimateBipartiteSBM)
myRef <- which.max(rowSums(t(sapply(initBipartite,function(m){m$nbBlocks}))))
collecTau_init <- initCollecTau(initBipartite,ref = myRef)
KRow <- initBipartite[[myRef]]$nbBlocks[1]
KCol <- initBipartite[[myRef]]$nbBlocks[2]
print(c(KRow,KCol))


###################################################### 
#---------------  Prior distribution
######################################################
hyperparamPrior_piCol <- setHyperparamPrior(M,KRow,KCol, emissionDist, model)

###########################################################################################
#------------------ VBEM  + SMC -VBEM 
###########################################################################################


#------------- variational estim 

estimOptionsVBEM <- list(maxIterVB = 1000,
                         maxIterVE = 100,
                         valStopCritVE = 10^-5,
                         valStopCritVB = 10^-5)

resEstimVBEM_piCol  <- VBEMBipartiteColSBM(collecNetworks,hyperparamPrior_piCol,collecTau_init,estimOptions = estimOptionsVBEM, emissionDist, model='piColBipartiteSBM')
LogLikMarg_VB_piCol <- sum(computeLogLikMarg_VB(collecNetworks,resEstimVBEM_piCol$collecTau,hyperparamPrior_piCol,emissionDist,model))



whereSaveResEstimVB <- paste0(getwd(),'/Simu/Res/',model,'/VB/resVB_',mySeed,'.Rdata')
#save(initBipartite,estimOptionsVBEM, resEstimVBEM_piCol,estimOptionsVBEM,hyperparamPrior_piCol,file=whereSaveResEstimVB)

#------------------ Set ApproxPost
hyperparamApproxPost_piCol <- resEstimVBEM_piCol$hyperparamPost
hyperparamApproxPost_piCol$collecTau <- resEstimVBEM_piCol$collecTau

#------------------  SMC  - VEM


HSample <- rParamZ(MC=1000, hyperparamApproxPost_piCol, emissionDist, model,nbNodes);
ZSample <-  HSample$ZSample
MI <- mutualInformationZ(list(HSample$ZSample[[3]]))
MI

estimOptionsSMC = list()
estimOptionsSMC$paramsMCMC <- list(nbIterMCMC=2)  
estimOptionsSMC$MC <- 1000
estimOptionsSMC$ESS.rate <- 0.9
estimOptionsSMC$cESS.rate <- 0.9
estimOptionsSMC$opSave <- FALSE
estimOptionsSMC$op.parallel  <- list(os =  .Platform$OS.type, mc.cores = 1); 
estimOptionsSMC$op.print<- TRUE
estimOptionsSMC$NB.iter.max  <- Inf # Inf
estimOptionsSMC$op.SMC.classic <- FALSE

resSMC_VB_piCol <- SMCColBipartiteSBM(data = mydata,hyperparamPrior_piCol,hyperparamApproxPost_piCol, emissionDist, model , estimOptionsSMC)

