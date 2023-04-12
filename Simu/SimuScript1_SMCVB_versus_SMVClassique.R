rm(list=ls())
setwd('~/WORK_ALL/RECHERCHE/TRAVAUX_RECHERCHE/Barbillon-Chabert/BayesianBipartiteColSBM/')

library(sbm)
library(gtools)
sapply(list.files(paste0(getwd(),'/Functions'),full.names = TRUE), source)


where_save_res <- paste0(getwd(),'/Simu/RES_SIMU1/')

###############################################################
######## MODELS 
###############################################################

emissionDist = 'bernoulli'
model = 'piColBipartiteSBM'

mySeed <- sample(1:999,1)
#mySeed <- 460
set.seed(mySeed)
###########################################################################################
############# simulation 
#############################################################################################
M = 3
KRow = 4
KCol = 3


#-----  block proportions simul piColSBM avec classes vides
if(model=='piColBipartiteSBM'){
  blockPropTrue <- list()
  blockPropTrue$row <-  rdirichlet(M,rep(1/(KRow-1),KRow))  #### emptying some blocks in certain netwokrs 
  blockPropTrue$col <- rdirichlet(M,rep(1/(KCol-1),KCol))   #### emptying some blocks in certain netwokros
  blockPropTrue$row[1,] <- rep(1/KRow,KRow) 
  blockPropTrue$col[1,] <- rep(1/KCol,KCol)
}

#-----  connectivity matrix
connectParamTrue <- list(mean = round( matrix(rbeta(KRow*KCol,1/1.1,1/1.1), KRow, KCol)  ,9))
oCol <- order(colSums(connectParamTrue$mean[,]),decreasing = TRUE)
oRow <- order(rowSums(connectParamTrue$mean[,]),decreasing = TRUE)
connectParamTrue$mean <- connectParamTrue$mean[oRow,oCol]
#-------  sizes of networks
nbNodes <- matrix(sample(10*c(6:10),M*2,replace = TRUE),M,2)
#nbNodes[1,] = c(150,200)
#nbNodes[5,] = c(30,20)

#---------  SIMULATION
mySampler <- lapply(1:M, function(m){sampleBipartiteSBM(nbNodes = nbNodes[m,],  list(row=blockPropTrue$row[m,],col=blockPropTrue$col[m,]),  connectParamTrue)})
collecNetworks <- lapply(mySampler,function(l){l$networkData})

save(collecNetworks,file='myTrialData.Rdata')
mydata <- list(collecNetworks = collecNetworks, M= M, nbNodes = nbNodes)

###########################################################################################
#------------------ VBEM 
###########################################################################################

#--------------- init CollecTau
initBipartite <- lapply(collecNetworks,estimateBipartiteSBM)
myRef <- which.max(rowSums(t(sapply(initBipartite,function(m){m$nbBlocks}))))
collecTau_init <- initCollecTau(initBipartite,ref = myRef)
KRow <- initBipartite[[myRef]]$nbBlocks[1]
KCol <- initBipartite[[myRef]]$nbBlocks[2]

print(c(KRow,KCol))

#---------------  Prior distribution

hyperparamPrior <- setHyperparamPrior(M,KRow,KCol, emissionDist, model)

#------------- variational estim 

estimOptionsVBEM <- list(maxIterVB = 1000,
                         maxIterVE = 100,
                         valStopCritVE = 10^-5,
                         valStopCritVB = 10^-5)

resEstimVBEM  <- VBEMBipartiteColSBM(collecNetworks,hyperparamPrior,collecTau_init,estimOptions = estimOptionsVBEM, emissionDist, model)

#------------------ Set ApproxPost
hyperparamApproxPost <- resEstimVBEM$hyperparamPost
hyperparamApproxPost$collecTau <- resEstimVBEM$collecTau



###########################################################################################
#------------------  SMC 
###########################################################################################
estimOptionsSMC = list()
estimOptionsSMC$paramsMCMC <- list(nbIterMCMC=2)  
estimOptionsSMC$MC <- 1000
estimOptionsSMC$ESS.rate <- 0.9
estimOptionsSMC$cESS.rate <- 0.9
estimOptionsSMC$opSave <- TRUE
estimOptionsSMC$op.parallel  <- list(os =  .Platform$OS.type, mc.cores = 1); 
estimOptionsSMC$op.print<- TRUE
estimOptionsSMC$NB.iter.max  <- Inf # Inf
estimOptionsSMC$op.SMC.classic <- FALSE
estimOptionsSMC$op.SMC.classic <- TRUE
resSMC <- SMCColBipartiteSBM(data = mydata,hyperparamPrior,hyperparamApproxPost = NULL, emissionDist, model, estimOptionsSMC)


save(mydata,resSMC,hyperparamApproxPost,file='myTrialSMCResults.Rdata')

estimOptionsSMC$op.SMC.classic <- TRUE
resSCM_classic <- SMCColBipartiteSBM(data = mydata,hyperparamPrior,hyperparamApproxPost, emissionDist, model, estimOptionsSMC)

plot(resSMC$alpha.vec,type='l')
