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

###########################################################################################
############# simulation 
#############################################################################################
M = 3
KRow = 4
KCol = 3


#-----  block proportions simul piColSBM avec classes vides
if(model=='piColBipartiteSBM'){
  blockPropTrue <- list()
  blockPropTrue$row <-  rdirichlet(M,rep(1/1.3,KRow))  #### emptying some blocks in certain netwokrs
  blockPropTrue$col <- rdirichlet(M,rep(1/1.3,KCol))   #### emptying some blocks in certain netwokros
  blockPropTrue$row[1,] <- rep(1/KRow,KRow)
  blockPropTrue$col[1,] <- rep(1/KCol,KCol)
  print(blockPropTrue)
}

#-----  connectivity matrix
connectParamTrue <- list(mean = round( matrix(rbeta(KRow*KCol,1/1.1,1/1.1), KRow, KCol)  ,9))
oCol <- order(colSums(connectParamTrue$mean[,]),decreasing = TRUE)
oRow <- order(rowSums(connectParamTrue$mean[,]),decreasing = TRUE)
connectParamTrue$mean <- connectParamTrue$mean[oRow,oCol]
#-------  sizes of networks
nbNodes <- matrix(sample(10*c(8:15),M*2,replace = TRUE),M,2)
nbNodes  <- 2*nbNodes
#nbNodes[1,] = c(150,200)
#nbNodes[5,] = c(30,20)

#---------  SIMULATION
mySampler <- lapply(1:M, function(m){sampleBipartiteSBM(nbNodes = nbNodes[m,],  list(row=blockPropTrue$row[m,],col=blockPropTrue$col[m,]),  connectParamTrue)})
collecNetworks <- lapply(mySampler,function(l){l$networkData})

mydata <- list(collecNetworks = collecNetworks, M= M, nbNodes = nbNodes)

save(connectParamTrue, blockPropTrue, model, emissionDist,mydata,file=paste0(where_save_res,'/data_simu1_2.Rdata'))

##############################################################
#---------- load data
##############################################################

load(file=paste0(where_save_res,'/data_simu1_2.Rdata'))
collecNetworks <- mydata$collecNetworks
M <- mydata$M
nbNodes <- mydata$nbNodes



###################################################### 
#---------------  Prior distribution
######################################################

KRow = 4
KCol = 3 
hyperparamPrior_piCol <- setHyperparamPrior(M,KRow,KCol, emissionDist, model='piColBipartiteSBM')




###########################################################################################
#------------------ VBEM  + SMC -VBEM 
###########################################################################################

#--------------- init CollecTau
initBipartite <- lapply(collecNetworks,estimateBipartiteSBM)
myRef <- which.max(rowSums(t(sapply(initBipartite,function(m){m$nbBlocks}))))
collecTau_init <- initCollecTau(initBipartite,ref = myRef)
KRow <- initBipartite[[myRef]]$nbBlocks[1]
KCol <- initBipartite[[myRef]]$nbBlocks[2]
print(c(KRow,KCol))

#------------- variational estim 

estimOptionsVBEM <- list(maxIterVB = 1000,
                         maxIterVE = 100,
                         valStopCritVE = 10^-5,
                         valStopCritVB = 10^-5)

resEstimVBEM_piCol  <- VBEMBipartiteColSBM(collecNetworks,hyperparamPrior_piCol,collecTau_init,estimOptions = estimOptionsVBEM, emissionDist, model='piColBipartiteSBM')
LogLikMarg_VB_piCol <- sum(computeLogLikMarg_VB(collecNetworks,resEstimVBEM_piCol$collecTau,hyperparamPrior_piCol,emissionDist,model))

#------------------ Set ApproxPost
hyperparamApproxPost_piCol <- resEstimVBEM_piCol$hyperparamPost
hyperparamApproxPost_piCol$collecTau <- resEstimVBEM_piCol$collecTau

#------------------  SMC 

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
resSMC_VB_piCol <- SMCColBipartiteSBM(data = mydata,hyperparamPrior_piCol,hyperparamApproxPost_piCol, emissionDist, model = 'piColBipartiteSBM', estimOptionsSMC)

save(model,hyperparamApproxPost_piCol, hyperparamPrior_piCol, resSMC_VB_piCol, estimOptionsSMC ,file=paste0(where_save_res,'/res_simu1_2_SMC_VB_piCol.Rdata'))



####################################################
#---------------------- SMC classic
####################################################

estimOptionsSMC$op.SMC.classic <- TRUE


resSMC_Classic <- SMCColBipartiteSBM(data = mydata,hyperparamPrior,hyperparamApproxPost = NULL, emissionDist, model, estimOptionsSMC)

save(estimOptionsSMC,hyperparamPrior,resSMC_Classic,file=paste0(where_save_res,'/res_simu1_1_SMC_Classic.Rdata'))


