rm(list=ls())
setwd('~/WORK_ALL/RECHERCHE/TRAVAUX_RECHERCHE/Barbillon-Chabert/BayesianBipartiteColSBM/')

library(sbm)
library(gtools)
sapply(list.files(paste0(getwd(),'/Functions'),full.names = TRUE), source)


where_save_res <- paste0(getwd(),'/Simu/RES_SIMU2/')

emissionDist = 'bernoulli'

###############################################################
######## LOAD data
###############################################################
load(paste0(getwd(),'/Simu/RES_SIMU1/data_simu1_1.Rdata'))

collecNetworks <- mydata$collecNetworks
M <- mydata$M
nbNodes <- mydata$nbNodes


###########################################################################################
############# Load res piCol
#############################################################################################

load(paste0(getwd(),'/Simu/RES_SIMU1/res_simu1_1_SMC_VB_piCol.Rdata'))

logMargLik_piColBipartiteSBM<- sum(resSMC_VB_piCol$vec.log.ratio.Z)


###########################################################################################
############# Load res iidCol
#############################################################################################
model = 'iidColBipartiteSBM'

KRow = 4
KCol = 3

###########################################################################################
#------------------ VBEM  + SMC -VBEM 
###########################################################################################


#------------ hyperparamPrior
hyperparamPrior_iidCol <- setHyperparamPrior(M,KRow,KCol, emissionDist, model =  'iidColBipartiteSBM')


#--------------- init CollecTau
initBipartite <- lapply(collecNetworks,estimateBipartiteSBM)
myRef <- which.max(rowSums(t(sapply(initBipartite,function(m){m$nbBlocks}))))
collecTau_init <- initCollecTau(initBipartite,ref = myRef)
KRow <- initBipartite[[myRef]]$nbBlocks[1]
KCol <- initBipartite[[myRef]]$nbBlocks[2]


#------------- variational estim 

estimOptionsVBEM <- list(maxIterVB = 1000,
                         maxIterVE = 100,
                         valStopCritVE = 10^-5,
                         valStopCritVB = 10^-5)

resEstimVBEM_iidCol  <- VBEMBipartiteColSBM(collecNetworks,hyperparamPrior_iidCol,collecTau_init,estimOptions = estimOptionsVBEM, emissionDist, model  = 'iidColBipartiteSBM')

#------------------ Set ApproxPost
hyperparamApproxPost_iidCol <- resEstimVBEM_iidCol$hyperparamPost
hyperparamApproxPost_iidCol$collecTau <- resEstimVBEM_iidCol$collecTau

#------------------  SMC 

resSMC_VB_iidCol <- SMCColBipartiteSBM(data = mydata,hyperparamPrior_iidCol,hyperparamApproxPost_iidCol, emissionDist, model  = 'iidColBipartiteSBM', estimOptionsSMC)

save(hyperparamApproxPost_iidCol, hyperparamPrior_iidCol, resSMC_VB_iidCol , estimOptionsSMC,file=paste0(where_save_res,'/res_simu1_1_SMC_VB_iidCol.Rdata'))

logMargLik_iidColBipartiteSBM <- sum(resSMC_VB_iidCol$vec.log.ratio.Z)


c(logMargLik_iidColBipartiteSBM, logMargLik_piColBipartiteSBM)

#################################################################################### 
#compar post of parameters
########################################################################################


par(mfrow=c(KRow,KCol))
for (k in 1:KRow){
  for (l in 1:KCol){
    #   if ((k ==1) & (l==1)){
    plot(density(resSMC_VB_piCol$HSample_end$connectParamSample[k,l,],weights = resSMC_VB_piCol$W.end),main='alpha')
    lines(density(resSMC_VB_iidCol$HSample_end$connectParamSample[k,l,],weights = resSMC_VB_iidCol$W.end),main='alpha',col='green')
     curve(dbeta(x,hyperparamApproxPost_piCol$connectParam$alpha[k,l],hyperparamApproxPost_piCol$connectParam$beta[k,l]),col='red',add=TRUE)
    abline(v = connectParamTrue$mean[k,l])
  }
}

#------------------------------------ 
par(mfrow=c(M,KRow))
for (m in 1:M){
  for (k in 1:KRow){
      plot(density(resSMC_VB_piCol$HSample_end$blockPropSample$row[m,k,],weights= resSMC_VB_piCol$W.end),main='pi_k',xlim = c(0,1))
      lines(density(resSMC_VB_iidCol$HSample_end$blockPropSample$row[k,],weights= resSMC_VB_iidCol$W.end),col='green')
      curve(dbeta(x,hyperparamApproxPost_piCol$blockProp$row[m,k], sum(hyperparamApproxPost_piCol$blockProp$row[m,])-hyperparamApproxPost_piCol$blockProp$row[m,k]),col='red',add=TRUE)
      abline(v = blockPropTrue$row[m,k])
    }  
  }

#------------------------------------ 
par(mfrow=c(M,KCol))
for (m in 1:M){
  for (l in 1:KCol){
    plot(density(resSMC_VB_piCol$HSample_end$blockPropSample$col[m,l,],weights= resSMC_VB_piCol$W.end),main='rho_l',xlim = c(0,1))
    lines(density(resSMC_VB_iidCol$HSample_end$blockPropSample$col[l,],weights= resSMC_VB_iidCol$W.end),col='green')
    curve(dbeta(x,hyperparamApproxPost_piCol$blockProp$col[m,l], sum(hyperparamApproxPost_piCol$blockProp$col[m,])-hyperparamApproxPost_piCol$blockProp$col[m,l]),col='red',add=TRUE)
    abline(v = blockPropTrue$col[m,l])
  }  
}



  

