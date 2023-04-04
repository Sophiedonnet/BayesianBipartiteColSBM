


#############################################################################
# MCMC kernel
#############################################################################


MCMC.Kernel <- function(data, H.mc.init, alpha.t, hyperparApproxPost, hyperparamPrior, emissionDist, model,paramsMCMC, opSave=FALSE){
  # rho = alpha.t
  
  ### 
  # HyperParms : Param?tres de la loi a priori
  # HyperParms.ApproxPost : param?tres de la loi a posteriori approch?e
  # alpha.t : param?tres de pond?ration post, post approch?e
  # B :  Number of it?rations
  
  ### variables latentes et param?tres
  
  M <- data$M
  collecNetworks <- data$collecNetworks
  nbNodes <- data$nbNodes
  KRow <- nrow(hyperparamPrior$connectParam$alpha)
  KCol <- ncol(hyperparamPrior$connectParam$alpha)
  H.mc <- H.mc.init
  
  for (iterMCMC in 1:paramsMCMC$B){
  
    ############# simulation of connectParam
    if (paramsMCMC$op.echan$connectParam == 1) {
      
      S  <- lapply(1:M,function(m){t(H.mc$Z[[m]]$row) %*% matrix(1,nbNodes[m,1],nbNodes[m,2]) %*% H.mc$Z[[m]]$col})
      SY  <- lapply(1:M,function(m){t(H.mc$Z[[m]]$row) %*% collecNetworks[[m]] %*% H.mc$Z[[m]]$col})
      alpha_sim <- alpha.t*(hyperparamPrior$connectParam$alpha +  Reduce(`+`, SY)) + (1-alpha.t)*hyperparamApproxPost$connectParam$alpha
      beta_sim <- alpha.t*(hyperparamPrior$connectParam$beta +  Reduce(`+`, SY) - (emissionDist == 'bernoulli') * Reduce(`+`, SY)) + (1-alpha.t)*hyperparamApproxPost$connectParam$beta
      if(emissionDist == 'bernoulli'){H.mc$connectParam <- matrix(rbeta(KRow*KCol,alpha_sim,beta_sim),KRow,KCol)}
      if(emissionDist == 'poisson'){H.mc$connectParam <- matrix(rgamma(KRow*KCol,alpha_sim,beta_sim),KRow,KCol)}
    }
  
    ############# simulation of blockParam
    if (paramsMCMC$op.echan$blockProp == 1) {
      
      S.row <- t(sapply(1:M, function(m){colSums(H.mc$Z[[m]]$row)}))
      S.col <- t(sapply(1:M, function(m){colSums(H.mc$Z[[m]]$col)}))
      
      if(model=='iidColBipartiteSBM'){
        e_row <- alpha.t*(colSums(S.row) + hyperparamPrior$blockProp$row) + (1 - alpha.t) * hyperparamApproxPost$blockProp$row
        H.mc$blockProp$row <- c(rdirichlet(1,e_row))
        e_col <- alpha.t*(colSums(S.col) + hyperparamPrior$blockProp$col) + (1 - alpha.t) * hyperparamApproxPost$blockProp$col
        H.mc$blockProp$col <- c(rdirichlet(1,e_col))
      }
      if(model=='piColBipartiteSBM'){
        e_row <- alpha.t*(S.row + hyperparamPrior$blockProp$row) + (1 - alpha.t) * hyperparamPrior$blockProp$row
        H.mc$blockProp$row <- t(sapply(1:M,function(m){rdirichlet(1,e_row[m,])}))
        e_col <- alpha.t*(S.col + hyperparamPrior$blockProp$col) + (1 - alpha.t) * hyperparamPrior$blockProp$col
        H.mc$blockProp$col <- t(sapply(1:M,function(m){rdirichlet(1,e_col[m,])}))
      }
    }
    ############# simulation of Z 
    if (paramsMCMC$op.echan$Znum == 1) {
      orderNetworks <-  sample(1:M,M,replace = FALSE)
      for (m in orderNetworks){ # for all networks
        
        if(model=='iidColBipartiteSBM'){
          piRow.m <- H.mc$blockProp$row
          piCol.m <-H.mc$blockProp$col
        }
        if(model=='piColBipartiteSBM'){
          piRow.m <-H.mc$blockProp$row[m,]
          piCol.m <-H.mc$blockProp$col[m,]
        }
        
        ############ Zrow
        logProbZRow <- collecNetworks[[m]]%*%H.mc$Z[[m]]$col %*%t(log(H.mc$connectParam))    # for bernoulli and poisson
        if(emissionDist == 'poisson'){
          logProbZRow <- logProbZRow  - matrix(1,nbNodes[m,1],nbNodes[m,2]) %*%  H.mc$Z[[m]]$col %*% t(H.mc$connectParam)
        }
        logProbZRow <- alpha.t * logProbZRow + alpha.t * matrix(log(piRow.m),nrow=nbNodes[m,1],ncol=KRow,byrow = TRUE)  
        logProbZRow  <- logProbZRow  + (1-alpha.t)*log(hyperparamApproxPost$collecTau[[m]]$row )
        ProbZRow <- fromBtoTau(logProbZRow, eps = 10^-10)
        for (i in 1:nbNodes[m,1]){H.mc$Z[[m]]$row[i,] = rmultinom(1,size=1,prob = ProbZRow[i,])}
        
        ############ Zcol 
        logProbZCol <-  t(collecNetworks[[m]])%*%  H.mc$Z[[m]]$row  %*%log(H.mc$connectParam) # for bernoulli and poisson
        if(emissionDist == 'poisson'){
          logProbZCol <- logProbZCol  - matrix(1,nbNodes[m,2],nbNodes[m,1]) %*%  H.mc$Z[[m]]$row   %*%H.mc$connectParam
          }
        logProbZCol <- alpha.t * logProbZCol + alpha.t * matrix(log(piCol.m),nrow=nbNodes[m,2],ncol=KCol,byrow = TRUE)  
        logProbZCol  <- logProbZCol  + (1-alpha.t)*log(hyperparamApproxPost$collecTau[[m]]$col )
        ProbZCol <- fromBtoTau(logProbZCol, eps = 10^-10)
        for (i in 1:nbNodes[m,2]){H.mc$Z[[m]]$col[i,] = rmultinom(1,size=1,prob = ProbZCol[i,])}
      }
    }
  }
  return(H.mc)
}
      
# cat('Fin MCM.Kernel')
