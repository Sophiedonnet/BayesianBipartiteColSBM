
#############################################################################
# Sampling in the prior distribution
#############################################################################
rParamZ <- function(MC, hyperparam, emissionDist, model ,nbNodes){
  
  KRow  <- nrow(hyperparam$connectParam$alpha)
  KCol  <- ncol(hyperparam$connectParam$alpha)
  M <- nrow(nbNodes)
  
  collecTau <- hyperparam$collecTau
  
  
  #--------- Simul of connectParamSample = array(dim = c(KRow, KCol,MC))
  if(emissionDist == 'bernoulli'){
    connectParamSample <- vapply(1:MC,function(mc){matrix(rbeta(KRow*KCol,hyperparam$connectParam$alpha,hyperparam$connectParam$beta),KRow,KCol)},matrix(0,KRow,KCol))
  }
  if(emissionDist == 'poisson'){
    connectParamSample <- vapply(1:MC,function(mc){matrix(rgamma(KRow*KCol,hyperparam$connectParam$alpha,hyperparam$connectParam$beta),KRow,KCol)},matrix(0,KRow,KCol))
  }
  
  #--------- Simul of blockPropSample
  blockPropSample <- vector(mode = "list", length = 2)
  names(blockPropSample) <- c('row','col')
 
  if(model == 'iidColBipartiteSBM'){
    #--------- Simul of blockPropSample$row  = array(dim = c(KRow,MC))
    blockPropSample$row <- vapply(1:MC,function(mc){rdirichlet(1,hyperparam$blockProp$row)},rep(0,KRow))
    #--------- Simul of blockPropSample$col  = array(dim = c(M, KCol,MC))
    blockPropSample$col <- vapply(1:MC,function(mc){rdirichlet(1,hyperparam$blockProp$col)},rep(0,KCol))
  }
  
  if(model == 'piColBipartiteSBM'){
    #--------- Simul of blockPropSample$row  = array(dim = c(KRow,MC))
    blockPropSample$row <- vapply(1:MC,function(mc){t(sapply(1:M,function(m){rdirichlet(1,hyperparam$blockProp$row[m,])}))},matrix(0,M,KRow))
    #--------- Simul of blockPropSample$col  = array(dim = c(M, KCol,MC))
    blockPropSample$col <- vapply(1:MC,function(mc){t(sapply(1:M,function(m){rdirichlet(1,hyperparam$blockProp$col[m,])}))},matrix(0,M,KCol))
  }
  
  
  # Sampling ZRow and Zcol
  
  if(is.null(collecTau)){
   
    
    ZSample <- lapply(1:M,function(m){
      Z.m <- vector(mode='list',length = 2); names(Z.m) = c('row','col')
      nRow.m <- nbNodes[m,1]
      nCol.m <- nbNodes[m,2]
      Z.m$row <- vapply(1:MC, function(mc){
        if(model == 'piColBipartiteSBM'){pi.row.m <- blockPropSample$row[m,,mc]}
        if(model == 'iidColBipartiteSBM'){pi.row.m <- blockPropSample$row[,mc]}
        t(rmultinom(nRow.m,size=1,prob = pi.row.m))},
        matrix(0,nRow.m,KRow)
        )
      
      Z.m$col <- vapply(1:MC, function(mc){
        if(model == 'piColBipartiteSBM'){pi.col.m <- blockPropSample$col[m,,mc]}
        if(model == 'iidColBipartiteSBM'){pi.col.m <- blockPropSample$col[,mc]}
        t(rmultinom(nCol.m,size=1,prob = pi.col.m))}
        ,matrix(0,nCol.m,KCol)
        )
      return(Z.m)
    })
  }else{
    ZSample <- lapply(1:M,function(m){
      Z.m <- vector(mode='list',length = 2); names(Z.m) = c('row','col')
      nRow.m <- nbNodes[m,1]
      nCol.m <- nbNodes[m,2]
      Z.m$row <- array(0,c(nRow.m,KRow,MC))
      Z.m$col <- array(0,c(nCol.m,KCol,MC))
      for (i in 1:nRow.m){Z.m$row[i,,] = rmultinom(MC,size=1,prob = collecTau[[m]]$row[i,])}
      for (j in 1:nCol.m){Z.m$col[j,,] = rmultinom(MC,size=1,prob = collecTau[[m]]$col[j,])}
      return(Z.m)
    })
  }
  
  HSample <- list(connectParamSample = connectParamSample, blockPropSample = blockPropSample,  ZSample = ZSample)
  return(HSample)
}



#############################################################################
# COND log lik
#############################################################################
condLogLik = function(data, H.mc,emissionDist){
  M <- data$M; 
  collecNetworks <- data$collecNetworks
  if(emissionDist == 'bernoulli'){
      v <- sapply(1:M,function(m){
        par.m <- H.mc$Z[[m]]$row%*% H.mc$connectParam %*% t(H.mc$Z[[m]]$col)
        res.m  <- sum(dbinom(collecNetworks[[m]],1,par.m,log = TRUE))
        }
        )
  }
  if(emissionDist == 'poisson'){
    v <- sapply(1:M,function(m){
      par.m <- H.mc$Z[[m]]$row%*% H.mc$connectParam %*% t(H.mc$Z[[m]]$col)
      res.m  <- sum(dpois(collecNetworks[[m]],par.m,log = TRUE))
    }
    )
  }
  res <- sum(v)
  return(res)
}

#############################################################################
#  COND log lik for a sample H.sample
#############################################################################

likelihood <- function(data, HSample,emissionDist){
  # cat('likelihood ')
 
  #------------------------------------
  condloglik.sample = vapply(1:MC, function(mc){
    H.mc <- list(connectParam = HSample$connectParamSample[,,mc])
    H.mc$Z <- lapply(HSample$ZSample,function(Zm){list(row = Zm$row[,,mc],col = Zm$col[,,mc])})
    return(condLogLik(data,H.mc,emissionDist))
    }
  ,1)
  #------------------------------------
  return(condloglik.sample) 
}




#################################################################################
# PRIOR log density
#################################################################################


#-----------------------------------------------------------------------
logDistConnectParam <- function(HSample,MC, hyperparam,emissionDist){
  # 
  if(emissionDist=='bernoulli'){
    res <- vapply(1:MC,function(mc){
      sum(dbeta(HSample$connectParamSample[,,mc],hyperparam$connectParam$alpha,hyperparam$connectParam$beta,log = TRUE))
    },1)
  }
  if(emissionDist=='poisson'){
    res <- vapply(1:MC,function(mc){
      sum(dgamma(HSample$connectParamSample[,,mc],hyperparam$connectParam$alpha,hyperparam$connectParam$beta,log = TRUE))
    },1)
  }
  
  return(res)
}

#-----------------------------------------------------------------------
logDirichletBlockProp <- function(HSample,MC, hyperparam,model){
 
  
  if (model=='piColBipartiteSBM'){
    res.row <- vapply(1:MC,function(mc){sum(log(ddirichlet(HSample$blockPropSample$row[,,mc],hyperparam$blockProp$row)))},0)
    res.col <- vapply(1:MC,function(mc){sum(log(ddirichlet(HSample$blockPropSample$col[,,mc],hyperparam$blockProp$col)))},0)
  }
  if(model =='iidColBipartiteSBM'){
    res.row <- vapply(1:MC,function(mc){log(ddirichlet(HSample$blockPropSample$row[,mc],hyperparam$blockProp$row))},0)
    res.col <- vapply(1:MC,function(mc){log(ddirichlet(HSample$blockPropSample$col[,mc],hyperparam$blockProp$col))},0)
    
  }
  return(res.row + res.col)
}


#-----------------------------------------------------------------------
logMultinomZ <- function(HSample, M, MC, hyperparam,model){

 res.m  = rep(0,MC)
 for (m in 1:M){
     Zsample.m.row <-HSample$ZSample[[m]]$row
     if(is.null(hyperparam$collectTau)){
       t.m.row <- apply(Zsample.m.row,c(2,3),sum)
       if(model == 'piColBipartiteSBM'){pi.row.m <- HSample$blockPropSample$row[m,, ]}
       if(model == 'iidColBipartiteSBM'){pi.row.m <- HSample$blockPropSample$row}
       u.m <-  apply(t.m.row*log(pi.row.m),c(2),sum)
     }else{
       u.m <- vapply(1:MC,function(mc){sum(Zsample.m.row[,,mc]*log(hyperparam$collecTau[[m]]$row))},1)
      }

     
     Zsample.m.col <-HSample$ZSample[[m]]$col
     if(is.null(hyperparam$collectTau)){
      t.m.col <- apply(Zsample.m.col,c(2,3),sum)
      if(model == 'piColBipartiteSBM'){pi.col.m <- HSample$blockPropSample$col[m,, ]}
      if(model == 'iidColBipartiteSBM'){pi.col.m <- HSample$blockPropSample$col}
      v.m <-  apply(t.m.col*log(pi.col.m),c(2),sum)
     }else{
       v.m <- vapply(1:MC,function(mc){sum(Zsample.m.col[,,mc]*log(hyperparam$collecTau[[m]]$col))},1)
     }
     
     res.m = res.m + u.m+v.m
 }
 return(res.m)
}

#-----------------------------------------------------------------------
logJointParamZ <- function(HSample,M, MC, hyperparam,emissionDist,model){
 
  
  a <- logDistConnectParam(HSample,MC, hyperparam,emissionDist)
  b <- logDirichletBlockProp(HSample,MC, hyperparam,model)
  c <- logMultinomZ(HSample, M, MC, hyperparam,model) 
  return(a+b+c)
}






#############################################################################
# MCMC kernel
#############################################################################


MCMC.Kernel <- function(data, H.mc, alpha.t, hyperparamApproxPost, hyperparamPrior, emissionDist, model,Parms.MCMC, op.save=FALSE, op.print=FALSE, op.SMC.classic = FALSE){
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
  KCol <- nrow(hyperparamPrior$connectParam$alpha)
  
 
  seqZnum = matrix(0,Parms.MCMC$B,n)
  seqpi = matrix(0,Parms.MCMC$B,K)
  seqalpha = array(0,c(Parms.MCMC$B,K,K))
  if (is.beta) {seqbeta = matrix(0,Parms.MCMC$B,p)}
  
  loglik <- cond.log.lik(data,H.mc,emissionDist)
  lprior.alpha <- logPrior.alpha(H, HyperParms.prior)
  lapproxPost.alpha <- ifelse(op.SMC.classic, lprior.alpha , logPrior.alpha(H,HyperParms.ApproxPost))
  if (is.beta) { 
    lprior.beta <- logPrior.beta(H, HyperParms.prior)
    lapproxPost.beta <- ifelse(op.SMC.classic, lprior.beta, logPrior.beta(H, HyperParms.ApproxPost))
  }
  
 
    
  
    
    ############# simulation of pi  | Z 
    
 
    
    if (Parms.MCMC$op.echan$pi == 1) {
      H.mc$blockProp = rdirichlet(1, alpha.t*(colSums(H$Z) + HyperParms.prior$e) + (1 - alpha.t) * HyperParms.ApproxPost$e)
    } 
    
    
    ############# simulation of connectParam (using adjusted Langevin metropolis Hastings) ????? 
    if (Parms.MCMC$op.echan$connectParam == 1) {
      
      S  <- lapply(1:M,function(m){t(H.mc$Z[[m]]$row) %*% matrix(1,nbNodes[m,1],nbNodes[m,2]) %*% H.mc$Z[[m]]$col})
      SY  <- lapply(1:M,function(m){t(H.mc$Z[[m]]$row) %*% collecNetworks[[m]] %*% H.mc$Z[[m]]$col})
      alpha_sim <- alpha.t*(hyperparamPrior$connectParam$alpha +  Reduce(`+`, SY)) + (1-alpha.t)*hyperparamApproxPost$connectParam$alpha
      beta_sim <- alpha.t*(hyperparamPrior$connectParam$beta +  Reduce(`+`, SY) - (emissionDist == 'bernoulli') * Reduce(`+`, SY)) + (1-alpha.t)*hyperparamApproxPost$connectParam$beta
      
      if(emissionDist == 'bernoulli'){H.mc$connectParam <- matrix(rbeta(KRow*KCol,alpha_sim,beta_sim),KRow,KCol)}
      if(emissionDist == 'poisson'){H.mc$connectParam <- matrix(rgamma(KRow*KCol,alpha_sim,beta_sim),KRow,KCol)}
    }
  ############# simulation of blockParam
  if (Parms.MCMC$op.echan$blockProp == 1) {
    
    S.row <- t(sapply(1:M, function(m){colSums(H.mc$Z[[m]]$row)}))
    S.col <- t(sapply(1:M, function(m){colSums(H.mc$Z[[m]]$col)}))
    
    if(model=='iidColBipartiteSBM'){
      e_row <- alpha.t*(colSums(S.row) + hyperparamPrior$blockProp$row) + (1 - alpha.t) * hyperparamPrior$blockProp$row
      H.mc$blockProp$row <- rdirichlet(1,e_row)
      e_col <- alpha.t*(colSums(S.col) + hyperparamPrior$blockProp$col) + (1 - alpha.t) * hyperparamPrior$blockProp$col
      H.mc$blockProp$col <- rdirichlet(1,e_col)
    }
    if(model=='piColBipartiteSBM'){
      e_row <- alpha.t*(S.row + hyperparamPrior$blockProp$row) + (1 - alpha.t) * hyperparamPrior$blockProp$row
      H.mc$blockProp$row <- t(sapply(1:M,function(m){rdirichlet(1,e_row[m,])}))
      e_col <- alpha.t*(S.col + hyperparamPrior$blockProp$col) + (1 - alpha.t) * hyperparamPrior$blockProp$col
      H.mc$blockProp$col <- t(sapply(1:M,function(m){rdirichlet(1,e_col[m,])}))
    }
  }
  ############# simulation of Z 
  if (Parms.MCMC$op.echan$Z == 1) {
    
    for (m in 1:M){
      
      
    }
    
  }
  
  
    
      })         
      
    }
    
  }  
    e_dim <- 
      H.mc$connectParam = rdirichlet
      H.mc$Z 
      H.mc$blockProp
  
      for (k in 1:K) {
        for (l in 1:k) {
      
          H.c <- H; 
          H.c$alpha <- as.matrix(H$alpha)
          H.c$alpha[k,l] <-  H$alpha[k,l] + Parms.MCMC$tau.alpha * sample(c(0.1,1,10),1)*rnorm(1)
          if (l != k) {H.c$alpha[l,k] = H.c$alpha[k,l]}
          
          loglik.c <- cond.log.lik(data, H.c)
          lprior.alpha.c <- logPrior.alpha(H.c, HyperParms.prior)
          lapproxPost.alpha.c <- ifelse(op.SMC.classic, lprior.alpha.c , logPrior.alpha(H.c,HyperParms.ApproxPost))
          
          #q.alpha.alpha.c <- q.alpha.c.alpha <- 0;
          T1.c <-  alpha.t * (loglik.c + lprior.alpha.c) +  (1 - alpha.t) *  lapproxPost.alpha.c
          T1   <-  alpha.t * (loglik   + lprior.alpha  ) +  (1 - alpha.t) *  lapproxPost.alpha
          
          log.prob.accept =  T1.c - T1
        
          if (log(runif(1)) < log.prob.accept) {
            H <- H.c 
            loglik <- loglik.c; 
            lprior.alpha <-  lprior.alpha.c
            lapproxPost.alpha <- lapproxPost.alpha.c  
            accep.alpha[k,l] <- accep.alpha[k,l] + 1
          }
        }
      }
    }
   
    ##############"" simulation of beta 
    if ((is.beta) & (Parms.MCMC$op.echan$beta == 1)) {
      
        H.c <- H; 
        H.c$beta <-   c(H$beta) + rmvnorm(1,rep(0,p),Parms.MCMC$tau * sample(c(1/10,1,10),1) * as.matrix(Sigma))
        
        loglik.c <- cond.log.lik(data, H.c)
        lprior.beta.c <- logPrior.beta(H.c,HyperParms.prior) 
        lapproxPost.beta.c <- ifelse(op.SMC.classic, lprior.beta.c , logPrior.beta(H.c,HyperParms.ApproxPost))
       
        q.beta.beta.c <- q.beta.c.beta <- 0;
        T1.c <-  alpha.t * (loglik.c + lprior.beta.c) + (1 - alpha.t) * lapproxPost.beta.c 
        T1   <-  alpha.t * (loglik  + lprior.beta) + (1 - alpha.t) * (lapproxPost.beta)
        
          
        log.prob.accept =  T1.c - T1 - (q.beta.c.beta - q.beta.beta.c)
        if (log(runif(1)) < log.prob.accept) {
          H <- H.c; 
          loglik <- loglik.c; 
          lprior.beta <-  lprior.beta.c
          lapproxPost.beta <- lapproxPost.beta.c  
          accep.beta = accep.beta + 1
      }
    }
    # cat('fin beta i= \n') 
    
    
    ##############"" simulation of Z 
    if (Parms.MCMC$op.echan$Znum == 1) {
      order_echan = sample(1:n,n,replace = FALSE)
      for (i in order_echan) {
       
        if (p == 1) {
          vec_tau = data$Y.mat[i, ] %*% as.matrix(H$alpha)[H$Znum, ] - colSums(exp(as.vector(data$X.array[i, -i, ] * H$beta) %o% rep(1, K) + as.matrix(H$alpha)[H$Znum[-i], ]))
        }
        if (p > 1) {
          vec_tau = data$Y.mat[i, ] %*% as.matrix(H$alpha)[H$Znum, ] - colSums(exp(as.vector(data$X.array[i, -i, ] %*% matrix(H$beta,ncol = 1)) %o% rep(1, K)  + as.matrix(H$alpha)[H$Znum[-i], ]))
        }
        
        if (p == 0) { 
          vec_tau = data$Y.mat[i, ] %*% as.matrix(H$alpha)[H$Znum, ] - colSums(exp(as.matrix(H$alpha)[H$Znum[-i], ]))
        }
        
        vec_tau <- vec_tau*0.5  #### (car matrice symétrique) 
        
       
        tau_i <- ifelse(rep(op.SMC.classic,K), alpha.t * vec_tau + log(H$pi), alpha.t * (vec_tau + log(H$pi)) + (1 - alpha.t)*log(HyperParms.ApproxPost$tau[i,]))
        tautmp = tau_i; if (sum(is.na(tau_i)) > 0) {print(tautmp); save(tautmp, file = 'TauTmp.Rdata')}
        tau_i = exp(tau_i - max(tau_i, na.rm = T))
        H$Z[i, ] = rmultinom(1, 1, tau_i)
        H$Znum <- H$Z %*% (1:K); 
    }
    loglik <- cond.log.lik(data,H)
    }
    
   if (op.save == TRUE) {
      seqZnum[b,] <- H$Znum
      seqpi[b,] <- H$pi
      if (is.beta) { seqbeta[b,] <- H$beta }
      seqalpha[b,,] <- H$alpha
    }
  }
  if (op.save == FALSE) { res <- H}
  if (op.save == TRUE) {
    res <- list(accep.alpha = accep.alpha,seqpi = seqpi,seqalpha = seqalpha,seqZnum = seqZnum)
    if (is.beta) { res$seqbeta = seqbeta }
  }
  return(res)
  
  
  # cat('Fin MCM.Kernel')
}



