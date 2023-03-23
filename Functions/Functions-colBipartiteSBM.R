
#############################################################################
# Sampling in the prior distribution
#############################################################################
rParamZ <- function(MC, hyperparam, emissionDist, model ,nbNodes, collecTau = NULL){
  
  KRow  <- nrow(hyperparam$connectParam$alpha)
  KCol  <- ncol(hyperparam$connectParam$alpha)
  M <- nrow(hyperparam$blockProp$row)
  
  
  
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
    blockPropSample$row <- vapply(1:MC,function(){rdirichlet(1,hyperparam$blockProp$row)},rep(0,KRow))
    #--------- Simul of blockPropSample$col  = array(dim = c(M, KCol,MC))
    blockPropSample$col <- vapply(1:MC,function(){rdirichlet(1,hyperparam$blockProp$col)},rep(0,KCol))
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
        t(rmultinom(nRow.m,size=1,prob = pi.row.m)) %*%(1:KRow)},
        rep(0,nRow.m)
        )
      
      Z.m$col <- vapply(1:MC, function(mc){
        if(model == 'piColBipartiteSBM'){pi.col.m <- blockPropSample$col[m,,mc]}
        if(model == 'iidColBipartiteSBM'){pi.col.m <- blockPropSample$col[,mc]}
        t(rmultinom(nCol.m,size=1,prob = pi.col.m)) %*%(1:KCol)}
        ,rep(0,nCol.m)
        )
      return(Z.m)
    })
  }else{
    ZSample <- lapply(1:M,function(m){
      Z.m <- vector(mode='list',length = 2); names(Z.m) = c('row','col')
      nRow.m <- nbNodes[m,1]
      nCol.m <- nbNodes[m,2]
      Z.m$row <- t(sapply(1:nRow.m,function(i){sample(1:KRow,MC,prob = collecTau[[m]]$row[i,],replace=TRUE)}))
      Z.m$col <- t(sapply(1:nCol.m,function(j){sample(1:KCol,MC,prob = collecTau[[m]]$col[j,],replace=TRUE)}))
      return(Z.m)
    })
  }
  
  H.sample <- list(connectParamSample = connectParamSample, blockPropSample = blockPropSample,  ZSample = ZSample)
  return(H.sample)
}



#############################################################################
# COND log lik
#############################################################################

cond.log.lik = function(data, H){
  
  A = as.matrix(H$alpha)[H$Znum, H$Znum]
  A.vec = F_Sym2Vec(A)
  
  is.beta = !is.null(data$X.mat) & (length(H$beta) > 0)
  
  
  if (is.beta) {
    beta = matrix(H$beta,ncol = 1)
    res = sum(dpois(data$Y.vec, exp(A.vec + data$X.mat %*%  beta), log = T))
  }else{
    res = sum(dpois(data$Y.vec, exp(A.vec), log = T))
  }
  return(res)
}

#############################################################################
#  COND log lik for a sample H.sample
#############################################################################

likelihood <- function(data, H.sample){
  # cat('likelihood ')
  K = ncol(as.matrix(H.sample$pi)); 
  
  is.beta = !is.null(data$X.mat) & (length(H.sample$beta[1]) > 0) 
  if (is.beta) {p <- dim(data$X.array)[3]} else {p <- 0 } 
   
  if (K == 1) {M = length(H.sample$pi)} else {M <- ifelse(is.vector(H.sample$pi), 1, nrow(H.sample$pi))}
  #------------------------------------
  condloglik.sample = vapply(1:M, function(m){
     if ((K == 1) & (p == 0)) {H.m = list(Znum =  H.sample$Znum[m,], alpha = H.sample$alpha[m], beta = c())}
     if ((K == 1) & (p == 1)) {H.m = list(Znum =  H.sample$Znum[m,], alpha = H.sample$alpha[m], beta = H.sample$beta[m])}
     if ((K == 1) & (p > 1)) {H.m = list(Znum =   H.sample$Znum[m,], alpha = H.sample$alpha[m], beta = H.sample$beta[m,])}
     if ((K > 1) & (p == 0)) {H.m = list(Znum =   H.sample$Znum[m,], alpha = H.sample$alpha[m,,], beta = c())}
     if ((K > 1) & (p == 1)) {H.m = list(Znum =   H.sample$Znum[m,], alpha = H.sample$alpha[m,,],beta = H.sample$beta[m])}
     if ((K > 1) & (p > 1)) {H.m = list(Znum =   H.sample$Znum[m,], alpha = H.sample$alpha[m,,],beta = H.sample$beta[m,])}
    return(cond.log.lik(data,H.m))
    }
  ,1)
  #------------------------------------
  return(condloglik.sample) 
}




#################################################################################
# PRIOR log density
#################################################################################

logNorm.alpha = function(alpha,malpha,Salpha){
  sum(dnorm(alpha[lower.tri(alpha,diag = TRUE)], mean = malpha[lower.tri(malpha,diag = TRUE)], sd = sqrt(Salpha[lower.tri(Salpha,diag = TRUE)]), log = T))
}

#-----------------------------------------------------------------------
logPrior.alpha <- function(H.sample, HyperParms){
  
  # 
  K = length(HyperParms$e); 
  if (K == 1) {M = length(H.sample$pi)}else{M <- ifelse(is.vector(H.sample$pi), 1, nrow(H.sample$pi))}
  if ( (K > 1) & (M == 1)) {H.sample$alpha = array(H.sample$alpha,c(1,K,K))}
  
  # cat('logprior ')
  
  if (K == 1) {
    alpha.logPrior = vapply(1:M, function(m){logNorm.alpha(H.sample$alpha[m], HyperParms$malpha,HyperParms$Salpha)}, 1)
  }else{
    alpha.logPrior = vapply(1:M, function(m){logNorm.alpha(H.sample$alpha[m,,], HyperParms$malpha,HyperParms$Salpha)}, 1)
  }
  return(alpha.logPrior)
}
#-----------------------------------------------------------------------
logPrior.pi <- function(H.sample, HyperParms){
  # cat('logprior ')
  K = length(HyperParms$e)
  if (K == 1) {M = length(H.sample$pi)}else{M <- ifelse(is.vector(H.sample$pi), 1, nrow(H.sample$pi))}
  
  pi.logPrior <- ifelse(rep(K == 1,M), rep(0, M), log(ddirichlet(H.sample$pi, rep(1, M) %*% t(HyperParms$e))) )
  
  return(pi.logPrior)
}
#-----------------------------------------------------------------------
logPrior.beta <- function(H.sample, HyperParms){
  # cat('logprior ')
 p = length(HyperParms$mbeta)
 K = length(HyperParms$e)
 #--------------------------------
 if (K == 1) {M = length(H.sample$pi)} else {M <- ifelse(is.vector(H.sample$pi), 1, nrow(H.sample$pi))}
 #-------------------------------------
 if (M == 1) {H.sample$beta <- matrix(H.sample$beta,nrow = 1)}
 if (p == 1) {beta.logPrior = vapply(1:M, function(m){dmvnorm(H.sample$beta[m], mean = HyperParms$mbeta, sigma = HyperParms$Sbeta, log = T)}, 1)
 }else{beta.logPrior = vapply(1:M, function(m){dmvnorm(H.sample$beta[m, ], mean = HyperParms$mbeta, sigma = HyperParms$Sbeta, log = T)}, 1)
  }
 
 return(beta.logPrior)
}
#-----------------------------------------------------------------------
logPrior.Z <- function(H.sample, HyperParms){

  if (K == 1) {M = length(H.sample$pi)} else {M <- ifelse(is.vector(H.sample$pi), 1, nrow(H.sample$pi))}      
  Z.logPrior = vapply(1:M, function(m){
    t.m <- c(table(H.sample$Znum[m, ]))
    sum(t.m*log(H.sample$pi[m, ]))},1)
 
  return(Z.logPrior)
}
#-----------------------------------------------------------------------
logPrior <- function(H.sample, HyperParms){
  is.beta <- (length(HyperParms.prior$mbeta)>0)
  
  logPrior.sample = logPrior.pi(H.sample, HyperParms) +
    logPrior.alpha(H.sample, HyperParms) + 
    logPrior.Z(H.sample, HyperParms)
  if(is.beta){logPrior.sample = logPrior.sample + logPrior.beta(H.sample, HyperParms)}
  return(logPrior.sample)
}


#############################################################################
# ABOUT de Approx posterior distribution
#############################################################################

#--------------- SAMPLING --------------------------- 
# rApproxPost <- function(M, HyperParms.ApproxPost){
#    
#    K = length(HyperParms.ApproxPost$e)
#    n = nrow(HyperParms.ApproxPost$tau)
#    
#    pi.sample = rdirichlet(M, HyperParms.ApproxPost$e)
#    alpha.sample = array(dim = c(M, K, K))
#    invisible(lapply(1:K, function(k){
#       lapply(k:K, function(l){
#          alpha.sample[, k, l] <<- HyperParms.ApproxPost$malpha[k, l] + 
#             sqrt(HyperParms.ApproxPost$Salpha[k, l]) * rnorm(M)
#          if (l > k) {alpha.sample[, l, k] <<- alpha.sample[, k, l]}
#       })}))
#    beta.sample = as.matrix(rmvnorm(M, mean = HyperParms.ApproxPost$mbeta, sigma = HyperParms.ApproxPost$Sbeta))
#    
#    # Sampling Z
#    Z.sample = array(dim = c(M, n, K)); Znum.sample = matrix(0, M, n)
#    invisible(lapply(1:n, function(i){
#       Z.sample[, i, ] <<- t(rmultinom(M, 1, HyperParms.ApproxPost$tau[i, ]))
#       if (K > 1) {
#          Znum.sample[, i] <<- Z.sample[, i, ] %*% (1:K)
#       }else{
#          Znum.sample[, i] <<- Z.sample[, i, ]
#       }
#    }))
# 
#    H.sample = list(pi = pi.sample, alpha = alpha.sample, beta = beta.sample, Znum = Znum.sample)
#    return(H.sample)
# }

rApproxPost <- function(M, HyperParms.ApproxPost){
  H.sample <- rPrior(M,HyperParms.ApproxPost)
  return(H.sample)
}
 
#---------------------------  Log approx post --------------------- 


# logApproxPost <- function(H.sample, HyperParms.ApproxPost){
#  
#   pi.logApproxPost = logPrior.pi(H.sample,HyperParms.ApproxPost)
#   alpha.logApproxPost = logPrior.alpha(H.sample,HyperParms.ApproxPost)
#   beta.logApproxPost = logPrior.beta(H.sample,HyperParms.ApproxPost)
#   
#   Z.logApproxPost = vapply(1:M, function(m){
#     Z.lAPtmp = 0
#     invisible(lapply(1:n, function(i){
#       Z.lAPtmp <<- Z.lAPtmp + log(HyperParms.ApproxPost$tau[i, H.sample$Znum[m, i]])
#     }))
#     return(Z.lAPtmp)}, 1)
#   #cat('*********************** \n')
#   logApproxPost.sample = pi.logApproxPost + alpha.logApproxPost + beta.logApproxPost + Z.logApproxPost
#   return(logApproxPost.sample)
# }
# 

logApproxPost <- function(M, HyperParms.ApproxPost){
  logApproxPost.sample <- logPrior(M,HyperParms.ApproxPost)
  return(logApproxPost.sample)
}



#############################################################################
# MCMC kernel
#############################################################################


MCMC.Kernel <- function(data, H, alpha.t, HyperParms.ApproxPost, HyperParms.prior, Parms.MCMC, op.save=FALSE, op.print=FALSE, op.SMC.classic = FALSE){
   # rho = alpha.t
   
  ### 
  # HyperParms : Param?tres de la loi a priori
  # HyperParms.ApproxPost : param?tres de la loi a posteriori approch?e
  # alpha.t : param?tres de pond?ration post, post approch?e
  # B :  Number of it?rations
  
  ### variables latentes et param?tres
 
  p <- length(HyperParms.ApproxPost$mbeta)
  is.beta <- (!is.null(data$X.array)) & (p > 0) & (length(H$beta) > 0); 
  if (is.beta == FALSE) {p <-  0; HyperParms.ApproxPost$mbeta = c()}
  
  
  
  #H <- H.0
  if (is.beta) {Sigma <- Parms.MCMC$Sigma[-1,-1]}
  
  n <- length(H$Znum); 
  K <- length(c(H$pi)); 
  if (K == 1) { Parms.MCMC$op.echan$pi <- Parms.MCMC$op.echan$Znum <- 0; H$alpha = matrix(H$alpha,1,1)}

  
  H$Z <- vapply(1:K, function(k){as.numeric(H$Znum == k)}, rep(1,n))
  
  seqZnum = matrix(0,Parms.MCMC$B,n)
  seqpi = matrix(0,Parms.MCMC$B,K)
  seqalpha = array(0,c(Parms.MCMC$B,K,K))
  if (is.beta) {seqbeta = matrix(0,Parms.MCMC$B,p)}
  
  loglik <- cond.log.lik(data,H)
  lprior.alpha <- logPrior.alpha(H, HyperParms.prior)
  lapproxPost.alpha <- ifelse(op.SMC.classic, lprior.alpha , logPrior.alpha(H,HyperParms.ApproxPost))
  if (is.beta) { 
    lprior.beta <- logPrior.beta(H, HyperParms.prior)
    lapproxPost.beta <- ifelse(op.SMC.classic, lprior.beta, logPrior.beta(H, HyperParms.ApproxPost))
  }
  
   
  
  accep.alpha = matrix(0,K,K) 
  if (is.beta) { accep.beta = 0; } 
  for (b in 1:Parms.MCMC$B) {
     # cat(b, '')
    if ((op.print > 0) & (b %% op.print == 0)) {
      if (is.beta) {print(c('MCMC Kernel : iteration',b,H$beta))
      }else{print(c('MCMC Kernel : iteration',b))}
    }
    
  
    
    ############# simulation of pi  | Z 
    
 
    
    if (Parms.MCMC$op.echan$pi == 1) {
      H$pi = rdirichlet(1, alpha.t*(colSums(H$Z) + HyperParms.prior$e) + (1 - alpha.t) * HyperParms.ApproxPost$e)
    } 
    
    
    ############# simulation of alpha (using adjusted Langevin metropolis Hastings) ????? 
    if (Parms.MCMC$op.echan$alpha == 1) {
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



