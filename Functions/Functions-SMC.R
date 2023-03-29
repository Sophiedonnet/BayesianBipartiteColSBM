library(plyr)

#-----------------------------------------------------------------
# Computes ESS from unnormalized log-weights
#-----------------------------------------------------------------


ComputecESS <- function(W.t.minus.1,w.increm.t) {
  cESS <- sum(W.t.minus.1  *  w.increm.t)^2 / sum(W.t.minus.1  *  w.increm.t^2)
  return(cESS)
}

#-----------------------------------------------------------------
# Dichotomie to find next alpha
#-----------------------------------------------------------------

FindAlpha.cESS <- function(W.t.minus.1, log.rho, cESS.rate, alpha.t.minus.1, tol=1e-4, op.save=FALSE) {
  
  cESS_vec <- c()
  alpha_vec <- c()
  
  threshold <- cESS.rate
  alpha.new = 1;
  cESS <- ComputecESS(W.t.minus.1,w.increm.t = exp((alpha.new - alpha.t.minus.1) * (log.rho - mean(log.rho))))
  if ( (cESS < threshold)|is.na(cESS)) {
    alpha.left <- alpha.t.minus.1;
    cESS.left <- ComputecESS(W.t.minus.1,w.increm.t = exp((alpha.left - alpha.t.minus.1)  *  (log.rho - mean(log.rho))))
    alpha.right <- 1;
    cESS.right <-  ComputecESS(W.t.minus.1,w.increm.t = exp((alpha.right - alpha.t.minus.1) * (log.rho - mean(log.rho))))
    
    while (is.na(cESS.right) & (alpha.right > alpha.t.minus.1)) {
      alpha.right <- alpha.right  *  0.99;
      cESS.right <-  ComputecESS(W.t.minus.1,w.increm.t = exp((alpha.right - alpha.t.minus.1)  *  (log.rho - mean(log.rho))))
    }
    
    alpha.new <- (alpha.left + alpha.right)/2; cESS <- ComputecESS(W.t.minus.1,w.increm.t=exp((alpha.new-alpha.t.minus.1) * (log.rho-mean(log.rho))))
    if (op.save == TRUE) {cESS_vec <- c(cESS_vec,cESS); alpha_vec <- c(alpha_vec,alpha.new)}
    diff <- 2 * tol; niter = 0
    while ((diff > tol) & (alpha.new > alpha.t.minus.1 + 1e-4) & (niter < 1e3)) {
      niter = niter + 1
      if (cESS > threshold) {alpha.left <- alpha.new; cESS.left <- cESS}else{alpha.right <- alpha.new; cESS.right <- cESS}
      alpha.new <- (alpha.left + alpha.right)/2; cESS <- ComputecESS(W.t.minus.1,w.increm.t=exp((alpha.new-alpha.t.minus.1) * (log.rho-mean(log.rho))))
      diff <- abs(cESS - threshold)
      # cat(alpha.new, cESS, diff, '\n')
      if (op.save == TRUE) {cESS_vec <- c(cESS_vec,cESS);alpha_vec<- c(alpha_vec,alpha.new)}
    }
    
  }
  
  output <-  list(alpha.new=alpha.new,cESS_vec=cESS_vec,alpha_vec=alpha_vec)
  return(output)
}

#-----------------------------------------------------------------
# Objects to sample
#-----------------------------------------------------------------

Func.sampling <- function(obj,Resample) {
  if (is.vector(obj)) {obj_2 <-  obj[Resample]}
  if (is.array(obj)) {
    if (length(dim(obj)) == 2) {obj_2 <- array(obj[Resample,], dim=dim(obj))}
    if (length(dim(obj)) == 3) {obj_2 <- array(obj[Resample,,], dim=dim(obj))}
    if (length(dim(obj)) == 4) {obj_2 <- array(obj[Resample,,,], dim=dim(obj))}
    if (length(dim(obj)) == 5) {obj_2 <- array(obj[Resample,,,,], dim=dim(obj))}
    if (length(dim(obj))>5) {print('erreur of array size')}
  }
  return(obj_2)
}


#-----------------------------------------------------------------
# Tools functions 
#-----------------------------------------------------------------

Func.search <- function(HSample,mc) {
  Nb_obj <- length(HSample)
  H.mc <-  list();
  for (o in 1:Nb_obj) {
    obj=HSample[[o]]; ##### echantillon de taille MC pour l'objet o. donc de taille   MCx(? ? ?)
    if (is.vector(obj)) {H.mc[[o]] = obj[mc]}
    if (is.array(obj)) {
      if (length(dim(obj)) == 2) {H.mc[[o]] = obj[mc,]}
      if (length(dim(obj)) == 3) {H.mc[[o]] = obj[mc,,]}
      if (length(dim(obj)) == 4) {H.mc[[o]] = obj[mc,,,]}
      if (length(dim(obj)) == 5) {H.mc[[o]] = obj[mc,,,,]}
      if (length(dim(obj))>5) {print('erreur of array size')}
    }
  }
  names(H.mc) = names(HSample)
  return(H.mc)
}

#------------------------------------
#   Computation of the log marginal likelihood
#-----------------------------------


estim.loglikmarg.U = function(RES_SMC) {
  U = RES_SMC$U
  niter = length(U)
  l <- sum(diff(RES_SMC$alpha.vec) * (U[1:(niter-1)] + U[2:niter]))/2
  return(l)
}

#-----------------------------------------------------------------
# SMC ############################################  SMC version 2
#-----------------------------------------------------------------

SMCColBipartiteSBM<- function(collecNetworks,hyperparamPrior,hyperparamApproxPost, emissionDist, model, estimOptionsSMC){
  
  
  Parms.MCMC.in.SMC <- estimOptionsSMC$Parms.MCMC.in.SMC
  MC <- estimOptionsSMC$MC
  ESS.rate <- estimOptionsSMC$ESS.rate
  cESS.rate <- estimOptionsSMC$cESS.rate
  op.save <- estimOptionsSMC$op.save #  FALSE
  op.parallel <- estimOptionsSMC$op.parallel 
  op.print <- estimOptionsSMC$op.print #TRUE
  NB.iter.max <- estimOptionsSMC$NB.iter.max # Inf
  op.SMC.classic <- estimOptionsSMC$typeSMC == 'classic'
  
  

  
  # ----------------- parameters of the SMC
  os  <- op.parallel$os
  #-------------------------------------------------------------------------------------
  #------------------   ITERATION 0  : simulation with the first init (VBEM or Prior)
  #-------------------------------------------------------------------------------------
  
  if (op.SMC.classic == FALSE) {
    HSample <- rParamZ(MC, hyperparamApproxPost, emissionDist, model,nbNodes); W.t <- rep(1/MC,MC)
  }else{
    HSample <- rParamZ(MC, hyperparamPrior, emissionDist, model,nbNodes); W.t <- rep(1/MC,MC)
    }
  RES <- list(HSample.0 = HSample,W.0 = W.t);
  alpha.t <- 0; 
  t <- 1;  
  alpha.vec <- c(alpha.t);
  vec.log.ratio.Z <- c(0); 
  KL <- c(0, NA);
  
  Vec.resampling = c(TRUE)

  if (op.save == TRUE) {
    RES.t <- list(); 
    RES.t[[t]] <- list(HSample = HSample,W.t = W.t,alpha.vec = alpha.vec,alpha.t = alpha.t)
  }
  
  logPrior.sample <- logPrior(HSample, hyperparamPrior)
  if (op.SMC.classic) {logApproxPost.sample <- logPrior.sample}else{logApproxPost.sample <- logApproxPost(HSample, hyperparamApproxPost)}
  condloglik.sample <- likelihood(collecNetworks, HSample)
  log.rho <- logPrior.sample  +  condloglik.sample - logApproxPost.sample
  
  r.Inf <- which(log.rho == -Inf)
  if (length(r.Inf) > 0) {log.rho[r.Inf] = min(log.rho[-r.Inf])/2}
  U <- c(sum(W.t * log.rho))
  Nb_obj <- length(HSample)
  
  #-------------------------------------------------------------------------------------
  #------------------   ITERATIONS suivantes
  #-------------------------------------------------------------------------------------
  while ((alpha.t < 1)  & (t < NB.iter.max)) {
    t <- t  + 1 #(first value  = 1)
    W.t.minus.1 <-  W.t
    
    #--- calculation of the new alpha
    res.find.alpha.new <- FindAlpha.cESS(W.t.minus.1,log.rho, cESS.rate,alpha.t, tol = 1e-4,op.save = FALSE)
    alpha.new <- res.find.alpha.new$alpha.new
    delta.alpha.t <- alpha.new - alpha.t;
    if (delta.alpha.t < 0) {browser()}
    alpha.t <- alpha.new
    if (alpha.t > 1) {alpha.t <- 1};
    alpha.vec <- c(alpha.vec,alpha.t);
    
    #---- Compute the weigths of the new particules. ----------------------------------------------------
    Anum = delta.alpha.t * mean(log.rho) ## nomalizing constant to prevent numerical problems
    w.increm.t <- exp(delta.alpha.t  *  log.rho - Anum)
    Ww.t <-  W.t.minus.1 * w.increm.t
    W.t <- Ww.t/sum(Ww.t)
    
    #---  compute the log ratio of the marginal likelihoods---------------------------------------
    #--- ratio.Z <- sum(W.t.minus.1 * w.increm.t) # Sum (W^(n-1)  * wtilde_n ). Cf calculs et Del Moral
    log.ratio.Z <- log(sum(W.t.minus.1 * w.increm.t)) + Anum
    vec.log.ratio.Z <- c(vec.log.ratio.Z,log.ratio.Z)
    
    #------ Resampling if necesary  -----------------------------------------------------------------
    ESS.t <- 1/sum((W.t)^2)
    if (is.nan(ESS.t)) {print("Error in ESS.t. Break"); t = t - 1; save.image(file = "fail.Rdata"); break}
    resampling.step <- (ESS.t < M * ESS.rate)
    Vec.resampling <- c(Vec.resampling,resampling.step)
    
    if (resampling.step == TRUE) {
      Resample <- sample(1:M, replace = T, prob = W.t)
      W.t <- rep(1/M,M)
      invisible(lapply(1:Nb_obj,function(o) {HSample[[o]] <<- Func.sampling(HSample[[o]],Resample)}))
    }
    
    # w.cumul.t <-  w.cumul.t*w.increm.t
    # if ((resampling.step) | (alpha.t == 1)) {
    #   ratio.Z.2 <- sum(w.cumul.t)
    #   vec.ratio.Z.2 <- c(vec.ratio.Z.2,ratio.Z.2)
    #   w.cumul.t <- W.t
    # }
    # 
    # 
    # 
    
    
    #--------  Propagate the particles with pi_t as invariant kernel -------------
    
    H.newsample <- HSample;
    f_MCMC <- function(m) {
      if ((m %% 500 == 0) & (model == 'SBMreg')) {cat(m, '')}
      # cat('\n m=', m, 'b=')
      H.m <- Func.search(HSample,m)
      GMMB.m <- MCMC.Kernel(data, H = H.m, alpha.t, hyperparam.ApproxPost, hyperparam.prior, op.save = FALSE, op.print = Parms.MCMC.in.SMC$B + 1, Parms.MCMC = Parms.MCMC.in.SMC,op.SMC.classic)
      return(GMMB.m)}
    
    
    if (os == "unix") {
      OUTPUT_MCMC <- mclapply(1:M, f_MCMC, mc.preschedule = TRUE, mc.cores = op.parallel$mc.cores)
    }else{
      OUTPUT_MCMC <- lapply(1:M,f_MCMC)
    }
    # cat('Fin MCMC \n')
    
    OUTPUT_MCMC <- do.call(cbind, OUTPUT_MCMC)
    for (o in 1:Nb_obj) {
      obj <- H.newsample[[o]];
      type.obj <- length(dim(obj)) - 1
      if (type.obj == 1) {obj <- matrix(unlist(OUTPUT_MCMC[o,]),nrow = M,byrow = TRUE)} ### recuperation pour des vecteurs
      if (type.obj == 2) {for (m in 1:M) {obj[m,,] <- OUTPUT_MCMC[o,m][[1]]}} ### recuperation pour des matrices
      if (type.obj == 3) {for (m in 1:M) {H.newsample[[o]][m,,,] <-  OUTPUT_MCMC[o,m][[1]]}}
      H.newsample[[o]] = obj
    }
    # cat(dim(HSample$pi), colMeans(HSample$pi), '/')
    HSample <- H.newsample
    # cat(dim(HSample$pi), colMeans(HSample$pi), '\n')
    # cat('Fin Output_MCMC \n')
    
    if (op.save == TRUE) {
      RES.t[[t]] <- list(HSample = HSample,W.t = W.t,alpha.vec = alpha.vec,alpha.t = alpha.t);
      RES.t[[t]]$w.increm.t <- w.increm.t;
      save(RES.t,file = 'encours.Rdata');
    }
    
    #------ Computation avec the log.rho to compute the w.increment.t  and U for marginal likelihood --------------------------------
    
    logPrior.sample <- logPrior(HSample, hyperparam.prior)
    if (op.SMC.classic == TRUE) {logApproxPost.sample <- logPrior.sample}else{logApproxPost.sample <- logApproxPost(HSample, hyperparam.ApproxPost)}
    condloglik.sample <- likelihood(data, HSample)
    log.rho <- logPrior.sample  +  condloglik.sample - logApproxPost.sample
    log.rho[which(log.rho == -Inf)] = -1e4
    log.rho[which(logApproxPost.sample == -Inf)] = -1e4
    log.rho[which(log.rho == Inf)] = -1e4
    
    U.t = sum(W.t * log.rho)
    U <- c(U,U.t)
    
    # cat('Fin log.rho \n')
    
    
    #----- Computation of the various KL   +  ussful constants
    KLalpha <- c(alpha.t * (sum(W.t * log.rho)) - sum(vec.log.ratio.Z), -(1 - alpha.t) * sum(W.t * log.rho) - sum(vec.log.ratio.Z))
    KL <- rbind(KL, KLalpha)
    
    
    # Moments
    if (model == 'LogReg') {
      MeanSample[[t]] <- colSums(matrix(W.t,ncol = M) %*% HSample$beta)
      VarSample[[t]] <- cov.wt(matrix(HSample$beta,nrow = M),wt = W.t)$cov
    }
    if (model == 'SBM') {
      theta.sample = cbind(HSample$pi[, -K], F_Array2Mat(HSample$gamma))
      MeanSample[[t]] = t(W.t)  %*%  theta.sample
      VarSample[[t]] = cov.wt(theta.sample, wt = W.t)$cov
    }
    
    # cat('Fin moments \n')
    
    if (op.print == TRUE) {print(c(t,alpha.t,ESS.t,sum(vec.log.ratio.Z)))}
    #if (op.print == TRUE) {print(c(t,alpha,ESS.t,U.t))}
  }### ----------------------------------------------------------- end  of algorithm
  
  
  
  log.Zalpha <- cumsum(vec.log.ratio.Z)
  KL[, 1] <- KL[, 1]
  log.Z <- log.Zalpha[length(log.Zalpha)]
  KL[, 2] <- KL[, 2]  +  log.Z
  if (is.element(model, c('LogReg', 'SBM'))) {
    KLgauss = matrix(0, t, 2)
    for (i in 1:t) {
      KLgauss[i, 1] = ComputeKLgauss(MeanSample[[i]], MeanSample[[1]], VarSample[[i]], VarSample[[1]])
      KLgauss[i, 2] = ComputeKLgauss(MeanSample[[i]], MeanSample[[t]], VarSample[[i]], VarSample[[t]])
      
    }
  }
  
  ##########
  if (op.save == TRUE) {RES$RES.t <- RES.t}
  RES$HSample_end <- HSample
  RES$W.end <- W.t
  RES$alpha.vec <- alpha.vec
  RES$vec.log.ratio.Z <- vec.log.ratio.Z;
  #RES$vec.ratio.Z.2 <- vec.ratio.Z.2;
  RES$KL <- KL
  RES$Vec.resampling <- Vec.resampling
  RES$U <- U
  
  
  if (is.element(model, c('LogReg', 'SBM'))) {
    RES$KLgauss <- KLgauss
    RES$MeanSample <- MeanSample
    RES$VarSample <- VarSample
  }
  return(RES)
}


mutualInformationZ  = function(Zsample){
  if (!is.matrix(Zsample)) {stop('Sample must be a Matrix with M rows if M is the size of the sample')}
  M <- nrow(Zsample)
  n <- ncol(Zsample)
  H.ind <- sum( sapply(1:n,function(i){freq_i = table(Zsample[,i])/M; return(sum(freq_i*log(freq_i)))}))
  U <- plyr::count(Zsample)$freq/M
  H.conj <- sum(U * log(U))
  return(H.conj  - H.ind)
}




