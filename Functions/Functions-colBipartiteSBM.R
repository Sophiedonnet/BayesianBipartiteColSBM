
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






