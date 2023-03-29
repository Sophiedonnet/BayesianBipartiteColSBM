checkSample = function(HSample,MC,M,KRow,KCol,hyperparam,emissionDist,model){
  
  blockPropSample <- HSample$blockPropSample
  ZSample <- HSample$ZSample
  connectParamSample <- HSample$connectParamSample  
  collecTau <- hyperparam$collecTau
  #---------------  connectParamSample
  print('--------------------dim of connectParamSample--------------------')
  print(sum(dim(connectParamSample)==c(KRow,KCol,MC))==3)
  expected_mean <- hyperparam$connectParam$alpha / (hyperparam$connectParam$beta + (emissionDist == 'bernoulli')*hyperparam$connectParam$alpha)
  print(cbind(c(expected_mean),c(apply(connectParamSample,c(1,2),mean))))
  
  #--------------- blockPropSample
  print('--------------------dim of blockPropSample--------------------')
  if(model == 'iidColBipartiteSBM'){
    print(dim(blockPropSample$row)==c(KRow,MC))
    print(dim(blockPropSample$col)==c(KCol,MC))
    
    print(cbind(hyperparam$blockProp$row/sum(hyperparam$blockProp$row),apply(blockPropSample$row,1,mean)))
    print(cbind(hyperparam$blockProp$col/sum(hyperparam$blockProp$col),apply(blockPropSample$col,1,mean)))
    
    
  }
  if(model == 'piColBipartiteSBM'){
    print(dim(blockPropSample$row)==c(M,KRow,MC))
    print(cbind(c(normByRow(hyperparam$blockProp$row)),c(apply(blockPropSample$row,c(1,2),mean))))
    
    print(dim(blockPropSample$col)==c(M,KCol,MC))
    print(cbind(c(normByRow(hyperparam$blockProp$col)),c(apply(blockPropSample$col,c(1,2),mean))))
    
  }
  
  
  #---------------   Z 
  print('--------------------dim of Z-----------------')
  dimZrow <- sapply(ZSample,function(Z.m){dim(Z.m$row)})
  dimZcol <- sapply(ZSample,function(Z.m){dim(Z.m$col)})
  
  print(sum(dimZrow[1,]==nbNodes[,1])==M)
  print(sum(dimZcol[1,]==nbNodes[,2])==M)
  
  print(sum(dimZrow[2,]==KRow)==M)
  print(sum(dimZcol[2,]==KCol)==M)
  
  print(sum(dimZrow[3,]==MC)==M)
  print(sum(dimZcol[3,]==MC)==M)
  
   
  ######################### 
  mrand <- sample(1:M,1)
  um <- apply(ZSample[[mrand]]$row,c(1,2),mean)
  i <- sample(1:nbNodes[mrand,1],1)
  print(um[i,])
  if(is.null(collecTau)){
    if(model == 'piColBipartiteSBM'){print(normByRow(hyperparam$blockProp$row)[mrand,])}
    if(model == 'iidColBipartiteSBM'){print(hyperparam$blockProp$row/sum(hyperparam$blockProp$row))}
  }
  
  if(!is.null(collecTau)){
    print(collecTau[[mrand]]$row[i,])
  }
  
  
  
  
  
}
