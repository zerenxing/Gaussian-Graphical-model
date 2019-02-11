generate_data <- function(n,p,initial_L,shared_L){
  X <- matrix(0,n,p)
  X[1,] <- rmvn.sparse(1,rep(0, p),initial_L)
  for(i in 2:n){
    X[i,] <- rou*X[i-1,] + rmvn.sparse(1,rep(0, p),shared_L)
  }
  return(X)
}

generateL_AR3<- function(p,lags=rep(1,3)){
  diags <- list(rep(1,p),rep(lags[1],p-1),rep(lags[2],p-2),rep(lags[3],p-3))
  true.L <- bandSparse(p,k=c(0,-1,-2,-3),diagonals = diags)
  true.L <- as.matrix(true.L)
  diag(true.L) <- 1
  return(true.L)
}


generateL_random <- function(p,scale=1){
  L <- matrix(1,p,p)
  L[upper.tri(L)] <- 0
  diag(L) <- 0
  for(i in 1:(p-1))
  {
    L[(i+1):p,i] <- L[(i+1):p,i]*rbinom(p-i,1,min(3/(p-i),1))
  }
  L <- L*scale
  diag(L) <- 1
  return(L)
}


Evaluation <- function(Adj1, Adj2){
  #Adj1 <- as.matrix(get.adjacency(graph1))
  #Adj2 <- as.matrix(get.adjacency(graph2))
  true.index <- which(Adj1==1)
  print(length(true.index))
  #false.index <- setdiff(which(upper.tri(Adj1)),true.index)
  false.index <- which(lower.tri(Adj1)&(Adj1==0))
  print(length(false.index))
  positive.index <- which(Adj2==1)
  print(length(positive.index))
  # negative.index <- setdiff(which(upper.tri(Adj2)),positive.index)
  negative.index <- which(lower.tri(Adj2)&(Adj2==0))
  print(length(negative.index))
  
  TP <- length(intersect(true.index,positive.index))
  FP <- length(intersect(false.index,positive.index))
  FN <- length(intersect(true.index,negative.index))
  TN <- length(intersect(false.index,negative.index))
  print(c(TP,FP,FN,TN))
  
  MCC.denom <- sqrt(TP+FP)*sqrt(TP+FN)*sqrt(TN+FP)*sqrt(TN+FN)
  if(MCC.denom==0) MCC.denom <- 1
  MCC <- (TP*TN-FP*FN)/MCC.denom
  if((TN+FP)==0) MCC <- 1
  
  Precision <- TP/(TP+FP)
  if((TP+FP)==0) Precision <- 1
  Recall <- TP/(TP+FN)
  if((TP+FN)==0) Recall <- 1
  Sensitivity <- Recall
  Specific <- TN/(TN+FP)
  if((TN+FP)==0) Specific <- 1
  
  #graph1.graph <- igraph.to.graphNEL(graph1)
  #graph2.graph <- igraph.to.graphNEL(graph2)
  #SHD <- shd(graph1.graph,graph2.graph)
  return(list(Precision=Precision,Recall=Recall,Sensitivity=Sensitivity,Specific=Specific,MCC=MCC,TP=TP,FP=FP,TN=TN,FN=FN))
}