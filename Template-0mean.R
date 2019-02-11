library(glmnet)
library(Matrix)
library(mvtnorm)
library(lassoshooting)
library(methods)
library(sparseMVN)
library(Rcpp)
library(pROC)

sourceCpp("GibbsC.cpp")
source("utils.R")

### Function to get diagonal submatrix
## U: orignal matrix
## ind: index of elements to get
set.seed(2017)
p <- PP
n <- NN
print(p)
rou <- ZZ
rho_est<-FALSE
q<- 1-6/(p-1)
tau <- 0.5
alpha1 <- 3
alpha2 <- 2
lambda2 <- 0.25
burnin<-20000
N<-100000
true_model <- TRUEL
threshs <- c(0.25,0.5,0.75)
generateL <- get(paste0("generateL_",true_model))
true.L <-generateL(p)

true.Adj <- true.L
true.Adj[true.Adj!=0]<-1
diag(true.Adj)<-0
d <- max(colSums(true.Adj)) + 1
param <- seq(0.1,0.5,by=0.05)
Omega <- true.L%*%t(true.L)
Omega.d <- as(Omega,"dgCMatrix")
L.d <- Cholesky(Omega.d)
xOmega <- (1-rou^2)*Omega
xOmega.d <- as(xOmega,"dgCMatrix")
xL.d <- Cholesky(xOmega.d)
error= NULL
for(iter in 1:10){
  set.seed(iter+1000)
  X <- generate_data(n,p,xL.d,L.d)
  if(rho_est){
    rou <- mean(apply(X,2, function(x) sum(x[1:(n-1)]*x[2:n])/sum(x[1:(n-1)]^2)))
  }
  res = X[2:n,]-rou*X[1:(n-1),]
  result<-sampling_DL(res,X,rou,q,tau,alpha1,alpha2,lambda2,burnin,N,0.5,true.L,rep(1,p))
  preds<-result$Lcounts[lower.tri(result$Lcounts)]/N
  labels <- true.Adj[lower.tri(true.Adj)]
  resultRoc <-roc(labels,preds,direction=">",auc=T)
  dir.create(file.path(paste0("results_",true_model)),showWarnings = FALSE)
  q <- round(q,4)
  pdf(file=paste0("results_",true_model,"/ROC_rho=",rou,"_p=",p,"_n=",n,"_iter=",iter,".pdf"))
  plot(resultRoc,print.auc=TRUE)
  dev.off()
  write.table(result$Lcounts,file=paste0("results_",true_model,"/Lcountsrho=",rou,"_p=",p,"_n=",n,"_iter=",iter,".txt"))
  graphs <- lapply(threshs,function(x){
    graph = ifelse(result$Lcounts>x*N,0,1)
    graph[upper.tri(graph)] <- 0
    diag(graph) <-0
    return(graph)
  })
  Lestimators <- lapply(graphs,function(G){
    Lhat = result$L
    Lhat[G<1] =0
    Lhat[G==1] = Lhat[G==1]/(N-result$Lcounts[G==1])
    diag(Lhat) <- 1
    return(Lhat)
  })
  Omega_est <- lapply(Lestimators,function(L) L%*%diag(as.vector(result$D))%*%t(L))
  
  resultThresh.est <- sapply(graphs, Evaluation, Adj1=true.Adj)
  colnames(resultThresh.est) <- threshs
  
  resultThresh.est<-cbind(resultThresh.est,resultRoc$auc)
  write.csv(resultThresh.est,file=paste0("results_",true_model,"/evaulation_rho=",rou,"_p=",p,"_n=",n,"_iter=",iter,".csv"))
  error_rho = result$rho-rou
  error_L = sapply(Lestimators, function(L) sqrt(sum((true.L-L)^2))/sqrt(sum(true.L^2)))
  error_Omega = sapply(Omega_est, function(O) sqrt(sum((O-Omega)^2))/sqrt(sum(Omega^2)))
  error_D = sqrt(sum((as.vector(result$D)-rep(1,p))^2)/p)
  dt= data.frame(error_rho,error_D,error_L,error_Omega,threshs)
  error= rbind(error,dt)
}
library(dplyr)

D_Rho_error<-error %>% filter(threshs==0.25) %>%summarise(errorD_mean=mean(error_D),errorD_sd=sd(error_D),
                                                         errorRho_mean=mean(error_rho),errorRho_sd=sd(error_rho))
write.csv(D_Rho_error,file=paste0("results_",true_model,"/D_Rho_errors_rho=",rou,"_p=",p,"_n=",n,"_N=",N,".csv"))
L_Omega_error<-error %>% group_by(threshs) %>% summarise(errorL_mean=mean(error_L),errorL_sd=sd(error_L),
                                                        errorOmega_mean=mean(error_Omega),errorOmega_sd=sd(error_Omega))
write.csv(L_Omega_error,file=paste0("results_",true_model,"/L_Omega_errors_rho=",rou,"_p=",p,"_n=",n,"_N=",N,".csv"))