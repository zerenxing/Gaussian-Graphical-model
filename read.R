
p_seq <- 50
n_seq <- c(200,400,800,1600,3200)
rho_seq <- c(0,0.5)
iters <- 10
N <- 1e05
true_model <- "random"
threshs <- c(0.25,0.5,0.75)
compile_evaluation<-function(p_seq,n_seq,rho_seq,iters,threshs,true_model){
  for(i in 1:length(p_seq)){
    for(rho in rho_seq){
      evaluation <- list()
      p=p_seq[i]
      for(j in 1:length(n_seq)){
        for(iter in 1:iters ){
          n=n_seq[j]
          evaluation_=read.csv(file=paste0("results_",true_model,"/evaulation_rho=",rho,"_p=",p,"_n=",n,"_iter=",iter,".csv"),row.names = 1)
          for(k in 1:length(threshs)){
            if(j==1 & iter==1){
              evaluation[[k]] <- matrix(0,nrow=nrow(evaluation_),ncol=length(n_seq))
              colnames(evaluation[[k]])<-paste0("n=",n_seq)
              rownames(evaluation[[k]])<-paste0("mean_",rownames(evaluation_))
            }
            evaluation[[k]][,j] <-evaluation[[k]][,j]+evaluation_[,k]/iters
          }
        }
      }
      for(k in 1:length(threshs)){
        write.csv(evaluation[[k]],file=paste0(true_model,"_evaluation_summary_p=",p,"_rho=",rho,"_thresh=",threshs[k],".csv"))
      }
    }
  }
}
library(reshape2)
library(ggplot2)
library(scales)  
compile_errors <- function(p_seq,n_seq,rho_seq,threshs,N,true_model){
  for(i in 1:length(p_seq)){
    for(rho in rho_seq){
      p=p_seq[i]
      for(j in 1:length(n_seq)){
        n=n_seq[j]
        D_Rho_error=read.csv(file=paste0("results_",true_model,"/D_Rho_errors_rho=",rho,"_p=",p,"_n=",n,"_N=",N,".csv"),row.names = 1)
        L_Omega_error = read.csv(file=paste0("results_",true_model,"/L_Omega_errors_rho=",rho,"_p=",p,"_n=",n,"_N=",N,".csv"),row.names = 1)
        if(j==1){
          D_Rho_errors<-NULL
          L_Omega_errors <- NULL
        }
        D_Rho_error$n <- n
        L_Omega_error$n <- n
        D_Rho_errors<-rbind(D_Rho_errors,D_Rho_error)
        L_Omega_errors <- rbind(L_Omega_errors,L_Omega_error)
      }
      L_Omega_errors_long<- reshape(L_Omega_errors, direction="long", 
                                    varying=list(c("errorL_mean","errorOmega_mean"), c("errorL_sd","errorOmega_sd")), 
                                    v.names=c("mean","sd"),timevar="parameter")
      L_Omega_errors_long$parameter = factor(L_Omega_errors_long$parameter,levels=c(1,2),labels=c("L","Omega"))
      L_Omega_errors_long$threshs<- factor(L_Omega_errors_long$threshs)
      write.csv(L_Omega_errors_long,paste0(true_model,"_L_Omega_errors_rho=",rho,"_p=",p,"_N=",N,".csv"))
      D_Rho_errors_long<- reshape(D_Rho_errors, direction="long", 
                                  varying=list(c("errorD_mean","errorRho_mean"), c("errorD_sd","errorRho_sd")), 
                                  v.names=c("mean","sd"),timevar="parameter")
      D_Rho_errors_long$parameter = factor(D_Rho_errors_long$parameter,levels=c(1,2),labels=c("D","Rho"))
      write.csv(D_Rho_errors_long,paste0(true_model,"_D_Rho_errors_rho=",rho,"_p=",p,"_N=",N,".csv"))
      p_L_Omega <- ggplot(L_Omega_errors_long,aes(x=as.factor(n),y=mean,colour=threshs,shape=parameter,
                                                  group=interaction(threshs,parameter)))+
        geom_point()+geom_line()+
        labs(title =paste0("Estimation error of ",true_model," L/Omega when p=",p), x = "n", y = "error")
      ggsave(paste0(true_model,"_L_Omega_errors_p=",p,"_rho=",rho,"_N=",N,".jpg"),p_L_Omega,width = 5,height=5)
      p_D_Rho <- ggplot(D_Rho_errors_long,aes(x=as.factor(n),y=mean,shape=parameter,group=parameter))+
        geom_point(color="red")+geom_line(color="red")+
        labs(title =paste0("Estimation error of ",true_model," D/Rho when p=",p), x = "n", y = "error")
      ggsave(paste0(true_model,"_D_Rho_errors_p=",p,"_rho=",rho,"_N=",N,".jpg"),p_D_Rho,width = 5,height=5)
    }
  }
}
compile_evaluation(p_seq,n_seq,rho_seq,iters,threshs,"AR3")
compile_evaluation(p_seq,n_seq,rho_seq,iters,threshs,"random")

compile_errors(p_seq,n_seq,rho_seq,threshs,N,"AR3")
compile_errors(p_seq,n_seq,rho_seq,threshs,N,"random")
