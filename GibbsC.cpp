#include <RcppArmadilloExtensions/sample.h>

using namespace Rcpp;
using namespace arma;


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
List GS_iterations(int burnin, int N, arma::vec Dstart,arma::mat Lstart,arma::mat Sstart,
                   double rho_start, mat X, int p, int n, double q,double tau, 
                   double alpha1, double alpha2, double lambda2){
  arma::vec Dcurrent = Dstart;
  arma::vec Destimator = arma::zeros<arma::vec>(p);
  arma::sp_mat Lcurrent(Lstart);
  arma::mat Lcounts = arma::zeros<arma::mat>(p,p);
  arma::mat Lestimator = arma::zeros<arma::mat>(p,p);
  arma::mat S =Sstart;
  double rho_current =rho_start;
  double rho_estimater = 0;
  double trace1;
  double proposal;
  double trace2;
  mat res_prop;
  mat Sprop;
  sp_mat matD;
  mat Omega;
  for(int iter=0;iter< (burnin+N);++iter){
    arma::mat zeroLs = arma::zeros<arma::mat>(p,p);
    for(int j=0;j<p-1;++j){
      int nzero=0;
      for(int i=1;i<p-j;++i){
        arma::rowvec Sj= S(i+j,arma::span(j,p-1));
        arma::vec Lj(Lcurrent(arma::span(j,p-1),j));
        double b = arma::as_scalar(sum(Sj*Lj)-Lj[i]*S((i+j),(i+j)));
        double d = S((i+j),(i+j))+1/pow(tau,2);
        
        double weightNorm1 = log(1-q)+pow(b,2)/(2*Dcurrent[j]*d)-log(tau)-log(d)/2;
        arma::vec weightNorm;
        
        weightNorm <<log(q)<<weightNorm1<<arma::endr;
        
        weightNorm = exp(weightNorm-max(weightNorm));
        /*if(i==2 and j==0)
        weightNorm.print();
        */
        //print(weightNorm)
        arma::vec choices;
        choices <<0<<1<<arma::endr;
        int choice = Rcpp::RcppArmadillo::sample(choices,1,FALSE,weightNorm)[0];
        if(choice){
          double temp = (randn<double>())*sqrt(Dcurrent[j]/d)-b/d;
          Lcurrent(i+j,j)= temp;
          
        }
        else{
          Lcurrent(i+j,j) = 0.0;
          zeroLs(i+j,j)+=1;
          nzero+=1;
        }
        //Lcurrent.print();
    }
      arma::vec Lj(Lcurrent(arma::span(j,p-1),j));
      double priorBeta = arma::as_scalar(Lj.t()*S(arma::span(j,p-1),arma::span(j,p-1))*Lj/2)+alpha2;
      double postBeta = arma::as_scalar(sum(pow(Lj.subvec(1,Lj.size()-1),2))/(2*pow(tau,2)))+priorBeta;
      double postAlpha = alpha1+(p-j-nzero+n-1)/2.0;
      double temp = randg<double>(arma::distr_param(postAlpha,1.0));
      Dcurrent[j] = postBeta/temp;
      //Rcout << postAlpha << " "<< postBeta <<"" <<Dcurrent[j] << endl;
  }
    Dcurrent[p] = (S(p-1,p-1)/2+alpha2)/randg<double>(arma::distr_param(alpha1+n/2,1.0));
    
    matD = sp_mat(diagmat(1/Dcurrent));
    proposal = (randn<double>())*0.1+rho_current;
    res_prop = X(span(1,n),span::all)-proposal*X(span(0,n-1),span::all);
    Sprop =res_prop.t()*res_prop;
    Omega =sp_mat(Lcurrent)*matD*sp_mat(Lcurrent.t());
    trace1 = trace(Omega*S);
    trace2 = trace(Omega*Sprop);
    if(log(randu<double>()) <(trace1-trace2)/2){
      S = Sprop;
      rho_current = proposal;
    }
    if(iter>=burnin){
      Lcounts = Lcounts+zeroLs;
      Lestimator = Lestimator+Lcurrent;
      Destimator = Destimator+Dcurrent/N;
      rho_estimater = rho_estimater+rho_current/N;
    }  
}
  return List::create(_["D"]=Destimator,_["L"]=Lestimator,_["Lcounts"]=Lcounts,_["rho"]=rho_estimater);
}

/*** R
sampling_DL <- function(res,X,rho,q,tau,alpha1,alpha2,lambda2,burnin,N,epsilon,L=NULL,D=NULL){
  p <- ncol(res)
  n <- nrow(res)
  S <- t(res)%*%res
  if(is.null(L)){
    ch <- chol(solve(S+epsilon*diag(p)))
    dd <- diag(ch)
    Lcurrent<- t(ch/dd)
#print(Lcurrent)
    Dcurrent <- 1/dd^2
  }else{
    Lcurrent=L
    Dcurrent=D
  }
  results <- GS_iterations(burnin,N,Dcurrent,Lcurrent,S,rho,X,p,n,q,tau,alpha1,alpha2,lambda2)
    return(results)
}

*/
